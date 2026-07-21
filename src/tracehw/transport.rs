//! Transport abstraction: the byte pipe between the controller and the
//! instrument.
//!
//! [`Transport`] is what the run controller ([`crate::tracehw::instrument`]) drives. The
//! only implementation today is [`PckReplay`], the **simulator** that replays a
//! recorded `.pck` session — it lets the whole control/acquisition stack run
//! with no hardware. A live serial/USB transport (Phase P1/M6) will implement
//! the same trait, so nothing above it changes.

use std::io::{ErrorKind, Read, Write};

use anyhow::{Context, Result};

use crate::tracehw::pck::{parse_front, FrontParse, Record, Session};
use crate::tracehw::protocol::Command;

#[cfg(feature = "serial")]
pub mod serial;

/// A bidirectional link to the instrument: send commands, poll inbound records.
pub trait Transport {
    /// Send a command to the instrument.
    fn send(&mut self, command: &Command) -> Result<()>;

    /// Poll the next inbound record. `Ok(None)` means the stream is exhausted
    /// (for a live link, that a record is not yet available).
    fn poll(&mut self) -> Result<Option<Record>>;
}

/// Replays a recorded `.pck` session as if it were a live instrument.
///
/// Records are handed out in file order, one per [`poll`](Transport::poll).
/// Commands are captured in [`sent`](PckReplay::sent) for inspection but do not
/// alter the replayed stream (a replay is inherently pre-recorded). This is the
/// no-hardware simulator behind `--simulate`.
pub struct PckReplay {
    records: std::vec::IntoIter<Record>,
    /// Commands the controller sent during the replay, in order.
    pub sent: Vec<Command>,
    /// The parsed session header (firmware, assay path, …).
    pub header: crate::tracehw::pck::Header,
}

impl PckReplay {
    /// Build a replay transport from an already-parsed session.
    pub fn new(session: Session) -> Self {
        PckReplay {
            header: session.header,
            records: session.records.into_iter(),
            sent: Vec::new(),
        }
    }

    /// Build a replay transport from a `.pck` file on disk.
    pub fn from_path(path: impl AsRef<std::path::Path>) -> Result<Self> {
        Ok(PckReplay::new(Session::read(path)?))
    }
}

impl Transport for PckReplay {
    fn send(&mut self, command: &Command) -> Result<()> {
        self.sent.push(command.clone());
        Ok(())
    }

    fn poll(&mut self) -> Result<Option<Record>> {
        Ok(self.records.next())
    }
}

/// Line-terminator appended to outbound commands (a serial Forth terminal
/// convention). **Unverified against hardware** — see the note on
/// [`StreamTransport`].
const COMMAND_TERMINATOR: &[u8] = b"\r";

/// A [`Transport`] over any byte stream (`Read + Write`) — the reusable core of
/// a **live serial/COM-port** link (milestone M6; the 2100 enumerates as a
/// virtual COM port — see `docs/bioanalyzer_protocol.md`).
///
/// The **receive** path is well-grounded: it deframes the same
/// `len | type | seq | payload | checksum` records the recorded `.pck` corpus
/// uses, resyncing on a bad length byte. The **send** path writes the command
/// text plus [`COMMAND_TERMINATOR`]; the exact outbound framing is not visible
/// in the corpus (only the instrument's echoes are), so it is a best guess and
/// **must be validated once hardware is available** before any real run (M7).
///
/// Binding this to a physical port only needs a `Read + Write` serial handle
/// (e.g. from the `serialport` crate) passed to [`StreamTransport::new`]; the
/// framing here is unchanged. It is deliberately not wired to real hardware yet.
pub struct StreamTransport<S> {
    stream: S,
    buf: Vec<u8>,
    scratch: [u8; 4096],
}

impl<S: Read + Write> StreamTransport<S> {
    /// Wrap a byte stream (a serial port handle, a TCP socket, or a test mock).
    pub fn new(stream: S) -> Self {
        StreamTransport {
            stream,
            buf: Vec::new(),
            scratch: [0u8; 4096],
        }
    }

    /// Pull whatever bytes are currently available into the internal buffer. A
    /// read timeout (a serial port with no data yet) is treated as zero bytes,
    /// not an error.
    fn fill(&mut self) -> Result<usize> {
        match self.stream.read(&mut self.scratch) {
            Ok(n) => {
                self.buf.extend_from_slice(&self.scratch[..n]);
                Ok(n)
            }
            Err(e) if matches!(e.kind(), ErrorKind::TimedOut | ErrorKind::WouldBlock) => Ok(0),
            Err(e) => Err(e).context("reading stream"),
        }
    }

    /// Extract one complete record from the buffer, resyncing past any leading
    /// garbage. Returns `None` if no complete record is buffered yet.
    fn take_record(&mut self) -> Option<Record> {
        loop {
            match parse_front(&self.buf) {
                FrontParse::Record(record, consumed) => {
                    self.buf.drain(..consumed);
                    return Some(record);
                }
                FrontParse::NeedMore => return None,
                FrontParse::Resync => {
                    // Bad length byte — drop one byte and retry.
                    self.buf.remove(0);
                }
            }
        }
    }
}

impl<S: Read + Write> Transport for StreamTransport<S> {
    fn send(&mut self, command: &Command) -> Result<()> {
        self.stream
            .write_all(command.wire().as_bytes())
            .context("writing command")?;
        self.stream
            .write_all(COMMAND_TERMINATOR)
            .context("writing command terminator")?;
        self.stream.flush().context("flushing command")
    }

    fn poll(&mut self) -> Result<Option<Record>> {
        if let Some(record) = self.take_record() {
            return Ok(Some(record));
        }
        // Nothing buffered — try to read more, then parse again.
        self.fill()?;
        Ok(self.take_record())
    }
}
