//! Live serial/COM-port transport (milestone M6). Enabled by the `serial`
//! Cargo feature.
//!
//! The 2100 Bioanalyzer enumerates as a virtual serial port (see
//! `docs/bioanalyzer_protocol.md`). [`open`] binds a port to the
//! [`StreamTransport`](super::StreamTransport) record framing — the same code
//! proven against the recorded corpus, now over a live byte pipe.
//!
//! **Hardware-validation status:** the connect + deframe (receive) path is
//! complete and compiles; it opens any serial device and reads framed records.
//! The default baud rate and the outbound command framing still need
//! confirmation on a physical instrument before a real run (M7).
//!
//! Live read loop: unlike the [`PckReplay`](super::PckReplay) simulator, a
//! serial [`poll`](super::Transport::poll) returns `Ok(None)` whenever no record
//! has arrived *yet* (a read timeout), not only at end-of-run. A live driver
//! must therefore keep polling until it sees [`Event::EndOfRun`](crate::tracehw::Event)
//! (or a deadline elapses) rather than stopping at the first `None`.

use std::time::{Duration, Instant};

use anyhow::{bail, Context, Result};
use serialport::SerialPort;

use super::{StreamTransport, Transport};
use crate::tracehw::instrument::{Acquisition, RunState};
use crate::tracehw::protocol::{self, Command};
use crate::tracehw::safety::{Action, RunApproval};

/// A live serial transport: [`StreamTransport`] over a boxed serial port.
pub type SerialTransport = StreamTransport<Box<dyn SerialPort>>;

/// Default connection parameters. Baud is a placeholder pending a reading from a
/// live instrument; the firmware is an 8-N-1 text terminal.
pub const DEFAULT_BAUD: u32 = 115_200;

/// List candidate serial ports on the host (best effort; may be empty when the
/// enumeration backend is unavailable).
pub fn list_ports() -> Vec<String> {
    serialport::available_ports()
        .map(|ports| ports.into_iter().map(|p| p.port_name).collect())
        .unwrap_or_default()
}

/// Open a serial port and wrap it as a [`SerialTransport`]. Uses 8-N-1 and a
/// short read timeout so [`poll`](super::Transport::poll) is non-blocking.
///
/// This performs a **read-only-safe** connect: opening the port and reading
/// framed records sends nothing to the instrument. Issuing commands (and thus a
/// run) goes through [`guarded_run`](crate::tracehw::guarded_run), which gates every
/// energizing action.
pub fn open(port: &str, baud: u32) -> Result<SerialTransport> {
    let handle = serialport::new(port, baud)
        .data_bits(serialport::DataBits::Eight)
        .parity(serialport::Parity::None)
        .stop_bits(serialport::StopBits::One)
        .flow_control(serialport::FlowControl::None)
        .timeout(Duration::from_millis(250))
        .open()
        .with_context(|| format!("opening serial port {port} at {baud} baud"))?;
    Ok(StreamTransport::new(handle))
}

/// Guarded **live** run over a serial transport (milestone M7). Same safety
/// gating as [`guarded_run`](crate::tracehw::guarded_run) — `AV.SETUP` and `START` each
/// require approval — but the acquisition loop keeps polling across read
/// timeouts until `END_OF_RUN` or `timeout` elapses (a live port yields no
/// records between samples), instead of stopping at the first `None`.
///
/// Not yet exercised against a physical instrument; the outbound command framing
/// and baud rate must be confirmed on hardware first.
pub fn guarded_run_live(
    transport: &mut SerialTransport,
    approval: &mut impl RunApproval,
    timeout: Duration,
) -> Result<Acquisition> {
    transport.send(&Command::AvClear)?;

    if !approval.approve(Action::Setup) {
        transport.send(&Command::Stop)?;
        bail!("run aborted: AV.SETUP was not approved");
    }
    transport.send(&Command::AvSetup)?;

    if !approval.approve(Action::Start) {
        transport.send(&Command::Stop)?;
        bail!("run aborted: START was not approved");
    }
    transport.send(&Command::Start)?;

    let mut acq = Acquisition::new();
    let deadline = Instant::now() + timeout;
    while acq.state != RunState::Ended && Instant::now() < deadline {
        match transport.poll()? {
            Some(record) => acq.apply(protocol::decode(&record)),
            None => std::thread::sleep(Duration::from_millis(10)),
        }
    }

    transport.send(&Command::Stop)?;
    if acq.state != RunState::Ended {
        bail!("run did not reach END_OF_RUN within {timeout:?}");
    }
    Ok(acq)
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use crate::tracehw::protocol::{decode, Event};
    use std::io::Write;

    fn record(type_id: u8, seq: u8, payload: &[u8]) -> Vec<u8> {
        let len = 3 + payload.len() + 1;
        let mut r = vec![len as u8, type_id, seq];
        r.extend_from_slice(payload);
        let ck = r.iter().fold(0u8, |a, b| a.wrapping_add(*b));
        r.push(ck);
        r
    }

    /// End-to-end over a real pseudo-terminal pair: bytes written to one end are
    /// deframed into records at the other, exercising the actual serial handle
    /// (not a mock). This is the M6 connect + receive path minus the instrument.
    #[test]
    fn deframes_over_a_real_pty() {
        let (mut host, mut device) = serialport::TTYPort::pair().expect("openpty");
        device
            .set_timeout(Duration::from_millis(200))
            .expect("timeout");

        // The "instrument" side writes framed records into the port.
        host.write_all(&record(0x01, 0x01, b"START  ok")).unwrap();
        host.write_all(&record(0x04, 0x6e, b"END_OF_RUN")).unwrap();
        host.flush().unwrap();

        let mut transport = StreamTransport::new(Box::new(device) as Box<dyn SerialPort>);
        let first = transport.poll().unwrap().expect("first record");
        assert_eq!(decode(&first), Event::Ack("START  ok".into()));
        let second = transport.poll().unwrap().expect("second record");
        assert_eq!(decode(&second), Event::EndOfRun);
    }
}
