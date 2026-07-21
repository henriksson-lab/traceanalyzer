//! Reader for 2100 Expert **`.pck`** recorded instrument-session captures.
//!
//! A `.pck` is the byte-for-byte log 2100 Expert writes of one instrument
//! session (its `data/packets/` folder). It is the reverse-engineering corpus
//! for the live protocol — decoding it needs no hardware. See
//! [`docs/pck_format.md`](../../docs/pck_format.md) for the full byte-level
//! notes; the essentials this reader relies on:
//!
//! * A fixed **2048-byte header** ([`HEADER_LEN`]). For assay runs it is
//!   UTF-16LE text, NUL-padded: line 0 = firmware version (`C.01.069`), line 1 =
//!   the absolute `.xsy` assay path, then a few assay fields. Diagnostics
//!   captures leave it non-text.
//! * The rest of the file is a flat stream of **records**, each:
//!   `len:u8 | type:u8 | seq:u8 | payload… | checksum:u8`, where `len` is the
//!   whole record length (so `payload` is `len - 4` bytes) and `checksum` is the
//!   **sum of every preceding byte of the record, mod 256**. `seq` increments
//!   per record type (wrapping). The stream ends with a short literal trailer
//!   (typically `b"\r\nok"`) that is not a record.
//!
//! Record `type` is overloaded and only partly understood (a *first cut* — see
//! [`type_label`]): `0x01` carries text command echoes / responses
//! (`AV.CLEAR  ok`, `START  ok`, `GET.ERROR …`) and the `END_OF_RUN` event;
//! `0x04` and `0x06` are the two high-rate per-sample streams (emitted 1:1
//! during acquisition); `0x02`/`0x03`/`0x05` are lower-rate telemetry. Turning
//! the sample payloads into calibrated fluorescence is deferred to protocol
//! Phase P2, so this reader exposes records faithfully rather than guessing.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

/// Fixed size of the leading header block, in bytes. Records begin here.
pub const HEADER_LEN: usize = 2048;

/// Minimum bytes for a well-formed record: `len | type | seq | checksum`.
const MIN_RECORD_LEN: usize = 4;

/// The 2048-byte session header.
#[derive(Debug, Clone)]
pub struct Header {
    /// The raw header bytes, verbatim (always [`HEADER_LEN`] long when present).
    pub raw: Vec<u8>,
    /// `true` when the header decoded as UTF-16LE text (assay runs); `false` for
    /// non-text headers (some diagnostics captures).
    pub is_text: bool,
    /// Firmware version string (header line 0), e.g. `"C.01.069"`.
    pub firmware_version: Option<String>,
    /// Absolute path to the `.xsy` assay used for the run (header line 1).
    pub assay_path: Option<String>,
    /// All decoded header text lines (includes the two fields above at 0/1).
    pub fields: Vec<String>,
}

impl Header {
    fn parse(raw: &[u8]) -> Header {
        // The text region is UTF-16LE, NUL-padded to `HEADER_LEN`. Every ASCII
        // char is `XX 00`, so a `0x0000` code unit only appears at the padding.
        let mut units = Vec::new();
        let mut i = 0;
        while i + 1 < raw.len() {
            let u = u16::from_le_bytes([raw[i], raw[i + 1]]);
            if u == 0 {
                break;
            }
            units.push(u);
            i += 2;
        }
        let text = String::from_utf16_lossy(&units);

        // Treat the header as text only if it is convincingly printable — this
        // rejects the binary headers of some diagnostics captures.
        let printable = text
            .chars()
            .filter(|c| *c == '\t' || !c.is_control())
            .count();
        let is_text = !text.is_empty() && printable * 100 >= text.chars().count() * 90;

        let fields: Vec<String> = if is_text {
            text.split("\r\n").map(str::to_string).collect()
        } else {
            Vec::new()
        };
        let firmware_version = fields.first().cloned().filter(|s| !s.is_empty());
        let assay_path = fields.get(1).cloned().filter(|s| !s.is_empty());
        Header {
            raw: raw.to_vec(),
            is_text,
            firmware_version,
            assay_path,
            fields,
        }
    }
}

/// One framed record from the session stream.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Byte offset of the record's `len` byte within the file (for diagnostics).
    pub offset: usize,
    /// Record type byte (see [`type_label`]).
    pub type_id: u8,
    /// Per-type sequence counter (wraps at 256).
    pub seq: u8,
    /// Record payload, excluding the trailing checksum byte.
    pub payload: Vec<u8>,
    /// The stored checksum byte.
    pub checksum: u8,
    /// Whether [`checksum`](Record::checksum) matches the recomputed sum.
    pub checksum_ok: bool,
}

impl Record {
    /// The payload decoded as text, if it is entirely printable ASCII (ignoring
    /// a possible trailing NUL). Command echoes and events like `START  ok` and
    /// `END_OF_RUN` come back as `Some`; binary sample records as `None`.
    pub fn as_text(&self) -> Option<String> {
        let bytes = match self.payload.split_last() {
            Some((&0, rest)) => rest, // tolerate a single trailing NUL
            _ => &self.payload[..],
        };
        if bytes.is_empty() {
            return None;
        }
        if bytes.iter().all(|&b| (0x20..0x7f).contains(&b)) {
            Some(String::from_utf8_lossy(bytes).into_owned())
        } else {
            None
        }
    }

    /// Whether this record carries a printable text line (as opposed to binary).
    pub fn is_text(&self) -> bool {
        self.as_text().is_some()
    }
}

/// Outcome of scanning to the end of the record stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamEnd {
    /// The stream ended cleanly on a record boundary; only the trailer remained.
    Clean,
    /// A record claimed more bytes than were left in the file (truncated log).
    Truncated,
}

/// A fully parsed `.pck` session.
#[derive(Debug, Clone)]
pub struct Session {
    /// The 2048-byte header.
    pub header: Header,
    /// All records, in file order.
    pub records: Vec<Record>,
    /// Leftover bytes after the last full record (the literal end trailer, and
    /// any truncated tail).
    pub trailer: Vec<u8>,
    /// How the record stream terminated.
    pub end: StreamEnd,
}

impl Session {
    /// Parse a `.pck` from its raw bytes.
    pub fn parse(bytes: &[u8]) -> Result<Session> {
        anyhow::ensure!(
            bytes.len() >= HEADER_LEN,
            "file is {} bytes, shorter than the {}-byte header",
            bytes.len(),
            HEADER_LEN
        );
        let header = Header::parse(&bytes[..HEADER_LEN]);

        let mut records = Vec::new();
        let mut pos = HEADER_LEN;
        let end = loop {
            let remaining = bytes.len() - pos;
            if remaining == 0 {
                break StreamEnd::Clean;
            }
            if remaining < MIN_RECORD_LEN {
                // Too short to be a record — the literal end trailer (e.g. b"\r\nok").
                break StreamEnd::Clean;
            }
            let len = bytes[pos] as usize;
            if len < MIN_RECORD_LEN {
                // Not a plausible record length: treat the rest as trailer.
                break StreamEnd::Clean;
            }
            if pos + len > bytes.len() {
                // A record that overruns the buffer is either the printable end
                // footer (e.g. b"\r\nok", whose first byte 0x0d reads as a bogus
                // length) or a genuinely truncated binary record. Distinguish by
                // whether the remaining bytes are all printable/whitespace.
                let rest = &bytes[pos..];
                let footer = rest
                    .iter()
                    .all(|&b| (0x20..0x7f).contains(&b) || matches!(b, b'\r' | b'\n' | b'\t'));
                break if footer {
                    StreamEnd::Clean
                } else {
                    StreamEnd::Truncated
                };
            }
            let rec = &bytes[pos..pos + len];
            let checksum = rec[len - 1];
            let computed = rec[..len - 1]
                .iter()
                .fold(0u8, |acc, &b| acc.wrapping_add(b));
            records.push(Record {
                offset: pos,
                type_id: rec[1],
                seq: rec[2],
                payload: rec[3..len - 1].to_vec(),
                checksum,
                checksum_ok: computed == checksum,
            });
            pos += len;
        };

        Ok(Session {
            header,
            records,
            trailer: bytes[pos..].to_vec(),
            end,
        })
    }

    /// Read and parse a `.pck` from a path.
    pub fn read(path: impl AsRef<Path>) -> Result<Session> {
        let path = path.as_ref();
        let bytes =
            fs::read(path).with_context(|| format!("reading .pck file {}", path.display()))?;
        Session::parse(&bytes).with_context(|| format!("parsing .pck file {}", path.display()))
    }

    /// Count of records per type id, ascending by type id.
    pub fn type_census(&self) -> Vec<(u8, usize)> {
        let mut counts: std::collections::BTreeMap<u8, usize> = std::collections::BTreeMap::new();
        for r in &self.records {
            *counts.entry(r.type_id).or_default() += 1;
        }
        counts.into_iter().collect()
    }

    /// Number of records whose stored checksum did not verify (0 for a healthy log).
    pub fn checksum_failures(&self) -> usize {
        self.records.iter().filter(|r| !r.checksum_ok).count()
    }

    /// The ordered text lines carried by the session — command echoes, responses
    /// and events such as `START  ok`, `GET.ERROR …`, `END_OF_RUN`. This is the
    /// human-readable command timeline of the run.
    pub fn text_timeline(&self) -> Vec<String> {
        self.records.iter().filter_map(Record::as_text).collect()
    }
}

/// Outcome of [`parse_front`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FrontParse {
    /// A complete record and the number of bytes it consumed.
    Record(Record, usize),
    /// Not enough bytes are buffered yet — read more and retry.
    NeedMore,
    /// The leading byte is not a valid record length; drop one byte and resync.
    Resync,
}

/// Try to parse one record from the front of a byte buffer (e.g. a live serial
/// stream), without the file-level header/trailer handling of [`Session::parse`].
pub fn parse_front(buf: &[u8]) -> FrontParse {
    if buf.is_empty() {
        return FrontParse::NeedMore;
    }
    let len = buf[0] as usize;
    if len < MIN_RECORD_LEN {
        return FrontParse::Resync;
    }
    if buf.len() < len {
        return FrontParse::NeedMore;
    }
    let rec = &buf[..len];
    let checksum = rec[len - 1];
    let computed = rec[..len - 1]
        .iter()
        .fold(0u8, |acc, &b| acc.wrapping_add(b));
    FrontParse::Record(
        Record {
            offset: 0,
            type_id: rec[1],
            seq: rec[2],
            payload: rec[3..len - 1].to_vec(),
            checksum,
            checksum_ok: computed == checksum,
        },
        len,
    )
}

/// Best-guess human label for a record type id. Types are only partly decoded
/// (protocol Phase P0/P2) — labels marked "?" are provisional.
pub fn type_label(type_id: u8) -> &'static str {
    match type_id {
        0x01 => "text (command echo / response / GET.ERROR)",
        0x02 => "telemetry? (low-rate)",
        0x03 => "aux? (low-rate)",
        0x04 => "sample stream A (high-rate; also END_OF_RUN event)",
        0x05 => "aux? (mid-rate)",
        0x06 => "sample stream B (high-rate)",
        _ => "unknown",
    }
}
