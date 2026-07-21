//! Protocol layer: the command set the instrument speaks, and decoding of
//! inbound records into semantic [`Event`]s.
//!
//! The wire protocol is a line-oriented text command set (Forth words the
//! `G2938C` firmware exposes) plus a stream of binary acquisition records. Both
//! directions are framed identically on disk — see [`crate::tracehw::pck`] and
//! [`docs/pck_format.md`](../../docs/pck_format.md). This module names the
//! commands seen in the recorded corpus and turns a raw [`Record`] into an
//! [`Event`]. Physical calibration of the sample payloads is protocol Phase P2;
//! until then [`Sample::raw`] is a provisional integer, not a fluorescence value.

use crate::tracehw::pck::Record;

/// A command sent PC → instrument. `wire()` renders the on-the-wire text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Command {
    /// `AV.CLEAR` — clear/prime the acquisition.
    AvClear,
    /// `AV.SETUP` — configure the run from the loaded assay.
    AvSetup,
    /// `START` — begin electrophoresis and acquisition.
    Start,
    /// `STOP` — end acquisition.
    Stop,
    /// `GET.ERROR` — poll the instrument's error status word.
    GetError,
    /// Any other raw command line (e.g. a script-download line).
    Raw(String),
}

impl Command {
    /// The exact text put on the wire (without the record framing).
    pub fn wire(&self) -> String {
        match self {
            Command::AvClear => "AV.CLEAR".into(),
            Command::AvSetup => "AV.SETUP".into(),
            Command::Start => "START".into(),
            Command::Stop => "STOP".into(),
            Command::GetError => "GET.ERROR".into(),
            Command::Raw(s) => s.clone(),
        }
    }
}

/// A provisional decode of one high-rate acquisition record. The `raw` value is
/// **not calibrated fluorescence** — decoding the `RawDataCu` encoding is
/// Phase P2. See [`docs/pck_format.md`](../../docs/pck_format.md).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sample {
    /// Source stream: `0x04` (stream A) or `0x06` (stream B).
    pub stream: u8,
    /// The per-sample counter byte, present on stream `0x06`.
    pub index: Option<u8>,
    /// Provisional raw integer pulled from the payload (uncalibrated).
    pub raw: i64,
}

/// A decoded inbound event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    /// A command acknowledgement / response line, e.g. `AV.CLEAR  ok`.
    Ack(String),
    /// A `GET.ERROR` response with its status word, e.g. `ERR16050015`.
    Error(String),
    /// The end-of-run event (`END_OF_RUN`).
    EndOfRun,
    /// A high-rate acquisition sample (provisional decode).
    Sample(Sample),
    /// A lower-rate telemetry / aux record we do not decode yet.
    Telemetry { type_id: u8 },
    /// A record whose checksum did not verify.
    Corrupt { type_id: u8 },
}

/// Little-endian unsigned integer from a byte slice (up to 8 bytes).
fn le(bytes: &[u8]) -> i64 {
    let mut v = 0i64;
    for (i, &b) in bytes.iter().enumerate() {
        v |= i64::from(b) << (8 * i);
    }
    v
}

/// Sign-extend an `n`-byte little-endian value.
fn le_signed(bytes: &[u8]) -> i64 {
    let n = bytes.len();
    let raw = le(bytes);
    let bits = 8 * n as u32;
    if bits < 64 && raw & (1 << (bits - 1)) != 0 {
        raw - (1 << bits)
    } else {
        raw
    }
}

/// Decode one framed record into a semantic [`Event`].
pub fn decode(record: &Record) -> Event {
    if !record.checksum_ok {
        return Event::Corrupt {
            type_id: record.type_id,
        };
    }
    if let Some(text) = record.as_text() {
        let trimmed = text.trim();
        if trimmed == "END_OF_RUN" {
            return Event::EndOfRun;
        }
        if let Some(rest) = trimmed.strip_prefix("GET.ERROR") {
            // e.g. "GET.ERROR ERR16050015 ok" -> "ERR16050015"
            let code = rest.split_whitespace().next().unwrap_or("").to_string();
            return Event::Error(code);
        }
        return Event::Ack(trimmed.to_string());
    }
    // Binary acquisition / telemetry records.
    let p = &record.payload;
    match record.type_id {
        // Stream A: 3 leading zero bytes + a 24-bit value (last 3 payload bytes).
        0x04 if p.len() >= 6 => Event::Sample(Sample {
            stream: 0x04,
            index: None,
            raw: le(&p[p.len() - 3..]),
        }),
        // Stream B: a 1-byte counter + a signed 24-bit value.
        0x06 if p.len() >= 4 => Event::Sample(Sample {
            stream: 0x06,
            index: Some(p[0]),
            raw: le_signed(&p[1..4]),
        }),
        other => Event::Telemetry { type_id: other },
    }
}
