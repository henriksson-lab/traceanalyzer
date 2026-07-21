//! Tests for the byte-stream framing transport (M6 receive/send core).
//!
//! These exercise the framing over an in-memory duplex mock — the same code
//! that will run over a real serial port, minus the port binding.

use std::io::{Cursor, Read, Write};

use traceanalyzer::tracehw::protocol::{decode, Event};
use traceanalyzer::tracehw::{Command, StreamTransport, Transport};

/// A `Read + Write` mock: reads from a fixed inbound buffer, records outbound
/// bytes for inspection.
struct Duplex {
    inbound: Cursor<Vec<u8>>,
    outbound: Vec<u8>,
}

impl Read for Duplex {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        self.inbound.read(out)
    }
}
impl Write for Duplex {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        self.outbound.extend_from_slice(data);
        Ok(data.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn record(type_id: u8, seq: u8, payload: &[u8]) -> Vec<u8> {
    let len = 3 + payload.len() + 1;
    let mut r = vec![len as u8, type_id, seq];
    r.extend_from_slice(payload);
    let ck = r.iter().fold(0u8, |a, b| a.wrapping_add(*b));
    r.push(ck);
    r
}

#[test]
fn deframes_records_from_a_stream() {
    let mut wire = Vec::new();
    wire.extend(record(0x01, 0x01, b"START  ok"));
    wire.extend(record(0x06, 0x00, &[0x11, 0x05, 0x00, 0x00]));
    wire.extend(record(0x04, 0x6e, b"END_OF_RUN"));

    let mut t = StreamTransport::new(Duplex {
        inbound: Cursor::new(wire),
        outbound: Vec::new(),
    });

    let mut events = Vec::new();
    while let Some(rec) = t.poll().expect("poll") {
        events.push(decode(&rec));
    }
    assert_eq!(events[0], Event::Ack("START  ok".into()));
    assert!(matches!(events[1], Event::Sample(_)));
    assert_eq!(events[2], Event::EndOfRun);
}

#[test]
fn resyncs_past_leading_garbage() {
    let mut wire = vec![0x00, 0x01, 0x02]; // junk with sub-minimal length bytes
    wire.extend(record(0x01, 0x01, b"AV.CLEAR  ok"));

    let mut t = StreamTransport::new(Duplex {
        inbound: Cursor::new(wire),
        outbound: Vec::new(),
    });
    let rec = t.poll().expect("poll").expect("record after resync");
    assert_eq!(decode(&rec), Event::Ack("AV.CLEAR  ok".into()));
}

#[test]
fn send_output_is_captured() {
    // Own the mock so we can inspect what was written.
    let mut mock = Duplex {
        inbound: Cursor::new(Vec::new()),
        outbound: Vec::new(),
    };
    {
        let mut t = StreamTransport::new(&mut mock);
        t.send(&Command::AvClear).unwrap();
        t.send(&Command::Start).unwrap();
    }
    assert_eq!(mock.outbound, b"AV.CLEAR\rSTART\r");
}
