//! Tests for the protocol/transport/instrument stack (M2/M3).

use traceanalyzer::tracehw::pck::Record;
use traceanalyzer::tracehw::protocol::{decode, Event, Sample};
use traceanalyzer::tracehw::{
    run_to_completion, Acquisition, Command, PckReplay, RunState, Session,
};

fn record(type_id: u8, seq: u8, payload: &[u8]) -> Vec<u8> {
    let len = 3 + payload.len() + 1;
    let mut r = vec![len as u8, type_id, seq];
    r.extend_from_slice(payload);
    let ck = r.iter().fold(0u8, |a, b| a.wrapping_add(*b));
    r.push(ck);
    r
}

fn header(lines: &[&str]) -> Vec<u8> {
    let mut h: Vec<u8> = lines
        .join("\r\n")
        .encode_utf16()
        .flat_map(u16::to_le_bytes)
        .collect();
    h.resize(traceanalyzer::tracehw::pck::HEADER_LEN, 0);
    h
}

/// A synthetic session that exercises the full handshake and both sample streams.
fn synthetic_session() -> Session {
    let mut bytes = header(&["C.01.069", "C:\\assays\\dsDNA\\DNA 1000 Series II.xsy"]);
    bytes.extend(record(0x01, 0x01, b"AV.CLEAR  ok"));
    bytes.extend(record(0x01, 0x02, b"AV.SETUP  ok"));
    bytes.extend(record(0x01, 0x03, b"START  ok"));
    bytes.extend(record(0x04, 0x00, &[0x00, 0x00, 0x00, 0x0a, 0x00, 0x00])); // A: 10
    bytes.extend(record(0x06, 0x00, &[0x11, 0x05, 0x00, 0x00])); // B: +5
    bytes.extend(record(0x06, 0x01, &[0x12, 0xfb, 0xff, 0xff])); // B: -5 (signed)
    bytes.extend(record(0x02, 0x00, &[0x03, 0x01, 0x00])); // telemetry (ignored)
    bytes.extend(record(0x04, 0x6e, b"END_OF_RUN"));
    bytes.extend(record(0x01, 0x04, b"STOP  ok"));
    bytes.extend_from_slice(b"\r\nok");
    Session::parse(&bytes).expect("parse")
}

#[test]
fn decode_classifies_records() {
    let s = synthetic_session();
    let events: Vec<Event> = s.records.iter().map(decode).collect();
    assert_eq!(events[0], Event::Ack("AV.CLEAR  ok".into()));
    assert_eq!(
        events[3],
        Event::Sample(Sample {
            stream: 0x04,
            index: None,
            raw: 10
        })
    );
    assert_eq!(
        events[4],
        Event::Sample(Sample {
            stream: 0x06,
            index: Some(0x11),
            raw: 5
        })
    );
    assert_eq!(
        events[5],
        Event::Sample(Sample {
            stream: 0x06,
            index: Some(0x12),
            raw: -5
        })
    );
    assert_eq!(events[6], Event::Telemetry { type_id: 0x02 });
    assert_eq!(events[7], Event::EndOfRun);
}

#[test]
fn corrupt_record_is_flagged() {
    // len=5, type=0x06, seq=0, payload=[0x11], checksum deliberately wrong.
    let raw = [5u8, 0x06, 0x00, 0x11, 0x00];
    let rec = Record {
        offset: 0,
        type_id: raw[1],
        seq: raw[2],
        payload: vec![raw[3]],
        checksum: raw[4],
        checksum_ok: false,
    };
    assert_eq!(decode(&rec), Event::Corrupt { type_id: 0x06 });
}

#[test]
fn run_to_completion_drives_state_machine() {
    let mut transport = PckReplay::new(synthetic_session());
    let acq = run_to_completion(&mut transport).expect("run");

    // Handshake commands were issued.
    assert_eq!(
        transport.sent,
        vec![
            Command::AvClear,
            Command::AvSetup,
            Command::Start,
            Command::Stop
        ]
    );

    assert_eq!(acq.state, RunState::Stopped);
    assert_eq!(acq.stream_a, vec![10]);
    assert_eq!(acq.stream_b, vec![5, -5]);
    assert_eq!(acq.corrupt, 0);

    // The state machine passes through each phase as acks arrive.
    let mut a = Acquisition::new();
    a.apply(Event::Ack("AV.CLEAR  ok".into()));
    assert_eq!(a.state, RunState::Cleared);
    a.apply(Event::Ack("AV.SETUP  ok".into()));
    assert_eq!(a.state, RunState::SetUp);
    a.apply(Event::Ack("START  ok".into()));
    assert_eq!(a.state, RunState::Running);
    a.apply(Event::EndOfRun);
    assert_eq!(a.state, RunState::Ended);
}

#[test]
fn produces_electropherogram() {
    let mut transport = PckReplay::new(synthetic_session());
    let acq = run_to_completion(&mut transport).expect("run");
    let ep = acq.to_electrophoresis("test", 0.5);
    assert_eq!(ep.samples.len(), 1);
    assert_eq!(ep.samples[0].fluorescence, vec![5.0, -5.0]);
    assert_eq!(ep.samples[0].time, vec![0.0, 0.5]); // sample_period_s honored
}

/// End-to-end on a real capture when `ext_software/` is present.
#[test]
fn simulates_a_real_capture_when_available() {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/ext_software/2100 expert/data/packets"
    );
    let Some(path) = std::fs::read_dir(dir).ok().and_then(|entries| {
        entries.flatten().map(|e| e.path()).find(|p| {
            p.extension().and_then(|e| e.to_str()) == Some("pck")
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("2100 expert_"))
        })
    }) else {
        eprintln!("skipping: no real captures present");
        return;
    };

    let mut transport = PckReplay::from_path(&path).expect("open");
    let acq = run_to_completion(&mut transport).expect("run");
    // Captures end either on the `STOP  ok` echo or a bare `ok` after
    // `END_OF_RUN`, so both terminal states are valid.
    assert!(
        matches!(acq.state, RunState::Ended | RunState::Stopped),
        "unexpected terminal state {:?}",
        acq.state
    );
    assert!(acq.sample_count() > 0, "no samples collected");
    // Streams A and B are emitted 1:1 during acquisition (the log can end
    // mid-pair, so allow an off-by-one at the boundary).
    let diff = (acq.stream_a.len() as i64 - acq.stream_b.len() as i64).abs();
    assert!(
        diff <= 1,
        "streams not paired: {} vs {}",
        acq.stream_a.len(),
        acq.stream_b.len()
    );
    assert_eq!(acq.corrupt, 0);
    eprintln!(
        "{}: {} samples, errors {:?}",
        path.display(),
        acq.sample_count(),
        acq.errors.first()
    );
}
