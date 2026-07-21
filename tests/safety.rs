//! Tests for the run-safety gating (M7). Exercised against the simulator; a
//! live run additionally needs hardware validation of the transport.

use traceanalyzer::tracehw::pck::HEADER_LEN;
use traceanalyzer::tracehw::safety::{guarded_run, Action, ApproveWith};
use traceanalyzer::tracehw::{ApproveAll, Command, DenyAll, PckReplay, RunState, Session};

fn record(type_id: u8, seq: u8, payload: &[u8]) -> Vec<u8> {
    let len = 3 + payload.len() + 1;
    let mut r = vec![len as u8, type_id, seq];
    r.extend_from_slice(payload);
    let ck = r.iter().fold(0u8, |a, b| a.wrapping_add(*b));
    r.push(ck);
    r
}

fn session() -> Session {
    let mut bytes = vec![0u8; HEADER_LEN];
    bytes.extend(record(0x01, 0x01, b"AV.CLEAR  ok"));
    bytes.extend(record(0x01, 0x02, b"AV.SETUP  ok"));
    bytes.extend(record(0x01, 0x03, b"START  ok"));
    bytes.extend(record(0x06, 0x00, &[0x11, 0x05, 0x00, 0x00]));
    bytes.extend(record(0x04, 0x6e, b"END_OF_RUN"));
    Session::parse(&bytes).unwrap()
}

#[test]
fn deny_all_aborts_before_setup() {
    let mut t = PckReplay::new(session());
    let err = guarded_run(&mut t, &mut DenyAll).unwrap_err();
    assert!(err.to_string().contains("AV.SETUP"), "{err}");
    // Only CLEAR then the safety STOP were sent — never SETUP or START.
    assert_eq!(t.sent, vec![Command::AvClear, Command::Stop]);
}

#[test]
fn refusing_start_aborts_after_setup() {
    let mut t = PckReplay::new(session());
    let mut approver = ApproveWith(|action| action != Action::Start);
    let err = guarded_run(&mut t, &mut approver).unwrap_err();
    assert!(err.to_string().contains("START"), "{err}");
    // SETUP happened (approved) but START never did.
    assert_eq!(
        t.sent,
        vec![Command::AvClear, Command::AvSetup, Command::Stop]
    );
}

#[test]
fn approve_all_runs_to_end() {
    let mut t = PckReplay::new(session());
    let acq = guarded_run(&mut t, &mut ApproveAll).unwrap();
    assert_eq!(acq.state, RunState::Ended);
    assert_eq!(acq.stream_b, vec![5]);
    assert_eq!(
        t.sent,
        vec![
            Command::AvClear,
            Command::AvSetup,
            Command::Start,
            Command::Stop
        ]
    );
}
