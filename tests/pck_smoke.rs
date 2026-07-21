//! Parser tests for the `.pck` recorded-session reader.
//!
//! The synthetic test builds a well-formed `.pck` in memory (so it runs
//! anywhere), exercising the header, framing, checksum verification and text
//! extraction. A second, data-gated test parses a real capture from
//! `ext_software/` when present — that folder is not committed, so the test
//! quietly skips when it is absent.

use traceanalyzer::tracehw::pck::{Session, StreamEnd, HEADER_LEN};

/// Build one record: `len | type | seq | payload | checksum`.
fn record(type_id: u8, seq: u8, payload: &[u8]) -> Vec<u8> {
    let len = 3 + payload.len() + 1;
    assert!(len <= u8::MAX as usize);
    let mut r = vec![len as u8, type_id, seq];
    r.extend_from_slice(payload);
    let ck = r.iter().fold(0u8, |acc, b| acc.wrapping_add(*b));
    r.push(ck);
    r
}

/// A 2048-byte UTF-16LE header from CRLF-joined text lines, NUL-padded.
fn header(lines: &[&str]) -> Vec<u8> {
    let text = lines.join("\r\n");
    let mut h: Vec<u8> = text.encode_utf16().flat_map(u16::to_le_bytes).collect();
    h.resize(HEADER_LEN, 0);
    h
}

#[test]
fn parses_synthetic_session() {
    let mut bytes = header(&[
        "C.01.069",
        "C:\\assays\\dsDNA\\DNA 1000 Series II.xsy",
        "1",
        "29",
        "0.050000 0.200000 1.000000 1.600000",
    ]);
    bytes.extend(record(0x01, 0x9e, b"AV.CLEAR  ok"));
    bytes.extend(record(0x01, 0x9f, b"AV.SETUP  ok"));
    bytes.extend(record(0x01, 0xa0, b"START  ok"));
    bytes.extend(record(0x04, 0x00, &[0x00, 0x00, 0x00, 0x09, 0xcd, 0xc0]));
    bytes.extend(record(0x06, 0x00, &[0x11, 0x00, 0x00, 0x00]));
    bytes.extend(record(0x04, 0x01, &[0x00, 0x00, 0x00, 0x04, 0xfc, 0x49]));
    bytes.extend(record(0x04, 0x6e, b"END_OF_RUN"));
    bytes.extend(record(0x01, 0xa1, b"STOP  ok"));
    bytes.extend_from_slice(b"\r\nok"); // literal end trailer

    let s = Session::parse(&bytes).expect("parse");

    // Header.
    assert!(s.header.is_text);
    assert_eq!(s.header.firmware_version.as_deref(), Some("C.01.069"));
    assert_eq!(
        s.header.assay_path.as_deref(),
        Some("C:\\assays\\dsDNA\\DNA 1000 Series II.xsy")
    );
    assert_eq!(s.header.fields.get(3).map(String::as_str), Some("29"));

    // Framing + checksums.
    assert_eq!(s.records.len(), 8);
    assert_eq!(s.checksum_failures(), 0);
    assert_eq!(s.end, StreamEnd::Clean);
    assert_eq!(s.trailer, b"\r\nok");

    // Census.
    let census = s.type_census();
    assert_eq!(census, vec![(0x01, 4), (0x04, 3), (0x06, 1)]);

    // Text timeline picks up commands, responses and the END_OF_RUN event, but
    // not the binary sample records.
    assert_eq!(
        s.text_timeline(),
        vec![
            "AV.CLEAR  ok",
            "AV.SETUP  ok",
            "START  ok",
            "END_OF_RUN",
            "STOP  ok",
        ]
    );

    // A binary sample record is not mistaken for text.
    let sample = s
        .records
        .iter()
        .find(|r| r.type_id == 0x04 && r.seq == 0x00);
    assert!(sample.unwrap().as_text().is_none());
}

#[test]
fn flags_a_corrupt_checksum() {
    let mut bytes = header(&["C.01.069"]);
    let mut rec = record(0x01, 0x01, b"START  ok");
    *rec.last_mut().unwrap() ^= 0xff; // corrupt the checksum
    bytes.extend(rec);

    let s = Session::parse(&bytes).expect("parse");
    assert_eq!(s.records.len(), 1);
    assert_eq!(s.checksum_failures(), 1);
    assert!(!s.records[0].checksum_ok);
}

#[test]
fn truncated_final_record_is_reported() {
    let mut bytes = header(&["C.01.069"]);
    // A record claiming length 20 but only a few bytes follow.
    bytes.extend_from_slice(&[20, 0x04, 0x00, 0x11, 0x22]);

    let s = Session::parse(&bytes).expect("parse");
    assert_eq!(s.records.len(), 0);
    assert_eq!(s.end, StreamEnd::Truncated);
}

/// Parse every real capture in `ext_software/` if that (uncommitted) folder is
/// present, asserting each frames cleanly with valid checksums.
#[test]
fn parses_real_captures_when_available() {
    let dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/ext_software/2100 expert/data/packets"
    );
    let entries = match std::fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => {
            eprintln!("skipping: {dir} not present");
            return;
        }
    };
    let mut checked = 0;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("pck") {
            continue;
        }
        let s = Session::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        // Strict framing/checksum checks apply to assay-run captures (named
        // "2100 expert_…"). Diagnostics captures ("HardwareDiagnosis_…") use a
        // different, not-yet-decoded layout and are only required to parse
        // without error.
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if name.starts_with("2100 expert_") {
            assert_eq!(
                s.checksum_failures(),
                0,
                "{}: checksum failures",
                path.display()
            );
            assert!(
                !s.records.is_empty(),
                "{}: no records parsed",
                path.display()
            );
            assert_eq!(
                s.header.firmware_version.as_deref(),
                Some("C.01.069"),
                "{}: unexpected firmware",
                path.display()
            );
            checked += 1;
        }
    }
    eprintln!("parsed {checked} real .pck captures");
    assert!(checked > 0, "packets folder present but held no .pck files");
}
