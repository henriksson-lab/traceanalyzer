//! Round-trip tests for `traceanalyzer::traceio::save::save_run` against real demo fixtures.
//! The demo `.xml.gz` files are gitignored; fetch with `bash scripts/fetch-testdata.sh`.
//! Tests skip gracefully when a fixture is absent (matching bioanalyzer_demo.rs).

use std::io::Read;
use std::path::PathBuf;

fn testdata(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join(name)
}

/// Load a gz fixture's XML, or `None` (with a skip note) if it is not present.
fn load_gz_xml(name: &str) -> Option<String> {
    let path = testdata(name);
    let raw = match std::fs::read(&path) {
        Ok(raw) => raw,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("skipping: fixture {} is not present", path.display());
            return None;
        }
        Err(e) => panic!("read fixture {}: {e}", path.display()),
    };
    let mut d = flate2::read::GzDecoder::new(&raw[..]);
    let mut s = String::new();
    d.read_to_string(&mut s).expect("gunzip fixture");
    Some(s)
}

/// A unique temp path so parallel test runs don't collide.
fn temp_xml(tag: &str) -> PathBuf {
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("traceio_save_{tag}_{pid}_{nanos}.xml"))
}

#[test]
fn round_trip_rename_from_gz_to_xml() {
    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let mut run = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse demo");
    assert!(run.samples.len() >= 2, "need >=2 samples for this test");

    // Remember an untouched sample and the original count.
    let original_count = run.samples.len();
    let other_idx = run.samples.len() - 1;
    let other_name = run.samples[other_idx].name.clone();
    assert_ne!(run.samples[0].name, "RENAMED_TEST_WELL");

    run.samples[0].name = "RENAMED_TEST_WELL".to_string();

    let src = testdata("demo_dna1000.xml.gz");
    let dst = temp_xml("roundtrip");
    traceanalyzer::traceio::save::save_run(&run, &src, &dst).expect("save_run");

    let reloaded_xml = std::fs::read_to_string(&dst).expect("read back");
    let reloaded = traceanalyzer::traceio::bioanalyzer::parse_xml(&reloaded_xml).expect("reparse");
    let _ = std::fs::remove_file(&dst);

    assert_eq!(
        reloaded.samples.len(),
        original_count,
        "sample count changed"
    );
    assert_eq!(
        reloaded.samples[0].name, "RENAMED_TEST_WELL",
        "rename did not persist"
    );
    assert_eq!(
        reloaded.samples[other_idx].name, other_name,
        "an unrelated sample name changed"
    );
}

#[test]
fn round_trip_xml_escaping() {
    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let mut run = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse demo");
    assert!(!run.samples.is_empty());

    let tricky = "A & B <x>";
    run.samples[0].name = tricky.to_string();

    let src = testdata("demo_dna1000.xml.gz");
    let dst = temp_xml("escape");
    traceanalyzer::traceio::save::save_run(&run, &src, &dst).expect("save_run");

    let reloaded_xml = std::fs::read_to_string(&dst).expect("read back");
    let reloaded = traceanalyzer::traceio::bioanalyzer::parse_xml(&reloaded_xml).expect("reparse");
    let _ = std::fs::remove_file(&dst);

    assert_eq!(
        reloaded.samples[0].name, tricky,
        "escaped name did not survive the round-trip exactly"
    );
}

#[test]
fn saving_as_xad_is_rejected() {
    // No fixture needed: the .xad rejection happens before src is read.
    let run = traceanalyzer::traceio::model::Electrophoresis {
        assay: Default::default(),
        ladder_peaks: vec![],
        regions: vec![],
        samples: vec![],
    };
    let err = traceanalyzer::traceio::save::save_run(
        &run,
        &PathBuf::from("whatever.xml"),
        &PathBuf::from("out.xad"),
    )
    .expect_err("save-as-.xad must error");
    let msg = format!("{err}");
    assert!(msg.contains(".xad"), "unhelpful error: {msg}");
}
