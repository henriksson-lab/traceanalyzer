//! Validation against real bioanalyzeR demo exports (jwfoley/bioanalyzeR, MIT).
//! The demo `.xml.gz` files are the exported-XML form of real runs, whose
//! schema is identical to the inner XML of native `.xad` files.

use std::io::Read;
use std::path::PathBuf;

fn testdata(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata")
        .join(name)
}

fn load_gz_xml(name: &str) -> String {
    let raw = std::fs::read(testdata(name)).expect("read fixture");
    let mut d = flate2::read::GzDecoder::new(&raw[..]);
    let mut s = String::new();
    d.read_to_string(&mut s).expect("gunzip fixture");
    s
}

#[test]
fn parses_dna1000_demo() {
    let xml = load_gz_xml("demo_dna1000.xml.gz");
    let run = traceio::bioanalyzer::parse_xml(&xml).expect("parse DNA 1000");

    assert_eq!(run.assay.length_unit, "bp");
    assert_eq!(run.assay.concentration_unit, "ng/µl");
    assert_eq!(run.assay.molarity_unit.as_deref(), Some("nM"));
    assert!(run.assay.has_upper_marker, "DNA 1000 has an upper marker");

    // Ladder is defined and starts at the 15 bp lower marker.
    assert!(
        run.ladder_peaks.len() >= 10,
        "expected a full ladder, got {}",
        run.ladder_peaks.len()
    );
    assert_eq!(run.ladder_peaks[0].size, 15.0);

    // Samples parsed, exactly one ladder well, traces populated.
    assert!(!run.samples.is_empty(), "no samples parsed");
    assert!(run.ladder_index().is_some(), "expected exactly one ladder well");

    let ladder = &run.samples[run.ladder_index().unwrap()];
    assert!(ladder.is_ladder);
    assert!(
        ladder.time.len() == ladder.fluorescence.len() && !ladder.fluorescence.is_empty(),
        "ladder trace time/fluorescence must be non-empty and equal length"
    );
    assert!(!ladder.peaks.is_empty(), "ladder should have called peaks");

    // Trace fluorescence must be finite and vary (not all zero).
    let (lo, hi) = ladder.fluorescence_range().expect("range");
    assert!(lo.is_finite() && hi.is_finite() && hi > lo, "trace is flat/NaN");

    // A DNA assay carries no RIN.
    assert!(run.samples.iter().all(|s| s.rin.is_none()));
}

#[test]
fn parses_rna_nano_demo_with_rin() {
    let xml = load_gz_xml("demo_rna_nano.xml.gz");
    let run = traceio::bioanalyzer::parse_xml(&xml).expect("parse RNA Nano");

    assert_eq!(run.assay.length_unit, "nt");
    assert!(!run.samples.is_empty());

    // At least one non-ladder RNA sample should report a RIN in [1, 10].
    let rin_count = run
        .samples
        .iter()
        .filter_map(|s| s.rin)
        .filter(|r| (1.0..=10.0).contains(r))
        .count();
    assert!(rin_count > 0, "expected at least one RIN value in [1,10]");
}

#[test]
fn calibrates_dna1000_to_ladder_range() {
    let xml = load_gz_xml("demo_dna1000.xml.gz");
    let mut run = traceio::bioanalyzer::parse_xml(&xml).expect("parse");
    traceio::calibration::calculate_length(&mut run, traceio::calibration::Method::Hyman)
        .expect("calibrate");

    // Every sample gets a per-point length vector aligned with its trace.
    for s in &run.samples {
        assert_eq!(s.length.len(), s.fluorescence.len(), "length/trace mismatch");
    }

    // The ladder well should span the DNA 1000 kit's markers: 15 bp .. 1500 bp.
    let ladder = &run.samples[run.ladder_index().unwrap()];
    let finite: Vec<f64> = ladder.length.iter().copied().filter(|v| v.is_finite()).collect();
    assert!(finite.len() > 100, "calibrated region too small");
    let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!((lo - 15.0).abs() < 1.0, "lower marker should be ~15 bp, got {lo}");
    assert!((hi - 1500.0).abs() < 10.0, "upper marker should be ~1500 bp, got {hi}");

    // Calibrated length must increase monotonically along the trace where finite
    // (larger fragments migrate later), i.e. the mobility model is monotone.
    let mut prev = f64::NEG_INFINITY;
    for &l in &ladder.length {
        if l.is_finite() {
            assert!(l >= prev - 1e-6, "length not monotone along trace");
            prev = l;
        }
    }
}

/// End-to-end native `.xad` decode. Runs only when a real sample is dropped in
/// at `testdata/sample.xad` (none is committed). This is how we will validate
/// the [`traceio::xad`] container unwrap once a file is available.
#[test]
fn decodes_native_xad_if_present() {
    let path = testdata("sample.xad");
    if !path.exists() {
        eprintln!("skipping: no {} present", path.display());
        return;
    }
    let run = traceio::xad::read_xad_file(&path).expect("decode native .xad");
    assert!(!run.samples.is_empty(), "native .xad produced no samples");
    assert!(
        run.samples.iter().any(|s| !s.fluorescence.is_empty()),
        "native .xad produced no traces"
    );
}
