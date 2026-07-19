//! Validation of the TapeStation reader against a real Agilent demo export
//! (a D1000 run from the MIT-licensed bioanalyzeR package). The fixtures are
//! fetched by `scripts/fetch-testdata.sh` into `testdata/tapestation/` and are
//! git-ignored, so the test skips cleanly when they are absent.

use std::path::PathBuf;

/// The D1000 demo XML (its `_Electropherogram.csv` sibling sits alongside), or
/// `None` when the fixtures have not been fetched.
fn d1000_xml() -> Option<PathBuf> {
    let p =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../testdata/tapestation/d1000.xml.gz");
    p.is_file().then_some(p)
}

#[test]
fn reads_d1000_export_with_traces_and_sizing() {
    let Some(xml) = d1000_xml() else {
        eprintln!("skipping: testdata/tapestation not fetched (run scripts/fetch-testdata.sh)");
        return;
    };
    let run = traceio::tapestation::read_tapestation(&xml).expect("read D1000 export");

    assert_eq!(run.assay.length_unit, "bp");
    assert_eq!(run.assay.assay_type, "DNA");
    assert!(run.assay.has_upper_marker);
    assert_eq!(run.samples.len(), 16, "D1000-Tubes-16 has 16 lanes");

    // Exactly one ladder, carrying the D1000 bands (25 bp .. 1500 bp markers).
    let ladders: Vec<&traceio::Sample> = run.samples.iter().filter(|s| s.is_ladder).collect();
    assert_eq!(ladders.len(), 1);
    let ladder = ladders[0];
    assert_eq!(ladder.name, "Ladder");
    assert_eq!(ladder.peaks.first().unwrap().observations, "Lower Marker");
    assert_eq!(ladder.peaks.last().unwrap().observations, "Upper Marker");
    assert!((ladder.peaks.first().unwrap().length - 25.0).abs() < 0.5);
    assert!((ladder.peaks.last().unwrap().length - 1500.0).abs() < 0.5);

    // Traces attached from the CSV and size-calibrated across the ladder span.
    assert!(ladder.fluorescence.len() > 500, "trace attached");
    assert_eq!(ladder.time.len(), ladder.fluorescence.len());
    let finite: Vec<f64> = ladder
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    assert!(finite.len() > 300, "sized across the ladder");
    let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (lo - 25.0).abs() < 5.0,
        "calibration starts ~25 bp, got {lo}"
    );
    assert!(
        (hi - 1500.0).abs() < 10.0,
        "calibration ends ~1500 bp, got {hi}"
    );

    // A "300bp fragment" sample sizes its main peak near 300 bp (from the XML).
    let frag = run
        .samples
        .iter()
        .find(|s| s.name.contains("300bp fragment"))
        .expect("a 300bp fragment lane");
    let main = frag
        .peaks
        .iter()
        .filter(|p| p.observations.is_empty()) // exclude the markers
        .max_by(|a, b| a.area.total_cmp(&b.area))
        .expect("a sample peak");
    assert!(
        (main.length - 300.0).abs() < 30.0,
        "300bp fragment sized {}",
        main.length
    );
    assert!(main.concentration.is_finite() && main.concentration > 0.0);
}
