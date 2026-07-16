//! Validation of the reverse-engineered Fragment Analyzer reader against a real
//! run. The run directory holds private instrument data (`fa_examples/`, which
//! is git-ignored), so every test skips cleanly when it is not present.

use std::path::PathBuf;

/// The reference FA run directory, or `None` if the private data is absent.
fn fa_run_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fa_examples/16-03-27");
    dir.is_dir().then_some(dir)
}

#[test]
fn reads_fa_run_traces_and_names() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    assert_eq!(run.assay.length_unit, "bp");
    assert_eq!(run.samples.len(), 12, "expected 12 capillaries");
    for (i, s) in run.samples.iter().enumerate() {
        assert_eq!(s.well_number, (i + 1) as i32);
        assert!(!s.name.is_empty(), "sample {i} has no name");
        assert!(s.name.starts_with(&format!("D{}", i + 1)), "well label in name: {}", s.name);
        assert_eq!(s.time.len(), s.fluorescence.len());
        assert!(s.fluorescence.len() > 1000, "trace too short");
    }
}

#[test]
fn fa_size_calibration_places_main_peak() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    // Calibration must span the 1..6000 bp ladder.
    let s0 = &run.samples[0];
    let finite: Vec<f64> = s0.length.iter().copied().filter(|v| v.is_finite()).collect();
    assert!(finite.len() > 500, "too few calibrated points");
    let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!((lo - 1.0).abs() < 5.0, "ladder should start ~1 bp, got {lo}");
    assert!((hi - 6000.0).abs() < 50.0, "ladder should end ~6000 bp, got {hi}");

    // D1's main sample peak is at ~294 bp (per the vendor Peak Table). Find the
    // brightest calibrated point in the 150..1000 bp window and check it lands
    // there — this validates .raw extraction and .PKS calibration together.
    let mut best_bp = f64::NAN;
    let mut best_f = f32::NEG_INFINITY;
    for i in 0..s0.length.len() {
        let bp = s0.length[i];
        if bp.is_finite() && (150.0..1000.0).contains(&bp) && s0.fluorescence[i] > best_f {
            best_f = s0.fluorescence[i];
            best_bp = bp;
        }
    }
    assert!((best_bp - 294.0).abs() < 45.0, "D1 main peak should be ~294 bp, got {best_bp}");
}

#[test]
fn fa_peaks_and_markers_from_pks() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    // D1: lower marker (1 bp), a ~294 bp sample peak, upper marker (6000 bp).
    let d1 = &run.samples[0];
    assert_eq!(d1.peaks.len(), 3, "D1 should have LM + sample + UM");
    let lm = &d1.peaks[0];
    let um = d1.peaks.last().unwrap();
    assert_eq!(lm.observations, "Lower Marker");
    assert_eq!(um.observations, "Upper Marker");
    assert!((lm.length - 1.0).abs() < 0.5, "LM sized 1 bp");
    assert!((um.length - 6000.0).abs() < 0.5, "UM sized 6000 bp");
    let sample = &d1.peaks[1];
    assert!((sample.length - 294.0).abs() < 20.0, "D1 sample ~294 bp, got {}", sample.length);
    assert!((sample.area - 77.4).abs() < 2.0, "D1 sample area ~77.4, got {}", sample.area);

    // Exactly one well is the size ladder (D12), with the full 16-point ladder.
    let ladders: Vec<&traceio::Sample> = run.samples.iter().filter(|s| s.is_ladder).collect();
    assert_eq!(ladders.len(), 1, "one ladder well expected");
    assert_eq!(ladders[0].peaks.len(), 16, "ladder should have 16 peaks");
    assert!((ladders[0].peaks.last().unwrap().length - 6000.0).abs() < 0.5);
}
