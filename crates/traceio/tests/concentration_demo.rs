//! Validation of per-point concentration/molarity against a real demo run
//! (jwfoley/bioanalyzeR DNA 1000 export, MIT).

use std::io::Read;
use std::path::PathBuf;

use traceio::calibration::{calculate_length, Method};
use traceio::concentration::{calculate_concentration, calculate_molarity};

fn testdata(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata")
        .join(name)
}

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

fn sum_finite(xs: &[f64]) -> f64 {
    xs.iter().copied().filter(|v| v.is_finite()).sum()
}

#[test]
fn concentration_and_molarity_on_dna1000() {
    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let mut run = traceio::bioanalyzer::parse_xml(&xml).expect("parse");

    calculate_length(&mut run, Method::Hyman).expect("length");
    calculate_concentration(&mut run).expect("concentration");
    calculate_molarity(&mut run).expect("molarity");

    // Per-point vectors are aligned with each sample's trace.
    for s in &run.samples {
        assert_eq!(s.concentration.len(), s.fluorescence.len(), "conc len");
        assert_eq!(s.molarity.len(), s.fluorescence.len(), "molarity len");
    }

    // The ladder's total concentration (over finite points) is positive/finite.
    let ladder = &run.samples[run.ladder_index().unwrap()];
    let total_conc = sum_finite(&ladder.concentration);
    assert!(
        total_conc.is_finite() && total_conc > 0.0,
        "ladder total concentration should be positive, got {total_conc}"
    );

    // Individual points can be slightly negative (baseline-subtracted signal
    // dips below zero), so check totals, not per-point sign: every sample's
    // summed molarity is positive and finite, with many finite points.
    let mut finite_molarity = 0usize;
    for s in &run.samples {
        for &m in &s.molarity {
            if m.is_finite() {
                finite_molarity += 1;
            }
        }
        let total_m = sum_finite(&s.molarity);
        assert!(
            total_m.is_finite() && total_m > 0.0,
            "sample {} total molarity should be positive, got {total_m}",
            s.well_number
        );
    }
    assert!(
        finite_molarity > 100,
        "expected many finite molarity points"
    );

    // Molarity must be NaN wherever length is NaN (MW undefined there).
    for s in &run.samples {
        for (m, l) in s.molarity.iter().zip(s.length.iter()) {
            if l.is_nan() {
                assert!(m.is_nan(), "molarity should be NaN where length is NaN");
            }
        }
    }
}
