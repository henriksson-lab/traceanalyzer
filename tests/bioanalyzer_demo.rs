//! Validation against real bioanalyzeR demo exports (jwfoley/bioanalyzeR, MIT).
//! The demo `.xml.gz` files are processed Bioanalyzer export XML. Native `.xad`
//! raw acquisition coverage is separate below and runs only when local `.xad`
//! fixtures are present.

use std::io::Read;
use std::path::PathBuf;

fn testdata(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
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

#[test]
fn parses_dna1000_demo() {
    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let run = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse DNA 1000");

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
    assert!(
        run.ladder_index().is_some(),
        "expected exactly one ladder well"
    );

    let ladder = &run.samples[run.ladder_index().unwrap()];
    assert!(ladder.is_ladder);
    assert!(
        ladder.time.len() == ladder.fluorescence.len() && !ladder.fluorescence.is_empty(),
        "ladder trace time/fluorescence must be non-empty and equal length"
    );
    assert!(!ladder.peaks.is_empty(), "ladder should have called peaks");

    // Trace fluorescence must be finite and vary (not all zero).
    let (lo, hi) = ladder.fluorescence_range().expect("range");
    assert!(
        lo.is_finite() && hi.is_finite() && hi > lo,
        "trace is flat/NaN"
    );

    // A DNA assay carries no RIN.
    assert!(run.samples.iter().all(|s| s.rin.is_none()));
}

#[test]
fn parses_rna_nano_demo_with_rin() {
    let Some(xml) = load_gz_xml("demo_rna_nano.xml.gz") else {
        return;
    };
    let run = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse RNA Nano");

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
    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let mut run = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse");
    traceanalyzer::traceio::calibration::calculate_length(
        &mut run,
        traceanalyzer::traceio::calibration::Method::Hyman,
    )
    .expect("calibrate");

    // Every sample gets a per-point length vector aligned with its trace.
    for s in &run.samples {
        assert_eq!(
            s.length.len(),
            s.fluorescence.len(),
            "length/trace mismatch"
        );
    }

    // The ladder well should span the DNA 1000 kit's markers: 15 bp .. 1500 bp.
    let ladder = &run.samples[run.ladder_index().unwrap()];
    let finite: Vec<f64> = ladder
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    assert!(finite.len() > 100, "calibrated region too small");
    let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (lo - 15.0).abs() < 1.0,
        "lower marker should be ~15 bp, got {lo}"
    );
    assert!(
        (hi - 1500.0).abs() < 10.0,
        "upper marker should be ~1500 bp, got {hi}"
    );

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

#[test]
fn marker_override_shifts_sizing() {
    use std::collections::HashMap;
    use traceanalyzer::traceio::calibration::{self, MarkerOverride, Method};

    let Some(xml) = load_gz_xml("demo_dna1000.xml.gz") else {
        return;
    };
    let base = traceanalyzer::traceio::bioanalyzer::parse_xml(&xml).expect("parse");

    // Pick a non-ladder sample and read its detected lower marker time.
    let idx = base
        .samples
        .iter()
        .position(|s| !s.is_ladder && !s.peaks.is_empty())
        .expect("a non-ladder sample with peaks");
    let (lower, _upper) = calibration::marker_times(&base, idx, None);
    let lower = lower.expect("detected lower marker");

    // Overriding with the *detected* time reproduces the automatic sizing.
    let mut same = base.clone();
    let mut ov_same = HashMap::new();
    ov_same.insert(
        idx,
        MarkerOverride {
            lower_time: Some(lower),
            upper_time: None,
        },
    );
    calibration::calculate_length_with(&mut same, Method::Hyman, &ov_same).expect("recalc");
    let mut auto = base.clone();
    calibration::calculate_length(&mut auto, Method::Hyman).expect("auto");
    let a: Vec<f64> = auto.samples[idx]
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    let b: Vec<f64> = same.samples[idx]
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    assert_eq!(
        a.len(),
        b.len(),
        "override-with-detected changed the calibrated span"
    );

    // Shifting the lower marker later must change the sizing measurably.
    let mut shifted = base.clone();
    let mut ov = HashMap::new();
    ov.insert(
        idx,
        MarkerOverride {
            lower_time: Some(lower + 3.0),
            upper_time: None,
        },
    );
    calibration::calculate_length_with(&mut shifted, Method::Hyman, &ov).expect("recalc shifted");
    let c: Vec<f64> = shifted.samples[idx]
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    let hi_auto = a.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let hi_shift = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (hi_auto - hi_shift).abs() > 1e-6,
        "marker shift did not change sizing"
    );

    // Other samples (no override) are unaffected by the per-sample override.
    let other = (0..base.samples.len()).find(|&i| i != idx && !base.samples[i].peaks.is_empty());
    if let Some(o) = other {
        let auto_o: Vec<f64> = auto.samples[o]
            .length
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .collect();
        let shift_o: Vec<f64> = shifted.samples[o]
            .length
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .collect();
        assert_eq!(auto_o.len(), shift_o.len(), "unrelated sample changed");
    }
}

/// Enumerate real `.xad` files to validate against: `testdata/sample.xad` plus
/// anything under a local `bioa_examples/` dir (private lab data, gitignored).
fn xad_samples() -> Vec<PathBuf> {
    let mut v = Vec::new();
    let one = testdata("sample.xad");
    if one.exists() {
        v.push(one);
    }
    let examples = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("bioa_examples");
    if let Ok(rd) = std::fs::read_dir(&examples) {
        for e in rd.flatten() {
            let p = e.path();
            if p.extension().map(|x| x == "xad").unwrap_or(false) {
                v.push(p);
            }
        }
    }
    v
}

/// End-to-end native `.xad` decode. Runs only when real `.xad` files are present
/// locally (none are committed). Validates the [`traceanalyzer::traceio::xad`] container unwrap.
///
/// Current `.xad` files carry raw acquisition + metadata, so we assert the
/// container decodes and yields assay info, the defined ladder, and sample
/// metadata — plus non-empty raw detector channels — rather than processed
/// per-well traces (see docs/xad_format.md §5).
#[test]
fn decodes_native_xad_if_present() {
    let samples = xad_samples();
    if samples.is_empty() {
        eprintln!("skipping: no real .xad files present");
        return;
    }
    for path in samples {
        let run = traceanalyzer::traceio::xad::read_xad_file(&path)
            .unwrap_or_else(|e| panic!("decode {}: {e}", path.display()));
        assert!(
            !run.samples.is_empty(),
            "no sample metadata in {}",
            path.display()
        );
        assert!(
            !run.assay.length_unit.is_empty(),
            "no assay info in {}",
            path.display()
        );
        assert!(
            !run.ladder_peaks.is_empty(),
            "no defined ladder in {}",
            path.display()
        );

        // Raw detector channels must be present and populated.
        let channels =
            traceanalyzer::traceio::xad::read_xad_raw_channels(&path).expect("raw channels");
        assert!(
            !channels.is_empty(),
            "no raw detector channels in {}",
            path.display()
        );
        assert!(
            channels.iter().any(|c| c.signal.len() > 1000),
            "raw detector channels are empty in {}",
            path.display()
        );
    }
}
