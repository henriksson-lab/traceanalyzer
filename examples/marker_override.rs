//! Verify manual marker overrides: shift a sample's lower marker and confirm
//! its per-point sizing (and hence length range) changes accordingly.
//!
//! Usage: cargo run --example marker_override -- <file> [idx]

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{bail, ensure};
use traceanalyzer::traceio::calibration::{self, MarkerOverride};

fn length_range(s: &traceanalyzer::traceio::Sample) -> (f64, f64) {
    let fin: Vec<f64> = s.length.iter().copied().filter(|v| v.is_finite()).collect();
    (
        fin.iter().cloned().fold(f64::INFINITY, f64::min),
        fin.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    )
}

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: marker_override <file> [idx]"));
    let idx: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(1);

    let mut run = traceanalyzer::loading::load(&path)?.run;
    if run.samples.is_empty() {
        bail!("{} contains no samples", path.display());
    }
    ensure!(
        idx < run.samples.len(),
        "sample index {idx} out of range; {} has {} samples (valid indices: 0..{})",
        path.display(),
        run.samples.len(),
        run.samples.len() - 1
    );

    let (lo, up) = calibration::marker_times(&run, idx, None);
    println!("sample {idx}: detected markers lower={lo:?} upper={up:?}");
    let (a0, a1) = length_range(&run.samples[idx]);
    println!("auto length range:       {a0:.1} .. {a1:.1}");

    // Shift the lower marker later by 2 s and recompute.
    let mut overrides = HashMap::new();
    overrides.insert(
        idx,
        MarkerOverride {
            lower_time: lo.map(|t| t + 2.0),
            upper_time: None,
        },
    );
    traceanalyzer::loading::recalibrate_with(&mut run, &overrides)?;
    let (b0, b1) = length_range(&run.samples[idx]);
    println!("override (+2s) range:    {b0:.1} .. {b1:.1}");

    assert!(
        (a0 - b0).abs() > 1e-6 || (a1 - b1).abs() > 1e-6,
        "override had no effect"
    );
    println!("OK: override changed the sizing");
    Ok(())
}
