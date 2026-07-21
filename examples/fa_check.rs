//! Load a Fragment Analyzer run via the reverse-engineered reader and report
//! what was decoded (samples, names, calibrated size range, trace stats).
//!
//! Usage: cargo run --example fa_check -- <run.raw | run_dir>

use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let path = PathBuf::from(
        std::env::args()
            .nth(1)
            .expect("usage: fa_check <run.raw|dir>"),
    );
    let run = traceanalyzer::traceio::fa::read_fa_run(&path)?;

    println!("assay: {} ({})", run.assay.assay_name, run.assay.assay_type);
    println!("file:  {}", run.assay.file_name);
    println!("samples: {}", run.samples.len());
    for s in &run.samples {
        let n = s.fluorescence.len();
        let fmin = s.fluorescence.iter().cloned().fold(f32::INFINITY, f32::min);
        let fmax = s
            .fluorescence
            .iter()
            .cloned()
            .fold(f32::NEG_INFINITY, f32::max);
        let finite: Vec<f64> = s.length.iter().copied().filter(|v| v.is_finite()).collect();
        let (lo, hi) = if finite.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            (
                finite.iter().cloned().fold(f64::INFINITY, f64::min),
                finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        };
        let ladder = if s.is_ladder { " [ladder]" } else { "" };
        println!(
            "  well {:>3}  {:24}  n={n}  fluor {fmin:.0}..{fmax:.0}  size {lo:.0}..{hi:.0} bp  {} peaks{ladder}",
            s.well_number, s.name, s.peaks.len()
        );
        for p in &s.peaks {
            let tag = if p.observations.is_empty() {
                "peak"
            } else {
                &p.observations
            };
            println!(
                "        {tag:14} {:.0} bp  (apex {:.0}, area {:.1})",
                p.length, p.time, p.area
            );
        }
    }
    Ok(())
}
