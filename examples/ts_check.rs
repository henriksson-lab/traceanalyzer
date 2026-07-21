//! Load a TapeStation export (metadata `.xml` or `_Electropherogram.csv`) and
//! report what was decoded.
//!
//! Usage: cargo run --example ts_check -- <file.xml | _Electropherogram.csv>

use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let path = PathBuf::from(std::env::args().nth(1).expect("usage: ts_check <xml|csv>"));
    let run = traceanalyzer::traceio::tapestation::read_tapestation(&path)?;

    println!(
        "assay: {} ({})  units: {} / {} / {}",
        run.assay.assay_name,
        run.assay.assay_type,
        run.assay.length_unit,
        run.assay.concentration_unit,
        run.assay.molarity_unit.as_deref().unwrap_or("-"),
    );
    println!("file: {}", run.assay.file_name);
    println!("samples: {}", run.samples.len());
    for s in &run.samples {
        let n = s.fluorescence.len();
        let finite: Vec<f64> = s.length.iter().copied().filter(|v| v.is_finite()).collect();
        let (lo, hi) = if finite.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            (
                finite.iter().cloned().fold(f64::INFINITY, f64::min),
                finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            )
        };
        let rin = s.rin.map(|v| format!(" RIN {v:.1}")).unwrap_or_default();
        let ladder = if s.is_ladder { " [ladder]" } else { "" };
        println!(
            "  well {:>3}  {:20}  trace {n}pts  size {lo:.0}..{hi:.0} {}  {} peaks{rin}{ladder}",
            s.well_number,
            s.name,
            run.assay.length_unit,
            s.peaks.len()
        );
        for p in &s.peaks {
            let tag = if p.observations.is_empty() {
                "peak"
            } else {
                &p.observations
            };
            println!(
                "        {tag:14} {:.0} {}  conc {:.2}  molarity {:.1}  area {:.3}",
                p.length, run.assay.length_unit, p.concentration, p.molarity, p.area
            );
        }
    }
    Ok(())
}
