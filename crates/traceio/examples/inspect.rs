//! Headless inspector: parse a Bioanalyzer file and print a summary.
//!
//! Usage: cargo run -p traceio --example inspect -- <file.xad | file.xml | file.xml.gz>

use std::io::Read;

fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).ok_or_else(|| {
        anyhow::anyhow!("usage: inspect <file.xad|file.xml|file.xml.gz>")
    })?;
    let mut run = load(&path)?;
    // Per-point ladder calibration (length in bp/nt for every trace point).
    if let Err(e) = traceio::calibration::calculate_length(
        &mut run,
        traceio::calibration::Method::Hyman,
    ) {
        eprintln!("(calibration skipped: {e})");
    }

    println!("File:        {}", run.assay.file_name);
    println!("Assay:       {} ({})", run.assay.assay_name, run.assay.assay_type);
    println!(
        "Units:       length={} conc={} molarity={:?}",
        run.assay.length_unit, run.assay.concentration_unit, run.assay.molarity_unit
    );
    println!("Ladder pks:  {}", run.ladder_peaks.len());
    println!("Regions:     {}", run.regions.len());
    println!("Samples:     {}", run.samples.len());
    println!();
    println!(
        "{:>4}  {:<24} {:>7} {:>6} {:>6}  {:>16}  {:>5}",
        "well", "name", "points", "peaks", "RIN", "size range", "ladder"
    );
    for s in &run.samples {
        println!(
            "{:>4}  {:<24} {:>7} {:>6} {:>6}  {:>16}  {:>5}",
            s.well_number,
            truncate(&s.name, 24),
            s.fluorescence.len(),
            s.peaks.len(),
            s.rin.map(|r| format!("{r:.1}")).unwrap_or_else(|| "-".into()),
            size_range(s, &run.assay.length_unit),
            if s.is_ladder { "yes" } else { "" },
        );
    }

    // For native .xad, also show the raw detector channels (the actual signal;
    // per-well processed traces are recomputed by 2100 Expert, see docs).
    if path.ends_with(".xad") {
        match traceio::xad::read_xad_raw_channels(std::path::Path::new(&path)) {
            Ok(chs) if !chs.is_empty() => {
                println!("\nRaw detector channels (whole-chip acquisition):");
                for c in &chs {
                    let secs = c.x_step * c.signal.len() as f64;
                    println!(
                        "  {:<18} {:>7} samples  @ {}s step  (~{:.0}s run)",
                        c.channel_id, c.signal.len(), c.x_step, secs
                    );
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn load(path: &str) -> anyhow::Result<traceio::Electrophoresis> {
    let p = std::path::Path::new(path);
    if path.ends_with(".xad") {
        traceio::xad::read_xad_file(p)
    } else if path.ends_with(".xml.gz") {
        let raw = std::fs::read(p)?;
        let mut d = flate2::read::GzDecoder::new(&raw[..]);
        let mut s = String::new();
        d.read_to_string(&mut s)?;
        traceio::bioanalyzer::parse_xml(&s)
    } else {
        let s = std::fs::read_to_string(p)?;
        traceio::bioanalyzer::parse_xml(&s)
    }
}

/// Min/max of the calibrated (finite) per-point length, e.g. "15–1500 bp".
fn size_range(s: &traceio::Sample, unit: &str) -> String {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &l in &s.length {
        if l.is_finite() {
            lo = lo.min(l);
            hi = hi.max(l);
        }
    }
    if lo.is_finite() {
        format!("{lo:.0}–{hi:.0} {unit}")
    } else {
        "-".into()
    }
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        s.chars().take(n - 1).collect::<String>() + "…"
    }
}
