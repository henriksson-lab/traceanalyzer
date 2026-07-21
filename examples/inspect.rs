//! Headless inspector: parse any supported electrophoresis file and print a summary.
//!
//! Usage: cargo run --example inspect -- <path>

fn main() -> anyhow::Result<()> {
    let path = std::env::args()
        .nth(1)
        .ok_or_else(|| anyhow::anyhow!("usage: inspect <path>"))?;
    let loaded = traceanalyzer::traceio::io::read_path(&path)?;
    for warning in &loaded.warnings {
        eprintln!("({warning})");
    }
    let run = &loaded.run;

    println!("Path:        {}", loaded.source.path.display());
    println!("Identity:    {}", loaded.source.identity.display());
    println!("Format:      {:?}", loaded.source.format);
    println!("Save:        {:?}", loaded.save_capabilities());
    println!("File:        {}", run.assay.file_name);
    println!(
        "Assay:       {} ({})",
        run.assay.assay_name, run.assay.assay_type
    );
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
            s.rin
                .map(|r| format!("{r:.1}"))
                .unwrap_or_else(|| "-".into()),
            size_range(s, &run.assay.length_unit),
            if s.is_ladder { "yes" } else { "" },
        );
    }

    if !loaded.raw_channels.is_empty() {
        println!("\nRaw detector channels (whole-chip acquisition):");
        for c in &loaded.raw_channels {
            let secs = c.x_step * c.signal.len() as f64;
            println!(
                "  {:<18} {:>7} samples  @ {}s step  (~{:.0}s run)",
                c.channel_id,
                c.signal.len(),
                c.x_step,
                secs
            );
        }
    }
    Ok(())
}

/// Min/max of the calibrated (finite) per-point length, e.g. "15–1500 bp".
fn size_range(s: &traceanalyzer::traceio::Sample, unit: &str) -> String {
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
