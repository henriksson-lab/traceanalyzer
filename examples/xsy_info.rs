//! Dump the contents of a 2100 Expert `.xsy` assay file.
//!
//! Usage: cargo run --example xsy_info -- <assay.xsy>

use anyhow::{Context, Result};

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .context("usage: xsy_info <assay.xsy>")?;
    let assay = traceanalyzer::traceio::xsy::read_xsy_file(&path)?;

    println!("name:        {}", assay.name);
    if !assay.comment.is_empty() {
        println!("comment:     {}", assay.comment);
    }
    println!(
        "variant:     {}",
        if assay.compressed {
            "compressed (Xceed)"
        } else {
            "plain"
        }
    );
    println!("size unit:   {}", assay.size_unit);
    if let Some(c) = assay.ladder_concentration {
        println!("ladder conc: {c}");
    }
    println!("ladder peaks: {}", assay.ladder_peaks.len());
    for p in &assay.ladder_peaks {
        println!(
            "  size {:>7} {}  areaB {:>8.2}  conc {:>6.2}",
            p.size, assay.size_unit, p.area_b, p.concentration
        );
    }
    println!("script values: {}", assay.script_values.len());
    if !assay.script_values.is_empty() {
        let head: Vec<String> = assay
            .script_values
            .iter()
            .take(16)
            .map(|v| v.to_string())
            .collect();
        println!("  first: {}", head.join(", "));
    }
    if !assay.sample_names.is_empty() {
        println!("preview samples: {}", assay.sample_names.join(", "));
    }
    Ok(())
}
