//! Print the peak/region table for one sample, for checking the Table tab data.
//!
//! Usage: cargo run -p traceanalyzer --example dump_table -- <file> [sample_index]

use std::path::PathBuf;

use traceanalyzer::{loading, table};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: dump_table <file> [idx]"));
    let idx: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(0);

    let run = loading::load(&path)?.run;
    let sample = &run.samples[idx];
    println!(
        "sample {idx}: {}  ({} peaks)",
        sample.name,
        sample.peaks.len()
    );
    println!("{}", table::HEADERS.join(" | "));
    for row in table::rows(&run, sample) {
        println!("{}", row.cells.join(" | "));
    }
    Ok(())
}
