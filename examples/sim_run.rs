//! Headless instrument simulator: replay a recorded `.pck` session through the
//! transport → protocol → run-state-machine stack and report the acquisition
//! it produces. This is the no-hardware `--simulate` path (milestone M3).
//!
//! Usage:
//!   cargo run --example sim_run -- <session.pck>

use anyhow::{Context, Result};
use traceanalyzer::tracehw::{run_to_completion, PckReplay};

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .context("usage: sim_run <session.pck>")?;

    let mut transport = PckReplay::from_path(&path)?;
    let assay = transport
        .header
        .assay_path
        .clone()
        .unwrap_or_else(|| "unknown assay".into());

    let acq = run_to_completion(&mut transport)?;

    println!("== simulated run: {path} ==");
    println!("assay:         {assay}");
    println!("commands sent: {:?}", transport.sent);
    println!("final state:   {:?}", acq.state);
    println!("acks:          {:?}", acq.acks);
    println!(
        "error polls:   {} ({})",
        acq.errors.len(),
        acq.errors.first().map(String::as_str).unwrap_or("-")
    );
    println!(
        "samples:       stream A {}, stream B {}",
        acq.stream_a.len(),
        acq.stream_b.len()
    );
    println!("corrupt recs:  {}", acq.corrupt);

    // Provisional trace (uncalibrated until Phase P2).
    let ep = acq.to_electrophoresis(&assay, 1.0);
    let s = &ep.samples[0];
    let range = s
        .fluorescence_range()
        .map(|(lo, hi)| format!("{lo:.0}..{hi:.0}"))
        .unwrap_or_else(|| "empty".into());
    println!(
        "electropherogram: {} points, provisional signal range {}",
        s.fluorescence.len(),
        range
    );
    Ok(())
}
