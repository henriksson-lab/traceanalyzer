//! Headless virtual-gel render, for display-free visual checking.
//!
//! Usage: cargo run --example render_gel -- <file> [out.ppm]

use std::io::Write;
use std::path::PathBuf;

use traceanalyzer::{gel, loading};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: render_gel <file> [out.ppm]"));
    let out = args.next().unwrap_or_else(|| "gel.ppm".into());

    let run = loading::load(&path)?.run;
    let (w, h) = gel::size(run.samples.len());
    let buf = gel::render(&run, w, h);

    let mut f = std::fs::File::create(&out)?;
    write!(f, "P6\n{w} {h}\n255\n")?;
    f.write_all(&buf)?;
    eprintln!("wrote {out}  ({} samples, {w}x{h}px)", run.samples.len());
    Ok(())
}
