//! Headless overview-grid render, for display-free visual checking.
//!
//! Usage: cargo run -p traceanalyzer --example render_overview -- <file> [out.ppm] [shared|per]

use std::io::Write;
use std::path::PathBuf;

use traceanalyzer::{loading, overview, plot};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: render_overview <file> [out.ppm] [shared|per]"));
    let out = args.next().unwrap_or_else(|| "overview.ppm".into());
    let shared = !matches!(args.next().as_deref(), Some("per"));

    let run = loading::load(&path)?.run;
    let layout = overview::layout(run.samples.len());
    let buf = overview::render(&run, plot::YMode::Fluorescence, shared, &layout);

    let mut f = std::fs::File::create(&out)?;
    write!(f, "P6\n{} {}\n255\n", layout.w, layout.h)?;
    f.write_all(&buf)?;
    eprintln!(
        "wrote {out}  ({} samples, {}x{} grid, {}x{}px, shared_y={shared})",
        layout.n, layout.rows, layout.cols, layout.w, layout.h
    );
    Ok(())
}
