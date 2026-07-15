//! Render a sample in marker-edit style: raw-time x-axis with the lower/upper
//! marker guide lines, for display-free checking of Phase H rendering.
//!
//! Usage: cargo run -p traceanalyzer --example render_markers -- <file> [out.ppm] [idx]

use std::io::Write;
use std::path::PathBuf;

use traceanalyzer::{loading, plot};
use traceio::calibration;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: render_markers <file> [out.ppm] [idx]"));
    let out = args.next().unwrap_or_else(|| "markers.ppm".into());
    let idx: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(1);

    let run = loading::load(&path)?.run;
    let s = &run.samples[idx];

    // Raw-time series (force_time = true) + effective marker lines.
    let series = plot::series(&run, s, plot::YMode::Fluorescence, true);
    let (lo, up) = calibration::marker_times(&run, idx, None);
    let markers: Vec<f64> = lo.into_iter().chain(up).collect();
    println!("sample {idx}: marker lines at raw times {markers:?}");

    let vp = plot::auto_viewport(&series);
    let buf = plot::render_overlay(&[&series], &vp, None, &markers, plot::PLOT_W, plot::PLOT_H);

    let mut f = std::fs::File::create(&out)?;
    write!(f, "P6\n{} {}\n255\n", plot::PLOT_W, plot::PLOT_H)?;
    f.write_all(&buf)?;
    eprintln!("wrote {out}");
    Ok(())
}
