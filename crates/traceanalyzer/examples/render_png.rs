//! Headless plot render — writes the electropherogram of one sample to a PPM,
//! for display-free visual checking of the plotters rendering.
//!
//! Usage: cargo run -p traceanalyzer --example render_png -- <file> [out.ppm] [sample_index] [y_mode]
//!   y_mode: fluorescence | concentration | molarity

use std::io::Write;
use std::path::PathBuf;

use traceanalyzer::{loading, plot};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: render_png <file> [out.ppm] [idx] [ymode]"));
    let out = args.next().unwrap_or_else(|| "plot.ppm".into());
    let idx: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(0);
    let y_mode = match args.next().as_deref() {
        Some("concentration") => plot::YMode::Concentration,
        Some("molarity") => plot::YMode::Molarity,
        _ => plot::YMode::Fluorescence,
    };

    let run = loading::load_and_calibrate(&path)?;
    let sample = run
        .samples
        .get(idx)
        .unwrap_or_else(|| panic!("sample {idx} out of range (have {})", run.samples.len()));

    let series = plot::series(&run, sample, y_mode);
    let vp = plot::auto_viewport(&series);
    let buf = plot::render_rgb(&series, &vp, plot::PLOT_W, plot::PLOT_H);

    let mut f = std::fs::File::create(&out)?;
    write!(f, "P6\n{} {}\n255\n", plot::PLOT_W, plot::PLOT_H)?;
    f.write_all(&buf)?;
    eprintln!(
        "wrote {out}  (sample {idx}: {}, {} points, {} peaks)",
        sample.name,
        series.xs.len(),
        series.peaks_x.len()
    );
    Ok(())
}
