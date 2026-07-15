//! Headless plot render — writes the electropherogram of one or more samples to
//! a PPM, for display-free visual checking of the plotters rendering (single
//! trace, or an overlay of several with a legend).
//!
//! Usage: cargo run -p traceanalyzer --example render_png -- <file> [out.ppm] [indices] [y_mode] [normalize]
//!   indices:   single (e.g. 0) or comma-separated for an overlay (e.g. 0,1,2)
//!   y_mode:    fluorescence | concentration | molarity
//!   normalize: "normalize" to peak-height-normalize overlaid traces

use std::io::Write;
use std::path::PathBuf;

use traceanalyzer::plot::Series;
use traceanalyzer::{loading, plot};

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let path = PathBuf::from(args.next().expect("usage: render_png <file> [out.ppm] [indices] [ymode] [normalize]"));
    let out = args.next().unwrap_or_else(|| "plot.ppm".into());
    let indices: Vec<usize> = args
        .next()
        .map(|s| s.split(',').filter_map(|p| p.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![0]);
    let y_mode = match args.next().as_deref() {
        Some("concentration") => plot::YMode::Concentration,
        Some("molarity") => plot::YMode::Molarity,
        _ => plot::YMode::Fluorescence,
    };
    let normalize = matches!(args.next().as_deref(), Some("normalize"));

    let loaded = loading::load(&path)?;
    let run = loaded.run;
    let overlay = indices.len() > 1;

    let series: Vec<Series> = indices
        .iter()
        .map(|&idx| {
            let sample = run
                .samples
                .get(idx)
                .unwrap_or_else(|| panic!("sample {idx} out of range (have {})", run.samples.len()));
            let s = plot::series(&run, sample, y_mode, false);
            if overlay && normalize { plot::normalized(&s) } else { s }
        })
        .collect();
    let refs: Vec<&Series> = series.iter().collect();

    // Optional 6th arg: highlight a data-x position (cross-highlight check).
    let highlight = args.next().and_then(|s| s.parse::<f64>().ok());

    let vp = plot::auto_viewport_multi(&refs);
    let buf = plot::render_overlay(&refs, &vp, highlight, &[], plot::PLOT_W, plot::PLOT_H);

    let mut f = std::fs::File::create(&out)?;
    write!(f, "P6\n{} {}\n255\n", plot::PLOT_W, plot::PLOT_H)?;
    f.write_all(&buf)?;
    eprintln!(
        "wrote {out}  (indices {indices:?}, {} series, y={:?}, normalize={normalize})",
        series.len(),
        y_mode,
    );
    Ok(())
}
