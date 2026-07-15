//! Overview tab: a small-multiples grid of every sample's trace, rendered as a
//! single `plotters` bitmap. Cells are laid out row-major; [`OverviewLayout`]
//! maps a click back to a sample index so the Detail tab can open it.

use plotters::prelude::*;

use crate::plot::{ensure_font, series, YMode, PLOT_W};
use traceio::Electrophoresis;

/// Grid geometry for the overview bitmap.
#[derive(Debug, Clone, Copy)]
pub struct OverviewLayout {
    pub cols: usize,
    pub rows: usize,
    pub n: usize,
    pub w: u32,
    pub h: u32,
}

const CELL_H: u32 = 200;

/// Choose a grid for `n` samples: up to 4 columns, as many rows as needed.
pub fn layout(n: usize) -> OverviewLayout {
    let cols = n.clamp(1, 4);
    let rows = n.div_ceil(cols).max(1);
    OverviewLayout {
        cols,
        rows,
        n,
        w: PLOT_W,
        h: rows as u32 * CELL_H,
    }
}

/// Map a fractional (0..1) click on the overview image to a sample index.
pub fn cell_at(layout: &OverviewLayout, fx: f64, fy: f64) -> Option<usize> {
    if layout.n == 0 {
        return None;
    }
    let col = (fx * layout.cols as f64).floor() as usize;
    let row = (fy * layout.rows as f64).floor() as usize;
    let col = col.min(layout.cols - 1);
    let row = row.min(layout.rows - 1);
    let idx = row * layout.cols + col;
    (idx < layout.n).then_some(idx)
}

/// Render the small-multiples grid into an RGB buffer (`w*h*3` bytes). With
/// `shared_y`, every cell uses the global y-range (comparable heights);
/// otherwise each cell auto-fits its own trace.
pub fn render(
    run: &Electrophoresis,
    y_mode: YMode,
    shared_y: bool,
    layout: &OverviewLayout,
) -> Vec<u8> {
    let mut buf = vec![255u8; (layout.w * layout.h * 3) as usize];
    if let Err(e) = draw(&mut buf, run, y_mode, shared_y, layout) {
        eprintln!("(overview render failed: {e})");
    }
    buf
}

fn draw(
    buf: &mut [u8],
    run: &Electrophoresis,
    y_mode: YMode,
    shared_y: bool,
    layout: &OverviewLayout,
) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font();
    let trace = RGBColor(31, 119, 180);
    let ladder_c = RGBColor(0, 158, 115);
    let grid = RGBColor(238, 238, 238);
    let text = RGBColor(60, 60, 60);

    // Pre-compute each sample's series once.
    let all: Vec<_> = run
        .samples
        .iter()
        .map(|s| series(run, s, y_mode, false))
        .collect();

    // Global y-range for shared-scale mode.
    let global_y = {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for s in &all {
            for &y in &s.ys {
                lo = lo.min(y);
                hi = hi.max(y);
            }
        }
        if lo.is_finite() {
            (lo, hi)
        } else {
            (0.0, 1.0)
        }
    };

    let root = BitMapBackend::with_buffer(buf, (layout.w, layout.h)).into_drawing_area();
    root.fill(&WHITE)?;
    let cells = root.split_evenly((layout.rows, layout.cols));

    for (i, cell) in cells.into_iter().enumerate() {
        let Some(sample) = run.samples.get(i) else {
            break;
        };
        let s = &all[i];

        // x-range and per-plot y-range.
        let (mut x0, mut x1) = (f64::INFINITY, f64::NEG_INFINITY);
        for &x in &s.xs {
            x0 = x0.min(x);
            x1 = x1.max(x);
        }
        if !x0.is_finite() {
            x0 = 0.0;
            x1 = 1.0;
        }
        let (y0, y1) = if shared_y {
            global_y
        } else {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for &y in &s.ys {
                lo = lo.min(y);
                hi = hi.max(y);
            }
            if lo.is_finite() {
                (lo, hi)
            } else {
                (0.0, 1.0)
            }
        };
        let pad = ((y1 - y0) * 0.08).max(f64::EPSILON);
        let x1 = x1.max(x0 + f64::EPSILON);
        let y1 = y1 + pad;

        let title = if sample.is_ladder {
            format!("★ {}", s.name)
        } else {
            s.name.clone()
        };
        let color = if sample.is_ladder { ladder_c } else { trace };

        let mut chart = ChartBuilder::on(&cell)
            .caption(&title, ("sans-serif", 14).into_font().color(&text))
            .margin(6)
            .x_label_area_size(18)
            .y_label_area_size(28)
            .build_cartesian_2d(x0..x1, (y0 - pad)..y1)?;
        chart
            .configure_mesh()
            .light_line_style(grid)
            .bold_line_style(WHITE)
            .label_style(("sans-serif", 10).into_font().color(&text))
            .x_labels(4)
            .y_labels(3)
            .draw()?;
        chart.draw_series(LineSeries::new(
            s.xs.iter().zip(&s.ys).map(|(&x, &y)| (x, y)),
            color.stroke_width(1),
        ))?;
    }
    root.present()?;
    Ok(())
}
