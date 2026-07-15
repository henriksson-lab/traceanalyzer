//! Virtual-gel view: each sample becomes a vertical lane whose darkness encodes
//! fluorescence along the common migration-time axis (well at top, molecules run
//! downward), mimicking a slab-gel image. Rendered as a single `plotters`
//! bitmap so it can share the Slint `Image` display path.

use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};

use crate::plot::{ensure_font, PLOT_H, PLOT_W};
use traceio::Electrophoresis;

// Chart insets (px), matching the ChartBuilder configuration in `draw`; used to
// map a fractional click back to a lane.
const LEFT: f64 = 66.0; // margin(10) + y_label_area(56)
const RIGHT: f64 = 10.0; // margin(10)

/// Which samples to draw and the overall image size for a gel render.
pub fn size(_n: usize) -> (u32, u32) {
    (PLOT_W, PLOT_H)
}

/// Original sample indices that become gel lanes (those with a trace), in order.
pub fn lane_indices(run: &Electrophoresis) -> Vec<usize> {
    run.samples
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.fluorescence.is_empty())
        .map(|(i, _)| i)
        .collect()
}

/// Map a fractional (0..1) x-click on the gel image to an original sample index.
pub fn lane_at(run: &Electrophoresis, fx: f64, w: u32) -> Option<usize> {
    let lanes = lane_indices(run);
    if lanes.is_empty() {
        return None;
    }
    let left = LEFT / w as f64;
    let span = (w as f64 - LEFT - RIGHT) / w as f64;
    let t = (fx - left) / span;
    if !(0.0..=1.0).contains(&t) {
        return None;
    }
    let pos = (t * lanes.len() as f64).floor() as usize;
    lanes.get(pos.min(lanes.len() - 1)).copied()
}

/// Render the virtual gel into an RGB buffer (`w*h*3` bytes). Lanes share one
/// intensity scale (global fluorescence max) so band darkness is comparable.
pub fn render(run: &Electrophoresis, w: u32, h: u32) -> Vec<u8> {
    let mut buf = vec![255u8; (w * h * 3) as usize];
    if let Err(e) = draw(&mut buf, run, w, h) {
        eprintln!("(gel render failed: {e})");
    }
    buf
}

fn draw(buf: &mut [u8], run: &Electrophoresis, w: u32, h: u32) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font();
    let text = RGBColor(40, 40, 90);

    let samples: Vec<_> = run.samples.iter().filter(|s| !s.fluorescence.is_empty()).collect();
    let n = samples.len();

    let root = BitMapBackend::with_buffer(buf, (w, h)).into_drawing_area();
    root.fill(&WHITE)?;
    if n == 0 {
        return Ok(());
    }

    // Common time axis and global intensity scale across all lanes.
    let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut f_max = f64::EPSILON;
    for s in &samples {
        for &t in &s.time {
            t_min = t_min.min(t);
            t_max = t_max.max(t);
        }
        for &f in &s.fluorescence {
            if (f as f64).is_finite() {
                f_max = f_max.max(f as f64);
            }
        }
    }
    if !t_min.is_finite() {
        return Ok(());
    }

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .margin_top(28) // room for rotated lane labels
        .x_label_area_size(0)
        .y_label_area_size(56)
        // x = lane index [0, n]; y = time, inverted so the well (min time) is up.
        .build_cartesian_2d(0f64..n as f64, t_max..t_min)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .y_desc("migration time (s)")
        .axis_desc_style(("sans-serif", 15).into_font().color(&text))
        .label_style(("sans-serif", 12).into_font().color(&text))
        .draw()?;

    // Each lane: a stack of grayscale rectangles between successive trace points.
    for (i, s) in samples.iter().enumerate() {
        let x0 = i as f64 + 0.06;
        let x1 = i as f64 + 0.94;
        let m = s.time.len().min(s.fluorescence.len());
        for k in 0..m.saturating_sub(1) {
            let (ta, tb) = (s.time[k], s.time[k + 1]);
            if !ta.is_finite() || !tb.is_finite() {
                continue;
            }
            let f = s.fluorescence[k] as f64;
            let norm = (f / f_max).clamp(0.0, 1.0);
            let g = (255.0 * (1.0 - norm)).round() as u8;
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x0, ta), (x1, tb)],
                RGBColor(g, g, g).filled(),
            )))?;
        }
        // Lane label (rotated) above the lane.
        let name = if s.name.is_empty() { format!("Well {}", s.well_number) } else { s.name.clone() };
        let style = ("sans-serif", 12)
            .into_font()
            .color(&text)
            .transform(FontTransform::Rotate270)
            .pos(Pos::new(HPos::Center, VPos::Bottom));
        root.draw(&Text::new(
            name,
            chart.backend_coord(&(i as f64 + 0.5, t_min)),
            style,
        ))?;
    }

    root.present()?;
    Ok(())
}
