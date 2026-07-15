//! Bitmap plot rendering with `plotters` (axes, gridlines, labels), returned as
//! an RGB buffer for display in a Slint `Image`.
//!
//! The rendered image is shown with `image-fit: fill`, so the widget stretches
//! it 1:1. Interaction code therefore works in **fractional** widget
//! coordinates and maps them to data via [`frac_to_data`], keeping overlays and
//! the bitmap in registration regardless of widget size.

use std::sync::OnceLock;

use plotters::prelude::*;
use plotters::style::register_font;
use traceio::{Electrophoresis, Sample};

/// Register the vendored font once, so `plotters` text works with the
/// `ab_glyph` backend (which bundles no system-font access).
fn ensure_font() {
    static FONT: OnceLock<()> = OnceLock::new();
    FONT.get_or_init(|| {
        let _ = register_font(
            "sans-serif",
            FontStyle::Normal,
            include_bytes!("../assets/DejaVuSans.ttf"),
        );
    });
}

/// Render resolution of the plot bitmap.
pub const PLOT_W: u32 = 1600;
pub const PLOT_H: u32 = 560;

// Chart insets (pixels at render resolution); also used as fractions for
// pixel→data mapping. Must match the ChartBuilder configuration below.
const MARGIN: f64 = 12.0;
const X_LABEL_AREA: f64 = 42.0;
const Y_LABEL_AREA: f64 = 66.0;

/// Quantity plotted on the y-axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum YMode {
    Fluorescence,
    Concentration,
    Molarity,
}

impl YMode {
    pub fn next(self) -> YMode {
        match self {
            YMode::Fluorescence => YMode::Concentration,
            YMode::Concentration => YMode::Molarity,
            YMode::Molarity => YMode::Fluorescence,
        }
    }
    pub fn label(self, run: &Electrophoresis) -> String {
        match self {
            YMode::Fluorescence => "fluorescence".into(),
            YMode::Concentration => format!("concentration ({})", run.assay.concentration_unit),
            YMode::Molarity => format!(
                "molarity ({})",
                run.assay.molarity_unit.as_deref().unwrap_or("")
            ),
        }
    }
}

/// A data window in (x, y) space.
#[derive(Debug, Clone, Copy)]
pub struct Viewport {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

/// Extracted plot data for one sample.
pub struct Series {
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub peaks_x: Vec<f64>,
    pub x_label: String,
    pub y_label: String,
}

/// Whether a sample should be plotted against calibrated length (bp/nt).
fn use_length(s: &Sample) -> bool {
    s.length.iter().filter(|v| v.is_finite()).count() >= 2
}

/// Build the (x, y) series for a sample under the given y-mode.
pub fn series(run: &Electrophoresis, s: &Sample, y_mode: YMode) -> Series {
    let use_len = use_length(s);
    let x_at = |i: usize| -> f64 {
        if use_len {
            s.length.get(i).copied().unwrap_or(f64::NAN)
        } else {
            s.time.get(i).copied().unwrap_or(f64::NAN)
        }
    };
    let y_at = |i: usize| -> f64 {
        match y_mode {
            YMode::Fluorescence => s.fluorescence.get(i).copied().map(|v| v as f64).unwrap_or(f64::NAN),
            YMode::Concentration => s.concentration.get(i).copied().unwrap_or(f64::NAN),
            YMode::Molarity => s.molarity.get(i).copied().unwrap_or(f64::NAN),
        }
    };

    let n = s.fluorescence.len();
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for i in 0..n {
        let (x, y) = (x_at(i), y_at(i));
        if x.is_finite() && y.is_finite() {
            xs.push(x);
            ys.push(y);
        }
    }
    let peaks_x = s
        .peaks
        .iter()
        .map(|p| if use_len { p.length } else { p.time })
        .filter(|v| v.is_finite())
        .collect();

    let x_label = if use_len {
        format!("size ({})", run.assay.length_unit)
    } else {
        "migration time (s)".to_string()
    };
    Series {
        xs,
        ys,
        peaks_x,
        x_label,
        y_label: y_mode.label(run),
    }
}

/// Build a series for a raw detector channel (time vs raw signal, no peaks).
pub fn raw_series(ch: &traceio::xad::RawChannel) -> Series {
    let mut xs = Vec::with_capacity(ch.signal.len());
    let mut ys = Vec::with_capacity(ch.signal.len());
    for (i, &v) in ch.signal.iter().enumerate() {
        let y = v as f64;
        if y.is_finite() {
            xs.push(ch.x_start + ch.x_step * i as f64);
            ys.push(y);
        }
    }
    Series {
        xs,
        ys,
        peaks_x: Vec::new(),
        x_label: "time (s)".to_string(),
        y_label: format!("{} (raw)", ch.channel_id),
    }
}

/// Auto-fit viewport around a series (with a little y-headroom).
pub fn auto_viewport(series: &Series) -> Viewport {
    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &x in &series.xs {
        x_min = x_min.min(x);
        x_max = x_max.max(x);
    }
    for &y in &series.ys {
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }
    if !x_min.is_finite() {
        x_min = 0.0;
        x_max = 1.0;
    }
    if !y_min.is_finite() {
        y_min = 0.0;
        y_max = 1.0;
    }
    let pad = ((y_max - y_min) * 0.05).max(f64::EPSILON);
    Viewport {
        x_min,
        x_max,
        y_min: y_min - pad,
        y_max: y_max + pad,
    }
}

/// Render one series into an RGB buffer (`w*h*3` bytes).
pub fn render_rgb(series: &Series, vp: &Viewport, w: u32, h: u32) -> Vec<u8> {
    let mut buf = vec![255u8; (w * h * 3) as usize];
    if let Err(e) = draw_into(&mut buf, series, vp, w, h) {
        eprintln!("(plot render failed: {e})");
    }
    buf
}

fn draw_into(
    buf: &mut [u8],
    series: &Series,
    vp: &Viewport,
    w: u32,
    h: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font();
    let trace = RGBColor(31, 119, 180);
    let peak = RGBColor(224, 130, 20);
    let grid = RGBColor(232, 232, 232);
    let text = RGBColor(60, 60, 60);

    let root = BitMapBackend::with_buffer(buf, (w, h)).into_drawing_area();
    root.fill(&WHITE)?;
    // Guard against a degenerate range.
    let (x0, x1) = (vp.x_min, vp.x_max.max(vp.x_min + f64::EPSILON));
    let (y0, y1) = (vp.y_min, vp.y_max.max(vp.y_min + f64::EPSILON));

    let mut chart = ChartBuilder::on(&root)
        .margin(MARGIN as i32)
        .x_label_area_size(X_LABEL_AREA as i32)
        .y_label_area_size(Y_LABEL_AREA as i32)
        .build_cartesian_2d(x0..x1, y0..y1)?;

    chart
        .configure_mesh()
        .x_desc(&series.x_label)
        .y_desc(&series.y_label)
        .axis_desc_style(("sans-serif", 16).into_font().color(&text))
        .label_style(("sans-serif", 13).into_font().color(&text))
        .light_line_style(grid)
        .draw()?;

    // Peak markers (behind the trace).
    for &px in &series.peaks_x {
        if px >= x0 && px <= x1 {
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(px, y0), (px, y1)],
                peak.mix(0.45),
            )))?;
        }
    }
    // Trace.
    chart.draw_series(LineSeries::new(
        series.xs.iter().zip(&series.ys).map(|(&x, &y)| (x, y)),
        trace.stroke_width(2),
    ))?;
    root.present()?;
    Ok(())
}

/// Map a fractional widget coordinate (0..1, origin top-left) to data space,
/// accounting for the chart insets. Used by zoom/pan and marker interactions.
pub fn frac_to_data(fx: f64, fy: f64, vp: &Viewport) -> (f64, f64) {
    let left = (MARGIN + Y_LABEL_AREA) / PLOT_W as f64;
    let right = (PLOT_W as f64 - MARGIN) / PLOT_W as f64;
    let top = MARGIN / PLOT_H as f64;
    let bottom = (PLOT_H as f64 - MARGIN - X_LABEL_AREA) / PLOT_H as f64;

    let tx = ((fx - left) / (right - left)).clamp(0.0, 1.0);
    let ty = ((fy - top) / (bottom - top)).clamp(0.0, 1.0);
    let data_x = vp.x_min + tx * (vp.x_max - vp.x_min);
    let data_y = vp.y_max - ty * (vp.y_max - vp.y_min); // y inverted
    (data_x, data_y)
}
