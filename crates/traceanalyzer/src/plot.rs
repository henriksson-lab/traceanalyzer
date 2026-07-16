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
use traceio::{Electrophoresis, Peak, Sample};

/// Register the vendored font once, so `plotters` text works with the
/// `ab_glyph` backend (which bundles no system-font access).
pub(crate) fn ensure_font() {
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
    /// All modes in dropdown order; index into this matches [`YMode::index`].
    pub const ALL: [YMode; 3] = [YMode::Fluorescence, YMode::Concentration, YMode::Molarity];

    pub fn next(self) -> YMode {
        match self {
            YMode::Fluorescence => YMode::Concentration,
            YMode::Concentration => YMode::Molarity,
            YMode::Molarity => YMode::Fluorescence,
        }
    }

    /// Position of this mode in [`YMode::ALL`].
    pub fn index(self) -> usize {
        match self {
            YMode::Fluorescence => 0,
            YMode::Concentration => 1,
            YMode::Molarity => 2,
        }
    }

    /// Mode for a dropdown index; out-of-range falls back to fluorescence.
    pub fn from_index(i: usize) -> YMode {
        YMode::ALL.get(i).copied().unwrap_or(YMode::Fluorescence)
    }
    /// True when the active run has data that can produce non-empty series for
    /// this mode. Fluorescence is the fallback because every supported run can
    /// render it; derived arrays are optional in some native formats.
    pub fn is_available(self, run: &Electrophoresis) -> bool {
        match self {
            YMode::Fluorescence => true,
            YMode::Concentration => has_finite_values(run, |s| &s.concentration),
            YMode::Molarity => has_finite_values(run, |s| &s.molarity),
        }
    }
    pub fn available_for(run: &Electrophoresis) -> Vec<YMode> {
        YMode::ALL
            .iter()
            .copied()
            .filter(|m| m.is_available(run))
            .collect()
    }
    pub fn from_available_index(run: &Electrophoresis, i: usize) -> YMode {
        YMode::available_for(run)
            .get(i)
            .copied()
            .unwrap_or(YMode::Fluorescence)
    }
    pub fn available_index(self, run: &Electrophoresis) -> usize {
        YMode::available_for(run)
            .iter()
            .position(|m| *m == self)
            .unwrap_or(0)
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

fn has_finite_values(run: &Electrophoresis, values: impl Fn(&Sample) -> &[f64]) -> bool {
    run.samples
        .iter()
        .any(|s| values(s).iter().any(|v| v.is_finite()))
}

/// A data window in (x, y) space.
#[derive(Debug, Clone, Copy)]
pub struct Viewport {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

/// Quantity plotted on the x-axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XAxis {
    Time,
    Length,
}

impl XAxis {
    pub fn label(self, run: &Electrophoresis) -> String {
        match self {
            XAxis::Time => "migration time (s)".to_string(),
            XAxis::Length => format!("size ({})", run.assay.length_unit),
        }
    }
}

/// Accepted x-axis selectors for [`series`]. This keeps the older bool call
/// sites working: `true` forces time, `false` chooses the sample default.
pub enum XAxisSelection {
    Default,
    Axis(XAxis),
}

impl From<XAxis> for XAxisSelection {
    fn from(axis: XAxis) -> Self {
        XAxisSelection::Axis(axis)
    }
}

impl From<bool> for XAxisSelection {
    fn from(force_time: bool) -> Self {
        if force_time {
            XAxisSelection::Axis(XAxis::Time)
        } else {
            XAxisSelection::Default
        }
    }
}

/// Extracted plot data for one sample.
pub struct Series {
    /// Legend label (sample/channel name); shown when overlaying.
    pub name: String,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub peaks_x: Vec<f64>,
    pub x_label: String,
    pub y_label: String,
}

/// Okabe–Ito colorblind-safe qualitative palette (yellow last: it is the
/// weakest as a thin line on white, so it is only reached when overlaying many
/// traces). Used to color overlaid samples and their legend entries.
pub const PALETTE: [(u8, u8, u8); 8] = [
    (0, 114, 178),   // blue
    (213, 94, 0),    // vermillion
    (0, 158, 115),   // bluish green
    (204, 121, 167), // reddish purple
    (86, 180, 233),  // sky blue
    (230, 159, 0),   // orange
    (0, 0, 0),       // black
    (204, 187, 0),   // dark yellow
];

/// Color for the i-th overlaid series.
pub fn palette_color(i: usize) -> RGBColor {
    let (r, g, b) = PALETTE[i % PALETTE.len()];
    RGBColor(r, g, b)
}

/// Whether a sample should be plotted against calibrated length (bp/nt).
/// Public so the table can label peak positions with the same x quantity.
pub fn use_length(s: &Sample) -> bool {
    s.length.iter().filter(|v| v.is_finite()).count() >= 2
}

pub fn default_x_axis(s: &Sample) -> XAxis {
    if use_length(s) {
        XAxis::Length
    } else {
        XAxis::Time
    }
}

/// Linearly interpolate a per-point trace array at a migration time.
pub fn value_at_time(times: &[f64], values: &[f64], target: f64) -> Option<f64> {
    if !target.is_finite() {
        return None;
    }
    let mut prev: Option<(f64, f64)> = None;
    for (&t, &v) in times.iter().zip(values) {
        if !t.is_finite() || !v.is_finite() {
            prev = None;
            continue;
        }
        if (t - target).abs() <= f64::EPSILON {
            return Some(v);
        }
        if let Some((pt, pv)) = prev {
            let between = (pt <= target && target <= t) || (t <= target && target <= pt);
            if between {
                let span = t - pt;
                if span.abs() <= f64::EPSILON {
                    return Some(v);
                }
                let f = (target - pt) / span;
                return Some(pv + f * (v - pv));
            }
        }
        prev = Some((t, v));
    }
    None
}

/// Index of the trace point nearest a reported peak/migration time.
pub fn nearest_time_index(s: &Sample, time: f64) -> Option<usize> {
    if s.time.is_empty() || !time.is_finite() {
        return None;
    }
    let mut best = None;
    let mut best_d = f64::INFINITY;
    for (i, &t) in s.time.iter().enumerate() {
        let d = (t - time).abs();
        if d < best_d {
            best_d = d;
            best = Some(i);
        }
    }
    best
}

/// Peak values derived from the calibrated per-point arrays at peak time.
#[derive(Debug, Clone, Copy)]
pub struct PeakPointValues {
    pub time: f64,
    pub length: f64,
    pub concentration: f64,
    pub molarity: f64,
}

impl PeakPointValues {
    pub fn x(self, axis: XAxis) -> f64 {
        match axis {
            XAxis::Time => self.time,
            XAxis::Length => self.length,
        }
    }
}

pub fn peak_point_values(s: &Sample, p: &Peak) -> PeakPointValues {
    let at = |xs: &[f64]| -> f64 { value_at_time(&s.time, xs, p.time).unwrap_or(f64::NAN) };
    PeakPointValues {
        time: p.time,
        length: at(&s.length),
        concentration: at(&s.concentration),
        molarity: at(&s.molarity),
    }
}

/// Build the (x, y) series for a sample under the given y-mode. When
/// `x_axis` is [`XAxis::Time`] the x-axis is raw migration time regardless of
/// calibration (used by marker-edit mode).
pub fn series(
    run: &Electrophoresis,
    s: &Sample,
    y_mode: YMode,
    x_axis: impl Into<XAxisSelection>,
) -> Series {
    let x_axis = match x_axis.into() {
        XAxisSelection::Default => default_x_axis(s),
        XAxisSelection::Axis(axis) => axis,
    };
    let x_at = |i: usize| -> f64 {
        match x_axis {
            XAxis::Time => s.time.get(i).copied().unwrap_or(f64::NAN),
            XAxis::Length => s.length.get(i).copied().unwrap_or(f64::NAN),
        }
    };
    let y_at = |i: usize| -> f64 {
        match y_mode {
            YMode::Fluorescence => s
                .fluorescence
                .get(i)
                .copied()
                .map(|v| v as f64)
                .unwrap_or(f64::NAN),
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
        .map(|p| peak_point_values(s, p).x(x_axis))
        .filter(|v| v.is_finite())
        .collect();

    let name = if s.name.is_empty() {
        format!("Well {}", s.well_number)
    } else {
        s.name.clone()
    };
    Series {
        name,
        xs,
        ys,
        peaks_x,
        x_label: x_axis.label(run),
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
        name: format!("{} (raw)", ch.channel_id),
        xs,
        ys,
        peaks_x: Vec::new(),
        x_label: "time (s)".to_string(),
        y_label: format!("{} (raw)", ch.channel_id),
    }
}

/// Return a copy of `s` with y-values scaled so its finite maximum is 1.0
/// (peak-height normalization for shape comparison across overlaid traces).
/// A non-positive or non-finite max leaves the data unchanged.
pub fn normalized(s: &Series) -> Series {
    let max =
        s.ys.iter()
            .copied()
            .filter(|v| v.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
    let scale = if max.is_finite() && max > 0.0 {
        1.0 / max
    } else {
        1.0
    };
    Series {
        name: s.name.clone(),
        xs: s.xs.clone(),
        ys: s.ys.iter().map(|v| v * scale).collect(),
        peaks_x: s.peaks_x.clone(),
        x_label: s.x_label.clone(),
        y_label: "normalized".to_string(),
    }
}

/// Auto-fit viewport around a series (with a little y-headroom).
pub fn auto_viewport(series: &Series) -> Viewport {
    auto_viewport_multi(&[series])
}

/// Auto-fit viewport covering the union of several series' extents.
pub fn auto_viewport_multi(series: &[&Series]) -> Viewport {
    viewport_impl(series, false)
}

/// Like [`auto_viewport_multi`], but the upper y-bound ignores an extreme narrow
/// spike so it doesn't squash the rest of the trace. Use for derived quantities
/// (concentration, molarity) where a tiny-size point can blow up numerically:
/// e.g. molarity = concentration / molecular_weight explodes as length → 0.
/// Broad real peaks (max within ~4× the 98th percentile) are never clipped.
pub fn auto_viewport_multi_robust(series: &[&Series]) -> Viewport {
    viewport_impl(series, true)
}

fn viewport_impl(series: &[&Series], robust: bool) -> Viewport {
    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut ys: Vec<f64> = Vec::new();
    for s in series {
        for &x in &s.xs {
            x_min = x_min.min(x);
            x_max = x_max.max(x);
        }
        for &y in &s.ys {
            y_min = y_min.min(y);
            y_max = y_max.max(y);
            if robust && y.is_finite() {
                ys.push(y);
            }
        }
    }
    if !x_min.is_finite() {
        x_min = 0.0;
        x_max = 1.0;
    }
    if !y_min.is_finite() {
        y_min = 0.0;
        y_max = 1.0;
    }
    if robust {
        y_max = robust_y_max(&mut ys, y_max);
    }
    let pad = ((y_max - y_min) * 0.05).max(f64::EPSILON);
    Viewport {
        x_min,
        x_max,
        y_min: y_min - pad,
        y_max: y_max + pad,
    }
}

/// Upper y-bound that discards an extreme narrow spike. A spike is detected when
/// the maximum exceeds `SPIKE_FACTOR` × the 98th percentile (broad peaks sit well
/// under that); the bound then becomes the largest value below the spike, so the
/// spike runs off the top of the plot while the rest of the trace fills it.
fn robust_y_max(ys: &mut [f64], plain_max: f64) -> f64 {
    const SPIKE_FACTOR: f64 = 4.0;
    // Too few points to characterize a distribution — keep the true max.
    if ys.len() < 20 {
        return plain_max;
    }
    ys.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p98 = ys[((ys.len() - 1) as f64 * 0.98).round() as usize];
    if p98 <= 0.0 || plain_max <= SPIKE_FACTOR * p98 {
        return plain_max; // no dominating spike
    }
    let threshold = SPIKE_FACTOR * p98;
    // Largest value at or below the spike threshold (ys is sorted ascending).
    ys.iter()
        .rev()
        .copied()
        .find(|&y| y <= threshold)
        .unwrap_or(plain_max)
}

/// Render one series into an RGB buffer (`w*h*3` bytes).
pub fn render_rgb(series: &Series, vp: &Viewport, w: u32, h: u32) -> Vec<u8> {
    render_overlay(&[series], vp, None, &[], w, h)
}

/// Render one or more series into an RGB buffer (`w*h*3` bytes). A single series
/// is drawn in the primary color with its peak markers; multiple series are
/// overlaid in distinct palette colors with a legend and no peak markers.
/// `highlight_x`, if set, draws a bold marker at that x (table→plot highlight).
/// `marker_xs` draws draggable marker guide lines (marker-edit mode).
pub fn render_overlay(
    series: &[&Series],
    vp: &Viewport,
    highlight_x: Option<f64>,
    marker_xs: &[f64],
    w: u32,
    h: u32,
) -> Vec<u8> {
    let mut buf = vec![255u8; (w * h * 3) as usize];
    if let Err(e) = draw_into(&mut buf, series, vp, highlight_x, marker_xs, w, h) {
        eprintln!("(plot render failed: {e})");
    }
    buf
}

fn draw_into(
    buf: &mut [u8],
    series: &[&Series],
    vp: &Viewport,
    highlight_x: Option<f64>,
    marker_xs: &[f64],
    w: u32,
    h: u32,
) -> Result<(), Box<dyn std::error::Error>> {
    ensure_font();
    let peak = RGBColor(224, 130, 20);
    let grid = RGBColor(232, 232, 232);
    let text = RGBColor(60, 60, 60);
    let overlay = series.len() > 1;

    let root = BitMapBackend::with_buffer(buf, (w, h)).into_drawing_area();
    root.fill(&WHITE)?;
    // Guard against a degenerate range.
    let (x0, x1) = (vp.x_min, vp.x_max.max(vp.x_min + f64::EPSILON));
    let (y0, y1) = (vp.y_min, vp.y_max.max(vp.y_min + f64::EPSILON));

    // Axis labels come from the first series (all overlaid series share units).
    let (x_label, y_label) = match series.first() {
        Some(s) => (s.x_label.as_str(), s.y_label.as_str()),
        None => ("", ""),
    };

    let mut chart = ChartBuilder::on(&root)
        .margin(MARGIN as i32)
        .x_label_area_size(X_LABEL_AREA as i32)
        .y_label_area_size(Y_LABEL_AREA as i32)
        .build_cartesian_2d(x0..x1, y0..y1)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .axis_desc_style(("sans-serif", 16).into_font().color(&text))
        .label_style(("sans-serif", 13).into_font().color(&text))
        .light_line_style(grid)
        .draw()?;

    // Peak markers (behind the trace); only for a single-sample view.
    if !overlay {
        for s in series {
            for &px in &s.peaks_x {
                if px >= x0 && px <= x1 {
                    chart.draw_series(std::iter::once(PathElement::new(
                        vec![(px, y0), (px, y1)],
                        peak.mix(0.45),
                    )))?;
                }
            }
        }
    }

    // Traces, one palette color each.
    for (i, s) in series.iter().enumerate() {
        let color = if overlay {
            palette_color(i)
        } else {
            RGBColor(31, 119, 180)
        };
        let line = chart.draw_series(LineSeries::new(
            s.xs.iter().zip(&s.ys).map(|(&x, &y)| (x, y)),
            color.stroke_width(2),
        ))?;
        if overlay {
            line.label(&s.name).legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 18, y)], color.stroke_width(3))
            });
        }
    }

    // Cross-highlight marker for a selected table row.
    if let Some(hx) = highlight_x {
        if hx >= x0 && hx <= x1 {
            let hl = RGBColor(214, 39, 40);
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(hx, y0), (hx, y1)],
                hl.stroke_width(2),
            )))?;
        }
    }

    // Draggable marker guide lines (marker-edit mode), in green.
    for &mx in marker_xs {
        if mx >= x0 && mx <= x1 {
            let mc = RGBColor(0, 158, 115);
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(mx, y0), (mx, y1)],
                mc.stroke_width(2),
            )))?;
        }
    }

    if overlay {
        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperRight)
            .background_style(WHITE.mix(0.85))
            .border_style(RGBColor(200, 200, 200))
            .label_font(("sans-serif", 14).into_font().color(&text))
            .draw()?;
    }
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

/// Inverse of [`frac_to_data`] on the x-axis: map a data x to a fractional
/// widget x (0..1). Used to hit-test draggable marker lines.
pub fn data_x_to_frac(data_x: f64, vp: &Viewport) -> f64 {
    let left = (MARGIN + Y_LABEL_AREA) / PLOT_W as f64;
    let right = (PLOT_W as f64 - MARGIN) / PLOT_W as f64;
    let span = vp.x_max - vp.x_min;
    let tx = if span.abs() > f64::EPSILON {
        (data_x - vp.x_min) / span
    } else {
        0.0
    };
    left + tx * (right - left)
}

#[cfg(test)]
mod tests {
    use super::*;
    use traceio::{AssayInfo, Electrophoresis, Sample};

    fn series_with_ys(ys: Vec<f64>) -> Series {
        let xs: Vec<f64> = (0..ys.len()).map(|i| i as f64).collect();
        Series {
            name: String::new(),
            xs,
            ys,
            peaks_x: Vec::new(),
            x_label: String::new(),
            y_label: String::new(),
        }
    }

    fn sample_with_derived(concentration: Vec<f64>, molarity: Vec<f64>) -> Sample {
        Sample {
            well_number: 1,
            name: "A1".to_string(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![1.0],
            fluorescence: vec![1.0],
            aligned_time: Vec::new(),
            length: vec![100.0],
            concentration,
            molarity,
            peaks: Vec::new(),
        }
    }

    fn run_with_sample(sample: Sample) -> Electrophoresis {
        Electrophoresis {
            assay: AssayInfo {
                concentration_unit: "ng/ul".to_string(),
                molarity_unit: Some("nM".to_string()),
                ..Default::default()
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![sample],
        }
    }

    #[test]
    fn y_mode_availability_omits_empty_derived_arrays() {
        let run = run_with_sample(sample_with_derived(Vec::new(), Vec::new()));

        assert_eq!(YMode::available_for(&run), vec![YMode::Fluorescence]);
        assert_eq!(YMode::from_available_index(&run, 1), YMode::Fluorescence);
        assert_eq!(YMode::Molarity.available_index(&run), 0);
    }

    #[test]
    fn y_mode_availability_preserves_present_derived_arrays() {
        let run = run_with_sample(sample_with_derived(vec![0.0], vec![f64::NAN, 2.0]));

        assert_eq!(
            YMode::available_for(&run),
            vec![YMode::Fluorescence, YMode::Concentration, YMode::Molarity]
        );
        assert_eq!(YMode::from_available_index(&run, 2), YMode::Molarity);
        assert_eq!(YMode::Molarity.available_index(&run), 2);
    }

    #[test]
    fn robust_viewport_clips_a_narrow_spike() {
        // 200 points near ~1.0, plus one huge spike: the marker-molarity case.
        let mut ys = vec![1.0; 200];
        ys[0] = 500.0;
        let s = series_with_ys(ys);
        let robust = auto_viewport_multi_robust(&[&s]);
        let plain = auto_viewport_multi(&[&s]);
        assert!(
            plain.y_max > 400.0,
            "plain keeps the spike: {}",
            plain.y_max
        );
        assert!(
            robust.y_max < 10.0,
            "robust drops the spike so the bulk fills the plot: {}",
            robust.y_max
        );
    }

    #[test]
    fn robust_viewport_preserves_a_broad_peak() {
        // A broad peak (many points near the top) is real signal, not a spike:
        // robust scaling must not clip it. Baseline 0, a wide plateau at 100.
        let mut ys = vec![0.0; 200];
        for y in ys.iter_mut().take(80) {
            *y = 100.0;
        }
        let s = series_with_ys(ys);
        let robust = auto_viewport_multi_robust(&[&s]);
        assert!(robust.y_max >= 100.0, "broad peak kept: {}", robust.y_max);
    }

    #[test]
    fn robust_viewport_noop_with_few_points() {
        let s = series_with_ys(vec![1.0, 2.0, 50.0]);
        let robust = auto_viewport_multi_robust(&[&s]);
        assert!(robust.y_max > 49.0, "too few points to judge: keep max");
    }
}
