//! Peak/region table rows for the Table tab. Builds display strings for the
//! focused sample's peaks (with computed height and % of total area) followed
//! by the run's smear regions. Each peak row carries its x-position so the plot
//! can cross-highlight it.

use crate::plot::use_length;
use traceio::{Electrophoresis, Peak, Sample};

/// One table row: the display cell strings plus, for peaks, the plot x-position
/// (calibrated length or migration time) used for cross-highlighting.
pub struct Row {
    pub cells: Vec<String>,
    /// x-position of this peak in plot space, or `None` for region rows.
    pub peak_x: Option<f64>,
}

/// Column headers, in order. Kept in sync with the `.slint` TableView columns.
pub const HEADERS: [&str; 9] =
    ["#", "size", "time (s)", "area", "height", "% total", "conc", "molarity", "note"];

/// Fluorescence height at a peak: the trace value nearest the peak's time.
fn peak_height(s: &Sample, p: &Peak) -> f64 {
    if s.time.is_empty() || !p.time.is_finite() {
        return f64::NAN;
    }
    // Nearest time index (traces are short; a linear scan is fine).
    let mut best = 0usize;
    let mut best_d = f64::INFINITY;
    for (i, &t) in s.time.iter().enumerate() {
        let d = (t - p.time).abs();
        if d < best_d {
            best_d = d;
            best = i;
        }
    }
    s.fluorescence.get(best).copied().map(|v| v as f64).unwrap_or(f64::NAN)
}

fn num(v: f64, decimals: usize) -> String {
    if v.is_finite() {
        format!("{v:.*}", decimals)
    } else {
        "–".to_string()
    }
}

/// Build the rows for one sample's peaks followed by the run's regions.
pub fn rows(run: &Electrophoresis, s: &Sample) -> Vec<Row> {
    let use_len = use_length(s);
    let total_area: f64 = s.peaks.iter().map(|p| p.area).filter(|v| v.is_finite()).sum();

    let mut out: Vec<Row> = s
        .peaks
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let height = peak_height(s, p);
            let pct = if total_area > 0.0 && p.area.is_finite() {
                p.area / total_area * 100.0
            } else {
                f64::NAN
            };
            let peak_x = if use_len { p.length } else { p.time };
            Row {
                cells: vec![
                    format!("{}", i + 1),
                    num(p.length, 0),
                    num(p.time, 1),
                    num(p.area, 1),
                    num(height, 1),
                    num(pct, 1),
                    num(p.concentration, 2),
                    num(p.molarity, 2),
                    p.observations.clone(),
                ],
                peak_x: peak_x.is_finite().then_some(peak_x),
            }
        })
        .collect();

    // Smear regions (size window only; other columns not applicable).
    for (i, r) in run.regions.iter().enumerate() {
        out.push(Row {
            cells: vec![
                format!("R{}", i + 1),
                format!("{:.0}–{:.0} {}", r.lower_length, r.upper_length, run.assay.length_unit),
                "–".into(),
                "–".into(),
                "–".into(),
                "–".into(),
                "–".into(),
                "–".into(),
                "region".into(),
            ],
            peak_x: None,
        });
    }
    out
}
