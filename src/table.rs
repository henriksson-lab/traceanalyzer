//! Peak/region table rows for the Table tab. Builds display strings for the
//! focused sample's peaks (with computed height and % of total area) followed
//! by the run's smear regions. Each peak row carries its x-position so the plot
//! can cross-highlight it.

use crate::plot::{self, XAxis};
use crate::traceio::{Electrophoresis, Peak, Sample};

/// One table row: the display cell strings plus, for peaks, the plot x-position
/// (calibrated length or migration time) used for cross-highlighting.
pub struct Row {
    pub cells: Vec<String>,
    /// x-position of this peak in plot space, or `None` for region rows.
    pub peak_x: Option<f64>,
}

/// Column headers, in order. Kept in sync with the `.slint` TableView columns.
pub const HEADERS: [&str; 9] = [
    "#", "size", "time (s)", "area", "height", "% total", "conc", "molarity", "note",
];

/// Fluorescence height at a peak: the trace value nearest the peak's time.
fn peak_height(s: &Sample, p: &Peak) -> f64 {
    let Some(best) = plot::nearest_time_index(s, p.time) else {
        return f64::NAN;
    };
    s.fluorescence
        .get(best)
        .copied()
        .map(|v| v as f64)
        .unwrap_or(f64::NAN)
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
    rows_with_axis(run, s, plot::default_x_axis(s))
}

/// Build rows using an explicit x-axis for peak cross-highlighting.
pub fn rows_with_axis(run: &Electrophoresis, s: &Sample, x_axis: XAxis) -> Vec<Row> {
    let total_area: f64 = s
        .peaks
        .iter()
        .map(|p| p.area)
        .filter(|v| v.is_finite())
        .sum();

    let mut out: Vec<Row> = s
        .peaks
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let height = peak_height(s, p);
            let values = plot::peak_point_values(s, p);
            let pct = if total_area > 0.0 && p.area.is_finite() {
                p.area / total_area * 100.0
            } else {
                f64::NAN
            };
            let length = finite_or(values.length, p.length);
            let concentration = finite_or(values.concentration, p.concentration);
            let molarity = finite_or(values.molarity, p.molarity);
            let peak_x = match x_axis {
                XAxis::Time => p.time,
                XAxis::Length => length,
            };
            Row {
                cells: vec![
                    format!("{}", i + 1),
                    num(length, 0),
                    num(p.time, 1),
                    num(p.area, 1),
                    num(height, 1),
                    num(pct, 1),
                    num(concentration, 2),
                    num(molarity, 2),
                    p.observations.clone(),
                ],
                peak_x: peak_x.is_finite().then_some(peak_x),
            }
        })
        .collect();

    // Smear regions (size window only; other columns not applicable). Prefer the
    // sample's own regions (TapeStation) and fall back to the run-level ones
    // (Bioanalyzer keeps regions at the assay level).
    let regions = if s.regions.is_empty() {
        &run.regions
    } else {
        &s.regions
    };
    for (i, r) in regions.iter().enumerate() {
        out.push(Row {
            cells: vec![
                format!("R{}", i + 1),
                format!(
                    "{:.0}–{:.0} {}",
                    r.lower_length, r.upper_length, run.assay.length_unit
                ),
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

fn finite_or(primary: f64, fallback: f64) -> f64 {
    if primary.is_finite() {
        primary
    } else {
        fallback
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traceio::{AssayInfo, Peak};

    fn sample_with_stale_peak_fields() -> Sample {
        Sample {
            well_number: 1,
            name: "A1".to_string(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![1.0, 2.0, 3.0],
            fluorescence: vec![4.0, 8.0, 4.0],
            aligned_time: Vec::new(),
            length: vec![100.0, 222.0, 300.0],
            concentration: vec![1.0, 5.25, 3.0],
            molarity: vec![10.0, 42.5, 30.0],
            peaks: vec![Peak {
                observations: String::new(),
                length: 999.0,
                time: 2.0,
                aligned_time: f64::NAN,
                start_time: f64::NAN,
                end_time: f64::NAN,
                aligned_start_time: f64::NAN,
                aligned_end_time: f64::NAN,
                area: 50.0,
                concentration: 99.0,
                molarity: 999.0,
            }],
            regions: Vec::new(),
        }
    }

    fn run() -> Electrophoresis {
        Electrophoresis {
            assay: AssayInfo {
                file_name: String::new(),
                creation_date: String::new(),
                assay_name: String::new(),
                assay_type: "DNA".to_string(),
                length_unit: "bp".to_string(),
                concentration_unit: "ng/ul".to_string(),
                molarity_unit: Some("nM".to_string()),
                has_upper_marker: false,
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: Vec::new(),
        }
    }

    #[test]
    fn peak_rows_use_recalibrated_trace_values_at_peak_time() {
        let sample = sample_with_stale_peak_fields();
        let rows = rows_with_axis(&run(), &sample, XAxis::Length);

        assert_eq!(rows[0].cells[1], "222");
        assert_eq!(rows[0].cells[6], "5.25");
        assert_eq!(rows[0].cells[7], "42.50");
        assert_eq!(rows[0].peak_x, Some(222.0));
    }

    #[test]
    fn peak_row_x_position_can_follow_time_axis() {
        let sample = sample_with_stale_peak_fields();
        let rows = rows_with_axis(&run(), &sample, XAxis::Time);

        assert_eq!(rows[0].peak_x, Some(2.0));
    }

    #[test]
    fn peak_rows_fall_back_to_parsed_values_without_calibrated_arrays() {
        let mut sample = sample_with_stale_peak_fields();
        sample.length.clear();
        sample.concentration.clear();
        sample.molarity.clear();

        let rows = rows_with_axis(&run(), &sample, XAxis::Length);

        assert_eq!(rows[0].cells[1], "999");
        assert_eq!(rows[0].cells[6], "99.00");
        assert_eq!(rows[0].cells[7], "999.00");
        assert_eq!(rows[0].peak_x, Some(999.0));
    }
}
