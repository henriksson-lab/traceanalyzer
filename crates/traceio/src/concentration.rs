//! Per-point concentration and molarity, ported from jwfoley/bioanalyzeR
//! (`calculate.concentration`, `calculate.molarity`, `molecular.weight`,
//! `integrate.peak`, MIT). Run these AFTER [`crate::calibration::calculate_length`],
//! which fills each sample's `aligned_time` and `length`.
//!
//! For the Bioanalyzer the integration x-variable is `aligned.time` (matching
//! bioanalyzeR's `get.x.name`). The mass coefficient is fit from the ladder's
//! non-marker peaks and then rescaled per sample by the ratio of the ladder's
//! marker area(s) to the sample's marker area(s), compensating for differences
//! in loading/detection.

use crate::model::{Electrophoresis, Sample};
use anyhow::{anyhow, Result};

/// Marker peak observation labels (from bioanalyzeR `electrophoresis.R`).
const LOWER_MARKER_NAMES: [&str; 2] = ["Lower Marker", "edited Lower Marker"];
const UPPER_MARKER_NAMES: [&str; 2] = ["Upper Marker", "edited Upper Marker"];

/// Estimated molecular weight (daltons) of a nucleic acid of the given length.
///
/// From bioanalyzeR `molecular.weight` (Thermo Fisher DNA/RNA MW conversions):
/// dsDNA = length·607.4 + 157.9, ssRNA = length·320.5 + 159.0. Unknown assay
/// types yield NaN.
pub fn molecular_weight(length: f64, assay_type: &str) -> f64 {
    if !length.is_finite() || length <= 0.0 {
        return f64::NAN;
    }
    match assay_type {
        "DNA" => length * 607.4 + 157.9,
        "RNA" => length * 320.5 + 159.0,
        _ => f64::NAN,
    }
}

/// Fill `concentration` for every sample of a Bioanalyzer run.
///
/// Requires exactly one ladder well and that `calculate_length` has already
/// run (so `aligned_time` is populated). Concentration is computed for every
/// point from its trapezoidal area; the first point of each sample is NaN.
pub fn calculate_concentration(run: &mut Electrophoresis) -> Result<()> {
    let ladder_idx = run
        .ladder_index()
        .ok_or_else(|| anyhow!("need exactly one ladder well for concentration"))?;

    // Defined ladder concentrations (the true values, from the assay setpoints).
    let ladder_concs: Vec<f64> = run.ladder_peaks.iter().map(|p| p.concentration).collect();
    if ladder_concs.is_empty() {
        return Err(anyhow!("no defined ladder concentrations"));
    }
    let lower_conc = *ladder_concs.first().unwrap();
    let upper_conc = *ladder_concs.last().unwrap();

    // Per-point trapezoidal area for every sample (first point of each = NaN).
    let areas: Vec<Vec<f64>> = run.samples.iter().map(per_point_area).collect();

    // --- ladder mass coefficient (single ladder shared by all samples) ------
    let ladder = &run.samples[ladder_idx];
    let ladder_area = &areas[ladder_idx];

    // Ladder peaks = those whose *reported* concentration is a defined ladder
    // concentration; markers are the lower/upper of those.
    let ladder_markers = marker_indices(ladder, lower_conc, upper_conc);
    let mut non_marker_ratios = Vec::new();
    for (i, p) in ladder.peaks.iter().enumerate() {
        if !conc_in_set(p.concentration, &ladder_concs) {
            continue; // not a recognized ladder peak
        }
        if ladder_markers.contains(&i) {
            continue; // markers are excluded from the mass-coefficient fit
        }
        let area = integrate_peak_area(ladder, p, ladder_area);
        // reported concentration of this ladder peak / its integrated area
        non_marker_ratios.push(p.concentration / area);
    }
    let ladder_mass_coefficient = mean_finite(&non_marker_ratios);

    // Ladder marker areas, ordered (lower before upper) for per-sample scaling.
    let ladder_marker_areas: Vec<f64> = ladder_markers
        .iter()
        .map(|&i| integrate_peak_area(ladder, &ladder.peaks[i], ladder_area))
        .collect();

    // --- per-sample mass coefficient and concentration ----------------------
    for (s, sample_area) in run.samples.iter_mut().zip(&areas) {
        let sample_markers = marker_indices(s, lower_conc, upper_conc);

        let mass_coefficient = if sample_markers.is_empty() {
            f64::NAN
        } else {
            let sample_marker_areas: Vec<f64> = sample_markers
                .iter()
                .map(|&i| integrate_peak_area(s, &s.peaks[i], sample_area))
                .collect();
            // Element-wise ladder/sample marker-area ratio, averaged. Both
            // marker lists are in ascending peak order (lower, then upper).
            let n = ladder_marker_areas.len().min(sample_marker_areas.len());
            let ratios: Vec<f64> = (0..n)
                .map(|k| ladder_marker_areas[k] / sample_marker_areas[k])
                .collect();
            ladder_mass_coefficient * mean_finite(&ratios)
        };

        s.concentration = sample_area.iter().map(|&a| mass_coefficient * a).collect();
    }

    Ok(())
}

/// Fill `molarity` for every sample: concentration / MW(length) · 1e6
/// (ng/µl → nmol/L or pg/µl → pmol/L). NaN wherever length/concentration is NaN.
pub fn calculate_molarity(run: &mut Electrophoresis) -> Result<()> {
    let assay_type = run.assay.assay_type.clone();
    for s in &mut run.samples {
        if s.concentration.is_empty() || s.length.is_empty() {
            s.molarity = Vec::new();
            continue;
        }
        s.molarity = s
            .concentration
            .iter()
            .zip(s.length.iter())
            .map(|(&c, &l)| {
                let mw = molecular_weight(l, &assay_type);
                if c.is_finite() && mw.is_finite() && mw > 0.0 {
                    c / mw * 1e6
                } else {
                    f64::NAN
                }
            })
            .collect();
    }
    Ok(())
}

/// Per-point trapezoidal area in aligned-time space, compensated for detector
/// dwell time (divided by aligned time). First point of the sample is NaN.
///
/// Ported from `calculate.concentration`:
///   area_i = (2·f_i − Δf_i)·Δx_i,  then /aligned_time_i   (x = aligned.time)
/// where Δf_i = f_i − f_{i-1}, Δx_i = x_i − x_{i-1}. Since (2·f_i − Δf_i) =
/// f_i + f_{i-1}, this is the trapezoid to the left of point i.
fn per_point_area(s: &Sample) -> Vec<f64> {
    let n = s.fluorescence.len();
    let mut area = vec![f64::NAN; n];
    if s.aligned_time.len() != n {
        return area; // not calibrated -> all NaN
    }
    for (i, slot) in area.iter_mut().enumerate().take(n).skip(1) {
        let f = s.fluorescence[i] as f64;
        let f_prev = s.fluorescence[i - 1] as f64;
        let dx = s.aligned_time[i] - s.aligned_time[i - 1];
        let at = s.aligned_time[i];
        if at.is_finite() && at > 0.0 && dx.is_finite() {
            *slot = (f + f_prev) * dx / at;
        }
    }
    area
}

/// Sum the per-point `area` over the points inside a peak, using the peak's
/// aligned-time boundaries (bioanalyzeR `in.peak` with x = aligned.time).
fn integrate_peak_area(s: &Sample, peak: &crate::model::Peak, area: &[f64]) -> f64 {
    let (lo, hi) = (peak.aligned_start_time, peak.aligned_end_time);
    if !lo.is_finite() || !hi.is_finite() {
        return f64::NAN;
    }
    let mut sum = 0.0;
    for (i, &at) in s.aligned_time.iter().enumerate() {
        if at >= lo && at <= hi {
            sum += area[i]; // NaN propagates, matching R sum() without na.rm
        }
    }
    sum
}

/// Indices of a sample's marker peaks: the lower marker (label + concentration
/// == `lower_conc`) and the upper marker (label + concentration == `upper_conc`)
/// if present, in ascending peak-index order.
fn marker_indices(s: &Sample, lower_conc: f64, upper_conc: f64) -> Vec<usize> {
    let mut out = Vec::new();
    for (i, p) in s.peaks.iter().enumerate() {
        let is_lower = LOWER_MARKER_NAMES.iter().any(|n| p.observations == *n)
            && approx_eq(p.concentration, lower_conc);
        let is_upper = UPPER_MARKER_NAMES.iter().any(|n| p.observations == *n)
            && approx_eq(p.concentration, upper_conc);
        if is_lower || is_upper {
            out.push(i);
        }
    }
    out
}

fn conc_in_set(c: f64, set: &[f64]) -> bool {
    set.iter().any(|&x| approx_eq(c, x))
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// Mean of the finite values (matching R's `mean(..., na.rm = TRUE)`); NaN if
/// there are none.
fn mean_finite(xs: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut count = 0usize;
    for &x in xs {
        if x.is_finite() {
            sum += x;
            count += 1;
        }
    }
    if count == 0 {
        f64::NAN
    } else {
        sum / count as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Peak;

    #[test]
    fn molecular_weight_matches_formulas() {
        assert!((molecular_weight(100.0, "DNA") - (100.0 * 607.4 + 157.9)).abs() < 1e-9);
        assert!((molecular_weight(100.0, "RNA") - (100.0 * 320.5 + 159.0)).abs() < 1e-9);
        assert!(molecular_weight(0.0, "DNA").is_nan());
        assert!(molecular_weight(-1.0, "RNA").is_nan());
        assert!(molecular_weight(100.0, "Protein").is_nan());
    }

    fn peak_between(lo: f64, hi: f64) -> Peak {
        Peak {
            observations: String::new(),
            length: f64::NAN,
            time: f64::NAN,
            aligned_time: f64::NAN,
            start_time: f64::NAN,
            end_time: f64::NAN,
            aligned_start_time: lo,
            aligned_end_time: hi,
            area: f64::NAN,
            concentration: f64::NAN,
            molarity: f64::NAN,
        }
    }

    #[test]
    fn integrate_peak_sums_area_within_bounds() {
        // Trace with aligned_time 1..=5 and per-point areas we control directly.
        let mut s = Sample {
            well_number: 1,
            name: String::new(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            fluorescence: vec![0.0; 5],
            aligned_time: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
        };
        s.peaks.push(peak_between(2.0, 4.0));
        let area = [10.0, 20.0, 30.0, 40.0, 50.0];
        // aligned_time 2,3,4 are within [2,4] -> 20+30+40 = 90
        let got = integrate_peak_area(&s, &s.peaks[0], &area);
        assert!((got - 90.0).abs() < 1e-9, "got {got}");
    }

    #[test]
    fn per_point_area_skips_non_positive_aligned_time() {
        let s = Sample {
            well_number: 1,
            name: String::new(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![0.0, 1.0, 2.0],
            fluorescence: vec![10.0, 20.0, 30.0],
            aligned_time: vec![-1.0, 0.0, 2.0],
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
        };

        let area = per_point_area(&s);

        assert!(area[0].is_nan());
        assert!(area[1].is_nan());
        assert!(area[2].is_finite());
    }
}
