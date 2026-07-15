//! Per-point ladder calibration: assign a molecular length (bp/nt) to every
//! trace point via the ladder standard curve.
//!
//! Ported from jwfoley/bioanalyzeR (`calculate.length`, MIT) and, underneath
//! it, R's own spline routines:
//!   - the FMM cubic spline from R `src/library/stats/src/splines.c`
//!     (`fmm_spline`, `spline_eval`),
//!   - the Hyman monotonicity filter and coefficient reconversion from R
//!     `src/library/stats/R/spline.R` (`hyman_filter`, `spl_coef_conv`).
//!
//! For the Bioanalyzer the mobility model is fit in *aligned-time* space, so
//! it is effectively recalibrated per sample against that sample's markers.
//! Each sample's per-point aligned time is derived from its marker peaks
//! (a linear map from raw `MigrationTime` to `AlignedMigrationTime`), matching
//! `read.bioanalyzer`.

use crate::model::{Electrophoresis, Sample};
use anyhow::{anyhow, Result};

/// Marker peak observation labels (from bioanalyzeR `electrophoresis.R`).
const LOWER_MARKER_NAMES: [&str; 2] = ["Lower Marker", "edited Lower Marker"];
const UPPER_MARKER_NAMES: [&str; 2] = ["Upper Marker", "edited Upper Marker"];

/// Mobility-model fitting method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// Hyman-filtered FMM cubic spline (bioanalyzeR default; monotone).
    Hyman,
    /// Piecewise-linear interpolation (`approxfun`); NA outside the ladder.
    Interpolation,
}

impl Default for Method {
    fn default() -> Self {
        Method::Hyman
    }
}

/// Fill `aligned_time` and `length` for every sample of a Bioanalyzer run.
///
/// Requires exactly one ladder well. Samples whose markers cannot be located
/// keep empty/NaN calibration but do not abort the run.
pub fn calculate_length(run: &mut Electrophoresis, method: Method) -> Result<()> {
    let ladder_idx = run
        .ladder_index()
        .ok_or_else(|| anyhow!("need exactly one ladder well to calibrate"))?;

    // Ladder concentrations that tag the lower/upper markers.
    let lower_conc = run.ladder_peaks.first().map(|p| p.concentration);
    let upper_conc = run.ladder_peaks.last().map(|p| p.concentration);

    // Standard curve points: ladder peaks' (aligned time -> length).
    let mut pts: Vec<(f64, f64)> = run.samples[ladder_idx]
        .peaks
        .iter()
        .filter(|p| p.aligned_time.is_finite() && p.length.is_finite())
        .map(|p| (p.aligned_time, p.length))
        .collect();
    let curve = StandardCurve::fit(&mut pts, method)?;

    let has_upper = run.assay.has_upper_marker;

    for s in &mut run.samples {
        // 1. Marker-based linear map: raw time -> aligned time, per point.
        let Some((coef, offset)) =
            alignment(s, has_upper, lower_conc, upper_conc)
        else {
            s.aligned_time = Vec::new();
            s.length = Vec::new();
            continue;
        };
        s.aligned_time = s.time.iter().map(|t| t * coef + offset).collect();

        // 2. Map aligned time -> length through the ladder curve; NA outside.
        s.length = s
            .aligned_time
            .iter()
            .map(|&at| curve.eval_in_range(at))
            .collect();
    }

    Ok(())
}

/// Compute a sample's (coefficient, offset) mapping raw time -> aligned time
/// from its lower (and, if the assay has one, upper) marker peaks.
fn alignment(
    s: &Sample,
    has_upper: bool,
    lower_conc: Option<f64>,
    upper_conc: Option<f64>,
) -> Option<(f64, f64)> {
    let lower = find_marker(s, &LOWER_MARKER_NAMES, lower_conc)?;
    if has_upper {
        let upper = find_marker(s, &UPPER_MARKER_NAMES, upper_conc)?;
        let dt = s.peaks[upper].time - s.peaks[lower].time;
        if dt == 0.0 {
            return None;
        }
        let coef = (s.peaks[upper].aligned_time - s.peaks[lower].aligned_time) / dt;
        let offset = s.peaks[lower].aligned_time - coef * s.peaks[lower].time;
        Some((coef, offset))
    } else {
        // Only a lower marker: scale through the origin.
        if s.peaks[lower].time == 0.0 {
            return None;
        }
        let coef = s.peaks[lower].aligned_time / s.peaks[lower].time;
        Some((coef, 0.0))
    }
}

/// Index of the marker peak whose label is in `names` and (when a reference
/// concentration is given) whose concentration matches it.
fn find_marker(s: &Sample, names: &[&str], conc: Option<f64>) -> Option<usize> {
    let mut found = None;
    for (i, p) in s.peaks.iter().enumerate() {
        let name_ok = names.iter().any(|n| p.observations == *n);
        let conc_ok = match conc {
            Some(c) => approx_eq(p.concentration, c),
            None => true,
        };
        if name_ok && conc_ok {
            if found.is_some() {
                return None; // conflicting markers -> refuse to guess
            }
            found = Some(i);
        }
    }
    found
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// A fitted length-vs-aligned-time standard curve.
struct StandardCurve {
    method: Method,
    x: Vec<f64>,
    y: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
}

impl StandardCurve {
    /// Fit from ladder points `(x = aligned time, y = length)`.
    fn fit(pts: &mut Vec<(f64, f64)>, method: Method) -> Result<StandardCurve> {
        // Sort by x and drop duplicate x (splines need strictly increasing x).
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        pts.dedup_by(|a, b| a.0 == b.0);
        if pts.len() < 2 {
            return Err(anyhow!("need >= 2 distinct ladder peaks to fit a curve"));
        }
        let x: Vec<f64> = pts.iter().map(|p| p.0).collect();
        let y: Vec<f64> = pts.iter().map(|p| p.1).collect();

        // splinefun(method = "hyman") requires monotone y.
        let increasing = y.windows(2).all(|w| w[1] >= w[0]);
        let decreasing = y.windows(2).all(|w| w[1] <= w[0]);
        if method == Method::Hyman && !(increasing || decreasing) {
            return Err(anyhow!(
                "ladder lengths are not monotone in migration order; \
                 cannot fit a Hyman spline"
            ));
        }

        let (b, c, d) = match method {
            Method::Hyman => {
                let (mut b, _, _) = fmm_spline(&x, &y);
                hyman_filter(&x, &y, &mut b);
                let (c, d) = spl_coef_conv(&x, &y, &b);
                (b, c, d)
            }
            // For linear interpolation the (b,c,d) polynomial is unused; eval
            // handles it directly.
            Method::Interpolation => (vec![0.0; x.len()], vec![0.0; x.len()], vec![0.0; x.len()]),
        };

        let x_min = *x.first().unwrap();
        let x_max = *x.last().unwrap();
        let y_min = y.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        Ok(StandardCurve {
            method,
            x,
            y,
            b,
            c,
            d,
            x_min,
            x_max,
            y_min,
            y_max,
        })
    }

    /// Evaluate the curve, returning NaN outside the ladder's interpolation
    /// range (matching `extrapolate = FALSE`): both the x-range and the
    /// resulting length-range must be within the ladder span.
    fn eval_in_range(&self, u: f64) -> f64 {
        if !(u >= self.x_min && u <= self.x_max) {
            return f64::NAN;
        }
        let v = self.eval(u);
        if v.is_finite() && v >= self.y_min && v <= self.y_max {
            v
        } else {
            f64::NAN
        }
    }

    fn eval(&self, u: f64) -> f64 {
        let n = self.x.len();
        // Locate segment: x[i] <= u <= x[i+1] (clamped to the ends).
        let i: usize = if u <= self.x[0] {
            0
        } else if u >= self.x[n - 1] {
            n - 1
        } else {
            let (mut lo, mut hi) = (0usize, n);
            while hi > lo + 1 {
                let mid = (lo + hi) / 2;
                if u < self.x[mid] {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            lo
        };
        match self.method {
            Method::Interpolation => {
                let j = i.min(n - 2);
                let t = (u - self.x[j]) / (self.x[j + 1] - self.x[j]);
                self.y[j] + t * (self.y[j + 1] - self.y[j])
            }
            Method::Hyman => {
                let dx = u - self.x[i];
                self.y[i] + dx * (self.b[i] + dx * (self.c[i] + dx * self.d[i]))
            }
        }
    }
}

// --- R spline internals, ported verbatim --------------------------------

/// FMM cubic spline (Forsythe, Malcolm & Moler), from R `splines.c`
/// `fmm_spline`. Returns first/second/third polynomial coefficients (b, c, d)
/// such that on `[x[i], x[i+1])`: `s(x) = y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]))`.
///
/// Ported keeping R's 1-based indexing (arrays sized `n+1`, index 0 unused).
fn fmm_spline(x: &[f64], y: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len();
    // 1-based scratch arrays.
    let mut xx = vec![0.0; n + 1];
    let mut yy = vec![0.0; n + 1];
    for i in 1..=n {
        xx[i] = x[i - 1];
        yy[i] = y[i - 1];
    }
    let mut b = vec![0.0; n + 1];
    let mut c = vec![0.0; n + 1];
    let mut d = vec![0.0; n + 1];

    if n < 3 {
        // n == 2 (n < 2 is rejected before fitting).
        let t = yy[2] - yy[1];
        b[1] = t / (xx[2] - xx[1]);
        b[2] = b[1];
        c[1] = 0.0;
        c[2] = 0.0;
        d[1] = 0.0;
        d[2] = 0.0;
        return (strip1(b), strip1(c), strip1(d));
    }

    let nm1 = n - 1;

    // Set up tridiagonal system: b = diagonal, d = offdiagonal, c = rhs.
    d[1] = xx[2] - xx[1];
    c[2] = (yy[2] - yy[1]) / d[1];
    for i in 2..n {
        d[i] = xx[i + 1] - xx[i];
        b[i] = 2.0 * (d[i - 1] + d[i]);
        c[i + 1] = (yy[i + 1] - yy[i]) / d[i];
        c[i] = c[i + 1] - c[i];
    }

    // End conditions: third derivatives from divided differences.
    b[1] = -d[1];
    b[n] = -d[nm1];
    c[1] = 0.0;
    c[n] = 0.0;
    if n > 3 {
        c[1] = c[3] / (xx[4] - xx[2]) - c[2] / (xx[3] - xx[1]);
        c[n] = c[nm1] / (xx[n] - xx[n - 2]) - c[n - 2] / (xx[nm1] - xx[n - 3]);
        c[1] = c[1] * d[1] * d[1] / (xx[4] - xx[1]);
        c[n] = -c[n] * d[nm1] * d[nm1] / (xx[n] - xx[n - 3]);
    }

    // Gaussian elimination.
    for i in 2..=n {
        let t = d[i - 1] / b[i - 1];
        b[i] -= t * d[i - 1];
        c[i] -= t * c[i - 1];
    }

    // Backward substitution.
    c[n] /= b[n];
    for i in (1..=nm1).rev() {
        c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
    }

    // Polynomial coefficients.
    b[n] = (yy[n] - yy[n - 1]) / d[nm1] + d[nm1] * (c[nm1] + 2.0 * c[n]);
    for i in 1..=nm1 {
        b[i] = (yy[i + 1] - yy[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
        d[i] = (c[i + 1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
    }
    c[n] = 3.0 * c[n];
    d[n] = d[nm1];

    (strip1(b), strip1(c), strip1(d))
}

/// Drop the unused 1-based sentinel slot, yielding a 0-based length-n vector.
fn strip1(mut v: Vec<f64>) -> Vec<f64> {
    v.remove(0);
    v
}

/// Hyman monotonicity filter on the first derivatives `b`, from R
/// `spline.R` `hyman_filter`. Operates in place on 0-based `b`.
fn hyman_filter(x: &[f64], y: &[f64], b: &mut [f64]) {
    let n = x.len();
    // ss[k] = (y[k+1]-y[k])/(x[k+1]-x[k]), length n-1.
    let ss: Vec<f64> = (0..n - 1)
        .map(|k| (y[k + 1] - y[k]) / (x[k + 1] - x[k]))
        .collect();
    // S0 = c(ss[0], ss);  S1 = c(ss, ss[n-2]).
    let s0 = |k: usize| if k == 0 { ss[0] } else { ss[k - 1] };
    let s1 = |k: usize| if k == n - 1 { ss[n - 2] } else { ss[k] };

    for k in 0..n {
        let t1 = s0(k).abs().min(s1(k).abs());
        let mut sig = b[k];
        if s0(k) * s1(k) > 0.0 {
            sig = s1(k);
        }
        if sig >= 0.0 {
            b[k] = b[k].max(0.0).min(3.0 * t1);
        } else {
            b[k] = b[k].min(0.0).max(-3.0 * t1);
        }
    }
}

/// Recompute `c` and `d` to be consistent with `y` and the (filtered) `b`,
/// from R `spline.R` `spl_coef_conv`. Returns 0-based `(c, d)` of length n.
fn spl_coef_conv(x: &[f64], y: &[f64], b: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|k| x[k + 1] - x[k]).collect();
    let yd: Vec<f64> = (0..n - 1).map(|k| -(y[k + 1] - y[k])).collect();

    let mut c = Vec::with_capacity(n);
    let mut d = Vec::with_capacity(n);
    for k in 0..n - 1 {
        let b0 = b[k];
        let b1 = b[k + 1];
        c.push(-(3.0 * yd[k] + (2.0 * b0 + b1) * h[k]) / (h[k] * h[k]));
        d.push((2.0 * yd[k] / h[k] + b0 + b1) / (h[k] * h[k]));
    }
    // Final c from the last segment; final d repeats the last interior d.
    let k = n - 2;
    let c_last = (3.0 * yd[k] + (b[k] + 2.0 * b[k + 1]) * h[k]) / (h[k] * h[k]);
    c.push(c_last);
    let d_last = d[n - 2];
    d.push(d_last);
    (c, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The spline must pass through every ladder node exactly (interpolation).
    #[test]
    fn hyman_spline_interpolates_nodes() {
        let mut pts = vec![
            (43.0, 15.0),
            (50.0, 25.0),
            (60.0, 50.0),
            (75.0, 100.0),
            (95.0, 200.0),
            (120.0, 500.0),
        ];
        let curve = StandardCurve::fit(&mut pts, Method::Hyman).unwrap();
        for &(x, y) in &pts {
            assert!((curve.eval(x) - y).abs() < 1e-6, "node ({x},{y}) not hit");
        }
        // Monotone between nodes.
        let mut prev = f64::NEG_INFINITY;
        let mut u = 43.0;
        while u <= 120.0 {
            let v = curve.eval(u);
            assert!(v >= prev - 1e-9, "not monotone at {u}");
            prev = v;
            u += 0.5;
        }
    }

    #[test]
    fn out_of_range_is_nan() {
        let mut pts = vec![(43.0, 15.0), (60.0, 50.0), (95.0, 200.0)];
        let curve = StandardCurve::fit(&mut pts, Method::Hyman).unwrap();
        assert!(curve.eval_in_range(10.0).is_nan());
        assert!(curve.eval_in_range(200.0).is_nan());
        assert!(curve.eval_in_range(60.0).is_finite());
    }
}
