//! The textual info line shown under the plot. (Plot geometry now lives in
//! [`crate::plot`], which renders via `plotters`.)

use traceio::xad::RawChannel;
use traceio::{Electrophoresis, Sample};

/// Info line for a raw detector channel (native `.xad`).
pub fn raw_info_line(ch: &RawChannel) -> String {
    let secs = ch.x_step * ch.signal.len() as f64;
    format!(
        "raw channel {}  ·  {} samples  ·  ~{:.0}s run  ·  x: time (s)  ·  processed per-well data is recomputed by 2100 Expert",
        ch.channel_id,
        ch.signal.len(),
        secs
    )
}

/// One-line summary of the selected sample: well, peaks, RIN, calibrated size
/// range, observations, and the current x-axis quantity.
pub fn info_line(run: &Electrophoresis, s: &Sample) -> String {
    let finite: Vec<f64> = s.length.iter().copied().filter(|v| v.is_finite()).collect();
    let use_len = finite.len() >= 2;

    let mut parts = vec![
        format!("Well {}", s.well_number),
        format!("{} peaks", s.peaks.len()),
    ];
    if let Some(rin) = s.rin {
        parts.push(format!("RIN {rin:.1}"));
    }
    if use_len {
        let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        parts.push(format!("{lo:.0}–{hi:.0} {}", run.assay.length_unit));
    }
    if !s.observations.is_empty() {
        parts.push(s.observations.clone());
    }
    let xlabel = if use_len {
        format!("size ({})", run.assay.length_unit)
    } else {
        "migration time (s)".to_string()
    };
    parts.push(format!("x: {xlabel}"));
    parts.join("   ·   ")
}
