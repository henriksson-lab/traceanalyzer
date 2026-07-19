//! GUI adapter around `traceio::io::read_path`.
//!
//! The parser crate owns format detection, instrument-specific analysis, and
//! native raw-channel extraction. This module keeps the older GUI-facing
//! `Loaded` shape and retains marker-override recalibration helpers.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use traceio::calibration::MarkerOverride;
use traceio::xad::RawChannel;
use traceio::Electrophoresis;

/// A loaded run plus any raw detector channels (populated for native `.xad`).
pub struct Loaded {
    pub run: Electrophoresis,
    pub raw_channels: Vec<RawChannel>,
    pub warning: Option<String>,
}

/// Load any supported electrophoresis path through the public `traceio` API.
pub fn load(path: &Path) -> Result<Loaded> {
    let loaded = traceio::io::read_path(path)?;

    Ok(Loaded {
        run: loaded.run,
        raw_channels: loaded.raw_channels,
        warning: warnings_to_gui_warning(loaded.warnings),
    })
}

fn warnings_to_gui_warning(warnings: Vec<String>) -> Option<String> {
    if warnings.is_empty() {
        None
    } else {
        for warning in &warnings {
            eprintln!("({warning})");
        }
        Some(warnings.join("\n"))
    }
}

/// Run per-point sizing, concentration and molarity. Failures are reported but
/// non-fatal (e.g. a native `.xad` has no per-well peaks to calibrate from).
pub fn calibrate(run: &mut Electrophoresis) {
    if let Err(e) = recalibrate_with(run, &HashMap::new()) {
        eprintln!("(calibration skipped: {e:#})");
    }
}

/// Re-run the full pipeline (sizing → concentration → molarity) applying manual
/// marker overrides. Used for live recompute when the user drags a marker.
pub fn recalibrate_with(
    run: &mut Electrophoresis,
    overrides: &HashMap<usize, MarkerOverride>,
) -> Result<()> {
    use traceio::{calibration, concentration};
    calibration::calculate_length_with(run, calibration::Method::Hyman, overrides)
        .context("sizing failed")?;
    concentration::calculate_concentration(run).context("concentration failed")?;
    concentration::calculate_molarity(run).context("molarity failed")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_traceio_warnings_leave_gui_warning_empty() {
        assert!(warnings_to_gui_warning(Vec::new()).is_none());
    }

    #[test]
    fn traceio_warnings_are_joined_for_gui_warning_slot() {
        let warning = warnings_to_gui_warning(vec![
            "Raw detector channels unavailable: synthetic decode failure".into(),
            "calibration skipped: sizing failed".into(),
        ])
        .expect("traceio warnings should be surfaced through the GUI warning slot");

        assert!(warning.contains("Raw detector channels unavailable"));
        assert!(warning.contains("synthetic decode failure"));
        assert!(warning.contains("calibration skipped"));
        assert!(warning.contains('\n'));
    }
}
