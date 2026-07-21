//! GUI adapter around `crate::traceio::io::read_path`.
//!
//! The parser crate owns format detection, instrument-specific analysis, and
//! native raw-channel extraction. This module keeps the older GUI-facing
//! `Loaded` shape and retains marker-override recalibration helpers.

use std::collections::HashMap;
use std::path::Path;

use crate::traceio::calibration::MarkerOverride;
use crate::traceio::xad::RawChannel;
use crate::traceio::Electrophoresis;
use anyhow::{Context, Result};

/// A loaded run plus any raw detector channels (populated for native `.xad`).
pub struct Loaded {
    pub run: Electrophoresis,
    pub source: crate::traceio::io::DetectedFormat,
    pub raw_channels: Vec<RawChannel>,
    pub fa_metadata: Option<crate::traceio::fa::FaMetadata>,
    pub warning: Option<String>,
}

/// Load any supported electrophoresis path through the public `traceio` API.
pub fn load(path: &Path) -> Result<Loaded> {
    let loaded = crate::traceio::io::read_path_with_metadata(path)?;
    let fa_metadata = loaded.fa_metadata();
    let loaded = loaded.loaded;

    Ok(Loaded {
        run: loaded.run,
        source: loaded.source,
        raw_channels: loaded.raw_channels,
        fa_metadata,
        warning: warnings_to_gui_warning(loaded.warnings),
    })
}

/// Run the no-hardware instrument simulator over a recorded `.pck` session and
/// return the acquisition as an [`Electrophoresis`] the viewer can display, the
/// same as any loaded file. This is the GUI end of milestone M5: `--simulate`.
///
/// The trace is the provisional (uncalibrated) acquisition signal — see
/// `docs/bioanalyzer_protocol.md`. The assay name is taken from the capture's
/// header path.
pub fn simulate_pck(path: &Path) -> Result<Electrophoresis> {
    let mut transport = crate::tracehw::PckReplay::from_path(path)
        .with_context(|| format!("opening recorded session {}", path.display()))?;
    let assay = transport
        .header
        .assay_path
        .as_deref()
        .map(|p| {
            p.rsplit(['\\', '/'])
                .next()
                .unwrap_or(p)
                .trim_end_matches(".xsy")
                .to_string()
        })
        .unwrap_or_else(|| "Simulated run".to_string());

    let acq = crate::tracehw::run_to_completion(&mut transport).context("running the simulator")?;
    if acq.sample_count() == 0 {
        anyhow::bail!("no samples were produced from {}", path.display());
    }
    Ok(acq.to_electrophoresis(&assay, 1.0))
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
    use crate::traceio::{calibration, concentration};
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
