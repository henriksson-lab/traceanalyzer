//! Loading electrophoresis runs from the supported file types, running the full
//! analysis pipeline (sizing → concentration → molarity), and — for native
//! `.xad` files — the raw detector channels.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use anyhow::Result;
use traceio::calibration::MarkerOverride;
use traceio::xad::RawChannel;
use traceio::Electrophoresis;

/// A loaded run plus any raw detector channels (populated for native `.xad`).
pub struct Loaded {
    pub run: Electrophoresis,
    pub raw_channels: Vec<RawChannel>,
}

/// Load a run from `.xad` / `.xml` / `.xml.gz`, calibrate it, and (for `.xad`)
/// read the raw detector channels.
pub fn load(path: &Path) -> Result<Loaded> {
    let mut run = parse(path)?;
    calibrate(&mut run);

    let raw_channels = if path.to_string_lossy().ends_with(".xad") {
        traceio::xad::read_xad_raw_channels(path).unwrap_or_default()
    } else {
        Vec::new()
    };

    Ok(Loaded { run, raw_channels })
}

/// Run per-point sizing, concentration and molarity. Failures are reported but
/// non-fatal (e.g. a native `.xad` has no per-well peaks to calibrate from).
pub fn calibrate(run: &mut Electrophoresis) {
    recalibrate_with(run, &HashMap::new());
}

/// Re-run the full pipeline (sizing → concentration → molarity) applying manual
/// marker overrides. Used for live recompute when the user drags a marker.
pub fn recalibrate_with(run: &mut Electrophoresis, overrides: &HashMap<usize, MarkerOverride>) {
    use traceio::{calibration, concentration};
    if let Err(e) = calibration::calculate_length_with(run, calibration::Method::Hyman, overrides) {
        eprintln!("(sizing skipped: {e})");
        return;
    }
    if let Err(e) = concentration::calculate_concentration(run) {
        eprintln!("(concentration skipped: {e})");
        return;
    }
    if let Err(e) = concentration::calculate_molarity(run) {
        eprintln!("(molarity skipped: {e})");
    }
}

fn parse(path: &Path) -> Result<Electrophoresis> {
    let name = path.to_string_lossy();
    if name.ends_with(".xad") {
        Ok(traceio::xad::read_xad_file(path)?)
    } else if name.ends_with(".xml.gz") {
        let raw = std::fs::read(path)?;
        let mut d = flate2::read::GzDecoder::new(&raw[..]);
        let mut s = String::new();
        d.read_to_string(&mut s)?;
        Ok(traceio::bioanalyzer::parse_xml(&s)?)
    } else {
        let s = std::fs::read_to_string(path)?;
        Ok(traceio::bioanalyzer::parse_xml(&s)?)
    }
}
