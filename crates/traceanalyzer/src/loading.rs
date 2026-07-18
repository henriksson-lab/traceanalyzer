//! Loading electrophoresis runs from the supported file types, running the full
//! analysis pipeline (sizing → concentration → molarity), and — for native
//! `.xad` files — the raw detector channels.

use std::collections::HashMap;
use std::io::Read;
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

/// Load a run from `.xad` / `.xml` / `.xml.gz`, calibrate it, and (for `.xad`)
/// read the raw detector channels.
pub fn load(path: &Path) -> Result<Loaded> {
    // Fragment Analyzer runs arrive already size-calibrated (the reader fills
    // per-point length itself), so they skip the marker-based `calibrate`.
    if traceio::fa::is_fa_path(path) {
        let run = traceio::fa::read_fa_run(path)?;
        return Ok(Loaded {
            run,
            raw_channels: Vec::new(),
            warning: None,
        });
    }

    // TapeStation exports (metadata .xml + _Electropherogram.csv) size
    // themselves against the ladder, so they also skip `calibrate`.
    if traceio::tapestation::is_tapestation_path(path) {
        let run = traceio::tapestation::read_tapestation(path)?;
        return Ok(Loaded {
            run,
            raw_channels: Vec::new(),
            warning: None,
        });
    }

    let mut run = parse(path)?;
    calibrate(&mut run);

    let (raw_channels, warning) = if is_xad_path(path) {
        read_xad_raw_channels_with_warning(path, traceio::xad::read_xad_raw_channels)
    } else {
        (Vec::new(), None)
    };

    Ok(Loaded {
        run,
        raw_channels,
        warning,
    })
}

fn is_xad_path(path: &Path) -> bool {
    path.file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .is_some_and(|n| n.ends_with(".xad"))
}

fn read_xad_raw_channels_with_warning(
    path: &Path,
    reader: impl FnOnce(&Path) -> Result<Vec<RawChannel>>,
) -> (Vec<RawChannel>, Option<String>) {
    match reader(path)
        .with_context(|| format!("reading raw detector channels from {}", path.display()))
    {
        Ok(channels) => (channels, None),
        Err(e) => {
            let msg = format!("Raw detector channels unavailable: {e:#}");
            eprintln!("({msg})");
            (Vec::new(), Some(msg))
        }
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

fn parse(path: &Path) -> Result<Electrophoresis> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
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

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::anyhow;

    #[test]
    fn xad_raw_channel_error_becomes_warning() {
        let path = Path::new("broken.xad");

        let (channels, warning) =
            read_xad_raw_channels_with_warning(path, |_| Err(anyhow!("synthetic decode failure")));

        assert!(channels.is_empty());
        let warning = warning.expect("raw-channel errors should be surfaced as warnings");
        assert!(warning.contains("Raw detector channels unavailable"));
        assert!(warning.contains("reading raw detector channels from broken.xad"));
        assert!(warning.contains("synthetic decode failure"));
    }

    #[test]
    fn non_xad_paths_do_not_request_raw_channels() {
        assert!(!is_xad_path(Path::new("run.xml")));
        assert!(!is_xad_path(Path::new("run.xml.gz")));
        assert!(is_xad_path(Path::new("run.xad")));
        assert!(is_xad_path(Path::new("RUN.XAD")));
    }
}
