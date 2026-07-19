//! Path-oriented loading API for library users.
//!
//! Format-specific modules remain available for callers that already know their
//! input type. This module provides the higher-level "read this path" behavior
//! used by applications.

use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::model::Electrophoresis;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TraceFormat {
    BioanalyzerXml,
    BioanalyzerXmlGz,
    BioanalyzerXad,
    TapeStationExport { entry: TapeStationEntry },
    FragmentAnalyzerRun { entry: FragmentAnalyzerEntry },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TapeStationEntry {
    Xml,
    Csv,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragmentAnalyzerEntry {
    Zip,
    Raw,
    Directory,
    Sidecar,
}

#[derive(Debug, Clone)]
pub struct DetectedFormat {
    pub path: PathBuf,
    pub identity: PathBuf,
    pub format: TraceFormat,
    pub capabilities: SourceCapabilities,
}

#[derive(Debug, Clone)]
pub struct LoadedRun {
    pub run: Electrophoresis,
    pub source: DetectedFormat,
    pub raw_channels: Vec<crate::xad::RawChannel>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourceCapabilities {
    pub can_rename: bool,
    pub can_save_in_place: bool,
    pub can_save_as_xml: bool,
    pub can_edit_markers: bool,
}

pub fn detect_format(path: impl AsRef<Path>) -> Result<Option<DetectedFormat>> {
    let path = path.as_ref();
    if crate::fa::is_fa_path(path) {
        let entry = if has_ext(path, "zip") {
            FragmentAnalyzerEntry::Zip
        } else if path.is_dir() {
            FragmentAnalyzerEntry::Directory
        } else if has_ext(path, "raw") {
            FragmentAnalyzerEntry::Raw
        } else {
            FragmentAnalyzerEntry::Sidecar
        };
        return Ok(Some(DetectedFormat {
            path: path.to_path_buf(),
            identity: crate::fa::run_identity(path),
            format: TraceFormat::FragmentAnalyzerRun { entry },
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: true,
                can_save_as_xml: false,
                can_edit_markers: false,
            },
        }));
    }

    if crate::tapestation::is_tapestation_path(path) {
        let entry = if file_name_lower(path).ends_with("_electropherogram.csv")
            || file_name_lower(path).ends_with("_electropherogram.csv.gz")
        {
            TapeStationEntry::Csv
        } else {
            TapeStationEntry::Xml
        };
        return Ok(Some(DetectedFormat {
            path: path.to_path_buf(),
            identity: crate::tapestation::run_identity(path).unwrap_or_else(|_| path.to_path_buf()),
            format: TraceFormat::TapeStationExport { entry },
            capabilities: SourceCapabilities {
                can_rename: false,
                can_save_in_place: false,
                can_save_as_xml: false,
                can_edit_markers: false,
            },
        }));
    }

    let lower = file_name_lower(path);
    let format = if lower.ends_with(".xad") {
        TraceFormat::BioanalyzerXad
    } else if lower.ends_with(".xml.gz") {
        TraceFormat::BioanalyzerXmlGz
    } else if lower.ends_with(".xml") {
        TraceFormat::BioanalyzerXml
    } else {
        return Ok(None);
    };
    let capabilities = match format {
        TraceFormat::BioanalyzerXml | TraceFormat::BioanalyzerXmlGz => SourceCapabilities {
            can_rename: true,
            can_save_in_place: true,
            can_save_as_xml: true,
            can_edit_markers: true,
        },
        TraceFormat::BioanalyzerXad => SourceCapabilities {
            can_rename: false,
            can_save_in_place: false,
            can_save_as_xml: true,
            can_edit_markers: false,
        },
        _ => unreachable!(),
    };
    Ok(Some(DetectedFormat {
        path: path.to_path_buf(),
        identity: path.to_path_buf(),
        format,
        capabilities,
    }))
}

pub fn read_path(path: impl AsRef<Path>) -> Result<LoadedRun> {
    let path = path.as_ref();
    let source = detect_format(path)?
        .ok_or_else(|| anyhow::anyhow!("unsupported electrophoresis path {}", path.display()))?;
    let mut warnings = Vec::new();
    let mut raw_channels = Vec::new();
    let mut run = match &source.format {
        TraceFormat::FragmentAnalyzerRun { .. } => crate::fa::read_fa_run(path)?,
        TraceFormat::TapeStationExport { .. } => crate::tapestation::read_tapestation(path)?,
        TraceFormat::BioanalyzerXad => {
            let run = crate::xad::read_xad_file(path)?;
            match crate::xad::read_xad_raw_channels(path)
                .with_context(|| format!("reading raw detector channels from {}", path.display()))
            {
                Ok(channels) => raw_channels = channels,
                Err(e) => warnings.push(format!("Raw detector channels unavailable: {e:#}")),
            }
            run
        }
        TraceFormat::BioanalyzerXml => {
            let xml = std::fs::read_to_string(path)
                .with_context(|| format!("reading {}", path.display()))?;
            crate::bioanalyzer::parse_xml(&xml)?
        }
        TraceFormat::BioanalyzerXmlGz => {
            let raw = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
            let mut d = flate2::read::GzDecoder::new(&raw[..]);
            let mut xml = String::new();
            d.read_to_string(&mut xml)
                .with_context(|| format!("gunzipping {}", path.display()))?;
            crate::bioanalyzer::parse_xml(&xml)?
        }
    };

    if matches!(
        source.format,
        TraceFormat::BioanalyzerXml | TraceFormat::BioanalyzerXmlGz | TraceFormat::BioanalyzerXad
    ) {
        if let Err(e) = run_bioanalyzer_pipeline(&mut run) {
            warnings.push(format!("calibration skipped: {e:#}"));
        }
    }

    Ok(LoadedRun {
        run,
        source,
        raw_channels,
        warnings,
    })
}

fn run_bioanalyzer_pipeline(run: &mut Electrophoresis) -> Result<()> {
    crate::calibration::calculate_length(run, crate::calibration::Method::Hyman)
        .context("sizing failed")?;
    crate::concentration::calculate_concentration(run).context("concentration failed")?;
    crate::concentration::calculate_molarity(run).context("molarity failed")?;
    Ok(())
}

fn has_ext(path: &Path, ext: &str) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case(ext))
}

fn file_name_lower(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_bioanalyzer_extensions_case_insensitively() {
        let xml = detect_format("RUN.XML").unwrap().unwrap();
        assert_eq!(xml.format, TraceFormat::BioanalyzerXml);
        assert!(xml.capabilities.can_save_as_xml);
        assert!(xml.capabilities.can_edit_markers);

        let gz = detect_format("RUN.XML.GZ").unwrap().unwrap();
        assert_eq!(gz.format, TraceFormat::BioanalyzerXmlGz);

        let xad = detect_format("RUN.XAD").unwrap().unwrap();
        assert_eq!(xad.format, TraceFormat::BioanalyzerXad);
        assert!(!xad.capabilities.can_save_in_place);
        assert!(xad.capabilities.can_save_as_xml);
    }

    #[test]
    fn detects_tapestation_csv_as_read_only_export() {
        let detected = detect_format("run_Electropherogram.csv").unwrap().unwrap();

        assert_eq!(
            detected.format,
            TraceFormat::TapeStationExport {
                entry: TapeStationEntry::Csv
            }
        );
        assert!(!detected.capabilities.can_rename);
        assert!(!detected.capabilities.can_save_in_place);
    }

    #[test]
    fn detects_fragment_analyzer_raw_identity_and_capabilities() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_io_detect_fa_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("RUN.RAW");
        std::fs::write(&raw, b"FA\0\0").unwrap();

        let detected = detect_format(&raw).unwrap().unwrap();

        assert_eq!(
            detected.format,
            TraceFormat::FragmentAnalyzerRun {
                entry: FragmentAnalyzerEntry::Raw
            }
        );
        assert_eq!(detected.identity, raw);
        assert!(detected.capabilities.can_rename);
        assert!(detected.capabilities.can_save_in_place);
        assert!(!detected.capabilities.can_edit_markers);
        std::fs::remove_dir_all(dir).unwrap();
    }
}
