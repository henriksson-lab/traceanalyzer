//! Path-oriented loading/saving API for library users.
//!
//! Format-specific modules remain available for callers that already know their
//! input type. This module provides the higher-level "read this path" and
//! "save this loaded run" behavior used by applications.

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

impl DetectedFormat {
    /// Detailed save support for this source format.
    pub fn save_capabilities(&self) -> SaveCapabilities {
        SaveCapabilities::for_source(self)
    }
}

#[derive(Debug, Clone)]
pub struct LoadedRun {
    pub run: Electrophoresis,
    pub source: DetectedFormat,
    pub raw_channels: Vec<crate::xad::RawChannel>,
    pub warnings: Vec<String>,
}

impl LoadedRun {
    pub fn new(
        run: Electrophoresis,
        source: DetectedFormat,
        raw_channels: Vec<crate::xad::RawChannel>,
        warnings: Vec<String>,
    ) -> Self {
        Self {
            run,
            source,
            raw_channels,
            warnings,
        }
    }

    /// Detailed save support for the source this run was loaded from.
    pub fn save_capabilities(&self) -> SaveCapabilities {
        self.source.save_capabilities()
    }
}

#[derive(Debug, Clone)]
pub struct LoadedRunWithMetadata {
    pub loaded: LoadedRun,
    pub fa_metadata: Option<crate::fa::FaMetadata>,
}

impl LoadedRunWithMetadata {
    pub fn new(loaded: LoadedRun, fa_metadata: Option<crate::fa::FaMetadata>) -> Self {
        Self {
            fa_metadata,
            loaded,
        }
    }

    /// Fragment Analyzer sidecar metadata, when loaded through a native FA path.
    pub fn fa_metadata(&self) -> Option<crate::fa::FaMetadata> {
        self.fa_metadata.clone()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SourceCapabilities {
    pub can_rename: bool,
    pub can_save_in_place: bool,
    pub can_save_as_xml: bool,
    pub can_edit_markers: bool,
}

/// Save support for a detected source.
///
/// `save_run` and `save_path` preserve unmodelled source data by using the
/// original source as a template. Today the supported persisted edit is the
/// per-sample name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SaveCapabilities {
    /// Sample names can be written back to a supported destination.
    pub can_save_sample_names: bool,
    /// The original source path can be used as the destination.
    pub can_save_in_place: bool,
    /// The run can be saved as plain Bioanalyzer-compatible XML.
    pub can_save_as_xml: bool,
    /// The run can be saved as gzip-compressed Bioanalyzer-compatible XML.
    pub can_save_as_xml_gz: bool,
}

impl SaveCapabilities {
    pub fn for_source(source: &DetectedFormat) -> Self {
        let mut capabilities = Self::for_format(&source.format);
        if matches!(source.format, TraceFormat::FragmentAnalyzerRun { .. })
            && !crate::fa::is_saveable_fa_destination_for_source(&source.path, &source.path)
        {
            capabilities.can_save_sample_names = false;
            capabilities.can_save_in_place = false;
        }
        capabilities
    }

    pub fn for_format(format: &TraceFormat) -> Self {
        match format {
            TraceFormat::BioanalyzerXml | TraceFormat::BioanalyzerXmlGz => Self {
                can_save_sample_names: true,
                can_save_in_place: true,
                can_save_as_xml: true,
                can_save_as_xml_gz: true,
            },
            TraceFormat::BioanalyzerXad => Self {
                can_save_sample_names: true,
                can_save_in_place: false,
                can_save_as_xml: true,
                can_save_as_xml_gz: true,
            },
            TraceFormat::FragmentAnalyzerRun { .. } => Self {
                can_save_sample_names: true,
                can_save_in_place: true,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            },
            TraceFormat::TapeStationExport { .. } => Self {
                can_save_sample_names: false,
                can_save_in_place: false,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            },
        }
    }
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
        let can_save_source_path = crate::fa::is_saveable_fa_destination_for_source(path, path);
        return Ok(Some(DetectedFormat {
            path: path.to_path_buf(),
            identity: crate::fa::run_identity(path),
            format: TraceFormat::FragmentAnalyzerRun { entry },
            capabilities: SourceCapabilities {
                can_rename: can_save_source_path,
                can_save_in_place: can_save_source_path,
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
            can_rename: true,
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
    Ok(read_path_with_metadata(path)?.loaded)
}

pub fn read_path_with_metadata(path: impl AsRef<Path>) -> Result<LoadedRunWithMetadata> {
    let path = path.as_ref();
    let source = detect_format(path)?
        .ok_or_else(|| anyhow::anyhow!("unsupported electrophoresis path {}", path.display()))?;
    let mut warnings = Vec::new();
    let mut raw_channels = Vec::new();
    let mut fa_metadata = None;
    let mut run = match &source.format {
        TraceFormat::FragmentAnalyzerRun { .. } => {
            let run = crate::fa::read_fa_run(path)?;
            match crate::fa::read_fa_metadata(path) {
                Ok(metadata) => {
                    warnings.extend(
                        metadata.warnings.iter().map(|warning| {
                            format!("Fragment Analyzer metadata warning: {warning}")
                        }),
                    );
                    fa_metadata = Some(metadata);
                }
                Err(e) => {
                    warnings.push(format!("Fragment Analyzer metadata unavailable: {e:#}"));
                }
            }
            run
        }
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

    Ok(LoadedRunWithMetadata::new(
        LoadedRun::new(run, source, raw_channels, warnings),
        fa_metadata,
    ))
}

/// Persist edits from a loaded run to `dst`, using the run's source path as the
/// template.
///
/// This is the path-oriented counterpart to [`read_path`]. It delegates to
/// [`crate::save::save_run`], so existing callers that already track source and
/// destination paths can keep using that lower-level API.
pub fn save_path(loaded: &LoadedRun, dst: impl AsRef<Path>) -> Result<()> {
    crate::save::save_path(loaded, dst)
}

/// Whether [`save_path`] supports writing this loaded run to `dst`.
///
/// This is a coarse preflight helper for UI/library gates; [`save_path`] remains
/// the authoritative operation and returns the detailed error if saving fails.
pub fn supports_save_path(loaded: &LoadedRun, dst: impl AsRef<Path>) -> bool {
    crate::save::supports_save_path(loaded, dst)
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
        assert_eq!(
            xml.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: true,
                can_save_in_place: true,
                can_save_as_xml: true,
                can_save_as_xml_gz: true,
            }
        );

        let gz = detect_format("RUN.XML.GZ").unwrap().unwrap();
        assert_eq!(gz.format, TraceFormat::BioanalyzerXmlGz);

        let xad = detect_format("RUN.XAD").unwrap().unwrap();
        assert_eq!(xad.format, TraceFormat::BioanalyzerXad);
        assert!(xad.capabilities.can_rename);
        assert!(!xad.capabilities.can_save_in_place);
        assert!(xad.capabilities.can_save_as_xml);
        assert_eq!(
            xad.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: true,
                can_save_in_place: false,
                can_save_as_xml: true,
                can_save_as_xml_gz: true,
            }
        );
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
        assert_eq!(
            detected.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: false,
                can_save_in_place: false,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            }
        );
    }

    #[test]
    fn tapestation_csv_in_fa_directory_is_not_detected_as_fa_sidecar() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_io_detect_ts_csv_in_fa_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("run.raw"), b"FA\0\0").unwrap();
        let csv = dir.join("run_Electropherogram.csv");
        std::fs::write(&csv, b"Sample,Time,Value\n").unwrap();

        let detected = detect_format(&csv).unwrap().unwrap();

        assert_eq!(
            detected.format,
            TraceFormat::TapeStationExport {
                entry: TapeStationEntry::Csv
            }
        );
        std::fs::remove_dir_all(dir).unwrap();
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
        std::fs::write(
            dir.join("RUN.txt"),
            "Capillary #: 1\nWell: D1\nSample ID: one\n",
        )
        .unwrap();

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
        assert_eq!(
            detected.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: true,
                can_save_in_place: true,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            }
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn arbitrary_fa_sibling_load_source_does_not_advertise_save_in_place() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_io_detect_fa_arbitrary_sibling_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        let notes = dir.join("notes.txt");
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(
            dir.join("run.txt"),
            "Capillary #: 1\nWell: D1\nSample ID: old\n",
        )
        .unwrap();
        std::fs::write(&notes, "not an FA sidecar").unwrap();

        let detected = detect_format(&notes).unwrap().unwrap();

        assert_eq!(
            detected.format,
            TraceFormat::FragmentAnalyzerRun {
                entry: FragmentAnalyzerEntry::Sidecar
            }
        );
        assert_eq!(detected.identity, raw);
        assert!(!detected.capabilities.can_rename);
        assert!(!detected.capabilities.can_save_in_place);
        assert_eq!(
            detected.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: false,
                can_save_in_place: false,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            }
        );

        let loaded = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            detected,
            vec![],
            vec![],
        );
        assert!(!supports_save_path(&loaded, &notes));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_zip_with_unpreserved_entry_metadata_does_not_advertise_save() {
        use std::io::Write;

        let zip_path = std::env::temp_dir().join(format!(
            "traceio_io_fa_zip_unsaveable_metadata_{}_{}.zip",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let out = std::fs::File::create(&zip_path).unwrap();
            let mut zip = zip::ZipWriter::new(out);
            let opts = zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated);
            zip.start_file("run.raw", opts).unwrap();
            zip.write_all(b"FA\0\0").unwrap();
            let mut txt_opts: zip::write::FullFileOptions<'_> = zip::write::FileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated);
            txt_opts
                .add_extra_data(0xcafe, vec![1, 2, 3].into_boxed_slice(), false)
                .unwrap();
            zip.start_file("run.txt", txt_opts).unwrap();
            zip.write_all(b"Capillary #: 1\nWell: D1\nSample ID: old\n")
                .unwrap();
            zip.finish().unwrap();
        }

        let detected = detect_format(&zip_path).unwrap().unwrap();

        assert!(!detected.capabilities.can_rename);
        assert!(!detected.capabilities.can_save_in_place);
        assert_eq!(
            detected.save_capabilities(),
            SaveCapabilities {
                can_save_sample_names: false,
                can_save_in_place: false,
                can_save_as_xml: false,
                can_save_as_xml_gz: false,
            }
        );

        let loaded = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            detected,
            vec![],
            vec![],
        );
        assert!(!supports_save_path(&loaded, &zip_path));
        std::fs::remove_file(zip_path).unwrap();
    }

    #[test]
    fn read_path_with_metadata_warns_for_bad_fa_current_without_dropping_metadata() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_io_fa_bad_current_warning_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        std::fs::write(&raw, synthetic_fa_raw(20, 3, &[5, 14])).unwrap();
        std::fs::write(
            dir.join("run.txt"),
            "Raw file: C:\\AATI\\Data\\run.raw\nOperator: admin\n\
             Capillary #: 1\nWell: D1\nSample ID: alpha\n\
             Capillary #: 2\nWell: D2\nSample ID: beta\n",
        )
        .unwrap();
        std::fs::write(dir.join("method.mthd"), "[Separation]\nKV=7.00\n").unwrap();
        std::fs::write(
            dir.join("run.current"),
            "Current(uA)\tVoltage(kV)\tPressure(PSI)\nnot-a-current-row\n",
        )
        .unwrap();

        let loaded = read_path_with_metadata(&raw).unwrap();
        let metadata = loaded.fa_metadata().unwrap();

        assert_eq!(metadata.run_header["Operator"], "admin");
        assert_eq!(metadata.method["Separation"]["KV"], "7.00");
        assert!(metadata.current.is_empty());
        assert!(
            loaded
                .loaded
                .warnings
                .iter()
                .any(|warning| warning.contains(".current ignored because it is malformed")),
            "got {:?}",
            loaded.loaded.warnings
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn loaded_run_constructor_has_no_fa_metadata() {
        let loaded = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            DetectedFormat {
                path: PathBuf::from("run.xml"),
                identity: PathBuf::from("run.xml"),
                format: TraceFormat::BioanalyzerXml,
                capabilities: SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: true,
                    can_edit_markers: true,
                },
            },
            vec![],
            vec![],
        );

        assert!(loaded.warnings.is_empty());
    }

    #[test]
    fn supports_save_path_checks_destination_shape() {
        let loaded = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            DetectedFormat {
                path: PathBuf::from("run.xml"),
                identity: PathBuf::from("run.xml"),
                format: TraceFormat::BioanalyzerXml,
                capabilities: SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: true,
                    can_edit_markers: true,
                },
            },
            vec![],
            vec![],
        );

        assert!(supports_save_path(&loaded, "out.xml"));
        assert!(supports_save_path(&loaded, "out.xml.gz"));
        assert!(!supports_save_path(&loaded, "out.xad"));
    }

    #[test]
    fn supports_save_path_allows_xad_sources_to_xml_but_not_xad_in_place() {
        let loaded = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            DetectedFormat {
                path: PathBuf::from("run.xad"),
                identity: PathBuf::from("run.xad"),
                format: TraceFormat::BioanalyzerXad,
                capabilities: SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: false,
                    can_save_as_xml: true,
                    can_edit_markers: false,
                },
            },
            vec![],
            vec![],
        );

        assert!(supports_save_path(&loaded, "out.xml"));
        assert!(supports_save_path(&loaded, "out.xml.gz"));
        assert!(!supports_save_path(&loaded, "run.xad"));
    }

    #[test]
    fn fa_metadata_is_owned_by_loaded_run_not_identity() {
        let source = DetectedFormat {
            path: PathBuf::from("run.raw"),
            identity: PathBuf::from("run.raw"),
            format: TraceFormat::FragmentAnalyzerRun {
                entry: FragmentAnalyzerEntry::Raw,
            },
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: true,
                can_save_as_xml: false,
                can_edit_markers: false,
            },
        };
        let with_metadata = LoadedRunWithMetadata::new(
            LoadedRun::new(
                Electrophoresis {
                    assay: Default::default(),
                    ladder_peaks: vec![],
                    regions: vec![],
                    samples: vec![],
                },
                source.clone(),
                vec![],
                vec![],
            ),
            Some(crate::fa::FaMetadata {
                run_header: [("Operator".to_string(), "first".to_string())]
                    .into_iter()
                    .collect(),
                ..Default::default()
            }),
        );
        assert_eq!(
            with_metadata
                .fa_metadata()
                .as_ref()
                .and_then(|metadata| metadata.run_header.get("Operator"))
                .map(String::as_str),
            Some("first")
        );

        let without_metadata = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            source,
            vec![],
            vec!["Fragment Analyzer metadata unavailable: test failure".to_string()],
        );

        assert_eq!(
            without_metadata.warnings,
            vec!["Fragment Analyzer metadata unavailable: test failure".to_string()]
        );
    }

    #[test]
    fn preexisting_loaded_run_does_not_gain_later_identity_metadata() {
        let source = DetectedFormat {
            path: PathBuf::from("run.raw"),
            identity: PathBuf::from("run.raw"),
            format: TraceFormat::FragmentAnalyzerRun {
                entry: FragmentAnalyzerEntry::Raw,
            },
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: true,
                can_save_as_xml: false,
                can_edit_markers: false,
            },
        };
        let old_manual = LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            source.clone(),
            vec![],
            vec![],
        );

        let later_loaded = LoadedRunWithMetadata::new(
            LoadedRun::new(
                Electrophoresis {
                    assay: Default::default(),
                    ladder_peaks: vec![],
                    regions: vec![],
                    samples: vec![],
                },
                source,
                vec![],
                vec![],
            ),
            Some(crate::fa::FaMetadata {
                run_header: [("Operator".to_string(), "later".to_string())]
                    .into_iter()
                    .collect(),
                ..Default::default()
            }),
        );

        assert!(old_manual.warnings.is_empty());
        assert_eq!(
            later_loaded
                .fa_metadata()
                .as_ref()
                .and_then(|metadata| metadata.run_header.get("Operator"))
                .map(String::as_str),
            Some("later")
        );
    }

    #[test]
    fn loaded_run_keeps_original_public_field_pattern() {
        let loaded = LoadedRun {
            run: Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![],
            },
            source: DetectedFormat {
                path: PathBuf::from("run.xml"),
                identity: PathBuf::from("run.xml"),
                format: TraceFormat::BioanalyzerXml,
                capabilities: SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: true,
                    can_edit_markers: true,
                },
            },
            raw_channels: vec![],
            warnings: vec![],
        };
        let LoadedRun {
            run,
            source,
            raw_channels,
            warnings,
        } = loaded;
        assert!(run.samples.is_empty());
        assert_eq!(source.format, TraceFormat::BioanalyzerXml);
        assert!(raw_channels.is_empty());
        assert!(warnings.is_empty());
    }

    fn put_u16(buf: &mut [u8], off: usize, value: u16) {
        buf[off..off + 2].copy_from_slice(&value.to_be_bytes());
    }

    fn synthetic_fa_raw(width: usize, scans: usize, columns: &[u16]) -> Vec<u8> {
        const DATA_START: usize = 0x7d0;
        let mut raw = vec![0u8; DATA_START + scans * width * 2];
        raw[..4].copy_from_slice(b"FA\0\0");
        put_u16(&mut raw, 0xff, width as u16);

        let table = 0x40;
        put_u16(&mut raw, table, 0);
        for (i, &col) in columns.iter().enumerate() {
            put_u16(&mut raw, table + 2 + i * 2, col);
        }
        put_u16(&mut raw, table + 2 + columns.len() * 2, 0);

        let mut off = DATA_START;
        for scan in 0..scans {
            for col in 0..width {
                let value = (scan * 100 + col) as u16;
                raw[off..off + 2].copy_from_slice(&value.to_be_bytes());
                off += 2;
            }
        }
        raw
    }
}
