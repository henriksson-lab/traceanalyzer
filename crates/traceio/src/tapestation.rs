//! Reader for Agilent **TapeStation** exported data (metadata XML + trace CSV).
//!
//! Ported from jwfoley/bioanalyzeR `R/tapestation.R` (MIT). The native
//! TapeStation Analysis project file (`.D1000`, `.HSD1000`, `.RNA`, …) is a
//! password-encrypted ZIP with no public key, so — like bioanalyzeR — this
//! reader consumes the **export**: *File → Export Data* in TapeStation Analysis
//! Software (v4.1+) writes two paired files per run:
//!
//! * `<name>.xml` — all metadata (samples, peaks, regions, units) but no trace,
//! * `<name>_Electropherogram.csv` — raw fluorescence only (one column per lane,
//!   one row per distance reading), Latin-1, with a header row.
//!
//! Either file opens the run; the reader derives the other by the
//! `_Electropherogram.csv` naming convention. Both may be gzip-compressed.
//! Sizing follows bioanalyzeR: each sample's lower/upper marker define a
//! marker-relative distance axis, and the ladder sample's peaks fit a monotone
//! spline from relative distance to length (bp/nt).

use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use roxmltree::{Document, Node};

use crate::calibration::StandardCurve;
use crate::model::{AssayInfo, Electrophoresis, Peak, Region, Sample};

/// Marker peak observation labels (shared with the Bioanalyzer model).
const LOWER_MARKER_NAMES: [&str; 2] = ["Lower Marker", "edited Lower Marker"];
const UPPER_MARKER_NAMES: [&str; 2] = ["Upper Marker", "edited Upper Marker"];

/// True if `path` is a TapeStation export entry point: an `_Electropherogram.csv`
/// (optionally gzipped), or an `.xml`/`.xml.gz` whose content looks like a
/// TapeStation metadata export (rather than a Bioanalyzer `<Chipset>` document).
pub fn is_tapestation_path(path: &Path) -> bool {
    let lower = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    if lower.ends_with("_electropherogram.csv") || lower.ends_with("_electropherogram.csv.gz") {
        return true;
    }
    if lower.ends_with(".xml") || lower.ends_with(".xml.gz") {
        if let Ok(text) = read_text(path) {
            return looks_like_tapestation_xml(&text);
        }
    }
    false
}

/// Cheap content check: TapeStation metadata has `<FileInformation>`/`<Samples>`
/// and lacks the Bioanalyzer `<Chipset>` root.
pub fn looks_like_tapestation_xml(xml: &str) -> bool {
    xml.contains("<FileInformation") && xml.contains("<Samples") && !xml.contains("<Chipset")
}

/// Read a TapeStation run from either exported file (the `.xml` or the
/// `_Electropherogram.csv`); the sibling is located by naming convention.
pub fn read_tapestation(path: &Path) -> Result<Electrophoresis> {
    let (xml_path, csv_path) = resolve_pair(path)?;
    let xml = read_text(&xml_path).with_context(|| format!("reading {}", xml_path.display()))?;
    let csv = csv_path
        .as_ref()
        .map(|p| read_text(p).with_context(|| format!("reading {}", p.display())))
        .transpose()?;

    let mut run = parse_xml(&xml)?;
    if let Some(csv) = csv.as_deref() {
        attach_traces(&mut run, csv);
    }
    calibrate(&mut run);
    Ok(run)
}

/// Canonical identity for a TapeStation export pair. Both the metadata XML and
/// `_Electropherogram.csv` entry point identify the metadata XML.
pub fn run_identity(path: &Path) -> Result<PathBuf> {
    resolve_pair(path).map(|(xml, _)| xml)
}

/// Given either export file, return `(xml_path, csv_path)`. The CSV is optional
/// (metadata still loads without a trace).
fn resolve_pair(path: &Path) -> Result<(PathBuf, Option<PathBuf>)> {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let lower = name.to_ascii_lowercase();
    if lower.ends_with("_electropherogram.csv") || lower.ends_with("_electropherogram.csv.gz") {
        // CSV given: derive the .xml stem by stripping the suffix.
        let gz = lower.ends_with(".gz");
        let base = &name[..name.len()
            - if gz {
                "_Electropherogram.csv.gz".len()
            } else {
                "_Electropherogram.csv".len()
            }];
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let xml = first_existing(dir, base, &[".xml", ".xml.gz"])?
            .ok_or_else(|| anyhow!("no metadata .xml next to {}", path.display()))?;
        return Ok((xml, Some(path.to_path_buf())));
    }
    // XML given: derive the electropherogram CSV.
    let gz = lower.ends_with(".xml.gz");
    if !gz && !lower.ends_with(".xml") {
        bail!(
            "unsupported TapeStation export extension for {}",
            path.display()
        );
    }
    let base = &name[..name.len() - if gz { ".xml.gz".len() } else { ".xml".len() }];
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    let csv = first_existing(
        dir,
        base,
        &["_Electropherogram.csv", "_Electropherogram.csv.gz"],
    )?;
    Ok((path.to_path_buf(), csv))
}

fn first_existing(dir: &Path, base: &str, suffixes: &[&str]) -> Result<Option<PathBuf>> {
    for suffix in suffixes {
        let exact = dir.join(format!("{base}{suffix}"));
        if exact.exists() {
            return Ok(Some(exact));
        }
    }

    let wanted = suffixes
        .iter()
        .map(|s| format!("{base}{s}").to_ascii_lowercase())
        .collect::<Vec<_>>();
    let mut matches = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let path = entry?.path();
        let Some(name) = path
            .file_name()
            .map(|n| n.to_string_lossy().to_ascii_lowercase())
        else {
            continue;
        };
        if wanted.iter().any(|w| w == &name) {
            matches.push(path);
        }
    }

    match matches.len() {
        0 => Ok(None),
        1 => Ok(matches.pop()),
        _ => Err(anyhow!(
            "ambiguous TapeStation sibling files for {base}: {}",
            matches
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

/// Read a possibly gzipped file as text. XML is UTF-8; the CSV is Latin-1 — both
/// decode losslessly here (Latin-1 maps each byte to U+00XX).
fn read_text(path: &Path) -> Result<String> {
    let raw = std::fs::read(path)?;
    let bytes = if path.to_string_lossy().to_ascii_lowercase().ends_with(".gz") {
        let mut d = flate2::read::GzDecoder::new(&raw[..]);
        let mut v = Vec::new();
        d.read_to_end(&mut v)?;
        v
    } else {
        raw
    };
    // Try UTF-8; fall back to Latin-1 (byte→codepoint) so CSV/XML both work.
    Ok(String::from_utf8(bytes.clone())
        .unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect()))
}

// --- XML metadata -----------------------------------------------------------

/// Parse the metadata XML into a run with samples/peaks/regions but no trace
/// (attached later from the CSV) and no per-point sizing (done in `calibrate`).
pub fn parse_xml(xml: &str) -> Result<Electrophoresis> {
    let doc = Document::parse(xml).context("parsing TapeStation XML")?;
    let root = doc.root_element();

    let file_info = child(root, "FileInformation");
    let assay_el = child(root, "Assay");
    let units = assay_el.and_then(|a| child(a, "Units"));

    let assay_name = file_info
        .and_then(|f| child_text(f, "Assay"))
        .unwrap_or_default();
    let length_unit = units
        .and_then(|u| child_text(u, "MolecularWeightUnit"))
        .unwrap_or_else(|| "bp".to_string());
    let concentration_unit = units
        .and_then(|u| child_text(u, "ConcentrationUnit"))
        .unwrap_or_default();
    // Molarity unit is in the XML (`MolarityUnit`); fall back to the bioanalyzeR
    // mapping from the concentration unit if it is absent.
    let molarity_unit = units
        .and_then(|u| child_text(u, "MolarityUnit"))
        .or_else(|| match concentration_unit.as_str() {
            "ng/µl" | "ng/ul" => Some("nM".to_string()),
            "pg/µl" | "pg/ul" => Some("pM".to_string()),
            _ => None,
        });
    let assay_type = if assay_name.to_ascii_uppercase().contains("RNA") {
        "RNA"
    } else {
        "DNA"
    };

    let samples: Vec<Sample> = child(root, "Samples")
        .map(|s| {
            s.children()
                .filter(Node::is_element)
                .map(parse_sample)
                .collect()
        })
        .unwrap_or_default();

    let has_upper_marker = samples
        .iter()
        .flat_map(|s| &s.peaks)
        .any(|p| UPPER_MARKER_NAMES.contains(&p.observations.as_str()));

    let assay = AssayInfo {
        file_name: file_info
            .and_then(|f| child_text(f, "FileName"))
            .unwrap_or_default(),
        creation_date: file_info
            .and_then(|f| child_text(f, "RunEndDate"))
            .unwrap_or_default(),
        assay_name,
        assay_type: assay_type.to_string(),
        length_unit,
        concentration_unit,
        molarity_unit,
        has_upper_marker,
    };

    Ok(Electrophoresis {
        assay,
        ladder_peaks: Vec::new(),
        regions: Vec::new(),
        samples,
    })
}

fn parse_sample(node: Node) -> Sample {
    let name = child_text(node, "Comment").unwrap_or_default();
    let observations = child_text(node, "Observations").unwrap_or_default();
    let well_label = child_text(node, "WellNumber").unwrap_or_default();
    // Integrity number: RINe (RNA) or DIN (gDNA); first finite one wins.
    let rine = child(node, "RNA")
        .map(|r| child_num(r, "RINe"))
        .unwrap_or(f64::NAN);
    let rin = [rine, child_num(node, "DIN")]
        .into_iter()
        .find(|v| v.is_finite());

    let peaks = child(node, "Peaks")
        .map(|p| {
            p.children()
                .filter(Node::is_element)
                .map(parse_peak)
                .collect()
        })
        .unwrap_or_default();
    // Per-sample smear-analysis regions.
    let regions: Vec<Region> = child(node, "Regions")
        .map(|r| {
            r.children()
                .filter(Node::is_element)
                .map(parse_region)
                .collect()
        })
        .unwrap_or_default();
    // ScreenTapeID groups which ladder calibrates which samples; stashed in
    // `category` (otherwise unused) so `calibrate` can size per tape.
    let screentape_id = child_text(node, "ScreenTapeID").unwrap_or_default();

    Sample {
        well_number: well_to_number(&well_label),
        name,
        category: screentape_id,
        is_ladder: is_ladder_observation(&observations),
        comment: String::new(),
        observations,
        rin,
        time: Vec::new(),
        fluorescence: Vec::new(),
        aligned_time: Vec::new(),
        length: Vec::new(),
        concentration: Vec::new(),
        molarity: Vec::new(),
        peaks,
        regions,
    }
}

fn parse_peak(node: Node) -> Peak {
    // RunDistance/FromPercent/ToPercent are stored as percent×100 (→ 0..1).
    let distance = child_num(node, "RunDistance") / 100.0;
    let lower = child_num(node, "FromPercent") / 100.0;
    let upper = child_num(node, "ToPercent") / 100.0;
    Peak {
        observations: child_text(node, "Observations").unwrap_or_default(),
        length: child_num(node, "Size"),
        time: distance,
        aligned_time: f64::NAN,
        start_time: lower,
        end_time: upper,
        aligned_start_time: f64::NAN,
        aligned_end_time: f64::NAN,
        area: child_num(node, "Area"),
        concentration: child_num(node, "CalibratedQuantity"),
        molarity: child_num(node, "Molarity"),
    }
}

fn parse_region(node: Node) -> Region {
    Region {
        lower_length: child_num(node, "From"),
        upper_length: child_num(node, "To"),
    }
}

/// Ladder detection (bioanalyzeR): `Observations` contains "Ladder" but not the
/// "run as sample" / "sizing changed" qualifiers that would make it a sample.
fn is_ladder_observation(obs: &str) -> bool {
    obs.contains("Ladder")
        && !obs.contains("Ladder run as sample")
        && !obs.contains("Ladder sizing changed")
}

// --- CSV trace --------------------------------------------------------------

/// Attach the electropherogram CSV columns to the samples (matched by order).
/// The distance axis is the reversed, normalized row index (bioanalyzeR):
/// `distance = rev(1..N) / N`, stored in each sample's `time` field.
fn attach_traces(run: &mut Electrophoresis, csv: &str) {
    let mut lines = csv.lines().filter(|l| !l.trim().is_empty());
    let Some(_header) = lines.next() else { return };
    let rows: Vec<Vec<f32>> = lines
        .map(|l| {
            l.split(',')
                .map(|c| c.trim().parse::<f32>().unwrap_or(f32::NAN))
                .collect()
        })
        .collect();
    let n = rows.len();
    if n == 0 {
        return;
    }
    let ncol = rows.iter().map(Vec::len).max().unwrap_or(0);
    let distance: Vec<f64> = (0..n).map(|i| (n - i) as f64 / n as f64).collect();

    for (col, sample) in run.samples.iter_mut().enumerate() {
        if col >= ncol {
            break;
        }
        sample.time = distance.clone();
        sample.fluorescence = rows
            .iter()
            .map(|r| r.get(col).copied().unwrap_or(f32::NAN))
            .collect();
    }
}

// --- sizing (distance → length) --------------------------------------------

/// Fill each sample's per-point `length` by mapping its marker-relative distance
/// through its tape's ladder spline. A run may hold several ScreenTapes, each
/// with its own ladder (`Sample::category` = `ScreenTapeID`); a sample is sized
/// by the ladder that shares its tape, or by the sole ladder when there is only
/// one. No ladder / no markers leaves `length` empty (trace on distance axis).
fn calibrate(run: &mut Electrophoresis) {
    use std::collections::HashMap;
    let has_upper = run.assay.has_upper_marker;

    // One relative-distance→size curve per tape that has a usable ladder.
    let mut curves: HashMap<String, StandardCurve> = HashMap::new();
    for s in &run.samples {
        if !s.is_ladder || curves.contains_key(&s.category) {
            continue;
        }
        if let Some(curve) = ladder_curve(s, has_upper) {
            curves.insert(s.category.clone(), curve);
        }
    }
    if curves.is_empty() {
        return;
    }
    // Fallback when a sample's tape has no ladder: the sole curve, if unambiguous.
    let single = (curves.len() == 1).then(|| curves.values().next().unwrap());

    for s in &mut run.samples {
        let Some(curve) = curves.get(&s.category).or(single) else {
            continue;
        };
        let Some((lo, up)) = marker_distances(s, has_upper) else {
            continue;
        };
        s.length = s
            .time
            .iter()
            .map(|&d| curve.eval_in_range(relative_distance(d, lo, up)))
            .collect();
    }
}

/// Fit a ladder sample's relative-distance→size spline.
fn ladder_curve(ladder: &Sample, has_upper: bool) -> Option<StandardCurve> {
    let (lo, up) = marker_distances(ladder, has_upper)?;
    let pts: Vec<(f64, f64)> = ladder
        .peaks
        .iter()
        .filter(|p| p.length.is_finite() && p.time.is_finite())
        .map(|p| (relative_distance(p.time, lo, up), p.length))
        .filter(|(rd, _)| rd.is_finite())
        .collect();
    StandardCurve::fit_hyman(&pts).ok()
}

/// A sample's `(lower_marker_distance, upper_marker_distance)`. When the assay
/// has no upper marker, the upper reference is 0 (bioanalyzeR convention).
fn marker_distances(s: &Sample, has_upper: bool) -> Option<(f64, f64)> {
    let lower = s
        .peaks
        .iter()
        .find(|p| LOWER_MARKER_NAMES.contains(&p.observations.as_str()))
        .map(|p| p.time)?;
    let upper = if has_upper {
        s.peaks
            .iter()
            .find(|p| UPPER_MARKER_NAMES.contains(&p.observations.as_str()))
            .map(|p| p.time)?
    } else {
        0.0
    };
    if (lower - upper).abs() < f64::EPSILON {
        return None;
    }
    Some((lower, upper))
}

fn relative_distance(d: f64, lower: f64, upper: f64) -> f64 {
    (d - upper) / (lower - upper)
}

// --- small helpers ----------------------------------------------------------

/// Map a well label like `A1`/`H12` to a 1-based number (row-major, 12/row).
fn well_to_number(well: &str) -> i32 {
    let well = well.trim();
    let mut chars = well.chars();
    let Some(row) = chars.next().map(|c| c.to_ascii_uppercase()) else {
        return 0;
    };
    if !row.is_ascii_alphabetic() {
        return well.parse().unwrap_or(0);
    }
    let col: i32 = chars.as_str().trim().parse().unwrap_or(0);
    (row as i32 - 'A' as i32) * 12 + col
}

fn child<'a, 'i>(node: Node<'a, 'i>, tag: &str) -> Option<Node<'a, 'i>> {
    node.children()
        .find(|c| c.is_element() && c.has_tag_name(tag))
}

fn child_text(node: Node, tag: &str) -> Option<String> {
    let c = child(node, tag)?;
    let text: String = c.children().filter_map(|t| t.text()).collect();
    let text = text.trim().to_string();
    (!text.is_empty()).then_some(text)
}

/// Parse a child element's text as `f64`, or NaN when absent/blank/non-numeric.
fn child_num(node: Node, tag: &str) -> f64 {
    child_text(node, tag)
        .and_then(|t| t.replace(',', ".").parse().ok())
        .unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    // A minimal but schema-faithful export: a ladder + one sample, sharing
    // markers at distance 0.845 (25 bp) and 0.15 (1500 bp), with a ladder band
    // at distance 0.50 (500 bp) so calibration is checkable.
    fn xml() -> String {
        let peak = |size: &str, dist: &str, obs: &str, area: &str, conc: &str, mol: &str| {
            format!(
                "<Peak><Size>{size}</Size><RunDistance>{dist}</RunDistance>\
                 <FromPercent>10</FromPercent><ToPercent>20</ToPercent>\
                 <Observations>{obs}</Observations><Area>{area}</Area>\
                 <CalibratedQuantity>{conc}</CalibratedQuantity><Molarity>{mol}</Molarity></Peak>"
            )
        };
        let ladder_peaks = format!(
            "{}{}{}",
            peak("25", "84.5", "Lower Marker", "1.0", "7.5", "462"),
            peak("500", "50.0", "", "0.3", "2.0", "6.4"),
            peak("1500", "15.0", "Upper Marker", "1.0", "6.5", "6.7"),
        );
        let sample_peaks = format!(
            "{}{}{}",
            peak("25", "84.5", "Lower Marker", "1.0", "8.0", "500"),
            peak("300", "60.0", "", "0.9", "5.0", "26"),
            peak("1500", "15.0", "Upper Marker", "1.0", "6.5", "6.7"),
        );
        format!(
            "<TapeStation>\
             <FileInformation><FileName>demo.D1000</FileName><RunEndDate>2020-01-01</RunEndDate>\
             <Assay>D1000</Assay></FileInformation>\
             <Assay><Units><MolecularWeightUnit>bp</MolecularWeightUnit>\
             <ConcentrationUnit>ng/µl</ConcentrationUnit><MolarityUnit>nmol/l</MolarityUnit></Units></Assay>\
             <Samples>\
             <Sample><WellNumber>A1</WellNumber><Comment>Ladder</Comment><Observations>Ladder</Observations>\
             <ScreenTapeID>TAPE-1</ScreenTapeID>\
             <RNA/><Peaks>{ladder_peaks}</Peaks></Sample>\
             <Sample><WellNumber>B1</WellNumber><Comment>Sample X</Comment><Observations></Observations>\
             <ScreenTapeID>TAPE-1</ScreenTapeID>\
             <RNA><RINe>8.4</RINe></RNA><Peaks>{sample_peaks}</Peaks>\
             <Regions><Region><From>200</From><To>400</To><AverageSize>300</AverageSize>\
             <Concentration>5.0</Concentration><Molarity>26</Molarity><PercentOfTotal>90</PercentOfTotal>\
             </Region></Regions></Sample>\
             </Samples></TapeStation>"
        )
    }

    // 10 rows → distances 1.0, 0.9, …, 0.1; row 5 has distance 0.50.
    fn csv() -> String {
        let mut s = String::from("A1: Ladder,B1: Sample X\n");
        for i in 0..10 {
            s.push_str(&format!("{},{}\n", i as f32, (i as f32) * 2.0));
        }
        s
    }

    #[test]
    fn parses_samples_peaks_and_units() {
        let mut run = parse_xml(&xml()).unwrap();
        attach_traces(&mut run, &csv());

        assert_eq!(run.assay.length_unit, "bp");
        assert_eq!(run.assay.concentration_unit, "ng/µl");
        assert_eq!(run.assay.molarity_unit.as_deref(), Some("nmol/l"));
        assert!(run.assay.has_upper_marker);
        assert_eq!(run.samples.len(), 2);

        let ladder = &run.samples[0];
        assert!(ladder.is_ladder);
        assert_eq!(ladder.name, "Ladder");
        assert_eq!(ladder.well_number, 1); // A1
        assert_eq!(ladder.peaks.len(), 3);
        assert_eq!(ladder.peaks[0].observations, "Lower Marker");
        assert_eq!(ladder.peaks[2].observations, "Upper Marker");
        assert!((ladder.peaks[1].length - 500.0).abs() < 1e-9);
        assert!((ladder.peaks[1].area - 0.3).abs() < 1e-9);
        assert!((ladder.peaks[1].concentration - 2.0).abs() < 1e-9);

        let s = &run.samples[1];
        assert!(!s.is_ladder);
        assert_eq!(s.name, "Sample X");
        assert_eq!(s.well_number, 13); // B1
        assert_eq!(s.rin, Some(8.4));
        assert_eq!(s.fluorescence.len(), 10);
        // distance axis is reversed/normalized: first row = 1.0, last = 0.1.
        assert!((s.time[0] - 1.0).abs() < 1e-9);

        // ScreenTapeID (→ category) and the per-sample smear region are parsed.
        assert_eq!(s.category, "TAPE-1");
        assert_eq!(s.regions.len(), 1);
        assert!((s.regions[0].lower_length - 200.0).abs() < 1e-9);
        assert!((s.regions[0].upper_length - 400.0).abs() < 1e-9);
    }

    #[test]
    fn calibration_maps_distance_to_size() {
        let mut run = parse_xml(&xml()).unwrap();
        attach_traces(&mut run, &csv());
        calibrate(&mut run);

        // Row 5 sits at distance 0.50 = the ladder's 500 bp band, so both the
        // ladder and the sample (shared markers → shared curve) size there to
        // ~500 bp. The Hyman spline passes exactly through the node.
        assert!(
            (run.samples[0].length[5] - 500.0).abs() < 1.0,
            "ladder: {}",
            run.samples[0].length[5]
        );
        assert!(
            (run.samples[1].length[5] - 500.0).abs() < 1.0,
            "sample: {}",
            run.samples[1].length[5]
        );

        // Points beyond the markers are uncalibrated (NaN), like the Bioanalyzer.
        assert!(
            run.samples[0].length[0].is_nan(),
            "distance 1.0 is past the lower marker"
        );
    }

    #[test]
    fn pair_resolution_accepts_case_variant_electropherogram_suffix() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_tapestation_pair_case_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let xml_path = dir.join("run.xml");
        let csv_path = dir.join("run_electropherogram.csv");
        std::fs::write(&xml_path, xml()).unwrap();
        std::fs::write(&csv_path, csv()).unwrap();

        let (_, found_csv) = resolve_pair(&xml_path).unwrap();

        assert_eq!(found_csv, Some(csv_path));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn pair_resolution_rejects_unsupported_extension() {
        let err = resolve_pair(Path::new("run.txt")).unwrap_err().to_string();

        assert!(
            err.contains("unsupported TapeStation export extension"),
            "got {err}"
        );
    }
}
