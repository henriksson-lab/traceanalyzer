//! Parser for the Agilent 2100 Bioanalyzer inner-XML data model.
//!
//! This XML is what "2100 Expert" writes via *File → Export to XML*, and it is
//! byte-for-byte the same document that lives (base64 + DEFLATE + UTF-16LE
//! compressed) inside a native `.xad` file — see [`crate::xad`]. So this parser
//! serves both the export and the native path.
//!
//! Logic ported from jwfoley/bioanalyzeR `R/bioanalyzer.R` (MIT).

use crate::model::*;
use anyhow::{anyhow, Context, Result};
use base64::Engine;
use roxmltree::{Document, Node};

/// Parse a Bioanalyzer XML document (already decompressed/decoded) into a run.
pub fn parse_xml(xml: &str) -> Result<Electrophoresis> {
    let doc = Document::parse(xml).context("parsing Bioanalyzer XML")?;
    let root = doc.root_element(); // <Chipset>
    let chip = descend(root, &["Chips", "Chip"])
        .ok_or_else(|| anyhow!("no Chips/Chip element (not a Bioanalyzer XML?)"))?;

    let assay = parse_assay_info(chip);
    let ladder_peaks = parse_ladder_peaks(chip);
    let regions = parse_regions(chip);

    let samples_node = descend(chip, &["Files", "File", "Samples"])
        .ok_or_else(|| anyhow!("no Files/File/Samples element"))?;

    let mut samples = Vec::new();
    for sample_node in children_named(samples_node, "Sample") {
        if let Some(sample) = parse_sample(sample_node)? {
            samples.push(sample);
        }
    }

    Ok(Electrophoresis {
        assay,
        ladder_peaks,
        regions,
        samples,
    })
}

fn parse_assay_info(chip: Node) -> AssayInfo {
    let assay_molecular = descend(
        chip,
        &["AssayBody", "DAAssaySetpoints", "DAMAssayInfoMolecular"],
    );
    let concentration_unit = assay_molecular
        .and_then(|n| child_text(n, "ConcentrationUnit"))
        .unwrap_or_default();
    // Hardcode molarity unit from concentration unit, matching bioanalyzeR:
    // otherwise the molecular-weight scales would be wrong.
    let molarity_unit = match concentration_unit.as_str() {
        "ng/µl" => Some("nM".to_string()),
        "pg/µl" => Some("pM".to_string()),
        _ => None,
    };
    let has_upper_marker = descend(
        chip,
        &[
            "AssayBody",
            "DASampleSetpoints",
            "DAMAlignment",
            "Channel",
            "AlignUpperMarker",
        ],
    )
    .and_then(|n| n.text())
    .map(|t| t.trim().eq_ignore_ascii_case("true"))
    .unwrap_or(false);

    AssayInfo {
        file_name: descend(chip, &["Files", "File", "FileInformation", "FileName"])
            .and_then(|n| n.text())
            .unwrap_or_default()
            .to_string(),
        creation_date: descend(chip, &["ChipInformation", "CreationDate"])
            .and_then(|n| n.text())
            .unwrap_or_default()
            .to_string(),
        assay_name: descend(chip, &["AssayHeader", "Title"])
            .and_then(|n| n.text())
            .unwrap_or_default()
            .to_string(),
        assay_type: descend(chip, &["AssayHeader", "Class"])
            .and_then(|n| n.text())
            .unwrap_or_default()
            .to_string(),
        length_unit: assay_molecular
            .and_then(|n| child_text(n, "SizeUnit"))
            .unwrap_or_default(),
        concentration_unit,
        molarity_unit,
        has_upper_marker,
    }
}

fn parse_ladder_peaks(chip: Node) -> Vec<LadderPeak> {
    let Some(ladder) = descend(
        chip,
        &[
            "AssayBody",
            "DAAssaySetpoints",
            "DAMAssayInfoMolecular",
            "LadderPeaks",
        ],
    ) else {
        return Vec::new();
    };
    children_named(ladder, "LadderPeak")
        .map(|p| LadderPeak {
            size: num(child_text(p, "Size")),
            area_a: num(child_text(p, "AreaAType")),
            area_b: num(child_text(p, "AreaBType")),
            concentration: num(child_text(p, "Concentration")),
        })
        .collect()
}

fn parse_regions(chip: Node) -> Vec<Region> {
    let Some(regions) = descend(
        chip,
        &[
            "AssayBody",
            "DASampleSetpoints",
            "DAMSmearAnalysis",
            "Channel",
            "RegionsMolecularSetpoints",
        ],
    ) else {
        return Vec::new();
    };
    regions
        .children()
        .filter(Node::is_element)
        .map(|r| Region {
            lower_length: num(child_text(r, "StartBasePair")),
            upper_length: num(child_text(r, "EndBasePair")),
        })
        .collect()
}

fn parse_sample(sample: Node) -> Result<Option<Sample>> {
    // Skip wells with no acquired data (matches HasData == "true" check).
    if child_text(sample, "HasData")
        .map(|t| !t.trim().eq_ignore_ascii_case("true"))
        .unwrap_or(true)
    {
        return Ok(None);
    }

    let well_number = child_text(sample, "WellNumber")
        .and_then(|t| t.trim().parse::<i32>().ok())
        .unwrap_or(-1);
    let name = child_text(sample, "Name")
        .unwrap_or_default()
        .trim()
        .to_string();
    let category = child_text(sample, "Category")
        .unwrap_or_default()
        .trim()
        .to_string();
    let is_ladder = category == "Ladder";
    let comment = child_text(sample, "Comment")
        .unwrap_or_default()
        .trim()
        .to_string();
    let observations = child_text(sample, "ResultFlagCommonLabel")
        .unwrap_or_default()
        .trim()
        .to_string();

    let rin = descend(sample, &["DAResultStructures", "DARRIN", "Channel", "RIN"])
        .and_then(|n| n.text())
        .and_then(|t| t.trim().parse::<f64>().ok());

    let peaks = descend(
        sample,
        &[
            "DAResultStructures",
            "DARIntegrator",
            "Channel",
            "PeaksMolecular",
        ],
    )
    .map(parse_peaks)
    .unwrap_or_default();

    let (time, fluorescence) = parse_signal(sample)?;

    Ok(Some(Sample {
        well_number,
        name,
        category,
        is_ladder,
        comment,
        observations,
        rin,
        time,
        fluorescence,
        aligned_time: Vec::new(),
        length: Vec::new(),
        concentration: Vec::new(),
        molarity: Vec::new(),
        peaks,
        regions: Vec::new(),
    }))
}

fn parse_peaks(peaks_molecular: Node) -> Vec<Peak> {
    children_named(peaks_molecular, "PeakMolecular")
        .map(|p| Peak {
            observations: child_text(p, "Observations")
                .unwrap_or_default()
                .trim()
                .to_string(),
            length: num(child_text(p, "FragmentSize")),
            time: num(child_text(p, "MigrationTime")),
            aligned_time: num(child_text(p, "AlignedMigrationTime")),
            start_time: num(child_text(p, "StartTime")),
            end_time: num(child_text(p, "EndTime")),
            aligned_start_time: num(child_text(p, "AlignedStartTime")),
            aligned_end_time: num(child_text(p, "AlignedEndTime")),
            area: num(child_text(p, "Area")),
            concentration: num(child_text(p, "Concentration")),
            molarity: num(child_text(p, "Molarity")),
        })
        .collect()
}

/// Read the processed fluorescence trace from the first detector channel.
///
/// Path: `DASignals/DetectorChannels/<first channel>/SignalData`, with the
/// signal itself in `ProcessedSignal` as base64 of little-endian float32.
/// Time axis: `t_i = XStart + XStep * i`.
fn parse_signal(sample: Node) -> Result<(Vec<f64>, Vec<f32>)> {
    let Some(channels) = descend(sample, &["DASignals", "DetectorChannels"]) else {
        return Ok((Vec::new(), Vec::new()));
    };
    let Some(channel) = channels.children().find(Node::is_element) else {
        return Ok((Vec::new(), Vec::new()));
    };
    let Some(signal) = child(channel, "SignalData") else {
        return Ok((Vec::new(), Vec::new()));
    };

    let n = child_text(signal, "NumberOfSamples")
        .and_then(|t| t.trim().parse::<usize>().ok())
        .unwrap_or(0);
    let x_start = num(child_text(signal, "XStart"));
    let x_step = num(child_text(signal, "XStep"));

    let fluorescence = match child_text(signal, "ProcessedSignal") {
        Some(b64) => decode_f32_le(&b64, n)?,
        None => Vec::new(),
    };

    let count = fluorescence.len();
    let time = (0..count).map(|i| x_start + x_step * i as f64).collect();
    Ok((time, fluorescence))
}

/// Decode a base64 string into up to `n` little-endian float32 values.
fn decode_f32_le(b64: &str, n: usize) -> Result<Vec<f32>> {
    let compact: String = b64.split_whitespace().collect();
    let bytes = base64::engine::general_purpose::STANDARD
        .decode(compact.as_bytes())
        .context("base64-decoding ProcessedSignal")?;
    let available = bytes.len() / 4;
    let count = if n == 0 { available } else { n.min(available) };
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let b = [
            bytes[i * 4],
            bytes[i * 4 + 1],
            bytes[i * 4 + 2],
            bytes[i * 4 + 3],
        ];
        out.push(f32::from_le_bytes(b));
    }
    Ok(out)
}

// --- small roxmltree navigation helpers ---------------------------------

/// First direct-child element with the given tag name.
fn child<'a, 'input>(node: Node<'a, 'input>, tag: &str) -> Option<Node<'a, 'input>> {
    node.children()
        .find(|c| c.is_element() && c.has_tag_name(tag))
}

/// Follow a chain of direct-child element names.
fn descend<'a, 'input>(node: Node<'a, 'input>, path: &[&str]) -> Option<Node<'a, 'input>> {
    let mut cur = node;
    for tag in path {
        cur = child(cur, tag)?;
    }
    Some(cur)
}

/// Iterator over direct-child elements with the given tag name.
fn children_named<'a, 'input>(
    node: Node<'a, 'input>,
    tag: &'a str,
) -> impl Iterator<Item = Node<'a, 'input>> {
    node.children()
        .filter(move |c| c.is_element() && c.has_tag_name(tag))
}

/// Concatenated text of a direct-child element.
fn child_text(node: Node, tag: &str) -> Option<String> {
    let c = child(node, tag)?;
    // roxmltree: element text is the first text child; concatenate to be safe.
    let mut s = String::new();
    for t in c.children().filter(|n| n.is_text()) {
        if let Some(txt) = t.text() {
            s.push_str(txt);
        }
    }
    Some(s)
}

/// Parse an optional numeric string, yielding NaN when absent/unparseable.
fn num(s: Option<String>) -> f64 {
    s.and_then(|t| t.trim().parse::<f64>().ok())
        .unwrap_or(f64::NAN)
}
