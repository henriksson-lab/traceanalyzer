//! Reader for Agilent 2100 Expert **`.xsy` assay** files (the per-chip run
//! method + analysis definition).
//!
//! An `.xsy` is the definition a user copies over to run a given chip kit. See
//! [`docs/xsy_format.md`](../../docs/xsy_format.md) for the container
//! notes. Two on-disk variants exist, both wrapping the same inner
//! `<Chipset>…` XML:
//!
//! * **Compressed** (most kits): an outer XML with a base64 `<compressed_data>`
//!   element. The blob is an **Xceed** container — a small header
//!   (`version`, sizes, the codec string `XceedSCO,1`) followed by a plain
//!   **raw-DEFLATE** payload that inflates to a UTF-16LE XML document. A second
//!   `Do not edit this preview` comment carries a readable base64 preview
//!   (assay title + sample table).
//! * **Plain** (some RNA kits): the inner `<Chipset>…` XML stored directly.
//!
//! This reader extracts the assay name, size unit and calibration **ladder**
//! (which feed [`crate::calibration`]/[`crate::concentration`]) and the raw
//! **run script** values (consumed later by the instrument controller). Deep
//! method fields are left in the XML for now.

use anyhow::{anyhow, bail, Context, Result};
use base64::Engine;
use std::io::Read;
use std::path::Path;

use crate::model::LadderPeak;

/// A parsed `.xsy` assay definition.
#[derive(Debug, Clone, Default)]
pub struct XsyAssay {
    /// Assay/method name, e.g. `"DNA 1000"`.
    pub name: String,
    /// Method comment/description, if present.
    pub comment: String,
    /// Molecular-size unit, e.g. `"bp"` or `"nt"`.
    pub size_unit: String,
    /// Ladder total concentration (`<LadderConcentration>`), if present.
    pub ladder_concentration: Option<f64>,
    /// Calibration ladder peaks (size + areas + concentration).
    pub ladder_peaks: Vec<LadderPeak>,
    /// Sample names from the readable preview table (compressed variant only).
    pub sample_names: Vec<String>,
    /// Decoded run-script values (comma-separated `<ScriptText>` payload). The
    /// column semantics are protocol Phase P2/P3; exposed raw for now.
    pub script_values: Vec<f64>,
    /// `true` if the file used the compressed/Xceed variant.
    pub compressed: bool,
}

/// Reduce an internal build-path `<Name>` to a display name (file stem).
fn clean_name(name: &str) -> String {
    name.rsplit(['\\', '/'])
        .next()
        .unwrap_or(name)
        .trim_end_matches(".xsy")
        .to_string()
}

/// Read and parse an `.xsy` assay file.
pub fn read_xsy_file(path: impl AsRef<Path>) -> Result<XsyAssay> {
    let path = path.as_ref();
    let bytes = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    parse_xsy(&bytes).with_context(|| format!("parsing {}", path.display()))
}

/// Parse an `.xsy` from its raw file bytes.
pub fn parse_xsy(bytes: &[u8]) -> Result<XsyAssay> {
    // The outer file is small text (UTF-8) when compressed, or the big inner XML
    // when plain. Decode leniently to inspect it.
    let outer = String::from_utf8_lossy(bytes);
    let compressed = outer.contains("<compressed_data");

    let (preview_title, sample_names) = extract_preview(&outer);

    let inner_xml = if compressed {
        decompress_body(&outer)?
    } else {
        decode_inner_bytes(bytes)?
    };

    let mut assay = parse_inner(&inner_xml)?;
    assay.compressed = compressed;
    assay.sample_names = sample_names;
    // The readable preview Title is the user-facing assay name; the inner XML's
    // first <Name> is an internal build path, so prefer the title when present,
    // and otherwise clean the path down to its file stem.
    assay.name = match preview_title {
        Some(title) => title,
        None => clean_name(&assay.name),
    };
    Ok(assay)
}

/// Extract, base64-decode and raw-inflate the `<compressed_data>` Xceed blob
/// into the inner XML string.
fn decompress_body(outer: &str) -> Result<String> {
    let start = outer
        .find("<compressed_data")
        .and_then(|i| outer[i..].find('>').map(|j| i + j + 1))
        .ok_or_else(|| anyhow!("no <compressed_data> element"))?;
    let end = outer[start..]
        .find("</compressed_data>")
        .map(|j| start + j)
        .ok_or_else(|| anyhow!("unterminated <compressed_data>"))?;
    let b64: String = outer[start..end].split_whitespace().collect();
    let raw = base64::engine::general_purpose::STANDARD
        .decode(b64.as_bytes())
        .context("base64-decoding compressed_data")?;

    if raw.len() < 20 {
        bail!("Xceed container too short ({} bytes)", raw.len());
    }
    let uncompressed_size = u32::from_le_bytes(raw[8..12].try_into().unwrap()) as usize;
    let codec_units = u32::from_le_bytes(raw[16..20].try_into().unwrap()) as usize;
    let base = 20 + codec_units * 2; // header + UTF-16 codec string ("XceedSCO,1")

    // The raw-DEFLATE payload begins a fixed distance past the codec string.
    // Scan a bounded window for the offset that inflates to exactly the declared
    // size (robust to small header variations across kits).
    for off in base..(base + 128).min(raw.len()) {
        if let Some(text) = try_inflate(&raw[off..], uncompressed_size) {
            return Ok(text);
        }
    }
    bail!("could not locate the DEFLATE stream in the Xceed container")
}

/// Raw-inflate `data`; return the UTF-16LE-decoded text iff it produced exactly
/// `expected` bytes (our validity check for a correct stream offset).
fn try_inflate(data: &[u8], expected: usize) -> Option<String> {
    let mut out = Vec::with_capacity(expected);
    if flate2::read::DeflateDecoder::new(data)
        .read_to_end(&mut out)
        .is_err()
    {
        return None;
    }
    if out.len() != expected {
        return None;
    }
    Some(decode_text(&out))
}

/// Decode a plain-variant `.xsy` (inner XML stored directly) to a string.
fn decode_inner_bytes(bytes: &[u8]) -> Result<String> {
    Ok(decode_text(bytes))
}

/// Decode bytes that are UTF-16LE (with or without BOM) or UTF-8 into a string.
fn decode_text(bytes: &[u8]) -> String {
    if bytes.starts_with(&[0xff, 0xfe]) {
        return encoding_rs::UTF_16LE.decode(&bytes[2..]).0.into_owned();
    }
    // Heuristic: lots of interleaved NULs ⇒ UTF-16LE without BOM.
    let sample = &bytes[..bytes.len().min(64)];
    let nul_ratio = sample.iter().filter(|&&b| b == 0).count() * 2;
    if nul_ratio >= sample.len() {
        encoding_rs::UTF_16LE.decode(bytes).0.into_owned()
    } else {
        String::from_utf8_lossy(bytes).into_owned()
    }
}

/// Parse the inner `<Chipset>…` XML for the fields we surface.
fn parse_inner(xml: &str) -> Result<XsyAssay> {
    let doc = roxmltree::Document::parse(xml).context("parsing inner assay XML")?;
    let root = doc.root_element();

    let text = |tag: &str| -> Option<String> {
        root.descendants()
            .find(|n| n.has_tag_name(tag))
            .and_then(|n| n.text())
            .map(|t| t.trim().to_string())
            .filter(|s| !s.is_empty())
    };

    let name = text("Name").unwrap_or_default();
    let comment = text("Comment").unwrap_or_default();
    let size_unit = text("SizeUnit").unwrap_or_default();
    let ladder_concentration = text("LadderConcentration").and_then(|s| s.parse().ok());

    // The document repeats the ladder across several method blocks. Take the
    // first complete, strictly-increasing run (it resets to the smallest size
    // when the next copy begins).
    let mut ladder_peaks: Vec<LadderPeak> = Vec::new();
    for peak in root
        .descendants()
        .filter(|n| n.has_tag_name("LadderPeak"))
        .filter_map(parse_ladder_peak)
    {
        if ladder_peaks
            .last()
            .is_some_and(|prev| peak.size <= prev.size)
        {
            break;
        }
        ladder_peaks.push(peak);
    }

    let script_values = root
        .descendants()
        .find(|n| n.has_tag_name("ScriptText"))
        .and_then(|n| n.text())
        .map(parse_script)
        .unwrap_or_default();

    Ok(XsyAssay {
        name,
        comment,
        size_unit,
        ladder_concentration,
        ladder_peaks,
        sample_names: Vec::new(),
        script_values,
        compressed: false,
    })
}

fn parse_ladder_peak(node: roxmltree::Node) -> Option<LadderPeak> {
    let child = |tag: &str| -> Option<f64> {
        node.children()
            .find(|c| c.has_tag_name(tag))
            .and_then(|c| c.text())
            .and_then(|t| t.trim().parse().ok())
    };
    Some(LadderPeak {
        size: child("Size")?,
        area_a: child("AreaAType").unwrap_or(f64::NAN),
        area_b: child("AreaBType").unwrap_or(f64::NAN),
        concentration: child("Concentration").unwrap_or(f64::NAN),
    })
}

/// Decode the base64 `<ScriptText>` payload (UTF-16LE CSV) into numbers.
fn parse_script(b64: &str) -> Vec<f64> {
    let cleaned: String = b64.split_whitespace().collect();
    let Ok(raw) = base64::engine::general_purpose::STANDARD.decode(cleaned.as_bytes()) else {
        return Vec::new();
    };
    decode_text(&raw)
        .split(',')
        .filter_map(|v| v.trim().parse().ok())
        .collect()
}

/// Pull the assay title and sample names out of the readable base64 preview
/// comment (compressed variant only). Best-effort: `(None, [])` for plain files.
///
/// The comment is `<!--Do not edit this preview infomation:<guid>:<base64>-->`,
/// where `<base64>` decodes to a small UTF-16 `<Preview>` XML with the title and
/// sample table.
fn extract_preview(outer: &str) -> (Option<String>, Vec<String>) {
    let Some(start) = outer.find("Do not edit this preview") else {
        return (None, Vec::new());
    };
    let tail = &outer[start..];
    let Some(end) = tail.find("-->") else {
        return (None, Vec::new());
    };
    // The base64 is the last ':'-delimited field before the comment close.
    let b64_field = tail[..end].rsplit(':').next().unwrap_or("");
    let cleaned: String = b64_field.split_whitespace().collect();
    let Ok(raw) = base64::engine::general_purpose::STANDARD.decode(cleaned.as_bytes()) else {
        return (None, Vec::new());
    };
    let xml = decode_text(&raw);
    // The preview payload has a short binary prefix before `<?xml`; start at the
    // first tag so roxmltree accepts it.
    let xml = match xml.find('<') {
        Some(i) => &xml[i..],
        None => return (None, Vec::new()),
    };
    let Ok(doc) = roxmltree::Document::parse(xml) else {
        return (None, Vec::new());
    };

    let title = doc
        .descendants()
        .find(|n| n.has_tag_name("Title"))
        .and_then(|n| n.text())
        .map(|t| t.trim().to_string())
        .filter(|s| !s.is_empty());

    // Names are the first <Cell><Value> of each <Row> in <SamplesInfo>.
    let names = doc
        .descendants()
        .filter(|n| n.has_tag_name("Row"))
        .filter_map(|row| {
            row.descendants()
                .find(|n| n.has_tag_name("Value"))
                .and_then(|n| n.text())
                .map(|t| t.trim().to_string())
        })
        .filter(|s| !s.is_empty())
        .collect();

    (title, names)
}
