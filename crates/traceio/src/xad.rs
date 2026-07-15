//! Native Agilent 2100 Bioanalyzer `.xad` container reader.
//!
//! A `.xad` is a text/XML wrapper whose analytical payload sits in a
//! `<compressed_data … dt:dt="bin.base64">…</compressed_data>` element. That
//! base64 decodes to a compressed blob which, once inflated, is a UTF-16LE XML
//! document identical to the *File → Export to XML* output that
//! [`crate::bioanalyzer::parse_xml`] consumes.
//!
//! Two container variants are handled — see `docs/xad_format.md` for the full
//! specification and how these findings differ from the original
//! grimbough/bioanalyzeR `readXAD.R`:
//!
//! * **Xceed** (observed in current 2100 Expert output): the decoded blob is a
//!   16-byte header (`[2, 79, uncompressed_len, compressed_len]` as LE u32) + a
//!   length-prefixed `"Xceed"` UTF-16 marker + padding, then a **raw DEFLATE**
//!   stream, then a 9-byte trailer. We locate and validate the stream from the
//!   header's self-described lengths.
//! * **Legacy** (original `readXAD.R`): strip 1 leading + 9 trailing bytes from
//!   the decoded blob and raw-inflate. Kept as a fallback; not validated here
//!   against a real legacy file.

use anyhow::{anyhow, Context, Result};
use base64::Engine;
use std::io::Read;

/// Bytes of trailer after the DEFLATE stream in both observed variants.
const FRAME_TAIL_BYTES: usize = 9;
/// Leading bytes stripped in the legacy (pre-Xceed) variant.
const LEGACY_FRAME_HEAD_BYTES: usize = 1;
/// `"Xceed"` in UTF-16LE, the marker identifying the Xceed container.
const XCEED_MARKER: &[u8] = b"X\x00c\x00e\x00e\x00d\x00";

/// Extract and decode the inner XML document from raw `.xad` file bytes.
pub fn extract_inner_xml(raw: &[u8]) -> Result<String> {
    let blob = decode_compressed_data(raw)?;
    let inflated = decompress_payload(&blob)?;

    // Inflated bytes are UTF-16LE XML (no BOM; starts directly with "<Chipset").
    let (decoded, _, had_errors) = encoding_rs::UTF_16LE.decode(&inflated);
    if had_errors {
        return Err(anyhow!("inflated .xad payload is not valid UTF-16LE"));
    }
    Ok(decoded.into_owned())
}

/// Convenience: read a `.xad` file and return the parsed run.
///
/// Note: current `.xad` files hold **raw acquisition data + sample metadata**,
/// not processed per-well results (see `docs/xad_format.md` §5). The returned
/// [`crate::model::Electrophoresis`] therefore has assay info, the defined
/// ladder, and sample metadata, but per-sample traces/peaks will be empty —
/// those are recomputed by 2100 Expert on open. Use [`read_xad_raw_channels`]
/// to get the raw detector electropherograms.
pub fn read_xad_file(path: &std::path::Path) -> Result<crate::model::Electrophoresis> {
    let raw = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let xml = extract_inner_xml(&raw)?;
    crate::bioanalyzer::parse_xml(&xml)
}

/// A raw detector channel from a `.xad`: the whole-chip continuous acquisition
/// (all wells in time order), before the software splits it per well.
#[derive(Debug, Clone)]
pub struct RawChannel {
    /// e.g. "BlueFluorescence", "RedFluorescence".
    pub channel_id: String,
    pub name: String,
    pub x_start: f64,
    pub x_step: f64,
    /// Raw signal samples (little-endian float32 in the file).
    pub signal: Vec<f32>,
}

/// Read the raw fluorescence detector channels from a native `.xad`.
///
/// Returns the `Chip/RawSignals/DetectorChannels/Channel` signals (the
/// electropherogram detectors, e.g. Blue/Red), skipping housekeeping channels
/// (temperature/voltage/current/position) and any empty channel.
pub fn read_xad_raw_channels(path: &std::path::Path) -> Result<Vec<RawChannel>> {
    let raw = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let xml = extract_inner_xml(&raw)?;
    let doc = roxmltree::Document::parse(&xml).context("parsing .xad inner XML")?;

    let mut channels = Vec::new();
    for sd in doc.descendants().filter(|n| n.has_tag_name("SignalData")) {
        // Only detector channels: parent chain …/DetectorChannels/Channel/SignalData.
        let is_detector = sd
            .parent()
            .and_then(|c| c.parent())
            .map(|dc| dc.has_tag_name("DetectorChannels"))
            .unwrap_or(false);
        if !is_detector {
            continue;
        }
        let signal = match child_text(sd, "RawSignal") {
            // NumberOfSamples is unreliable in raw files; decode all floats.
            Some(b64) => decode_f32_le(&b64)?,
            None => continue,
        };
        if signal.is_empty() {
            continue;
        }
        channels.push(RawChannel {
            channel_id: child_text(sd, "ChannelID").unwrap_or_default(),
            name: child_text(sd, "Name").unwrap_or_default(),
            x_start: child_text(sd, "XStart")
                .and_then(|t| t.trim().parse().ok())
                .unwrap_or(0.0),
            x_step: child_text(sd, "XStep")
                .and_then(|t| t.trim().parse().ok())
                .unwrap_or(1.0),
            signal,
        });
    }
    Ok(channels)
}

/// Concatenated text of a direct-child element.
fn child_text(node: roxmltree::Node, tag: &str) -> Option<String> {
    let c = node.children().find(|n| n.is_element() && n.has_tag_name(tag))?;
    Some(c.children().filter_map(|t| t.text()).collect())
}

/// Decode a base64 string into little-endian float32 values (all available).
fn decode_f32_le(b64: &str) -> Result<Vec<f32>> {
    let compact: String = b64.split_whitespace().collect();
    let bytes = base64::engine::general_purpose::STANDARD
        .decode(compact.as_bytes())
        .context("base64-decoding RawSignal")?;
    Ok(bytes
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect())
}

/// Locate the `<compressed_data …>…</compressed_data>` element and base64-decode
/// its contents. Falls back to the legacy line-slicing scheme if the element is
/// absent.
fn decode_compressed_data(raw: &[u8]) -> Result<Vec<u8>> {
    let b64 = match find_between(raw, b"<compressed_data", b"</compressed_data>") {
        Some(inner) => {
            // Skip past the opening tag's '>' to the base64 body.
            let gt = inner
                .iter()
                .position(|&b| b == b'>')
                .ok_or_else(|| anyhow!("malformed <compressed_data> tag"))?;
            &inner[gt + 1..]
        }
        None => return legacy_line_extract(raw),
    };
    // base64 ignoring embedded whitespace/newlines.
    let compact: Vec<u8> = b64
        .iter()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect();
    base64::engine::general_purpose::STANDARD
        .decode(&compact)
        .context("base64-decoding <compressed_data>")
}

/// Decompress the decoded blob, dispatching on the container variant.
fn decompress_payload(blob: &[u8]) -> Result<Vec<u8>> {
    if is_xceed(blob) {
        decompress_xceed(blob)
    } else {
        // Legacy: buffer[2 : len-9] in readXAD.R (1-indexed) -> raw inflate.
        if blob.len() <= LEGACY_FRAME_HEAD_BYTES + FRAME_TAIL_BYTES {
            return Err(anyhow!("decoded payload too short for legacy framing"));
        }
        let framed = &blob[LEGACY_FRAME_HEAD_BYTES..blob.len() - FRAME_TAIL_BYTES];
        raw_inflate(framed).context("inflating legacy .xad DEFLATE stream")
    }
}

/// True if the blob carries the Xceed header (`"Xceed"` marker at offset 20).
fn is_xceed(blob: &[u8]) -> bool {
    blob.len() >= 30 && &blob[20..20 + XCEED_MARKER.len()] == XCEED_MARKER
}

/// Xceed variant: use the header's self-described lengths to locate the raw
/// DEFLATE stream, then validate the inflated size.
fn decompress_xceed(blob: &[u8]) -> Result<Vec<u8>> {
    let uncompressed_len = read_u32_le(blob, 8)? as usize;
    let compressed_len = read_u32_le(blob, 12)? as usize;

    // The DEFLATE stream is `compressed_len` bytes ending `FRAME_TAIL_BYTES`
    // before the end of the blob.
    let end = blob
        .len()
        .checked_sub(FRAME_TAIL_BYTES)
        .ok_or_else(|| anyhow!("Xceed blob shorter than its trailer"))?;
    let start = end
        .checked_sub(compressed_len)
        .ok_or_else(|| anyhow!("Xceed compressed_len {compressed_len} exceeds blob"))?;
    if start < 16 {
        return Err(anyhow!("Xceed DEFLATE stream overlaps the header"));
    }

    let inflated = raw_inflate(&blob[start..end]).context("inflating Xceed DEFLATE stream")?;
    if inflated.len() != uncompressed_len {
        return Err(anyhow!(
            "Xceed size mismatch: header says {uncompressed_len} bytes, inflated {}",
            inflated.len()
        ));
    }
    Ok(inflated)
}

/// Legacy `readXAD.R` extraction: locate the payload between the 5th and 6th
/// lines that contain a '<', then trim to the base64 body. Preserved for older
/// files that do not use the `<compressed_data>` element.
fn legacy_line_extract(raw: &[u8]) -> Result<Vec<u8>> {
    const B64_MARKER: &str = "Oy9";
    const TAIL_TRIM_CHARS: usize = 18;
    let text = String::from_utf8_lossy(raw);
    let lines: Vec<&str> = text.lines().collect();
    let tag_lines: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, l)| l.contains('<'))
        .map(|(i, _)| i)
        .collect();
    if tag_lines.len() < 6 {
        return Err(anyhow!(
            "no <compressed_data> element and too few tag lines ({}) for legacy extraction",
            tag_lines.len()
        ));
    }
    let block = &lines[tag_lines[4]..=tag_lines[5]];
    if block.len() < 2 {
        return Err(anyhow!("legacy compressed block too short"));
    }
    let mut payload: Vec<String> = block[1..].iter().map(|s| s.to_string()).collect();
    let marker_pos = payload[0]
        .find(B64_MARKER)
        .ok_or_else(|| anyhow!("legacy base64 marker {B64_MARKER:?} not found"))?;
    payload[0] = payload[0][marker_pos.saturating_sub(1)..].to_string();
    let last = payload.len() - 1;
    let keep_to = payload[last].len().saturating_sub(TAIL_TRIM_CHARS);
    payload[last] = payload[last][..keep_to].to_string();
    let b64: String = payload.concat();
    base64::engine::general_purpose::STANDARD
        .decode(b64.as_bytes())
        .context("base64-decoding legacy .xad block")
}

/// Find the bytes strictly between the first `start` marker and the following
/// `end` marker (the returned slice starts at `start`, i.e. includes the opening
/// tag, so the caller can skip past its attributes).
fn find_between<'a>(hay: &'a [u8], start: &[u8], end: &[u8]) -> Option<&'a [u8]> {
    let s = find_sub(hay, start)?;
    let e = find_sub(&hay[s..], end)? + s;
    Some(&hay[s..e])
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

fn read_u32_le(b: &[u8], off: usize) -> Result<u32> {
    let bytes = b
        .get(off..off + 4)
        .ok_or_else(|| anyhow!("truncated header at offset {off}"))?;
    Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
}

fn raw_inflate(data: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut decoder = flate2::read::DeflateDecoder::new(data);
    decoder.read_to_end(&mut out).context("raw DEFLATE decode")?;
    Ok(out)
}
