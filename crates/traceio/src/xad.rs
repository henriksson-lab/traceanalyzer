//! Native Agilent 2100 Bioanalyzer `.xad` container reader.
//!
//! A `.xad` file is a line-oriented text/XML wrapper. The main analytical
//! payload sits in one section as base64 text which, once decoded, is a raw
//! DEFLATE stream framed by a 1-byte header and 9-byte trailer; inflating it
//! yields a UTF-16LE XML document identical to the "Export to XML" output that
//! [`crate::bioanalyzer::parse_xml`] consumes.
//!
//! This is a direct port of grimbough/bioanalyzeR `R/readXAD.R`
//! (`extractCompressed`). That code, and therefore this port, relies on a set
//! of file-position magic constants that were derived from specific sample
//! files. They are collected here as named constants and MUST be validated
//! against real `.xad` files from several 2100 Expert versions (B.01.x and
//! B.02.x); the framing may not be version-stable, and one known variant needs
//! the [`B64_MARKER`] changed from "Oy9" to "Ox9".
//!
//! Status: UNVALIDATED end-to-end — no real `.xad` sample was available at
//! authoring time. The inner-XML parser it feeds *is* validated (against real
//! exported XML). See the crate tests.

use anyhow::{anyhow, Context, Result};
use base64::Engine;
use std::io::Read;

/// The compressed payload lives between the 5th and 6th tag-bearing lines
/// (1-indexed in R: `tagLines[5]` .. `tagLines[6]`).
const TAG_LINE_START: usize = 5;
const TAG_LINE_END: usize = 6;
/// Marker anchoring the true start of the base64 on the first payload line.
const B64_MARKER: &str = "Oy9";
/// Characters of closing markup to drop from the last payload line.
const TAIL_TRIM_CHARS: usize = 18;
/// Bytes of framing to drop after base64-decoding (1 leading, 9 trailing).
const FRAME_HEAD_BYTES: usize = 1;
const FRAME_TAIL_BYTES: usize = 9;

/// Extract and decode the inner XML document from raw `.xad` file bytes.
pub fn extract_inner_xml(raw: &[u8]) -> Result<String> {
    // Read as text lines. The wrapper + base64 are ASCII; lossy UTF-8 is fine
    // for locating tags and slicing base64 (the true payload is UTF-16LE only
    // *after* inflation, handled below).
    let text = String::from_utf8_lossy(raw);
    let lines: Vec<&str> = text.lines().collect();

    // 1-indexed indices of lines that contain a '<' (R: which(grepl("<"))).
    let tag_lines: Vec<usize> = lines
        .iter()
        .enumerate()
        .filter(|(_, l)| l.contains('<'))
        .map(|(i, _)| i + 1)
        .collect();

    if tag_lines.len() < TAG_LINE_END {
        return Err(anyhow!(
            "too few tag lines in .xad ({}, need >= {})",
            tag_lines.len(),
            TAG_LINE_END
        ));
    }

    // Block = originalLines[tagLines[5] : tagLines[6]] (inclusive, 1-indexed).
    let start = tag_lines[TAG_LINE_START - 1]; // 1-indexed line number
    let end = tag_lines[TAG_LINE_END - 1];
    let block: Vec<&str> = lines[start - 1..end].to_vec(); // slice to 0-indexed

    // readXAD.R writes compressedData[2:end] — i.e. it drops the first line of
    // the block and uses the rest.
    if block.len() < 2 {
        return Err(anyhow!("compressed block too short"));
    }
    let mut payload: Vec<String> = block[1..].iter().map(|s| s.to_string()).collect();

    // First payload line: cut everything before the B64 marker (keeping the one
    // character immediately preceding it, matching `regexpr(...)-1`).
    let first = &payload[0];
    let marker_pos = first
        .find(B64_MARKER)
        .ok_or_else(|| anyhow!("base64 marker {:?} not found", B64_MARKER))?;
    let keep_from = marker_pos.saturating_sub(1);
    payload[0] = first[keep_from..].to_string();

    // Last payload line: drop the trailing closing-markup characters.
    let last_idx = payload.len() - 1;
    let last = &payload[last_idx];
    let keep_to = last.len().saturating_sub(TAIL_TRIM_CHARS);
    payload[last_idx] = last[..keep_to].to_string();

    // Concatenate and base64-decode (whitespace/newlines are irrelevant).
    let b64: String = payload.concat();
    let mut bytes = base64::engine::general_purpose::STANDARD
        .decode(b64.as_bytes())
        .context("base64-decoding .xad compressed block")?;

    // Strip the framing bytes: buffer[2 : (length - 9)] in R (1-indexed).
    if bytes.len() <= FRAME_HEAD_BYTES + FRAME_TAIL_BYTES {
        return Err(anyhow!("decoded payload too short to hold framing"));
    }
    let framed = &bytes[FRAME_HEAD_BYTES..bytes.len() - FRAME_TAIL_BYTES];

    // Raw DEFLATE inflate (tinf_uncompress in the reference is raw inflate).
    let inflated = raw_inflate(framed).context("inflating .xad DEFLATE stream")?;
    bytes.clear();

    // Inflated bytes are UTF-16LE text.
    let (decoded, _, had_errors) = encoding_rs::UTF_16LE.decode(&inflated);
    if had_errors {
        return Err(anyhow!("inflated .xad payload is not valid UTF-16LE"));
    }
    Ok(decoded.into_owned())
}

/// Convenience: read a `.xad` file and return the parsed run.
pub fn read_xad_file(path: &std::path::Path) -> Result<crate::model::Electrophoresis> {
    let raw = std::fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    let xml = extract_inner_xml(&raw)?;
    crate::bioanalyzer::parse_xml(&xml)
}

fn raw_inflate(data: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    let mut decoder = flate2::read::DeflateDecoder::new(data);
    decoder.read_to_end(&mut out).context("raw DEFLATE decode")?;
    Ok(out)
}
