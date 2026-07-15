//! Writer for user edits back to Agilent Bioanalyzer files.
//!
//! The reader side ([`crate::bioanalyzer`], [`crate::xad`]) is feature-rich; this
//! module does the inverse for the one edit the GUI currently supports: renaming
//! wells. Rather than re-serialize the whole document (which would drop anything
//! we don't model), we **patch in place** — locate each `<Sample>` by its
//! `<WellNumber>` and rewrite only the inner text of its direct-child `<Name>`,
//! leaving every other byte untouched.
//!
//! Supported source encodings (by extension): `.xml` (UTF-8), `.xml.gz` (gzip)
//! and native `.xad` (UTF-16LE inner XML via [`crate::xad::extract_inner_xml`]).
//! Output encoding follows the `dst` extension (`.xml` or `.xml.gz`); writing a
//! real `.xad` container back is out of scope.

use crate::model::Electrophoresis;
use anyhow::{anyhow, Context, Result};
use roxmltree::{Document, Node};
use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::Path;

/// File encodings we can read from / write to.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// Plain UTF-8 XML (`.xml`).
    Xml,
    /// Gzip-compressed UTF-8 XML (`.xml.gz`).
    XmlGz,
    /// Native Bioanalyzer container (`.xad`) — readable, not writable.
    Xad,
}

/// Persist the run to `dst`, using `src` (the originally-loaded file) as the
/// template so nothing but the user's edits changes. Currently the only edit
/// supported is per-well renames: for each sample, the `<Name>` element of the
/// matching `<Sample>` (matched by `<WellNumber>`) is rewritten to
/// `run.samples[i].name`.
pub fn save_run(run: &Electrophoresis, src: &Path, dst: &Path) -> Result<()> {
    let dst_kind = classify(dst)
        .ok_or_else(|| anyhow!("unsupported output extension for {}", dst.display()))?;
    if dst_kind == Kind::Xad {
        return Err(anyhow!(
            "saving back to a native .xad container is not yet supported; \
             save as .xml (or .xml.gz) instead"
        ));
    }

    // Decode the template into a UTF-8 XML string we can splice.
    let xml = read_xml(src)?;
    let patched = patch_names(&xml, run)?;

    write_xml(dst, dst_kind, &patched)
}

/// Read `src` and return its XML as a UTF-8 string, decoding per extension.
fn read_xml(src: &Path) -> Result<String> {
    let kind = classify(src)
        .ok_or_else(|| anyhow!("unsupported source extension for {}", src.display()))?;
    let raw = std::fs::read(src).with_context(|| format!("reading {}", src.display()))?;
    match kind {
        Kind::Xml => {
            String::from_utf8(raw).with_context(|| format!("{} is not valid UTF-8", src.display()))
        }
        Kind::XmlGz => {
            let mut d = flate2::read::GzDecoder::new(&raw[..]);
            let mut s = String::new();
            d.read_to_string(&mut s)
                .with_context(|| format!("gunzipping {}", src.display()))?;
            Ok(s)
        }
        // The native container's inner XML is UTF-16LE; extract_inner_xml
        // returns it already decoded to a Rust (UTF-8) String.
        Kind::Xad => crate::xad::extract_inner_xml(&raw)
            .with_context(|| format!("extracting inner XML from {}", src.display())),
    }
}

/// Write `xml` to `dst` in the requested encoding.
fn write_xml(dst: &Path, kind: Kind, xml: &str) -> Result<()> {
    match kind {
        Kind::Xml => {
            std::fs::write(dst, xml.as_bytes())
                .with_context(|| format!("writing {}", dst.display()))?;
        }
        Kind::XmlGz => {
            let file =
                std::fs::File::create(dst).with_context(|| format!("creating {}", dst.display()))?;
            let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            enc.write_all(xml.as_bytes())
                .with_context(|| format!("gzip-writing {}", dst.display()))?;
            enc.finish()
                .with_context(|| format!("finishing gzip {}", dst.display()))?;
        }
        Kind::Xad => unreachable!("save-as-.xad is rejected earlier"),
    }
    Ok(())
}

/// Classify a path by its (case-insensitive) extension. Returns `None` for
/// anything we neither read nor write.
fn classify(path: &Path) -> Option<Kind> {
    let name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    if name.ends_with(".xml.gz") {
        Some(Kind::XmlGz)
    } else if name.ends_with(".xad") {
        Some(Kind::Xad)
    } else if name.ends_with(".xml") {
        Some(Kind::Xml)
    } else {
        None
    }
}

/// A single byte-range replacement in the source string.
struct Edit {
    start: usize,
    end: usize,
    text: String,
}

/// Rewrite the `<Name>` text of every `<Sample>` whose `<WellNumber>` matches a
/// sample in `run`, leaving all other bytes intact.
fn patch_names(xml: &str, run: &Electrophoresis) -> Result<String> {
    // well number -> desired name.
    let wanted: HashMap<i32, &str> = run
        .samples
        .iter()
        .map(|s| (s.well_number, s.name.as_str()))
        .collect();

    // Collect edits with the parsed document borrowing `xml`, then drop the
    // document so we can mutate the owned string.
    let edits = {
        let doc = Document::parse(xml).context("parsing template XML")?;
        collect_name_edits(&doc, &wanted)
    };

    // Apply from the end so earlier byte offsets stay valid.
    let mut out = xml.to_string();
    let mut edits = edits;
    edits.sort_by(|a, b| b.start.cmp(&a.start));
    for e in edits {
        out.replace_range(e.start..e.end, &e.text);
    }
    Ok(out)
}

/// Find, for each matching `<Sample>`, the byte range whose replacement renames
/// its direct-child `<Name>`.
fn collect_name_edits(doc: &Document, wanted: &HashMap<i32, &str>) -> Vec<Edit> {
    let mut edits = Vec::new();
    // A well is a `<Sample>` element directly under a `<Samples>` element (the
    // detector channels also carry `<Name>`, but they are nested deeper).
    for sample in doc.descendants().filter(|n| {
        n.has_tag_name("Sample")
            && n.parent().map(|p| p.has_tag_name("Samples")).unwrap_or(false)
    }) {
        let Some(well) = child_text(sample, "WellNumber")
            .and_then(|t| t.trim().parse::<i32>().ok())
        else {
            continue;
        };
        let Some(new_name) = wanted.get(&well) else {
            continue;
        };
        let Some(name_elem) = child(sample, "Name") else {
            continue; // nothing to rename
        };
        if let Some(edit) = name_edit(name_elem, escape_xml(new_name)) {
            edits.push(edit);
        }
    }
    edits
}

/// Build the edit that replaces `name_elem`'s inner text with `escaped`.
fn name_edit(name_elem: Node, escaped: String) -> Option<Edit> {
    // Preferred: replace the existing text node in place (robust to attributes).
    if let Some(tc) = name_elem.children().find(|c| c.is_text()) {
        let r = tc.range();
        return Some(Edit {
            start: r.start,
            end: r.end,
            text: escaped,
        });
    }

    // Empty (`<Name></Name>`) or self-closing (`<Name/>`) element: fall back to
    // the element's own range and insert/replace accordingly.
    let r = name_elem.range();
    let slice = name_elem.document().input_text().get(r.start..r.end)?;
    let gt = slice.find('>')?; // end of the opening tag
    let self_closing = slice.as_bytes().get(gt.wrapping_sub(1)) == Some(&b'/');
    if self_closing {
        // Replace `<Name/>` wholesale with `<Name>text</Name>`.
        Some(Edit {
            start: r.start,
            end: r.end,
            text: format!("<Name>{escaped}</Name>"),
        })
    } else {
        // Insert between `<Name>` and `</Name>` (inner span is empty).
        let inner_start = r.start + gt + 1;
        let inner_end = r.end - "</Name>".len();
        Some(Edit {
            start: inner_start,
            end: inner_end,
            text: escaped,
        })
    }
}

/// Escape the XML text-content metacharacters (`&`, `<`, `>`). `&` must go first.
fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

// --- small roxmltree navigation helpers (mirroring bioanalyzer.rs) ----------

/// First direct-child element with the given tag name.
fn child<'a, 'input>(node: Node<'a, 'input>, tag: &str) -> Option<Node<'a, 'input>> {
    node.children()
        .find(|c| c.is_element() && c.has_tag_name(tag))
}

/// Concatenated text of a direct-child element.
fn child_text(node: Node, tag: &str) -> Option<String> {
    let c = child(node, tag)?;
    Some(c.children().filter_map(|t| t.text()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escape_orders_ampersand_first() {
        assert_eq!(escape_xml("A & B <x>"), "A &amp; B &lt;x&gt;");
        assert_eq!(escape_xml("plain"), "plain");
    }

    #[test]
    fn renames_only_the_direct_child_name() {
        // Minimal document mirroring the real nesting: a detector-channel <Name>
        // lives deeper inside the same <Sample> and must NOT be touched.
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample key=\"0\"><Name>Ladder 1</Name><WellNumber>1</WellNumber>\
            <DASignals><DetectorChannels><Channel><SignalData>\
            <Name>Blue Fluorescence</Name></SignalData></Channel></DetectorChannels></DASignals>\
            </Sample>\
            <Sample key=\"1\"><Name>old</Name><WellNumber>2</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";

        let mut run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![],
        };
        run.samples.push(make_sample(1, "New Ladder & Co <1>"));
        run.samples.push(make_sample(2, "second"));

        let out = patch_names(xml, &run).unwrap();
        assert!(out.contains("<Name>New Ladder &amp; Co &lt;1&gt;</Name>"));
        assert!(out.contains("<Name>second</Name>"));
        // The detector-channel name is untouched.
        assert!(out.contains("<Name>Blue Fluorescence</Name>"));
        assert!(!out.contains("<Name>Ladder 1</Name>"));
        assert!(!out.contains("<Name>old</Name>"));
    }

    #[test]
    fn fills_an_empty_name_element() {
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name></Name><WellNumber>3</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(3, "filled")],
        };
        let out = patch_names(xml, &run).unwrap();
        assert!(out.contains("<Name>filled</Name>"), "got: {out}");
    }

    fn make_sample(well: i32, name: &str) -> crate::model::Sample {
        crate::model::Sample {
            well_number: well,
            name: name.to_string(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![],
            fluorescence: vec![],
            aligned_time: vec![],
            length: vec![],
            concentration: vec![],
            molarity: vec![],
            peaks: vec![],
        }
    }
}
