//! Writer for user edits back to instrument files.
//!
//! The reader side ([`crate::bioanalyzer`], [`crate::xad`]) is feature-rich; this
//! module does the inverse for the one edit the GUI currently supports: renaming
//! wells. Rather than re-serialize the whole document (which would drop anything
//! we don't model), we **patch in place** — locate each `<Sample>` by its
//! `<WellNumber>` and rewrite only the inner text of its direct-child `<Name>`,
//! leaving every other byte untouched.
//!
//! Bioanalyzer source encodings (by extension): `.xml` (UTF-8), `.xml.gz`
//! (gzip), and native `.xad` (decoded through the container's inner XML).
//! Output encoding follows the `dst` extension (`.xml` or `.xml.gz`); writing a
//! real `.xad` container back is out of scope. Fragment Analyzer `.raw` runs save
//! in place by patching the `.txt` sidecar's `Sample ID:` values.

use crate::model::Electrophoresis;
use anyhow::{anyhow, Context, Result};
use roxmltree::{Document, Node};
use std::collections::{HashMap, HashSet};
use std::io::{Read, Write};
use std::path::Path;

/// Persist edits from a [`crate::io::LoadedRun`] to `dst`, using its source path
/// as the template.
///
/// This is the path-oriented counterpart to [`crate::io::read_path`]. Use
/// [`save_run`] when the source path and run are tracked separately.
pub fn save_path(loaded: &crate::io::LoadedRun, dst: impl AsRef<Path>) -> Result<()> {
    save_run(&loaded.run, &loaded.source.path, dst.as_ref())
}

/// Coarse preflight for [`save_path`].
///
/// Returns `true` only for destination shapes that the high-level save API is
/// expected to support for this loaded source. The actual save can still fail on
/// malformed inputs, missing sidecars, permissions, or unmatched sample names.
pub fn supports_save_path(loaded: &crate::io::LoadedRun, dst: impl AsRef<Path>) -> bool {
    let dst = dst.as_ref();
    match loaded.source.format {
        crate::io::TraceFormat::BioanalyzerXml | crate::io::TraceFormat::BioanalyzerXmlGz => {
            matches!(classify(dst), Some(Kind::Xml | Kind::XmlGz))
        }
        crate::io::TraceFormat::BioanalyzerXad => {
            matches!(classify(dst), Some(Kind::Xml | Kind::XmlGz))
        }
        crate::io::TraceFormat::FragmentAnalyzerRun { .. } => {
            crate::fa::is_saveable_fa_destination_for_source(&loaded.source.path, dst)
        }
        crate::io::TraceFormat::TapeStationExport { .. } => false,
    }
}

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
    if crate::fa::is_fa_path(src) {
        if !same_path(src, dst) && !crate::fa::is_saveable_fa_destination_for_source(src, dst) {
            return Err(anyhow!(
                "Fragment Analyzer runs can only be saved in place or to their .txt sidecar"
            ));
        }
        return crate::fa::save_txt_names(dst, run);
    }
    if crate::tapestation::is_tapestation_path(src) {
        return Err(anyhow!(
            "TapeStation exports are read-only; sample renames cannot be saved back"
        ));
    }

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
    // TapeStation exports share the `.xml` extension but a different schema (the
    // Bioanalyzer patcher would silently no-op); they are read-only for now.
    if crate::tapestation::looks_like_tapestation_xml(&xml) {
        return Err(anyhow!(
            "TapeStation exports are read-only; sample renames cannot be saved back"
        ));
    }
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

fn same_path(a: &Path, b: &Path) -> bool {
    let comparable =
        |path: &Path| std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
    comparable(a) == comparable(b)
}

/// Write `xml` to `dst` in the requested encoding.
fn write_xml(dst: &Path, kind: Kind, xml: &str) -> Result<()> {
    let mut tmp = TempOutput::create_next_to(dst)?;
    match kind {
        Kind::Xml => {
            tmp.file
                .as_mut()
                .expect("temporary output file should be present while writing")
                .write_all(xml.as_bytes())
                .with_context(|| format!("writing {}", tmp.path.display()))?;
        }
        Kind::XmlGz => {
            let mut enc = flate2::write::GzEncoder::new(
                tmp.file
                    .as_mut()
                    .expect("temporary output file should be present while writing"),
                flate2::Compression::default(),
            );
            enc.write_all(xml.as_bytes())
                .with_context(|| format!("gzip-writing {}", tmp.path.display()))?;
            enc.finish()
                .with_context(|| format!("finishing gzip {}", tmp.path.display()))?;
        }
        Kind::Xad => unreachable!("save-as-.xad is rejected earlier"),
    }
    tmp.persist(dst)?;
    Ok(())
}

struct TempOutput {
    path: std::path::PathBuf,
    file: Option<std::fs::File>,
}

impl TempOutput {
    fn create_next_to(dst: &Path) -> Result<Self> {
        let parent = destination_parent(dst);
        let name = dst
            .file_name()
            .ok_or_else(|| anyhow!("destination has no file name: {}", dst.display()))?
            .to_string_lossy();
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        for attempt in 0..100u32 {
            let path = parent.join(format!(
                ".{name}.tmp-{}-{nonce}-{attempt}",
                std::process::id()
            ));
            match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
            {
                Ok(file) => {
                    return Ok(Self {
                        path,
                        file: Some(file),
                    });
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => {
                    return Err(e).with_context(|| {
                        format!("creating temporary output next to {}", dst.display())
                    });
                }
            }
        }
        Err(anyhow!(
            "could not create a unique temporary output next to {}",
            dst.display()
        ))
    }

    fn persist(mut self, dst: &Path) -> Result<()> {
        let file = self
            .file
            .take()
            .expect("temporary output file should be present before persist");
        preserve_existing_permissions(&self.path, dst)?;
        file.sync_all()
            .with_context(|| format!("fsyncing temporary output {}", self.path.display()))?;
        drop(file);
        replace_file(&self.path, dst).with_context(|| {
            format!(
                "replacing {} with temporary output {}",
                dst.display(),
                self.path.display()
            )
        })?;
        let parent = destination_parent(dst);
        sync_destination_dir(parent)?;
        Ok(())
    }
}

fn preserve_existing_permissions(tmp: &Path, dst: &Path) -> Result<()> {
    match std::fs::metadata(dst) {
        Ok(meta) => std::fs::set_permissions(tmp, meta.permissions()).with_context(|| {
            format!(
                "preserving permissions from {} on temporary file {}",
                dst.display(),
                tmp.display()
            )
        }),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e).with_context(|| format!("reading permissions for {}", dst.display())),
    }
}

#[cfg(unix)]
fn sync_destination_dir(parent: &Path) -> Result<()> {
    std::fs::File::open(parent)
        .and_then(|dir| dir.sync_all())
        .with_context(|| format!("fsyncing destination directory {}", parent.display()))
}

#[cfg(not(unix))]
fn sync_destination_dir(_parent: &Path) -> Result<()> {
    Ok(())
}

#[cfg(windows)]
fn replace_file(tmp: &Path, dst: &Path) -> std::io::Result<()> {
    use std::os::windows::ffi::OsStrExt;

    const MOVEFILE_REPLACE_EXISTING: u32 = 0x1;
    const MOVEFILE_WRITE_THROUGH: u32 = 0x8;

    #[link(name = "kernel32")]
    extern "system" {
        fn MoveFileExW(
            lpExistingFileName: *const u16,
            lpNewFileName: *const u16,
            dwFlags: u32,
        ) -> i32;
    }

    fn wide(path: &Path) -> Vec<u16> {
        path.as_os_str()
            .encode_wide()
            .chain(std::iter::once(0))
            .collect()
    }

    let tmp = wide(tmp);
    let dst = wide(dst);
    let ok = unsafe {
        MoveFileExW(
            tmp.as_ptr(),
            dst.as_ptr(),
            MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH,
        )
    };
    if ok == 0 {
        Err(std::io::Error::last_os_error())
    } else {
        Ok(())
    }
}

#[cfg(not(windows))]
fn replace_file(tmp: &Path, dst: &Path) -> std::io::Result<()> {
    std::fs::rename(tmp, dst)
}

impl Drop for TempOutput {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

fn destination_parent(dst: &Path) -> &Path {
    dst.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
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
    reject_duplicate_bioanalyzer_wells(run)?;

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
    if !run.samples.is_empty() {
        if edits.is_empty() {
            return Err(anyhow!(
                "template XML did not contain any matching Bioanalyzer sample names to patch"
            ));
        }
        if edits.len() != wanted.len() {
            return Err(anyhow!(
                "template XML matched {} of {} Bioanalyzer sample wells",
                edits.len(),
                wanted.len()
            ));
        }
    }

    // Apply from the end so earlier byte offsets stay valid.
    let mut out = xml.to_string();
    let mut edits = edits;
    edits.sort_by(|a, b| b.start.cmp(&a.start));
    for e in edits {
        out.replace_range(e.start..e.end, &e.text);
    }
    Ok(out)
}

fn reject_duplicate_bioanalyzer_wells(run: &Electrophoresis) -> Result<()> {
    let mut seen = HashSet::with_capacity(run.samples.len());
    if run
        .samples
        .iter()
        .any(|sample| !seen.insert(sample.well_number))
    {
        return Err(anyhow!(
            "loaded Bioanalyzer run has duplicate sample well numbers; refusing partial XML save"
        ));
    }
    Ok(())
}

/// Find, for each matching `<Sample>`, the byte range whose replacement renames
/// its direct-child `<Name>`.
fn collect_name_edits(doc: &Document, wanted: &HashMap<i32, &str>) -> Vec<Edit> {
    let mut edits = Vec::new();
    // A well is a `<Sample>` element directly under a `<Samples>` element (the
    // detector channels also carry `<Name>`, but they are nested deeper).
    for sample in doc.descendants().filter(|n| {
        n.has_tag_name("Sample")
            && n.parent()
                .map(|p| p.has_tag_name("Samples"))
                .unwrap_or(false)
    }) {
        let Some(well) =
            child_text(sample, "WellNumber").and_then(|t| t.trim().parse::<i32>().ok())
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
    let self_closing = slice[..gt].trim_end().ends_with('/');
    if self_closing {
        // Replace `<Name/>` / `<Name attr="..."/>` while preserving attributes.
        let opening = slice[..gt].trim_end().trim_end_matches('/').trim_end();
        Some(Edit {
            start: r.start,
            end: r.end,
            text: format!("{opening}>{escaped}</Name>"),
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

    #[test]
    fn fills_self_closing_name_elements_with_spacing_or_attributes() {
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name /><WellNumber>3</WellNumber></Sample>\
            <Sample><Name kind=\"sample\" /><WellNumber>4</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(3, "first"), make_sample(4, "second")],
        };

        let out = patch_names(xml, &run).unwrap();

        assert!(out.contains("<Name>first</Name>"), "got: {out}");
        assert!(
            out.contains("<Name kind=\"sample\">second</Name>"),
            "attributes on empty Name are preserved: {out}"
        );
    }

    #[test]
    fn bioanalyzer_patch_errors_when_no_sample_names_match() {
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>9</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(7, "new")],
        };

        let err = patch_names(xml, &run).unwrap_err().to_string();

        assert!(err.contains("did not contain any matching"), "got: {err}");
    }

    #[test]
    fn bioanalyzer_patch_errors_when_only_some_sample_names_match() {
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "new"), make_sample(2, "missing")],
        };

        let err = patch_names(xml, &run).unwrap_err().to_string();

        assert!(err.contains("matched 1 of 2"), "got: {err}");
    }

    #[test]
    fn bioanalyzer_patch_rejects_duplicate_sample_wells() {
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "new"), make_sample(1, "duplicate")],
        };

        let err = patch_names(xml, &run).unwrap_err().to_string();

        assert!(err.contains("duplicate sample well numbers"), "got: {err}");
        assert!(err.contains("refusing partial XML save"), "got: {err}");
    }

    #[test]
    fn bioanalyzer_duplicate_wells_do_not_overwrite_in_place() {
        let dir = unique_temp_dir("traceio_bioanalyzer_duplicate_save");
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xml");
        let original = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        std::fs::write(&src, original).unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "new"), make_sample(1, "duplicate")],
        };

        let err = save_run(&run, &src, &src).unwrap_err().to_string();

        assert!(err.contains("duplicate sample well numbers"), "got: {err}");
        assert_eq!(std::fs::read_to_string(&src).unwrap(), original);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn bioanalyzer_xml_save_in_place_uses_replacement_output() {
        let dir = unique_temp_dir("traceio_bioanalyzer_xml_in_place");
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xml");
        std::fs::write(
            &src,
            "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>",
        )
        .unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "new")],
        };

        save_run(&run, &src, &src).unwrap();

        let patched = std::fs::read_to_string(&src).unwrap();
        assert!(patched.contains("<Name>new</Name>"), "got: {patched}");
        assert!(!patched.contains("<Name>old</Name>"));
        assert_no_bioanalyzer_temp_files(&dir);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn bioanalyzer_xml_gz_save_in_place_uses_replacement_output() {
        let dir = unique_temp_dir("traceio_bioanalyzer_xml_gz_in_place");
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xml.gz");
        let xml = "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>";
        write_gz_xml(&src, xml);
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "new")],
        };

        save_run(&run, &src, &src).unwrap();

        let patched = read_gz_xml(&src);
        assert!(patched.contains("<Name>new</Name>"), "got: {patched}");
        assert!(!patched.contains("<Name>old</Name>"));
        assert_no_bioanalyzer_temp_files(&dir);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_raw_save_patches_txt_sidecar_in_place() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_save_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        let txt = dir.join("run.txt");
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(
            &txt,
            "Capillary #: 1\nWell: D1\nSample ID: old\nCapillary #: 2\nWell: D2\nSample ID: old2\n",
        )
        .unwrap();

        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "D1: renamed"), make_sample(2, "plain")],
        };

        save_run(&run, &raw, &raw).unwrap();

        let patched = std::fs::read_to_string(&txt).unwrap();
        assert!(patched.contains("Sample ID: renamed"));
        assert!(patched.contains("Sample ID: plain"));
        assert!(!patched.contains("Sample ID: old\n"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_save_allows_exact_sidecar_source_as_in_place_destination() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_save_source_identity_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        let pks = dir.join("run.PKS");
        let txt = dir.join("run.txt");
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(&pks, b"sidecar").unwrap();
        std::fs::write(&txt, "Capillary #: 1\nWell: D1\nSample ID: old\n").unwrap();

        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "D1: renamed")],
        };

        save_run(&run, &pks, &pks).unwrap();

        let patched = std::fs::read_to_string(&txt).unwrap();
        assert!(patched.contains("Sample ID: renamed"), "got {patched}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_save_rejects_non_txt_sidecar_destination() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_save_pks_destination_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        let pks = dir.join("run.PKS");
        let txt = dir.join("run.txt");
        let original = "Capillary #: 1\nWell: D1\nSample ID: old\n";
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(&pks, b"sidecar").unwrap();
        std::fs::write(&txt, original).unwrap();

        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "D1: renamed")],
        };
        let loaded = crate::io::LoadedRun::new(
            run.clone(),
            crate::io::DetectedFormat {
                path: raw.clone(),
                identity: raw.clone(),
                format: crate::io::TraceFormat::FragmentAnalyzerRun {
                    entry: crate::io::FragmentAnalyzerEntry::Raw,
                },
                capabilities: crate::io::SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: false,
                    can_edit_markers: false,
                },
            },
            vec![],
            vec![],
        );

        assert!(!supports_save_path(&loaded, &pks));
        let err = save_run(&run, &raw, &pks).unwrap_err().to_string();

        assert!(
            err.contains("Fragment Analyzer runs can only be saved in place"),
            "got {err}"
        );
        assert_eq!(std::fs::read_to_string(&txt).unwrap(), original);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_raw_only_source_is_not_save_supported_and_errors_clearly() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_raw_only_save_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        std::fs::write(&raw, b"FA\0\0").unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "renamed")],
        };
        let loaded = crate::io::LoadedRun::new(
            run.clone(),
            crate::io::DetectedFormat {
                path: raw.clone(),
                identity: raw.clone(),
                format: crate::io::TraceFormat::FragmentAnalyzerRun {
                    entry: crate::io::FragmentAnalyzerEntry::Raw,
                },
                capabilities: crate::io::SourceCapabilities {
                    can_rename: false,
                    can_save_in_place: false,
                    can_save_as_xml: false,
                    can_edit_markers: false,
                },
            },
            vec![],
            vec![],
        );

        assert!(!supports_save_path(&loaded, &raw));
        let err = save_run(&run, &raw, &raw).unwrap_err().to_string();
        assert!(err.contains("FA .txt sidecar not found"), "got {err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_save_rejects_arbitrary_txt_destination_without_patching_sidecar() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_save_arbitrary_destination_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        let txt = dir.join("run.txt");
        let arbitrary = dir.join("notes.txt");
        let original = "Capillary #: 1\nWell: D1\nSample ID: old\n";
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(&txt, original).unwrap();
        std::fs::write(&arbitrary, "not an FA sidecar").unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "D1: renamed")],
        };
        let loaded = crate::io::LoadedRun::new(
            run.clone(),
            crate::io::DetectedFormat {
                path: raw.clone(),
                identity: raw.clone(),
                format: crate::io::TraceFormat::FragmentAnalyzerRun {
                    entry: crate::io::FragmentAnalyzerEntry::Raw,
                },
                capabilities: crate::io::SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: false,
                    can_edit_markers: false,
                },
            },
            vec![],
            vec![],
        );

        assert!(!supports_save_path(&loaded, &arbitrary));
        let err = save_run(&run, &raw, &arbitrary).unwrap_err().to_string();

        assert!(
            err.contains("Fragment Analyzer runs can only be saved in place"),
            "got {err}"
        );
        assert_eq!(std::fs::read_to_string(&txt).unwrap(), original);
        assert_eq!(
            std::fs::read_to_string(&arbitrary).unwrap(),
            "not an FA sidecar"
        );
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn fa_zip_without_txt_is_not_save_supported_and_errors_clearly() {
        use std::io::Write;

        let zip_path = std::env::temp_dir().join(format!(
            "traceio_fa_zip_no_txt_save_{}_{}.zip",
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
            zip.finish().unwrap();
        }
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "renamed")],
        };
        let loaded = crate::io::LoadedRun::new(
            run.clone(),
            crate::io::DetectedFormat {
                path: zip_path.clone(),
                identity: zip_path.clone(),
                format: crate::io::TraceFormat::FragmentAnalyzerRun {
                    entry: crate::io::FragmentAnalyzerEntry::Zip,
                },
                capabilities: crate::io::SourceCapabilities {
                    can_rename: false,
                    can_save_in_place: false,
                    can_save_as_xml: false,
                    can_edit_markers: false,
                },
            },
            vec![],
            vec![],
        );

        assert!(!supports_save_path(&loaded, &zip_path));
        let err = save_run(&run, &zip_path, &zip_path)
            .unwrap_err()
            .to_string();
        assert!(err.contains("has no .txt entry to rename"), "got {err}");
        std::fs::remove_file(zip_path).unwrap();
    }

    #[test]
    fn tapestation_csv_save_is_rejected_as_read_only() {
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "renamed")],
        };

        let err = save_run(
            &run,
            Path::new("run_Electropherogram.csv"),
            Path::new("run.xml"),
        )
        .unwrap_err()
        .to_string();

        assert!(
            err.contains("TapeStation exports are read-only"),
            "got {err}"
        );
    }

    #[test]
    fn xad_source_saves_as_xml_but_not_in_place_xad() {
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "renamed")],
        };
        let dir = unique_temp_dir("traceio_xad_save_as_xml");
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xad");
        let dst = dir.join("run.xml");
        std::fs::write(
            &src,
            synthetic_xad(
                "<Chipset><Chips><Chip><Files><File><Samples>\
                 <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
                 </Samples></File></Files></Chip></Chips></Chipset>",
            ),
        )
        .unwrap();

        save_run(&run, &src, &dst).unwrap();

        let saved = std::fs::read_to_string(&dst).unwrap();
        assert!(saved.contains("<Name>renamed</Name>"), "got {saved}");
        assert!(!saved.contains("<Name>old</Name>"));

        let err = save_run(&run, &src, &src).unwrap_err().to_string();

        assert!(err.contains("native .xad container"), "got {err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn xad_source_saves_as_xml_gz() {
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: vec![],
            regions: vec![],
            samples: vec![make_sample(1, "gz renamed")],
        };
        let dir = unique_temp_dir("traceio_xad_save_as_xml_gz");
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xad");
        let dst = dir.join("run.xml.gz");
        std::fs::write(
            &src,
            synthetic_xad(
                "<Chipset><Chips><Chip><Files><File><Samples>\
                 <Sample><Name>old</Name><WellNumber>1</WellNumber></Sample>\
                 </Samples></File></Files></Chip></Chips></Chipset>",
            ),
        )
        .unwrap();

        save_run(&run, &src, &dst).unwrap();

        let saved = read_gz_xml(&dst);
        assert!(saved.contains("<Name>gz renamed</Name>"), "got {saved}");
        assert!(!saved.contains("<Name>old</Name>"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn save_path_uses_loaded_source_template() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_save_path_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("run.xml");
        let dst = dir.join("renamed.xml");
        std::fs::write(
            &src,
            "<Chipset><Chips><Chip><Files><File><Samples>\
            <Sample><Name>old</Name><WellNumber>7</WellNumber></Sample>\
            </Samples></File></Files></Chip></Chips></Chipset>",
        )
        .unwrap();
        let loaded = crate::io::LoadedRun::new(
            Electrophoresis {
                assay: Default::default(),
                ladder_peaks: vec![],
                regions: vec![],
                samples: vec![make_sample(7, "new")],
            },
            crate::io::DetectedFormat {
                path: src.clone(),
                identity: src,
                format: crate::io::TraceFormat::BioanalyzerXml,
                capabilities: crate::io::SourceCapabilities {
                    can_rename: true,
                    can_save_in_place: true,
                    can_save_as_xml: true,
                    can_edit_markers: true,
                },
            },
            vec![],
            vec![],
        );

        save_path(&loaded, &dst).unwrap();

        let patched = std::fs::read_to_string(&dst).unwrap();
        assert!(patched.contains("<Name>new</Name>"), "got: {patched}");
        assert!(!patched.contains("<Name>old</Name>"));
        std::fs::remove_dir_all(dir).unwrap();
    }

    fn unique_temp_dir(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "{tag}_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    fn write_gz_xml(path: &Path, xml: &str) {
        let file = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        enc.write_all(xml.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    fn read_gz_xml(path: &Path) -> String {
        let raw = std::fs::read(path).unwrap();
        let mut dec = flate2::read::GzDecoder::new(&raw[..]);
        let mut xml = String::new();
        dec.read_to_string(&mut xml).unwrap();
        xml
    }

    fn synthetic_xad(xml: &str) -> Vec<u8> {
        use base64::Engine;

        let mut utf16 = Vec::with_capacity(xml.len() * 2);
        for unit in xml.encode_utf16() {
            utf16.extend_from_slice(&unit.to_le_bytes());
        }

        let mut deflater =
            flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
        deflater.write_all(&utf16).unwrap();
        let deflated = deflater.finish().unwrap();

        let mut framed = Vec::with_capacity(1 + deflated.len() + 9);
        framed.push(0);
        framed.extend_from_slice(&deflated);
        framed.extend_from_slice(&[0; 9]);

        format!(
            "<root><compressed_data dt:dt=\"bin.base64\">{}</compressed_data></root>",
            base64::engine::general_purpose::STANDARD.encode(framed)
        )
        .into_bytes()
    }

    fn assert_no_bioanalyzer_temp_files(dir: &Path) {
        let leftovers: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .filter(|name| name.contains(".tmp-"))
            .collect();
        assert!(leftovers.is_empty(), "leftover temp files: {leftovers:?}");
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
            regions: vec![],
        }
    }
}
