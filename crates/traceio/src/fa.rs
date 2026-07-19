//! Native Agilent / Advanced Analytical (AATI) **Fragment Analyzer** run reader.
//!
//! Reverse-engineered from a real run directory; see `docs/fa_format.md` for the
//! full byte-level notes. A run is a directory of sibling files sharing one
//! timestamp stem, e.g. `2025 11 19 16H 03M.{raw,PKS,txt,...}`. The pieces this
//! reader uses:
//!
//! * **`.raw`** — the CCD acquisition. `FA\0\0` magic, all values **big-endian
//!   `u16`**. A short header carries the CCD line width (pixels) and a table of
//!   the capillary centre columns; the payload at [`DATA_START`] is
//!   `scans × width` intensities. One capillary's electropherogram is its centre
//!   pixel column read across every scan (scan index ≈ migration time in s).
//! * **`.PKS`** — LabVIEW-flattened peak/analysis data. Used here only for the
//!   size-calibration anchor **times** (a strictly-increasing `u16` run), paired
//!   with the standard size ladder to map scan → base pairs.
//! * **`.txt`** — human-readable capillary → well → sample-name list.
//!
//! Fragment Analyzer runs arrive already size-calibrated, so unlike the
//! Bioanalyzer path this reader fills each sample's `length` directly by
//! interpolating the ladder and does **not** use [`crate::calibration`].

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

use crate::model::{AssayInfo, Electrophoresis, LadderPeak, Peak, Sample};

/// Byte offset where the `.raw` CCD payload begins (after the header).
const DATA_START: usize = 0x7d0;
/// `.raw` magic bytes.
const RAW_MAGIC: &[u8] = b"FA\0\0";
/// Half-width (in pixels) of the window averaged around each capillary centre.
const CAP_WINDOW: usize = 4;

/// Standard Fragment Analyzer dsDNA size ladder (bp), paired 1:1 with the
/// calibration anchor times read from `.PKS`. Matches the "1–6000 bp" reagent
/// kit ladder used by the reference run. If a run's anchor-time count differs
/// from this, calibration is skipped (traces stay on the time axis).
const LADDER_BP: [f64; 16] = [
    1.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1200.0, 1500.0,
    2000.0, 3000.0, 6000.0,
];

/// Mass concentration setpoints for the standard 1-6000 bp FA dsDNA ladder.
/// These are used when `.PKS` does not expose trustworthy quantity columns:
/// concentration is computed from the ladder peak areas and marker scaling.
const LADDER_CONC: [f64; 16] = [
    0.1824, 1.6124, 1.6549, 1.5459, 1.5098, 4.0840, 1.2509, 1.2085, 1.2342, 1.1831, 4.0657, 1.3504,
    1.2080, 1.1280, 1.0403, 0.0309,
];

/// One capillary's identity from the `.txt` sidecar.
struct CapInfo {
    well: String,
    sample_id: String,
    description: String,
}

/// Fragment Analyzer run metadata that does not fit the normalized
/// [`Electrophoresis`] model.
#[derive(Debug, Clone, Default)]
pub struct FaMetadata {
    /// Header fields from the `.txt` sidecar before the capillary list.
    pub run_header: BTreeMap<String, String>,
    /// Parsed `method.mthd` sections.
    pub method: IniDocument,
    /// Parsed `.ANAI` analysis settings.
    pub analysis: IniDocument,
    /// Current/voltage/pressure samples from `.current`, one row per scan when
    /// present.
    pub current: Vec<FaCurrentPoint>,
    /// Raw per-frame timing/status lines from `Timing.txt`.
    pub timing: Vec<String>,
    /// Raw exposure/status lines from `ExpTime.txt`.
    pub exposure: Vec<String>,
    /// Bytes from `CameraImage.bmp`; observed files are PNG despite the name.
    pub camera_image: Option<Vec<u8>>,
}

/// INI-like document, preserving section and key order within sorted maps.
pub type IniDocument = BTreeMap<String, BTreeMap<String, String>>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FaCurrentPoint {
    pub current_ua: f64,
    pub voltage_kv: f64,
    pub pressure_psi: f64,
}

/// Read a Fragment Analyzer run. `path` may be the `.raw` file itself or the run
/// directory (in which case the single `.raw` inside is used).
pub fn read_fa_run(path: &Path) -> Result<Electrophoresis> {
    if has_zip_extension(path) {
        return read_fa_zip(path);
    }
    // Filesystem run: a `.raw` file (or a directory holding one) plus siblings.
    let raw_path = resolve_raw_path(path)?;
    let dir = raw_path.parent().unwrap_or_else(|| Path::new("."));

    let raw =
        std::fs::read(&raw_path).with_context(|| format!("reading {}", raw_path.display()))?;
    let pks = sibling_path(dir, &raw_path, "PKS")?
        .map(|p| std::fs::read(&p).with_context(|| format!("reading {}", p.display())))
        .transpose()?;
    let txt = sibling_path(dir, &raw_path, "txt")?
        .map(|p| std::fs::read(&p).with_context(|| format!("reading {}", p.display())))
        .transpose()?
        .map(|b| String::from_utf8_lossy(&b).into_owned());
    let file_name = raw_path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    build_fa_run(&raw, pks.as_deref(), txt.as_deref(), file_name)
}

/// Read FA-specific method, analysis, and QC sidecars. Missing optional files
/// are left empty; malformed present files return an error with their file name.
pub fn read_fa_metadata(path: &Path) -> Result<FaMetadata> {
    if has_zip_extension(path) {
        return read_fa_metadata_zip(path);
    }

    let raw_path = resolve_raw_path(path)?;
    let dir = raw_path.parent().unwrap_or_else(|| Path::new("."));
    let text_sibling = |ext: &str| -> Result<Option<String>> {
        sibling_path(dir, &raw_path, ext)?
            .map(|p| read_lossy_text(&p))
            .transpose()
    };
    let txt = text_sibling("txt")?;
    let method = dir.join("method.mthd");
    let anai = text_sibling("ANAI")?;
    let current = text_sibling("current")?;
    let timing = dir.join("Timing.txt");
    let exposure = dir.join("ExpTime.txt");
    let camera = dir.join("CameraImage.bmp");

    Ok(FaMetadata {
        run_header: txt.as_deref().map(parse_txt_header).unwrap_or_default(),
        method: if method.exists() {
            parse_ini(&read_lossy_text(&method)?)
        } else {
            IniDocument::default()
        },
        analysis: anai.as_deref().map(parse_ini).unwrap_or_default(),
        current: current
            .as_deref()
            .map(parse_current_log)
            .transpose()?
            .unwrap_or_default(),
        timing: if timing.exists() {
            read_lines_lossy(&timing)?
        } else {
            Vec::new()
        },
        exposure: if exposure.exists() {
            read_lines_lossy(&exposure)?
        } else {
            Vec::new()
        },
        camera_image: if camera.exists() {
            Some(std::fs::read(&camera).with_context(|| format!("reading {}", camera.display()))?)
        } else {
            None
        },
    })
}

/// Build an [`Electrophoresis`] from the raw bytes of the three FA files this
/// reader uses, independent of whether they came from disk or a zip entry.
fn build_fa_run(
    raw: &[u8],
    pks: Option<&[u8]>,
    txt: Option<&str>,
    file_name: String,
) -> Result<Electrophoresis> {
    let caps = txt.map(parse_txt).unwrap_or_default();
    let traces = parse_raw(raw, caps.len())?;

    // Size calibration and peaks both come from `.PKS`.
    let anchor_times = pks.and_then(|d| pks_anchor_times(d, traces.scans));
    let calib = anchor_times.filter(|t| t.len() == LADDER_BP.len());
    let peaks = pks
        .and_then(|d| pks_peaks(d, traces.scans, calib.as_deref()))
        .unwrap_or_default();

    let mut run = build_run(traces, &caps, calib.as_deref(), &peaks, file_name)?;
    compute_quantities(&mut run);
    Ok(run)
}

/// Read a Fragment Analyzer run packaged as a single `.zip` — the whole run
/// folder (or just its files) compressed into one archive. Locates the `.raw`
/// entry by its `FA\0\0` magic and reads the sibling `.PKS`/`.txt` entries that
/// share its stem (directory prefix included, so nested zips work).
fn read_fa_zip(path: &Path) -> Result<Electrophoresis> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut zip = zip::ZipArchive::new(std::io::BufReader::new(file))
        .with_context(|| format!("reading zip {}", path.display()))?;

    // Entry names, and the `.raw` entry identified by its magic bytes.
    let mut names = Vec::with_capacity(zip.len());
    let mut raw_name: Option<String> = None;
    for i in 0..zip.len() {
        let mut e = zip.by_index(i)?;
        let name = e.name().to_string();
        if raw_name.is_none() {
            if let Some(name) = fa_raw_entry_name(&mut e) {
                raw_name = Some(name.clone());
            }
        }
        names.push(name);
    }
    let raw_name = raw_name.ok_or_else(|| {
        anyhow!(
            "{} is not a Fragment Analyzer run (no .raw entry with FA magic)",
            path.display()
        )
    })?;

    let stem = strip_known_ext(&raw_name);
    let pks_name = find_sibling(&names, stem, "PKS")?;
    let txt_name = find_sibling(&names, stem, "txt")?;

    let raw = read_zip_entry(&mut zip, &raw_name)?;
    let pks = pks_name.map(|n| read_zip_entry(&mut zip, &n)).transpose()?;
    let txt = txt_name
        .map(|n| read_zip_entry(&mut zip, &n))
        .transpose()?
        .map(|b| String::from_utf8_lossy(&b).into_owned());
    let file_name = basename(&raw_name).to_string();
    build_fa_run(&raw, pks.as_deref(), txt.as_deref(), file_name)
}

fn read_fa_metadata_zip(path: &Path) -> Result<FaMetadata> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut zip = zip::ZipArchive::new(std::io::BufReader::new(file))
        .with_context(|| format!("reading zip {}", path.display()))?;

    let mut names = Vec::with_capacity(zip.len());
    let mut raw_name: Option<String> = None;
    for i in 0..zip.len() {
        let mut e = zip.by_index(i)?;
        let name = e.name().to_string();
        if raw_name.is_none() {
            if let Some(name) = fa_raw_entry_name(&mut e) {
                raw_name = Some(name.clone());
            }
        }
        names.push(name);
    }
    let raw_name = raw_name.ok_or_else(|| {
        anyhow!(
            "{} is not a Fragment Analyzer run (no .raw entry with FA magic)",
            path.display()
        )
    })?;
    let stem = strip_known_ext(&raw_name);
    let prefix = raw_name
        .rsplit_once('/')
        .map(|(prefix, _)| format!("{prefix}/"))
        .unwrap_or_default();

    let txt = find_sibling(&names, stem, "txt")?
        .map(|n| read_zip_entry(&mut zip, &n))
        .transpose()?
        .map(lossy_string);
    let anai = find_sibling(&names, stem, "ANAI")?
        .map(|n| read_zip_entry(&mut zip, &n))
        .transpose()?
        .map(lossy_string);
    let current = find_sibling(&names, stem, "current")?
        .map(|n| read_zip_entry(&mut zip, &n))
        .transpose()?
        .map(lossy_string);
    let method_name = find_named_entry(&names, &prefix, "method.mthd");
    let timing_name = find_named_entry(&names, &prefix, "Timing.txt");
    let exposure_name = find_named_entry(&names, &prefix, "ExpTime.txt");
    let camera_name = find_named_entry(&names, &prefix, "CameraImage.bmp");

    Ok(FaMetadata {
        run_header: txt.as_deref().map(parse_txt_header).unwrap_or_default(),
        method: method_name
            .map(|n| read_zip_entry(&mut zip, &n).map(lossy_string))
            .transpose()?
            .as_deref()
            .map(parse_ini)
            .unwrap_or_default(),
        analysis: anai.as_deref().map(parse_ini).unwrap_or_default(),
        current: current
            .as_deref()
            .map(parse_current_log)
            .transpose()?
            .unwrap_or_default(),
        timing: timing_name
            .map(|n| read_zip_entry(&mut zip, &n).map(lossy_string))
            .transpose()?
            .as_deref()
            .map(lines)
            .unwrap_or_default(),
        exposure: exposure_name
            .map(|n| read_zip_entry(&mut zip, &n).map(lossy_string))
            .transpose()?
            .as_deref()
            .map(lines)
            .unwrap_or_default(),
        camera_image: camera_name
            .map(|n| read_zip_entry(&mut zip, &n))
            .transpose()?,
    })
}

fn read_zip_entry<R: std::io::Read + std::io::Seek>(
    zip: &mut zip::ZipArchive<R>,
    name: &str,
) -> Result<Vec<u8>> {
    use std::io::Read;
    let mut e = zip
        .by_name(name)
        .with_context(|| format!("zip entry {name}"))?;
    let mut v = Vec::with_capacity(e.size() as usize);
    e.read_to_end(&mut v)?;
    Ok(v)
}

fn fa_raw_entry_name(entry: &mut zip::read::ZipFile<'_>) -> Option<String> {
    use std::io::Read;
    if !entry.is_file() || !name_has_ext(entry.name(), "raw") {
        return None;
    }
    let name = entry.name().to_string();
    let mut magic = [0u8; 4];
    (entry.read_exact(&mut magic).is_ok() && magic == RAW_MAGIC).then_some(name)
}

fn find_named_entry(names: &[String], prefix: &str, name: &str) -> Option<String> {
    let exact = format!("{prefix}{name}");
    names
        .iter()
        .find(|n| n.eq_ignore_ascii_case(&exact))
        .cloned()
        .or_else(|| {
            names
                .iter()
                .find(|n| basename(n).eq_ignore_ascii_case(name))
                .cloned()
        })
}

fn read_lossy_text(path: &Path) -> Result<String> {
    std::fs::read(path)
        .map(lossy_string)
        .with_context(|| format!("reading {}", path.display()))
}

fn read_lines_lossy(path: &Path) -> Result<Vec<String>> {
    read_lossy_text(path).map(|s| lines(&s))
}

fn lossy_string(bytes: Vec<u8>) -> String {
    String::from_utf8_lossy(&bytes).into_owned()
}

fn lines(text: &str) -> Vec<String> {
    text.lines()
        .map(|l| l.trim_end_matches('\r').to_string())
        .collect()
}

fn parse_txt_header(text: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.starts_with("Capillary #:") || line.starts_with("====") {
            break;
        }
        if let Some((key, value)) = line.split_once(':') {
            out.insert(key.trim().to_string(), value.trim().to_string());
        }
    }
    out
}

fn parse_ini(text: &str) -> IniDocument {
    let mut out = IniDocument::new();
    let mut section = String::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with(';') || line.starts_with('#') {
            continue;
        }
        if let Some(name) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
            section = name.trim().to_string();
            out.entry(section.clone()).or_default();
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        out.entry(section.clone()).or_default().insert(
            key.trim().to_string(),
            value.trim().trim_matches('"').to_string(),
        );
    }
    out
}

fn parse_current_log(text: &str) -> Result<Vec<FaCurrentPoint>> {
    let mut rows = Vec::new();
    for (idx, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if idx == 0 && line.to_ascii_lowercase().contains("current") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 3 {
            bail!("FA .current line {} has {} columns", idx + 1, cols.len());
        }
        rows.push(FaCurrentPoint {
            current_ua: parse_decimal(cols[0])?,
            voltage_kv: parse_decimal(cols[1])?,
            pressure_psi: parse_decimal(cols[2])?,
        });
    }
    Ok(rows)
}

fn parse_decimal(s: &str) -> Result<f64> {
    s.trim()
        .replace(',', ".")
        .parse::<f64>()
        .with_context(|| format!("parsing numeric value {s:?}"))
}

/// Save user-renamed FA sample IDs back to the `.txt` sidecar. The `.raw`
/// payload is immutable; only `Sample ID:` values are rewritten, preserving the
/// rest of the human-readable sidecar. For a zipped run the archive is rewritten
/// in place with only its `.txt` entry patched.
pub fn save_txt_names(path: &Path, run: &Electrophoresis) -> Result<()> {
    if has_zip_extension(path) {
        return save_txt_names_zip(path, run);
    }
    let raw_path = resolve_raw_path(path)?;
    let dir = raw_path.parent().unwrap_or_else(|| Path::new("."));
    let txt_path = sibling_path(dir, &raw_path, "txt")?
        .ok_or_else(|| anyhow!("FA .txt sidecar not found next to {}", raw_path.display()))?;
    let text = std::fs::read_to_string(&txt_path)
        .or_else(|_| std::fs::read(&txt_path).map(|b| String::from_utf8_lossy(&b).into_owned()))
        .with_context(|| format!("reading {}", txt_path.display()))?;
    let (patched, changed) = patch_txt_names(&text, run);
    if changed == 0 {
        bail!(
            "FA .txt sidecar {} did not contain any Sample ID entries",
            txt_path.display()
        );
    }
    write_file_atomic(&txt_path, patched.as_bytes())?;
    Ok(())
}

/// Rewrite a zipped FA run in place, patching only the `.txt` entry's sample
/// IDs. Other entries are read and written into a fresh archive; the new archive
/// is written to a unique temp file and swapped in.
fn save_txt_names_zip(path: &Path, run: &Electrophoresis) -> Result<()> {
    use std::io::{Read, Write};

    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut zip = zip::ZipArchive::new(std::io::BufReader::new(file))
        .with_context(|| format!("reading zip {}", path.display()))?;

    // Locate the `.txt` entry sharing the `.raw` entry's stem.
    let names: Vec<String> = (0..zip.len())
        .map(|i| zip.by_index(i).map(|e| e.name().to_string()))
        .collect::<std::result::Result<_, _>>()?;
    let raw_name = {
        let mut found = None;
        for i in 0..zip.len() {
            let mut entry = zip.by_index(i)?;
            if let Some(name) = fa_raw_entry_name(&mut entry) {
                found = Some(name);
                break;
            }
        }
        found.ok_or_else(|| anyhow!("zip {} has no FA .raw entry", path.display()))?
    };
    let txt_name = find_sibling(&names, strip_known_ext(&raw_name), "txt")?
        .ok_or_else(|| anyhow!("zip {} has no .txt entry to rename", path.display()))?;

    let mut entries = Vec::with_capacity(zip.len());
    let mut old_text = None;
    for i in 0..zip.len() {
        let mut entry = zip.by_index(i)?;
        let name = entry.name().to_string();
        let options = zip_entry_options(&entry);
        let is_dir = entry.is_dir();
        let mut data = Vec::with_capacity(entry.size() as usize);
        if !is_dir {
            entry
                .read_to_end(&mut data)
                .with_context(|| format!("reading zip entry {name}"))?;
        }
        if name == txt_name {
            old_text = Some(String::from_utf8_lossy(&data).into_owned());
        }
        entries.push(ZipEntryData {
            name,
            data,
            options,
            is_dir,
        });
    }
    let old_text = old_text.ok_or_else(|| anyhow!("zip entry {txt_name} disappeared"))?;
    let (patched, changed) = patch_txt_names(&old_text, run);
    if changed == 0 {
        bail!("FA .txt entry {txt_name} did not contain any Sample ID entries");
    }

    // Write a fresh archive next to the target, then replace it.
    let (tmp, out) = create_unique_temp_file(path)?;
    {
        let mut writer = zip::ZipWriter::new(std::io::BufWriter::new(out));
        for entry in entries {
            if entry.is_dir {
                writer.add_directory(&entry.name, entry.options)?;
            } else if entry.name == txt_name {
                writer.start_file(&entry.name, entry.options)?;
                writer.write_all(patched.as_bytes())?;
            } else {
                writer.start_file(&entry.name, entry.options)?;
                writer.write_all(&entry.data)?;
            }
        }
        writer.finish()?;
    }
    drop(zip); // close the original before replacing it (matters on Windows)
    replace_file(&tmp, path)
        .with_context(|| format!("replacing {} with updated archive", path.display()))?;
    Ok(())
}

/// True if `path` looks like a Fragment Analyzer run entry point: a single
/// `.zip` packaging the run (the recommended form), a `.raw` file with the
/// `FA\0\0` magic, a directory containing one, or any other sibling file inside
/// such a directory (so dropping a run's `.PKS`/`.txt`/etc. opens the run).
pub fn is_fa_path(path: &Path) -> bool {
    if path.is_dir() {
        return matches!(find_raw_in_dir(path), Ok(Some(_)));
    }
    if has_zip_extension(path) {
        return zip_has_fa_raw(path);
    }
    if has_raw_extension(path) {
        let mut buf = [0u8; 4];
        return matches!(read_prefix(path, &mut buf), Ok(())) && buf == RAW_MAGIC;
    }
    // Forgiving: a plain member of an FA run folder. Never steal a Bioanalyzer
    // file — it has its own reader and could legitimately sit in the same dir.
    if is_bioanalyzer_ext(path) {
        return false;
    }
    parent_dir(path).is_some_and(|p| matches!(find_raw_in_dir(p), Ok(Some(_))))
}

/// Canonical identity of an FA run for any path `is_fa_path` accepts: the `.zip`
/// itself for a zipped run, otherwise the run's `.raw` file. Lets the app tell
/// when two entry points (the folder, the `.raw`, or a sibling like `.PKS` /
/// `.txt`) refer to the same run, so a multi-file drop opens it only once.
pub fn run_identity(path: &Path) -> PathBuf {
    if has_zip_extension(path) {
        return path.to_path_buf();
    }
    resolve_raw_path(path).unwrap_or_else(|_| path.to_path_buf())
}

/// The path's parent directory, or `None` for a bare filename (empty parent).
fn parent_dir(path: &Path) -> Option<&Path> {
    path.parent().filter(|p| !p.as_os_str().is_empty())
}

struct ZipEntryData {
    name: String,
    data: Vec<u8>,
    options: zip::write::SimpleFileOptions,
    is_dir: bool,
}

fn zip_entry_options(entry: &zip::read::ZipFile<'_>) -> zip::write::SimpleFileOptions {
    let mut options =
        zip::write::SimpleFileOptions::default().compression_method(entry.compression());
    if let Some(time) = entry.last_modified() {
        options = options.last_modified_time(time);
    }
    if let Some(mode) = entry.unix_mode() {
        options = options.unix_permissions(mode);
    }
    options
}

fn sibling_path(dir: &Path, raw_path: &Path, ext: &str) -> Result<Option<PathBuf>> {
    let raw_stem = raw_path.file_stem().ok_or_else(|| {
        anyhow!(
            "FA run: cannot derive file stem from {}",
            raw_path.display()
        )
    })?;
    let exact_name = PathBuf::from(format!("{}.{}", raw_stem.to_string_lossy(), ext));
    let exact = dir.join(&exact_name);
    if exact.exists() {
        return Ok(Some(exact));
    }

    let raw_stem = raw_stem.to_string_lossy();
    let mut matches = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let path = entry?.path();
        let stem_matches = path
            .file_stem()
            .map(|s| s.to_string_lossy().eq_ignore_ascii_case(&raw_stem))
            .unwrap_or(false);
        let ext_matches = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case(ext));
        if stem_matches && ext_matches {
            matches.push(path);
        }
    }
    match matches.len() {
        0 => Ok(None),
        1 => Ok(matches.pop()),
        _ => Err(anyhow!(
            "ambiguous FA .{ext} sidecars next to {}: {}",
            raw_path.display(),
            matches
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

fn create_unique_temp_file(dst: &Path) -> Result<(PathBuf, std::fs::File)> {
    let dir = dst.parent().unwrap_or_else(|| Path::new("."));
    let name = dst
        .file_name()
        .map(|n| n.to_string_lossy())
        .ok_or_else(|| anyhow!("cannot derive temp name from {}", dst.display()))?;
    for attempt in 0..100 {
        let candidate = dir.join(format!(
            ".{name}.tmp-{}-{}-{attempt}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_nanos()
        ));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => return Ok((candidate, file)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => {
                return Err(e).with_context(|| format!("creating {}", candidate.display()));
            }
        }
    }
    Err(anyhow!(
        "could not allocate unique temp path next to {}",
        dst.display()
    ))
}

fn write_file_atomic(dst: &Path, contents: &[u8]) -> Result<()> {
    use std::io::Write;

    let (tmp, mut out) = create_unique_temp_file(dst)?;
    out.write_all(contents)
        .with_context(|| format!("writing {}", tmp.display()))?;
    out.sync_all()
        .with_context(|| format!("syncing {}", tmp.display()))?;
    drop(out);
    replace_file(&tmp, dst)?;
    Ok(())
}

fn replace_file(tmp: &Path, dst: &Path) -> Result<()> {
    match std::fs::rename(tmp, dst) {
        Ok(()) => Ok(()),
        Err(first) if dst.exists() => {
            let backup = dst.with_extension(format!(
                "replace-bak-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos()
            ));
            std::fs::rename(dst, &backup).with_context(|| {
                format!(
                    "moving {} aside after replace failed with {first}",
                    dst.display()
                )
            })?;
            match std::fs::rename(tmp, dst) {
                Ok(()) => {
                    let _ = std::fs::remove_file(&backup);
                    Ok(())
                }
                Err(e) => {
                    let _ = std::fs::rename(&backup, dst);
                    Err(e).with_context(|| {
                        format!(
                            "renaming {} to {} after moving original aside",
                            tmp.display(),
                            dst.display()
                        )
                    })
                }
            }
        }
        Err(e) => {
            Err(e).with_context(|| format!("renaming {} to {}", tmp.display(), dst.display()))
        }
    }
}

/// True for the Bioanalyzer extensions handled by another reader.
fn is_bioanalyzer_ext(path: &Path) -> bool {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    name.ends_with(".xad") || name.ends_with(".xml") || name.ends_with(".gz")
}

/// True if `path` is a zip archive containing a `.raw` entry with the FA magic.
fn zip_has_fa_raw(path: &Path) -> bool {
    use std::io::Read;
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let Ok(mut zip) = zip::ZipArchive::new(std::io::BufReader::new(file)) else {
        return false;
    };
    for i in 0..zip.len() {
        let Ok(mut e) = zip.by_index(i) else { continue };
        if e.is_file() && name_has_ext(e.name(), "raw") {
            let mut magic = [0u8; 4];
            if e.read_exact(&mut magic).is_ok() && magic == RAW_MAGIC {
                return true;
            }
        }
    }
    false
}

fn has_zip_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("zip"))
}

/// Case-insensitive extension test for a zip-entry name (no `Path` needed).
fn name_has_ext(name: &str, ext: &str) -> bool {
    std::path::Path::new(name)
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case(ext))
}

/// Strip a trailing `.<ext>` (any extension) from a zip-entry name, keeping any
/// directory prefix, so `dir/run.raw` → `dir/run`.
fn strip_known_ext(name: &str) -> &str {
    match name.rfind('.') {
        Some(dot) if !name[dot + 1..].contains('/') => &name[..dot],
        _ => name,
    }
}

/// Find the entry sharing `stem` with the given extension. Prefer the exact
/// sidecar path, then accept one unambiguous case-insensitive fallback.
fn find_sibling(names: &[String], stem: &str, ext: &str) -> Result<Option<String>> {
    let exact = format!("{stem}.{ext}");
    if names.iter().any(|n| n == &exact) {
        return Ok(Some(exact));
    }

    let mut matches = names
        .iter()
        .filter(|n| strip_known_ext(n).eq_ignore_ascii_case(stem) && name_has_ext(n, ext))
        .cloned()
        .collect::<Vec<_>>();
    match matches.len() {
        0 => Ok(None),
        1 => Ok(matches.pop()),
        _ => Err(anyhow!(
            "ambiguous FA .{ext} sidecars in zip for stem {stem}: {}",
            matches.join(", ")
        )),
    }
}

fn basename(name: &str) -> &str {
    name.rsplit('/').next().unwrap_or(name)
}

fn resolve_raw_path(path: &Path) -> Result<PathBuf> {
    if path.is_dir() {
        return find_raw_in_dir(path)?
            .ok_or_else(|| anyhow!("no .raw file found in FA run dir {}", path.display()));
    }
    if has_raw_extension(path) {
        return Ok(path.to_path_buf());
    }
    // A sibling member of the run (e.g. `.PKS`/`.txt`): use its folder's `.raw`.
    if let Some(parent) = parent_dir(path) {
        if let Some(raw) = find_raw_in_dir(parent)? {
            return Ok(raw);
        }
    }
    Ok(path.to_path_buf())
}

fn find_raw_in_dir(dir: &Path) -> Result<Option<PathBuf>> {
    let mut matches = Vec::new();
    for entry in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let path = entry?.path();
        if !has_raw_extension(&path) {
            continue;
        }
        let mut buf = [0u8; 4];
        if read_prefix(&path, &mut buf).is_ok() && buf == RAW_MAGIC {
            matches.push(path);
        }
    }

    match matches.len() {
        0 => Ok(None),
        1 => Ok(matches.pop()),
        _ => Err(anyhow!(
            "ambiguous FA .raw files in {}: {}",
            dir.display(),
            matches
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

fn has_raw_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("raw"))
}

fn read_prefix(path: &Path, buf: &mut [u8]) -> Result<()> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    f.read_exact(buf)?;
    Ok(())
}

/// Decoded `.raw` payload: the CCD grid plus the capillary centre columns.
struct RawTraces {
    scans: usize,
    width: usize,
    creation_date: String,
    /// `scans * width` big-endian u16 intensities.
    data: Vec<u16>,
    /// Centre pixel column for each capillary (in acquisition order).
    columns: Vec<usize>,
}

impl RawTraces {
    /// Averaged intensity window around a capillary centre for one scan.
    fn value(&self, cap: usize, scan: usize) -> f32 {
        let c = self.columns[cap];
        let lo = c.saturating_sub(CAP_WINDOW);
        let hi = (c + CAP_WINDOW).min(self.width - 1);
        let base = scan * self.width;
        let mut sum = 0u32;
        for p in lo..=hi {
            sum += self.data[base + p] as u32;
        }
        sum as f32 / (hi - lo + 1) as f32
    }
}

fn parse_raw(raw: &[u8], nwells_hint: usize) -> Result<RawTraces> {
    if raw.len() < DATA_START || &raw[..4] != RAW_MAGIC {
        bail!("not a Fragment Analyzer .raw file (bad magic)");
    }
    // CCD line width (pixels): big-endian u16 in the header.
    let width = be_u16(raw, 0xff) as usize;
    if !(16..=8192).contains(&width) {
        bail!("FA .raw: implausible CCD width {width}");
    }
    let body = raw.len() - DATA_START;
    if !body.is_multiple_of(2 * width) {
        bail!("FA .raw: payload {body} not a multiple of line width {width}");
    }
    let scans = body / (2 * width);
    if scans < 2 {
        bail!("FA .raw: too few scans ({scans})");
    }

    let mut data = Vec::with_capacity(scans * width);
    let mut o = DATA_START;
    for _ in 0..scans * width {
        data.push(u16::from_be_bytes([raw[o], raw[o + 1]]));
        o += 2;
    }

    let columns = capillary_columns(raw, width, nwells_hint)?;
    let creation_date = raw_header_datetime(raw);
    Ok(RawTraces {
        scans,
        width,
        creation_date,
        data,
        columns,
    })
}

/// Locate the capillary centre columns in the `.raw` header: a run of
/// strictly-increasing big-endian `u16` values, all in `1..width`, flanked by
/// zeros. Prefers a run whose length matches `nwells_hint` when that is > 0.
fn capillary_columns(raw: &[u8], width: usize, nwells_hint: usize) -> Result<Vec<usize>> {
    let end = DATA_START.min(raw.len());
    let mut best: Option<Vec<usize>> = None;
    let mut i = 4;
    while i + 4 <= end {
        // Require a zero separator before the run.
        if be_u16(raw, i) != 0 {
            i += 2;
            continue;
        }
        let mut j = i + 2;
        let mut run = Vec::new();
        let mut prev = 0u16;
        while j + 2 <= end {
            let v = be_u16(raw, j);
            if v == 0 || v as usize >= width || v <= prev {
                break;
            }
            run.push(v as usize);
            prev = v;
            j += 2;
        }
        // Followed by a zero, and a plausible length.
        let terminated = j + 2 <= end && be_u16(raw, j) == 0;
        if terminated && run.len() >= 2 {
            let matches_hint = nwells_hint > 0 && run.len() == nwells_hint;
            if matches_hint {
                return Ok(run);
            }
            if best.as_ref().is_none_or(|b| run.len() > b.len()) {
                best = Some(run);
            }
        }
        i = j.max(i + 2);
    }
    best.ok_or_else(|| anyhow!("FA .raw: could not locate the capillary column table"))
}

fn raw_header_datetime(raw: &[u8]) -> String {
    let date = ascii_cstr_field(raw, 0x64, 0x30);
    let time = ascii_cstr_field(raw, 0x96, 0x20);
    match (date.as_deref(), time.as_deref()) {
        (Some(date), Some(time)) => format!("{date} {time}"),
        (Some(date), None) => date.to_string(),
        (None, Some(time)) => time.to_string(),
        (None, None) => String::new(),
    }
}

fn ascii_cstr_field(raw: &[u8], start: usize, len: usize) -> Option<String> {
    let bytes = raw.get(start..start + len)?;
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let s = std::str::from_utf8(&bytes[..end]).ok()?.trim();
    (!s.is_empty() && s.chars().all(|c| c.is_ascii() && !c.is_control())).then(|| s.to_string())
}

/// Parse the `.txt` sidecar text into per-capillary metadata, in order.
fn parse_txt(text: &str) -> Vec<CapInfo> {
    let mut caps = Vec::new();
    let (mut well, mut sample, mut description) = (String::new(), String::new(), String::new());
    let mut in_cap = false;
    let flush = |caps: &mut Vec<CapInfo>,
                 w: &mut String,
                 s: &mut String,
                 d: &mut String,
                 in_cap: &mut bool| {
        if *in_cap {
            caps.push(CapInfo {
                well: std::mem::take(w),
                sample_id: std::mem::take(s),
                description: std::mem::take(d),
            });
        }
        *in_cap = false;
    };
    for line in text.lines() {
        let line = line.trim();
        if let Some(_rest) = line.strip_prefix("Capillary #:") {
            flush(
                &mut caps,
                &mut well,
                &mut sample,
                &mut description,
                &mut in_cap,
            );
            in_cap = true;
        } else if let Some(rest) = line.strip_prefix("Well:") {
            well = rest.trim().to_string();
        } else if let Some(rest) = line.strip_prefix("Sample ID:") {
            sample = rest.trim().to_string();
        } else if let Some(rest) = line.strip_prefix("Description:") {
            description = rest.trim().to_string();
        }
    }
    flush(
        &mut caps,
        &mut well,
        &mut sample,
        &mut description,
        &mut in_cap,
    );
    caps
}

fn patch_txt_names(text: &str, run: &Electrophoresis) -> (String, usize) {
    let mut out = String::with_capacity(text.len());
    let mut cap_idx: Option<usize> = None;
    let mut current_well = String::new();
    let mut changed = 0usize;

    for line in text.split_inclusive('\n') {
        let (body, ending) = line_body_and_ending(line);
        let trimmed = body.trim_start();
        if trimmed.starts_with("Capillary #:") {
            cap_idx = Some(cap_idx.map_or(0, |i| i + 1));
            current_well.clear();
            out.push_str(line);
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("Well:") {
            current_well = rest.trim().to_string();
            out.push_str(line);
            continue;
        }
        if trimmed.starts_with("Sample ID:") {
            if let Some(i) = cap_idx {
                if let Some(sample) = run.samples.get(i) {
                    let prefix_len = body.find("Sample ID:").unwrap_or(0);
                    let id = sidecar_sample_id(&sample.name, &current_well);
                    out.push_str(&body[..prefix_len]);
                    out.push_str("Sample ID: ");
                    out.push_str(&id);
                    out.push_str(ending);
                    changed += 1;
                    continue;
                }
            }
        }
        out.push_str(line);
    }

    (out, changed)
}

fn line_body_and_ending(line: &str) -> (&str, &str) {
    if let Some(body) = line.strip_suffix("\r\n") {
        (body, "\r\n")
    } else if let Some(body) = line.strip_suffix('\n') {
        (body, "\n")
    } else {
        (line, "")
    }
}

fn sidecar_sample_id(display_name: &str, well: &str) -> String {
    let name = display_name.trim();
    if well.is_empty() {
        return name.to_string();
    }
    if name == well {
        return String::new();
    }
    let Some(rest) = name.strip_prefix(well) else {
        return name.to_string();
    };
    rest.strip_prefix(':')
        .map(str::trim)
        .unwrap_or(name)
        .to_string()
}

/// Extract the size-calibration anchor times (scan numbers) from `.PKS`: the
/// longest strictly-increasing big-endian `u16` run whose values are all within
/// `(0, scans]`.
fn pks_anchor_times(pks: &[u8], scans: usize) -> Option<Vec<f64>> {
    let mut best: Vec<u16> = Vec::new();
    let mut i = 0usize;
    while i + 2 <= pks.len() {
        let mut j = i;
        let mut run: Vec<u16> = Vec::new();
        let mut prev = 0u16;
        while j + 2 <= pks.len() {
            let v = be_u16(pks, j);
            if v == 0 || v <= prev || v as usize > scans {
                break;
            }
            run.push(v);
            prev = v;
            j += 2;
        }
        if run.len() > best.len() {
            best = run;
        }
        i = j.max(i + 2);
    }
    // The run is preceded by a LabVIEW-style count; when captured as the first
    // element (it is smaller than the first time), drop it.
    if best.len() >= 2 && best[0] as usize == best.len() - 1 {
        best.remove(0);
    }
    (best.len() >= 8).then(|| best.into_iter().map(|v| v as f64).collect())
}

/// Parse the per-well peak table from `.PKS`. Layout (big-endian), verified by
/// reverse engineering (see `docs/fa_format.md`):
///
/// ```text
/// u32 nwells
/// per well:
///   20-byte summary  { u16 lm_apex, f32 lm_rfu, f32 lm_raw_area,
///                      u16 um_apex, f32 um_rfu, f32 um_raw_area }
///   u32 npeaks;  npeaks × 26-byte record   (copy 1, used here)
///   u32 npeaks;  npeaks × 26-byte record   (copy 2, an aligned duplicate)
///   8-byte trailer
/// ```
/// Each 26-byte record is `u16 start, u16 apex, u16 end` then 5×`f32`, all
/// identified against the ProSize CSV:
/// `[raw_area, rfu/height, baseline_a, baseline_b, corrected_area]`. Only the
/// corrected area (and geometry) is used. Concentration/molarity are *not*
/// stored anywhere in the native run (verified by an exhaustive value search) —
/// ProSize computes them from area + the size standard, so the reader does too.
/// Peaks are kept within the well's `[lm_apex, um_apex]` window; the lm/um peaks
/// are labelled and assigned the ladder's end sizes (markers are 1 bp / 6000 bp
/// by definition), samples are sized from the calibration. Returns `None` on any
/// framing inconsistency, so a malformed `.PKS` leaves peaks empty (never fails
/// the load).
fn pks_peaks(pks: &[u8], scans: usize, calib: Option<&[f64]>) -> Option<Vec<Vec<Peak>>> {
    const REC: usize = 26;
    let u16a = |o: usize| -> Option<u16> { (o + 2 <= pks.len()).then(|| be_u16(pks, o)) };
    let u32a = |o: usize| -> Option<u32> {
        (o + 4 <= pks.len())
            .then(|| u32::from_be_bytes([pks[o], pks[o + 1], pks[o + 2], pks[o + 3]]))
    };
    let f32a = |o: usize| -> Option<f32> {
        (o + 4 <= pks.len())
            .then(|| f32::from_be_bytes([pks[o], pks[o + 1], pks[o + 2], pks[o + 3]]))
    };

    let mut pos = 0usize;
    let nwells = u32a(pos)? as usize;
    if nwells == 0 || nwells > 1024 {
        return None;
    }
    pos += 4;

    let mut out = Vec::with_capacity(nwells);
    for _ in 0..nwells {
        // Summary { u16 lm_apex, f32 lm_rfu, f32 _, u16 um_apex, f32 um_rfu, f32 _ }.
        // The lower/upper marker apex times bound the reported peaks.
        let lm_apex = u16a(pos)? as usize;
        let um_apex = u16a(pos + 10)? as usize;
        pos += 20;

        let mut chosen: Vec<Peak> = Vec::new();
        for copy in 0..2 {
            let npk = u32a(pos)? as usize;
            if npk > 4096 {
                return None;
            }
            pos += 4;
            for _ in 0..npk {
                let start = u16a(pos)? as f64;
                let apex = u16a(pos + 2)? as usize;
                let end = u16a(pos + 4)? as f64;
                if apex == 0 || apex > scans {
                    return None;
                }
                // Record: [u16 start,apex,end] then f32
                // [raw_area, rfu, baseline_a, baseline_b, corrected_area].
                let area = f32a(pos + 22)?;
                pos += REC;
                // Build peaks only from the first copy, and only markers/samples
                // inside the marker window.
                if copy == 0 && apex >= lm_apex && apex <= um_apex {
                    let (observations, length) = if apex == lm_apex {
                        ("Lower Marker".to_string(), LADDER_BP[0])
                    } else if apex == um_apex {
                        ("Upper Marker".to_string(), LADDER_BP[LADDER_BP.len() - 1])
                    } else {
                        let len = calib.map_or(f64::NAN, |t| interp_bp(apex as f64, t, &LADDER_BP));
                        (String::new(), len)
                    };
                    let concentration = if apex == lm_apex {
                        LADDER_CONC[0]
                    } else if apex == um_apex {
                        LADDER_CONC[LADDER_CONC.len() - 1]
                    } else {
                        f64::NAN
                    };
                    chosen.push(Peak {
                        observations,
                        length,
                        time: apex as f64,
                        aligned_time: apex as f64,
                        start_time: start,
                        end_time: end,
                        aligned_start_time: start,
                        aligned_end_time: end,
                        area: area as f64,
                        concentration,
                        molarity: f64::NAN,
                    });
                }
            }
        }
        pos += 8; // trailer
        if chosen.len() >= 8 {
            for p in &mut chosen {
                if let Some(idx) = ladder_index_for_length(p.length) {
                    p.concentration = LADDER_CONC[idx];
                }
            }
        }
        out.push(chosen);
    }
    Some(out)
}

fn build_run(
    traces: RawTraces,
    caps: &[CapInfo],
    calib: Option<&[f64]>,
    peaks: &[Vec<Peak>],
    file_name: String,
) -> Result<Electrophoresis> {
    let ncap = traces.columns.len();
    // Per-scan length (bp) from the anchor times ↔ ladder, or empty if no calib.
    let lengths: Vec<f64> = match calib {
        Some(times) => (0..traces.scans)
            .map(|s| interp_bp(s as f64, times, &LADDER_BP))
            .collect(),
        None => Vec::new(),
    };

    let mut samples = Vec::with_capacity(ncap);
    for cap in 0..ncap {
        let time: Vec<f64> = (0..traces.scans).map(|s| s as f64).collect();
        let mut fluorescence: Vec<f32> = (0..traces.scans).map(|s| traces.value(cap, s)).collect();
        // Subtract the per-capillary dark baseline so traces start near zero,
        // matching the vendor's baseline-corrected export (the .raw carries a
        // large CCD DC offset). Uses a low percentile to resist outliers.
        let baseline = low_percentile(&fluorescence, 0.05);
        for v in &mut fluorescence {
            *v -= baseline;
        }
        let info = caps.get(cap);
        // Prefer the well label (e.g. "D1") as the display name; fall back to the
        // sample id. Well number is the capillary position (1-based).
        let name = match info {
            Some(c) if !c.well.is_empty() && !c.sample_id.is_empty() => {
                format!("{}: {}", c.well, c.sample_id)
            }
            Some(c) if !c.sample_id.is_empty() => c.sample_id.clone(),
            Some(c) => c.well.clone(),
            None => String::new(),
        };
        let well_number = (cap + 1) as i32;
        let well_peaks = peaks.get(cap).cloned().unwrap_or_default();
        // A well whose peaks span most of the ladder is the size standard.
        let is_ladder = well_peaks.len() >= 8;
        samples.push(Sample {
            well_number,
            name,
            category: String::new(),
            is_ladder,
            comment: info.map(|c| c.description.clone()).unwrap_or_default(),
            observations: String::new(),
            rin: None,
            time,
            fluorescence,
            aligned_time: (0..traces.scans).map(|s| s as f64).collect(),
            length: lengths.clone(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: well_peaks,
            regions: Vec::new(),
        });
    }

    let assay = AssayInfo {
        file_name,
        creation_date: traces.creation_date,
        assay_name: "Fragment Analyzer".to_string(),
        assay_type: "DNA".to_string(),
        length_unit: "bp".to_string(),
        concentration_unit: "ng/µl".to_string(),
        molarity_unit: Some("nmol/L".to_string()),
        has_upper_marker: true,
    };
    Ok(Electrophoresis {
        assay,
        ladder_peaks: standard_ladder_peaks(calib.is_some()),
        regions: Vec::new(),
        samples,
    })
}

fn standard_ladder_peaks(enabled: bool) -> Vec<LadderPeak> {
    if !enabled {
        return Vec::new();
    }
    LADDER_BP
        .iter()
        .zip(LADDER_CONC)
        .map(|(&size, concentration)| LadderPeak {
            size,
            area_a: f64::NAN,
            area_b: f64::NAN,
            concentration,
        })
        .collect()
}

fn ladder_index_for_length(length: f64) -> Option<usize> {
    LADDER_BP
        .iter()
        .position(|&bp| (length - bp).abs() <= bp.abs().max(1.0) * 1e-6)
}

fn compute_quantities(run: &mut Electrophoresis) {
    if crate::concentration::calculate_concentration(run).is_ok()
        && crate::concentration::calculate_molarity(run).is_ok()
    {
        fill_peak_quantities_from_points(run);
    }
}

fn fill_peak_quantities_from_points(run: &mut Electrophoresis) {
    for sample in &mut run.samples {
        for peak in &mut sample.peaks {
            let Some(i) = nearest_time_index(&sample.time, peak.time) else {
                continue;
            };
            if let Some(c) = sample
                .concentration
                .get(i)
                .copied()
                .filter(|v| v.is_finite())
            {
                peak.concentration = c;
            }
            if let Some(m) = sample.molarity.get(i).copied().filter(|v| v.is_finite()) {
                peak.molarity = m;
            }
        }
    }
}

fn nearest_time_index(time: &[f64], target: f64) -> Option<usize> {
    let mut best = None;
    let mut best_d = f64::INFINITY;
    for (i, &t) in time.iter().enumerate() {
        let d = (t - target).abs();
        if d < best_d {
            best = Some(i);
            best_d = d;
        }
    }
    best
}

/// Linear interpolation of scan → bp over the (times, bp) anchor points.
/// Returns NaN outside the calibrated range (matching the Bioanalyzer path).
fn interp_bp(scan: f64, times: &[f64], bp: &[f64]) -> f64 {
    if scan < times[0] || scan > times[times.len() - 1] {
        return f64::NAN;
    }
    for w in 0..times.len() - 1 {
        if scan >= times[w] && scan <= times[w + 1] {
            let f = (scan - times[w]) / (times[w + 1] - times[w]);
            return bp[w] + f * (bp[w + 1] - bp[w]);
        }
    }
    f64::NAN
}

fn be_u16(b: &[u8], off: usize) -> u16 {
    u16::from_be_bytes([b[off], b[off + 1]])
}

/// Value at the given fraction (0..1) of the sorted trace — a robust baseline
/// estimate (`frac = 0.05` ≈ the noise floor, resistant to a few low outliers).
fn low_percentile(v: &[f32], frac: f32) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s: Vec<f32> = v.to_vec();
    s.sort_by(|a, b| a.total_cmp(b));
    let idx = ((s.len() as f32 * frac) as usize).min(s.len() - 1);
    s[idx]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn put_u16(buf: &mut [u8], off: usize, value: u16) {
        buf[off..off + 2].copy_from_slice(&value.to_be_bytes());
    }

    fn synthetic_raw(width: usize, scans: usize, columns: &[u16]) -> Vec<u8> {
        let mut raw = vec![0u8; DATA_START + scans * width * 2];
        raw[..4].copy_from_slice(RAW_MAGIC);
        put_u16(&mut raw, 0xff, width as u16);
        raw[0x64..0x64 + b"Wednesday, November 19, 2025".len()]
            .copy_from_slice(b"Wednesday, November 19, 2025");
        raw[0x96..0x96 + b"4:19:23 PM".len()].copy_from_slice(b"4:19:23 PM");

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

    #[test]
    fn raw_parser_validates_magic_width_payload_and_capillary_table() {
        let raw = synthetic_raw(20, 3, &[5, 14]);
        let parsed = parse_raw(&raw, 2).expect("synthetic FA raw");

        assert_eq!(parsed.width, 20);
        assert_eq!(parsed.scans, 3);
        assert_eq!(
            parsed.creation_date,
            "Wednesday, November 19, 2025 4:19:23 PM"
        );
        assert_eq!(parsed.columns, vec![5, 14]);
        assert_eq!(parsed.data.len(), 60);
        assert_eq!(parsed.value(0, 2), 205.0);
        assert_eq!(parsed.value(1, 2), 214.0);

        let mut bad_magic = raw.clone();
        bad_magic[0] = b'X';
        assert!(parse_raw(&bad_magic, 2).is_err());

        let mut bad_width = raw.clone();
        put_u16(&mut bad_width, 0xff, 15);
        assert!(parse_raw(&bad_width, 2).is_err());

        let truncated = &raw[..raw.len() - 2];
        assert!(parse_raw(truncated, 2).is_err());
    }

    #[test]
    fn pks_anchor_count_is_removed_and_interpolation_bounds_are_nan() {
        let mut pks = vec![0xaa, 0xbb];
        for value in [
            16u16, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170,
        ] {
            pks.extend_from_slice(&value.to_be_bytes());
        }
        pks.extend_from_slice(&0u16.to_be_bytes());

        let anchors = pks_anchor_times(&pks, 200).expect("anchor times");
        assert_eq!(anchors.len(), LADDER_BP.len());
        assert_eq!(anchors[0], 20.0);
        assert_eq!(anchors[15], 170.0);

        assert!(interp_bp(19.0, &anchors, &LADDER_BP).is_nan());
        assert_eq!(interp_bp(20.0, &anchors, &LADDER_BP), 1.0);
        assert!((interp_bp(25.0, &anchors, &LADDER_BP) - 50.5).abs() < f64::EPSILON);
        assert_eq!(interp_bp(170.0, &anchors, &LADDER_BP), 6000.0);
        assert!(interp_bp(171.0, &anchors, &LADDER_BP).is_nan());
    }

    #[test]
    fn txt_sidecar_patch_rewrites_sample_ids_and_strips_well_prefix() {
        let text = "Capillary #: 1\r\nWell: D1\r\nSample ID: old one\r\n\
                    Capillary #: 2\r\nWell: D2\r\nSample ID: old two\r\n";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![make_sample(1, "D1: new one"), make_sample(2, "custom two")],
        };

        let (patched, changed) = patch_txt_names(text, &run);

        assert_eq!(changed, 2);
        assert!(patched.contains("Sample ID: new one\r\n"));
        assert!(patched.contains("Sample ID: custom two\r\n"));
        assert!(!patched.contains("old one"));
        assert!(!patched.contains("old two"));
    }

    #[test]
    fn txt_sidecar_patch_preserves_blank_sample_id_for_well_only_name() {
        let text = "Capillary #: 1\nWell: D1\nSample ID: \n";
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![make_sample(1, "D1")],
        };

        let (patched, changed) = patch_txt_names(text, &run);

        assert_eq!(changed, 1);
        assert_eq!(patched, "Capillary #: 1\nWell: D1\nSample ID: \n");
    }

    #[test]
    fn txt_sidecar_parse_keeps_sample_descriptions() {
        let caps = parse_txt(
            "Capillary #: 1\nWell: D1\nSample ID: alpha\nDescription: note one\n\
             Capillary #: 2\nWell: D2\nSample ID: beta\nDescription: note two\n",
        );

        assert_eq!(caps.len(), 2);
        assert_eq!(caps[0].well, "D1");
        assert_eq!(caps[0].sample_id, "alpha");
        assert_eq!(caps[0].description, "note one");
        assert_eq!(caps[1].description, "note two");
    }

    #[test]
    fn fa_path_detection_accepts_uppercase_raw_extension() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_upper_raw_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("RUN.RAW");
        std::fs::write(&raw, b"FA\0\0").unwrap();

        assert!(is_fa_path(&raw));
        assert!(is_fa_path(&dir));

        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn filesystem_sidecars_are_found_case_insensitively() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_sidecar_case_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("RUN.RAW");
        std::fs::write(&raw, synthetic_raw(20, 3, &[5, 14])).unwrap();
        std::fs::write(
            dir.join("run.TxT"),
            "Capillary #: 1\nWell: D1\nSample ID: alpha\n\
             Capillary #: 2\nWell: D2\nSample ID: beta\n",
        )
        .unwrap();

        let run = read_fa_run(&raw).unwrap();

        assert_eq!(run.samples[0].name, "D1: alpha");
        assert_eq!(run.samples[1].name, "D2: beta");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn build_run_maps_raw_creation_date_and_txt_descriptions() {
        let raw = synthetic_raw(20, 3, &[5, 14]);
        let txt = "Capillary #: 1\nWell: D1\nSample ID: alpha\nDescription: first note\n\
                   Capillary #: 2\nWell: D2\nSample ID: beta\nDescription: second note\n";

        let run = build_fa_run(&raw, None, Some(txt), "run.raw".to_string()).unwrap();

        assert_eq!(
            run.assay.creation_date,
            "Wednesday, November 19, 2025 4:19:23 PM"
        );
        assert_eq!(run.samples[0].comment, "first note");
        assert_eq!(run.samples[1].comment, "second note");
    }

    #[test]
    fn metadata_parses_ini_current_logs_and_header() {
        let txt = "Raw file: C:\\AATI\\Data\\run.raw\nOperator: admin\nNotes: ok\n\
                   ============================================================\n\
                   Capillary #: 1\nWell: D1\nSample ID: alpha\n";
        let method = "[Separation]\nKV=7.00\nTime minutes=25.00\n";
        let anai = "[Size Calibration]\nCalibration Ladder = \"1,100,200\"\nColumn = 11\n";
        let current = "Current(uA)\tVoltage(kV)\tPressure(PSI)\n54.0\t6.99\t-0.007\n";

        let header = parse_txt_header(txt);
        let method = parse_ini(method);
        let anai = parse_ini(anai);
        let current = parse_current_log(current).unwrap();

        assert_eq!(header["Operator"], "admin");
        assert_eq!(header["Raw file"], "C:\\AATI\\Data\\run.raw");
        assert_eq!(method["Separation"]["KV"], "7.00");
        assert_eq!(anai["Size Calibration"]["Calibration Ladder"], "1,100,200");
        assert_eq!(
            current[0],
            FaCurrentPoint {
                current_ua: 54.0,
                voltage_kv: 6.99,
                pressure_psi: -0.007,
            }
        );
    }

    #[test]
    fn filesystem_metadata_reads_known_sidecars() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_metadata_fs_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        std::fs::write(&raw, synthetic_raw(20, 3, &[5, 14])).unwrap();
        std::fs::write(
            dir.join("run.txt"),
            "Raw file: C:\\AATI\\Data\\run.raw\nOperator: admin\n\
             Capillary #: 1\nWell: D1\nSample ID: alpha\n",
        )
        .unwrap();
        std::fs::write(dir.join("method.mthd"), "[Separation]\nKV=7.00\n").unwrap();
        std::fs::write(dir.join("run.ANAI"), "[General]\nMode = 2\n").unwrap();
        std::fs::write(
            dir.join("run.current"),
            "Current(uA)\tVoltage(kV)\tPressure(PSI)\n54.0\t7.00\t-0.01\n",
        )
        .unwrap();
        std::fs::write(dir.join("Timing.txt"), "Frame 1: Start\n").unwrap();
        std::fs::write(dir.join("ExpTime.txt"), "845ms, 50\n").unwrap();
        std::fs::write(dir.join("CameraImage.bmp"), b"\x89PNG\r\n").unwrap();

        let meta = read_fa_metadata(&raw).unwrap();

        assert_eq!(meta.run_header["Operator"], "admin");
        assert_eq!(meta.method["Separation"]["KV"], "7.00");
        assert_eq!(meta.analysis["General"]["Mode"], "2");
        assert_eq!(meta.current.len(), 1);
        assert_eq!(meta.timing, vec!["Frame 1: Start"]);
        assert_eq!(meta.exposure, vec!["845ms, 50"]);
        assert_eq!(meta.camera_image.as_deref(), Some(&b"\x89PNG\r\n"[..]));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn filesystem_raw_ambiguity_is_rejected() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_raw_ambiguous_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.raw"), b"FA\0\0").unwrap();
        std::fs::write(dir.join("b.RAW"), b"FA\0\0").unwrap();

        let err = resolve_raw_path(&dir).unwrap_err().to_string();

        assert!(err.contains("ambiguous FA .raw files"), "got {err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn filesystem_sidecar_save_uses_case_insensitive_unique_match() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_sidecar_save_case_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("RUN.RAW");
        let txt = dir.join("run.TxT");
        std::fs::write(&raw, synthetic_raw(20, 3, &[5, 14])).unwrap();
        std::fs::write(&txt, "Capillary #: 1\nWell: D1\nSample ID: old\n").unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![make_sample(1, "D1: renamed")],
        };

        save_txt_names(&raw, &run).unwrap();

        let patched = std::fs::read_to_string(&txt).unwrap();
        assert!(patched.contains("Sample ID: renamed"), "got {patched}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn filesystem_sidecar_ambiguity_is_rejected() {
        let dir = std::env::temp_dir().join(format!(
            "traceio_fa_sidecar_ambiguous_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("RUN.RAW");
        std::fs::write(&raw, synthetic_raw(20, 3, &[5, 14])).unwrap();
        std::fs::write(
            dir.join("run.TxT"),
            "Capillary #: 1\nWell: D1\nSample ID: one\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("Run.txt"),
            "Capillary #: 1\nWell: D1\nSample ID: two\n",
        )
        .unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![make_sample(1, "D1: renamed")],
        };

        let err = save_txt_names(&raw, &run).unwrap_err().to_string();

        assert!(err.contains("ambiguous FA .txt sidecars"), "got {err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn zip_read_targets_magic_raw_sidecar_not_first_raw_name() {
        use std::io::Write;

        let zip_path = std::env::temp_dir().join(format!(
            "traceio_fa_decoy_raw_read_{}_{}.zip",
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
            zip.start_file("decoy.raw", opts).unwrap();
            zip.write_all(b"not-fragment-analyzer").unwrap();
            zip.start_file("decoy.txt", opts).unwrap();
            zip.write_all(b"Capillary #: 1\nWell: Z1\nSample ID: decoy\n")
                .unwrap();
            zip.start_file("RUN.RAW", opts).unwrap();
            zip.write_all(&synthetic_raw(20, 3, &[5, 14])).unwrap();
            zip.start_file("run.TxT", opts).unwrap();
            zip.write_all(b"Capillary #: 1\nWell: D1\nSample ID: real\n")
                .unwrap();
            zip.finish().unwrap();
        }

        let run = read_fa_run(&zip_path).unwrap();

        assert_eq!(run.samples[0].name, "D1: real");
        std::fs::remove_file(zip_path).unwrap();
    }

    #[test]
    fn zip_metadata_reads_known_sidecars() {
        use std::io::Write;

        let zip_path = std::env::temp_dir().join(format!(
            "traceio_fa_metadata_zip_{}_{}.zip",
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
            zip.start_file("folder/run.raw", opts).unwrap();
            zip.write_all(&synthetic_raw(20, 3, &[5, 14])).unwrap();
            zip.start_file("folder/run.txt", opts).unwrap();
            zip.write_all(b"Raw file: C:\\AATI\\Data\\run.raw\nOperator: admin\n")
                .unwrap();
            zip.start_file("folder/run.ANAI", opts).unwrap();
            zip.write_all(b"[General]\nMode = 2\n").unwrap();
            zip.start_file("folder/run.current", opts).unwrap();
            zip.write_all(b"Current(uA)\tVoltage(kV)\tPressure(PSI)\n54.0\t7.00\t-0.01\n")
                .unwrap();
            zip.start_file("folder/method.mthd", opts).unwrap();
            zip.write_all(b"[Separation]\nKV=7.00\n").unwrap();
            zip.start_file("folder/Timing.txt", opts).unwrap();
            zip.write_all(b"Frame 1: Start\n").unwrap();
            zip.start_file("folder/ExpTime.txt", opts).unwrap();
            zip.write_all(b"845ms, 50\n").unwrap();
            zip.start_file("folder/CameraImage.bmp", opts).unwrap();
            zip.write_all(b"\x89PNG\r\n").unwrap();
            zip.finish().unwrap();
        }

        let meta = read_fa_metadata(&zip_path).unwrap();

        assert_eq!(meta.run_header["Operator"], "admin");
        assert_eq!(meta.method["Separation"]["KV"], "7.00");
        assert_eq!(meta.analysis["General"]["Mode"], "2");
        assert_eq!(meta.current.len(), 1);
        assert_eq!(meta.timing, vec!["Frame 1: Start"]);
        assert_eq!(meta.exposure, vec!["845ms, 50"]);
        assert_eq!(meta.camera_image.as_deref(), Some(&b"\x89PNG\r\n"[..]));
        std::fs::remove_file(zip_path).unwrap();
    }

    #[test]
    fn zip_save_targets_magic_raw_sidecar_not_first_raw_name() {
        use std::io::{Read, Write};

        let zip_path = std::env::temp_dir().join(format!(
            "traceio_fa_decoy_raw_{}_{}.zip",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let stale_tmp = zip_path.with_extension("zip.tmp");
        {
            let out = std::fs::File::create(&zip_path).unwrap();
            let mut zip = zip::ZipWriter::new(out);
            let opts = zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Deflated);
            zip.start_file("decoy.raw", opts).unwrap();
            zip.write_all(b"not-fragment-analyzer").unwrap();
            zip.start_file("decoy.txt", opts).unwrap();
            zip.write_all(b"Capillary #: 1\nWell: Z1\nSample ID: decoy\n")
                .unwrap();
            zip.start_file("real.raw", opts).unwrap();
            zip.write_all(b"FA\0\0").unwrap();
            zip.start_file("real.txt", opts).unwrap();
            zip.write_all(b"Capillary #: 1\nWell: D1\nSample ID: old\n")
                .unwrap();
            zip.finish().unwrap();
        }
        std::fs::write(&stale_tmp, b"do not overwrite").unwrap();
        let run = Electrophoresis {
            assay: Default::default(),
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![make_sample(1, "D1: renamed")],
        };

        save_txt_names(&zip_path, &run).unwrap();

        let file = std::fs::File::open(&zip_path).unwrap();
        let mut zip = zip::ZipArchive::new(file).unwrap();
        let mut decoy = String::new();
        zip.by_name("decoy.txt")
            .unwrap()
            .read_to_string(&mut decoy)
            .unwrap();
        let mut real = String::new();
        zip.by_name("real.txt")
            .unwrap()
            .read_to_string(&mut real)
            .unwrap();
        assert!(decoy.contains("Sample ID: decoy"), "got {decoy}");
        assert!(real.contains("Sample ID: renamed"), "got {real}");
        assert_eq!(std::fs::read(&stale_tmp).unwrap(), b"do not overwrite");
        std::fs::remove_file(zip_path).unwrap();
        std::fs::remove_file(stale_tmp).unwrap();
    }

    fn make_sample(well: i32, name: &str) -> Sample {
        Sample {
            well_number: well,
            name: name.to_string(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: Vec::new(),
            fluorescence: Vec::new(),
            aligned_time: Vec::new(),
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
            regions: Vec::new(),
        }
    }
}
