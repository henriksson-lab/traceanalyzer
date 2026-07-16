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

use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};

use crate::model::{AssayInfo, Electrophoresis, Peak, Sample};

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

/// One capillary's identity from the `.txt` sidecar.
struct CapInfo {
    well: String,
    sample_id: String,
}

/// Read a Fragment Analyzer run. `path` may be the `.raw` file itself or the run
/// directory (in which case the single `.raw` inside is used).
pub fn read_fa_run(path: &Path) -> Result<Electrophoresis> {
    let raw_path = resolve_raw_path(path)?;
    let stem = raw_path
        .file_stem()
        .ok_or_else(|| anyhow!("FA run: cannot derive file stem from {}", raw_path.display()))?;
    let dir = raw_path.parent().unwrap_or_else(|| Path::new("."));
    let sibling = |ext: &str| dir.join(format!("{}.{ext}", stem.to_string_lossy()));

    let caps = read_txt(&sibling("txt")).unwrap_or_default();
    let raw = std::fs::read(&raw_path).with_context(|| format!("reading {}", raw_path.display()))?;
    let traces = parse_raw(&raw, caps.len())?;

    // Size calibration and peaks both come from `.PKS`.
    let pks = std::fs::read(sibling("PKS")).ok();
    let anchor_times = pks.as_deref().and_then(|d| pks_anchor_times(d, traces.scans));
    let calib = anchor_times.filter(|t| t.len() == LADDER_BP.len());
    let peaks = pks
        .as_deref()
        .and_then(|d| pks_peaks(d, traces.scans, calib.as_deref()))
        .unwrap_or_default();

    let file_name = raw_path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    build_run(traces, &caps, calib.as_deref(), &peaks, file_name)
}

/// True if `path` looks like a Fragment Analyzer run entry point (a `.raw` file
/// with the `FA\0\0` magic, or a directory containing one).
pub fn is_fa_path(path: &Path) -> bool {
    if path.is_dir() {
        return find_raw_in_dir(path).is_some();
    }
    if path.extension().and_then(|e| e.to_str()) != Some("raw") {
        return false;
    }
    let mut buf = [0u8; 4];
    matches!(read_prefix(path, &mut buf), Ok(())) && buf == RAW_MAGIC
}

fn resolve_raw_path(path: &Path) -> Result<PathBuf> {
    if path.is_dir() {
        return find_raw_in_dir(path)
            .ok_or_else(|| anyhow!("no .raw file found in FA run dir {}", path.display()));
    }
    Ok(path.to_path_buf())
}

fn find_raw_in_dir(dir: &Path) -> Option<PathBuf> {
    let mut buf = [0u8; 4];
    std::fs::read_dir(dir).ok()?.flatten().map(|e| e.path()).find(|p| {
        p.extension().and_then(|e| e.to_str()) == Some("raw")
            && read_prefix(p, &mut buf).is_ok()
            && buf == RAW_MAGIC
    })
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
    Ok(RawTraces { scans, width, data, columns })
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

/// Parse the `.txt` sidecar into per-capillary (well, sample id), in order.
fn read_txt(path: &Path) -> Result<Vec<CapInfo>> {
    let text = std::fs::read_to_string(path)
        .or_else(|_| std::fs::read(path).map(|b| String::from_utf8_lossy(&b).into_owned()))
        .with_context(|| format!("reading {}", path.display()))?;
    let mut caps = Vec::new();
    let (mut well, mut sample) = (String::new(), String::new());
    let mut in_cap = false;
    let flush = |caps: &mut Vec<CapInfo>, w: &mut String, s: &mut String, in_cap: &mut bool| {
        if *in_cap {
            caps.push(CapInfo { well: std::mem::take(w), sample_id: std::mem::take(s) });
        }
        *in_cap = false;
    };
    for line in text.lines() {
        let line = line.trim();
        if let Some(_rest) = line.strip_prefix("Capillary #:") {
            flush(&mut caps, &mut well, &mut sample, &mut in_cap);
            in_cap = true;
        } else if let Some(rest) = line.strip_prefix("Well:") {
            well = rest.trim().to_string();
        } else if let Some(rest) = line.strip_prefix("Sample ID:") {
            sample = rest.trim().to_string();
        }
    }
    flush(&mut caps, &mut well, &mut sample, &mut in_cap);
    Ok(caps)
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
///   20-byte summary  { u16 lm_apex, f32 lm_rfu, f32 _, u16 um_apex, f32 um_rfu, f32 _ }
///   u32 npeaks;  npeaks × 26-byte record   (copy 1, used here)
///   u32 npeaks;  npeaks × 26-byte record   (copy 2, an aligned duplicate)
///   8-byte trailer
/// ```
/// Each 26-byte record is `u16 start, u16 apex, u16 end` then 5×`f32`
/// (`[_, rfu, _, _, corrected_area]`). Peaks are kept within the well's
/// `[lm_apex, um_apex]` window; the lm/um peaks are labelled and assigned the
/// ladder's end sizes (markers are 1 bp / 6000 bp by definition), samples are
/// sized from the calibration. Returns `None` on any framing inconsistency, so
/// a malformed `.PKS` simply leaves peaks empty without failing the load.
fn pks_peaks(pks: &[u8], scans: usize, calib: Option<&[f64]>) -> Option<Vec<Vec<Peak>>> {
    const REC: usize = 26;
    let u16a = |o: usize| -> Option<u16> { (o + 2 <= pks.len()).then(|| be_u16(pks, o)) };
    let u32a = |o: usize| -> Option<u32> {
        (o + 4 <= pks.len()).then(|| u32::from_be_bytes([pks[o], pks[o + 1], pks[o + 2], pks[o + 3]]))
    };
    let f32a = |o: usize| -> Option<f32> {
        (o + 4 <= pks.len()).then(|| f32::from_be_bytes([pks[o], pks[o + 1], pks[o + 2], pks[o + 3]]))
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
                // Record: [u16 start, apex, end] then f32 [_, rfu, _, _, area].
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
                    chosen.push(Peak {
                        observations,
                        length,
                        time: apex as f64,
                        aligned_time: f64::NAN,
                        start_time: start,
                        end_time: end,
                        aligned_start_time: f64::NAN,
                        aligned_end_time: f64::NAN,
                        area: area as f64,
                        concentration: f64::NAN,
                        molarity: f64::NAN,
                    });
                }
            }
        }
        pos += 8; // trailer
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
        Some(times) => (0..traces.scans).map(|s| interp_bp(s as f64, times, &LADDER_BP)).collect(),
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
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time,
            fluorescence,
            aligned_time: Vec::new(),
            length: lengths.clone(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: well_peaks,
        });
    }

    let assay = AssayInfo {
        file_name,
        creation_date: String::new(),
        assay_name: "Fragment Analyzer".to_string(),
        assay_type: "DNA".to_string(),
        length_unit: "bp".to_string(),
        concentration_unit: "ng/µl".to_string(),
        molarity_unit: Some("nmol/L".to_string()),
        has_upper_marker: true,
    };
    Ok(Electrophoresis { assay, ladder_peaks: Vec::new(), regions: Vec::new(), samples })
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
