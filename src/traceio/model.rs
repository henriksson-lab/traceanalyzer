//! Common in-memory data model for automated-electrophoresis runs.
//!
//! This mirrors the `electrophoresis` object of the R package `bioanalyzeR`
//! (jwfoley/bioanalyzeR, MIT) closely enough that its parsing logic ports
//! one-to-one, but is intended to be shared across all instrument families
//! (Bioanalyzer, TapeStation, Fragment Analyzer).

/// One row of the assay's calibration ladder (known size + concentration).
///
/// From the Bioanalyzer `<LadderPeak>` element:
/// `Size`, `AreaAType`, `AreaBType`, `Concentration`.
#[derive(Debug, Clone)]
pub struct LadderPeak {
    pub size: f64,
    pub area_a: f64,
    pub area_b: f64,
    pub concentration: f64,
}

/// A detected/reported peak within one sample.
///
/// Fields follow `<PeakMolecular>`. Missing numeric fields are stored as NaN.
#[derive(Debug, Clone)]
pub struct Peak {
    pub observations: String,
    pub length: f64,
    pub time: f64,
    pub aligned_time: f64,
    pub start_time: f64,
    pub end_time: f64,
    pub aligned_start_time: f64,
    pub aligned_end_time: f64,
    pub area: f64,
    pub concentration: f64,
    pub molarity: f64,
}

/// A user/assay-defined smear region (integrated size window).
#[derive(Debug, Clone)]
pub struct Region {
    pub lower_length: f64,
    pub upper_length: f64,
}

/// One well/lane: metadata, the raw fluorescence trace, and its peaks.
#[derive(Debug, Clone)]
pub struct Sample {
    pub well_number: i32,
    pub name: String,
    pub category: String,
    pub is_ladder: bool,
    pub comment: String,
    pub observations: String,
    /// RNA Integrity Number, present only for RNA assays.
    pub rin: Option<f64>,
    /// Trace time axis (same length as `fluorescence`).
    pub time: Vec<f64>,
    /// Processed fluorescence signal (stored as float32 in the source file).
    pub fluorescence: Vec<f32>,
    /// Marker-aligned time per point (filled by calibration; empty until then).
    pub aligned_time: Vec<f64>,
    /// Molecular length (bp/nt) per point from the ladder standard curve
    /// (filled by calibration; NaN outside the ladder's interpolation range).
    pub length: Vec<f64>,
    /// Mass concentration per point (e.g. ng/µl or pg/µl), from
    /// `concentration::calculate_concentration`. Empty until then; the first
    /// point of each sample is NaN (no trapezoid to its left).
    pub concentration: Vec<f64>,
    /// Molarity per point (nM or pM), from `concentration::calculate_molarity`.
    /// Empty until then; NaN wherever `length` or `concentration` is NaN.
    pub molarity: Vec<f64>,
    pub peaks: Vec<Peak>,
    /// Per-sample smear/analysis regions (TapeStation). Empty for instruments
    /// whose regions are assay-level (Bioanalyzer keeps those in
    /// [`Electrophoresis::regions`]).
    pub regions: Vec<Region>,
}

impl Sample {
    /// Min/max of the fluorescence trace, or None if empty.
    pub fn fluorescence_range(&self) -> Option<(f32, f32)> {
        let mut it = self.fluorescence.iter().copied().filter(|v| v.is_finite());
        let first = it.next()?;
        let (mut lo, mut hi) = (first, first);
        for v in it {
            lo = lo.min(v);
            hi = hi.max(v);
        }
        Some((lo, hi))
    }
}

/// Per-run assay metadata.
#[derive(Debug, Clone, Default)]
pub struct AssayInfo {
    pub file_name: String,
    pub creation_date: String,
    pub assay_name: String,
    pub assay_type: String,
    /// Molecular-length unit, e.g. "bp" or "nt".
    pub length_unit: String,
    /// Concentration unit, e.g. "ng/µl".
    pub concentration_unit: String,
    /// Derived molarity unit ("nM"/"pM"), if determinable.
    pub molarity_unit: Option<String>,
    pub has_upper_marker: bool,
}

/// A fully parsed electrophoresis run.
#[derive(Debug, Clone)]
pub struct Electrophoresis {
    pub assay: AssayInfo,
    pub ladder_peaks: Vec<LadderPeak>,
    pub regions: Vec<Region>,
    pub samples: Vec<Sample>,
}

impl Electrophoresis {
    /// Index of the ladder sample, if exactly one is present.
    pub fn ladder_index(&self) -> Option<usize> {
        let ladders: Vec<usize> = self
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| s.is_ladder)
            .map(|(i, _)| i)
            .collect();
        match ladders.as_slice() {
            [i] => Some(*i),
            _ => None,
        }
    }
}
