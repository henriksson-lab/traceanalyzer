//! `traceio` — readers for Agilent automated-electrophoresis data.
//!
//! Coverage status:
//! - **Bioanalyzer 2100**: native `.xad` ([`xad`]) and exported XML
//!   ([`bioanalyzer`]). Inner-XML parsing validated against real demo files;
//!   the native `.xad` container unwrap is a faithful but not-yet-validated
//!   port of grimbough/bioanalyzeR (needs real `.xad` samples).
//! - **TapeStation**: exported reader ([`tapestation`]) — metadata `.xml` +
//!   `_Electropherogram.csv` (a port of jwfoley/bioanalyzeR `tapestation.R`),
//!   with per-sample regions, multi-tape (`ScreenTapeID`) ladder grouping, and
//!   marker-relative sizing. Validated against Agilent demo exports across every
//!   assay (D1000/D5000, HS variants, cfDNA, gDNA, RNA). The native project file
//!   (`.D1000`/`.RNA`/…) is a password-encrypted ZIP with no public key, so —
//!   like bioanalyzeR — only the export is read (per-point concentration/molarity
//!   are not in the export and are not fabricated; see `docs/tapestation_format.md`).
//! - **Fragment Analyzer (AATI)**: native run reader ([`fa`]) — reverse-engineered
//!   from a real run directory (see `docs/fa_format.md`). Reads the `.raw` CCD
//!   acquisition into 12 per-capillary electropherograms, plus the `.PKS`
//!   size-calibration anchors and per-well peak table (sizes, areas, lower/upper
//!   markers, ladder detection). Computes concentration/molarity from the
//!   standard FA ladder setpoints and marker-area scaling. Validated against the
//!   ProSize CSV export.
//!
//! Most applications should start with [`io::detect_format`], [`io::read_path`],
//! and [`io::save_path`]. The format modules remain available for callers that
//! already know their input type, and [`save::save_run`] remains available when
//! the source path and run are tracked separately.

pub mod bioanalyzer;
pub mod calibration;
pub mod concentration;
pub mod fa;
pub mod io;
pub mod model;
pub mod save;
pub mod tapestation;
pub mod xad;

pub use model::*;
