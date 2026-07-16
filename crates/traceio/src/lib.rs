//! `traceio` — readers for Agilent automated-electrophoresis data.
//!
//! Coverage status:
//! - **Bioanalyzer 2100**: native `.xad` ([`xad`]) and exported XML
//!   ([`bioanalyzer`]). Inner-XML parsing validated against real demo files;
//!   the native `.xad` container unwrap is a faithful but not-yet-validated
//!   port of grimbough/bioanalyzeR (needs real `.xad` samples).
//! - **TapeStation**: TODO — native container is an encrypted ZIP (blocked);
//!   plan is XML + unaligned-CSV export.
//! - **Fragment Analyzer (AATI)**: native run reader ([`fa`]) — reverse-engineered
//!   from a real run directory (see `docs/fa_format.md`). Reads the `.raw` CCD
//!   acquisition into 12 per-capillary electropherograms, plus the `.PKS`
//!   size-calibration anchors and per-well peak table (sizes, areas, lower/upper
//!   markers, ladder detection). Computes concentration/molarity from the
//!   standard FA ladder setpoints and marker-area scaling. Validated against the
//!   ProSize CSV export.

pub mod bioanalyzer;
pub mod calibration;
pub mod concentration;
pub mod fa;
pub mod model;
pub mod save;
pub mod xad;

pub use model::*;
