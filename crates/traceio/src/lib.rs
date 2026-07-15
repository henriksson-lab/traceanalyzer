//! `traceio` — readers for Agilent automated-electrophoresis data.
//!
//! Coverage status:
//! - **Bioanalyzer 2100**: native `.xad` ([`xad`]) and exported XML
//!   ([`bioanalyzer`]). Inner-XML parsing validated against real demo files;
//!   the native `.xad` container unwrap is a faithful but not-yet-validated
//!   port of grimbough/bioanalyzeR (needs real `.xad` samples).
//! - **TapeStation**: TODO — native container is an encrypted ZIP (blocked);
//!   plan is XML + unaligned-CSV export.
//! - **Fragment Analyzer / ProSize**: TODO — native `.raw`/`.db3` is likely
//!   SQLite (to be confirmed on a sample); plan is CSV export first.

pub mod bioanalyzer;
pub mod calibration;
pub mod concentration;
pub mod model;
pub mod xad;

pub use model::*;
