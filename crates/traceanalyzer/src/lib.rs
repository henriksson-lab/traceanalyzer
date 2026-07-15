//! Non-UI logic for the traceanalyzer viewer, exposed as a library so it can be
//! exercised headlessly (e.g. the `render_png` example) and unit-tested. The
//! Slint UI itself lives in the binary (`main.rs`).

pub mod loading;
pub mod plot;
pub mod render;
pub mod state;
