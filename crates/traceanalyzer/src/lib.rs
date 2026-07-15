//! Non-UI logic for the Trace analyzer viewer, exposed as a library so it can be
//! exercised headlessly (e.g. the `render_png` example) and unit-tested. The
//! Slint UI itself lives in the binary (`main.rs`).

pub mod gel;
pub mod loading;
pub mod overview;
pub mod plot;
pub mod render;
pub mod state;
pub mod table;
