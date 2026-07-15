//! Application state shared across Slint callbacks.

use traceio::xad::RawChannel;
use traceio::Electrophoresis;

use crate::plot::{Viewport, YMode};

pub struct AppState {
    pub run: Electrophoresis,
    /// Raw detector channels (populated for native `.xad`; see `raw_mode`).
    pub raw_channels: Vec<RawChannel>,
    /// Index of the currently focused entry (a sample, or a raw channel in
    /// raw mode).
    pub selected: usize,
    pub y_mode: YMode,
    /// Zoom/pan window; `None` means auto-fit the current entry.
    pub viewport: Option<Viewport>,
    /// Last load error, shown in the UI.
    pub error: Option<String>,
}

impl AppState {
    pub fn new(run: Electrophoresis, raw_channels: Vec<RawChannel>) -> Self {
        AppState {
            run,
            raw_channels,
            selected: 0,
            y_mode: YMode::Fluorescence,
            viewport: None,
            error: None,
        }
    }

    /// True when there is no processed per-well signal but raw detector channels
    /// are available (a native `.xad`): the UI lists/plots the raw channels.
    pub fn raw_mode(&self) -> bool {
        let no_processed = self.run.samples.iter().all(|s| s.fluorescence.is_empty());
        no_processed && !self.raw_channels.is_empty()
    }

    /// Header summary line for the window.
    pub fn title(&self) -> String {
        let mode = if self.raw_mode() { "  [raw acquisition]" } else { "" };
        format!(
            "{}  —  {} ({}),  {} samples{mode}",
            self.run.assay.file_name,
            self.run.assay.assay_name,
            self.run.assay.assay_type,
            self.run.samples.len()
        )
    }

    /// Number of selectable entries in the current mode.
    pub fn entry_count(&self) -> usize {
        if self.raw_mode() {
            self.raw_channels.len()
        } else {
            self.run.samples.len()
        }
    }

    /// Display labels for the entry list.
    pub fn entry_labels(&self) -> Vec<String> {
        if self.raw_mode() {
            return self
                .raw_channels
                .iter()
                .map(|c| format!("{} (raw)", c.channel_id))
                .collect();
        }
        self.run
            .samples
            .iter()
            .map(|s| {
                let label = if s.name.is_empty() {
                    format!("Well {}", s.well_number)
                } else {
                    format!("{}: {}", s.well_number, s.name)
                };
                if s.is_ladder {
                    format!("{label}  [ladder]")
                } else {
                    label
                }
            })
            .collect()
    }
}
