//! Application state shared across Slint callbacks.

use std::collections::HashMap;

use traceio::calibration::MarkerOverride;
use traceio::xad::RawChannel;
use traceio::Electrophoresis;

use crate::plot::{Viewport, YMode};

/// Which marker line is being dragged in marker-edit mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Marker {
    Lower,
    Upper,
}

/// Marker drag state, including the override to restore if validation or
/// recalibration fails.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MarkerDrag {
    pub sample_idx: usize,
    pub marker: Marker,
    pub original_override: Option<MarkerOverride>,
}

pub struct AppState {
    pub run: Electrophoresis,
    /// Raw detector channels (populated for native `.xad`; see `raw_mode`).
    pub raw_channels: Vec<RawChannel>,
    /// Selected entry indices, in click order. The last element is the primary
    /// (focused) entry; always non-empty. In raw mode this holds a single
    /// channel index. Multiple entries overlay in the Detail plot.
    pub selection: Vec<usize>,
    /// Anchor for shift-range selection (index of the last plain/ctrl click).
    pub anchor: usize,
    /// Peak-height-normalize overlaid traces (shape comparison).
    pub normalize: bool,
    /// Overview tab: share one y-scale across all mini-plots (vs. per-plot fit).
    pub overview_shared_y: bool,
    /// Overview tab: show the virtual gel instead of the small-multiples grid.
    pub overview_gel: bool,
    pub y_mode: YMode,
    /// Zoom/pan window; `None` means auto-fit the current entry.
    pub viewport: Option<Viewport>,
    /// x-position of a cross-highlighted peak (table→plot), if any.
    pub highlight_x: Option<f64>,
    /// Plot x-position for each current table row (`None` for region rows),
    /// parallel to the rows pushed to the Table tab.
    pub table_peak_x: Vec<Option<f64>>,
    /// Marker-edit mode: plot in raw migration time with draggable lower/upper
    /// marker lines.
    pub marker_edit: bool,
    /// Manual marker overrides, per sample index (empty = fully automatic).
    pub overrides: HashMap<usize, MarkerOverride>,
    /// Marker currently being dragged, if any (transient).
    pub grabbed: Option<MarkerDrag>,
    /// Last load error, shown in the UI.
    pub error: Option<String>,
}

impl AppState {
    pub fn new(run: Electrophoresis, raw_channels: Vec<RawChannel>, error: Option<String>) -> Self {
        AppState {
            run,
            raw_channels,
            selection: vec![0],
            anchor: 0,
            normalize: false,
            overview_shared_y: true,
            overview_gel: false,
            y_mode: YMode::Fluorescence,
            viewport: None,
            highlight_x: None,
            table_peak_x: Vec::new(),
            marker_edit: false,
            overrides: HashMap::new(),
            grabbed: None,
            error,
        }
    }

    /// The primary (focused) entry index: drives the info line, viewport
    /// anchoring, and raw-channel display.
    pub fn primary(&self) -> usize {
        *self.selection.last().unwrap_or(&0)
    }

    /// True if entry `i` is in the current selection.
    pub fn is_selected(&self, i: usize) -> bool {
        self.selection.contains(&i)
    }

    /// Per-entry selection flags for the list highlight.
    pub fn selection_flags(&self) -> Vec<bool> {
        (0..self.entry_count())
            .map(|i| self.is_selected(i))
            .collect()
    }

    /// Apply a list click with modifier keys, updating the selection set.
    /// Returns true if the selection changed (caller resets the viewport).
    /// Ctrl toggles one entry (keeping at least one); Shift selects the range
    /// from the anchor; a plain click selects only that entry. Raw mode is
    /// always single-select (overlay is meaningful only for samples).
    pub fn select_click(&mut self, idx: usize, ctrl: bool, shift: bool) -> bool {
        if idx >= self.entry_count() {
            return false;
        }
        let before = self.selection.clone();
        if self.raw_mode() {
            self.selection = vec![idx];
            self.anchor = idx;
        } else if ctrl {
            if let Some(pos) = self.selection.iter().position(|&x| x == idx) {
                if self.selection.len() > 1 {
                    self.selection.remove(pos);
                }
            } else {
                self.selection.push(idx);
            }
            self.anchor = idx;
        } else if shift {
            let (lo, hi) = if self.anchor <= idx {
                (self.anchor, idx)
            } else {
                (idx, self.anchor)
            };
            self.selection = (lo..=hi).collect();
            // Keep `idx` as primary so the info line follows the click.
            if self.primary() != idx {
                self.selection.retain(|&x| x != idx);
                self.selection.push(idx);
            }
        } else {
            self.selection = vec![idx];
            self.anchor = idx;
        }
        let changed = self.selection != before;
        if changed {
            self.highlight_x = None; // stale peak highlight from the old sample
        }
        changed
    }

    /// True when there is no processed per-well signal but raw detector channels
    /// are available (a native `.xad`): the UI lists/plots the raw channels.
    pub fn raw_mode(&self) -> bool {
        let no_processed = self.run.samples.iter().all(|s| s.fluorescence.is_empty());
        no_processed && !self.raw_channels.is_empty()
    }

    /// Header summary line for the window.
    pub fn title(&self) -> String {
        let mode = if self.raw_mode() {
            "  [raw acquisition]"
        } else {
            ""
        };
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
