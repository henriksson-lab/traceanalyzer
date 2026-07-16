//! Application state shared across Slint callbacks.

use std::collections::HashMap;
use std::path::PathBuf;

use traceio::calibration::MarkerOverride;
use traceio::xad::RawChannel;
use traceio::Electrophoresis;

use crate::plot::{Viewport, YMode};

/// Parallel columns describing the visible rows of the well tree: row 0 is the
/// file, rows 1.. are its wells (omitted when the file node is collapsed). All
/// vectors have the same length.
pub struct TreeRows {
    pub labels: Vec<String>,
    pub is_file: Vec<bool>,
    /// Sample index for well rows; `-1` for the file row.
    pub well_index: Vec<i32>,
    pub selected: Vec<bool>,
    /// Visible-row index of the primary well, or `-1`.
    pub primary_row: i32,
}

/// Last path component of an instrument file path. The path may use Windows
/// (`\`) or Unix (`/`) separators, so split on both; falls back to the whole
/// string when there is no separator.
fn file_base_name(path: &str) -> String {
    path.rsplit(|c| c == '\\' || c == '/')
        .next()
        .unwrap_or(path)
        .to_string()
}

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
    /// Overview tab: include ladder wells in the all-sample overview.
    pub overview_show_ladders: bool,
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
    /// Whether the file node in the well tree is expanded.
    pub expanded: bool,
    /// Path the run was loaded from (target for File → Save).
    pub source_path: Option<PathBuf>,
    /// True when the run has in-memory edits (e.g. renamed wells) not yet saved.
    pub dirty: bool,
    /// Last load error, shown in the UI.
    pub error: Option<String>,
}

impl AppState {
    pub fn new(
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
        error: Option<String>,
    ) -> Self {
        AppState {
            run,
            raw_channels,
            selection: vec![0],
            anchor: 0,
            normalize: false,
            overview_shared_y: true,
            overview_gel: false,
            overview_show_ladders: false,
            y_mode: YMode::Fluorescence,
            viewport: None,
            highlight_x: None,
            table_peak_x: Vec::new(),
            marker_edit: false,
            overrides: HashMap::new(),
            grabbed: None,
            expanded: true,
            source_path,
            dirty: false,
            error,
        }
    }

    /// Visible rows of the well tree (see [`TreeRows`]).
    /// Full path of the loaded file as reported by the instrument (may be an
    /// absolute Windows path); empty if none. Shown as the tree file-node tooltip.
    pub fn file_path(&self) -> String {
        self.run.assay.file_name.clone()
    }

    pub fn tree_rows(&self) -> TreeRows {
        let mut file_label = if self.run.assay.file_name.is_empty() {
            "run".to_string()
        } else {
            file_base_name(&self.run.assay.file_name)
        };
        // Flag unsaved edits on the file node only (not the wells below it).
        if self.dirty {
            file_label.push_str(" (mod)");
        }
        let mut t = TreeRows {
            labels: vec![file_label],
            is_file: vec![true],
            well_index: vec![-1],
            selected: vec![false],
            primary_row: -1,
        };
        if self.expanded {
            for (i, label) in self.entry_labels().into_iter().enumerate() {
                t.labels.push(label);
                t.is_file.push(false);
                t.well_index.push(i as i32);
                t.selected.push(self.is_selected(i));
                if i == self.primary() {
                    t.primary_row = (t.labels.len() - 1) as i32;
                }
            }
        }
        t
    }

    /// Whether the primary-selected entry is a renameable well (not a raw channel).
    pub fn can_rename(&self) -> bool {
        !self.raw_mode() && self.primary() < self.run.samples.len()
    }

    /// Current name of the primary-selected well (empty in raw mode).
    pub fn primary_name(&self) -> String {
        if self.raw_mode() {
            return String::new();
        }
        self.run
            .samples
            .get(self.primary())
            .map(|s| s.name.clone())
            .unwrap_or_default()
    }

    /// Rename the primary-selected well in memory; marks the run dirty.
    /// Returns true if the name actually changed.
    pub fn rename_primary(&mut self, name: &str) -> bool {
        if self.raw_mode() {
            return false;
        }
        let idx = self.primary();
        if let Some(s) = self.run.samples.get_mut(idx) {
            let name = name.trim();
            if !name.is_empty() && s.name != name {
                s.name = name.to_string();
                self.dirty = true;
                return true;
            }
        }
        false
    }

    /// Whether there is a loaded file that File → Save can write to.
    pub fn can_save(&self) -> bool {
        self.source_path.is_some()
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

#[cfg(test)]
mod tests {
    use super::*;
    use traceio::AssayInfo;

    fn state_with_file(name: &str) -> AppState {
        let run = Electrophoresis {
            assay: AssayInfo {
                file_name: name.to_string(),
                ..Default::default()
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: Vec::new(),
        };
        AppState::new(run, Vec::new(), None, None)
    }

    #[test]
    fn tree_file_label_shows_basename_and_mod_marker_when_dirty() {
        let mut st = state_with_file(r"C:\Program Files\demo\Demo DNA 1000.xad");
        // Clean: basename only, no marker.
        assert_eq!(st.tree_rows().labels[0], "Demo DNA 1000.xad");
        // Dirty: same basename plus the " (mod)" marker.
        st.dirty = true;
        assert_eq!(st.tree_rows().labels[0], "Demo DNA 1000.xad (mod)");
    }
}
