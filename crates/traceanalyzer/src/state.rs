//! Application state shared across Slint callbacks.
//!
//! The app is multi-file: [`AppState`] holds a `Vec<OpenFile>` and an `active`
//! index. Everything per-file (the run, its raw channels, selection, viewport,
//! marker overrides, dirty flag, tree expand state) lives on [`OpenFile`]; view
//! preferences that apply across files (y-mode, normalization, overview toggles)
//! stay on `AppState`. A single open file is just a list of length one — there
//! is no separate single-file path. `active == None` means no file is loaded.

use std::collections::HashMap;
use std::path::PathBuf;

use traceio::calibration::MarkerOverride;
use traceio::xad::RawChannel;
use traceio::Electrophoresis;

use crate::plot::{Viewport, YMode};

/// Parallel columns describing the visible rows of the multi-file well tree.
/// Rows span every open file: each file contributes a file-node row followed by
/// its well rows (omitted when that file is collapsed). All vectors share length.
pub struct TreeRows {
    pub labels: Vec<String>,
    pub is_file: Vec<bool>,
    /// Index of the file each row belongs to.
    pub file_index: Vec<i32>,
    /// Per-row: the owning file's expanded state (meaningful on file rows).
    pub file_expanded: Vec<bool>,
    /// Per-row: full file path (file rows only; empty on well rows) for tooltips.
    pub file_path: Vec<String>,
    /// Selection highlight (the active file's selected wells only).
    pub selected: Vec<bool>,
    /// Visible-row index of the active file's primary well, or `-1`.
    pub primary_row: i32,
}

/// Last path component of an instrument file path. The path may use Windows
/// (`\`) or Unix (`/`) separators, so split on both; falls back to the whole
/// string when there is no separator.
fn file_base_name(path: &str) -> String {
    path.rsplit(['\\', '/'])
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

/// One open file and all of its per-file view state.
pub struct OpenFile {
    pub run: Electrophoresis,
    /// Raw detector channels (populated for native `.xad`; see `raw_mode`).
    pub raw_channels: Vec<RawChannel>,
    /// Path the run was loaded from (target for File → Save).
    pub source_path: Option<PathBuf>,
    /// True when the run has in-memory edits (e.g. renamed wells) not yet saved.
    pub dirty: bool,
    /// Whether this file's node in the well tree is expanded.
    pub expanded: bool,
    /// Selected entry indices, in click order; the last is the primary (focused)
    /// entry. In raw mode this holds a single channel index.
    pub selection: Vec<usize>,
    /// Anchor for shift-range selection (index of the last plain/ctrl click).
    pub anchor: usize,
    /// Zoom/pan window; `None` means auto-fit the current entry.
    pub viewport: Option<Viewport>,
    /// Manual marker overrides, per sample index (empty = fully automatic).
    pub overrides: HashMap<usize, MarkerOverride>,
}

impl OpenFile {
    pub fn new(
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
    ) -> Self {
        OpenFile {
            run,
            raw_channels,
            source_path,
            dirty: false,
            expanded: true,
            selection: vec![0],
            anchor: 0,
            viewport: None,
            overrides: HashMap::new(),
        }
    }

    /// True when there is no processed per-well signal but raw detector channels
    /// are available (a native `.xad`): the UI lists/plots the raw channels.
    fn raw_mode(&self) -> bool {
        let no_processed = self.run.samples.iter().all(|s| s.fluorescence.is_empty());
        no_processed && !self.raw_channels.is_empty()
    }

    /// True for native Fragment Analyzer runs. FA saves are in-place `.txt`
    /// sidecar patches next to the immutable `.raw` acquisition file.
    fn fragment_analyzer_mode(&self) -> bool {
        self.run.assay.assay_name == "Fragment Analyzer"
            || self
                .source_path
                .as_ref()
                .and_then(|p| p.extension())
                .and_then(|e| e.to_str())
                .is_some_and(|e| e.eq_ignore_ascii_case("raw"))
    }

    /// Number of selectable entries in the current mode.
    fn entry_count(&self) -> usize {
        if self.raw_mode() {
            self.raw_channels.len()
        } else {
            self.run.samples.len()
        }
    }

    /// The primary (focused) entry index.
    fn primary(&self) -> usize {
        *self.selection.last().unwrap_or(&0)
    }

    /// Display labels for the entry list.
    fn entry_labels(&self) -> Vec<String> {
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

    /// File-node label: basename (or "run"), plus " (mod)" when unsaved.
    fn node_label(&self) -> String {
        let mut label = if self.run.assay.file_name.is_empty() {
            "run".to_string()
        } else {
            file_base_name(&self.run.assay.file_name)
        };
        if self.dirty {
            label.push_str(" (mod)");
        }
        label
    }
}

/// Which file (and, for a well row, which well) a visible tree row maps to.
#[derive(Clone, Copy)]
struct RowRef {
    file: usize,
    well: Option<usize>,
}

pub struct AppState {
    /// Open files. Empty ⇒ no file loaded (`active == None`).
    pub files: Vec<OpenFile>,
    /// Index of the active file, or `None` when none are open.
    pub active: Option<usize>,
    /// Peak-height-normalize overlaid traces (shape comparison).
    pub normalize: bool,
    /// Overview tab: share one y-scale across all mini-plots (vs. per-plot fit).
    pub overview_shared_y: bool,
    /// Overview tab: show the virtual gel instead of the small-multiples grid.
    pub overview_gel: bool,
    /// Overview tab: include ladder wells in the all-sample overview.
    pub overview_show_ladders: bool,
    pub y_mode: YMode,
    /// Marker-edit mode: plot in raw migration time with draggable marker lines.
    pub marker_edit: bool,
    /// x-position of a cross-highlighted peak (table→plot), if any.
    pub highlight_x: Option<f64>,
    /// Plot x-position for each current table row (`None` for region rows),
    /// parallel to the rows pushed to the table.
    pub table_peak_x: Vec<Option<f64>>,
    /// Marker currently being dragged, if any (transient).
    pub grabbed: Option<MarkerDrag>,
    /// Last load/save error, shown in the UI.
    pub error: Option<String>,
}

impl AppState {
    /// An empty state with no files loaded (the "No file loaded" screen). Files
    /// are added with [`AppState::add_file`].
    pub fn empty() -> Self {
        AppState {
            files: Vec::new(),
            active: None,
            normalize: false,
            // Per-plot y-scale by default: each overview cell auto-fits, so files
            // with very different signal magnitudes stay readable side by side.
            overview_shared_y: false,
            overview_gel: false,
            overview_show_ladders: false,
            y_mode: YMode::Fluorescence,
            marker_edit: false,
            highlight_x: None,
            table_peak_x: Vec::new(),
            grabbed: None,
            error: None,
        }
    }

    pub fn new(
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
        error: Option<String>,
    ) -> Self {
        let mut s = Self::empty();
        s.files.push(OpenFile::new(run, raw_channels, source_path));
        s.active = Some(0);
        s.error = error;
        s
    }

    // ---- active file access ------------------------------------------------

    pub fn active_file(&self) -> Option<&OpenFile> {
        self.active.and_then(|i| self.files.get(i))
    }

    pub fn active_file_mut(&mut self) -> Option<&mut OpenFile> {
        match self.active {
            Some(i) => self.files.get_mut(i),
            None => None,
        }
    }

    /// Add a file and make it active.
    pub fn add_file(
        &mut self,
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
    ) {
        self.files
            .push(OpenFile::new(run, raw_channels, source_path));
        self.active = Some(self.files.len() - 1);
        self.reset_transient();
    }

    /// Index of an already-open file whose source matches `path` (compared
    /// canonically), if any — so re-opening a run just re-activates it.
    pub fn find_file_by_source(&self, path: &std::path::Path) -> Option<usize> {
        let target = canonical(path);
        self.files
            .iter()
            .position(|f| f.source_path.as_deref().map(canonical).as_ref() == Some(&target))
    }

    /// Make `idx` the active file (no-op if out of range), preserving that
    /// file's own selection/viewport.
    pub fn activate_file(&mut self, idx: usize) {
        if idx < self.files.len() {
            self.active = Some(idx);
            self.reset_transient();
        }
    }

    /// Remove a file, fixing up the active index (⇒ `None` when the list empties).
    pub fn close_file(&mut self, idx: usize) {
        if idx >= self.files.len() {
            return;
        }
        self.files.remove(idx);
        self.active = if self.files.is_empty() {
            None
        } else {
            let a = self.active.unwrap_or(0);
            let a = if a > idx { a - 1 } else { a };
            Some(a.min(self.files.len() - 1))
        };
        self.reset_transient();
    }

    /// Clear transient view state that must not survive a file switch/close.
    fn reset_transient(&mut self) {
        self.highlight_x = None;
        self.grabbed = None;
        self.table_peak_x.clear();
    }

    // ---- pervasive accessors delegating to the active file -----------------
    // These unwrap and so must only be called when a file is active; main.rs
    // guards its entry points (empty-state renders before touching these).

    pub fn run(&self) -> &Electrophoresis {
        &self.active_file().expect("no active file").run
    }
    pub fn set_run(&mut self, run: Electrophoresis) {
        if let Some(f) = self.active_file_mut() {
            f.run = run;
        }
    }
    pub fn raw_channels(&self) -> &[RawChannel] {
        self.active_file()
            .map(|f| f.raw_channels.as_slice())
            .unwrap_or(&[])
    }
    pub fn selection(&self) -> &[usize] {
        self.active_file()
            .map(|f| f.selection.as_slice())
            .unwrap_or(&[])
    }
    pub fn viewport(&self) -> Option<Viewport> {
        self.active_file().and_then(|f| f.viewport)
    }
    pub fn set_viewport(&mut self, vp: Option<Viewport>) {
        if let Some(f) = self.active_file_mut() {
            f.viewport = vp;
        }
    }
    pub fn overrides(&self) -> &HashMap<usize, MarkerOverride> {
        &self.active_file().expect("no active file").overrides
    }
    pub fn overrides_mut(&mut self) -> &mut HashMap<usize, MarkerOverride> {
        &mut self.active_file_mut().expect("no active file").overrides
    }
    pub fn source_path(&self) -> Option<PathBuf> {
        self.active_file().and_then(|f| f.source_path.clone())
    }
    pub fn set_source_path(&mut self, p: Option<PathBuf>) {
        if let Some(f) = self.active_file_mut() {
            f.source_path = p;
        }
    }
    pub fn is_dirty(&self) -> bool {
        self.active_file().is_some_and(|f| f.dirty)
    }
    pub fn set_dirty(&mut self, dirty: bool) {
        if let Some(f) = self.active_file_mut() {
            f.dirty = dirty;
        }
    }

    /// Full path of the active file (may be an absolute Windows path); empty if
    /// none. Shown as the tree file-node tooltip.
    pub fn file_path(&self) -> String {
        self.active_file()
            .map(|f| f.run.assay.file_name.clone())
            .unwrap_or_default()
    }

    // ---- tree layout -------------------------------------------------------

    /// Flat list of visible rows across all files, in display order. The single
    /// source of truth for both `tree_rows` and row→(file, well) resolution.
    fn row_layout(&self) -> Vec<RowRef> {
        let mut rows = Vec::new();
        for (fi, f) in self.files.iter().enumerate() {
            rows.push(RowRef {
                file: fi,
                well: None,
            });
            if f.expanded {
                for w in 0..f.entry_count() {
                    rows.push(RowRef {
                        file: fi,
                        well: Some(w),
                    });
                }
            }
        }
        rows
    }

    pub fn tree_rows(&self) -> TreeRows {
        let mut t = TreeRows {
            labels: Vec::new(),
            is_file: Vec::new(),
            file_index: Vec::new(),
            file_expanded: Vec::new(),
            file_path: Vec::new(),
            selected: Vec::new(),
            primary_row: -1,
        };
        for (fi, f) in self.files.iter().enumerate() {
            t.labels.push(f.node_label());
            t.is_file.push(true);
            t.file_index.push(fi as i32);
            t.file_expanded.push(f.expanded);
            t.file_path.push(f.run.assay.file_name.clone());
            t.selected.push(false);
            if !f.expanded {
                continue;
            }
            let is_active = self.active == Some(fi);
            for (w, label) in f.entry_labels().into_iter().enumerate() {
                t.labels.push(label);
                t.is_file.push(false);
                t.file_index.push(fi as i32);
                t.file_expanded.push(false);
                t.file_path.push(String::new());
                t.selected.push(is_active && f.selection.contains(&w));
                if is_active && w == f.primary() {
                    t.primary_row = (t.labels.len() - 1) as i32;
                }
            }
        }
        t
    }

    fn row_ref(&self, row: usize) -> Option<RowRef> {
        self.row_layout().get(row).copied()
    }

    /// The file index a visible row belongs to (for closing/expanding).
    pub fn row_file(&self, row: usize) -> Option<usize> {
        self.row_ref(row).map(|r| r.file)
    }

    /// Click a visible tree row: activate its file and, for a well row, apply the
    /// selection click (resetting that file's viewport when the selection moves).
    /// Returns whether anything changed.
    pub fn select_row(&mut self, row: usize, ctrl: bool, shift: bool) -> bool {
        let Some(rr) = self.row_ref(row) else {
            return false;
        };
        let switched = self.active != Some(rr.file);
        if switched {
            self.active = Some(rr.file);
            self.reset_transient();
        }
        match rr.well {
            Some(w) => {
                let changed = self.select_click(w, ctrl, shift);
                if changed {
                    if let Some(f) = self.active_file_mut() {
                        f.viewport = None;
                    }
                }
                switched || changed
            }
            None => switched,
        }
    }

    /// Make `file_idx` active and plain-select well `well_idx` in it. Used by an
    /// Overview click, which can land in any open file (not just the active one).
    pub fn activate_and_select(&mut self, file_idx: usize, well_idx: usize) {
        if file_idx >= self.files.len() {
            return;
        }
        if self.active != Some(file_idx) {
            self.active = Some(file_idx);
            self.reset_transient();
        }
        self.select_click(well_idx, false, false);
        if let Some(f) = self.active_file_mut() {
            f.viewport = None;
        }
    }

    /// Toggle the expanded state of the file owning a visible row.
    pub fn toggle_expand_row(&mut self, row: usize) {
        if let Some(rr) = self.row_ref(row) {
            if let Some(f) = self.files.get_mut(rr.file) {
                f.expanded ^= true;
            }
        }
    }

    // ---- active-file queries (empty-safe) ----------------------------------

    /// Whether the primary-selected entry is a renameable well.
    pub fn can_rename(&self) -> bool {
        self.active_file()
            .is_some_and(|f| !f.raw_mode() && f.primary() < f.run.samples.len())
    }

    /// Current name of the primary-selected well (empty in raw mode / no file).
    pub fn primary_name(&self) -> String {
        self.active_file()
            .map(|f| {
                if f.raw_mode() {
                    String::new()
                } else {
                    f.run
                        .samples
                        .get(f.primary())
                        .map(|s| s.name.clone())
                        .unwrap_or_default()
                }
            })
            .unwrap_or_default()
    }

    /// Rename the active file's primary-selected well; marks that file dirty.
    /// Returns true if the name actually changed.
    pub fn rename_primary(&mut self, name: &str) -> bool {
        let Some(f) = self.active_file_mut() else {
            return false;
        };
        if f.raw_mode() {
            return false;
        }
        let idx = f.primary();
        if let Some(s) = f.run.samples.get_mut(idx) {
            let name = name.trim();
            if !name.is_empty() && s.name != name {
                s.name = name.to_string();
                f.dirty = true;
                return true;
            }
        }
        false
    }

    /// Whether the active file has a source path that File → Save can write to.
    pub fn can_save(&self) -> bool {
        self.active_file().is_some_and(|f| f.source_path.is_some())
    }

    /// True when the active file is a native Fragment Analyzer run.
    pub fn fragment_analyzer_mode(&self) -> bool {
        self.active_file()
            .is_some_and(|f| f.fragment_analyzer_mode())
    }

    /// The active file's primary (focused) entry index; 0 when no file is open.
    pub fn primary(&self) -> usize {
        self.active_file().map_or(0, |f| f.primary())
    }

    /// True if entry `i` is in the active file's selection.
    pub fn is_selected(&self, i: usize) -> bool {
        self.active_file().is_some_and(|f| f.selection.contains(&i))
    }

    /// Per-entry selection flags for the list highlight.
    pub fn selection_flags(&self) -> Vec<bool> {
        (0..self.entry_count())
            .map(|i| self.is_selected(i))
            .collect()
    }

    /// Apply a list click with modifier keys to the active file, updating its
    /// selection set. Returns true if the selection changed (caller resets the
    /// viewport). Ctrl toggles one entry (keeping at least one); Shift selects
    /// the range from the anchor; a plain click selects only that entry. Raw
    /// mode is always single-select (overlay is meaningful only for samples).
    pub fn select_click(&mut self, idx: usize, ctrl: bool, shift: bool) -> bool {
        let raw = self.raw_mode();
        if idx >= self.entry_count() {
            return false;
        }
        let changed = {
            let Some(f) = self.active_file_mut() else {
                return false;
            };
            let before = f.selection.clone();
            if raw {
                f.selection = vec![idx];
                f.anchor = idx;
            } else if ctrl {
                if let Some(pos) = f.selection.iter().position(|&x| x == idx) {
                    if f.selection.len() > 1 {
                        f.selection.remove(pos);
                    }
                } else {
                    f.selection.push(idx);
                }
                f.anchor = idx;
            } else if shift {
                let (lo, hi) = if f.anchor <= idx {
                    (f.anchor, idx)
                } else {
                    (idx, f.anchor)
                };
                f.selection = (lo..=hi).collect();
                // Keep `idx` as primary so the info line follows the click.
                if f.selection.last() != Some(&idx) {
                    f.selection.retain(|&x| x != idx);
                    f.selection.push(idx);
                }
            } else {
                f.selection = vec![idx];
                f.anchor = idx;
            }
            f.selection != before
        };
        if changed {
            self.highlight_x = None; // stale peak highlight from the old sample
        }
        changed
    }

    /// True when the active file is in raw-channel mode.
    pub fn raw_mode(&self) -> bool {
        self.active_file().is_some_and(|f| f.raw_mode())
    }

    /// Header summary line for the info area.
    pub fn title(&self) -> String {
        match self.active_file() {
            None => "No file loaded".to_string(),
            Some(f) => {
                let mode = if f.raw_mode() {
                    "  [raw acquisition]"
                } else {
                    ""
                };
                format!(
                    "{}  —  {} ({}),  {} samples{mode}",
                    f.run.assay.file_name,
                    f.run.assay.assay_name,
                    f.run.assay.assay_type,
                    f.run.samples.len()
                )
            }
        }
    }

    /// Number of selectable entries in the active file's current mode.
    pub fn entry_count(&self) -> usize {
        self.active_file().map_or(0, |f| f.entry_count())
    }

    /// Display labels for the active file's entry list.
    pub fn entry_labels(&self) -> Vec<String> {
        self.active_file()
            .map(|f| f.entry_labels())
            .unwrap_or_default()
    }
}

/// Canonicalize a path for identity comparison, falling back to the path as-is
/// when it cannot be resolved (e.g. it no longer exists).
fn canonical(path: &std::path::Path) -> PathBuf {
    std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf())
}

#[cfg(test)]
mod tests {
    use super::*;
    use traceio::{AssayInfo, Sample};

    fn sample(well: i32, name: &str) -> Sample {
        Sample {
            well_number: well,
            name: name.to_string(),
            category: String::new(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: Vec::new(),
            fluorescence: Vec::new(),
            aligned_time: Vec::new(),
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
        }
    }

    fn run_named(name: &str, wells: i32) -> Electrophoresis {
        Electrophoresis {
            assay: AssayInfo {
                file_name: name.to_string(),
                ..Default::default()
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: (0..wells).map(|w| sample(w + 1, "")).collect(),
        }
    }

    fn state_with_file(name: &str) -> AppState {
        AppState::new(run_named(name, 0), Vec::new(), None, None)
    }

    fn fa_state_with_source(source: &str) -> AppState {
        let mut run = run_named(source, 1);
        run.assay.assay_name = "Fragment Analyzer".to_string();
        AppState::new(run, Vec::new(), Some(PathBuf::from(source)), None)
    }

    #[test]
    fn tree_file_label_shows_basename_and_mod_marker_when_dirty() {
        let mut st = state_with_file(r"C:\Program Files\demo\Demo DNA 1000.xad");
        assert_eq!(st.tree_rows().labels[0], "Demo DNA 1000.xad");
        st.active_file_mut().unwrap().dirty = true;
        assert_eq!(st.tree_rows().labels[0], "Demo DNA 1000.xad (mod)");
    }

    #[test]
    fn tree_rows_span_multiple_files_with_per_file_mod() {
        let mut st = state_with_file("a.xml");
        st.add_file(run_named("b.xml", 0), Vec::new(), None);
        let t = st.tree_rows();
        assert_eq!(t.labels, vec!["a.xml", "b.xml"]);
        assert_eq!(t.is_file, vec![true, true]);
        assert_eq!(t.file_index, vec![0, 1]);
        // (mod) is per file: dirtying the second must not touch the first.
        st.files[1].dirty = true;
        let t = st.tree_rows();
        assert_eq!(t.labels[0], "a.xml");
        assert_eq!(t.labels[1], "b.xml (mod)");
    }

    #[test]
    fn row_mapping_across_files_with_mixed_expand() {
        // File 0: 2 wells (expanded); File 1: 1 well (expanded).
        let mut st = AppState::new(run_named("a", 2), Vec::new(), None, None);
        st.add_file(run_named("b", 1), Vec::new(), None);
        // Rows: [f0, a-w0, a-w1, f1, b-w0]
        assert_eq!(st.tree_rows().labels.len(), 5);
        assert_eq!(st.row_file(0), Some(0));
        assert_eq!(st.row_file(2), Some(0));
        assert_eq!(st.row_file(3), Some(1));
        assert_eq!(st.row_file(4), Some(1));
        // Clicking file 1's well activates file 1 and selects that well.
        st.select_row(4, false, false);
        assert_eq!(st.active, Some(1));
        assert_eq!(st.selection(), &[0]);
        // Collapsing file 0 shifts file 1 up: rows [f0, f1, b-w0].
        st.files[0].expanded = false;
        assert_eq!(st.tree_rows().labels.len(), 3);
        assert_eq!(st.row_file(1), Some(1)); // file 1's node
        assert_eq!(st.row_file(2), Some(1)); // file 1's well
    }

    #[test]
    fn close_file_fixes_up_active_index() {
        let mut st = state_with_file("a");
        st.add_file(run_named("b", 0), Vec::new(), None); // active = 1
        st.add_file(run_named("c", 0), Vec::new(), None); // active = 2
        assert_eq!(st.active, Some(2));
        st.close_file(2); // removed the active last file → clamp to 1
        assert_eq!(st.active, Some(1));
        st.close_file(0); // removing below active shifts it down → 0
        assert_eq!(st.active, Some(0));
        st.close_file(0); // last one gone → no active file
        assert_eq!(st.active, None);
        assert!(st.files.is_empty());
    }

    #[test]
    fn bioanalyzer_wells_can_be_renamed_and_saved() {
        let mut st = AppState::new(
            run_named("run.xml", 1),
            Vec::new(),
            Some(PathBuf::from("run.xml")),
            None,
        );

        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
        assert_eq!(st.primary_name(), "A1");
    }

    #[test]
    fn fragment_analyzer_runs_can_be_renamed_and_saved_in_place() {
        let mut st = fa_state_with_source("run.raw");

        assert!(st.fragment_analyzer_mode());
        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
        assert_eq!(st.run().samples[0].name, "A1");
    }

    #[test]
    fn raw_extension_source_is_treated_as_fragment_analyzer_saveable() {
        let mut st = AppState::new(
            run_named("run.raw", 1),
            Vec::new(),
            Some(PathBuf::from("run.RAW")),
            None,
        );

        assert!(st.fragment_analyzer_mode());
        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
    }

    #[test]
    fn find_file_by_source_dedups_and_activates() {
        let mut st = AppState::empty();
        st.add_file(run_named("a", 1), Vec::new(), Some(PathBuf::from("/runs/a.raw")));
        st.add_file(run_named("b", 1), Vec::new(), Some(PathBuf::from("/runs/b.raw")));

        assert_eq!(st.find_file_by_source(std::path::Path::new("/runs/a.raw")), Some(0));
        assert_eq!(st.find_file_by_source(std::path::Path::new("/runs/b.raw")), Some(1));
        assert_eq!(st.find_file_by_source(std::path::Path::new("/runs/c.raw")), None);

        st.activate_file(0);
        assert_eq!(st.active, Some(0));
    }
}
