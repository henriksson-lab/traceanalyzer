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

use traceio::calibration::{marker_times, MarkerOverride};
use traceio::io::{DetectedFormat, SaveCapabilities, SourceCapabilities};
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
    path.rsplit(['\\', '/']).next().unwrap_or(path).to_string()
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
    /// Optional human-readable sidecar metadata/QC summary for native formats
    /// that expose it outside the normalized trace model.
    pub metadata_text: Option<String>,
    /// Compact source-specific key/value summary shown in the Metadata tab.
    pub metadata_rows: Vec<(String, String)>,
    /// Path the run was loaded from (target for File → Save).
    pub source_path: Option<PathBuf>,
    /// User-facing path that was opened or last saved. For multi-file formats
    /// this can differ from the canonical save/dedup identity in `source_path`.
    pub opened_path: Option<PathBuf>,
    /// Source capabilities from `traceio` detection, when this file came through
    /// the path-oriented loader.
    pub source_capabilities: Option<SourceCapabilities>,
    /// Save capabilities from `traceio` detection, when available.
    pub save_capabilities: Option<SaveCapabilities>,
    /// True when the run has saveable in-memory edits (e.g. renamed wells).
    /// Marker overrides are session-only calibration/view state and must not
    /// make Save imply persistence until the save layer supports them.
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
            metadata_text: None,
            metadata_rows: Vec::new(),
            opened_path: source_path.clone(),
            source_path,
            source_capabilities: None,
            save_capabilities: None,
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

    fn bioanalyzer_xad_mode(&self) -> bool {
        self.source_path
            .as_ref()
            .and_then(|p| p.extension())
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("xad"))
    }

    /// True for TapeStation exports. They are intentionally read-only: the XML
    /// is a derived export, and native project files are encrypted.
    fn tapestation_mode(&self) -> bool {
        self.source_path
            .as_ref()
            .and_then(|p| p.file_name())
            .map(|n| n.to_string_lossy().to_ascii_lowercase())
            .is_some_and(|n| {
                n.ends_with("_electropherogram.csv") || n.ends_with("_electropherogram.csv.gz")
            })
            || is_tapestation_assay_name(&self.run.assay.assay_name)
            || self
                .run
                .assay
                .file_name
                .rsplit_once('.')
                .map(|(_, ext)| {
                    matches!(
                        ext.to_ascii_lowercase().as_str(),
                        "d1000"
                            | "hsd1000"
                            | "d5000"
                            | "hsd5000"
                            | "cfdna"
                            | "gdna"
                            | "rna"
                            | "hsrna"
                    )
                })
                .unwrap_or(false)
            || self.run.samples.iter().any(|s| !s.regions.is_empty())
    }

    fn has_processed_sample_data(&self) -> bool {
        self.run.samples.iter().any(sample_has_processed_data)
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
        let source_label = self
            .opened_path
            .as_ref()
            .or(self.source_path.as_ref())
            .and_then(|path| path.file_name())
            .map(|name| name.to_string_lossy().into_owned());
        let mut label = if let Some(name) = source_label {
            name
        } else if self.run.assay.file_name.is_empty() {
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
    /// Last non-error user-visible status message, shown until an error replaces
    /// it or a new action updates it.
    pub status: Option<String>,
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
            status: None,
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

    pub fn add_file_with_capabilities(
        &mut self,
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
        source_capabilities: SourceCapabilities,
        save_capabilities: SaveCapabilities,
    ) {
        self.files
            .push(OpenFile::new(run, raw_channels, source_path));
        if let Some(file) = self.files.last_mut() {
            file.source_capabilities = Some(source_capabilities);
            file.save_capabilities = Some(save_capabilities);
        }
        self.active = Some(self.files.len() - 1);
        self.reset_transient();
    }

    pub fn replace_active_file_with_capabilities(
        &mut self,
        run: Electrophoresis,
        raw_channels: Vec<RawChannel>,
        source_path: Option<PathBuf>,
        source_capabilities: SourceCapabilities,
        save_capabilities: SaveCapabilities,
    ) {
        let Some(idx) = self.active else {
            return;
        };
        if idx >= self.files.len() {
            return;
        }
        let old = &self.files[idx];
        let expanded = old.expanded;
        let opened_path = old.opened_path.clone();
        let old_selection = old.selection.clone();
        let old_anchor = old.anchor;
        let viewport = old.viewport;
        let overrides = old.overrides.clone();
        let mut file = OpenFile::new(run, raw_channels, source_path);
        let entry_count = file.entry_count();
        file.expanded = expanded;
        file.opened_path = opened_path;
        file.selection = old_selection
            .into_iter()
            .filter(|&idx| idx < entry_count)
            .collect();
        if file.selection.is_empty() && entry_count > 0 {
            file.selection.push(0);
        }
        file.anchor = old_anchor.min(entry_count.saturating_sub(1));
        file.viewport = viewport;
        let sample_count = file.run.samples.len();
        file.overrides = overrides
            .into_iter()
            .filter(|(idx, _)| *idx < sample_count)
            .collect();
        file.source_capabilities = Some(source_capabilities);
        file.save_capabilities = Some(save_capabilities);
        self.files[idx] = file;
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
            self.reset_no_active_file();
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

    /// Drop file-scoped view controls when there is no active file to apply
    /// them to. Overview toggles remain application preferences.
    fn reset_no_active_file(&mut self) {
        self.normalize = false;
        self.y_mode = YMode::Fluorescence;
        self.marker_edit = false;
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
    pub fn metadata_text(&self) -> &str {
        self.active_file()
            .and_then(|f| f.metadata_text.as_deref())
            .unwrap_or("")
    }
    pub fn set_active_metadata_text(&mut self, text: Option<String>) {
        if let Some(f) = self.active_file_mut() {
            f.metadata_text = text;
        }
    }
    pub fn metadata_rows(&self) -> &[(String, String)] {
        self.active_file()
            .map(|f| f.metadata_rows.as_slice())
            .unwrap_or(&[])
    }
    pub fn set_active_metadata_rows(&mut self, rows: Vec<(String, String)>) {
        if let Some(f) = self.active_file_mut() {
            f.metadata_rows = rows;
        }
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
    pub fn has_marker_overrides(&self) -> bool {
        self.active_file()
            .is_some_and(|f| file_has_effective_marker_overrides(f))
    }
    pub fn file_has_marker_overrides(&self, idx: usize) -> bool {
        self.files
            .get(idx)
            .is_some_and(file_has_effective_marker_overrides)
    }
    pub fn file_sample_has_marker_override(&self, file_idx: usize, sample_idx: usize) -> bool {
        self.files
            .get(file_idx)
            .is_some_and(|f| sample_has_effective_marker_override(f, sample_idx))
    }
    pub fn source_path(&self) -> Option<PathBuf> {
        self.active_file().and_then(|f| f.source_path.clone())
    }
    pub fn opened_path(&self) -> Option<PathBuf> {
        self.active_file().and_then(|f| f.opened_path.clone())
    }
    pub fn set_active_opened_path(&mut self, p: Option<PathBuf>) {
        if let Some(f) = self.active_file_mut() {
            f.opened_path = p;
        }
    }
    pub fn set_source_path(&mut self, p: Option<PathBuf>) {
        if let Some(f) = self.active_file_mut() {
            f.source_path = p;
        }
    }
    pub fn adopt_detected_source(&mut self, detected: DetectedFormat) {
        if let Some(f) = self.active_file_mut() {
            let save_capabilities = detected.save_capabilities();
            f.opened_path = Some(detected.path);
            f.source_path = Some(detected.identity);
            f.source_capabilities = Some(detected.capabilities);
            f.save_capabilities = Some(save_capabilities);
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

    /// True if a session-only marker override changed its effective marker
    /// positions since `before`. Empty overrides are equivalent to no override.
    pub fn marker_override_changed(&self, idx: usize, before: Option<MarkerOverride>) -> bool {
        let Some(f) = self.active_file() else {
            return false;
        };
        let before = effective_marker_times(&f.run, idx, before);
        let after = effective_marker_times(&f.run, idx, f.overrides.get(&idx).copied());
        before != after
    }

    /// Full path of the active file (may be an absolute Windows path); empty if
    /// none. Shown as the tree file-node tooltip.
    pub fn file_path(&self) -> String {
        self.active_file()
            .and_then(|f| {
                f.source_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .or_else(|| {
                        (!f.run.assay.file_name.is_empty()).then(|| f.run.assay.file_name.clone())
                    })
            })
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
            t.file_path.push(
                f.source_path
                    .as_ref()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|| f.run.assay.file_name.clone()),
            );
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
            f.expanded = true;
            f.viewport = None;
        }
    }

    /// Activate and toggle the expanded state of the file owning a visible row.
    /// Returns whether the active file changed.
    pub fn toggle_expand_row(&mut self, row: usize) -> bool {
        let Some(rr) = self.row_ref(row) else {
            return false;
        };
        let switched = self.active != Some(rr.file);
        if switched {
            self.active = Some(rr.file);
            self.reset_transient();
        }
        if let Some(f) = self.files.get_mut(rr.file) {
            f.expanded ^= true;
        }
        switched
    }

    // ---- active-file queries (empty-safe) ----------------------------------

    /// Whether the primary-selected entry is a renameable well.
    pub fn can_rename(&self) -> bool {
        self.active_file().is_some_and(|f| {
            if f.raw_mode() || f.primary() >= f.run.samples.len() {
                return false;
            }
            if let Some(capabilities) = f.source_capabilities {
                return capabilities.can_rename;
            }
            !f.bioanalyzer_xad_mode()
                && !f.tapestation_mode()
                && (!f.fragment_analyzer_mode()
                    || f.source_path
                        .as_deref()
                        .is_some_and(traceio::fa::has_saveable_txt_sidecar))
        })
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
        if !self.can_rename() {
            return false;
        }
        let Some(f) = self.active_file_mut() else {
            return false;
        };
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

    /// Whether File → Save can persist the active file. Most formats save in
    /// place; native `.xad` sources use the Save As fallback because the source
    /// container is not rewritable.
    pub fn can_save(&self) -> bool {
        self.active_file().is_some_and(|f| {
            if f.source_path.is_none() {
                return false;
            }
            if let Some(capabilities) = f.save_capabilities {
                return capabilities.can_save_in_place
                    || (f.bioanalyzer_xad_mode()
                        && (capabilities.can_save_as_xml || capabilities.can_save_as_xml_gz));
            }
            !f.bioanalyzer_xad_mode()
                && !f.tapestation_mode()
                && (!f.fragment_analyzer_mode()
                    || f.source_path
                        .as_deref()
                        .is_some_and(traceio::fa::has_saveable_txt_sidecar))
        })
    }

    /// Whether the active file can be saved to a new Bioanalyzer XML file.
    pub fn can_save_as(&self) -> bool {
        if self.has_marker_overrides() {
            return false;
        }
        self.active_file().is_some_and(|f| {
            if f.source_path.is_none() {
                return false;
            }
            if let Some(capabilities) = f.save_capabilities {
                return capabilities.can_save_as_xml || capabilities.can_save_as_xml_gz;
            }
            if f.raw_mode() {
                return false;
            }
            !f.bioanalyzer_xad_mode() && !f.fragment_analyzer_mode() && !f.tapestation_mode()
        })
    }

    /// Whether marker override editing is meaningful for the active file.
    pub fn can_edit_markers(&self) -> bool {
        self.active_file().is_some_and(|f| {
            if f.raw_mode()
                || f.selection.len() != 1
                || !f.selection.iter().all(|&idx| idx < f.run.samples.len())
            {
                return false;
            }
            let source_allows = if let Some(capabilities) = f.source_capabilities {
                capabilities.can_edit_markers
            } else {
                !f.bioanalyzer_xad_mode() && !f.fragment_analyzer_mode() && !f.tapestation_mode()
            };
            source_allows && sample_has_usable_marker_times(f, f.selection[0])
        })
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
            if !self.can_edit_markers() {
                self.marker_edit = false;
                self.grabbed = None;
            }
        }
        changed
    }

    /// True when the active file is in raw-channel mode.
    pub fn raw_mode(&self) -> bool {
        self.active_file().is_some_and(|f| f.raw_mode())
    }

    /// True when the active run has at least one sample with processed trace
    /// arrays, as opposed to metadata-only sample shells.
    pub fn has_processed_sample_data(&self) -> bool {
        self.active_file()
            .is_some_and(|f| f.has_processed_sample_data())
    }

    pub fn sample_has_processed_data(&self, idx: usize) -> bool {
        self.active_file()
            .and_then(|f| f.run.samples.get(idx))
            .is_some_and(sample_has_processed_data)
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

fn is_tapestation_assay_name(name: &str) -> bool {
    let normalized: String = name
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_uppercase())
        .collect();
    matches!(
        normalized.as_str(),
        "D1000" | "HSD1000" | "D5000" | "HSD5000" | "CFDNA" | "GDNA" | "RNA" | "HSRNA"
    )
}

fn sample_has_processed_data(sample: &traceio::Sample) -> bool {
    !sample.fluorescence.is_empty()
        && (!sample.time.is_empty() || !sample.aligned_time.is_empty() || !sample.length.is_empty())
}

fn sample_has_usable_marker_times(file: &OpenFile, idx: usize) -> bool {
    let Some(sample) = file.run.samples.get(idx) else {
        return false;
    };
    if !sample_has_processed_data(sample) {
        return false;
    }
    let ov = file.overrides.get(&idx).filter(|ov| !ov.is_empty());
    let (lower, upper) = marker_times(&file.run, idx, ov);
    let Some(lower) = lower.filter(|time| time.is_finite()) else {
        return false;
    };
    if file.run.assay.has_upper_marker {
        let Some(upper) = upper.filter(|time| time.is_finite()) else {
            return false;
        };
        lower < upper
    } else {
        true
    }
}

fn effective_marker_times(
    run: &Electrophoresis,
    idx: usize,
    ov: Option<MarkerOverride>,
) -> (Option<f64>, Option<f64>) {
    let ov = ov.filter(|ov| !ov.is_empty());
    marker_times(run, idx, ov.as_ref())
}

fn file_has_effective_marker_overrides(file: &OpenFile) -> bool {
    file.overrides
        .iter()
        .any(|(&idx, _)| sample_has_effective_marker_override(file, idx))
}

fn sample_has_effective_marker_override(file: &OpenFile, idx: usize) -> bool {
    let Some(override_state) = file.overrides.get(&idx).copied() else {
        return false;
    };
    if override_state.is_empty() {
        return false;
    }
    effective_marker_times(&file.run, idx, None)
        != effective_marker_times(&file.run, idx, Some(override_state))
}

#[cfg(test)]
mod tests {
    use super::*;
    use traceio::io::TraceFormat;
    use traceio::{AssayInfo, LadderPeak, Peak, Sample};

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
            regions: Vec::new(),
        }
    }

    fn processed_sample(well: i32, name: &str) -> Sample {
        let mut sample = sample(well, name);
        sample.time = vec![0.0, 1.0, 9.0, 10.0];
        sample.fluorescence = vec![1.0, 2.0, 2.0, 1.0];
        sample
    }

    fn marker_peak(name: &str, time: f64, concentration: f64) -> Peak {
        Peak {
            observations: name.to_string(),
            length: f64::NAN,
            time,
            aligned_time: time,
            start_time: f64::NAN,
            end_time: f64::NAN,
            aligned_start_time: f64::NAN,
            aligned_end_time: f64::NAN,
            area: f64::NAN,
            concentration,
            molarity: f64::NAN,
        }
    }

    fn state_with_detected_markers() -> AppState {
        let mut run = processed_run_named("run.xml", 2);
        run.assay.has_upper_marker = true;
        run.ladder_peaks = vec![
            LadderPeak {
                size: 15.0,
                area_a: f64::NAN,
                area_b: f64::NAN,
                concentration: 1.0,
            },
            LadderPeak {
                size: 1500.0,
                area_a: f64::NAN,
                area_b: f64::NAN,
                concentration: 2.0,
            },
        ];
        for sample in &mut run.samples {
            sample.peaks = vec![
                marker_peak("Lower Marker", 1.0, 1.0),
                marker_peak("Upper Marker", 9.0, 2.0),
            ];
        }
        AppState::new(run, Vec::new(), Some(PathBuf::from("run.xml")), None)
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

    fn processed_run_named(name: &str, wells: i32) -> Electrophoresis {
        let mut run = run_named(name, 0);
        run.samples = (0..wells).map(|w| processed_sample(w + 1, "")).collect();
        run
    }

    fn state_with_file(name: &str) -> AppState {
        AppState::new(run_named(name, 0), Vec::new(), None, None)
    }

    fn fa_state_with_source(source: &str) -> AppState {
        let mut run = processed_run_named(source, 1);
        run.assay.assay_name = "Fragment Analyzer".to_string();
        AppState::new(run, Vec::new(), Some(PathBuf::from(source)), None)
    }

    fn tapestation_state_with_source(source: &str) -> AppState {
        let mut run = run_named("demo.D1000", 1);
        run.samples[0].regions.push(traceio::Region {
            lower_length: 200.0,
            upper_length: 400.0,
        });
        AppState::new(run, Vec::new(), Some(PathBuf::from(source)), None)
    }

    fn raw_channel() -> RawChannel {
        RawChannel {
            channel_id: "BlueFluorescence".to_string(),
            name: "Blue".to_string(),
            x_start: 0.0,
            x_step: 1.0,
            signal: vec![1.0, 2.0],
        }
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
    fn marker_override_change_is_session_only_and_not_save_dirty() {
        let mut st = state_with_detected_markers();
        let before = st.overrides().get(&0).copied();

        st.overrides_mut().insert(
            0,
            MarkerOverride {
                lower_time: Some(2.0),
                upper_time: None,
            },
        );

        assert!(st.marker_override_changed(0, before));
        assert!(!st.is_dirty());
    }

    #[test]
    fn unchanged_marker_override_does_not_mark_save_dirty() {
        let mut st = state_with_detected_markers();
        let ov = MarkerOverride {
            lower_time: Some(1.0),
            upper_time: None,
        };
        st.overrides_mut().insert(0, ov);

        assert!(!st.marker_override_changed(0, Some(ov)));
        assert!(!st.is_dirty());
    }

    #[test]
    fn empty_marker_override_is_equivalent_to_absent_for_change_tracking() {
        let mut st = state_with_detected_markers();
        let before = st.overrides().get(&0).copied();

        st.overrides_mut().insert(0, MarkerOverride::default());

        assert!(!st.marker_override_changed(0, before));
        assert!(!st.is_dirty());
    }

    #[test]
    fn marker_override_equal_to_detected_marker_is_not_changed() {
        let mut st = state_with_detected_markers();
        let before = st.overrides().get(&0).copied();

        st.overrides_mut().insert(
            0,
            MarkerOverride {
                lower_time: Some(1.0),
                upper_time: Some(9.0),
            },
        );

        assert!(!st.marker_override_changed(0, before));
        assert!(!st.has_marker_overrides());
        assert!(st.can_save_as());
        assert!(!st.is_dirty());
    }

    #[test]
    fn marker_edit_requires_single_selected_sample() {
        let mut st = state_with_detected_markers();
        assert!(st.can_edit_markers());

        assert!(st.select_click(1, true, false));
        assert!(!st.can_edit_markers());

        st.marker_edit = true;
        assert!(st.select_click(0, false, false));
        assert!(st.can_edit_markers());
        assert!(st.marker_edit);
    }

    #[test]
    fn marker_edit_requires_processed_selected_sample() {
        let mut st = state_with_detected_markers();
        st.active_file_mut().unwrap().run.samples[0]
            .fluorescence
            .clear();

        assert!(!st.can_edit_markers());

        st.active_file_mut().unwrap().run.samples[0].fluorescence = vec![1.0, 2.0, 2.0, 1.0];
        assert!(st.can_edit_markers());
    }

    #[test]
    fn marker_edit_requires_finite_lower_marker() {
        let mut st = state_with_detected_markers();
        st.active_file_mut().unwrap().run.samples[0].peaks[0].time = f64::NAN;

        assert!(!st.can_edit_markers());

        st.active_file_mut().unwrap().run.samples[0].peaks[0].time = 1.0;
        assert!(st.can_edit_markers());
    }

    #[test]
    fn marker_edit_requires_upper_marker_when_assay_has_upper_marker() {
        let mut st = state_with_detected_markers();
        st.active_file_mut().unwrap().run.samples[0].peaks.pop();

        assert!(!st.can_edit_markers());

        st.active_file_mut().unwrap().run.assay.has_upper_marker = false;
        assert!(st.can_edit_markers());
    }

    #[test]
    fn marker_edit_requires_lower_before_upper_marker() {
        let mut st = state_with_detected_markers();
        st.active_file_mut().unwrap().run.samples[0].peaks[1].time = 0.5;

        assert!(!st.can_edit_markers());
    }

    #[test]
    fn multi_selection_turns_off_marker_edit_mode() {
        let mut st = state_with_detected_markers();
        st.marker_edit = true;

        assert!(st.select_click(1, true, false));

        assert!(!st.can_edit_markers());
        assert!(!st.marker_edit);
        assert_eq!(st.grabbed, None);
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
    fn file_row_toggle_activates_that_file() {
        let mut st = AppState::new(run_named("a", 1), Vec::new(), None, None);
        st.add_file(run_named("b", 1), Vec::new(), None);
        st.active = Some(0);

        assert!(st.toggle_expand_row(2));

        assert_eq!(st.active, Some(1));
        assert!(!st.files[1].expanded);
    }

    #[test]
    fn overview_selection_expands_collapsed_file() {
        let mut st = AppState::new(run_named("a", 1), Vec::new(), None, None);
        st.add_file(run_named("b", 1), Vec::new(), None);
        st.files[0].expanded = false;
        st.active = Some(1);

        st.activate_and_select(0, 0);

        assert_eq!(st.active, Some(0));
        assert!(st.files[0].expanded);
        assert_eq!(st.selection(), &[0]);
        assert_eq!(st.tree_rows().primary_row, 1);
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
    fn close_inactive_file_preserves_active_logical_file() {
        let mut st = state_with_file("a");
        st.add_file(run_named("b", 0), Vec::new(), None);
        st.add_file(run_named("c", 0), Vec::new(), None);
        st.add_file(run_named("d", 0), Vec::new(), None);
        st.activate_file(2);
        assert_eq!(st.active_file().unwrap().run.assay.file_name, "c");

        st.close_file(0);
        assert_eq!(st.active, Some(1));
        assert_eq!(st.active_file().unwrap().run.assay.file_name, "c");

        st.close_file(2);
        assert_eq!(st.active, Some(1));
        assert_eq!(st.active_file().unwrap().run.assay.file_name, "c");
    }

    #[test]
    fn close_last_file_clears_file_scoped_view_state() {
        let mut st = state_with_file("a");
        st.normalize = true;
        st.y_mode = YMode::Concentration;
        st.marker_edit = true;
        st.highlight_x = Some(42.0);
        st.table_peak_x = vec![Some(42.0)];
        st.overview_shared_y = true;

        st.close_file(0);

        assert_eq!(st.active, None);
        assert!(!st.normalize);
        assert_eq!(st.y_mode, YMode::Fluorescence);
        assert!(!st.marker_edit);
        assert_eq!(st.highlight_x, None);
        assert!(st.table_peak_x.is_empty());
        assert!(st.overview_shared_y);
    }

    #[test]
    fn bioanalyzer_wells_can_be_renamed_and_saved() {
        let mut st = state_with_detected_markers();

        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(st.can_save_as());
        assert!(st.can_edit_markers());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
        assert_eq!(st.primary_name(), "A1");
    }

    #[test]
    fn marker_overrides_disable_save_as_without_dirtying_file() {
        let mut st = state_with_detected_markers();

        st.overrides_mut().insert(
            0,
            MarkerOverride {
                lower_time: Some(2.0),
                upper_time: None,
            },
        );

        assert!(st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.is_dirty());
    }

    #[test]
    fn xad_sources_are_read_only_even_with_processed_samples() {
        let mut st = AppState::new(
            processed_run_named("run.xad", 1),
            Vec::new(),
            Some(PathBuf::from("run.xad")),
            None,
        );

        assert!(!st.can_rename());
        assert!(!st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.can_edit_markers());
        assert!(!st.rename_primary("A1"));
        assert!(!st.is_dirty());
    }

    #[test]
    fn raw_xad_can_save_as_when_detected_capabilities_allow_it() {
        let detected = traceio::io::DetectedFormat {
            path: PathBuf::from("run.xad"),
            identity: PathBuf::from("run.xad"),
            format: TraceFormat::BioanalyzerXad,
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: false,
                can_save_as_xml: true,
                can_edit_markers: false,
            },
        };
        let mut st = AppState::empty();
        st.add_file_with_capabilities(
            run_named("run.xad", 1),
            vec![raw_channel()],
            Some(detected.identity.clone()),
            detected.capabilities,
            detected.save_capabilities(),
        );

        assert!(st.raw_mode());
        assert!(st.can_save_as());
    }

    #[test]
    fn detected_xad_can_use_save_menu_to_fall_back_to_save_as() {
        let detected = traceio::io::DetectedFormat {
            path: PathBuf::from("run.xad"),
            identity: PathBuf::from("run.xad"),
            format: TraceFormat::BioanalyzerXad,
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: false,
                can_save_as_xml: true,
                can_edit_markers: false,
            },
        };
        let mut st = AppState::empty();
        st.add_file_with_capabilities(
            processed_run_named("run.xad", 1),
            Vec::new(),
            Some(detected.identity.clone()),
            detected.capabilities,
            detected.save_capabilities(),
        );

        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(st.can_save_as());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
    }

    #[test]
    fn adopt_detected_source_updates_active_identity_and_capabilities() {
        let xad = traceio::io::DetectedFormat {
            path: PathBuf::from("run.xad"),
            identity: PathBuf::from("run.xad"),
            format: TraceFormat::BioanalyzerXad,
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: false,
                can_save_as_xml: true,
                can_edit_markers: false,
            },
        };
        let xml = traceio::io::DetectedFormat {
            path: PathBuf::from("saved.xml"),
            identity: PathBuf::from("/canonical/saved.xml"),
            format: TraceFormat::BioanalyzerXml,
            capabilities: SourceCapabilities {
                can_rename: true,
                can_save_in_place: true,
                can_save_as_xml: true,
                can_edit_markers: true,
            },
        };
        let mut st = AppState::empty();
        st.add_file_with_capabilities(
            processed_run_named("run.xad", 1),
            Vec::new(),
            Some(xad.identity.clone()),
            xad.capabilities,
            xad.save_capabilities(),
        );

        st.adopt_detected_source(xml);

        assert_eq!(
            st.source_path(),
            Some(PathBuf::from("/canonical/saved.xml"))
        );
        assert!(st.can_save());
        assert!(st.can_save_as());
        let file = st.active_file().unwrap();
        assert!(file.source_capabilities.unwrap().can_edit_markers);
        assert!(file.save_capabilities.unwrap().can_save_in_place);
    }

    #[test]
    fn replace_active_file_preserves_session_view_state() {
        let mut st = state_with_detected_markers();
        st.files[0].expanded = false;
        st.files[0].selection = vec![1];
        st.files[0].anchor = 1;
        st.files[0].viewport = Some(Viewport {
            x_min: 1.0,
            x_max: 9.0,
            y_min: 0.0,
            y_max: 10.0,
        });
        st.files[0].overrides.insert(
            1,
            MarkerOverride {
                lower_time: Some(2.0),
                upper_time: Some(8.0),
            },
        );

        st.replace_active_file_with_capabilities(
            processed_run_named("saved.xml", 2),
            Vec::new(),
            Some(PathBuf::from("saved.xml")),
            SourceCapabilities {
                can_rename: true,
                can_save_in_place: true,
                can_save_as_xml: true,
                can_edit_markers: true,
            },
            SaveCapabilities {
                can_save_sample_names: true,
                can_save_in_place: true,
                can_save_as_xml: true,
                can_save_as_xml_gz: true,
            },
        );

        let file = st.active_file().unwrap();
        assert!(!file.expanded);
        assert_eq!(file.selection, vec![1]);
        assert_eq!(file.anchor, 1);
        assert_eq!(file.viewport.unwrap().x_min, 1.0);
        assert_eq!(file.overrides.get(&1).unwrap().lower_time, Some(2.0));
        assert!(!file.dirty);
    }

    #[test]
    fn file_path_prefers_loaded_source_path() {
        let st = AppState::new(
            run_named("exported-name.xml", 0),
            Vec::new(),
            Some(PathBuf::from("/actual/source.xml")),
            None,
        );

        assert_eq!(st.file_path(), "/actual/source.xml");
        assert_eq!(st.tree_rows().file_path[0], "/actual/source.xml");
    }

    #[test]
    fn tree_label_prefers_opened_path_over_canonical_identity() {
        let mut st = AppState::new(
            run_named("canonical.raw", 1),
            Vec::new(),
            Some(PathBuf::from("/run/canonical.raw")),
            None,
        );
        st.set_active_opened_path(Some(PathBuf::from("/run/run.fa.zip")));

        assert_eq!(st.tree_rows().labels[0], "run.fa.zip");
        assert_eq!(st.opened_path(), Some(PathBuf::from("/run/run.fa.zip")));
        assert_eq!(st.source_path(), Some(PathBuf::from("/run/canonical.raw")));
    }

    #[test]
    fn fragment_analyzer_runs_can_be_renamed_and_saved_in_place() {
        let dir = unique_temp_dir("traceanalyzer_fa_saveable");
        std::fs::create_dir_all(&dir).unwrap();
        let raw = dir.join("run.raw");
        std::fs::write(&raw, b"FA\0\0").unwrap();
        std::fs::write(
            dir.join("run.txt"),
            b"Capillary #: 1\nWell: 1\nSample ID: old\n",
        )
        .unwrap();
        let mut st = fa_state_with_source(raw.to_str().unwrap());

        assert!(st.fragment_analyzer_mode());
        assert!(st.can_rename());
        assert!(st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.can_edit_markers());
        assert!(st.rename_primary("A1"));
        assert!(st.is_dirty());
        assert_eq!(st.run().samples[0].name, "A1");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn raw_extension_source_without_txt_sidecar_is_read_only() {
        let mut st = AppState::new(
            processed_run_named("run.raw", 1),
            Vec::new(),
            Some(PathBuf::from("run.RAW")),
            None,
        );

        assert!(st.fragment_analyzer_mode());
        assert!(!st.can_rename());
        assert!(!st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.rename_primary("A1"));
        assert!(!st.is_dirty());
    }

    #[test]
    fn processed_sample_data_requires_real_trace_arrays() {
        let st = AppState::new(
            run_named("metadata-only.xml", 1),
            Vec::new(),
            Some(PathBuf::from("metadata-only.xml")),
            None,
        );
        assert!(!st.raw_mode());
        assert!(!st.has_processed_sample_data());
        assert!(!st.sample_has_processed_data(0));

        let st = AppState::new(
            processed_run_named("processed.xml", 1),
            Vec::new(),
            Some(PathBuf::from("processed.xml")),
            None,
        );
        assert!(st.has_processed_sample_data());
        assert!(st.sample_has_processed_data(0));
    }

    #[test]
    fn tapestation_exports_are_read_only() {
        let mut st = tapestation_state_with_source("run.xml");

        assert!(!st.can_rename());
        assert!(!st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.can_edit_markers());
        assert!(!st.rename_primary("A1"));
        assert!(!st.is_dirty());
    }

    #[test]
    fn tapestation_assay_name_is_read_only_even_without_regions() {
        let mut run = run_named("export.xml", 1);
        run.assay.assay_name = "D1000".to_string();
        let mut st = AppState::new(run, Vec::new(), Some(PathBuf::from("export.xml")), None);

        assert!(!st.can_rename());
        assert!(!st.can_save());
        assert!(!st.can_save_as());
        assert!(!st.can_edit_markers());
        assert!(!st.rename_primary("A1"));
        assert!(!st.is_dirty());
    }

    #[test]
    fn find_file_by_source_dedups_and_activates() {
        let mut st = AppState::empty();
        st.add_file(
            run_named("a", 1),
            Vec::new(),
            Some(PathBuf::from("/runs/a.raw")),
        );
        st.add_file(
            run_named("b", 1),
            Vec::new(),
            Some(PathBuf::from("/runs/b.raw")),
        );

        assert_eq!(
            st.find_file_by_source(std::path::Path::new("/runs/a.raw")),
            Some(0)
        );
        assert_eq!(
            st.find_file_by_source(std::path::Path::new("/runs/b.raw")),
            Some(1)
        );
        assert_eq!(
            st.find_file_by_source(std::path::Path::new("/runs/c.raw")),
            None
        );

        st.activate_file(0);
        assert_eq!(st.active, Some(0));
    }

    fn unique_temp_dir(prefix: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "{prefix}_{}_{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }
}
