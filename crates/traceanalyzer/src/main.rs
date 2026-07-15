//! Slint viewer for automated-electrophoresis runs.
//!
//! Usage: cargo run -p traceanalyzer -- <file.xad | file.xml | file.xml.gz>
//! With no argument it loads the bundled DNA 1000 demo, if present.

use std::cell::RefCell;
use std::path::PathBuf;
use std::rc::Rc;

use slint::winit_030::winit::event::WindowEvent;
use slint::winit_030::{EventResult, WinitWindowAccessor};
use slint::{
    Image, ModelRc, Rgb8Pixel, SharedPixelBuffer, SharedString, StandardListViewItem, VecModel,
};

use traceanalyzer::plot::{self, Series, Viewport, XAxis};
use traceanalyzer::state::{AppState, Marker, MarkerDrag};
use traceanalyzer::{gel, loading, overview, render, table};
use traceio::calibration::{marker_times, MarkerOverride};

slint::include_modules!();

type SharedState = Rc<RefCell<AppState>>;

fn main() -> anyhow::Result<()> {
    // Let Slint choose the backend. Explicitly probing winit is tempting for
    // drag-drop, but a stale DISPLAY can leave the event loop half-initialized
    // and prevent any backend from starting. The winit file-drop hook below is
    // active when Slint's selected backend is winit.

    let path = match std::env::args().nth(1) {
        Some(p) => PathBuf::from(p),
        None => {
            let demo = PathBuf::from(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../testdata/demo_dna1000.xml.gz"
            ));
            if !demo.exists() {
                anyhow::bail!("no file given and demo fixture not found");
            }
            demo
        }
    };
    let loaded = loading::load(&path)?;
    let state: SharedState = Rc::new(RefCell::new(AppState::new(
        loaded.run,
        loaded.raw_channels,
        Some(path.clone()),
        loaded.warning,
    )));

    let ui = AppWindow::new()?;
    let table = {
        let mut s = state.borrow_mut();
        refresh_all(&ui, &mut s)
    };
    table.apply(&ui);

    // Open a different file via native dialog.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_open_file(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Electrophoresis", &["xad", "xml", "gz"])
                .pick_file()
            {
                reload(&ui, &st, &path);
            }
        });
    }

    // External file drag-and-drop: reload on the last dropped file when the
    // selected Slint backend is winit. Otherwise the Open… dialog remains the
    // portable path.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.window().on_winit_window_event(move |_win, ev| {
            if let WindowEvent::DroppedFile(path) = ev {
                if let Some(ui) = ui_weak.upgrade() {
                    reload(&ui, &st, path);
                }
            }
            EventResult::Propagate
        });
    }

    // Selection (with ctrl/shift modifiers) -> update selection set, reset
    // viewport to auto-fit, and re-render.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_select(move |idx, ctrl, shift| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                if s.select_click(idx as usize, ctrl, shift) {
                    s.viewport = None; // auto-fit the new selection
                }
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            refresh_selection(&ui, &st.borrow());
            show_selected(&ui, &st.borrow());
        });
    }

    // File → Quit.
    ui.on_quit(|| {
        let _ = slint::quit_event_loop();
    });

    // Help → About.
    ui.on_about(|| {
        rfd::MessageDialog::new()
            .set_level(rfd::MessageLevel::Info)
            .set_title("About traceanalyzer")
            .set_description(about_text())
            .set_buttons(rfd::MessageButtons::Ok)
            .show();
    });

    // File → Save (write edits back to the source file).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_save_file(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            let dst = st.borrow().source_path.clone();
            match dst {
                // A native .xad cannot be rewritten in place; offer Save As.
                Some(p) if p.extension().and_then(|e| e.to_str()) == Some("xad") => {
                    save_as_dialog(&ui, &st)
                }
                Some(p) => do_save(&ui, &st, p),
                None => save_as_dialog(&ui, &st),
            }
        });
    }

    // File → Save As….
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_save_file_as(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            save_as_dialog(&ui, &st);
        });
    }

    // Expand/collapse the file node in the well tree.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_file_expand(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().expanded ^= true;
            refresh_tree(&ui, &st.borrow());
        });
    }

    // Rename the selected well in memory (persist later via File → Save).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_rename_well(move |name| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let changed = st.borrow_mut().rename_primary(name.as_str());
            if changed {
                let s = st.borrow();
                refresh_tree(&ui, &s);
                ui.set_window_title(SharedString::from(window_title(&s)));
            }
        });
    }

    // Toggle peak-height normalization for overlaid traces.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_normalize(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                s.normalize = !s.normalize;
                s.viewport = None; // y-range changes with normalization
            }
            show_selected(&ui, &st.borrow());
        });
    }

    // Table row selected -> cross-highlight that peak in the plot.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_table_row_changed(move |row| {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                s.highlight_x = usize::try_from(row)
                    .ok()
                    .and_then(|r| s.table_peak_x.get(r).copied().flatten());
            }
            show_selected(&ui, &st.borrow());
        });
    }

    // Plot clicked -> select the nearest peak's table row and highlight it.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_plot_click(move |fx, fy| {
            let Some(ui) = ui_weak.upgrade() else { return };
            if st.borrow().marker_edit {
                return; // clicks grab markers, not peaks, in marker-edit mode
            }
            let row = {
                let mut s = st.borrow_mut();
                let vp = current_viewport(&s);
                let (data_x, _) = plot::frac_to_data(fx as f64, fy as f64, &vp);
                match nearest_peak_row(&s, data_x, &vp) {
                    Some((row, x)) => {
                        s.highlight_x = Some(x);
                        Some(row)
                    }
                    None => None,
                }
            };
            if let Some(row) = row {
                ui.set_table_current_row(row as i32);
                show_selected(&ui, &st.borrow());
            }
        });
    }

    // Overview: click a mini-plot to open it in Detail.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_overview_click(move |fx, fy| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let idx = {
                let s = st.borrow();
                if s.raw_mode() {
                    None
                } else if s.overview_gel {
                    let (w, _) = gel::size(s.run.samples.len());
                    gel::lane_at(&s.run, fx as f64, w)
                } else {
                    let layout = overview::layout(s.run.samples.len());
                    overview::cell_at(&layout, fx as f64, fy as f64)
                }
            };
            if let Some(idx) = idx {
                let table = {
                    let mut s = st.borrow_mut();
                    s.select_click(idx, false, false);
                    s.viewport = None;
                    build_table_refresh(&mut s)
                };
                table.apply(&ui);
                refresh_selection(&ui, &st.borrow());
                show_selected(&ui, &st.borrow());
                ui.set_active_tab(0); // jump to Detail
            }
        });
    }

    // Overview: toggle shared vs. per-plot y-scale.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_overview_yscale(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().overview_shared_y ^= true;
            refresh_overview(&ui, &st.borrow());
        });
    }

    // Overview: toggle traces grid vs. virtual gel.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_overview_gel(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().overview_gel ^= true;
            refresh_overview(&ui, &st.borrow());
        });
    }

    // Zoom around a fractional x/y anchor (wheel).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_zoom(move |fx, fy, factor| {
            let Some(ui) = ui_weak.upgrade() else { return };
            zoom(&st, fx as f64, fy as f64, factor as f64);
            show_selected(&ui, &st.borrow());
        });
    }
    // Pan by a fractional delta (drag) — or drag a grabbed marker line.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_pan(move |dfx, dfy| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let grabbed = st.borrow().grabbed;
            if let Some(grabbed) = grabbed {
                drag_marker(&ui, &st, grabbed, dfx as f64);
            } else {
                pan(&st, dfx as f64, dfy as f64);
                show_selected(&ui, &st.borrow());
            }
        });
    }
    // Marker-edit: pointer down may grab a marker line near the cursor.
    {
        let st = state.clone();
        ui.on_plot_press(move |fx, _fy| {
            let mut s = st.borrow_mut();
            s.grabbed = None;
            if !s.marker_edit || s.raw_mode() || s.selection.len() != 1 {
                return;
            }
            let vp = current_viewport(&s);
            let idx = s.primary();
            let ov = s.overrides.get(&idx).copied();
            let (lo, up) = marker_times(&s.run, idx, ov.as_ref());
            let tol = 0.02_f64; // fractional grab tolerance
            let mut best: Option<(Marker, f64)> = None;
            for (m, t) in [(Marker::Lower, lo), (Marker::Upper, up)] {
                if let Some(t) = t {
                    let d = (plot::data_x_to_frac(t, &vp) - fx as f64).abs();
                    if d <= tol && best.is_none_or(|b| d < b.1) {
                        best = Some((m, d));
                    }
                }
            }
            s.grabbed = best.map(|b| MarkerDrag {
                sample_idx: idx,
                marker: b.0,
                original_override: ov,
            });
        });
    }
    // Marker-edit: pointer up releases any grabbed marker (and refreshes views).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_plot_release(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                release_marker_drag(&mut s).map(|()| build_table_refresh(&mut s))
            };
            if let Some(table) = table {
                table.apply(&ui);
                ui.set_error_text(SharedString::from(
                    st.borrow().error.clone().unwrap_or_default(),
                ));
                show_selected(&ui, &st.borrow());
                refresh_overview(&ui, &st.borrow());
            }
        });
    }
    // Toggle marker-edit mode (raw-time axis + draggable marker lines).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_marker_edit(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                s.marker_edit ^= true;
                s.viewport = None; // x-axis space changes
                s.grabbed = None;
                s.highlight_x = None;
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            show_selected(&ui, &st.borrow());
        });
    }
    // Reset the focused sample's markers to automatic detection.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_reset_markers(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                let idx = s.primary();
                let previous = s.overrides.remove(&idx);
                s.grabbed = None;
                match commit_recalibration(&mut s) {
                    Ok(()) => {
                        s.error = None;
                        s.viewport = None;
                    }
                    Err(e) => {
                        restore_override(&mut s, idx, previous);
                        s.error = Some(format!("Could not reset markers: {e:#}"));
                    }
                }
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            ui.set_error_text(SharedString::from(
                st.borrow().error.clone().unwrap_or_default(),
            ));
            show_selected(&ui, &st.borrow());
            refresh_overview(&ui, &st.borrow());
        });
    }
    // Reset zoom/pan to auto-fit.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_reset_view(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().viewport = None;
            show_selected(&ui, &st.borrow());
        });
    }
    // Cycle the y-axis quantity.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_cycle_y_mode(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                s.y_mode = s.y_mode.next();
                s.viewport = None; // y-range changes with the quantity
                refresh_overview(&ui, &s);
            }
            show_selected(&ui, &st.borrow());
        });
    }

    show_selected(&ui, &state.borrow());
    ui.run()?;
    Ok(())
}

/// Load a new file into the shared state and refresh the whole UI.
fn reload(ui: &AppWindow, state: &SharedState, path: &std::path::Path) {
    match loading::load(path) {
        Ok(loaded) => {
            let table = {
                let mut s = state.borrow_mut();
                *s = AppState::new(
                    loaded.run,
                    loaded.raw_channels,
                    Some(path.to_path_buf()),
                    loaded.warning,
                );
                refresh_all(ui, &mut s)
            };
            table.apply(ui);
            show_selected(ui, &state.borrow());
        }
        Err(e) => {
            ui.set_error_text(SharedString::from(format!(
                "Could not open {}: {e}",
                path.display()
            )));
        }
    }
}

/// Push run-level data (title, entry list, error) into the UI and build a table update.
fn refresh_all(ui: &AppWindow, st: &mut AppState) -> TableRefresh {
    ui.set_assay_title(SharedString::from(st.title()));
    ui.set_window_title(SharedString::from(window_title(st)));
    ui.set_error_text(SharedString::from(st.error.clone().unwrap_or_default()));
    refresh_tree(ui, st);
    refresh_overview(ui, st);
    build_table_refresh(st)
}

/// Native-window title: app name plus the loaded file, with a `*` when there are
/// unsaved edits.
fn window_title(st: &AppState) -> String {
    let base = if st.run.assay.file_name.is_empty() {
        "traceanalyzer".to_string()
    } else {
        format!("traceanalyzer — {}", st.run.assay.file_name)
    };
    if st.dirty {
        format!("{base} *")
    } else {
        base
    }
}

/// Push the well tree (rows, selection, expansion) and the rename/save enable
/// state into the UI. Also sets the rename field to the selected well's name.
fn refresh_tree(ui: &AppWindow, st: &AppState) {
    let t = st.tree_rows();
    let labels: Vec<SharedString> = t.labels.into_iter().map(SharedString::from).collect();
    ui.set_tree_labels(ModelRc::from(Rc::new(VecModel::from(labels))));
    ui.set_tree_is_file(ModelRc::from(Rc::new(VecModel::from(t.is_file))));
    ui.set_tree_well_index(ModelRc::from(Rc::new(VecModel::from(t.well_index))));
    ui.set_tree_selected(ModelRc::from(Rc::new(VecModel::from(t.selected))));
    ui.set_tree_primary_row(t.primary_row);
    ui.set_tree_expanded(st.expanded);
    ui.set_entry_count(st.entry_count() as i32);
    ui.set_current_index(st.primary() as i32);
    ui.set_can_rename(st.can_rename());
    ui.set_can_save(st.can_save());
    ui.set_rename_text(SharedString::from(st.primary_name()));
}

/// Render the Overview tab (small-multiples grid or virtual gel) into the UI.
fn refresh_overview(ui: &AppWindow, st: &AppState) {
    ui.set_overview_shared_y(st.overview_shared_y);
    ui.set_overview_gel(st.overview_gel);
    if st.raw_mode() {
        let w = plot::PLOT_W;
        let h = plot::PLOT_H;
        let series: Vec<Series> = st.raw_channels.iter().map(plot::raw_series).collect();
        let refs: Vec<&Series> = series.iter().collect();
        let vp = plot::auto_viewport_multi(&refs);
        let buf = plot::render_overlay(&refs, &vp, None, &[], w, h);
        ui.set_overview_image_width(w as i32);
        ui.set_overview_image_height(h as i32);
        ui.set_overview_image(rgb_to_image(&buf, w, h));
    } else if st.overview_gel {
        let (w, h) = gel::size(st.run.samples.len());
        let buf = gel::render(&st.run, w, h);
        ui.set_overview_image_width(w as i32);
        ui.set_overview_image_height(h as i32);
        ui.set_overview_image(rgb_to_image(&buf, w, h));
    } else {
        let layout = overview::layout(st.run.samples.len());
        let buf = overview::render(&st.run, st.y_mode, st.overview_shared_y, &layout);
        ui.set_overview_image_width(layout.w as i32);
        ui.set_overview_image_height(layout.h as i32);
        ui.set_overview_image(rgb_to_image(&buf, layout.w, layout.h));
    }
}

struct TableRefresh {
    rows: ModelRc<ModelRc<StandardListViewItem>>,
    current_row: i32,
}

impl TableRefresh {
    fn empty() -> Self {
        Self {
            rows: ModelRc::from(Rc::new(VecModel::<ModelRc<StandardListViewItem>>::default())),
            current_row: -1,
        }
    }

    fn apply(self, ui: &AppWindow) {
        ui.set_table_rows(self.rows);
        ui.set_table_current_row(self.current_row);
    }
}

/// Build the peak/region table for the focused sample.
/// Also records each row's plot x-position (for cross-highlighting) into `st`.
fn build_table_refresh(st: &mut AppState) -> TableRefresh {
    if st.raw_mode() || st.run.samples.get(st.primary()).is_none() {
        st.table_peak_x.clear();
        return TableRefresh::empty();
    }
    let x_axis = table_x_axis(st);
    let rows = table::rows_with_axis(&st.run, &st.run.samples[st.primary()], x_axis);
    st.table_peak_x = rows.iter().map(|r| r.peak_x).collect();

    let model_rows: Vec<ModelRc<StandardListViewItem>> = rows
        .iter()
        .map(|r| {
            let cells: Vec<StandardListViewItem> = r
                .cells
                .iter()
                .map(|c| StandardListViewItem::from(SharedString::from(c.as_str())))
                .collect();
            ModelRc::from(Rc::new(VecModel::from(cells)))
        })
        .collect();
    TableRefresh {
        rows: ModelRc::from(Rc::new(VecModel::from(model_rows))),
        current_row: -1,
    }
}

/// Push the current selection (tree highlight, primary index, rename field) into
/// the UI.
fn refresh_selection(ui: &AppWindow, st: &AppState) {
    refresh_tree(ui, st);
}

/// Build the plot series for the current selection (raw channel, single sample,
/// or several overlaid samples), applying normalization when requested.
fn selected_series(st: &AppState) -> Vec<Series> {
    if st.raw_mode() {
        return vec![plot::raw_series(&st.raw_channels[st.primary()])];
    }
    let overlay = st.selection.len() > 1;
    st.selection
        .iter()
        .filter_map(|&i| st.run.samples.get(i))
        .map(|s| {
            let x_axis = sample_x_axis(st, s, overlay);
            let series = plot::series(&st.run, s, st.y_mode, x_axis);
            if overlay && st.normalize {
                plot::normalized(&series)
            } else {
                series
            }
        })
        .collect()
}

fn sample_x_axis(st: &AppState, s: &traceio::Sample, overlay: bool) -> XAxis {
    if st.marker_edit && !overlay {
        XAxis::Time
    } else {
        plot::default_x_axis(s)
    }
}

fn table_x_axis(st: &AppState) -> XAxis {
    let overlay = st.selection.len() > 1;
    sample_x_axis(st, &st.run.samples[st.primary()], overlay)
}

/// Effective marker x-positions (raw times) to draw in marker-edit mode.
fn marker_lines(st: &AppState) -> Vec<f64> {
    if !st.marker_edit || st.raw_mode() || st.selection.len() != 1 {
        return Vec::new();
    }
    let idx = st.primary();
    let ov = st.overrides.get(&idx);
    let (lo, up) = traceio::calibration::marker_times(&st.run, idx, ov);
    lo.into_iter().chain(up).collect()
}

/// Render the current selection into the plot.
fn show_selected(ui: &AppWindow, st: &AppState) {
    if st.primary() >= st.entry_count() {
        return;
    }
    ui.set_current_index(st.primary() as i32);

    let series = selected_series(st);
    let refs: Vec<&Series> = series.iter().collect();
    let vp = st
        .viewport
        .unwrap_or_else(|| plot::auto_viewport_multi(&refs));
    // Highlight the selected peak only in the single-sample view.
    let highlight = if st.selection.len() == 1 {
        st.highlight_x
    } else {
        None
    };
    let markers = marker_lines(st);
    let buf = plot::render_overlay(&refs, &vp, highlight, &markers, plot::PLOT_W, plot::PLOT_H);
    ui.set_plot_image(rgb_to_image(&buf, plot::PLOT_W, plot::PLOT_H));

    let info = if st.raw_mode() {
        render::raw_info_line(&st.raw_channels[st.primary()])
    } else if st.selection.len() > 1 {
        format!("{} samples overlaid", st.selection.len())
    } else {
        render::info_line(&st.run, &st.run.samples[st.primary()])
    };
    ui.set_sample_info(SharedString::from(info));
    ui.set_y_mode_label(SharedString::from(st.y_mode.label(&st.run)));
    ui.set_normalize_on(st.normalize);
    ui.set_marker_edit(st.marker_edit);
}

/// Text shown in the Help → About dialog.
fn about_text() -> String {
    format!(
        "traceanalyzer {}\n\nOpen-source post-measurement analysis for automated-electrophoresis runs (Agilent Bioanalyzer, TapeStation, Fragment Analyzer).\n\n© {}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_AUTHORS"),
    )
}

/// Prompt for a destination and save the run there (used by Save As, and by Save
/// when there is no writable source path — e.g. a native `.xad`).
fn save_as_dialog(ui: &AppWindow, st: &SharedState) {
    let start_name = st
        .borrow()
        .source_path
        .as_ref()
        .and_then(|p| p.file_stem())
        .and_then(|n| n.to_str())
        .map(|n| format!("{n}.xml"))
        .unwrap_or_else(|| "run.xml".to_string());

    if let Some(dst) = rfd::FileDialog::new()
        .add_filter("Bioanalyzer XML", &["xml"])
        .add_filter("Bioanalyzer XML (gzip)", &["gz"])
        .set_file_name(start_name)
        .save_file()
    {
        do_save(ui, st, dst);
    }
}

/// Write the current run to `dst`, using the loaded source file as the template,
/// then adopt `dst` as the new source and clear the dirty flag.
fn do_save(ui: &AppWindow, st: &SharedState, dst: std::path::PathBuf) {
    let result = {
        let s = st.borrow();
        match s.source_path.clone() {
            Some(src) => traceio::save::save_run(&s.run, &src, &dst),
            None => Err(anyhow::anyhow!("no source file to save from")),
        }
    };
    match result {
        Ok(()) => {
            {
                let mut s = st.borrow_mut();
                s.dirty = false;
                s.source_path = Some(dst);
            }
            let s = st.borrow();
            ui.set_can_save(s.can_save());
            ui.set_window_title(SharedString::from(window_title(&s)));
            ui.set_error_text(SharedString::default());
        }
        Err(e) => ui.set_error_text(SharedString::from(format!("Save failed: {e:#}"))),
    }
}

fn rgb_to_image(buf: &[u8], w: u32, h: u32) -> Image {
    let mut pixbuf = SharedPixelBuffer::<Rgb8Pixel>::new(w, h);
    pixbuf.make_mut_bytes().copy_from_slice(buf);
    Image::from_rgb8(pixbuf)
}

/// Ensure a concrete viewport exists (materializing the auto-fit), returning it.
fn current_viewport(st: &AppState) -> Viewport {
    if let Some(vp) = st.viewport {
        return vp;
    }
    let series = selected_series(st);
    let refs: Vec<&Series> = series.iter().collect();
    plot::auto_viewport_multi(&refs)
}

/// Find the table row whose peak is nearest `data_x`, within a small tolerance
/// (3% of the visible x-span), returning `(row_index, peak_x)`.
fn nearest_peak_row(st: &AppState, data_x: f64, vp: &Viewport) -> Option<(usize, f64)> {
    let tol = (vp.x_max - vp.x_min).abs() * 0.03;
    st.table_peak_x
        .iter()
        .enumerate()
        .filter_map(|(i, px)| px.map(|x| (i, x)))
        .map(|(i, x)| (i, x, (x - data_x).abs()))
        .filter(|&(_, _, d)| d <= tol)
        .min_by(|a, b| a.2.total_cmp(&b.2))
        .map(|(i, x, _)| (i, x))
}

fn zoom(state: &SharedState, fx: f64, fy: f64, factor: f64) {
    let mut st = state.borrow_mut();
    if st.primary() >= st.entry_count() {
        return;
    }
    let vp = current_viewport(&st);
    let (ax, ay) = plot::frac_to_data(fx, fy, &vp);
    // Zoom towards the anchor: scale spans by `factor`, keep anchor fixed.
    let nx0 = ax - (ax - vp.x_min) * factor;
    let nx1 = ax + (vp.x_max - ax) * factor;
    let ny0 = ay - (ay - vp.y_min) * factor;
    let ny1 = ay + (vp.y_max - ay) * factor;
    st.viewport = Some(Viewport {
        x_min: nx0,
        x_max: nx1,
        y_min: ny0,
        y_max: ny1,
    });
}

fn pan(state: &SharedState, dfx: f64, dfy: f64) {
    let mut st = state.borrow_mut();
    if st.primary() >= st.entry_count() {
        return;
    }
    let vp = current_viewport(&st);
    // Fractional drag across the plot area translates to a data shift.
    let dx = -dfx * (vp.x_max - vp.x_min);
    let dy = dfy * (vp.y_max - vp.y_min); // screen-down is data-up
    st.viewport = Some(Viewport {
        x_min: vp.x_min + dx,
        x_max: vp.x_max + dx,
        y_min: vp.y_min + dy,
        y_max: vp.y_max + dy,
    });
}

/// Move a grabbed marker by a fractional-x drag delta, keeping marker overrides,
/// calibrated trace arrays, table rows, and overview in sync during the drag.
fn drag_marker(ui: &AppWindow, state: &SharedState, drag: MarkerDrag, dfx: f64) {
    let table = {
        let mut st = state.borrow_mut();
        let idx = drag.sample_idx;
        if idx >= st.run.samples.len() || st.raw_mode() {
            return;
        }
        let vp = current_viewport(&st);
        let span = vp.x_max - vp.x_min;
        let ov = st.overrides.get(&idx).copied();
        let (lo, up) = marker_times(&st.run, idx, ov.as_ref());
        let cur = match drag.marker {
            Marker::Lower => lo,
            Marker::Upper => up,
        };
        let Some(cur) = cur else { return };
        let requested = cur + dfx * span;
        let Some(new_time) = valid_marker_time(&st, idx, drag.marker, requested) else {
            return;
        };

        let entry = st.overrides.entry(idx).or_default();
        match drag.marker {
            Marker::Lower => entry.lower_time = Some(new_time),
            Marker::Upper => entry.upper_time = Some(new_time),
        }
        match validate_marker_state(&st, idx, st.overrides.get(&idx).copied())
            .and_then(|()| commit_recalibration(&mut st))
        {
            Ok(()) => {
                st.error = None;
            }
            Err(e) => {
                let restore_msg =
                    restore_override_and_recalibrate(&mut st, idx, drag.original_override);
                st.error = Some(format!(
                    "Could not apply marker override: {e:#}{restore_msg}"
                ));
            }
        }
        build_table_refresh(&mut st)
    };
    table.apply(ui);
    ui.set_error_text(SharedString::from(
        state.borrow().error.clone().unwrap_or_default(),
    ));
    show_selected(ui, &state.borrow());
    refresh_overview(ui, &state.borrow());
}

fn release_marker_drag(st: &mut AppState) -> Option<()> {
    let drag = st.grabbed.take()?;
    let idx = drag.sample_idx;
    let current = st.overrides.get(&idx).copied();
    match validate_marker_state(st, idx, current).and_then(|()| commit_recalibration(st)) {
        Ok(()) => {
            st.error = None;
        }
        Err(e) => {
            let restore_msg = restore_override_and_recalibrate(st, idx, drag.original_override);
            st.error = Some(format!(
                "Could not apply marker override: {e:#}{restore_msg}"
            ));
        }
    }
    Some(())
}

fn commit_recalibration(st: &mut AppState) -> anyhow::Result<()> {
    let mut run = st.run.clone();
    let ov = st.overrides.clone();
    loading::recalibrate_with(&mut run, &ov)?;
    st.run = run;
    st.highlight_x = None;
    Ok(())
}

fn restore_override(st: &mut AppState, idx: usize, previous: Option<MarkerOverride>) {
    if let Some(previous) = previous.filter(|ov| !ov.is_empty()) {
        st.overrides.insert(idx, previous);
    } else {
        st.overrides.remove(&idx);
    }
}

fn restore_override_and_recalibrate(
    st: &mut AppState,
    idx: usize,
    previous: Option<MarkerOverride>,
) -> String {
    restore_override(st, idx, previous);
    commit_recalibration(st)
        .err()
        .map(|restore_err| format!("; restoring previous sizing also failed: {restore_err:#}"))
        .unwrap_or_default()
}

fn validate_marker_state(
    st: &AppState,
    idx: usize,
    ov: Option<MarkerOverride>,
) -> anyhow::Result<()> {
    let (lo, up) = marker_times(&st.run, idx, ov.as_ref());
    for (name, time) in [("lower", lo), ("upper", up)] {
        if let Some(time) = time {
            if !time.is_finite() {
                anyhow::bail!("{name} marker time is not finite");
            }
        }
    }
    if let (Some(lo), Some(up)) = (lo, up) {
        if lo >= up {
            anyhow::bail!("lower marker must be before upper marker");
        }
    }

    let Some(sample) = st.run.samples.get(idx) else {
        anyhow::bail!("sample {idx} is no longer available");
    };
    let mut times = sample.time.iter().copied().filter(|t| t.is_finite());
    let Some(first) = times.next() else {
        return Ok(());
    };
    let (mut min_t, mut max_t) = (first, first);
    for t in times {
        min_t = min_t.min(t);
        max_t = max_t.max(t);
    }
    for (name, time) in [("lower", lo), ("upper", up)] {
        if let Some(time) = time {
            if time < min_t || time > max_t {
                anyhow::bail!("{name} marker must stay within the sample trace");
            }
        }
    }
    Ok(())
}

fn valid_marker_time(st: &AppState, idx: usize, marker: Marker, requested: f64) -> Option<f64> {
    if !requested.is_finite() {
        return None;
    }
    let sample = st.run.samples.get(idx)?;
    let mut times = sample.time.iter().copied().filter(|t| t.is_finite());
    let first = times.next()?;
    let (mut min_t, mut max_t) = (first, first);
    for t in times {
        min_t = min_t.min(t);
        max_t = max_t.max(t);
    }
    if max_t <= min_t {
        return None;
    }

    let ov = st.overrides.get(&idx).copied();
    let (lo, up) = marker_times(&st.run, idx, ov.as_ref());
    let gap = ((max_t - min_t) * 1e-6).max(f64::EPSILON);
    let mut t = requested.clamp(min_t, max_t);
    match marker {
        Marker::Lower => {
            if let Some(up) = up.filter(|v| v.is_finite()) {
                t = t.min(up - gap);
            }
        }
        Marker::Upper => {
            if let Some(lo) = lo.filter(|v| v.is_finite()) {
                t = t.max(lo + gap);
            }
        }
    }
    (t >= min_t && t <= max_t).then_some(t)
}
