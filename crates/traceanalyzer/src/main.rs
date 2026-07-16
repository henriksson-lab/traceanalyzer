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

use traceanalyzer::plot::{self, Series, Viewport, XAxis, YMode};
use traceanalyzer::state::{AppState, Marker, MarkerDrag};
use traceanalyzer::{gel, loading, overview, render, table};
use traceio::calibration::{marker_times, MarkerOverride};
use traceio::Electrophoresis;

slint::include_modules!();

type SharedState = Rc<RefCell<AppState>>;

fn main() -> anyhow::Result<()> {
    // Let Slint choose the backend. Explicitly probing winit is tempting for
    // drag-drop, but a stale DISPLAY can leave the event loop half-initialized
    // and prevent any backend from starting. The winit file-drop hook below is
    // active when Slint's selected backend is winit.

    // Every command-line argument is a file to open (`traceanalyzer a.xad b.xml
    // c.xml.gz`). With no arguments, fall back to the bundled DNA 1000 demo.
    let args: Vec<PathBuf> = std::env::args().skip(1).map(PathBuf::from).collect();
    let paths: Vec<PathBuf> = if args.is_empty() {
        let demo = PathBuf::from(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../testdata/demo_dna1000.xml.gz"
        ));
        if demo.exists() {
            vec![demo]
        } else {
            Vec::new()
        }
    } else {
        args
    };

    // Load each file; a failure to open one is a non-fatal warning (shown in the
    // UI) rather than aborting — the other files still open.
    let mut app = AppState::empty();
    let mut warnings: Vec<String> = Vec::new();
    for path in &paths {
        match loading::load(path) {
            Ok(loaded) => {
                let source = if traceio::fa::is_fa_path(path) {
                    traceio::fa::run_identity(path)
                } else {
                    path.clone()
                };
                if let Some(idx) = app.find_file_by_source(&source) {
                    app.active = Some(idx);
                    continue;
                }
                app.add_file(loaded.run, loaded.raw_channels, Some(source));
                if let Some(w) = loaded.warning {
                    warnings.push(w);
                }
            }
            Err(e) => warnings.push(format!("Could not open {}: {e:#}", path.display())),
        }
    }
    // `add_file` leaves the last-added file active; show the first one instead.
    if !app.files.is_empty() {
        app.active = Some(0);
    }
    if !warnings.is_empty() {
        app.error = Some(warnings.join("\n"));
    }
    let state: SharedState = Rc::new(RefCell::new(app));

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
                .add_filter("Electrophoresis", &["xad", "xml", "gz", "zip", "raw"])
                .pick_file()
            {
                open_added_file(&ui, &st, &path);
            }
        });
    }

    // Open an unzipped Fragment Analyzer run folder directly.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_open_folder(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            if let Some(path) = rfd::FileDialog::new().pick_folder() {
                open_added_file(&ui, &st, &path);
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
                    open_added_file(&ui, &st, path);
                }
            }
            EventResult::Propagate
        });
    }

    // Keyboard selection within the active file (well index + ctrl/shift).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_select(move |idx, ctrl, shift| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                if s.select_click(idx as usize, ctrl, shift) {
                    s.set_viewport(None); // auto-fit the new selection
                }
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            refresh_selection(&ui, &st.borrow());
            show_selected(&ui, &st.borrow());
        });
    }

    // Mouse click on a visible tree row (row index): activates its file and, for
    // a well row, selects the well. A file switch refreshes the whole UI.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_select_row(move |row, ctrl, shift| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                s.select_row(row as usize, ctrl, shift);
                refresh_all(&ui, &mut s)
            };
            table.apply(&ui);
            show_selected(&ui, &st.borrow());
        });
    }

    // Close the file owning a visible row (its [x] button). Prompts to save that
    // file's unsaved edits first; Cancel keeps it open.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_close_file(move |row| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let Some(file_idx) = st.borrow().row_file(row as usize) else {
                return;
            };
            if confirm_close_file(&ui, &st, file_idx) {
                st.borrow_mut().close_file(file_idx);
            }
            // Refresh regardless: the confirm step may have switched the active
            // file, and a close changes the tree/overview/detail.
            let table = {
                let mut s = st.borrow_mut();
                refresh_all(&ui, &mut s)
            };
            table.apply(&ui);
            show_selected(&ui, &st.borrow());
        });
    }

    // File → Quit (prompt to save unsaved edits first).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_quit(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            if confirm_exit(&ui, &st) {
                let _ = slint::quit_event_loop();
            }
        });
    }

    // Window close button (X): same save prompt; veto the close on Cancel.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.window().on_close_requested(move || {
            let Some(ui) = ui_weak.upgrade() else {
                return slint::CloseRequestResponse::HideWindow;
            };
            if confirm_exit(&ui, &st) {
                let _ = slint::quit_event_loop();
                slint::CloseRequestResponse::HideWindow
            } else {
                slint::CloseRequestResponse::KeepWindowShown
            }
        });
    }

    // Help → About.
    ui.on_about(|| {
        rfd::MessageDialog::new()
            .set_level(rfd::MessageLevel::Info)
            .set_title("About Trace analyzer")
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
            save_current(&ui, &st);
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

    // Expand/collapse a file node in the well tree (by visible row index).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_file_expand(move |row| {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().toggle_expand_row(row as usize);
            refresh_tree(&ui, &st.borrow());
        });
    }

    // Rename the selected well in memory (persist later via File → Save), then
    // advance to the next well so a list can be renamed top-to-bottom.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_rename_well(move |name| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                s.rename_primary(name.as_str());
                let next = s.primary() + 1;
                if next < s.entry_count() {
                    s.select_click(next, false, false);
                    s.set_viewport(None); // auto-fit the newly focused well
                }
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            let s = st.borrow();
            refresh_selection(&ui, &s);
            show_selected(&ui, &s);
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
                s.set_viewport(None); // y-range changes with normalization
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
            let hit = {
                let s = st.borrow();
                if s.active_file().is_none() || s.raw_mode() {
                    None
                } else {
                    let (run, origin) = overview_combined(&s);
                    let idx = if s.overview_gel {
                        let (w, _) = gel::size(run.samples.len());
                        gel::lane_at(&run, fx as f64, w)
                    } else {
                        let layout = overview::layout(run.samples.len());
                        overview::cell_at(&layout, fx as f64, fy as f64)
                    };
                    idx.and_then(|i| origin.get(i).copied())
                }
            };
            if let Some((file_idx, well_idx)) = hit {
                let table = {
                    let mut s = st.borrow_mut();
                    s.activate_and_select(file_idx, well_idx);
                    ensure_y_mode_available(&mut s);
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

    // Overview: include/exclude ladder wells.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_overview_ladders(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            st.borrow_mut().overview_show_ladders ^= true;
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
            if !s.marker_edit || s.raw_mode() || s.selection().len() != 1 {
                return;
            }
            let vp = current_viewport(&s);
            let idx = s.primary();
            let ov = s.overrides().get(&idx).copied();
            let (lo, up) = marker_times(s.run(), idx, ov.as_ref());
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
                s.set_viewport(None); // x-axis space changes
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
            if st.borrow().active_file().is_none() {
                return;
            }
            let table = {
                let mut s = st.borrow_mut();
                let idx = s.primary();
                let previous = s.overrides_mut().remove(&idx);
                s.grabbed = None;
                match commit_recalibration(&mut s) {
                    Ok(()) => {
                        s.error = None;
                        s.set_viewport(None);
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
            st.borrow_mut().set_viewport(None);
            show_selected(&ui, &st.borrow());
        });
    }
    // Select the y-axis quantity from the dropdown.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_select_y_mode(move |idx| {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                s.y_mode = YMode::from_available_index(s.run(), idx as usize);
                s.set_viewport(None); // y-range changes with the quantity
                refresh_overview(&ui, &s);
            }
            show_selected(&ui, &st.borrow());
        });
    }

    show_selected(&ui, &state.borrow());
    ui.run()?;
    Ok(())
}

/// Open a file as a new entry in the multi-file list, make it active, and
/// refresh the whole UI. Existing open files are kept.
fn open_added_file(ui: &AppWindow, state: &SharedState, path: &std::path::Path) {
    // Canonical run identity: for a Fragment Analyzer run every entry point (the
    // `.zip`, the `.raw`, the folder, or a sibling like `.PKS`/`.txt`) maps to
    // one identity, so a multi-file drop opens the run just once. Also drives
    // save targeting and the tree tooltip path.
    let source = if traceio::fa::is_fa_path(path) {
        traceio::fa::run_identity(path)
    } else {
        path.to_path_buf()
    };

    // Already open? Re-activate it instead of loading a duplicate.
    if let Some(idx) = state.borrow().find_file_by_source(&source) {
        let table = {
            let mut s = state.borrow_mut();
            s.activate_file(idx);
            refresh_all(ui, &mut s)
        };
        table.apply(ui);
        show_selected(ui, &state.borrow());
        return;
    }

    match loading::load(path) {
        Ok(loaded) => {
            let table = {
                let mut s = state.borrow_mut();
                s.add_file(loaded.run, loaded.raw_channels, Some(source));
                s.error = loaded.warning;
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
    ensure_y_mode_available(st);
    ui.set_assay_title(SharedString::from(st.title()));
    ui.set_error_text(SharedString::from(st.error.clone().unwrap_or_default()));
    refresh_tree(ui, st);
    refresh_overview(ui, st);
    build_table_refresh(st)
}

/// Clamp global y-mode state after switching files or loading data whose
/// derived arrays are absent (for example some Fragment Analyzer imports).
fn ensure_y_mode_available(st: &mut AppState) {
    if st.active_file().is_some() && !st.y_mode.is_available(st.run()) {
        st.y_mode = YMode::Fluorescence;
        st.set_viewport(None);
    }
}

/// Push the well tree (rows, selection, expansion) and the rename/save enable
/// state into the UI. Also sets the rename field to the selected well's name.
fn refresh_tree(ui: &AppWindow, st: &AppState) {
    let t = st.tree_rows();
    let labels: Vec<SharedString> = t.labels.into_iter().map(SharedString::from).collect();
    let file_path: Vec<SharedString> = t.file_path.into_iter().map(SharedString::from).collect();
    ui.set_tree_labels(ModelRc::from(Rc::new(VecModel::from(labels))));
    ui.set_tree_is_file(ModelRc::from(Rc::new(VecModel::from(t.is_file))));
    ui.set_tree_file_index(ModelRc::from(Rc::new(VecModel::from(t.file_index))));
    ui.set_tree_file_expanded(ModelRc::from(Rc::new(VecModel::from(t.file_expanded))));
    ui.set_tree_file_path(ModelRc::from(Rc::new(VecModel::from(file_path))));
    ui.set_tree_selected(ModelRc::from(Rc::new(VecModel::from(t.selected))));
    ui.set_tree_primary_row(t.primary_row);
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
    ui.set_overview_show_ladders(st.overview_show_ladders);
    if st.active_file().is_none() {
        ui.set_overview_image_width(plot::PLOT_W as i32);
        ui.set_overview_image_height(plot::PLOT_H as i32);
        ui.set_overview_image(blank_plot_image());
        return;
    }
    if st.raw_mode() {
        let w = plot::PLOT_W;
        let h = plot::PLOT_H;
        let series: Vec<Series> = st.raw_channels().iter().map(plot::raw_series).collect();
        let refs: Vec<&Series> = series.iter().collect();
        let vp = plot::auto_viewport_multi(&refs);
        let buf = plot::render_overlay(&refs, &vp, None, &[], w, h);
        ui.set_overview_image_width(w as i32);
        ui.set_overview_image_height(h as i32);
        ui.set_overview_image(rgb_to_image(&buf, w, h));
    } else if st.overview_gel {
        let (run, _) = overview_combined(st);
        let (w, h) = gel::size(run.samples.len());
        let buf = gel::render(&run, w, h);
        ui.set_overview_image_width(w as i32);
        ui.set_overview_image_height(h as i32);
        ui.set_overview_image(rgb_to_image(&buf, w, h));
    } else {
        let (run, _) = overview_combined(st);
        let layout = overview::layout(run.samples.len());
        let buf = overview::render(&run, st.y_mode, st.overview_shared_y, &layout);
        ui.set_overview_image_width(layout.w as i32);
        ui.set_overview_image_height(layout.h as i32);
        ui.set_overview_image(rgb_to_image(&buf, layout.w, layout.h));
    }
}

/// Every open file's overview samples flattened into one run for the grid/gel,
/// with a parallel map from each combined sample index to its `(file, well)`
/// origin. Ladder wells are dropped unless the "show ladders" toggle is on; when
/// more than one file is open each sample title is prefixed with its file name so
/// cells stay distinguishable.
fn overview_combined(st: &AppState) -> (Electrophoresis, Vec<(usize, usize)>) {
    let multi = st.files.len() > 1;
    // Use the active file's run as the template (assay units, ladder), then swap
    // in the flattened sample set.
    let mut combined = st.run().clone();
    combined.samples = Vec::new();
    let mut origin: Vec<(usize, usize)> = Vec::new();
    for (fi, f) in st.files.iter().enumerate() {
        let tag = file_base_name(&f.run.assay.file_name);
        for (wi, sample) in f.run.samples.iter().enumerate() {
            if !st.overview_show_ladders && sample.is_ladder {
                continue;
            }
            let mut sample = sample.clone();
            if multi && !tag.is_empty() {
                let base = if sample.name.is_empty() {
                    format!("Well {}", sample.well_number)
                } else {
                    sample.name.clone()
                };
                sample.name = format!("{tag} · {base}");
            }
            combined.samples.push(sample);
            origin.push((fi, wi));
        }
    }
    (combined, origin)
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
    if st.active_file().is_none() || st.raw_mode() || st.run().samples.get(st.primary()).is_none() {
        st.table_peak_x.clear();
        return TableRefresh::empty();
    }
    let x_axis = table_x_axis(st);
    let rows = table::rows_with_axis(st.run(), &st.run().samples[st.primary()], x_axis);
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
        return vec![plot::raw_series(&st.raw_channels()[st.primary()])];
    }
    let overlay = st.selection().len() > 1;
    st.selection()
        .iter()
        .filter_map(|&i| st.run().samples.get(i))
        .map(|s| {
            let x_axis = sample_x_axis(st, s, overlay);
            let series = plot::series(st.run(), s, st.y_mode, x_axis);
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
    let overlay = st.selection().len() > 1;
    sample_x_axis(st, &st.run().samples[st.primary()], overlay)
}

/// Effective marker x-positions (raw times) to draw in marker-edit mode.
fn marker_lines(st: &AppState) -> Vec<f64> {
    if !st.marker_edit || st.raw_mode() || st.selection().len() != 1 {
        return Vec::new();
    }
    let idx = st.primary();
    let ov = st.overrides().get(&idx);
    let (lo, up) = traceio::calibration::marker_times(st.run(), idx, ov);
    lo.into_iter().chain(up).collect()
}

/// Render the current selection into the plot.
fn show_selected(ui: &AppWindow, st: &AppState) {
    if st.active_file().is_none() {
        ui.set_plot_image(blank_plot_image());
        ui.set_sample_info(SharedString::default());
        return;
    }
    if st.primary() >= st.entry_count() {
        return;
    }
    ui.set_current_index(st.primary() as i32);

    let series = selected_series(st);
    let refs: Vec<&Series> = series.iter().collect();
    let vp = st
        .viewport()
        .unwrap_or_else(|| auto_fit_viewport(st, &refs));
    // Highlight the selected peak only in the single-sample view.
    let highlight = if st.selection().len() == 1 {
        st.highlight_x
    } else {
        None
    };
    let markers = marker_lines(st);
    let buf = plot::render_overlay(&refs, &vp, highlight, &markers, plot::PLOT_W, plot::PLOT_H);
    ui.set_plot_image(rgb_to_image(&buf, plot::PLOT_W, plot::PLOT_H));

    let info = if st.raw_mode() {
        render::raw_info_line(&st.raw_channels()[st.primary()])
    } else if st.selection().len() > 1 {
        format!("{} samples overlaid", st.selection().len())
    } else {
        render::info_line(st.run(), &st.run().samples[st.primary()])
    };
    ui.set_sample_info(SharedString::from(info));
    let y_options: Vec<SharedString> = YMode::ALL
        .into_iter()
        .filter(|m| m.is_available(st.run()))
        .map(|m| SharedString::from(m.label(st.run())))
        .collect();
    ui.set_y_mode_options(ModelRc::from(Rc::new(VecModel::from(y_options))));
    ui.set_y_mode_index(st.y_mode.available_index(st.run()) as i32);
    ui.set_normalize_on(st.normalize);
    ui.set_marker_edit(st.marker_edit);
}

/// Text shown in the Help → About dialog.
fn about_text() -> String {
    format!(
        "Trace analyzer {}\n\nOpen-source post-measurement analysis for automated-electrophoresis runs (Agilent Bioanalyzer, TapeStation, Fragment Analyzer).\n\n© {}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_AUTHORS"),
    )
}

/// Prompt for a destination and save the run there (used by Save As, and by Save
/// when there is no writable source path — e.g. a native `.xad`).
/// Save the current run using the same policy as File → Save: rewrite the source
/// file in place when possible, otherwise (no source, or a native `.xad` that
/// can't be rewritten) prompt with Save As. Clears `dirty` on success.
fn save_current(ui: &AppWindow, st: &SharedState) {
    let dst = st.borrow().source_path();
    match dst {
        // A native .xad cannot be rewritten in place; offer Save As.
        Some(p) if p.extension().and_then(|e| e.to_str()) == Some("xad") => save_as_dialog(ui, st),
        Some(p) => do_save(ui, st, p),
        None => save_as_dialog(ui, st),
    }
}

/// Ask Save / Don't Save / Cancel for the ACTIVE file's unsaved edits. Returns
/// `true` if it is OK to proceed (edits saved or explicitly discarded), `false`
/// to abort (Cancel, dismissed, or a failed/aborted save). Assumes the caller
/// has made the file in question active. No borrow held across the dialog.
fn prompt_unsaved_active(ui: &AppWindow, st: &SharedState) -> bool {
    let description = {
        let s = st.borrow();
        let name = file_base_name(&s.file_path());
        if name.is_empty() {
            "This run has unsaved changes.\n\nSave them before closing?".to_string()
        } else {
            format!("“{name}” has unsaved changes.\n\nSave them before closing?")
        }
    };
    let choice = rfd::MessageDialog::new()
        .set_level(rfd::MessageLevel::Warning)
        .set_title("Unsaved changes")
        .set_description(description)
        .set_buttons(rfd::MessageButtons::YesNoCancel) // Yes = Save, No = Don't Save
        .show();
    match choice {
        rfd::MessageDialogResult::Yes => {
            // Attempt the save; only proceed if it actually cleared `dirty`
            // (a failed write or a cancelled Save-As dialog leaves it set).
            save_current(ui, st);
            !st.borrow().is_dirty()
        }
        rfd::MessageDialogResult::No => true, // discard edits
        _ => false,                           // Cancel / dismissed
    }
}

/// Before closing one file, prompt to save its unsaved edits. Makes that file
/// active (so Save / Save-As target it) and returns whether it may be closed.
fn confirm_close_file(ui: &AppWindow, st: &SharedState, file_idx: usize) -> bool {
    {
        let mut s = st.borrow_mut();
        if file_idx < s.files.len() {
            s.active = Some(file_idx);
        }
    }
    if !st.borrow().is_dirty() {
        return true;
    }
    prompt_unsaved_active(ui, st)
}

/// On exit, prompt per dirty file (Save / Don't Save / Cancel). Any Cancel or a
/// failed save aborts the whole quit and keeps the window open.
fn confirm_exit(ui: &AppWindow, st: &SharedState) -> bool {
    let dirty: Vec<usize> = {
        let s = st.borrow();
        (0..s.files.len()).filter(|&i| s.files[i].dirty).collect()
    };
    for idx in dirty {
        {
            let mut s = st.borrow_mut();
            if idx >= s.files.len() || !s.files[idx].dirty {
                continue; // already saved, or vanished
            }
            s.active = Some(idx); // target this file for the prompt/save
        }
        if !prompt_unsaved_active(ui, st) {
            return false;
        }
    }
    true
}

/// Last path component of an instrument file path (Windows `\` or Unix `/`).
fn file_base_name(path: &str) -> String {
    path.rsplit(['\\', '/']).next().unwrap_or(path).to_string()
}

fn save_as_dialog(ui: &AppWindow, st: &SharedState) {
    if st.borrow().fragment_analyzer_mode() {
        ui.set_error_text(SharedString::from(
            "Fragment Analyzer runs can only be saved in place by updating the .txt sidecar",
        ));
        return;
    }

    let start_name = st
        .borrow()
        .source_path()
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
        match s.source_path() {
            Some(src) => traceio::save::save_run(s.run(), &src, &dst),
            None => Err(anyhow::anyhow!("no source file to save from")),
        }
    };
    match result {
        Ok(()) => {
            {
                let mut s = st.borrow_mut();
                s.set_dirty(false);
                s.set_source_path(Some(dst));
            }
            let s = st.borrow();
            ui.set_can_save(s.can_save());
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

/// A plot-sized image with just empty axes, shown when no file is open.
fn blank_plot_image() -> Image {
    let refs: Vec<&Series> = Vec::new();
    let vp = plot::auto_viewport_multi(&refs);
    let buf = plot::render_overlay(&refs, &vp, None, &[], plot::PLOT_W, plot::PLOT_H);
    rgb_to_image(&buf, plot::PLOT_W, plot::PLOT_H)
}

/// Ensure a concrete viewport exists (materializing the auto-fit), returning it.
fn current_viewport(st: &AppState) -> Viewport {
    if st.active_file().is_none() {
        let refs: Vec<&Series> = Vec::new();
        return plot::auto_viewport_multi(&refs);
    }
    if let Some(vp) = st.viewport() {
        return vp;
    }
    let series = selected_series(st);
    let refs: Vec<&Series> = series.iter().collect();
    auto_fit_viewport(st, &refs)
}

/// Auto-fit viewport for the current selection, using robust y-scaling for the
/// derived quantities (concentration/molarity) so a numeric spike at tiny sizes
/// doesn't squash the trace. Keeps display and interaction viewports in sync.
fn auto_fit_viewport(st: &AppState, refs: &[&Series]) -> Viewport {
    let robust_y = !st.raw_mode() && matches!(st.y_mode, YMode::Concentration | YMode::Molarity);
    if robust_y {
        plot::auto_viewport_multi_robust(refs)
    } else {
        plot::auto_viewport_multi(refs)
    }
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
    st.set_viewport(Some(Viewport {
        x_min: nx0,
        x_max: nx1,
        y_min: ny0,
        y_max: ny1,
    }));
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
    st.set_viewport(Some(Viewport {
        x_min: vp.x_min + dx,
        x_max: vp.x_max + dx,
        y_min: vp.y_min + dy,
        y_max: vp.y_max + dy,
    }));
}

/// Move a grabbed marker by a fractional-x drag delta, keeping marker overrides,
/// calibrated trace arrays, table rows, and overview in sync during the drag.
fn drag_marker(ui: &AppWindow, state: &SharedState, drag: MarkerDrag, dfx: f64) {
    let table = {
        let mut st = state.borrow_mut();
        if st.active_file().is_none() {
            return;
        }
        let idx = drag.sample_idx;
        if idx >= st.run().samples.len() || st.raw_mode() {
            return;
        }
        let vp = current_viewport(&st);
        let span = vp.x_max - vp.x_min;
        let ov = st.overrides().get(&idx).copied();
        let (lo, up) = marker_times(st.run(), idx, ov.as_ref());
        let cur = match drag.marker {
            Marker::Lower => lo,
            Marker::Upper => up,
        };
        let Some(cur) = cur else { return };
        let requested = cur + dfx * span;
        let Some(new_time) = valid_marker_time(&st, idx, drag.marker, requested) else {
            return;
        };

        let entry = st.overrides_mut().entry(idx).or_default();
        match drag.marker {
            Marker::Lower => entry.lower_time = Some(new_time),
            Marker::Upper => entry.upper_time = Some(new_time),
        }
        match validate_marker_state(&st, idx, st.overrides().get(&idx).copied())
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
    let current = st.overrides().get(&idx).copied();
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
    let mut run = st.run().clone();
    let ov = st.overrides().clone();
    loading::recalibrate_with(&mut run, &ov)?;
    st.set_run(run);
    st.highlight_x = None;
    Ok(())
}

fn restore_override(st: &mut AppState, idx: usize, previous: Option<MarkerOverride>) {
    if let Some(previous) = previous.filter(|ov| !ov.is_empty()) {
        st.overrides_mut().insert(idx, previous);
    } else {
        st.overrides_mut().remove(&idx);
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
    let (lo, up) = marker_times(st.run(), idx, ov.as_ref());
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

    let Some(sample) = st.run().samples.get(idx) else {
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
    let sample = st.run().samples.get(idx)?;
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

    let ov = st.overrides().get(&idx).copied();
    let (lo, up) = marker_times(st.run(), idx, ov.as_ref());
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
