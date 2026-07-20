//! Slint viewer for automated-electrophoresis runs.
//!
//! Usage: cargo run -p traceanalyzer -- <file.xad | file.xml | file.xml.gz>
//! With no argument it loads the bundled DNA 1000 demo, if present.

use std::cell::RefCell;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::rc::Rc;

use slint::winit_030::winit::event::WindowEvent;
use slint::winit_030::{EventResult, WinitWindowAccessor};
use slint::{
    Image, ModelRc, Rgb8Pixel, SharedPixelBuffer, SharedString, StandardListViewItem, VecModel,
};

use traceanalyzer::plot::{self, Series, Viewport, XAxis, YMode};
use traceanalyzer::state::{AppState, Marker, MarkerDrag};
use traceanalyzer::{export, gel, loading, overview, render, table};
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
                let source = loaded.source.identity.clone();
                if let Some(idx) = app.find_file_by_source(&source) {
                    app.active = Some(idx);
                    continue;
                }
                let source_capabilities = loaded.source.capabilities;
                let save_capabilities = loaded.source.save_capabilities();
                let fa_metadata = loaded.fa_metadata.clone();
                app.add_file_with_capabilities(
                    loaded.run,
                    loaded.raw_channels,
                    Some(source),
                    source_capabilities,
                    save_capabilities,
                );
                app.set_active_opened_path(Some(path.to_path_buf()));
                if let Some(metadata) = loaded_fa_metadata(path, fa_metadata) {
                    app.set_active_metadata_text(metadata.text);
                    app.set_active_metadata_rows(metadata.rows);
                    if let Some(warning) = metadata.warning {
                        warnings.push(warning);
                    }
                }
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
                .add_filter(
                    "All supported runs",
                    &["xad", "xml", "gz", "fa.zip", "zip", "raw", "csv"],
                )
                .add_filter("Bioanalyzer", &["xad", "xml", "gz"])
                .add_filter("TapeStation export", &["xml", "gz", "csv"])
                .add_filter("Fragment Analyzer", &["fa.zip", "zip", "raw"])
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

    // File -> Export Peak Table CSV...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_peak_table(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_peak_table_dialog(&ui, &st);
        });
    }

    // File -> Export All Peak Tables CSV...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_run_peak_table(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_run_peak_table_dialog(&ui, &st);
        });
    }

    // File -> Export Trace Data CSV...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_trace_data(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_trace_data_dialog(&ui, &st);
        });
    }

    // File -> Export Raw Detector Channels CSV...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_raw_channels(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_raw_channels_dialog(&ui, &st);
        });
    }

    // File -> Export Metadata/QC CSV...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_metadata_qc(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_metadata_qc_dialog(&ui, &st);
        });
    }

    // File -> Export Plot PNG...
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_export_plot_png(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            export_plot_png_dialog(&ui, &st);
        });
    }

    // Expand/collapse a file node in the well tree (by visible row index).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_toggle_file_expand(move |row| {
            let Some(ui) = ui_weak.upgrade() else { return };
            let table = {
                let mut s = st.borrow_mut();
                s.toggle_expand_row(row as usize);
                refresh_all(&ui, &mut s)
            };
            table.apply(&ui);
            show_selected(&ui, &st.borrow());
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
                let renamed = s.rename_primary(name.as_str());
                let next = s.primary() + 1;
                if renamed && next < s.entry_count() {
                    s.select_click(next, false, false);
                    s.set_viewport(None); // auto-fit the newly focused well
                }
                build_table_refresh(&mut s)
            };
            table.apply(&ui);
            let s = st.borrow();
            refresh_selection(&ui, &s);
            refresh_metadata(&ui, &s);
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
                if s.active_file().is_none() {
                    return;
                }
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
                    refresh_all(&ui, &mut s)
                };
                table.apply(&ui);
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
            refresh_tree(&ui, &st.borrow());
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
            refresh_tree(&ui, &st.borrow());
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
            refresh_tree(&ui, &st.borrow());
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
                refresh_tree(&ui, &st.borrow());
                ui.set_error_text(SharedString::from(
                    st.borrow().error.clone().unwrap_or_default(),
                ));
                show_selected(&ui, &st.borrow());
                refresh_metadata(&ui, &st.borrow());
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
                if !s.can_edit_markers() {
                    return;
                }
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
            if !st.borrow().can_edit_markers() {
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
            refresh_tree(&ui, &st.borrow());
            ui.set_error_text(SharedString::from(
                st.borrow().error.clone().unwrap_or_default(),
            ));
            show_selected(&ui, &st.borrow());
            refresh_metadata(&ui, &st.borrow());
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
                if s.active_file().is_none() {
                    return;
                }
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
    if !is_supported_open_path(path) {
        set_open_error(
            ui,
            state,
            format!(
                "Unsupported input: {}. Open a Bioanalyzer .xad/.xml/.xml.gz file, a TapeStation exported .xml or _Electropherogram.csv file, or a Fragment Analyzer .fa.zip/.raw/run folder.",
                path.display()
            ),
        );
        return;
    }

    match loading::load(path) {
        Ok(loaded) => {
            let source = loaded.source.identity.clone();
            // Already open? Re-activate it instead of loading a duplicate.
            if let Some(idx) = state.borrow().find_file_by_source(&source) {
                let table = {
                    let mut s = state.borrow_mut();
                    s.activate_file(idx);
                    s.error = None;
                    s.status = Some(format!("Already open: {}", source.display()));
                    refresh_all(ui, &mut s)
                };
                table.apply(ui);
                show_selected(ui, &state.borrow());
                return;
            }
            let table = {
                let mut s = state.borrow_mut();
                let source_capabilities = loaded.source.capabilities;
                let save_capabilities = loaded.source.save_capabilities();
                let fa_metadata = loaded.fa_metadata.clone();
                s.add_file_with_capabilities(
                    loaded.run,
                    loaded.raw_channels,
                    Some(source),
                    source_capabilities,
                    save_capabilities,
                );
                s.set_active_opened_path(Some(path.to_path_buf()));
                let metadata_warning = loaded_fa_metadata(path, fa_metadata).and_then(|metadata| {
                    s.set_active_metadata_text(metadata.text);
                    s.set_active_metadata_rows(metadata.rows);
                    metadata.warning
                });
                s.error = join_optional_messages([loaded.warning, metadata_warning]);
                s.status = s
                    .error
                    .is_none()
                    .then(|| format!("Opened {}", path.display()));
                refresh_all(ui, &mut s)
            };
            table.apply(ui);
            show_selected(ui, &state.borrow());
        }
        Err(e) => {
            set_open_error(
                ui,
                state,
                format!("Could not open {}: {e:#}", path.display()),
            );
        }
    }
}

fn set_open_error(ui: &AppWindow, state: &SharedState, message: String) {
    {
        let mut s = state.borrow_mut();
        s.error = Some(message.clone());
        s.status = None;
    }
    ui.set_error_text(SharedString::from(message));
    ui.set_status_text(SharedString::default());
}

fn is_supported_open_path(path: &std::path::Path) -> bool {
    traceio::io::detect_format(path).is_ok_and(|detected| detected.is_some())
}

struct LoadedFaMetadata {
    text: Option<String>,
    rows: Vec<(String, String)>,
    warning: Option<String>,
}

fn loaded_fa_metadata(
    path: &std::path::Path,
    metadata: Option<traceio::fa::FaMetadata>,
) -> Option<LoadedFaMetadata> {
    if !traceio::fa::is_fa_path(path) {
        return None;
    }
    metadata.map(|metadata| LoadedFaMetadata {
        text: format_fa_metadata(&metadata),
        rows: summarize_fa_metadata(&metadata),
        warning: None,
    })
}

fn summarize_fa_metadata(metadata: &traceio::fa::FaMetadata) -> Vec<(String, String)> {
    let mut rows = Vec::new();
    for (key, value) in &metadata.run_header {
        if !value.trim().is_empty() {
            rows.push((format!("Run header: {key}"), value.clone()));
        }
    }
    if !metadata.method.is_empty() {
        rows.push((
            "Method sections".to_string(),
            metadata.method.len().to_string(),
        ));
    }
    if !metadata.analysis.is_empty() {
        rows.push((
            "Analysis sections".to_string(),
            metadata.analysis.len().to_string(),
        ));
    }
    if !metadata.current.is_empty() {
        rows.push((
            "Current log points".to_string(),
            metadata.current.len().to_string(),
        ));
        push_fa_current_range(
            &mut rows,
            "Current range",
            metadata.current.iter().map(|p| p.current_ua),
            "uA",
        );
        push_fa_current_range(
            &mut rows,
            "Voltage range",
            metadata.current.iter().map(|p| p.voltage_kv),
            "kV",
        );
        push_fa_current_range(
            &mut rows,
            "Pressure range",
            metadata.current.iter().map(|p| p.pressure_psi),
            "psi",
        );
    }
    if !metadata.timing.is_empty() {
        rows.push((
            "Timing lines".to_string(),
            metadata.timing.len().to_string(),
        ));
    }
    if !metadata.exposure.is_empty() {
        rows.push((
            "Exposure lines".to_string(),
            metadata.exposure.len().to_string(),
        ));
    }
    if let Some(image) = &metadata.camera_image {
        rows.push(("Camera image bytes".to_string(), image.len().to_string()));
    }
    rows
}

fn push_fa_current_range(
    rows: &mut Vec<(String, String)>,
    label: &str,
    values: impl Iterator<Item = f64>,
    unit: &str,
) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for value in values.filter(|v| v.is_finite()) {
        min = min.min(value);
        max = max.max(value);
    }
    if min.is_finite() && max.is_finite() {
        rows.push((label.to_string(), format!("{min:.3} to {max:.3} {unit}")));
    }
}

fn format_fa_metadata(metadata: &traceio::fa::FaMetadata) -> Option<String> {
    let mut out = String::new();

    if !metadata.run_header.is_empty() {
        out.push_str("Run header\n");
        for (key, value) in &metadata.run_header {
            let _ = writeln!(out, "{key}: {value}");
        }
        out.push('\n');
    }

    write_ini_summary(&mut out, "Method", &metadata.method);
    write_ini_summary(&mut out, "Analysis", &metadata.analysis);

    if !metadata.current.is_empty() {
        let mut current_min = f64::INFINITY;
        let mut current_max = f64::NEG_INFINITY;
        let mut voltage_min = f64::INFINITY;
        let mut voltage_max = f64::NEG_INFINITY;
        let mut pressure_min = f64::INFINITY;
        let mut pressure_max = f64::NEG_INFINITY;
        for point in &metadata.current {
            current_min = current_min.min(point.current_ua);
            current_max = current_max.max(point.current_ua);
            voltage_min = voltage_min.min(point.voltage_kv);
            voltage_max = voltage_max.max(point.voltage_kv);
            pressure_min = pressure_min.min(point.pressure_psi);
            pressure_max = pressure_max.max(point.pressure_psi);
        }
        out.push_str("Current log\n");
        let _ = writeln!(out, "Points: {}", metadata.current.len());
        let _ = writeln!(out, "Current: {current_min:.3} to {current_max:.3} uA");
        let _ = writeln!(out, "Voltage: {voltage_min:.3} to {voltage_max:.3} kV");
        let _ = writeln!(out, "Pressure: {pressure_min:.3} to {pressure_max:.3} psi");
        out.push('\n');
    }

    if !metadata.timing.is_empty() {
        out.push_str("Timing\n");
        let _ = writeln!(out, "Lines: {}", metadata.timing.len());
        for line in metadata.timing.iter().take(12) {
            let _ = writeln!(out, "{line}");
        }
        if metadata.timing.len() > 12 {
            let _ = writeln!(out, "... {} more lines", metadata.timing.len() - 12);
        }
        out.push('\n');
    }

    if !metadata.exposure.is_empty() {
        out.push_str("Exposure\n");
        let _ = writeln!(out, "Lines: {}", metadata.exposure.len());
        for line in metadata.exposure.iter().take(12) {
            let _ = writeln!(out, "{line}");
        }
        if metadata.exposure.len() > 12 {
            let _ = writeln!(out, "... {} more lines", metadata.exposure.len() - 12);
        }
        out.push('\n');
    }

    if let Some(image) = &metadata.camera_image {
        out.push_str("Camera image\n");
        let _ = writeln!(out, "Bytes: {}", image.len());
    }

    let trimmed = out.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}

fn write_ini_summary(out: &mut String, title: &str, document: &traceio::fa::IniDocument) {
    if document.is_empty() {
        return;
    }

    let _ = writeln!(out, "{title}");
    for (section, values) in document {
        let _ = writeln!(out, "[{section}]");
        for (key, value) in values {
            let _ = writeln!(out, "{key}: {value}");
        }
    }
    out.push('\n');
}

fn join_optional_messages(messages: impl IntoIterator<Item = Option<String>>) -> Option<String> {
    let joined = messages
        .into_iter()
        .flatten()
        .filter(|message| !message.is_empty())
        .collect::<Vec<_>>()
        .join("\n");
    (!joined.is_empty()).then_some(joined)
}

/// Push run-level data (title, entry list, error) into the UI and build a table update.
fn refresh_all(ui: &AppWindow, st: &mut AppState) -> TableRefresh {
    ensure_y_mode_available(st);
    ensure_marker_edit_available(st);
    ui.set_assay_title(SharedString::from(st.title()));
    ui.set_error_text(SharedString::from(st.error.clone().unwrap_or_default()));
    ui.set_status_text(SharedString::from(st.status.clone().unwrap_or_default()));
    ui.set_metadata_text(SharedString::from(st.metadata_text()));
    refresh_metadata(ui, st);
    refresh_tree(ui, st);
    refresh_overview(ui, st);
    build_table_refresh(st)
}

fn refresh_metadata(ui: &AppWindow, st: &AppState) {
    let (summary, run_rows) = run_metadata_rows(st);
    ui.set_metadata_run_summary(SharedString::from(summary));
    ui.set_metadata_run_rows(kv_rows_model(&run_rows));

    let source_title = if st.active_file().is_none() {
        ""
    } else if st.fragment_analyzer_mode() {
        "Fragment Analyzer sidecar"
    } else {
        "Source-specific metadata"
    };
    ui.set_metadata_source_title(SharedString::from(source_title));
    ui.set_metadata_source_rows(kv_rows_model(st.metadata_rows()));
}

fn run_metadata_rows(st: &AppState) -> (String, Vec<(String, String)>) {
    let Some(file) = st.active_file() else {
        return (String::new(), Vec::new());
    };
    let run = &file.run;
    let sample_count = run.samples.len();
    let ladder_count = run.samples.iter().filter(|s| s.is_ladder).count();
    let peak_count: usize = run.samples.iter().map(|s| s.peaks.len()).sum();
    let sample_region_count: usize = run.samples.iter().map(|s| s.regions.len()).sum();
    let trace_point_count: usize = run.samples.iter().map(|s| s.fluorescence.len()).sum();
    let raw_mode = st.raw_mode();

    let mode = if raw_mode {
        "raw acquisition"
    } else {
        "processed samples"
    };
    let summary = format!(
        "{} samples, {} ladders, {} peaks ({mode})",
        sample_count, ladder_count, peak_count
    );

    let mut rows = Vec::new();
    push_metadata_row(&mut rows, "File name", &run.assay.file_name);
    if let Some(source) = file.opened_path.as_ref().or(file.source_path.as_ref()) {
        rows.push(("Loaded from".to_string(), source.display().to_string()));
    }
    push_metadata_row(&mut rows, "Created", &run.assay.creation_date);
    push_metadata_row(&mut rows, "Assay name", &run.assay.assay_name);
    push_metadata_row(&mut rows, "Assay type", &run.assay.assay_type);
    push_metadata_row(&mut rows, "Length unit", &run.assay.length_unit);
    push_metadata_row(
        &mut rows,
        "Concentration unit",
        &run.assay.concentration_unit,
    );
    if let Some(unit) = &run.assay.molarity_unit {
        push_metadata_row(&mut rows, "Molarity unit", unit);
    }
    rows.push(("Samples".to_string(), sample_count.to_string()));
    rows.push(("Ladder samples".to_string(), ladder_count.to_string()));
    rows.push(("Reported peaks".to_string(), peak_count.to_string()));
    rows.push((
        "Ladder peaks".to_string(),
        run.ladder_peaks.len().to_string(),
    ));
    rows.push((
        "Analysis regions".to_string(),
        (run.regions.len() + sample_region_count).to_string(),
    ));
    if raw_mode {
        rows.push((
            "Raw channels".to_string(),
            file.raw_channels.len().to_string(),
        ));
    } else {
        rows.push(("Trace points".to_string(), trace_point_count.to_string()));
    }
    rows.push((
        "Upper marker".to_string(),
        if run.assay.has_upper_marker {
            "yes"
        } else {
            "no"
        }
        .to_string(),
    ));
    rows.push((
        "Unsaved edits".to_string(),
        if file.dirty { "yes" } else { "no" }.to_string(),
    ));
    rows.push((
        "Session marker overrides".to_string(),
        if st.has_marker_overrides() {
            "yes; recalibration uses manual marker positions that Save/Save As does not persist"
        } else {
            "no"
        }
        .to_string(),
    ));
    (summary, rows)
}

fn push_metadata_row(rows: &mut Vec<(String, String)>, label: &str, value: &str) {
    let value = value.trim();
    if !value.is_empty() {
        rows.push((label.to_string(), value.to_string()));
    }
}

/// Clamp global y-mode state after switching files or loading data whose
/// derived arrays are absent (for example some Fragment Analyzer imports).
fn ensure_y_mode_available(st: &mut AppState) {
    if st.active_file().is_some() && !st.y_mode.is_available(st.run()) {
        st.y_mode = YMode::Fluorescence;
        st.set_viewport(None);
    }
}

fn ensure_marker_edit_available(st: &mut AppState) {
    if st.marker_edit && !st.can_edit_markers() {
        st.marker_edit = false;
        st.grabbed = None;
        st.highlight_x = None;
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
    ui.set_can_save(st.is_dirty() && st.can_save());
    ui.set_can_save_as(st.can_save_as());
    ui.set_can_export_table(can_export_peak_table(st));
    ui.set_can_export_run_table(can_export_run_peak_table(st));
    ui.set_can_export_traces(can_export_trace_data(st));
    ui.set_can_export_raw_channels(can_export_raw_channels(st));
    ui.set_can_export_metadata(can_export_metadata_qc(st));
    ui.set_can_export_plot(can_export_plot(st));
    ui.set_can_export_detail_plot(can_export_plot_for_tab(st, ExportPlotTab::Detail));
    ui.set_can_export_overview_plot(can_export_plot_for_tab(st, ExportPlotTab::Overview));
    ui.set_can_edit_markers(st.can_edit_markers());
    ui.set_rename_text(SharedString::from(st.primary_name()));
}

/// Render the Overview tab (small-multiples grid or virtual gel) into the UI.
fn refresh_overview(ui: &AppWindow, st: &AppState) {
    ui.set_overview_shared_y(st.overview_shared_y);
    ui.set_overview_gel(st.overview_gel);
    ui.set_overview_show_ladders(st.overview_show_ladders);
    let (buf, w, h) = render_overview_bitmap(st).unwrap_or_else(blank_plot_bitmap);
    ui.set_overview_image_width(w as i32);
    ui.set_overview_image_height(h as i32);
    ui.set_overview_image(rgb_to_image(&buf, w, h));
}

fn render_overview_bitmap(st: &AppState) -> Option<(Vec<u8>, u32, u32)> {
    st.active_file()?;

    if st.raw_mode() {
        let w = plot::PLOT_W;
        let h = plot::PLOT_H;
        let series: Vec<Series> = st.raw_channels().iter().map(plot::raw_series).collect();
        if !series.iter().any(series_has_points) {
            return None;
        }
        let refs: Vec<&Series> = series.iter().collect();
        let vp = plot::auto_viewport_multi(&refs);
        let buf = plot::render_overlay(&refs, &vp, None, &[], w, h);
        Some((buf, w, h))
    } else if st.overview_gel {
        let (run, _) = overview_combined(st);
        if !run.samples.iter().any(sample_has_gel_trace_data) {
            return None;
        }
        let (w, h) = gel::size(run.samples.len());
        let buf = gel::render(&run, w, h);
        Some((buf, w, h))
    } else {
        let (run, _) = overview_combined(st);
        if !run
            .samples
            .iter()
            .any(|sample| sample_has_overview_trace_data(&run, sample, st.y_mode))
        {
            return None;
        }
        let layout = overview::layout(run.samples.len());
        let buf = overview::render(&run, st.y_mode, st.overview_shared_y, &layout);
        Some((buf, layout.w, layout.h))
    }
}

fn series_has_points(series: &Series) -> bool {
    !series.xs.is_empty() && !series.ys.is_empty()
}

fn sample_has_overview_trace_data(
    run: &Electrophoresis,
    sample: &traceio::Sample,
    y_mode: YMode,
) -> bool {
    series_has_points(&plot::series(run, sample, y_mode, false))
}

fn sample_has_gel_trace_data(sample: &traceio::Sample) -> bool {
    sample.time.iter().any(|time| time.is_finite())
        && sample
            .fluorescence
            .iter()
            .any(|fluorescence| (*fluorescence as f64).is_finite())
}

/// Every open file's overview samples flattened into one run for the grid/gel,
/// with a parallel map from each combined sample index to its `(file, well)`
/// origin. Ladder wells are dropped unless the "show ladders" toggle is on; when
/// more than one file is open each sample title is prefixed with its file name so
/// cells stay distinguishable.
fn overview_combined(st: &AppState) -> (Electrophoresis, Vec<(usize, usize)>) {
    let multi = st.files.len() > 1;
    let file_labels = overview_file_labels(st);
    // Use the active file's run as the template (assay units, ladder), then swap
    // in the flattened sample set.
    let mut combined = st.run().clone();
    combined.samples = Vec::new();
    let mut origin: Vec<(usize, usize)> = Vec::new();
    for (fi, f) in st.files.iter().enumerate() {
        let tag = file_labels.get(fi).map(String::as_str).unwrap_or_default();
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

fn overview_file_labels(st: &AppState) -> Vec<String> {
    let rows = st.tree_rows();
    rows.labels
        .into_iter()
        .zip(rows.is_file)
        .zip(rows.file_index)
        .filter_map(|((label, is_file), file_index)| {
            if is_file {
                Some((file_index as usize, label))
            } else {
                None
            }
        })
        .fold(
            vec![String::new(); st.files.len()],
            |mut labels, (idx, label)| {
                if let Some(slot) = labels.get_mut(idx) {
                    *slot = label;
                }
                labels
            },
        )
}

fn overview_marker_overrides_active(st: &AppState) -> bool {
    let (_, origin) = overview_combined(st);
    origin
        .iter()
        .any(|&(file_idx, sample_idx)| st.file_sample_has_marker_override(file_idx, sample_idx))
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

fn kv_rows_model(rows: &[(String, String)]) -> ModelRc<ModelRc<StandardListViewItem>> {
    let model_rows: Vec<ModelRc<StandardListViewItem>> = rows
        .iter()
        .map(|(key, value)| {
            let cells = vec![
                StandardListViewItem::from(SharedString::from(key.as_str())),
                StandardListViewItem::from(SharedString::from(value.as_str())),
            ];
            ModelRc::from(Rc::new(VecModel::from(cells)))
        })
        .collect();
    ModelRc::from(Rc::new(VecModel::from(model_rows)))
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
    if st.active_file().is_none() || st.primary() >= st.entry_count() {
        clear_detail(ui);
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

fn clear_detail(ui: &AppWindow) {
    ui.set_plot_image(blank_plot_image());
    ui.set_sample_info(SharedString::default());
    ui.set_y_mode_options(ModelRc::from(Rc::new(VecModel::from(
        Vec::<SharedString>::new(),
    ))));
    ui.set_y_mode_index(-1);
    ui.set_current_index(-1);
    ui.set_normalize_on(false);
    ui.set_marker_edit(false);
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
        Some(p) if is_xad_path(&p) => save_as_dialog(ui, st),
        Some(p) => do_save(ui, st, p),
        None => save_as_dialog(ui, st),
    }
}

fn is_xad_path(path: &std::path::Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("xad"))
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

/// Before closing one file, prompt to save its unsaved edits. Temporarily makes
/// that file active so Save / Save-As target it, then restores the previous
/// active file before the caller removes an inactive file.
fn confirm_close_file(ui: &AppWindow, st: &SharedState, file_idx: usize) -> bool {
    let previous_active = st.borrow().active;
    let switched_active = previous_active != Some(file_idx);
    {
        let mut s = st.borrow_mut();
        if file_idx >= s.files.len() {
            return false;
        }
        s.active = Some(file_idx);
    }
    if !st.borrow().is_dirty() {
        if switched_active {
            st.borrow_mut().active = previous_active;
        }
        return true;
    }
    let confirmed = prompt_unsaved_active(ui, st);
    if switched_active {
        st.borrow_mut().active = previous_active;
    }
    confirmed
}

/// On exit, prompt per dirty file (Save / Don't Save / Cancel). Any Cancel or a
/// failed save aborts the whole quit and keeps the window open.
fn confirm_exit(ui: &AppWindow, st: &SharedState) -> bool {
    let previous_active = st.borrow().active;
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
            restore_active_file_ui(ui, st, previous_active);
            return false;
        }
    }
    true
}

fn restore_active_file_ui(ui: &AppWindow, st: &SharedState, active: Option<usize>) {
    let table = {
        let mut s = st.borrow_mut();
        s.active = active;
        refresh_all(ui, &mut s)
    };
    table.apply(ui);
    show_selected(ui, &st.borrow());
}

/// Last path component of an instrument file path (Windows `\` or Unix `/`).
fn file_base_name(path: &str) -> String {
    path.rsplit(['\\', '/']).next().unwrap_or(path).to_string()
}

fn save_as_dialog(ui: &AppWindow, st: &SharedState) {
    let blocking_error = {
        let s = st.borrow();
        save_as_blocking_error(&s)
    };
    if let Some(message) = blocking_error {
        set_status_error(ui, st, message);
        return;
    }

    let start_name = save_as_start_name(st.borrow().source_path().as_deref());

    if let Some(dst) = rfd::FileDialog::new()
        .add_filter("Bioanalyzer XML", SAVE_AS_XML_EXTENSIONS)
        .set_file_name(start_name)
        .save_file()
    {
        do_save(ui, st, dst);
    }
}

const SAVE_AS_XML_EXTENSIONS: &[&str] = &["xml"];

fn save_as_start_name(source: Option<&std::path::Path>) -> String {
    let Some(source) = source else {
        return "run.xml".to_string();
    };
    let Some(name) = source.file_name().and_then(|n| n.to_str()) else {
        return "run.xml".to_string();
    };
    for suffix in [".xml.gz", ".xml", ".xad"] {
        if let Some(stem) = name.strip_suffix(suffix) {
            return format!("{stem}.xml");
        }
    }
    source
        .file_stem()
        .and_then(|n| n.to_str())
        .map(|n| format!("{n}.xml"))
        .unwrap_or_else(|| "run.xml".to_string())
}

fn save_as_blocking_error(st: &AppState) -> Option<String> {
    if st.fragment_analyzer_mode() {
        Some(
            "Fragment Analyzer runs can only be saved in place by updating the .txt sidecar"
                .to_string(),
        )
    } else if st.has_marker_overrides() {
        Some(
            "Marker adjustments are session-only; reset markers before Save As, or export CSV/PNG with marker provenance instead."
                .to_string(),
        )
    } else {
        None
    }
}

fn export_peak_table_dialog(ui: &AppWindow, st: &SharedState) {
    if !can_export_peak_table(&st.borrow()) {
        set_status_error(ui, st, "No peak table is available to export".to_string());
        return;
    }

    let start_name = {
        let s = st.borrow();
        format!("{}_peaks.csv", export_stem(&s))
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_peak_table(&s) {
            Err(anyhow::anyhow!("no peak table is available"))
        } else {
            let sample = &s.run().samples[s.primary()];
            export::write_peak_table_csv(&dst, s.run(), sample, table_x_axis(&s))
        }
    };
    finish_export(ui, st, "Export peak table", &dst, result);
}

fn export_run_peak_table_dialog(ui: &AppWindow, st: &SharedState) {
    if !can_export_run_peak_table(&st.borrow()) {
        set_status_error(
            ui,
            st,
            "No processed samples are available to export".to_string(),
        );
        return;
    }

    let start_name = {
        let s = st.borrow();
        format!("{}_all_peaks.csv", export_run_stem(&s))
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_run_peak_table(&s) {
            Err(anyhow::anyhow!("no processed samples are available"))
        } else {
            export::write_run_peak_table_csv(&dst, s.run(), table_x_axis(&s))
        }
    };
    finish_export(ui, st, "Export all peak tables", &dst, result);
}

fn export_trace_data_dialog(ui: &AppWindow, st: &SharedState) {
    if !can_export_trace_data(&st.borrow()) {
        set_status_error(
            ui,
            st,
            "No processed trace data is available to export".to_string(),
        );
        return;
    }

    let start_name = {
        let s = st.borrow();
        format!("{}_traces.csv", export_stem(&s))
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_trace_data(&s) {
            Err(anyhow::anyhow!("no processed trace data is available"))
        } else {
            let samples: Vec<export::TraceSample<'_>> = if s.selection().len() > 1 {
                s.selection()
                    .iter()
                    .copied()
                    .filter(|&sample_index| s.sample_has_processed_data(sample_index))
                    .filter_map(|sample_index| {
                        s.run()
                            .samples
                            .get(sample_index)
                            .map(|sample| export::TraceSample {
                                sample_index,
                                sample,
                            })
                    })
                    .collect()
            } else {
                s.run()
                    .samples
                    .get(s.primary())
                    .map(|sample| export::TraceSample {
                        sample_index: s.primary(),
                        sample,
                    })
                    .into_iter()
                    .collect()
            };
            export::write_trace_data_csv(&dst, &samples, export::DEFAULT_TRACE_COLUMNS)
        }
    };
    finish_export(ui, st, "Export trace data", &dst, result);
}

fn export_raw_channels_dialog(ui: &AppWindow, st: &SharedState) {
    if !can_export_raw_channels(&st.borrow()) {
        set_status_error(
            ui,
            st,
            "No raw detector channels are available to export".to_string(),
        );
        return;
    }

    let start_name = {
        let s = st.borrow();
        format!("{}_raw_channels.csv", export_run_stem(&s))
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_raw_channels(&s) {
            Err(anyhow::anyhow!("no raw detector channels are available"))
        } else {
            let channels = raw_channels_for_export(&s);
            export::write_raw_channels_csv(&dst, &channels)
        }
    };
    finish_export(ui, st, "Export raw detector channels", &dst, result);
}

fn raw_channels_for_export(st: &AppState) -> Vec<&traceio::xad::RawChannel> {
    st.raw_channels().iter().collect()
}

fn export_metadata_qc_dialog(ui: &AppWindow, st: &SharedState) {
    if !can_export_metadata_qc(&st.borrow()) {
        set_status_error(ui, st, "No run metadata is available to export".to_string());
        return;
    }

    let start_name = {
        let s = st.borrow();
        format!("{}_metadata_qc.csv", export_run_stem(&s))
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_metadata_qc(&s) {
            Err(anyhow::anyhow!("no run metadata is available"))
        } else {
            let notes = metadata_export_notes(&s);
            export::write_metadata_qc_csv_with_notes_and_provenance(
                &dst,
                s.run(),
                Some(notes.as_str()),
                s.has_marker_overrides(),
            )
        }
    };
    finish_export(ui, st, "Export metadata/QC", &dst, result);
}

fn export_plot_png_dialog(ui: &AppWindow, st: &SharedState) {
    let Some(tab) = ExportPlotTab::from_tab_index(ui.get_active_tab()) else {
        set_status_error(ui, st, "No plot is available on this tab".to_string());
        return;
    };
    if !can_export_plot_for_tab(&st.borrow(), tab) {
        set_status_error(ui, st, "No plot is available to export".to_string());
        return;
    }

    let start_name = {
        let s = st.borrow();
        let suffix = match tab {
            ExportPlotTab::Detail => "plot",
            ExportPlotTab::Overview => "overview",
        };
        let stem = if tab == ExportPlotTab::Overview {
            export_run_stem(&s)
        } else {
            export_stem(&s)
        };
        format!("{stem}_{suffix}.png")
    };
    let Some(dst) = rfd::FileDialog::new()
        .add_filter("PNG image", &["png"])
        .set_file_name(start_name)
        .save_file()
    else {
        return;
    };

    let result = {
        let s = st.borrow();
        if !can_export_plot_for_tab(&s, tab) {
            Err(anyhow::anyhow!("no plot is available"))
        } else if tab == ExportPlotTab::Overview {
            match render_overview_bitmap(&s) {
                Some((buf, w, h)) => export::write_rgb_png(&dst, &buf, w, h),
                None => Err(anyhow::anyhow!("no overview is available")),
            }
        } else {
            let series = selected_series(&s);
            let refs: Vec<&Series> = series.iter().collect();
            let vp = s.viewport().unwrap_or_else(|| auto_fit_viewport(&s, &refs));
            let highlight = if s.selection().len() == 1 {
                s.highlight_x
            } else {
                None
            };
            let markers = marker_lines(&s);
            let buf =
                plot::render_overlay(&refs, &vp, highlight, &markers, plot::PLOT_W, plot::PLOT_H);
            export::write_rgb_png(&dst, &buf, plot::PLOT_W, plot::PLOT_H)
        }
    };
    let marker_note = {
        let s = st.borrow();
        match tab {
            ExportPlotTab::Detail => s.has_marker_overrides(),
            ExportPlotTab::Overview => overview_marker_overrides_active(&s),
        }
    };
    finish_export_with_marker_note(ui, st, "Export plot", &dst, result, marker_note);
}

fn can_export_peak_table(st: &AppState) -> bool {
    st.active_file().is_some()
        && !st.raw_mode()
        && st
            .run()
            .samples
            .get(st.primary())
            .is_some_and(|sample| sample_has_table_rows(st.run(), sample))
}

fn can_export_run_peak_table(st: &AppState) -> bool {
    st.active_file().is_some()
        && !st.raw_mode()
        && st
            .run()
            .samples
            .iter()
            .any(|sample| sample_has_table_rows(st.run(), sample))
}

fn can_export_trace_data(st: &AppState) -> bool {
    st.active_file().is_some()
        && !st.raw_mode()
        && st
            .selection()
            .iter()
            .any(|&sample_index| st.sample_has_processed_data(sample_index))
}

fn can_export_raw_channels(st: &AppState) -> bool {
    st.active_file().is_some() && st.raw_mode() && !st.raw_channels().is_empty()
}

fn can_export_metadata_qc(st: &AppState) -> bool {
    st.active_file().is_some()
}

fn can_export_plot(st: &AppState) -> bool {
    can_export_plot_for_tab(st, ExportPlotTab::Detail)
        || can_export_plot_for_tab(st, ExportPlotTab::Overview)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExportPlotTab {
    Detail,
    Overview,
}

impl ExportPlotTab {
    fn from_tab_index(idx: i32) -> Option<Self> {
        match idx {
            0 => Some(Self::Detail),
            1 => Some(Self::Overview),
            _ => None,
        }
    }
}

fn can_export_plot_for_tab(st: &AppState, tab: ExportPlotTab) -> bool {
    if st.active_file().is_none() {
        return false;
    }
    match tab {
        ExportPlotTab::Detail => {
            st.primary() < st.entry_count()
                && (st.raw_mode() || st.sample_has_processed_data(st.primary()))
        }
        ExportPlotTab::Overview => render_overview_bitmap(st).is_some(),
    }
}

fn sample_has_table_rows(run: &traceio::Electrophoresis, sample: &traceio::Sample) -> bool {
    !sample.peaks.is_empty() || !sample.regions.is_empty() || !run.regions.is_empty()
}

fn finish_export(
    ui: &AppWindow,
    st: &SharedState,
    operation: &str,
    dst: &std::path::Path,
    result: anyhow::Result<()>,
) {
    let marker_note = st.borrow().has_marker_overrides();
    finish_export_with_marker_note(ui, st, operation, dst, result, marker_note);
}

fn finish_export_with_marker_note(
    ui: &AppWindow,
    st: &SharedState,
    operation: &str,
    dst: &std::path::Path,
    result: anyhow::Result<()>,
    marker_note: bool,
) {
    match result {
        Ok(()) => {
            let message = format!(
                "{}: {}{}",
                operation.replacen("Export", "Exported", 1),
                dst.display(),
                if marker_note {
                    " (session marker overrides active)"
                } else {
                    ""
                }
            );
            {
                let mut s = st.borrow_mut();
                s.error = None;
                s.status = Some(message.clone());
            }
            ui.set_error_text(SharedString::default());
            ui.set_status_text(SharedString::from(message));
        }
        Err(e) => set_status_error(ui, st, format!("{operation} failed: {e:#}")),
    }
}

fn metadata_export_notes(st: &AppState) -> String {
    let mut notes = String::new();
    if st.has_marker_overrides() {
        notes.push_str(
            "Session marker overrides are active; calibrated lengths, concentrations, and molarity values were recomputed in this GUI session and are not persisted in the source file.\n",
        );
    }
    if !st.metadata_text().trim().is_empty() {
        notes.push_str(st.metadata_text());
    }
    notes
}

fn set_status_error(ui: &AppWindow, st: &SharedState, message: String) {
    {
        let mut s = st.borrow_mut();
        set_status_error_state(&mut s, message.clone());
    }
    ui.set_error_text(SharedString::from(message));
    ui.set_status_text(SharedString::default());
}

fn set_status_error_state(st: &mut AppState, message: String) {
    st.error = Some(message);
    st.status = None;
}

fn export_stem(st: &AppState) -> String {
    let mut parts = Vec::new();
    parts.push(export_run_stem(st));

    if st.selection().len() > 1 {
        parts.push("overlay".to_string());
    } else if st.raw_mode() {
        if let Some(ch) = st.raw_channels().get(st.primary()) {
            parts.push(ch.channel_id.to_string());
        }
    } else if let Some(sample) = st.run().samples.get(st.primary()) {
        let label = if sample.name.trim().is_empty() {
            format!("well_{}", sample.well_number)
        } else {
            format!("well_{}_{}", sample.well_number, sample.name.trim())
        };
        parts.push(label);
    }

    sanitize_export_stem(&parts.join("_"))
}

fn export_run_stem(st: &AppState) -> String {
    let run_part = st
        .source_path()
        .as_ref()
        .and_then(|p| p.file_name())
        .and_then(|n| n.to_str())
        .map(|n| strip_known_extension(n).to_string())
        .or_else(|| {
            let name = file_base_name(&st.run().assay.file_name);
            (!name.is_empty()).then_some(strip_known_extension(&name).to_string())
        })
        .unwrap_or_else(|| "run".to_string());
    sanitize_export_stem(&run_part)
}

fn strip_known_extension(name: &str) -> &str {
    for suffix in [
        ".xml.gz", ".csv.gz", ".fa.zip", ".xml", ".csv", ".xad", ".zip", ".raw",
    ] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            return stripped;
        }
    }
    name
}

fn sanitize_export_stem(raw: &str) -> String {
    let sanitized: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let trimmed = sanitized.trim_matches(['_', '.']);
    if trimmed.is_empty() {
        "traceanalyzer".to_string()
    } else {
        trimmed.chars().take(96).collect()
    }
}

/// Write the current run to `dst`, using the loaded source file as the template,
/// then adopt `dst` as the new source and clear the dirty flag.
fn do_save(ui: &AppWindow, st: &SharedState, dst: std::path::PathBuf) {
    struct SavedRun {
        loaded: loading::Loaded,
        metadata: Option<LoadedFaMetadata>,
    }

    let result = (|| -> anyhow::Result<SavedRun> {
        let s = st.borrow();
        match s.source_path() {
            Some(src) => {
                if save_destination_conflict(&s, &dst).is_some() {
                    anyhow::bail!(
                        "destination {} is already open; close that file before saving over it",
                        dst.display()
                    );
                }
                traceio::save::save_run(s.run(), &src, &dst)?;
                let loaded = loading::load(&dst)?;
                let metadata = loaded_fa_metadata(&dst, loaded.fa_metadata.clone());
                Ok(SavedRun { loaded, metadata })
            }
            None => Err(anyhow::anyhow!("no source file to save from")),
        }
    })();
    match result {
        Ok(saved) => {
            let marker_note = st.borrow().has_marker_overrides();
            let table = {
                let mut s = st.borrow_mut();
                let source = saved.loaded.source.identity.clone();
                let source_capabilities = saved.loaded.source.capabilities;
                let save_capabilities = saved.loaded.source.save_capabilities();
                s.replace_active_file_with_capabilities(
                    saved.loaded.run,
                    saved.loaded.raw_channels,
                    Some(source),
                    source_capabilities,
                    save_capabilities,
                );
                s.set_active_opened_path(Some(dst.clone()));
                let metadata_warning = saved.metadata.and_then(|metadata| {
                    s.set_active_metadata_text(metadata.text);
                    s.set_active_metadata_rows(metadata.rows);
                    metadata.warning
                });
                s.error = join_optional_messages([saved.loaded.warning, metadata_warning]);
                if s.has_marker_overrides() {
                    if let Err(e) = commit_recalibration(&mut s) {
                        let recalibration_warning = Some(format!(
                            "Could not reapply marker overrides after save: {e:#}"
                        ));
                        s.error = join_optional_messages([s.error.take(), recalibration_warning]);
                    }
                }
                s.status = Some(if marker_note {
                    "Saved sample-name changes; marker adjustments remain session-only".to_string()
                } else {
                    "Saved changes".to_string()
                });
                refresh_all(ui, &mut s)
            };
            table.apply(ui);
            let s = st.borrow();
            ui.set_can_save(s.is_dirty() && s.can_save());
            ui.set_can_save_as(s.can_save_as());
            ui.set_error_text(SharedString::from(s.error.clone().unwrap_or_default()));
            ui.set_status_text(SharedString::from(s.status.clone().unwrap_or_default()));
            refresh_tree(ui, &s);
            show_selected(ui, &s);
        }
        Err(e) => {
            let message = format!("Save failed: {e:#}");
            {
                let mut s = st.borrow_mut();
                s.error = Some(message.clone());
                s.status = None;
            }
            ui.set_error_text(SharedString::from(message));
            ui.set_status_text(SharedString::default());
        }
    }
}

fn save_destination_conflict(st: &AppState, dst: &std::path::Path) -> Option<usize> {
    st.find_file_by_source(dst)
        .filter(|&existing| Some(existing) != st.active)
}

fn rgb_to_image(buf: &[u8], w: u32, h: u32) -> Image {
    let mut pixbuf = SharedPixelBuffer::<Rgb8Pixel>::new(w, h);
    pixbuf.make_mut_bytes().copy_from_slice(buf);
    Image::from_rgb8(pixbuf)
}

/// A plot-sized image with just empty axes, shown when no file is open.
fn blank_plot_image() -> Image {
    let (buf, w, h) = blank_plot_bitmap();
    rgb_to_image(&buf, w, h)
}

fn blank_plot_bitmap() -> (Vec<u8>, u32, u32) {
    let refs: Vec<&Series> = Vec::new();
    let vp = plot::auto_viewport_multi(&refs);
    let buf = plot::render_overlay(&refs, &vp, None, &[], plot::PLOT_W, plot::PLOT_H);
    (buf, plot::PLOT_W, plot::PLOT_H)
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
        if !st.can_edit_markers() {
            return;
        }
        let idx = drag.sample_idx;
        if idx >= st.run().samples.len() || st.raw_mode() {
            return;
        }
        let vp = current_viewport(&st);
        let span = vp.x_max - vp.x_min;
        let previous = st.overrides().get(&idx).copied();
        let (lo, up) = marker_times(st.run(), idx, previous.as_ref());
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
    refresh_tree(ui, &state.borrow());
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

#[cfg(test)]
mod tests {
    use super::*;
    use traceio::{AssayInfo, Sample};

    fn sample(well_number: i32, is_ladder: bool) -> Sample {
        Sample {
            well_number,
            name: format!("S{well_number}"),
            category: String::new(),
            is_ladder,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time: vec![0.0, 1.0],
            fluorescence: vec![1.0, 2.0],
            aligned_time: Vec::new(),
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
            regions: Vec::new(),
        }
    }

    fn metadata_only_sample(well_number: i32) -> Sample {
        Sample {
            time: Vec::new(),
            fluorescence: Vec::new(),
            ..sample(well_number, false)
        }
    }

    fn run(name: &str, samples: Vec<Sample>) -> Electrophoresis {
        Electrophoresis {
            assay: AssayInfo {
                file_name: name.to_string(),
                ..Default::default()
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples,
        }
    }

    #[test]
    fn overview_export_is_unavailable_when_no_samples_are_visible() {
        let mut st = AppState::new(
            run("ladders.xml", vec![sample(1, true)]),
            Vec::new(),
            None,
            None,
        );

        assert!(!can_export_plot_for_tab(&st, ExportPlotTab::Overview));

        st.overview_show_ladders = true;
        assert!(can_export_plot_for_tab(&st, ExportPlotTab::Overview));
    }

    #[test]
    fn overview_export_is_unavailable_for_metadata_only_samples() {
        let mut st = AppState::new(
            run("metadata.xml", vec![metadata_only_sample(1)]),
            Vec::new(),
            None,
            None,
        );

        assert!(!can_export_plot_for_tab(&st, ExportPlotTab::Overview));

        st.overview_gel = true;
        assert!(!can_export_plot_for_tab(&st, ExportPlotTab::Overview));
    }

    #[test]
    fn overview_export_remains_available_when_any_sample_has_trace_data() {
        let mut st = AppState::new(
            run("mixed.xml", vec![metadata_only_sample(1), sample(2, false)]),
            Vec::new(),
            None,
            None,
        );

        assert!(can_export_plot_for_tab(&st, ExportPlotTab::Overview));

        st.overview_gel = true;
        assert!(can_export_plot_for_tab(&st, ExportPlotTab::Overview));
    }

    #[test]
    fn multi_file_overview_labels_use_opened_file_labels() {
        let mut st = AppState::new(
            run("duplicate-assay.xml", vec![sample(1, false)]),
            Vec::new(),
            Some(PathBuf::from("/canonical/a.raw")),
            None,
        );
        st.set_active_opened_path(Some(PathBuf::from("/opened/first.fa.zip")));
        st.add_file(
            run("duplicate-assay.xml", vec![sample(1, false)]),
            Vec::new(),
            Some(PathBuf::from("/canonical/b.raw")),
        );
        st.set_active_opened_path(Some(PathBuf::from("/opened/second.fa.zip")));

        let (combined, _) = overview_combined(&st);

        assert_eq!(combined.samples[0].name, "first.fa.zip · S1");
        assert_eq!(combined.samples[1].name, "second.fa.zip · S1");
    }

    #[test]
    fn overview_marker_provenance_checks_all_contributing_files() {
        let mut st = AppState::new(run("a.xml", vec![sample(1, false)]), Vec::new(), None, None);
        st.add_file(run("b.xml", vec![sample(1, false)]), Vec::new(), None);
        st.active = Some(0);
        st.files[1].overrides.insert(
            0,
            MarkerOverride {
                lower_time: Some(0.2),
                upper_time: None,
            },
        );

        assert!(!st.has_marker_overrides());
        assert!(overview_marker_overrides_active(&st));
    }

    #[test]
    fn save_destination_conflict_rejects_other_open_file() {
        let mut st = AppState::new(
            run("a.xml", vec![sample(1, false)]),
            Vec::new(),
            Some(PathBuf::from("/tmp/a.xml")),
            None,
        );
        st.add_file(
            run("b.xml", vec![sample(1, false)]),
            Vec::new(),
            Some(PathBuf::from("/tmp/b.xml")),
        );
        st.active = Some(0);

        assert_eq!(
            save_destination_conflict(&st, PathBuf::from("/tmp/b.xml").as_path()),
            Some(1)
        );
        assert_eq!(
            save_destination_conflict(&st, PathBuf::from("/tmp/a.xml").as_path()),
            None
        );
    }

    #[test]
    fn overview_marker_provenance_ignores_hidden_sample_overrides() {
        let mut st = AppState::new(
            run("run.xml", vec![sample(1, true), sample(2, false)]),
            Vec::new(),
            None,
            None,
        );
        st.overview_show_ladders = false;
        st.overrides_mut().insert(
            0,
            MarkerOverride {
                lower_time: Some(0.2),
                upper_time: None,
            },
        );

        assert!(!overview_marker_overrides_active(&st));

        st.overrides_mut().insert(
            1,
            MarkerOverride {
                lower_time: Some(0.3),
                upper_time: None,
            },
        );
        assert!(overview_marker_overrides_active(&st));
    }

    #[test]
    fn save_as_gui_filter_only_offers_plain_xml_extension() {
        assert_eq!(SAVE_AS_XML_EXTENSIONS, &["xml"]);
    }

    #[test]
    fn save_as_start_name_does_not_double_xml_gz_source_suffix() {
        assert_eq!(
            save_as_start_name(Some(PathBuf::from("run.xml.gz").as_path())),
            "run.xml"
        );
        assert_eq!(
            save_as_start_name(Some(PathBuf::from("run.xad").as_path())),
            "run.xml"
        );
    }

    #[test]
    fn save_as_fragment_analyzer_reports_state_backed_error() {
        let mut run = run("run.raw", vec![sample(1, false)]);
        run.assay.assay_name = "Fragment Analyzer".to_string();
        let mut st = AppState::new(run, Vec::new(), Some(PathBuf::from("run.raw")), None);
        st.status = Some("Ready".to_string());
        let message =
            "Fragment Analyzer runs can only be saved in place by updating the .txt sidecar"
                .to_string();

        assert_eq!(save_as_blocking_error(&st), Some(message.clone()));

        set_status_error_state(&mut st, message.clone());
        assert_eq!(st.error, Some(message));
        assert_eq!(st.status, None);
    }

    #[test]
    fn save_as_marker_overrides_report_state_backed_error() {
        let mut st = AppState::new(
            run("run.xml", vec![sample(1, false)]),
            Vec::new(),
            Some(PathBuf::from("run.xml")),
            None,
        );
        st.overrides_mut().insert(
            0,
            MarkerOverride {
                lower_time: Some(0.2),
                upper_time: None,
            },
        );

        assert_eq!(
            save_as_blocking_error(&st),
            Some(
                "Marker adjustments are session-only; reset markers before Save As, or export CSV/PNG with marker provenance instead."
                    .to_string()
            )
        );
    }

    #[test]
    fn export_plot_is_unavailable_on_metadata_tab() {
        assert_eq!(ExportPlotTab::from_tab_index(2), None);
    }
}
