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

use traceanalyzer::plot::{self, Series, Viewport};
use traceanalyzer::state::{AppState, Marker};
use traceanalyzer::{gel, loading, overview, render, table};
use traceio::calibration::marker_times;

slint::include_modules!();

type SharedState = Rc<RefCell<AppState>>;

fn main() -> anyhow::Result<()> {
    // Select the winit backend explicitly so the external file-drop hook below
    // (registered via `on_winit_window_event`) is available.
    if let Err(e) = slint::BackendSelector::new().backend_name("winit".into()).select() {
        eprintln!("(winit backend unavailable, drag-drop disabled: {e})");
    }

    let path = match std::env::args().nth(1) {
        Some(p) => PathBuf::from(p),
        None => {
            let demo =
                PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../../testdata/demo_dna1000.xml.gz"));
            if !demo.exists() {
                anyhow::bail!("no file given and demo fixture not found");
            }
            demo
        }
    };
    let loaded = loading::load(&path)?;
    let state: SharedState = Rc::new(RefCell::new(AppState::new(loaded.run, loaded.raw_channels)));

    let ui = AppWindow::new()?;
    refresh_all(&ui, &mut state.borrow_mut());

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

    // External file drag-and-drop: reload on the last dropped file. Falls back
    // silently to the Open… dialog when the winit backend isn't active.
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
            {
                let mut s = st.borrow_mut();
                if s.select_click(idx as usize, ctrl, shift) {
                    s.viewport = None; // auto-fit the new selection
                }
                refresh_table(&ui, &mut s);
            }
            refresh_selection(&ui, &st.borrow());
            show_selected(&ui, &st.borrow());
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
                s.highlight_x = usize::try_from(row).ok().and_then(|r| s.table_peak_x.get(r).copied().flatten());
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
                if s.overview_gel {
                    let (w, _) = gel::size(s.run.samples.len());
                    gel::lane_at(&s.run, fx as f64, w)
                } else {
                    let layout = overview::layout(s.entry_count());
                    overview::cell_at(&layout, fx as f64, fy as f64)
                }
            };
            if let Some(idx) = idx {
                {
                    let mut s = st.borrow_mut();
                    s.select_click(idx, false, false);
                    s.viewport = None;
                    refresh_table(&ui, &mut s);
                }
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
            if let Some(marker) = grabbed {
                drag_marker(&ui, &st, marker, dfx as f64);
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
            s.grabbed = best.map(|b| b.0);
        });
    }
    // Marker-edit: pointer up releases any grabbed marker (and refreshes views).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_plot_release(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            let was = st.borrow().grabbed.is_some();
            st.borrow_mut().grabbed = None;
            if was {
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
            {
                let mut s = st.borrow_mut();
                s.marker_edit ^= true;
                s.viewport = None; // x-axis space changes
                s.grabbed = None;
            }
            show_selected(&ui, &st.borrow());
        });
    }
    // Reset the focused sample's markers to automatic detection.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_reset_markers(move || {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                let idx = s.primary();
                s.overrides.remove(&idx);
                let ov = s.overrides.clone();
                loading::recalibrate_with(&mut s.run, &ov);
                s.viewport = None;
                refresh_table(&ui, &mut s);
            }
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
            {
                let mut s = state.borrow_mut();
                *s = AppState::new(loaded.run, loaded.raw_channels);
                refresh_all(ui, &mut s);
            }
            ui.set_error_text(SharedString::new());
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

/// Push run-level data (title, entry list, error, table) into the UI.
fn refresh_all(ui: &AppWindow, st: &mut AppState) {
    ui.set_assay_title(SharedString::from(st.title()));
    let names: Vec<SharedString> = st.entry_labels().into_iter().map(SharedString::from).collect();
    ui.set_sample_names(ModelRc::from(Rc::new(VecModel::from(names))));
    ui.set_error_text(SharedString::from(st.error.clone().unwrap_or_default()));
    refresh_selection(ui, st);
    refresh_overview(ui, st);
    refresh_table(ui, st);
}

/// Render the Overview tab (small-multiples grid or virtual gel) into the UI.
fn refresh_overview(ui: &AppWindow, st: &AppState) {
    ui.set_overview_shared_y(st.overview_shared_y);
    ui.set_overview_gel(st.overview_gel);
    if st.overview_gel {
        let (w, h) = gel::size(st.run.samples.len());
        let buf = gel::render(&st.run, w, h);
        ui.set_overview_image(rgb_to_image(&buf, w, h));
    } else {
        let layout = overview::layout(st.run.samples.len());
        let buf = overview::render(&st.run, st.y_mode, st.overview_shared_y, &layout);
        ui.set_overview_image(rgb_to_image(&buf, layout.w, layout.h));
    }
}

/// Build the peak/region table for the focused sample and push it to the UI.
/// Also records each row's plot x-position (for cross-highlighting) into `st`.
fn refresh_table(ui: &AppWindow, st: &mut AppState) {
    if st.raw_mode() || st.run.samples.get(st.primary()).is_none() {
        st.table_peak_x.clear();
        ui.set_table_rows(ModelRc::from(Rc::new(VecModel::<ModelRc<StandardListViewItem>>::default())));
        ui.set_table_current_row(-1);
        return;
    }
    let rows = table::rows(&st.run, &st.run.samples[st.primary()]);
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
    ui.set_table_rows(ModelRc::from(Rc::new(VecModel::from(model_rows))));
    ui.set_table_current_row(-1);
}

/// Push the current selection (primary index + per-row flags) into the UI.
fn refresh_selection(ui: &AppWindow, st: &AppState) {
    ui.set_current_index(st.primary() as i32);
    let flags: Vec<bool> = st.selection_flags();
    ui.set_selected_flags(ModelRc::from(Rc::new(VecModel::from(flags))));
}

/// Build the plot series for the current selection (raw channel, single sample,
/// or several overlaid samples), applying normalization when requested.
fn selected_series(st: &AppState) -> Vec<Series> {
    if st.raw_mode() {
        return vec![plot::raw_series(&st.raw_channels[st.primary()])];
    }
    let overlay = st.selection.len() > 1;
    // Marker-edit mode forces the raw-time x-axis (single sample only).
    let force_time = st.marker_edit && !overlay;
    st.selection
        .iter()
        .filter_map(|&i| st.run.samples.get(i))
        .map(|s| {
            let series = plot::series(&st.run, s, st.y_mode, force_time);
            if overlay && st.normalize {
                plot::normalized(&series)
            } else {
                series
            }
        })
        .collect()
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
    let vp = st.viewport.unwrap_or_else(|| plot::auto_viewport_multi(&refs));
    // Highlight the selected peak only in the single-sample view.
    let highlight = if st.selection.len() == 1 { st.highlight_x } else { None };
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

/// Move a grabbed marker by a fractional-x drag delta: update the override raw
/// time, re-run sizing, and refresh the plot + table live. The viewport (in
/// raw-time space) is kept fixed so the trace does not jump while dragging.
fn drag_marker(ui: &AppWindow, state: &SharedState, marker: Marker, dfx: f64) {
    let mut st = state.borrow_mut();
    if st.primary() >= st.entry_count() {
        return;
    }
    let idx = st.primary();
    let vp = current_viewport(&st);
    let span = vp.x_max - vp.x_min;
    let ov = st.overrides.get(&idx).copied();
    let (lo, up) = marker_times(&st.run, idx, ov.as_ref());
    let cur = match marker {
        Marker::Lower => lo,
        Marker::Upper => up,
    };
    let Some(cur) = cur else { return };
    let new_time = cur + dfx * span;

    let entry = st.overrides.entry(idx).or_default();
    match marker {
        Marker::Lower => entry.lower_time = Some(new_time),
        Marker::Upper => entry.upper_time = Some(new_time),
    }
    let ov = st.overrides.clone();
    loading::recalibrate_with(&mut st.run, &ov);
    refresh_table(ui, &mut st);
    drop(st);
    show_selected(ui, &state.borrow());
}
