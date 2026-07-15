//! Slint viewer for automated-electrophoresis runs.
//!
//! Usage: cargo run -p traceanalyzer -- <file.xad | file.xml | file.xml.gz>
//! With no argument it loads the bundled DNA 1000 demo, if present.

use std::cell::RefCell;
use std::path::PathBuf;
use std::rc::Rc;

use slint::{Image, ModelRc, Rgb8Pixel, SharedPixelBuffer, SharedString, VecModel};

use traceanalyzer::plot::{self, Series, Viewport};
use traceanalyzer::state::AppState;
use traceanalyzer::{loading, render};

slint::include_modules!();

type SharedState = Rc<RefCell<AppState>>;

fn main() -> anyhow::Result<()> {
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
    refresh_all(&ui, &state.borrow());

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

    // Selection -> reset viewport to auto-fit and re-render.
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_select(move |idx| {
            let Some(ui) = ui_weak.upgrade() else { return };
            {
                let mut s = st.borrow_mut();
                s.selected = idx as usize;
                s.viewport = None; // auto-fit new sample
            }
            show_selected(&ui, &st.borrow());
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
    // Pan by a fractional delta (drag).
    {
        let ui_weak = ui.as_weak();
        let st = state.clone();
        ui.on_pan(move |dfx, dfy| {
            let Some(ui) = ui_weak.upgrade() else { return };
            pan(&st, dfx as f64, dfy as f64);
            show_selected(&ui, &st.borrow());
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
            *state.borrow_mut() = AppState::new(loaded.run, loaded.raw_channels);
            ui.set_error_text(SharedString::new());
            refresh_all(ui, &state.borrow());
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

/// Push run-level data (title, entry list, error) into the UI.
fn refresh_all(ui: &AppWindow, st: &AppState) {
    ui.set_assay_title(SharedString::from(st.title()));
    let names: Vec<SharedString> = st.entry_labels().into_iter().map(SharedString::from).collect();
    ui.set_sample_names(ModelRc::from(Rc::new(VecModel::from(names))));
    ui.set_error_text(SharedString::from(st.error.clone().unwrap_or_default()));
}

/// Render the currently selected entry (sample or raw channel) into the plot.
fn show_selected(ui: &AppWindow, st: &AppState) {
    if st.selected >= st.entry_count() {
        return;
    }
    ui.set_current_index(st.selected as i32);

    let (series, info) = if st.raw_mode() {
        let ch = &st.raw_channels[st.selected];
        (plot::raw_series(ch), render::raw_info_line(ch))
    } else {
        let sample = &st.run.samples[st.selected];
        (
            plot::series(&st.run, sample, st.y_mode),
            render::info_line(&st.run, sample),
        )
    };

    let vp = st.viewport.unwrap_or_else(|| plot::auto_viewport(&series));
    let buf = plot::render_rgb(&series, &vp, plot::PLOT_W, plot::PLOT_H);
    ui.set_plot_image(rgb_to_image(&buf, plot::PLOT_W, plot::PLOT_H));
    ui.set_sample_info(SharedString::from(info));
    ui.set_y_mode_label(SharedString::from(st.y_mode.label(&st.run)));
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
    let series: Series = if st.raw_mode() {
        plot::raw_series(&st.raw_channels[st.selected])
    } else {
        plot::series(&st.run, &st.run.samples[st.selected], st.y_mode)
    };
    plot::auto_viewport(&series)
}

fn zoom(state: &SharedState, fx: f64, fy: f64, factor: f64) {
    let mut st = state.borrow_mut();
    if st.selected >= st.entry_count() {
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
    if st.selected >= st.entry_count() {
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
