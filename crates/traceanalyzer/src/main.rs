//! Minimal Slint viewer for a Bioanalyzer run: sample list + electropherogram.
//!
//! Usage: cargo run -p traceanalyzer -- <file.xad | file.xml | file.xml.gz>
//! With no argument it loads the bundled DNA 1000 demo, if present.

use std::io::Read;
use std::rc::Rc;

use slint::{ModelRc, SharedString, VecModel};
use traceio::{Electrophoresis, Sample};

slint::include_modules!();

// Plot viewbox (must match ui/app.slint).
const VB_W: f64 = 1000.0;
const VB_H: f64 = 400.0;
const MARGIN: f64 = 12.0;

fn main() -> anyhow::Result<()> {
    let mut run = match std::env::args().nth(1) {
        Some(path) => load(&path)?,
        None => {
            let demo = concat!(env!("CARGO_MANIFEST_DIR"), "/../../testdata/demo_dna1000.xml.gz");
            if std::path::Path::new(demo).exists() {
                load(demo)?
            } else {
                anyhow::bail!("no file given and demo fixture not found");
            }
        }
    };
    // Per-point ladder calibration so the x-axis can be shown in bp/nt.
    if let Err(e) = traceio::calibration::calculate_length(
        &mut run,
        traceio::calibration::Method::Hyman,
    ) {
        eprintln!("(calibration skipped: {e})");
    }
    let run = Rc::new(run);

    let ui = AppWindow::new()?;
    ui.set_assay_title(SharedString::from(format!(
        "{}  —  {} ({}),  {} samples",
        run.assay.file_name, run.assay.assay_name, run.assay.assay_type, run.samples.len()
    )));

    let names: Vec<SharedString> = run
        .samples
        .iter()
        .map(|s| {
            let label = if s.name.is_empty() {
                format!("Well {}", s.well_number)
            } else {
                format!("{}: {}", s.well_number, s.name)
            };
            SharedString::from(if s.is_ladder {
                format!("{label}  [ladder]")
            } else {
                label
            })
        })
        .collect();
    ui.set_sample_names(ModelRc::from(Rc::new(VecModel::from(names))));

    // Wire selection -> recompute plot.
    {
        let ui_weak = ui.as_weak();
        let run = run.clone();
        ui.on_select(move |idx| {
            if let Some(ui) = ui_weak.upgrade() {
                show_sample(&ui, &run, idx as usize);
            }
        });
    }

    // Initial view.
    if !run.samples.is_empty() {
        show_sample(&ui, &run, 0);
    }

    ui.run()?;
    Ok(())
}

fn show_sample(ui: &AppWindow, run: &Electrophoresis, idx: usize) {
    let Some(sample) = run.samples.get(idx) else {
        return;
    };
    ui.set_current_index(idx as i32);
    let axis = XAxis::for_sample(sample, run);
    ui.set_trace_commands(SharedString::from(trace_path(sample, &axis)));
    ui.set_peak_commands(SharedString::from(peak_path(sample, &axis)));
    ui.set_sample_info(SharedString::from(sample_info(run, sample, &axis)));
}

/// Which quantity to plot along x: calibrated length if available, else time.
struct XAxis {
    use_length: bool,
    min: f64,
    max: f64,
    label: String,
}

impl XAxis {
    fn for_sample(s: &Sample, run: &Electrophoresis) -> XAxis {
        // Prefer calibrated length (bp/nt) when the sample has finite values.
        let finite_len: Vec<f64> = s.length.iter().copied().filter(|v| v.is_finite()).collect();
        if finite_len.len() >= 2 {
            let min = finite_len.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = finite_len.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            XAxis {
                use_length: true,
                min,
                max,
                label: format!("size ({})", run.assay.length_unit),
            }
        } else {
            let (min, max) = (
                *s.time.first().unwrap_or(&0.0),
                *s.time.last().unwrap_or(&1.0),
            );
            XAxis {
                use_length: false,
                min,
                max,
                label: "migration time (s)".to_string(),
            }
        }
    }

    /// x-value for trace point `i`, or None if it should be skipped.
    fn point_x(&self, s: &Sample, i: usize) -> Option<f64> {
        let v = if self.use_length {
            *s.length.get(i)?
        } else {
            *s.time.get(i)?
        };
        v.is_finite().then_some(v)
    }

    fn map_x(&self, v: f64) -> f64 {
        let span = (self.max - self.min).max(f64::EPSILON);
        MARGIN + (v - self.min) / span * (VB_W - 2.0 * MARGIN)
    }
}

/// Build an SVG polyline for the fluorescence trace in the plot viewbox.
fn trace_path(s: &Sample, axis: &XAxis) -> String {
    if s.fluorescence.is_empty() {
        return String::new();
    }
    let (lo, hi) = s.fluorescence_range().unwrap_or((0.0, 1.0));
    let y_span = ((hi - lo) as f64).max(f64::EPSILON);
    // flip y: higher fluorescence -> smaller y (towards top)
    let map_y = |f: f32| VB_H - MARGIN - (f as f64 - lo as f64) / y_span * (VB_H - 2.0 * MARGIN);

    // Downsample to keep the path string small.
    let n = s.fluorescence.len();
    let step = (n / 1500).max(1);

    let mut cmd = String::with_capacity(n / step * 16);
    let mut started = false;
    let mut i = 0;
    while i < n {
        if let Some(xv) = axis.point_x(s, i) {
            let x = axis.map_x(xv);
            let y = map_y(s.fluorescence[i]);
            if !started {
                cmd.push_str(&format!("M {:.1} {:.1} ", x, y));
                started = true;
            } else {
                cmd.push_str(&format!("L {:.1} {:.1} ", x, y));
            }
        }
        i += step;
    }
    cmd
}

/// Vertical ticks at the bottom marking called peak positions.
fn peak_path(s: &Sample, axis: &XAxis) -> String {
    if s.peaks.is_empty() {
        return String::new();
    }
    let mut cmd = String::new();
    for p in &s.peaks {
        let v = if axis.use_length { p.length } else { p.time };
        if !v.is_finite() || v < axis.min || v > axis.max {
            continue;
        }
        let x = axis.map_x(v);
        cmd.push_str(&format!("M {:.1} {:.1} L {:.1} {:.1} ", x, VB_H - MARGIN, x, VB_H - MARGIN - 24.0));
    }
    cmd
}

fn sample_info(run: &Electrophoresis, s: &Sample, axis: &XAxis) -> String {
    let mut parts = vec![
        format!("Well {}", s.well_number),
        format!("{} peaks", s.peaks.len()),
    ];
    if let Some(rin) = s.rin {
        parts.push(format!("RIN {rin:.1}"));
    }
    if axis.use_length {
        parts.push(format!("{:.0}–{:.0} {}", axis.min, axis.max, run.assay.length_unit));
    }
    if !s.observations.is_empty() {
        parts.push(s.observations.clone());
    }
    parts.push(format!("x-axis: {};  y: fluorescence", axis.label));
    parts.join("   ·   ")
}

fn load(path: &str) -> anyhow::Result<Electrophoresis> {
    let p = std::path::Path::new(path);
    if path.ends_with(".xad") {
        Ok(traceio::xad::read_xad_file(p)?)
    } else if path.ends_with(".xml.gz") {
        let raw = std::fs::read(p)?;
        let mut d = flate2::read::GzDecoder::new(&raw[..]);
        let mut s = String::new();
        d.read_to_string(&mut s)?;
        Ok(traceio::bioanalyzer::parse_xml(&s)?)
    } else {
        let s = std::fs::read_to_string(p)?;
        Ok(traceio::bioanalyzer::parse_xml(&s)?)
    }
}
