# Trace analyzer GUI — implementation plan

Living checklist for building out the Slint viewer. Tick items as they land.
Legend: `[ ]` todo · `[~]` in progress · `[x]` done.

## Decisions
- **Rendering:** hybrid — trace/axes/gel rendered to a bitmap with `plotters`,
  shown in a Slint `Image`; interactive elements (marker lines, crosshair,
  selection, hover) drawn as Slint elements on top. One shared pixel↔data
  transform, defined in Phase A, used by both.
- **File dialog:** `rfd` (native picker).
- Format-reader (`traceio`) changes are confined to Phase H (manual markers);
  all other phases are GUI-only (`src`).
- **Headless/CI:** bitmap render tests and examples run without a display. GUI
  smoke runs need a real display or Xvfb (`xvfb-run -a cargo run --
  testdata/demo_dna1000.xml.gz`). Fetch gitignored fixtures
  first with `bash scripts/fetch-testdata.sh` for parser/file-loading coverage;
  the headless render smoke test does not need downloaded fixtures.
- **System packages:** Slint/winit/rfd Linux builds may require X11/Wayland and
  GTK development packages, depending on the distribution image.

## New dependencies (added as phases need them)
- [x] `plotters` (bitmap + ab_glyph, vendored DejaVuSans font) — Phase B
- [x] `rfd` — Phase C
- [x] `slint` feature `unstable-winit-030` (file-drop hook) — Phase C

## Packaging/test coverage
- [x] Repository collapsed to one package; `traceio` now lives under
      `src/traceio`, so `cargo package` has no internal path-crate registry
      dependency.
- [x] Native `.xad` loading surfaces raw-channel decode errors in the GUI instead
      of silently treating them as no raw channels.
- [x] Headless render coverage asserts dimensions and nonblank plot output.

## Feature → phase map
| Requested feature | Phase |
|---|---|
| File-open dialog | C |
| Drag-drop of files | C |
| Drawn axes (ticks/gridlines/labels) | B |
| Zoom / pan | B |
| Concentration/molarity display | B + F |
| Native `.xad` raw-channel view | C |
| Multi-select + overlapped traces | D |
| Overview tab (all traces) | E |
| Peak / region table | F |
| Gel / virtual-gel view | G |
| Manual low/high markers | H |

---

## Phase A — Foundation refactor
- [x] `AppState` struct (run, raw channels, selection set, marker overrides,
      viewport, axis/y-mode, colors) behind `Rc<RefCell<…>>`.
- [x] Split GUI into modules: `state.rs`, `render.rs`, UI wiring in `main.rs`.
- [x] Slint shell with `TabWidget`: **Detail** | **Overview**. Detail contains
      the sample list, plot, and peak table; Overview is the all-sample view.
- [x] Run full pipeline on load (`calculate_length` → `calculate_concentration`
      → `calculate_molarity`).
- [x] Define the shared pixel↔data transform used by renderer + overlays.
- *Deliverable:* today's behavior, restructured; no user-visible change.

## Phase B — Real plot (axes, zoom/pan, y-mode)
- [x] `render.rs`: plotters bitmap of the trace with "nice" ticks, gridlines,
      numeric axis labels (x = size/time, y = fluorescence/conc/molarity).
- [x] Viewport state: wheel zoom, drag-pan, drag-box zoom, reset/auto-fit.
- [x] Y-axis mode toggle: fluorescence ↔ concentration ↔ molarity.

## Phase C — File loading UX
- [x] Open dialog via `rfd` (`.xad` / `.xml` / `.xml.gz`).
- [x] Drag-drop of files (winit-backend `on_winit_window_event` →
      `WindowEvent::DroppedFile`; needs `unstable-winit-030` feature and is
      active when Slint selects winit; fallback = dialog).
- [x] `.xad` raw-channel view (Blue/Red detectors via `read_xad_raw_channels`),
      labeled "raw acquisition".
- [x] In-window error surfacing (no crashes on bad files).

## Phase D — Multi-select & overlay
- [x] Multi-select sample list (ctrl/⌘-click toggles, shift-click range; primary
      = last click, drives info line + viewport).
- [x] Overlay selected traces in one axes; per-sample Okabe–Ito color + legend
      (`render_overlay`; peak markers only in single-sample view).
- [x] Normalize toggle (peak-height-normalized) for shape comparison.

## Phase E — Overview tab
- [x] Small-multiples grid of all traces (mini-plots), labeled, ladder flagged
      (★ + green) — `overview.rs`, one plotters bitmap, up to 4 columns.
- [x] Shared/per-plot y-scale toggle.
- [x] "Show ladders" toggle, disabled by default, filters ladder wells from
      Overview traces/gel and click hit-testing.
- [x] Click a mini-plot → select it + jump to Detail (fractional click →
      `overview::cell_at`).

## Phase F — Peak / region table
- [x] `StandardTableView` under the Detail trace: peak
      #/size/time/area/height/%total/conc/molarity + note, followed by
      smear-region rows (`table.rs`; height = trace value at peak time, %total =
      area fraction).
- [x] Row↔plot cross-highlight (row→bold red marker line; plot click→nearest
      peak selects its row, within 3% of the visible x-span).

## Phase G — Gel / virtual-gel view
- [x] Fluorescence→grayscale lanes side-by-side, shared intensity scale, common
      migration-time axis (well at top), rotated lane labels (`gel.rs`).
- [x] "View: traces / gel" toggle in the Overview tab; clicking a lane opens it
      in Detail (`gel::lane_at`).

## Phase H — Manual low/high markers
- [x] `traceio`: `MarkerOverride` + `calculate_length_with(overrides)` +
      `marker_times()`. Overridden samples re-derive alignment (raw→aligned
      time) from user-placed raw marker times, keeping the ladder's canonical
      aligned targets; sizing/conc/molarity re-run via `loading::recalibrate_with`.
      Covered by `marker_override_shifts_sizing` test.
- [x] GUI: "Scale to markers" is enabled by default; disabling it switches the
      Detail plot to raw-time and draws green lower/upper lines. Drag a line
      (grab within 2% of the cursor) to move it with live transactional
      recalibration/table refresh; "Reset markers" restores automatic detection
      and recalibrates immediately.

## Risks
- Slint external file drag-drop isn't first-class → winit-backend handling;
  fallback is the open dialog.
- Keep the pixel↔data transform single-sourced (Phase A) so overlays and the
  plotters image stay in registration under zoom/pan.
