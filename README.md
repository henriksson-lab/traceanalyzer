# traceanalyzer

Open-source, **post-measurement analysis** for automated-electrophoresis
instruments — a replacement for the vendor software that ships with the Agilent
**2100 Bioanalyzer**, Agilent **TapeStation**, and (bonus) Agilent/AATI
**Fragment Analyzer**. No hardware control; this reads and analyses saved runs.

## Status

| Platform | Native format | Reader status |
|---|---|---|
| **Bioanalyzer 2100** | `.xad` | native container **decode validated on real files** (Xceed-framed DEFLATE); raw detector signals + metadata + defined ladder extracted. Export-XML parser (processed per-well traces, peaks, sizing) **validated on real data**. |
| TapeStation | `.D1000`/… (encrypted ZIP) | not started — native is blocked; plan is XML + unaligned-CSV export |
| Fragment Analyzer | `.raw`/`.db3` (likely SQLite) | not started — plan is ProSize CSV first, then probe the SQLite schema |

This is an early prototype focused on the Bioanalyzer reader.

> **Key finding:** a native `.xad` stores **raw acquisition data + sample
> metadata only** — 2100 Expert recomputes processed per-well traces, peaks,
> sizing and RIN when the file is opened (only exports capture those). So
> reading a `.xad` today yields the raw detector electropherograms + metadata;
> reproducing the software's numbers needs the processing pipeline. Full
> details and the container spec: [`docs/xad_format.md`](docs/xad_format.md).

## Layout

- `crates/traceio` — GUI-free library: file readers + the shared
  `Electrophoresis` data model. Fast to build and test.
  - `model.rs` — common data model (assay info, samples, traces, peaks, ladder).
  - `bioanalyzer.rs` — parser for the Bioanalyzer inner XML.
  - `xad.rs` — native `.xad` container unwrap (base64 → DEFLATE → UTF-16LE).
  - `calibration.rs` — per-point ladder calibration (marker alignment +
    Hyman-filtered FMM spline) giving a size (bp/nt) for every trace point.
  - `concentration.rs` — per-point concentration + molarity (trapezoidal area,
    marker-ratio mass coefficients, molecular weight).
- `docs/xad_format.md` — reverse-engineered `.xad` container + schema spec.
- `scripts/fetch-testdata.sh` — download the demo fixtures (they aren't committed).
- `crates/traceanalyzer` — Slint GUI: sample list + electropherogram plot.
- `testdata/` — real demo runs (from jwfoley/bioanalyzeR, MIT).

## Try it

```sh
# One-time fixture setup for parser tests and file-loading examples:
scripts/fetch-testdata.sh

# Headless summary of a run (no display needed):
cargo run -p traceio --example inspect -- testdata/demo_dna1000.xml.gz
cargo run -p traceio --example inspect -- testdata/demo_rna_nano.xml.gz

# GUI viewer (needs a display):
cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz

# Parser tests (validate against real demo files):
cargo test -p traceio

# Headless render smoke test (no display or downloaded fixtures needed):
cargo test -p traceanalyzer --test headless_render
```

The `inspect`/GUI loaders accept `.xad` (native), `.xml`, and `.xml.gz`.

### Fixtures, display, and packaging prerequisites

- `scripts/fetch-testdata.sh` downloads the gitignored demo fixtures into
  `testdata/`. The `traceio` integration tests and file-loading examples
  expect `demo_dna1000.xml.gz` and `demo_rna_nano.xml.gz` there.
- `traceanalyzer` uses Slint plus the `winit` backend. Running the GUI needs a
  working display (`DISPLAY` or a Wayland session). In headless Linux CI, run GUI
  smoke commands under Xvfb, for example:
  `xvfb-run -a cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz`.
- Headless render examples/tests use the bitmap backend and vendored
  DejaVuSans font, so they do not require a display or system fonts.
- Linux builds may need the native libraries used by Slint/winit and file
  dialogs, commonly X11/Wayland development packages plus GTK development
  packages for `rfd` on distributions that do not install them by default.
- Packaging can be checked with `cargo package -p traceio`. The GUI crate's
  local `traceio` path dependency is versioned for registry packaging, but
  `cargo package -p traceanalyzer` requires `traceio 0.1.0` to be published in
  the target registry first; even `--no-verify` still resolves that registry
  dependency while preparing the package.

## The native `.xad` format

See [`docs/xad_format.md`](docs/xad_format.md) for the full reverse-engineered
specification. In brief: a `.xad` is a text/XML wrapper whose `<compressed_data>`
element holds base64 → an **Xceed-framed** blob → a **raw DEFLATE** stream →
**UTF-16LE XML**. `xad.rs` locates the element by name, parses the Xceed header's
self-describing lengths, inflates, and validates. A legacy 1-byte/9-byte framing
(original grimbough `readXAD.R`) is kept as a fallback.

The inner XML is the **raw acquisition** layout (`Chip/RawSignals/…` +
per-sample metadata), **not** the processed export — so `read_xad_file` returns
assay info + defined ladder + sample metadata, and `read_xad_raw_channels`
returns the raw detector electropherograms (Blue/Red). The
`bioanalyzer::parse_xml` element paths (`ProcessedSignal`, `PeakMolecular`,
`DARRIN`, …) apply to the **export** XML, which does contain processed per-well
results.

The `decodes_native_xad_if_present` test validates the container decode against
any real `.xad` under `bioa_examples/` (private, gitignored) or `testdata/sample.xad`.

## Per-point ladder calibration

`calibration::calculate_length` assigns a molecular length (bp/nt) to every
trace point, ported from bioanalyzeR `calculate.length`:

1. Each sample's marker peaks give a linear map from raw `MigrationTime` to
   marker-aligned time (so every sample is normalised to the ladder).
2. The ladder well's peaks fit a standard curve, aligned-time → length. The
   default is R's `splinefun(method = "hyman")` — an FMM cubic spline with the
   Hyman monotonicity filter — ported verbatim from R's `splines.c`/`spline.R`.
   Linear interpolation is also available.
3. Every point is mapped through the curve; points outside the ladder's range
   are left `NaN` (no extrapolation).

Validated against the demo runs: the DNA 1000 ladder calibrates to its 15–1500
bp markers, RNA 6000 Nano to ~25–5600 nt, monotone throughout. The GUI plots on
this bp/nt axis when a run is calibrated.

## Next steps

1. **`.xad` processing pipeline** — the big one: go from the raw `.xad` detector
   signal to processed per-well results (split the continuous acquisition into
   wells, baseline-subtract, detect peaks, align markers, size against the
   ladder), so `.xad` alone reproduces what the software shows. The sizing
   ([`calibration`]) and concentration/molarity ([`concentration`]) stages
   already exist and plug in once per-well traces + peaks are produced.
2. Gel-like (virtual-gel) rendering.
3. TapeStation (XML + CSV export) and Fragment Analyzer (ProSize CSV) readers
   into the same `Electrophoresis` model.

[`calibration`]: crates/traceio/src/calibration.rs
[`concentration`]: crates/traceio/src/concentration.rs

## Credits

Format knowledge and demo data derive from the MIT-licensed R packages
`jwfoley/bioanalyzeR` and `grimbough/bioanalyzeR`.
