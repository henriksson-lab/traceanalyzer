# Architecture

How `traceanalyzer` is organised and how the analysis works. For build/test
instructions see [`development.md`](development.md); for the native container
format, [`xad_format.md`](xad_format.md).

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
- `crates/traceanalyzer` — Slint GUI: well tree + electropherogram plot.
- `docs/xad_format.md` — reverse-engineered `.xad` container + schema spec.
- `scripts/fetch-testdata.sh` — download the demo fixtures (they aren't committed).
- `testdata/` — real demo runs (from jwfoley/bioanalyzeR, MIT).

## The native `.xad` format

See [`xad_format.md`](xad_format.md) for the full reverse-engineered
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

## Roadmap

1. **`.xad` processing pipeline** — the big one: go from the raw `.xad` detector
   signal to processed per-well results (split the continuous acquisition into
   wells, baseline-subtract, detect peaks, align markers, size against the
   ladder), so `.xad` alone reproduces what the software shows. The sizing
   ([`calibration`]) and concentration/molarity ([`concentration`]) stages
   already exist and plug in once per-well traces + peaks are produced.
2. Gel-like (virtual-gel) rendering.
3. TapeStation (XML + CSV export) and Fragment Analyzer (ProSize CSV) readers
   into the same `Electrophoresis` model.

[`calibration`]: ../crates/traceio/src/calibration.rs
[`concentration`]: ../crates/traceio/src/concentration.rs
