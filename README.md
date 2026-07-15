# traceanalyzer

Open-source, **post-measurement analysis** for automated-electrophoresis
instruments — a replacement for the vendor software that ships with the Agilent
**2100 Bioanalyzer**, Agilent **TapeStation**, and (bonus) Agilent/AATI
**Fragment Analyzer**. No hardware control; this reads and analyses saved runs.

## Status

| Platform | Native format | Reader status |
|---|---|---|
| **Bioanalyzer 2100** | `.xad` | inner-XML parser **validated on real data**; native `.xad` container unwrap ported, **awaiting a real `.xad` to validate** |
| TapeStation | `.D1000`/… (encrypted ZIP) | not started — native is blocked; plan is XML + unaligned-CSV export |
| Fragment Analyzer | `.raw`/`.db3` (likely SQLite) | not started — plan is ProSize CSV first, then probe the SQLite schema |

This is an early prototype focused on the Bioanalyzer reader.

## Layout

- `crates/traceio` — GUI-free library: file readers + the shared
  `Electrophoresis` data model. Fast to build and test.
  - `model.rs` — common data model (assay info, samples, traces, peaks, ladder).
  - `bioanalyzer.rs` — parser for the Bioanalyzer inner XML.
  - `xad.rs` — native `.xad` container unwrap (base64 → DEFLATE → UTF-16LE).
  - `calibration.rs` — per-point ladder calibration (marker alignment +
    Hyman-filtered FMM spline) giving a size (bp/nt) for every trace point.
- `crates/traceanalyzer` — Slint GUI: sample list + electropherogram plot.
- `testdata/` — real demo runs (from jwfoley/bioanalyzeR, MIT).

## Try it

```sh
# Headless summary of a run (no display needed):
cargo run -p traceio --example inspect -- testdata/demo_dna1000.xml.gz
cargo run -p traceio --example inspect -- testdata/demo_rna_nano.xml.gz

# GUI viewer (needs a display):
cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz

# Tests (validate the parser against real demo files):
cargo test -p traceio
```

The `inspect`/GUI loaders accept `.xad` (native), `.xml`, and `.xml.gz`.

## The native `.xad` format

A `.xad` is a line-oriented text/XML wrapper. The analytical payload is one
section of **base64** text that decodes to a **raw DEFLATE** stream framed by a
1-byte header and 9-byte trailer; inflating it yields a **UTF-16LE XML**
document. That inner XML is identical to the *File → Export to XML* output, so
one parser (`bioanalyzer::parse_xml`) serves both the native and export paths.

Inside the inner XML:
- Root `Chipset` → `Chips` → `Chip`.
- Raw trace: `…/Samples/Sample/DASignals/DetectorChannels/<ch>/SignalData`,
  with `ProcessedSignal` = base64 of little-endian **float32**, time axis
  `t_i = XStart + XStep·i`.
- Peaks: `…/DAResultStructures/DARIntegrator/Channel/PeaksMolecular/PeakMolecular`.
- RIN (RNA only): `…/DAResultStructures/DARRIN/Channel/RIN`.
- Ladder: `…/DAAssaySetpoints/DAMAssayInfoMolecular/LadderPeaks/LadderPeak`.

The container unwrap in `xad.rs` is a faithful port of grimbough/bioanalyzeR's
`readXAD.R` and depends on file-position magic constants (documented at the top
of that file) derived from specific samples. **To validate it, drop a real
`.xad` at `testdata/sample.xad`** and run `cargo test -p traceio` — the
`decodes_native_xad_if_present` test will exercise the full path.

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

1. Obtain real `.xad` files (several 2100 Expert versions) and validate/robustify
   the container unwrap.
2. **Concentration & molarity** per point (bioanalyzeR `calculate.concentration`
   / `calculate.molarity`: trapezoidal area, marker-ratio mass coefficients).
3. Gel-like (virtual-gel) rendering.
4. TapeStation (XML + CSV export) and Fragment Analyzer (ProSize CSV) readers
   into the same `Electrophoresis` model.

## Credits

Format knowledge and demo data derive from the MIT-licensed R packages
`jwfoley/bioanalyzeR` and `grimbough/bioanalyzeR`.
