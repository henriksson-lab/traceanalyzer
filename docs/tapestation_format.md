# Agilent TapeStation export format

Notes for the TapeStation reader ([`src/traceio/tapestation.rs`](../src/traceio/tapestation.rs)),
ported from jwfoley/bioanalyzeR `R/tapestation.R` (MIT) and validated against
Agilent's own demo exports (the bioanalyzeR `inst/extdata/tapestation` fixtures)
plus Benchling's `allotropy` reference XMLs.

## Native file: not readable

TapeStation Analysis project files use assay-named extensions — `.D1000`,
`.HSD1000`, `.D5000`, `.HSD5000`, `.cfDNA`, `.gDNA`, `.RNA`, `.HSRNA` — and are
**password-encrypted ZIP archives with no public key** (Agilent declines to
disclose it; see bioanalyzeR issue #23, closed `wontfix`). There is no way to
read them directly, so — like bioanalyzeR — we read the **export** instead.

## Exported files (what we read)

*File → Export Data* in TapeStation Analysis Software (v4.1+) writes a pair that
share a stem:

- **`<name>.xml`** — all metadata (samples, peaks, regions, units); no trace.
- **`<name>_Electropherogram.csv`** — raw fluorescence only; no metadata.

Either opens the run; the reader derives the sibling by the
`_Electropherogram.csv` naming convention. Both may be gzip-compressed. Detection
is by that suffix, or (for `.xml`) by content — TapeStation XML has
`<FileInformation>`/`<Samples>` and no Bioanalyzer `<Chipset>` root.

### CSV (trace)

Latin-1, a header row, then one **column per lane** and one **row per distance
reading**. Header cells are lane labels (`A1: Ladder,B1: …`); we don't use them
(names come from the XML). The distance axis is the reversed, normalized row
index — `distance = rev(1..N)/N` — stored per sample as `time` (first row = 1.0).
Aligned exports are exactly 760 rows; we use the unaligned export.

### XML (metadata)

Flat schema (distinct from Bioanalyzer's — a separate parser, same target model):

- `FileInformation` → `FileName`, `RunEndDate` (creation date), `Assay` (kit name)
- `Assay/Units` → `MolecularWeightUnit` (bp/nt), `ConcentrationUnit`, `MolarityUnit`
- `Samples` → per `Sample`: `WellNumber` (`A1`…), `Comment` (= display name),
  `Observations` (ladder/marker detection), `RNA/RINe` or `DIN` (integrity → `rin`),
  `ScreenTapeID`, `Peaks`, `Regions`
- per `Peak`: `Size` → length, `RunDistance` → distance (percent ÷100),
  `FromPercent`/`ToPercent` → boundary distances (÷100), `Area`, `Height`,
  `CalibratedQuantity` → concentration, `Molarity`. Missing numeric fields are
  the literal `-` (parse to NaN).
- per `Region` (per-sample smear analysis): `From`/`To` (bp), plus vendor
  summary fields such as `AverageSize`, `Concentration`, `Molarity`, and
  `PercentOfTotal`. The shared model currently stores and surfaces only the
  region bounds.

### Sizing (distance → length)

Following bioanalyzeR: each sample's **Lower/Upper Marker** peaks (by
`Observations`) define a marker-relative distance
`rel = (distance − upper) / (lower − upper)` (upper = 0 when the assay has no
upper marker). The **ladder** sample (detected via `Observations` containing
"Ladder") fits a monotone Hyman spline from `rel` → `Size` (reusing
[`calibration::StandardCurve`]), and every sample's per-point `length` is that
spline evaluated at its own `rel`. Points outside the marker range get NaN, so
uncalibrated traces still plot on the distance axis.

## Model mapping & limitations

`traceanalyzer::traceio::io::read_path` dispatches via `tapestation::is_tapestation_path` and,
like the Fragment Analyzer path, fills `length` itself and **skips** the Bioanalyzer
`calibrate`. Peaks carry size/area/concentration/molarity and marker labels from
the XML. **Per-sample region bounds** are parsed into `Sample::regions` and
surfaced in the Table tab. **Multi-tape** runs are handled: each `ScreenTapeID` (kept in
`Sample::category`) is sized by its own ladder, with the sole ladder as a
fallback. TapeStation exports are **read-only** — a rename cannot be written back
(the native file is encrypted, the XML is a derived artifact), so the GUI
disables rename/save for these runs.

**Per-point concentration/molarity are intentionally not computed.** The export
gives concentration/molarity only per *peak* (surfaced in the Table tab), never
per trace point, and the Bioanalyzer per-point algorithm ([`concentration`])
does not transfer: it locates markers by matching fixed ladder concentrations,
whereas TapeStation marker concentrations vary per sample. Rather than fabricate
a misleading concentration trace, the y-axis for TapeStation offers fluorescence
only.

[`calibration::StandardCurve`]: ../src/traceio/calibration.rs
[`concentration`]: ../src/traceio/concentration.rs
