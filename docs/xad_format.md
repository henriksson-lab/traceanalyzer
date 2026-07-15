# Agilent 2100 Bioanalyzer `.xad` file format

Reverse-engineered specification of the native `.xad` "chip data" file written
by Agilent **2100 Expert** software, as implemented by `crates/traceio/src/xad.rs`.

This documents the format as observed in **real files from current 2100 Expert
output** (High Sensitivity DNA assay, 2026). Where it differs from the original
reverse engineering in grimbough/bioanalyzeR `readXAD.R`, that is called out
explicitly under [Differences vs. the original parser](#differences-vs-the-original-readxadr-parser).

## 1. Overview

A `.xad` is a text/XML wrapper file. Its analytical payload is a single
`<compressed_data>` element holding a base64-encoded, DEFLATE-compressed
**UTF-16LE XML** document. Decompressing it yields the "inner XML".

Crucially — see [§5](#5-key-finding-xad-holds-raw-acquisition-not-results) —
the inner XML of these files is the **raw acquisition** layout: raw detector
signals plus per-sample *metadata*. It does **not** contain the processed
per-well electropherograms, called peaks, sizing, concentration, or RIN. 2100
Expert recomputes all of those when the file is opened; only exports/reports
(`File → Export to XML`, PDF) capture the processed results.

## 2. Outer container

The file is line-oriented text. Reading it as UTF-8 (lossy) is sufficient to
locate the payload. Structure:

```
<?xml version="1.0" encoding="UTF-8"?>
<!--Do not edit this comment tag:<uuid>:10:<preview blob>...
    ...<base64 of a small UTF-16LE "Preview"/GelImage document>-->
<!--Do not edit this header information:<uuid>:...-->
<compressed_data xmlns:dt="urn:schemas-microsoft-com:datatypes" dt:dt="bin.base64">
<base64 ... spanning thousands of lines ...>
</compressed_data>
```

- The first `<!--…comment tag…-->` holds a base64 **preview/gel thumbnail**
  (base64 → UTF-16LE → a tiny `<Preview>/<GelImage>` XML). Not required for
  analysis; ignored by the reader.
- The **main payload** is the `<compressed_data … dt:dt="bin.base64">` element
  (a Microsoft XML "data island"). Everything between the opening tag's `>` and
  `</compressed_data>`, with whitespace/newlines removed, is the base64 blob.

The reader extracts the payload by locating the `<compressed_data>` element by
name (robust to line layout), **not** by line index.

## 3. The compressed blob (Xceed container)

base64-decoding the `<compressed_data>` body yields a blob produced by the
**Xceed** compression library (a .NET component). Layout:

| offset | size | contents |
|-------:|-----:|----------|
| 0  | 4  | `u32` = **2** — format marker/version (constant) |
| 4  | 4  | `u32` = **79** (0x4F) — constant; purpose unconfirmed (likely an Xceed method id) |
| 8  | 4  | `u32` **uncompressed length** — byte length of the inflated UTF-16LE XML |
| 12 | 4  | `u32` **compressed length** — byte length of the DEFLATE stream |
| 16 | 4  | `u32` = **10** — byte length of the marker string that follows |
| 20 | 10 | `"Xceed"` in UTF-16LE |
| 30 | 12 | `"SCO,10"` in UTF-16LE — Xceed method/level string |
| 42 | 34 | zero padding |
| 76 | *compressed length* | **raw DEFLATE stream** (no zlib/gzip header) |
| end−9 | 9 | trailer (observed all-zero / small; not needed) |

Both sample files place the DEFLATE stream at offset **76**, and satisfy:

```
deflate_start = blob.len() − 9 − compressed_length      # = 76
inflated.len() == uncompressed_length                    # self-validating
```

The reader locates the stream via those **header-derived lengths** (not a
hard-coded 76), raw-inflates it (`flate2` `DeflateDecoder`, i.e. `wbits = -15`),
and asserts the inflated length equals the header's uncompressed length.

The inflated bytes are **UTF-16LE** XML with no BOM, beginning directly with
`<Chipset …>` (no `<?xml?>` prolog).

## 4. Inner XML schema (raw acquisition)

Root `<Chipset>` → `<Chips>` → `<Chip>`. Two relevant branches:

### 4a. `Chip/RawSignals/` — raw instrument acquisition

Each channel is a `…/SignalData` with scalar fields
(`ChannelID`, `Name`, `XStart`, `XStep`, `NumberOfSamples`, `UnitX`, `UnitY`,
`AlignmentBias`, `AlignmentScale`, …) and a `RawSignal` element:

- `RawSignal` = **base64 of little-endian float32** — the raw signal samples.
- ⚠️ `NumberOfSamples` is **0/unreliable** in the raw file; decode the actual
  count from the base64 length (`bytes/4`).
- Time axis: `t_i = XStart + XStep·i` (detector `XStep` ≈ 0.05 s).

Channels observed:

| path under `RawSignals` | meaning |
|---|---|
| `BoardTemperature/SignalData`, `ChipTemperature/SignalData` | housekeeping |
| `Voltages/Voltage/SignalData`, `Currents/Current/SignalData` | housekeeping |
| `HorizontalPosition`, `VerticalPosition`, `DebugSignal` | housekeeping |
| **`DetectorChannels/Channel/SignalData`** | **the fluorescence detectors** |

The detector channels are the electropherogram signal, e.g. `BlueFluorescence`
and `RedFluorescence`, each a single **whole-chip continuous acquisition**
(~38k samples covering *all* wells in time order, before the software splits
them per well).

### 4b. `Chip/Files/File/Samples/Sample` — per-sample metadata

Same path as the export. But in the `.xad` each `<Sample>` contains **only
metadata**: `Index`, `HasData`, `Category` (e.g. `"Ladder"`), `Name`,
`Comment`, `WellNumber`, `DASampleSetpoints`, review/workflow status, etc.
There is **no `DASignals`**, no per-well `ProcessedSignal`, and no
`DAResultStructures`/`PeakMolecular` under the samples.

Assay setpoints **are** present: `Chip/AssayBody/DAAssaySetpoints/
DAMAssayInfoMolecular` (`SizeUnit`, `ConcentrationUnit`, `LadderPeaks` — the
*defined* ladder), and `DAMAlignment`/`AlignUpperMarker`.

## 5. KEY FINDING: `.xad` holds raw acquisition, not results

For these files, the native `.xad` and the `Export to XML` output are **not the
same document**:

| element | native `.xad` inner XML | `Export to XML` |
|---|---:|---:|
| `ProcessedSignal` (per-well) | 0 | 27 |
| `PeakMolecular` (called peaks) | 0 | 225 |
| per-sample `DASignals`/`DetectorChannels` | 0 | 13 |
| `RawSignal` (raw detector) | 39 | 66 |
| `Sample` (metadata) | 12 | 12 |

The inner XML diverges from the export early (in a `<Method>`/`ScriptLines`
section) and is ~30% smaller. **2100 Expert derives ProcessedSignal, calls
peaks, aligns markers, sizes against the ladder and computes RIN on open** — the
`.xad` persists only the raw signals and setpoints.

**Implication for an open-source replacement:** decoding the container gets you
the raw detector electropherograms + sample metadata + assay setpoints. To
reproduce the numbers the software shows (and that appear in exports), we must
implement the processing pipeline ourselves: per-well splitting from the
continuous acquisition, baseline subtraction, marker alignment, peak detection,
ladder sizing, region/smear analysis, and the quality metrics. `traceio`
already has the ladder **sizing** ([`calibration`]) and **concentration/
molarity** ([`concentration`]) steps, which operate once per-well traces + peaks
exist — i.e. they currently consume the **export** XML; wiring them to raw
`.xad` data needs the front half of that pipeline.

## 6. Differences vs. the original `readXAD.R` parser

grimbough/bioanalyzeR `readXAD.R` was the starting point. What we found differs:

1. **Payload location.** `readXAD.R` slices the payload by *line index* — it
   takes the block between the 5th and 6th lines containing `<`, trims a fixed
   72-char prefix / 18-char suffix, and anchors the base64 on a literal `"Oy9"`
   marker (one user had to patch this to `"Ox9"`). That is brittle and **fails
   on these files** (they have only 4 such lines). We instead locate the
   `<compressed_data>` element **by name**, which is robust across layouts.
2. **Container framing.** `readXAD.R` strips **1 leading + 9 trailing** bytes
   from the decoded blob and raw-inflates. These files instead carry a **76-byte
   Xceed header** (with self-describing uncompressed/compressed lengths) before
   the DEFLATE stream, plus the same 9-byte trailer. The reader detects the
   `"Xceed"` marker and uses the header; it keeps the 1-byte/9-byte legacy path
   as a fallback for older files.
3. **Contents.** The inner XML `readXAD.R` + jwfoley `bioanalyzer.R` consume has
   per-well `ProcessedSignal` and `PeakMolecular` (i.e. processed results).
   These real `.xad` files contain **neither** — only raw signals + metadata
   (§5). So on current files, native `.xad` is a raw-data source, and the
   processed per-well data must be computed or read from an export.

## 7. Decoder recipe (summary)

```text
1. Read file bytes.
2. Find "<compressed_data" … ">"  … "</compressed_data>"; take the base64 body.
3. Strip ASCII whitespace; base64-decode  ->  blob.
4. If blob[20..30] == "Xceed" (UTF-16LE):           # Xceed container
     uncompressed = u32le(blob, 8);  compressed = u32le(blob, 12)
     end   = blob.len() - 9
     start = end - compressed
     inflated = raw_inflate(blob[start..end])        # DEFLATE, wbits = -15
     assert inflated.len() == uncompressed
   else:                                             # legacy readXAD.R
     inflated = raw_inflate(blob[1 .. blob.len()-9])
5. inner_xml = decode_utf16le(inflated)              # starts with "<Chipset"
```

## 8. Provenance

Reverse-engineered from real 2100 Expert `.xad` files (High Sensitivity DNA
assay) together with their `Export to XML` counterparts, cross-checked with
grimbough/bioanalyzeR (`readXAD.R`) and jwfoley/bioanalyzeR (`bioanalyzer.R`),
both MIT. Sample files are private lab data and are not committed to the repo.
