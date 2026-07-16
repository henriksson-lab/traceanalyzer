# Fragment Analyzer native run format (reverse-engineered)

Notes from reverse-engineering an Advanced Analytical / Agilent **Fragment
Analyzer** run directory (software version `1.2.0.11`). No public spec or parser
exists for the binary files; the vendor's ProSize **CSV export** (in
`Exported_data/`) was used as ground truth. Implemented in
[`crates/traceio/src/fa.rs`](../crates/traceio/src/fa.rs).

A run is a directory of sibling files sharing a timestamp stem, e.g.
`2025 11 19 16H 03M.{raw,raw2D,PKS,ANNT,GANNT,ANAI,RQN,current,txt}` plus
`method.mthd`, `CameraImage.bmp`, `ExpTime.txt`, `Timing.txt`, and an
`Exported_data/` folder of CSV/PNG exports.

The FA software is written in **LabVIEW**, so most binaries are big-endian and
several use LabVIEW *flattened data* conventions (4-byte big-endian length/count
prefixes; `.ANNT`/`.GANNT` are LabVIEW *Flatten To XML*).

## `.raw` — CCD acquisition (the electropherograms)

The centerpiece. All multi-byte values are **big-endian**.

| Offset | Meaning |
|---|---|
| `0x00` | magic `46 41 00 00` (`"FA\0\0"`) |
| `0x33` | version string (`"1.2.0.11"`) |
| `0x6c` | creation date string; `0x97` time string |
| `0xff` | **CCD line width** in pixels, `u16` (observed `701`) |
| ~`0x3e8` | **capillary centre-column table**: `[8×u16 zero][N×u16 columns, strictly increasing][8×u16 zero]`. `N` = capillary count (12). A second copy sits ~`0x5dc`. |
| `0x7d0` (`DATA_START`) | start of the CCD payload |

Payload = `scans × width` `u16` intensities, **row-major by scan**. The scan
count is derived: `scans = (file_size − 0x7d0) / (2 × width)` (observed `1501`,
matching `Timing.txt`'s frame count and the `1501` stored in `.PKS`).

**Extracting one capillary's electropherogram:** each capillary images to a
fixed pixel column (the centre-column table). Capillary *c*'s trace is that
column read across every scan — i.e. `value[scan] = data[scan*width + col_c]`,
averaged over a small `±CAP_WINDOW` pixel window to reduce noise. The scan index
is the migration-time axis (scan ≈ second). Capillary order = well order
(capillary 1 → well D1, per the `.txt` sidecar).

Validation: per-pixel-column variance across scans peaks exactly at the
centre-column positions; extracted, size-calibrated traces correlate 0.82–0.95
with the ProSize `Electropherogram.csv` for wells with real sample (empty wells
correlate low, as expected — only the co-injected markers are present).

## `.PKS` — peaks + size calibration (partially decoded)

LabVIEW-flattened, big-endian. Starts with `u32` count `= 12` (wells), then
per-well peak records mixing `u16` (scan times) and `f32` (RFU, area, conc, …).

**Size calibration (decoded):** a block `[u16 count=16][16×u16 anchor times]
[u16 0][u16 scans=1501]` (observed at `0xbe0`). The 16 anchor **times** are scan
numbers; they pair 1:1 with the standard ladder base-pair values
`[1,100,200,…,1000,1200,1500,2000,3000,6000]` (from the size-standard kit).
`fa.rs` locates this block as the longest strictly-increasing `u16` run within
`(0, scans]` (dropping the leading count element), then linearly interpolates
`scan → bp`. Points outside `[time₀, time₁₅]` get `NaN` length (uncalibrated).

**Peak table (decoded).** Layout, big-endian, from file start:

```text
u32 nwells
per well:
  20-byte summary  { u16 lm_apex, f32 lm_rfu, f32 _, u16 um_apex, f32 um_rfu, f32 _ }
  u32 npeaks;  npeaks × 26-byte record   (copy 1 — used)
  u32 npeaks;  npeaks × 26-byte record   (copy 2 — an aligned duplicate)
  8-byte trailer
```

Each **26-byte record** is `u16 start_scan, u16 apex_scan, u16 end_scan` then
5×`f32` `[_, RFU/height, _, _, corrected_area]` (the `f32`s sit on odd 2-byte
offsets; the three unlabelled ones, and where concentration/molarity live, are
still unknown). The per-well **summary** gives the lower/upper marker apex scan
times, which bound the reported peaks (records outside `[lm_apex, um_apex]`, e.g.
a sub-marker injection artifact, are dropped — matching the vendor table). The
lower/upper-marker peaks are assigned the ladder end sizes (1 bp / 6000 bp by
definition); other peaks are sized from the calibration via their apex scan.

Verified against the ProSize `Peak Table.csv`: e.g. D1 → LM `1 bp, area 23.0`,
sample `293 bp, area 77.4` (CSV 294 bp / 77.431), UM `6000 bp`; the ladder well
(D12) yields all 16 points `1,100,…,6000 bp`. `fa.rs` parses copy 1, labels the
markers, and flags the ladder well (`is_ladder` when ≥ 8 peaks). Any framing
inconsistency makes the parser return no peaks rather than fail the load.
Concentration/molarity are not populated (not yet located in `.PKS`).

## Other files (not used yet)

- **`.raw2D`** (`"## #"` magic) — a small 2-D CCD image / gel snapshot.
- **`.ANNT` / `.GANNT`** — LabVIEW Flatten-To-XML annotation trees (a 4-byte
  length prefix then `<Array>/<Cluster>/<NumElts>/<Dimsize>/<Val>` XML).
- **`.ANAI`** — INI-style per-capillary analysis settings.
- **`.current`** — TSV log of Current/Voltage/Pressure during the run.
- **`.txt`** — capillary → well → sample-name (used for names/wells).

## Model mapping

Fragment Analyzer runs arrive already size-calibrated, so the FA reader fills
each `Sample::length` directly (scan→bp interpolation) and **skips** the
Bioanalyzer marker-based `calibration` path. `loading::load` dispatches to it via
`fa::is_fa_path` (a `.raw` file with the `FA\0\0` magic, or a directory holding
one). Peaks (with lower/upper marker labels and ladder-well detection) are
populated from `.PKS`; per-point and per-peak concentration/molarity are not.
