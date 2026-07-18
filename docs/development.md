# Development

Building, testing, and packaging `traceanalyzer`. For how the code is organised
and how the analysis works, see [`architecture.md`](architecture.md); for the
native container format, [`xad_format.md`](xad_format.md).

## Building and running

```sh
# One-time fixture setup for parser tests and file-loading examples:
scripts/fetch-testdata.sh

# Headless summary of a run (no display needed):
cargo run -p traceio --example inspect -- testdata/demo_dna1000.xml.gz
cargo run -p traceio --example inspect -- testdata/demo_rna_nano.xml.gz

# GUI viewer (needs a display):
cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz

# macOS app bundle:
make osx-app
make osx-app-universal
```

The `inspect`/GUI loaders accept Bioanalyzer `.xad`/`.xml`/`.xml.gz`,
TapeStation exported XML/`_Electropherogram.csv` pairs, and Fragment Analyzer
`.raw` files, run directories, or zipped runs.

## Testing

```sh
# Parser/analysis tests (validate against real demo files — need fixtures):
cargo test -p traceio

# Headless render smoke test (no display or downloaded fixtures needed):
cargo test -p traceanalyzer --test headless_render
```

The `decodes_native_xad_if_present` test validates the container decode against
any real `.xad` under `bioa_examples/` (private, gitignored) or
`testdata/sample.xad`.

## Fixtures, display, and packaging prerequisites

- `scripts/fetch-testdata.sh` downloads the gitignored demo fixtures into
  `testdata/`. The `traceio` integration tests and file-loading examples
  expect `demo_dna1000.xml.gz` and `demo_rna_nano.xml.gz` there. The fixtures
  are real demo runs (DNA 1000 and Eukaryote Total RNA Nano) from the
  MIT-licensed [jwfoley/bioanalyzeR](https://github.com/jwfoley/bioanalyzeR).
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
