# Development

Building, testing, and packaging `traceanalyzer`. For how the code is organised
and how the analysis works, see [`architecture.md`](architecture.md); for the
native container format, [`xad_format.md`](xad_format.md).

## Building and running

```sh
# One-time fixture setup for parser tests and file-loading examples:
bash scripts/fetch-testdata.sh

# Headless summary of a run (no display needed):
cargo run -p traceio --example inspect -- testdata/demo_dna1000.xml.gz
cargo run -p traceio --example inspect -- testdata/demo_rna_nano.xml.gz
cargo run -p traceio --example inspect -- testdata/tapestation/d1000.xml.gz

# GUI viewer (needs a display):
cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz

# macOS app bundle:
make osx-app
make osx-app-universal
```

The `inspect` example uses `traceio::io::read_path`, which accepts Bioanalyzer
`.xad`/`.xml`/`.xml.gz`, TapeStation exported XML/`_Electropherogram.csv` pairs,
and Fragment Analyzer `.fa.zip`/`.zip` archives, `.raw` files, or run
directories.

Library callers can use the path-oriented API without knowing the instrument
format ahead of time:

```rust
fn main() -> anyhow::Result<()> {
    let detected = traceio::io::detect_format("testdata/demo_dna1000.xml.gz")?
        .expect("supported electrophoresis file");
    println!("{:?}", detected.save_capabilities());

    let loaded = traceio::io::read_path("testdata/demo_dna1000.xml.gz")?;
    if traceio::io::supports_save_path(&loaded, "out.xml") {
        traceio::io::save_path(&loaded, "out.xml")?;
    }

    Ok(())
}
```

The lower-level `traceio::save::save_run(&run, &src, &dst)` API is still
available for applications that store the model and original source path
separately.

## Testing

```sh
# Parser/analysis tests (validate against real demo files — need fixtures):
cargo test -p traceio

# GUI crate tests, including the headless render smoke test:
cargo test -p traceanalyzer
```

The `decodes_native_xad_if_present` test validates the container decode against
any real `.xad` under `bioa_examples/` (private, gitignored) or
`testdata/sample.xad`.

## Fixtures, display, and packaging prerequisites

- `scripts/fetch-testdata.sh` downloads the gitignored demo fixtures into
  `testdata/`. The `traceio` integration tests and file-loading examples
  expect `demo_dna1000.xml.gz`, `demo_rna_nano.xml.gz`, and the matched
  TapeStation `tapestation/d1000.xml.gz` plus
  `tapestation/d1000_Electropherogram.csv.gz` pair there. The TapeStation pair
  is required for full `traceio` coverage of exported TapeStation XML/CSV
  loading. The fixtures are real demo runs (DNA 1000, Eukaryote Total RNA Nano,
  and TapeStation D1000) from the MIT-licensed
  [jwfoley/bioanalyzeR](https://github.com/jwfoley/bioanalyzeR). Run it as
  `bash scripts/fetch-testdata.sh` so the command works even when a source
  archive did not preserve executable bits.
- `traceanalyzer` uses Slint plus the `winit` backend. Running the GUI needs a
  working display (`DISPLAY` or a Wayland session). In headless Linux CI, run GUI
  smoke commands under Xvfb, for example:
  `xvfb-run -a cargo run -p traceanalyzer -- testdata/demo_dna1000.xml.gz`.
- Headless render tests use the bitmap backend and vendored DejaVuSans font, so
  that coverage does not require a display or system fonts.
- Linux builds may need the native libraries used by Slint/winit and file
  dialogs, commonly X11/Wayland development packages plus GTK development
  packages for `rfd` on distributions that do not install them by default.
- `make deb` derives the Debian architecture from `dpkg --print-architecture`
  and versioned ELF runtime dependencies with `dpkg-shlibdeps`; packagers can
  override `DEB_ARCH`, `DEB_DEPENDS` for extra non-ELF dependencies, and
  `DEB_RECOMMENDS`.
- Packaging can be checked with `cargo package -p traceio`. The GUI crate's
  local `traceio` path dependency is versioned for registry packaging, but
  `cargo package -p traceanalyzer` requires `traceio 0.1.0` to be published in
  the target registry first; even `--no-verify` still resolves that registry
  dependency while preparing the package.
