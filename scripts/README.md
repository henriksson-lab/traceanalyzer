# scripts

## fetch-testdata.sh

Downloads the demo fixtures used by the `traceio` integration tests and
file-loading examples into `testdata/`, so the binary blobs don't have to be
committed to git. Headless renderer tests in `traceanalyzer` use synthetic data
and do not need these downloads.

```sh
bash scripts/fetch-testdata.sh
cargo test -p traceio
```

The fixtures are real demo runs (Bioanalyzer DNA 1000 and Eukaryote Total RNA
Nano, plus a TapeStation D1000 XML/electropherogram CSV export pair) from the
MIT-licensed [jwfoley/bioanalyzeR](https://github.com/jwfoley/bioanalyzeR)
package. The script is idempotent: it skips existing fixtures only after
validating that they are non-empty gzip files, and replaces stale or corrupt
cached files.

## build-macos-app.sh

Compatibility wrapper around `make osx-app`. It builds the canonical
`target/osx/Trace analyzer.app` bundle, then copies it under `dist/` at the repo
root. macOS only.

```sh
bash scripts/build-macos-app.sh          # build dist/Trace analyzer.app
bash scripts/build-macos-app.sh --dmg    # also build dist/trace-analyzer-<version>.dmg
```

The Makefile owns the bundle layout, icon generation, and `Info.plist` metadata,
so the script stays in sync with CI and `make osx-app`.

The bundle is **unsigned and un-notarized**, so Gatekeeper will warn on first
launch; see the commented `codesign --deep --sign` note in the script for
future signing. The `dist/` output directory is git-ignored.
