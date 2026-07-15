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

The fixtures are real demo runs (DNA 1000 and Eukaryote Total RNA Nano) from the
MIT-licensed [jwfoley/bioanalyzeR](https://github.com/jwfoley/bioanalyzeR)
package. The script is idempotent — it skips any fixture already present — and
verifies each download is a valid gzip file.

## build-macos-app.sh

Packages the release `traceanalyzer` binary into a macOS `.app` bundle (and,
optionally, a `.dmg`) under `dist/` at the repo root. macOS only.

```sh
bash scripts/build-macos-app.sh          # build dist/traceanalyzer.app
bash scripts/build-macos-app.sh --dmg    # also build dist/traceanalyzer-<version>.dmg
```

The script runs `cargo build --release -p traceanalyzer`, assembles
`traceanalyzer.app/Contents/{MacOS,Resources}`, and writes an `Info.plist`
declaring `.xad`/`.xml`/`.xml.gz` as openable document types. It is idempotent
(the bundle is rebuilt from scratch each run) and validates the generated
plist with `plutil -lint`.

The bundle is **unsigned and un-notarized**, so Gatekeeper will warn on first
launch; see the commented `codesign --deep --sign` note in the script for
future signing. There is no app icon yet — add `Resources/traceanalyzer.icns`
(the script picks up `scripts/traceanalyzer.icns` automatically) to include one.
The `dist/` output directory is git-ignored.
