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
