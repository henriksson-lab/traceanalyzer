# scripts

## fetch-testdata.sh

Downloads the test fixtures that `cargo test` needs into `testdata/`, so the
binary blobs don't have to be committed to git.

```sh
bash scripts/fetch-testdata.sh
cargo test
```

The fixtures are real demo runs (DNA 1000 and Eukaryote Total RNA Nano) from the
MIT-licensed [jwfoley/bioanalyzeR](https://github.com/jwfoley/bioanalyzeR)
package. The script is idempotent — it skips any fixture already present — and
verifies each download is a valid gzip file.
