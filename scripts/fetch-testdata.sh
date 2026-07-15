#!/usr/bin/env bash
#
# Fetch the test fixtures used by `cargo test`.
#
# These are real demo runs from the MIT-licensed jwfoley/bioanalyzeR package.
# They are downloaded rather than committed, so binary blobs stay out of git.
# Idempotent: an existing non-empty file is left untouched.
#
set -euo pipefail

# Resolve repo root from this script's location so it works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
dest="$repo_root/testdata"

mkdir -p "$dest"

base="https://raw.githubusercontent.com/jwfoley/bioanalyzeR/master/inst/extdata/bioanalyzer"

# local name  <TAB>  remote filename (URL-encoded)
fixtures=(
  "demo_dna1000.xml.gz|Demo%20DNA%201000%20Series%20II.xml.gz"
  "demo_rna_nano.xml.gz|Demo%20Eukaryote%20Total%20RNA%20Nano%20Series%20II.xml.gz"
)

is_gzip() {
  # true if the first two bytes are the gzip magic 1f 8b
  [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]
}

for entry in "${fixtures[@]}"; do
  local_name="${entry%%|*}"
  remote_name="${entry##*|}"
  out="$dest/$local_name"

  if [ -s "$out" ]; then
    echo "skip   $local_name (already present)"
    continue
  fi

  echo "fetch  $local_name"
  tmp="$out.tmp.$$"
  if ! curl -fSL --retry 3 "$base/$remote_name" -o "$tmp"; then
    rm -f "$tmp"
    echo "ERROR: download failed for $local_name" >&2
    exit 1
  fi
  if ! is_gzip "$tmp"; then
    rm -f "$tmp"
    echo "ERROR: $local_name is not a gzip file (bad download?)" >&2
    exit 1
  fi
  mv "$tmp" "$out"
  echo "ok     $local_name ($(wc -c < "$out") bytes)"
done

echo "testdata ready in $dest"
