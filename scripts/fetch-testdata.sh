#!/usr/bin/env bash
#
# Fetch the test fixtures used by `cargo test`.
#
# These are real demo runs from the MIT-licensed jwfoley/bioanalyzeR package.
# They are downloaded rather than committed, so binary blobs stay out of git.
# Idempotent: an existing non-empty, valid gzip file is left untouched.
#
set -euo pipefail

# Resolve repo root from this script's location so it works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
dest="$repo_root/testdata"

mkdir -p "$dest"

repo="https://raw.githubusercontent.com/jwfoley/bioanalyzeR/master/inst/extdata"

# each entry: "<local path under testdata/>|<remote path under $repo, URL-encoded>"
fixtures=(
  # Bioanalyzer exported XML.
  "demo_dna1000.xml.gz|bioanalyzer/Demo%20DNA%201000%20Series%20II.xml.gz"
  "demo_rna_nano.xml.gz|bioanalyzer/Demo%20Eukaryote%20Total%20RNA%20Nano%20Series%20II.xml.gz"
  # TapeStation exported metadata XML + electropherogram CSV (a matched pair;
  # the reader derives the CSV from the XML stem, so keep these names in sync).
  "tapestation/d1000.xml.gz|tapestation/D1000-Tubes-16-D1000.xml.gz"
  "tapestation/d1000_Electropherogram.csv.gz|tapestation/D1000-Tubes-16-D1000_Electropherogram.csv.gz"
)

is_gzip() {
  # true if the first two bytes are the gzip magic 1f 8b
  [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]
}

validate_gzip() {
  is_gzip "$1" && gzip -t "$1" >/dev/null 2>&1
}

for entry in "${fixtures[@]}"; do
  local_name="${entry%%|*}"
  remote_name="${entry##*|}"
  out="$dest/$local_name"
  mkdir -p "$(dirname "$out")"

  if [ -s "$out" ] && validate_gzip "$out"; then
    echo "skip   $local_name (already present)"
    continue
  fi

  if [ -e "$out" ]; then
    echo "stale  $local_name (missing or invalid gzip; replacing)" >&2
    rm -f "$out"
  fi

  echo "fetch  $local_name"
  tmp="$out.tmp.$$"
  if ! curl -fSL --retry 3 "$repo/$remote_name" -o "$tmp"; then
    rm -f "$tmp"
    echo "ERROR: download failed for $local_name" >&2
    exit 1
  fi
  if ! validate_gzip "$tmp"; then
    rm -f "$tmp"
    echo "ERROR: $local_name is not a valid gzip file (bad download?)" >&2
    exit 1
  fi
  mv "$tmp" "$out"
  echo "ok     $local_name ($(wc -c < "$out") bytes)"
done

echo "testdata ready in $dest"
