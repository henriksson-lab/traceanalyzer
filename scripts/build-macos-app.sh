#!/usr/bin/env bash
#
# Compatibility wrapper for the Makefile macOS bundle target.
#
# Output is copied to `dist/` at the repo root (git-ignored). The canonical app
# bundle is `target/osx/Trace analyzer.app`.
#
# Usage:
#   bash scripts/build-macos-app.sh            # build the .app
#   bash scripts/build-macos-app.sh --dmg      # also build a .dmg (needs hdiutil)
#
set -euo pipefail

# Resolve repo root from this script's location so it works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

app_name="Trace analyzer"
version="$(sed -n 's/^version = "\(.*\)"/\1/p' "$repo_root/Cargo.toml" | head -n 1)"

make_dmg=false
for arg in "$@"; do
  case "$arg" in
    --dmg) make_dmg=true ;;
    -h|--help)
      echo "Usage: bash scripts/build-macos-app.sh [--dmg]"
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $arg" >&2
      exit 1
      ;;
  esac
done

if [ "$(uname -s)" != "Darwin" ]; then
  echo "ERROR: this script builds a macOS .app and must run on macOS (Darwin)." >&2
  exit 1
fi

dist_dir="$repo_root/dist"
app_dir="$repo_root/target/osx/$app_name.app"
dist_app_dir="$dist_dir/$app_name.app"

echo "==> Building $app_name.app"
( cd "$repo_root" && make osx-app )

echo "==> Copying app bundle to dist/"
rm -rf "$dist_app_dir"
mkdir -p "$dist_dir"
ditto "$app_dir" "$dist_app_dir"

# Blank the extended-attribute quarantine flag if present (best-effort).
xattr -cr "$dist_app_dir" 2>/dev/null || true

# --- Code signing (not performed) -----------------------------------------
# This bundle is unsigned and un-notarized; Gatekeeper will warn on first run.
# To sign for distribution in the future:
#   codesign --deep --force --options runtime \
#     --sign "Developer ID Application: <NAME> (<TEAMID>)" "$app_dir"
# followed by notarization via `xcrun notarytool submit`.

if command -v plutil >/dev/null 2>&1; then
  plutil -lint "$dist_app_dir/Contents/Info.plist" >/dev/null
  echo "    Info.plist validated (plutil -lint)"
fi

echo "ok     $dist_app_dir"

if [ "$make_dmg" = true ]; then
  if command -v hdiutil >/dev/null 2>&1; then
    dmg_path="$dist_dir/trace-analyzer-$version.dmg"
    echo "==> Building $app_name-$version.dmg"
    rm -f "$dmg_path"
    hdiutil create \
      -volname "$app_name" \
      -srcfolder "$dist_app_dir" \
      -ov -format UDZO \
      "$dmg_path" >/dev/null
    echo "ok     $dmg_path"
  else
    echo "WARN: --dmg requested but hdiutil not found; skipping .dmg" >&2
  fi
fi

echo "done   output in $dist_dir"
