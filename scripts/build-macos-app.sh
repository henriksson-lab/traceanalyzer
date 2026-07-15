#!/usr/bin/env bash
#
# Package the release `traceanalyzer` binary into a macOS `.app` bundle,
# and optionally a distributable `.dmg`.
#
# Output goes to `dist/` at the repo root (git-ignored). The script is
# idempotent: it rebuilds the bundle from scratch on each run.
#
# Usage:
#   bash scripts/build-macos-app.sh            # build the .app
#   bash scripts/build-macos-app.sh --dmg      # also build a .dmg (needs hdiutil)
#
set -euo pipefail

# Resolve repo root from this script's location so it works from any CWD.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

app_name="traceanalyzer"
bundle_id="org.traceanalyzer.app"
version="0.1.0"   # keep in sync with [workspace.package] version in Cargo.toml
min_macos="11.0"

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
app_dir="$dist_dir/$app_name.app"
macos_dir="$app_dir/Contents/MacOS"
resources_dir="$app_dir/Contents/Resources"
binary="$repo_root/target/release/$app_name"

echo "==> Building release binary ($app_name)"
( cd "$repo_root" && cargo build --release -p "$app_name" )

if [ ! -x "$binary" ]; then
  echo "ERROR: expected binary not found at $binary" >&2
  exit 1
fi

echo "==> Assembling $app_name.app"
rm -rf "$app_dir"
mkdir -p "$macos_dir" "$resources_dir"
cp "$binary" "$macos_dir/$app_name"
chmod +x "$macos_dir/$app_name"

# --- App icon -------------------------------------------------------------
# There is currently no app icon. If/when one is added, drop an .icns file at
# Resources/$app_name.icns and uncomment the CFBundleIconFile entry in the
# Info.plist below. Do NOT fabricate binary icon assets here.
# TODO: add Resources/traceanalyzer.icns and reference it via CFBundleIconFile.
icon_plist_entry=""
if [ -f "$repo_root/scripts/$app_name.icns" ]; then
  cp "$repo_root/scripts/$app_name.icns" "$resources_dir/$app_name.icns"
  icon_plist_entry="	<key>CFBundleIconFile</key>
	<string>$app_name.icns</string>
"
  echo "    bundled icon: $app_name.icns"
else
  echo "    no icon found (skipping; see TODO in script)"
fi

echo "==> Writing Info.plist"
cat > "$app_dir/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
	<key>CFBundleName</key>
	<string>$app_name</string>
	<key>CFBundleDisplayName</key>
	<string>$app_name</string>
	<key>CFBundleIdentifier</key>
	<string>$bundle_id</string>
	<key>CFBundleExecutable</key>
	<string>$app_name</string>
	<key>CFBundleVersion</key>
	<string>$version</string>
	<key>CFBundleShortVersionString</key>
	<string>$version</string>
	<key>CFBundlePackageType</key>
	<string>APPL</string>
	<key>CFBundleInfoDictionaryVersion</key>
	<string>6.0</string>
	<key>LSMinimumSystemVersion</key>
	<string>$min_macos</string>
	<key>NSHighResolutionCapable</key>
	<true/>
${icon_plist_entry}	<key>CFBundleDocumentTypes</key>
	<array>
		<dict>
			<key>CFBundleTypeName</key>
			<string>Electrophoresis trace</string>
			<key>CFBundleTypeRole</key>
			<string>Viewer</string>
			<key>LSHandlerRank</key>
			<string>Owner</string>
			<key>CFBundleTypeExtensions</key>
			<array>
				<string>xad</string>
				<string>xml</string>
				<string>xml.gz</string>
			</array>
		</dict>
	</array>
</dict>
</plist>
PLIST

# Blank the extended-attribute quarantine flag if present (best-effort).
xattr -cr "$app_dir" 2>/dev/null || true

# --- Code signing (not performed) -----------------------------------------
# This bundle is unsigned and un-notarized; Gatekeeper will warn on first run.
# To sign for distribution in the future:
#   codesign --deep --force --options runtime \
#     --sign "Developer ID Application: <NAME> (<TEAMID>)" "$app_dir"
# followed by notarization via `xcrun notarytool submit`.

if command -v plutil >/dev/null 2>&1; then
  plutil -lint "$app_dir/Contents/Info.plist" >/dev/null
  echo "    Info.plist validated (plutil -lint)"
fi

echo "ok     $app_dir"

if [ "$make_dmg" = true ]; then
  if command -v hdiutil >/dev/null 2>&1; then
    dmg_path="$dist_dir/$app_name-$version.dmg"
    echo "==> Building $app_name-$version.dmg"
    rm -f "$dmg_path"
    hdiutil create \
      -volname "$app_name" \
      -srcfolder "$app_dir" \
      -ov -format UDZO \
      "$dmg_path" >/dev/null
    echo "ok     $dmg_path"
  else
    echo "WARN: --dmg requested but hdiutil not found; skipping .dmg" >&2
  fi
fi

echo "done   output in $dist_dir"
