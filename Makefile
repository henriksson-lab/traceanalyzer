APP_NAME := Trace analyzer
APP_VERSION := $(shell sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n 1)
APP_BUNDLE := target/osx/$(APP_NAME).app
APP_BINARY := target/release/traceanalyzer
APP_UNIVERSAL_BINARY := target/osx/traceanalyzer-universal
APP_EXE := $(APP_BUNDLE)/Contents/MacOS/traceanalyzer
APP_PLIST := $(APP_BUNDLE)/Contents/Info.plist
# App icon: assets/icon.svg is the source of truth; assets/icon-1024.png is the
# committed 1024px master it renders to. The .icns is built from that master at
# packaging time with macOS built-ins only (sips + iconutil) — no SVG rasterizer.
APP_ICON_SRC := assets/icon-1024.png
APP_ICONSET := target/osx/AppIcon.iconset
APP_ICNS := $(APP_BUNDLE)/Contents/Resources/AppIcon.icns

# --- Linux install (freedesktop layout) -------------------------------------
# Standard prefix vars; packagers override DESTDIR (staging) and PREFIX.
PREFIX ?= /usr/local
DESTDIR ?=
BINDIR := $(DESTDIR)$(PREFIX)/bin
DATADIR := $(DESTDIR)$(PREFIX)/share
ICONDIR := $(DATADIR)/icons/hicolor
# The window's Wayland app_id / X11 WM_CLASS; the icon PNG/SVG and the .desktop
# StartupWMClass are all keyed to this name so the desktop finds the icon.
LINUX_APP_ID := traceanalyzer

.PHONY: osx-app osx-app-universal osx-bundle clean-osx-app install uninstall

# Install the release binary, the .desktop entry, and the icons into a
# freedesktop layout. Uses the scalable SVG (any size, no rasterizer needed)
# plus the committed 256px PNG as a fallback for themes that ignore SVG.
install:
	cargo build -p traceanalyzer --release
	install -Dm755 "$(APP_BINARY)" "$(BINDIR)/$(LINUX_APP_ID)"
	install -Dm644 packaging/traceanalyzer.desktop \
		"$(DATADIR)/applications/$(LINUX_APP_ID).desktop"
	install -Dm644 assets/icon.svg \
		"$(ICONDIR)/scalable/apps/$(LINUX_APP_ID).svg"
	install -Dm644 crates/traceanalyzer/assets/window-icon.png \
		"$(ICONDIR)/256x256/apps/$(LINUX_APP_ID).png"
	@printf 'Installed to %s. If not staging (DESTDIR empty), refresh caches:\n' "$(DESTDIR)$(PREFIX)"
	@printf '  update-desktop-database %s/applications\n' "$(DATADIR)"
	@printf '  gtk-update-icon-cache %s\n' "$(ICONDIR)"

uninstall:
	rm -f "$(BINDIR)/$(LINUX_APP_ID)" \
		"$(DATADIR)/applications/$(LINUX_APP_ID).desktop" \
		"$(ICONDIR)/scalable/apps/$(LINUX_APP_ID).svg" \
		"$(ICONDIR)/256x256/apps/$(LINUX_APP_ID).png"

# --- macOS .app bundle ------------------------------------------------------
osx-app:
	cargo build -p traceanalyzer --release
	$(MAKE) osx-bundle APP_BUNDLE_BINARY="$(APP_BINARY)"

osx-app-universal:
	rustup target add x86_64-apple-darwin aarch64-apple-darwin
	cargo build -p traceanalyzer --release --target x86_64-apple-darwin
	cargo build -p traceanalyzer --release --target aarch64-apple-darwin
	mkdir -p target/osx
	lipo -create \
		target/x86_64-apple-darwin/release/traceanalyzer \
		target/aarch64-apple-darwin/release/traceanalyzer \
		-output "$(APP_UNIVERSAL_BINARY)"
	lipo -verify_arch x86_64 arm64 "$(APP_UNIVERSAL_BINARY)"
	$(MAKE) osx-bundle APP_BUNDLE_BINARY="$(APP_UNIVERSAL_BINARY)"

osx-bundle:
	mkdir -p "$(APP_BUNDLE)/Contents/MacOS" "$(APP_BUNDLE)/Contents/Resources" "$(APP_ICONSET)"
	cp "$(APP_BUNDLE_BINARY)" "$(APP_EXE)"
	chmod +x "$(APP_EXE)"
	set -e; for pair in "16 icon_16x16" "32 icon_16x16@2x" "32 icon_32x32" \
		"64 icon_32x32@2x" "128 icon_128x128" "256 icon_128x128@2x" \
		"256 icon_256x256" "512 icon_256x256@2x" "512 icon_512x512" \
		"1024 icon_512x512@2x"; do \
		set -- $$pair; \
		sips -z $$1 $$1 "$(APP_ICON_SRC)" --out "$(APP_ICONSET)/$$2.png" >/dev/null; \
	done
	iconutil -c icns "$(APP_ICONSET)" -o "$(APP_ICNS)"
	printf '%s\n' \
		'<?xml version="1.0" encoding="UTF-8"?>' \
		'<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">' \
		'<plist version="1.0">' \
		'<dict>' \
		'  <key>CFBundleDevelopmentRegion</key>' \
		'  <string>en</string>' \
		'  <key>CFBundleDisplayName</key>' \
		'  <string>$(APP_NAME)</string>' \
		'  <key>CFBundleExecutable</key>' \
		'  <string>traceanalyzer</string>' \
		'  <key>CFBundleIconFile</key>' \
		'  <string>AppIcon</string>' \
		'  <key>CFBundleIdentifier</key>' \
		'  <string>org.traceanalyzer.TraceAnalyzer</string>' \
		'  <key>CFBundleInfoDictionaryVersion</key>' \
		'  <string>6.0</string>' \
		'  <key>CFBundleName</key>' \
		'  <string>$(APP_NAME)</string>' \
		'  <key>CFBundlePackageType</key>' \
		'  <string>APPL</string>' \
		'  <key>CFBundleShortVersionString</key>' \
		'  <string>$(APP_VERSION)</string>' \
		'  <key>CFBundleVersion</key>' \
		'  <string>$(APP_VERSION)</string>' \
		'  <key>LSMinimumSystemVersion</key>' \
		'  <string>11.0</string>' \
		'  <key>NSHighResolutionCapable</key>' \
		'  <true/>' \
		'</dict>' \
		'</plist>' \
		> "$(APP_PLIST)"
	@printf 'Built %s\n' "$(APP_BUNDLE)"

clean-osx-app:
	rm -rf "$(APP_BUNDLE)" "$(APP_ICONSET)" "$(APP_UNIVERSAL_BINARY)"
