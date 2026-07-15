APP_NAME := Trace analyzer
APP_BUNDLE := target/osx/$(APP_NAME).app
APP_BINARY := target/release/traceanalyzer
APP_EXE := $(APP_BUNDLE)/Contents/MacOS/traceanalyzer
APP_PLIST := $(APP_BUNDLE)/Contents/Info.plist

.PHONY: osx-app clean-osx-app

osx-app:
	cargo build -p traceanalyzer --release
	mkdir -p "$(APP_BUNDLE)/Contents/MacOS"
	cp "$(APP_BINARY)" "$(APP_EXE)"
	chmod +x "$(APP_EXE)"
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
		'  <key>CFBundleIdentifier</key>' \
		'  <string>org.traceanalyzer.TraceAnalyzer</string>' \
		'  <key>CFBundleInfoDictionaryVersion</key>' \
		'  <string>6.0</string>' \
		'  <key>CFBundleName</key>' \
		'  <string>$(APP_NAME)</string>' \
		'  <key>CFBundlePackageType</key>' \
		'  <string>APPL</string>' \
		'  <key>CFBundleShortVersionString</key>' \
		'  <string>0.1.0</string>' \
		'  <key>CFBundleVersion</key>' \
		'  <string>0.1.0</string>' \
		'  <key>LSMinimumSystemVersion</key>' \
		'  <string>11.0</string>' \
		'  <key>NSHighResolutionCapable</key>' \
		'  <true/>' \
		'</dict>' \
		'</plist>' \
		> "$(APP_PLIST)"
	@printf 'Built %s\n' "$(APP_BUNDLE)"

clean-osx-app:
	rm -rf "$(APP_BUNDLE)"
