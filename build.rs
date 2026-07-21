fn main() {
    slint_build::compile("src/ui/app.slint").expect("compile Slint UI");

    // On Windows, embed the app icon into the .exe so Explorer/taskbar show it.
    // Gated to Windows hosts (where winresource + its toolchain are available);
    // a non-fatal warning keeps other setups building if the toolchain is absent.
    #[cfg(windows)]
    {
        println!("cargo:rerun-if-changed=assets/icon.ico");
        let mut res = winresource::WindowsResource::new();
        res.set_icon("assets/icon.ico");
        if let Err(e) = res.compile() {
            println!("cargo:warning=could not embed Windows icon: {e}");
        }
    }
}
