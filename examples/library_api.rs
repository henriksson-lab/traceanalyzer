//! Minimal library-consumer flow: detect, read, inspect metadata, and optionally save.
//!
//! Usage:
//!   cargo run --example library_api -- <input> [save-output]

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let input = args
        .next()
        .ok_or_else(|| anyhow::anyhow!("usage: library_api <input> [save-output]"))?;
    let output = args.next();

    let detected = traceanalyzer::traceio::io::detect_format(&input)?
        .ok_or_else(|| anyhow::anyhow!("unsupported electrophoresis path: {input}"))?;
    println!("detected: {:?}", detected.format);
    println!("identity: {}", detected.identity.display());
    println!("save: {:?}", detected.save_capabilities());

    let loaded = traceanalyzer::traceio::io::read_path(&input)?;
    println!("samples: {}", loaded.run.samples.len());
    println!("warnings: {}", loaded.warnings.len());

    if let Some(output) = output {
        if !traceanalyzer::traceio::io::supports_save_path(&loaded, &output) {
            anyhow::bail!(
                "detected source cannot be saved to this destination: {:?} -> {output}",
                loaded.source.format
            );
        }
        traceanalyzer::traceio::io::save_path(&loaded, &output)?;
        println!("saved: {output}");
    }

    Ok(())
}
