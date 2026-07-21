//! Tests for the `.xsy` assay reader. Data-gated on the uncommitted
//! `ext_software/` vendor install; skips cleanly when it is absent.

use std::path::PathBuf;

fn assays_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/ext_software/2100 expert/assays"
    ));
    dir.is_dir().then_some(dir)
}

#[test]
fn reads_compressed_dna1000() {
    let Some(dir) = assays_dir() else {
        eprintln!("skipping: ext_software/ not present");
        return;
    };
    let assay =
        traceanalyzer::traceio::xsy::read_xsy_file(dir.join("dsDNA/DNA 1000 Series II.xsy"))
            .expect("read");

    assert!(assay.compressed, "DNA 1000 is the Xceed-compressed variant");
    assert!(assay.name.contains("DNA 1000"), "name was {:?}", assay.name);
    assert_eq!(assay.size_unit, "bp");

    // DNA 1000 ladder: 13 peaks, 15 bp … 1500 bp, strictly increasing.
    assert_eq!(assay.ladder_peaks.len(), 13);
    assert_eq!(assay.ladder_peaks.first().unwrap().size, 15.0);
    assert_eq!(assay.ladder_peaks.last().unwrap().size, 1500.0);
    assert!(assay.ladder_peaks.windows(2).all(|w| w[1].size > w[0].size));

    // The run script decoded to a non-trivial list of numeric values.
    assert!(assay.script_values.len() > 100);
}

#[test]
fn reads_plain_rna_nano() {
    let Some(dir) = assays_dir() else {
        eprintln!("skipping: ext_software/ not present");
        return;
    };
    let assay = traceanalyzer::traceio::xsy::read_xsy_file(
        dir.join("RNA/Eukaryote Total RNA Nano Series II.xsy"),
    )
    .expect("read");

    assert!(
        !assay.compressed,
        "this RNA kit stores the method uncompressed"
    );
    assert!(assay.name.contains("RNA"), "name was {:?}", assay.name);
    assert_eq!(assay.size_unit, "nt");
    assert!(!assay.ladder_peaks.is_empty());
    assert!(assay.ladder_peaks.windows(2).all(|w| w[1].size > w[0].size));
}

/// Every shipped assay must parse, frame cleanly and yield a plausible ladder.
#[test]
fn all_shipped_assays_parse() {
    let Some(dir) = assays_dir() else {
        eprintln!("skipping: ext_software/ not present");
        return;
    };
    let mut count = 0;
    let mut stack = vec![dir];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d).unwrap().flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().and_then(|e| e.to_str()) == Some("xsy") {
                let assay = traceanalyzer::traceio::xsy::read_xsy_file(&path)
                    .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
                assert!(!assay.name.is_empty(), "{}: empty name", path.display());
                count += 1;
            }
        }
    }
    eprintln!("parsed {count} shipped .xsy assays");
    assert!(count > 0);
}
