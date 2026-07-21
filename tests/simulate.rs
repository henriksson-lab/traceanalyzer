//! Headless test of the `--simulate` GUI path (milestone M5): the simulator
//! turns a recorded `.pck` capture into a viewable run. Data-gated on the
//! uncommitted `ext_software/` install.

use std::path::PathBuf;

#[test]
fn simulate_pck_produces_a_viewable_run() {
    let dir = PathBuf::from(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/ext_software/2100 expert/data/packets"
    ));
    let Some(pck) = std::fs::read_dir(&dir).ok().and_then(|entries| {
        entries.flatten().map(|e| e.path()).find(|p| {
            p.extension().and_then(|e| e.to_str()) == Some("pck")
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("2100 expert_"))
        })
    }) else {
        eprintln!("skipping: no captures present");
        return;
    };

    let run = traceanalyzer::loading::simulate_pck(&pck).expect("simulate");
    assert_eq!(run.samples.len(), 1, "one continuous acquisition trace");
    let trace = &run.samples[0];
    assert!(!trace.fluorescence.is_empty(), "trace should have points");
    assert_eq!(
        trace.time.len(),
        trace.fluorescence.len(),
        "time and signal axes align"
    );
    assert!(!run.assay.assay_name.is_empty(), "assay name from header");
    eprintln!(
        "{}: {} points, assay {:?}",
        pck.display(),
        trace.fluorescence.len(),
        run.assay.assay_name
    );
}
