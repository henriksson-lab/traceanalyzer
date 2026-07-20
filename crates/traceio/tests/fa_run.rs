//! Validation of the reverse-engineered Fragment Analyzer reader against a real
//! run. The run directory holds private instrument data (`fa_examples/`, which
//! is git-ignored), so every test skips cleanly when it is not present.

use std::io::Write;
use std::path::{Path, PathBuf};

/// The reference FA run directory, or `None` if the private data is absent.
fn fa_run_dir() -> Option<PathBuf> {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fa_examples/16-03-27");
    dir.is_dir().then_some(dir)
}

/// Zip the run directory's top-level files (flat entries) into `dest`.
fn zip_run_dir(dir: &Path, dest: &Path) {
    let out = std::fs::File::create(dest).expect("create zip");
    let mut w = zip::ZipWriter::new(out);
    let opts = zip::write::SimpleFileOptions::default()
        .compression_method(zip::CompressionMethod::Deflated);
    for entry in std::fs::read_dir(dir).expect("read run dir").flatten() {
        let p = entry.path();
        if p.is_file() {
            let name = p.file_name().unwrap().to_string_lossy().into_owned();
            w.start_file(name, opts).expect("start entry");
            w.write_all(&std::fs::read(&p).expect("read file"))
                .expect("write entry");
        }
    }
    w.finish().expect("finish zip");
}

#[test]
fn reads_fa_run_traces_and_names() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    assert_eq!(run.assay.length_unit, "bp");
    assert_eq!(run.samples.len(), 12, "expected 12 capillaries");
    for (i, s) in run.samples.iter().enumerate() {
        assert_eq!(s.well_number, (i + 1) as i32);
        assert!(!s.name.is_empty(), "sample {i} has no name");
        assert!(
            s.name.starts_with(&format!("D{}", i + 1)),
            "well label in name: {}",
            s.name
        );
        assert_eq!(s.time.len(), s.fluorescence.len());
        assert!(s.fluorescence.len() > 1000, "trace too short");
    }
}

#[test]
fn fa_size_calibration_places_main_peak() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    // Calibration must span the 1..6000 bp ladder.
    let s0 = &run.samples[0];
    let finite: Vec<f64> = s0
        .length
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    assert!(finite.len() > 500, "too few calibrated points");
    let lo = finite.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = finite.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert!(
        (lo - 1.0).abs() < 5.0,
        "ladder should start ~1 bp, got {lo}"
    );
    assert!(
        (hi - 6000.0).abs() < 50.0,
        "ladder should end ~6000 bp, got {hi}"
    );

    // D1's main sample peak is at ~294 bp (per the vendor Peak Table). Find the
    // brightest calibrated point in the 150..1000 bp window and check it lands
    // there — this validates .raw extraction and .PKS calibration together.
    let mut best_bp = f64::NAN;
    let mut best_f = f32::NEG_INFINITY;
    for i in 0..s0.length.len() {
        let bp = s0.length[i];
        if bp.is_finite() && (150.0..1000.0).contains(&bp) && s0.fluorescence[i] > best_f {
            best_f = s0.fluorescence[i];
            best_bp = bp;
        }
    }
    assert!(
        (best_bp - 294.0).abs() < 45.0,
        "D1 main peak should be ~294 bp, got {best_bp}"
    );
}

#[test]
fn fa_peaks_and_markers_from_pks() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    // D1: lower marker (1 bp), a ~294 bp sample peak, upper marker (6000 bp).
    let d1 = &run.samples[0];
    assert_eq!(d1.peaks.len(), 3, "D1 should have LM + sample + UM");
    let lm = &d1.peaks[0];
    let um = d1.peaks.last().unwrap();
    assert_eq!(lm.observations, "Lower Marker");
    assert_eq!(um.observations, "Upper Marker");
    assert!((lm.length - 1.0).abs() < 0.5, "LM sized 1 bp");
    assert!((um.length - 6000.0).abs() < 0.5, "UM sized 6000 bp");
    let sample = &d1.peaks[1];
    assert!(
        (sample.length - 294.0).abs() < 20.0,
        "D1 sample ~294 bp, got {}",
        sample.length
    );
    assert!(
        (sample.area - 77.4).abs() < 2.0,
        "D1 sample area ~77.4, got {}",
        sample.area
    );

    // Exactly one well is the size ladder (D12), with the full 16-point ladder.
    let ladders: Vec<&traceio::Sample> = run.samples.iter().filter(|s| s.is_ladder).collect();
    assert_eq!(ladders.len(), 1, "one ladder well expected");
    assert_eq!(ladders[0].peaks.len(), 16, "ladder should have 16 peaks");
    assert!((ladders[0].peaks.last().unwrap().length - 6000.0).abs() < 0.5);
}

#[test]
fn fa_concentration_and_molarity_are_computed() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let run = traceio::fa::read_fa_run(&dir).expect("read FA run");

    let d1 = &run.samples[0];
    assert_eq!(d1.concentration.len(), d1.fluorescence.len());
    assert_eq!(d1.molarity.len(), d1.fluorescence.len());
    assert!(d1.concentration.iter().any(|v| v.is_finite() && *v > 0.0));
    assert!(d1.molarity.iter().any(|v| v.is_finite() && *v > 0.0));

    let sample_peak = &d1.peaks[1];
    assert!(
        sample_peak.concentration.is_finite() && sample_peak.concentration > 0.0,
        "D1 sample concentration should be finite, got {}",
        sample_peak.concentration
    );
    assert!(
        sample_peak.molarity.is_finite() && sample_peak.molarity > 0.0,
        "D1 sample molarity should be finite, got {}",
        sample_peak.molarity
    );
}

#[test]
fn path_api_reads_fa_metadata_from_fa_zip() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let zip_path = std::env::temp_dir().join(format!(
        "traceio_fa_path_api_{}_{}.fa.zip",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    zip_run_dir(&dir, &zip_path);

    let loaded = traceio::io::read_path_with_metadata(&zip_path).expect("path API reads FA zip");

    assert_eq!(loaded.loaded.source.identity, zip_path);
    assert!(matches!(
        loaded.loaded.source.format,
        traceio::io::TraceFormat::FragmentAnalyzerRun {
            entry: traceio::io::FragmentAnalyzerEntry::Zip
        }
    ));
    assert!(loaded.fa_metadata().is_some());
    assert!(!loaded.loaded.run.samples.is_empty());
    std::fs::remove_file(zip_path).unwrap();
}

#[test]
fn path_api_saves_fa_zip_sidecar_names() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let zip_path = std::env::temp_dir().join(format!(
        "traceio_fa_path_save_{}_{}.fa.zip",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    zip_run_dir(&dir, &zip_path);

    let mut loaded = traceio::io::read_path(&zip_path).expect("read path API");
    loaded.run.samples[0].name = "D1: public api rename".to_string();

    assert!(traceio::io::supports_save_path(&loaded, &zip_path));
    traceio::io::save_path(&loaded, &zip_path).expect("save path API");

    let reloaded = traceio::io::read_path(&zip_path).expect("reload saved FA zip");
    assert_eq!(reloaded.run.samples[0].name, "D1: public api rename");
    std::fs::remove_file(zip_path).unwrap();
}

#[test]
fn reads_fa_run_from_zip() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let zip_path = std::env::temp_dir().join("traceio_fa_read.zip");
    zip_run_dir(&dir, &zip_path);

    // A zipped run must decode identically to the folder.
    let from_dir = traceio::fa::read_fa_run(&dir).expect("read dir");
    let from_zip = traceio::fa::read_fa_run(&zip_path).expect("read zip");
    assert!(
        traceio::fa::is_fa_path(&zip_path),
        "zip should be recognized as FA"
    );
    assert_eq!(from_zip.samples.len(), from_dir.samples.len());
    for (a, b) in from_zip.samples.iter().zip(&from_dir.samples) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.peaks.len(), b.peaks.len());
        assert_eq!(a.fluorescence, b.fluorescence);
    }
    let _ = std::fs::remove_file(&zip_path);
}

#[test]
fn saves_renamed_names_into_zip() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    let zip_path = std::env::temp_dir().join("traceio_fa_save.zip");
    zip_run_dir(&dir, &zip_path);

    let mut run = traceio::fa::read_fa_run(&zip_path).expect("read zip");
    let orig_trace = run.samples[0].fluorescence.clone();
    run.samples[0].name = "D1: RENAMED_SAMPLE".to_string();
    traceio::fa::save_txt_names(&zip_path, &run).expect("save into zip");

    // Reload: the rename persists and the .raw-derived trace is untouched.
    let reloaded = traceio::fa::read_fa_run(&zip_path).expect("reload zip");
    assert!(
        reloaded.samples[0].name.contains("RENAMED_SAMPLE"),
        "got {}",
        reloaded.samples[0].name
    );
    assert_eq!(
        reloaded.samples[0].fluorescence, orig_trace,
        "traces must survive a save"
    );
    assert_eq!(reloaded.samples.len(), 12);
    let _ = std::fs::remove_file(&zip_path);
}

#[test]
fn opens_from_any_run_member_file() {
    let Some(dir) = fa_run_dir() else {
        eprintln!("skipping: fa_examples run not present");
        return;
    };
    // Locate the run's .raw and two sibling members.
    let file_with_ext = |ext: &str| -> PathBuf {
        std::fs::read_dir(&dir)
            .unwrap()
            .flatten()
            .map(|e| e.path())
            .find(|p| {
                p.extension()
                    .and_then(|e| e.to_str())
                    .is_some_and(|e| e.eq_ignore_ascii_case(ext))
            })
            .unwrap_or_else(|| panic!("no .{ext} in run dir"))
    };
    let raw = file_with_ext("raw");
    let pks = file_with_ext("PKS");
    let txt = file_with_ext("txt");

    // Any run member (or the folder) is recognized as FA...
    assert!(traceio::fa::is_fa_path(&pks), "a .PKS member should be FA");
    assert!(traceio::fa::is_fa_path(&txt), "a .txt member should be FA");
    assert!(traceio::fa::is_fa_path(&dir), "the run folder should be FA");

    // ...and all map to one identity (the .raw), so a multi-file drop dedups.
    let id = traceio::fa::run_identity(&raw);
    assert_eq!(traceio::fa::run_identity(&pks), id);
    assert_eq!(traceio::fa::run_identity(&txt), id);
    assert_eq!(traceio::fa::run_identity(&dir), id);

    // Opening via a member decodes the whole run.
    let via_member = traceio::fa::read_fa_run(&pks).expect("open via .PKS");
    assert_eq!(via_member.samples.len(), 12);

    // A Bioanalyzer file sitting in the folder must NOT be hijacked as FA.
    let stray_xml = dir.join("decoy.xml");
    std::fs::write(&stray_xml, b"<Chipset/>").unwrap();
    let hijacked = traceio::fa::is_fa_path(&stray_xml);
    let _ = std::fs::remove_file(&stray_xml);
    assert!(!hijacked, ".xml must be left to the Bioanalyzer reader");
}
