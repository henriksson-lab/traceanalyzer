use traceanalyzer::plot::{self, YMode};
use traceio::{AssayInfo, Electrophoresis, Sample};

fn demo_run() -> Electrophoresis {
    let sample = Sample {
        well_number: 1,
        name: "A1".to_string(),
        category: "Sample".to_string(),
        is_ladder: false,
        comment: String::new(),
        observations: String::new(),
        rin: None,
        time: vec![0.0, 1.0, 2.0, 3.0, 4.0],
        fluorescence: vec![0.0, 4.0, 16.0, 4.0, 0.0],
        aligned_time: Vec::new(),
        length: Vec::new(),
        concentration: Vec::new(),
        molarity: Vec::new(),
        peaks: Vec::new(),
    };

    Electrophoresis {
        assay: AssayInfo {
            file_name: "synthetic".to_string(),
            creation_date: String::new(),
            assay_name: "Synthetic".to_string(),
            assay_type: "DNA".to_string(),
            length_unit: "bp".to_string(),
            concentration_unit: "ng/ul".to_string(),
            molarity_unit: Some("nM".to_string()),
            has_upper_marker: false,
        },
        ladder_peaks: Vec::new(),
        regions: Vec::new(),
        samples: vec![sample],
    }
}

fn count_primary_trace_pixels(pixels: &[u8]) -> usize {
    pixels
        .chunks_exact(3)
        .filter(|rgb| {
            let [r, g, b] = **rgb else { unreachable!() };
            r < 80 && (90..=160).contains(&g) && b > 140
        })
        .count()
}

#[test]
fn render_rgb_has_expected_dimensions_and_draws_trace_series() {
    let run = demo_run();
    let series = plot::series(&run, &run.samples[0], YMode::Fluorescence, false);
    let viewport = plot::auto_viewport(&series);
    let width = 320;
    let height = 160;

    let pixels = plot::render_rgb(&series, &viewport, width, height);

    assert_eq!(pixels.len(), (width * height * 3) as usize);
    assert!(
        count_primary_trace_pixels(&pixels) > 20,
        "rendered plot should contain the blue trace series, not only axes/grid/text"
    );
}
