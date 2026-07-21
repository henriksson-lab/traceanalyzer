//! File export helpers for the GUI. These operate on already-loaded, in-memory
//! data and rendered plot buffers; instrument parser internals stay out of this
//! layer.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::traceio::xad::RawChannel;
use crate::traceio::{Electrophoresis, Sample};
use anyhow::{bail, Context};

use crate::plot::XAxis;
use crate::table;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceColumn {
    Time,
    Fluorescence,
    AlignedTime,
    Length,
    Concentration,
    Molarity,
}

impl TraceColumn {
    fn header(self) -> &'static str {
        match self {
            Self::Time => "time_s",
            Self::Fluorescence => "fluorescence",
            Self::AlignedTime => "aligned_time_s",
            Self::Length => "length",
            Self::Concentration => "concentration",
            Self::Molarity => "molarity",
        }
    }

    fn len(self, sample: &Sample) -> usize {
        match self {
            Self::Time => sample.time.len(),
            Self::Fluorescence => sample.fluorescence.len(),
            Self::AlignedTime => sample.aligned_time.len(),
            Self::Length => sample.length.len(),
            Self::Concentration => sample.concentration.len(),
            Self::Molarity => sample.molarity.len(),
        }
    }

    fn value(self, sample: &Sample, index: usize) -> Option<f64> {
        match self {
            Self::Time => sample.time.get(index).copied(),
            Self::Fluorescence => sample.fluorescence.get(index).map(|v| *v as f64),
            Self::AlignedTime => sample.aligned_time.get(index).copied(),
            Self::Length => sample.length.get(index).copied(),
            Self::Concentration => sample.concentration.get(index).copied(),
            Self::Molarity => sample.molarity.get(index).copied(),
        }
        .filter(|v| v.is_finite())
    }
}

pub const DEFAULT_TRACE_COLUMNS: &[TraceColumn] = &[
    TraceColumn::Time,
    TraceColumn::Fluorescence,
    TraceColumn::AlignedTime,
    TraceColumn::Length,
    TraceColumn::Concentration,
    TraceColumn::Molarity,
];

#[derive(Debug, Clone, Copy)]
pub struct TraceSample<'a> {
    /// Zero-based sample position in the source run.
    pub sample_index: usize,
    pub sample: &'a Sample,
}

/// Write long-format raw detector channel data from a native `.xad`.
pub fn write_raw_channels_csv(dst: &Path, channels: &[&RawChannel]) -> anyhow::Result<()> {
    if channels.is_empty() {
        bail!("raw channel export needs at least one channel");
    }

    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let mut out = BufWriter::new(file);

    write_csv_record(
        &mut out,
        [
            "channel_index",
            "channel_id",
            "channel_name",
            "point_index",
            "time_s",
            "signal",
        ],
    )?;

    for (channel_index, channel) in channels.iter().enumerate() {
        let prefix = [
            (channel_index + 1).to_string(),
            channel.channel_id.clone(),
            channel.name.clone(),
        ];
        for (point_index, signal) in channel.signal.iter().enumerate() {
            let time = channel.x_start + channel.x_step * point_index as f64;
            write_csv_record(
                &mut out,
                [
                    prefix[0].clone(),
                    prefix[1].clone(),
                    prefix[2].clone(),
                    (point_index + 1).to_string(),
                    csv_num(time),
                    csv_num(*signal as f64),
                ],
            )?;
        }
    }

    finish_writer(out, dst)
}

/// Write the focused sample's peak/region table as CSV.
pub fn write_peak_table_csv(
    dst: &Path,
    run: &Electrophoresis,
    sample: &Sample,
    x_axis: XAxis,
) -> anyhow::Result<()> {
    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let mut out = BufWriter::new(file);

    write_csv_record(&mut out, table::HEADERS)?;
    for row in table::rows_with_axis(run, sample, x_axis) {
        write_csv_record(&mut out, &row.cells)?;
    }
    out.flush()
        .with_context(|| format!("could not finish writing {}", dst.display()))?;
    Ok(())
}

/// Write every sample's peak/region table rows into one CSV.
pub fn write_run_peak_table_csv(
    dst: &Path,
    run: &Electrophoresis,
    x_axis: XAxis,
) -> anyhow::Result<()> {
    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let mut out = BufWriter::new(file);

    let mut headers = sample_headers();
    headers.extend(table::HEADERS.iter().map(|h| h.to_string()));
    write_csv_record(&mut out, &headers)?;

    for (sample_index, sample) in run.samples.iter().enumerate() {
        let prefix = sample_prefix(sample_index, sample);
        for row in table::rows_with_axis(run, sample, x_axis) {
            let mut fields = prefix.clone();
            fields.extend(row.cells);
            write_csv_record(&mut out, &fields)?;
        }
    }

    finish_writer(out, dst)
}

/// Write long-format trace data for selected samples and trace columns.
pub fn write_trace_data_csv(
    dst: &Path,
    samples: &[TraceSample<'_>],
    columns: &[TraceColumn],
) -> anyhow::Result<()> {
    if columns.is_empty() {
        bail!("trace export needs at least one column");
    }

    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let mut out = BufWriter::new(file);

    let mut headers = sample_headers();
    headers.push("point_index".to_string());
    headers.extend(columns.iter().map(|c| c.header().to_string()));
    write_csv_record(&mut out, &headers)?;

    for row in samples {
        let sample = row.sample;
        let max_len = columns.iter().map(|c| c.len(sample)).max().unwrap_or(0);
        let prefix = sample_prefix(row.sample_index, sample);
        for point_index in 0..max_len {
            let mut fields = prefix.clone();
            fields.push((point_index + 1).to_string());
            fields.extend(columns.iter().map(|c| {
                c.value(sample, point_index)
                    .map(csv_num)
                    .unwrap_or_default()
            }));
            write_csv_record(&mut out, &fields)?;
        }
    }

    finish_writer(out, dst)
}

/// Write normalized run metadata plus basic per-sample QC/summary metrics.
pub fn write_metadata_qc_csv(dst: &Path, run: &Electrophoresis) -> anyhow::Result<()> {
    write_metadata_qc_csv_with_notes(dst, run, None)
}

/// Write normalized metadata/QC plus optional source-specific notes, such as a
/// Fragment Analyzer sidecar summary shown in the GUI Metadata tab.
pub fn write_metadata_qc_csv_with_notes(
    dst: &Path,
    run: &Electrophoresis,
    notes: Option<&str>,
) -> anyhow::Result<()> {
    write_metadata_qc_csv_with_notes_and_provenance(dst, run, notes, false)
}

/// Write normalized metadata/QC plus optional source-specific notes and GUI
/// session provenance.
pub fn write_metadata_qc_csv_with_notes_and_provenance(
    dst: &Path,
    run: &Electrophoresis,
    notes: Option<&str>,
    marker_overrides_active: bool,
) -> anyhow::Result<()> {
    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let mut out = BufWriter::new(file);

    write_csv_record(
        &mut out,
        [
            "section",
            "sample_index",
            "well",
            "sample_name",
            "key",
            "value",
        ],
    )?;

    write_run_metric(&mut out, "file_name", &run.assay.file_name)?;
    write_run_metric(&mut out, "creation_date", &run.assay.creation_date)?;
    write_run_metric(&mut out, "assay_name", &run.assay.assay_name)?;
    write_run_metric(&mut out, "assay_type", &run.assay.assay_type)?;
    write_run_metric(&mut out, "length_unit", &run.assay.length_unit)?;
    write_run_metric(
        &mut out,
        "concentration_unit",
        &run.assay.concentration_unit,
    )?;
    write_run_metric(
        &mut out,
        "molarity_unit",
        run.assay.molarity_unit.as_deref().unwrap_or(""),
    )?;
    write_run_metric(
        &mut out,
        "has_upper_marker",
        &run.assay.has_upper_marker.to_string(),
    )?;
    write_run_metric(
        &mut out,
        "session_marker_overrides_active",
        &marker_overrides_active.to_string(),
    )?;
    if marker_overrides_active {
        write_run_metric(
            &mut out,
            "session_marker_overrides_provenance",
            "Manual marker overrides were active in the GUI session; recalibrated values may differ from source-file marker detection and overrides are not persisted by Save/Save As.",
        )?;
    }
    write_run_metric(&mut out, "sample_count", &run.samples.len().to_string())?;
    write_run_metric(
        &mut out,
        "ladder_peak_count",
        &run.ladder_peaks.len().to_string(),
    )?;
    write_run_metric(&mut out, "run_region_count", &run.regions.len().to_string())?;

    for (sample_index, sample) in run.samples.iter().enumerate() {
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "is_ladder",
            sample.is_ladder,
        )?;
        write_sample_metric(&mut out, sample_index, sample, "category", &sample.category)?;
        write_sample_metric(&mut out, sample_index, sample, "comment", &sample.comment)?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "observations",
            &sample.observations,
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "rin",
            sample.rin.map(csv_num).unwrap_or_default(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "trace_points",
            sample.fluorescence.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "aligned_time_points",
            sample.aligned_time.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "length_points",
            sample.length.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "concentration_points",
            sample.concentration.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "molarity_points",
            sample.molarity.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "peak_count",
            sample.peaks.len(),
        )?;
        write_sample_metric(
            &mut out,
            sample_index,
            sample,
            "region_count",
            sample.regions.len(),
        )?;
        if let Some((lo, hi)) = sample.fluorescence_range() {
            write_sample_metric(
                &mut out,
                sample_index,
                sample,
                "fluorescence_min",
                csv_num(lo as f64),
            )?;
            write_sample_metric(
                &mut out,
                sample_index,
                sample,
                "fluorescence_max",
                csv_num(hi as f64),
            )?;
        }
        if let Some((lo, hi)) = finite_range(&sample.time) {
            write_sample_metric(&mut out, sample_index, sample, "time_min_s", csv_num(lo))?;
            write_sample_metric(&mut out, sample_index, sample, "time_max_s", csv_num(hi))?;
        }
    }

    if let Some(notes) = notes.filter(|notes| !notes.trim().is_empty()) {
        for (line_index, line) in notes.lines().enumerate() {
            write_csv_record(
                &mut out,
                [
                    "source_metadata".to_string(),
                    String::new(),
                    String::new(),
                    String::new(),
                    (line_index + 1).to_string(),
                    line.to_string(),
                ],
            )?;
        }
    }

    finish_writer(out, dst)
}

/// Encode an RGB buffer (`width * height * 3` bytes) as an 8-bit PNG.
pub fn write_rgb_png(dst: &Path, rgb: &[u8], width: u32, height: u32) -> anyhow::Result<()> {
    let expected = width as usize * height as usize * 3;
    if rgb.len() != expected {
        bail!(
            "plot buffer has {} bytes, expected {} for {}x{} RGB",
            rgb.len(),
            expected,
            width,
            height
        );
    }

    let file = File::create(dst).with_context(|| format!("could not create {}", dst.display()))?;
    let writer = BufWriter::new(file);
    let mut encoder = png::Encoder::new(writer, width, height);
    encoder.set_color(png::ColorType::Rgb);
    encoder.set_depth(png::BitDepth::Eight);
    let mut png = encoder
        .write_header()
        .with_context(|| format!("could not write PNG header to {}", dst.display()))?;
    png.write_image_data(rgb)
        .with_context(|| format!("could not write PNG data to {}", dst.display()))?;
    Ok(())
}

fn finish_writer(mut out: BufWriter<File>, dst: &Path) -> anyhow::Result<()> {
    out.flush()
        .with_context(|| format!("could not finish writing {}", dst.display()))?;
    Ok(())
}

fn sample_headers() -> Vec<String> {
    [
        "sample_index",
        "well",
        "sample_name",
        "category",
        "is_ladder",
    ]
    .iter()
    .map(|h| h.to_string())
    .collect()
}

fn sample_prefix(sample_index: usize, sample: &Sample) -> Vec<String> {
    vec![
        (sample_index + 1).to_string(),
        sample.well_number.to_string(),
        sample.name.clone(),
        sample.category.clone(),
        sample.is_ladder.to_string(),
    ]
}

fn csv_num(v: f64) -> String {
    if v.is_finite() {
        v.to_string()
    } else {
        String::new()
    }
}

fn finite_range(values: &[f64]) -> Option<(f64, f64)> {
    let mut it = values.iter().copied().filter(|v| v.is_finite());
    let first = it.next()?;
    let (mut lo, mut hi) = (first, first);
    for v in it {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    Some((lo, hi))
}

fn write_run_metric<W: Write>(out: &mut W, key: &str, value: &str) -> std::io::Result<()> {
    write_csv_record(out, ["run", "", "", "", key, value])
}

fn write_sample_metric<W, V>(
    out: &mut W,
    sample_index: usize,
    sample: &Sample,
    key: &str,
    value: V,
) -> std::io::Result<()>
where
    W: Write,
    V: ToString,
{
    write_csv_record(
        out,
        [
            "sample".to_string(),
            (sample_index + 1).to_string(),
            sample.well_number.to_string(),
            sample.name.clone(),
            key.to_string(),
            value.to_string(),
        ],
    )
}

fn write_csv_record<W, I, S>(out: &mut W, fields: I) -> std::io::Result<()>
where
    W: Write,
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let mut first = true;
    for field in fields {
        if first {
            first = false;
        } else {
            out.write_all(b",")?;
        }
        write_csv_field(out, field.as_ref())?;
    }
    out.write_all(b"\n")
}

fn write_csv_field<W: Write>(out: &mut W, field: &str) -> std::io::Result<()> {
    let needs_quotes = field
        .as_bytes()
        .iter()
        .any(|b| matches!(b, b',' | b'"' | b'\n' | b'\r'));
    if !needs_quotes {
        out.write_all(field.as_bytes())?;
        return Ok(());
    }

    out.write_all(b"\"")?;
    for b in field.bytes() {
        if b == b'"' {
            out.write_all(b"\"\"")?;
        } else {
            out.write_all(&[b])?;
        }
    }
    out.write_all(b"\"")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traceio::xad::RawChannel;
    use crate::traceio::{AssayInfo, Peak, Region};

    fn sample(well_number: i32, name: &str) -> Sample {
        Sample {
            well_number,
            name: name.to_string(),
            category: "Sample".to_string(),
            is_ladder: false,
            comment: "note".to_string(),
            observations: String::new(),
            rin: Some(8.5),
            time: vec![0.0, 1.0, 2.0],
            fluorescence: vec![10.0, 20.0, 15.0],
            aligned_time: vec![0.1, 1.1, 2.1],
            length: vec![100.0, 200.0, f64::NAN],
            concentration: vec![f64::NAN, 2.5, 3.5],
            molarity: vec![f64::NAN, 4.5, 5.5],
            peaks: vec![Peak {
                observations: "main".to_string(),
                length: 200.0,
                time: 1.0,
                aligned_time: 1.1,
                start_time: 0.8,
                end_time: 1.2,
                aligned_start_time: 0.9,
                aligned_end_time: 1.3,
                area: 50.0,
                concentration: 2.5,
                molarity: 4.5,
            }],
            regions: vec![Region {
                lower_length: 150.0,
                upper_length: 250.0,
            }],
        }
    }

    fn run() -> Electrophoresis {
        Electrophoresis {
            assay: AssayInfo {
                file_name: "demo.xml".to_string(),
                creation_date: "2024-01-02".to_string(),
                assay_name: "DNA".to_string(),
                assay_type: "DNA".to_string(),
                length_unit: "bp".to_string(),
                concentration_unit: "ng/ul".to_string(),
                molarity_unit: Some("nM".to_string()),
                has_upper_marker: true,
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![sample(1, "A1"), sample(2, "B1, quoted")],
        }
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "traceanalyzer_export_test_{}_{}_{}.csv",
            name,
            std::process::id(),
            std::thread::current().name().unwrap_or("thread")
        ));
        path
    }

    fn raw_channel(channel_id: &str, name: &str, x_start: f64, x_step: f64) -> RawChannel {
        RawChannel {
            channel_id: channel_id.to_string(),
            name: name.to_string(),
            x_start,
            x_step,
            signal: vec![1.0, f32::NAN, 3.5],
        }
    }

    fn trace_sample(sample_index: usize, sample: &Sample) -> TraceSample<'_> {
        TraceSample {
            sample_index,
            sample,
        }
    }

    #[test]
    fn csv_fields_escape_commas_quotes_and_newlines() {
        let mut out = Vec::new();
        write_csv_record(&mut out, ["plain", "a,b", "say \"yes\"", "two\nlines"]).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "plain,\"a,b\",\"say \"\"yes\"\"\",\"two\nlines\"\n"
        );
    }

    #[test]
    fn run_peak_table_includes_sample_identity_for_all_samples() {
        let dst = temp_path("run_peaks");
        write_run_peak_table_csv(&dst, &run(), XAxis::Length).unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text.starts_with(
            "sample_index,well,sample_name,category,is_ladder,#,size,time (s),area,height,% total,conc,molarity,note\n"
        ));
        assert!(text.contains("1,1,A1,Sample,false,1,200,1.0,50.0,20.0,100.0,2.50,4.50,main\n"));
        assert!(text.contains("2,2,\"B1, quoted\",Sample,false,R1,150"));
    }

    #[test]
    fn trace_data_export_writes_selected_columns_and_blank_missing_values() {
        let run = run();
        let samples = [&run.samples[0]];
        let dst = temp_path("trace_data");

        write_trace_data_csv(
            &dst,
            &[trace_sample(0, samples[0])],
            &[
                TraceColumn::Time,
                TraceColumn::Fluorescence,
                TraceColumn::Length,
            ],
        )
        .unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert_eq!(
            text,
            "sample_index,well,sample_name,category,is_ladder,point_index,time_s,fluorescence,length\n\
             1,1,A1,Sample,false,1,0,10,100\n\
             1,1,A1,Sample,false,2,1,20,200\n\
             1,1,A1,Sample,false,3,2,15,\n"
        );
    }

    #[test]
    fn trace_data_export_preserves_original_sample_indices() {
        let run = run();
        let dst = temp_path("trace_original_indices");

        write_trace_data_csv(
            &dst,
            &[
                trace_sample(1, &run.samples[1]),
                trace_sample(0, &run.samples[0]),
            ],
            &[TraceColumn::Time],
        )
        .unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text
            .starts_with("sample_index,well,sample_name,category,is_ladder,point_index,time_s\n"));
        assert!(text.contains("2,2,\"B1, quoted\",Sample,false,1,0\n"));
        assert!(text.contains("1,1,A1,Sample,false,1,0\n"));
    }

    #[test]
    fn raw_channel_export_writes_time_and_blank_non_finite_signal() {
        let dst = temp_path("raw_channels");
        let ch = raw_channel("BlueFluorescence", "Blue, detector", 0.5, 0.25);

        write_raw_channels_csv(&dst, &[&ch]).unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert_eq!(
            text,
            "channel_index,channel_id,channel_name,point_index,time_s,signal\n\
             1,BlueFluorescence,\"Blue, detector\",1,0.5,1\n\
             1,BlueFluorescence,\"Blue, detector\",2,0.75,\n\
             1,BlueFluorescence,\"Blue, detector\",3,1,3.5\n"
        );
    }

    #[test]
    fn raw_channel_export_writes_all_supplied_channels() {
        let dst = temp_path("raw_all_channels");
        let blue = raw_channel("BlueFluorescence", "Blue detector", 0.5, 0.25);
        let red = raw_channel("RedFluorescence", "Red detector", 1.0, 0.5);

        write_raw_channels_csv(&dst, &[&blue, &red]).unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text.contains("1,BlueFluorescence,Blue detector,1,0.5,1\n"));
        assert!(text.contains("2,RedFluorescence,Red detector,1,1,1\n"));
        assert!(text.contains("2,RedFluorescence,Red detector,3,2,3.5\n"));
    }

    #[test]
    fn raw_channel_export_rejects_empty_channel_selection() {
        let dst = temp_path("raw_empty_channels");
        let err = write_raw_channels_csv(&dst, &[]).unwrap_err();
        assert!(err.to_string().contains("at least one channel"));
    }

    #[test]
    fn metadata_qc_export_includes_run_and_sample_metrics() {
        let dst = temp_path("metadata_qc");
        write_metadata_qc_csv(&dst, &run()).unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text.contains("run,,,,file_name,demo.xml\n"));
        assert!(text.contains("run,,,,has_upper_marker,true\n"));
        assert!(text.contains("sample,1,1,A1,rin,8.5\n"));
        assert!(text.contains("sample,1,1,A1,fluorescence_min,10\n"));
        assert!(text.contains("sample,1,1,A1,time_max_s,2\n"));
    }

    #[test]
    fn metadata_qc_export_can_include_source_metadata_notes() {
        let dst = temp_path("metadata_qc_notes");
        write_metadata_qc_csv_with_notes(&dst, &run(), Some("Method\n[Separation]\nKV: 7.00"))
            .unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text.contains("source_metadata,,,,1,Method\n"));
        assert!(text.contains("source_metadata,,,,3,KV: 7.00\n"));
    }

    #[test]
    fn metadata_qc_export_flags_session_marker_overrides() {
        let dst = temp_path("metadata_qc_marker_overrides");
        write_metadata_qc_csv_with_notes_and_provenance(&dst, &run(), None, true).unwrap();

        let text = std::fs::read_to_string(&dst).unwrap();
        let _ = std::fs::remove_file(&dst);

        assert!(text.contains("run,,,,session_marker_overrides_active,true\n"));
        assert!(text.contains("run,,,,session_marker_overrides_provenance,"));
    }

    #[test]
    fn trace_data_export_rejects_empty_column_selection() {
        let dst = temp_path("trace_empty_columns");
        let run = run();
        let err = write_trace_data_csv(&dst, &[trace_sample(0, &run.samples[0])], &[]).unwrap_err();
        assert!(err.to_string().contains("at least one column"));
    }
}
