//! File export helpers for the GUI. These operate on already-loaded, in-memory
//! data and rendered plot buffers; instrument parser internals stay out of this
//! layer.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context};
use traceio::{Electrophoresis, Sample};

use crate::plot::XAxis;
use crate::table;

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

    #[test]
    fn csv_fields_escape_commas_quotes_and_newlines() {
        let mut out = Vec::new();
        write_csv_record(&mut out, ["plain", "a,b", "say \"yes\"", "two\nlines"]).unwrap();

        assert_eq!(
            String::from_utf8(out).unwrap(),
            "plain,\"a,b\",\"say \"\"yes\"\"\",\"two\nlines\"\n"
        );
    }
}
