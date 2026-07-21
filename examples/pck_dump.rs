//! Decode a 2100 Expert `.pck` recorded session into a human timeline.
//!
//! Usage:
//!   cargo run --example pck_dump -- <session.pck> [--samples N]
//!
//! Prints the header, a record-type census, checksum health, and an event
//! timeline (text commands/responses/events, with runs of high-rate sample
//! records collapsed). With `--samples N` it also dumps the first `N` records of
//! each high-rate stream as a provisional (Phase-P0) numeric decode.

use anyhow::{bail, Context, Result};
use traceanalyzer::tracehw::pck::{type_label, Session, StreamEnd};

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let mut path: Option<String> = None;
    let mut samples = 0usize;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--samples" => {
                samples = args
                    .next()
                    .context("--samples needs a number")?
                    .parse()
                    .context("--samples value must be an integer")?;
            }
            _ if path.is_none() => path = Some(a),
            other => bail!("unexpected argument: {other}"),
        }
    }
    let path = path.context("usage: pck_dump <session.pck> [--samples N]")?;

    let session = Session::read(&path)?;

    println!("== {path} ==");
    print_header(&session);
    print_census(&session);
    print_timeline(&session);
    if samples > 0 {
        print_samples(&session, samples);
    }
    Ok(())
}

fn print_header(s: &Session) {
    println!("\n-- header ({} bytes) --", s.header.raw.len());
    if !s.header.is_text {
        println!("  (non-text header — likely a diagnostics capture)");
        return;
    }
    if let Some(fw) = &s.header.firmware_version {
        println!("  firmware: {fw}");
    }
    if let Some(assay) = &s.header.assay_path {
        println!("  assay:    {assay}");
    }
    for (i, f) in s.header.fields.iter().enumerate().skip(2) {
        if !f.is_empty() {
            println!("  field[{i}]: {f}");
        }
    }
}

fn print_census(s: &Session) {
    println!("\n-- record census ({} records) --", s.records.len());
    for (type_id, count) in s.type_census() {
        println!("  0x{type_id:02x}  {count:>7}  {}", type_label(type_id));
    }
    let bad = s.checksum_failures();
    println!(
        "  checksums: {} ({} failed)",
        if bad == 0 { "all OK" } else { "PROBLEM" },
        bad
    );
    match s.end {
        StreamEnd::Clean => {
            let t = String::from_utf8_lossy(&s.trailer);
            println!(
                "  stream end: clean, trailer {:?}",
                t.trim_end_matches('\0')
            );
        }
        StreamEnd::Truncated => println!("  stream end: TRUNCATED (log cut off mid-record)"),
    }
}

/// Print the ordered events, collapsing runs of high-rate binary sample records
/// (types 0x04/0x06 and their aux companions) into a single summary line.
fn print_timeline(s: &Session) {
    println!("\n-- timeline --");
    let mut run = 0usize;
    let flush = |run: &mut usize| {
        if *run > 0 {
            println!("  … {} sample/telemetry records …", run);
            *run = 0;
        }
    };
    for r in &s.records {
        match r.as_text() {
            Some(text) => {
                flush(&mut run);
                let flag = if r.checksum_ok {
                    ""
                } else {
                    "  [BAD CHECKSUM]"
                };
                println!("  t0x{:02x} seq{:02x}  {}{}", r.type_id, r.seq, text, flag);
            }
            None => run += 1,
        }
    }
    flush(&mut run);
}

/// Provisional numeric decode of the two high-rate streams. The exact meaning of
/// these bytes is protocol Phase P2 — this just surfaces the raw structure.
fn print_samples(s: &Session, n: usize) {
    for type_id in [0x04u8, 0x06u8] {
        let recs: Vec<_> = s
            .records
            .iter()
            .filter(|r| r.type_id == type_id && r.as_text().is_none())
            .take(n)
            .collect();
        if recs.is_empty() {
            continue;
        }
        println!(
            "\n-- first {} records of stream 0x{type_id:02x} --",
            recs.len()
        );
        for r in &recs {
            let hex: Vec<String> = r.payload.iter().map(|b| format!("{b:02x}")).collect();
            // First-cut fields: a 24-bit little-endian value from the last 3
            // payload bytes (see docs/pck_format.md).
            let le24 = if r.payload.len() >= 3 {
                let p = &r.payload[r.payload.len() - 3..];
                Some(u32::from(p[0]) | (u32::from(p[1]) << 8) | (u32::from(p[2]) << 16))
            } else {
                None
            };
            println!(
                "  seq{:02x}  [{}]  le24_tail={}",
                r.seq,
                hex.join(" "),
                le24.map_or("-".into(), |v| v.to_string())
            );
        }
    }
}
