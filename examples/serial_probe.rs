//! Live serial connect probe (milestones M6/M7). Requires the `serial` feature:
//!
//!   # list candidate ports:
//!   cargo run --features serial --example serial_probe
//!
//!   # read-only connect: open a port and print framed records for ~5 s:
//!   cargo run --features serial --example serial_probe -- <port> [baud]
//!
//!   # guarded run (prompts before each energizing step): add --run
//!   cargo run --features serial --example serial_probe -- <port> <baud> --run
//!
//! Without a physical instrument this connects to whatever serial device is
//! given and simply reads nothing — the connect + framing path is what it
//! exercises.

#[cfg(not(feature = "serial"))]
fn main() {
    eprintln!("rebuild with `--features serial` to use this example");
}

#[cfg(feature = "serial")]
fn main() -> anyhow::Result<()> {
    use std::io::Write;
    use std::time::{Duration, Instant};
    use traceanalyzer::tracehw::protocol::decode;
    use traceanalyzer::tracehw::safety::{Action, ApproveWith};
    use traceanalyzer::tracehw::transport::serial;
    use traceanalyzer::tracehw::Transport;

    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        let ports = serial::list_ports();
        if ports.is_empty() {
            println!("no serial ports found");
        } else {
            println!("serial ports:");
            for p in ports {
                println!("  {p}");
            }
        }
        println!("\npass a port name to connect: serial_probe <port> [baud] [--run]");
        return Ok(());
    }

    let port = &args[0];
    let baud = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(serial::DEFAULT_BAUD);
    let run = args.iter().any(|a| a == "--run");

    let mut transport = serial::open(port, baud)?;
    println!("connected to {port} at {baud} baud");

    if run {
        // Console approver embodying the safety gate: the operator must confirm
        // each energizing action.
        let mut approver = ApproveWith(|action: Action| {
            print!(
                "About to {} — this energizes the chip. Type 'yes' to proceed: ",
                match action {
                    Action::Setup => "AV.SETUP",
                    Action::Start => "START (apply voltage)",
                }
            );
            std::io::stdout().flush().ok();
            let mut line = String::new();
            std::io::stdin().read_line(&mut line).ok();
            line.trim().eq_ignore_ascii_case("yes")
        });
        let acq =
            serial::guarded_run_live(&mut transport, &mut approver, Duration::from_secs(3600))?;
        println!(
            "run finished: state {:?}, {} samples",
            acq.state,
            acq.sample_count()
        );
        return Ok(());
    }

    // Read-only probe: print framed records for a few seconds.
    println!("reading (read-only) for 5 s…");
    let deadline = Instant::now() + Duration::from_secs(5);
    let mut count = 0;
    while Instant::now() < deadline {
        match transport.poll()? {
            Some(rec) => {
                count += 1;
                if count <= 20 {
                    println!("  {:?}", decode(&rec));
                }
            }
            None => std::thread::sleep(Duration::from_millis(20)),
        }
    }
    println!("received {count} records");
    Ok(())
}
