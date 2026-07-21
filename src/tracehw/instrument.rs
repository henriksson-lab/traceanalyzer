//! Run control: a state machine that consumes protocol [`Event`]s and grows an
//! acquisition, plus a driver that pumps a [`Transport`] to completion.
//!
//! [`Acquisition::apply`] takes one event at a time, so the same type backs both
//! the headless simulator (M3) and the live GUI panel (M5), which will feed it
//! events as they arrive and re-render the growing trace. Converting the
//! collected samples into an [`Electrophoresis`] is [`Acquisition::to_electrophoresis`];
//! the fluorescence values are provisional until protocol Phase P2.

use crate::traceio::{AssayInfo, Electrophoresis, Sample};
use anyhow::Result;

use crate::tracehw::protocol::{self, Command, Event};
use crate::tracehw::transport::Transport;

/// Where a run is in its lifecycle. Advances on the acknowledged commands and
/// the `END_OF_RUN` event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunState {
    /// Before any command has been acknowledged.
    Idle,
    /// `AV.CLEAR` acknowledged.
    Cleared,
    /// `AV.SETUP` acknowledged.
    SetUp,
    /// `START` acknowledged — acquisition in progress.
    Running,
    /// `END_OF_RUN` seen — acquisition finished, awaiting `STOP`.
    Ended,
    /// `STOP` acknowledged — run complete.
    Stopped,
}

/// A growing acquisition assembled from protocol events.
#[derive(Debug, Clone)]
pub struct Acquisition {
    /// Current lifecycle state.
    pub state: RunState,
    /// Provisional raw values from stream A (record type `0x04`).
    pub stream_a: Vec<i64>,
    /// Provisional raw values from stream B (record type `0x06`).
    pub stream_b: Vec<i64>,
    /// Command acknowledgements / responses seen, in order.
    pub acks: Vec<String>,
    /// `GET.ERROR` status words seen, in order.
    pub errors: Vec<String>,
    /// Count of records that failed their checksum.
    pub corrupt: usize,
}

impl Default for Acquisition {
    fn default() -> Self {
        Acquisition {
            state: RunState::Idle,
            stream_a: Vec::new(),
            stream_b: Vec::new(),
            acks: Vec::new(),
            errors: Vec::new(),
            corrupt: 0,
        }
    }
}

impl Acquisition {
    /// A fresh, idle acquisition.
    pub fn new() -> Self {
        Self::default()
    }

    /// Fold one protocol event into the acquisition, advancing state and
    /// appending samples as appropriate.
    pub fn apply(&mut self, event: Event) {
        match event {
            Event::Ack(line) => {
                self.state = match () {
                    _ if line.starts_with("AV.CLEAR") => RunState::Cleared,
                    _ if line.starts_with("AV.SETUP") => RunState::SetUp,
                    _ if line.starts_with("START") => RunState::Running,
                    _ if line.starts_with("STOP") => RunState::Stopped,
                    _ => self.state,
                };
                self.acks.push(line);
            }
            Event::Error(code) => {
                if !code.is_empty() {
                    self.errors.push(code);
                }
            }
            Event::EndOfRun => {
                if self.state == RunState::Running {
                    self.state = RunState::Ended;
                }
            }
            Event::Sample(s) => match s.stream {
                0x04 => self.stream_a.push(s.raw),
                0x06 => self.stream_b.push(s.raw),
                _ => {}
            },
            Event::Telemetry { .. } => {}
            Event::Corrupt { .. } => self.corrupt += 1,
        }
    }

    /// Number of samples collected on the primary (stream B) trace so far.
    pub fn sample_count(&self) -> usize {
        self.stream_b.len()
    }

    /// Build a single-trace [`Electrophoresis`] from the collected samples.
    ///
    /// The continuous acquisition is returned as one [`Sample`] (well 1), the
    /// way the raw `.xad` reader represents an un-split run. `sample_period_s`
    /// sets seconds-per-point on the time axis (pass `1.0` to fall back to a
    /// sample-index axis). Fluorescence is the provisional stream-B value —
    /// uncalibrated until Phase P2.
    pub fn to_electrophoresis(&self, assay_name: &str, sample_period_s: f64) -> Electrophoresis {
        let fluorescence: Vec<f32> = self.stream_b.iter().map(|&v| v as f32).collect();
        let time: Vec<f64> = (0..fluorescence.len())
            .map(|i| i as f64 * sample_period_s)
            .collect();

        let sample = Sample {
            well_number: 1,
            name: "Live acquisition".to_string(),
            category: "Sample".to_string(),
            is_ladder: false,
            comment: String::new(),
            observations: String::new(),
            rin: None,
            time,
            fluorescence,
            aligned_time: Vec::new(),
            length: Vec::new(),
            concentration: Vec::new(),
            molarity: Vec::new(),
            peaks: Vec::new(),
            regions: Vec::new(),
        };

        Electrophoresis {
            assay: AssayInfo {
                assay_name: assay_name.to_string(),
                assay_type: "Bioanalyzer live (simulated)".to_string(),
                ..AssayInfo::default()
            },
            ladder_peaks: Vec::new(),
            regions: Vec::new(),
            samples: vec![sample],
        }
    }
}

/// Drive a transport to completion, folding every inbound record into an
/// [`Acquisition`]. For the [`PckReplay`](crate::tracehw::transport::PckReplay) simulator
/// this runs the whole recorded session; for a live transport it would loop
/// until `STOP`/disconnect. The controller sends the run handshake first.
pub fn run_to_completion(transport: &mut impl Transport) -> Result<Acquisition> {
    // Issue the run handshake. For a replay these are captured but do not gate
    // the (pre-recorded) inbound stream; a live transport would await each ack.
    for cmd in [Command::AvClear, Command::AvSetup, Command::Start] {
        transport.send(&cmd)?;
    }

    let mut acq = Acquisition::new();
    while let Some(record) = transport.poll()? {
        acq.apply(protocol::decode(&record));
    }

    transport.send(&Command::Stop)?;
    Ok(acq)
}
