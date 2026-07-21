//! `tracehw` — Bioanalyzer 2100 instrument protocol.
//!
//! This module holds everything specific to *talking to* the instrument, kept
//! separate from the format readers in `traceio` and the GUI in
//! `traceanalyzer`. It is built out in the phases of
//! [`docs/hardware_control_plan.md`](../../docs/hardware_control_plan.md).
//!
//! Landed so far:
//!
//! * [`pck`] — reader for 2100 Expert's recorded-session `.pck` captures
//!   (protocol Phase P0). The reverse-engineering corpus for the wire protocol;
//!   needs no hardware.
//! * [`protocol`] — the command set and record→[`Event`](protocol::Event) decode.
//! * [`transport`] — the [`Transport`](transport::Transport) trait and the
//!   [`PckReplay`](transport::PckReplay) `.pck`-replay simulator.
//! * [`instrument`] — the run state machine ([`Acquisition`](instrument::Acquisition))
//!   and a driver that turns a replayed session into an
//!   [`Electrophoresis`](crate::traceio::Electrophoresis).
//! * [`safety`] — approval-gated [`guarded_run`](safety::guarded_run) (plan §5).
//!
//! The transport is a virtual **serial/COM port** (Phase P1/M6 finding).
//! [`StreamTransport`](transport::StreamTransport) frames records over any byte
//! stream; with the `serial` feature, [`transport::serial`] opens a real port
//! (round-trip tested over a pseudo-terminal) and
//! [`serial::guarded_run_live`](transport::serial) drives a gated live run. The
//! only step left is validating the outbound framing and baud against a physical
//! 2100 — pure hardware.

pub mod instrument;
pub mod pck;
pub mod protocol;
pub mod safety;
pub mod transport;

pub use instrument::{run_to_completion, Acquisition, RunState};
pub use pck::{Header, Record, Session, StreamEnd};
pub use protocol::{Command, Event};
pub use safety::{guarded_run, Action, ApproveAll, DenyAll, RunApproval};
pub use transport::{PckReplay, StreamTransport, Transport};
