//! Safety gating for real runs (plan §5).
//!
//! Driving the instrument applies **high voltage, a laser and a heater** to a
//! physical chip — every energizing action is outward-facing and irreversible.
//! [`guarded_run`] therefore requires explicit [`RunApproval`] before each such
//! action and aborts (sending `STOP`) the moment approval is withheld. The safe
//! default is [`DenyAll`]; [`ApproveAll`] is for the simulator and tests only.
//!
//! **Status:** this gating layer is implemented and tested against the
//! [`PckReplay`](crate::tracehw::transport::PckReplay) simulator. Running it against a
//! live serial transport (M7) additionally needs the outbound command framing
//! validated on hardware — see [`StreamTransport`](crate::tracehw::transport::StreamTransport).
//! Firmware download is intentionally out of scope.

use anyhow::{bail, Result};

use crate::tracehw::instrument::{Acquisition, RunState};
use crate::tracehw::protocol::{self, Command};
use crate::tracehw::transport::Transport;

/// A physical action the controller is about to take.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    /// `AV.SETUP` — loads/energizes the run configuration.
    Setup,
    /// `START` — begins electrophoresis (applies voltage). The critical gate.
    Start,
}

/// Approves (or refuses) each energizing action. Implementors typically prompt
/// the user; the run aborts if any action is refused.
pub trait RunApproval {
    /// Return `true` to allow `action`, `false` to abort the run.
    fn approve(&mut self, action: Action) -> bool;
}

/// Refuses every action — the safe default when no approver is wired up.
pub struct DenyAll;
impl RunApproval for DenyAll {
    fn approve(&mut self, _action: Action) -> bool {
        false
    }
}

/// Approves every action. **Simulator/tests only** — never for real hardware.
pub struct ApproveAll;
impl RunApproval for ApproveAll {
    fn approve(&mut self, _action: Action) -> bool {
        true
    }
}

/// Approves from a closure (e.g. a GUI confirmation dialog).
pub struct ApproveWith<F>(pub F);
impl<F: FnMut(Action) -> bool> RunApproval for ApproveWith<F> {
    fn approve(&mut self, action: Action) -> bool {
        (self.0)(action)
    }
}

/// Run the guarded handshake over `transport`, gating each energizing step
/// through `approval`. `AV.CLEAR` (non-energizing prep) runs first; `AV.SETUP`
/// and `START` each require approval. If either is refused, `STOP` is sent and
/// the run aborts with an error. On approval it drives acquisition until
/// `END_OF_RUN`, then sends `STOP`.
///
/// Works against the simulator today; a live run also needs the transport's
/// outbound framing hardware-validated (M7).
pub fn guarded_run(
    transport: &mut impl Transport,
    approval: &mut impl RunApproval,
) -> Result<Acquisition> {
    transport.send(&Command::AvClear)?;

    if !approval.approve(Action::Setup) {
        transport.send(&Command::Stop)?;
        bail!("run aborted: AV.SETUP was not approved");
    }
    transport.send(&Command::AvSetup)?;

    if !approval.approve(Action::Start) {
        transport.send(&Command::Stop)?;
        bail!("run aborted: START was not approved");
    }
    transport.send(&Command::Start)?;

    let mut acq = Acquisition::new();
    while let Some(record) = transport.poll()? {
        acq.apply(protocol::decode(&record));
        if acq.state == RunState::Ended {
            break;
        }
    }

    transport.send(&Command::Stop)?;
    Ok(acq)
}
