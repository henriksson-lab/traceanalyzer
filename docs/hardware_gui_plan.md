# Instrument-control GUI — implementation plan

A live **Instrument** view for driving a run (real or simulated) from the
Slint GUI, built on the existing `tracehw` backend. Companion to
`docs/gui_plan.md` (file-viewer plan) and the `tracehw` module docs.

## Decisions

- **Placement:** a dedicated **"Instrument" view mode**, toggled from the top
  bar like Plot / Overview / Gel. The full content area becomes the instrument
  panel; it persists across a run, so the user can flip to an already-open
  file's Plot and back without interrupting acquisition.
- **v1 transport scope:** wire the whole UX against the **`.pck`-replay
  simulator** (`PckReplay`) only. No hardware, builds **without** the `serial`
  feature, fully headless-testable. The serial transport is a drop-in
  follow-up behind the existing feature flag (same `Transport` trait).
- **Non-blocking:** acquisition runs on a **background thread**; the UI thread
  only renders periodic snapshots. A multi-minute run must never freeze the UI.
- **Safety is non-negotiable:** every energizing action (`AV.SETUP`, `START`)
  is gated by an explicit modal approval, per `tracehw::safety`. Default =
  abort.
- **Provisional signal:** the live trace is raw stream-B, uncalibrated until
  protocol Phase P2 — always shown behind a "PROVISIONAL" banner.

## What already exists (backend — do not rebuild)

| Piece | Where | Role in this view |
|---|---|---|
| `Acquisition` + `RunState` | `tracehw/instrument.rs` | Folds `Event`s via `apply()`; live counters `sample_count()`, `acks`, `errors`, `corrupt`; `to_electrophoresis()` for the growing trace. |
| `Transport` trait | `tracehw/transport.rs` | `PckReplay` (sim) and `serial::SerialTransport` are interchangeable. |
| `RunApproval` / `ApproveWith` | `tracehw/safety.rs` | `ApproveWith(closure)` is the GUI-dialog bridge; `guarded_run` gates Setup/Start. |
| `serial::list_ports / open / guarded_run_live` | `tracehw/transport/serial.rs` | The follow-up serial path. `open` is read-only-safe. |
| `simulate_pck` + `simulate_added_file` | `loading.rs`, `main.rs` | Today's synchronous path; the new view is the incremental version of it. |

The `Acquisition` docstring already names this panel as its intended second
consumer ("the live GUI panel, which will feed it events as they arrive and
re-render the growing trace") — the state machine takes events one at a time
precisely so it can be pumped incrementally.

## Layout

```
[ Plot ] [ Overview ] [ Gel ] [ ●Instrument ]
────────────────────────────────────────────────────────────────────
 Transport: [Simulator (.pck) ▾]   Port: [—]  Baud:[115200] ↻
 [ Connect ]        ● Connected · read-only
────────────────────────────────────────────────────────────────────
 Lifecycle:  Idle ─ Cleared ─ SetUp ─▶ RUNNING ─ Ended ─ Stopped
   [ ▶ Start Run ]   [ ■ Stop / Abort ]
──────────────────────────────────────────┬─────────────────────────
  ⚠ PROVISIONAL — uncalibrated raw signal  │  Samples:   12 480
        (live electropherogram,            │  Elapsed:   0:42
         redraws every ~200 ms)            │  Corrupt:   0
     ╱╲      ╱╲___                          │  Errors:    —
 ___╱  ╲____╱      ╲___                     │  ── Log ──
  x = elapsed time · y = raw stream-B       │  → AV.CLEAR ack
                                           │  → START    ack
──────────────────────────────────────────┴─────────────────────────
  Run finished ·  [ Open in analyzer ]   [ Save session (.pck)… ]
```

### Regions

1. **Connection bar.** Transport dropdown (v1: "Simulator" enabled; "Serial
   port…" shown but disabled with tooltip *"rebuild with `--features serial`"*).
   Port dropdown + Refresh (↻) + baud, all inert in v1. **Connect / Disconnect**
   toggle. Status pill: `Disconnected → Connected · read-only → Running → Ended
   → Error`. Connect sends nothing to the instrument.
2. **Lifecycle stepper.** A breadcrumb reflecting `RunState`, so the user always
   sees where the run is. The current state is highlighted.
3. **Run controls.** `▶ Start Run` (enabled only when Connected & Idle);
   `■ Stop / Abort` (enabled while Running). Start launches the guarded
   handshake on the background thread.
4. **Live plot.** Reuses `plot.rs` unchanged — the snapshot is a real
   `Electrophoresis`. A persistent "PROVISIONAL / uncalibrated" banner.
5. **Telemetry sidebar.** `sample_count()`, elapsed time, `corrupt`, `errors`,
   and a scrolling command/ack log (`acks`).
6. **Post-run bar.** On `Ended → Stopped`: **Open in analyzer** (hands the final
   `to_electrophoresis()` to `add_file`, exactly like `simulate_added_file`, so
   the normal calibration pipeline runs) and **Save session…** (write `.pck`).

### Safety approval modal

When the handshake reaches `Action::Setup`, then `Action::Start`:

```
┌ ⚠  Confirm energizing action ─────────────────────┐
│ START begins electrophoresis: it applies HIGH      │
│ VOLTAGE, a LASER and a HEATER to the physical      │
│ chip. This is irreversible.                        │
│                                                    │
│                        [ Abort ]  [ Approve ]      │
└────────────────────────────────────────────────────┘
```

Wording tracks `safety.rs`. **Abort is the default focus.** Refusal sends
`STOP` and ends the run cleanly. (Against the simulator the actions are
harmless, but the modal still fires — the UX is identical to a live run, which
is the point of testing it on the simulator.)

## Threading & state model

- **New module `instrument_view.rs`** (in `traceanalyzer`, headlessly testable,
  same as the other view modules).
- A background thread owns the `Transport` + `Acquisition` and runs a guarded
  poll loop (the `guarded_run` shape, but yielding snapshots). It communicates
  with the UI thread over channels:
  - **snapshot channel** (thread → UI): periodic `AcqSnapshot { state, trace,
    sample_count, corrupt, errors, acks }`, drained by a Slint timer (~5 Hz)
    which re-renders. Cheap: clone counters + the stream-B slice.
  - **approval channel** (thread → UI → thread): the `ApproveWith` closure
    parks, posts an approval request via `slint::invoke_from_event_loop` to show
    the modal, and blocks on a reply channel for Approve/Abort.
  - **stop flag** (`AtomicBool`): Stop/Abort sets it; the poll loop checks it
    each iteration and sends `STOP`.
- `AppState` gains an `Option<InstrumentSession>` (connection state, current
  `RunState`, latest snapshot, log) and a `view_mode` that includes `Instrument`.
- Disconnect / window-close joins the thread and drops the transport.

## Slint surface (new callbacks/properties)

Callbacks: `instrument-connect`, `instrument-disconnect`, `instrument-start`,
`instrument-stop`, `instrument-approve(action, bool)`, `instrument-open-result`,
`instrument-save-session`.
Properties: `instrument-state` (enum→string), `instrument-connected`,
`instrument-transport-list`, `instrument-selected-transport`, `port-list`,
`baud`, `sample-count`, `elapsed`, `corrupt-count`, `error-text`, `log-lines`,
plus the shared live-plot image the other views already use.

## Phases

- **Phase I0 — view scaffold.** Add the `Instrument` view mode + top-bar toggle;
  empty panel with the connection bar and disabled controls. No backend yet.
- **Phase I1 — simulator connect.** Transport dropdown (Simulator only) →
  Connect wraps a `PckReplay`; status pill; Disconnect. Read-only, no run.
- **Phase I2 — background acquisition + live plot.** Start Run spawns the poll
  thread; snapshot channel + Slint timer drive incremental `plot.rs` redraw;
  lifecycle stepper + telemetry sidebar update live; Stop flag.
- **Phase I3 — safety modal.** `ApproveWith` → `invoke_from_event_loop` modal;
  Abort-default; refusal path sends `STOP`. Verify Setup then Start both gate.
- **Phase I4 — post-run.** Open-in-analyzer (→ `add_file`, calibration runs) and
  Save-session (`.pck` writer, if not already present).
- **Phase I5 — serial follow-up (feature-gated).** Enable the Serial transport
  option under `--features serial`: `list_ports` populates the dropdown,
  `serial::open` + `guarded_run_live` replace `PckReplay`. Gated on validating
  outbound framing/baud against a physical 2100 (per `serial.rs` / `safety.rs`).

## Test coverage

- `examples/render_instrument.rs` — render the Instrument panel (disconnected,
  running-with-trace, finished) to PNG headlessly, like the other `render_*`
  examples.
- Headless test (`tests/`) driving I1–I4 against a bundled/synthetic `.pck`:
  connect → start → auto-approve → run to `Ended` → assert `sample_count > 0`
  and that Open-in-analyzer yields a calibratable `Electrophoresis`.
- Approval-refusal test: deny `Start` → assert run aborts and `STOP` sent (mirror
  of the existing `safety.rs` simulator tests, but through the view layer).
- All of the above run with **no display and no `serial` feature**.

## Risks / open questions

- **Snapshot cost.** Cloning the whole stream-B each tick is O(n) in samples;
  for long runs, switch to appending deltas or capping redraw resolution. Start
  simple (full clone at 5 Hz), optimize if a long `.pck` stutters.
- **Elapsed-time axis.** `to_electrophoresis` currently takes a fixed
  `sample_period_s` (sim passes `1.0` → index axis). Real seconds-per-point
  comes from telemetry/assay metadata (Phase P2); until then the live x-axis is
  labeled "sample index", not seconds.
- **Simulator pacing.** `PckReplay` yields all records immediately, so a
  simulated run finishes instantly — fine for tests, but for a *demonstration*
  of the live UX we may want an optional throttle so the trace visibly grows.
- **CLAUDE.md is stale** — it still says "No hardware control"; update it once
  this lands.
