# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md`. Updated at the
end of every session.

Last reviewed: 2026-10-01 (overnight session, the user away; no questions asked, decisions
recorded under "Needs the user"). P0 to P3 are implemented. P4 is done for the
synchronous-machine set and, **new, for the grid-forming converters** (GFM_VSM, GFM_DROOP,
GFM_VOC). P5 passes on base, vsm, droop and voc; P6 passes on base, vsm and voc, and on
droop everything but the residual (1.4e-12 against the 1e-12 gate: a Float64 floor, see
"Needs the user" 1). P7's PHPS gates pass on the base case; the converter cases need parity
pack v6 (generated and bound locally this session, **not published**; see P7 in the table).
Route B (structure-search LMI) now runs to validated, interval-checked results. The P0
Windows test-harness error is fixed. The PowerFactory driver is as before. **Full suite
with pack v6 (2026-10-01): 2916 passes, 4 failures, 2 broken** (the GFL skips). The four
failures are the droop P6 residual and the three converter IDA gates, all explained under
"Needs the user" A and F; nothing else fails.

**Committed** on `master`: everything through P10 plus the governor and terminal passivity
diagnostics. `origin/master` is at `795d6a3` (P10); the three later diagnostic commits are
not pushed yet. **Everything from this session is uncommitted** (the user said not to commit
without asking): Route B (`src/ph/structure.jl`, `ext/PorthosJuMPExt`, `scripts/storage_search.jl`),
the converter models (`src/components/converters/`), converter initialisation
(`src/init/components.jl`, `src/init/equilibrium.jl`), the generator changes
(`parity/generate/components.py`, `generate_pack.py`), `Artifacts.toml` (bound to v6),
tests, `parity/README.md`, this file. **Do not push before pack v6 is published**: CI would
fail to download it.

## Status by phase

| Phase | State | Gate |
|---|---|---|
| P0 environment + parity pack | done | Pack v5 published; pack v6 generated and bound locally (not published). The tamper test no longer uses `chmod`: its copies are written as new files (the `EINVAL` on Windows came from the artifact's permissions). |
| P1 input / output | done | All 60 JSON under `cases/` load, validate and round-trip; all 200 PHPS case files also load (optional test, runs when PHPS_Opt is present). |
| P2 network | done | Y-bus (power-flow and DAE), load admittances, Norton stamps and fault shunts equal PHPS **bit for bit** on all 6 parity cases (gate asks for 1e-12). |
| P3 power flow | done | `V`/`theta` match PHPS within 4e-15 (gate 1e-8), same iteration counts; the case `v0`/`a0` within 1e-6 / 2e-6 relative to the slack (gate as reworded, see "Decided"). |
| P4 components, items 1 to 3 | done | GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD: at 200 samples each, all 43 state-dependent branch sites agree with PHPS and each goes both ways; rhs, both output sets, injection, H and grad H within rtol 1e-12 (rhs bit-identical at 1190 of 1200 samples); parameters identical to PHPS's; contract entry identical. `rhs!` / `outputs!` are type-stable and allocation-free; ForwardDiff Jacobians checked. |
| P4 item 5 (grid-forming converters) | done | GFM_VSM_PHTRUE, GFM_DROOP_PHTRUE, GFM_VOC_PHTRUE against pack v6 (local): branches agree, rhs bit-identical at 199/200, 200/200, 200/200 samples; parameters and contracts identical. Each type also has one **parameter variant** (pf_frame on, VSM field controller, droop adaptive droop and Zv, VOC current-reference limiter) so the dormant branches are checked too. The 60-step virtual-impedance bisection is a root solve (`vi_bisection`): its comparisons are not branch sites, its ForwardDiff derivative is the implicit-function one (value bit for bit). One VOC site (`pvoc_mode = 1`) has no PHPS Python reference (py_codegen lacks `tanh`) and is recorded as unreached. Allocation-free, type-stable, Jacobians checked (unit tests). |
| P4 item 4 (GFL, GFL_ZIF) | waits | for the GFL reservoir model (control discussion) |
| P5 assembly | done for base and the GFM cases | Against PHPS's **compiled C++** `dae_residual`: base, vsm, droop, voc within the 1e-12 gate, all but 1228 / 1104 / 454 / 326 values bit-identical (the rest differ in the last bits: sin/cos of different math libraries). Layout, wiring, Y-bus, loads, COI and faults identical. gfl: `@test_skip` until GFL is ported. |
| P6 initialisation | done for base, vsm, voc; droop residual 1.4e-12 | Base as before (residual 2.7e-13, 9.0e-10 from PHPS). Converters: PHPS keeps each converter's Norton source u_out at its first-pass value and recaptures the set-points (p_set, u_set, q_set, PSET_REF, V_nom) at the final bus voltage (checked on the pack: u_out unchanged to 1.6e-14, no polish). Porthos does the same: in the Newton solve each converter's set-point rows are replaced by "Norton current = first-pass value", the set-points are recaptured after it (`converter_init`), and the lag states are settled on their own equations in Float64 (`lag_states`, `_settle_lags!`). vsm: residual 2.1e-13, voc 4.3e-13, droop **1.4e-12 (fails the 1e-12 gate)**; all three 3.2e-9 from PHPS (allowed about 2e-5), first pass within 1e-12, every initialised parameter within 3.2e-9 of PHPS's. |
| P7 simulation | done for base; converters: BDF1 and CSV pass, IDA gate fails | Pack v4, `test/parity/p7_sim.jl`, all from PHPS's initialised state and parameters. **BDF1** vs PHPS's compiled BDF1: max 1.7e-11 over the whole 6 s (gate 1e-9), the same 268 non-converged steps, the first at the same time. **IDA at 1e-10** vs PHPS IDA at 1e-10: max 1.2e-7 before the first sliding mode (gate 1e-6), 4.1e-3 after it over 15 s (reported). The sliding-mode start comes from the pack: the first run of at least 10 consecutive non-converged BDF1 steps on one equation (IEEEG3_11.xp, 239 steps from 1.324 s; limiter crossings cost 1 to 3). **CSV**: header identical, the 20 complete rows (states and observables) within 2.5e-10 (gate 1e-9, chosen by the agent: the trajectory gate). **IDA production**: 1.7e-2 from PHPS's (reported; step sequences differ). PowerFactory item dropped and LineFault / rk4 unsupported (user decisions). **Converter cases (pack v6, local):** BDF1 vs PHPS's compiled BDF1: vsm 1.2e-11, droop 1.6e-11, voc 6.3e-11 (gate 1e-9), the same non-converged steps (127, 523, 9); CSV header identical, rows within 2.7e-10. **IDA at 1e-10 fails its 1e-6 gate**: 9.6e-3 (vsm), 2.9e-4 (droop), 6.1e-3 (voc) before the detected sliding mode (vsm 2.269 s, droop 1.4165 s, voc none). IDA agrees with PHPS to about 1e-7 through the fault and clearing, then departs at an isolated limiter event: vsm 1.355 s and droop 1.405 s (IEEEG3_11 pilot valve), voc 1.955 s (IEEEG1_10 valve). Diagnosis: Porthos's own IDA at 1e-10 and 1e-11 separates at the same time in voc (4.3e-4 at 1.950 s), and in vsm IDA at 1e-11 fails at 1.3365 s (h = 3e-11); so the IDA solution there depends on tolerance and step placement, as the base case's after its sliding mode; BDF1's 1e-11 agreement shows the models are the same. See "Needs the user" F. |
| P10 PH audits | parity gate done; phase partial (base) | `src/ph/` (storage, reduced field and exact Jacobian, physical projection, shifted-storage audit, port models), `test/parity/p10_ph.jl` against pack v5, `scripts/ph_audit.jl`. Shifted storage at PHPS's point: rank 54/171, 41 positive eigenvalues of sym(SA), max 12.882693018 (PHPS 12.882693020, which differenced; Porthos's Jacobian is exact), min within 3.4e-10, exact counterexamples along the top eigenvector within 1e-9; held states = the ten IEEET1 `xr` (Tr = 0). Governor ports: all nine IEEEG1 at speed -> Tm match PHPS's KYP and exact-arithmetic records to 1e-15 (crossing 1.934717940050879 rad/s); KYP infeasibility certified in the frequency domain. Route B structure search (`src/ph/structure.jl`, `decay_margin` in the JuMP extension, `scripts/storage_search.jl`): decay-margin SDP per pattern with interval-checked dual certificates and Lyapunov checks (results under "Route B results"). Open: the port-residual audit. |
| P8, P9, P11, P12 | not started | |

## Needs the user

New from the overnight session (2026-10-01):

A. **P6 droop residual: 1.42e-12 against the 1e-12 gate.** Everything else in P6 passes on
   droop (3.2e-9 from PHPS, allowed 2.6e-5; parameters within 3.2e-9). The rows left are the
   measured-current lags (I - I_meas)/Tmeas_i: I = (u_out - V)/(j Zseries) carries the
   rounding of u_out divided by Zseries = 0.0143 (1.5e-14 at |I| = 7 pu), then 1/Tmeas_i = 200
   (1.4e-12 is one to two ulps of I). In the droop, u_out depends on I_meas through
   r_vi = Zseries (loop gain 1), so no floating-point I_meas zeroes both rows; Newton plus the
   lag settling reaches 1.42e-12 (vsm 2.1e-13 and voc 4.3e-13 pass). The gate was not
   loosened and the test fails. Options: keep it failing; or state the gate per row relative
   to the row's scale (e.g. |f_i| <= 1e-12 max(1, |a_i|/T_i)); or accept 2e-12. Your call.
B. **Publish parity pack v6?** Generated this session from PHPS ba11ea1 (same docs-only
   dirt), bound locally in `Artifacts.toml`; tarball `parity/dist/parity_pack-v6.tar.gz`
   (hashes under "Decided"). Publishing it (GitHub release `parity-pack-v6`, as v3 to v5) is
   outward-facing and waits for your approval; until then do not push `master`.
C. **P4 gate reading for the converters** (agent's choices, please confirm): the
   comparisons inside the 60-step virtual-impedance bisection are treated as a root solve,
   not branch sites (marked `in_loop` in the pack, not matched; the regime switch after the
   loop is matched); one parameter variant per converter type exercises the branches the
   case parameters keep fixed (stronger than exempting them); the VOC `pvoc_mode = 1` site
   is recorded as unreached because PHPS's Python kernel translation (`py_codegen`) has no
   `tanh`, so PHPS itself cannot run that mode outside C++ (a PHPS finding; no case uses it).
D. **Route B results** (below) need the method discussion: a strictly decaying quadratic
   storage needs machine-to-machine couplings across units; which physical storage terms
   those should become is a control question.
E. **Commit**: nothing from this session is committed.
F. **P7 IDA gate on the converter cases** (fails, not loosened): the gate's "until the first
   limiter enters a sliding mode", detected as 10 or more consecutive non-converged BDF1
   steps on one equation, does not cover the isolated limiter events (IEEEG3 pilot valve,
   IEEEG1 valve) after which IDA at 1e-10 already depends on its step placement (evidence in
   the P7 row). Options: gate until the first limiter event after the fault is cleared (a
   rule to agree on), or add event detection on the switching surfaces (roadmap 2.3) so IDA
   steps onto the kinks, or keep IDA reported-only on these cases with BDF1 as the exact
   check (as after the sliding mode on the base case).

Earlier items:

1. **Clarabel** is left out of the JuMP extension. Clarabel pins TimerOutputs 0.5, and the
   SciML stack (NonlinearSolveBase via Sundials and NonlinearSolve) needs TimerOutputs 1.x,
   so they cannot resolve together. SCS and COSMO are in. Revisit when Clarabel updates, or
   decide whether it matters for P10 and II.1.
2. **PHPS Hamiltonian at 50 Hz** (for information, no parity case affected): PHPS's
   `hamiltonian` / `grad_hamiltonian` for GENROU and GENSAL compute the base frequency as
   `2*pi*params.get('fn', 60)`, but machine params never contain `fn`, so they use 60 Hz
   even in a 50 Hz case. Porthos uses the machine's `omega_b`, which is the same number at
   60 Hz (all parity cases) and the correct one at 50 Hz.
3. The P6 check of the initialised parameters uses an absolute tolerance of 1e-8 chosen by
   the agent (the roadmap gives none); the actual differences are below 8.4e-10.

## Decided (2026-10-01, overnight, by the agent; see "Needs the user" to revisit)

- **Parity pack v6 generated and bound locally, not published**: `parity/pack-v6`, tarball
  `parity/dist/parity_pack-v6.tar.gz`, sha256
  `400974ee39cb16ac45c9628916facbf21583fa225bad0af4d863d7219b050136`, tree hash
  `c99fe023705c11cc10260a9dcbc608421daa3cb8`. Against v5 it adds 21 files
  (`components/GFM_*.json`, `sim/{vsm,droop,voc}/*`); the only changed files are the three
  `sim/base/*.json`, by PHPS's wall time. The component seeds are now stable per type
  (`SEED_ORDER` in `generate/components.py`): a first build used the sorted-list index and
  shifted the seeds of the IEEE types, so it was discarded and the pack regenerated.
  `parity/README.md` has its v6 row, marked not uploaded; update it after publishing.
- **Converter initialisation follows PHPS's invariant** (the Norton source stays at its first
  pass, set-points are recaptured; checked on the pack), and lag states are settled on their
  own equations after the Newton solve (general mechanism: `lag_states` per model).
- **Route B uses the trace-normalised decay margin**, not "Q >= I, rate <= -eps I", and
  patterns are mapped from the physical coordinates (the reference-angle fix).
- **P0**: the tamper test writes fresh copies instead of `chmod` (the `EINVAL`).

## Decided (2026-09-29)

- **P10 audit tooling first** (user: "let's follow your next step"), before any model rework:
  the tools judge every candidate storage. Findings beyond PHPS's records, from
  `scripts/ph_audit.jl` on Porthos's own equilibrium: the IEEEG1 speed -> Tm port has all
  zeros in the left half-plane; its non-passivity comes from relative degree 2 (6 poles, 4
  zeros), so no parameter choice makes that port positive real. IEEEG3 has a zero at
  +1.333 1/s and Re H < 0 on 1.3189 to 14.8932 rad/s (PHPS's letter: +1.33, 1.32 to 14.9).
  **KYP infeasibility** is certified by a frequency with Re H(jw) < 0 (then the positive-real
  LMI has no solution, a theorem), not by a solver status; the JuMP LMI solver comes with
  the storage search. Parity pack v5 adds `records/model_review/audit_controller_kyp.json`
  and `audit_governor_nonpassivity_exact.json` (generator: `MODEL_REVIEW_RECORDS`).
- **Parity pack v5 published** (user approved), as v3 and v4: REST API with the stored git
  credential, download checked against both bound hashes.
- **PowerFactory driver** (user: build the PowerFactory tooling first, in Porthos, without
  copying PHPS_Opt's 102 MB `pf/` folder). `pf/` is a small Python package, standard
  library only, run with PowerFactory 2022 SP1's Python 3.10: `py -3.10 pf/run.py inspect`
  (model and load flow to JSON) and `py -3.10 pf/run.py simulate <scenario>` (the
  scenario's BusFaults as `EvtShc` in a fresh copy "Porthos" of study case "Base"; full
  precision CSV and `run.json` with the column map). Julia: `src/io/powerfactory.jl`
  (`pf_simulate`, `pf_inspect`, `read_pf_results`, `pf_compare` with PHPS's metrics) and
  `scripts/pf_compare_fault.jl` (one command). `AGENTS.md` and the roadmap (2.5, P7) now
  name `pf/` as the second Python exception (it drives PowerFactory, computes nothing) and
  the only way to run PowerFactory. What was learned about PowerFactory 2022 is in
  `pf/README.md` (engine vs window, the INI workaround, time units in ms, full-precision
  export, the column map).
- **Result on the base case**: load flow within 2.2e-9 pu / 1.8e-6 deg; bus-16 fault
  (150 ms), Porthos IDA against PowerFactory: rotor angles 0.38 to 0.65 deg RMS (max
  1.32 deg), speeds below 9.5e-5 pu RMS, P 0.015 to 0.035 pu RMS; every metric smaller than
  PHPS's recorded comparison (0.48 to 0.79 deg RMS, max 1.68 deg).

- **Parity pack v4 is published**: release `parity-pack-v4` with `parity_pack-v4.tar.gz`
  (sha256 `d055ebf55875f72d6f5c1e6632b9be5f98c7a549e68cd2dd7c681733f36399c0`, tree hash
  `2d5a8c246703fe42d90791076935381a385f6b98`, both checked on the downloaded file), made
  like v3 (REST API, stored git credential).
- **P7 is closed on the PHPS gates** (user, 2026-09-29): the PowerFactory item is dropped.
  The user validated PHPS against PowerFactory, so parity with PHPS carries it over;
  PowerFactory 24 is not on this machine, and the raw PF traces are gone from PHPS_Opt.
  LineFault and rk4 stay unsupported for now (their only scenario,
  `IEEE39Bus/mid_line_fault.json`, points at a system file that has never existed; it is not
  the most severe case anyway). `docs/ROADMAP.md` P7 updated.
- **Control-related work waits for a method discussion** (user, 2026-09-29): the three
  energy-reservoir models (governors and GFL) and the control tools (PH audits, KYP /
  passivity, ROA certificates, design) need careful thought about the method first. The
  software tools come first.
- **Parity pack v3 is published**: GitHub release `parity-pack-v3` on Rasoolpey/Porthos.jl
  with `parity_pack-v3.tar.gz` (sha256
  `b2ceaed1c8b44ed95b1180fc33e2a771b5463f1d08b8c4c91729cb32e335398e`, tree hash
  `c7353f05cb079034198de322f9dcfae7f9feaf68`, both checked on the downloaded file), created
  through the GitHub REST API with the stored git credential (no `gh` CLI on this machine).
  CI and fresh clones can fetch it. A future pack (v4, ...) needs its own release the same
  way; the tarball must be uploaded exactly as `bind_parity_pack.jl` wrote it
  (`archive_artifact` output is not byte-reproducible).
- **GFL-ZIF is deferred.** PHPS at ba11ea1 cannot initialise any of the three GFL-ZIF system
  files (coupled network solve does not converge at t = 0), so `GFL_ZIF_PHTRUE` has no
  reference yet. The user will come back to the energy-reservoir model of the GFL converter
  later; until then GFL_ZIF_PHTRUE is not ported. (It is not a P5 to P7 parity case.)
- **P3 gate wording changed** (user approval): bus voltages match PHPS within 1e-8, and the
  case `v0` within 1e-6 and `a0` within 2e-6 relative to the slack bus (the stored values
  are rounded to 6 decimals, and the converter cases keep PowerFactory's angle reference).
  `docs/ROADMAP.md`, `README.md` and `test/parity/p3_powerflow.jl` updated; no more
  `@test_broken` there.
- The repository has one branch, `master`; commits go there directly (commit and push only
  when asked). Pushing `master` runs CI.
- **IDA gate reworded** (user: both integrators' results are accurate within their
  limits): IDA against PHPS IDA at `rtol = atol = 1e-10` must agree within 1e-6 until the
  first limiter enters a sliding mode; afterwards the difference is reported, not gated.
  Reason: where the IEEEG3 pilot valve chatters (bus-16 fault, from 1.3 s) the IDA solution
  depends on step placement, not tolerance (Porthos at 1e-9 vs 1e-10 differs by 4e-2 there;
  at 1e-11 IDA fails). BDF1 stays the exact check. Possible later improvements (not now):
  event detection on switching surfaces (roadmap 2.3), finer steps, or retuned limiter
  parameters. `docs/ROADMAP.md` updated.
- **Scripts must be one command** (user preference): `julia --project=. scripts\<name>.jl`
  from the repository root, no extra environments or environment variables. Optional
  tools such as CairoMakie are installed into Julia's default environment on first use.

## Environment on this machine

- **PowerFactory 2022 SP1** at `C:\Program Files\DIgSILENT\PowerFactory 2022 SP1`, Python
  module for 3.10 (Python 3.10.11 is installed: `py -3.10`). Project
  "39 Bus New England System", study case "Base". **Close the PowerFactory window** before
  any `pf/` run (the engine cannot start while it is open: exit code 4002). Compare a
  scenario: `julia --project=. scripts/pf_compare_fault.jl [scenario.json]` (about 5 s of
  PowerFactory).

One-command setup: `powershell -ExecutionPolicy Bypass -File scripts\setup.ps1 [-RunTests]`
(juliaup, Julia 1.12 as a juliaup override for this directory, `Pkg.instantiate`, the parity
pack, the generator venv from `parity/generate/requirements.txt`, MSYS2 UCRT64 g++ and
SUNDIALS at `C:\msys64`). Safe to re-run.

- Julia 1.12.7 via juliaup (`winget` id 9NJNWW8PVKMN), default channel 1.12. `julia` is on
  PATH (WindowsApps alias).
- `Project.toml` compat `julia = "1.12"`; `Manifest.toml` resolved with 1.12.7; CI runs
  1.12 on ubuntu and windows.
- PHPS runs here: venv at `parity/generate/.venv` (Python 3.13.15, pinned in
  `parity/generate/requirements.txt`; git-ignored). Run PHPS with `PYTHONDONTWRITEBYTECODE=1`
  and `PYTHONIOENCODING=utf-8` from `PHPS_Opt/phps` so its tree stays untouched. PHPS's full
  initialisation runs in pure Python: `DiracRunner(path, output_dir=<temp dir>).build(
  solver="scipy")` stops before C++ generation (always pass `output_dir`, or it writes into
  PHPS_Opt); with `solver="bdf1"` or `"ida"` it also writes and compiles the C++ kernel.
- **g++ 16.2.0 and SUNDIALS 7.5.0** from MSYS2 UCRT64 at `C:\msys64\ucrt64` (installed with
  winget `MSYS2.MSYS2` and pacman), which is where PHPS looks for them on Windows. So PHPS's
  compiled BDF1 and IDA runs work here (needed for the P7 references).
- PHPS_Opt is at `ba11ea1`, with uncommitted docs-only changes (recorded in the pack
  manifest). The generator refuses to run if PHPS inputs are dirty.
- Run a fault and plot it: `julia --project=. scripts/run_base_fault.jl` (BDF1 and IDA of
  the base-case bus-16 fault; results in `outputs/IEEE39Bus_PF_base_bus16_150ms_{bdf1,ida}/`,
  figures in `outputs/figures/`). The user ran it successfully (2026-09-29); the figures
  then got a distinct colour per machine and legends outside the axes. BDF1 takes about
  80 s (the faithful PHPS scheme: a finite-difference Jacobian of 282 residual calls and a
  dense LU per Newton iteration; PHPS's C++ takes 67 s); IDA about 5 s (PHPS: 6.8 s).
- **Start-up time** (the user's complaint about compiling on every run): `src/precompile.jl`
  is a PrecompileTools workload (load, equilibrium, BDF1 and IDA steps with a fault switch,
  the CSV / JLD2 / run.json writers). Precompiling Porthos takes about 20 s after a source
  change; a fresh process then reaches the simulation in about 2 s (before: about 22 s of
  compilation). A PackageCompiler system image remains possible if the 2 s ever matter.
- Tests: `julia --project=. test/runtests.jl` (about 40 s), or `Pkg.test()`, which is
  slower (it precompiles a fresh test environment). One gate alone:
  `julia --project=. -e 'using Porthos, Test; const ROOT = pwd(); include("test/parity/common.jl"); include("test/parity/p4_components.jl")'`.
  To test against a local pack before binding it: `PORTHOS_PARITY_PACK=parity/pack-v4`.
- Regenerate the pack: see `parity/README.md` (new baseline = new version, new reason).

## What exists

```
src/Porthos.jl            module, exports
src/io/expr.jl            safe expression reader (C order, Float64, env of named values)
src/io/json.jl            read/write JSON keeping key order and int/float; json_identical
src/io/schema.jl          JSONSchema validation (cases/schema/*.schema.json)
src/io/case.jl            Case + typed tables with PHPS defaults; ComponentSpec, Wire
src/io/scenario.jl        Scenario, SolverSettings, BusFault/LineFault/OtherEvent (PHPS DAE defaults)
src/io/contracts.jl       ContractSet / ContractEntry / DomainClause
src/io/params.jl          MODEL_TYPES (PHTRUE set + COMPLEXLOAD), machine-base normalisation,
                          component_params + per-type constructor defaults (type_defaults!)
src/io/parity.jl          parity pack: artifact lookup, SHA-256 verification, readers
                          (pack_json; pack_binary for the sim .bin files)
src/network/ybus.jl       Network (sorted bus ids), ybus / ybus_pf / ybus_dae, Norton stamps,
                          load admittances; CPython complex division for bit parity
src/network/events.jl     fault admittance, fault shunts, with_fault, LineFault line split
src/powerflow/newton.jl   exact port of PHPS Newton-Raphson; skip_pf_solve overrides
src/components/primitives.jl   branch primitives with recorders (gt/ge/lt/le, clamp_mode,
                          nonwindup, outband_relax, guard_min, select); NoModes / ModeLog
src/components/interface.jl    AbstractComponent, rhs!/outputs!/step_outputs!/modes,
                          hamiltonian, injection, norton_admittance, contract,
                          build_component / with_params, COMPONENT_CONSTRUCTORS
src/components/machines/  genrou.jl, gensal.jl
src/components/exciters/  ieeet1.jl
src/components/governors/ ieeeg1.jl, ieeeg3.jl
src/components/loads/     complexload.jl
src/components/converters/ common.jl (vi_bisection, the virtual-impedance divider), gfm_vsm.jl,
                          gfm_droop.jl, gfm_voc.jl
src/assembly/wiring.jl    InputSource, resolve_wiring (PHPS wire semantics + post-init refresh)
src/assembly/dae.jl       DAESystem, assemble (Y-bus at the power-flow voltages, PHPS's C++
                          constant rounding), dae_residual! (port of the C++ dae_residual)
src/assembly/sparsity.jl  jacobian_pattern (structural, valid in every limiter mode)
src/init/components.jl    init_from_phasor (machines), init_from_targets (exciter, governors),
                          converter_init / converter_current (converters),
                          first_pass (PHPS Initializer.run: power split, machine links)
src/init/equilibrium.jl   solve_equilibrium: first pass -> Gauss-Newton on the full DAE,
                          slack set-point closing the power balance, Vini fixed point,
                          reservoir references; component_io
src/components/dispatch.jl     concrete-type dispatch over the models (the residual takes
                          5.6 us and does not allocate on IEEE-39)
src/components/observables.jl  PHPS's logged signals per model (delta_deg, Te, H_steam, ...)
src/assembly/dae.jl       also DAEWorkspace + the allocation-free residual with per-fault flags
src/sim/bdf1.jl           simulate_bdf1: port of PHPS's solve_bdf1; SimResult
src/sim/ida.jl            simulate_ida (Sundials.jl), consistent_voltages! (IDACalcIC-like)
src/sim/results.jl        csv_columns / csv_row / write_results_csv / write_results_jld2,
                          run_metadata (records the settings the integrator used), simulate
src/precompile.jl         PrecompileTools workload for the simulation path
scripts/run_base_fault.jl simulate the base-case bus-16 fault and plot it (one command)
src/io/powerfactory.jl    PowerFactory results reader, pf_simulate / pf_inspect (run pf/),
                          pf_machine_map (case _pf keys), pf_compare (PHPS's metrics)
pf/                       PowerFactory driver (Python 3.10, stdlib): run.py, porthospf/
                          (session, inspect_model, simulate), config.json, README.md
scripts/pf_compare_fault.jl  PowerFactory and Porthos on one scenario, compared (one command)
src/ph/storage.jl         storage components, H / grad H / Hess H, solve_network, reduced_field,
                          reduced_jacobian (exact: f_x - f_V g_V^-1 g_x)
src/ph/audit.jl           PhysicalProjection (contract reservoirs, held states, COI section),
                          shifted_storage_audit
src/ph/ports.jl           PortModel (a component linearised at one port), transfer,
                          real_part_crossings, passivity_certificate, port_zeros,
                          frequency_response (Hessenberg), loop_port_model (the rest of
                          the system seen from one component, exact)
scripts/ph_audit.jl       the audits at Porthos's equilibrium, JSON report (one command)
scripts/governor_passivity.jl  governor shortage vs the rest's excess, per governor
src/ph/dissipativity.jl   MultiPortModel, open_loops_model (several loops at once), multiport_margin,
                          port_margin, loop_margin, kyp_riccati (storage from the KYP Riccati
                          equation), port_storage, rest_storage
scripts/joint_governor_storage.jl  joint storage at the governor ports (negative result)
src/ph/terminal.jl        route A: sync_jacobian, TerminalModel / NetworkModel,
                          terminal_models (the cut at the machine terminals), terminal_margins
scripts/terminal_passivity.jl  route A test, per unit and per layer (negative result)
src/ph/structure.jl       route B: state_groups, storage_pattern, section_pattern,
                          reference_section, pow2_scaling, lyapunov_check, margin_residuals,
                          verified_min_eig (interval), pattern_certificate, verified_lyapunov,
                          rank_couplings, add_coupling!, couple_states!; decay_margin and
                          structured_lyapunov in ext/PorthosJuMPExt
scripts/storage_search.jl route B search (margins, hypotheses, dual-guided greedy, L1;
                          Hypatia, capped at 300 s per solve)
parity/generate/sim.py    sim section: PHPS's compiled BDF1 / IDA runs (5 ms grid, binary)
ext/                      PorthosMakieExt (CairoMakie); PorthosJuMPExt (JuMP), with the
                          structure-search SDP currently in progress
parity/generate/generate_pack.py   pack generator (network, powerflow, records, components)
parity/generate/components.py      components section: PHPS init, instrumented kernels, sampling
parity/generate/dae.py             dae section: PHPS's C++ kernel compiled with a residual harness
scripts/bind_parity_pack.jl        tarball + Artifacts.toml binding (hash taken from the tarball)
scripts/setup.ps1                  one-command install
test/unit/, test/parity/p0..p7, p10   unit tests and the P0 to P7 and P10 gates (common.jl:
                                   phps_init_params, phps_initial_state)
```

## Findings about PHPS that later phases depend on

- **Only the PHTRUE models** (plus COMPLEXLOAD) are ported; the retired models stay in
  PHPS_Opt as reference. Five PHTRUE classes inherit their equations from retired parents
  (GENROU, GENSAL, IEEET1, IEEEG1, IEEEG3); each Porthos model is the flattened result.
- **DAE, not ODE.** Porthos ports PHPS's DAE path (`src/dirac/`: `DiracCompiler`,
  `DiracRunner`; 68 uses in tools and studies). The Kron-reduced ODE path
  (`SimulationRunner`, `compiler.get_z_bus_kron`) is legacy (`tools/run_simulation.py`
  only) and is not ported.
  - `LineFault` exists only in the ODE runner; `DiracRunner._inject_events` ignores it.
    Porthos has the topology (`split_line_for_fault`), but PHPS has no DAE reference for the
    one LineFault scenario. Decide at P7.
  - The one `rk4` scenario cannot run on the DAE path in PHPS either. Decide at P7.
- **The Y-bus PHPS simulates with is not the one from the case `v0`.** Before building the
  DAE, `DiracRunner.build` overwrites every bus `v0`/`a0` with the solved power-flow
  voltages, and every load's `V0` and `Vini` too. So the simulation Y-bus has the PQ loads
  at `(P - jQ)/V_pf^2`. The P2 gate compares the case-`v0` matrix (`DiracCompiler.build`
  alone), which Porthos matches bit for bit. **P5 must build `ybus_dae` with the power-flow
  voltages** (add a `v0` override to `ybus_dae` / `load_admittances`) and the pack's `dae`
  section should record the matrix `runner.dae_compiler.Y_full` after `runner.build`.
- **Initialisation writes parameters.** PHPS's initialisation adds or changes: GENROU /
  GENSAL `Efd0`, `Tm0`; IEEET1 `PFD_REF`, `Vref`; IEEEG1 / IEEEG3 `PM_REF`, `Pref`;
  COMPLEXLOAD `V0`, `Vini` (power-flow voltage, later `Vini` = DAE-consistent voltage).
  The pack records them per instance (`init_set`); the P4 gate takes them from the pack.
  P6 must reproduce them. The kernels bake some into literals (`PFD_REF`, `PM_REF`), so
  they must be set before a simulation starts.
- **PHPS's DAE is its compiled C++**, not `py_solver` (which differs slightly: unrounded
  COI weights, a network solve instead of KCL residuals). `DiracCompiler._emit_dae_residual`
  writes the network constants with `%.10e` (Y-bus, load arrays, slack references, fault
  G/B and times: 11 significant digits) and the COI weights with `%.6f`; Porthos rounds the
  same way (`assemble(...; phps_rounding = true)`, the default). The residual is
  `res = [ydot - f; g]` with KCL `g_i = I_inj,i - (Y V)_i - Y_fault,i V_i - I_freq,i`, dense
  sums over j in order; a slack bus without a machine gets `V - V_ref` (none in the parity
  cases). COI: members are generators with an `omega` state (weight 2H, or Ta), every rotor
  angle gets `- omega_b (omega_COI - 1)`, `delta_COI' = omega_b (omega_COI - 1)`, `omega_b`
  from the first member.
- **How Porthos initialises (P6)**, replacing PHPS's multi-pass refinement: the first pass
  (a port of `Initializer.run`: machines from `I = conj(S/V)` at their power-flow phasors,
  exciters and governors from the machine targets, loads with `V0 = Vini = |V_pf|`), then a
  Gauss-Newton solve (QR least squares, rows equilibrated) of the full DAE. Unknowns: all
  states except `delta_COI`, the reservoir levels and the slack machine's rotor angle
  (rotation gauge), plus all bus voltages. Equations: all rows except the reservoir rows
  and the `delta_COI` row (the frequency is left free). Around it, the slack machine's
  set-point (its governor's `Pref`, or its `Tm0`) is adjusted by a scalar Newton step until
  `omega_COI = 1`, the way the power-flow slack absorbs the losses, and `Vini` is set to the
  solved `|V|` until it no longer changes. Finally each reservoir reference is set to the
  power it supplies. On the base case PHPS's whole refinement moves the state by at most
  6.3e-12 from its first pass, which is why this reproduces PHPS.
- **PHPS's equilibrium is not tight**: at the initialised `(x*, V*)` its own residual is
  up to 3e-10 in f (8e-9 in the converter cases) and 1e-8 in g. The P6 gate allows for this
  (`||J^-1|| r_PHPS`).
- **Stale outputs in the output pass**: in PHPS's residual, an output kernel that reads
  inputs sees the outputs of components later in the order from the previous call. None of
  the synchronous-machine set's output kernels read component outputs (only bus voltages),
  so this does not matter yet; check it for the converters.
- **Order of evaluation in PHPS's DAE:** all `out` kernels first, then all `step` kernels;
  a step kernel overwrites some outputs (machines: Pe, Qe, id, iq, It, i_fd). Porthos:
  `outputs!` = the out kernel, `step_outputs!` = after the step.
- **COMPLEXLOAD constants are rounded** to 13 significant digits (`%.12e`, exponents
  `%.6f`) when PHPS writes them into C++; Porthos rounds the same way (`_c12e`, `_c6f`).
- DAE bus fault: `Y_f = (r - jx)/(r^2 + x^2)` with `z^2 >= 1e-20`; a missing `x` means
  `1e-5` (bolted, PowerFactory-like), a missing `r` means 0.
- DAE Y-bus = lines + shunts + PQ loads as `(P - jQ)/v0^2` + generator Norton
  `1/(ra + j xd'')` for every generator except GFL and GFL_ZIF. VSM, droop and VOC get
  `ra = 0`, `xd'' = Zseries` (default 0.10) from their constructors. `load_G`/`load_B` for
  the residual are overridden by COMPLEXLOAD `P0`, `Q0`, `V0`.
- Machine params (GENROU/GENSAL) are converted from the Sn base to the system base at load
  time, with an "already normalised" heuristic, a `D = 2 Sn/Sbase` default when `D` is
  absent, and an `xd'' <= xl` repair (`src/io/params.jl`).
- PHPS `YBusBuilder` defaults a missing line `x` to 0.001 (`system_graph` uses 0.01); Porthos
  follows the Y-bus.
- Power flow: no Q-limits in PHPS. PHPS's power flow on IEEE-39 converges in 4 iterations to
  2e-11.
- The DAE state count on the base case is n_diff = 203 (including delta_COI) and
  n_alg = 78 (Vd, Vq per bus). GENROU_1 (bus 39, Sn = 10000) has no exciter or governor.
- How PHPS's Python solver assembles the right-hand side (`src/dirac/py_solver.py`
  `PyDAESolver.rhs`): network solve for V (Norton currents + voltage-dependent load
  corrections, slack buses pinned), machine dq frames, output pass, step pass, then the COI
  correction `dxdt[delta] -= omega_b (omega_coi - 1)` with `dxdt[delta_COI] =
  omega_b (omega_coi - 1)` (weights from `coi_weight`). Wiring expressions come from
  `dae_compiler.wiring_map` (e.g. `CONST:0.982000`, `BUS_31.Vterm`, `GENROU_2.i_fd`).

## Deferred from P4 (by design)

- `initialize(c, targets)`: belongs with P6, where PHPS's initialisation chain
  (`init_from_phasor`, `init_from_targets`, the rebalancing passes, the Pref sync and the
  equilibrium polish) is ported and checked end to end.
- `observables(c)`: belongs with P7 (CSV columns such as `delta_deg`, `Te`, `H_steam`).
- Interval methods for the branch primitives (`_gt`, ... returning the decided branch or
  throwing `UndecidedBranch`, with margins): P11.
- The Norton rule still lives in `network/ybus.jl` (`norton_stamps`); the components'
  `norton_admittance` is tested equal to it. Switch `ybus_dae` to the components with the
  converters.
- The residual is allocation-free now (`DAEWorkspace`) but still uses the dense Y-bus
  (39 buses: fine). Larger networks need a sparse KCL with the same summation order.
- id_dq / iq_dq refresh: PHPS overwrites these outputs between the output and step passes
  with 6-decimal coefficients. Nothing in the parity cases reads them; `assemble` refuses a
  case that does.

## Next steps

**Start here next conversation (2026-10-01):** go through "Needs the user" A to F with the
user (P6 droop residual, publish pack v6, P4 reading for the converters, Route B method,
commit, P7 IDA gate on the converter cases). Then, software tools: P8 reports (item 2
below), or event detection on the switching surfaces (roadmap 2.3), which would also answer
F. Route B: the L1 stage did not converge in its cap (see "Route B results", step 6).

The user's direction (2026-09-29): the models must be properly defined before the
stability tooling. The governor, exciter and GFL reservoir models are not passive and need
rework, so that `H` can serve as the Lyapunov function and the ROA builds on it (sources:
`PHPS_Opt/study/presentation/response_to_reviewer_component_models.tex`, the supervisor's
review; `PHPS_Opt/phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` 0.2b and 0.2c, Target B).
The PowerFactory driver was built first, so that reworked models can be checked against
PowerFactory. The agent's proposal, **waiting for the user's answers**:

- Diagnosis: (1) the reservoirs are one-way coupled, give no storage along plant
  directions (H_s rank 54/171; H_gfl rank 1 of 13) and are constant-power supplies;
  (2) the controller layers are not passive at their ports (IEEEG1 proved for any storage,
  IEEEG3 right-half-plane zero, IEEET1 open); (3) the GFL states have no storage.
- Principle for the rework: each reservoir a physical store coupled both ways, fed by a
  constant-effort source (pressure, head, dc voltage behind a loss) instead of a
  constant-power one, with the controller only modulating the power path (valve,
  modulation index). IEEEG1: boiler / steam chest pressure. IEEEG3: penstock water column
  plus head (the RHP zero becomes physical storage). IEEET1: dc supply behind the exact
  E_fd i_fd port, plus AVR storage (maybe joint with the machine). GFL: dc link with dc
  voltage control. Controller states still need storage (roadmap B2, B3).
- Keep the PHPS models as the parity-anchored set; add the reworked ones as new component
  types; validate each against its original (a limit in which it reduces to it, and the
  bounded deviation on the fault cases; now also against PowerFactory with `pf/`).
- **User's answer (2026-09-29): keep the dynamics the same if at all possible**, that is
  the storage-only route (PHPS roadmap B1 to B3 and B4 items 1 to 3): new storage terms,
  joint storages (machine-governor, machine-exciter), a different supply rate or port
  pairing, and different state realisations of the same validated equations. Trajectories,
  PHPS parity and the PowerFactory match stay exact. Known limit: at the *current* ports,
  IEEEG1 is proved non-passive for any storage and IEEEG3 has a right-half-plane zero, and
  no realisation changes an input-output property, so these need a joint storage with the
  machine or another port pairing. The physical two-way reservoirs above would change the
  dynamics and are the fallback only.
- Remaining questions: (1) answered above;
  (2) what counts as validation for a reworked model? (3) do the physical choices match
  what the supervisor wants? (4) start the P10 audit tooling now?
- **Governor test done (2026-09-29, user: "let's do the first test")**:
  `scripts/governor_passivity.jl`. For each governor it opens the loop exactly at the
  governor (`loop_port_model`; the rest of the system as input Tm -> output omega on the
  common-angle section; closing it again reproduces the 171 eigenvalues to 1e-14) and
  compares the governor's shortage nu with the rest's excess rho = min Re(1/G_rest(jw))
  (grid and the limit at infinity, -c'Ab/(c'b)^2). Results on the base case:
  - IEEEG1_2 to IEEEG1_9 (8 sets): rho - nu = +1.38 to +3.79, **passes**. A joint
    governor-plus-rest storage exists with the classical split. Next: construct it and check
    it with `shifted_storage_audit`.
  - IEEEG1_10 (G 09, bus 38): **fails, and not because of the governor.** The rest of the
    system seen from G 09's speed port has a lightly damped zero pair at -0.0367 +- 1.4213i;
    near 1.55 rad/s Re(1/G_rest) = -1470, so the rest is not passive at that port (even the
    frequency-wise margin is -1467). Next: find which part of the rest causes it (exciters,
    other governors, loads), e.g. by opening further loops.
  - IEEEG3_11 (hydro): the index test fails (nu = 50.7 against rho = 2.1), but the
    frequency-wise margin is positive at every frequency (+2.105, at infinity). A joint
    storage exists, but its split needs a frequency-dependent multiplier.
- **Joint governor storage attempted (2026-09-29): it does not exist at the Tm ports.**
  `scripts/joint_governor_storage.jl`, `src/ph/dissipativity.jl` (`open_loops_model`,
  `multiport_margin`, `port_margin`, `loop_margin`, and exact Riccati storages
  `kyp_riccati` / `port_storage` / `rest_storage`, no SDP). The eight governors that pass
  one at a time fail together:
  - constant split (rest pays n_k w^2): fails at 0.048 rad/s, the slow system-wide
    frequency mode, for every eta scanned; without these governors the rest is too weakly
    damped there.
  - even the most generous frequency-wise diagonal split (each governor pays its own
    Re H(jw)) fails: -0.0243 at 6.59 rad/s (-0.027 with all ten). So **no storage split port
    by port at the governors' Tm ports exists, static or dynamic**. The individual passes
    relied on the other governors being inside the rest.
  - the rotor-plus-governor cut cannot help: every machine has D = 0 (as in PowerFactory,
    dpu = 0), so the swing damping lives in the damper windings.
  - at low frequency the lossy network makes the differential speed directions
    non-passive, as Herm(K/jw) with K non-symmetric (the PES paper's path dependence of
    U_net). The eta*Tm^2 supply term is the tool for that.
  - **Next: route A first; if it does not work, route B** (user, 2026-09-29). See "Plan:
    storage that keeps the dynamics" below.
- P10 tooling: **done** (see the status table). Next, in this order: (a) the port-residual
  audit (reconstruct each port's power independently and compare with grad H' f, grouped
  by component and connection; PHPS work package 1 item 3); (b) the storage search that
  keeps the dynamics (PHPS roadmap B1 to B3, B4 items 1 to 3): the network potential U_net,
  then the structure-search LMI for controller and cross-term blocks (JuMP extension,
  block-sparse), then joint machine-governor storage (passivity indices: IEEEG1 shortage
  1.01 to 3.37 against the machine's excess at the same port). Every candidate is judged
  by `shifted_storage_audit`. Ask the user about the method before (b).

### Plan: storage that keeps the dynamics (A first, then B)

Both routes look for a storage V, positive definite on the physical coordinates, with
dV/dt <= 0 along the **unchanged** dynamics (PHPS parity and the PowerFactory match stay
exact). Every candidate is judged by `shifted_storage_audit` and the P10 checks (P > 0 and
sym(PA) < 0 on the section, then the exact nonlinear dH/dt).

**Route A: cut at the machine terminals (physical units). Try first.**
- Units: each machine with its own controllers (governor, and the exciter). The storage is
  the machine's physical energy (rotor, fluxes, damper windings: the machines have D = 0,
  so all swing damping is here) plus governor and exciter terms plus cross-terms *inside
  the unit only*.
- Port: the terminal power V*I (d and q), real electrical power, which cancels exactly at the
  network. Each unit must satisfy the shifted dissipation inequality
  dV_unit/dt <= (power delivered at its terminal), around the operating point.
- Network: Y = N + D. The lossless part N gives the network potential U_net (restoring
  energy for the rotor angles); D gives dissipation. The loads (constant-power parts) and the
  line losses are handled on the network side (incremental supply for each load law).
- Steps: (1) port model of a unit at its terminal (2-port, d and q, in the rotating frame;
  the tools `port_model` / `open_loops_model` need a multi-input, multi-output version);
  (2) frequency test per unit (is it dissipative at its terminal, and with what shortage or
  excess?); (3) network side: N, D, loads; (4) the storages (Riccati, as for the governors)
  and the joint check.
- Expected risk: a unit may still fail at its terminal (watch GENROU_10 / IEEEG1_10, whose
  rest was non-passive near 1.55 rad/s), and the constant-power loads and losses must be
  paid on the network side.
- Result if it works: a Lyapunov function that is a sum of per-machine energies plus the
  network potential; every term has a meaning (the PH picture of review comment 3).

**Route A result (2026-09-29): fails with the incremental terminal power dV'dI.**
`src/ph/terminal.jl` (`sync_jacobian`: the COI correction removed, synchronous frame;
`terminal_models`: 11 units and the network with 19 loads from the exact Jacobian, every
coupling outside the terminals checked absent; `terminal_margins`), `scripts/terminal_passivity.jl`.
The cut is exact: units and network closed again give the 171 section eigenvalues to 1e-12
plus the common rotation (0). But:
- every unit is non-passive at its terminal, and so is the bare machine (Tm, Efd fixed):
  at DC it is a constant-power source (the rotor angle adjusts until Pe = Tm), with an
  indefinite Herm(Y(0)) (eigenvalues +-5.5 to +-26), and near its swing frequency
  (6.9 to 11.3 rad/s, D = 0) it is non-passive too (-13.5 to -217);
- the exciters make it much worse at DC (a high-gain AVR is a regulated voltage source):
  -43 to -3422; the governors add almost nothing;
- the network with its loads is non-passive at DC too (-18.5);
- the whole cut: -3420 at DC. So synchronous machines are not incrementally passive in
  (V, I) terms, whatever the controllers: a known fact, which is why the power-system
  passivity literature uses other port pairs (frequency / active power, voltage magnitude /
  reactive power).
- Possible variant A' (not tried): the same physical units with polar ports (P, omega) and
  (Q, |V|); the lossy network is then the obstacle at low frequency (skew part of the
  synchronising matrix, seen in the governor test).

**Route B: structure-search LMI (coupled storage). If A fails, or to diagnose it.**
- One quadratic storage over all 171 section coordinates, V = z'Pz, P > 0,
  A'P + PA <= -eps I. The trusted physical energy blocks (kinetic, magnetic) are fixed; the
  controller blocks and cross-terms are free within an allowed pattern
  (machine-governor, machine-exciter, network neighbours).
- Objective: the sparsest set of cross-terms (weighted L1). The output lists which couplings
  are needed and how strong, the direct answer to "what storage is missing". Full coupling
  always has a solution (the local V_P that PHPS certified), so the information is the
  sparsity.
- Needs a block-sparse SDP (JuMP extension; PHPS's dense 342-dimensional SCS solve did not
  finish in 25 min). Local and quadratic; the cross-terms are mathematical unless each gets
  a physical derivation, then turned into nonlinear storage terms (roadmap B3).
- Use with A: run B restricted to A's pattern (cross-terms only inside each machine unit,
  plus U_net). Feasible: A works, and B gives the blocks to build. Infeasible: B shows which
  extra coupling is missing.

**Route B results so far (2026-09-30), `scripts/storage_search.jl`, `src/ph/structure.jl`,
`ext/PorthosJuMPExt` (`structured_lyapunov`).**
- Set-up: V = z'Qz on the common-angle section in reference-angle coordinates (plain
  physical states; the reference is GENROU_1.delta, the largest COI weight; the other
  coordinate is recovered from l'x = 0), with Q >= I and As'Q + QAs <= -eps I (eps = 1e-3).
  Patterns (`storage_pattern`): `:unit` (each machine with its governor and exciter as one
  dense block, each load its own block, all rotor angles coupled; 1236 free entries),
  `:component` (each component alone, angles coupled; 541 entries), `:full`.
- **Numerics, learned the hard way:** the unscaled problem is badly conditioned (the dense
  Lyapunov solution has condition 1.36e8; its diagonal runs from 0.01 on the exciter states
  to 4.4e5 on GENROU_1.omega; max |As| = 4250 from fast exciter and damper time constants
  against the slowest mode at -0.027). SCS (first order) with tight tolerances did not
  finish in about 30 minutes (and printed nothing: output to a file needs `flush`); this
  is probably also why PHPS's SCS run never finished. Fix: a diagonal rescaling from the
  dense solution (condition 3.7e3, max |As| 86; it keeps every pattern) and the
  interior-point solver Hypatia (pure Julia, installed into the default environment; no
  conflict like Clarabel's), each stage capped at 300 s. Whole run: about 7 minutes.
- Results:
  - 0. dense (Lyapunov equation): Lyapunov, as expected (PHPS's V_P).
  - 1. `:unit` pattern: solver ALMOST_INFEASIBLE (88 s); the returned Q is about 0,
    so not Lyapunov.
  - 2. `:component` pattern: INFEASIBLE (21 s).
  - 3. `:unit` pattern with L1 on the cross-component entries: ALMOST_INFEASIBLE (304 s,
    the time cap).
- Reading (softened after review, 2026-09-30): **no block-local quadratic storage was
  found.** Hypatia reported the `:component` pattern infeasible (a solver certificate, not
  yet validated independently) and the larger `:unit` pattern ALMOST_INFEASIBLE, which is an
  unconfirmed numerical result, not infeasibility. Even a certified infeasibility of the unit
  pattern would prove only that *some* off-pattern entry is needed in a strictly decaying
  quadratic storage; it would not identify which. Couplings between units (machine speeds,
  fluxes of electrically neighbouring machines) are **hypotheses to test**, not findings.
  Route A is consistent with this, but the two failures are not logically equivalent (A: a
  passivity split with the V*I supply; B: a pattern-restricted quadratic Lyapunov function).
  Part of "almost infeasible" may be the normalisation: Q >= I with a fixed eps is
  homogeneous, and solvers struggle near the boundary.

**Route B plan (agreed 2026-09-30, after review; start the next conversation here).**
1. **Soften and document** (done in this file): the wording above; the steps below.
2. **Solver records and output:**
   - `structured_lyapunov` (ext/PorthosJuMPExt) returns and every report stores the JuMP
     and Hypatia versions (Hypatia 0.11.0 is installed in the default environment, outside
     the project Manifest, so its version must be recorded), primal and dual status, raw
     status, iterations, primal and dual residuals, solve time, and the dual certificate
     (the dual matrices).
   - `scripts/storage_search.jl`: print and flush each stage label *before* the solve; state
     in the report that the L1 penalty acts in the scaled coordinates (T = diag(Q0)^(-1/2)),
     so sparsity rankings depend on that scaling.
3. **Validate the `:component` infeasibility independently:** check Hypatia's dual
   certificate in Porthos (the dual matrix positive semidefinite, and the certificate
   conditions for the pattern), first with Float64 eigenvalues, later with interval
   arithmetic (the roadmap's certificate discipline).
4. **Decay-margin SDP instead of binary feasibility:** maximise t subject to
   I <= Q <= kappa I and As'Q + QAs <= -t I on the pattern (the upper bound keeps it bounded
   and removes the arbitrary eps). It returns the best margin per pattern (how far a pattern
   is from working) and its dual.
5. **Rank the missing inter-unit blocks with the dual:** the dual matrix (or the worst-decay
   eigenvector) at the forbidden positions says which added couplings would raise the margin
   most. Add the highest-ranked blocks, and test the hypotheses along the way (machine speeds
   coupled across units; flux states of electrical neighbours).
6. **Only then the L1 search** (fewest cross-terms), on the first larger pattern that is
   feasible.
Keep every stage bounded (time limit, progress printed and flushed), and say how long it
should take before starting it.

**Route B results (2026-10-01, plan steps 2 to 5 done; step 6 pending).** `julia --project=.
scripts/storage_search.jl` (about 30 min; report `outputs/ph_audit/system_phtrue/storage_search.json`
with the JuMP 1.31.2 / Hypatia 0.11.0 versions, statuses, iterations, times and residuals).
- **Method** (`src/ph/structure.jl`): instead of "Q >= I, rate <= -eps I" (homogeneous, which
  made the solver report ALMOST_INFEASIBLE), the decay margin
  gamma(S) = min lambda_max(As'P + PAs) over P on the pattern S, P >= 0, tr P = 1 (compact;
  S carries a strict quadratic Lyapunov function iff gamma(S) < 0). Its dual is the
  certificate: for any Z = LL' and Y equal to As Z + Z As' on S, gamma(S) >= lambda_min(Y)/tr Z;
  `pattern_certificate` checks it in interval arithmetic (`verified_min_eig`: congruence by
  approximate eigenvectors, then Gershgorin); `verified_lyapunov` checks a primal P the same
  way. Rigorous for the Float64 Jacobian; the scaling is by powers of two (`pow2_scaling`),
  so the scaled matrix is exact. Tested on an analytic 2x2 case (a hand-made certificate).
- **Pattern fix:** the reference angle is -l'z/l_ref, so terms with it become couplings with
  every rotor angle; the earlier patterns, written on the section coordinates, missed these
  for GENROU_1's own states. Patterns are now written on the physical coordinates and mapped
  (`section_pattern`).
- **Numbers** (scaled coordinates; the dense Lyapunov solution gives -4.07e-3, and
  gamma(full) >= 2 max Re eig = -5.4e-2):
  - `:component` and `:unit`: gamma = 0 to solver precision (in [-2e-10, 6.7e-9] and
    [-6.9e-11, 2.4e-9], rigorous lower bounds). No certified infeasibility (that needs an
    exact Z with P_S(As Z + Z As') = 0), but the best block-local margin is at most 2.4e-9.
    So a block-local quadratic storage cannot decay strictly by any useful margin.
  - unit + speed-speed, + speed-angle, + machine-angle across units: still 0.
  - **unit + all machine states coupled across machines (3053 entries): -2.03e-6, verified
    Lyapunov** (rigorous).
  - greedy from `:unit`, adding the 4 inter-component blocks the dual ranks highest per
    step: 0 for 6 steps, then verified Lyapunov from step 7 (2233 entries, -8.4e-7) to step
    10 (2686, -1.8e-5). The blocks it picks: machine-machine pairs (GENROU_7-GENROU_8,
    GENROU_4-GENROU_7, GENROU_8-GENROU_10, GENROU_9-GENROU_11; whether these are electrical
    neighbours is not checked), the reference machine GENROU_1 (bus 39, no
    controllers, D = 0, the largest inertia) with most machines, exciter-machine across
    units, and last the hydro governor IEEEG3_11 with machines.
  - reading: strict decay needs cross-unit coupling of machine states, above all of
    GENROU_1; the margins are 2 to 3 orders below the dense storage's. Hypotheses to take
    to the method discussion, not findings about physics yet.
- **Step 6 (L1): no result yet.** The rerun (2026-10-01, same numbers as above, deterministic)
  ran L1 on the first verified pattern, unit + machine-machine (3053 entries), keeping half its
  margin (rate 1.0e-6): Hypatia stopped at the 300 s cap (TIME_LIMIT) and the iterate is not a
  verified Lyapunov function, so the block list it printed (40 of 55 inter-unit machine blocks,
  led by GENROU_7-GENROU_8, GENROU_9-GENROU_11, GENROU_1 with most machines) is not a finding.
  Next: run L1 on the smaller greedy-7 pattern (2233 entries), or restrict the L1 weights to
  the machine-machine blocks, with a longer cap agreed with the user (per the bounded-run
  rule: probe first, say how long).

Software tools, as before:

1. **Grid-forming converters** (P4 item 5): **ported (2026-10-01)**; P7 on their cases
   with pack v6 (see the P7 row). The stale-output order is not an issue for them (their
   output kernels read only their own states). Still open: switch `ybus_dae` to the
   components' `norton_admittance` (it uses the case parameters, which agree). Earlier text: `GFM_VSM_PHTRUE`, `GFM_DROOP_PHTRUE`,
   `GFM_VOC_PHTRUE`, ported as PHPS has them at `ba11ea1`. Then P5 to P7 on the vsm, droop
   and voc cases: pack v5 with their `components` samples and `sim` runs. **GFL waits**
   (user, 2026-09-29): the GFL converter needs its own energy-reservoir model, like the
   governors and exciters have, and those reservoir models have problems to fix; that is
   control work for the method discussion (below). `GFL_PHTRUE` and `GFL_ZIF_PHTRUE` are
   ported after it. With them: switch
   `ybus_dae` to the components' `norton_admittance`, and check the stale-output order of
   the output pass (see "Findings").
2. **P8 reports** (PorthosMakieExt): the scenario `plots` blocks, PHPS's figure style
   (`tools/figstyle.py`), a summary per run; gate: figures regenerate from records alone.
   Also PHPS tools the roadmap does not list yet: `plot_results.py`, `plot_system_graph.py`,
   `print_system_model.py`, `export_system_pdf.py` (add them to P8 if wanted).
3. **P9 simulation studies**, each gated on a recorded anchor: CCT by bisection and fault
   series (`tools/sweep_critical_clearing.py`, `study/src/cct_*`), frequency indices
   (RoCoF, nadir: `tools/frequency_indices.py`, `mode_rocof_*`), voltage (static Q-V, snap,
   LVRT recovery: `tools/voltage_*.py`, `qv_*`). Then the energy-function studies
   (transient energy CCT, analytic and UEP walls, SOS wall). Needs threaded scenario sweeps
   (roadmap 2.6).
4. **Performance and infrastructure**: `bench/` timing suite; threads for sweeps; a fast
   BDF1 option (sparse AD Jacobian, sparse LU) next to the PHPS-exact one; sparse KCL for
   larger networks; optionally FBDF / Rodas5P with a mass matrix; later events (line trip,
   load step).
5. **P12**: `python/porthos/` wrappers (juliacall) with the lint rule; Documenter pages;
   the fresh-clone gate (one command reproduces the parity suite).

Control-related, waiting for the method discussion (do not start alone):

- the energy-reservoir models: a GFL reservoir built like the governors' and exciters'
  (then GFL_PHTRUE and GFL_ZIF_PHTRUE are ported), and fixes to the governor and exciter
  reservoirs' problems;
- the remaining P10 port-residual audit and validation of the structure-search SDP
  (Clarabel remains an open solver choice; see "Needs the user" 1);
- P11 ROA certificates (interval primitives for the limiters come with it);
- the P9 reservoir, design, OPF and control-mode studies, and the energy-margin /
  reservoir-work study;
- Part II.
