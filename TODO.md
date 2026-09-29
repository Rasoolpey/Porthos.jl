# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md`. Updated at the
end of every session.

Last session: 2026-09-29 (second session that day). P0 to P3 are done. P4 is done for the
synchronous-machine set (GENROU, GENSAL, IEEET1, IEEEG1, IEEEG3, COMPLEXLOAD), P5 and P6 for
the base case, and P7's PHPS gates pass on the base case against parity pack v4
(published). P7 is closed on the PHPS gates (PowerFactory item dropped, LineFault / rk4
left unsupported; see "Decided"). The converters come after the P5 to P7 slice. The test suite passes (1460 tests, about 2 minutes: the P7 gate runs the full
6 s BDF1); the only tests marked broken are the P5 and P6 converter cases.

**Committed and pushed** on `master` (2026-09-29): everything up to and including this
session's P6, P7, precompile workload and pack v4 binding.

## Status by phase

| Phase | State | Gate |
|---|---|---|
| P0 environment + parity pack | done | Pack v4 is bound and published (release `parity-pack-v4`); it loads and every hash checks, locally and from the download. |
| P1 input / output | done | All 60 JSON under `cases/` load, validate and round-trip; all 200 PHPS case files also load (optional test, runs when PHPS_Opt is present). |
| P2 network | done | Y-bus (power-flow and DAE), load admittances, Norton stamps and fault shunts equal PHPS **bit for bit** on all 6 parity cases (gate asks for 1e-12). |
| P3 power flow | done | `V`/`theta` match PHPS within 4e-15 (gate 1e-8), same iteration counts; the case `v0`/`a0` within 1e-6 / 2e-6 relative to the slack (gate as reworded, see "Decided"). |
| P4 components, items 1 to 3 | done | GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD: at 200 samples each, all 43 state-dependent branch sites agree with PHPS and each goes both ways; rhs, both output sets, injection, H and grad H within rtol 1e-12 (rhs bit-identical at 1190 of 1200 samples); parameters identical to PHPS's; contract entry identical. `rhs!` / `outputs!` are type-stable and allocation-free; ForwardDiff Jacobians checked. |
| P4 items 4 and 5 (converters) | after P7 | per the roadmap order of work |
| P5 assembly | done for base | Against PHPS's **compiled C++** `dae_residual` (the function IDA and BDF1 solve): f and g at the equilibrium bit-identical, at 50 random states fault off and on within 1.3e-14 (gate 1e-12; 1228 of 28100 values differ in the last bits, from sin/cos of different math libraries). Layout, wiring, Y-bus, loads, COI and faults identical. Structural Jacobian pattern checked against ForwardDiff. gfl, vsm, droop, voc: `@test_skip` until their models are ported. |
| P6 initialisation | done for base | Residual at Porthos's equilibrium 2.7e-13 (gate 1e-12). Distance to PHPS's equilibrium 9.0e-10, within the allowance 7.2e-7 = `||J^+||` (67) times the residual at PHPS's point (1.07e-8, PHPS's own residual). All initialised parameters (Efd0, Tm0, Vref, Pref, PFD_REF, PM_REF, V0, Vini) within 8.4e-10 of PHPS's. The first pass reproduces PHPS's pre-refinement state within 1e-12. Converter cases pending. |
| P7 simulation | done for base | Pack v4, `test/parity/p7_sim.jl`, all from PHPS's initialised state and parameters. **BDF1** vs PHPS's compiled BDF1: max 1.7e-11 over the whole 6 s (gate 1e-9), the same 268 non-converged steps, the first at the same time. **IDA at 1e-10** vs PHPS IDA at 1e-10: max 1.2e-7 before the first sliding mode (gate 1e-6), 4.1e-3 after it over 15 s (reported). The sliding-mode start comes from the pack: the first run of at least 10 consecutive non-converged BDF1 steps on one equation (IEEEG3_11.xp, 239 steps from 1.324 s; limiter crossings cost 1 to 3). **CSV**: header identical, the 20 complete rows (states and observables) within 2.5e-10 (gate 1e-9, chosen by the agent: the trajectory gate). **IDA production**: 1.7e-2 from PHPS's (reported; step sequences differ). PowerFactory item dropped and LineFault / rk4 unsupported (user decisions). |
| P8 to P12 | not started | |

## Needs the user

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

## Decided (2026-09-29)

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
src/assembly/wiring.jl    InputSource, resolve_wiring (PHPS wire semantics + post-init refresh)
src/assembly/dae.jl       DAESystem, assemble (Y-bus at the power-flow voltages, PHPS's C++
                          constant rounding), dae_residual! (port of the C++ dae_residual)
src/assembly/sparsity.jl  jacobian_pattern (structural, valid in every limiter mode)
src/init/components.jl    init_from_phasor (machines), init_from_targets (exciter, governors),
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
parity/generate/sim.py    sim section: PHPS's compiled BDF1 / IDA runs (5 ms grid, binary)
ext/                      PorthosMakieExt (CairoMakie), PorthosJuMPExt (JuMP): empty stubs
parity/generate/generate_pack.py   pack generator (network, powerflow, records, components)
parity/generate/components.py      components section: PHPS init, instrumented kernels, sampling
parity/generate/dae.py             dae section: PHPS's C++ kernel compiled with a residual harness
scripts/bind_parity_pack.jl        tarball + Artifacts.toml binding (hash taken from the tarball)
scripts/setup.ps1                  one-command install
test/unit/, test/parity/p0..p7     unit tests and the P0 to P7 gates (common.jl:
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

Software tools first (user, 2026-09-29); the control-related items at the end wait for a
method discussion with the user.

1. **Grid-forming converters** (P4 item 5): `GFM_VSM_PHTRUE`, `GFM_DROOP_PHTRUE`,
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
- P10 PH audits (shifted storage, KYP / passivity with JuMP; Clarabel, "Needs the user" 1);
- P11 ROA certificates (interval primitives for the limiters come with it);
- the P9 reservoir, design, OPF and control-mode studies, and the energy-margin /
  reservoir-work study;
- Part II.
