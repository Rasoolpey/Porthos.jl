# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md`. Updated at the
end of every session.

Last session: 2026-09-29. P0 to P3 are done. P4 is done for the synchronous-machine set
(GENROU, GENSAL, IEEET1, IEEEG1, IEEEG3, COMPLEXLOAD) and P5 for the base case; the
converters come after the P5 to P7 slice. The test suite passes (`1313 pass, 10 broken`:
6 are the open P3 question below, 4 are the P5 converter cases, pending their models).

## Status by phase

| Phase | State | Gate |
|---|---|---|
| P0 environment + parity pack | done locally | Pack v2 loads and every hash checks locally. **CI cannot fetch the pack until the tarball is uploaded** (see "Needs the user"). |
| P1 input / output | done | All 60 JSON under `cases/` load, validate and round-trip; all 200 PHPS case files also load (optional test, runs when PHPS_Opt is present). |
| P2 network | done | Y-bus (power-flow and DAE), load admittances, Norton stamps and fault shunts equal PHPS **bit for bit** on all 6 parity cases (gate asks for 1e-12). |
| P3 power flow | done, one clause open | `V`/`theta` match PHPS within 4e-15 (gate 1e-8), same iteration counts. The "case v0/a0 within 1e-8" clause cannot hold (below); marked `@test_broken`, not loosened. |
| P4 components, items 1 to 3 | done | GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD: at 200 samples each, all 43 state-dependent branch sites agree with PHPS and each goes both ways; rhs, both output sets, injection, H and grad H within rtol 1e-12 (rhs bit-identical at 1190 of 1200 samples); parameters identical to PHPS's; contract entry identical. `rhs!` / `outputs!` are type-stable and allocation-free; ForwardDiff Jacobians checked. |
| P4 items 4 and 5 (converters) | after P7 | per the roadmap order of work |
| P5 assembly | done for base | Against PHPS's **compiled C++** `dae_residual` (the function IDA and BDF1 solve): f and g at the equilibrium bit-identical, at 50 random states fault off and on within 1.3e-14 (gate 1e-12; 1228 of 28100 values differ in the last bits, from sin/cos of different math libraries). Layout, wiring, Y-bus, loads, COI and faults identical. Structural Jacobian pattern checked against ForwardDiff. gfl, vsm, droop, voc: `@test_skip` until their models are ported. |
| P6 initialisation | next | |
| P7 to P12 | not started | |

## Needs the user

1. **Upload pack v3** (outward-facing, so not done by the agent): create the GitHub release
   `parity-pack-v3` on Rasoolpey/Porthos.jl and attach `parity/dist/parity_pack-v3.tar.gz`
   (sha256 `b2ceaed1c8b44ed95b1180fc33e2a771b5463f1d08b8c4c91729cb32e335398e`, bound in
   `Artifacts.toml`). Until then CI fails at the parity gates. The tarball is git-ignored
   and must not be regenerated before upload: `archive_artifact` output is not
   byte-reproducible, and `Artifacts.toml` binds this exact file. (Packs v1 and v2 were
   superseded before they were uploaded; nothing to do for them.)
2. **GFL-ZIF has no PHPS DAE reference.** PHPS at ba11ea1 cannot initialise any of the three
   GFL-ZIF system files (`DiracRunner.build`: "coupled network solve did not converge at
   t=0, max |KCL| 5.7e-3 after 40 iterations"). GFL-ZIF is not one of the roadmap's P5 to P7
   parity cases, so those gates are unaffected, but `GFL_ZIF_PHTRUE` is in the P4 active set
   and its component samples come from an initialised case. Options: fix the case or the
   model in PHPS (a new PHPS commit, then a new pack), sample GFL_ZIF_PHTRUE on another case,
   or drop it from the active set.
3. **P3 gate wording** (AGENTS.md: "if a gate looks wrong, explain why and ask"). "Bus
   voltages match ... the case v0/a0 within 1e-8" cannot hold, even for PHPS:
   - the cases store PowerFactory voltages to 6 decimals, so PHPS differs by up to
     4.6e-7 in magnitude and 8.6e-7 rad in angle on the base case;
   - in the five converter cases, `a0` keeps PowerFactory's angle reference (bus 31 = 0)
     but the slack is bus 39 with `a0 = 0`, a uniform shift of 0.1755 rad.

   Suggested wording: "match PHPS within 1e-8; match the case v0 within 1e-6, and a0 within
   2e-6 relative to the slack bus". The test suite already checks this as a supplementary
   test; the literal clause stays `@test_broken` until the roadmap is changed.
4. **Clarabel** is left out of the JuMP extension. Clarabel pins TimerOutputs 0.5, and the
   SciML stack (NonlinearSolveBase via Sundials and NonlinearSolve) needs TimerOutputs 1.x,
   so they cannot resolve together. SCS and COSMO are in. Revisit when Clarabel updates, or
   decide whether it matters for P10 and II.1.
5. **PHPS Hamiltonian at 50 Hz** (for information, no parity case affected): PHPS's
   `hamiltonian` / `grad_hamiltonian` for GENROU and GENSAL compute the base frequency as
   `2*pi*params.get('fn', 60)`, but machine params never contain `fn`, so they use 60 Hz
   even in a 50 Hz case. Porthos uses the machine's `omega_b`, which is the same number at
   60 Hz (all parity cases) and the correct one at 50 Hz.
6. Commit when you are ready (the agent does not commit without being asked). Work is on
   branch `p0-p3-io-network-powerflow`.

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
- Tests: `julia --project=. test/runtests.jl` (about 30 s), or `Pkg.test()`, which is
  slower (it precompiles a fresh test environment). One gate alone:
  `julia --project=. -e 'using Porthos, Test; const ROOT = pwd(); include("test/parity/common.jl"); include("test/parity/p4_components.jl")'`.
  To test against a local pack before binding it: `PORTHOS_PARITY_PACK=parity/pack-v2`.
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
ext/                      PorthosMakieExt (CairoMakie), PorthosJuMPExt (JuMP): empty stubs
parity/generate/generate_pack.py   pack generator (network, powerflow, records, components)
parity/generate/components.py      components section: PHPS init, instrumented kernels, sampling
parity/generate/dae.py             dae section: PHPS's C++ kernel compiled with a residual harness
scripts/bind_parity_pack.jl        tarball + Artifacts.toml binding (hash taken from the tarball)
scripts/setup.ps1                  one-command install
test/unit/, test/parity/p0..p5     unit tests and the P0 to P5 gates
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
- The residual allocates (buffers per call) and uses the dense Y-bus. That is fine for
  parity; P7 needs a preallocated, sparse (CSR, same summation order) version for speed.
- id_dq / iq_dq refresh: PHPS overwrites these outputs between the output and step passes
  with 6-decimal coefficients. Nothing in the parity cases reads them; `assemble` refuses a
  case that does.

## Next steps (P6, then P7 on the base case)

1. P6, initialisation. Port PHPS's chain as `DiracRunner.build` runs it: power flow ->
   `Initializer.run` (each generator's `init_from_phasor`, exciters' and governors'
   `init_from_targets`, the targets dictionary) -> the v0/a0 and load V0/Vini overwrite ->
   DAE-consistent voltages (`_compute_dae_consistent_voltages`, `CoupledNetwork`) ->
   `_rebalance_te_for_fullbus`, `_sync_voltages_to_states`, the closure loop -> Pref sync
   (`pref_from_rest_state`) -> `_equilibrium_polish` (Newton on f = 0). Then replace it by one
   Newton solve of the full DAE for `(x*, V*)` (roadmap). Gate: residual at most 1e-12, and
   `x*`, `V*` within 1e-10 of the pack's `dae/<case>.json` equilibrium, or within
   `||J^-1|| r_PHPS` where PHPS's residual is larger. P6 must also reproduce every
   `init_set` parameter (the P4 and P5 gates take them from the pack until then).
2. P7, simulation: BDF1 first (PHPS's fixed-step scheme in `_emit_main`), then IDA via
   Sundials.jl; events, results writer (`simulation_results.csv` with PHPS's columns),
   `run.json`. Reference trajectories: PHPS's compiled runs (g++ and SUNDIALS are
   installed), at production tolerances, `rtol = atol = 1e-10`, and BDF1 `dt = 5e-4`.
3. Then the converters (P4 items 4 and 5) and P5 to P7 on the GFL, VSM, droop and VOC cases.
