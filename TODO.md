# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md` (its Part I banner:
parity is frozen). Updated at the end of every session.

Last reviewed: 2026-10-01.

## Goal

A **rigorous region of attraction (ROA)** for the power system, built on a **storage
function that truly decreases**, ideally the physical port-Hamiltonian energy extended by
well-founded terms. Two layers (decided 2026-10-01):

1. **Now:** certify the ROA with the local quadratic candidate `V_P` (the Lyapunov equation
   on the common-angle section). It exists, covers every physical state, and is the test bed
   for the rigorous interval machinery (P11).
2. **Target:** a nonlinear composite storage
   `V_ext = H_machines,s + U_net,s + V_controllers + V_cross`, certified by the same pipeline
   and then the primary ROA certificate. Call it `H_ext` only where every added term has a
   physical or port-Hamiltonian derivation; otherwise call it `V_ext` and say which parts are
   physical and which mathematical. Do not treat a `V_P` certificate as the answer to the
   physical-storage question.

The user's standing preference: keep the component dynamics unchanged if at all possible
(storage-only route: new storage terms, joint storages, other supply rates or ports);
physical model changes (two-way reservoirs) are the fallback. Control-method choices are
discussed with the user before they are implemented.

## Decisions (2026-10-01, user with a reviewer)

- **Candidate-independent certificate pipeline.** A Lyapunov candidate provides: its value,
  gradient, interval evaluation, equilibrium shift, physical-coordinate projection and a
  fingerprint. First implementation `QuadraticCandidate(P)`; second
  `ExtendedStorageCandidate(H_ext)`; both go through the same proof gates.
- **Network: hybrid.** Keep the bus voltages and the network explicit when deriving storage,
  port cancellation and physical meaning. In the proof, validate a unique KCL branch and
  eliminate the voltages by the interval implicit-function argument. The dense Schur
  complement `f_x - f_V g_V^-1 g_x` must not be mistaken for physical cross-machine storage.
- **Limiters: single mode first.** The certified ROA must stay inside one decided branch of
  every limiter and domain clause (the branch primitives raise `UndecidedBranch` otherwise).
  Multi-mode or boundary certificates only once the smooth-region result is useful.
- **Parity is frozen** as a regression guardrail (pack v6); PHPS is not ground truth; no
  more parity refinement (see the roadmap banner). PHPS's certificate numbers (the P11 gate)
  are sanity references only.
- Route B (the structure-search LMI) is a **diagnostic** for missing storage blocks, not a
  source of physical terms.

## Next steps, in order

1. ~~Record the certificate decision~~ and ~~add the roadmap override~~ (done 2026-10-01).
2. **Publish pack v6** (the one loose end, below), with the user's go-ahead, and verify the
   download; this restores CI. Then no more parity-pack work.
3. **Generic Lyapunov-candidate interface with the quadratic `V_P`** (the first code task; a
   new `src/roa/`). The candidate API as decided above; `QuadraticCandidate(P)` with `P`
   from the Lyapunov equation on the section (`physical_projection`, `reference_section` and
   the dense `lyap` solution in `scripts/storage_search.jl` show how), shifted to the
   equilibrium from `solve_equilibrium`; interval evaluation of `V` and of `grad V' f`
   through the generic component code (the models are generic in the number type; interval
   methods for the branch primitives, returning the decided branch or throwing
   `UndecidedBranch`, belong here); the KCL branch and the voltage elimination by the
   interval implicit-function argument; the single-mode containment audit with the contract
   `domain_clauses`. Roadmap P11 lists the pieces (equilibrium Krawczyk, first-order and
   centered hulls, definiteness test, mean-value gate, records, `ROACheck`). Probe small and
   bound every long computation (memory: bounded long runs).
4. **Nonlinear port-power residual audit** (the missing P10 audit): for each component,
   reconstruct independently its supply (port power), internal storage derivative,
   dissipation and the network cancellation, and check `grad H' f = supply - dissipation`
   at random states and along trajectories. This gives Route A' reliable signs and units.
5. **Polar ports and explicit `U_net`** (Route A'): derive the conjugate polar supply from
   the network energy balance (do not assume (P, omega) and (Q, |V|) are the right pairs);
   test the lossless network first, then quantify the passivity shortage from the
   conductances and the active (constant-power) loads. The cheapest remaining test of
   whether unchanged dynamics admit a physically structured `H_ext`.
6. **Construct `V_ext`**: machine energies + `U_net` + controller terms + cross terms, each
   with a stated origin. Use Route B only to diagnose missing blocks, after checking that its
   blocks survive a change of reference angle (or an orthonormal COI basis) and the explicit
   `U_net`. Controllers that cannot be passive alone (IEEEG1: relative degree 2; IEEEG3:
   right-half-plane zero) need a joint machine-governor storage or another port.
7. **Certify `V_P` and `V_ext`** through the same single-mode pipeline and compare the sets.

## Loose ends (no thinking needed; step 2)

- **Parity pack v6 is not published**, but `master` (pushed, `78ae66f`) binds it, so CI fails
  until the release `parity-pack-v6` exists. Tarball `parity/dist/parity_pack-v6.tar.gz`,
  sha256 `400974ee39cb16ac45c9628916facbf21583fa225bad0af4d863d7219b050136`, tree hash
  `c99fe023705c11cc10260a9dcbc608421daa3cb8`. Publishing is outward-facing: ask the user,
  then upload as for v3 to v5 (GitHub REST API with the stored git credential; there is no
  `gh` on this machine) and verify both hashes on the downloaded file.
- Committed on 2026-10-01 after a green local suite (2926 passes, 2 intentional GFL skips,
  no failures): converter IDA comparison report-only in P7 (BDF1 is the check there); the
  droop P6 residual exception (2e-12 on its two measured-current rows, every other row
  1e-12); `vi_bisection` through the primitives with `NoModes()`; VOC `pvoc_mode = 1`
  rejected at construction (no reference); the `docs/ROADMAP.md` banner; this file. Not
  pushed: `origin/master` is at `78ae66f`.
- GFL / GFL_ZIF wait for the GFL reservoir model (part of the modelling work).
- Later software: P8 reports, P9 studies (CCT sweeps etc.), P12 Python wrappers and docs.

## Evidence so far: storage and passivity (base case, IEEE-39, at the equilibrium)

The port is done for everything the stability work needs: I/O, network, power flow, the
synchronous-machine set (GENROU, GENSAL, IEEET1, IEEEG1, IEEEG3, COMPLEXLOAD), the
grid-forming converters (GFM_VSM, GFM_DROOP, GFM_VOC), assembly, initialisation, BDF1 and IDA
simulation, the PowerFactory driver, and the P10 port-Hamiltonian audit tools. Sources for
the modelling work: the supervisor's review
(`PHPS_Opt/study/presentation/response_to_reviewer_component_models.tex`),
`PHPS_Opt/phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` (0.2b, 0.2c, Target B), and roadmap
Part II.2 (B1 to B7).

1. **The assembled storage is not a Lyapunov function.** `shifted_storage_audit`: Hess H has
   rank 54 of 171 on the common-angle section; sym(SA) has 41 positive eigenvalues (up to
   12.9). The reservoirs are one-way accounts fed at constant power (no storage along plant
   directions); controller states carry no storage.
2. **Controller ports.** IEEEG1 at speed -> Tm is non-passive for any storage (relative
   degree 2: 6 poles, 4 zeros; Re G(jw) changes sign at 1.9347 rad/s). IEEEG3 has a zero at
   +1.333 1/s (Re H < 0 on 1.32 to 14.9 rad/s). IEEET1 not settled. No realisation changes
   these input-output facts, so they need a joint storage or another port pairing.
3. **Governor loops** (`scripts/governor_passivity.jl`): one at a time, 8 IEEEG1 pass against
   the rest of the system; IEEEG1_10's rest is itself non-passive near 1.55 rad/s (a lightly
   damped zero pair at -0.037 +- 1.421i); IEEEG3_11 passes only with a frequency-dependent
   split. **Jointly, no storage split port by port at the Tm ports exists**, static or
   dynamic (`scripts/joint_governor_storage.jl`): fails at 0.048 rad/s (the slow system mode)
   and at 6.59 rad/s. All machines have D = 0; the swing damping is in the damper windings.
   At low frequency the lossy network makes differential speed directions non-passive
   (the path dependence of U_net).
4. **Route A, the cut at the machine terminals with the power V*I**
   (`src/ph/terminal.jl`, `scripts/terminal_passivity.jl`): every unit, and the bare machine,
   is non-passive at its terminal (a constant-power source at DC, indefinite Herm Y(0),
   eigenvalues +-5.5 to +-26; non-passive near its swing frequency; the AVR makes DC much
   worse, -43 to -3422); the network with its loads is non-passive at DC (-18.5). Synchronous
   machines are not incrementally passive in (V, I) terms. Not tried: A', the same units
   with polar ports (P, omega) and (Q, |V|).
5. **Route B, quadratic storage on the reduced dynamics** (`src/ph/structure.jl`,
   `scripts/storage_search.jl`, decay margin gamma(S) = min lambda_max(As'P + PAs) over P on
   a pattern, P >= 0, tr P = 1; dual certificates and Lyapunov checks in interval
   arithmetic):
   - each component alone, or each unit (machine + governor + exciter) alone: margin 0 (at
     most 2.4e-9; the full storage reaches -4.1e-3). A block-local quadratic storage does
     not decay by any useful margin.
   - adding speed-speed, speed-angle or machine-angle couplings across units: still 0.
   - coupling all machine states across machines: a verified Lyapunov function, margin
     -2.0e-6. A dual-guided greedy finds verified patterns from 2233 free entries (-8.4e-7)
     to 2686 (-1.8e-5); it picks machine-machine blocks and above all GENROU_1 (bus 39: the
     reference angle, no controllers, D = 0, the largest inertia). The sparsest-set (L1)
     stage did not converge in its 300 s cap.
   - **Caution before reading it physically**: As = f_x - f_V g_V^-1 g_x eliminates the
     network, so cross-machine terms may be network energy, not joint machine storage; and
     GENROU_1's role may depend on it being the eliminated reference angle.

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
- PHPS (read-only reference, only for regenerating the parity pack): venv at
  `parity/generate/.venv`, g++ and SUNDIALS at `C:\msys64\ucrt64`; details in
  `parity/README.md`.
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
- Tests: `julia --project=. test/runtests.jl` (about 12 min, mostly the P7 runs of four
  cases). One file alone: `julia --project=. -e 'using Porthos, Test; const ROOT = pwd();
  include("test/parity/common.jl"); include("test/unit/ph.jl")'`.

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

