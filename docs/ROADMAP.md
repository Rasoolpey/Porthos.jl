# Porthos.jl roadmap: an all-Julia port-Hamiltonian power system simulator and stability toolbox

Status: plan, 2026-09-29. Porthos.jl (Port-Hamiltonian Operation and Stability) is a new
repository. The current PHPS repository (Python/C++, PHPS_Opt) is not modified; it is the
**frozen parity reference** that every Julia stage is checked against. In this document
"PHPS" always means that reference, and "Porthos" the Julia implementation.

---

## 0. Scope and working rules

**What goes in Julia: everything that computes.** That covers:

- reading the case and scenario JSON (network, components, parameters, wiring, events, solver
  settings);
- building the system (Y-bus, component models, wiring, state layout);
- power flow and initialisation to an equilibrium;
- time-domain simulation with events, and saving results;
- all stability studies (CCT, indices, voltage and frequency studies, energy methods);
- port-Hamiltonian audits;
- the rigorous ROA certificates;
- design optimisation;
- reports and figures.

**What Python may do, and nothing more:**

1. create or edit JSON files (case variants, scenario sweeps);
2. call Julia entry points to run a study that already exists in Julia, and print the path
   of the record it wrote.

When a run starts, no numerics happen in Python: no NumPy computation, no plotting, no
validation. The Python wrappers are a few lines each (`juliacall`, or a subprocess to a Julia
script).

**Order.** Part I builds the Julia implementation to parity with PHPS. Part II then continues
the modelling and stability programme (larger certified ROA, making `H` a Lyapunov function,
converter physics) in Julia only.

**Rules carried over from PHPS:**

- **ROA/stability certificates are proof-based:** Lyapunov and energy functions with interval or
  validated enclosures. They use no simulation-based boundaries and no CCT sweeps as ROA
  evidence.
- **Certificate discipline:** outward rounding, same-box identity between gates, fingerprints of
  case, model and candidate, and `verified_valid_level` only when every gate passes.
  Optimisation (JuMP) output is a candidate only.
- **PowerFactory stays the external reference** for model validation; PHPS stays the internal
  reference for parity.

---

## 1. The PHPS reference, and what to take from it

Everything from the previous implementation stays in PHPS_Opt:

- local checkout: **`C:\Users\em18736\Documents\PHPS_Opt`**
- GitHub: [Rasoolpey/PHPS_Opt](https://github.com/Rasoolpey/PHPS_Opt)

It is read and never modified. Nothing is copied in advance. A file comes into this
repository only when a phase needs it, unchanged, with the PHPS commit it came from recorded.
Small inputs go into git; large data goes into content-hashed artifacts bound in
`Artifacts.toml` (tarballs on a GitHub release).

| Item | Path in PHPS_Opt | Use | Brought into Porthos |
|---|---|---|---|
| Case and scenario JSON | `phps/cases/**` (system files `system*.json`, scenario files with `solver`/`events`/`output`/`plots`) | inputs, unchanged format | at P1, into `cases/` (git) |
| Port contracts | `phps/src/components/model_port_contracts.json` (schema `2.1-reservoir-corrected`, with `domain_clauses`) | component contracts and the fail-closed domain audit | at P1, into `contracts/` (git) |
| PowerFactory references | `study/pf_reference/**` (about 81 MB) | external validation targets | at P7, as an artifact |
| Certificate records | `study/pf_reference/certificates/*.json` | ROA parity targets | at P0, inside the parity pack |
| Model code | `phps/src/` (`components/`, `ybus.py`, `powerflow.py`, `initialization.py`, `runner.py`, `roa/`, `certification.py`) | what each phase ports | read in place |
| Julia interval code | `phps/julia/` (`ROAIntervalJets.jl`, `SecondOrderJets.jl`) | starting point for `src/roa/` | ported at P11 |
| Tests | `phps/tests/` (including the analytic 1-D and 2-D certificate cases) | test ideas and analytic cases | read in place |
| Study scripts | `study/src/`, `phps/tools/` | the source of each P9 study | read in place |
| Documentation | `phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md`, `phps/HOW_TO_ADD_PH_COMPONENT.md`, `phps/README.md`, `study/presentation/*.tex` | model equations, proofs, Part II plans | read in place |

**Parity pack.** Generate one set of reference numbers from PHPS for a fixed set of cases.
The generator scripts live in this repository under `parity/generate/`, so PHPS_Opt stays
unmodified; they run once against a PHPS checkout at a recorded commit. They are the only
place in Porthos where PHPS code runs, and they are exempt from the Python lint rule (2.5)
because their job is to compute the reference. The pack must contain:

- for each active component type: `f`, the outputs, the Norton injection, `H` and `grad H` at 200
  random states and inputs, with the inputs recorded and the side of every limiter recorded.
  States are drawn in a physically meaningful box around the equilibrium, and samples within
  rounding distance of a switching surface are redrawn;
- for IEEE-39 base, GFL, VSM, droop and VOC: the Y-bus, power-flow voltages, the equilibrium
  `x*` and `V*` with the PHPS DAE residual at `x*`, state names in kernel order, and `f`, `g`
  at 50 random states;
- for the bus-16 fault (150 ms) on the same cases: `simulation_results.csv` with IDA at the
  production tolerances, IDA at tight tolerances (`rtol = atol = 1e-10`), and BDF1
  (`dt = 5e-4`), with the Newton tolerance of each run and how PHPS handles limiter switching;
- the certificate records (`2.69e-12`, and `5.16e-10` with the centered hull) and the 83
  contract-clause rows;
- the PH audit numbers of P10 and the study anchors listed in P9 below;
- a `MANIFEST.json` with the PHPS commit, tool versions, seeds and a SHA-256 per file.

The parity pack is read-only and is the regression baseline. A new baseline is a new artifact
with a new hash and a written reason in `parity/README.md`.

Component types to port. The active set comes first, because it is what the IEEE-39 study uses:

- machines: `GENROU_PHTRUE`, `GENSAL_PHTRUE`
- exciter: `IEEET1_PHTRUE`
- governors: `IEEEG1_PHTRUE`, `IEEEG3_PHTRUE`
- load: `COMPLEXLOAD` (contract key `ComplexLoad`)
- converters: `GFL_PHTRUE`, `GFL_ZIF_PHTRUE`, `GFM_VSM_PHTRUE`, `GFM_DROOP_PHTRUE`, `GFM_VOC_PHTRUE`
- network: `PiLine`, `Transformer2W`, shunts
- not ported for now: the older models in PHPS `components/retired/` (`GENROU_PHS`,
  `IEEEX1_PHS`, `TGOV1_PHS`, `IEEEST_PHS`, ...). They stay in PHPS_Opt as reference. Five
  PHTRUE classes are Python subclasses of retired ones (GENROU_PHTRUE of `GenRouPHS`,
  GENSAL_PHTRUE of `GenSal`, IEEET1_PHTRUE, IEEEG1_PHTRUE and IEEEG3_PHTRUE of their `_PHS`
  parents); each is ported as one flat model whose equations are the inherited code plus the
  PHTRUE changes, and the parent is read only to find out what the PHTRUE model computes.

Event types in use: `BusFault` (about 1,650 scenarios) and `LineFault` (one). Solver methods in
use: `ida` (about 1,390 scenarios), `bdf1` (about 200) and `rk4` (one).

---

## 2. Architecture

### 2.1 Repository layout

The annotated tree is in the README ("Repository layout"). In short:

```
Porthos.jl/
  Project.toml, Manifest.toml         pinned environment (Manifest committed from P0)
  Artifacts.toml                      content hashes of the parity pack and PF references
  src/
    Porthos.jl                        top module, public API
    io/          JSON schema + loaders (case, scenario, contracts), parameter expressions,
                 results writers (CSV compatible with PHPS, Arrow/JLD2), run metadata + hashes
    network/     buses, lines, transformers, shunts, Y-bus, fault and topology events
    components/  interface.jl, primitives.jl (2.3), and one file per model type in
                 machines/, exciters/, governors/, loads/, converters/
    assembly/    wiring graph from `connections`, state layout and names, DAE residual,
                 sparsity pattern, COI reference, one-way (reservoir) coordinates
    powerflow/   Newton-Raphson AC power flow
    init/        equilibrium: power flow -> component init -> full DAE Newton
    sim/         DAE solvers (IDA via Sundials.jl, fixed-step BDF1, OrdinaryDiffEq mass-matrix
                 methods), event callbacks, re-initialisation of algebraic states
    studies/     CCT, indices, voltage and frequency studies, energy methods, sweeps
    ph/          Hamiltonians, PH structure checks, shifted-storage audit
    roa/         interval core, jets, certificates, containment audit, independent check
    precompile.jl  PrecompileTools workloads
  ext/
    PorthosMakieExt/  report: summaries, figures (CairoMakie), scenario `plots` blocks
    PorthosJuMPExt/   design (conditioning, hull-aware LMI, PH-informed storage, OPF-type)
                      and the KYP / passivity programs
  cases/                              copied inputs (2.4), with schema/
  contracts/                          model_port_contracts.json
  parity/                             pack description; generate/ builds the pack (section 1)
  test/                               unit/, parity/ (one file per gate), lint/
  bench/                              timing suite for the targets in 2.6
  python/                             thin wrappers only (2.5)
  scripts/                            Julia command-line entry points
  docs/                               this roadmap, Documenter pages
```

The report and design code sit in package extensions so that `using Porthos` does not load
CairoMakie, JuMP or the SDP solvers. `src/` declares the plotting and design functions; the
extensions implement them. This keeps the per-run load time low, which is the point of the
rewrite.

### 2.2 One model code for simulation and proof

The central design decision: each component is written once, as plain Julia functions that
are generic in the number type. The same `f` then runs with:

- `Float64`, for simulation and power flow;
- `ForwardDiff.Dual` or sparse automatic differentiation, for Jacobians and Newton;
- interval types and first- and second-order interval jets, for the ROA certificates;
- symbolic types, where exact structure is wanted (Symbolics.jl), for documentation and checks.

This removes the current split between the C++ kernel, the SymPy parser and the generated
Julia code, and with it every parity layer between them.

Component interface (abstract type `AbstractComponent`):

| Function | Meaning |
|---|---|
| `states(c)`, `params(c)`, `ports(c)` | declared state names, parameters and input/output ports (as in `PowerComponent`) |
| `rhs!(dx, c, x, u, p)` | differential equations, generic in the element type |
| `outputs!(y, c, x, u, p)` | port outputs (for example `Efd`, `Tm`, currents) |
| `injection(c, x, V, p)` | Norton or current injection into KCL, and whether it contributes admittance |
| `initialize(c, targets)` | states from power-flow targets |
| `hamiltonian(c, x)`, `grad_hamiltonian(c, x)` | declared storage |
| `observables(c)` | derived signals written to results (`delta_deg`, `H_steam`, ...) |
| `contract(c)` | the port-contract entry, including `domain_clauses` |
| `modes(c, x, u)` | the limiter and branch decisions (2.3) |

### 2.3 Limiters and branches are explicit

Every limiter, clamp, anti-windup or non-windup test and guard is written with a small set of
primitives (`clamp_mode`, `nonwindup`, `select(cond, a, b)`, `guard_min`), not with bare `if`.
Each primitive:

- evaluates normally on `Float64`;
- on intervals, returns the decided branch and its margin, or raises "undecided". This is what
  the Stage 2.3 containment audit needs.
- reports its switching surface, for event detection in the simulator. For the P7 parity
  gate the simulator handles limiter switching exactly as PHPS does (recorded in the parity
  pack); root-finding on switching surfaces is an option added after the gate passes.

### 2.4 Inputs and outputs

- **Inputs:** the current JSON formats are kept unchanged:
  - system files: `config`, `Bus`, `PQ`, `PV`, `Slack`, `Line`, `Shunt`, `components`,
    `connections`;
  - scenario files: `system`, `solver`, `events`, `output`, `plots`.

  A JSON schema is written for each, and the loader validates against it. Parameter strings
  such as `"2.0 * M_PI * 60.0"` are parsed by a small safe expression reader, never `eval`.
- **Outputs:** `simulation_results.csv` with the same column names as PHPS, so the existing
  study scripts and PowerFactory comparisons still apply. There is also a binary copy
  (Arrow or JLD2) and a `run.json` with the case, scenario, package-source and Manifest hashes,
  solver settings and wall time.

### 2.5 The Python layer

`python/porthos/` imports `juliacall`, activates the project and exposes, for example,
`simulate(scenario_path)`, `cct(case, bus, ...)` and `certify_roa(case, ...)`. Each returns
the record path. A lint test enforces the rule: modules under `python/` may not import NumPy,
SciPy or matplotlib, except in the JSON-editing helpers. `parity/generate/` is outside
`python/` and is not covered by the rule (section 1).

### 2.6 Performance measures

- `PrecompileTools` workloads for loading, simulation and certificates, with an optional
  PackageCompiler system image, so there is no compile delay per run.
- A persistent session for interactive work and sweeps.
- Type-stable, allocation-free right-hand sides, with sparse Jacobians coloured from the
  component structure.
- Threads or `Distributed` for scenario sweeps and per-row certificate loops.

Timings of the current pipeline, which these measures target (2026-09-29, IEEE-39):

| Stage | Current time |
|---|---|
| SymPy DAE build from the C++ kernels | about 28 s per run |
| Equilibrium Krawczyk (Python) | about 26 s |
| First-order hull | 15.5 s (Python), 0.4 s (Julia) |
| Centered hull | about 6 s (Julia) |
| Structured SDP (cvxpy/SCS) | not finished in 25 min |

---

## Part I: Julia implementation

Each phase ends with a **parity gate** against the parity pack. A phase is done only when its
gate passes in CI.

Numerical comparisons use `|a - b| <= atol + rtol * |b|`, with `atol` set per quantity from
its scale, so that values near zero or near a cancellation do not fail a purely relative test.

**Order of work: a thin end-to-end slice first.** P2 to P7 are done first for the IEEE-39 base
case with the synchronous-machine set only (P4 items 1 to 3). Assembly, state ordering,
initialisation and event bugs then show up before the converters are written. The converters
(P4 items 4 and 5) follow, and the P5 to P7 gates are re-run on the GFL, VSM, droop and VOC
cases.

### P0. Environment and parity pack
- Pin the Julia version and project dependencies:
  - core: JSON3, StructTypes, SparseArrays, ForwardDiff, DifferentiationInterface with
    SparseConnectivityTracer and SparseMatrixColorings (sparse Jacobians; SparseDiffTools is
    deprecated), Sundials.jl, OrdinaryDiffEqBDF and OrdinaryDiffEqRosenbrock (not the full
    OrdinaryDiffEq), NonlinearSolve, IntervalArithmetic, Arrow, JLD2 and PrecompileTools;
  - extensions (weak dependencies): CairoMakie; JuMP with SCS, COSMO and Clarabel.
- Commit the Manifest.
- Generate the parity pack (section 1), bind it in `Artifacts.toml`, and import it.
- **Gate:** CI runs; the parity-pack loader reads every file and checks every hash.

### P1. Input and output
- Case and scenario loaders with schema validation, the contract loader and the parameter
  expression reader.
- **Gate:** every JSON under `cases/` loads and round-trips.

### P2. Network
- Y-bus from lines, transformers (tap, phase) and shunts; the bus map; fault admittances.
- **Gate:** Y-bus equal to PHPS within `1e-12` on all parity cases.

### P3. Power flow
- Newton-Raphson with PV/PQ/slack handling and Q-limit options as in PHPS.
- **Gate:** bus voltages match PHPS and the case `v0`/`a0` within `1e-8`.

### P4. Components
Port one type at a time, active set first:

1. `GENROU_PHTRUE`, `GENSAL_PHTRUE`
2. `IEEET1_PHTRUE`, `IEEEG1_PHTRUE`, `IEEEG3_PHTRUE`
3. `COMPLEXLOAD`
4. `GFL_PHTRUE`, `GFL_ZIF_PHTRUE`
5. `GFM_VSM_PHTRUE`, `GFM_DROOP_PHTRUE`, `GFM_VOC_PHTRUE`

For each type, implement `rhs!`, `outputs!`, `injection`, `initialize`, `hamiltonian`,
`observables`, `contract` and `modes`, including the corrected constant-power reservoir law.

- **Gate (per type):** at the 200 random states of the parity pack, the limiter branches agree
  first (both sides exercised); then `rhs`, outputs, injection, `H` and `grad H` agree with
  `rtol = 1e-12` and the per-quantity `atol`. The contract entry is identical.

### P5. Assembly
- Wiring graph from `connections`, port resolution, state layout and names in PHPS order (so
  CSV columns and records align), the COI reference angle, and one-way reservoir coordinates.
- The DAE residual `F(x, V)` / `G(x, V) = 0` (KCL), with its sparsity pattern.
- **Gate:** `f` and `g` equal the parity pack at 50 random states per case within `1e-12`.

### P6. Initialisation
- Power flow, then component `initialize`, then a Newton solve of the full DAE for `(x*, V*)`,
  replacing PHPS's multi-pass refinement.
- **Gate:** on all parity cases, the Porthos DAE residual is at most `1e-12`, and `x*` and `V*`
  match PHPS within `1e-10` or, where the PHPS residual is larger, within the distance that
  residual allows (`‖J⁻¹‖ · r_PHPS`, computed in the test). The single Newton solve may be more
  accurate than PHPS's multi-pass refinement, and the gate must not fail for that reason.

### P7. Simulation
- IDA (Sundials.jl) with the PHPS tolerances, fixed-step BDF1 at `dt = 5e-4`, and optionally
  OrdinaryDiffEq's FBDF or Rodas5P with a mass matrix.
- Events: `BusFault` (fault admittance on a time window), `LineFault`, and line trip and load
  step for later. Algebraic states are re-initialised after each switch.
- Results writer and `run.json`.
- **Gate:** bus-16 fault trajectories:
  - BDF1 against PHPS BDF1: at most `1e-9` (same scheme, same Newton tolerance, same switching
    times, all taken from the parity pack);
  - IDA against PHPS IDA, both at `rtol = atol = 1e-10`: at most `1e-6` over 15 s. At production
    tolerances the step sequences differ, so two IDA runs are not expected to agree that
    closely;
  - at production tolerances, against the PowerFactory references, with the same metrics and
    thresholds the current validation uses.

### P8. Reports
- The scenario `plots` blocks, figure style, and a summary per run (CairoMakie).
- **Gate:** figures regenerate from records alone.

### P9. Studies
Port in this order. Each has a gate that reproduces a recorded number or table.

| Study | PHPS source today | Gate |
|---|---|---|
| CCT by fault series and bisection | `study/src/cct_*`, `phps_cct_triple.py`, `build_cct_series_cases.py` | recorded CCT anchors (bus-16 PF anchor and the sweep tables) |
| Transient energy and gate wall | `tools/transient_energy_cct.py`, `analytic_walls.py`, `uep_reduced_wall.py`, `study/src/audit_gate_wall39.py` | recorded wall levels and rankings |
| Frequency indices (RoCoF, nadir) | `tools/frequency_indices.py`, `study/src/mode_rocof_*` | recorded index tables |
| Voltage: static Q-V, snap, LVRT recovery | `tools/voltage_stability.py`, `voltage_snap.py`, `voltage_recovery.py`, `study/src/qv_*` | recorded Q-V margins (for example `6.055 -> 12.383 pu`) |
| Energy margin and reservoir work | `tools/compute_energy_margin.py`, reservoir audits | recorded `D_r` and energy records |
| Stability-index and design studies (reservoir distribution, OPF comparators, control modes) | `study/src/reservoir_*`, `opt_*`, `opf_only_*`, `control_mode_*` | recorded study summaries |

### P10. Port-Hamiltonian audits
- Hamiltonians and gradients, the shifted-storage audit (Hessian rank, `sym(SA)`, exact
  nonlinear `dH/dt`), the physical projection and common-angle quotient, and the KYP and
  passivity tests (JuMP).
- **Gate:** reproduce rank `54/171`, 41 positive eigenvalues and max `+12.882693`; IEEEG1
  `Re(-G(jw))` crossing at `1.934718 rad/s` and KYP infeasibility on all nine sets.

### P11. ROA certificate pipeline
Components:

- the interval core (IntervalArithmetic.jl) with first- and second-order sparse jets, generic
  over the same component code;
- the candidate `V_P` (Lyapunov equation on the quotient);
- the equilibrium Krawczyk, KCL branch and first-order and centered hulls;
- the interval definiteness test and the mean-value gate with bisection;
- the containment audit, using the `modes` primitives and the fail-closed `domain_clauses`;
- records and fingerprints;
- `ROACheck`, an independent second check: it recomputes digests and runs a different
  definiteness method (interval `LDL^T`).

**Gate:**

- the first-order hull reproduces `derivative_verified_level = verified_valid_level = 2.69e-12`
  for both candidate `P`;
- the centered hull reproduces decay at `5.16e-10`;
- the same 83 clauses pass.

### P12. Python wrappers and retirement of the old pipeline
- `python/porthos/` exposes the Julia entry points, and the lint rule of 2.5 applies.
- All studies are run from Julia. The old PHPS remains only as the source of the parity pack.
- **Gate:** a fresh clone on a clean machine reproduces the parity suite with one command.

---

## Part II: Modelling and stability programme (after Part I)

The details are in `phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` in PHPS_Opt (§0.2,
§0.2b, §0.2c). It is carried out in Porthos only.

### II.1 Larger certified ROA (Target A)
1. Stage 2.3 containment at the centered level `5.16e-10`, which makes it a certified ROA; then
   the conditioning-SDP `P`.
2. Structured hull-aware LMI for `P`, solved with JuMP and sparsity (block or chordal) instead of
   dense 342-dimensional cones; re-certify the result.
3. Fault reach without simulation: prove `dot(V_P) <= a V_P + b` on the fault-on field on
   `{V_P <= c_f}`, with the fault-on branch, modes and domains certified. The comparison
   solution then gives a certified CCT lower bound.
4. Falsification evidence (gate 6), once the certified set is large enough to matter.

### II.2 Making `H` itself a Lyapunov function (Target B; answers review comments 3 and 4)
- **B1:** add the network potential `U_net` to the storage, which gives curvature along the
  rotor angles.
- **B2:** a PH-informed storage LMI. It fixes the physical energy blocks, frees the controller
  and cross-term blocks, and finds the sparsest set of cross-terms that makes the storage work.
- **B3:** turn the required blocks into nonlinear storage terms (controller integrators and
  lags, machine-governor and machine-exciter couplings, the PLL potential `K (1 - cos theta)`),
  each with a stated origin.
- **B4:** repair the decay. In order: passivity indices with a joint machine-governor storage; a
  general quadratic supply (KYP/IQC); a new port; a physical steam or hydraulic extension; a
  controller redesign only as a last resort. IEEEG3 (right-half-plane zero) starts at the port or
  extension step. IEEET1 needs a full-matrix and joint machine-exciter storage.
- **B5:** close the network port with `Y = N + D`; incremental supply rates for constant-power
  loads; the internal-reactance closure.
- **B6:** converter layers: VSM, droop and VOC storage and terminal inequalities; GFL PLL, PI,
  lag and filter storage.
- **B7:** certify the non-quadratic `H_ext` with the same gates (outer box from a quadratic lower
  bound, centered enclosure of `grad H_ext^T f`).

### II.3 Physical converter models (future work)
Add the dc link `(1/2) C_dc v_dc^2`, the LCL filter and the inner current loop as real PH storage,
replacing the reservoir account. This is a new model class with EMT-type time scales and needs
its own validation.

---

## 3. Testing and validation discipline

- **Unit tests:** per component, per primitive (limiter decisions on intervals) and per solver.
- **Parity tests:** every P-gate above, in CI, against the read-only parity pack.
- **External validation:** the PowerFactory comparisons with the current metrics.
- **Certificates:** analytic test systems (the 1-D and 2-D cases already in the PHPS tests) plus
  the IEEE-39 records.
- **Records:** every run and study writes `run.json` with case, scenario, source and Manifest
  hashes, so any number can be traced to what produced it.

## 4. Risks and how they are handled

| Risk | Handling |
|---|---|
| Limiter semantics differ from the C++ kernel | explicit primitives (2.3) and per-type parity at states that exercise both branches |
| Solver differences hide model errors | compare BDF1 fixed-step first (same scheme, tight tolerance), IDA second |
| Generic code loses type stability and gets slow | `@code_warntype` and allocation tests in CI on the RHS |
| Interval evaluation of generic models is slow | the sparse jets already written for PHPS (`ROAIntervalJets.jl`, `SecondOrderJets.jl`) move into `roa/` |
| Parity baseline drifts | the parity pack is generated once, hashed, bound as an artifact and read-only |
| Parity gates fail on rounding, not on model errors | mixed `atol`/`rtol` tests; samples kept away from switching surfaces; branch agreement checked before values; tight-tolerance IDA reference |
| Heavy dependencies make every run slow to start | CairoMakie and JuMP with its solvers as package extensions; OrdinaryDiffEq sublibraries only; PrecompileTools and an optional system image |
| Bugs in assembly or initialisation found late | the thin end-to-end slice (Part I, order of work) before the converters |
| Python creeps back into computation | the lint rule in 2.5, plus review |
