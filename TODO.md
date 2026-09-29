# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md`. Updated at the
end of every session.

Last session: 2026-09-29. Phases P0 to P3 implemented; the test suite passes
(`736 pass, 6 broken`; the broken ones are the open P3 question below).

## Status by phase

| Phase | State | Gate |
|---|---|---|
| P0 environment + parity pack | done locally | Pack v1 loads and every hash checks locally. **CI cannot fetch the pack until the tarball is uploaded** (see "Needs the user"). |
| P1 input / output | done | All 60 JSON under `cases/` load, validate and round-trip; all 200 PHPS case files also load (optional test, runs when PHPS_Opt is present). |
| P2 network | done | Y-bus (power-flow and DAE), load admittances, Norton stamps and fault shunts equal PHPS **bit for bit** on all 6 parity cases (gate asks for 1e-12). |
| P3 power flow | done, one clause open | `V`/`theta` match PHPS within 4e-15 (gate 1e-8), same iteration counts. The "case v0/a0 within 1e-8" clause cannot hold (below); marked `@test_broken`, not loosened. |
| P4 components | next | |
| P5 to P12 | not started | |

## Needs the user

1. **Upload pack v1** (outward-facing, so not done by the agent): create the GitHub release
   `parity-pack-v1` on Rasoolpey/Porthos.jl and attach `parity/dist/parity_pack-v1.tar.gz`
   (sha256 `fa5fd638f8c0b11b367b98e617240af862ba97f4284031240e65c5e3392deaa0`, bound in
   `Artifacts.toml`). Until then CI fails at the parity gates. The tarball is git-ignored
   and must not be regenerated before upload: `archive_artifact` output is not
   byte-reproducible, and `Artifacts.toml` binds this exact file.
2. **P3 gate wording** (AGENTS.md: "if a gate looks wrong, explain why and ask"). "Bus
   voltages match ... the case v0/a0 within 1e-8" cannot hold, even for PHPS:
   - the cases store PowerFactory voltages to 6 decimals, so PHPS differs by up to
     4.6e-7 in magnitude and 8.6e-7 rad in angle on the base case;
   - in the five converter cases, `a0` keeps PowerFactory's angle reference (bus 31 = 0)
     but the slack is bus 39 with `a0 = 0`, a uniform shift of 0.1755 rad.

   Suggested wording: "match PHPS within 1e-8; match the case v0 within 1e-6, and a0 within
   2e-6 relative to the slack bus". The test suite already checks this as a supplementary
   test; the literal clause stays `@test_broken` until the roadmap is changed.
3. **Clarabel** is left out of the JuMP extension. Clarabel pins TimerOutputs 0.5, and the
   SciML stack (NonlinearSolveBase via Sundials and NonlinearSolve) needs TimerOutputs 1.x,
   so they cannot resolve together. SCS and COSMO are in. Revisit when Clarabel updates, or
   decide whether it matters for P10 and II.1.
4. Commit when you are ready (the agent does not commit without being asked). Everything is
   uncommitted: `Project.toml`, `Manifest.toml` (to be committed per P0), `Artifacts.toml`,
   `src/`, `ext/`, `test/`, `cases/`, `contracts/`, `parity/`, `scripts/`, CI and this file.

## Environment on this machine

One-command setup: `powershell -ExecutionPolicy Bypass -File scripts\setup.ps1 [-RunTests]`
(juliaup, Julia 1.12 as a juliaup override for this directory, `Pkg.instantiate`, the parity
pack, the generator venv from `parity/generate/requirements.txt`). Safe to re-run. It
does not install a C++ compiler or SUNDIALS yet (see below); add them to the script when
P7 needs them.

- Julia 1.12.7 via juliaup (`winget` id 9NJNWW8PVKMN), default channel 1.12. `julia` is on
  PATH (WindowsApps alias).
- `Project.toml` compat `julia = "1.12"`; `Manifest.toml` resolved with 1.12.7; CI runs
  1.12 on ubuntu and windows.
- PHPS runs here: venv at `parity/generate/.venv` (Python 3.13.15, PHPS
  `requirements.txt`; git-ignored). Run PHPS with `PYTHONDONTWRITEBYTECODE=1` from
  `PHPS_Opt/phps` so its tree stays untouched. **No g++ and no SUNDIALS yet**: needed at
  P7 for PHPS's compiled DAE runs (IDA and BDF1 reference CSVs). Options: MSYS2
  (`mingw-w64-x86_64-gcc`, `mingw-w64-x86_64-sundials`), or PHPS's pure-Python solver
  `src/dirac/py_solver.py`, which says it is "identical to the C++ BDF-1".
- PHPS_Opt is at `ba11ea1`, with uncommitted docs-only changes (recorded in the pack
  manifest). The generator refuses to run if PHPS inputs are dirty.
- Tests: `julia --project=. test/runtests.jl` (about 15 s), or `Pkg.test()`, which is
  slower (it precompiles a fresh test environment).

## What exists

```
src/Porthos.jl            module, exports
src/io/expr.jl            safe expression reader (C order, Float64, env of named values)
src/io/json.jl            read/write JSON keeping key order and int/float; json_identical
src/io/schema.jl          JSONSchema validation (cases/schema/*.schema.json)
src/io/case.jl            Case + typed tables with PHPS defaults; ComponentSpec, Wire
src/io/scenario.jl        Scenario, SolverSettings, BusFault/LineFault/OtherEvent (PHPS DAE defaults)
src/io/contracts.jl       ContractSet / ContractEntry / DomainClause
src/io/params.jl          PHPS machine-base normalisation (_normalise_genrou_params) + ctor defaults
src/io/parity.jl          parity pack: artifact lookup, SHA-256 verification, readers
src/network/ybus.jl       Network (sorted bus ids), ybus / ybus_pf / ybus_dae, Norton stamps,
                          load admittances; CPython complex division for bit parity
src/network/events.jl     fault admittance, fault shunts, with_fault, LineFault line split
src/powerflow/newton.jl   exact port of PHPS Newton-Raphson; skip_pf_solve overrides
ext/                      PorthosMakieExt (CairoMakie), PorthosJuMPExt (JuMP): empty stubs
parity/generate/generate_pack.py   pack generator (sections network, powerflow, records)
scripts/bind_parity_pack.jl        tarball + Artifacts.toml binding (hash taken from the tarball)
test/unit/, test/parity/p0..p3     unit tests and the P0 to P3 gates
```

## Findings about PHPS that later phases depend on

- **DAE, not ODE.** Porthos ports PHPS's DAE path (`src/dirac/`: `DiracCompiler`,
  `DiracRunner`; 68 uses in tools and studies). The Kron-reduced ODE path
  (`SimulationRunner`, `compiler.get_z_bus_kron`) is legacy (`tools/run_simulation.py`
  only) and is not ported.
  - `LineFault` exists only in the ODE runner; `DiracRunner._inject_events` ignores it.
    Porthos has the topology (`split_line_for_fault`), but PHPS has no DAE reference for the
    one LineFault scenario. Decide at P7.
  - The one `rk4` scenario cannot run on the DAE path in PHPS either. Decide at P7.
- DAE bus fault: `Y_f = (r - jx)/(r^2 + x^2)` with `z^2 >= 1e-20`; a missing `x` means
  `1e-5` (bolted, PowerFactory-like), a missing `r` means 0.
- DAE Y-bus = lines + shunts + PQ loads as `(P - jQ)/v0^2` (case v0) + generator Norton
  `1/(ra + j xd'')` for every `component_role == "generator"` except GFL and GFL_ZIF.
  VSM, droop and VOC get `ra = 0`, `xd'' = Zseries` (default 0.10) from their constructors.
  `load_G`/`load_B` for the residual are overridden by COMPLEXLOAD `P0`, `Q0`, `V0`.
- Machine params (GENROU/GENSAL families) are converted from the Sn base to the system base
  at load time, with an "already normalised" heuristic, a `D = 2 Sn/Sbase` default when `D`
  is absent, and an `xd'' <= xl` repair. Ported in `src/io/params.jl`; P4 must use
  `component_params`, not the raw JSON.
- PHPS `YBusBuilder` defaults a missing line `x` to 0.001 (`system_graph` uses 0.01); Porthos
  follows the Y-bus.
- Power flow: no Q-limits in PHPS. The unknowns are the angles of non-slack buses, then the
  magnitudes of PQ buses. PHPS's power flow on IEEE-39 converges in 4 iterations to 2e-11.
- The DAE state count on the base case is n_diff = 203 (including delta_COI) and
  n_alg = 78 (Vd, Vq per bus).

## Next steps (P4, then the P5 to P7 slice)

Roadmap "Order of work": the synchronous-machine set first (P4 items 1 to 3), then P5 to
P7 on the base case, then the converters.

1. Extend `generate_pack.py` with a `components` section (pack v2, with a reason in
   `parity/README.md`): for GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE,
   IEEEG3_PHTRUE and COMPLEXLOAD, evaluate `f`, outputs, the Norton injection, `H` and
   `grad H` at 200 random states and inputs (seeded, box around the equilibrium, samples
   near switching surfaces redrawn), recording the limiter sides. Find out how PHPS
   evaluates a single component in Python (C++ snippets via `src/dirac/py_codegen.py` and
   `py_solver.py`, or the symbolic PHS in `get_symbolic_phs`), and use the same code the
   DAE runs.
2. `src/components/interface.jl` (`AbstractComponent`, the roadmap 2.2 functions) and
   `primitives.jl` (`clamp_mode`, `nonwindup`, `select`, `guard_min`, with interval
   decisions and switching surfaces).
3. Port the six types in `src/components/{machines,exciters,governors,loads}/`: generic
   number type, no bare `if` on state, non-allocating `rhs!`, with an allocation test each.
   Move the Norton rule from `network/ybus.jl` into each model's `injection`.
4. Then P5 (assembly: wiring from `connections`, state order = PHPS `state_offsets`, COI,
   reservoirs, the residual `F`/`G`), P6 and P7 on the base case.
