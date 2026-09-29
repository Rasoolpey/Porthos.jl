# Working rules for agents in Porthos.jl

Porthos.jl is the all-Julia successor of PHPS. Your job is to build it, phase by phase.

**Read first:**

- [TODO.md](TODO.md): where the work stands, what needs the user, and the next steps.
  Update it at the end of every session.
- [docs/ROADMAP.md](docs/ROADMAP.md): the plan. Phases P0 to P12 each end with a parity gate;
  Part II follows. Work through the phases in order. "Order of work" at the start of Part I
  says how P2 to P7 are sequenced.
- [README.md](README.md): the design, the architecture and the target repository layout.

## The PHPS reference

Everything from the previous implementation is in
**`C:\Users\em18736\Documents\PHPS_Opt`** (GitHub: Rasoolpey/PHPS_Opt). It is the frozen
parity reference.

| What | Where in PHPS_Opt |
|---|---|
| Case and scenario JSON (the format Porthos reads unchanged) | `phps/cases/` |
| Port contracts, with `domain_clauses` | `phps/src/components/model_port_contracts.json` |
| Component models to port | `phps/src/components/` (`generators/`, `exciters/`, `governors/`, `loads/`, `renewables/`, `network/`) |
| Y-bus, power flow, initialisation, run pipeline | `phps/src/ybus.py`, `powerflow.py`, `initialization.py`, `runner.py` |
| ROA certificates and the physical-state projection | `phps/src/roa/`, `phps/src/certification.py` |
| Existing Julia interval code (sparse jets) | `phps/julia/` |
| Tests, including the analytic 1-D and 2-D certificate cases | `phps/tests/` |
| Stability studies (the sources for P9) | `study/src/`, `phps/tools/` |
| PowerFactory references and certificate records | `study/pf_reference/`, `study/pf_reference/certificates/` |
| Model equations, proofs, Part II plans | `phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md`, `phps/HOW_TO_ADD_PH_COMPONENT.md`, `study/presentation/*.tex` |
| PowerFactory API scripts | `pf/api/` |

Rules for using it:

- Read PHPS_Opt; never modify it. The only PHPS code that runs is the parity-pack generator
  in `parity/generate/` (roadmap section 1), which records the PHPS commit it ran against.
- Don't copy data in advance. Bring a file in only when a phase needs it, unchanged, and note
  the PHPS commit it came from. Large data (PowerFactory references, the parity pack) goes into
  artifacts, not git.
- Never start `PowerFactory.exe` (or its launcher) from the shell. Run and inspect
  PowerFactory only through Porthos's driver in `pf/` (`py -3.10 pf/run.py ...`, or
  `pf_simulate` / `pf_inspect` from Julia; see `pf/README.md`), which works in a copy of the
  study case. `PHPS_Opt\pf\api` is reference only. The engine cannot start while the
  PowerFactory window is open: ask the user to close it.

## Rules for this repository

- Julia does all computing. Python may only edit JSON and call Porthos entry points (roadmap
  section 0 and 2.5). Two exceptions, both computing nothing themselves:
  `parity/generate/` (runs PHPS to make the parity pack) and `pf/` (drives PowerFactory,
  whose API is Python only, and writes what PowerFactory computed).
- Component code:
  - is generic in the number type;
  - uses the limiter primitives, never a bare `if` on a state-dependent value;
  - has an `rhs!` that does not allocate (roadmap 2.2 and 2.3).
- Certificates are proofs, under the certificate discipline in roadmap section 0. Optimiser
  output is only a candidate. Simulation and sampling are never ROA evidence.
- A phase is done only when its parity gate passes. Never loosen a gate just to make it pass.
  If a gate looks wrong, explain why and ask.
- The repository is public. Never commit credentials, unpublished reviewer correspondence or
  local data dumps.
- Do not commit unless the user explicitly asks.
