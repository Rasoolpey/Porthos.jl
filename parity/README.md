# Parity pack

The parity pack holds the reference numbers from PHPS that every Porthos parity gate is
tested against (roadmap section 1). It is read-only: it is generated once, content-hashed,
and bound as the lazy artifact `parity_pack` in [`Artifacts.toml`](../Artifacts.toml). The
test suite installs it and `Porthos.load_parity_pack()` checks the SHA-256 of every file
against `MANIFEST.json` and rejects files the manifest does not list.

A new baseline is a new artifact with a new hash, and a new entry below saying why.

## Generating

`generate/generate_pack.py` is the only place in Porthos where PHPS code runs. It imports
PHPS from a read-only checkout, refuses to run when PHPS inputs (`phps/src`, `phps/cases`,
`phps/julia`, `study/pf_reference`) have uncommitted changes, and records the PHPS commit
and any other dirty paths.

```
# once: the generator venv, pinned in generate/requirements.txt (git-ignored);
# scripts/setup.ps1 does this
py -3.13 -m venv parity/generate/.venv
parity/generate/.venv/Scripts/python.exe -m pip install -r parity/generate/requirements.txt

# generate into a new directory, then bind it
parity/generate/.venv/Scripts/python.exe parity/generate/generate_pack.py \
    --phps <PHPS_Opt> --out parity/pack --reason "..."
julia --project scripts/bind_parity_pack.jl parity/pack v<N>
```

The bind script writes `parity/dist/parity_pack-v<N>.tar.gz` and binds it with the URL of
the GitHub release `parity-pack-v<N>`. Upload the tarball to that release so CI and fresh
clones can fetch it. To test against a local pack directory before binding, set
`PORTHOS_PARITY_PACK=<dir>`.

## Layout (format 1)

| File | Contents | Gate |
|---|---|---|
| `MANIFEST.json` | pack format, reason, PHPS commit and dirty paths, tool versions, seeds, SHA-256 per file | P0 |
| `cases.json` | parity cases: name, system file, bus-16 fault scenario (paths relative to `phps/cases`, same under `cases/`) | all |
| `network/<case>.json` | bus order; power-flow Y-bus (`YBusBuilder.build(include_loads=False)`); DAE Y-bus (`DiracCompiler.Y_full`, after the scenario's events are injected as `DiracRunner` does); per-bus load `G`, `B`, `P`, `Q`, `kpf`, `kqf`; Norton stamps (`ra`, `xd''` as PHPS used them); bus-fault shunts; state counts | P2 |
| `powerflow/<case>.json` | bus types, `P`/`Q` specifications, start point, solved `V`/`theta` (`PowerFlow.solve`, tol 1e-6, 20 iterations), iterations and final mismatch, computed `P`/`Q`, the case `v0`/`a0` | P3 |
| `records/certificates/*.json` | the certificate records from `study/pf_reference/certificates/`, byte-for-byte from the commit (the `.npz` arrays are not included yet) | P11 |
| `dae/<case>.json` | per P5 case (`base`, `gfl`, `vsm`, `droop`, `voc`; module `generate/dae.py`): PHPS's compiled C++ residual (`dae_residual`, the function IDA and BDF1 solve) at the initialised equilibrium `(x*, V*)` and at 50 random states, fault off and fault on (`f = -res[:n_diff]`, `g = res[n_diff:]` with `ydot = 0`); the state layout (component order, offsets, names, `delta_COI`), the resolved wiring as C++ expressions, the simulation Y-bus (`Y_full`, loads at the power-flow voltages), load arrays, slack buses, COI members and weights, `omega_b`, the fault shunts, every component's initialised parameters, and the SHA-256 of the generated kernel | P5, P6 |
| `components/<TYPE>.json` | per model type (generator module `generate/components.py`): state, input and output names; every instance's parameters as PHPS uses them (`params_used`), the ones its initialisation added or changed (`init_set`), and its equilibrium `x*`, `u*`; the branch sites of the `step` and `out` kernels (condition text, parameter-only or state-dependent) and their coverage; 200 samples with `x`, `u`, `dxdt`, the outputs after the `out` kernel and after the `step` kernel, `H`, `grad H` and the outcome of every branch site | P4 |

Complex matrices are stored as row-major `{"re": [[...]], "im": [[...]]}`. Floats are
written with Python's `repr`, so Julia reads the exact binary64 values.

Parity cases: `base` (IEEE39Bus_PF/system_phtrue), `gfl`, `gfl_zif`, `vsm`, `droop`, `voc`,
each with its bus-16 fault scenario.

## Baselines

| Version | Artifact tree hash | PHPS commit | Sections | Reason |
|---|---|---|---|---|
| v1 | `6e7d6bf617fa0d8c226616c740339d574a37c7a4` | `ba11ea1827c344ace9dd2734215b5f5b4dbbe822` (dirty: docs only, `PHPSjl_ROADMAP.md` deleted, `phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` modified) | network, powerflow, records | Initial baseline for P0 to P3. Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0, Windows 11. Superseded by v2 before it was uploaded. |
| v2 | `035b093fb4306b9d186be6b48ae82e67f357b960` | same as v1 | network, powerflow, records, components | P4: adds the `components` section for the synchronous-machine set (GENROU_PHTRUE, GENSAL_PHTRUE, IEEET1_PHTRUE, IEEEG1_PHTRUE, IEEEG3_PHTRUE, COMPLEXLOAD, sampled on the base case). The 24 files of v1 are byte-identical in v2. Same tool versions. Superseded by v3 before it was uploaded. |
| v3 | `c7353f05cb079034198de322f9dcfae7f9feaf68` | same as v1 | network, powerflow, records, components, dae | P5: adds the `dae` section. The C++ kernels were compiled with g++ 16.2.0 (MSYS2 UCRT64, `-O3`, as PHPS does). The 30 files of v2 are byte-identical in v3. |

How the `dae` samples are made: `DiracRunner.build(solver="bdf1")` initialises the case
and writes PHPS's production C++ kernel; the kernel is compiled unchanged in one
translation unit with a small harness (PHPS's `main` renamed) that reads points from a
binary file, sets the fault flags, calls `dae_residual` twice (the output pass reads some
outputs of later components from the previous call) and writes the residual back; the
harness echoes every point and the generator checks it was read bit for bit. `gfl_zif` is
not a P5 case: PHPS at ba11ea1 cannot initialise any of its three system files (the coupled
network solve does not converge at t = 0).

How the `components` samples are made: PHPS initialises the case in pure Python
(`DiracRunner.build` with a Python solver stops before C++ generation), each component's C++
kernels are translated to Python by PHPS's own `py_codegen`, and every condition is wrapped
so its outcome is recorded. Each sample is checked to be bit-identical to PHPS's unmodified
`make_step_func` / `make_out_func`. States and inputs are drawn around the equilibrium with
wide draws on limited quantities; a sample whose branch outcomes change under a 1e-9
relative perturbation is redrawn; the generator fails if any state-dependent site is not
exercised both ways.

Planned sections for later baselines: `components` for the converters (P4 items 4 and 5),
`init` (P6: the initialisation chain's intermediate values, if the `dae` equilibrium is not
enough), `sim` (P7: bus-16 fault CSVs for IDA production, IDA 1e-10 and BDF1), `ph` (P10),
`studies` (P9).
