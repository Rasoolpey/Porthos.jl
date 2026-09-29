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
# once: a venv with PHPS's requirements (git-ignored)
py -3.13 -m venv parity/generate/.venv
parity/generate/.venv/Scripts/python.exe -m pip install -r <PHPS_Opt>/phps/requirements.txt

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

Complex matrices are stored as row-major `{"re": [[...]], "im": [[...]]}`. Floats are
written with Python's `repr`, so Julia reads the exact binary64 values.

Parity cases: `base` (IEEE39Bus_PF/system_phtrue), `gfl`, `gfl_zif`, `vsm`, `droop`, `voc`,
each with its bus-16 fault scenario.

## Baselines

| Version | Artifact tree hash | PHPS commit | Sections | Reason |
|---|---|---|---|---|
| v1 | `6e7d6bf617fa0d8c226616c740339d574a37c7a4` | `ba11ea1827c344ace9dd2734215b5f5b4dbbe822` (dirty: docs only, `PHPSjl_ROADMAP.md` deleted, `phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` modified) | network, powerflow, records | Initial baseline for P0 to P3. Python 3.13.15, NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0, Windows 11. |

Planned sections for later baselines: `components` (P4: `f`, outputs, injection, `H`,
`grad H` at 200 random states with limiter sides), `dae` (P5: `f`, `g` at 50 states, state
names), `init` (P6: `x*`, `V*`, residual), `sim` (P7: bus-16 fault CSVs for IDA production,
IDA 1e-10 and BDF1), `ph` (P10), `studies` (P9).
