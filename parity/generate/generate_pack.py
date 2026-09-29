"""Generate the Porthos parity pack from a PHPS checkout.

This is the only place in Porthos where PHPS code runs (roadmap section 1). It imports
PHPS from a read-only checkout, evaluates the reference quantities for a fixed set of
cases, and writes them as JSON with a MANIFEST.json that records the PHPS commit, the
tool versions and a SHA-256 per file.

The pack is built in sections, so later phases can extend it:

    network     Y-bus (power-flow and DAE variants), load admittances, Norton stamps,
                fault shunts                                           (P2)
    powerflow   Newton-Raphson bus voltages and specifications         (P3)
    records     certificate and PH-audit records copied unchanged      (P10, P11 targets)
    components  model kernels, outputs, H and grad H at random states and inputs, with
                branch outcomes (see components.py)                    (P4)
    dae         PHPS's compiled DAE residual f, g at the equilibrium and at random states,
                fault off and on, with the system layout (see dae.py)  (P5, P6)
    sim         PHPS's compiled BDF1 and IDA runs of the bus-16 fault (see sim.py)  (P7)

A new or regenerated pack is a new baseline: it gets a new artifact hash and a written
reason in parity/README.md.

Usage (from the Porthos.jl root, with the venv made from PHPS's requirements.txt):

    parity/generate/.venv/Scripts/python.exe parity/generate/generate_pack.py \
        --phps C:/Users/em18736/Documents/PHPS_Opt --out parity/pack

Floats are written with Python's repr (shortest round-trip), so Julia reads the exact
binary64 values PHPS computed.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import hashlib
import io
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import types
from pathlib import Path

PACK_FORMAT = 1

# Parity cases: (name, system file, bus-16 fault scenario), paths relative to phps/cases.
PARITY_CASES = [
    ("base", "IEEE39Bus_PF/system_phtrue.json", "IEEE39Bus_PF/bus_fault_bus16_150ms.json"),
    ("gfl", "IEEE39Bus_PF_gfl/system_gfl.json", "IEEE39Bus_PF_gfl/bus_fault_gfl_bus16.json"),
    ("gfl_zif", "IEEE39Bus_PF_gfl-zif/system_gfl_zif.json",
     "IEEE39Bus_PF_gfl-zif/bus_fault_gfl_zif_bus16_40ms_dt250us_bolted.json"),
    ("vsm", "IEEE39Bus_PF_gfm-vsm/system_gfm_vsm.json",
     "IEEE39Bus_PF_gfm-vsm/bus_fault_gfm_vsm_bus16.json"),
    ("droop", "IEEE39Bus_PF_gfm-droop/system_gfm_droop.json",
     "IEEE39Bus_PF_gfm-droop/bus_fault_gfm_droop_bus16.json"),
    ("voc", "IEEE39Bus_PF_gfm-voc/system_gfm_voc.json",
     "IEEE39Bus_PF_gfm-voc/bus_fault_gfm_voc_bus16.json"),
]

# Cases of the P5 to P7 gates (roadmap section 1). gfl_zif is not among them: PHPS at
# ba11ea1 cannot initialise it (coupled network solve does not converge at t = 0).
DAE_CASES = ("base", "gfl", "vsm", "droop", "voc")

# Cases with reference simulations: the base case and the grid-forming converter cases (gfl
# follows with its model, after its reservoir rework).
SIM_CASES = ("base", "vsm", "droop", "voc")

SECTIONS = ("network", "powerflow", "records", "components", "dae", "sim")

# PHPS paths whose modification would change the reference numbers. The generator refuses
# to run if any of them is dirty; other dirty paths (documentation) are recorded.
PHPS_INPUT_PREFIXES = ("phps/src/", "phps/cases/", "phps/julia/", "study/pf_reference/")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _git(phps: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(phps), *args], check=True,
                          capture_output=True, text=True).stdout


def phps_provenance(phps: Path) -> dict:
    commit = _git(phps, "rev-parse", "HEAD").strip()
    status = _git(phps, "status", "--porcelain", "--untracked-files=all")
    dirty = []
    for line in status.splitlines():
        path = line[3:].strip().strip('"')
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        dirty.append({"status": line[:2].strip(), "path": path})
    blocking = [d["path"] for d in dirty
                if d["path"].startswith(PHPS_INPUT_PREFIXES) and "__pycache__" not in d["path"]]
    if blocking:
        raise SystemExit("PHPS inputs are modified; commit or stash them first:\n  "
                         + "\n  ".join(blocking))
    remote = _git(phps, "config", "--get", "remote.origin.url").strip()
    return {"commit": commit, "remote": remote, "dirty_paths": dirty}


def _check_finite(obj, where="") -> None:
    if isinstance(obj, float):
        if not math.isfinite(obj):
            raise ValueError(f"non-finite value at {where}")
    elif isinstance(obj, dict):
        for k, v in obj.items():
            _check_finite(v, f"{where}.{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            _check_finite(v, f"{where}[{i}]")


def write_json(path: Path, obj) -> None:
    _check_finite(obj, path.name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(obj, fh, indent=1, allow_nan=False)
        fh.write("\n")


def fl(x) -> float:
    return float(x)


def vec(a) -> list:
    return [float(v) for v in a]


def cmat(M) -> dict:
    """Dense complex matrix as row-major real and imaginary parts."""
    return {"re": [[float(v.real) for v in row] for row in M],
            "im": [[float(v.imag) for v in row] for row in M]}


@contextlib.contextmanager
def quiet():
    """Silence PHPS's progress prints (they go to stdout)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        yield buf


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def section_network(phps_root: Path, name: str, system_rel: str, scenario_rel: str) -> dict:
    from src.ybus import YBusBuilder
    from src.dirac.dae_compiler import DiracCompiler
    from src.dirac.dae_runner import DiracRunner

    cases = phps_root / "phps" / "cases"
    scenario = json.loads((cases / scenario_rel).read_text(encoding="utf-8"))
    system_path = cases / system_rel

    with quiet():
        # Power-flow Y-bus: lines and shunts only (PowerFlow uses build(include_loads=False)).
        dc = DiracCompiler(str(system_path))
        y_pf = YBusBuilder(dc.base_compiler.data).build(include_loads=False).copy()
        # Inject the scenario events exactly as the DAE runner does, then build the DAE
        # structure (full Y-bus with constant-Z loads and Norton stamps, fault shunts).
        runner_view = types.SimpleNamespace(base_compiler=dc.base_compiler)
        DiracRunner._inject_events(runner_view, scenario.get("events", []))
        dc.build()

    norton = []
    for comp in dc.components:
        if comp.component_role != "generator" or not comp.contributes_norton_admittance:
            continue
        norton.append({
            "name": comp.name,
            "type": type(comp).__name__,
            "bus": comp.params.get("bus"),
            "ra": fl(comp.params.get("ra", 0.0)),
            "xd_pp": fl(comp.params.get("xd_double_prime", comp.params.get("xd1", 0.2))),
        })

    faults = [{
        "bus": int(f["bus_id"]), "bus_index": int(f["bus_idx"]) + 1,
        "t_start": fl(f["t_start"]), "t_end": fl(f["t_end"]),
        "g": fl(f["g"]), "b": fl(f["b"]),
    } for f in dc.fault_events_dae]

    return {
        "case": name,
        "system": system_rel,
        "scenario": scenario_rel,
        "bus_ids": [int(b) for b in dc.bus_indices],
        "Y_pf": cmat(y_pf),
        "Y_dae": cmat(dc.Y_full),
        "load_G": vec(dc.load_G), "load_B": vec(dc.load_B),
        "load_P": vec(dc.load_P), "load_Q": vec(dc.load_Q),
        "load_kpf": vec(dc.load_kpf), "load_kqf": vec(dc.load_kqf),
        "norton": norton,
        "excluded_dyn_lines": [ln.get("idx") for ln in getattr(dc, "_excluded_dyn_lines", [])],
        "faults": faults,
        "n_diff": int(dc.n_diff), "n_alg": int(dc.n_alg),
    }


def section_powerflow(phps_root: Path, name: str, system_rel: str) -> dict:
    from src.system_graph import build_system_graph
    from src.powerflow import PowerFlow

    system_path = phps_root / "phps" / "cases" / system_rel
    with quiet():
        graph = build_system_graph(str(system_path))
    data = graph.raw_data
    skip = bool(data.get("config", {}).get("skip_pf_solve", False))

    tol, max_iter = 1e-6, 20          # PowerFlow.solve defaults, as Initializer calls it
    pf = PowerFlow(data)
    V_start, th_start = pf.V.copy(), pf.theta.copy()
    with quiet() as buf:
        if skip:
            converged = bool(pf.load_bus_overrides())
        else:
            converged = bool(pf.solve(tol=tol, max_iter=max_iter))
    log = buf.getvalue()
    m = re.search(r"Converged in (\d+) iterations\. Norm: (\S+)", log)
    S = pf.calculate_power()

    buses = sorted(data["Bus"], key=lambda b: b["idx"])
    return {
        "case": name,
        "system": system_rel,
        "skip_pf_solve": skip,
        "tol": tol, "max_iter": max_iter,
        "converged": converged,
        "iterations": int(m.group(1)) if m else None,
        "final_mismatch": float(m.group(2)) if m else None,
        "bus_ids": [int(b) for b in pf.buses],
        "bus_types": [int(t) for t in pf.bus_types],
        "P_spec": vec(pf.P_spec), "Q_spec": vec(pf.Q_spec),
        "V_start": vec(V_start), "theta_start": vec(th_start),
        "V": vec(pf.V), "theta": vec(pf.theta),
        "P_calc": vec(S.real), "Q_calc": vec(S.imag),
        "case_v0": [fl(b.get("v0", 1.0)) for b in buses],
        "case_a0": [fl(b.get("a0", 0.0)) for b in buses],
    }


# PH-audit records of the component models review (the P10 targets), from
# study/pf_reference/model_review/.
MODEL_REVIEW_RECORDS = ("audit_controller_kyp.json", "audit_governor_nonpassivity_exact.json")


def section_records(phps_root: Path, out: Path) -> list:
    """Copy records unchanged: every certificate record (the P11 targets, with the shifted-
    storage audits of P10) and the model-review PH audits in MODEL_REVIEW_RECORDS (P10)."""
    base = phps_root / "study" / "pf_reference"
    sources = [("certificates", src) for src in sorted((base / "certificates").glob("*.json"))]
    sources += [("model_review", base / "model_review" / name) for name in MODEL_REVIEW_RECORDS]
    copied = []
    for sub, src in sources:
        dst = out / "records" / sub / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Take the committed blob, so the pack matches the recorded commit byte for byte.
        rel = src.relative_to(phps_root).as_posix()
        blob = subprocess.run(["git", "-C", str(phps_root), "show", f"HEAD:{rel}"],
                              check=True, capture_output=True).stdout
        dst.write_bytes(blob)
        copied.append({"file": f"records/{sub}/{src.name}", "phps_path": rel})
    return copied


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--phps", required=True, type=Path, help="PHPS_Opt checkout (read only)")
    ap.add_argument("--out", required=True, type=Path, help="output directory (must not exist)")
    ap.add_argument("--sections", default=",".join(SECTIONS))
    ap.add_argument("--reason", default="initial baseline",
                    help="why this baseline exists (recorded in the manifest)")
    args = ap.parse_args()

    phps_root = args.phps.resolve()
    out = args.out.resolve()
    sections = [s for s in args.sections.split(",") if s]
    for s in sections:
        if s not in SECTIONS:
            raise SystemExit(f"unknown section {s!r}; known: {SECTIONS}")
    if out.exists():
        raise SystemExit(f"{out} exists; a pack is never regenerated in place")

    provenance = phps_provenance(phps_root)
    # PHPS resolves its own paths relative to phps/ (imports as `src.*`).
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(phps_root / "phps"))
    os.chdir(phps_root / "phps")

    out.mkdir(parents=True)
    cases_index = []
    for name, system_rel, scenario_rel in PARITY_CASES:
        entry = {"name": name, "system": system_rel, "scenario": scenario_rel}
        if "network" in sections:
            write_json(out / "network" / f"{name}.json",
                       section_network(phps_root, name, system_rel, scenario_rel))
        if "powerflow" in sections:
            write_json(out / "powerflow" / f"{name}.json",
                       section_powerflow(phps_root, name, system_rel))
        cases_index.append(entry)
        print(f"  {name}: done", file=sys.stderr)
    write_json(out / "cases.json", cases_index)

    records = section_records(phps_root, out) if "records" in sections else []

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    if "dae" in sections:
        from dae import section_dae
        for k, (name, system_rel, scenario_rel) in enumerate(PARITY_CASES):
            if name not in DAE_CASES:
                continue
            write_json(out / "dae" / f"{name}.json",
                       section_dae(phps_root, name, system_rel, scenario_rel, 20260930 + k))
            print(f"  dae: {name} done", file=sys.stderr)

    if "sim" in sections:
        from sim import section_sim
        for name, system_rel, scenario_rel in PARITY_CASES:
            if name not in SIM_CASES:
                continue
            section_sim(phps_root, name, system_rel, scenario_rel, out)
            print(f"  sim: {name} done", file=sys.stderr)

    if "components" in sections:
        from components import section_components
        case_by_name = {n: s for n, s, _ in PARITY_CASES}
        for ctype, rec in section_components(phps_root, case_by_name).items():
            write_json(out / "components" / f"{ctype}.json", rec)
            print(f"  components: {ctype} done", file=sys.stderr)

    import numpy, scipy, sympy
    files = {}
    for p in sorted(out.rglob("*")):
        if p.is_file():
            files[p.relative_to(out).as_posix()] = sha256_file(p)
    manifest = {
        "pack_format": PACK_FORMAT,
        "created_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "reason": args.reason,
        "phps": provenance,
        "generator": {
            "script": "parity/generate/generate_pack.py",
            "sha256": sha256_file(Path(__file__)),
        },
        "tools": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": numpy.__version__, "scipy": scipy.__version__, "sympy": sympy.__version__,
        },
        "seeds": {},
        "sections": sections,
        "cases": cases_index,
        "records": records,
        "files": files,
    }
    write_json(out / "MANIFEST.json", manifest)
    print(f"parity pack written to {out} ({len(files)} files)", file=sys.stderr)


if __name__ == "__main__":
    main()
