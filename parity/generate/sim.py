"""Parity-pack section `sim` (roadmap P7).

PHPS's compiled simulations of the bus-16 fault scenario of a parity case:

  bdf1       fixed-step BDF1, dt = 5e-4, over the scenario duration
  ida_prod   IDA at PHPS's production tolerances (scenario settings), over the duration
  ida_tight  IDA at rtol = atol = 1e-10, over 15 s

For each run: `<run>.bin` holds the state and bus-voltage columns (the columns Porthos's
state vector has, delta_COI excepted) on a 5 ms grid, little-endian float64, row-major
(time by column); `<run>.json` holds the column names, the times, 20 complete CSV rows
(all columns, for the observables), PHPS's solver log figures (Newton failures of BDF1 with
their times and equations, IDA step statistics) and the run settings. The runs start from
PHPS's initialised equilibrium, which the `dae` section records.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import re
import tempfile
import time
from pathlib import Path

import numpy as np

GRID = 0.005
N_FULL_ROWS = 20


def _run(phps_root: Path, system_rel: str, events, solver: str, dt: float, duration: float,
         log_dt: float, rtol=None, atol=None) -> tuple[Path, str, float]:
    from src.dirac.dae_runner import DiracRunner

    out = Path(tempfile.mkdtemp())
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        r = DiracRunner(str(phps_root / "phps" / "cases" / system_rel), output_dir=str(out),
                        events=events)
        r.build(dt=dt, duration=duration, solver=solver, rtol=rtol, atol=atol, log_dt=log_dt)
        t0 = time.time()
        r.run()
        seconds = time.time() - t0
    return out, log.getvalue(), seconds


def _record(out: Path, text: str, seconds: float, settings: dict, dest: Path, label: str):
    with open(out / "simulation_results.csv", newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    data = np.array(body, dtype=float)
    t = data[:, 0]
    keep_cols = [j for j, h in enumerate(header)
                 if (h.startswith("Vd_Bus") or h.startswith("Vq_Bus")
                     or (h.count(".") == 1 and j > 0))]
    # states only among the dotted columns: the component state names come from the kernel;
    # observables are dotted too, so keep every dotted column and let Porthos pick by name
    k = np.round(t / GRID)
    grid = np.where(np.abs(t - k * GRID) <= 1e-9)[0]
    block = data[np.ix_(grid, keep_cols)]
    dest.mkdir(parents=True, exist_ok=True)
    block.astype("<f8").tofile(dest / f"{label}.bin")
    full = np.linspace(0, len(body) - 1, N_FULL_ROWS).astype(int)
    meta = {
        "run": label,
        "settings": settings,
        "run_seconds": seconds,
        "columns": [header[j] for j in keep_cols],
        "times": [float(v) for v in t[grid]],
        "shape": [len(grid), len(keep_cols)],
        "csv_header": header,
        "full_rows": [[float(v) for v in data[i]] for i in full],
        "n_records": len(body),
    }
    m = re.search(r"Newton did not reach tol on (\d+) step", text)
    meta["nonconverged_steps"] = int(m.group(1)) if m else 0
    bad = out / "nonconverged_steps.csv"
    if bad.exists():
        with open(bad, newline="") as fh:
            meta["nonconverged"] = [{"t": float(r["t"]), "row": r["row_name"],
                                     "res_norm": float(r["res_norm"])}
                                    for r in csv.DictReader(fh)]
    m = re.search(r"Steps: (\d+)\s+Residual evals: (\d+)\s+Jacobian evals: (\d+)", text)
    if m:
        meta["ida_steps"], meta["ida_res_evals"], meta["ida_jac_evals"] = map(int, m.groups())
    return meta


def section_sim(phps_root: Path, name: str, system_rel: str, scenario_rel: str, out: Path):
    scenario = json.loads((phps_root / "phps" / "cases" / scenario_rel).read_text(encoding="utf-8"))
    events = scenario.get("events", [])
    duration = float(scenario["solver"]["duration"])
    runs = [
        ("bdf1", dict(solver="bdf1", dt=5e-4, duration=duration, log_dt=1e-3)),
        ("ida_prod", dict(solver="ida", dt=5e-4, duration=duration, log_dt=1e-3,
                          rtol=scenario["solver"].get("rtol"),
                          atol=scenario["solver"].get("atol"))),
        ("ida_tight", dict(solver="ida", dt=5e-4, duration=15.0, log_dt=1e-3,
                           rtol=1e-10, atol=1e-10)),
    ]
    index = {}
    for label, s in runs:
        run_out, text, seconds = _run(phps_root, system_rel, events, **s)
        index[label] = _record(run_out, text, seconds, s, out / "sim" / name, label)
    for label, meta in index.items():
        (out / "sim" / name / f"{label}.json").write_text(
            json.dumps(meta, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return list(index)
