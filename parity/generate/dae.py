"""Parity-pack section `dae` (roadmap P5, and the P6 equilibrium).

For each parity case with its bus-16 fault scenario, PHPS initialises the system and
generates its production C++ DAE kernel (`DiracRunner.build`). The kernel's
`dae_residual(y, ydot, res, t)` is then compiled UNCHANGED together with a small harness
(appended to the same translation unit, with PHPS's `main` renamed) and evaluated at:

  - the initialised equilibrium (x*, V*), which is also the P6 reference;
  - 50 random states and bus voltages around it, fault off and fault on.

With ydot = 0 the residual is res = [-f(x, V); g(x, V)], so f = -res[:n_diff] and
g = res[n_diff:]. The section also records what Porthos needs to build the same system
before its own initialisation exists (P6): the state layout, the resolved wiring, the
network constants exactly as the C++ has them, and every component's initialised
parameters.

Compiler: g++ from MSYS2 UCRT64 (the location PHPS itself searches on Windows), with PHPS's
flags (-O3).
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

N_SAMPLES = 50
UCRT64 = Path(r"C:\msys64\ucrt64")

HARNESS = r"""
#undef main
#include <cstdio>
// Porthos parity harness: evaluate PHPS's dae_residual at points read from a binary file
// (int32 count, then per point: int32 fault flag, float64 t, N_TOTAL float64 y).
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    FILE* fi = fopen(argv[1], "rb");
    FILE* fo = fopen(argv[2], "wb");
    if (!fi || !fo) return 3;
    int npts = 0;
    if (fread(&npts, sizeof(int), 1, fi) != 1) return 4;
    static double y[N_TOTAL], yd[N_TOTAL], res[N_TOTAL];
    for (int k = 0; k < npts; ++k) {
        int flag = 0; double t = 0.0;
        if (fread(&flag, sizeof(int), 1, fi) != 1) return 5;
        if (fread(&t, sizeof(double), 1, fi) != 1) return 5;
        if (fread(y, sizeof(double), N_TOTAL, fi) != (size_t)N_TOTAL) return 6;
        for (int f = 0; f < N_FAULTS; ++f) fault_active[f] = flag;
        for (int i = 0; i < N_TOTAL; ++i) yd[i] = 0.0;
        // twice: the output pass reads some outputs of components later in the order
        // from the previous call; the second call is the settled value
        dae_residual(y, yd, res, t);
        dae_residual(y, yd, res, t);
        fwrite(y, sizeof(double), N_TOTAL, fo);
        fwrite(res, sizeof(double), N_TOTAL, fo);
    }
    fclose(fi);
    fclose(fo);
    return 0;
}
"""


def _compile(source: Path, binary: Path) -> None:
    env = dict(os.environ)
    env["PATH"] = str(UCRT64 / "bin") + os.pathsep + r"C:\msys64\usr\bin" + os.pathsep + env["PATH"]
    cmd = [str(UCRT64 / "bin" / "g++.exe"), "-O3", f"-I{UCRT64 / 'include'}", str(source),
           "-o", str(binary), f"-L{UCRT64 / 'lib'}", "-lm"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError("harness compilation failed:\n" + r.stderr[-4000:])


def _run_residual(binary: Path, points, flags) -> list:
    tmp = Path(tempfile.mkdtemp())
    fin, fout = tmp / "in.bin", tmp / "out.bin"
    n = len(points[0])
    with open(fin, "wb") as fh:
        fh.write(np.int32(len(points)).tobytes())
        for y, flag in zip(points, flags):
            fh.write(np.int32(flag).tobytes())
            fh.write(np.float64(0.0).tobytes())
            fh.write(np.asarray(y, dtype="<f8").tobytes())
    env = dict(os.environ)
    env["PATH"] = str(UCRT64 / "bin") + os.pathsep + env["PATH"]
    r = subprocess.run([str(binary), str(fin), str(fout)], capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"harness failed ({r.returncode}): {r.stderr}")
    vals = np.fromfile(fout, dtype="<f8").reshape(len(points), 2, n)
    out = []
    for k, y in enumerate(points):
        if not np.array_equal(vals[k, 0].view(np.int64), np.asarray(y, float).view(np.int64)):
            raise RuntimeError("harness did not read the points back exactly")
        out.append(vals[k, 1].tolist())
    return out


def _num_params(params: dict) -> dict:
    return {k: float(v) for k, v in params.items()
            if isinstance(v, (int, float)) and not isinstance(v, bool)
            and math.isfinite(float(v))}


def section_dae(phps_root: Path, name: str, system_rel: str, scenario_rel: str,
                seed: int) -> dict:
    from src.dirac.dae_runner import DiracRunner
    from src.dirac.frame import coi_weight

    cases = phps_root / "phps" / "cases"
    scenario = json.loads((cases / scenario_rel).read_text(encoding="utf-8"))
    solver = scenario.get("solver", {})
    out_dir = Path(tempfile.mkdtemp())
    with contextlib.redirect_stdout(io.StringIO()):
        runner = DiracRunner(str(cases / system_rel), output_dir=str(out_dir),
                             events=scenario.get("events", []))
        runner.build(dt=float(solver.get("dt", 0.0005)),
                     duration=float(solver.get("duration", 10.0)), solver="bdf1",
                     log_dt=solver.get("log_dt"))
    dc = runner.dae_compiler
    bc = runner.base_compiler

    # the production kernel, unchanged, plus the harness
    src = Path(runner.source_path).read_text(encoding="utf-8")
    harness_src = out_dir / "phps_residual_harness.cpp"
    harness_src.write_text("#define main phps_embedded_main\n" + src + HARNESS,
                           encoding="utf-8")
    harness_bin = out_dir / "phps_residual_harness.exe"
    _compile(harness_src, harness_bin)

    n_diff, n_bus = dc.n_diff, dc.n_bus
    x_star = np.array(runner.x0, dtype=float)
    Vd_star = np.array(runner.Vd_init, dtype=float)
    Vq_star = np.array(runner.Vq_init, dtype=float)
    y_star = np.concatenate([x_star, np.column_stack([Vd_star, Vq_star]).reshape(-1)])

    rng = np.random.default_rng(seed)
    points = []
    for _ in range(N_SAMPLES):
        x = x_star * (1.0 + rng.uniform(-0.05, 0.05, x_star.size)) \
            + rng.uniform(-0.02, 0.02, x_star.size)
        Vd = Vd_star + rng.uniform(-0.05, 0.05, n_bus)
        Vq = Vq_star + rng.uniform(-0.05, 0.05, n_bus)
        points.append(np.concatenate([x, np.column_stack([Vd, Vq]).reshape(-1)]))

    has_fault = len(dc.fault_events_dae) > 0
    all_pts = [y_star] + points + ([y_star] + points if has_fault else [])
    flags = [0] * (1 + N_SAMPLES) + ([1] * (1 + N_SAMPLES) if has_fault else [])
    res = _run_residual(harness_bin, all_pts, flags)

    def split(r):
        r = np.asarray(r)
        return [float(-v) for v in r[:n_diff]], [float(v) for v in r[n_diff:]]

    samples = []
    for k, y in enumerate(points):
        f, g = split(res[1 + k])
        s = {"x": [float(v) for v in y[:n_diff]], "V": [float(v) for v in y[n_diff:]],
             "f": f, "g": g}
        if has_fault:
            ff, gf = split(res[2 + N_SAMPLES + k])
            s["f_fault"], s["g_fault"] = ff, gf
        samples.append(s)
    f0, g0 = split(res[0])
    equilibrium = {"x": [float(v) for v in x_star],
                   "Vd": [float(v) for v in Vd_star], "Vq": [float(v) for v in Vq_star],
                   "f": f0, "g": g0,
                   "x_initial": [float(v) for v in runner.x0_initial]}
    if has_fault:
        equilibrium["f_fault"], equilibrium["g_fault"] = split(res[1 + N_SAMPLES])

    gen_bus_ids = set(dc.gen_bus_map.values())
    slack = []
    for bus_id, V_ref in dc.slack_V.items():
        if bus_id in dc.bus_map and bus_id not in gen_bus_ids:
            slack.append({"bus": int(bus_id), "index": int(dc.bus_map[bus_id]) + 1,
                          "Vd": float(V_ref.real), "Vq": float(V_ref.imag)})
    coi = []
    for comp in dc.components:
        if comp.component_role == "generator" and "omega" in comp.state_schema:
            coi.append({"component": comp.name, "weight": float(coi_weight(comp))})
    wiring = {}
    for comp in dc.components:
        for pname, _, _ in comp.port_schema["in"]:
            wiring[f"{comp.name}.{pname}"] = dc.wiring_map.get((comp.name, pname))

    return {
        "case": name, "system": system_rel, "scenario": scenario_rel, "seed": seed,
        "n_diff": int(n_diff), "n_bus": int(n_bus), "n_alg": int(dc.n_alg),
        "delta_coi_index": int(dc.delta_coi_idx) + 1,
        "bus_ids": [int(b) for b in dc.bus_indices],
        "components": [{"name": c.name, "type": bc.graph.raw_data["components"][c.name]["type"],
                        "offset": int(dc.state_offsets[c.name]) + 1,
                        "states": list(c.state_schema),
                        "inputs": [p[0] for p in c.port_schema["in"]],
                        "outputs": [p[0] for p in c.port_schema["out"]],
                        "role": c.component_role,
                        "params": _num_params(c.params)}
                       for c in dc.components],
        "wiring": wiring,
        "Y_full": {"re": [[float(v) for v in row] for row in dc.Y_full.real],
                   "im": [[float(v) for v in row] for row in dc.Y_full.imag]},
        "load_G": [float(v) for v in dc.load_G], "load_B": [float(v) for v in dc.load_B],
        "load_P": [float(v) for v in dc.load_P], "load_Q": [float(v) for v in dc.load_Q],
        "load_kpf": [float(v) for v in dc.load_kpf], "load_kqf": [float(v) for v in dc.load_kqf],
        "slack": slack,
        "coi": coi,
        # the COI base frequency: the first swing source's omega_b, as the C++ emits it
        "omega_b_sys": next((str(c.params.get("omega_b", "2.0 * M_PI * 60.0"))
                             for c in dc.components
                             if c.component_role == "generator" and "omega" in c.state_schema),
                            None),
        "faults": [{"bus": int(f["bus_id"]), "index": int(f["bus_idx"]) + 1,
                    "t_start": float(f["t_start"]), "t_end": float(f["t_end"]),
                    "g": float(f["g"]), "b": float(f["b"])} for f in dc.fault_events_dae],
        "equilibrium": equilibrium,
        "samples": samples,
        "kernel_sha256": hashlib.sha256(src.encode("utf-8")).hexdigest(),
    }
