"""Parity-pack section `components` (roadmap P4).

For each component type, evaluate PHPS's own model code at random states and inputs:
the state derivatives, the outputs (the `out` kernel, then the `step` kernel, in the order
PHPS's DAE evaluates them), the Hamiltonian and its gradient, and which way every branch
site in the kernels went.

The model code is PHPS's: each kernel is the component's C++ (`get_cpp_step_code`,
`get_cpp_compute_outputs_code`) with its parameters prepended, translated to Python by
PHPS's own `src.dirac.py_codegen` (the translation PHPS's Python solvers run). The only
change made here is that every condition is wrapped in `_br(site, cond)`, which records
the outcome and returns `cond` unchanged. Every sample is checked to give bit-identical
results to PHPS's unmodified `make_step_func` / `make_out_func`.

Components are taken after PHPS's initialisation of the case (`DiracRunner.build` with a
Python solver, which stops before C++ generation), because initialisation sets parameters
the kernels use (reservoir references, load V0 and Vini from the power flow, governor Pref,
machine Efd0, ...). Every parameter that initialisation added or changed, found by
comparing with a snapshot taken before it ran, is recorded as `init_set`.

Sampling: states and inputs are drawn around the initialised equilibrium, with wide draws
on limited quantities so both sides of every limiter occur. A sample whose branch outcomes
change under a relative perturbation of 1e-9 is within rounding distance of a switching
surface and is redrawn.

Conditions inside a loop (the converters' 60-step bisection for the virtual-impedance
current) belong to a numerical root solve: their late outcomes sit on the root by
construction and flip under any perturbation. They are recorded and marked `in_loop`, and
left out of the redraw test; Porthos does not match them (the solve's result is compared
through the outputs).

Parameter variants (`VARIANTS`): some branches cannot switch with a case's parameters (the
pf_frame reference-frame input, the VSM field controller, the droop's adaptive droop and
virtual impedance, the passivity-based VOC). A variant is a copy of the type's first case
instance with a few parameters overridden (recorded as `variant_of` and `overrides`); it is
sampled like the others, so those branches are exercised and checked too.
"""

from __future__ import annotations

import contextlib
import copy
import io
import math
import re
import tempfile

import numpy as np

N_SAMPLES = 200
PERTURB_REL = 1e-9
N_PERTURB = 4

# component type -> parity case it is sampled from (the synchronous-machine set first,
# roadmap "Order of work"; converters follow in a later baseline)
COMPONENT_SOURCES = {
    "GENROU_PHTRUE": "base",
    "GENSAL_PHTRUE": "base",
    "IEEET1_PHTRUE": "base",
    "IEEEG1_PHTRUE": "base",
    "IEEEG3_PHTRUE": "base",
    "COMPLEXLOAD": "base",
    # grid-forming converters (P4 item 5), each on its own case
    "GFM_VSM_PHTRUE": "vsm",
    "GFM_DROOP_PHTRUE": "droop",
    "GFM_VOC_PHTRUE": "voc",
}

# seed offsets: the order in which the types joined the pack (v2: the sorted synchronous-
# machine set; v6: the grid-forming converters). Append new types at the end.
SEED_ORDER = ["COMPLEXLOAD", "GENROU_PHTRUE", "GENSAL_PHTRUE", "IEEEG1_PHTRUE",
              "IEEEG3_PHTRUE", "IEEET1_PHTRUE", "GFM_DROOP_PHTRUE", "GFM_VOC_PHTRUE",
              "GFM_VSM_PHTRUE"]

# type -> [(suffix, parameter overrides)]: copies of the first instance that exercise the
# branches the case's parameters keep fixed
VARIANTS = {
    "GFM_VSM_PHTRUE": [("modes", {"pf_frame": 1.0, "K_field": 0.5, "Dq": 5.0, "n_vi": 2.0,
                                  "D_direct": 2.0})],
    "GFM_DROOP_PHTRUE": [("modes", {"pf_frame": 1.0, "adapt_droop": 1.0, "adapt_vi": 1.0,
                                    "vi_mode": 1.0, "x_vi": 0.01})],
    "GFM_VOC_PHTRUE": [("modes", {"pf_frame": 1.0, "i_ref_max": 5.0})],
}

# type -> {site text: reason}: state-dependent sites that no sample can reach, recorded as
# `unreached` (the coverage check then requires that they were never reached)
UNREACHED = {
    # pvoc_mode = 1 (Kong et al.) calls tanh, which PHPS's Python translation of the
    # kernels (src/dirac/py_codegen at ba11ea1) does not provide: PHPS can only run that
    # mode compiled, so there is no Python reference for it. No parity case uses it.
    "GFM_VOC_PHTRUE": {"vsq > 1.0e-6": "pvoc_mode = 1 only; PHPS's py_codegen has no tanh; "
                                       "unsupported in Porthos (rejected at construction)"},
}

_KEYWORDS = {"and", "or", "not", "True", "False", "math", "abs", "min", "max"}
_IDENT = re.compile(r"\b[A-Za-z_][A-Za-z_0-9]*\b")
_NUMERIC = re.compile(r"^[\s0-9eE+\-*/().]*$")


# ---------------------------------------------------------------------------
# Kernel instrumentation
# ---------------------------------------------------------------------------

def _find_matching_paren(s: str, start: int) -> int:
    depth = 0
    for i in range(start, len(s)):
        if s[i] == "(":
            depth += 1
        elif s[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError(f"unbalanced parentheses in {s!r}")


class Kernel:
    """One translated kernel (`step` or `out`) with recorded branch sites."""

    def __init__(self, comp, mode: str):
        from src.dirac.py_codegen import _build_kernel_code, translate_cpp_kernel

        self.mode = mode
        body = translate_cpp_kernel(_build_kernel_code(comp, mode))
        self.sites: list[dict] = []
        # Names bound to constants: every assignment to them (anywhere in the kernel,
        # conditional or not) involves only literals and other constants. Fixpoint over
        # the whole kernel, so `v = x[1]; if v > HI: v = HI` leaves v state-dependent.
        assigns = []
        for line in body.split("\n"):
            m = re.match(r"^\s*([A-Za-z_]\w*)\s*=\s*(?!=)(.+)$", line)
            if m:
                assigns.append((m.group(1), m.group(2)))
        assigned = {n for n, _ in assigns}
        nonconst = set()
        changed = True
        while changed:
            changed = False
            for n, rhs in assigns:
                if n in nonconst:
                    continue
                ids = {t for t in _IDENT.findall(rhs)} - _KEYWORDS
                if ("x[" in rhs or "inputs[" in rhs or "outputs[" in rhs
                        or ids & nonconst or not ids <= assigned):
                    nonconst.add(n)
                    changed = True
        const = assigned - nonconst
        lines = []
        loop_indent = None
        for line in body.split("\n"):
            indent = len(line) - len(line.lstrip())
            if loop_indent is not None and line.strip() and indent <= loop_indent:
                loop_indent = None
            self._in_loop = loop_indent is not None
            lines.append(self._instrument(line, const))
            if loop_indent is None and re.match(r"^\s*for .+:\s*$", line):
                loop_indent = indent
        self.source = "\n".join(lines)
        args = "x, dxdt, inputs, outputs, t" if mode == "step" else "x, inputs, outputs, t"
        src = "import math as math\ndef _k(%s, _br):\n" % args
        src += "\n".join("    " + ln for ln in self.source.split("\n")) + "\n"
        ns: dict = {}
        exec(compile(src, f"<{mode}_{comp.name}>", "exec"), ns)
        self._fn = ns["_k"]

    def _site(self, cond: str, const: set) -> int:
        ids = {t for t in _IDENT.findall(cond)} - _KEYWORDS
        state_dep = ("x[" in cond or "inputs[" in cond or "outputs[" in cond
                     or not ids <= const)
        site = {"id": len(self.sites), "text": cond.strip(), "param_only": not state_dep}
        if getattr(self, "_in_loop", False):
            site["in_loop"] = True
        self.sites.append(site)
        return len(self.sites) - 1

    def _instrument(self, line: str, const: set) -> str:
        m = re.match(r"^(\s*)(if|elif) (.+):\s*$", line)
        if m:
            k = self._site(m.group(3), const)
            return f"{m.group(1)}{m.group(2)} _br({k}, ({m.group(3)})):"
        # ternaries: `A if (COND) else B`
        out, i = "", 0
        while True:
            j = line.find(" if (", i)
            if j < 0:
                out += line[i:]
                break
            p = j + len(" if ")
            q = _find_matching_paren(line, p)
            cond = line[p + 1:q]
            if not line[q + 1:].lstrip().startswith("else"):
                raise ValueError(f"unexpected conditional form: {line!r}")
            k = self._site(cond, const)
            out += line[i:j] + f" if _br({k}, ({cond}))"
            i = q + 1
        return out

    def __call__(self, *args):
        rec: list = []

        def _br(k, c):
            rec.append((k, bool(c)))
            return c

        self._fn(*args, _br)
        return rec


# ---------------------------------------------------------------------------
# Sampling boxes
# ---------------------------------------------------------------------------

def _mix(rng, p_wide, local, wide):
    return wide() if rng.random() < p_wide else local()


def _reservoir(rng, center, C):
    # mostly near the equilibrium level; sometimes near empty, where the guard p >= 1e-6
    # (x < 1e-6 C) switches
    if rng.random() < 0.15:
        return rng.uniform(0.0, 2e-6 * C)
    return center * (1.0 + rng.uniform(-0.1, 0.1))


def _rel(rng, c, r, a=0.0):
    return c * (1.0 + rng.uniform(-r, r)) + rng.uniform(-a, a)


def sample_point(ctype: str, comp, x_star, u_star, rng):
    p = comp.params
    s = list(comp.state_schema)
    ins = [q[0] for q in comp.port_schema["in"]]
    x = np.array(x_star, dtype=float)
    u = np.array(u_star, dtype=float)
    X = dict(zip(s, x))
    U = dict(zip(ins, u))

    if ctype in ("GENROU_PHTRUE", "GENSAL_PHTRUE"):
        X["delta"] += rng.uniform(-0.5, 0.5)
        X["omega"] = rng.uniform(0.95, 1.05)
        for k in s[2:]:
            X[k] = _rel(rng, X[k], 0.2, 0.05)
        if ctype == "GENSAL_PHTRUE" and rng.random() < 0.05:
            # saturation site `psi_ag > A_sat && psi_ag > 1e-6`: a (nearly) dead terminal
            U["Vd"], U["Vq"] = rng.uniform(-1e-7, 1e-7, 2)
        else:
            U["Vd"] += rng.uniform(-0.2, 0.2)
            U["Vq"] += rng.uniform(-0.2, 0.2)
        U["Tm"] = _rel(rng, U["Tm"], 0.3, 0.1)
        U["Efd"] = _rel(rng, U["Efd"], 0.3, 0.1)

    elif ctype == "IEEET1_PHTRUE":
        lo, hi = float(p["VRMIN"]), float(p["VRMAX"])
        X["xr"] = _rel(rng, X["xr"], 0.2, 0.05)
        X["xa"] = _mix(rng, 0.5, lambda: _rel(rng, X["xa"], 0.2, 0.1),
                       lambda: rng.uniform(lo - 0.5, hi + 0.5))
        X["xe"] = _rel(rng, X["xe"], 0.3, 0.1)
        X["xf"] = _rel(rng, X["xf"], 0.3, 0.1)
        X["x_field"] = _reservoir(rng, X["x_field"], float(p["C_FIELD"]))
        U["Vterm"] = _rel(rng, U["Vterm"], 0.3)
        U["Vref"] = _rel(rng, U["Vref"], 0.1)
        U["upss"] = rng.uniform(-0.05, 0.05)
        U["i_fd"] = _rel(rng, U["i_fd"], 0.3, 0.001)

    elif ctype == "IEEEG1_PHTRUE":
        lo, hi = float(p["PMIN"]), float(p["PMAX"])
        span = hi - lo
        for k in ("x0", "x2", "x3", "x4", "x5"):
            X[k] = _rel(rng, X[k], 0.3, 0.1)
        X["x1"] = _mix(rng, 0.5, lambda: _rel(rng, X["x1"], 0.2, 0.1),
                       lambda: rng.uniform(lo - 0.2 * span - 0.1, hi + 0.2 * span + 0.1))
        X["x_steam"] = _reservoir(rng, X["x_steam"], float(p["C_TANK"]))
        U["omega"] = rng.uniform(0.95, 1.05)
        U["Pref"] = _rel(rng, U["Pref"], 0.1, 0.05)
        U["u_agc"] = rng.uniform(-0.05, 0.05)

    elif ctype == "IEEEG3_PHTRUE":
        rb = float(p["R_base"])
        gmin, gmax = float(p["PMIN"]) * rb, float(p["PMAX"]) * rb
        rdn, rup = float(p["UC"]) * rb, float(p["UO"]) * rb
        X["xp"] = _mix(rng, 0.5, lambda: _rel(rng, X["xp"], 0.2, 0.2 * rup),
                       lambda: rng.uniform(1.5 * rdn, 1.5 * rup))
        X["xr"] = _rel(rng, X["xr"], 0.3, 0.1)
        span = gmax - gmin
        X["at"] = _mix(rng, 0.5, lambda: _rel(rng, X["at"], 0.2, 0.1),
                       lambda: rng.uniform(gmin - 0.2 * span, gmax + 0.2 * span))
        X["x1"] = _rel(rng, X["x1"], 0.3, 0.1)
        X["x_water"] = _reservoir(rng, X["x_water"], float(p["C_TANK"]))
        U["omega"] = rng.uniform(0.95, 1.05)
        U["Pref"] = _rel(rng, U["Pref"], 0.3, 0.05)
        U["u_agc"] = rng.uniform(-0.05, 0.05)

    elif ctype == "COMPLEXLOAD":
        # |V| over all four PowerFactory voltage bands, angle near the equilibrium
        ang = math.atan2(U["Vq"], U["Vd"]) + rng.uniform(-0.3, 0.3)
        # 5 %: a collapsed bus, below the |V|^2 > 1e-8 clip of the correction current
        vm = rng.uniform(0.0, 5e-5) if rng.random() < 0.05 else rng.uniform(0.0, 1.4)
        U["Vd"], U["Vq"] = vm * math.cos(ang), vm * math.sin(ang)
        if "z" in X:
            vini = float(p.get("Vini", p.get("V0", 1.0)))
            X["z"] = _mix(rng, 0.1, lambda: X["z"] + rng.uniform(-0.3, 0.3),
                          lambda: rng.uniform(-1.3 * vini, -0.9 * vini))
    elif ctype in ("GFM_VSM_PHTRUE", "GFM_DROOP_PHTRUE"):
        # mostly local draws (virtual impedance inert), sometimes wide ones (it engages)
        wide = rng.random() < 0.4
        X["theta"] += rng.uniform(-0.5, 0.5) if wide else rng.uniform(-0.02, 0.02)
        X["omega"] = rng.uniform(0.95, 1.05)
        X["u_mag"] = _rel(rng, X["u_mag"], 0.2 if wide else 0.02)
        X["x_tank"] = _reservoir(rng, X["x_tank"], float(p["C_TANK"]))
        for k in ("Vd_meas", "Vq_meas"):
            X[k] += rng.uniform(-0.4, 0.4) if wide else rng.uniform(-0.02, 0.02)
        U["Vd"] += rng.uniform(-0.2, 0.2)
        U["Vq"] += rng.uniform(-0.2, 0.2)
        # the pf_frame speed input: unwired (0) or live
        U["omega_ref"] = 0.0 if rng.random() < 0.3 else rng.uniform(0.95, 1.05)
        if ctype == "GFM_VSM_PHTRUE":
            X["omega_f"] = rng.uniform(0.95, 1.05)
            for k in ("Id_meas", "Iq_meas"):
                X[k] = _rel(rng, X[k], 0.3, 0.5)
            # the field flux near either projection bound, or near its equilibrium
            lo, hi, bw = float(p["uf_min"]), float(p["uf_max"]), float(p["uf_band"])
            r = rng.random()
            X["u_field"] = (rng.uniform(hi - 2 * bw, hi + bw) if r < 0.3 else
                            rng.uniform(lo - bw, lo + 2 * bw) if r < 0.6 else
                            _rel(rng, X["u_field"], 0.2))
        else:
            X["q_lpf"] = _rel(rng, X["q_lpf"], 0.3, 0.5)
            # measured current: near the equilibrium, or up to beyond the overcurrent
            # threshold i_lim and the converter limit i_con_lim
            if rng.random() < 0.5:
                Id, Iq = _rel(rng, X["Id_meas"], 0.3, 0.5), _rel(rng, X["Iq_meas"], 0.3, 0.5)
            else:
                m = rng.uniform(0.0, 1.4 * max(float(p["i_lim"]), float(p["i_con_lim"])))
                a = rng.uniform(-math.pi, math.pi)
                Id, Iq = m * math.cos(a), m * math.sin(a)
            X["Id_meas"], X["Iq_meas"] = Id, Iq
            for k, c in (("Id_lpf2", Id), ("Iq_lpf2", Iq), ("Id_lpfz", Id), ("Iq_lpfz", Iq)):
                X[k] = c + rng.uniform(-1.0, 1.0)
            # available capacity: around S_min, above a1/-a0 (Rv < 0), or huge (mp_e < 1e-9)
            r = rng.random()
            X["Sa"] = (rng.uniform(0.0, 1.5) if r < 0.6 else rng.uniform(1.5, 6.0) if r < 0.95
                       else rng.uniform(1e8, 2e8))
            U["S_avail"] = rng.uniform(-0.5, 1.5)

    elif ctype == "GFM_VOC_PHTRUE":
        v = complex(X["v_alpha"], X["v_beta"])
        r = rng.random()
        if r < 0.05:        # |v| < 1e-3: below the pVOC guard v^2 > 1e-6
            a = rng.uniform(-math.pi, math.pi)
            v = rng.uniform(0.0, 1e-3) * complex(math.cos(a), math.sin(a))
        elif r < 0.45:      # wide: the virtual impedance engages
            v *= rng.uniform(0.6, 1.4) * complex(math.cos(a := rng.uniform(-0.5, 0.5)), math.sin(a))
        else:
            v *= rng.uniform(0.98, 1.02) * complex(math.cos(a := rng.uniform(-0.02, 0.02)), math.sin(a))
        X["v_alpha"], X["v_beta"] = v.real, v.imag
        X["x_tank"] = _reservoir(rng, X["x_tank"], float(p["C_TANK"]))
        wide = r < 0.45
        for k in ("Vd_meas", "Vq_meas"):
            X[k] += rng.uniform(-0.4, 0.4) if wide else rng.uniform(-0.02, 0.02)
        for k in ("Id_meas", "Iq_meas", "Id_lpf2", "Iq_lpf2", "p_pvoc", "q_pvoc"):
            X[k] = _rel(rng, X[k], 0.3, 0.5)
        if rng.random() < 0.1:   # a (nearly) bolted terminal: |V| < 0.01 pu
            m, a = rng.uniform(0.0, 0.02), rng.uniform(-math.pi, math.pi)
            U["Vd"], U["Vq"] = m * math.cos(a), m * math.sin(a)
        else:
            U["Vd"] += rng.uniform(-0.2, 0.2)
            U["Vq"] += rng.uniform(-0.2, 0.2)
        for k in ("omega_ref", "omega_aux"):
            U[k] = 0.0 if rng.random() < 0.3 else rng.uniform(0.95, 1.05)

    else:
        raise ValueError(f"no sampling box for {ctype}")

    return np.array([X[k] for k in s]), np.array([U[k] for k in ins])


# ---------------------------------------------------------------------------
# Section
# ---------------------------------------------------------------------------

def _num_params(params: dict) -> dict:
    out = {}
    for k, v in params.items():
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)) and math.isfinite(float(v)):
            out[k] = float(v)
    return out


def build_initialised(phps_root, system_rel):
    """Run PHPS's initialisation (Python path, no C++) and return the runner, the solver
    helper and the equilibrium inputs of every component."""
    from src.dirac.dae_runner import DiracRunner
    from src.dirac.py_solver import PyDAESolver

    system_path = phps_root / "phps" / "cases" / system_rel
    with contextlib.redirect_stdout(io.StringIO()):
        runner = DiracRunner(str(system_path), output_dir=tempfile.mkdtemp())
        # parameters as the constructors left them, before initialisation touches them
        pre_init = {c.name: _num_params(dict(c.params))
                    for c in runner.base_compiler.components}
        runner.build(solver="scipy")
        solver = PyDAESolver(runner)
    x0 = runner.x0.copy()

    # Equilibrium inputs: the first half of PyDAESolver.rhs at x0 (network solve, dq frames,
    # output pass), then each component's gathered inputs.
    Vd, Vq = solver._solve_network(x0, 0.0)
    Vterm = np.sqrt(Vd ** 2 + Vq ** 2)
    vd_dq, vq_dq = {}, {}
    for comp in solver.components:
        if comp.name in solver._gen_dq_params:
            off = solver.state_offsets[comp.name]
            vd, vq, _, _ = solver._compute_dq_frame(comp.name, x0, Vd, Vq, off)
            vd_dq[comp.name], vq_dq[comp.name] = vd, vq
    for comp in solver.components:
        off = solver.state_offsets[comp.name]
        n = len(comp.state_schema)
        inp = solver._gather_inputs(comp, Vd, Vq, Vterm, vd_dq, vq_dq)
        solver._out_fns[comp.name](x0[off:off + n], inp, solver._outputs[comp.name], 0.0)
    u_star = {}
    for comp in solver.components:
        u_star[comp.name] = solver._gather_inputs(comp, Vd, Vq, Vterm, vd_dq, vq_dq).copy()
    return runner, solver, x0, u_star, pre_init


def _evaluate(kstep, kout, comp, x, u):
    n_out = len(comp.port_schema["out"])
    outputs = np.zeros(n_out)
    br_out = kout(x, u, outputs, 0.0)
    out_pass1 = outputs.copy()
    dxdt = np.zeros(len(comp.state_schema))
    br_step = kstep(x, dxdt, u, outputs, 0.0)
    return dxdt, out_pass1, outputs.copy(), br_out, br_step


def _same_bits(a, b) -> bool:
    return np.array_equal(np.asarray(a, float).view(np.int64), np.asarray(b, float).view(np.int64))


def section_components(phps_root, case_by_name: dict, seed: int = 20260929) -> dict:
    """Return {type: record} for every type in COMPONENT_SOURCES."""
    from src.dirac.py_codegen import make_out_func, make_step_func

    built = {}
    records = {}
    for ctype, case_name in sorted(COMPONENT_SOURCES.items()):
        # the seed offset is the type's place in SEED_ORDER, so adding a type does not
        # change the samples of the others
        t_index = SEED_ORDER.index(ctype)
        if case_name not in built:
            built[case_name] = build_initialised(phps_root, case_by_name[case_name])
        runner, solver, x0, u_star, pre_init = built[case_name]
        case_json = runner.base_compiler.graph.raw_data["components"]
        comps = [c for c in solver.components if case_json[c.name]["type"] == ctype]
        if not comps:
            raise RuntimeError(f"no {ctype} in case {case_name}")

        instances = {}
        kernels = {}
        for c in comps:
            off = solver.state_offsets[c.name]
            n = len(c.state_schema)
            used = _num_params(c.params)
            before = pre_init[c.name]
            # everything PHPS's initialisation added or changed
            init_set = {k: v for k, v in used.items() if k not in before or before[k] != v}
            instances[c.name] = {
                "params_used": used,
                "init_set": init_set,
                "x_star": [float(v) for v in x0[off:off + n]],
                "u_star": [float(v) for v in u_star[c.name]],
            }
            kernels[c.name] = (Kernel(c, "step"), Kernel(c, "out"),
                               make_step_func(c), make_out_func(c))

        # parameter variants of the first instance (VARIANTS)
        variants = []
        for suffix, overrides in VARIANTS.get(ctype, []):
            base = comps[0]
            v = copy.deepcopy(base)
            v.name = f"{base.name}__{suffix}"
            v.params.update(overrides)
            ib = instances[base.name]
            instances[v.name] = {
                "variant_of": base.name,
                "overrides": dict(overrides),
                "params_used": _num_params(v.params),
                "init_set": ib["init_set"],
                "x_star": ib["x_star"],
                "u_star": ib["u_star"],
            }
            kernels[v.name] = (Kernel(v, "step"), Kernel(v, "out"),
                               make_step_func(v), make_out_func(v))
            variants.append(v)
        # the variants take about half of the samples
        order = comps + variants * len(comps)

        first = kernels[comps[0].name]
        for c in comps[1:] + variants:
            for kidx in (0, 1):
                a = [s["text"] for s in first[kidx].sites]
                b = [s["text"] for s in kernels[c.name][kidx].sites]
                if a != b:
                    raise RuntimeError(f"{ctype}: {c.name} has different branch sites "
                                       f"from {comps[0].name}; sample them separately")
        rng = np.random.default_rng(seed + t_index)
        samples = []
        redrawn = 0
        while len(samples) < N_SAMPLES:
            c = order[len(samples) % len(order)]
            kstep, kout, ref_step, ref_out = kernels[c.name]
            inst = instances[c.name]
            x, u = sample_point(ctype, c, inst["x_star"], inst["u_star"], rng)
            dxdt, o1, o2, br_out, br_step = _evaluate(kstep, kout, c, x, u)

            # redraw samples within rounding distance of a switching surface (conditions
            # inside a loop are a root solve's, not switching surfaces)
            def switching(br, kern):
                loop = {st["id"] for st in kern.sites if st.get("in_loop")}
                return [b for b in br if b[0] not in loop]
            near = False
            for _ in range(N_PERTURB):
                xp = x * (1.0 + PERTURB_REL * rng.uniform(-1, 1, x.size)) \
                    + PERTURB_REL * rng.uniform(-1, 1, x.size)
                up = u * (1.0 + PERTURB_REL * rng.uniform(-1, 1, u.size)) \
                    + PERTURB_REL * rng.uniform(-1, 1, u.size)
                _, _, _, bo, bs = _evaluate(kstep, kout, c, xp, up)
                if (switching(bo, kout) != switching(br_out, kout)
                        or switching(bs, kstep) != switching(br_step, kstep)):
                    near = True
                    break
            if near:
                redrawn += 1
                continue

            # the instrumented kernels must reproduce PHPS's own functions bit for bit
            r_out = np.zeros(len(c.port_schema["out"]))
            ref_out(x, u, r_out, 0.0)
            r1 = r_out.copy()
            r_dx = np.zeros(len(c.state_schema))
            ref_step(x, r_dx, u, r_out, 0.0)
            if not (_same_bits(r1, o1) and _same_bits(r_out, o2) and _same_bits(r_dx, dxdt)):
                raise RuntimeError(f"instrumented kernel differs from PHPS for {c.name}")

            samples.append({
                "component": c.name,
                "x": [float(v) for v in x],
                "u": [float(v) for v in u],
                "dxdt": [float(v) for v in dxdt],
                "outputs_out": [float(v) for v in o1],
                "outputs_step": [float(v) for v in o2],
                "H": float(c.hamiltonian(np.array(x))),
                "grad_H": [float(v) for v in c.grad_hamiltonian(np.array(x))],
                "branches_out": [[k, b] for k, b in br_out],
                "branches_step": [[k, b] for k, b in br_step],
            })

        # coverage: every state-dependent site must have gone both ways
        coverage = {}
        for mode, kidx in (("out", 1), ("step", 0)):
            sites = first[kidx].sites
            seen = {s["id"]: set() for s in sites}
            for smp in samples:
                for k, b in smp[f"branches_{mode}"]:
                    seen[k].add(b)
            coverage[mode] = {str(k): sorted(v) for k, v in seen.items()}
            for s in sites:
                reason = UNREACHED.get(ctype, {}).get(s["text"])
                if reason is not None:
                    if seen[s["id"]]:
                        raise RuntimeError(f"{ctype} {mode} site {s['id']} ({s['text']}) "
                                           "is listed as unreached but was reached")
                    s["unreached"] = reason
                    continue
                if not s["param_only"] and seen[s["id"]] != {False, True}:
                    raise RuntimeError(
                        f"{ctype} {mode} site {s['id']} ({s['text']}) not exercised both "
                        f"ways: {sorted(seen[s['id']])}; widen its sampling box")

        records[ctype] = {
            "type": ctype,
            "case": case_name,
            "phps_class": type(comps[0]).__name__,
            "seed": seed + t_index,
            "n_samples": N_SAMPLES,
            "redrawn_near_switching": redrawn,
            "states": list(comps[0].state_schema),
            "inputs": [q[0] for q in comps[0].port_schema["in"]],
            "outputs": [q[0] for q in comps[0].port_schema["out"]],
            "sites": {"step": first[0].sites, "out": first[1].sites},
            "site_coverage": coverage,
            "instances": instances,
            "samples": samples,
        }
    return records
