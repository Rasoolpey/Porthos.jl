"""RMS simulation of a Porthos scenario in PowerFactory.

The scenario JSON (the format Porthos reads) gives the events and the duration; its
`BusFault` events become PowerFactory short-circuit events (`EvtShc`) on the terminal
named by `bus_name` (config), with the fault impedance converted from per unit on the
case's MVA base to ohms on the terminal's nominal voltage. Other event types are refused.

The run happens in the working study case (a fresh copy of the user's; see `session`):
its own events are taken out of service, the fault events are created, the result
variables are registered, and PowerFactory's CSV export and a `run.json` are written.
`run.json` maps each registered (object, variable) to its result column, so the reader
does not depend on the header texts PowerFactory writes.
"""

import json
import os
import time

from .session import PFError, name_of, output_lines, safe_get, utc_now, write_json


def request_from_scenario(path, cfg):
    """The PowerFactory run a Porthos scenario asks for."""
    with open(path, encoding="utf-8") as fh:
        sc = json.load(fh)
    system = os.path.normpath(os.path.join(os.path.dirname(path), sc["system"]))
    with open(system, encoding="utf-8") as fh:
        mva_base = float(json.load(fh).get("config", {}).get("mva_base", 100.0))
    events = []
    for e in sc.get("events", []):
        if e.get("type") != "BusFault":
            raise PFError(f"event type {e.get('type')!r} is not supported (BusFault only)")
        t0 = float(e["t_start"])
        t1 = (t0 + float(e["t_duration"]) if "t_duration" in e else
              float(e["t_end"]) if "t_end" in e else t0 + 0.1)
        events.append({"type": "BusFault", "bus": int(e["bus"]), "t_start": t0, "t_end": t1,
                       "r_pu": float(e.get("r", 0.0)), "x_pu": float(e.get("x", 1e-5))})
    return {"scenario": os.path.abspath(path), "system": system, "mva_base": mva_base,
            "t_stop": float(sc["solver"]["duration"]), "events": events}


def _terminal(app, name):
    t = next((b for b in app.GetCalcRelevantObjects("*.ElmTerm") if b.loc_name == name), None)
    if t is None:
        raise PFError(f"no terminal {name!r}")
    return t


def _set_events(app, cfg, request):
    evt = app.GetFromStudyCase("IntEvt")
    if evt is None:
        raise PFError("the study case has no event folder (IntEvt)")
    for e in evt.GetContents("*"):          # copies of the base case's events
        e.SetAttribute("outserv", 1)
    made = []
    for k, ev in enumerate(request["events"]):
        term = _terminal(app, cfg["bus_name"].format(ev["bus"]))
        zb = float(safe_get(term, "uknom")) ** 2 / request["mva_base"]      # ohm
        a = evt.CreateObject("EvtShc", f"porthos_{k:02d}_apply")
        a.SetAttribute("p_target", term)
        a.SetAttribute("i_shc", 0)                    # three-phase short circuit
        a.SetAttribute("time", ev["t_start"])
        a.SetAttribute("R_f", ev["r_pu"] * zb)
        a.SetAttribute("X_f", ev["x_pu"] * zb)
        c = evt.CreateObject("EvtShc", f"porthos_{k:02d}_clear")
        c.SetAttribute("p_target", term)
        c.SetAttribute("i_shc", 4)                    # clear
        c.SetAttribute("time", ev["t_end"])
        for o in (a, c):
            o.SetAttribute("outserv", 0)
        made.append({"bus": term.loc_name, "apply_s": safe_get(a, "time"),
                     "clear_s": safe_get(c, "time"), "R_f_ohm": safe_get(a, "R_f"),
                     "X_f_ohm": safe_get(a, "X_f")})
    live = [e for e in evt.GetContents("*") if not safe_get(e, "outserv")]
    if len(live) != 2 * len(request["events"]):
        raise PFError(f"expected {2 * len(request['events'])} live events, found {len(live)}")
    return made


def _notable(lines, limit=50):
    """The error and warning lines of an output-window capture."""
    return [l for l in lines if "error" in l.lower() or "warn" in l.lower()][:limit]


def run(app, cfg, study_case, request, out_dir, dt_ms=None):
    """Run the request in the active (working) study case; returns the path of run.json."""
    rms = cfg["rms"]
    dt_ms = float(dt_ms if dt_ms is not None else rms["dt_ms"])
    events = _set_events(app, cfg, request)

    inc = app.GetFromStudyCase("ComInc")
    sim = app.GetFromStudyCase("ComSim")
    res = study_case.CreateObject("ElmRes", "Porthos results")
    inc.SetAttribute("iopt_sim", "rms")
    inc.SetAttribute("iopt_net", "sym")
    inc.SetAttribute("iopt_adapt", 0)                  # fixed step
    inc.SetAttribute("dtgrd", dt_ms)                   # ms in PowerFactory 2022
    inc.SetAttribute("tstart", float(rms["t_start_ms"]))
    inc.SetAttribute("p_resvar", res)
    app.ClearOutputWindow()
    rc_inc = inc.Execute()
    msgs_inc = output_lines(app)
    if rc_inc != 0:
        raise PFError("initial conditions failed (ComInc rc %s):\n  %s"
                      % (rc_inc, "\n  ".join(_notable(msgs_inc) or msgs_inc[-20:])))

    missing = []
    for cls, variables in cfg["variables"].items():
        for obj in app.GetCalcRelevantObjects("*." + cls):
            if safe_get(obj, "outserv"):
                continue
            for var in variables:
                try:
                    ok = res.AddVariable(obj, var) == 0
                except Exception:  # noqa: BLE001
                    ok = False
                if not ok:
                    missing.append({"object": obj.loc_name, "class": cls, "variable": var})
    res.InitialiseWriting()

    sim.SetAttribute("tstop", request["t_stop"])
    app.ClearOutputWindow()
    t0 = time.time()
    rc_sim = sim.Execute()
    wall = time.time() - t0
    msgs_sim = output_lines(app)

    # the column map, from the result object itself (column k is CSV column k + 2, after
    # the time column; the CSV's order is not the registration order)
    res.Load()
    columns = []
    for k in range(res.GetNumberOfColumns()):
        obj = res.GetObject(k)
        columns.append({"column": k, "object": name_of(obj),
                        "class": obj.GetClassName() if obj is not None else None,
                        "variable": res.GetVariable(k), "unit": res.GetUnit(k),
                        "description": res.GetDescription(k)})
    n_rows = res.GetNumberOfRows()
    res.Release()

    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "pf_results.csv")
    if os.path.exists(csv_path):
        os.remove(csv_path)
    exp = app.GetFromStudyCase("ComRes")
    exp.SetAttribute("pResult", res)
    exp.SetAttribute("iopt_exp", 6)                    # CSV
    exp.SetAttribute("iopt_sep", 0)                    # own separators, not the system's
    exp.SetAttribute("col_Sep", ",")
    exp.SetAttribute("dec_Sep", ".")
    exp.SetAttribute("iopt_honly", 0)
    exp.SetAttribute("iopt_csel", 0)                   # every registered variable
    exp.SetAttribute("numberFormat", 1)                # scientific, full double precision
    exp.SetAttribute("numberPrecisionScientific", 15)
    exp.SetAttribute("f_name", csv_path)
    rc_exp = exp.Execute()
    for _ in range(50):                                # the file is written asynchronously
        if os.path.exists(csv_path):
            break
        time.sleep(0.1)
    if not os.path.exists(csv_path):
        raise PFError(f"the result export wrote no file (ComRes rc {rc_exp})")

    record = {
        "created_utc": utc_now(),
        "powerfactory_dir": cfg["powerfactory_dir"],
        "project": cfg["project"],
        "study_case": study_case.loc_name,
        "base_study_case": cfg["study_case"],
        "request": request,
        "events": events,
        "settings": {"simulation": safe_get(inc, "iopt_sim"), "network": safe_get(inc, "iopt_net"),
                     "dt_ms": safe_get(inc, "dtgrd"), "t_start_ms": safe_get(inc, "tstart"),
                     "t_stop": safe_get(sim, "tstop"), "adaptive": safe_get(inc, "iopt_adapt")},
        "rc": {"ComInc": rc_inc, "ComSim": rc_sim, "ComRes": rc_exp},
        "messages": {"ComInc": _notable(msgs_inc), "ComSim": _notable(msgs_sim)},
        "wall_time_s": wall,
        "csv": os.path.basename(csv_path),
        "rows": n_rows,
        "columns": columns,
        "not_recorded": missing,
    }
    return write_json(os.path.join(out_dir, "run.json"), record)
