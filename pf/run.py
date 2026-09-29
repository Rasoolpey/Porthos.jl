"""Porthos PowerFactory entry point. Close the PowerFactory window first.

    py -3.10 pf/run.py inspect [--out DIR]
        Dump the model of the configured study case (network, machines, controllers, load
        flow) to DIR/model.json (default outputs/pf/model/).

    py -3.10 pf/run.py simulate SCENARIO.json [--out DIR] [--dt-ms 1.0]
        RMS simulation of a Porthos scenario's bus faults in a fresh copy of the study case;
        writes DIR/pf_results.csv and DIR/run.json (default outputs/pf/<scenario name>/).

Settings: pf/config.json (project, study case, PowerFactory folder, recorded variables),
overridden by pf/config.local.json if it exists.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from porthospf import inspect_model, session, simulate  # noqa: E402
from porthospf.session import PFError  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pf/run.py", description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("inspect", help="dump the PowerFactory model to model.json")
    p.add_argument("--out", default=None)
    p = sub.add_parser("simulate", help="RMS run of a Porthos scenario's bus faults")
    p.add_argument("scenario")
    p.add_argument("--out", default=None)
    p.add_argument("--dt-ms", type=float, default=None)
    a = ap.parse_args(argv)

    cfg = session.load_config()
    t0 = time.time()
    try:
        _, app = session.start(cfg)
        if a.cmd == "inspect":
            case, _ = session.activate(app, cfg, working=False)
            out = a.out or os.path.join(session.OUT_DIR, "model")
            model = inspect_model.dump(app)
            model.update({"created_utc": session.utc_now(), "project": cfg["project"],
                          "study_case": case.loc_name,
                          "powerfactory_dir": cfg["powerfactory_dir"]})
            path = session.write_json(os.path.join(out, "model.json"), model)
        else:
            request = simulate.request_from_scenario(a.scenario, cfg)
            case, _ = session.activate(app, cfg, working=True)
            name = os.path.splitext(os.path.basename(a.scenario))[0]
            out = a.out or os.path.join(session.OUT_DIR, name)
            path = simulate.run(app, cfg, case, request, out, dt_ms=a.dt_ms)
    except PFError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print(f"{a.cmd}: {path}  ({time.time() - t0:.1f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
