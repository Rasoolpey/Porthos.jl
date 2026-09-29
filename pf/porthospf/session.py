"""PowerFactory session: configuration, engine start, project and study case.

PowerFactory runs as an engine (no window) started by `GetApplicationExt`. Two facts about
PowerFactory 2022 SP1 on the build machine shape this module:

- The engine cannot start while the PowerFactory window is open (exit code 4002: the
  licence or the database is in use). Close the window first.
- The installed `PowerFactory.ini` has a `type` key under `[license]` that the 2022 engine
  rejects ("Unknown parameter"). The engine is started with a copy of the INI without the
  keys listed in `ini_drop` (`pf/config.json`), written to `outputs/pf/engine.ini`. The
  installation is never modified.

Every run works in a fresh copy of the configured study case (`working_study_case`,
default "Porthos"), made at the start of the process, so the user's study case is never
changed. Network-element parameters are shared by all study cases, so any change to them
must go through `Restore`.
"""

import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PF_DIR = os.path.join(ROOT, "pf")
OUT_DIR = os.path.join(ROOT, "outputs", "pf")


class PFError(RuntimeError):
    """A PowerFactory step failed; the message says what to do."""


def load_config():
    """`pf/config.json`, updated by `pf/config.local.json` if present (machine-specific)."""
    with open(os.path.join(PF_DIR, "config.json"), encoding="utf-8") as fh:
        cfg = json.load(fh)
    local = os.path.join(PF_DIR, "config.local.json")
    if os.path.exists(local):
        with open(local, encoding="utf-8") as fh:
            cfg.update(json.load(fh))
    return cfg


def safe_get(obj, attr, default=None):
    """An attribute, or `default` when it does not exist or cannot be read yet."""
    if obj is None:
        return default
    try:
        return obj.GetAttribute(attr)
    except Exception:  # noqa: BLE001 - PowerFactory raises plain exceptions
        return default


def name_of(obj):
    return getattr(obj, "loc_name", None) if obj is not None else None


class Restore:
    """Record every attribute change and undo them in reverse order."""

    def __init__(self):
        self._log = []

    def set(self, obj, attr, value):
        self._log.append((obj, attr, safe_get(obj, attr)))
        obj.SetAttribute(attr, value)

    def restore(self):
        for obj, attr, old in reversed(self._log):
            try:
                obj.SetAttribute(attr, old)
            except Exception:  # noqa: BLE001
                pass
        self._log = []


def _engine_ini(cfg):
    """A copy of the installed INI without the keys the engine rejects, or None."""
    src = os.path.join(cfg["powerfactory_dir"], "PowerFactory.ini")
    drop = {k.lower(): {x.lower() for x in v} for k, v in cfg.get("ini_drop", {}).items()}
    if not os.path.exists(src) or not drop:
        return None
    section, out, dropped = None, [], []
    with open(src, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("[") and s.endswith("]"):
                section = s[1:-1].strip().lower()
            elif "=" in s and section in drop and s.split("=", 1)[0].strip().lower() in drop[section]:
                dropped.append(f"[{section}] {s.split('=', 1)[0].strip()}")
                continue
            out.append(line)
    if not dropped:
        return None
    os.makedirs(OUT_DIR, exist_ok=True)
    dst = os.path.join(OUT_DIR, "engine.ini")
    if " " in dst:
        raise PFError(f"the engine INI path must not contain spaces: {dst}")
    with open(dst, "w", encoding="utf-8") as fh:
        fh.writelines(out)
    return dst


def start(cfg):
    """Start the PowerFactory engine. Returns (powerfactory module, application)."""
    want = tuple(int(x) for x in cfg["python_version"].split("."))
    if sys.version_info[:2] != want:
        raise PFError(f"run with Python {cfg['python_version']} (PowerFactory's module is built "
                      f"for it): py -{cfg['python_version']} pf/run.py ...")
    pydir = os.path.join(cfg["powerfactory_dir"], "Python", cfg["python_version"])
    if not os.path.isdir(pydir):
        raise PFError(f"no PowerFactory Python module at {pydir} (set powerfactory_dir in "
                      f"pf/config.local.json)")
    if pydir not in sys.path:
        sys.path.insert(0, pydir)
    import powerfactory as pf  # noqa: E402 - path set just above

    ini = _engine_ini(cfg)
    args = f"/ini {ini}" if ini else None
    try:
        app = pf.GetApplicationExt(None, None, args)
    except Exception as e:  # noqa: BLE001 - pf.ExitError
        msg = str(e)
        hint = (" Close the PowerFactory window and run again: the engine cannot start while "
                "it is open." if "4002" in msg else "")
        raise PFError(f"PowerFactory did not start ({msg}).{hint}") from e
    if app is None:
        raise PFError("PowerFactory did not start (no application; see the messages above)")
    return pf, app


def activate(app, cfg, working=True):
    """Activate the project and the study case; with `working`, a fresh copy of it.

    Returns (study case object, base study case name).
    """
    if app.ActivateProject(cfg["project"]) != 0:
        raise PFError(f"cannot activate project {cfg['project']!r}")
    folder = app.GetProjectFolder("study")
    cases = {c.loc_name: c for c in folder.GetContents("*.IntCase")}
    base = cases.get(cfg["study_case"])
    if base is None:
        raise PFError(f"no study case {cfg['study_case']!r}; there are {sorted(cases)}")
    if not working:
        base.Activate()
        return base, base.loc_name
    wname = cfg["working_study_case"]
    if wname == cfg["study_case"]:
        raise PFError("working_study_case must differ from study_case")
    base.Activate()
    old = cases.get(wname)
    if old is not None:
        old.Delete()
    wc = folder.AddCopy(base, wname)
    if wc is None:
        raise PFError(f"cannot copy study case {cfg['study_case']!r} to {wname!r}")
    wc.Activate()
    return wc, base.loc_name


def output_lines(app):
    """The messages in PowerFactory's output window."""
    ow = app.GetOutputWindow()
    return [str(x) for x in (ow.GetContent() or [])] if ow else []


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, allow_nan=True)
        fh.write("\n")
    return path
