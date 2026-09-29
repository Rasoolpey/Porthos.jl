# PowerFactory driver

`pf/` runs DIgSILENT PowerFactory for Porthos: it dumps the PowerFactory model and runs the
RMS simulation of a Porthos scenario, so that Porthos can be compared with PowerFactory on
the same case. It is the one place besides `parity/generate/` where Python does more than
edit JSON: PowerFactory's API is Python only. The rule stays that **Python computes
nothing**. `pf/` only operates PowerFactory and writes what PowerFactory computed (its own
CSV export) plus JSON records. Every comparison is done in Julia (`src/io/powerfactory.jl`,
`scripts/pf_compare_fault.jl`).

The code is written for PowerFactory 2022 SP1 and uses Python's standard library only.
No virtual environment is needed.

## Use

Close the PowerFactory window first. The API starts PowerFactory as a separate engine,
which cannot run while the window is open (exit code 4002).

```
julia --project=. scripts/pf_compare_fault.jl [scenario.json]    # PowerFactory + Porthos, compared
py -3.10 pf/run.py inspect                                       # model and load flow -> outputs/pf/model/model.json
py -3.10 pf/run.py simulate cases/IEEE39Bus_PF/bus_fault_bus16_150ms.json
```

From Julia: `pf_simulate(scenario)` (returns `PFResults`), `pf_inspect()`,
`read_pf_results(dir)`, `pf_signal(r, "G 04", "s:speed")`, `pf_compare(...)`.

Outputs go to `outputs/pf/<scenario name>/` (git-ignored):

| File | Contents |
|---|---|
| `pf_results.csv` | PowerFactory's CSV export, full precision (scientific, 15 digits): time, then one column per recorded variable; two header rows (object, description) |
| `run.json` | project, study cases, the request (scenario, events, stop time), the PowerFactory events as created (ohmic fault impedance), solver settings, return codes, error and warning messages, wall time, and the **column map** (object, class, variable, unit per column), read from the result object |
| `model.json` | `inspect`: buses, lines, transformers, loads, machines (element and type data), composite models, DSL controllers (parameter names and values), load flow |
| `compare.json`, `pf_*.png` | written by `scripts/pf_compare_fault.jl` |

## Configuration

`pf/config.json`: PowerFactory folder and Python version, project
(`39 Bus New England System`), study case (`Base`), the working study case (`Porthos`), the
terminal naming (`Bus {:02d}`), the RMS step (1 ms) and start time, and the recorded
variables per class. Machine-specific overrides go into `pf/config.local.json`
(git-ignored), which updates the keys it contains.

## What the driver does to the PowerFactory project

- `simulate` copies the configured study case to `Porthos` (deleting the previous copy)
  and works only there. It takes the copy's events out of service, creates its own
  `EvtShc` pairs, and adds a result object. Your study case is not changed.
- `inspect` activates your study case and runs a load flow, which changes no input data.
- Network-element parameters are shared by all study cases. Nothing changes them now. Any
  future tool that does must use `session.Restore`, which undoes every change in reverse
  order.

## PowerFactory 2022 SP1: facts learned

- **Engine start.** `GetApplicationExt` starts an engine, it does not attach to an open
  window, so close the window first.
- **INI workaround.** The installed `PowerFactory.ini` has `type = server` under
  `[license]`, and the 2022 engine rejects that key ("Unknown parameter"). The driver starts
  the engine with `/ini <copy>`, a copy without the keys in `ini_drop`, written to
  `outputs/pf/engine.ini` (git-ignored: it holds the licence container id). The copy's path
  must not contain spaces (the argument is passed unquoted). The installation is not
  modified.
- **Time units.** `ComInc.dtgrd` and `ComInc.tstart` are in **milliseconds** (defaults
  10 ms and -100 ms); `ComSim.tstop` is in seconds.
- **Result access.** The result object (`ElmRes`) belongs to `ComInc.p_resvar`; `ComSim`
  has no `p_resvar` of its own. After `res.Load()`, `res.GetObject(k)`,
  `GetVariable(k)`, `GetUnit(k)` and `FindColumn(obj, var)` work. `GetColumnValues` fails,
  and `GetValue(row, col)` costs about 44 us per value.
- **Export.** The CSV export (`ComRes`, `iopt_exp = 6`) writes 6 decimals by default. Set
  `numberFormat = 1` and `numberPrecisionScientific = 15` for full precision. Its column
  order is the result object's, not the registration order. The second header row holds
  descriptions ("u, Magnitude in p.u."), not variable names, hence the column map in
  `run.json`.
- **Event times.** At each event PowerFactory writes two rows with the same time, before
  and after the event.
- **The 39-bus project.** The ten `ElmSym` machines are G 01 to G 10. G 02 at bus 31 is the
  reference machine (`ip_ctrl = 1`); G 01 is "Rest of U.S.A. / Canada" (10 000 MVA, no
  controllers); G 05 has two parallel units (`ngnum = 2`), which Porthos models as two
  machines. The controllers are IEEET1 AVRs, IEEEG1 governors (IEEEG3 at G 10), and
  `pss_CONV` stabilisers, all out of service. Loads use the General Load Type
  (`kpu = kqu = 2`).
- **Traps inherited from the PHPS scripts (`PHPS_Opt/pf/api`, PowerFactory 2024).**
  PowerFactory binds a study case's events at initialisation. An event deleted after a
  simulation in the same process can stay live. `Delete()` returns an error code instead of
  raising. `CreateObject` with an existing name makes `name(1)`. The driver avoids all four
  by working in a fresh copy and running one simulation per process.

## Result on the base case (2026-09-29)

`scripts/pf_compare_fault.jl` on the bus-16 fault (150 ms), PowerFactory 2022 SP1 against
Porthos IDA at production tolerances:
- **Load flow:** Porthos matches PowerFactory within 2.2e-9 pu and 1.8e-6 deg.
- **Rotor angles** (relative to G 01): 0.38 to 0.65 deg RMS, at most 1.32 deg.
- **Speeds:** RMS below 9.5e-5 pu.
- **Terminal P:** 0.015 to 0.035 pu RMS (0.15 pu for the 10 000 MVA equivalent G 01).

Every metric is smaller than PHPS's recorded comparison against PowerFactory 2024
(`PHPS_Opt/study/pf_reference/model_review/pf_compare/pf_compare_bus16_v4_x1e-5.json`:
0.48 to 0.79 deg RMS, at most 1.68 deg).
