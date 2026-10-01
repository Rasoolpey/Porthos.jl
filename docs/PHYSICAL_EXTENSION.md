# Physical extension: an electromagnetic dq benchmark with two-way storage

Design note, 2026-09-30. Status: agreed design, not yet implemented. Context: `TODO.md`,
"Rule-4 decision" and "Physical-extension design decision".

## Why

On the reduced RMS IEEE-39 model the storage completion `V_completion` certifies a strict
local Lyapunov function, but only with a broad repair: 202 of 263 component-pair groups
stay significant, the repair reshapes the machines' own Hessian blocks, and cond(P) is about
1e5. The obstruction already appears on the integrable-loss equivalent, where every network
supply is exact, so it comes from the controllers and the algebraic reduction, not only from
the losses. The DC1A, IEEEG1 and IEEEG3 blocks are reduced transfer models (IEEE 421.5-2016;
IEEE PES-TR1): their states carry no identified physical energy. The extension gives every
controller power path a physical, two-way storage and keeps the fast electromagnetic states,
so that the plant Hamiltonian is physical by construction.  The regulator may add a separate
controller storage, which is labelled as such rather than counted as electromagnetic energy.

## Scope

A separate balanced, averaged electromagnetic model in the synchronous dq frame (stiff, not
switching EMT). The frozen RMS case models, their parity and the `V_P` certificate are
untouched. Every element is a port-Hamiltonian component, `x' = (J - R) grad H + g u`, and the
interconnection is power preserving (Kirchhoff/Dirac, as in Gernandt et al. 2021).

| Element | States | Storage `H` | Ports, dissipation |
|---|---|---|---|
| Synchronous machine | stator fluxes `psi_d, psi_q`; field `lambda_f`; dampers; `omega` | `lambda' L^-1 lambda / 2 + J omega^2 / 2` | terminal `(v_dq, i_dq)`, shaft `(T_m, omega)`, field `(v_f, i_f)`, `i_f = dH/d lambda_f`; winding resistances |
| Excitation | DC-link voltage `v_dc` | `C_dc v_dc^2 / 2` | ideal averaged converter `v_f = m v_dc`, `i_dc = m i_f` (power preserving); DC and field losses explicit; modulation `m` from a simple passive AVR |
| Hydro prime mover | penstock flow `q`, head `h`, gate-servo state | `L_h q^2 / 2 + C_h (h - h*)^2 / 2 + H_servo` | turbine shaft power, friction dissipation; upstream reservoir head |
| Steam prime mover | servo state; steam-chest and reheater pressure or mass | servo storage + shifted (convex) availability, quadratic near equilibrium | enthalpy-flow ports to shaft power; boiler boundary |
| π line | series current `i_dq`, shunt voltages `v_dq` | `L |i|^2 / 2 + C |v|^2 / 2` | `R`, `G` dissipation; the `omega L`, `omega C` terms skew |
| Load | RLC branch states | inductor and capacitor energy | resistive dissipation |

Rules that keep the energy honest:
- **No double count.** The field-winding magnetic energy belongs to the machine Hamiltonian;
  the excitation system adds only the DC link and the converter.
- **No open supplies.** The DC source, the upstream reservoir head and the boiler are either
  modelled or declared fixed chemostatted boundaries; every attraction claim is then stated
  conditional on those boundaries (exergetic bookkeeping as in Lohmayer, Kotyczka and
  Leyendecker 2021).
- **No guessed energies on lag states.** Physical parameters (inductances, capacitances,
  penstock and surge-tank data, steam volumes) come from stated data; the RMS time constants
  are matched only as the reduced limit, through an explicit parameter map.
- **Rotating DC exciter later,** only with its winding, armature-reaction and shaft data
  (DC1A's `xe, TE, KE` do not determine that energy).

## Gates

Build and pass in this order; do not scale up past a failed gate.
1. One hydro unit to an infinite bus (hydro first: the IEEEG3 mode, `0.027 /s`, is the
   spectral bottleneck of the IEEE-39 study).
2. One steam unit to an infinite bus.
3. Two machines, steam and hydro, to expose the swing and AVR inter-unit coupling.
4. IEEE-39 only after gates 1 to 3.

At every gate:
- exact component and total power balances along the flow (the identity pattern of
  `strain_balance`, errors at round-off);
- a convex equilibrium-shifted Hamiltonian (Bregman Hessian positive on the quotient);
- local decay with only own-unit or nearest-neighbour terms (the completion LMI and
  `generalized_decay_rate`; interval check by `verified_lyapunov`);
- recovery of the RMS model when the fast derivatives are set to zero (stator, line and load
  stamps equal to the Y-bus model at nominal frequency).
If gate 3 again needs broad cross-unit terms, record that the legacy controller structure,
not the algebraic network elimination, is the obstruction.

## Inputs fixed for the first pilot

### Voltage regulator

Derive the converter's **incremental** passive output from the model before choosing a
regulator.  With the shifted plant Hamiltonian `H_s`, modulation input `u = m - m_star`, and
converter input vector field `g_m`, use

```text
y_m       = g_m' grad(H_s)
xi'       = y_m
u_command = -k_p y_m - k_i xi
m         = m_star + u_command
H_AVR     = k_i xi^2 / 2.
```

The ideal converter is power preserving for every `m`, so `g_m' grad(H) = 0` for the raw
physical energy.  The output above is a shifted-energy control output, not another physical
power port; it can be nonzero because `grad(H_s) = grad(H) - grad(H_star)`.  The base system
at `m_star` must first satisfy the corresponding shifted dissipation identity.  While the
modulation limit is inactive,
`y_m u_command + H_AVR' = -k_p y_m^2`; the AVR therefore adds damping without an
unaccounted supply.  Choose `m_star` from the desired terminal-voltage equilibrium.  Do not
replace `y_m` by `V_ref - |V_t|` unless the corresponding matching equation is derived and
checked.  The first certificate stays inside the modulation limits.  A later limiter proof
may use a projected integrator with anti-windup, but its dissipation inequality must be proved
on every admitted branch.

### DC source

Use a fixed DC source as a declared chemostatted boundary for gates 1--3, connected to the DC
link through a finite source resistance:

```text
C_dc v_dc' = (v_s - v_dc) / R_s - m i_f.
```

Here `v_s` is fixed, `R_s > 0` supplies explicit dissipation, and the attraction claim is
conditional on `v_s` remaining fixed.  This is preferable to an ideal voltage clamp because
the DC-link state and its energy remain in the model.  Add a finite source model only when a
study asks about source depletion or source-side controls.

### Parameter provenance

- Use the chosen IEEE-39 machine's electrical base, operating point and GENROU/GENSAL
  reactance and time-constant data.  Derive the full winding inductance and resistance data by
  an explicit map, then check positive definiteness, the open-circuit time constants and the
  quasi-steady reduction.
- Use an adjacent IEEE-39 branch for the pilot π line.  Convert its series and shunt stamp at
  nominal frequency to `R`, `L` and `C`, and verify that eliminating the electromagnetic states
  reproduces that stamp.
- Replace the selected bus demand by an RLC branch fitted to the same equilibrium `P` and `Q`.
  Record that this is a passive surrogate, not an identification of the original ComplexLoad.
- The IEEE-39/IEEE governor records do not identify penstock geometry, surge-tank area, steam
  volumes or thermodynamic states.  Take one complete hydro data set from Gil-Gonzalez et al.
  for gate 1 and one complete published steam-turbine data set cited by IEEE PES-TR1 for gate 2.
  Label both as benchmark parameters.  Match the existing RMS time constants only after the
  physical models pass their energy balances; do not infer missing energy coefficients from
  the RMS lags.

These choices make gate 1 reproducible without pretending that the IEEE-39 dynamic record
contains plant geometry that it does not contain.

## References

- IEEE Std 421.5-2016, *Recommended Practice for Excitation System Models for Power System
  Stability Studies*.
- IEEE PES-TR1 (2013), *Dynamic Models for Turbine-Governors in Power System Studies*,
  [PES Resource Center](https://resourcecenter.ieee-pes.org/publications/technical-reports/PESTR1.html).
- Gil-González, Garces, Fosso and Escobar-Mejía, "Passivity-based control of power systems
  considering hydro-turbine with surge tank", IEEE Trans. Power Systems 35 (2020) 2002-2011,
  [IEEE Xplore](https://ieeexplore.ieee.org/document/8877767/).
- Lohmayer, Kotyczka and Leyendecker, "Exergetic port-Hamiltonian systems: modelling basics",
  Math. Comput. Model. Dyn. Syst. 27 (2021) 489-521, [arXiv:2008.04091](https://arxiv.org/abs/2008.04091).
- Gernandt, Haller, Reis and van der Schaft, "Port-Hamiltonian formulation of nonlinear
  electrical circuits" (2021), [DOI](https://doi.org/10.1016/j.geomphys.2020.103959).
