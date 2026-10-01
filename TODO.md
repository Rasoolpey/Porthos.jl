# Porthos.jl progress and hand-off

Read this first in a new session, then `docs/ROADMAP.md` (its Part I banner: parity is
frozen) and `README.md`. Updated at the end of every session.

Last reviewed: 2026-10-01 (the storage completion is rejected as the primary physical
`H_ext` under rule 4 and retained as the benchmark `V_completion`; both closing diagnostics
are done: the decay-rate bisection and the backward block elimination, whose fixed-point
pattern keeps 31 of 40 units and 95 % of the repair entries; next is gate 1 of the
small-system physical extension, `docs/PHYSICAL_EXTENSION.md`).

## Start here (hand-off for the next conversation)

**Where things stand.** The simulator port is done and parity is frozen (pack v6 published).
The ROA pipeline (`src/roa/`) certifies the quadratic scaffold `V_P` on IEEE-39 (level
`1.714e-10`, `ROACheck` passes). Step 5a showed that PHPS's declared GENROU storage `Q`
leaves a non-exact stator exchange; `cdc429b` showed that GENROU admits the
Nishino-Chakrabortty-Ishizaki strain-energy structure with a metric `M`. **Step 5b (this
session, `src/ph/strain.jl`, `scripts/strain_energy_check.jl`, about 2.5 minutes):**
- the structure holds for all 11 rotors, the GENSAL included (`strain_rotor`,
  `strain_metric`; GENSAL's saturation enters as a field input and is inactive here);
- the identity `dU_ext/dt = -z'Mz' + Efd B_f'Mz' + sum Te delta' + residual` is exact to
  2e-14 along the flow (equilibrium, random states, the fault trajectory), with every
  residual term in closed form. With `dU_net/d delta = Te = Pe/omega` exactly. On the lossless
  variant the residual is only the speed voltage and the ComplexLoad filter lag; on the case
  network the conductance term dominates (up to 125 against a rotor dissipation of 26);
- the joint candidate `U_ext + sum omega_b H (omega - 1)^2` (torque-consistent kinetic
  storage) has a Bregman Hessian **positive semidefinite on every machine and network
  coordinate** (the 13 negative directions of 5a's `H_w` are gone); its only negative
  directions are the 16 ComplexLoad filter states, fixed by a small gauge `5 dz^2/2`;
- **local decay fails because of the network losses, not the machines.** Supply-shifted
  candidate `S'` and the exact decomposition of its rate into second-order terms
  (`strain_decay_forms`, closes to 4e-15): with the plant alone (controller states frozen)
  13 to 17 of 83 directions create energy (up to 7 to 15 /s). The creators are the
  equilibrium loss currents times the KCL-branch curvature (`loss_curvature`, up to 86 /s
  alone), the incremental conductance one-form (13 to 17 /s), the ComplexLoad active part
  (21 to 26 /s) and the load lag (energy-creating at every inductive load); the speed voltage
  is minor (0.9 /s), the COI coupling 3.5 /s. Details in 5b below.

**Method decision (2026-09-30): keep the present plant and complete the storage.** Do not
reclassify the active loads as unclosed external supplies: that would prove an open-system
dissipativity statement, not attraction of the closed IEEE-39 model. Keep the physical
strain energy as a named core and add the smallest exact, reference-invariant scalar repair
needed for positivity and decay. Solve the controller and loss problems together, in this
order:
1. Add quadratic curvature on the four feedback states of each IEEET1 (not its one-way
   reservoir) and fix the physical cross-term `-(Efd-Efd*) B_f^T M(z-z*)`. Determine the
   remaining exciter blocks/cross-terms from a structured local LMI. Treat each governor in
   the same joint machine-controller LMI; do not require a standalone passive governor,
   which the port tests have already ruled out.
2. With all controller states included, rerun the local decay decomposition. Only then use
   Route B to rank free network/load cross-terms and solve a constrained completion:
   preserve the Hessian of the physical core, allow sparse scalar corrections on the
   identified load/network/controller blocks, and require `P > 0`, `sym(PA) < 0` on the
   quotient. Check every selected block under other reference machines or the COI basis.
3. Lift the successful quadratic blocks to explicit nonlinear scalar terms and keep the
   report split as `H_ext = H_physical + H_controller + H_repair`. Certify the full scalar;
   do not call `H_repair` physical energy.
4. If the constrained completion is infeasible or needs a nearly dense repair, return to a
   physical model extension: dynamic network/bus energy states and passive load models.
   Dynamic π lines alone are not expected to fix this case because the dominant conductance
   is in load shunts. A load port is useful only when it is closed by a specified passive
   load subsystem and its storage; leaving its supply open does not establish an ROA.

**Completion result (2026-09-30; `src/ph/sink.jl`, `src/ph/completion.jl`,
`structured_completion` in the JuMP extension, `scripts/storage_completion.jl`, about 15
minutes, report `outputs/completion/storage_completion.txt`; tests in `test/unit/ph.jl`).**
- *Stage-1 plant (user's choice):* the plain `lossless_variant` has no equilibrium at the
  IEEE-39 dispatch (no active consumption, while Tm sums to about 58 pu). Stage 1 uses the
  **integrable-loss equivalent** `integrable_loss_variant(sys, x*, V*)`: G and the ComplexLoad
  active part removed, each bus's equilibrium loss power (36 buses, 61.4 pu) drawn by a
  certificate-only `ConstantPowerSink` with the exact potential `P theta + Q ln|V|`. Same
  equilibrium (3e-13), identity exact (7e-14), residual only speed voltage and load lag.
  `certificate_only` exempts the sinks from the contract lookup; the case model is unchanged.
- *Physical core there:* the Hessian of S' has one negative machine direction (-2.2; the
  constant-power sinks and the frozen-z reactive loads soften the network) and 18 load ones.
  A per-load filter-error term `kappa (|V| - Vini - z)^2 / 2` (kappa >= 5) removes the machine
  one; a diagonal `c (z - z*)^2 / 2` the load ones. These two are the load repair basis.
- *LMI* (`structured_completion`: `P = Hfix + F + sum k_j B_j`, max `t` with `P >= tI`,
  `-(A'P + PA) >= tI`, power-of-two Lyapunov scaling; Hypatia): with only the own-unit
  controller blocks (exciter 3x3 on xa, xe, xf, since xr is held with TR = 0; governor blocks;
  their cross-blocks with the machine), the field cross-terms fixed and the load terms:
  **infeasible**, t = -0.069. The binding dual is on omega and delta (the decay constraint
  binds, tr Z_R = 1); Route B ranks angle-speed couplings between machines first. Per-machine
  `eps_k Ddelta_k Domega_k` terms: -0.065; the free swing block: -0.037. A dual-guided greedy
  search (three component-pair blocks per round) reaches t > 0 after 15 rounds and 39 blocks
  (exciter-exciter between units, machine-machine, the big machine with exciters): interval-
  verified, t = 1.3e-4. Its reference-invariant closure (`invariant_closure`: 2806 free
  entries, 19 % of the symmetric matrix) gives t = 5.6e-3, interval-verified.
- *Stage 2, real plant:* the same closed pattern certifies immediately, **no extra
  loss-specific block needed**: t = 1.19e-3, interval `P_min >= 0.0039..0.0045`,
  `decay_min >= 0.0012`, cond(P) 3.7e5 to 4.2e5.
- *Lift:* `storage_completion` / `completion_energy` give the explicit nonlinear
  `H_ext = H_physical + H_controller + H_repair` (S' with the field cross-terms; the
  controller quadratic; inter-unit quadratics and the load terms through `|V(x)|`); its
  Hessian equals the LMI matrix to 8e-16, gradient 4e-13.
- *Size of the repair (the rule-4 question):* 2070 repair entries (14 % of the matrix) after
  the invariant closure; a group-sparsity solve (SOC per component pair, 263 groups, hit its
  20-minute cap at a feasible iterate) leaves 202 groups above 1e-2 of the largest, the
  largest being the machines' **own** blocks (the rotor physical Hessian is reshaped, not
  only coupled); pruning at 1e-3 / 1e-2 re-closes to the same pattern. In scaled Frobenius
  norm (max-margin solution) physical 681, controller 1.7e3, repair 1.75e3. So: feasible and
  interval-certified locally, but neither sparse nor small, the margin is small (t ~ 1e-3)
  and cond(P) ~ 4e5. The repair was already needed on the integrable-loss plant, so the
  obstruction is the machine-controller-network coupling (AVRs sensing |V|, the swing
  couplings), not the losses.

**Rule-4 decision (2026-09-30): reject this completion as the primary `H_ext`.** Keep it as
an exact local Lyapunov/existence witness and a diagnostic benchmark, but do not spend the
ROA-certificate effort on it yet. The decisive evidence is structural rather than the raw
14 % entry count: 202 of 263 component-pair groups remain significant, the free repair
reshapes the machines' own Hessian blocks, both the controller and repair norms exceed the
physical core, and the verified margin is small with condition number about `4e5`. It no
longer supports the intended statement "physical strain energy plus a modest correction."

Before changing the plant, close the negative result with two capped diagnostics, not an
open-ended sparsity run:
1. Compute the generalized decay rate of the present solution and maximize a *normalized*
   rate `mu` by bisection on `-(A'P+PA) >= 2mu P`, with explicit condition-number caps in the
   scaled coordinates. The present absolute `t` is scaling dependent.
2. Starting from the feasible pattern, use block backward elimination at a fixed useful
   `mu` and condition cap. Record an infeasible dual certificate when removing a group fails.
   Stop when a complete pass removes nothing; this decides indispensability more directly
   than waiting for the group-SOC objective to converge.

**Diagnostics so far (2026-09-30; `completion_rate_feasibility` and
`generalized_decay_rate`, the benchmark is named `V_completion`).**
- *Normalised decay rate* `mu` (largest with `-(A'P+PA) >= 2 mu P`, coordinate-invariant):
  the max-margin `V_completion` has only `mu = 5.2e-4 /s` (cond 3.7e5). Bisection over the same
  pattern (geometric on [5e-4, 0.027]; 0.027 = -max Re eig(A), the IEEEG3 hydro-governor mode,
  bounds every quadratic Lyapunov function): uncapped `mu* in [0.0120, 0.0128]` (the upper
  end certified: converged negative margin), cond 6.3e5; cap `cond <= 1e6`: `mu >= 0.0120`
  (cond 1.4e5); cap `1e5`: `mu >= 0.0113` (cond 1.0e5). For the capped runs the upper ends are
  time limits (900 s), not certificates. So the pattern supports about 44 % of the spectral
  bound with cond about 1e5.
- *Backward block elimination* at `mu = 0.006`, cap `1e5` (40 units: the 39 recorded pairs and
  the swing block; each test re-closes the pattern). A test is: the uncapped max-margin
  problem at rate `mu` (a converged `t < 0` gives a saved solver-level infeasibility
  certificate with its duals, even without a cap); if feasible but worse conditioned than
  the cap, the capped feasibility problem decides; "undecided" means feasible below the cap
  only at a rate just under `0.99 mu`, or the capped solve hit its time limit. All verdicts
  are conditional on
  this pattern, `mu = 0.006` and `cond(P) <= 1e5`; they are Float64 solver results (the
  saved duals are not yet interval-checked). Feasibility is monotone in the allowed blocks,
  so a block indispensable against an accepted set stays indispensable when more blocks are
  removed.
  **Phase 1, single removals (done, 40 of 40):** 8 blocks indispensable with every other
  block present: `IEEET1_10 x IEEET1_11`, `GENROU_1 x IEEET1_10`, `GENROU_1 x GENROU_11`,
  `GENROU_11 x GENROU_9`, `GENROU_6 x GENROU_7`, `GENROU_5 x GENROU_7`,
  `GENROU_10 x GENROU_10` and the swing block; 31 removable alone; 1 undecided
  (`GENROU_8 x GENROU_8`, rate 0.005936). Each removal takes out only 9 to 24 entries after
  re-closure.
  **Phase 2, batched removal of the 32 candidates (done, 2026-10-01; one worker, about 15
  minutes per solve):** removing all 32 together fails (666 entries), as do both halves; the
  blocks substitute for one another. **Jointly removable, with rate and cap kept: 9 blocks**
  (`GENROU_10 x GENROU_11`, `IEEET1_6 x IEEET1_9`, `IEEET1_5 x IEEET1_9`,
  `GENROU_4 x GENROU_7`, `IEEET1_5 x IEEET1_6`, `IEEET1_7 x IEEET1_8`, `IEEET1_4 x IEEET1_5`,
  `GENROU_11 x IEEET1_9`, `IEEET1_10 x IEEET1_4`). After the invariant re-closure they remove
  only **109 of the 2070 repair entries (5 %)**: the other blocks re-cover most of their
  angle entries.
  **Phase 3, final pass (done):** every block neither accepted nor already indispensable
  against an accepted subset was retested alone against the final 9 (capped solves allowed
  1200 s); the pass removed nothing, so the procedure is at its fixed point. Against the
  final accepted set: **indispensable, the 8 of phase 1 plus 21 more** (`IEEET1_4 x IEEET1_6`,
  `GENROU_11 x GENROU_8`, `IEEET1_4 x IEEET1_8`, `IEEET1_10 x IEEET1_9`,
  `GENROU_10 x GENROU_9`, `GENROU_1 x IEEET1_5`, `GENROU_1 x IEEET1_6`,
  `GENROU_11 x GENROU_2`, `IEEET1_10 x IEEET1_7`, `GENROU_10 x GENROU_8`,
  `GENROU_2 x GENROU_3`, `GENROU_11 x GENROU_3`, `GENROU_10 x GENROU_4`,
  `GENROU_4 x GENROU_8`, `GENROU_11 x GENROU_7`, `IEEET1_10 x IEEET1_2`,
  `IEEET1_2 x IEEET1_5`, `GENROU_8 x GENROU_8`, `GENROU_3 x GENROU_4`, `IEEET1_4 x IEEET1_7`,
  `GENROU_7 x GENROU_9`; the last four were undecided against fewer accepted blocks);
  **undecided, kept in the pattern: 2** (`GENROU_3 x GENROU_8`, `GENROU_10 x GENROU_3`:
  feasible uncapped only above the cap; the capped solve hits its time limit).
  **Result:** the fixed-point pattern keeps 31 of the 40 units, 2697 free entries
  (736 own-unit controller, **1961 repair, 13.3 % of the symmetric matrix**): 29
  solver-indispensable blocks and 2 unresolved ones. It is a fixed point under the
  procedure's acceptance rules and time limits, for this ordering; strict irreducibility is
  not established for the 2 unresolved blocks, and it is not necessarily a globally sparsest
  pattern.
  It keeps inter-unit exciter-exciter, machine-machine and big-machine-exciter couplings,
  machines' own blocks and the swing block; every indispensable verdict comes with a
  converged negative margin. This supports the rule-4 conclusion: broad mathematical
  coupling is required by this reduced linear model and certificate class at this rate and
  conditioning, rather than being introduced solely by the greedy search. It does not make
  those blocks physical energy or prove that the plant has a corresponding two-component
  energy store. The study is closed; the method is in `completion_rate_feasibility`,
  `generalized_decay_rate` and `invariant_closure` (the driver scripts were session
  scratch work).

Then move to a physical extension, but treat it as a machine-network-load-controller model,
not a π-line-only change. Retain dynamic line/bus electromagnetic energy, close shunt and
ComplexLoad power through specified passive load dynamics, and give the exciter field and
governor actuator/steam paths physical two-way storage ports. Re-run the own-unit passivity
and compositional-storage gates before an IEEE-39 model is built. Stage 1 already shows that
network lifting alone may leave the AVR/swing inter-unit obstruction.

**Physical-extension design decision (2026-09-30; one-page design note in
`docs/PHYSICAL_EXTENSION.md`).** Build this as a separate balanced
electromagnetic dq model; do not modify the frozen RMS case models. The IEEE DC1A, IEEEG1 and
IEEEG3 blocks are reduced stability-study transfer models: their time constants do not by
themselves identify inductances, capacitances, moving masses or thermodynamic storage. Do not
attach guessed physical Hamiltonians to their existing lag states.
1. *Excitation:* the generator field-winding magnetic energy belongs to the full machine
   Hamiltonian and must not be counted twice. Close its physical `(v_f, i_f)` port through an
   averaged ideal converter and a DC-link capacitor,
   `H_dc = C_dc v_dc^2/2`, with the modulation interconnection power preserving and field/DC
   resistances explicit dissipation. Close the source side with an explicit supply model or
   declare a fixed-voltage source as a chemostatted boundary and make the resulting ROA
   conditional on it; never leave an unspecified power supply. Use a simple passive AVR for
   the first structural test.
   A rotating dc-commutator exciter is a later DC1A-specific variant only if its winding,
   armature-reaction and shaft parameters are supplied; the DC1A `xe, TE, KE` block alone is
   not enough to reconstruct that energy.
2. *Prime mover:* make two pilots. For IEEEG3/hydro, use gate-servo hydraulic compliance,
   penstock water-column kinetic energy and surge/head-tank gravitational/compliance energy;
   start here because its mode is the spectral bottleneck. For IEEEG1/steam, use servo
   hydraulic storage plus steam-chest/reheater **availability (exergy)** in pressure/mass
   states, with enthalpy-flow ports and turbine shaft power. The existing lag states may be
   matched as a reduced limit, not declared to be those energies without a parameter map.
   Treat upstream head and boiler chemical potential as explicit chemostatted boundaries (or
   model their finite reservoirs); state that condition in every attraction claim.
3. *Electrical test system:* full winding/stator-flux synchronous machine, dynamic dq π line
   (series inductor and shunt capacitor states), and a passive RLC load. This is a stiff
   averaged electromagnetic model, not a switching EMT model. Prove the component power
   balances and the total Dirac interconnection first; then show that setting the fast
   derivatives to zero recovers the RMS stator, line and load stamps.
4. Gate order: one hydro unit to an infinite bus; one steam unit; a two-machine steam/hydro
   system; only then an IEEE-39 extension. At each gate require a convex shifted Hamiltonian,
   local decay with own-unit or nearest-neighbour terms, and the singular-limit comparison.
   If broad cross-unit completion returns in the two-machine model, record that the legacy
   controllers rather than algebraic network elimination are the obstruction.
5. Inputs for the first pilot are fixed in `docs/PHYSICAL_EXTENSION.md`. Derive the passive
   converter output `y_m = g_m' grad(H_s)` and use passive PI damping
   (`xi' = y_m`, `m-m_star = -k_p y_m-k_i xi`) while the modulation limit is inactive; set
   `m_star` from the voltage target rather than assuming voltage error is power conjugate.
   Feed the DC link from a fixed, declared DC chemostat through a finite source resistance.
   Use IEEE-39 data for the electrical operating point and an adjacent line/RLC surrogate,
   but use complete published benchmark data for hydro and steam storage parameters that the
   RMS records do not identify. Every map and quasi-steady reduction is a gate, not an
   assumption.

**Upcoming steps (agreed 2026-10-01).**
1. ~~Finish the elimination study, record it, commit and push~~ (done 2026-10-01).
2. Next session: gate 1, one hydro unit to an infinite bus in the electromagnetic dq model
   (`docs/PHYSICAL_EXTENSION.md`). Prove, in this order:
   1. the equilibrium closes with `m = m_star`, the fixed DC boundary and `R_s > 0`;
   2. the raw physical-energy balance closes exactly, including the source power and every
      resistive and hydraulic loss;
   3. the shifted unregulated identity has no unexplained residual and identifies its signed
      dissipation. **Stop here if it has a positive residual outside the span of `y_m`:** the
      proposed AVR cannot repair that through the modulation alone;
   4. for `u = m - m_star`: `dH_s/dt(m_star + u) - dH_s/dt(m_star) = u y_m`, `y_m = g_m' grad H_s`;
   5. with the PI controller: `dH_s/dt + dH_AVR/dt = dH_s/dt|_{m_star} - k_p y_m^2`;
   6. the certified region stays strictly inside the modulation limits (limiter branches
      afterwards).
3. Then gate 2 (one steam unit), gate 3 (steam and hydro, two machines), and IEEE-39 only
   after the first three pass (`docs/PHYSICAL_EXTENSION.md`, "Gates").

Reuse: `strain_balance`, `strain_decay_forms`, `completion_problem`, `completion_pattern`,
`structured_completion` (max margin, L1 or group sparsity), `storage_completion`,
`quotient_hessian`, `integrable_loss_variant`. Probe small and bound long runs (a quotient
Hessian takes about 50 s, a completion LMI 2 to 6 minutes, a sparsity solve over 20).

**Working rules to keep:** commits only when the user asks, with the owner's Git identity and
no AI attribution trailer (`Co-Authored-By`, `Signed-off-by` or similar); do not change
component dynamics, runtime outputs or the IEEET1 reservoir wiring; label certificate-only
terms as such; method choices on controllers go to the user first. CI results are not a
priority for the user. Standing project rules: Julia does all computing (Python only edits
JSON, runs `parity/generate/` and drives PowerFactory in `pf/`); the PHPS reference at
`C:\Users\em18736\Documents\PHPS_Opt` is read-only; PowerFactory runs only through `pf/`
with its window closed; component code stays generic in the number type, uses the branch
primitives and an allocation-free `rhs!`; certificates are proofs, never samples; the
repository is public (no credentials, unpublished reviewer correspondence or data dumps).

## Goal

A **rigorous region of attraction (ROA)** for the power system, built on a **storage
function that truly decreases**, ideally the physical port-Hamiltonian energy extended by
well-founded terms. Two layers (decided 2026-10-01):

1. **Now:** certify the ROA with the local quadratic candidate `V_P` (the Lyapunov equation
   on the common-angle section). It exists, covers every physical state, and is the test bed
   for the rigorous interval machinery (P11).
2. **Target:** a nonlinear composite storage
   `V_ext = H_machines,s + U_net,s + V_controllers + V_cross`, certified by the same pipeline
   and then the primary ROA certificate. Call it `H_ext` only where every added term has a
   physical or port-Hamiltonian derivation; otherwise call it `V_ext` and say which parts are
   physical and which mathematical. Do not treat a `V_P` certificate as the answer to the
   physical-storage question.

The user's standing preference: keep the component dynamics unchanged if at all possible
(storage-only route: new storage terms, joint storages, other supply rates or ports);
physical model changes (two-way reservoirs) are the fallback. Control-method choices are
discussed with the user before they are implemented. Commit attribution stays with the
repository owner's configured identity: never add AI `Co-Authored-By` or similar trailers.

## Decisions (2026-10-01, user with a reviewer)

- **Candidate-independent certificate pipeline.** A Lyapunov candidate provides: its value,
  gradient, interval evaluation, equilibrium shift, physical-coordinate projection and a
  fingerprint. First implementation `QuadraticCandidate(P)`; second
  `ExtendedStorageCandidate(H_ext)`; both go through the same proof gates.
- **Network: hybrid.** Keep the bus voltages and the network explicit when deriving storage,
  port cancellation and physical meaning. In the proof, validate a unique KCL branch and
  eliminate the voltages by the interval implicit-function argument. The dense Schur
  complement `f_x - f_V g_V^-1 g_x` must not be mistaken for physical cross-machine storage.
- **Limiters: single mode first.** The certified ROA must stay inside one decided branch of
  every limiter and domain clause (the branch primitives raise `UndecidedBranch` otherwise).
  Multi-mode or boundary certificates only once the smooth-region result is useful.
- **Parity is frozen** as a regression guardrail (pack v6); PHPS is not ground truth; no
  more parity refinement (see the roadmap banner). PHPS's certificate numbers (the P11 gate)
  are sanity references only.
- Route B (the structure-search LMI) is a **diagnostic** for missing storage blocks, not a
  source of physical terms.

## Next steps, in order

1. ~~Record the certificate decision~~ and ~~add the roadmap override~~ (done 2026-10-01).
2. ~~**Publish pack v6**~~ (done 2026-09-30: release `parity-pack-v6`, the downloaded file's
   sha256 and `Tar.tree_hash` match `Artifacts.toml`). No more parity-pack work.
3. ~~**Generic Lyapunov-candidate interface with the quadratic `V_P`**~~ (done 2026-09-30,
   commit `cb2a777`; `V_P` stays scaffolding that validates the machinery, the main effort
   moves to steps 4 to 6;
   see "ROA certificate pipeline" below). `V_P` certified at `1.71e-10` on IEEE-39 with every
   gate; `ROACheck` passes. Open inside this step, none blocking step 4:
   - speed: the coordinate-by-coordinate centered hull takes about 25 s per level (171
     nested-dual Jacobians over the whole model); structure-aware sparse second-order jets
     would cut it, and larger inner chunks compile far too slowly (tried: chunk 16 did not
     finish in 9 minutes);
   - the certified set is tiny in physical terms (box half-widths up to about 2e-4 in scaled
     coordinates); growing it is Target A (roadmap II.1: conditioning or hull-aware `P`,
     then fault reach), not part of the Target B track;
   - the `H_ext` candidate needs its own `sublevel_half_widths` (a quadratic lower bound),
     `gradient_matrix_hull` (an interval Hessian of `H_ext`) and `positivity_proof`; the
     gates themselves are shared.
   Review closures before the record is final (review of 2026-09-30):
   - ~~claim~~ (done): the record now states attraction of the retained physical quotient,
     conditional on the excluded states (reservoirs, delta_COI) staying in the ranges over
     which they are proved not to feed back (widest `x0 * [2^-k, 2^k]`, k in 60/40/20/10;
     the reservoirs get `2^10`, the next span straddles their own `guard_min`), with a bound
     on each one's drift rate over the certified box; convergence of the reservoirs and a
     full-state ROA are explicitly not claimed. Global nonfeedback (for every reservoir
     value) is not proved;
   - ~~rotation guard~~ (done): `rotation_action` per model type (machines: `delta` and the
     network-frame V / I pairs; ComplexLoad: its V / I pair; exciters and governors:
     invariant) and `check_rotation_symmetry` on the wiring, injections and fixed-voltage
     buses; model types without a declaration (the GFM converters, VOC's Cartesian states
     in particular) are rejected;
   - ~~analytic tests~~ (done): `AnalyticModel`; 1-D `x' = -x + x^3` (all three hulls
     certify exactly `c < 1/6`: pass 0.16, fail 0.17 and at the true boundary 0.5), 2-D
     radial `x' = -x (1 - |x|^2)` (pass 0.05, fail 0.5), a clamp (certified inside one mode,
     `UndecidedBranch` across it), an unstable case rejected. They exposed and fixed two
     bugs: the centre of the centered hull now needs the equilibrium voltages inside the
     branch's uniqueness box `X` and is evaluated on `eq.V cap V`, and `verified_min_eig`
     failed on 1 x 1 matrices;
   - ~~README status~~ (done); full suite with the final default: see the loose ends;
   - ~~clean record~~ (done): `outputs/roa/IEEE39Bus_PF/certificate_quadratic.json` and
     `roa_check.json` regenerated from commit `cb2a777` (`porthos_src_modified = false`):
     level `1.714e-10`, decay bound `-3.44e-3`, 83 clauses, 21 excluded states with ranges and
     drift bounds, ROACheck passed. `outputs/` is not in git; rerun the script to reproduce.
4. ~~**Nonlinear port-power residual audit**~~ (done 2026-09-30; `src/ph/power.jl`,
   `scripts/port_power_audit.jl`, about 10 s, report `outputs/ph_audit/<case>/port_power.json`;
   tests in `test/unit/ph.jl`). At the equilibrium, 16 random KCL-consistent states and 601
   points of the bus-16 fault trajectory (IDA):
   - identities: network balance (injected = `V'GV` + fault + frequency loads + `V . g`)
     2.8e-13, machine terminal (`V . It = V . I_norton - Re(y_n)|V|^2`) 3e-14, `Pe` = terminal
     power 7e-15, ComplexLoad declared P = drawn power 7e-11;
   - one-way reservoirs (IEEET1, IEEEG1, IEEEG3): residual exactly 0 (lossless accounts);
   - machines, against the contract ports `Tm`, `Efd*i_fd`, `-(Vd*Id+Vq*Iq)`: with the
     declared shifted kinetic storage `H (omega-1)^2` every unit "creates" energy (residual
     down to -19.9 in the fault), because `Tm` is conjugate to `H omega^2`, not to the shifted
     form; with `H omega^2` the mechanical side (`Tm - Pe - d(H omega^2)/dt`, D = 0) and the
     terminal side close to 1e-14, and **the whole residual is in the magnetic block**
     (`Efd i_fd - dH_mag/dt`): GENSAL is dissipative everywhere ([0.0022, 0.021]); **GENROU is
     not** ([-0.022, 0.054]; negative for GENROU_2, 4, 7, 8, 9, near the equilibrium at 1e-2
     perturbations as well as in the transient). So GENROU's declared circuit storage, with
     the field port as its only electrical supply, is not a storage function: its rotor
     circuits exchange energy with the stator through a path no port accounts for (the
     contract's own status is `M1 rotor_field_pass_terminal_network_open`). This needs a
     method discussion before step 6: which storage or port closes GENROU's magnetic block
     (e.g. an explicit stator-exchange port, or a different magnetic storage in the same
     coordinates); GENSAL shows the structure that works.
   - ComplexLoad's contract port is prose, not an expression, so its port is not evaluated
     (its power identity is checked on the network side).
   ~~Review before pushing~~ (done 2026-09-30): the trajectory voltages are now re-solved on
   the healthy or fault-on KCL branch (`solve_network(...; faults_on)`), so every sample
   satisfies KCL to round-off (`|V . g| <= 3e-13`, was `1.14e-3` from IDA's tolerance). It
   does not change GENROU's `-2.18e-2` magnetic residual, which also appears in the
   KCL-consistent near-equilibrium samples.
   **GENROU method decision:** keep the dynamics unchanged and first derive the exact
   quadratic rotor-circuit identity. With `z = (Eq', psi_d, Ed', psi_q)`, write
   `zdot = A z + B_f Efd + B_s [id, iq]` and
   `H_mag = z' Q z / 2`; the existing `i_fd` already checks the field colocation
   `B_f' Q z`. Declare the missing bilinear term `[id, iq]' B_s' Q z` as an internal
   stator-exchange port and prove that the remaining quadratic term is nonnegative rotor
   loss. In step 5, require the opposite stator-exchange supply to appear in the explicit
   stator/network energy identity, so it cancels on interconnection. Only if the completed
   identity leaves an indefinite rotor-loss term should a different `Q` be sought, constrained
   by the same field and stator port colocation. Keep raw physical power (`H omega^2`, `Tm`)
   separate from the equilibrium-shifted/Bregman storage and its incremental supply used for
   Lyapunov decay.
   **Result (2026-09-30): the route works with the existing `Q`.** `rotor_structure` (in
   `genrou.jl`) restates the rotor as `z' = A z + B_f Efd + B_s [id, iq]`; the audit checks it
   against the right-hand side at every sample (5e-15), the field colocation
   `i_fd = B_f' Q z` (2e-17) and the balance `dH_mag/dt = Efd i_fd + [id, iq]' B_s' Q z -
   rotor loss` (2.5e-16). The rotor loss matrix `-sym(QA)` is proved positive definite for
   every GENROU unit (rigorous `lambda_min >= 4.5e-4`), so the rotor loss is nonnegative
   (sampled [0.011, 0.24]) and the negative magnetic residual is exactly the stator exchange
   (sign-indefinite, [-0.029, 0.21]): no other `Q` is needed. Step 5 tests how much of this
   exchange can be represented by the stator/network energy identity; the result is in 5a.
5a. **Step 5.1 finding and decision (2026-09-30; the coordinate decision in `374d9b5`, the
   exactness gate and the review wording below in the commit after it): the
   reviewer's route cannot close in the declared coordinates.** At the equilibrium the stator exchange
   `p_s = [id, iq]' B_s' Q z` is nonzero (0.02 to 0.10 per GENROU unit, 1 to 20 times the
   field power) while the mechanical and terminal balances are exact (`Tm = Pe`), and every
   storage has zero rate there: no stator/network storage or lossless port can supply it. It
   is not caused by PHPS's simplified `dEq'/dt` (no Sauer-Pai correction; PHPS
   `genrou_phs.py` line 316): with the Sauer-Pai equation `p_s` at steady state is the same.
   It is not the neglected transformer power `i . dpsi''/dt` either (correlation -0.6 to 0.4).
   Cause: `H_mag(z)` is the magnetic energy at zero stator current, so `Qz` is not the rotor
   current when stator current flows (the dampers show "current" and loss in steady state).
   **Rotor-current coordinates close it exactly, with the same `Q` and dynamics:**
   `w = z - D i`, `D = -A^-1 B_s`, so `z' = A w + B_f Efd`. At the equilibrium the damper
   components of `Qw` are 0 (1e-18) and field power `Efd B_f' Q w` equals the loss
   `-w' sym(QA) w` (1e-15), on all 10 units; the loss matrix is the one already proved
   positive definite. Then `d/dt (w'Qw/2) = Efd B_f'Qw - w'(-sym QA)w - w'QD di/dt`: the
   stator exchange becomes a transformer-type port `(di/dt, -D'Qw)` that vanishes at every
   equilibrium. **Decision:** adopt `w` as the current-corrected rotor coordinate of a
   *joint machine-network candidate*; `Qw`, rather than `w`, is the rotor-current/coenergy
   variable. Do not replace the component-local `hamiltonian(c, x)`, which remains the
   open-stator intrinsic audit, because the new storage depends on the algebraic current.
   Since `i` itself depends on the states and KCL voltages, call `w` a derived variable until
   the reduced map has been proved locally full-rank; positivity of `Q` alone does not prove
   positivity of the equilibrium-shifted storage after this composition.
   Define the new conjugate flow `i_fd_energy = B_f'Qw` in the energy/audit API only. Keep
   the existing runtime `i_fd = B_f'Qz` and the IEEET1 one-way reservoir unchanged for the
   frozen model; that legacy account is excluded from `V_ext` and must not be presented as
   the field balance of the new joint storage. Replacing or rewiring it would require a new
   model version after the hardware field-current base is validated.

   Continue step 5 through the KCL branch: compute
   `di/dt = (i_x - i_V g_V^-1 g_x) xdot` (or its section-coordinate equivalent), first for
   the lossless internal-node network and one smooth mode. Do not assume a scalar `U_net`
   exists: test whether the complementary one-form is exact (reference-angle invariant and
   symmetric Jacobian/curl zero). If it is exact, integrate it and verify cancellation of
   `-w'QD di/dt`; if not, retain the term as a physically motivated exchange/passivity
   shortage rather than calling it network energy or scalar storage. Derive the polar
   conjugate variables only after this identity. Build the equilibrium-shifted/Bregman form
   of the resulting joint storage on the certified KCL branch before using it as a Lyapunov
   candidate.
   **Result of the exactness gate (2026-09-30; `src/ph/exchange.jl`,
   `scripts/exchange_one_form.jl`, about 30 s; tests in `test/unit/ph.jl`): the one-form is
   not exact, so no scalar `U_net` cancels the exchange.**
   - `w = z - D i` verified on all 10 GENROU units at the equilibrium (damper currents 1e-17,
     field power = loss 4e-15); `i_fd_energy = B_f'Qw` differs from the runtime `i_fd` by 1
     to 55 % (the runtime output and the IEEET1 reservoir are unchanged).
   - `di/dt` through the KCL branch: `kcl_solve` iterates with the Float64 `g_V` at the
     converged point (zero contraction on the dual parts), so duals carry the exact
     implicit-function derivatives; they match finite differences to 1.4e-8.
   - Reference angle: the stator currents are invariant under a common rotation (1e-14), so
     the one-form `a . d eta`, `a = (dI/d eta)' D'Qw`, is well defined on the quotient.
   - Exactness: the Jacobian of `a` is not symmetric: relative curl 0.109 on the case
     network, 0.113 on the lossless variant (`lossless_variant`: every conductance, active
     load and ComplexLoad `P0` removed); antisymmetric/symmetric part 6.8 % (Frobenius);
     path integrals over a 1e-2 path differ by 2e-6 of 5e-3. The curl has a closed form,
     `(dI/d eta)' D'Q (dz/d eta) - transpose` (matches the measured Jacobian to 1.4e-16):
     writing `a . d eta = -d(I'D'QD I / 2) + (D'Qz)' dI`, the second part is not exact
     because the rotor states and the stator currents vary independently on the KCL branch.
     This is intrinsic to that one-form for the declared `Q`, not to network losses; the
     strain-energy metric found below changes the storage and removes this armature-reaction
     obstruction.
   - **Decision at this stage, superseded for the physical candidate below:** if the declared
     `Q` were retained, its non-exact one-form would remain an explicit, sign-indefinite
     shortage in `V_ext`'s decay identity. Do not call the one-form `V_cross`: it is not a
     scalar storage. Splitting off the exact scalar `F = I'D'QD I / 2` gives the gauge family
     `H_alpha = H_w - alpha F`; changing `alpha` redistributes an exact differential but
     cannot change the curl. Use `alpha = 0` (the raw nonnegative rotor energy `H_w`) as the
     baseline, and test `alpha = 1` and any later structured scalar cross-terms only through
     the equilibrium-shifted/Bregman Hessian on the quotient
     and the decay calculation. In particular, `H_w - F = z'Qz/2 - z'QD I` is not positive
     merely because `F` is exact. A coordinate change, including polar coordinates, cannot
     remove the curl; only an independently derived physical supply may cancel it. If the
     shortage is not dominated by proved losses, use structured exact scalar cross-terms
     (guided by Route B) or a dynamic extension, and state which construction was used.
   **Calculations 1 to 3 (2026-09-30; `src/ph/joint_storage.jl`,
   `scripts/joint_storage_check.jl`, about 4 minutes; fast tests in `test/unit/ph.jl`).**
   On the KCL branch `H_alpha` is a function of the section coordinates through
   `w = z - D I(eta)`, so its Bregman Hessian at the equilibrium is
   `W'QW + sum_j (Qw*)_j Hess w_j`, with the field current `(Qw*)_fd` nonzero.
   - Quotient Hessians (171 coordinates; inertia +/-/0): `H_w` (alpha = 0) 69/13/89, most
     negative -0.083; `H_w - F` (alpha = 1) 68/14/89, -0.053; with the shifted kinetic
     storage `sum H (omega - 1)^2`: 71/11/89 (-0.072) and 78/4/89 (-0.045). **Every form is
     indefinite**: the stator-current curvature makes the composed storage non-convex, and
     alpha = 1 reduces but does not remove the negative directions. The 89 zero directions
     are what the machine storage does not see (controllers, network).
   - Exchange against dissipation, at linear order around the equilibrium (quadratic forms
     on the section: rotor loss `W'LW`, exchange `-sym(W'QD J_I A)`, incremental network
     conductance loss `J_V'G J_V`): the exchange vanishes where rotor + network dissipation
     vanish (1.8e-12), but on the rest **it exceeds that dissipation by up to 927 times**
     (largest generalized eigenvalue; scale-invariant). Proved rotor and network losses do
     not dominate the shortage locally **for the declared `Q`**. This comparison is obsolete
     for the physical strain-energy candidate and must not be carried over to `M`. Controllers
     carry no storage yet, so their dissipation was not in the comparison.
   - **Literature check and revised candidate contribution:**
     `docs/LITERATURE_FINDINGS.md` records 30 core papers and six cross-field references. The
     non-closed one-form is a valid diagnosis of PHPS's `Q`, but the reproduced
     Nishino-Chakrabortty-Ishizaki construction shows that GENROU admits a different convex
     strain energy. The candidate contribution is now the four-state GENROU extension, the
     precise residuals of the detailed model, and an interval-certified joint `H_ext`, subject
     to a broader novelty search.
   - **Revised method decision (2026-09-30): physical strain energy first; structured scalar
     repair second; dynamic extension last.** Verify `U_rot + U_B`, then close or bound the
     speed-voltage, field-controller, load and conductance terms. Use Route B only to select
     sparse, rotation-invariant exact scalar terms for residual curvature or decay. Require a
     positive Bregman Hessian and useful local decay margin before the nonlinear interval
     proof. If no static candidate survives, test Krasovskii/Brayton-Moser rate storage and
     only then a certificate-only dynamic-supply/IQC extension.
5. **Polar ports and the network balance** (Route A'). **Step 1 done (2026-09-30;
   `src/ph/polar.jl`, `scripts/polar_balance.jl`, about 1 minute; tests in `test/unit/ph.jl`)**,
   at the equilibrium, 8 random KCL-consistent states and the healthy part of the bus-16
   fault (70 samples):
   - **Network, lossless part: exact.** `U_B = -1/2 sum B_ij (Vd_i Vd_j + Vq_i Vq_j)` has
     `dU_B/dt = sum_i (P_i^B theta_i' + Q_i^B d ln|V_i|/dt)` (2.8e-13 against rates up to
     210): the conjugate pairs derived from the identity are (P, theta) and (Q, ln|V|).
   - **Conductance part: not exact.** `sum (P^G theta' + Q^G d ln|V|/dt) = -Im(V'^H G V)`
     (1e-13); its curl in Cartesian voltages is `2 (G kron J)`, large on IEEE-39 (`G` up to
     65, mostly the constant-impedance loads): transfer and shunt conductances make the
     network balance path-dependent, as known for structure-preserving energy functions.
   - **Machines, at the internal EMF node** `E'' = omega (-psi_q'', psi_d'')` behind `j x''`
     (ra = 0, xd'' = xq''): the network delivers `(P_int phi' + Q_int d ln|E''|/dt)/omega_b`,
     which is exactly `omega` times the transformer power `i . dpsi''/dt / omega_b` plus
     `omega' i . psi''/omega_b` (1.4e-17). **It does not match the rotor exchange
     `w'QD di/dt`**: the exchange is 5 to 50 times larger (rms 0.05 to 0.10 against 0.002
     to 0.011) and only weakly correlated (0.16 to 0.59), so the mismatch is essentially the
     whole exchange. This rules out cancellation of the exchange generated by PHPS's declared
     magnetic `Q`; the two-axis result below shows that changing to the network-compatible
     rotor strain energy makes armature reaction an exact gradient instead.
   **Two-axis reproduction, first result (2026-09-30; `rotor_gradient_metric` in
   `src/ph/joint_storage.jl`, test in `test/unit/ph.jl`; commit `cdc429b`).** Source: the arXiv
   version, Ishizaki, Nishino, Chakrabortty, arXiv:2304.00987v2 (the TAC 2026 paper's
   preprint). Their mechanism: on a lossless network the two-axis flux dynamics are a
   *gradient flow* of the strain energy, `tau E' = -(X - X') dU/dE + V_fd`, with
   `P = dU/d delta`; the Bregman form of `U` is the EI storage, the flux terms give the
   dissipation `-tau |E'|^2/(X - X')`, and losslessness is needed for the gradient to exist.
   The GENROU question is then algebraic: per axis, a symmetric metric `M` with
   `M B_s = -Psi''^T` (the sign fixed by the internal-EMF network potential,
   `grad_z U_net = Psi''^T i` for `E'' = (-psi_q'', psi_d'')` behind `j x''`) and `M A`
   symmetric. **It exists on all 10 GENROU units with `M > 0` and `U_rot = -z^T M A z / 2`
   strictly convex** (smallest eigenvalues of `M` 0.21 to 11, of `-sym(MA)` 0.87 to 39;
   `M B_s + Psi''^T` at 1e-16), and every machine has `ra = 0`, `xd'' = xq''`, `D = 0`. So
   GENROU does admit the strain-energy structure; the non-exact exchange of 5a came from
   PHPS's declared `Q`, which differs from `U_rot` (shape 0.7 to 20 %). What does not carry
   over: (1) the speed voltage (`E'' = omega psi''` and `Te = Pe/omega` in the model, no
   `omega` in theirs), a sign-indefinite `(omega - 1)` residual; (2) conductances (their
   Theorem 1); (3) the ComplexLoads' voltage-dependent corrections, which need their own
   potential; (4) the field port becomes `(Efd, B_f^T M z')` (a rate, like `V_fd E_q'` in
   theirs), not `(Efd, i_fd)`: with constant `Efd` it drops out of the Bregman form, with an
   AVR it is a method choice.
   **Field-port decision:** keep runtime `i_fd`, the IEEET1 interface and the plant dynamics
   unchanged. Define `y_fd^U = B_f^T M z'` only for the strain-energy audit and certificate.
   In the Bregman identity the incremental supply is `(Efd - Efd*) y_fd^U`. With
   `r = B_f^T M (z-z*)`, first test the exact scalar machine-exciter cross-term
   `-(Efd-Efd*)r`, which replaces that supply by `-Efd' r`; accept it only if the joint
   Hessian is positive and the exciter dynamics give decay. Do not modify runtime `i_fd` to
   force this identity.
   Next: verify the full identity
   `d/dt (U_rot + U_net) = -z'^T M z' + Efd B_f^T M z' + P_int delta'` on the lossless
   variant, quantifying the speed-voltage residual. The present `M > 0` and convexity checks
   are Float64 candidate checks with large margins, not yet interval certificates.
   **Option (user, 2026-09-30): dynamic π lines as port-Hamiltonian subsystems.** The
   interconnection can be lossless even when the lines contain resistance: restoring each
   π line's inductor and capacitor states makes the line a PH subsystem (storage
   `L|i|^2/2 + C|v|^2/2`, the dq-frame `omega L`, `omega C` terms skew), terminal powers cancel
   through the power-preserving Kirchhoff interconnection (a Dirac structure), and `R` and
   `G` appear as explicit nonnegative, incrementally passive dissipation. This explains the
   paper's losslessness condition: it concerns the quasi-static phasor network, where
   eliminating the R-L dynamics in the rotating frame turns resistance into transfer
   conductance, whose one-form is non-exact (the `2 (G kron J)` curl of `polar_balance`).
   Caution: the standard static π admittance in Porthos's algebraic phasor Y-bus and a
   dynamic π-line model are not the same system. The construction is established in the
   literature: Fiaz, Zonetti, Ortega, Scherpen, and van der Schaft (2013) assemble generators,
   static loads, and physical π lines into one PH network through power-preserving graph
   interconnections; Gernandt et al. (2021) give the corresponding circuit/Dirac formulation;
   and Gernandt and Hinsen (2024) prove a power balance for lossy telegraph-line networks.
   Porthos still has to derive the synchronous-dq π equations and show that their nominal-
   frequency steady state gives exactly its Y-bus stamp. Consequences: a consistent
   dynamic network needs the machines' stator flux transients too (the `dpsi/dt` stator
   terms GENROU's quasi-static stator drops, the transformer term behind 5a); the result is
   an EMT-like, stiff model with millisecond time constants, not the frozen RMS plant, so a
   certificate for it transfers to the RMS model only through a separate
   singular-perturbation argument (the user's physical-model-change fallback, not the
   storage-only route); constant-power ComplexLoads stay non-passive either way. Plan: keep
   the RMS route (5b) as the main line; as a side experiment after 5b, test the idea on a
   single machine behind a dynamic π line to an infinite bus and on a two-machine case
   (full-order machine with stator transients), checking incremental passivity (Bregman
   Hessian and decay) with and without line resistance, as a reference for how much of the
   RMS conductance residual is an artefact of the quasi-static reduction.
5b. **The joint strain-energy identity and its local decay (2026-09-30; `src/ph/strain.jl`,
   `scripts/strain_energy_check.jl`, about 2.5 minutes, report in
   `outputs/strain/strain_energy_check.txt`; tests in `test/unit/ph.jl`).**
   - *Storage.* `U_ext = U_rot + U_net`: `U_rot = -z'MAz/2` per rotor (`strain_rotor`,
     `strain_metric`: the 10 GENROU and the GENSAL, whose one-state q axis has
     `M = Tq0''/(xq - xq'')`; its saturation is read as the field input `u_f = Efd - Sat (xd - xl)`);
     `U_net = U_B(V) + sum_k [|E1_k|^2/(2x'') - E1_k . V_k/x''] + sum_L U_L` with the
     nominal-speed EMF `E1 = (-psi_q'', psi_d'')`, and the ComplexLoad reactive potential
     `U_L = Q_act(z) ln v - Q0 v^2/(2V0^2) - Phi(z)`, `Phi = Q_act (ln s - 1/kqv)`,
     `s = Vini + z` (valid in the load law's middle band `udmin < |V| < udmax`; samples
     outside it are skipped: 10 of the healthy trajectory records on the case).
   - *Identity* (`strain_balance`, closed form of every term, error 2e-14 relative):
     `dU_ext/dt = -z'Mz' + Efd B_f'Mz' + sum Te delta'` + GENSAL saturation
     + speed voltage `(i1 - i) . Psi z' + (omega - 1) E1 . V'/x''` + conductances
     `[-jGV] . V'` + ComplexLoad active part `[-jGpV] . V'` + load lag
     `(ln v - ln s) dQ_act/dz z'`. `dU_net/d delta_k = Te_k = Pe_k/omega_k` exactly. The lag has
     the sign of `Q0` (`t1 z' = v - s`): it creates energy at the 18 inductive loads (of 19).
     On the lossless variant the residual is the speed voltage and the lag only (as
     predicted); on the case network the conductance term dominates.
   - *Kinetic storage.* With `Te = Pe/omega` and `Tm/omega`, `K = sum omega_b H (omega - 1)^2`
     pairs exactly with `Te delta'`: `dK/dt + sum Te delta' = sum omega_b (omega - 1) Tm/omega
     - omega_b (omega_coi - 1) sum Te` (D = 0 on every machine). With constant `Tm` this leaves
     the dissipation `-omega_b Tm (omega - 1)^2/omega` and a COI coupling term.
   - *Hessians* (quotient, 171 coordinates, at the equilibrium; `omega* = 1`, `Tm* = Te*`):
     the Bregman form of `U_ext + K` has no negative direction on the 64 machine coordinates
     and 16 on the 19 ComplexLoad filter states (the `-Phi` gauge gives `U_zz = -Q'/s`); a
     gauge `c_L dz^2/2` with `c_L >= 5` makes it positive definite on the plant (min 0.51).
     The 88 zero directions are the controllers (no storage yet).
   - *Decay.* `S' = U_ext + K - sum u_f* B_f'Mz - sum Tm* delta - c*'V`,
     `c* = -j(GV* + Gp*V*)`, removes every exact first-order supply (gradient 3e-13 at the
     equilibrium); its rate is a sum of second-order terms whose quadratic forms
     `strain_decay_forms` returns (sum = `sym(H A)` to 4e-15). The Bregman form of
     `U_ext + K` is `S' + C` with `C = Hess(c*'V(eta))`. On the plant (controller states
     frozen, load gauge 5), against the candidate's own Hessian: family `S' + alpha C`,
     alpha < 0.3 not positive (one machine direction), alpha = 0.5 max rate 7.3 /s (13 of 83
     directions creating), alpha = 1 max 15 /s (17). Each term alone (max rate, alpha = 1):
     `loss_curvature` 86, load lag 90 (5 with a larger gauge), active load 22, conductance
     17, COI 3.5, speed voltage 0.9; rotor dissipation, the `Tm/omega` damping and the frozen
     controllers are nonpositive. So the machines' part of the storage works; the local
     obstruction is the network's loss currents (the constant-impedance load conductance and
     the ComplexLoad active power) and the load lag. For comparison, the machine block with
     frozen loads is itself unstable (+0.24 /s), so the loads cannot be left out.
   - *Next construction:* item 4 and the loss repair are one constrained storage-completion
     problem. Give each IEEET1's four feedback states a positive quadratic storage, impose
     the exact field cross-term, include the governors jointly with their machines, and only
     then add Route-B-guided network/load scalar terms. The final function must expose its
     physical, controller and mathematical-repair pieces separately. The `M > 0` and
     convexity checks are Float64, not interval certificates.
6. **Construct `V_ext`**: machine energies + any independently derived exact `U_net` +
   controller terms + scalar cross terms, each
   with a stated origin. Use Route B only to diagnose missing blocks, after checking that its
   blocks survive a change of reference angle (or an orthonormal COI basis) and any exact
   `U_net`. Controllers that cannot be passive alone (IEEEG1: relative degree 2; IEEEG3:
   right-half-plane zero) need a joint machine-governor storage or another port.
7. **Certify `V_P` and `V_ext`** through the same single-mode pipeline and compare the sets.

## ROA certificate pipeline (`src/roa/`, 2026-09-30)

Run: `julia --project=. scripts/certify_roa.jl [scenario.json]` (about 8 minutes). Records:
`outputs/roa/<case>/certificate_quadratic.json` (with `P`, the coordinates and every tried
level) and `roa_check.json`.

- **Theorem** (`certificate_claim`): attraction of the retained physical quotient (171
  coordinates on IEEE-39, modulo the common rotation) to the enclosed equilibrium, while the
  excluded reservoirs and the monitor stay in their recorded feedback-free ranges; see
  `excluded_ranges` in the record (with drift-rate bounds, about 2e-3 per second at most).
- **Coordinates** (`SectionModel`, `section_model(eq)`): the physical coordinates of
  `physical_projection` on the common-angle section, all but one reference rotor angle
  (`GENROU_1.delta`), scaled by powers of two (171 on IEEE-39); the bus voltages stay
  explicit unknowns. The field is projected along the common rotation of every rotor angle:
  PHPS's rounded COI weights do not sum exactly to the rounded total, so the unprojected
  section drifts (the equilibrium is a relative equilibrium). The rotation symmetry is
  declared per model type (`rotation_action`) and checked on the wiring
  (`check_rotation_symmetry`); undeclared types are rejected.
- **Models** (`AbstractSectionModel`): `SectionModel` and `AnalyticModel` (closed-form test
  systems through the same gates).
- **Candidate interface** (`LyapunovCandidate`): `candidate_value`, `candidate_gradient`,
  `positivity_proof`, `sublevel_half_widths`, `gradient_matrix_hull`, `candidate_fingerprint`,
  `candidate_record`; the decay gate is `xi' sym(N'M) xi < 0` with `grad V = N xi` and
  `h = M xi`, so the same gate serves `V_P` (`N = 2P`) and a later `H_ext`.
  `QuadraticCandidate` / `quadratic_candidate(m)`: `A'P + PA = -I` at the equilibrium.
- **Proof steps**: interval methods for the branch primitives (decided or `UndecidedBranch`;
  duals compare by value); `ProofFailure` for "not proved"; equilibrium by Krawczyk on
  `[h; g]`; KCL branch over the box by the parametric Krawczyk test; first-order hull
  `h_eta + h_V Dv` with `Dv` by a verified interval solve; centered hull
  `J(eta*) + [-R, R]` with `R` from nested duals (coordinate by coordinate by default, or one
  interval direction); definiteness by Weyl + Collatz-Wielandt (and the eigenvector
  Gershgorin test, whichever is lower); containment audit on the same boxes (evaluation
  decorations, every branch decided, held states exactly constant, reservoirs and the
  monitor proved not to feed back over their recorded ranges, every contract domain clause in
  interval arithmetic; the expression reader is now generic in the number type).
- **Results on IEEE-39**: `verified_valid_level` `1.71e-10` (centered, per coordinate),
  `1.29e-11` (centered, one direction), `1.88e-12` (first order); 83 domain clauses and 242
  branch sites; above about `1e-6` the ComplexLoad voltage branch (`udmax`) is undecided,
  which is where single-mode certificates end in any case. `ROACheck` recomputes the case,
  contract, system, model and candidate digests and re-proves `P > 0` and the decay by
  interval Cholesky (preconditioned by the inverse Cholesky factor for `P`).
- The PHPS P11 numbers (`2.69e-12`, `5.16e-10`) are in PHPS's own normalisation of `P`;
  Porthos's levels are of the same order and are not meant to match.

## Loose ends

- Parity pack v6 is published (release `parity-pack-v6`, 2026-09-30). Verify a pack's tree
  hash with `Tar.tree_hash` on the decompressed tarball, as `scripts/bind_parity_pack.jl`
  does; hashing a GNU-tar extraction on Windows gives a different (wrong) hash.
- The ROA pipeline and the review closures are committed as one commit (2026-09-30) after a
  green full suite (3043 passes, 2 intentional GFL skips, no failures, 14.5 min).
- Committed on 2026-09-30 after a green full suite (3211 passes, 2 intentional GFL skips,
  no failures, 79 min under load): step 5b, the storage completion, the diagnostics tools
  and this file; not pushed. Before that commit, `origin/master` was at
  `4c84d7d` (the reviewer's documentation revisions). What it contains: step 5b
  (`src/ph/strain.jl`, `scripts/strain_energy_check.jl`, the exports, the shared per-axis
  metric solver in `rotor_gradient_metric`, its tests, this file) and the storage completion
  (`src/ph/sink.jl`, `src/ph/completion.jl`, `certificate_only` in `src/ph/audit.jl`,
  `structured_completion` in the JuMP extension, `scripts/storage_completion.jl`, tests). The
  full suite includes the new strain, variant and completion testsets.
- Committed on 2026-10-01 after a green local suite (2926 passes, 2 intentional GFL skips,
  no failures): converter IDA comparison report-only in P7 (BDF1 is the check there); the
  droop P6 residual exception (2e-12 on its two measured-current rows, every other row
  1e-12); `vi_bisection` through the primitives with `NoModes()`; VOC `pvoc_mode = 1`
  rejected at construction (no reference); the `docs/ROADMAP.md` banner; this file. Not
  pushed: `origin/master` is at `78ae66f`.
- GFL / GFL_ZIF wait for the GFL reservoir model (part of the modelling work).
- Later software: P8 reports, P9 studies (CCT sweeps etc.), P12 Python wrappers and docs.

## Evidence so far: storage and passivity (base case, IEEE-39, at the equilibrium)

The port is done for everything the stability work needs: I/O, network, power flow, the
synchronous-machine set (GENROU, GENSAL, IEEET1, IEEEG1, IEEEG3, COMPLEXLOAD), the
grid-forming converters (GFM_VSM, GFM_DROOP, GFM_VOC), assembly, initialisation, BDF1 and IDA
simulation, the PowerFactory driver, and the P10 port-Hamiltonian audit tools. Sources for
the modelling work: the supervisor's review
(`PHPS_Opt/study/presentation/response_to_reviewer_component_models.tex`),
`PHPS_Opt/phps/PHPS_nonlinear_PH_Lyapunov_ROA_roadmap.md` (0.2b, 0.2c, Target B), and roadmap
Part II.2 (B1 to B7).

1. **The assembled storage is not a Lyapunov function.** `shifted_storage_audit`: Hess H has
   rank 54 of 171 on the common-angle section; sym(SA) has 41 positive eigenvalues (up to
   12.9). The reservoirs are one-way accounts fed at constant power (no storage along plant
   directions); controller states carry no storage.
2. **Controller ports.** IEEEG1 at speed -> Tm is non-passive for any storage (relative
   degree 2: 6 poles, 4 zeros; Re G(jw) changes sign at 1.9347 rad/s). IEEEG3 has a zero at
   +1.333 1/s (Re H < 0 on 1.32 to 14.9 rad/s). IEEET1 not settled. No realisation changes
   these input-output facts, so they need a joint storage or another port pairing.
3. **Governor loops** (`scripts/governor_passivity.jl`): one at a time, 8 IEEEG1 pass against
   the rest of the system; IEEEG1_10's rest is itself non-passive near 1.55 rad/s (a lightly
   damped zero pair at -0.037 +- 1.421i); IEEEG3_11 passes only with a frequency-dependent
   split. **Jointly, no storage split port by port at the Tm ports exists**, static or
   dynamic (`scripts/joint_governor_storage.jl`): fails at 0.048 rad/s (the slow system mode)
   and at 6.59 rad/s. All machines have D = 0; the swing damping is in the damper windings.
   At low frequency the lossy network makes differential speed directions non-passive
   (the path dependence of U_net).
4. **Route A, the cut at the machine terminals with the power V*I**
   (`src/ph/terminal.jl`, `scripts/terminal_passivity.jl`): every unit, and the bare machine,
   is non-passive at its terminal (a constant-power source at DC, indefinite Herm Y(0),
   eigenvalues +-5.5 to +-26; non-passive near its swing frequency; the AVR makes DC much
   worse, -43 to -3422); the network with its loads is non-passive at DC (-18.5). Synchronous
   machines are not incrementally passive in (V, I) terms. Not tried: A', the same units
   with polar ports (P, omega) and (Q, |V|).
5. **Route B, quadratic storage on the reduced dynamics** (`src/ph/structure.jl`,
   `scripts/storage_search.jl`, decay margin gamma(S) = min lambda_max(As'P + PAs) over P on
   a pattern, P >= 0, tr P = 1; dual certificates and Lyapunov checks in interval
   arithmetic):
   - each component alone, or each unit (machine + governor + exciter) alone: margin 0 (at
     most 2.4e-9; the full storage reaches -4.1e-3). A block-local quadratic storage does
     not decay by any useful margin.
   - adding speed-speed, speed-angle or machine-angle couplings across units: still 0.
   - coupling all machine states across machines: a verified Lyapunov function, margin
     -2.0e-6. A dual-guided greedy finds verified patterns from 2233 free entries (-8.4e-7)
     to 2686 (-1.8e-5); it picks machine-machine blocks and above all GENROU_1 (bus 39: the
     reference angle, no controllers, D = 0, the largest inertia). The sparsest-set (L1)
     stage did not converge in its 300 s cap.
   - **Caution before reading it physically**: As = f_x - f_V g_V^-1 g_x eliminates the
     network, so cross-machine terms may be network energy, not joint machine storage; and
     GENROU_1's role may depend on it being the eliminated reference angle.

## Environment on this machine

- **PowerFactory 2022 SP1** at `C:\Program Files\DIgSILENT\PowerFactory 2022 SP1`, Python
  module for 3.10 (Python 3.10.11 is installed: `py -3.10`). Project
  "39 Bus New England System", study case "Base". **Close the PowerFactory window** before
  any `pf/` run (the engine cannot start while it is open: exit code 4002). Compare a
  scenario: `julia --project=. scripts/pf_compare_fault.jl [scenario.json]` (about 5 s of
  PowerFactory).

One-command setup: `powershell -ExecutionPolicy Bypass -File scripts\setup.ps1 [-RunTests]`
(juliaup, Julia 1.12 as a juliaup override for this directory, `Pkg.instantiate`, the parity
pack, the generator venv from `parity/generate/requirements.txt`, MSYS2 UCRT64 g++ and
SUNDIALS at `C:\msys64`). Safe to re-run.

- Julia 1.12.7 via juliaup (`winget` id 9NJNWW8PVKMN), default channel 1.12. `julia` is on
  PATH (WindowsApps alias).
- `Project.toml` compat `julia = "1.12"`; `Manifest.toml` resolved with 1.12.7; CI runs
  1.12 on ubuntu and windows.
- PHPS (read-only reference, only for regenerating the parity pack): venv at
  `parity/generate/.venv`, g++ and SUNDIALS at `C:\msys64\ucrt64`; details in
  `parity/README.md`.
- Run a fault and plot it: `julia --project=. scripts/run_base_fault.jl` (BDF1 and IDA of
  the base-case bus-16 fault; results in `outputs/IEEE39Bus_PF_base_bus16_150ms_{bdf1,ida}/`,
  figures in `outputs/figures/`). The user ran it successfully (2026-09-29); the figures
  then got a distinct colour per machine and legends outside the axes. BDF1 takes about
  80 s (the faithful PHPS scheme: a finite-difference Jacobian of 282 residual calls and a
  dense LU per Newton iteration; PHPS's C++ takes 67 s); IDA about 5 s (PHPS: 6.8 s).
- **Start-up time** (the user's complaint about compiling on every run): `src/precompile.jl`
  is a PrecompileTools workload (load, equilibrium, BDF1 and IDA steps with a fault switch,
  the CSV / JLD2 / run.json writers). Precompiling Porthos takes about 20 s after a source
  change; a fresh process then reaches the simulation in about 2 s (before: about 22 s of
  compilation). A PackageCompiler system image remains possible if the 2 s ever matter.
- Tests: `julia --project=. test/runtests.jl` (about 15 min, mostly the P7 runs of four
  cases). One file alone: `julia --project=. -e 'using Porthos, Test; const ROOT = pwd();
  include("test/parity/common.jl"); include("test/unit/ph.jl")'`.

## What exists

```
src/Porthos.jl            module, exports
src/io/expr.jl            safe expression reader (C order, Float64, env of named values)
src/io/json.jl            read/write JSON keeping key order and int/float; json_identical
src/io/schema.jl          JSONSchema validation (cases/schema/*.schema.json)
src/io/case.jl            Case + typed tables with PHPS defaults; ComponentSpec, Wire
src/io/scenario.jl        Scenario, SolverSettings, BusFault/LineFault/OtherEvent (PHPS DAE defaults)
src/io/contracts.jl       ContractSet / ContractEntry / DomainClause
src/io/params.jl          MODEL_TYPES (PHTRUE set + COMPLEXLOAD), machine-base normalisation,
                          component_params + per-type constructor defaults (type_defaults!)
src/io/parity.jl          parity pack: artifact lookup, SHA-256 verification, readers
                          (pack_json; pack_binary for the sim .bin files)
src/network/ybus.jl       Network (sorted bus ids), ybus / ybus_pf / ybus_dae, Norton stamps,
                          load admittances; CPython complex division for bit parity
src/network/events.jl     fault admittance, fault shunts, with_fault, LineFault line split
src/powerflow/newton.jl   exact port of PHPS Newton-Raphson; skip_pf_solve overrides
src/components/primitives.jl   branch primitives with recorders (gt/ge/lt/le, clamp_mode,
                          nonwindup, outband_relax, guard_min, select); NoModes / ModeLog
src/components/interface.jl    AbstractComponent, rhs!/outputs!/step_outputs!/modes,
                          hamiltonian, injection, norton_admittance, contract,
                          build_component / with_params, COMPONENT_CONSTRUCTORS
src/components/machines/  genrou.jl, gensal.jl
src/components/exciters/  ieeet1.jl
src/components/governors/ ieeeg1.jl, ieeeg3.jl
src/components/loads/     complexload.jl
src/components/converters/ common.jl (vi_bisection, the virtual-impedance divider), gfm_vsm.jl,
                          gfm_droop.jl, gfm_voc.jl
src/assembly/wiring.jl    InputSource, resolve_wiring (PHPS wire semantics + post-init refresh)
src/assembly/dae.jl       DAESystem, assemble (Y-bus at the power-flow voltages, PHPS's C++
                          constant rounding), dae_residual! (port of the C++ dae_residual)
src/assembly/sparsity.jl  jacobian_pattern (structural, valid in every limiter mode)
src/init/components.jl    init_from_phasor (machines), init_from_targets (exciter, governors),
                          converter_init / converter_current (converters),
                          first_pass (PHPS Initializer.run: power split, machine links)
src/init/equilibrium.jl   solve_equilibrium: first pass -> Gauss-Newton on the full DAE,
                          slack set-point closing the power balance, Vini fixed point,
                          reservoir references; component_io
src/components/dispatch.jl     concrete-type dispatch over the models (the residual takes
                          5.6 us and does not allocate on IEEE-39)
src/components/observables.jl  PHPS's logged signals per model (delta_deg, Te, H_steam, ...)
src/assembly/dae.jl       also DAEWorkspace + the allocation-free residual with per-fault flags
src/sim/bdf1.jl           simulate_bdf1: port of PHPS's solve_bdf1; SimResult
src/sim/ida.jl            simulate_ida (Sundials.jl), consistent_voltages! (IDACalcIC-like)
src/sim/results.jl        csv_columns / csv_row / write_results_csv / write_results_jld2,
                          run_metadata (records the settings the integrator used), simulate
src/precompile.jl         PrecompileTools workload for the simulation path
scripts/run_base_fault.jl simulate the base-case bus-16 fault and plot it (one command)
src/io/powerfactory.jl    PowerFactory results reader, pf_simulate / pf_inspect (run pf/),
                          pf_machine_map (case _pf keys), pf_compare (PHPS's metrics)
pf/                       PowerFactory driver (Python 3.10, stdlib): run.py, porthospf/
                          (session, inspect_model, simulate), config.json, README.md
scripts/pf_compare_fault.jl  PowerFactory and Porthos on one scenario, compared (one command)
src/ph/storage.jl         storage components, H / grad H / Hess H, solve_network, reduced_field,
                          reduced_jacobian (exact: f_x - f_V g_V^-1 g_x)
src/ph/audit.jl           PhysicalProjection (contract reservoirs, held states, COI section),
                          shifted_storage_audit
src/ph/ports.jl           PortModel (a component linearised at one port), transfer,
                          real_part_crossings, passivity_certificate, port_zeros,
                          frequency_response (Hessenberg), loop_port_model (the rest of
                          the system seen from one component, exact)
scripts/ph_audit.jl       the audits at Porthos's equilibrium, JSON report (one command)
scripts/governor_passivity.jl  governor shortage vs the rest's excess, per governor
src/ph/dissipativity.jl   MultiPortModel, open_loops_model (several loops at once), multiport_margin,
                          port_margin, loop_margin, kyp_riccati (storage from the KYP Riccati
                          equation), port_storage, rest_storage
scripts/joint_governor_storage.jl  joint storage at the governor ports (negative result)
src/ph/terminal.jl        route A: sync_jacobian, TerminalModel / NetworkModel,
                          terminal_models (the cut at the machine terminals), terminal_margins
scripts/terminal_passivity.jl  route A test, per unit and per layer (negative result)
src/ph/structure.jl       route B: state_groups, storage_pattern, section_pattern,
                          reference_section, pow2_scaling, lyapunov_check, margin_residuals,
                          verified_min_eig (interval), pattern_certificate, verified_lyapunov,
                          rank_couplings, add_coupling!, couple_states!; decay_margin and
                          structured_lyapunov in ext/PorthosJuMPExt
scripts/storage_search.jl route B search (margins, hypotheses, dual-guided greedy, L1;
                          Hypatia, capped at 300 s per solve)
src/roa/interval.jl       interval core: branch primitives on intervals, ProofFailure,
                          interval_solve, verified_max_eig, weyl_max_eig, ellipsoid boxes
src/roa/section.jl        SectionModel, section_model, lift, section_field (projected),
                          state_field, section_residual, section_jacobian
src/roa/candidate.jl      LyapunovCandidate interface, QuadraticCandidate, fingerprints
                          (system_digest, model_fingerprint)
src/roa/enclosure.jl      krawczyk, enclose_equilibrium, enclose_kcl_branch, jacobian_hull,
                          centered_hull (per coordinate or one direction)
src/roa/containment.jl    containment_audit (categories + contract domain clauses)
src/roa/certificate.jl    certify_level, certify_roa (candidate or scenario path), records
src/roa/check.jl          ROACheck: roa_check, interval_cholesky, cholesky_positive_definite
scripts/certify_roa.jl    V_P certificate + ROACheck on IEEE-39 (one command)
src/ph/power.jl           port-power residual audit: component_power, network_power,
                          power_audit_samples, port_power_audit
scripts/port_power_audit.jl  the audit at the equilibrium, random states and the fault
                          trajectory (one command)
src/ph/exchange.jl        GENROU stator exchange through the KCL branch: current_correction,
                          kcl_solve (implicit-function duals), rotor_current_coordinates,
                          exchange_one_form, one_form_exactness, one_form_path_test,
                          exchange_curl, lossless_variant
scripts/exchange_one_form.jl  the exactness gate on the case and its lossless variant
src/ph/joint_storage.jl   rotor_energy (gauge family H_alpha), shifted_kinetic_energy,
                          quotient_hessian (Bregman Hessians), exchange_dissipation_forms
scripts/joint_storage_check.jl  Hessians and exchange-vs-dissipation (one command)
src/ph/polar.jl           network_potential (U_B), polar_balance (lossless identity,
                          conductance part, internal-EMF supply vs exchange), conductance_curl
scripts/polar_balance.jl  the polar network balance (one command)
src/ph/strain.jl          joint strain energy: strain_rotor / strain_metric (GENROU, GENSAL),
                          strain_energy (U_ext), strain_balance (the identity term by term),
                          strain_reference / shifted_strain_energy (S'), strain_decay_forms
scripts/strain_energy_check.jl  step 5b: metric, identity, Hessians, decay attribution (one command)
src/ph/sink.jl            ConstantPowerSink (certificate-only), integrable_loss_variant
src/ph/completion.jl      completion_problem, completion_pattern, invariant_closure, StorageCompletion,
                          storage_completion, completion_energy (the lifted H_ext), completion_hessian_check
scripts/storage_completion.jl  stage 1, stage 2, interval checks and the lift of H_ext (one command)
parity/generate/sim.py    sim section: PHPS's compiled BDF1 / IDA runs (5 ms grid, binary)
ext/                      PorthosMakieExt (CairoMakie); PorthosJuMPExt (JuMP), with the
                          structure-search SDP currently in progress
parity/generate/generate_pack.py   pack generator (network, powerflow, records, components)
parity/generate/components.py      components section: PHPS init, instrumented kernels, sampling
parity/generate/dae.py             dae section: PHPS's C++ kernel compiled with a residual harness
scripts/bind_parity_pack.jl        tarball + Artifacts.toml binding (hash taken from the tarball)
scripts/setup.ps1                  one-command install
test/unit/, test/parity/p0..p7, p10   unit tests (roa.jl: the pipeline, with sampled checks
                                   of the enclosures) and the P0 to P7 and P10 gates
                                   (common.jl: phps_init_params, phps_initial_state)
```

