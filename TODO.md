# Porthos.jl progress and hand-off

Read this first in a new session, then `AGENTS.md` and `docs/ROADMAP.md` (its Part I banner:
parity is frozen). Updated at the end of every session.

Last reviewed: 2026-09-30 (step 3, the ROA certificate pipeline for `V_P`, done).

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
discussed with the user before they are implemented.

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
     This is intrinsic to the machine coupling, not to the network losses.
   - **Decision:** retain the non-exact one-form as an explicit, sign-indefinite exchange
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
     not dominate the shortage locally. Controllers carry no storage yet, so their
     dissipation is not in this comparison.
   - So, by the decision above, the shortage needs structured exact scalar cross-terms
     (guided by Route B) or a dynamic extension, and the negative curvature of `H_alpha`
     needs terms that dominate it; both are method choices for the user and the reviewer.
     Step 5's polar network balance can proceed independently (it cannot remove the curl).
   - **Literature check and candidate contribution:** `docs/LITERATURE_FINDINGS.md` records 30
     core papers and six cross-field references. The issue is a non-closed work/supply one-form
     (a circulatory force in mechanics, non-integrable differential supply in control). The literature has
     general remedies, but the review found no work combining the closed-form obstruction for
     a detailed GENROU multimachine DAE with construction and interval ROA certification of a
     repairing scalar `H_ext`; this is a candidate research contribution, subject to a broader
     novelty search. The implementation order is: reproduce the recent two-axis strain-energy
     result; search for a structured exact scalar `H_ext`; then test Krasovskii/Brayton-Moser
     rate storage and, if needed, a certificate-only dynamic-supply/IQC extension. Nonzero curl
     cannot be removed by a coordinate change, and cyclo-dissipativity alone is not an ROA
     certificate.
5. **Polar ports and the network balance** (Route A'): derive the conjugate polar supply from
   the network energy balance (do not assume (P, omega) and (Q, |V|) are the right pairs);
   test the lossless network first, then quantify the passivity shortage from the
   conductances and the active (constant-power) loads. The cheapest remaining test of
   whether unchanged dynamics admit a physically structured `H_ext`.
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
- Step 4 (the port-power audit) and the review corrections (KCL re-solve of the trajectory
  samples, GENROU's rotor port split with the rotor-loss proof) are committed on `master`
  but not pushed; the full suite passed with them (3065 passes, 2 intentional GFL skips, no
  failures). Push when the user says so. The first CI run with pack v6 (`f53bfa0`) was still running when
  the session moved on; the user checks it.
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

