# Literature findings: GENROU strain energy, storage choice and Lyapunov construction

Date: 2026-09-30

## The initial non-exact result and its resolution

For the current-corrected GENROU rotor coordinate

```text
w = z - D i,
```

the magnetic balance contains the negative of the complementary work one-form

```text
omega_ex = w' Q D dI,
```

that is, `dot(H_w)` contains `-omega_ex/dt`. On the KCL branch its exterior derivative is
nonzero. Equivalently, the Jacobian of its
coefficient vector is not symmetric. This is invariant under a smooth change of coordinates,
so polar coordinates cannot turn it into the differential of a scalar network energy. The
tested split

```text
omega_ex = -d(1/2 I' D' Q D I) + (D'Qz)' dI
```

isolates an exact gauge term and a non-exact remainder. The latter must appear in the decay
identity **if PHPS's declared `Q` is retained**.

That qualifier is decisive. Nonzero exterior derivative is invariant under coordinate changes,
so polar variables cannot repair this particular one-form. It is not invariant under changing
the storage itself. Applying the strain-energy construction of Ishizaki, Nishino and
Chakrabortty produces a different positive metric `M`. On all ten IEEE-39 GENROU units,

```text
M B_s = -Psi''^T,       M A = A' M,
U_rot(z) = -1/2 z' M A z.
```

The computed `M` and `-sym(MA)` are positive definite with comfortable Float64 margins, and
the armature-reaction term is exactly the gradient of the lossless internal-EMF network
potential. Thus the machine model admits a scalar strain energy. The non-closed form diagnosed
a mismatch between PHPS's magnetic storage and the network-compatible strain energy; it was
not an intrinsic obstruction of GENROU.

The most useful names for this issue are:

- **non-closed (non-exact) work or supply one-form** in differential geometry and mechanics;
- **inexact Pfaffian form** in classical thermodynamics and differential equations;
- **circulatory/nonconservative positional force** in mechanical stability;
- **non-integrable differential feedback or differential supply** in nonlinear control;
- **passivity shortage** or **equilibrium-independent passive-short behavior** when the
  remainder is bounded by a quadratic supply rate.

These names describe related formulations, not interchangeable theorems. In particular,
cyclo-dissipativity permits an indefinite storage and therefore does not by itself produce an
ROA certificate.

## Candidate research contribution

The potentially important contribution is now the extension and certification of a known
strain-energy mechanism, together with a precise diagnosis of an inherited storage choice.
This review did not find a work that establishes all of the following for a detailed GENROU
multimachine DAE:

1. show that the declared PHPS magnetic storage produces a non-closed exchange and identify
   its closed-form curl;
2. derive a positive gradient metric `M` for all four GENROU rotor-circuit states, including
   the damper circuits, with the sign fixed by the internal-EMF network potential;
3. prove the lossless armature-reaction interconnection closes with the resulting convex
   `U_rot`, rather than with the declared `Q`;
4. isolate the residuals caused by speed-scaled voltage, conductances, voltage-dependent loads
   and a dynamic field voltage;
5. construct a joint scalar Lyapunov function `H_ext` for those residuals; and
6. certify a nonlinear ROA for `H_ext` with interval arithmetic, including the algebraic KCL
   branch and limiter modes.

Items 1 and 2 are established numerically with algebraic residual checks; the network identity
in item 3 is the next full-flow check. Items 4 to 6 remain open. Until the identities receive
validated bounds, a broader novelty search is completed, and the nonlinear certificate passes,
this should be described as a **candidate contribution**, not a priority claim.

## What has worked in the literature

### 1. Quantify and dominate the shortage

Shifted passivity, equilibrium-independent dissipativity, and passive-short methods retain a
scalar storage but add a quadratic shortage term to its supply rate. Stability follows when
the interconnection or other components contribute enough excess dissipation. This remains a
fallback framework for terms that survive the strain-energy identity.

For Porthos, define the scalar rate along the reduced field,

```text
sigma_ex(eta) = omega_ex(eta)[f_eta(eta)],
```

and seek a structured quadratic bound rather than one global scalar index. For PHPS's declared
`Q`, rotor plus network conductance loss was exceeded by a factor of 927, ruling out plain
loss dominance for that candidate. The gradient metric `M` supersedes this comparison for the
physical strain-energy route; shortage bounds should be recomputed only for its residual terms.

### 2. Use an exact scalar Lyapunov function that need not be physical energy

The Lyapunov-functions-family and structured LMI literature treats physical energy as a seed,
then searches for exact scalar cross-terms. This is the most direct route to the stated goal:
let `H` denote the Lyapunov function, while recording which pieces are physical storage and
which are certificate terms. Nonzero curl forbids adding the one-form itself to `H`; it does
not forbid another scalar `H_ext` whose derivative dominates the resulting rate.

### 3. Lift from state space to tangent or path space

Differential passivity, differential Lyapunov functions, and control contraction metrics put
the storage on the tangent bundle. A positive differential metric can be integrated along a
curve, usually a minimizing geodesic, to obtain an incremental distance or energy. Manchester
and Slotine explicitly address non-integrable differential feedback by a path/geodesic
construction; their later amendment states an additional integrability condition needed by a
stronger feedback claim.

This is a genuine solution to non-integrability, but it changes the certificate. It does not
make `omega_ex` an endpoint potential. A rigorous ROA implementation would have to certify the
metric and the geodesic/path construction.

### 4. Differentiate the ports or use a rate-based storage

Brayton-Moser power shaping and Krasovskii passivity replace the ordinary energy/passive map
with a mixed potential, a differentiated port, or a storage based on the vector field such as
`f(x)' M(x) f(x) / 2`. This can yield a scalar Lyapunov function even when a chosen work
one-form is not integrable. It is now a fallback if the strain-energy and joint-controller
construction leaves an undominated residual.

### 5. Add certificate dynamics

Dynamic supply rates and IQCs pass the troublesome signals through an auxiliary stable system
and include its state in the storage. Path dependence is represented as memory. The plant
dynamics need not change, but the certificate gains states and terminal-cost conditions. This
is a strong fallback if no useful static `H_ext` exists.

### 6. Separate reversible structure from dissipation

Dirac, GENERIC, and metriplectic formulations keep skew/circulatory exchange separate from the
scalar energy and symmetric dissipation. Helmholtz-Hodge decomposition similarly separates
exact, coexact, and harmonic pieces. These frameworks give the right bookkeeping and can
identify the maximal exact part, but they do not make the coexact remainder a storage.

## Direct implications for Porthos

1. Use `U_rot = -z'MAz/2`, not PHPS's declared magnetic `Q`, as the GENROU rotor strain-energy
   candidate. Prove the metric and convexity with validated arithmetic before certification.
2. Combine it with the exact lossless network potential and verify the complete flow identity.
   The old curl remains a valid negative result only for the discarded `Q`.
3. Keep runtime `i_fd` and the component dynamics unchanged. For the strain energy define the
   audit/certificate flow `y_fd^U = B_f'M zdot`. In the Bregman identity its supply is
   `(Efd-Efd*) y_fd^U` and vanishes when `Efd` is constant.
4. For a dynamic AVR, set `r = B_f'M(z-z*)` and test the exact scalar cross-term
   `-(Efd-Efd*)r`. It converts the field supply to `-Efd_dot r`; a joint machine-exciter
   storage must then prove positivity and decay.
5. Quantify the speed-voltage, conductance and ComplexLoad residuals separately. Do not fold
   them back into the obsolete 927-times shortage of the declared `Q`.
6. Use structured scalar terms guided by Route B only for residuals that remain after this
   physical strain-energy construction. Use Krasovskii or dynamic-supply extensions as
   fallbacks.

## Dynamic π lines: lifting the lossy algebraic network

The user's observation supplies a second, physically stronger route around the conductance
obstruction. A physical π line is an R-L series branch with shunt capacitors (and, when
present, shunt conductances). If its inductor fluxes and capacitor charges are retained as
states, it is a port-Hamiltonian subsystem. Kirchhoff interconnection of its terminal
voltage-current ports is power preserving, so terminal powers cancel when the network is
assembled. The total line balance has the form

`H_line_dot = p_from + p_to - i_series' R i_series - v_shunt' G v_shunt`.

Thus the **interconnection is lossless**, while a line with `R` or `G` is passive and
dissipative rather than lossless. In a synchronous dq frame the nominal-frequency rotation
terms belong to the skew interconnection matrix and do not change the stored energy.

This is an established construction. Fiaz, Zonetti, Ortega, Scherpen, and van der Schaft
build a complete PH power-network model from physical π-line components and a
power-preserving graph interconnection. Gernandt, Haller, Reis, and van der Schaft give the
general PH/Dirac construction for nonlinear RLC circuits, including resistive relations.
Gernandt and Hinsen derive power balance and passivity for networks of lossy distributed
transmission lines governed by the telegraph equations. Structure-preserving discretizations
of those equations give finite-dimensional PH line models.

This does not turn Porthos's present algebraic phasor Y-bus into a dynamic PH network by
renaming it. The Y-bus has already eliminated the line's electromagnetic states. In that
reduced model, resistance appears as transfer conductance and the polar work one-form has the
nonzero curl measured by `polar_balance`. The dynamic π model changes the plant and adds fast,
stiff states; transferring an ROA result back to the RMS DAE would require a proved
singular-perturbation or invariant-slow-manifold argument.

The proposed check is therefore concrete:

1. derive one three-phase π line in synchronous dq coordinates and prove its PH balance;
2. set its electromagnetic derivatives to zero at nominal frequency and recover exactly the
   series and shunt stamps used by Porthos's Y-bus;
3. interconnect one full-order machine, one π line, and an infinite bus and verify the joint
   Bregman balance with `R = G = 0` and with losses;
4. repeat on two machines, then compare the slow eigenstructure with the current RMS DAE;
5. use this either as the physical plant for a new certificate or as the lifted model in a
   rigorous reduction theorem. Do not claim that its certificate already applies to the
   existing 171-state quotient.

## Recommended experiment order

1. **Full strain-energy identity:** verify `U_rot + U_B` along the lossless GENROU flow and
   measure the speed-voltage residual.
2. **Field and load closure:** derive the joint machine-exciter cross-term and the
   voltage-dependent load potential; quantify the conductance shortage separately.
3. **Structured static gate:** use Route B to add exact scalar terms only where the completed
   physical identity still lacks curvature or decay.
4. **Nonlinear certification gate:** Bregman-shift the successful candidate and run the same
   interval positivity, branch-containment, KCL-uniqueness, and decay proof used for `V_P`.
5. **Rate-storage gate:** test `S_K = f_eta' M f_eta / 2` and a Brayton-Moser mixed-potential
   variant in which `dI/dt` is a port variable.
6. **Dynamic extension gate:** add the smallest stable filter that gives a feasible dynamic
   supply/IQC certificate; keep those states certificate-only.

## Core records

The records below were checked through Scite metadata and, where available, indexed full text.
Links point to DOI records or open preprints.

| # | Record | What it contributes here |
|---:|---|---|
| 1 | Willems, “Dissipative dynamical systems Part I: General theory” (1972), [DOI](https://doi.org/10.1007/BF00276493) | Available storage and the general supply-rate framework. |
| 2 | van der Schaft, “Cyclo-Dissipativity Revisited” (2021), [DOI](https://doi.org/10.1109/TAC.2020.3013941) | Closed-cycle characterization when storage need not be bounded below. Useful diagnosis, insufficient alone for ROA. |
| 3 | van der Schaft and Jeltsema, “Limits to Energy Conversion” (2022), [DOI](https://doi.org/10.1109/TAC.2021.3075652) | Cyclo-passivity and partial Legendre transforms for physical systems. |
| 4 | Ortega, García-Canseco, and Stanković, “A cyclo-dissipativity condition for power factor improvement in electrical circuits” (2006), [DOI](https://doi.org/10.1109/ACC.2006.1656596) | Electrical example of a useful cycle inequality without ordinary passivity. |
| 5 | Forni, Sepulchre, and van der Schaft, “On differential passivity of physical systems” (2013), [DOI](https://doi.org/10.1109/CDC.2013.6760930) | Tangent-bundle storage and differential ports for physical systems. |
| 6 | Forni and Sepulchre, “A Differential Lyapunov Framework for Contraction Analysis” (2014), [DOI](https://doi.org/10.1109/TAC.2013.2285771) | Finsler/differential Lyapunov functions integrated into path distances. |
| 7 | Manchester and Slotine, “Control Contraction Metrics: Convex and Intrinsic Criteria for Nonlinear Feedback Design” (2017), [DOI](https://doi.org/10.1109/TAC.2017.2668380) | Convex differential metric; path/geodesic construction for non-integrable differential control. |
| 8 | Manchester and Chaffey, amendment to the preceding paper (2017), [arXiv](https://arxiv.org/abs/1711.08128) | Clarifies the additional integrability condition and limits of the original feedback claim. |
| 9 | Kosaraju, Kawano, and Scherpen, “Krasovskii’s Passivity” (2019), [DOI](https://doi.org/10.1016/j.ifacol.2019.12.005) | Builds storage from system rates/vector fields and changes the passive map. |
| 10 | Kawano, Kosaraju, and Scherpen, “Krasovskii and Shifted Passivity-Based Control” (2021), [DOI](https://doi.org/10.1109/TAC.2020.3040252) | Necessary/sufficient Krasovskii conditions and relations to differential, incremental, and shifted passivity. |
| 11 | Ortega, Jeltsema, and Scherpen, “Power Shaping: A New Paradigm for Stabilization of Nonlinear RLC Circuits” (2003), [DOI](https://doi.org/10.1109/TAC.2003.817918) | Differentiated ports and power-like storage through Brayton-Moser structure. |
| 12 | Favache, Dochain, and Winkin, “Power-shaping control: Writing the system dynamics into the Brayton–Moser form” (2011), [DOI](https://doi.org/10.1016/j.sysconle.2011.04.021) | General PDE/sign conditions for obtaining a mixed-potential representation. |
| 13 | Blankenstein, “Power balancing for a new class of non-linear systems and stabilization of RLC circuits” (2005), [DOI](https://doi.org/10.1080/00207170500036191) | Extends power shaping beyond the original circuit class. |
| 14 | Jeltsema and van der Schaft, “Pseudo-gradient and Lagrangian boundary control system formulation of electromagnetic fields” (2007), [DOI](https://doi.org/10.1088/1751-8113/40/38/013) | Pseudo-gradient/mixed-potential treatment of Maxwell curl equations. |
| 15 | Kosaraju, Pasumarthy, and Jeltsema, “Alternative Passive Maps for Infinite-Dimensional Systems Using Mixed-Potential Functions” (2015), [DOI](https://doi.org/10.1016/j.ifacol.2015.10.205) | Shows that changing the passive map can recover a useful dissipativity identity. |
| 16 | Monshizadeh, Monshizadeh, and Ortega, “Conditions on shifted passivity of port-Hamiltonian systems” (2019), [DOI](https://doi.org/10.1016/j.sysconle.2018.10.010) | Bregman-shifted Hamiltonian, monotonicity tests, and passivity excess/shortage; includes a synchronous-generator example. |
| 17 | Simpson-Porco, “Equilibrium-Independent Dissipativity With Quadratic Supply Rates” (2019), [DOI](https://doi.org/10.1109/TAC.2018.2838664) | General equilibrium-independent quadratic supply rates and a modified Hill-Moylan characterization. |
| 18 | Wu, van der Schaft, and Chen, “Stabilization of Port-Hamiltonian Systems Based on Shifted Passivity via Feedback” (2021), [DOI](https://doi.org/10.1109/TAC.2020.3005156) | Alternative storage and feedback passivation subject to explicit integrability conditions. |
| 19 | Sharf, Jain, and Zelazo, “Geometric Method for Passivation and Cooperative Control of Equilibrium-Independent Passive-Short Systems” (2021), [DOI](https://doi.org/10.1109/TAC.2020.3043390) | Geometric treatment and compensation of equilibrium-independent passivity shortage. |
| 20 | Yang, Liu, and Wang, “Distributed Stability Conditions for Power Systems With Heterogeneous Nonlinear Bus Dynamics” (2020), [DOI](https://doi.org/10.1109/TPWRS.2019.2951202) | Network stability from local passivity indices; directly relevant to shortage domination. |
| 21 | Chellaboina, Haddad, and Kamath, “Dynamic Dissipativity Theory for Stability of Nonlinear Feedback Dynamical Systems” (2005), [DOI](https://doi.org/10.1109/CDC.2005.1582912) | Dynamic operators in the supply-rate framework. |
| 22 | Megretski and Rantzer, “System analysis via integral quadratic constraints” (1997), [DOI](https://doi.org/10.1109/9.587335) | Foundational dynamic-multiplier/IQC framework. |
| 23 | Seiler, “Stability Analysis With Dissipation Inequalities and Integral Quadratic Constraints” (2015), [DOI](https://doi.org/10.1109/TAC.2014.2361004) | Connects IQCs to state-space dissipation inequalities. |
| 24 | Scherer and Veenman, “Stability analysis by dynamic dissipation inequalities” (2018), [DOI](https://doi.org/10.1016/j.sysconle.2018.08.005) | Finite-horizon dynamic multipliers and terminal costs, suitable for invariant-set analysis. |
| 25 | Khong, Chen, and Lanzon, “Feedback stability analysis via dissipativity with dynamic supply rates” (2025), [DOI](https://doi.org/10.1016/j.automatica.2024.112000) | Nonlinear dynamic supply rates with auxiliary systems; closest modern dynamic-extension result. |
| 26 | Caliskan and Tabuada, “Compositional Transient Stability Analysis of Multimachine Power Networks” (2014), [DOI](https://doi.org/10.1109/TCNS.2014.2304868) | Compositional Lyapunov construction for multimachine networks, including network losses. |
| 27 | Vu and Turitsyn, “Lyapunov Functions Family Approach to Transient Stability Assessment” (2016), [DOI](https://doi.org/10.1109/TPWRS.2015.2425885) | Replaces a unique physical energy with an SDP-generated family of scalar Lyapunov functions. |
| 28 | Stegink, De Persis, and van der Schaft, “A Unifying Energy-Based Approach to Stability of Power Grids With Market Dynamics” (2017), [DOI](https://doi.org/10.1109/TAC.2016.2613901) | Port-Hamiltonian interconnection and shifted/Bregman energy for a coupled grid model. |
| 29 | Ishizaki and Chakrabortty, “Necessity of Lossless Transmission and Convexity of Potential Energy Function for Equilibrium Independent Passivity of Power Systems” (2021), [DOI](https://doi.org/10.1109/CDC45484.2021.9683357) | Establishes losslessness and convex potential energy as decisive conditions in a power-system model. |
| 30 | Nishino, Chakrabortty, and Ishizaki, “A Necessary and Sufficient Condition for Equilibrium-Independent Passivity of Power Systems With Two-Axis Generators” (2026), [DOI](https://doi.org/10.1109/TAC.2025.3609489), [open preprint](https://arxiv.org/abs/2304.00987) | Most direct result: lossless transmission is necessary and strain-energy convexity characterizes the EI-passive equilibrium set; its preprint supplied the equations reproduced for GENROU. |
| 31 | Fiaz, Zonetti, Ortega, Scherpen, and van der Schaft, “A port-Hamiltonian approach to power network modeling and analysis” (2013), [DOI](https://doi.org/10.1016/j.ejcon.2013.09.002) | Direct precedent for assembling generators, static loads, and physical π transmission-line elements through power-preserving graph interconnections. |
| 32 | Gernandt, Haller, Reis, and van der Schaft, “Port-Hamiltonian formulation of nonlinear electrical circuits” (2021), [DOI](https://doi.org/10.1016/j.geomphys.2020.103959), [open preprint](https://arxiv.org/abs/2004.10821) | Gives a compositional PH/Dirac formulation of inductors, capacitors, resistive relations, sources, and Kirchhoff interconnection. |
| 33 | Gernandt and Hinsen, “Stability and passivity for a class of distributed port-Hamiltonian networks” (2024), [open preprint](https://arxiv.org/abs/2212.02792) | Proves power balance, passivity, and stability results for networks of lossy telegraph-equation transmission lines under Kirchhoff-type interconnection. |
| 34 | Šešlija, Scherpen, and van der Schaft, “Explicit simplicial discretization of distributed-parameter port-Hamiltonian systems” (2014), [DOI](https://doi.org/10.1016/j.automatica.2013.11.020), [open preprint](https://arxiv.org/abs/1201.5764) | Shows how to discretize a transmission-line PH system while preserving its Dirac/interconnection structure and boundary power ports. |

## Cross-field records that clarify the geometry

- Bulatović, “On the stability and instability criteria for circulatory systems: A review”
  (2020), [DOI](https://doi.org/10.2298/TAM201021013B): mechanical systems with
  nonconservative positional/circulatory forces, where ordinary energy arguments fail.
- Anderson and Thompson, *The Inverse Problem of the Calculus of Variations for Ordinary
  Differential Equations* (1992), [DOI](https://doi.org/10.1090/memo/0473): Helmholtz
  conditions and multipliers for deciding whether equations admit a scalar variational origin.
- Bhatia, Norgard, and Pascucci, “The Helmholtz-Hodge Decomposition—A Survey” (2013),
  [DOI](https://doi.org/10.1109/TVCG.2012.316): exact/coexact/harmonic decomposition and its
  dependence on metric and boundary conditions.
- Courant, “Dirac manifolds” (1990),
  [DOI](https://doi.org/10.1090/S0002-9947-1990-0998124-1): geometric structure for skew
  interconnection that is not itself stored energy.
- Morrison, “A paradigm for joined Hamiltonian and dissipative systems” (1986),
  [DOI](https://doi.org/10.1016/0167-2789(86)90209-5): metriplectic split into antisymmetric
  Hamiltonian and symmetric dissipative brackets.
- Grmela and Öttinger, “Dynamics and thermodynamics of complex fluids. I” (1997),
  [DOI](https://doi.org/10.1103/PhysRevE.56.6620): GENERIC decomposition of reversible and
  irreversible dynamics with degeneracy conditions.

## Reading priority

Read in this order for the immediate implementation decision:

1. Nishino, Chakrabortty, and Ishizaki (2026), then Ishizaki and Chakrabortty (2021).
2. Monshizadeh et al. (2019), Simpson-Porco (2019), Yang et al. (2020), and Sharf et al.
   (2021) for shortage and shifted storage.
3. Kawano et al. (2021) and Ortega et al. (2003) for a `dI/dt`/rate-storage alternative.
4. Khong et al. (2025), Scherer and Veenman (2018), and Seiler (2015) if a dynamic
   certificate is needed.
5. Forni and Sepulchre (2014) and Manchester and Slotine (2017) if a differential/geodesic
   certificate becomes necessary.

This is a structured literature map, not a claim that all thirty papers solve the identical
GENROU problem. The first twenty-five supply the main mathematical remedies; records 26 to 30
are the closest power-system applications and model-specific tests.
