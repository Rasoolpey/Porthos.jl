# Literature findings: non-exact GENROU exchange and Lyapunov construction

Date: 2026-09-30

## The mathematical issue

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
identity unless another physical subsystem supplies the opposite one-form.

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

The potentially important contribution is the combination of a model-specific obstruction
and a constructive, rigorous response. The literature contains the general mathematical
ingredients below, but this review did not find a work that establishes all of the following
for a detailed GENROU multimachine DAE:

1. derive the current-corrected rotor coordinate `w = z - D i` without changing the plant;
2. identify the complementary stator exchange as a one-form on the KCL quotient;
3. prove its non-exactness analytically through its closed-form curl and numerically on both
   lossy and lossless networks;
4. show that the obstruction persists independently of reference angle and network losses;
5. quantify why the obvious physical candidate fails: its Bregman Hessian is indefinite and
   its exchange shortage exceeds the presently proved rotor and network losses;
6. construct a scalar Lyapunov function `H_ext` that handles the remainder honestly; and
7. certify a nonlinear ROA for that `H_ext` with interval arithmetic, including the algebraic
   KCL branch and limiter modes.

Items 1 to 5 are findings already obtained in Porthos. Items 6 and 7 are the open constructive
part. Until a broader systematic novelty search and a complete proof are finished, this should
be described as a **candidate contribution**, not as the first solution in the literature.

The intended claim is not that physical energy is useless. It is that, for this detailed model,
the physical rotor energy generates a non-closed exchange on the reduced state space. The
correct final object is therefore a scalar Lyapunov function `H_ext`: physical storage provides
its interpretable core, and exact cross-terms, rate storage, or certificate-only dynamics repair
its positivity and decay. Every added term remains subject to the same rigorous ROA checks.

## What has worked in the literature

### 1. Quantify and dominate the shortage

Shifted passivity, equilibrium-independent dissipativity, and passive-short methods retain a
scalar storage but add a quadratic shortage term to its supply rate. Stability follows when
the interconnection or other components contribute enough excess dissipation. This is the
closest static framework to Porthos's present identity.

For Porthos, define the scalar rate along the reduced field,

```text
sigma_ex(eta) = omega_ex(eta)[f_eta(eta)],
```

and seek a structured quadratic bound rather than one global scalar index. The current local
calculation says rotor plus network conductance loss can be exceeded by a factor of 927, so a
plain loss-dominance claim is already ruled out. Controller terms and structured cross-terms
can still change the combined inequality.

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
`f(x)' M(x) f(x) / 2`. This can yield a scalar Lyapunov function even when the original work
one-form is not integrable. The appearance of `dI/dt` in GENROU's corrected balance makes this
family particularly relevant.

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

1. Keep `H_w` as the interpretable rotor energy and keep the non-exact exchange in `dot(H)`.
   Do not name the one-form `V_cross` or integrate it along an arbitrary path.
2. Use `H` for the final scalar Lyapunov function. It may contain physical storage, Bregman
   shifts, and exact certificate cross-terms; every term must have a well-defined scalar value
   and pass the quotient Hessian and interval decay checks.
3. Reproduce the Nishino-Ishizaki two-axis strain-energy construction on a reduced Porthos
   model. Their result says losslessness and convexity characterize equilibrium-independent
   passivity for that model class. The Porthos curl survives removal of losses, so either their
   global electromagnetic subsystem uses a different state/port split or the extra GENROU
   rotor circuits break the two-axis structure. This comparison can distinguish the two.
4. Treat the passive-short calculation as a diagnostic. The factor 927 rules out the present
   rotor-plus-network losses as a sufficient bound; repeat it only after controller storage or
   structured exact cross-terms have been added.
5. Search next for a static `H_ext` with structured quadratic cross-terms, seeded by the full
   Lyapunov matrix and constrained by rotation symmetry. Certify positivity and decay with the
   existing interval pipeline.
6. In parallel, test a Krasovskii/Brayton-Moser candidate using `dI/dt`. If the static search
   cannot produce a useful certified set, try a low-order dynamic-supply/IQC extension.
7. Continue the polar network audit because it can reveal useful conjugate variables and
   exact terms. It cannot remove the measured curl.

## Recommended experiment order

1. **Two-axis reproduction gate:** implement the strain energy from the 2021/2026
   Ishizaki-Nishino line on the lossless reduced model; compare its Hessian and supply identity
   with the GENROU audit.
2. **Structured static gate:** solve for exact quadratic cross-terms in `H_ext`, using the
   Route-B matrix as a seed and imposing quotient invariance and a useful decay margin.
3. **Nonlinear certification gate:** Bregman-shift the successful candidate and run the same
   interval positivity, branch-containment, KCL-uniqueness, and decay proof used for `V_P`.
4. **Rate-storage gate:** test `S_K = f_eta' M f_eta / 2` and a Brayton-Moser mixed-potential
   variant in which `dI/dt` is a port variable.
5. **Dynamic extension gate:** add the smallest stable filter that gives a feasible dynamic
   supply/IQC certificate; keep those states certificate-only.

## Thirty core records

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
| 30 | Nishino, Chakrabortty, and Ishizaki, “A Necessary and Sufficient Condition for Equilibrium-Independent Passivity of Power Systems With Two-Axis Generators” (2026), [DOI](https://doi.org/10.1109/TAC.2025.3609489) | Most direct current result: lossless transmission is necessary and strain-energy convexity characterizes the EI-passive equilibrium set for two-axis generators. |

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
