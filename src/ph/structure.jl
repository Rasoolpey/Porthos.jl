# Route B: structure search for a quadratic storage (TODO.md, "Plan: storage that keeps the
# dynamics"; PHPS roadmap 0.2c B2).
#
# Coordinates: the physical coordinates (`PhysicalProjection.keep`) on the common-angle
# section, parametrised by all of them except one reference rotor angle, which the COI
# constraint l'x = 0 determines (`reference_section`). These are plain states, so a pattern
# on them reads component by component. A local storage V = z' Q z must satisfy
#   Q >= I  and  As' Q + Q As <= -eps I,
# with As the exact reduced Jacobian in these coordinates. Q is restricted to an allowed
# pattern (`storage_pattern`), and a solver (the JuMP extension: `structured_lyapunov`) looks
# for one, possibly with an L1 objective on chosen entries, so that the answer shows which
# cross-terms are needed. The solver's Q is a candidate; `lyapunov_check` checks it.

"""
    state_groups(sys, projection) -> Vector{NamedTuple}

For each physical coordinate: the component, its role (`:machine`, `:governor`,
`:exciter`, `:load`, `:other`), the unit (the machine it belongs to: itself, or the machine
it feeds), and whether it is a rotor angle.
"""
function state_groups(sys::DAESystem, projection::PhysicalProjection)
    owner = Dict{Int,Int}()
    for m in sys.coi_members
        owner[m] = m
        for s in sys.sources[m]
            s.kind == SRC_OUTPUT && (owner[s.index] = m)
        end
    end
    comp_of = zeros(Int, sys.n_diff)
    for k in eachindex(sys.comps), i in _state_range(sys, k)
        comp_of[i] = k
    end
    out = NamedTuple[]
    for i in projection.keep
        k = comp_of[i]
        c = sys.comps[k]
        role = component_role(c)
        r = role === :generator ? :machine : role === :governor ? :governor :
            role === :exciter ? :exciter : role === :load ? :load : :other
        s = state_names(c)[i - sys.offsets[k] + 1]
        u = get(owner, k, 0)
        push!(out, (component = name(c), role = r, unit = u == 0 ? name(c) : name(sys.comps[u]),
                    angle = r === :machine && s in ("delta", "theta"), state = s))
    end
    return out
end

"""
    storage_pattern(groups, kind) -> BitMatrix

Allowed nonzero entries of `P` (symmetric):
- `:full` every entry;
- `:unit` within each unit (machine, governor, exciter together), within each other
  component (the loads), and among all rotor angles;
- `:component` within each component only, and among all rotor angles.
"""
function storage_pattern(groups::AbstractVector, kind::Symbol)
    n = length(groups)
    kind === :full && return trues(n, n)
    kind in (:unit, :component) || throw(ArgumentError("pattern $kind: use :full, :unit or :component"))
    key = kind === :unit ? (g -> g.unit) : (g -> g.component)
    return BitMatrix([key(groups[i]) == key(groups[j]) || (groups[i].angle && groups[j].angle)
                      for i in 1:n, j in 1:n])
end

"""
    section_pattern(mask, U) -> BitMatrix

The support, in section coordinates, of `Q = U' P U` for `P` on the physical-coordinate
pattern `mask` (`U` from `reference_section`). The reference angle is `-l_z'z / l_ref`, so a
term coupling a state with the reference angle becomes couplings of that state with every
rotor angle; a pattern written directly on the section coordinates would miss them (for the
reference machine's own states in particular). The support is a superset of the image
subspace (its entries are free, not tied to the COI weights).
"""
section_pattern(mask::AbstractMatrix{Bool}, U::AbstractMatrix) =
    BitMatrix(abs.(U)' * Float64.(mask) * abs.(U) .> 0)

"""
    reference_section(sys, x, V, projection) -> (z, U, As, ref)

The section coordinates: `z` indexes the physical coordinates kept (all but the reference
angle `ref`, the COI source with the largest weight); `U` maps them to the physical
coordinates (identity, plus the row `-l_z'/l_ref` for the reference angle); `As` is the exact
reduced Jacobian in these coordinates (`E A U`, `E` selecting `z`), valid because the section
is invariant.
"""
function reference_section(sys::DAESystem, x::AbstractVector, V::AbstractVector,
                           projection::PhysicalProjection)
    keep = projection.keep
    l = projection.coi[keep]
    ref = argmax(l)
    z = [k for k in eachindex(keep) if k != ref]
    U = zeros(length(keep), length(z))
    for (j, k) in enumerate(z)
        U[k, j] = 1.0
    end
    U[ref, :] .= -l[z] ./ l[ref]
    A = reduced_jacobian(sys, x, V)[keep, keep]
    maximum(abs, l' * A) <= 1e-9 * max(1.0, maximum(abs, A)) ||
        error("the common-angle section is not invariant")
    return z, U, A[z, :] * U, ref
end

"""
    lyapunov_check(Q, As) -> NamedTuple

The smallest eigenvalue of `Q` and the largest of `As'Q + QAs` (a local Lyapunov function
of the section dynamics when the first is > 0 and the second < 0).
"""
function lyapunov_check(Q::AbstractMatrix, As::AbstractMatrix)
    all(isfinite, Q) || return (Q_min = NaN, rate_max = NaN, lyapunov = false)
    qmin = eigmin(Symmetric((Q + Q') / 2))
    rmax = eigmax(Symmetric(As' * Q + Q * As))
    return (Q_min = qmin, rate_max = rmax, lyapunov = qmin > 0 && rmax < 0)
end

"""
    structured_lyapunov(As, mask; weights = nothing, eps = 1e-3, optimizer, silent = true)
        -> (Q, status)

Solve for `Q` on the pattern `mask` with `Q >= I` and `As'Q + QAs <= -eps I`, minimising
`sum(weights .* |Q|)` (a feasibility problem when `weights` is `nothing`). Implemented in the
JuMP extension (`using JuMP` and a conic solver, e.g. SCS). The result is a candidate: check
it with `lyapunov_check`.
"""
function structured_lyapunov end

"""
    structured_completion(A, Hfix, mask, basis; optimizer, mu = 0, nonnegative, silent = true)
        -> (P, t, coefficients, F, record, dual_P, dual_rate)

The constrained storage completion: `P = Hfix + F + sum_j k_j basis[j]` with `Hfix` fixed (the
physical core and its fixed cross-terms), `F` free on the pattern `mask` (controller blocks)
and scalar coefficients `k_j` on the structured repair matrices `basis[j]` (`k_j >= 0` where
`nonnegative[j]`); maximises `t` subject to `P >= t I` and `-(A'P + PA) - 2 mu P >= t I`.
`t > 0` is a strict quadratic Lyapunov function with decay rate `mu` in these coordinates
(scale them first, e.g. by powers of two). With `weights` (on the free entries) and
`basis_weights`, it instead minimises the weighted L1 norm of the free entries and the
coefficients subject to `t >= t_min`: the sparsest completion with that margin. Implemented
in the JuMP extension; the result is a candidate, to be checked (`lyapunov_check`,
`verified_lyapunov`).
"""
function structured_completion end

"""
    completion_rate_feasibility(A, Hfix, mask, basis; optimizer, mu, cond_cap, nonnegative)
        -> (feasible, gamma, P, record, dual_P, dual_cond, dual_rate)

Whether the completion `P = Hfix + F + sum k_j basis[j]` (as in `structured_completion`)
can reach the normalised decay rate `mu` with a bounded condition number: maximise `gamma`
subject to `gamma I <= P <= cond_cap gamma I` and `-(A'P + PA) >= 2 mu P`. `feasible` when
`gamma > 0`; otherwise the solver's duals are the (Float64) infeasibility certificate.
`generalized_decay_rate` gives the rate of a given `P`.
"""
function completion_rate_feasibility end

"""
    generalized_decay_rate(P, A) -> Real

The largest `mu` with `-(A'P + PA) >= 2 mu P` for `P > 0`: half the smallest eigenvalue of
`-(A'P + PA)` against `P`. Invariant under a change of coordinates (a congruence of `P`
with the similarity of `A`).
"""
function generalized_decay_rate(P::AbstractMatrix, A::AbstractMatrix)
    L = cholesky(Symmetric((P + P') / 2)).L
    R = -(A' * P + P * A)
    return eigmin(Symmetric(L \ R / L')) / 2
end

# ---------------------------------------------------------------------------------------
# Decay margin of a pattern, its dual certificate, and the ranking of missing couplings
# (TODO.md, "Route B plan", steps 3 to 5).
#
# For a pattern S (symmetric, diagonal included) the decay margin is
#   gamma(S) = min { lambda_max(As'P + P As) : P = P' on S, P >= 0, tr P = 1 }.
# S carries a strict quadratic Lyapunov function (P > 0, As'P + PAs < 0) exactly when
# gamma(S) < 0 (such a P is automatically definite). The problem is compact and homogeneous,
# so it needs no epsilon and no upper bound on P, and it says how far a pattern is from
# working. Its dual gives the certificate: for any Z = L L' (L arbitrary) and any symmetric Y
# that agrees with As Z + Z As' on S,
#   gamma(S) >= lambda_min(Y) / tr Z,
# because tr(Z) lambda_max(As'P + PAs) >= <Z, As'P + PAs> = <As Z + Z As', P> = <Y, P>
# >= lambda_min(Y) tr P for every admissible P (P vanishes off S). So Y > 0 proves that no
# quadratic Lyapunov function with pattern S exists. The solver's duals supply Z and the
# entries of Y off S; `pattern_certificate` checks the bound in interval arithmetic.
# Off S, G = As Z + Z As' - M (M the dual of P >= 0) is the derivative of the margin with
# respect to freeing those entries: `rank_couplings` ranks the missing blocks by it.
# A diagonal congruence keeps every pattern and the sign of gamma; `pow2_scaling` picks
# powers of two, so the scaled matrix is exactly T^-1 As T in floating point.

"""
    pow2_scaling(d) -> Diagonal

The diagonal scaling `T` with `T_ii` the power of two nearest to `1/sqrt(d_i)` (e.g. `d` the
diagonal of the dense Lyapunov solution). Scaling by powers of two is exact in floating point.
"""
pow2_scaling(d::AbstractVector) = Diagonal([exp2(round(-log2(v) / 2)) for v in d])

"""
    decay_margin(As, mask; optimizer, weights = nothing, rate = nothing, silent = true)
        -> NamedTuple

The decay margin `gamma(S)` of the pattern `mask` (see the comment above `pow2_scaling`):
minimise `lambda` over `P` on the pattern with `P >= 0`, `tr P = 1`, `lambda I - (As'P + PAs)
>= 0`. With `weights` (and `rate > 0`), minimise instead `sum(weights .* |P|)` subject to
`As'P + PAs <= -rate I` (the fewest weighted entries that keep a margin `rate`).
Returns `P`, `gamma` (the computed `lambda_max(As'P + PAs)`), the duals `Z` (of the rate
constraint) and `M` (of `P >= 0`), and `record` (solver statuses, time, iterations, and
residuals computed here, not taken from the solver). Implemented in the JuMP extension.
The result is a candidate; `pattern_certificate` and `verified_lyapunov` judge it.
"""
function decay_margin end

"""
    margin_residuals(As, mask, P, Z, M) -> Dict

Primal and dual residuals of a `decay_margin` solution, computed in Float64 from the returned
matrices: `P` off the pattern, `lambda_min(P)`, `|tr P - 1|`; `lambda_min(Z)`, `|tr Z - 1|`,
`lambda_min(M)`; the dual stationarity on the pattern (`As Z + Z As' - M` a multiple of the
identity on S) and the duality gap.
"""
function margin_residuals(As::AbstractMatrix, mask::AbstractMatrix{Bool}, P::AbstractMatrix,
                          Z::AbstractMatrix, M::AbstractMatrix)
    n = size(As, 1)
    G = As * Z + Z * As' - M
    nu = sum(G[i, i] for i in 1:n) / n
    stat = maximum(abs(G[i, j] - (i == j ? nu : 0.0)) for i in 1:n, j in 1:n if mask[i, j])
    lam = eigmax(Symmetric(As' * P + P * As))
    return Dict{String,Any}(
        "P_off_pattern" => maximum((abs(P[i, j]) for i in 1:n, j in 1:n if !mask[i, j]); init = 0.0),
        "P_min_eig" => eigmin(Symmetric(P)), "P_trace_error" => abs(tr(P) - 1),
        "Z_min_eig" => eigmin(Symmetric(Z)), "Z_trace_error" => abs(tr(Z) - 1),
        "M_min_eig" => eigmin(Symmetric(M)), "dual_stationarity" => stat,
        "dual_value" => nu, "primal_value" => lam, "duality_gap" => lam - nu)
end

"""
    verified_min_eig(Y) -> Float64

A rigorous lower bound `c` on the smallest eigenvalue of every symmetric matrix in the
interval matrix `Y` (`-Inf` when the check fails). Method: `X` the approximate eigenvectors of
`mid(Y)`; if the interval matrix `X'(Y - cI)X` is strictly diagonally dominant with a
positive diagonal (Gershgorin), then `X'(Y - cI)X` is positive definite, `X` is nonsingular,
and by Sylvester's law of inertia `Y - cI > 0`.
"""
function verified_min_eig(Y::AbstractMatrix{<:IntervalArithmetic.Interval})
    I_ = IntervalArithmetic
    m = I_.mid.(Y)
    F = eigen(Symmetric((m + m') / 2))
    X = I_.interval.(F.vectors)
    Xt = permutedims(X)
    l1 = F.values[1]
    s = maximum(abs, F.values)
    for d in (1e-6, 1e-4, 1e-2, 1e-1, 0.5)
        c = l1 - d * max(abs(l1), 1e-12 * s)
        Yc = copy(Y)
        for i in axes(Yc, 1)
            Yc[i, i] -= I_.interval(c)
        end
        B = Xt * Yc * X
        ok = all(axes(B, 1)) do i
            I_.inf(B[i, i]) > sum((I_.mag(B[i, j]) for j in axes(B, 2) if j != i); init = 0.0)
        end
        ok && return c
    end
    return -Inf
end

"""
    pattern_certificate(As, mask, Z, M) -> NamedTuple

Checks the dual certificate of `decay_margin` for the pattern `mask`: `L = V sqrt(max(D, 0))`
from the eigen-decomposition of `Z` (so `Z = LL'` exactly, as a real matrix), `Y` equal to
`As Z + Z As'` on the pattern and to `M` off it. Returns the Float64 bound `lambda_min(Y) /
tr Z` and the rigorous one (`L L'`, `As Z + Z As'` and the trace in interval arithmetic,
`verified_min_eig`), for the matrix `As` as given. `infeasible = true` means the rigorous
bound is positive: no quadratic Lyapunov function with this pattern exists for `As`.
"""
function pattern_certificate(As::AbstractMatrix, mask::AbstractMatrix{Bool},
                             Z::AbstractMatrix, M::AbstractMatrix)
    I_ = IntervalArithmetic
    n = size(As, 1)
    F = eigen(Symmetric((Z + Z') / 2))
    L = F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
    Zf = L * L'
    Wf = As * Zf + Zf * As'
    Yf = [mask[i, j] ? Wf[i, j] : (M[i, j] + M[j, i]) / 2 for i in 1:n, j in 1:n]
    float_bound = eigmin(Symmetric(Yf)) / tr(Zf)
    Li = I_.interval.(L)
    Ai = I_.interval.(As)
    Zi = Li * permutedims(Li)
    Wi = Ai * Zi + Zi * permutedims(Ai)
    Yi = [mask[i, j] ? Wi[i, j] : I_.interval((M[i, j] + M[j, i]) / 2) for i in 1:n, j in 1:n]
    trZ = sum(Li .^ 2)
    c = verified_min_eig(Yi)
    verified = c > 0 ? c / I_.sup(trZ) : (isfinite(c) ? c / I_.inf(trZ) : -Inf)
    return (float_bound = float_bound, verified_bound = verified, infeasible = verified > 0)
end

"""
    verified_lyapunov(P, As) -> NamedTuple

Rigorous check of a candidate `P` (a strict Lyapunov function of `As` when both bounds are
positive): lower bounds on `lambda_min(P)` and on `lambda_min(-(As'P + PAs))`, in interval
arithmetic (`verified_min_eig`), for the matrices as given.
"""
function verified_lyapunov(P::AbstractMatrix, As::AbstractMatrix)
    I_ = IntervalArithmetic
    Pi = I_.interval.((P + P') / 2)
    Ai = I_.interval.(As)
    R = -(permutedims(Ai) * Pi + Pi * Ai)
    p = verified_min_eig(Pi)
    r = verified_min_eig(R)
    return (P_min = p, decay_min = r, lyapunov = p > 0 && r > 0)
end

"""
    rank_couplings(G, mask, groups; by = :component) -> Vector{NamedTuple}

The blocks missing from the pattern, ranked by the Frobenius norm of `G` over their
off-pattern entries (`G = As Z + Z As' - M` from `decay_margin`'s duals: the rate at which
freeing those entries lowers the margin). `by = :component` groups coordinates by component,
`:unit` by unit, `:state` lists single entries.
"""
function rank_couplings(G::AbstractMatrix, mask::AbstractMatrix{Bool}, groups::AbstractVector;
                        by::Symbol = :component)
    key = by === :component ? (i -> groups[i].component) :
          by === :unit ? (i -> groups[i].unit) :
          by === :state ? (i -> groups[i].component * "." * groups[i].state) :
          throw(ArgumentError("by = $by: use :component, :unit or :state"))
    acc = Dict{Tuple{String,String},Float64}()
    cnt = Dict{Tuple{String,String},Int}()
    n = size(G, 1)
    for j in 1:n, i in 1:j
        mask[i, j] && continue
        a, b = key(i), key(j)
        k = a <= b ? (a, b) : (b, a)
        acc[k] = get(acc, k, 0.0) + (i == j ? 1 : 2) * G[i, j]^2
        cnt[k] = get(cnt, k, 0) + 1
    end
    out = [(a = k[1], b = k[2], norm = sqrt(v), entries = cnt[k]) for (k, v) in acc]
    return sort!(out; by = r -> -r.norm)
end

"""
    couple_states!(mask, groups, f, g) -> mask

Frees every entry between a coordinate whose group satisfies `f` and one whose group
satisfies `g`, symmetrically (e.g. `f = g -> g.state == "omega"`, `g = g -> g.angle`: all
speed-angle cross-terms).
"""
function couple_states!(mask::AbstractMatrix{Bool}, groups::AbstractVector, f, g)
    for i in eachindex(groups), j in eachindex(groups)
        if f(groups[i]) && g(groups[j])
            mask[i, j] = true
            mask[j, i] = true
        end
    end
    return mask
end

"""
    add_coupling!(mask, groups, a, b; by = :component) -> mask

Frees every entry between the coordinates of `a` and those of `b` (component or unit names,
per `by`), symmetrically.
"""
function add_coupling!(mask::AbstractMatrix{Bool}, groups::AbstractVector, a::AbstractString,
                       b::AbstractString; by::Symbol = :component)
    key = by === :component ? (g -> g.component) : by === :unit ? (g -> g.unit) :
          throw(ArgumentError("by = $by: use :component or :unit"))
    ia = findall(g -> key(g) == a, groups)
    ib = findall(g -> key(g) == b, groups)
    (isempty(ia) || isempty(ib)) && throw(ArgumentError("no coordinates for $a or $b"))
    for i in ia, j in ib
        mask[i, j] = true
        mask[j, i] = true
    end
    return mask
end
