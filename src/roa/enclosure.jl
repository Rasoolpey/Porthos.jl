# Validated enclosures on the section model (roadmap P11): the equilibrium, the KCL branch
# over a box of section coordinates, and the first-order hull of the reduced Jacobian.
# Interval Jacobians come from ForwardDiff over IntervalArithmetic intervals, through the
# same component code that simulates; the branch primitives decide on the box or throw.

"""
    EquilibriumEnclosure

A box `eta x V` containing exactly one zero of `[h; g]` (the relative equilibrium on the
section and its bus voltages), proved by the Krawczyk test.
"""
struct EquilibriumEnclosure
    eta::Vector{Ival}
    V::Vector{Ival}
    record::Dict{String,Any}
end

"""
    krawczyk(F, J, what, center, r; maxiter = 12) -> (K, X, record)

The Krawczyk test for `F(w) = 0` around `center`: with `C` an approximate inverse of the
Jacobian at the centre and `X = center + [-r, r]`,
`K(X) = center - C F(center) + (I - C J(X)) (X - center)`. `K` in the interior of `X` proves
that `X` contains exactly one zero, and it lies in `K`. On failure the radius is inflated to
twice the reach of `K` (epsilon inflation), at most `maxiter` times. `F` maps an interval
vector to an interval vector; `J` maps an interval box to an interval Jacobian; `what` names
the problem in errors.
"""
function krawczyk(F, J, what::AbstractString, center::AbstractVector{Float64},
                  r::AbstractVector{Float64}; C::AbstractMatrix{Float64}, maxiter::Integer = 12)
    n = length(center)
    Ci = ival.(C)
    Fc = F(ival.(center))
    well_defined(Fc) || fail("$what: the residual at the centre is not well defined")
    base = ival.(center) .- Ci * Fc
    Id = ival.(Matrix{Float64}(I, n, n))
    r = copy(r)
    for it in 1:maxiter
        X = ibox(center, r)
        JX = J(X)
        well_defined(JX) || fail("$what: the Jacobian is not well defined on the box")
        K = base .+ (Id .- Ci * JX) * (X .- ival.(center))
        if strictly_inside(K, X)
            return K, X, Dict{String,Any}("iterations" => it, "radius_max" => maximum(r),
                                          "enclosure_width_max" => maximum(IA.diam, K))
        end
        reach = [max(abs(IA.inf(k) - c), abs(IA.sup(k) - c)) for (k, c) in zip(K, center)]
        all(isfinite, reach) || fail("$what: the Krawczyk operator is unbounded")
        r = max.(2 .* reach, 2 .* r)
    end
    fail("$what: the Krawczyk test did not pass in $maxiter inflations")
end

"""
    enclose_equilibrium(m::AbstractSectionModel; newton = 4) -> EquilibriumEnclosure

The equilibrium on the section: a few Newton steps from `(0, V0)` in Float64, then the
Krawczyk test on `[h; g]` in the unknowns `[eta; V]`.
"""
function enclose_equilibrium(m::AbstractSectionModel; newton::Integer = 4)
    n = neta(m)
    F(w) = section_residual(m, w)
    w = vcat(zeros(n), m.V0)
    for _ in 1:newton
        w .-= ForwardDiff.jacobian(F, w) \ F(w)
    end
    C = inv(ForwardDiff.jacobian(F, w))
    r0 = 4 .* abs.(C * F(w)) .+ 8 .* eps.(max.(abs.(w), 1e-3))
    K, X, rec = krawczyk(F, X -> ForwardDiff.jacobian(F, X), "equilibrium", w, r0; C)
    rec["newton_residual_inf"] = maximum(abs, F(w))
    rec["eta_width_max"] = maximum(IA.diam, K[1:n])
    rec["V_width_max"] = maximum(IA.diam, K[n+1:end])
    return EquilibriumEnclosure(K[1:n], K[n+1:end], rec)
end

"""
    KCLBranch

For every `eta` in the box `eta`, KCL `g(eta, V) = 0` has exactly one solution `V(eta)` in
the box `X`, and it lies in `V` (subset of `X`); `g_V` is nonsingular on `eta x X`.
"""
struct KCLBranch
    eta::Vector{Ival}
    X::Vector{Ival}
    V::Vector{Ival}
    record::Dict{String,Any}
end

"""
    enclose_kcl_branch(m, E, Vc) -> KCLBranch

The KCL branch over the section box `E`, centred at the voltages `Vc` (Float64): the
parametric Krawczyk test on `V -> g(E, V)`.
"""
function enclose_kcl_branch(m::AbstractSectionModel, E::AbstractVector{<:IA.Interval},
                            Vc::AbstractVector{Float64})
    G(V) = section_field(m, E, V)[2]
    Gf(V) = section_field(m, IA.mid.(E), V)[2]
    C = inv(ForwardDiff.jacobian(Gf, Vc))
    spread = [IA.mag(v) for v in ival.(C) * G(ival.(Vc))]
    r0 = 2 .* spread .+ 8 .* eps.(max.(abs.(Vc), 1e-3))
    K, X, rec = krawczyk(G, X -> ForwardDiff.jacobian(G, X), "KCL branch", Vc, r0; C)
    return KCLBranch(collect(E), X, K, rec)
end

"""
    jacobian_hull(m, E, V) -> (M, record)

An interval matrix containing the Jacobian of the network-reduced section field
`eta -> h(eta, v(eta))` at every `eta` in `E` with `v(eta)` in `V` on the KCL branch:
`h_eta + h_V Dv` with `Dv` enclosing `-g_V^-1 g_eta` over `E x V` (`interval_solve`, which
also proves `g_V` nonsingular there). The integral of the Jacobian along any segment in `E`
(the mean-value matrix) lies in the same hull.
"""
function jacobian_hull(m::AbstractSectionModel, E::AbstractVector{<:IA.Interval},
                       V::AbstractVector{<:IA.Interval})
    M, rec, _ = _reduced_jacobian_hull(m, E, V)
    return M, rec
end

function _reduced_jacobian_hull(m::AbstractSectionModel, E, V)
    n = neta(m)
    J = ForwardDiff.jacobian(w -> section_residual(m, w), vcat(collect(E), collect(V)))
    well_defined(J) || fail("the section Jacobian is not well defined on the box")
    h_eta, h_V = J[1:n, 1:n], J[1:n, n+1:end]
    g_eta, g_V = J[n+1:end, 1:n], J[n+1:end, n+1:end]
    Dv, beta = interval_solve(g_V, .-g_eta)
    M = h_eta .+ h_V * Dv
    rec = Dict{String,Any}("g_V_preconditioned_norm" => beta,
                           "hull_width_max" => maximum(IA.diam, M),
                           "Dv_width_max" => maximum(IA.diam, Dv))
    return M, rec, Dv
end

# Tag of the inner (directional) dual of the centered form; created before any outer
# jacobian tag, so ForwardDiff orders the nesting correctly.
_centered_direction(x) = x
const _CenteredTag = typeof(ForwardDiff.Tag(_centered_direction, Ival))

"""
    centered_hull(m, eq, E, V, Xi; second_order = :coordinates, chunk = 1, X = V) -> (M, record)

A tighter enclosure of the mean-value matrix `M(xi) = int_0^1 J(eta* + t xi) dt` for every
equilibrium `eta*` in `eq` and every `xi` in the box `Xi` (with `E` containing `eq.eta + Xi`
and `V` the KCL-branch root box over `E`): the centered form
`M in J(eta*) + [-R, R]`, where `R` bounds half the variation `|d/ds J(p + s xi)|` at every
`p` in `E` along every `xi` in `Xi` (the factor 1/2 is `int_0^1 t dt`). Along a direction
the reduced Jacobian changes by `D = dh_eta + dh_V Dv + h_V dDv`,
`dDv = -g_V^-1 (dg_eta + dg_V Dv)`, with the KCL branch moving as `v' = Dv xi`, and the
derivatives come from nested duals over the box. `second_order = :coordinates` (default)
bounds `R = (1/2) sum_k |xi_k| max |dJ/deta_k|` coordinate by coordinate (171 nested-dual
Jacobians on IEEE-39, about 25 s; `chunk` directions per evaluation, but larger chunks
compile far too slowly); `:direction` uses the single interval direction `[Xi; Dv Xi]` (one
evaluation, about 3 times wider on IEEE-39). The result is intersected entrywise with the first-order hull over `E x V`,
which also encloses `M`.
`J(eta*)` is the first-order hull on the thin equilibrium boxes. The equilibrium voltage box
must lie in the branch's uniqueness box `X` (`KCLBranch.X`): then the equilibrium's voltages
are the branch root `v(eta*)`, which lies in `V` too, and the centre is evaluated on
`eq.V cap V`.
"""
function centered_hull(m::AbstractSectionModel, eq::EquilibriumEnclosure, E::AbstractVector{<:IA.Interval},
                       V::AbstractVector{<:IA.Interval}, Xi::AbstractVector{<:IA.Interval};
                       second_order::Symbol = :coordinates, chunk::Integer = 1,
                       X::AbstractVector{<:IA.Interval} = V)
    n = neta(m)
    all(IA.inf(x) <= IA.inf(e) && IA.sup(e) <= IA.sup(x) for (e, x) in zip(eq.V, X)) ||
        fail("the equilibrium voltages are not inside the KCL branch's uniqueness box")
    M1, rec, Dv = _reduced_jacobian_hull(m, E, V)
    Mc, _, _ = _reduced_jacobian_hull(m, eq.eta, iintersect(eq.V, V))
    R = second_order === :coordinates ? _coordinate_radius(m, E, V, Dv, Xi, chunk) :
        second_order === :direction ? _direction_radius(m, E, V, Dv, Xi) :
        throw(ArgumentError("second_order: :coordinates or :direction"))
    rec["second_order"] = string(second_order)
    C = [c + ival(-r, r) for (c, r) in zip(Mc, R)]
    M = iintersect(C, M1)
    rec["centre_hull_width_max"] = maximum(IA.diam, Mc)
    rec["centered_radius_max"] = maximum(R)
    rec["first_order_hull_width_max"] = rec["hull_width_max"]
    rec["hull_width_max"] = maximum(IA.diam, M)
    rec["entries_from_first_order_hull"] = count(IA.diam.(M1) .< IA.diam.(C))
    return M, rec
end

# R = |D| / 2 with D the directional derivative along the single interval direction
# [Xi; Dv Xi] (one nested-dual Jacobian). Cheap, but the interval tangent forgets that the
# voltage perturbation is tied to xi, which loses the network's cancellations.
function _direction_radius(m::AbstractSectionModel, E, V, Dv, Xi)
    n = neta(m)
    T = vcat(collect(Xi), Dv * collect(Xi))
    w = [ForwardDiff.Dual{_CenteredTag}(b, t) for (b, t) in zip(vcat(collect(E), collect(V)), T)]
    J = ForwardDiff.jacobian(z -> section_residual(m, z), w)
    well_defined(J) || fail("the second-order jets are not well defined on the box")
    dJ = map(y -> ForwardDiff.partials(y)[1], J)
    Jv = map(ForwardDiff.value, J)
    h_V, g_V = Jv[1:n, n+1:end], Jv[n+1:end, n+1:end]
    dh_eta, dh_V = dJ[1:n, 1:n], dJ[1:n, n+1:end]
    dg_eta, dg_V = dJ[n+1:end, 1:n], dJ[n+1:end, n+1:end]
    dDv, _ = interval_solve(g_V, .-(dg_eta .+ dg_V * Dv))
    D = dh_eta .+ dh_V * Dv .+ h_V * dDv
    return [IA.sup(ival(IA.mag(d)) / 2) for d in D]
end

# R = (1/2) sum_k w_k |dJ_red/deta_k| over the box, coordinate by coordinate (PHPS's
# Hessian form, with the voltage tangent of direction k the column Dv[:, k]). The
# derivatives come in chunks of `chunk` directions (nested duals); only their weighted
# magnitudes are kept, and the implicit part is bounded without solving per direction:
#   |dDv_k| <= (I - |R|)^-1 |C| |dg_eta_k + dg_V_k Dv|,  R = I - C g_V, beta = ||R||_inf < 1,
#   (I - |R|)^-1 y <= y + beta/(1 - beta) ||y||_inf 1 (Neumann series), column by column,
# so sum_k w_k |dJ_red_k| <= A1 + A2 |Dv| + |h_V| Y with the weighted sums
# A1, A2, G1, G2 of |dh_eta_k|, |dh_V_k|, |dg_eta_k|, |dg_V_k| and
# Y = (I - |R|)^-1 |C| (G1 + G2 |Dv|). Every product and sum is outward rounded.
function _coordinate_radius(m::AbstractSectionModel, E, V, Dv, Xi, chunk::Integer)
    n = neta(m)
    nv = nvolt(m)
    base = vcat(collect(E), collect(V))
    w = [IA.sup(ival(IA.mag(x))) for x in Xi]
    Z = ival(0.0)
    A1, A2 = fill(Z, n, n), fill(Z, n, nv)
    G1, G2 = fill(Z, nv, n), fill(Z, nv, nv)
    h_V = g_V = nothing
    for k0 in 1:chunk:n
        ks = k0:min(k0 + chunk - 1, n)
        q = length(ks)
        seeds = [ntuple(j -> i <= n ? ival(i == ks[j] ? 1.0 : 0.0) : Dv[i - n, ks[j]], q)
                 for i in 1:n+nv]
        x = [ForwardDiff.Dual{_CenteredTag}(base[i], seeds[i]...) for i in 1:n+nv]
        J = ForwardDiff.jacobian(z -> section_residual(m, z), x)
        well_defined(J) || fail("the second-order jets are not well defined on the box")
        if h_V === nothing
            Jv = map(ForwardDiff.value, J)
            h_V, g_V = Jv[1:n, n+1:end], Jv[n+1:end, n+1:end]
        end
        for (j, k) in enumerate(ks)
            wk = ival(w[k])
            for c in 1:n+nv, r in 1:n+nv
                d = ForwardDiff.partials(J[r, c])[j]
                (IA.inf(d) == 0 && IA.sup(d) == 0) && continue
                t = wk * ival(IA.mag(d))
                if r <= n
                    c <= n ? (A1[r, c] += t) : (A2[r, c - n] += t)
                else
                    c <= n ? (G1[r - n, c] += t) : (G2[r - n, c - n] += t)
                end
            end
        end
    end
    C = inv(IA.mid.(g_V))
    Rm = IA.mag.(ival.(Matrix{Float64}(I, nv, nv)) .- ival.(C) * g_V)
    beta = inorm_inf_up(ival.(Rm))
    beta < 1 || fail("centered hull: ||I - C g_V|| <= $beta is not below 1")
    aDv = ival.(IA.mag.(Dv))
    Y0 = ival.(abs.(C)) * (G1 .+ G2 * aDv)
    f = ival(beta) / (ival(1.0) - ival(beta))
    Y = copy(Y0)
    for j in 1:n
        s = f * ival(maximum(IA.sup, view(Y0, :, j)))
        for i in 1:nv
            Y[i, j] += s
        end
    end
    S = A1 .+ A2 * aDv .+ ival.(IA.mag.(h_V)) * Y
    return [IA.sup(s / 2) for s in S]
end
