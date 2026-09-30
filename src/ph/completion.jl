# The storage completion H_ext = H_physical + H_controller + H_repair (TODO.md, method
# decision of 2026-09-30).
#
# H_physical  the supply-shifted strain energy S' (`shifted_strain_energy`: rotor strain
#             energy, the torque-consistent kinetic storage, the reactive-load and network
#             potentials) plus the fixed field cross-terms -(Efd - Efd*) B_f'M (z - z*), which
#             cancel each machine's incremental field supply.
# H_controller  quadratic storage on the feedback states of each exciter and governor and
#             their cross-blocks with their own machine (the one-way reservoirs excluded),
#             determined by the completion LMI, not required to be passive alone.
# H_repair    reference-invariant scalar terms: quadratic blocks between units, the swing
#             block (the angle and speed deviations of every machine), and per ComplexLoad
#             the filter-error term kappa/2 (|V| - Vini - z)^2 and c/2 (z - z*)^2; optionally
#             eps_k (delta_k - delta_k*)(omega_k - 1). Not physical energy.
#
# The quadratic parts live on the common-angle section (state deviations of the kept
# coordinates, the reference angle eliminated through the COI constraint, which the flow
# preserves). A block pattern is invariant under a change of angle parametrization (another
# reference machine, an orthonormal COI basis) exactly when, for every non-angle coordinate,
# its free angle entries are all or none; `invariant_closure` enforces that.
#
# `completion_problem` builds the LMI data of `structured_completion` (the JuMP extension)
# on a section model with `scale = :none`; `storage_completion` turns a solution into the
# explicit nonlinear scalar `completion_energy`, whose Hessian at the equilibrium is the
# LMI's matrix (`completion_hessian_check`).

# section coordinates of a component's kept states
_section_coords(m::SectionModel, coord::Dict{Int,Int}, k::Integer) =
    [coord[i] for i in _state_range(m.sys, k) if haskey(coord, i)]

# the controller feeding input `nm` of machine component k, or 0
function _feeding(sys::DAESystem, k::Integer, nm::AbstractString)
    j = findfirst(==(nm), input_names(sys.comps[k]))
    j === nothing && return 0
    s = sys.sources[k][j]
    return s.kind === SRC_OUTPUT ? s.index : 0
end

"""
    completion_problem(m; forms = strain_decay_forms(m)) -> NamedTuple

The data of the storage-completion LMI on the section model `m` (`scale = :none`, at an
equilibrium): `A` (section Jacobian), `H` (Hessian of S'), `Hfix` (`H` plus the field
cross-terms), `own` (the own-unit controller pattern: exciter and governor blocks and their
cross-blocks with their machine), `swing` (every angle and speed pair), `basis` (records with
`kind` = `:load_filter`, `:load_z` or `:delta_omega`, the component index `k` and the
`matrix`), `names`, `units` (machine, exciter and governor component indices) and `field`
(per machine: its exciter's `xe` state, `M B_f`, the rotor state indices).
"""
function completion_problem(m::SectionModel; forms = nothing)
    sys = m.sys
    r = forms === nothing ? strain_decay_forms(m) : forms
    H, A = r.H, r.A
    n = neta(m)
    coord = Dict(m.keep[m.z[j]] => j for j in 1:n)
    Hfix = copy(H)
    own = falses(n, n)
    units = NamedTuple[]
    field = NamedTuple[]
    for (k, c) in enumerate(sys.comps)
        sr = strain_rotor(c)
        sr === nothing && continue
        mc = _section_coords(m, coord, k)
        ke, kg = _feeding(sys, k, "Efd"), _feeding(sys, k, "Tm")
        for kc in (ke, kg)
            kc == 0 && continue
            ctrl = _section_coords(m, coord, kc)
            own[ctrl, ctrl] .= true
            own[ctrl, mc] .= true
            own[mc, ctrl] .= true
        end
        push!(units, (machine = k, exciter = ke, governor = kg))
        ke == 0 && continue
        je = findfirst(==("xe"), state_names(sys.comps[ke]))
        je === nothing && continue
        xe = _state_range(sys, ke)[je]
        MBf = strain_metric(sr).M * sr.Bf
        zs = collect(_state_range(sys, k)[sr.states])
        e = coord[xe]
        for (i, zi) in enumerate(zs)
            Hfix[e, coord[zi]] -= MBf[i]
            Hfix[coord[zi], e] -= MBf[i]
        end
        push!(field, (machine = k, xe = xe, MBf = MBf, z = zs))
    end
    names = m.names
    sw = findall(nm -> any(endswith(nm, s) for s in (".delta", ".omega")) &&
                       strain_rotor(sys.comps[findfirst(c -> name(c) == split(nm, ".")[1], sys.comps)]) !== nothing, names)
    swing = falses(n, n)
    swing[sw, sw] .= true
    # repair basis: load-filter error and filter curvature per ComplexLoad, delta-omega per machine
    Vc = m.V0
    lk = [(k, sys.inj[k][1]) for (k, c) in enumerate(sys.comps) if c isa COMPLEXLOAD && c.p.t1 > 0.0]
    Jv = ForwardDiff.jacobian(e -> (x = lift(m, e); V = kcl_solve(sys, x, Vc);
                                    [sqrt(V[2b-1]^2 + V[2b]^2) for (_, b) in lk]), zeros(n))
    Lx = ForwardDiff.jacobian(e -> lift(m, e), zeros(n))
    basis = NamedTuple[]
    for (l, (k, _)) in enumerate(lk)
        iz = coord[first(_state_range(sys, k))]
        w = Jv[l, :]
        w[iz] -= 1.0
        push!(basis, (kind = :load_filter, k = k, matrix = w * w', nonnegative = true))
        E = zeros(n, n)
        E[iz, iz] = 1.0
        push!(basis, (kind = :load_z, k = k, matrix = E, nonnegative = true))
    end
    for u in units
        r_ = _state_range(sys, u.machine)
        ud, uw = Lx[r_[1], :], Lx[r_[2], :]
        push!(basis, (kind = :delta_omega, k = u.machine, matrix = ud * uw' + uw * ud', nonnegative = false))
    end
    return (A = A, H = H, Hfix = Hfix, own = own, swing = swing, basis = basis, names = names,
            units = units, field = field, m = m)
end

"""
    completion_pattern(prob, pairs) -> BitMatrix

The free pattern `own .| swing` plus the full blocks between the listed component pairs
(`(name_a, name_b)`), closed with `invariant_closure`.
"""
function completion_pattern(prob, pairs)
    mask = prob.own .| prob.swing
    comp(nm) = split(nm, ".")[1]
    for (a, b) in pairs
        ia = findall(nm -> comp(nm) == a, prob.names)
        ib = findall(nm -> comp(nm) == b, prob.names)
        (isempty(ia) || isempty(ib)) && throw(ArgumentError("no section coordinates for $a or $b"))
        mask[ia, ib] .= true
        mask[ib, ia] .= true
    end
    return invariant_closure(mask, prob.names)
end

"""
    invariant_closure(mask, names) -> BitMatrix

The smallest superset of `mask` that is invariant under a change of the angle
parametrization of the common-angle section: every non-angle coordinate with a free angle
entry gets all angle entries, and the angle-angle block is free as soon as any entry is.
"""
function invariant_closure(mask::AbstractMatrix{Bool}, names::AbstractVector{<:AbstractString})
    isang = [endswith(nm, ".delta") for nm in names]
    cl = BitMatrix(mask .| mask')
    if any(cl[isang, isang])
        cl[isang, isang] .= true
    end
    for y in findall(.!isang)
        if any(cl[isang, y])
            cl[isang, y] .= true
            cl[y, isang] .= true
        end
    end
    return cl
end

"""
    StorageCompletion

An explicit nonlinear `H_ext` built by `storage_completion`: the physical core (S' with the
fixed field cross-terms), the controller quadratic `Fc` and the repair quadratic `Fr` (section
coordinates) and the repair basis terms with their coefficients.
"""
struct StorageCompletion
    m::SectionModel
    ref::NamedTuple
    field::Vector{NamedTuple}
    Fc::Matrix{Float64}
    Fr::Matrix{Float64}
    terms::Vector{NamedTuple}        # (kind, k, coefficient)
end

"""
    storage_completion(prob, sol; T = I) -> StorageCompletion

From a solution `sol` of `structured_completion` on `prob` (in coordinates scaled by the
diagonal `T`: `A_s = T^-1 A T`, `H_s = T H T`), the explicit storage: the free part split
into controller (`prob.own`) and repair entries, the basis coefficients attached to their
terms.
"""
function storage_completion(prob, sol; T = I, basis_index = eachindex(prob.basis))
    Ti = inv(T)
    F = Ti * sol.F * Ti
    Fc = F .* prob.own
    Fr = F .* .!prob.own
    m = prob.m
    ref = strain_reference(m.sys, m.x0, m.V0)
    terms = [(kind = prob.basis[j].kind, k = prob.basis[j].k, coefficient = sol.coefficients[i])
             for (i, j) in enumerate(basis_index)]
    return StorageCompletion(m, ref, prob.field, Fc, Fr, terms)
end

"""
    completion_energy(sc, eta; V0 = sc.m.V0) -> NamedTuple

`(physical, controller, repair, total)` of the explicit storage at section coordinates `eta`
(voltages on the KCL branch through `kcl_solve` from `V0`, converged near `lift(m, eta)`),
generic in the number type.
"""
function completion_energy(sc::StorageCompletion, eta::AbstractVector; V0::AbstractVector = sc.m.V0)
    m = sc.m
    sys = m.sys
    x = lift(m, eta)
    V = kcl_solve(sys, x, V0)
    phys = shifted_strain_energy(sys, x, V, sc.ref)
    for f in sc.field
        dz = x[f.z] .- m.x0[f.z]
        phys -= (x[f.xe] - m.x0[f.xe]) * dot(f.MBf, dz)
    end
    ctrl = dot(eta, sc.Fc * eta) / 2
    rep = dot(eta, sc.Fr * eta) / 2
    for t in sc.terms
        t.coefficient == 0 && continue
        if t.kind === :load_filter
            c = sys.comps[t.k]
            b = sys.inj[t.k][1]
            z = x[first(_state_range(sys, t.k))]
            rep += t.coefficient * (sqrt(V[2b-1]^2 + V[2b]^2) - c.p.Vini_c - z)^2 / 2
        elseif t.kind === :load_z
            i = first(_state_range(sys, t.k))
            rep += t.coefficient * (x[i] - m.x0[i])^2 / 2
        elseif t.kind === :delta_omega
            r_ = _state_range(sys, t.k)
            rep += t.coefficient * (x[r_[1]] - m.x0[r_[1]]) * (x[r_[2]] - m.x0[r_[2]])
        end
    end
    return (physical = phys, controller = ctrl, repair = rep, total = phys + ctrl + rep)
end

"""
    completion_hessian_check(sc, P) -> NamedTuple

The Hessian of `completion_energy(sc, .).total` at the equilibrium against `P` (the LMI's
matrix in unscaled section coordinates): the relative difference, and the gradient norm.
"""
function completion_hessian_check(sc::StorageCompletion, P::AbstractMatrix)
    n = neta(sc.m)
    f(e) = completion_energy(sc, e).total
    g = ForwardDiff.gradient(f, zeros(n))
    Hs = ForwardDiff.hessian(f, zeros(n))
    return (relative = maximum(abs, Hs .- P) / maximum(abs, P), gradient_norm = norm(g))
end
