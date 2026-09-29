# Physical-state projection and the local shifted-storage audit (roadmap P10; PHPS
# `certification.py` and `study/src/audit_hs_decay39.py`).
#
# The audit asks whether the assembled storage H can serve as a Lyapunov function near the
# equilibrium x* of the network-reduced field f:
#   - positivity: the Hessian S of H at x* on the physical coordinates;
#   - decay: dH_s/dt = (grad H(x) - grad H(x*))' f(x) has the quadratic part z' M z,
#     M = sym(S A), A = Df(x*), on the invariant section of the common angle;
#   - and, along the top eigenvector of M, the exact nonlinear dH_s/dt (a positive value is
#     a concrete counterexample to monotone decay).
# It is a diagnostic that can falsify a candidate storage, never a proof.

"""
    PhysicalProjection

The physical coordinates of the full state: the one-way reservoir states (from the port
contracts), the identically held states and the `delta_COI` monitor are removed (`keep`);
the COI covector `coi` (weight on each COI source's angle state) defines the invariant
section, whose orthonormal basis in `keep` coordinates is `basis`.
"""
struct PhysicalProjection
    names::Vector{String}
    reservoir::Vector{Int}
    held::Vector{Int}
    monitor::Vector{Int}
    keep::Vector{Int}
    coi::Vector{Float64}
    basis::Matrix{Float64}
end

Base.show(io::IO, p::PhysicalProjection) =
    print(io, "PhysicalProjection(", length(p.names), " states: ", length(p.reservoir),
          " reservoirs, ", length(p.held), " held, ", length(p.monitor), " monitor; ",
          length(p.keep), " kept, section ", size(p.basis, 2), ")")

"""
    reservoir_states(sys; contracts = load_contracts()) -> Vector{Int}

The one-way accounting states: each component's reservoir state as its port contract
declares it. A component type without a contract is an error (fail closed).
"""
function reservoir_states(sys::DAESystem; contracts::ContractSet = default_contracts())
    out = Int[]
    for (k, c) in enumerate(sys.comps)
        entry = contract(contracts, model_type(c))
        res = get(entry.raw, :reservoir, nothing)
        res isa JSON3.Object || continue
        s = String(res[:state])
        j = findfirst(==(s), state_names(c))
        j === nothing && error("$(name(c)): contract reservoir state $s is not a state")
        push!(out, sys.offsets[k] + j - 1)
    end
    return out
end

"""
    physical_projection(sys, x, V; contracts, probes = 3, scale = 1e-3, seed = 7)

Held states are those whose reduced-field row is exactly zero at `probes` states
`x + scale * randn` (PHPS's rule). The COI weights are `sys.coi_weights`.
"""
function physical_projection(sys::DAESystem, x::AbstractVector, V::AbstractVector;
                             contracts::ContractSet = default_contracts(), probes::Integer = 3,
                             scale::Real = 1e-3, seed::Integer = 7)
    n = sys.n_diff
    res = reservoir_states(sys; contracts)
    mon = [sys.delta_coi]
    rng = Random.Xoshiro(seed)
    rows = [reduced_field(sys, x .+ scale .* randn(rng, n), V)[1] for _ in 1:probes]
    held = [i for i in 1:n if !(i in res) && !(i in mon) && all(r -> r[i] == 0.0, rows)]
    keep = [i for i in 1:n if !(i in res) && !(i in mon) && !(i in held)]
    coi = zeros(n)
    for (m, w) in zip(sys.coi_members, sys.coi_weights)
        c = sys.comps[m]
        j = findfirst(in(("delta", "theta")), state_names(c))
        j === nothing && error("$(name(c)): COI source without a delta/theta state")
        coi[sys.offsets[m] + j - 1] += w
    end
    basis = nullspace(reshape(coi[keep], 1, :))
    return PhysicalProjection(state_names(sys), res, held, mon, keep, coi, basis)
end

"""
    shifted_storage_audit(sys, x, V; projection, eps, rank_tol = 1e-9, rate_tol = 1e-12)
        -> Dict

The local audit of the assembled storage at the equilibrium `(x, V)`:

- `storage_rank`: rank of `S = Hess H` on the kept coordinates (eigenvalues above
  `rank_tol * max(1, max|S|)`), and the eigenvalues of `S` on the section;
- `sym_SA`: the eigenvalues of `M = U' sym(S A) U` on the section (`U` = section basis,
  `A` the exact reduced Jacobian), the number above `rate_tol * max(1, max|eig|)`;
- `exact`: `H_s` and `dH_s/dt` at `x + sign * eps * v` along the top eigenvector `v` of `M`
  (lifted, normalised), with `dH_s/dt` split by storage component.
"""
function shifted_storage_audit(sys::DAESystem, x::AbstractVector, V::AbstractVector;
                               projection::PhysicalProjection = physical_projection(sys, x, V),
                               eps = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4), rank_tol::Real = 1e-9,
                               rate_tol::Real = 1e-12)
    p = projection
    keep, U = p.keep, p.basis
    A = reduced_jacobian(sys, x, V)[keep, keep]
    S = hessian_total_hamiltonian(sys, x)[keep, keep]
    eS = eigvals(Symmetric(S))
    rank = count(>(rank_tol * max(1.0, maximum(abs, S))), abs.(eS))
    Ssec = eigvals(Symmetric(U' * S * U))
    M = Symmetric(U' * (0.5 .* (S * A .+ A' * S)) * U)
    w, W = eigen(M)
    npos = count(>(rate_tol * max(1.0, maximum(abs, w))), w)
    zero_curv = [p.names[i] for (l, i) in enumerate(keep) if abs(S[l, l]) <= rank_tol * max(1.0, maximum(abs, S))]

    g0 = grad_total_hamiltonian(sys, x)
    H0 = total_hamiltonian(sys, x)
    v = U * W[:, end]
    v ./= norm(v)
    comps = storage_components(sys)
    exact = Dict{String,Any}[]
    for e in eps, sgn in (1.0, -1.0)
        xe = collect(float(x))
        xe[keep] .+= sgn * e .* v
        f, _ = reduced_field(sys, xe, V)
        dg = grad_total_hamiltonian(sys, xe) .- g0
        per = Dict(name(sys.comps[k]) => sum(dg[_state_range(sys, k)] .* f[_state_range(sys, k)])
                   for k in comps)
        push!(exact, Dict{String,Any}("eps" => e, "sign" => sgn,
                                      "Hs" => total_hamiltonian(sys, xe) - H0 - dot(g0, xe .- x),
                                      "dHs_dt" => dot(dg, f), "by_component" => per))
    end
    return Dict{String,Any}(
        "n_states" => sys.n_diff,
        "reservoir_states" => p.names[p.reservoir],
        "held_states" => p.names[p.held],
        "monitor_states" => p.names[p.monitor],
        "n_keep" => length(keep), "n_section" => size(U, 2),
        "storage_components" => [name(sys.comps[k]) for k in comps],
        "storage_rank" => rank,
        "storage_section_eig_min" => first(Ssec), "storage_section_eig_max" => last(Ssec),
        "zero_curvature_states" => zero_curv,
        "sym_SA_eig_min" => first(w), "sym_SA_eig_max" => last(w),
        "sym_SA_positive" => npos,
        "exact_along_top_eigenvector" => exact,
        "top_eigenvector" => Dict(p.names[keep[l]] => v[l] for l in eachindex(v) if abs(v[l]) > 1e-3),
    )
end
