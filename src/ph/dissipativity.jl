# Dissipativity across ports: open several component loops at once, test whether the rest
# of the system can pay the components' passivity shortages, and construct storages
# (roadmap P10 / Part II B4 item 1: joint storage by passivity indices).
#
# Supply rates, for a component k with port input u_k (e.g. speed) and output y_k (e.g. Tm):
#   component: dV_k/dt <= -y_k u_k + d_k u_k^2      (shortage d_k, H_k = -G_k)
#   rest:      dV_r/dt <= sum_k y_k u_k - n_k u_k^2  (the rest's input is y_k, output u_k)
# With n_k >= d_k the sum is <= -sum (n_k - d_k) u_k^2 <= 0: V_r + sum V_k is a joint storage,
# and the dynamics are untouched (the split is only in the certificate).

"""
    MultiPortModel

The rest of the system with several components removed: `x' = A x + B y`, `u = C x`, where
`y` stacks the removed components' outputs (the rest's inputs) and `u` their inputs (the
rest's outputs), on the common-angle section of the remaining states.
"""
struct MultiPortModel
    components::Vector{String}
    input::String
    output::String
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
end

Base.show(io::IO, m::MultiPortModel) =
    print(io, "MultiPortModel(rest without ", length(m.components), " components: ",
          m.output, " -> ", m.input, ", ", size(m.A, 1), " states)")

"""
    open_loops_model(sys, x, V, ks; input, output, projection, contracts)
        -> (rest::MultiPortModel, ports::Vector{PortModel}, T::Matrix)

Open the loops at the components `ks` at once (as `loop_port_model` does for one): the rest
of the system from their outputs to their inputs, and each component's own port model.
`T` maps the coordinates `[rest section; component states...]` to the physical (`keep`)
coordinates, so that `T' A T` is the whole system in those coordinates.
"""
function open_loops_model(sys::DAESystem, x::AbstractVector, V::AbstractVector,
                          ks::AbstractVector{<:Integer}; input::AbstractString,
                          output::AbstractString,
                          projection::PhysicalProjection = physical_projection(sys, x, V),
                          contracts::ContractSet = default_contracts())
    ins, _ = component_io(sys, x, V)
    keep = projection.keep
    pos = Dict(i => l for (l, i) in enumerate(keep))
    A = reduced_jacobian(sys, x, V)[keep, keep]
    tol = 1e-9 * max(1.0, maximum(abs, A))
    ports = PortModel[]
    Ss = Vector{Int}[]
    for k in ks
        c = sys.comps[k]
        r = _state_range(sys, k)
        m = port_model(c, x[r], ins[k]; input, output, contracts)
        S = [pos[i] for i in r if haskey(pos, i)]
        length(S) == length(m.states) ||
            error("$(name(c)): its kept states differ from its port model's states")
        push!(ports, m)
        push!(Ss, S)
    end
    allS = reduce(vcat, Ss)
    R = setdiff(1:length(keep), allS)
    for (j, Sj) in enumerate(Ss), (k, Sk) in enumerate(Ss)
        j != k && maximum(abs, A[Sj, Sk]) > tol &&
            error("$(ports[j].component) and $(ports[k].component) are coupled directly")
    end
    Bs, Cs = Vector{Float64}[], Vector{Float64}[]
    for (m, S) in zip(ports, Ss)
        b = A[R, S] * m.C / dot(m.C, m.C)
        cr = vec(m.B' * A[S, R]) / dot(m.B, m.B)
        maximum(abs, A[R, S] .- b * m.C') <= tol ||
            error("$(m.component) acts on the rest through more than its output $output")
        maximum(abs, A[S, R] .- m.B * cr') <= tol ||
            error("$(m.component) is driven by the rest through more than its input $input")
        push!(Bs, b)
        push!(Cs, cr)
    end
    l = projection.coi[keep][R]
    maximum(abs, l' * A[R, R]) <= tol || error("the common-angle section is not invariant")
    U = nullspace(reshape(l, 1, :))
    n, nS = size(U, 2), length(allS)
    T = zeros(length(keep), n + nS)
    T[R, 1:n] .= U
    for (j, i) in enumerate(allS)
        T[i, n + j] = 1.0
    end
    rest = MultiPortModel([m.component for m in ports], String(input), String(output),
                          U' * A[R, R] * U, U' * reduce(hcat, Bs), permutedims(reduce(hcat, Cs)) * U)
    return rest, ports, T
end

"""`G(s) = C (sI - A)^{-1} B` of a multi-port model (square: one row and column per port)."""
transfer(m::MultiPortModel, s::Number) = m.C * ((s * I - m.A) \ complex.(m.B))

"""
    multiport_margin(rest, n; eta = 0.0, ws) -> NamedTuple

The frequency condition for a storage of the rest with supply
`sum_k y_k u_k - n_k u_k^2 + eta y_k^2` (`u = G y`): the smallest eigenvalue over `ws` of
`Herm(G(jw)) - G(jw)' diag(n) G(jw) + eta I`. `passes` when it is positive (then
`rest_storage(rest, n; eta)` exists, by the KYP lemma in its strict form). As `w -> infinity`
it tends to `eta`. With `eta = 0` it is the classical constant split; `eta > 0` lets the
components pay for directions the rest cannot dissipate (the differential speed directions
at low frequency, which a lossy network makes non-passive).
"""
function multiport_margin(rest::MultiPortModel, n::AbstractVector; eta::Real = 0.0,
                          ws = exp10.(range(-3, 3; length = 20001)))
    F = hessenberg(rest.A)
    B = complex.(rest.B)
    N = Diagonal(n)
    lmin, wmin = Inf, NaN
    for w in ws
        G = rest.C * ((F - (im * w) * I) \ B)
        l = eigmin(Hermitian((G + G') / 2 - G' * N * G + eta * I))
        l < lmin && ((lmin, wmin) = (l, w))
    end
    return (margin = lmin, at = wmin, passes = lmin > 0)
end

"""
    port_margin(m, d; eta = 0.0, ws) -> NamedTuple

The frequency condition for `port_storage(m, d; eta)`: the smallest over `ws` of
`Re H(jw) + d - eta |H(jw)|^2`, `H = -G`, and its value as `w -> infinity` (`d`).
"""
function port_margin(m::PortModel, d::Real; eta::Real = 0.0,
                     ws = exp10.(range(-3, 3; length = 20001)))
    H = -frequency_response(m, ws)
    v = real.(H) .+ d .- eta .* abs2.(H)
    k = argmin(v)
    return (margin = min(v[k], d), at = v[k] < d ? ws[k] : Inf, passes = min(v[k], d) > 0)
end

"""
    kyp_riccati(A, B, Q, S, R; strict = 0.0) -> (P, residual)

A quadratic storage `V = x' P x` from the KYP inequality

    K(P) = [A'P + PA + Q + strict*I   PB - S;  (PB - S)'   -R] <= 0,   R > 0,

i.e. `dV/dt <= s(x, w) - strict |x|^2` for the supply `s = -x'Qx + 2 x'S w + w'R w` of
`x' = A x + B w`. `P` is the stabilising solution of the Riccati equation (the Schur complement
of `K` equal to zero), from the stable invariant subspace of its Hamiltonian matrix; it exists
when the frequency condition of the supply holds strictly. `residual` is the largest
eigenvalue of `K(P)`, which must be <= 0 up to rounding.
"""
function kyp_riccati(A::AbstractMatrix, B::AbstractMatrix, Q::AbstractMatrix, S::AbstractMatrix,
                     R::AbstractMatrix; strict::Real = 0.0)
    n = size(A, 1)
    Ri = inv(Symmetric(Matrix(R)))
    A1 = A - B * Ri * S'
    Ham = [A1 B * Ri * B'; -(Q + S * Ri * S' + strict * I) -A1']
    F = schur(Ham)
    ordschur!(F, real.(F.values) .< 0)
    count(<(0), real.(F.values)) == n ||
        error("kyp_riccati: eigenvalues on the imaginary axis (the frequency condition is not strict)")
    X = F.Z[:, 1:n]
    P = X[n + 1:end, :] / X[1:n, :]
    P = (P + P') / 2
    K = [A' * P + P * A + Q + strict * I  P * B - S; (P * B - S)' -R]
    return P, eigmax(Symmetric(K))
end

"""
    port_storage(m, d; eta = 0.0, strict = 0.0) -> (P, residual)

The storage of a component at its port with `dV/dt <= u y + d u^2 - eta y^2`, `y = -(C x)`
(a governor: `u` the speed, `y = -Tm`; the port model must have `D = 0`).
"""
function port_storage(m::PortModel, d::Real; eta::Real = 0.0, strict::Real = 0.0)
    m.D == 0 || throw(ArgumentError("port_storage: the port model has a feedthrough"))
    c = m.C
    return kyp_riccati(m.A, reshape(m.B, :, 1), eta .* (c * c'), reshape(-c ./ 2, :, 1),
                       fill(Float64(d), 1, 1); strict)
end

"""
    rest_storage(rest, n; eta, strict = 0.0) -> (P, residual)

The storage of the rest of the system with `dV/dt <= sum_k y_k u_k - n_k u_k^2 + eta y_k^2`
(`y` its inputs, the components' outputs; `u = C z` its outputs). `eta > 0` makes the
Riccati equation regular; the components' storages pay it back (`port_storage(...; eta)`),
so it cancels in the sum.
"""
rest_storage(rest::MultiPortModel, n::AbstractVector; eta::Real, strict::Real = 0.0) =
    kyp_riccati(rest.A, rest.B, rest.C' * Diagonal(n) * rest.C, rest.C' ./ 2,
                Matrix(eta * I, length(n), length(n)); strict)

"""
    loop_margin(rest, ports; ws) -> NamedTuple

The most generous diagonal split at the ports: each component pays exactly its own
`Re H_k(jw)` (`H_k = -G_k`) at every frequency. The loop admits a storage split port by port
(with any static or frequency-dependent diagonal multiplier) only if
`Herm(G(jw)) + G(jw)' diag(Re H_k(jw)) G(jw)` is positive definite at every frequency;
returns its smallest eigenvalue over `ws` and where. A negative value rules the port-by-port
split out.
"""
function loop_margin(rest::MultiPortModel, ports::AbstractVector{PortModel};
                     ws = exp10.(range(-3, 3; length = 20001)))
    F = hessenberg(rest.A)
    B = complex.(rest.B)
    Hs = [-frequency_response(m, ws) for m in ports]
    lmin, wmin = Inf, NaN
    for (j, w) in enumerate(ws)
        G = rest.C * ((F - (im * w) * I) \ B)
        l = eigmin(Hermitian((G + G') / 2 + G' * Diagonal([real(H[j]) for H in Hs]) * G))
        l < lmin && ((lmin, wmin) = (l, w))
    end
    return (margin = lmin, at = wmin, passes = lmin > 0)
end
