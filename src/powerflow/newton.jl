# Newton-Raphson AC power flow in polar coordinates.
#
# Port of PHPS powerflow.py (commit ba11ea1): the same bus classification, specifications,
# unknown ordering (angles of all non-slack buses, then magnitudes of PQ buses), Jacobian,
# convergence test (max |mismatch| < tol, checked before each update) and defaults
# (tol = 1e-6, 20 iterations). PHPS applies no reactive-power limits, so neither does this.

@enum BusType::Int8 PQ_BUS = 0 PV_BUS = 1 SLACK_BUS = 2

"""
    PowerFlowSpec

Power-flow data in network bus order: bus types, net injection specifications
(generation minus load) and the starting point (flat, with PV/slack magnitudes and the slack
angle from the case).
"""
struct PowerFlowSpec
    net::Network
    Y::Matrix{ComplexF64}
    types::Vector{BusType}
    P_spec::Vector{Float64}
    Q_spec::Vector{Float64}
    V0::Vector{Float64}
    theta0::Vector{Float64}
end

function PowerFlowSpec(case::Case; net::Network = Network(case))
    n = nbus(net)
    Y = Matrix(ybus_pf(case; net))
    types = fill(PQ_BUS, n)
    V = ones(n)
    th = zeros(n)
    P = zeros(n)
    Q = zeros(n)
    for s in case.slack
        haskey(net.index, s.bus) || continue
        i = net.index[s.bus]
        types[i] = SLACK_BUS
        V[i] = s.v0
        th[i] = s.a0
    end
    for g in case.pv
        haskey(net.index, g.bus) || continue
        i = net.index[g.bus]
        types[i] == SLACK_BUS && continue
        types[i] = PV_BUS
        V[i] = g.v0
        P[i] += g.p0
    end
    for l in case.pq
        haskey(net.index, l.bus) || continue
        i = net.index[l.bus]
        P[i] -= l.p0
        Q[i] -= l.q0
    end
    return PowerFlowSpec(net, Y, types, P, Q, V, th)
end

"""
    PowerFlowResult
"""
struct PowerFlowResult
    spec::PowerFlowSpec
    V::Vector{Float64}
    theta::Vector{Float64}
    converged::Bool
    iterations::Int          # Newton updates taken (the PHPS "Converged in N iterations")
    mismatch::Float64        # max |mismatch| at the returned point (NaN if not evaluated)
    from_case::Bool          # true if the voltages came from Bus v0/a0 (skip_pf_solve)
end

Base.show(io::IO, r::PowerFlowResult) =
    print(io, "PowerFlowResult(", nbus(r.spec.net), " buses, ",
          r.from_case ? "from case v0/a0" :
          r.converged ? "converged in $(r.iterations) iterations, mismatch $(r.mismatch)" :
          "NOT converged", ")")

"""
    bus_power(Y, V, theta) -> Vector{Complex}

Complex power injected at each bus, `S = V .* conj(Y * V)`.
"""
function bus_power(Y::AbstractMatrix, V::AbstractVector, theta::AbstractVector)
    Vc = V .* cis.(theta)
    return Vc .* conj.(Y * Vc)
end

"""
    solve_powerflow(case; tol = 1e-6, max_iter = 20) -> PowerFlowResult

Solve the power flow. If the case sets `config.skip_pf_solve`, the voltages are taken from
`Bus` `v0`/`a0` instead (PHPS `load_bus_overrides`).
"""
function solve_powerflow(case::Case; tol::Real = 1e-6, max_iter::Integer = 20,
                         net::Network = Network(case))
    spec = PowerFlowSpec(case; net)
    if case.config.skip_pf_solve
        V, th = _bus_overrides(case, spec)
        return PowerFlowResult(spec, V, th, true, 0, NaN, true)
    end
    return solve_powerflow(spec; tol, max_iter)
end

# PHPS matches Bus entries by name against "BUS<id>", "<id>" and the id itself.
function _bus_overrides(case::Case, spec::PowerFlowSpec)
    V, th = copy(spec.V0), copy(spec.theta0)
    byname = Dict{String,Int}()
    for (id, i) in spec.net.index
        byname["BUS$id"] = i
        byname[string(id)] = i
    end
    for (k, b) in enumerate(case.buses)
        i = get(byname, b.name, nothing)
        i === nothing && continue
        raw = case.raw[:Bus][k]
        haskey(raw, :v0) && (V[i] = b.v0)
        haskey(raw, :a0) && (th[i] = b.a0)
    end
    return V, th
end

function solve_powerflow(spec::PowerFlowSpec; tol::Real = 1e-6, max_iter::Integer = 20)
    n = nbus(spec.net)
    Y = spec.Y
    V = copy(spec.V0)
    th = copy(spec.theta0)

    theta_idx = [i for i in 1:n if spec.types[i] != SLACK_BUS]
    v_idx = [i for i in 1:n if spec.types[i] == PQ_BUS]
    nt, nv = length(theta_idx), length(v_idx)
    trow = zeros(Int, n)               # position of bus i among the unknowns, 0 if none
    vrow = zeros(Int, n)
    for (k, i) in enumerate(theta_idx)
        trow[i] = k
    end
    for (k, i) in enumerate(v_idx)
        vrow[i] = nt + k
    end
    G, B = real(Y), imag(Y)
    m = zeros(nt + nv)
    J = zeros(nt + nv, nt + nv)

    for it in 0:(max_iter - 1)
        S = bus_power(Y, V, th)
        Pc, Qc = real(S), imag(S)
        for (k, i) in enumerate(theta_idx)
            m[k] = spec.P_spec[i] - Pc[i]
        end
        for (k, i) in enumerate(v_idx)
            m[nt + k] = spec.Q_spec[i] - Qc[i]
        end
        # PHPS leaves the loop without success when there is nothing to solve.
        isempty(m) && break
        norm_inf = maximum(abs, m)
        if norm_inf < tol
            return PowerFlowResult(spec, V, th, true, it, norm_inf, false)
        end

        fill!(J, 0.0)
        for i in 1:n
            if trow[i] > 0
                r = trow[i]
                for j in 1:n
                    c = trow[j]
                    c > 0 || continue
                    if i == j
                        J[r, c] = -Qc[i] - V[i]^2 * B[i, i]
                    else
                        tij = th[i] - th[j]
                        J[r, c] = V[i] * V[j] * (G[i, j] * sin(tij) - B[i, j] * cos(tij))
                    end
                end
                for j in 1:n
                    c = vrow[j]
                    c > 0 || continue
                    if i == j
                        J[r, c] = Pc[i] / V[i] + V[i] * G[i, i]
                    else
                        tij = th[i] - th[j]
                        J[r, c] = V[i] * (G[i, j] * cos(tij) + B[i, j] * sin(tij))
                    end
                end
            end
            if vrow[i] > 0
                r = vrow[i]
                for j in 1:n
                    c = trow[j]
                    c > 0 || continue
                    if i == j
                        J[r, c] = Pc[i] - V[i]^2 * G[i, i]
                    else
                        tij = th[i] - th[j]
                        J[r, c] = -V[i] * V[j] * (G[i, j] * cos(tij) + B[i, j] * sin(tij))
                    end
                end
                for j in 1:n
                    c = vrow[j]
                    c > 0 || continue
                    if i == j
                        J[r, c] = Qc[i] / V[i] - V[i] * B[i, i]
                    else
                        tij = th[i] - th[j]
                        J[r, c] = V[i] * (G[i, j] * sin(tij) - B[i, j] * cos(tij))
                    end
                end
            end
        end

        dx = J \ m
        for (k, i) in enumerate(theta_idx)
            th[i] += dx[k]
        end
        for (k, i) in enumerate(v_idx)
            V[i] += dx[nt + k]
        end
    end
    return PowerFlowResult(spec, V, th, false, max_iter, NaN, false)
end
