# Fixed-step BDF1 (backward Euler) for the DAE (roadmap P7).
#
# Port of PHPS's `solve_bdf1` (DiracCompiler._emit_main at ba11ea1), so that the P7 gate can
# compare "the same scheme, the same Newton tolerance, the same switching times":
#
#   step k -> t_new = (k+1) dt; fault q is on for t_start_q + dt/4 < t_new <= t_end_q + dt/4
#   Newton on G(y) = F(t_new, y, (y - y_old)/dt), F = [ydot - f(x, V); g(x, V)]:
#     stop when ||G||_2 < newton_tol (1e-8), at most 40 iterations;
#     Jacobian by forward differences, h_j = 1e-7 (1 + |y_j|);
#     dense LU with partial pivoting (PHPS's lu_solve);
#     trust region: the whole step is scaled so no bus voltage moves more than
#     dv_max (nit + 1) (dv_max = 0.15; <= 0 disables it);
#     backtracking line search (up to 5 halvings; one more halving if none improved);
#     a step whose Newton did not converge is accepted and counted.
#   The state is logged every round(log_dt/dt) steps, before the step, and at the end.
#   The run stops if |V| at the first bus exceeds 5 pu or is NaN (checked every 1 s).

"""
    SimResult

A simulation record: logged times `t`, logged states `Y` (one column per time,
`[x; Vd_1, Vq_1, ...]`), the system, solver statistics, and the settings the integrator
actually used (defaults filled in), for `run.json`.
"""
struct SimResult
    sys::DAESystem
    t::Vector{Float64}
    Y::Matrix{Float64}
    method::Symbol
    nonconverged::Int
    first_nonconverged_t::Float64
    stopped_early::Bool
    wall_time::Float64
    settings::Dict{String,Any}
end

Base.show(io::IO, r::SimResult) =
    print(io, "SimResult(", r.method, ", ", length(r.t), " records to t = ",
          isempty(r.t) ? 0.0 : r.t[end], ", ", r.nonconverged, " non-converged steps",
          r.stopped_early ? ", STOPPED EARLY" : "", ", ", round(r.wall_time; digits = 2), " s)")

# F(t, y, ydot) = [ydot - f; g] into res
function _dae_F!(res, sys::DAESystem, ws::DAEWorkspace, y, ydot, active)
    nd = sys.n_diff
    n = length(y)
    f = view(res, 1:nd)
    g = view(res, nd + 1:n)
    dae_residual!(f, g, sys, ws, view(y, 1:nd), view(y, nd + 1:n), active)
    for i in 1:nd
        res[i] = ydot[i] - res[i]
    end
    return res
end

# ||r||_2 summed in index order, as PHPS's C++ does
function _norm2_seq(r)
    s = 0.0
    for v in r
        s += v * v
    end
    return sqrt(s)
end

# PHPS lu_solve: Gaussian elimination with partial pivoting, in place, on A (element (i, j)
# is stored at At[j, i], so PHPS's row operations run down contiguous columns) and b.
function _phps_lu_solve!(At::Matrix{Float64}, b::Vector{Float64})
    n = length(b)
    @inbounds for k in 1:n
        pivot = k
        pmax = abs(At[k, k])
        for i in k + 1:n
            if abs(At[k, i]) > pmax
                pmax = abs(At[k, i])
                pivot = i
            end
        end
        if pivot != k
            for j in 1:n
                At[j, k], At[j, pivot] = At[j, pivot], At[j, k]
            end
            b[k], b[pivot] = b[pivot], b[k]
        end
        akk = At[k, k]
        abs(akk) < 1e-30 && (akk = 1e-30)
        for i in k + 1:n
            factor = At[k, i] / akk
            for j in k + 1:n
                At[j, i] -= factor * At[j, k]
            end
            b[i] -= factor * b[k]
        end
    end
    @inbounds for i in n:-1:1
        s = b[i]
        for j in i + 1:n
            s -= At[j, i] * b[j]
        end
        aii = At[i, i]
        abs(aii) < 1e-30 && (aii = 1e-30)
        b[i] = s / aii
    end
    return b
end

"""
    simulate_bdf1(sys, y0; dt, duration, log_dt = dt, dv_max = 0.15, newton_tol = 1e-8,
                  max_newton = 40, eps_fd = 1e-7) -> SimResult

Integrate the DAE from `y0 = [x0; V0]` with PHPS's fixed-step BDF1 (see the file header).
"""
function simulate_bdf1(sys::DAESystem, y0::AbstractVector{Float64}; dt::Real, duration::Real,
                       log_dt::Real = dt, dv_max::Real = 0.15, newton_tol::Real = 1e-8,
                       max_newton::Integer = 40, eps_fd::Real = 1e-7)
    wall = time()
    nd = sys.n_diff
    n = length(y0)
    n == nd + nalg(sys) || throw(DimensionMismatch("y0 has length $n, expected $(nd + nalg(sys))"))
    ws = DAEWorkspace(sys)
    y = collect(Float64, y0)
    y_old = similar(y)
    ydot = similar(y)
    res = similar(y)
    res_pert = similar(y)
    y_pert = similar(y)
    ydot_pert = similar(y)
    dy = similar(y)
    Jt = zeros(n, n)                    # Jt[j, i] = dG_i / dy_j
    active = falses(length(sys.faults))
    steps = round(Int, duration / dt)
    log_every = max(1, round(Int, log_dt / dt))
    check_every = trunc(Int, 1.0 / dt)

    ts = Float64[]
    cols = Vector{Vector{Float64}}()
    t = 0.0
    last_logged = -1.0
    nonconv = 0
    first_bad = -1.0
    stopped = false
    event_eps = 0.25 * dt

    for step in 0:steps - 1
        if step % log_every == 0
            push!(ts, t)
            push!(cols, copy(y))
            last_logged = t
        end
        y_old .= y
        t_new = (step + 1) * dt
        for (q, fs) in enumerate(sys.faults)
            active[q] = t_new > fs.t_start + event_eps && t_new <= fs.t_end + event_eps
        end

        res_norm = 0.0
        for nit in 0:max_newton - 1
            @. ydot = (y - y_old) / dt
            _dae_F!(res, sys, ws, y, ydot, active)
            res_norm = _norm2_seq(res)
            res_norm < newton_tol && break

            # forward-difference Jacobian
            for j in 1:n
                y_pert .= y
                h = eps_fd * (1.0 + abs(y[j]))
                y_pert[j] += h
                @. ydot_pert = (y_pert - y_old) / dt
                _dae_F!(res_pert, sys, ws, y_pert, ydot_pert, active)
                for i in 1:n
                    Jt[j, i] = (res_pert[i] - res[i]) / h
                end
            end
            @. dy = -res
            _phps_lu_solve!(Jt, dy)

            # trust region on the bus voltages
            if dv_max > 0.0
                radius = dv_max * (nit + 1)
                vstep = 0.0
                for i in nd + 1:n
                    abs(dy[i]) > vstep && (vstep = abs(dy[i]))
                end
                if vstep > radius
                    s = radius / vstep
                    dy .*= s
                end
            end
            # backtracking line search
            alpha = 1.0
            improved = false
            for ls in 0:4
                @. y_pert = y + alpha * dy
                @. ydot_pert = (y_pert - y_old) / dt
                _dae_F!(res_pert, sys, ws, y_pert, ydot_pert, active)
                if _norm2_seq(res_pert) < res_norm
                    improved = true
                    break
                end
                alpha *= 0.5
            end
            improved || (alpha *= 0.5)
            @. y += alpha * dy
        end
        if res_norm >= newton_tol
            nonconv += 1
            first_bad < 0 && (first_bad = t_new)
        end
        t = t_new

        if step % check_every == 0
            V1 = hypot(y[nd + 1], y[nd + 2])
            if V1 > 5.0 || isnan(V1)
                stopped = true
                break
            end
        end
    end
    if t > last_logged + 1e-12
        push!(ts, t)
        push!(cols, copy(y))
    end
    settings = Dict{String,Any}("dt" => dt, "duration" => duration, "log_dt" => log_dt,
                                "dv_max" => dv_max, "newton_tol" => newton_tol,
                                "max_newton" => max_newton, "eps_fd" => eps_fd)
    return SimResult(sys, ts, reduce(hcat, cols), :bdf1, nonconv, first_bad, stopped,
                     time() - wall, settings)
end
