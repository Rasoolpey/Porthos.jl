# Variable-order BDF with SUNDIALS IDA (roadmap P7).
#
# Follows PHPS's IDA driver (DiracCompiler._emit_main_ida at ba11ea1): dense linear solver,
# differential / algebraic ID vector, maximum order 3 and a maximum step of min(dt, 5 ms)
# for systems with more than 30 unknowns, tolerances rtol = 1e-4, atol = 1e-6 for such
# systems (1e-6 / 1e-8 otherwise) unless the scenario sets them. Integration runs in
# segments between the fault switching times; at each switch the fault flags change and the
# initial condition is made consistent the way IDACalcIC(IDA_YA_YDP_INIT) does it: the
# differential states are held, the bus voltages solved from KCL, and the state
# derivatives set to f. As in PHPS, the Newton for the voltages starts, at a fault
# clearing, from the voltages saved just before the fault was applied.

"""
    consistent_voltages!(V, sys, x, active; ws = DAEWorkspace(sys), tol = 1e-12, maxiter = 50)

Solve KCL `g(x, V) = 0` for the bus voltages with the states `x` held (in place, from the
guess in `V`). Returns the final max-abs residual.
"""
function consistent_voltages!(V, sys::DAESystem, x, active::AbstractVector{Bool};
                              tol::Real = 1e-12, maxiter::Integer = 50)
    nd = sys.n_diff
    f = zeros(nd)
    g = zeros(length(V))
    ws = DAEWorkspace(sys)
    G(v) = begin
        wsd = DAEWorkspace(sys, eltype(v))
        gg = zeros(eltype(v), length(v))
        dae_residual!(zeros(eltype(v), nd), gg, sys, wsd, x, v, active)
        gg
    end
    r = Inf
    for _ in 1:maxiter
        dae_residual!(f, g, sys, ws, x, V, active)
        r = maximum(abs, g)
        r <= tol && break
        J = ForwardDiff.jacobian(G, V)
        V .-= J \ g
    end
    return r
end

"""
    simulate_ida(sys, y0; dt, duration, log_dt = dt, rtol = nothing, atol = nothing)
        -> SimResult

Integrate the DAE from `y0 = [x0; V0]` with SUNDIALS IDA, restarting at every fault switch
(see the file header). Records every `log_dt` and at the end.
"""
function simulate_ida(sys::DAESystem, y0::AbstractVector{Float64}; dt::Real, duration::Real,
                      log_dt::Real = dt, rtol = nothing, atol = nothing,
                      maxiters::Integer = 5_000_000)
    wall = time()
    nd = sys.n_diff
    n = length(y0)
    big = n > 30
    rtol = rtol === nothing ? (big ? 1e-4 : 1e-6) : Float64(rtol)
    atol = atol === nothing ? (big ? 1e-6 : 1e-8) : Float64(atol)
    dtmax = big ? min(dt, 0.005) : dt

    # segment ends: the fault switching times inside (0, T), then T
    ev = Float64[]
    for fs in sys.faults
        0.0 < fs.t_start < duration && push!(ev, fs.t_start)
        0.0 < fs.t_end < duration && push!(ev, fs.t_end)
    end
    seg_ends = vcat(sort!(unique(ev)), Float64(duration))

    active = falses(length(sys.faults))
    ws = DAEWorkspace(sys)
    function F!(res, du, u, p, t)
        _dae_F!(res, sys, ws, u, du, active)
        return nothing
    end
    diffvars = [trues(nd); falses(n - nd)]

    y = collect(Float64, y0)
    ts = [0.0]
    cols = [copy(y)]
    klog = 1                                   # next log time: klog * log_dt
    t0 = 0.0
    pre_fault_V = nothing
    nsteps = 0
    stopped = false
    for t_end in seg_ends
        t_end <= t0 + 1e-12 && continue
        had = any(active)
        for (q, fs) in enumerate(sys.faults)
            active[q] = t0 >= fs.t_start && t0 < fs.t_end
        end
        has = any(active)
        V = y[nd + 1:end]
        if !had && has
            pre_fault_V = copy(V)
        elseif had && !has && pre_fault_V !== nothing
            V .= pre_fault_V
        end
        # consistent initial condition for this segment (IDACalcIC, YA_YDP)
        x = y[1:nd]
        consistent_voltages!(V, sys, x, active)
        y[nd + 1:end] .= V
        f = zeros(nd)
        dae_residual!(f, zeros(n - nd), sys, ws, x, V, active)
        du0 = [f; zeros(n - nd)]

        prob = SciMLBase.DAEProblem(F!, du0, copy(y), (t0, t_end); differential_vars = diffvars)
        alg = Sundials.IDA(linear_solver = :Dense, max_order = big ? 3 : 5)
        # log on PHPS's grid: every log_dt from the start, plus the segment end
        saves = Float64[]
        while klog * log_dt <= t_end + 1e-12
            push!(saves, min(klog * log_dt, t_end))
            klog += 1
        end
        sol = SciMLBase.solve(prob, alg; reltol = rtol, abstol = atol, dtmax = dtmax,
                              saveat = saves, save_start = true, save_end = true,
                              maxiters = maxiters, initializealg = SciMLBase.NoInit())
        # Sundials.jl can return the start time without its state when save_start = false,
        # so it is saved, and times and states are checked to pair up.
        length(sol.t) == length(sol.u) ||
            error("IDA returned $(length(sol.t)) times but $(length(sol.u)) states")
        # record exactly the planned log times (saveat points; the end of the segment is
        # among them when it falls on the grid)
        for tt in saves
            tt > ts[end] + 1e-12 || continue
            k = findfirst(s -> abs(s - tt) <= 1e-12, sol.t)
            k === nothing && continue
            push!(ts, tt)
            push!(cols, copy(sol.u[k]))
        end
        nsteps += length(sol.t)
        if !SciMLBase.successful_retcode(sol)
            stopped = true
            break
        end
        y .= sol.u[end]
        t0 = t_end
        V1 = hypot(y[nd + 1], y[nd + 2])
        (V1 > 5.0 || isnan(V1)) && (stopped = true; break)
    end
    settings = Dict{String,Any}("dt" => dt, "duration" => duration, "log_dt" => log_dt,
                                "rtol" => rtol, "atol" => atol, "dtmax" => dtmax,
                                "max_order" => big ? 3 : 5, "linear_solver" => "dense",
                                "internal_steps" => nsteps)
    return SimResult(sys, ts, reduce(hcat, cols), :ida, 0, -1.0, stopped, time() - wall,
                     settings)
end
