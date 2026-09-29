# Equilibrium (roadmap P6): power flow -> component initialisation -> one Newton solve of the
# full DAE for (x*, V*), in place of PHPS's multi-pass refinement.
#
# Unknowns: every state except delta_COI (free: its derivative is zero at any value), the
# reservoir levels (one-way accounts: nothing depends on them) and, when no bus has a fixed
# voltage, one rotor angle (the network equations are invariant under a common rotation of
# all angles; PHPS's reference is the power-flow slack angle); plus every bus voltage.
# Equations: every residual row except the reservoir rows. One row is dependent (the
# rotation), so each step is a least-squares (QR) solve of a consistent system.
#
# Parameters that fix the operating point come from the power flow via `first_pass`. After
# the solve, as in PHPS: load Vini is set to the solved |V|, governor Pref to its rest-state
# value, and the solve is repeated; finally each reservoir reference is set to the power it
# actually supplies, so the reservoirs are at rest.

"""
    EquilibriumResult

`sys` (assembled with the initialised parameters), the equilibrium `x`, `V`, the initial
guess `x_initial` from the first pass, the max-abs residual, and the Newton iterations.
"""
struct EquilibriumResult
    sys::DAESystem
    x::Vector{Float64}
    V::Vector{Float64}
    x_initial::Vector{Float64}
    init_params::Dict{String,Dict{String,Float64}}
    residual::Float64
    iterations::Int
end

Base.show(io::IO, r::EquilibriumResult) =
    print(io, "EquilibriumResult(", r.sys.n_diff, " states, residual ", r.residual, ", ",
          r.iterations, " Newton iterations)")

const RESERVOIR_STATES = ("x_field", "x_steam", "x_water")

# the state index held fixed for the rotation gauge, or 0 when a bus voltage is fixed
function _gauge_state(sys::DAESystem)
    any(sys.slack) && return 0
    slack_buses = Set(s.bus for s in sys.case.slack)
    for pass in (true, false), (k, c) in enumerate(sys.comps)
        component_role(c) === :generator || continue
        pass && !(bus(c) in slack_buses) && continue
        state_names(c)[1] == "delta" && return sys.offsets[k]
    end
    return 0
end

function _reservoir_states(sys::DAESystem)
    out = Int[]
    for (k, c) in enumerate(sys.comps), (j, s) in enumerate(state_names(c))
        s in RESERVOIR_STATES && push!(out, sys.offsets[k] + j - 1)
    end
    return out
end

"""
    component_io(sys, x, V) -> (inputs, outputs)

Every component's inputs and outputs at `(x, V)`, as the residual evaluates them (after the
step pass).
"""
function component_io(sys::DAESystem, x, V)
    nb = nbus(sys)
    Vd = [V[2i - 1] for i in 1:nb]
    Vq = [V[2i] for i in 1:nb]
    T = promote_type(eltype(x), eltype(V))
    outs = [zeros(T, noutputs(c)) for c in sys.comps]
    ins = [zeros(T, ninputs(c)) for c in sys.comps]
    for (k, c) in enumerate(sys.comps)
        o = sys.offsets[k]
        _gather!(ins[k], k, T, Vd, Vq, x, sys, outs)
        _outputs!(outs[k], c, view(x, o:o + nstates(c) - 1), ins[k], params(c), NoModes())
    end
    for (k, c) in enumerate(sys.comps)
        o = sys.offsets[k]
        n = nstates(c)
        _gather!(ins[k], k, T, Vd, Vq, x, sys, outs)
        _step!(zeros(T, n), outs[k], c, view(x, o:o + n - 1), ins[k], params(c), NoModes())
    end
    return ins, outs
end

# Inner solve: every residual row except the reservoir rows and the delta_COI row, so the
# common frequency is free (it settles where the set-points balance the power).
function _newton!(y, sys::DAESystem; tol, maxiter)
    nd = sys.n_diff
    res = _reservoir_states(sys)
    fixed = Set(res)
    push!(fixed, sys.delta_coi)
    g0 = _gauge_state(sys)
    g0 > 0 && push!(fixed, g0)
    cols = [j for j in 1:length(y) if !(j in fixed)]
    rows = [r for r in 1:length(y) if !(r in res) && r != sys.delta_coi]
    F(z) = vcat(dae_residual(sys, view(z, 1:nd), view(z, nd + 1:length(z)))...)
    it = 0
    r = F(y)[rows]
    while maximum(abs, r) > tol && it < maxiter
        J = ForwardDiff.jacobian(F, y)[rows, cols]
        # equilibrate the rows so equations with small coefficients (e.g. a swing row of a
        # large-inertia machine, scaled by 1/2H) weigh as much as the others
        s = [1.0 / max(maximum(abs, view(J, i, :)), 1e-300) for i in axes(J, 1)]
        y[cols] .-= qr(s .* J, ColumnNorm()) \ (s .* r)
        it += 1
        r_new = F(y)[rows]
        # stop when rounding, not the linearisation, limits the residual
        maximum(abs, r_new) >= 0.5 * maximum(abs, r) && (r = r_new; break)
        r = r_new
    end
    return it
end

# The set-point that closes the power balance (like the power-flow slack): the Pref of the
# governor of the machine at the slack bus, or that machine's Tm0 if it has no governor.
function _slack_setpoint(sys::DAESystem)
    slack_buses = Set(s.bus for s in sys.case.slack)
    k = findfirst(c -> component_role(c) === :generator && bus(c) in slack_buses, sys.comps)
    k === nothing && (k = findfirst(c -> component_role(c) === :generator, sys.comps))
    k === nothing && return nothing
    gen = name(sys.comps[k])
    for (j, c) in enumerate(sys.comps)
        component_role(c) === :governor || continue
        s = sys.sources[j][findfirst(==("omega"), input_names(c))]
        s.kind === SRC_OUTPUT && s.index == k && return (name(c), "Pref")
    end
    return (gen, "Tm0")
end

_coi_omega(sys::DAESystem, x) = isempty(sys.coi_members) ? 1.0 :
    sum(w * x[sys.offsets[k] + 1] for (w, k) in zip(sys.coi_weights, sys.coi_members)) /
    sys.coi_total

"""
    solve_equilibrium(case[, scenario]; tol = 1e-13, maxiter = 20, phps_rounding = true)
        -> EquilibriumResult

Initialise the case (power flow, component initialisation) and solve the full DAE for its
equilibrium `(x*, V*)`. The set-points come from the power flow; the slack machine's
set-point (its governor's `Pref`, or its `Tm0`) is adjusted so the common frequency is
exactly nominal, as the power-flow slack absorbs the losses.
"""
function solve_equilibrium(case::Case, scenario::Union{Nothing,Scenario} = nothing;
                           tol::Real = 1e-13, maxiter::Integer = 20,
                           phps_rounding::Bool = true)
    pf = solve_powerflow(case)
    x_init, init, _ = first_pass(case; pf)
    sys = assemble(case, scenario; init_params = init, phps_rounding)
    nb = nbus(sys)
    nd = sys.n_diff
    V = zeros(2nb)
    for i in 1:nb
        V[2i - 1] = pf.V[i] * cos(pf.theta[i])
        V[2i] = pf.V[i] * sin(pf.theta[i])
    end
    y = vcat(copy(x_init), V)
    iters = _newton!(y, sys; tol, maxiter)

    slack = _slack_setpoint(sys)
    rebuild() = (sys = assemble(case, scenario; init_params = init, phps_rounding))
    function set_vini!()
        change = 0.0
        for c in sys.comps
            component_role(c) === :load || continue
            i = sys.net.index[bus(c)]
            v = hypot(y[nd + 2i - 1], y[nd + 2i])
            change = max(change, abs(init[name(c)]["Vini"] - v))
            init[name(c)]["Vini"] = v
        end
        return change
    end
    # frequency error as a function of the slack set-point
    eps() = _coi_omega(sys, y) - 1.0
    dedp = 0.0
    if slack !== nothing
        cname, key = slack
        P0 = init[cname][key]
        e0 = eps()
        h = 1e-6 * max(abs(P0), 1.0)
        init[cname][key] = P0 + h
        rebuild()
        yh = copy(y)
        iters += _newton!(yh, sys; tol, maxiter)
        dedp = (_coi_omega(sys, yh) - 1.0 - e0) / h
        init[cname][key] = P0
        rebuild()
    end
    for outer in 1:20
        vchange = set_vini!()
        e = slack === nothing ? 0.0 : eps()
        if slack !== nothing && dedp != 0.0 && abs(e) > 0.0
            init[slack[1]][slack[2]] -= e / dedp
        end
        rebuild()
        iters += _newton!(y, sys; tol, maxiter)
        (vchange == 0.0 && abs(eps()) * sys.omega_b <= tol) && break
    end
    # reservoir references: the power each reservoir supplies at the equilibrium
    x, Vs = y[1:nd], y[nd + 1:end]
    ins, outs = component_io(sys, x, Vs)
    for (k, c) in enumerate(sys.comps)
        o = sys.offsets[k]
        if c isa IEEET1_PHTRUE
            init[name(c)]["PFD_REF"] = x[o + 2] * ins[k][4]      # Efd * i_fd
        elseif c isa IEEEG1_PHTRUE || c isa IEEEG3_PHTRUE
            init[name(c)]["PM_REF"] = outs[k][1]                 # Tm
        end
    end
    rebuild()
    f, g = dae_residual(sys, x, Vs)
    res = max(maximum(abs, f), maximum(abs, g))
    return EquilibriumResult(sys, x, Vs, x_init, init, res, iters)
end
