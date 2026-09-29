# The assembled DAE (roadmap P5).
#
#   x' = f(x, V)          component states, then the centre-of-inertia angle delta_COI
#   0  = g(x, V)          KCL at every bus: I_inj(x, V) - (Y + Y_fault) V - I_freq
#                         (a fixed-voltage slack bus without a machine: V - V_ref)
#
# Port of the residual PHPS compiles for IDA and BDF1 (DiracCompiler._emit_dae_residual at
# ba11ea1), evaluated in the same order:
#
#   1. every component's output kernel, in case order; the currents of the components with
#      `Id`/`Iq` outputs are summed per bus (Norton sources and load corrections);
#   2. every component's step kernel, in case order (a step may rewrite some of its own
#      outputs, which later components then read);
#   3. the centre-of-inertia frame: omega_COI = sum(2H_k omega_k) / sum(2H_k) over the swing
#      sources, each rotor angle rotates with omega_b (omega_COI - 1), and delta_COI tracks it;
#   4. KCL with the dense Y-bus, the active fault shunts and the frequency-dependent loads.
#
# PHPS writes the network constants into its C++ with printf("%.10e") (11 significant
# digits) and the COI weights with "%.6f". With `phps_rounding = true` (the default, for
# parity) the system uses the same rounded constants.

"""
    DAESystem

An assembled system: components, state layout, wiring, network constants and events.
States `x` are in PHPS order (components in case order, then `delta_COI`); bus voltages
are `V = [Vd_1, Vq_1, Vd_2, Vq_2, ...]` in network order.
"""
struct DAESystem
    case::Case
    net::Network
    comps::Vector{AbstractComponent}
    offsets::Vector{Int}                 # first state of each component (1-based)
    n_diff::Int                          # component states + delta_COI
    delta_coi::Int
    sources::Vector{Vector{InputSource}}
    inj::Vector{NTuple{3,Int}}           # (bus index, Id output, Iq output), zeros if none
    G::Matrix{Float64}                   # Y-bus, real and imaginary parts
    B::Matrix{Float64}
    load::LoadAdmittances
    slack::Vector{Bool}                  # fixed-voltage buses (slack without a machine)
    Vd_ref::Vector{Float64}
    Vq_ref::Vector{Float64}
    coi_members::Vector{Int}             # component indices of the swing sources
    coi_weights::Vector{Float64}
    coi_total::Float64
    omega_b::Float64
    faults::Vector{FaultShunt}
    pf::PowerFlowResult
end

nbus(sys::DAESystem) = nbus(sys.net)
"""Number of algebraic variables (2 per bus)."""
nalg(sys::DAESystem) = 2 * nbus(sys)

Base.show(io::IO, s::DAESystem) = print(io, "DAESystem(", length(s.comps), " components, ",
                                        s.n_diff, " states, ", nbus(s), " buses, ",
                                        length(s.faults), " faults)")

"""
    state_names(sys) -> Vector{String}

`"COMPONENT.state"` for every differential state, then `"delta_COI"`.
"""
function state_names(sys::DAESystem)
    out = String[]
    for c in sys.comps, s in state_names(c)
        push!(out, name(c) * "." * s)
    end
    push!(out, "delta_COI")
    return out
end

_r10(x::Float64, on::Bool) = on ? parse(Float64, Printf.@sprintf("%.10e", x)) : x

"""
    assemble(case[, scenario]; init_params = Dict(), phps_rounding = true) -> DAESystem

Build the DAE of a case, with the bus faults of `scenario`. `init_params` maps component
names to parameters that initialisation sets (until P6, these come from PHPS). The
constant-impedance loads use the solved power-flow voltages, as PHPS's simulation does.
"""
function assemble(case::Case, scenario::Union{Nothing,Scenario} = nothing;
                  init_params::AbstractDict = Dict{String,Any}(),
                  phps_rounding::Bool = true)
    net = Network(case)
    pf = solve_powerflow(case; net)

    pre = [component_params(case, s) for s in case.components]
    comps = AbstractComponent[]
    for (s, pd) in zip(case.components, pre)
        c = build_component(s.type, s.name, copy(pd))
        haskey(init_params, s.name) && (c = with_params(c, init_params[s.name]))
        push!(comps, c)
    end

    offsets = Int[]
    n = 0
    for c in comps
        push!(offsets, n + 1)
        n += nstates(c)
    end
    delta_coi = n + 1
    n_diff = n + 1

    sources = resolve_wiring(case, comps, net; pre_params = pre)
    for (k, c) in enumerate(comps), s in sources[k]
        if s.kind === SRC_OUTPUT && output_names(comps[s.index])[s.port] in ("id_dq", "iq_dq")
            throw(ArgumentError("$(name(c)) reads $(name(comps[s.index])).$(output_names(comps[s.index])[s.port]); " *
                                "PHPS refreshes id_dq / iq_dq with 6-decimal coefficients, not ported yet"))
        end
    end

    inj = NTuple{3,Int}[]
    for c in comps
        b = bus(c)
        outs = output_names(c)
        if b !== nothing && haskey(net.index, b) && "Id" in outs && "Iq" in outs
            push!(inj, (net.index[b], findfirst(==("Id"), outs), findfirst(==("Iq"), outs)))
        else
            push!(inj, (0, 0, 0))
        end
    end

    Y = Matrix(ybus_dae(case; net, v0 = pf.V))
    G = [_r10(real(y), phps_rounding) for y in Y]
    B = [_r10(imag(y), phps_rounding) for y in Y]
    loads = Dict(name(c) => param_dict(c) for c in comps if component_role(c) === :load)
    la = load_admittances(case; net, v0 = pf.V, load_params = loads)
    r(v) = [_r10(x, phps_rounding) for x in v]
    la = LoadAdmittances(r(la.G), r(la.B), r(la.P), r(la.Q), r(la.kpf), r(la.kqf),
                         la.has_complex_loads)

    nb = nbus(net)
    slack = falses(nb)
    Vd_ref, Vq_ref = zeros(nb), zeros(nb)
    gen_buses = Set(bus(c) for c in comps if component_role(c) === :generator)
    for s in case.slack
        haskey(net.index, s.bus) || continue
        s.bus in gen_buses && continue          # a machine at the slack: KCL
        i = net.index[s.bus]
        slack[i] = true
        V = s.v0 * cis(s.a0)
        Vd_ref[i], Vq_ref[i] = _r10(real(V), phps_rounding), _r10(imag(V), phps_rounding)
    end

    members = [k for (k, c) in enumerate(comps)
               if component_role(c) === :generator && "omega" in state_names(c)]
    weights = Float64[]
    for k in members
        pd = param_dict(comps[k])
        w = haskey(pd, "H") ? 2.0 * param_value(pd["H"]) :
            haskey(pd, "Ta") ? param_value(pd["Ta"]) :
            throw(ArgumentError("$(name(comps[k])): no inertia (H or Ta) for the COI weight"))
        push!(weights, phps_rounding ? _c6f(w) : w)
    end
    total = isempty(weights) ? 0.0 : foldl(+, weights)
    omega_b = isempty(members) ? 2.0 * π * 60.0 :
              param_value(get(param_dict(comps[members[1]]), "omega_b", "2.0 * M_PI * 60.0"))

    faults = FaultShunt[]
    if scenario !== nothing
        any(e -> e isa LineFault, scenario.events) &&
            throw(ArgumentError("LineFault is not supported by PHPS's DAE path (see TODO.md)"))
        for f in fault_shunts(net, scenario.events)
            push!(faults, FaultShunt(f.bus, f.index, _r10(f.t_start, phps_rounding),
                                     _r10(f.t_end, phps_rounding), _r10(f.g, phps_rounding),
                                     _r10(f.b, phps_rounding)))
        end
    end

    return DAESystem(case, net, comps, offsets, n_diff, delta_coi, sources, inj, G, B, la,
                     collect(slack), Vd_ref, Vq_ref, members, weights, total, omega_b, faults, pf)
end

# -- residual ---------------------------------------------------------------------------

@inline function _input(s::InputSource, T, Vd, Vq, x, sys::DAESystem, outs)
    k = s.kind
    k === SRC_ZERO && return zero(T)
    k === SRC_CONST && return T(s.value)
    k === SRC_VD && return T(Vd[s.index])
    k === SRC_VQ && return T(Vq[s.index])
    k === SRC_VTERM && return T(sqrt(Vd[s.index] * Vd[s.index] + Vq[s.index] * Vq[s.index]))
    k === SRC_OUTPUT && return T(outs[s.index][s.port])
    # dq-frame bus voltage of a machine (angle = its first state)
    c = sys.comps[s.index]
    i = sys.net.index[bus(c)]
    d = x[sys.offsets[s.index]]
    k === SRC_DQ_VD && return T(Vd[i] * sin(d) - Vq[i] * cos(d))
    return T(Vd[i] * cos(d) + Vq[i] * sin(d))
end

function _gather!(u, k::Int, T, Vd, Vq, x, sys::DAESystem, outs)
    for (j, s) in enumerate(sys.sources[k])
        u[j] = _input(s, T, Vd, Vq, x, sys, outs)
    end
    return u
end

"""
    dae_residual!(f, g, sys, x, V; faults_on = false, t = 0.0) -> (f, g)

`f = x'` (differential right-hand side, length `n_diff`) and the algebraic residual `g`
(length `2 nbus`, `[KCL_d, KCL_q]` per bus) at states `x` and bus voltages
`V = [Vd_1, Vq_1, ...]`. With `faults_on`, every fault shunt of the scenario is applied.
PHPS's residual is `[x' - f; g]`.
"""
function dae_residual!(f, g, sys::DAESystem, x, V; faults_on::Bool = false, t::Real = 0.0)
    T = promote_type(eltype(x), eltype(V))
    nb = nbus(sys)
    Vd = [V[2i - 1] for i in 1:nb]
    Vq = [V[2i] for i in 1:nb]
    comps = sys.comps
    outs = [zeros(T, noutputs(c)) for c in comps]
    ins = [zeros(T, ninputs(c)) for c in comps]
    Id_inj = zeros(T, nb)
    Iq_inj = zeros(T, nb)

    # 1. outputs and injections
    for (k, c) in enumerate(comps)
        o = sys.offsets[k]
        xc = view(x, o:o + nstates(c) - 1)
        _gather!(ins[k], k, T, Vd, Vq, x, sys, outs)
        _outputs!(outs[k], c, xc, ins[k], params(c), NoModes())
        b, jd, jq = sys.inj[k]
        if b > 0
            Id_inj[b] += outs[k][jd]
            Iq_inj[b] += outs[k][jq]
        end
    end
    # 2. dynamics
    for (k, c) in enumerate(comps)
        o = sys.offsets[k]
        n = nstates(c)
        xc = view(x, o:o + n - 1)
        _gather!(ins[k], k, T, Vd, Vq, x, sys, outs)
        _step!(view(f, o:o + n - 1), outs[k], c, xc, ins[k], params(c), NoModes())
    end
    # 3. centre-of-inertia frame
    m = sys.coi_members
    if length(m) > 1
        acc = sys.coi_weights[1] * x[sys.offsets[m[1]] + 1]
        for q in 2:length(m)
            acc += sys.coi_weights[q] * x[sys.offsets[m[q]] + 1]
        end
        coi_omega = acc / sys.coi_total
        for k in m
            f[sys.offsets[k]] -= sys.omega_b * (coi_omega - 1.0)
        end
        f[sys.delta_coi] = sys.omega_b * (coi_omega - 1.0)
    else
        coi_omega = isempty(m) ? one(T) : T(x[sys.offsets[m[1]] + 1])
        f[sys.delta_coi] = zero(T)
    end
    # 4. KCL
    Yf_g = zeros(nb)
    Yf_b = zeros(nb)
    if faults_on
        for fs in sys.faults
            Yf_g[fs.index] += fs.g
            Yf_b[fs.index] += fs.b
        end
    end
    la = sys.load
    for i in 1:nb
        if sys.slack[i]
            g[2i - 1] = Vd[i] - sys.Vd_ref[i]
            g[2i] = Vq[i] - sys.Vq_ref[i]
            continue
        end
        Id_ybus = zero(T)
        Iq_ybus = zero(T)
        for j in 1:nb
            Gij = sys.G[i, j]
            Bij = sys.B[i, j]
            Id_ybus += Gij * Vd[j] - Bij * Vq[j]
            Iq_ybus += Gij * Vq[j] + Bij * Vd[j]
        end
        Id_ybus += Yf_g[i] * Vd[i] - Yf_b[i] * Vq[i]
        Iq_ybus += Yf_g[i] * Vq[i] + Yf_b[i] * Vd[i]
        if la.kpf[i] != 0.0 || la.kqf[i] != 0.0
            dw = coi_omega - 1.0
            dP = la.kpf[i] * dw * la.G[i]
            dQ = la.kqf[i] * dw * la.B[i]
            Id_ybus += dP * Vd[i] - dQ * Vq[i]
            Iq_ybus += dP * Vq[i] + dQ * Vd[i]
        end
        g[2i - 1] = Id_inj[i] - Id_ybus
        g[2i] = Iq_inj[i] - Iq_ybus
    end
    return f, g
end

"""
    dae_residual(sys, x, V; faults_on = false) -> (f, g)
"""
function dae_residual(sys::DAESystem, x, V; kwargs...)
    T = promote_type(eltype(x), eltype(V))
    return dae_residual!(zeros(T, sys.n_diff), zeros(T, nalg(sys)), sys, x, V; kwargs...)
end
