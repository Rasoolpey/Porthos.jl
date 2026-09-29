# Component initialisation from the power flow (PHPS Initializer.run, the first pass).
#
# Machines are initialised from their terminal phasors (V, I = conj(S/V)); they return the
# targets (Efd, Tm, Vt, i_fd, ...) their exciter and governor are then initialised from.
# Initialisation also sets parameters: Efd0 / Tm0 on the machines, Vref and the field
# reservoir reference on the exciter, Pref and the steam / water reservoir reference on the
# governors, V0 / Vini on the loads. Ports of the PHPS methods named in each function.

"""
    MachineTargets

Steady-state quantities of a machine at its power-flow operating point, used to
initialise its exciter and governor.
"""
struct MachineTargets
    Efd::Float64
    Tm::Float64
    Vt::Float64
    i_fd::Float64
    vd::Float64
    vq::Float64
    id::Float64
    iq::Float64
end

# PHPS GenRouPHS.init_from_phasor / GenSal.init_from_phasor: rotor angle from the q-axis EMF
# V + (ra + j xq) I, then the dq quantities.
function _dq_at_phasor(ra, xq, V::ComplexF64, I::ComplexF64)
    Eq_phasor = V + complex(ra, xq) * I
    delta = angle(Eq_phasor)
    dq_factor = cis(-(delta - π / 2))
    Vdq = V * dq_factor
    Idq = I * dq_factor
    return delta, real(Vdq), imag(Vdq), real(Idq), imag(Idq)
end

"""
    init_from_phasor(c, V, I) -> (x, targets::MachineTargets)

Machine states at terminal voltage `V` and current `I` (network frame, pu).
"""
function init_from_phasor(c::GENROU_PHTRUE, V::ComplexF64, I::ComplexF64)
    p = c.p
    ra = p.ra
    delta, vd, vq, id, iq = _dq_at_phasor(ra, p.xq, V, I)
    psi_d_pp = vq + ra * iq + p.xd_double_prime * id
    psi_q_pp = -(vd + ra * id - p.xq_double_prime * iq)
    # settled damper: psi_d = Eq' - (xd' - xl) id, so Eq' = psi_d'' + (xd' - xd'') id
    Eq_p = psi_d_pp + (p.xd_prime - p.xd_double_prime) * id
    psi_d = Eq_p - (p.xd_prime - p.xl) * id
    Ed_p = -psi_q_pp - (p.xq_prime - p.xq_double_prime) * iq
    psi_q = -Ed_p - (p.xq_prime - p.xl) * iq
    Efd = Eq_p + (p.xd - p.xd_prime) * id
    Tm = vd * id + vq * iq
    x = [delta, 1.0, Eq_p, psi_d, Ed_p, psi_q]
    return x, MachineTargets(Efd, Tm, abs(V), _genrou_ifd(p, Eq_p, psi_d), vd, vq, id, iq)
end

function init_from_phasor(c::GENSAL_PHTRUE, V::ComplexF64, I::ComplexF64)
    p = c.p
    ra = p.ra
    delta, vd, vq, id, iq = _dq_at_phasor(ra, p.xq, V, I)
    psi_d_pp = vq + ra * iq + p.xd_double_prime * id
    psi_q_pp = -(p.xq - p.xq_double_prime) * iq
    Eq_p = psi_d_pp + (p.xd_prime - p.xd_double_prime) * id
    psi_d = Eq_p - (p.xd_prime - p.xl) * id
    Efd = Eq_p + (p.xd - p.xd_prime) * id          # saturation ignored at init, as in PHPS
    Tm = vd * id + vq * iq
    x = [delta, 1.0, Eq_p, psi_d, psi_q_pp]
    return x, MachineTargets(Efd, Tm, abs(V), _gensal_ifd(p, Eq_p, psi_d), vd, vq, id, iq)
end

# -- grid-forming converters -------------------------------------------------------------
#
# PHPS init_from_phasor of the converters: the internal EMF behind the series (and virtual)
# impedance at the operating point, the measurement states settled at V and I, the
# reservoir at P_STAR C_TANK, and the set-points captured there (p_set, u_set, q_set, ...;
# PowerFactory's p_set_eff / u_set_eff). They are also PHPS's `rebalance_for_bus_voltage`
# at a new bus voltage, with I = (u_out - V) / (j Zseries) from the existing Norton source.

const ConverterComponent = Union{GFM_VSM_PHTRUE,GFM_DROOP_PHTRUE,GFM_VOC_PHTRUE}

"""
    converter_init(c, V, I) -> (x, params::Dict, targets::MachineTargets)

A grid-forming converter's states at terminal voltage `V` and injected current `I`, and the
parameters PHPS's initialisation sets (`Efd0`, `Tm0` and the model's set-points).
"""
function converter_init(c::GFM_VSM_PHTRUE, V::ComplexF64, I::ComplexF64)
    p = c.p
    d = param_dict(c)
    Zs = p.Zseries
    n_vi_init = p.n_vi < 1.0 ? 1.0 : p.n_vi
    # the adaptive boost at the load-flow current, so init is exact in overcurrent too
    ov_i = max(abs(I) - p.i_lim, 0.0)
    r_eff_i = p.r_vi + p.kpr * ov_i
    x_eff_i = p.x_vi + p.kpx * ov_i
    Z_tot_i = complex(r_eff_i, Zs + x_eff_i)
    E = V + n_vi_init * Z_tot_i * I
    theta0 = angle(E)
    u_mag0 = abs(E)
    P_LF = real(V * conj(I))
    params = Dict{String,Float64}("p_set" => P_LF, "u_set" => u_mag0, "PSET_REF" => P_LF,
                                  "Efd0" => u_mag0, "Tm0" => P_LF)
    p.v_set <= 1e-6 && (params["v_set"] = abs(V))
    u_field0 = p.f_set > 1e-6 ? u_mag0 / p.f_set : u_mag0
    x_tank0 = _p(d, "P_STAR", c.name) * p.C_TANK
    x = [theta0, p.f_set, p.f_set, u_mag0, x_tank0, real(V), imag(V), real(I), imag(I), u_field0]
    return x, params, _converter_targets(u_mag0, P_LF, V)
end

function converter_init(c::GFM_DROOP_PHTRUE, V::ComplexF64, I::ComplexF64)
    p = c.p
    d = param_dict(c)
    Zs = p.Zseries
    # the virtual-impedance drop present at the operating point: R always, X only in the
    # paper Zv mode (the legacy X path is a high-pass, zero at steady state)
    paper_zv = p.vi_mode > 0.5
    Rv = 0.0
    if p.adapt_vi > 0.5
        S_min = max(p.S_min, 1.0e-3)
        Sa0 = max(1.0, S_min)
        Rv = max(p.a1 / Sa0 + p.a0, 0.0)
    end
    r_vi_i = p.r_vi + Rv
    x_vi_i = paper_zv ? p.x_vi + Rv : 0.0
    E = V + complex(r_vi_i, Zs + x_vi_i) * I
    theta0 = angle(E)
    u_mag0 = abs(E)
    S = V * conj(I)
    P_LF, Q_LF = real(S), imag(S)
    params = Dict{String,Float64}("p_set" => P_LF, "q_set" => Q_LF, "u_set" => u_mag0,
                                  "Efd0" => u_mag0, "Tm0" => P_LF)
    x_tank0 = _p(d, "P_STAR", c.name) * p.C_TANK
    x = [theta0, p.f_set, Q_LF, u_mag0, x_tank0, real(V), imag(V), real(I), imag(I), real(I),
         imag(I), 1.0, real(I), imag(I)]
    return x, params, _converter_targets(u_mag0, P_LF, V)
end

function converter_init(c::GFM_VOC_PHTRUE, V::ComplexF64, I::ComplexF64)
    p = c.p
    d = param_dict(c)
    v0 = V + complex(p.r_vi, p.Zseries) * I
    S = V * conj(I)
    P_LF, Q_LF = real(S), imag(S)
    params = Dict{String,Float64}("p_set" => P_LF, "q_set" => Q_LF, "PSET_REF" => P_LF,
                                  "Efd0" => abs(v0), "Tm0" => P_LF)
    _p(d, "vnom_from_lf", c.name) > 0.5 && (params["V_nom"] = abs(v0))
    x_tank0 = _p(d, "P_STAR", c.name) * p.C_TANK
    x = [real(v0), imag(v0), x_tank0, real(V), imag(V), real(I), imag(I), real(I), imag(I),
         real(S), imag(S)]
    return x, params, _converter_targets(abs(v0), P_LF, V)
end

_converter_targets(Efd, Tm, V) = MachineTargets(Efd, Tm, abs(V), 0.0, real(V), imag(V), 0.0, 0.0)

"""
    converter_current(c, x, V) -> ComplexF64

The current a converter injects at bus voltage `V`: (u_out - V) / (j Zseries), with u_out
its Norton source (PHPS rebalance_for_bus_voltage).
"""
function converter_current(c::ConverterComponent, x, V::ComplexF64)
    u_out = complex(_converter_uout(c, x)...)
    return _py_cdiv(u_out - V, complex(0.0, c.p.Zseries))
end

_converter_uout(c::GFM_VSM_PHTRUE, x) = _vsm_uout(NoModes(), c.p, x[1], _vsm_umag(c.p, x), x[6], x[7])
_converter_uout(c::GFM_DROOP_PHTRUE, x) = _droop_uout(NoModes(), c.p, x)
_converter_uout(c::GFM_VOC_PHTRUE, x) = _voc_uout(NoModes(), c.p, x)

"""
    init_from_targets(c, t::MachineTargets) -> (x, params::Dict)

Exciter or governor states at its machine's steady state, and the parameters
initialisation sets (set-points and reservoir references).
"""
function init_from_targets(c::IEEET1_PHTRUE, t::MachineTargets)
    p = c.p
    # PHPS Ieeet1PHS.init_from_targets
    Se = p.SAT_A * exp(p.SAT_B * t.Efd)
    Vr = max(p.VRMIN, min(p.VRMAX, (p.KE + Se) * t.Efd))
    Vref = p.KA > 1e-6 ? t.Vt + Vr / p.KA : t.Vt
    # Ieeet1PHTrue: field tank at p = F_STAR, reference = delivered field power
    x = [t.Vt, Vr, t.Efd, t.Efd, p.F_STAR * p.C_FIELD]
    return x, Dict("Vref" => Vref, "PFD_REF" => t.Efd * t.i_fd)
end

function init_from_targets(c::IEEEG1_PHTRUE, t::MachineTargets)
    p = c.p
    # PHPS Ieeeg1PHS.init_from_targets (defaults PMIN -999, PMAX 999 only if unset)
    Tm = max(p.PMIN, min(p.PMAX, t.Tm))
    x = [Tm, Tm, Tm, Tm, Tm, Tm, p.P_STAR * p.C_TANK]
    return x, Dict("Pref" => 1.0 + Tm / p.K, "PM_REF" => Tm)
end

function init_from_targets(c::IEEEG3_PHTRUE, t::MachineTargets)
    p = c.p
    # PHPS Ieeeg3PHS._equilibrium_states: xp = 0, xr = at = x1 = Tm/a23 (gate clamped)
    gate = max(p.PMIN * p.R_base, min(p.PMAX * p.R_base, t.Tm / p.a23))
    x = [0.0, gate, gate, gate, p.P_STAR * p.C_TANK]
    return x, Dict("Pref" => p.Sigma * gate, "PM_REF" => t.Tm)
end

"""
    first_pass(case; pf = solve_powerflow(case)) -> (x, init_params, targets)

PHPS `Initializer.run` for the ported models: machines from their power-flow phasors (power
split among the machines of a bus as PHPS does), then exciters and governors from their
machine's targets, loads at rest with `V0 = Vini =` the power-flow voltage. `x` has the
state layout of [`assemble`](@ref) (`delta_COI = 0`); `init_params` maps component names to
the parameters initialisation sets.
"""
function first_pass(case::Case; pf::PowerFlowResult = solve_powerflow(case))
    net = pf.spec.net
    comps = [build_component(case, s) for s in case.components]
    index = Dict(name(c) => k for (k, c) in enumerate(comps))
    offsets = cumsum([1; [nstates(c) for c in comps]])
    x = zeros(offsets[end])                               # + delta_COI
    init = Dict{String,Dict{String,Float64}}()
    targets = Dict{String,MachineTargets}()

    S = bus_power(Matrix(ybus_pf(case; net)), pf.V, pf.theta)
    # machines, grouped by bus in order of first appearance
    gens_by_bus = Dict{Int,Vector{Int}}()
    order = Int[]
    for (k, c) in enumerate(comps)
        component_role(c) === :generator || continue
        b = bus(c)
        haskey(gens_by_bus, b) || (gens_by_bus[b] = Int[]; push!(order, b))
        push!(gens_by_bus[b], k)
    end
    for b in order
        i = net.index[b]
        Vp = pf.V[i] * cis(pf.theta[i])
        P_load = 0.0
        Q_load = 0.0
        for l in case.pq
            l.bus == b && (P_load += l.p0; Q_load += l.q0)
        end
        S_total = complex(real(S[i]) + P_load, imag(S[i]) + Q_load)
        ks = gens_by_bus[b]
        S_gens = _split_machine_power(comps, ks, S_total)
        for (k, Sg) in zip(ks, S_gens)
            I = conj(_py_cdiv(Sg, Vp))
            if comps[k] isa ConverterComponent
                xk, pk, t = converter_init(comps[k], Vp, I)
                init[name(comps[k])] = pk
            else
                xk, t = init_from_phasor(comps[k], Vp, I)
                init[name(comps[k])] = Dict("Efd0" => t.Efd, "Tm0" => t.Tm)
            end
            x[offsets[k]:offsets[k] + length(xk) - 1] .= xk
            targets[name(comps[k])] = t
        end
    end
    # exciters and governors, from their machine
    for (k, c) in enumerate(comps)
        role = component_role(c)
        role in (:exciter, :governor) || continue
        g = _machine_of(case, comps, index, k)
        (g === nothing || !haskey(targets, name(comps[g]))) && continue
        xk, pk = init_from_targets(c, targets[name(comps[g])])
        x[offsets[k]:offsets[k] + length(xk) - 1] .= xk
        init[name(c)] = pk
    end
    # loads: at rest (z = 0) at their power-flow voltage
    for c in comps
        component_role(c) === :load || continue
        V0 = pf.V[net.index[bus(c)]]
        init[name(c)] = Dict("V0" => V0, "Vini" => V0)
    end
    return x, init, targets
end

# PHPS: one machine takes the bus total (or its p_override / q_override); several machines
# take their p0 / q0 and share the remainder equally.
function _split_machine_power(comps, ks, S_total::ComplexF64)
    if length(ks) == 1
        pd = param_dict(comps[ks[1]])
        P = haskey(pd, "p_override") ? param_value(pd["p_override"]) : real(S_total)
        Q = haskey(pd, "q_override") ? param_value(pd["q_override"]) : imag(S_total)
        return (haskey(pd, "p_override") || haskey(pd, "q_override")) ? [complex(P, Q)] : [S_total]
    end
    S = fill(0.0im, length(ks))
    remaining = S_total
    free = Int[]
    for (j, k) in enumerate(ks)
        pd = param_dict(comps[k])
        if haskey(pd, "p0")
            q = haskey(pd, "q0") ? param_value(pd["q0"]) : 0.0
            S[j] = complex(param_value(pd["p0"]), q)
            remaining -= S[j]
        else
            push!(free, j)
        end
    end
    if isempty(free)
        S[1] += remaining
    else
        share = remaining / length(free)
        for j in free
            S[j] = share
        end
    end
    return S
end

# The machine an exciter or governor belongs to: the machine whose Efd is wired from the
# exciter / whose omega feeds the governor; else the `syn` parameter (PHPS
# Initializer._get_generator_for_comp).
function _machine_of(case::Case, comps, index, k::Int)
    c = comps[k]
    nm = name(c)
    if component_role(c) === :exciter
        for w in case.connections
            src = split(w.from, '.'; limit = 2)
            (length(src) == 2 && src[1] == nm) || continue
            dst = split(w.to, '.'; limit = 2)
            g = get(index, dst[1], 0)
            g > 0 && component_role(comps[g]) === :generator && return g
        end
    else
        for w in case.connections
            startswith(w.to, nm * ".omega") || continue
            g = get(index, split(w.from, '.')[1], 0)
            g > 0 && component_role(comps[g]) === :generator && return g
        end
    end
    syn = get(param_dict(c), "syn", nothing)
    syn === nothing && return nothing
    return get(index, string(syn), nothing)
end
