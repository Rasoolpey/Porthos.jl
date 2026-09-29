# COMPLEXLOAD (contract key ComplexLoad): PowerFactory General Load (ElmLod / TypLod) with a
# t1-filtered voltage dependence and the out-of-range multiplier k(V).
#
# Port of PHPS ComplexLoadPHS (phps/src/components/loads/complex_load_phs.py at ba11ea1).
# The Y-bus already carries the constant-impedance base (P0 - jQ0)/V0^2, so the model
# injects only the difference between the PowerFactory law and that base.
#
# PHPS writes its constants into the C++ kernel with printf("%.12e") (exponents "%.6f"), so
# the kernel computes with those rounded values; this port rounds the same way.
#
# States  [z] when t1 > 0 (filtered voltage deviation), none otherwise
# Inputs  [Vd, Vq]
# Outputs [Id, Iq, Pload, Qload]
#
# Branch sites (output kernel): the voltage bands |V| < udmin/2, |V| <= udmin,
# |V| >= udmax (evaluated as an if/elseif chain), the guard of the voltage ratio, and the
# |V|^2 > 1e-8 clip of the correction current.

struct ComplexLoadParams
    t1::Float64          # > 0: filtered (one state)
    has_band::Bool       # udmin > 0: out-of-range multiplier active
    t1_c::Float64        # constants as written into PHPS's C++
    Vini_c::Float64
    udmin_c::Float64
    udmax_c::Float64
    P0_c::Float64
    Q0_c::Float64
    kpv_c::Float64
    kqv_c::Float64
    G_load_c::Float64
    Q_base_c::Float64
end

struct COMPLEXLOAD <: AbstractComponent
    name::String
    bus::Int
    p::ComplexLoadParams
    params::ParamDict
end

const _COMPLEXLOAD_INPUTS = ["Vd", "Vq"]
const _COMPLEXLOAD_OUTPUTS = ["Id", "Iq", "Pload", "Qload"]

model_type(::COMPLEXLOAD) = "COMPLEXLOAD"
component_role(::COMPLEXLOAD) = :load
state_names(c::COMPLEXLOAD) = c.p.t1 > 0.0 ? ["z"] : String[]
input_names(::COMPLEXLOAD) = _COMPLEXLOAD_INPUTS
output_names(::COMPLEXLOAD) = _COMPLEXLOAD_OUTPUTS
bus(c::COMPLEXLOAD) = c.bus

const _COMPLEXLOAD_DEFAULTS = ("P0" => 0.0, "Q0" => 0.0, "V0" => 1.0, "kpf" => 0.0,
                               "kqf" => 0.0, "kpv" => 2.0, "kqv" => 2.0, "t1" => 0.0,
                               "udmin" => 0.0, "udmax" => 0.0)

function type_defaults!(::Val{:COMPLEXLOAD}, p::ParamDict)
    for (k, v) in _COMPLEXLOAD_DEFAULTS
        get!(p, k, v)
    end
    # Vini: voltage at RMS initialisation (set by initialisation; defaults to V0)
    get!(p, "Vini", p["V0"])
    return p
end

function COMPLEXLOAD(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    V0 = f("V0")
    V02 = max(V0 * V0, 1e-12)
    P0, Q0 = f("P0"), f("Q0")
    udmin = f("udmin")
    p = ComplexLoadParams(f("t1"), udmin > 0.0, _c12e(f("t1")), _c12e(f("Vini")), _c12e(udmin),
                          _c12e(f("udmax")), _c12e(P0), _c12e(Q0), _c6f(f("kpv")),
                          _c6f(f("kqv")), _c12e(P0 / V02), _c12e(Q0 / V02))
    return COMPLEXLOAD(name, Int(f("bus")), p, d)
end
COMPONENT_CONSTRUCTORS["COMPLEXLOAD"] = COMPLEXLOAD

function _step!(dx, y, c::COMPLEXLOAD, x, u, p::ComplexLoadParams, rec::ModeRecorder)
    if p.t1 > 0.0
        Vd, Vq = u[1], u[2]
        Vmag_s = sqrt(Vd * Vd + Vq * Vq)
        dx[1] = ((Vmag_s - p.Vini_c) - x[1]) / p.t1_c
    end
    return dx
end

# The PowerFactory law and the correction current: (Id, Iq, P, Q).
@inline function _load_law(rec::ModeRecorder, p::ComplexLoadParams, x, u)
    Vd, Vq = u[1], u[2]
    T = promote_type(eltype(x), eltype(u))
    Vmag2 = Vd * Vd + Vq * Vq
    Vmag = sqrt(Vmag2)
    Vini_val = p.Vini_c
    # PowerFactory out-of-range multiplier k(V), C1-continuous
    kV = if !p.has_band
        one(T)
    elseif lt(rec, Vmag, 0.5 * p.udmin_c)
        T(2.0 * Vmag2 / (p.udmin_c * p.udmin_c))
    elseif le(rec, Vmag, p.udmin_c)
        _d = (Vmag - p.udmin_c) / p.udmin_c
        T(1.0 - 2.0 * _d * _d)
    elseif decide(rec, p.udmax_c > p.udmin_c && _ge(Vmag, p.udmax_c))
        _d = Vmag - p.udmax_c
        T(1.0 + _d * _d)
    else
        one(T)
    end
    # filtered voltage ratio 1 + z/Vini (z = |V| - Vini without the filter)
    z = p.t1 > 0.0 ? T(x[1]) : T(Vmag - Vini_val)
    vr = guard_min(rec, 1.0 + z / Vini_val, 1.0e-6)
    P_act = p.P0_c * kV * vr^p.kpv_c
    Q_act = p.Q0_c * kV * vr^p.kqv_c
    # correction relative to the constant-Z base already in the Y-bus
    Vm2 = select(rec, _gt(Vmag2, 1.0e-8), T(Vmag2), T(1.0e-8))
    Gp = P_act / Vm2 - p.G_load_c
    Bq = Q_act / Vm2 - p.Q_base_c
    return -(Gp * Vd + Bq * Vq), -(Gp * Vq - Bq * Vd), P_act, Q_act
end

function _outputs!(y, c::COMPLEXLOAD, x, u, p::ComplexLoadParams, rec::ModeRecorder)
    Id, Iq, P_act, Q_act = _load_law(rec, p, x, u)
    y[3] = P_act
    y[4] = Q_act
    y[1] = Id
    y[2] = Iq
    return y
end

hamiltonian(c::COMPLEXLOAD, x, p::ComplexLoadParams) = 0.0      # dissipative port
grad_hamiltonian!(g, c::COMPLEXLOAD, x, p::ComplexLoadParams) = fill!(g, 0.0)

"""Correction current into the bus for bus voltage `V = (Vd, Vq)` (depends on V)."""
function injection(c::COMPLEXLOAD, x, V, p::ComplexLoadParams = c.p)
    Id, Iq, _, _ = _load_law(NoModes(), p, x, V)
    return Id, Iq
end
