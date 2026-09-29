# IEEEG3_PHTRUE: IEEE hydro governor (PowerFactory pcu_IEEEG3) with a headrace-reservoir
# account.
#
# Port of the equations PHPS runs for IEEEG3_PHTRUE (phps/src/components/governors/
# ieeeg3_phtrue.py at ba11ea1): the four states of its retired parent Ieeeg3PHS, unchanged,
# plus the one-way headrace drained by the turbine power (corrected constant-power law).
# Expressions keep PHPS's order of operations.
#
# States  [xp, xr, at, x1, x_water]
# Inputs  [omega, Pref, u_agc]
# Outputs [Tm]
#
# Branch sites: gate position clamp (also in the output kernel), pilot-valve rate clamp,
# non-windup and out-of-band relaxation of the pilot valve and of the gate, the headrace
# guard. Gate position and rate limits are carried to the system base by R_base.

struct Ieeeg3Params
    Tg::Float64
    Tp::Float64
    Sigma::Float64
    Delta::Float64
    Tr::Float64
    a23::Float64
    UO::Float64
    UC::Float64
    PMAX::Float64
    PMIN::Float64
    R_base::Float64
    Ta_t::Float64        # a11 Tw, water-column lag
    Tb_t::Float64        # (a11 - a13 a21 / a23) Tw, lead (negative)
    C_TANK::Float64
    P_STAR::Float64
    PM_REF::Float64
end

struct IEEEG3_PHTRUE <: AbstractComponent
    name::String
    p::Ieeeg3Params
    params::ParamDict
end

const _IEEEG3_STATES = ["xp", "xr", "at", "x1", "x_water"]

model_type(::IEEEG3_PHTRUE) = "IEEEG3_PHTRUE"
component_role(::IEEEG3_PHTRUE) = :governor
state_names(::IEEEG3_PHTRUE) = _IEEEG3_STATES
input_names(::IEEEG3_PHTRUE) = _GOVERNOR_INPUTS
output_names(::IEEEG3_PHTRUE) = _GOVERNOR_OUTPUTS

const _IEEEG3_DEFAULTS = ("Tg" => 0.05, "Tp" => 0.04, "Sigma" => 0.04, "Delta" => 0.2,
                          "Tr" => 10.0, "a11" => 0.5, "a13" => 1.0, "a21" => 1.5,
                          "a23" => 1.0, "Tw" => 0.75, "UO" => 0.1, "UC" => -0.1,
                          "PMAX" => 1.0, "PMIN" => 0.0, "R_base" => 1.0)

function type_defaults!(::Val{:IEEEG3_PHTRUE}, p::ParamDict)
    for (k, v) in _IEEEG3_DEFAULTS
        get!(p, k, v)
    end
    get!(p, "C_TANK", 1000.0)
    get!(p, "P_STAR", 1.0)
    get!(p, "PM_REF", 0.0)
    return p
end

function IEEEG3_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    a11, a13, a21, a23, Tw = f("a11"), f("a13"), f("a21"), f("a23"), f("Tw")
    abs(a23) < 1e-12 && throw(ArgumentError("$name: a23 must be non-zero (turbine gain)"))
    Ta = a11 * Tw
    Ta <= 0.0 && throw(ArgumentError("$name: a11*Tw must be > 0 (got a11=$a11, Tw=$Tw)"))
    Tb = (a11 - a13 * a21 / a23) * Tw
    p = Ieeeg3Params(f("Tg"), f("Tp"), f("Sigma"), f("Delta"), f("Tr"), a23, f("UO"), f("UC"),
                     f("PMAX"), f("PMIN"), f("R_base"), Ta, Tb, f("C_TANK"), f("P_STAR"),
                     f("PM_REF"))
    return IEEEG3_PHTRUE(name, p, d)
end
COMPONENT_CONSTRUCTORS["IEEEG3_PHTRUE"] = IEEEG3_PHTRUE

# Clamped gate and the one turbine output built from it (PHPS _cpp_clamped_gate_and_tm).
@inline function _ieeeg3_gate_tm(rec::ModeRecorder, p::Ieeeg3Params, x)
    GMAX = p.PMAX * p.R_base
    GMIN = p.PMIN * p.R_base
    at_lim = clamp_mode(rec, x[3], GMIN, GMAX)
    dx1_t = (at_lim - x[4]) / p.Ta_t
    Tm_out = p.a23 * (x[4] + p.Tb_t * dx1_t)
    return at_lim, dx1_t, Tm_out, GMIN, GMAX
end

function _outputs!(y, c::IEEEG3_PHTRUE, x, u, p::Ieeeg3Params, rec::ModeRecorder)
    _, _, Tm_out, _, _ = _ieeeg3_gate_tm(rec, p, x)
    y[1] = Tm_out
    return y
end

function _step!(dx, y, c::IEEEG3_PHTRUE, x, u, p::Ieeeg3Params, rec::ModeRecorder)
    omega, Pref, u_agc = u[1], u[2], u[3]
    at_lim, dx1_t, Tm_out, GMIN, GMAX = _ieeeg3_gate_tm(rec, p, x)
    RUP = p.UO * p.R_base
    RDN = p.UC * p.R_base

    # governor summing junction
    droop_perm = p.Sigma * at_lim
    droop_tran = p.Delta * (at_lim - x[2])
    e = Pref + u_agc - p.R_base * (omega - 1.0) - droop_perm - droop_tran

    # pilot valve (Tp): its state is the commanded gate rate; non-windup limited state
    xp_lim = clamp_mode(rec, x[1], RDN, RUP)
    dxp = (e / p.Tg - xp_lim) / p.Tp
    dxp = nonwindup(rec, xp_lim, dxp, RDN, RUP)
    dxp = outband_relax(rec, dxp, x[1], RDN, RUP, p.Tp)
    dx[1] = dxp

    # dashpot / transient-droop reset
    dx[2] = (at_lim - x[2]) / p.Tr

    # gate servo integrator, position-limited
    gate_rate = nonwindup(rec, at_lim, xp_lim, GMIN, GMAX)
    gate_rate = outband_relax(rec, gate_rate, x[3], GMIN, GMAX, p.Tg)
    dx[3] = gate_rate

    # water column (waterhammer) lag
    dx[4] = dx1_t

    # headrace reservoir (one-way; does not feed back)
    h_head = guard_min(rec, x[5] / p.C_TANK, 1e-6)
    dx[5] = (p.PM_REF - Tm_out) / h_head
    return dx
end

hamiltonian(c::IEEEG3_PHTRUE, x, p::Ieeeg3Params) = 0.5 * (x[5] * x[5]) / p.C_TANK

function grad_hamiltonian!(g, c::IEEEG3_PHTRUE, x, p::Ieeeg3Params)
    fill!(g, 0.0)
    g[5] = x[5] / p.C_TANK
    return g
end
