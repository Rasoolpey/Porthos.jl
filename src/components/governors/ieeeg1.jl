# IEEEG1_PHTRUE: IEEE steam-turbine governor with a steam-reservoir account.
#
# Port of the equations PHPS runs for IEEEG1_PHTRUE (phps/src/components/governors/
# ieeeg1_phtrue.py at ba11ea1): the six states of its retired parent Ieeeg1PHS, unchanged,
# plus the one-way steam tank drained by the turbine power (corrected constant-power law).
# Expressions keep PHPS's order of operations.
#
# States  [x0, x1, x2, x3, x4, x5, x_steam]
# Inputs  [omega, Pref, u_agc]
# Outputs [Tm]
#
# Branch sites (step): valve position clamp, servo rate limits UO/UC, non-windup at
# PMAX/PMIN, out-of-band relaxation of the raw valve state, the tank pressure guard.
# The droop speed path is scaled to the system base by R_base = Pnom/Sbase.

struct Ieeeg1Params
    K::Float64
    T1::Float64
    T2::Float64
    T3::Float64
    T4::Float64
    T5::Float64
    T6::Float64
    T7::Float64
    K1::Float64
    K3::Float64
    K5::Float64
    K7::Float64
    PMAX::Float64
    PMIN::Float64
    UO::Float64
    UC::Float64
    R_base::Float64
    C_TANK::Float64
    P_STAR::Float64
    PM_REF::Float64
end

struct IEEEG1_PHTRUE <: AbstractComponent
    name::String
    p::Ieeeg1Params
    params::ParamDict
end

const _IEEEG1_STATES = ["x0", "x1", "x2", "x3", "x4", "x5", "x_steam"]
const _GOVERNOR_INPUTS = ["omega", "Pref", "u_agc"]
const _GOVERNOR_OUTPUTS = ["Tm"]

model_type(::IEEEG1_PHTRUE) = "IEEEG1_PHTRUE"
component_role(::IEEEG1_PHTRUE) = :governor
state_names(::IEEEG1_PHTRUE) = _IEEEG1_STATES
input_names(::IEEEG1_PHTRUE) = _GOVERNOR_INPUTS
output_names(::IEEEG1_PHTRUE) = _GOVERNOR_OUTPUTS

function type_defaults!(::Val{:IEEEG1_PHTRUE}, p::ParamDict)
    get!(p, "C_TANK", 1000.0)
    get!(p, "P_STAR", 1.0)
    get!(p, "PM_REF", 0.0)
    return p
end

function IEEEG1_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    p = Ieeeg1Params(f("K"), f("T1"), f("T2"), f("T3"), f("T4"), f("T5"), f("T6"), f("T7"),
                     f("K1"), f("K3"), f("K5"), f("K7"), f("PMAX"), f("PMIN"), f("UO"),
                     f("UC"), _p(d, "R_base", name, 1.0), f("C_TANK"), f("P_STAR"),
                     f("PM_REF"))
    return IEEEG1_PHTRUE(name, p, d)
end
COMPONENT_CONSTRUCTORS["IEEEG1_PHTRUE"] = IEEEG1_PHTRUE

@inline _ieeeg1_tm(p::Ieeeg1Params, x) = p.K1 * x[3] + p.K3 * x[4] + p.K5 * x[5] + p.K7 * x[6]

function _outputs!(y, c::IEEEG1_PHTRUE, x, u, p::Ieeeg1Params, rec::ModeRecorder)
    y[1] = _ieeeg1_tm(p, x)
    return y
end

function _step!(dx, y, c::IEEEG1_PHTRUE, x, u, p::Ieeeg1Params, rec::ModeRecorder)
    omega, Pref, u_agc = u[1], u[2], u[3]
    err = Pref + u_agc - 1.0 - p.R_base * (omega - 1.0)

    # lead-lag governor K (1 + s T2) / (1 + s T1)
    dx[1] = (p.K * err - x[1]) / p.T1
    GV = x[1] * (1.0 - p.T2 / p.T1) + p.K * (p.T2 / p.T1) * err

    # servo: position-limited, rate-limited, leak-proof
    x1_lim = clamp_mode(rec, x[2], p.PMIN, p.PMAX)
    srv_rate = (GV - x1_lim) / p.T3
    srv_rate = clamp_mode(rec, srv_rate, p.UC, p.UO)
    srv_rate = nonwindup(rec, x1_lim, srv_rate, p.PMIN, p.PMAX)
    srv_rate = outband_relax(rec, srv_rate, x[2], p.PMIN, p.PMAX, p.T3)
    dx[2] = srv_rate

    # cascaded turbine stages, fed by the limited valve
    dx[3] = (x1_lim - x[3]) / p.T4
    dx[4] = (x[3] - x[4]) / p.T5
    dx[5] = (x[4] - x[5]) / p.T6
    dx[6] = (x[5] - x[6]) / p.T7

    # steam reservoir (one-way; does not feed back)
    Tm_out = _ieeeg1_tm(p, x)
    p_tank = guard_min(rec, x[7] / p.C_TANK, 1e-6)
    dx[7] = (p.PM_REF - Tm_out) / p_tank
    return dx
end

hamiltonian(c::IEEEG1_PHTRUE, x, p::Ieeeg1Params) = 0.5 * (x[7] * x[7]) / p.C_TANK

function grad_hamiltonian!(g, c::IEEEG1_PHTRUE, x, p::Ieeeg1Params)
    fill!(g, 0.0)
    g[7] = x[7] / p.C_TANK
    return g
end
