# GENSAL_PHTRUE: fifth-order salient-pole synchronous machine (no q-axis transient winding)
# with the Efd-collocated field flow i_fd exposed.
#
# Port of the equations PHPS runs for GENSAL_PHTRUE (phps/src/components/generators/
# gensal_phtrue.py at ba11ea1): the kernels of its retired parent GenSal with the PHTRUE
# "V4" substitutions applied (subtransient EMF omega*psi'', swing Tm/omega), plus
# outputs[9] = i_fd. Expressions keep PHPS's order of operations.
#
# States  [delta, omega, E_q_prime, psi_d, psi_q_pp]
# Inputs  [Vd, Vq, Tm, Efd]
# Outputs [Id, Iq, omega, Pe, Qe, id_dq, iq_dq, It_Re, It_Im, i_fd]
#
# Branch sites (step): the saturation test `psi_ag > A_sat && psi_ag > 1e-6`.

struct GensalParams
    H::Float64
    D::Float64
    ra::Float64
    xd::Float64
    xq::Float64
    xd_prime::Float64
    xd_double_prime::Float64
    xq_double_prime::Float64
    Td0_prime::Float64
    Td0_double_prime::Float64
    Tq0_double_prime::Float64
    xl::Float64
    A_sat::Float64
    B_sat::Float64
    omega_b::Float64
end

struct GENSAL_PHTRUE <: AbstractComponent
    name::String
    bus::Int
    p::GensalParams
    params::ParamDict
end

const _GENSAL_STATES = ["delta", "omega", "E_q_prime", "psi_d", "psi_q_pp"]

model_type(::GENSAL_PHTRUE) = "GENSAL_PHTRUE"
state_names(::GENSAL_PHTRUE) = _GENSAL_STATES
input_names(::GENSAL_PHTRUE) = _MACHINE_INPUTS
output_names(::GENSAL_PHTRUE) = _MACHINE_OUTPUTS
bus(c::GENSAL_PHTRUE) = c.bus

"""
    gensal_sat_coeffs(S10, S12) -> (A, B)

Scaled-quadratic saturation `Sat(x) = B (x - A)^2 / x` through `(1.0, S10)` and
`(1.2, S12)` (PHPS GenSal._sat_coeffs).
"""
function gensal_sat_coeffs(S10::Float64, S12::Float64)
    (S10 <= 0.0 || S12 <= 0.0) && return 0.0, 0.0
    ratio = sqrt(1.2 * S12 / S10)
    abs(ratio - 1.0) < 1e-10 && return 0.0, 0.0
    u = 0.2 / (ratio - 1.0)
    return 1.0 - u, S10 / (u * u)
end

function type_defaults!(::Val{:GENSAL_PHTRUE}, p::ParamDict)
    get!(p, "S10", 0.0)
    get!(p, "S12", 0.0)
    p["A_sat"], p["B_sat"] = gensal_sat_coeffs(param_value(p["S10"]), param_value(p["S12"]))
    return p
end

function GENSAL_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    abs(f("xd_double_prime") - f("xq_double_prime")) > 1e-9 &&
        throw(ArgumentError("$name: GENSAL_PHTRUE requires xd_double_prime == " *
                            "xq_double_prime; subtransient saliency is not modelled"))
    p = GensalParams(f("H"), f("D"), f("ra"), f("xd"), f("xq"), f("xd_prime"),
                     f("xd_double_prime"), f("xq_double_prime"), f("Td0_prime"),
                     f("Td0_double_prime"), f("Tq0_double_prime"), f("xl"), f("A_sat"),
                     f("B_sat"), f("omega_b"))
    return GENSAL_PHTRUE(name, Int(f("bus")), p, d)
end
COMPONENT_CONSTRUCTORS["GENSAL_PHTRUE"] = GENSAL_PHTRUE

@inline function _gensal_norton(p::GensalParams, x)
    delta = x[1]
    sin_d = sin(delta)
    cos_d = cos(delta)
    Eq_p, psi_d, psi_q_pp = x[3], x[4], x[5]
    k_d = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    psi_d_pp = Eq_p * k_d + psi_d * (1.0 - k_d)
    det = p.ra * p.ra + p.xd_double_prime * p.xq_double_prime
    id_no = x[2] * (-p.ra * psi_q_pp + p.xq_double_prime * psi_d_pp) / det
    iq_no = x[2] * (p.xd_double_prime * psi_q_pp + p.ra * psi_d_pp) / det
    I_Re = id_no * sin_d + iq_no * cos_d
    I_Im = -id_no * cos_d + iq_no * sin_d
    return I_Re, I_Im, id_no, iq_no
end

@inline _gensal_ifd(p, Eq_p, psi_d) =
    (Eq_p / (p.xd - p.xd_prime)
     + (p.xd_prime - p.xd_double_prime) * (Eq_p - psi_d)
       / ((p.xd_prime - p.xl) * (p.xd_prime - p.xl))) / (p.omega_b * p.Td0_prime)

function _outputs!(y, c::GENSAL_PHTRUE, x, u, p::GensalParams, rec::ModeRecorder)
    I_Re, I_Im, id_no, iq_no = _gensal_norton(p, x)
    y[1] = I_Re
    y[2] = I_Im
    y[3] = x[2]
    y[4] = 0.0
    y[5] = 0.0
    y[6] = id_no
    y[7] = iq_no
    y[8] = I_Re
    y[9] = I_Im
    y[10] = _gensal_ifd(p, x[3], x[4])
    return y
end

function _step!(dx, y, c::GENSAL_PHTRUE, x, u, p::GensalParams, rec::ModeRecorder)
    V_Re, V_Im, Tm, Efd = u[1], u[2], u[3], u[4]
    delta, omega, Eq_p, psi_d, psi_q_pp = x[1], x[2], x[3], x[4], x[5]
    sin_d = sin(delta)
    cos_d = cos(delta)
    vd = V_Re * sin_d - V_Im * cos_d
    vq = V_Re * cos_d + V_Im * sin_d
    k_d = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    psi_d_pp = Eq_p * k_d + psi_d * (1.0 - k_d)
    rhs_d = vd + omega * psi_q_pp
    rhs_q = vq - omega * psi_d_pp
    det = p.ra * p.ra + p.xd_double_prime * p.xq_double_prime
    id = (-p.ra * rhs_d - p.xq_double_prime * rhs_q) / det
    iq = (p.xd_double_prime * rhs_d - p.ra * rhs_q) / det
    Te = psi_d_pp * iq - psi_q_pp * id

    # additive saturation on the air-gap flux (scaled quadratic)
    psi_d_term = vq + p.ra * iq
    psi_q_term = vd + p.ra * id
    psi_ag = sqrt(psi_d_term * psi_d_term + psi_q_term * psi_q_term)
    sat_val = if decide(rec, _gt(psi_ag, p.A_sat) && _gt(psi_ag, 1e-6))
        diff = psi_ag - p.A_sat
        p.B_sat * diff * diff / psi_ag
    else
        zero(psi_ag)
    end
    Efd_net = Efd - sat_val * (p.xd - p.xl)

    dx[1] = p.omega_b * (omega - 1.0)
    dx[2] = (Tm / omega - Te - p.D * (omega - 1.0)) / (2.0 * p.H)
    dx[3] = (Efd_net - Eq_p - (p.xd - p.xd_prime) * id) / p.Td0_prime
    dx[4] = (-psi_d - (p.xd_prime - p.xl) * id + Eq_p) / p.Td0_double_prime
    dx[5] = (-(psi_q_pp) - (p.xq - p.xq_double_prime) * iq) / p.Tq0_double_prime

    if y !== nothing
        y[4] = vd * id + vq * iq
        y[5] = vq * id - vd * iq
        y[6] = id
        y[7] = iq
        y[8] = id * sin_d + iq * cos_d
        y[9] = -id * cos_d + iq * sin_d
        y[10] = _gensal_ifd(p, x[3], x[4])
    end
    return dx
end

"""Kinetic plus reduced-order rotor magnetic storage (PHPS GenSal.hamiltonian)."""
function hamiltonian(c::GENSAL_PHTRUE, x, p::GensalParams)
    omega, Eq_p, psi_d, psi_q_pp = x[2], x[3], x[4], x[5]
    wb = p.omega_b
    kd = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    Eq_pp = kd * Eq_p + (1.0 - kd) * psi_d
    h_mech = p.H * (omega - 1.0)^2
    h_d = (Eq_p^2 / (p.xd - p.xd_prime)
           + (Eq_p - Eq_pp)^2 / (p.xd_prime - p.xd_double_prime)) / (2.0 * wb)
    h_q = (psi_q_pp^2 / (p.xq - p.xq_double_prime)) / (2.0 * wb)
    return h_mech + h_d + h_q
end

function grad_hamiltonian!(g, c::GENSAL_PHTRUE, x, p::GensalParams)
    omega, Eq_p, psi_d, psi_q_pp = x[2], x[3], x[4], x[5]
    wb = p.omega_b
    qd_mut = (p.xd_prime - p.xd_double_prime) / (p.xd_prime - p.xl)^2
    g[1] = 0.0
    g[2] = 2.0 * p.H * (omega - 1.0)
    g[3] = (Eq_p / (p.xd - p.xd_prime) + qd_mut * (Eq_p - psi_d)) / wb
    g[4] = -qd_mut * (Eq_p - psi_d) / wb
    g[5] = psi_q_pp / (wb * (p.xq - p.xq_double_prime))
    return g
end

"""Norton current (network frame), independent of the bus voltage."""
function injection(c::GENSAL_PHTRUE, x, V, p::GensalParams = c.p)
    I_Re, I_Im, _, _ = _gensal_norton(p, x)
    return I_Re, I_Im
end

norton_admittance(c::GENSAL_PHTRUE) = _norton_y(c.p.ra, c.p.xd_double_prime)
