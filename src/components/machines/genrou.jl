# GENROU_PHTRUE: sixth-order round-rotor synchronous machine with the Efd-collocated field
# flow i_fd exposed for the exciter's field-supply account.
#
# Port of the equations PHPS runs for GENROU_PHTRUE (phps/src/components/generators/
# genrou_phtrue.py at ba11ea1). Its Python class inherits the dynamics from the retired
# GenRouPHS; this file is the flattened model: the parent's step and output kernels plus
# the PHTRUE additions (Pe = terminal power, outputs[9] = i_fd). Expressions keep PHPS's
# order of operations.
#
# States  [delta, omega, E_q_prime, psi_d, E_d_prime, psi_q]
# Inputs  [Vd, Vq, Tm, Efd]              (bus voltage in the network RI frame)
# Outputs [Id, Iq, omega, Pe, Qe, id_dq, iq_dq, It_Re, It_Im, i_fd]
#
# Speed convention: subtransient EMF omega*psi'' (PowerFactory partially neglected speed
# variation), swing (Tm/omega - Te - D (omega-1)) / 2H with Te = P/omega.

struct GenrouParams
    H::Float64
    D::Float64
    ra::Float64
    xd::Float64
    xq::Float64
    xd_prime::Float64
    xq_prime::Float64
    xd_double_prime::Float64
    xq_double_prime::Float64
    Td0_prime::Float64
    Tq0_prime::Float64
    Td0_double_prime::Float64
    Tq0_double_prime::Float64
    xl::Float64
    omega_b::Float64
end

struct GENROU_PHTRUE <: AbstractComponent
    name::String
    bus::Int
    p::GenrouParams
    params::ParamDict
end

const _GENROU_STATES = ["delta", "omega", "E_q_prime", "psi_d", "E_d_prime", "psi_q"]
const _MACHINE_INPUTS = ["Vd", "Vq", "Tm", "Efd"]
const _MACHINE_OUTPUTS = ["Id", "Iq", "omega", "Pe", "Qe", "id_dq", "iq_dq", "It_Re", "It_Im",
                          "i_fd"]

model_type(::GENROU_PHTRUE) = "GENROU_PHTRUE"
component_role(::GENROU_PHTRUE) = :generator
state_names(::GENROU_PHTRUE) = _GENROU_STATES
input_names(::GENROU_PHTRUE) = _MACHINE_INPUTS
output_names(::GENROU_PHTRUE) = _MACHINE_OUTPUTS
bus(c::GENROU_PHTRUE) = c.bus

function GENROU_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    for k in ("Td0_prime", "Tq0_prime", "Td0_double_prime", "Tq0_double_prime")
        f(k) > 0.0 || throw(ArgumentError("$name: GENROU_PHTRUE requires $k > 0; use a " *
                                          "reduced-order machine model when a rotor circuit " *
                                          "time constant is zero"))
    end
    p = GenrouParams(f("H"), f("D"), f("ra"), f("xd"), f("xq"), f("xd_prime"), f("xq_prime"),
                     f("xd_double_prime"), f("xq_double_prime"), f("Td0_prime"),
                     f("Tq0_prime"), f("Td0_double_prime"), f("Tq0_double_prime"), f("xl"),
                     f("omega_b"))
    return GENROU_PHTRUE(name, Int(f("bus")), p, d)
end
COMPONENT_CONSTRUCTORS["GENROU_PHTRUE"] = GENROU_PHTRUE

# i_fd = [Eq'/(xd-xd') + (xd'-xd'')(Eq'-psi_d)/(xd'-xl)^2] / (omega_b Td0')
@inline _genrou_ifd(p, Eq_p, psi_d) =
    (Eq_p / (p.xd - p.xd_prime)
     + (p.xd_prime - p.xd_double_prime) * (Eq_p - psi_d)
       / ((p.xd_prime - p.xl) * (p.xd_prime - p.xl))) / (p.omega_b * p.Td0_prime)

# Norton current of the subtransient EMF omega*psi'' behind ra + j x'' (network frame).
@inline function _genrou_norton(p::GenrouParams, x)
    delta = x[1]
    sin_d = sin(delta)
    cos_d = cos(delta)
    Eq_p, psi_d, Ed_p, psi_q = x[3], x[4], x[5], x[6]
    k_d = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    k_q = (p.xq_double_prime - p.xl) / (p.xq_prime - p.xl)
    psi_d_pp = Eq_p * k_d + psi_d * (1.0 - k_d)
    psi_q_pp = -Ed_p * k_q + psi_q * (1.0 - k_q)
    det = p.ra * p.ra + p.xd_double_prime * p.xq_double_prime
    omega = x[2]
    id_no = omega * (-p.ra * psi_q_pp + p.xq_double_prime * psi_d_pp) / det
    iq_no = omega * (p.xd_double_prime * psi_q_pp + p.ra * psi_d_pp) / det
    I_Re = id_no * sin_d + iq_no * cos_d
    I_Im = -id_no * cos_d + iq_no * sin_d
    return I_Re, I_Im, id_no, iq_no
end

function _outputs!(y, c::GENROU_PHTRUE, x, u, p::GenrouParams, rec::ModeRecorder)
    I_Re, I_Im, id_no, iq_no = _genrou_norton(p, x)
    y[1] = I_Re
    y[2] = I_Im
    y[3] = x[2]
    y[4] = 0.0
    y[5] = 0.0
    y[6] = id_no
    y[7] = iq_no
    y[8] = I_Re
    y[9] = I_Im
    y[10] = _genrou_ifd(p, x[3], x[4])
    return y
end

function _step!(dx, y, c::GENROU_PHTRUE, x, u, p::GenrouParams, rec::ModeRecorder)
    V_Re, V_Im, Tm, Efd = u[1], u[2], u[3], u[4]
    delta, omega, Eq_p, psi_d, Ed_p, psi_q = x[1], x[2], x[3], x[4], x[5], x[6]
    sin_d = sin(delta)
    cos_d = cos(delta)
    # Park transform, network -> dq
    vd = V_Re * sin_d - V_Im * cos_d
    vq = V_Re * cos_d + V_Im * sin_d
    k_d = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    k_q = (p.xq_double_prime - p.xl) / (p.xq_prime - p.xl)
    psi_d_pp = Eq_p * k_d + psi_d * (1.0 - k_d)
    psi_q_pp = -Ed_p * k_q + psi_q * (1.0 - k_q)
    # stator algebraic equations, solved for (id, iq)
    rhs_d = vd + omega * psi_q_pp
    rhs_q = vq - omega * psi_d_pp
    det_s = p.ra * p.ra + p.xd_double_prime * p.xq_double_prime
    id = (-p.ra * rhs_d - p.xq_double_prime * rhs_q) / det_s
    iq = (p.xd_double_prime * rhs_d - p.ra * rhs_q) / det_s
    Te = (vd * id + vq * iq) / omega

    inv_2H = 1.0 / (2.0 * p.H)
    dx[1] = p.omega_b * (omega - 1.0)
    dx[2] = inv_2H * (Tm / omega - Te - p.D * (omega - 1.0))
    dx[3] = (Efd - Eq_p - (p.xd - p.xd_prime) * id) / p.Td0_prime
    dx[4] = (Eq_p - psi_d - (p.xd_prime - p.xl) * id) / p.Td0_double_prime
    dx[5] = (-Ed_p + (p.xq - p.xq_prime) * iq) / p.Tq0_prime
    dx[6] = (-Ed_p - psi_q - (p.xq_prime - p.xl) * iq) / p.Tq0_double_prime

    if y !== nothing
        y[4] = vd * id + vq * iq          # Pe: terminal active power (PHTRUE)
        y[5] = vq * id - vd * iq          # Qe
        y[6] = id
        y[7] = iq
        y[8] = id * sin_d + iq * cos_d    # It_Re
        y[9] = -id * cos_d + iq * sin_d   # It_Im
        y[10] = _genrou_ifd(p, x[3], x[4])
    end
    return dx
end

"""
Circuit storage in (E', psi) coordinates: `H (omega-1)^2` plus the d- and q-axis magnetic
energies (see PHPS GenRouPHS.hamiltonian).
"""
function hamiltonian(c::GENROU_PHTRUE, x, p::GenrouParams)
    omega, Eq_p, psi_d, Ed_p, psi_q = x[2], x[3], x[4], x[5], x[6]
    wb = p.omega_b
    kd = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    kq = (p.xq_double_prime - p.xl) / (p.xq_prime - p.xl)
    Eq_pp = kd * Eq_p + (1.0 - kd) * psi_d
    Ed_pp = kq * Ed_p - (1.0 - kq) * psi_q
    h_mech = p.H * (omega - 1.0)^2
    h_d = (Eq_p^2 / (p.xd - p.xd_prime)
           + (Eq_p - Eq_pp)^2 / (p.xd_prime - p.xd_double_prime)) / (2.0 * wb)
    h_q = (Ed_p^2 / (p.xq - p.xq_prime)
           + (Ed_p - Ed_pp)^2 / (p.xq_prime - p.xq_double_prime)) / (2.0 * wb)
    return h_mech + h_d + h_q
end

function grad_hamiltonian!(g, c::GENROU_PHTRUE, x, p::GenrouParams)
    omega, Eq_p, psi_d, Ed_p, psi_q = x[2], x[3], x[4], x[5], x[6]
    wb = p.omega_b
    qd_mut = (p.xd_prime - p.xd_double_prime) / (p.xd_prime - p.xl)^2
    qq_mut = (p.xq_prime - p.xq_double_prime) / (p.xq_prime - p.xl)^2
    g[1] = 0.0
    g[2] = 2.0 * p.H * (omega - 1.0)
    g[3] = (Eq_p / (p.xd - p.xd_prime) + qd_mut * (Eq_p - psi_d)) / wb
    g[4] = -qd_mut * (Eq_p - psi_d) / wb
    g[5] = (Ed_p / (p.xq - p.xq_prime) + qq_mut * (Ed_p + psi_q)) / wb
    g[6] = qq_mut * (Ed_p + psi_q) / wb
    return g
end

"""
    rotor_structure(c) -> NamedTuple or nothing

The rotor circuits of a machine as a linear port system, restating its step kernel:
`z' = A z + B_f Efd + B_s [id, iq]` for the rotor states `z = x[states]`, with the magnetic
storage `H_mag = z' Q z / 2` (the declared storage minus `H (omega - 1)^2`). `nothing` for
models whose rotor is not linear in these terms (GENSAL: saturation). The port-power audit
checks this form against the model's right-hand side at every sample.
"""
rotor_structure(::AbstractComponent) = nothing
function rotor_structure(c::GENROU_PHTRUE)
    p = c.p
    A = [-1/p.Td0_prime 0.0 0.0 0.0;
         1/p.Td0_double_prime -1/p.Td0_double_prime 0.0 0.0;
         0.0 0.0 -1/p.Tq0_prime 0.0;
         0.0 0.0 -1/p.Tq0_double_prime -1/p.Tq0_double_prime]
    Bf = [1/p.Td0_prime, 0.0, 0.0, 0.0]
    Bs = [-(p.xd - p.xd_prime)/p.Td0_prime 0.0;
          -(p.xd_prime - p.xl)/p.Td0_double_prime 0.0;
          0.0 (p.xq - p.xq_prime)/p.Tq0_prime;
          0.0 -(p.xq_prime - p.xl)/p.Tq0_double_prime]
    x = zeros(6)
    x[2] = 1.0
    Q = ForwardDiff.jacobian(z -> (xz = Vector{eltype(z)}(x); xz[3:6] .= z;
                                   grad_hamiltonian(c, xz)[3:6]), zeros(4))
    return (states = 3:6, A = A, Bf = Bf, Bs = Bs, Q = (Q + Q') / 2, currents = ("id_dq", "iq_dq"))
end

"""Norton current (network frame), independent of the bus voltage."""
function injection(c::GENROU_PHTRUE, x, V, p::GenrouParams = c.p)
    I_Re, I_Im, _, _ = _genrou_norton(p, x)
    return I_Re, I_Im
end

norton_admittance(c::GENROU_PHTRUE) = _norton_y(c.p.ra, c.p.xd_double_prime)

# PHPS YBusBuilder.add_generator_impedance: 1/(ra + j xd''), tiny z replaced by 1e-4 j.
function _norton_y(ra::Float64, xd_pp::Float64)
    z = complex(ra, xd_pp)
    abs(z) < 1e-6 && (z = complex(0.0, 0.0001))
    return _py_cdiv(1.0, z)
end
