# Observables: the derived signals PHPS logs after each component's states in
# simulation_results.csv (the `observables` property of each PHPS model, in its order, with
# its C++ expressions). `u` and `y` are the component's inputs and outputs (after the step
# kernel) at the logged state.

"""
    observable_names(c) -> Vector{String}

Names of the component's logged signals, in PHPS's column order.
"""
function observable_names end

"""
    observable_values!(v, c, x, u, y) -> v

Values of [`observable_names`](@ref) at states `x`, inputs `u` and outputs `y`.
"""
function observable_values! end

# PHPS writes the rotor angle in degrees with this value of pi
const _PHPS_PI = 3.14159265359

# -- GENROU_PHTRUE ----------------------------------------------------------------------
observable_names(::GENROU_PHTRUE) =
    ["delta_deg", "Te", "Tm_in", "Pe", "Qe", "V_term", "Eq_p", "omega", "H_total", "i_fd"]

function observable_values!(v, c::GENROU_PHTRUE, x, u, y)
    p = c.p
    v[1] = x[1] * 180.0 / _PHPS_PI
    v[2] = y[4] / x[2]                        # Te = P / omega (PHTRUE)
    v[3] = u[3]                               # Tm_in
    v[4] = y[4]                               # Pe = terminal active power
    v[5] = y[5]                               # Qe
    v[6] = sqrt(u[1] * u[1] + u[2] * u[2])    # V_term
    v[7] = x[3]                               # Eq_p
    v[8] = x[2]                               # omega
    v[9] = hamiltonian(c, x, p)               # H_total
    v[10] = _genrou_ifd(p, x[3], x[4])        # i_fd
    return v
end

# -- GENSAL_PHTRUE ----------------------------------------------------------------------
observable_names(::GENSAL_PHTRUE) = ["delta_deg", "Te", "Pe", "Qe", "V_term", "H_total", "i_fd"]

function observable_values!(v, c::GENSAL_PHTRUE, x, u, y)
    p = c.p
    v[1] = x[1] * 180.0 / _PHPS_PI
    v[2] = (y[4] + p.ra * (y[6] * y[6] + y[7] * y[7])) / x[2]   # Te = P_ag / omega
    v[3] = y[4]
    v[4] = y[5]
    v[5] = sqrt(u[1] * u[1] + u[2] * u[2])
    v[6] = hamiltonian(c, x, p)
    v[7] = _gensal_ifd(p, x[3], x[4])
    return v
end

# -- IEEET1_PHTRUE ----------------------------------------------------------------------
observable_names(::IEEET1_PHTRUE) = ["Efd", "Vr", "H_field", "p_field"]

function observable_values!(v, c::IEEET1_PHTRUE, x, u, y)
    p = c.p
    v[1] = x[3]
    v[2] = x[2] > p.VRMAX ? p.VRMAX : (x[2] < p.VRMIN ? p.VRMIN : x[2])   # limited Vr
    v[3] = 0.5 * x[5] * x[5] / p.C_FIELD
    v[4] = x[5] / p.C_FIELD
    return v
end

# -- IEEEG1_PHTRUE ----------------------------------------------------------------------
observable_names(::IEEEG1_PHTRUE) = ["Tm", "Valve", "H_steam", "p_tank"]

function observable_values!(v, c::IEEEG1_PHTRUE, x, u, y)
    p = c.p
    v[1] = _ieeeg1_tm(p, x)
    v[2] = x[2]
    v[3] = 0.5 * x[7] * x[7] / p.C_TANK
    v[4] = x[7] / p.C_TANK
    return v
end

# -- IEEEG3_PHTRUE ----------------------------------------------------------------------
observable_names(::IEEEG3_PHTRUE) =
    ["Tm", "Gate", "GateState", "GateRate", "Dashpot", "H_water", "h_head"]

function observable_values!(v, c::IEEEG3_PHTRUE, x, u, y)
    p = c.p
    gmin, gmax = p.PMIN * p.R_base, p.PMAX * p.R_base
    rdn, rup = p.UC * p.R_base, p.UO * p.R_base
    gate = min(max(x[3], gmin), gmax)
    v[1] = p.a23 * (x[4] + p.Tb_t * ((gate - x[4]) / p.Ta_t))
    v[2] = gate
    v[3] = x[3]
    v[4] = min(max(x[1], rdn), rup)
    v[5] = p.Delta * (gate - x[2])
    v[6] = 0.5 * x[5] * x[5] / p.C_TANK
    v[7] = x[5] / p.C_TANK
    return v
end

# -- COMPLEXLOAD -------------------------------------------------------------------------
observable_names(::COMPLEXLOAD) = ["Pload", "Qload"]

function observable_values!(v, c::COMPLEXLOAD, x, u, y)
    v[1] = y[3]
    v[2] = y[4]
    return v
end

# -- grid-forming converters -------------------------------------------------------------
# (their C++ expressions use this value of pi)
const _PHPS_PI_GFM = 3.14159265358979

observable_names(::GFM_VSM_PHTRUE) =
    ["theta_deg", "omega", "u_mag", "i_mag", "H_kin", "H_damp", "H_tank"]

function observable_values!(v, c::GFM_VSM_PHTRUE, x, u, y)
    p = c.p
    v[1] = x[1] * 180.0 / _PHPS_PI_GFM
    v[2] = x[2]
    v[3] = x[4]
    v[4] = sqrt(x[8] * x[8] + x[9] * x[9])
    v[5] = 0.5 * p.Ta * (x[2] - p.f_set) * (x[2] - p.f_set)
    v[6] = 0.5 * p.Dp / p.omega_c * (x[3] - p.f_set) * (x[3] - p.f_set)
    v[7] = 0.5 * x[5] * x[5] / p.C_TANK
    return v
end

observable_names(::GFM_DROOP_PHTRUE) =
    ["theta_deg", "omega", "q_lpf", "u_mag", "i_mag", "H_kin", "H_tank"]

function observable_values!(v, c::GFM_DROOP_PHTRUE, x, u, y)
    p = c.p
    v[1] = x[1] * 180.0 / _PHPS_PI_GFM
    v[2] = x[2]
    v[3] = x[3]
    v[4] = x[4]
    v[5] = sqrt(x[8] * x[8] + x[9] * x[9])
    v[6] = 0.5 * p.Ta * (x[2] - p.f_set) * (x[2] - p.f_set)
    v[7] = 0.5 * x[5] * x[5] / p.C_TANK
    return v
end

observable_names(::GFM_VOC_PHTRUE) =
    ["theta_deg", "v_mag", "chi", "i_mag", "H_osc", "H_orbit", "H_tank"]

function observable_values!(v, c::GFM_VOC_PHTRUE, x, u, y)
    p = c.p
    s = x[1] * x[1] + x[2] * x[2]
    v[1] = atan(x[2], x[1]) * 180.0 / _PHPS_PI_GFM
    v[2] = sqrt(s)
    v[3] = p.xi * (p.V_nom * p.V_nom - s)
    v[4] = sqrt(x[6] * x[6] + x[7] * x[7])
    v[5] = 0.5 * s / p.eta
    v[6] = 0.25 * (s - p.V_nom * p.V_nom) * (s - p.V_nom * p.V_nom) / (p.eta * p.xi)
    v[7] = 0.5 * x[3] * x[3] / p.C_TANK
    return v
end
