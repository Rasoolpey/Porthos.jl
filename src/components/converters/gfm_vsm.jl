# GFM_VSM_PHTRUE: grid-forming virtual synchronous machine (PHPS mirror of PowerFactory's
# gfm_vsm_pf DLL, VI-only build), with a one-way source reservoir.
#
# Port of phps/src/components/renewables/gfm_vsm_phtrue.py at ba11ea1. Expressions keep
# PHPS's order of operations.
#
# An EMF u_mag at angle theta behind Zseries (Norton admittance 1/(j Zseries) in the Y-bus):
#   Ta omega'  = (p_set - p_mea) - Dp (omega_f - f_set) - D_direct (omega - f_set)
#   omega_f'   = omega_c (omega - omega_f)                (damping low-pass)
#   theta'     = (omega - f_set + (w_ref - 1)) omega_n
#   u_mag'     = (u_set - u_mag) / Tv                     (K_field = 0; the Natarajan-Weiss
#                                                         field controller when K_field > 0)
# with measured voltage and current lags, the predictive virtual-impedance divider
# (`_vi_divider`) and the reservoir x_tank' = (p_set - p_mea) / max(x_tank/C, 1e-6).
#
# States  [theta, omega, omega_f, u_mag, x_tank, Vd_meas, Vq_meas, Id_meas, Iq_meas, u_field]
# Inputs  [Vd, Vq, omega_ref]
# Outputs [Id, Iq, omega]   (Norton current u_out / (j Zseries), network frame)

struct GfmVsmParams
    Ta::Float64
    Dp::Float64
    omega_c::Float64
    omega_n::Float64
    f_set::Float64
    Tv::Float64
    Zseries::Float64
    ra::Float64
    xd_double_prime::Float64
    D_direct::Float64
    i_lim::Float64
    r_vi::Float64
    x_vi::Float64
    kpr::Float64
    kpx::Float64
    p_set::Float64
    u_set::Float64
    q_set::Float64
    K_field::Float64
    Dq::Float64
    v_set::Float64
    uf_min::Float64
    uf_max::Float64
    uf_band::Float64
    n_vi::Float64
    pf_frame::Float64
    C_TANK::Float64
    P_STAR::Float64
    Tmeas::Float64
    Tmeas_i::Float64
end

struct GFM_VSM_PHTRUE <: AbstractComponent
    name::String
    bus::Int
    p::GfmVsmParams
    params::ParamDict
end

const _GFM_VSM_STATES = ["theta", "omega", "omega_f", "u_mag", "x_tank", "Vd_meas", "Vq_meas",
                         "Id_meas", "Iq_meas", "u_field"]
const _GFM_VSM_INPUTS = ["Vd", "Vq", "omega_ref"]

model_type(::GFM_VSM_PHTRUE) = "GFM_VSM_PHTRUE"
component_role(::GFM_VSM_PHTRUE) = :generator
state_names(::GFM_VSM_PHTRUE) = _GFM_VSM_STATES
input_names(::GFM_VSM_PHTRUE) = _GFM_VSM_INPUTS
output_names(::GFM_VSM_PHTRUE) = _CONVERTER_OUTPUTS
bus(c::GFM_VSM_PHTRUE) = c.bus

# GfmVsmPHTrue.__init__ defaults, in its order (ra and xd_double_prime: component_params)
const _GFM_VSM_DEFAULTS = ("C_TANK" => 1000.0, "P_STAR" => 1.0, "PSET_REF" => 0.0,
                           "p_set" => 0.5, "u_set" => 1.0, "i_con_lim" => 1.5, "r_vi" => 0.0,
                           "x_vi" => 0.0, "i_lim" => 1.1, "kpr" => 1.0, "kpx" => 5.0,
                           "D_direct" => 0.0, "Tmeas" => 0.002, "Tmeas_i" => 0.005,
                           "q_set" => 0.0, "K_field" => 0.0, "Dq" => 0.0, "v_set" => 0.0,
                           "uf_min" => 0.1, "uf_max" => 3.0, "uf_band" => 1.0e-3,
                           "n_vi" => 1.0, "pf_frame" => 0.0)

function type_defaults!(::Val{:GFM_VSM_PHTRUE}, p::ParamDict)
    for (k, v) in _GFM_VSM_DEFAULTS
        get!(p, k, v)
    end
    return p
end

function GfmVsmParams(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    return GfmVsmParams(f("Ta"), f("Dp"), f("omega_c"), f("omega_n"), f("f_set"), f("Tv"),
                        f("Zseries"), f("ra"), f("xd_double_prime"), f("D_direct"), f("i_lim"),
                        f("r_vi"), f("x_vi"), f("kpr"), f("kpx"), f("p_set"), f("u_set"),
                        f("q_set"), f("K_field"), f("Dq"), f("v_set"), f("uf_min"),
                        f("uf_max"), f("uf_band"), f("n_vi"), f("pf_frame"), f("C_TANK"),
                        f("P_STAR"), f("Tmeas"), f("Tmeas_i"))
end

GFM_VSM_PHTRUE(name::String, d::ParamDict) =
    GFM_VSM_PHTRUE(name, Int(_p(d, "bus", name)), GfmVsmParams(name, d), d)
COMPONENT_CONSTRUCTORS["GFM_VSM_PHTRUE"] = GFM_VSM_PHTRUE

# PHPS _limiter_chain_cpp (VSM): the Mod 3 virtual inductor (a parameter switch), then the
# predictive divider on the measured voltage.
@inline function _vsm_uout(rec::ModeRecorder, p::GfmVsmParams, theta, u_mag, Vdm, Vqm)
    sin_t = sin(theta)
    cos_t = cos(theta)
    E_Re = u_mag * cos_t
    E_Im = u_mag * sin_t
    if p.n_vi > 1.0
        E_Re = Vdm + (E_Re - Vdm) / p.n_vi
        E_Im = Vqm + (E_Im - Vqm) / p.n_vi
    end
    return _vi_divider(rec, E_Re, E_Im, Vdm, Vqm, p.Zseries, p.i_lim, p.r_vi, p.x_vi, p.kpr,
                       p.kpx)
end

# e = Mf if omega when the field controller is on (K_field > 0), else the LPF state
@inline _vsm_umag(p::GfmVsmParams, x) = p.K_field > 0.0 ? x[10] * x[2] : x[4]

function _outputs!(y, c::GFM_VSM_PHTRUE, x, u, p::GfmVsmParams, rec::ModeRecorder)
    uout_Re, uout_Im = _vsm_uout(rec, p, x[1], _vsm_umag(p, x), x[6], x[7])
    y[1] = uout_Im / p.Zseries
    y[2] = -uout_Re / p.Zseries
    y[3] = x[2]
    return y
end

function _step!(dx, y, c::GFM_VSM_PHTRUE, x, u, p::GfmVsmParams, rec::ModeRecorder)
    theta, omega, omega_f, u_field = x[1], x[2], x[3], x[10]
    u_mag = _vsm_umag(p, x)
    Vdm, Vqm, Idm, Iqm = x[6], x[7], x[8], x[9]
    V_Re, V_Im = u[1], u[2]
    w_ref = _wref(rec, p.pf_frame, u[3])
    uout_Re, uout_Im = _vsm_uout(rec, p, theta, u_mag, Vdm, Vqm)

    # actual injected current (u_out - V) / (j Zseries) and delivered power
    dR = uout_Re - V_Re
    dI = uout_Im - V_Im
    I_Re = dI / p.Zseries
    I_Im = -dR / p.Zseries
    p_mea = V_Re * I_Re + V_Im * I_Im

    # VSM controller: swing, damping low-pass, angle
    d_omega = ((p.p_set - p_mea) - p.Dp * (omega_f - p.f_set)
               - p.D_direct * (omega - p.f_set)) / p.Ta
    d_omega_f = p.omega_c * (omega - omega_f)
    d_theta = (omega - p.f_set + (w_ref - 1.0)) * p.omega_n

    # Natarajan-Weiss field controller with its smoothstep projection (K_field > 0 only)
    q_mea_vsm = V_Im * I_Re - V_Re * I_Im
    v_mag_t = sqrt(V_Re * V_Re + V_Im * V_Im)
    d_ufield = zero(p_mea)
    if p.K_field > 0.0
        vset_eff = p.v_set > 1.0e-6 ? p.v_set : 1.0
        w_field = (p.q_set - q_mea_vsm + p.Dq * (vset_eff - v_mag_t)) / p.K_field
        bw = p.uf_band > 1.0e-9 ? p.uf_band : 1.0e-9
        t_hi = (p.uf_max - u_field) / bw
        t_lo = (u_field - p.uf_min) / bw
        t_hi = lt(rec, t_hi, 0.0) ? zero(t_hi) : t_hi
        t_hi = gt(rec, t_hi, 1.0) ? one(t_hi) : t_hi
        t_lo = lt(rec, t_lo, 0.0) ? zero(t_lo) : t_lo
        t_lo = gt(rec, t_lo, 1.0) ? one(t_lo) : t_lo
        gate_hi = t_hi * t_hi * (3.0 - 2.0 * t_hi)
        gate_lo = t_lo * t_lo * (3.0 - 2.0 * t_lo)
        d_ufield = gt(rec, w_field, 0.0) ? w_field * gate_hi : w_field * gate_lo
    end
    d_umag = p.K_field > 0.0 ? d_ufield * omega + u_field * d_omega : (p.u_set - u_mag) / p.Tv

    # measurement lags
    Tm = _tguard(p.Tmeas, 1e-5)
    Tmi = _tguard(p.Tmeas_i, 1e-5)

    dx[1] = d_theta
    dx[2] = d_omega
    dx[3] = d_omega_f
    dx[4] = d_umag
    p_res = guard_min(rec, x[5] / p.C_TANK, 1.0e-6)
    dx[5] = (p.p_set - p_mea) / p_res
    dx[6] = (V_Re - Vdm) / Tm
    dx[7] = (V_Im - Vqm) / Tm
    dx[8] = (I_Re - Idm) / Tmi
    dx[9] = (I_Im - Iqm) / Tmi
    dx[10] = d_ufield
    return dx
end

function hamiltonian(c::GFM_VSM_PHTRUE, x, p::GfmVsmParams)
    H_damp = p.Dp <= 0.0 ? zero(eltype(x)) : 0.5 * p.Dp / p.omega_c * (x[3] - p.f_set)^2
    H_field = zero(eltype(x))
    if p.K_field > 0.0
        uf_ref = p.u_set / (p.f_set > 1e-6 ? p.f_set : 1.0)
        H_field = 0.5 * p.K_field * (x[10] - uf_ref)^2
    end
    return 0.5 * p.Ta * (x[2] - p.f_set)^2 + H_damp + 0.5 * x[5] * x[5] / p.C_TANK + H_field
end

function grad_hamiltonian!(g, c::GFM_VSM_PHTRUE, x, p::GfmVsmParams)
    fill!(g, 0.0)
    g[2] = p.Ta * (x[2] - p.f_set)
    p.Dp > 0.0 && (g[3] = p.Dp / p.omega_c * (x[3] - p.f_set))
    g[5] = x[5] / p.C_TANK
    if p.K_field > 0.0
        uf_ref = p.u_set / (p.f_set > 1e-6 ? p.f_set : 1.0)
        g[10] = p.K_field * (x[10] - uf_ref)
    end
    return g
end

"""Norton current (network frame): state-only, independent of the bus voltage."""
function injection(c::GFM_VSM_PHTRUE, x, V, p::GfmVsmParams = c.p)
    uout_Re, uout_Im = _vsm_uout(NoModes(), p, x[1], _vsm_umag(p, x), x[6], x[7])
    return uout_Im / p.Zseries, -uout_Re / p.Zseries
end

norton_admittance(c::GFM_VSM_PHTRUE) = _norton_y(c.p.ra, c.p.xd_double_prime)

lag_states(c::GFM_VSM_PHTRUE) =
    [(j, _tguard(j <= 7 ? c.p.Tmeas : c.p.Tmeas_i, 1e-5)) for j in 6:9]
