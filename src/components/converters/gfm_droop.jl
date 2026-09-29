# GFM_DROOP_PHTRUE: grid-forming droop converter (PHPS mirror of PowerFactory's gfm_droop_pf
# DLL), with a one-way source reservoir.
#
# Port of phps/src/components/renewables/gfm_droop_phtrue.py at ba11ea1 (a subclass of
# GfmVsmPHTrue: the VSM defaults apply, the output chain, dynamics and states are its own).
# Expressions keep PHPS's order of operations.
#
# D'Arco form of the P-f droop, Ta = Tf_p/mp:
#   Ta_e omega' = (p_set - p_mea) - (omega - f_set)/mp_e   (+ the mp_dot correction with the
#                                                           Li et al. adaptive droop)
#   q_lpf'      = (q_mea - q_lpf)/Tf_q
#   u_mag'      = (u_set - mq_e (q_lpf - q_set) - u_mag)/Tv
#   theta'      = (omega - f_set + (w_ref - 1)) omega_n
# Output chain: the split-path virtual impedance (R on the measured current, X on its
# high-pass), the overcurrent boost and the two-stage voltage-source limiter.
#
# States  [theta, omega, q_lpf, u_mag, x_tank, Vd_meas, Vq_meas, Id_meas, Iq_meas, Id_lpf2,
#          Iq_lpf2, Sa, Id_lpfz, Iq_lpfz]
# Inputs  [Vd, Vq, omega_ref, S_avail]
# Outputs [Id, Iq, omega]

struct GfmDroopParams
    mp::Float64
    mq::Float64
    Tf_p::Float64
    Tf_q::Float64
    omega_n::Float64
    f_set::Float64
    Tv::Float64
    Zseries::Float64
    ra::Float64
    xd_double_prime::Float64
    r_vi::Float64
    x_vi::Float64
    i_lim::Float64
    kpr::Float64
    kpx::Float64
    i_con_lim::Float64
    p_set::Float64
    q_set::Float64
    u_set::Float64
    Ta::Float64
    adapt_droop::Float64
    dw_droop::Float64
    dv_droop::Float64
    S_min::Float64
    T_sa::Float64
    adapt_vi::Float64
    a1::Float64
    a0::Float64
    vi_mode::Float64
    T_vz::Float64
    T_vi_filt::Float64
    pf_frame::Float64
    C_TANK::Float64
    Tmeas::Float64
    Tmeas_i::Float64
end

struct GFM_DROOP_PHTRUE <: AbstractComponent
    name::String
    bus::Int
    p::GfmDroopParams
    params::ParamDict
end

const _GFM_DROOP_STATES = ["theta", "omega", "q_lpf", "u_mag", "x_tank", "Vd_meas", "Vq_meas",
                           "Id_meas", "Iq_meas", "Id_lpf2", "Iq_lpf2", "Sa", "Id_lpfz",
                           "Iq_lpfz"]
const _GFM_DROOP_INPUTS = ["Vd", "Vq", "omega_ref", "S_avail"]

model_type(::GFM_DROOP_PHTRUE) = "GFM_DROOP_PHTRUE"
component_role(::GFM_DROOP_PHTRUE) = :generator
state_names(::GFM_DROOP_PHTRUE) = _GFM_DROOP_STATES
input_names(::GFM_DROOP_PHTRUE) = _GFM_DROOP_INPUTS
output_names(::GFM_DROOP_PHTRUE) = _CONVERTER_OUTPUTS
bus(c::GFM_DROOP_PHTRUE) = c.bus

# GfmDroopPHTrue.__init__ defaults (its own first, then GfmVsmPHTrue's)
function type_defaults!(::Val{:GFM_DROOP_PHTRUE}, p::ParamDict)
    for (k, v) in ("mp" => 0.05, "mq" => 0.05, "Tf_p" => 0.1, "Tf_q" => 0.1, "q_set" => 0.0,
                   "omega_n" => 376.991, "f_set" => 1.0, "Tv" => 0.02, "T_vi_filt" => 0.005)
        get!(p, k, v)
    end
    get!(p, "Ta", param_value(p["Tf_p"]) / param_value(p["mp"]))
    for (k, v) in ("adapt_droop" => 0.0, "dw_droop" => 0.05, "dv_droop" => 0.05,
                   "S_min" => 0.05, "T_sa" => 0.02, "adapt_vi" => 0.0, "a1" => 0.036,
                   "a0" => -0.0115, "vi_mode" => 0.0, "T_vz" => 0.001)
        get!(p, k, v)
    end
    get!(p, "Dp", 1.0 / param_value(p["mp"]))
    return type_defaults!(Val(:GFM_VSM_PHTRUE), p)
end

function GFM_DROOP_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    p = GfmDroopParams(f("mp"), f("mq"), f("Tf_p"), f("Tf_q"), f("omega_n"), f("f_set"),
                       f("Tv"), f("Zseries"), f("ra"), f("xd_double_prime"), f("r_vi"),
                       f("x_vi"), f("i_lim"), f("kpr"), f("kpx"), f("i_con_lim"), f("p_set"),
                       f("q_set"), f("u_set"), f("Ta"), f("adapt_droop"), f("dw_droop"),
                       f("dv_droop"), f("S_min"), f("T_sa"), f("adapt_vi"), f("a1"), f("a0"),
                       f("vi_mode"), f("T_vz"), f("T_vi_filt"), f("pf_frame"), f("C_TANK"),
                       f("Tmeas"), f("Tmeas_i"))
    return GFM_DROOP_PHTRUE(name, Int(f("bus")), p, d)
end
COMPONENT_CONSTRUCTORS["GFM_DROOP_PHTRUE"] = GFM_DROOP_PHTRUE

# PHPS _limiter_chain_cpp (droop): virtual impedance, overcurrent boost, then Stage 1
# (|u_con - V_meas| clipped to Zseries i_con_lim) and Stage 2 (measured |i| clamp).
@inline function _droop_uout(rec::ModeRecorder, p::GfmDroopParams, x)
    T = eltype(x)
    theta, u_mag = x[1], x[4]
    Vdm, Vqm, Idm, Iqm, Il2d, Il2q = x[6], x[7], x[8], x[9], x[10], x[11]
    sin_t = sin(theta)
    cos_t = cos(theta)
    E_Re = u_mag * cos_t
    E_Im = u_mag * sin_t

    ih_Re = Idm - Il2d
    ih_Im = Iqm - Il2q
    r_vi_u = convert(T, p.r_vi)
    x_vi_u = convert(T, p.x_vi)
    if p.adapt_vi > 0.5
        Sa_l = x[12]
        Smn = p.S_min > 1e-3 ? p.S_min : 1e-3
        Sa_gl = gt(rec, Sa_l, Smn) ? Sa_l : convert(T, Smn)
        Rv = p.a1 / Sa_gl + p.a0
        Rv = lt(rec, Rv, 0.0) ? zero(T) : Rv
        r_vi_u = p.r_vi + Rv
        x_vi_u = p.x_vi + Rv
    end
    if p.vi_mode > 0.5
        Izd = x[13]
        Izq = x[14]
        du_Re = r_vi_u * Izd - x_vi_u * Izq
        du_Im = r_vi_u * Izq + x_vi_u * Izd
    else
        du_Re = r_vi_u * Idm - x_vi_u * ih_Im
        du_Im = r_vi_u * Iqm + x_vi_u * ih_Re
    end
    imag_m = sqrt(Idm * Idm + Iqm * Iqm)
    if gt(rec, imag_m, p.i_lim)
        dri = p.kpr * (imag_m - p.i_lim)
        dxi = p.kpx * (imag_m - p.i_lim)
        du_Re += dri * Idm - dxi * Iqm
        du_Im += dri * Iqm + dxi * Idm
    end
    ucon_Re = E_Re - du_Re
    ucon_Im = E_Im - du_Im

    # Stage 1
    dvr = ucon_Re - Vdm
    dvi = ucon_Im - Vqm
    dv_mag = sqrt(dvr * dvr + dvi * dvi)
    dv_max = p.Zseries * p.i_con_lim
    if decide(rec, _gt(dv_mag, dv_max) && _gt(dv_mag, 1.0e-9))
        sc1 = dv_max / dv_mag
        ucon_Re = Vdm + dvr * sc1
        ucon_Im = Vqm + dvi * sc1
    end
    # Stage 2
    if decide(rec, _gt(imag_m, p.i_con_lim) && _gt(imag_m, 1.0e-9))
        sc2 = p.i_con_lim / imag_m
        ucon_Re = Vdm + (ucon_Re - Vdm) * sc2
        ucon_Im = Vqm + (ucon_Im - Vqm) * sc2
    end
    return ucon_Re, ucon_Im
end

function _outputs!(y, c::GFM_DROOP_PHTRUE, x, u, p::GfmDroopParams, rec::ModeRecorder)
    uout_Re, uout_Im = _droop_uout(rec, p, x)
    y[1] = uout_Im / p.Zseries
    y[2] = -uout_Re / p.Zseries
    y[3] = x[2]
    return y
end

function _step!(dx, y, c::GFM_DROOP_PHTRUE, x, u, p::GfmDroopParams, rec::ModeRecorder)
    omega, q_lpf, u_mag = x[2], x[3], x[4]
    Vdm, Vqm, Idm, Iqm, Il2d, Il2q = x[6], x[7], x[8], x[9], x[10], x[11]
    V_Re, V_Im = u[1], u[2]
    w_ref = _wref(rec, p.pf_frame, u[3])
    # available capacity S_a: an unwired input (0) is guarded to 1 (rated)
    Sa_in_g = le(rec, u[4], 0.0) ? one(u[4]) : u[4]
    Tsa_g = p.T_sa > 1e-4 ? p.T_sa : 1e-4
    uout_Re, uout_Im = _droop_uout(rec, p, x)

    dR = uout_Re - V_Re
    dI = uout_Im - V_Im
    I_Re = dI / p.Zseries
    I_Im = -dR / p.Zseries
    p_mea = V_Re * I_Re + V_Im * I_Im
    q_mea = V_Im * I_Re - V_Re * I_Im

    # Li et al. adaptive droop (parameter switch), on the guarded S_a state
    Sa_s = x[12]
    Smin = p.S_min > 1e-3 ? p.S_min : 1e-3
    Sa_g = gt(rec, Sa_s, Smin) ? Sa_s : convert(typeof(Sa_s), Smin)
    mp_e = p.adapt_droop > 0.5 ? p.dw_droop / Sa_g : convert(typeof(Sa_g), p.mp)
    mq_e = p.adapt_droop > 0.5 ? p.dv_droop / Sa_g : convert(typeof(Sa_g), p.mq)
    Ta_e = p.Tf_p / (gt(rec, mp_e, 1e-9) ? mp_e : convert(typeof(mp_e), 1e-9))

    # droop dynamics (D'Arco form) with the mp_dot correction of the adaptive droop
    d_Sa = (Sa_in_g - Sa_s) / Tsa_g
    d_omega = ((p.p_set - p_mea) - (omega - p.f_set) / mp_e) / Ta_e
    if decide(rec, p.adapt_droop > 0.5 && _gt(Sa_s, Smin))
        d_omega += (d_Sa / Sa_g) * (p.f_set - omega)
    end
    d_qlpf = (q_mea - q_lpf) / p.Tf_q
    u_mag_set = p.u_set - mq_e * (q_lpf - p.q_set)
    d_umag = (u_mag_set - u_mag) / p.Tv
    d_theta = (omega - p.f_set + (w_ref - 1.0)) * p.omega_n

    Tm = _tguard(p.Tmeas, 1e-5)
    Tmi = _tguard(p.Tmeas_i, 1e-5)
    Tv2 = _tguard(p.T_vi_filt, 1e-5)

    dx[1] = d_theta
    dx[2] = d_omega
    dx[3] = d_qlpf
    dx[4] = d_umag
    p_res = guard_min(rec, x[5] / p.C_TANK, 1.0e-6)
    dx[5] = (p.p_set - p_mea) / p_res
    dx[6] = (V_Re - Vdm) / Tm
    dx[7] = (V_Im - Vqm) / Tm
    dx[8] = (I_Re - Idm) / Tmi
    dx[9] = (I_Im - Iqm) / Tmi
    dx[10] = (Idm - Il2d) / Tv2
    dx[11] = (Iqm - Il2q) / Tv2
    dx[12] = d_Sa
    Tvz_g = _tguard(p.T_vz, 1e-5)
    dx[13] = (I_Re - x[13]) / Tvz_g
    dx[14] = (I_Im - x[14]) / Tvz_g
    return dx
end

# PHPS _Ta_eff: the parameter Ta, or Tf_p/mp_eff with the adaptive droop
function _droop_ta_eff(p::GfmDroopParams, x)
    p.adapt_droop <= 0.5 && return p.Ta
    S_min = max(p.S_min, 1.0e-3)
    Sa = max(x[12], S_min)
    mp_e = p.dw_droop / Sa
    return p.Tf_p / max(mp_e, 1.0e-9)
end

hamiltonian(c::GFM_DROOP_PHTRUE, x, p::GfmDroopParams) =
    0.5 * _droop_ta_eff(p, x) * (x[2] - p.f_set)^2 + 0.5 * x[5] * x[5] / p.C_TANK

function grad_hamiltonian!(g, c::GFM_DROOP_PHTRUE, x, p::GfmDroopParams)
    fill!(g, 0.0)
    g[2] = _droop_ta_eff(p, x) * (x[2] - p.f_set)
    g[5] = x[5] / p.C_TANK
    if p.adapt_droop > 0.5
        S_min = max(p.S_min, 1.0e-3)
        dw = max(p.dw_droop, 1.0e-12)
        x[12] > S_min && (g[12] = 0.5 * (p.Tf_p / dw) * (x[2] - p.f_set)^2)
    end
    return g
end

"""Norton current (network frame): state-only, independent of the bus voltage."""
function injection(c::GFM_DROOP_PHTRUE, x, V, p::GfmDroopParams = c.p)
    uout_Re, uout_Im = _droop_uout(NoModes(), p, x)
    return uout_Im / p.Zseries, -uout_Re / p.Zseries
end

norton_admittance(c::GFM_DROOP_PHTRUE) = _norton_y(c.p.ra, c.p.xd_double_prime)

function lag_states(c::GFM_DROOP_PHTRUE)
    p = c.p
    Tm, Tmi = _tguard(p.Tmeas, 1e-5), _tguard(p.Tmeas_i, 1e-5)
    Tv2, Tvz = _tguard(p.T_vi_filt, 1e-5), _tguard(p.T_vz, 1e-5)
    Tsa = p.T_sa > 1e-4 ? p.T_sa : 1e-4
    return [(6, Tm), (7, Tm), (8, Tmi), (9, Tmi), (10, Tv2), (11, Tv2), (12, Tsa), (13, Tvz),
            (14, Tvz)]
end
