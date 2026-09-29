# GFM_VOC_PHTRUE: grid-forming Andronov-Hopf virtual-oscillator converter (PHPS mirror of
# PowerFactory's gfm_voc_pf DLL, VI-only v1.1 build), with a one-way source reservoir.
#
# Port of phps/src/components/renewables/gfm_voc_phtrue.py at ba11ea1. Expressions keep
# PHPS's order of operations.
#
# The state is the internal EMF phasor v = v_alpha + j v_beta (no swing state, no angle
# integrator, so not a centre-of-inertia member):
#   v' = chi v + eta J (i* - i),  chi = xi (V_nom^2 - |v|^2),  J = e^{j phi},
#   i* = (p_set - j q_set) / conj(V)            (terminal voltage)
# or, with pvoc_mode = 1, the passivity-based VOC of Kong et al. (2022). Output chain: the
# predictive virtual-impedance divider (`_vi_divider`) on v.
#
# States  [v_alpha, v_beta, x_tank, Vd_meas, Vq_meas, Id_meas, Iq_meas, Id_lpf2, Iq_lpf2,
#          p_pvoc, q_pvoc]
# Inputs  [Vd, Vq, omega_ref, omega_aux]
# Outputs [Id, Iq]

struct GfmVocParams
    xi::Float64
    eta::Float64
    V_nom::Float64
    omega_n::Float64
    phi::Float64
    Zseries::Float64
    ra::Float64
    xd_double_prime::Float64
    p_set::Float64
    q_set::Float64
    r_vi::Float64
    x_vi::Float64
    i_lim::Float64
    kpr::Float64
    kpx::Float64
    i_ref_max::Float64
    pvoc_mode::Float64
    xi2::Float64
    xi3::Float64
    xi2_db::Float64
    T_pq::Float64
    Tmeas::Float64
    Tmeas_i::Float64
    T_vi_filt::Float64
    pf_frame::Float64
    coi_h_ref::Float64
    coi_h_aux::Float64
    C_TANK::Float64
end

struct GFM_VOC_PHTRUE <: AbstractComponent
    name::String
    bus::Int
    p::GfmVocParams
    params::ParamDict
end

const _GFM_VOC_STATES = ["v_alpha", "v_beta", "x_tank", "Vd_meas", "Vq_meas", "Id_meas",
                         "Iq_meas", "Id_lpf2", "Iq_lpf2", "p_pvoc", "q_pvoc"]
const _GFM_VOC_INPUTS = ["Vd", "Vq", "omega_ref", "omega_aux"]
const _GFM_VOC_OUTPUTS = ["Id", "Iq"]

model_type(::GFM_VOC_PHTRUE) = "GFM_VOC_PHTRUE"
component_role(::GFM_VOC_PHTRUE) = :generator
state_names(::GFM_VOC_PHTRUE) = _GFM_VOC_STATES
input_names(::GFM_VOC_PHTRUE) = _GFM_VOC_INPUTS
output_names(::GFM_VOC_PHTRUE) = _GFM_VOC_OUTPUTS
bus(c::GFM_VOC_PHTRUE) = c.bus

# GfmVocPHTrue.__init__ defaults, in its order (ra and xd_double_prime: component_params)
const _GFM_VOC_DEFAULTS = ("C_TANK" => 1000.0, "P_STAR" => 1.0, "PSET_REF" => 0.0,
                           "xi" => 50.0, "eta" => 20.0, "V_nom" => 1.0, "omega_n" => 376.991,
                           "phi" => 1.5707963, "p_set" => 0.5, "q_set" => 0.0, "r_vi" => 0.1,
                           "x_vi" => 0.0, "i_lim" => 10.0, "kpr" => 1.0, "kpx" => 5.0,
                           "i_con_lim" => 1.5, "i_ref_max" => 0.0, "pvoc_mode" => 0.0,
                           "xi2" => 20.0, "xi3" => 20.0, "xi2_db" => 1.0e-4, "T_pq" => 0.005,
                           "Tmeas" => 0.002, "Tmeas_i" => 0.005, "T_vi_filt" => 0.005,
                           "vnom_from_lf" => 0.0, "pf_frame" => 0.0, "coi_h_ref" => 32.0,
                           "coi_h_aux" => 10.0)

function type_defaults!(::Val{:GFM_VOC_PHTRUE}, p::ParamDict)
    for (k, v) in _GFM_VOC_DEFAULTS
        get!(p, k, v)
    end
    return p
end

function GFM_VOC_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    p = GfmVocParams(f("xi"), f("eta"), f("V_nom"), f("omega_n"), f("phi"), f("Zseries"),
                     f("ra"), f("xd_double_prime"), f("p_set"), f("q_set"), f("r_vi"),
                     f("x_vi"), f("i_lim"), f("kpr"), f("kpx"), f("i_ref_max"), f("pvoc_mode"),
                     f("xi2"), f("xi3"), f("xi2_db"), f("T_pq"), f("Tmeas"), f("Tmeas_i"),
                     f("T_vi_filt"), f("pf_frame"), f("coi_h_ref"), f("coi_h_aux"),
                     f("C_TANK"))
    return GFM_VOC_PHTRUE(name, Int(f("bus")), p, d)
end
COMPONENT_CONSTRUCTORS["GFM_VOC_PHTRUE"] = GFM_VOC_PHTRUE

@inline _voc_uout(rec::ModeRecorder, p::GfmVocParams, x) =
    _vi_divider(rec, x[1], x[2], x[4], x[5], p.Zseries, p.i_lim, p.r_vi, p.x_vi, p.kpr, p.kpx)

function _outputs!(y, c::GFM_VOC_PHTRUE, x, u, p::GfmVocParams, rec::ModeRecorder)
    uout_Re, uout_Im = _voc_uout(rec, p, x)
    y[1] = uout_Im / p.Zseries
    y[2] = -uout_Re / p.Zseries
    return y
end

function _step!(dx, y, c::GFM_VOC_PHTRUE, x, u, p::GfmVocParams, rec::ModeRecorder)
    v_alpha, v_beta = x[1], x[2]
    Vdm, Vqm = x[4], x[5]
    p_pvoc, q_pvoc = x[10], x[11]
    Idm, Iqm, Il2d, Il2q = x[6], x[7], x[8], x[9]
    V_Re, V_Im = u[1], u[2]

    # pf_frame: the reference and auxiliary machine speeds give the COI; unwired inputs (0)
    # are guarded to 1 (the feedforward is then exactly 0)
    w_ref = p.pf_frame > 0.5 ? u[3] : one(u[3])
    w_aux = p.pf_frame > 0.5 ? u[4] : one(u[4])
    w_ref = lt(rec, w_ref, 0.5) ? one(w_ref) : w_ref
    w_aux = lt(rec, w_aux, 0.5) ? one(w_aux) : w_aux
    w_coi = (p.coi_h_ref * w_ref + p.coi_h_aux * w_aux) / (p.coi_h_ref + p.coi_h_aux)
    uout_Re, uout_Im = _voc_uout(rec, p, x)

    dR = uout_Re - V_Re
    dI = uout_Im - V_Im
    I_Re = dI / p.Zseries
    I_Im = -dR / p.Zseries
    p_mea = V_Re * I_Re + V_Im * I_Im

    # bolted-terminal gate: |V| < 0.01 pu reads as zero measured current
    vsq_t = V_Re * V_Re + V_Im * V_Im
    zero_i = lt(rec, vsq_t, 1.0e-4)
    I_msr_Re = zero_i ? zero(I_Re) : I_Re
    I_msr_Im = zero_i ? zero(I_Im) : I_Im

    # current reference on the terminal voltage, optional circular limiter
    iref_a = zero(vsq_t)
    iref_b = zero(vsq_t)
    if gt(rec, vsq_t, 1.0e-4)
        iref_a = (V_Re * p.p_set + V_Im * p.q_set) / vsq_t
        iref_b = (V_Im * p.p_set - V_Re * p.q_set) / vsq_t
    end
    if p.i_ref_max > 0.0
        irm = sqrt(iref_a * iref_a + iref_b * iref_b)
        if decide(rec, _gt(irm, p.i_ref_max) && _gt(irm, 1.0e-12))
            scir = p.i_ref_max / irm
            iref_a = iref_a * scir
            iref_b = iref_b * scir
        end
    end

    # AHO core (RMS: the frame carries the rotation)
    vsq = v_alpha * v_alpha + v_beta * v_beta
    chi = p.xi * (p.V_nom * p.V_nom - vsq)
    e_a = iref_a - I_msr_Re
    e_b = iref_b - I_msr_Im
    cph = cos(p.phi)
    sph = sin(p.phi)
    Je_a = cph * e_a - sph * e_b
    Je_b = sph * e_a + cph * e_b
    if p.pvoc_mode < 0.5
        d_va = chi * v_alpha + p.eta * Je_a
        d_vb = chi * v_beta + p.eta * Je_b
    else
        # passivity-based VOC (Kong et al. 2022, eq. 6-8 and the smoothed eq. 28)
        Vr2 = p.V_nom * p.V_nom
        v2s = gt(rec, vsq, 1.0e-6) ? vsq : convert(typeof(vsq), 1.0e-6)
        q_err = p.q_set / Vr2 - q_pvoc / v2s
        p_err = p.p_set / Vr2 - p_pvoc / v2s
        s_sw = q_err * (vsq - Vr2)
        dbw = p.xi2_db > 1.0e-12 ? p.xi2_db : 1.0e-12
        xi2_eff = -abs(p.xi2) * tanh(s_sw / dbw)
        oc11 = chi + xi2_eff * q_err
        oc21 = p.xi3 * p_err
        d_va = oc11 * v_alpha - oc21 * v_beta
        d_vb = oc21 * v_alpha + oc11 * v_beta
    end
    # pf_frame feedforward (exactly 0 in the physical frame)
    d_va = d_va - (w_ref - w_coi) * p.omega_n * v_beta
    d_vb = d_vb + (w_ref - w_coi) * p.omega_n * v_alpha

    Tm = _tguard(p.Tmeas, 1e-5)
    Tmi = _tguard(p.Tmeas_i, 1e-5)
    Tv2 = _tguard(p.T_vi_filt, 1e-5)

    dx[1] = d_va
    dx[2] = d_vb
    p_res = guard_min(rec, x[3] / p.C_TANK, 1.0e-6)
    dx[3] = (p.p_set - p_mea) / p_res
    dx[4] = (V_Re - Vdm) / Tm
    dx[5] = (V_Im - Vqm) / Tm
    dx[6] = (I_msr_Re - Idm) / Tmi
    dx[7] = (I_msr_Im - Iqm) / Tmi
    dx[8] = (Idm - Il2d) / Tv2
    dx[9] = (Iqm - Il2q) / Tv2
    Tpq = _tguard(p.T_pq, 1e-5)
    q_mea = V_Im * I_Re - V_Re * I_Im
    dx[10] = (p_mea - p_pvoc) / Tpq
    dx[11] = (q_mea - q_pvoc) / Tpq
    return dx
end

hamiltonian(c::GFM_VOC_PHTRUE, x, p::GfmVocParams) =
    0.5 * (x[1] * x[1] + x[2] * x[2]) / p.eta + 0.5 * x[3] * x[3] / p.C_TANK

function grad_hamiltonian!(g, c::GFM_VOC_PHTRUE, x, p::GfmVocParams)
    fill!(g, 0.0)
    g[1] = x[1] / p.eta
    g[2] = x[2] / p.eta
    g[3] = x[3] / p.C_TANK
    return g
end

"""Norton current (network frame): state-only, independent of the bus voltage."""
function injection(c::GFM_VOC_PHTRUE, x, V, p::GfmVocParams = c.p)
    uout_Re, uout_Im = _voc_uout(NoModes(), p, x)
    return uout_Im / p.Zseries, -uout_Re / p.Zseries
end

norton_admittance(c::GFM_VOC_PHTRUE) = _norton_y(c.p.ra, c.p.xd_double_prime)

function lag_states(c::GFM_VOC_PHTRUE)
    p = c.p
    Tm, Tmi = _tguard(p.Tmeas, 1e-5), _tguard(p.Tmeas_i, 1e-5)
    Tv2, Tpq = _tguard(p.T_vi_filt, 1e-5), _tguard(p.T_pq, 1e-5)
    return [(4, Tm), (5, Tm), (6, Tmi), (7, Tmi), (8, Tv2), (9, Tv2), (10, Tpq), (11, Tpq)]
end
