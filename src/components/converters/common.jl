# Shared pieces of the grid-forming converter models (GFM_VSM_PHTRUE, GFM_DROOP_PHTRUE,
# GFM_VOC_PHTRUE): an internal EMF behind the series reactance Zseries, Norton-injected into
# the Y-bus (`xd_double_prime = Zseries`, `ra = 0`), with a one-way source reservoir.

const _CONVERTER_OUTPUTS = ["Id", "Iq", "omega"]

"""
    vi_bisection(du_mag, Zseries, i_lim, r_vi, x_vi, kpx, kpr) -> i

PHPS's 60-step bisection for the current magnitude of the predictive virtual impedance: the
root of the monotone equation i |Z_tot(i)| = du_mag, Z_tot(i) = r_eff(i) + j(Zseries +
x_eff(i)), searched on [0, du_mag / Zseries], with the comparisons of PHPS's C++
(`mid > i_lim` for the boost, `g_sq < du_sq` for the side). A numerical root solve, not a
limiter: its comparisons are not branch sites (the late ones sit on the root by
construction, so they are not recorded or matched against PHPS); the regime switch is the
`ipred_mag > i_lim` site after it. The comparisons go through the branch primitives with
the non-recording `NoModes()` recorder; a validated root enclosure for interval arguments
belongs to P11.
"""
@inline function vi_bisection(du_mag, Zseries, i_lim, r_vi, x_vi, kpx, kpr)
    lo = zero(du_mag)
    hi = du_mag / Zseries
    du_sq = du_mag * du_mag
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        ov_k = gt(NoModes(), mid, i_lim) ? mid - i_lim : zero(mid)
        re_k = r_vi + kpr * ov_k
        xe_k = Zseries + x_vi + kpx * ov_k
        g_sq = mid * mid * (re_k * re_k + xe_k * xe_k)
        if lt(NoModes(), g_sq, du_sq)
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

# Derivatives of the root: differentiating through the bisection gives the secant slope
# i/du_mag, right only while the virtual impedance is off. The root's value is the
# bisection's (bit for bit), its derivative the implicit-function one of
# F(i, du) = i^2 |Z_tot(i)|^2 - du^2 = 0 (PHPS's squared form): di/ddu = 2 du / F_i.
function vi_bisection(du_mag::ForwardDiff.Dual{T}, Zseries, i_lim, r_vi, x_vi, kpx,
                      kpr) where {T}
    d = ForwardDiff.value(du_mag)
    i = vi_bisection(d, Zseries, i_lim, r_vi, x_vi, kpx, kpr)
    on = gt(NoModes(), i, i_lim)
    ov = on ? i - i_lim : zero(i)
    re = r_vi + kpr * ov
    xe = Zseries + x_vi + kpx * ov
    F_i = 2 * i * (re * re + xe * xe) + (on ? i * i * (2 * re * kpr + 2 * xe * kpx) : zero(i))
    didd = iszero(F_i) ? 1 / hypot(r_vi, Zseries + x_vi) + zero(i) : 2 * d / F_i
    return ForwardDiff.Dual{T}(i, didd * ForwardDiff.partials(du_mag))
end

# PHPS's predictive impedance-divider virtual impedance (the only fault-current limiter of
# the VI-only VSM and VOC builds): the converter is the EMF E behind
# Z_tot = r_eff + j(Zseries + x_eff), with r_eff = r_vi + kpr max(0, |i| - i_lim) and
# x_eff = x_vi + kpx max(0, |i| - i_lim); |i_pred| from `vi_bisection`. Returns the
# commanded EMF u_out = E - (r_eff + j x_eff) i_pred. One recorded branch site (VI engaged).
@inline function _vi_divider(rec::ModeRecorder, E_Re, E_Im, Vdm, Vqm, Zseries, i_lim, r_vi,
                             x_vi, kpr, kpx)
    du_Re = E_Re - Vdm
    du_Im = E_Im - Vqm
    du_mag = sqrt(du_Re * du_Re + du_Im * du_Im)
    ipred_mag = vi_bisection(du_mag, Zseries, i_lim, r_vi, x_vi, kpx, kpr)
    ov_vi = gt(rec, ipred_mag, i_lim) ? ipred_mag - i_lim : zero(ipred_mag)
    r_eff = r_vi + kpr * ov_vi
    x_eff = x_vi + kpx * ov_vi
    xt_vi = Zseries + x_eff
    den_vi = r_eff * r_eff + xt_vi * xt_vi
    ipred_Re = (du_Re * r_eff + du_Im * xt_vi) / den_vi
    ipred_Im = (du_Im * r_eff - du_Re * xt_vi) / den_vi
    uout_Re = E_Re - (r_eff * ipred_Re - x_eff * ipred_Im)
    uout_Im = E_Im - (r_eff * ipred_Im + x_eff * ipred_Re)
    return uout_Re, uout_Im
end

# PHPS guards `(T > 1e-5) ? T : 1e-5` on measurement time constants (parameters only)
@inline _tguard(T::Float64, floor::Float64) = T > floor ? T : floor

# The reference-machine speed input of the pf_frame validation mode: read only when
# pf_frame > 0.5, and an unwired input (0) is guarded back to 1 (one recorded site).
@inline function _wref(rec::ModeRecorder, pf_frame::Float64, w)
    w_ref = pf_frame > 0.5 ? w : one(w)
    return lt(rec, w_ref, 0.5) ? one(w_ref) : w_ref
end
