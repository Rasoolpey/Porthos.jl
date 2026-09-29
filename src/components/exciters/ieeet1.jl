# IEEET1_PHTRUE: IEEE Type-1 (DC1A) exciter with a field-supply reservoir account.
#
# Port of the equations PHPS runs for IEEET1_PHTRUE (phps/src/components/exciters/
# ieeet1_phtrue.py at ba11ea1): the four control states of its retired parent Ieeet1PHS,
# unchanged, plus the one-way field-supply tank drained by Efd*i_fd (corrected
# constant-power law). Expressions keep PHPS's order of operations.
#
# States  [xr, xa, xe, xf, x_field]
# Inputs  [Vterm, Vref, upss, i_fd]
# Outputs [Efd]
#
# Branch sites (step): regulator non-windup at VRMAX / VRMIN, the clamp of the regulator
# output seen by the exciter, the reservoir pressure guard. The TR, TF, TA, TE tests are
# parameter conditions, not branches.

struct Ieeet1Params
    TR::Float64
    KA::Float64
    TA::Float64
    KE::Float64
    TE::Float64
    KF::Float64
    TF::Float64
    VRMAX::Float64
    VRMIN::Float64
    SAT_A::Float64
    SAT_B::Float64
    C_FIELD::Float64
    F_STAR::Float64
    PFD_REF::Float64
end

struct IEEET1_PHTRUE <: AbstractComponent
    name::String
    p::Ieeet1Params
    params::ParamDict
end

const _IEEET1_STATES = ["xr", "xa", "xe", "xf", "x_field"]
const _IEEET1_INPUTS = ["Vterm", "Vref", "upss", "i_fd"]
const _IEEET1_OUTPUTS = ["Efd"]

model_type(::IEEET1_PHTRUE) = "IEEET1_PHTRUE"
state_names(::IEEET1_PHTRUE) = _IEEET1_STATES
input_names(::IEEET1_PHTRUE) = _IEEET1_INPUTS
output_names(::IEEET1_PHTRUE) = _IEEET1_OUTPUTS

"""
    ieeet1_sat_coeffs(E1, SE1, E2, SE2) -> (A, B)

PowerFactory avr_IEEET1 exponential saturation `Se(E) = A exp(B E)` through
`(E1, SE1)` and `(E2, SE2)`; `(0, 0)` when saturation is disabled.
"""
function ieeet1_sat_coeffs(E1, SE1, E2, SE2)
    (E1 <= 0 || E2 <= 0 || SE1 <= 0 || SE2 <= 0 || abs(E2 - E1) < 1e-10) && return 0.0, 0.0
    B = log(SE2 / SE1) / (E2 - E1)
    return SE2 / exp(B * E2), B
end

function type_defaults!(::Val{:IEEET1_PHTRUE}, p::ParamDict)
    get!(p, "C_FIELD", 1000.0)
    get!(p, "F_STAR", 1.0)
    get!(p, "PFD_REF", 0.0)
    if !haskey(p, "SAT_A") || !haskey(p, "SAT_B")
        # PHPS Ieeet1PHS.__init__: fit from E1/SE1/E2/SE2 (and drop them)
        g(k, alt = nothing) = haskey(p, k) ? param_value(pop!(p, k)) :
                              alt !== nothing && haskey(p, alt) ? param_value(pop!(p, alt)) : 0.0
        E1 = g("E1"); SE1 = g("SE1", "Se1"); E2 = g("E2"); SE2 = g("SE2", "Se2")
        p["SAT_A"], p["SAT_B"] = ieeet1_sat_coeffs(E1, SE1, E2, SE2)
    end
    return p
end

function IEEET1_PHTRUE(name::String, d::ParamDict)
    f(k) = _p(d, k, name)
    p = Ieeet1Params(f("TR"), f("KA"), f("TA"), f("KE"), f("TE"), f("KF"), f("TF"),
                     f("VRMAX"), f("VRMIN"), f("SAT_A"), f("SAT_B"), f("C_FIELD"),
                     f("F_STAR"), f("PFD_REF"))
    return IEEET1_PHTRUE(name, p, d)
end
COMPONENT_CONSTRUCTORS["IEEET1_PHTRUE"] = IEEET1_PHTRUE

function _outputs!(y, c::IEEET1_PHTRUE, x, u, p::Ieeet1Params, rec::ModeRecorder)
    y[1] = x[3]                               # Efd = xe
    return y
end

function _step!(dx, y, c::IEEET1_PHTRUE, x, u, p::Ieeet1Params, rec::ModeRecorder)
    T = promote_type(eltype(x), eltype(u))
    xr, xa, xe, xf = x[1], x[2], x[3], x[4]
    Vterm, Vref, upss = u[1], u[2], u[3]

    # 1. transducer; TR = 0 is an exact passthrough (parameter condition)
    Vc = p.TR > 1.0e-9 ? T(xr) : T(Vterm)
    dx[1] = p.TR > 1.0e-9 ? (Vterm - xr) / p.TR : zero(T)
    # 2. rate-feedback washout
    TF_eff = p.TF > 1e-4 ? p.TF : 1e-4
    Vfb = p.KF * (xe - xf) / TF_eff
    # 3. voltage error (+ PSS)
    Verr = Vref - Vc - Vfb + upss
    # 4. regulator Ka/Ta, non-windup at VRMAX/VRMIN
    TA_eff = p.TA > 1e-4 ? p.TA : 1e-4
    dxa = (p.KA * Verr - xa) / TA_eff
    dxa = nonwindup(rec, xa, dxa, p.VRMIN, p.VRMAX)
    dx[2] = dxa
    # 5. exciter Ke/Te with exponential saturation; sees the limited regulator output
    Se = p.SAT_A * exp(p.SAT_B * xe)
    TE_eff = p.TE > 1e-4 ? p.TE : 1e-4
    xa_lim = clamp_mode(rec, xa, p.VRMIN, p.VRMAX)
    dx[3] = (xa_lim - (p.KE + Se) * xe) / TE_eff
    # 6. rate-feedback state
    dx[4] = (xe - xf) / TF_eff

    # field-supply reservoir (one-way; does not feed back)
    x_field = x[5]
    i_fd = u[4]
    _pfld = guard_min(rec, x_field / p.C_FIELD, 1e-6)
    dx[5] = (p.PFD_REF - xe * i_fd) / _pfld
    return dx
end

hamiltonian(c::IEEET1_PHTRUE, x, p::Ieeet1Params) = 0.5 * (x[5] * x[5]) / p.C_FIELD

function grad_hamiltonian!(g, c::IEEET1_PHTRUE, x, p::Ieeet1Params)
    fill!(g, 0.0)
    g[5] = x[5] / p.C_FIELD
    return g
end
