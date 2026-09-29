# Component parameter processing done by PHPS at load time.
#
# PHPS (json_compat.instantiate_components) converts synchronous-machine parameters from the
# machine base Sn to the system base mva_base when a new-format case is loaded, and the model
# constructors add defaults. The numbers every later stage uses are the processed ones, so
# Porthos reproduces the processing exactly, including its heuristics. Source:
# phps/src/json_compat.py `_normalise_genrou_params` (PHPS commit ba11ea1).

"""
Component model types Porthos implements: the final PHTRUE models of PHPS plus
`COMPLEXLOAD`. The older models in PHPS `components/retired/` are not ported for now; they
stay in PHPS_Opt as reference. Cases that use them still load (`load_case` does not depend
on model types), but every model-level function rejects them.
"""
const MODEL_TYPES = Set(["GENROU_PHTRUE", "GENSAL_PHTRUE", "IEEET1_PHTRUE", "IEEEG1_PHTRUE",
                         "IEEEG3_PHTRUE", "COMPLEXLOAD", "GFL_PHTRUE", "GFL_ZIF_PHTRUE",
                         "GFM_VSM_PHTRUE", "GFM_DROOP_PHTRUE", "GFM_VOC_PHTRUE"])

"""A component whose model type Porthos does not implement."""
struct UnsupportedModelError <: Exception
    component::String
    type::String
end
Base.showerror(io::IO, e::UnsupportedModelError) =
    print(io, "UnsupportedModelError: component ", e.component, " has model type ", e.type,
          ", which Porthos does not implement (supported: ",
          join(sort!(collect(MODEL_TYPES)), ", "), ")")

"""
    check_model_type(spec)

Throw [`UnsupportedModelError`](@ref) unless `spec.type` is in [`MODEL_TYPES`](@ref).
"""
check_model_type(spec::ComponentSpec) =
    spec.type in MODEL_TYPES ? nothing : throw(UnsupportedModelError(spec.name, spec.type))

"""Machine types whose parameters PHPS converts from machine base to system base."""
const MACHINE_BASE_TYPES = Set(["GENROU_PHTRUE", "GENSAL_PHTRUE"])

const _GENROU_ALIASES = ("xd1" => "xd_prime", "xq1" => "xq_prime",
                         "xd2" => "xd_double_prime", "xq2" => "xq_double_prime",
                         "Td10" => "Td0_prime", "Tq10" => "Tq0_prime",
                         "Td20" => "Td0_double_prime", "Tq20" => "Tq0_double_prime")

"""
    ParamDict

Processed parameters of one component: numbers as `Float64`, everything else (names,
flags, strings that are not expressions) as read.
"""
const ParamDict = Dict{String,Any}

_pv(x) = x isa Real && !(x isa Bool) ? Float64(x) : x
_fget(d::ParamDict, k, default) = haskey(d, k) ? param_value(d[k]) : Float64(default)

function _kfd_scale!(out::ParamDict, xd_val, xl_val, xd1_val)
    Xad = xd_val - xl_val
    if Xad > 1e-6 && (xd1_val - xl_val) > 1e-6
        denom = Xad - (xd1_val - xl_val)
        if abs(denom) > 1e-6
            Xfl = (Xad * (xd1_val - xl_val)) / denom
            out["Kfd_scale"] = (Xad + Xfl) / Xad
        else
            out["Kfd_scale"] = 1.0
        end
    else
        out["Kfd_scale"] = 1.0
    end
    return out
end

_fn_string(fn::Float64) = "2.0 * M_PI * " * _py_float_str(fn)

# Python's str(float) for the values that occur here (e.g. 60.0 -> "60.0", 50.0 -> "50.0").
function _py_float_str(x::Float64)
    isinteger(x) && abs(x) < 1e16 && return string(Int(x)) * ".0"
    return repr(x)
end

"""
    normalise_machine_params(params, mva_base, fn) -> ParamDict

Port of PHPS `_normalise_genrou_params`: rename the short aliases, scale impedances by
`mva_base/Sn` and inertia and damping by `Sn/mva_base`, apply the default `D = 2 Sn/mva_base`
when `D` is absent, repair `xd'' <= xl`, and compute `Kfd_scale`. Parameters that are
already on the system base (detected as PHPS does) are not scaled again.
"""
function normalise_machine_params(params::ParamDict, mva_base::Float64, fn::Float64)
    out = copy(params)
    _Sn = _fget(params, "Sn", mva_base)
    _Z = mva_base / _Sn
    already = false
    if abs(_Z - 1.0) > 0.01 && haskey(params, "xd1") && haskey(params, "xd_prime")
        _xd1 = param_value(params["xd1"])
        _xdp = param_value(params["xd_prime"])
        if _xd1 > 0 && abs(_xdp - _xd1 * _Z) / abs(_xdp) < 0.01
            already = true
        end
    end

    if already
        for (k_old, k_new) in _GENROU_ALIASES
            if haskey(params, k_old) && !haskey(params, k_new)
                out[k_new] = params[k_old]
            end
        end
        out["omega_b"] = _fn_string(fn)
        if haskey(params, "M") && !haskey(params, "H")
            out["H"] = param_value(params["M"]) / 2.0
        end
        xd_pp = _fget(out, "xd_double_prime", 0.2)
        xd_p = _fget(out, "xd_prime", 0.3)
        xq_pp = _fget(out, "xq_double_prime", 0.2)
        xq_p = _fget(out, "xq_prime", 0.55)
        xl = _fget(out, "xl", 0.15)
        xd_pp <= xl && (out["xd_double_prime"] = xl + 0.5 * (xd_p - xl))
        xq_pp <= xl && (out["xq_double_prime"] = xl + 0.5 * (xq_p - xl))
        return _kfd_scale!(out, _fget(out, "xd", 1.8), _fget(out, "xl", 0.15),
                           _fget(out, "xd_prime", 0.3))
    end

    for (k_old, k_new) in _GENROU_ALIASES
        haskey(params, k_old) && (out[k_new] = params[k_old])
    end
    haskey(params, "M") && (out["H"] = param_value(params["M"]) / 2.0)
    out["omega_b"] = _fn_string(fn)

    Sn = _fget(out, "Sn", mva_base)
    Z_scale = mva_base / Sn
    M_scale = Sn / mva_base
    for key in ("ra", "xd", "xq", "xd_prime", "xq_prime", "xd_double_prime",
                "xq_double_prime", "xl")
        haskey(out, key) && (out[key] = param_value(out[key]) * Z_scale)
    end
    if haskey(out, "M")
        out["M"] = param_value(out["M"]) * M_scale
        out["H"] = out["M"] / 2.0
    elseif haskey(out, "H")
        out["H"] = param_value(out["H"]) * M_scale
    end
    # An explicit D is honoured (including 0); the 2 pu machine-base floor applies only
    # when the case declares no D.
    out["D"] = haskey(params, "D") ? param_value(params["D"]) * M_scale : 2.0 * M_scale

    xd_pp = haskey(out, "xd_double_prime") ? param_value(out["xd_double_prime"]) :
            _fget(out, "xd2", 0.2)
    xq_pp = haskey(out, "xq_double_prime") ? param_value(out["xq_double_prime"]) :
            _fget(out, "xq2", 0.2)
    xd_p = haskey(out, "xd_prime") ? param_value(out["xd_prime"]) : _fget(out, "xd1", 0.3)
    xq_p = haskey(out, "xq_prime") ? param_value(out["xq_prime"]) : _fget(out, "xq1", 0.55)
    xl = _fget(out, "xl", 0.15)
    xd_pp <= xl && (out["xd_double_prime"] = xl + 0.5 * (xd_p - xl))
    xq_pp <= xl && (out["xq_double_prime"] = xl + 0.5 * (xq_p - xl))

    xd1_val = haskey(out, "xd_prime") ? param_value(out["xd_prime"]) : _fget(out, "xd1", 0.3)
    return _kfd_scale!(out, _fget(out, "xd", 1.8), _fget(out, "xl", 0.15), xd1_val)
end

"""
    component_params(case, spec) -> ParamDict

The parameters PHPS uses for a component: the case values, normalised to the system base
for synchronous machines (unless `_params_normalized` is set), plus the defaults and derived
values its constructor adds ([`type_defaults!`](@ref)). Parameters set by initialisation
(reservoir references, load `Vini`, `Pref`/`Vref`) keep their defaults here.
"""
function component_params(case::Case, spec::ComponentSpec)
    check_model_type(spec)
    p = ParamDict(string(k) => _pv(v) for (k, v) in spec.params)
    if spec.type in MACHINE_BASE_TYPES && !(get(p, "_params_normalized", false) === true)
        p = normalise_machine_params(p, case.config.mva_base, case.config.fn)
        p["_params_normalized"] = true
    elseif spec.type in ("GFM_VSM_PHTRUE", "GFM_DROOP_PHTRUE", "GFM_VOC_PHTRUE")
        # gfm_vsm_phtrue.py / gfm_voc_phtrue.py constructors (droop inherits from VSM).
        get!(p, "ra", 0.0)
        get!(p, "xd_double_prime", haskey(p, "Zseries") ? p["Zseries"] : 0.10)
    elseif spec.type in ("GFL_PHTRUE", "GFL_ZIF_PHTRUE")
        get!(p, "ra", 0.0)
    end
    return type_defaults!(Val(Symbol(spec.type)), p)
end

"""
    type_defaults!(::Val{type}, p::ParamDict) -> p

The defaults and derived parameters a PHPS model constructor adds (e.g. reservoir
capacities, saturation coefficients). Each model file adds its method.
"""
type_defaults!(::Val, p::ParamDict) = p
