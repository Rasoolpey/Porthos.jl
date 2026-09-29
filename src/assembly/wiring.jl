# Wiring: where each component input comes from.
#
# Port of PHPS SystemCompiler._resolve_wiring / _wire_src_to_placeholder and the post-
# initialisation refresh of DiracCompiler._refresh_control_params (PHPS ba11ea1):
#
#   CONST:<v>            a constant
#   PARAM:<comp>.<key>   a parameter value, read when the wiring is built (0 if missing)
#   BUS_<id>.Vd|Vq|Vterm the bus voltage (network RI frame); other signals give 0
#   DQ_<gen>.Vd|Vq       the bus voltage in the machine's dq frame
#   <comp>.<port>        an output of another component (0 if the port does not exist)
#
# Wires are applied in file order; a later wire to the same input replaces an earlier one.
# An input without a wire is 0. After initialisation, an exciter's Vref and a governor's
# Pref become their initialised set-points, and a machine's Tm / Efd become Tm0 / Efd0
# unless they are driven by another component's output.

@enum SourceKind::UInt8 SRC_ZERO SRC_CONST SRC_VD SRC_VQ SRC_VTERM SRC_OUTPUT SRC_DQ_VD SRC_DQ_VQ

"""
    InputSource

Where one component input is read from: `kind`, and the constant `value`, the bus or
component `index`, and the output `port` (1-based), as the kind needs.
"""
struct InputSource
    kind::SourceKind
    value::Float64
    index::Int
    port::Int
end
InputSource(kind::SourceKind) = InputSource(kind, 0.0, 0, 0)
const_source(v::Real) = InputSource(SRC_CONST, Float64(v), 0, 0)

function Base.show(io::IO, s::InputSource)
    k = s.kind
    k === SRC_ZERO ? print(io, "0") :
    k === SRC_CONST ? print(io, "const(", s.value, ")") :
    k === SRC_OUTPUT ? print(io, "output(", s.index, ", ", s.port, ")") :
    print(io, lowercase(string(k)[5:end]), "(", s.index, ")")
end

# a numeric literal as C++ reads it (CONST: text is pasted into the generated code)
function _cpp_literal(text::AbstractString)
    v = tryparse(Float64, strip(text))
    return v === nothing ? parse_param_expr(text) : v
end

function _wire_source(src::String, comps, index::Dict{String,Int}, net::Network,
                      pre_params::Vector{ParamDict})
    if startswith(src, "DQ_")
        rest = src[4:end]
        occursin('.', rest) || return nothing
        gen, sig = split(rest, '.'; limit = 2)
        k = get(index, gen, 0)
        k == 0 && throw(ArgumentError("wire source $src: unknown machine $gen"))
        sig in ("Vd", "Vd_dq") && return InputSource(SRC_DQ_VD, 0.0, k, 0)
        sig in ("Vq", "Vq_dq") && return InputSource(SRC_DQ_VQ, 0.0, k, 0)
        return nothing
    elseif startswith(src, "CONST:")
        return const_source(_cpp_literal(src[7:end]))
    elseif startswith(src, "PARAM:")
        rest = src[7:end]
        occursin('.', rest) || return const_source(0.0)
        cname, key = split(rest, '.'; limit = 2)
        k = get(index, cname, 0)
        (k == 0 || !haskey(pre_params[k], key)) && return const_source(0.0)
        return const_source(param_value(pre_params[k][key]))
    elseif startswith(src, "BUS_")
        parts = split(src, '.')
        length(parts) >= 2 || return const_source(0.0)
        bus_id = tryparse(Int, parts[1][5:end])
        bus_id === nothing && return const_source(0.0)
        haskey(net.index, bus_id) || return const_source(0.0)
        i = net.index[bus_id]
        sig = parts[2]
        sig == "Vd" && return InputSource(SRC_VD, 0.0, i, 0)
        sig == "Vq" && return InputSource(SRC_VQ, 0.0, i, 0)
        sig == "Vterm" && return InputSource(SRC_VTERM, 0.0, i, 0)
        return const_source(0.0)
    else
        occursin('.', src) || return const_source(0.0)
        cname = split(src, '.')[1]
        port = split(src, '.'; limit = 2)[2]
        k = get(index, cname, 0)
        k == 0 && return const_source(0.0)
        j = findfirst(==(port), output_names(comps[k]))
        j === nothing && return const_source(0.0)
        return InputSource(SRC_OUTPUT, 0.0, k, j)
    end
end

"""
    resolve_wiring(case, comps, net; pre_params) -> Vector{Vector{InputSource}}

The source of every input of every component (`comps` in case order). `pre_params` are
the parameter dictionaries before initialisation (read by `PARAM:` wires); the refresh
after initialisation uses the components' own (initialised) parameters.
"""
function resolve_wiring(case::Case, comps::Vector{<:AbstractComponent}, net::Network;
                        pre_params::Vector{ParamDict} = [param_dict(c) for c in comps])
    index = Dict(name(c) => k for (k, c) in enumerate(comps))
    src = [fill(InputSource(SRC_ZERO), ninputs(c)) for c in comps]
    for w in case.connections
        occursin('.', w.to) || continue
        dcomp, dport = split(w.to, '.'; limit = 2)
        k = get(index, dcomp, 0)
        k == 0 && continue
        j = findfirst(==(dport), input_names(comps[k]))
        j === nothing && continue
        s = _wire_source(w.from, comps, index, net, pre_params)
        s === nothing || (src[k][j] = s)
    end
    # refresh after initialisation (PHPS DiracCompiler._refresh_control_params)
    for (k, c) in enumerate(comps)
        role = component_role(c)
        pd = param_dict(c)
        ins = input_names(c)
        if role === :exciter && "Vref" in ins && haskey(pd, "Vref")
            src[k][findfirst(==("Vref"), ins)] = const_source(param_value(pd["Vref"]))
        elseif role === :governor && "Pref" in ins && haskey(pd, "Pref")
            src[k][findfirst(==("Pref"), ins)] = const_source(param_value(pd["Pref"]))
        elseif role === :generator
            for (port, key) in (("Tm", "Tm0"), ("Efd", "Efd0"))
                (port in ins && haskey(pd, key)) || continue
                j = findfirst(==(port), ins)
                src[k][j].kind === SRC_OUTPUT && continue
                src[k][j] = const_source(param_value(pd[key]))
            end
        end
    end
    return src
end
