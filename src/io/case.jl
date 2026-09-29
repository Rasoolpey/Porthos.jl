# System case loader: the PHPS system JSON (config, Bus, PQ, PV, Slack, Line, Shunt,
# components, connections), unchanged format.
#
# Typed tables carry the values PHPS computes with, including PHPS's defaults for missing
# keys (taken from ybus.py, powerflow.py and system_graph.py). The parsed document is kept
# as `raw`, so provenance keys (`_pf`, `rating_src`, ...) and everything not yet interpreted
# survive a write.

"""Case-wide settings from the `config` block."""
struct CaseConfig
    mva_base::Float64
    fn::Float64
    omega_b::Float64          # evaluated from its expression, e.g. "2.0 * M_PI * 60.0"
    skip_pf_solve::Bool       # PHPS: use Bus v0/a0 instead of solving the power flow
    kpf::Float64              # global frequency dependence of loads (fallback)
    kqf::Float64
end

struct BusData
    idx::Int
    name::String
    Vn::Float64
    v0::Float64
    a0::Float64
end

"""A `PQ` entry: constant power consumption (positive = drawn from the network)."""
struct LoadData
    idx::String
    bus::Int
    p0::Float64
    q0::Float64
end

"""A `PV` or `Slack` entry."""
struct SourceData
    idx::String
    bus::Int
    p0::Float64
    q0::Float64
    v0::Float64
    a0::Float64               # Slack angle; 0 for PV
    pmax::Float64
    pmin::Float64
    qmax::Float64
    qmin::Float64
end

"""A `Line` entry: pi line or two-winding transformer (tap and phase shift on `bus1`)."""
struct LineData
    idx::String
    bus1::Int
    bus2::Int
    r::Float64
    x::Float64
    b::Float64                # total line charging
    tap::Float64
    phi::Float64              # phase shift [rad]
end

struct ShuntData
    idx::String
    bus::Int
    g::Float64
    b::Float64
end

"""A dynamic component: its name, model type and parameters as written in the case."""
struct ComponentSpec
    name::String
    type::String
    params::JSON3.Object
end

Base.:(==)(a::ComponentSpec, b::ComponentSpec) =
    a.name == b.name && a.type == b.type && json_identical(a.params, b.params)

"""A signal wire `from -> to` from the `connections` list (`"COMP.port"` strings)."""
struct Wire
    from::String
    to::String
end

"""
    Case

A loaded system file. Tables keep file order; `raw` is the parsed document.
"""
struct Case
    path::String
    raw::JSON3.Object
    description::String
    config::CaseConfig
    buses::Vector{BusData}
    pq::Vector{LoadData}
    pv::Vector{SourceData}
    slack::Vector{SourceData}
    lines::Vector{LineData}
    shunts::Vector{ShuntData}
    components::Vector{ComponentSpec}
    connections::Vector{Wire}
end

function Base.:(==)(a::Case, b::Case)
    for f in fieldnames(Case)
        f === :path && continue
        if f === :raw
            json_identical(a.raw, b.raw) || return false
        else
            getfield(a, f) == getfield(b, f) || return false
        end
    end
    return true
end

Base.show(io::IO, c::Case) = print(io, "Case(\"", basename(c.path), "\": ", length(c.buses),
                                   " buses, ", length(c.lines), " lines, ",
                                   length(c.components), " components)")

# -- field readers ------------------------------------------------------------------------

_num(o, k::Symbol, default) = haskey(o, k) ? param_value(o[k]) : Float64(default)
_str(o, k::Symbol, default) = haskey(o, k) ? string(o[k]) : String(default)
_int(o, k::Symbol) = Int(o[k])

function _source(o, is_slack::Bool)
    SourceData(_str(o, :idx, string(o[:bus])), _int(o, :bus),
               _num(o, :p0, 0.0), _num(o, :q0, 0.0), _num(o, :v0, 1.0),
               is_slack ? _num(o, :a0, 0.0) : 0.0,
               _num(o, :pmax, Inf), _num(o, :pmin, -Inf),
               _num(o, :qmax, Inf), _num(o, :qmin, -Inf))
end

_table(raw, key::Symbol) = haskey(raw, key) && raw[key] !== nothing ? raw[key] : ()

"""
    load_case(path; validate = true) -> Case

Read a system file. With `validate`, the document is first checked against
`cases/schema/system.schema.json`.
"""
function load_case(path::AbstractString; validate::Bool = true)
    raw = read_json(path)
    raw isa JSON3.Object || throw(SchemaError(String(path), :system, "top level is not an object"))
    validate && validate_json(:system, raw; file = path)
    return Case(String(path), raw)
end

function Case(path::String, raw::JSON3.Object)
    cfg = haskey(raw, :config) ? raw[:config] : JSON3.Object()
    fn = _num(cfg, :fn, 60.0)
    omega_b = haskey(cfg, :omega_b) ? param_value(cfg[:omega_b]) : 2.0 * π * fn
    config = CaseConfig(_num(cfg, :mva_base, 100.0), fn, omega_b,
                        Bool(get(cfg, :skip_pf_solve, false)),
                        _num(cfg, :kpf, 0.0), _num(cfg, :kqf, 0.0))

    buses = [BusData(_int(b, :idx), _str(b, :name, string(b[:idx])), _num(b, :Vn, 1.0),
                     _num(b, :v0, 1.0), _num(b, :a0, 0.0)) for b in _table(raw, :Bus)]
    pq = [LoadData(_str(p, :idx, string(p[:bus])), _int(p, :bus), _num(p, :p0, 0.0),
                   _num(p, :q0, 0.0)) for p in _table(raw, :PQ)]
    pv = [_source(p, false) for p in _table(raw, :PV)]
    slack = [_source(s, true) for s in _table(raw, :Slack)]
    # Defaults as in PHPS YBusBuilder (x = 0.001, not the 0.01 of system_graph.LineData).
    lines = [LineData(_str(l, :idx, ""), _int(l, :bus1), _int(l, :bus2), _num(l, :r, 0.0),
                      _num(l, :x, 0.001), _num(l, :b, 0.0), _num(l, :tap, 1.0),
                      _num(l, :phi, 0.0)) for l in _table(raw, :Line)]
    shunts = [ShuntData(_str(s, :idx, ""), _int(s, :bus), _num(s, :g, 0.0), _num(s, :b, 0.0))
              for s in _table(raw, :Shunt)]
    components = [ComponentSpec(string(name), string(c[:type]), c[:params])
                  for (name, c) in raw[:components]]
    connections = [Wire(string(w[:from]), string(w[:to])) for w in raw[:connections]]

    ids = Set{Int}()
    for b in buses
        b.idx in ids && throw(JSONPathError(path, "Bus", "duplicate bus idx $(b.idx)"))
        push!(ids, b.idx)
    end

    return Case(path, raw, _str(raw, :description, ""), config, buses, pq, pv, slack, lines,
                shunts, components, connections)
end

"""
    write_case(path, case)

Write the case document. Reading it back gives an identical [`Case`](@ref).
"""
write_case(path::AbstractString, case::Case) = write_json(path, case.raw)

"""
    component(case, name) -> ComponentSpec
"""
function component(case::Case, name::AbstractString)
    i = findfirst(c -> c.name == name, case.components)
    i === nothing && throw(KeyError(name))
    return case.components[i]
end

"""
    param(spec, key[, default]) -> Float64

Numeric parameter of a component (numbers or expression strings).
"""
param(c::ComponentSpec, key::Union{Symbol,AbstractString}) = param_value(c.params[Symbol(key)])
param(c::ComponentSpec, key::Union{Symbol,AbstractString}, default) =
    haskey(c.params, Symbol(key)) ? param_value(c.params[Symbol(key)]) : Float64(default)

hasparam(c::ComponentSpec, key::Union{Symbol,AbstractString}) = haskey(c.params, Symbol(key))
