# Scenario loader: the PHPS scenario JSON (system, solver, events, output, plots).
#
# Event defaults follow the PHPS DAE pipeline (dirac/dae_compiler.py `_parse_bus_faults`
# and runner.py for LineFault): a BusFault without `x` is a bolted fault with r = 0,
# x = 1e-5; a window without `t_duration` or `t_end` lasts 0.1 s.

abstract type AbstractEvent end

"""Three-phase shunt fault `1/(r + jx)` at a bus on `[t_start, t_end)`."""
struct BusFault <: AbstractEvent
    bus::Int
    r::Float64
    x::Float64
    t_start::Float64
    t_end::Float64
end

"""Fault at fraction `distance` along a line (from `bus1`), through `1/(r + jx)`."""
struct LineFault <: AbstractEvent
    line_idx::String
    distance::Float64
    r::Float64
    x::Float64
    t_start::Float64
    t_end::Float64
end

"""An event type Porthos reads but does not interpret yet (kept verbatim)."""
struct OtherEvent <: AbstractEvent
    type::String
    raw::JSON3.Object
end
Base.:(==)(a::OtherEvent, b::OtherEvent) = a.type == b.type && json_identical(a.raw, b.raw)

struct SolverSettings
    method::Symbol                  # :ida, :bdf1 or :rk4
    dt::Float64
    duration::Float64
    log_dt::Union{Nothing,Float64}
    rtol::Union{Nothing,Float64}
    atol::Union{Nothing,Float64}
end

"""
    Scenario

A loaded scenario file. `system_path` is resolved relative to the scenario file.
"""
struct Scenario
    path::String
    raw::JSON3.Object
    description::String
    system_path::String
    solver::SolverSettings
    events::Vector{AbstractEvent}
    output_dir::Union{Nothing,String}
end

function Base.:(==)(a::Scenario, b::Scenario)
    json_identical(a.raw, b.raw) && a.description == b.description &&
        a.system_path == b.system_path && a.solver == b.solver && a.events == b.events &&
        a.output_dir == b.output_dir
end

Base.show(io::IO, s::Scenario) = print(io, "Scenario(\"", basename(s.path), "\": ",
                                       s.solver.method, ", ", length(s.events), " events)")

function _window(e)
    t0 = param_value(e[:t_start])
    t1 = haskey(e, :t_duration) ? t0 + param_value(e[:t_duration]) :
         haskey(e, :t_end) ? param_value(e[:t_end]) : t0 + 0.1
    return t0, t1
end

function _event(e)
    type = string(e[:type])
    if type == "BusFault"
        t0, t1 = _window(e)
        return BusFault(Int(e[:bus]), _num(e, :r, 0.0), _num(e, :x, 1e-5), t0, t1)
    elseif type == "LineFault"
        t0, t1 = _window(e)
        return LineFault(string(e[:line_idx]), _num(e, :distance, 0.5), _num(e, :r, 0.0),
                         _num(e, :x, 0.001), t0, t1)
    end
    return OtherEvent(type, e)
end

_optnum(o, k::Symbol) = haskey(o, k) && o[k] !== nothing ? param_value(o[k]) : nothing

"""
    load_scenario(path; validate = true) -> Scenario
"""
function load_scenario(path::AbstractString; validate::Bool = true)
    raw = read_json(path)
    raw isa JSON3.Object ||
        throw(SchemaError(String(path), :scenario, "top level is not an object"))
    validate && validate_json(:scenario, raw; file = path)
    return Scenario(String(path), raw)
end

function Scenario(path::String, raw::JSON3.Object)
    s = raw[:solver]
    solver = SolverSettings(Symbol(s[:method]), param_value(s[:dt]), param_value(s[:duration]),
                            _optnum(s, :log_dt), _optnum(s, :rtol), _optnum(s, :atol))
    events = AbstractEvent[_event(e) for e in _table(raw, :events)]
    out = haskey(raw, :output) && haskey(raw[:output], :directory) ?
          string(raw[:output][:directory]) : nothing
    system_path = normpath(joinpath(dirname(path), string(raw[:system])))
    return Scenario(path, raw, _str(raw, :description, ""), system_path, solver, events, out)
end

"""
    write_scenario(path, scenario)

Write the scenario document. Reading it back gives an identical [`Scenario`](@ref) (the
`system` entry is written as it was, so keep the file next to its system file).
"""
write_scenario(path::AbstractString, sc::Scenario) = write_json(path, sc.raw)

"""
    load_json_input(path) -> Case or Scenario

Load a system or scenario file, telling them apart by their top-level keys.
"""
function load_json_input(path::AbstractString; validate::Bool = true)
    raw = read_json(path)
    raw isa JSON3.Object || throw(ArgumentError("$path: top level is not a JSON object"))
    if haskey(raw, :Bus)
        validate && validate_json(:system, raw; file = path)
        return Case(String(path), raw)
    elseif haskey(raw, :system)
        validate && validate_json(:scenario, raw; file = path)
        return Scenario(String(path), raw)
    end
    throw(ArgumentError("$path is neither a system file (Bus) nor a scenario (system)"))
end
