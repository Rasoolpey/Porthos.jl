# Component interface (roadmap 2.2).
#
# Every model is a subtype of `AbstractComponent` holding its name, its typed parameters
# `p` and the processed parameter dictionary it was built from. The model functions are
# generic in the number type of the states and inputs (Float64, ForwardDiff.Dual, and the
# interval types of P11), take the parameters explicitly, and do not allocate.
#
# A model implements
#
#   state_names(c), input_names(c), output_names(c)
#   _step!(dx, y, c, x, u, p, rec)   state derivatives, and (y !== nothing) the outputs PHPS
#                                    writes during its step
#   _outputs!(y, c, x, u, p, rec)    the output kernel (what the network and the other
#                                    components read)
#   hamiltonian(c, x, p), grad_hamiltonian!(g, c, x, p)
#
# and, when it is attached to a bus, `injection` and `norton_admittance`. Everything else
# (`rhs!`, `outputs!`, `step_outputs!`, `modes`, `contract`, ...) is derived here.

"""Supertype of every component model."""
abstract type AbstractComponent end

"""Model type name as used in case files (e.g. `"GENROU_PHTRUE"`)."""
function model_type end
"""State names, in PHPS order."""
function state_names end
"""Input port names, in PHPS order."""
function input_names end
"""Output port names, in PHPS order."""
function output_names end

nstates(c::AbstractComponent) = length(state_names(c))
ninputs(c::AbstractComponent) = length(input_names(c))
noutputs(c::AbstractComponent) = length(output_names(c))

"""Typed parameters of a component."""
params(c::AbstractComponent) = c.p
"""The processed parameter dictionary the component was built from."""
param_dict(c::AbstractComponent) = c.params
"""Instance name."""
name(c::AbstractComponent) = c.name

"""
    ports(c) -> (in = [...], out = [...])
"""
ports(c::AbstractComponent) = (in = input_names(c), out = output_names(c))

"""
    rhs!(dx, c, x, u[, p]) -> dx

State derivatives. Generic in the element types of `x` and `u`; does not allocate.
"""
rhs!(dx, c::AbstractComponent, x, u, p = params(c)) = (_step!(dx, nothing, c, x, u, p, NoModes()); dx)

"""
    outputs!(y, c, x, u[, p]) -> y

The output kernel: the port outputs the network (KCL) and the other components read.
"""
outputs!(y, c::AbstractComponent, x, u, p = params(c)) = (_outputs!(y, c, x, u, p, NoModes()); y)

"""
    step_outputs!(y, c, x, u[, p]) -> y

The outputs after PHPS's step kernel has run: `outputs!` first, then the entries the step
overwrites (terminal power, dq currents, ...). These are what PHPS logs.
"""
function step_outputs!(y, c::AbstractComponent, x, u, p = params(c))
    _outputs!(y, c, x, u, p, NoModes())
    dx = similar(x, promote_type(eltype(x), eltype(u)), nstates(c))
    _step!(dx, y, c, x, u, p, NoModes())
    return y
end

"""
    modes(c, x, u[, p]; kernel = :step) -> Vector{Bool}

The state-dependent branch decisions of the step (or `:out`) kernel, in evaluation order.
"""
function modes(c::AbstractComponent, x, u, p = params(c); kernel::Symbol = :step)
    rec = ModeLog()
    T = promote_type(eltype(x), eltype(u))
    if kernel === :step
        _step!(zeros(T, nstates(c)), nothing, c, x, u, p, rec)
    elseif kernel === :out
        _outputs!(zeros(T, noutputs(c)), c, x, u, p, rec)
    else
        throw(ArgumentError("kernel must be :step or :out"))
    end
    return rec.decisions
end

"""
    hamiltonian(c, x[, p])

Declared storage `H(x)`.
"""
function hamiltonian end
hamiltonian(c::AbstractComponent, x) = hamiltonian(c, x, params(c))

"""
    grad_hamiltonian!(g, c, x[, p]) -> g
"""
function grad_hamiltonian! end
grad_hamiltonian!(g, c::AbstractComponent, x) = grad_hamiltonian!(g, c, x, params(c))

"""
    grad_hamiltonian(c, x[, p]) -> Vector
"""
grad_hamiltonian(c::AbstractComponent, x, p = params(c)) =
    grad_hamiltonian!(zeros(eltype(x), nstates(c)), c, x, p)

"""Bus id of a network-attached component, `nothing` otherwise."""
bus(c::AbstractComponent) = nothing

"""
    component_role(c) -> Symbol

`:generator`, `:exciter`, `:governor` or `:load`, as PHPS's `component_role` (it decides the
wiring refresh after initialisation and the centre-of-inertia members).
"""
function component_role end

"""
    injection(c, x, V[, p]) -> (I_re, I_im)

Current the component injects into its bus, in the network (RI) frame, for bus voltage
`V = (Vd, Vq)`. For a Norton source this is the Norton current; its admittance is
[`norton_admittance`](@ref) and sits in the Y-bus.
"""
function injection end

"""
    norton_admittance(c) -> Union{Nothing, ComplexF64}

Admittance the component adds to the Y-bus diagonal at its bus, or `nothing`.
"""
norton_admittance(c::AbstractComponent) = nothing

"""
    contract(c) -> ContractEntry

The component's port-contract entry (contracts/model_port_contracts.json).
"""
contract(c::AbstractComponent) = contract(default_contracts(), model_type(c))

const _CONTRACTS = Ref{Union{Nothing,ContractSet}}(nothing)
const _CONTRACTS_LOCK = ReentrantLock()

"""The contract set in `contracts/`, loaded once."""
function default_contracts()
    lock(_CONTRACTS_LOCK) do
        _CONTRACTS[] === nothing && (_CONTRACTS[] = load_contracts())
        return _CONTRACTS[]
    end
end

Base.show(io::IO, c::AbstractComponent) = print(io, model_type(c), "(\"", name(c), "\")")

# -- construction ----------------------------------------------------------------------

"""Component constructors by model type: `(name, params::ParamDict) -> component`."""
const COMPONENT_CONSTRUCTORS = Dict{String,Any}()

"""
    build_component(case, spec) -> AbstractComponent
    build_component(type, name, params::ParamDict) -> AbstractComponent

Build a model from a case entry, with the parameters PHPS would use
([`component_params`](@ref)).
"""
build_component(case::Case, spec::ComponentSpec) =
    build_component(spec.type, spec.name, component_params(case, spec))

function build_component(type::AbstractString, name::AbstractString, params::ParamDict)
    haskey(COMPONENT_CONSTRUCTORS, type) ||
        throw(UnsupportedModelError(String(name), String(type)))
    return COMPONENT_CONSTRUCTORS[type](String(name), params)
end

"""
    with_params(c, overrides) -> component

A copy of `c` with some parameters replaced (e.g. the references set by initialisation).
"""
with_params(c::AbstractComponent, overrides::AbstractDict) =
    build_component(model_type(c), name(c), merge(param_dict(c), ParamDict(string(k) => v
                                                                          for (k, v) in overrides)))

# Parameter access for constructors: the processed value, or an error naming the component.
function _p(d::ParamDict, key::String, cname::String)
    haskey(d, key) || throw(ArgumentError("component $cname: missing parameter $key"))
    return param_value(d[key])
end
_p(d::ParamDict, key::String, cname::String, default) =
    haskey(d, key) ? param_value(d[key]) : Float64(default)

# Round as PHPS does when it writes a constant into C++ with printf("%.12e") / ("%.6f").
_c12e(x::Float64) = parse(Float64, Printf.@sprintf("%.12e", x))
_c6f(x::Float64) = parse(Float64, Printf.@sprintf("%.6f", x))
