# Port-contract loader: contracts/model_port_contracts.json (schema 2.1-reservoir-corrected),
# copied unchanged from PHPS. Each model type has a status, storage and domain statements,
# its ports, an optional reservoir, and fail-closed `domain_clauses`. The clauses are read
# here; they are evaluated by the containment audit (P11).

const CONTRACTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "contracts",
                                         "model_port_contracts.json"))
const CONTRACT_SCHEMA_VERSION = "2.1-reservoir-corrected"

"""
    DomainClause

One fail-closed domain clause. `check` is one of `state_bounds`, `parameter_conditions`,
`nonbinding_accounting_state`, `single_smooth_mode` or `bus_voltage_squared_above`. Bounds
and conditions are expression strings over the component parameters.
"""
struct DomainClause
    id::String
    text::String
    check::Symbol
    state::Union{Nothing,String}
    lower::Union{Nothing,String}
    upper::Union{Nothing,String}
    conditions::Vector{String}
    threshold::Union{Nothing,Float64}
    raw::JSON3.Object
end
Base.:(==)(a::DomainClause, b::DomainClause) = json_identical(a.raw, b.raw)

const DOMAIN_CHECKS = (:state_bounds, :parameter_conditions, :nonbinding_accounting_state,
                       :single_smooth_mode, :bus_voltage_squared_above)

struct ContractEntry
    model::String
    status::String
    domain_clauses::Vector{DomainClause}
    raw::JSON3.Object
end
Base.:(==)(a::ContractEntry, b::ContractEntry) = a.model == b.model && json_identical(a.raw, b.raw)

struct ContractSet
    path::String
    schema_version::String
    power_convention::String
    raw::JSON3.Object
    entries::Dict{String,ContractEntry}
end

Base.show(io::IO, c::ContractSet) =
    print(io, "ContractSet(", c.schema_version, ", ", length(c.entries), " models)")

_optstr(o, k::Symbol) = haskey(o, k) && o[k] !== nothing ? string(o[k]) : nothing

function DomainClause(o::JSON3.Object)
    check = Symbol(o[:check])
    check in DOMAIN_CHECKS || throw(ArgumentError("unknown domain clause check '$check'"))
    DomainClause(string(o[:id]), string(get(o, :text, "")), check, _optstr(o, :state),
                 _optstr(o, :lower), _optstr(o, :upper),
                 String[string(c) for c in get(o, :conditions, ())],
                 haskey(o, :threshold) ? param_value(o[:threshold]) : nothing, o)
end

"""
    load_contracts(path = CONTRACTS_PATH) -> ContractSet

Read the port contracts. The schema version must be $(CONTRACT_SCHEMA_VERSION).
"""
function load_contracts(path::AbstractString = CONTRACTS_PATH)
    raw = read_json(path)
    ver = string(raw[:schema_version])
    ver == CONTRACT_SCHEMA_VERSION ||
        throw(ArgumentError("$path: contract schema $ver, expected $CONTRACT_SCHEMA_VERSION"))
    entries = Dict{String,ContractEntry}()
    for (model, c) in raw[:contracts]
        clauses = DomainClause[DomainClause(d) for d in get(c, :domain_clauses, ())]
        entries[string(model)] = ContractEntry(string(model), string(get(c, :status, "")),
                                               clauses, c)
    end
    return ContractSet(String(path), ver, string(get(raw, :power_convention, "")), raw, entries)
end

"""
    contract_key(type) -> String

Contract key of a component type. Only `COMPLEXLOAD` differs (`ComplexLoad`).
"""
contract_key(type::AbstractString) = type == "COMPLEXLOAD" ? "ComplexLoad" : String(type)

"""
    contract(set, type) -> ContractEntry
"""
contract(set::ContractSet, type::AbstractString) = set.entries[contract_key(type)]
