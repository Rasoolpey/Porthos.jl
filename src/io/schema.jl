# JSON Schema validation of case and scenario files (cases/schema/*.schema.json).

const SCHEMA_DIR = normpath(joinpath(@__DIR__, "..", "..", "cases", "schema"))

const _SCHEMAS = Dict{Symbol,JSONSchema.Schema}()
const _SCHEMA_LOCK = ReentrantLock()

function _schema(kind::Symbol)
    lock(_SCHEMA_LOCK) do
        get!(_SCHEMAS, kind) do
            file = joinpath(SCHEMA_DIR, "$(kind).schema.json")
            JSONSchema.Schema(JSON3.read(read(file, String), Dict{String,Any}))
        end
    end
end

"""
    SchemaError

A case or scenario file does not satisfy its JSON schema.
"""
struct SchemaError <: Exception
    file::String
    kind::Symbol
    issue::String
end
Base.showerror(io::IO, e::SchemaError) =
    print(io, "SchemaError: ", e.file, " is not a valid ", e.kind, " file:\n", e.issue)

"""
    validate_json(kind, obj; file = "<memory>")

Validate a parsed JSON value against the `:system` or `:scenario` schema. Throws
[`SchemaError`](@ref) on failure and returns `nothing` otherwise.
"""
function validate_json(kind::Symbol, obj; file::AbstractString = "<memory>")
    # JSONSchema.jl works on plain Dict/Vector trees.
    plain = JSON3.read(JSON3.write(obj), Dict{String,Any})
    issue = JSONSchema.validate(_schema(kind), plain)
    issue === nothing || throw(SchemaError(String(file), kind, sprint(show, issue)))
    return nothing
end
