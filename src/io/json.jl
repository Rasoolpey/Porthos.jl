# JSON reading and writing that preserves what PHPS wrote: key order, the integer / float
# distinction, and exact binary64 values (JSON3 parses with correct rounding and writes the
# shortest representation that reads back to the same Float64).

"""
    read_json(path) -> JSON3.Object / JSON3.Array

Parse a JSON file. The result is immutable and keeps key order.
"""
read_json(path::AbstractString) = JSON3.read(read(path, String))

"""
    write_json(path, obj; indent = 1)

Write `obj` as pretty-printed JSON with LF line endings (PHPS uses `indent=1`).
"""
function write_json(path::AbstractString, obj; indent::Integer = 1)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        write_json(io, obj; indent)
    end
    return path
end

function write_json(io::IO, obj; indent::Integer = 1)
    _write_pretty(io, obj, 0, indent)
    write(io, '\n')
    return nothing
end

_write_pretty(io::IO, x, level, indent) = JSON3.write(io, x)

function _write_pretty(io::IO, x::AbstractDict, level, indent)
    if isempty(x)
        write(io, "{}")
        return
    end
    write(io, "{\n")
    first = true
    for (k, v) in x
        first || write(io, ",\n")
        first = false
        write(io, ' '^(indent * (level + 1)))
        JSON3.write(io, string(k))
        write(io, ": ")
        _write_pretty(io, v, level + 1, indent)
    end
    write(io, '\n', ' '^(indent * level), '}')
end

function _write_pretty(io::IO, x::AbstractVector, level, indent)
    if isempty(x)
        write(io, "[]")
        return
    end
    write(io, "[\n")
    for (i, v) in enumerate(x)
        i > 1 && write(io, ",\n")
        write(io, ' '^(indent * (level + 1)))
        _write_pretty(io, v, level + 1, indent)
    end
    write(io, '\n', ' '^(indent * level), ']')
end

"""
    json_identical(a, b) -> Bool

Structural identity of two parsed JSON values: same key order, same value types (`1` and
`1.0` differ) and the same values (`isequal`, so `-0.0 != 0.0`).
"""
function json_identical(a::AbstractDict, b::AbstractDict)
    length(a) == length(b) || return false
    for ((ka, va), (kb, vb)) in zip(a, b)
        string(ka) == string(kb) || return false
        json_identical(va, vb) || return false
    end
    return true
end

function json_identical(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return false
    return all(json_identical(x, y) for (x, y) in zip(a, b))
end

json_identical(a::AbstractString, b::AbstractString) = a == b
json_identical(a::Nothing, b::Nothing) = true
json_identical(a::Bool, b::Bool) = a == b
json_identical(a::Number, b::Number) = typeof(a) == typeof(b) && isequal(a, b)
json_identical(a, b) = false

"""
    JSONPathError

A required key is missing, or a value has the wrong type, while reading a case.
"""
struct JSONPathError <: Exception
    file::String
    path::String
    msg::String
end
Base.showerror(io::IO, e::JSONPathError) = print(io, "JSONPathError: ", e.file, " at ", e.path,
                                                 ": ", e.msg)
