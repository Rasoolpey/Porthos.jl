# Parity pack: read-only reference numbers generated once from PHPS by
# parity/generate/generate_pack.py and bound as the lazy artifact `parity_pack` in
# Artifacts.toml. Loading verifies every file against the SHA-256 in MANIFEST.json and
# rejects files the manifest does not list.

const PARITY_ARTIFACT = "parity_pack"
const PARITY_PACK_FORMAT = 1
const ARTIFACTS_TOML = normpath(joinpath(@__DIR__, "..", "..", "Artifacts.toml"))

"""
    ParityPack

A verified parity pack: its directory and manifest.
"""
struct ParityPack
    dir::String
    manifest::JSON3.Object
end

Base.show(io::IO, p::ParityPack) =
    print(io, "ParityPack(", length(p.manifest[:files]), " files, PHPS ",
          first(string(p.manifest[:phps][:commit]), 10), ")")

struct ParityPackError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ParityPackError) = print(io, "ParityPackError: ", e.msg)

"""
    parity_pack_dir() -> String

Directory of the parity pack: `ENV["PORTHOS_PARITY_PACK"]` if set, otherwise the installed
`$(PARITY_ARTIFACT)` artifact. The artifact is lazy; install it with
`Pkg.Artifacts.ensure_artifact_installed("$(PARITY_ARTIFACT)", "$(basename(ARTIFACTS_TOML))")`
(the test suite does this).
"""
function parity_pack_dir()
    env = get(ENV, "PORTHOS_PARITY_PACK", "")
    isempty(env) || return env
    isfile(ARTIFACTS_TOML) || throw(ParityPackError("no Artifacts.toml at $ARTIFACTS_TOML"))
    h = Artifacts.artifact_hash(PARITY_ARTIFACT, ARTIFACTS_TOML)
    h === nothing && throw(ParityPackError("Artifacts.toml has no '$PARITY_ARTIFACT' entry"))
    Artifacts.artifact_exists(h) || throw(ParityPackError(
        "the parity pack artifact $(bytes2hex(h.bytes)) is not installed; run " *
        "Pkg.Artifacts.ensure_artifact_installed(\"$PARITY_ARTIFACT\", \"$ARTIFACTS_TOML\")"))
    return Artifacts.artifact_path(h)
end

function _sha256_file(path::AbstractString)
    open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

"""
    verify_parity_pack(dir) -> JSON3.Object

Check the pack format, that every file listed in `MANIFEST.json` exists with the recorded
SHA-256, and that the directory holds no unlisted file. Returns the manifest.
"""
function verify_parity_pack(dir::AbstractString)
    mpath = joinpath(dir, "MANIFEST.json")
    isfile(mpath) || throw(ParityPackError("no MANIFEST.json in $dir"))
    manifest = read_json(mpath)
    fmt = Int(get(manifest, :pack_format, 0))
    fmt == PARITY_PACK_FORMAT ||
        throw(ParityPackError("pack format $fmt, this Porthos reads $PARITY_PACK_FORMAT"))
    listed = Set{String}()
    for (rel, sha) in manifest[:files]
        rel = string(rel)
        push!(listed, rel)
        f = joinpath(dir, rel)
        isfile(f) || throw(ParityPackError("missing file $rel"))
        got = _sha256_file(f)
        got == string(sha) || throw(ParityPackError("hash mismatch for $rel: $got != $sha"))
    end
    for (root, _, files) in walkdir(dir), f in files
        rel = replace(relpath(joinpath(root, f), dir), '\\' => '/')
        rel == "MANIFEST.json" && continue
        rel in listed || throw(ParityPackError("file not in the manifest: $rel"))
    end
    return manifest
end

"""
    load_parity_pack(dir = parity_pack_dir()) -> ParityPack

Load and verify the parity pack.
"""
load_parity_pack(dir::AbstractString = parity_pack_dir()) =
    ParityPack(String(dir), verify_parity_pack(dir))

"""
    pack_json(pack, relpath) -> JSON3 value
"""
function pack_json(pack::ParityPack, rel::AbstractString)
    haskey(pack.manifest[:files], Symbol(rel)) ||
        throw(ParityPackError("$rel is not in the parity pack"))
    return read_json(joinpath(pack.dir, rel))
end

"""
    pack_cases(pack) -> Vector of (name, system, scenario) named tuples

The parity cases, with paths relative to the PHPS `phps/cases` directory (and therefore
to the Porthos `cases/` directory, which holds the same files).
"""
pack_cases(pack::ParityPack) =
    [(name = string(c[:name]), system = string(c[:system]), scenario = string(c[:scenario]))
     for c in pack.manifest[:cases]]

has_section(pack::ParityPack, s::AbstractString) = s in string.(pack.manifest[:sections])

"""
    pack_matrix(obj) -> Matrix{ComplexF64}

A complex matrix stored as row-major `{"re": [[...]], "im": [[...]]}`.
"""
function pack_matrix(obj)
    re, im = obj[:re], obj[:im]
    n, m = length(re), length(re[1])
    M = Matrix{ComplexF64}(undef, n, m)
    for i in 1:n, j in 1:m
        M[i, j] = complex(Float64(re[i][j]), Float64(im[i][j]))
    end
    return M
end

"""
    pack_binary(pack, relpath, (rows, cols)) -> Matrix{Float64}

A little-endian float64 array stored row-major (the `sim` section's `.bin` files: one row
per time, one column per signal), returned as a `rows x cols` matrix.
"""
function pack_binary(pack::ParityPack, rel::AbstractString, dims::Tuple{Integer,Integer})
    haskey(pack.manifest[:files], Symbol(rel)) ||
        throw(ParityPackError("$rel is not in the parity pack"))
    rows, cols = dims
    raw = read(joinpath(pack.dir, rel))
    length(raw) == 8 * rows * cols ||
        throw(ParityPackError("$rel has $(length(raw)) bytes, expected 8 x $rows x $cols"))
    v = ltoh.(reinterpret(Float64, raw))
    return permutedims(reshape(v, cols, rows))
end

"""
    pack_vector(obj) -> Vector{Float64}
"""
pack_vector(obj) = Float64[Float64(v) for v in obj]
