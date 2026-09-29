# Shared helpers for the parity gates.

using Pkg

"""
The verified parity pack. `PORTHOS_PARITY_PACK` points at a local pack directory
(used while generating a new baseline); otherwise the lazy artifact is installed from
Artifacts.toml. A missing pack is an error, never a skip: every gate needs it.
"""
function test_parity_pack()
    if isempty(get(ENV, "PORTHOS_PARITY_PACK", ""))
        Pkg.Artifacts.ensure_artifact_installed(Porthos.PARITY_ARTIFACT, Porthos.ARTIFACTS_TOML)
    end
    return load_parity_pack()
end

const PACK = test_parity_pack()

"""`|a - b| <= atol + rtol*|b|` elementwise (roadmap Part I)."""
close_to(a, b; atol, rtol = 0.0) = all(abs.(a .- b) .<= atol .+ rtol .* abs.(b))

"""Number of entries whose binary64 components differ."""
bitdiff(A::AbstractArray{ComplexF64}, B::AbstractArray{ComplexF64}) =
    count(i -> reinterpret(UInt128, A[i]) != reinterpret(UInt128, B[i]), eachindex(A, B))
bitdiff(A::AbstractArray{Float64}, B::AbstractArray{Float64}) =
    count(i -> reinterpret(UInt64, A[i]) != reinterpret(UInt64, B[i]), eachindex(A, B))

case_path(rel) = joinpath(ROOT, "cases", rel)
