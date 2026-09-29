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

"""
    phps_init_params(case, d) -> Dict(component name => Dict(param => value))

The parameters of the pack's `dae` record `d` that differ from Porthos's own processing of
the case: those PHPS's initialisation sets (the P5 gate checks that nothing else differs).
"""
function phps_init_params(case, d)
    init = Dict{String,Dict{String,Float64}}()
    for c in d[:components]
        nm = string(c[:name])
        pd = param_dict(build_component(case, component(case, nm)))
        diff = Dict{String,Float64}()
        for (k, v) in c[:params]
            key = string(k)
            if !haskey(pd, key) || !(pd[key] isa Real) || Float64(pd[key]) !== Float64(v)
                diff[key] = Float64(v)
            end
        end
        init[nm] = diff
    end
    return init
end

"""PHPS's initialised state `[x; Vd_1, Vq_1, ...]` from the pack's `dae` record `d`."""
phps_initial_state(d) =
    vcat(Float64.(d[:equilibrium][:x]),
         vec(permutedims(hcat(Float64.(d[:equilibrium][:Vd]), Float64.(d[:equilibrium][:Vq])))))
