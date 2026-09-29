# Closed-form systems for the certificate pipeline (tests of the theorem on problems with
# known answers, independent of any power-system model): eta' = h(eta, v), 0 = g(eta, v),
# with the same proof steps as a `SectionModel`. The containment audit has no contract
# clauses here: it checks that the field is well defined and every branch primitive decided
# on the box.

"""
    AnalyticModel(name, h, g, n, V0; description = "")

A closed-form model for the certificate pipeline: `h(eta, v)` (length `n`) and `g(eta, v)`
(length `length(V0)`, at least 1), both generic in the number type; `V0` is the algebraic
unknown near the equilibrium `eta = 0`. Branches must go through the branch primitives
(e.g. `Porthos.clamp_mode(NoModes(), ...)`).
"""
struct AnalyticModel{H,G} <: AbstractSectionModel
    name::String
    h::H
    g::G
    n::Int
    V0::Vector{Float64}
    description::String
end

function AnalyticModel(name::AbstractString, h, g, n::Integer, V0::AbstractVector;
                       description::AbstractString = "")
    isempty(V0) && throw(ArgumentError("an AnalyticModel needs at least one algebraic unknown"))
    return AnalyticModel(String(name), h, g, Int(n), collect(Float64, V0), String(description))
end

Base.show(io::IO, m::AnalyticModel) = print(io, "AnalyticModel(", m.name, ", ", m.n, " + ", length(m.V0), ")")

neta(m::AnalyticModel) = m.n
nvolt(m::AnalyticModel) = length(m.V0)
section_field(m::AnalyticModel, eta::AbstractVector, V::AbstractVector; logs = nothing) =
    (m.h(eta, V), m.g(eta, V))

model_fingerprint(m::AnalyticModel) =
    _digest("analytic" => m.name, "description" => m.description, "n" => m.n, "V0" => m.V0)

function containment_audit(m::AnalyticModel, E::AbstractVector{<:IA.Interval},
                           V::AbstractVector{<:IA.Interval}; kwargs...)
    verdict = Dict("evaluation" => true, "limiter_mode" => true)
    records = Dict{String,Any}[]
    try
        h, g = section_field(m, E, V)
        ok = well_defined(h) && well_defined(g)
        verdict["evaluation"] = ok
        push!(records, Dict{String,Any}("category" => "evaluation", "id" => "field_well_defined",
                                        "component" => m.name, "passed" => ok, "detail" => ""))
    catch e
        e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
        verdict["limiter_mode"] = false
        verdict["evaluation"] = false
        push!(records, Dict{String,Any}("category" => "limiter_mode", "id" => "all_branches_decided",
                                        "component" => m.name, "passed" => false,
                                        "detail" => sprint(showerror, e)))
    end
    return ContainmentAudit(verdict, records, Dict{String,Any}(), Dict{String,Any}(), box_digest(E, V))
end

section_record(m::AnalyticModel) = Dict{String,Any}(
    "fingerprint" => model_fingerprint(m), "analytic" => m.name, "description" => m.description,
    "n" => m.n, "V0" => m.V0)

certificate_claim(m::AnalyticModel) =
    "Omega_c = {(eta, v): V_cand(eta - eta*) <= c, v = v(eta)} with c = verified_valid_level, " *
    "on the certified branch of g = 0: dV_cand/dt < 0 except at eta*, one smooth mode; so " *
    "Omega_c is forward invariant and inside the region of attraction of the enclosed " *
    "equilibrium eta* of $(m.name)"

model_assumptions(m::AnalyticModel) = ["h and g as written, in exact real arithmetic"]
