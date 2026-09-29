# The coordinates of the ROA certificate (roadmap P11; TODO.md, decisions of 2026-10-01).
#
# The certificate is stated for the network-reduced dynamics on the common-angle section,
# with the bus voltages kept explicit ("network: hybrid"): the proof carries V as unknowns
# on a certified KCL branch instead of eliminating them symbolically.
#
# Coordinates. `keep` are the physical coordinates of `physical_projection` (the one-way
# reservoirs, the identically held states and the delta_COI monitor removed). The section
# fixes the COI-weighted angle `l'x_keep = l'x0_keep`, and is parametrised, as in
# `reference_section`, by every kept coordinate except one reference rotor angle, scaled by
# powers of two (exact in floating point):
#
#   x[keep[z_j]] = x0[keep[z_j]] + s_j eta_j,
#   x[keep[ref]] = x0[keep[ref]] - sum_j l_{z_j} s_j eta_j / l_ref,
#   every other state at its value in x0.
#
# Field. The model is invariant under a common rotation of every rotor angle and every bus
# voltage phasor (declared per model type by `rotation_action`, checked on the wiring by
# `check_rotation_symmetry`).
# The quotient by this symmetry is the section, and its dynamics are the field projected
# along the rotation direction `a` (1 on every rotor angle):
#
#   eta' = h(eta, V) = [(I - a l' / (l'a)) f_keep(x(eta), V)]_z ./ s,    0 = g(x(eta), V).
#
# The projection makes the section exactly invariant. It matters: PHPS's rounded COI weights
# do not sum exactly to the rounded total the model divides by, so the unprojected field
# drifts off the section (the equilibrium is a relative equilibrium rotating at about
# 1e-10 rad/s).

"""
    AbstractSectionModel

A system the certificate pipeline works on: section coordinates `eta` (`neta`), algebraic
unknowns `V` (`nvolt`) with a reference value `V0`, and `section_field(m, eta, V) -> (h, g)`
generic in the number type. `SectionModel` (a power system on its common-angle section) and
`AnalyticModel` (closed-form test systems) implement it.
"""
abstract type AbstractSectionModel end

"""
    SectionModel

The network-reduced dynamics on the common-angle section in scaled coordinates `eta`, with
the bus voltages explicit (see the comment at the top of `src/roa/section.jl`). Built by
`section_model`.
"""
struct SectionModel <: AbstractSectionModel
    sys::DAESystem
    x0::Vector{Float64}          # full state at the reference point (eta = 0)
    V0::Vector{Float64}          # bus voltages there
    projection::PhysicalProjection
    keep::Vector{Int}            # physical coordinates (state indices)
    z::Vector{Int}               # section coordinates (positions in keep)
    ref::Int                     # the reference rotor angle (position in keep)
    l::Vector{Float64}           # COI weights on keep (nonzero on the rotor angles)
    angle::BitVector             # rotor angles on keep: the rotation direction a
    la::Float64                  # l'a in Float64 ...
    la_iv::Ival                  # ... and enclosed
    scale::Vector{Float64}       # powers of two
    names::Vector{String}        # the state behind each eta coordinate
    init_params::Dict{String,Dict{String,Float64}}   # set by initialisation (for the record)
end

Base.show(io::IO, m::SectionModel) =
    print(io, "SectionModel(", length(m.z), " coordinates, ", nalg(m.sys),
          " voltages; reference angle ", state_names(m.sys)[m.keep[m.ref]], ")")

neta(m::SectionModel) = length(m.z)
nvolt(m::SectionModel) = nalg(m.sys)

# -- the rotation action --------------------------------------------------------------------
#
# The common-angle quotient needs the model to be equivariant under the rotation R(theta):
# every bus voltage phasor and every network-frame current turns by theta, every rotor angle
# shifts by theta, and everything else is unchanged. This is a property of each model's
# equations, declared per model type below after reading its code, and of the wiring and the
# network, which `check_rotation_symmetry` verifies structurally. A model type without a
# declaration is rejected: in particular the grid-forming converters (GFM_VOC's internal
# Cartesian voltage and current states would turn as vectors, not shift as an angle) until
# their action is written down and reviewed.

"""
    rotation_action(c) -> NamedTuple or nothing

How component `c` transforms under the common rotation: `angles`, the states shifted by
theta; `inputs` and `outputs`, the (d, q) pairs in the network frame, turned by theta; every
other state, input and output is invariant. `nothing` when the model type declares no action
(`section_model` then rejects the system).
"""
rotation_action(::AbstractComponent) = nothing
# GENROU / GENSAL: the Park transform uses (V_Re, V_Im) and delta only through
# delta - angle(V); the Norton current and the terminal current are network-frame pairs;
# omega, Pe, Qe, id_dq, iq_dq and i_fd are frame-independent.
rotation_action(::Union{GENROU_PHTRUE,GENSAL_PHTRUE}) =
    (angles = ["delta"], inputs = [("Vd", "Vq")], outputs = [("Id", "Iq"), ("It_Re", "It_Im")])
# COMPLEXLOAD: its state filters |V|; its current is -(Gp - j Bq) V with Gp, Bq functions of
# |V| and the state, so it turns with V; P and Q are invariant.
rotation_action(::COMPLEXLOAD) = (angles = String[], inputs = [("Vd", "Vq")], outputs = [("Id", "Iq")])
# Exciters and governors see only |V|, speeds, powers and set-points.
rotation_action(::Union{IEEET1_PHTRUE,IEEEG1_PHTRUE,IEEEG3_PHTRUE}) =
    (angles = String[], inputs = Tuple{String,String}[], outputs = Tuple{String,String}[])

"""
    check_rotation_symmetry(sys) -> Vector{Int}

Verifies that the assembled system is equivariant under the common rotation, given the
per-model `rotation_action`s: no fixed-voltage bus; every network-frame input pair is wired
to a bus voltage pair or to a network-frame output pair, and every invariant input to an
invariant source (a constant, |V|, a machine's dq-frame voltage, an invariant output);
every KCL injection is a network-frame output pair. The Y-bus, the fault shunts and the
frequency-dependent loads (through the COI speed) are linear or invariant. Returns the
state indices of the rotating angles. Throws `ArgumentError` otherwise.
"""
function check_rotation_symmetry(sys::DAESystem)
    any(sys.slack) && throw(ArgumentError(
        "a fixed-voltage bus breaks the rotation symmetry the common-angle section needs"))
    acts = map(sys.comps) do c
        a = rotation_action(c)
        a === nothing && throw(ArgumentError(
            "$(name(c)): no rotation action is declared for $(model_type(c)); the common-angle " *
            "section is not available for this model type yet"))
        a
    end
    function idx(names, s, c)
        j = findfirst(==(s), names)
        j === nothing && throw(ArgumentError("$(name(c)): its rotation action names $s, which it does not have"))
        return j
    end
    # (component, output port) -> :d / :q, for the network-frame outputs
    frame_out = Dict{Tuple{Int,Int},Symbol}()
    for (k, c) in enumerate(sys.comps), (d, q) in acts[k].outputs
        frame_out[(k, idx(output_names(c), d, c))] = :d
        frame_out[(k, idx(output_names(c), q, c))] = :q
    end
    rotating = Int[]
    for (k, c) in enumerate(sys.comps)
        for s in acts[k].angles
            push!(rotating, sys.offsets[k] + idx(state_names(c), s, c) - 1)
        end
        role = Dict{Int,Tuple{Symbol,Int}}()          # input -> (:d/:q, pair number)
        for (p, (d, q)) in enumerate(acts[k].inputs)
            role[idx(input_names(c), d, c)] = (:d, p)
            role[idx(input_names(c), q, c)] = (:q, p)
        end
        pair_src = Dict{Int,Any}()
        for (j, s) in enumerate(sys.sources[k])
            where_ = "$(name(c)).$(input_names(c)[j])"
            if haskey(role, j)
                ax, p = role[j]
                ok = (ax === :d && s.kind === SRC_VD) || (ax === :q && s.kind === SRC_VQ) ||
                     (s.kind === SRC_OUTPUT && get(frame_out, (s.index, s.port), nothing) === ax)
                ok || throw(ArgumentError("$where_ must read the $(ax) part of a network-frame pair"))
                key = s.kind === SRC_OUTPUT ? (:out, s.index) : (:bus, s.index)
                haskey(pair_src, p) ? (pair_src[p] == key ||
                    throw(ArgumentError("$where_: its network-frame pair reads two different sources"))) :
                    (pair_src[p] = key)
            else
                ok = s.kind in (SRC_ZERO, SRC_CONST, SRC_VTERM) ||
                     (s.kind in (SRC_DQ_VD, SRC_DQ_VQ) && "delta" in acts[s.index].angles &&
                      state_names(sys.comps[s.index])[1] == "delta") ||
                     (s.kind === SRC_OUTPUT && !haskey(frame_out, (s.index, s.port)))
                ok || throw(ArgumentError("$where_ is frame-independent but reads a network-frame signal"))
            end
        end
        b, jd, jq = sys.inj[k]
        b > 0 && (get(frame_out, (k, jd), nothing) === :d && get(frame_out, (k, jq), nothing) === :q ||
                  throw(ArgumentError("$(name(c)): its KCL injection is not a network-frame output pair")))
    end
    return rotating
end

"""
    section_model(sys, x0, V0; projection, scale = :lyapunov) -> SectionModel

The section model at the reference point `(x0, V0)` (normally the equilibrium from
`solve_equilibrium`). `scale`: a vector of powers of two, `:none`, or `:lyapunov` (the
diagonal of the dense Lyapunov solution of the unscaled section Jacobian, rounded to powers
of two by `pow2_scaling`, which makes the candidate `P` well conditioned).

Fails closed unless the rotation symmetry holds (`check_rotation_symmetry`: every model type
declares its `rotation_action`, and the wiring and network respect it) and the rotating
angles are exactly the COI sources' angles, all physical coordinates.
"""
function section_model(sys::DAESystem, x0::AbstractVector, V0::AbstractVector;
                       projection::PhysicalProjection = physical_projection(sys, x0, V0),
                       scale = :lyapunov,
                       init_params::AbstractDict = Dict{String,Dict{String,Float64}}())
    ip = Dict{String,Dict{String,Float64}}(string(k) => Dict{String,Float64}(string(q) => Float64(v) for (q, v) in d)
                                           for (k, d) in init_params)
    keep = projection.keep
    names = state_names(sys)
    rotating = check_rotation_symmetry(sys)
    angle = BitVector([i in rotating for i in keep])
    all(i -> i in keep, rotating) || throw(ArgumentError(
        "a rotating angle state is not a physical coordinate: " *
        join(names[filter(i -> !(i in keep), rotating)], ", ")))
    l = projection.coi[keep]
    for (k, i) in enumerate(keep)
        angle[k] == (l[k] != 0) || throw(ArgumentError(
            "$(names[i]): the rotating angles must be exactly the COI sources' angles"))
    end
    count(angle) >= 2 || throw(ArgumentError("the common-angle section needs at least two rotor angles"))
    ref = argmax(l)
    z = [k for k in eachindex(keep) if k != ref]
    la = sum(l[angle])
    la_iv = sum(ival.(l[angle]))
    m = SectionModel(sys, collect(Float64, x0), collect(Float64, V0), projection, keep, z, ref, l,
                     angle, la, la_iv, ones(length(z)), names[keep[z]], ip)
    s = scale === :none ? ones(length(z)) :
        scale === :lyapunov ? _lyapunov_scale(m) :
        scale isa AbstractVector ? collect(Float64, scale) :
        throw(ArgumentError("scale: a vector of powers of two, :none or :lyapunov"))
    length(s) == length(z) || throw(DimensionMismatch("scale has $(length(s)) entries, need $(length(z))"))
    all(v -> v > 0 && exp2(round(log2(v))) == v, s) ||
        throw(ArgumentError("the scale must be powers of two (exact in floating point)"))
    return SectionModel(sys, m.x0, m.V0, projection, keep, z, ref, l, angle, la, la_iv, s, m.names, ip)
end

"""
    section_model(eq::EquilibriumResult; kwargs...) -> SectionModel

The section model at the equilibrium of `solve_equilibrium` (its system, state, voltages and
initialisation parameters).
"""
section_model(eq::EquilibriumResult; kwargs...) =
    section_model(eq.sys, eq.x, eq.V; init_params = eq.init_params, kwargs...)

function _lyapunov_scale(m::SectionModel)
    A = section_jacobian(m)
    Q0 = lyap(Matrix(A'), Matrix(1.0I, size(A)...))
    return diag(pow2_scaling(diag(Q0)))
end

"""
    lift(m, eta) -> x

The full state at section coordinates `eta` (in `eta`'s number type).
"""
function lift(m::SectionModel, eta::AbstractVector{T}) where {T}
    x = Vector{T}(undef, length(m.x0))
    for i in eachindex(x)
        x[i] = m.x0[i]
    end
    acc = zero(T)
    for (j, k) in enumerate(m.z)
        d = m.scale[j] * eta[j]
        i = m.keep[k]
        x[i] = m.x0[i] + d
        m.l[k] != 0 && (acc += m.l[k] * d)
    end
    iref = m.keep[m.ref]
    x[iref] = m.x0[iref] - acc / m.l[m.ref]
    return x
end

"""
    section_coordinates(m, x) -> eta

The section coordinates of a full state on the section (the inverse of `lift` there; the
states outside `keep` and the reference angle are ignored).
"""
section_coordinates(m::SectionModel, x::AbstractVector) =
    [(x[m.keep[k]] - m.x0[m.keep[k]]) / m.scale[j] for (j, k) in enumerate(m.z)]

"""
    section_field(m, eta, V) -> (h, g)

The projected section field `h` and the KCL residual `g` (see the file comment), in the
promoted number type of `eta` and `V` (Float64, intervals, duals over either). With
`logs = (out = [...], step = [...])`, per-component `ModeLog`s record the branch decisions.
"""
function section_field(m::SectionModel, eta::AbstractVector, V::AbstractVector; logs = nothing)
    T = promote_type(eltype(eta), eltype(V))
    f, g = state_field(m, lift(m, convert(AbstractVector{T}, eta)), V; logs)
    return project_field(m, f), g
end

"""
    state_field(m, x, V; logs = nothing) -> (f, g)

The unprojected field and KCL residual at a full state `x` (no fault).
"""
function state_field(m::SectionModel, x::AbstractVector, V::AbstractVector; logs = nothing)
    T = promote_type(eltype(x), eltype(V))
    sys = m.sys
    f, g = zeros(T, sys.n_diff), zeros(T, nalg(sys))
    dae_residual!(f, g, sys, DAEWorkspace(sys, T), x, V, fill(false, length(sys.faults)), logs)
    return f, g
end

"""
    project_field(m, f) -> h

The section field from the full field: `[(I - a l'/(l'a)) f_keep]_z ./ s`.
"""
function project_field(m::SectionModel, f::AbstractVector{T}) where {T}
    lf = zero(T)
    for (k, i) in enumerate(m.keep)
        m.l[k] != 0 && (lf += m.l[k] * f[i])
    end
    w = lf / (is_interval_type(T) ? m.la_iv : m.la)
    h = Vector{T}(undef, length(m.z))
    for (j, k) in enumerate(m.z)
        fk = f[m.keep[k]]
        h[j] = (m.angle[k] ? fk - w : fk) / m.scale[j]
    end
    return h
end

"""
    section_residual(m, w) -> [h; g]

`section_field` on the stacked unknowns `w = [eta; V]`.
"""
function section_residual(m::AbstractSectionModel, w::AbstractVector)
    n = neta(m)
    h, g = section_field(m, view(w, 1:n), view(w, n + 1:length(w)))
    return vcat(h, g)
end

"""
    section_jacobian(m[, eta, V]) -> Matrix

The Jacobian of the network-reduced section field `eta -> h(eta, V(eta))` at `(eta, V)`
with `g(eta, V) = 0` (default: `eta = 0`, `V = V0`), in Float64:
`h_eta - h_V g_V^-1 g_eta`. For constructing candidates; the proof uses `jacobian_hull`.
"""
function section_jacobian(m::AbstractSectionModel, eta::AbstractVector = zeros(neta(m)),
                          V::AbstractVector = m.V0)
    n = neta(m)
    J = ForwardDiff.jacobian(w -> section_residual(m, w), vcat(collect(float(eta)), collect(float(V))))
    return J[1:n, 1:n] - J[1:n, n+1:end] * (J[n+1:end, n+1:end] \ J[n+1:end, 1:n])
end
