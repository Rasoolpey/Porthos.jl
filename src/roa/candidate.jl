# Lyapunov candidates (TODO.md, decision of 2026-10-01: a candidate-independent certificate
# pipeline). A candidate is a function V on the section coordinates of a `SectionModel`,
# centred at the true equilibrium eta* (enclosed by `enclose_equilibrium`, never known
# exactly): it is evaluated at the offset xi = eta - eta*. Every candidate goes through the
# same proof gates (`certify_level`); only these methods differ.
#
#   candidate_model(c)                 the SectionModel (coordinates, physical projection)
#   candidate_value(c, xi)             V(eta* + xi), generic in the number type
#   candidate_gradient(c, xi)          grad V(eta* + xi)
#   positivity_proof(c)                V > 0 off the equilibrium, proved
#   sublevel_half_widths(c, level)     w with {V(eta* + xi) <= level} inside |xi| <= w, rigorous
#   gradient_matrix_hull(c, xi_box)    an interval matrix N with grad V(eta* + xi) = N xi for
#                                      some N in it, for every xi in the box (V has its
#                                      minimum at eta*, so grad V(eta*) = 0)
#   candidate_fingerprint(c)           SHA-256 of everything that defines V and its coordinates
#   candidate_record(c)                a Dict for the certificate record
#
# The first implementation is the quadratic V_P = xi' P xi; the extended storage H_ext
# (Target B) is the second.

"""
    LyapunovCandidate

Abstract type of the Lyapunov candidates of the ROA pipeline; see the comment at the top of
`src/roa/candidate.jl` for the methods a candidate provides.
"""
abstract type LyapunovCandidate end

function candidate_model end
function candidate_value end
function candidate_gradient end
function positivity_proof end
function sublevel_half_widths end
function gradient_matrix_hull end
function candidate_fingerprint end
function candidate_record end

"""
    lift(c::LyapunovCandidate, eta) -> x

The full physical state at section coordinates `eta`.
"""
lift(c::LyapunovCandidate, eta::AbstractVector) = lift(candidate_model(c), eta)

# SHA-256 over named arrays, with their type and shape.
function _digest(items::Pair...)
    ctx = SHA.SHA256_CTX()
    for (k, v) in items
        SHA.update!(ctx, Vector{UInt8}(string(k)))
        a = v isa AbstractArray ? collect(v) : [v]
        SHA.update!(ctx, Vector{UInt8}(string(eltype(a), size(a))))
        if eltype(a) <: Union{Float64,Int,Bool}
            SHA.update!(ctx, collect(reinterpret(UInt8, vec(a))))
        else
            SHA.update!(ctx, Vector{UInt8}(join(string.(a), "\n")))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

"""
    system_digest(sys::DAESystem) -> String

SHA-256 of the healthy-network model: component types, names and parameters (after
initialisation), wiring, Y-bus, load admittances, fixed-voltage buses and the COI constants.
The fault shunts are left out (the certificate is for the healthy network).
"""
function system_digest(sys::DAESystem)
    la = sys.load
    return _digest("types" => model_type.(sys.comps), "names" => name.(sys.comps),
                   "params" => [repr(params(c)) for c in sys.comps],
                   "sources" => [repr(s) for s in sys.sources], "inj" => repr(sys.inj),
                   "G" => sys.G, "B" => sys.B,
                   "load" => hcat(la.G, la.B, la.P, la.Q, la.kpf, la.kqf),
                   "slack" => sys.slack, "Vd_ref" => sys.Vd_ref, "Vq_ref" => sys.Vq_ref,
                   "coi_members" => sys.coi_members, "coi_weights" => sys.coi_weights,
                   "coi_total" => sys.coi_total, "omega_b" => sys.omega_b)
end

"""
    model_fingerprint(m::SectionModel) -> String

SHA-256 of the model (`system_digest`) and of the section coordinates: reference point, kept
states, section, scale and names.
"""
model_fingerprint(m::SectionModel) =
    _digest("system" => system_digest(m.sys), "x0" => m.x0, "V0" => m.V0, "keep" => m.keep,
            "z" => m.z, "ref" => m.ref, "l" => m.l, "angle" => collect(m.angle),
            "scale" => m.scale, "names" => m.names)

# -- the quadratic candidate V_P ----------------------------------------------------------

"""
    QuadraticCandidate

`V(eta* + xi) = xi' P xi` on the section coordinates of `model`, `P` exactly symmetric.
Built by `quadratic_candidate`.
"""
struct QuadraticCandidate <: LyapunovCandidate
    model::AbstractSectionModel
    P::Matrix{Float64}
    construction::Dict{String,Any}
end

Base.show(io::IO, c::QuadraticCandidate) =
    print(io, "QuadraticCandidate(", size(c.P, 1), " coordinates, ", c.construction["method"], ")")

"""
    quadratic_candidate(m::AbstractSectionModel; P = nothing) -> QuadraticCandidate

The local quadratic candidate `V_P`. By default `P` solves `A'P + PA = -I` with `A` the
Jacobian of the section field at the reference point (`section_jacobian`); a given `P` (for
instance from an optimiser) is only symmetrised. Either way `P` is a candidate: its
positivity and every decay claim are proved later.
"""
function quadratic_candidate(m::AbstractSectionModel; P::Union{Nothing,AbstractMatrix} = nothing)
    A = section_jacobian(m)
    n = size(A, 1)
    abscissa = maximum(real, eigvals(A))
    if P === nothing
        abscissa < 0 || throw(ArgumentError(
            "the section Jacobian is not Hurwitz (max Re eig = $abscissa): no quadratic Lyapunov function"))
        P0 = lyap(Matrix(A'), Matrix(1.0I, n, n))
        method = "A'P + PA = -I at the reference point"
    else
        size(P) == (n, n) || throw(DimensionMismatch("P must be $n x $n"))
        P0 = Matrix{Float64}(P)
        method = "given"
    end
    Ps = (P0 + P0') / 2                     # exactly symmetric: (a + b)/2 == (b + a)/2
    Q = -(A' * Ps + Ps * A)
    ep = eigvals(Symmetric(Ps))
    eq = eigvals(Symmetric((Q + Q') / 2))
    construction = Dict{String,Any}(
        "method" => method, "dimension" => n, "spectral_abscissa" => abscissa,
        "P_eig_min" => first(ep), "P_eig_max" => last(ep), "P_condition" => last(ep) / first(ep),
        "decay_matrix_eig_min" => first(eq), "decay_matrix_eig_max" => last(eq),
        "note" => "floating-point construction data, not proof inputs")
    return QuadraticCandidate(m, Ps, construction)
end

candidate_model(c::QuadraticCandidate) = c.model
candidate_value(c::QuadraticCandidate, xi::AbstractVector) = dot(xi, c.P * xi)
candidate_gradient(c::QuadraticCandidate, xi::AbstractVector) = 2 .* (c.P * xi)

"""
    positivity_proof(c::QuadraticCandidate) -> NamedTuple

`P > 0` proved by `verified_min_eig` (a rigorous lower bound on its smallest eigenvalue).
"""
function positivity_proof(c::QuadraticCandidate)
    lo = verified_min_eig(ival.(c.P))
    return (proved = lo > 0, lower_bound = lo,
            method = "lambda_min(P) >= $lo: Gershgorin on X'(P - cI)X, X approximate eigenvectors (Sylvester inertia)")
end

"""
    sublevel_half_widths(c::QuadraticCandidate, level) -> Vector{Float64}

`|xi_i| <= sqrt(level (P^-1)_ii)` on `{xi' P xi <= level}`, with `(P^-1)_ii` bounded above
by an interval enclosure of `P^-1` and the square root rounded up.
"""
function sublevel_half_widths(c::QuadraticCandidate, level::Real)
    d = get!(() -> verified_inverse_diagonal(c.P), c.construction, "inverse_diagonal_upper")
    return ellipsoid_half_widths(level, d)
end

"""
    gradient_matrix_hull(c::QuadraticCandidate, xi_box) -> Matrix{Interval}

`grad V(eta* + xi) = 2P xi` exactly: the thin matrix `2P`.
"""
gradient_matrix_hull(c::QuadraticCandidate, xi_box::AbstractVector) = ival.(2 .* c.P)

candidate_fingerprint(c::QuadraticCandidate) =
    _digest("kind" => "quadratic", "model" => model_fingerprint(c.model), "P" => c.P)

candidate_record(c::QuadraticCandidate) = Dict{String,Any}(
    "kind" => "quadratic V_P = xi' P xi",
    "fingerprint" => candidate_fingerprint(c),
    "P_columns" => [c.P[:, j] for j in axes(c.P, 2)],
    "construction" => Dict(k => v for (k, v) in c.construction if k != "inverse_diagonal_upper"))
