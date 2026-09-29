# Interval core of the ROA certificates (roadmap P11).
#
# Everything here is outward rounded: IntervalArithmetic.jl does the arithmetic, and a
# Float64 enters a proof only as an exact constant (a matrix entry, a power-of-two scale)
# or as the bound of an interval. Floating-point results that are not enclosed (an
# approximate inverse used as a preconditioner, approximate eigenvectors) only steer a
# method whose conclusion is then checked in interval arithmetic.

const IA = IntervalArithmetic
const Ival = IA.Interval{Float64}

"""
    ProofFailure

Thrown when a proof step does not go through (a Krawczyk test that does not contract, a
matrix not proved nonsingular, an evaluation that is not well defined on its box). It
means "not proved", never "false"; any other exception is a bug.
"""
struct ProofFailure <: Exception
    msg::String
end
Base.showerror(io::IO, e::ProofFailure) = print(io, "ProofFailure: ", e.msg)
fail(msg::AbstractString) = throw(ProofFailure(msg))

# -- branch primitives on intervals (roadmap 2.3) ------------------------------------------
#
# A comparison of intervals is decided when it has the same outcome for every pair of
# points: `a > b` is true when inf(a) > sup(b) and false when sup(a) <= inf(b). Otherwise
# the branch is undecided and the primitive throws, so an interval evaluation either runs
# in one smooth mode on the whole box or fails.

_undecided(op, a, b) = throw(UndecidedBranch("$(a) $(op) $(b) is not decided on the box"))

function _gt_val(a::IA.Interval, b::IA.Interval)
    IA.inf(a) > IA.sup(b) && return true
    IA.sup(a) <= IA.inf(b) && return false
    _undecided(">", a, b)
end
function _ge_val(a::IA.Interval, b::IA.Interval)
    IA.inf(a) >= IA.sup(b) && return true
    IA.sup(a) < IA.inf(b) && return false
    _undecided(">=", a, b)
end
_lt_val(a::IA.Interval, b::IA.Interval) = _gt_val(b, a)
_le_val(a::IA.Interval, b::IA.Interval) = _ge_val(b, a)
for op in (:_gt_val, :_ge_val, :_lt_val, :_le_val)
    @eval $op(a::IA.Interval, b::Real) = $op(a, IA.interval(b))
    @eval $op(a::Real, b::IA.Interval) = $op(IA.interval(a), b)
end

"""
    is_interval_type(T) -> Bool

Whether `T` is an interval, or a (nested) ForwardDiff dual over intervals.
"""
is_interval_type(::Type) = false
is_interval_type(::Type{<:IA.Interval}) = true
is_interval_type(::Type{ForwardDiff.Dual{G,V,N}}) where {G,V,N} = is_interval_type(V)

# -- small helpers -------------------------------------------------------------------------

ival(x::Real) = IA.interval(x)
ival(lo::Real, hi::Real) = IA.interval(lo, hi)
ival(x::IA.Interval) = x

"""Midpoints and radii of an interval array."""
imid(X::AbstractArray{<:IA.Interval}) = IA.mid.(X)

"""Outward-rounded `x + [-r, r]` for Float64 vectors (r >= 0)."""
function ibox(center::AbstractVector{<:Real}, r::AbstractVector{<:Real})
    return [ival(c) + ival(-ri, ri) for (c, ri) in zip(center, r)]
end

"""
    strictly_inside(K, X) -> Bool

`K` lies in the interior of `X`, entry by entry.
"""
strictly_inside(K::AbstractArray{<:IA.Interval}, X::AbstractArray{<:IA.Interval}) =
    all(IA.inf(x) < IA.inf(k) && IA.sup(k) < IA.sup(x) for (k, x) in zip(K, X))

"""Entrywise intersection (errors when empty: the arguments must both enclose one set)."""
function iintersect(A::AbstractArray{<:IA.Interval}, B::AbstractArray{<:IA.Interval})
    return map(A, B) do a, b
        lo, hi = max(IA.inf(a), IA.inf(b)), min(IA.sup(a), IA.sup(b))
        lo <= hi || error("iintersect: empty intersection")
        ival(lo, hi)
    end
end

"""Upper bound of the infinity norm of an interval matrix."""
inorm_inf_up(A::AbstractMatrix{<:IA.Interval}) =
    maximum(IA.sup(sum(ival.(IA.mag.(view(A, i, :))))) for i in axes(A, 1))

"""
    well_defined(y) -> Bool

Every entry is a nonempty, bounded interval whose decoration is at least `dac` (defined
and continuous on the whole box: no division by an interval containing zero, no square
root of a negative part). Duals are checked through their value and partials.
"""
well_defined(y::IA.Interval) = IA.decoration(y) >= IA.dac && isfinite(IA.inf(y)) &&
                               isfinite(IA.sup(y))
well_defined(y::ForwardDiff.Dual) = well_defined(ForwardDiff.value(y)) &&
                                    all(well_defined, ForwardDiff.partials(y))
well_defined(y::AbstractArray) = all(well_defined, y)

# -- linear systems -----------------------------------------------------------------------

"""
    interval_solve(A, B; iterations = 3) -> (X, beta)

An enclosure `X` of `{A^-1 B : A in [A], B in [B]}` (column by column) for interval
matrices, and the bound `beta >= ||I - C A||_inf` (`C` an approximate inverse of mid(A)).
`beta < 1` proves every `A` in `[A]` nonsingular; otherwise this throws. Method: with
`R = I - C[A]`, every solution satisfies `x = C b + R x`, so `|x| <= |C b| / (1 - beta)`
bounds it, and `X <- (C[B] + R X) cap X` refines the bound (Krawczyk's operator for linear
systems).
"""
function interval_solve(A::AbstractMatrix{<:IA.Interval}, B::AbstractMatrix{<:IA.Interval};
                        iterations::Integer = 3)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("A must be square"))
    C = try
        inv(imid(A))
    catch e
        e isa SingularException ? fail("interval_solve: the midpoint matrix is singular") : rethrow()
    end
    all(isfinite, C) || fail("interval_solve: the midpoint matrix is singular")
    Ci = ival.(C)
    R = ival.(Matrix{Float64}(I, n, n)) .- Ci * A
    beta = inorm_inf_up(R)
    beta < 1 || fail("interval_solve: ||I - C A|| <= $beta is not below 1 (A not proved nonsingular)")
    Z = Ci * B
    X = similar(Z)
    for j in axes(Z, 2)
        zmax = maximum(IA.mag, view(Z, :, j))
        r = IA.sup(ival(zmax) / (ival(1.0) - ival(beta)))
        for i in 1:n
            X[i, j] = ival(-r, r)
        end
    end
    for _ in 1:iterations
        X = iintersect(Z .+ R * X, X)
    end
    return X, beta
end
interval_solve(A::AbstractMatrix{<:IA.Interval}, b::AbstractVector{<:IA.Interval}; kwargs...) =
    (r = interval_solve(A, reshape(b, :, 1); kwargs...); (vec(r[1]), r[2]))

# -- symmetric definiteness ---------------------------------------------------------------

"""
    verified_max_eig(S) -> Float64

A rigorous upper bound on the largest eigenvalue of every symmetric matrix in the interval
matrix `S` (`Inf` when the check fails); `verified_min_eig` applied to `-S`.
"""
verified_max_eig(S::AbstractMatrix{<:IA.Interval}) = -verified_min_eig(.-S)

"""
    weyl_max_eig(S) -> NamedTuple

A rigorous upper bound on the largest eigenvalue of every symmetric matrix in the interval
matrix `S`, for wide `S`: by Weyl, `lambda_max(Sc + D) <= lambda_max(Sc) + ||D||_2` with
`Sc` the (Float64, symmetric) midpoint and `|D| <= Rad` entrywise, and
`||D||_2 <= rho(|D|) <= rho(Rad)` (Perron-Frobenius). `lambda_max(Sc)` is bounded by
`verified_max_eig` on the thin matrix, `rho(Rad)` by the Collatz-Wielandt bound
`max_i (Rad v)_i / v_i` for a positive vector `v` near the Perron vector (outward rounded).
Returns the bound and its two parts.
"""
function weyl_max_eig(S::AbstractMatrix{<:IA.Interval})
    Sc = IA.mid.(S)
    Sc = (Sc + Sc') / 2                                    # exactly symmetric
    Rad = [IA.mag(s - ival(c)) for (s, c) in zip(S, Sc)]   # |S - Sc| <= Rad, rounded up
    Rad = max.(Rad, Rad')
    lam = verified_max_eig(ival.(Sc))
    v = ones(size(Rad, 1))
    for _ in 1:50                                          # power iteration (steers only)
        w = Rad * v
        nw = maximum(w)
        nw > 0 || break
        v = w ./ nw
    end
    v .= max.(v, 1e-3 * maximum(v))
    Rv = ival.(Rad) * ival.(v)
    rho = maximum(IA.sup(Rv[i] / ival(v[i])) for i in eachindex(v))
    return (bound = IA.sup(ival(lam) + ival(rho)), midpoint = lam, radius = rho)
end

"""
    verified_inverse_diagonal(P) -> Vector{Float64}

Rigorous upper bounds on the diagonal of `P^-1` for a Float64 matrix `P` (which must be
proved positive definite separately): the diagonal of an interval enclosure of `P^-1`.
"""
function verified_inverse_diagonal(P::AbstractMatrix{Float64})
    n = size(P, 1)
    X, _ = interval_solve(ival.(P), ival.(Matrix{Float64}(I, n, n)))
    return [IA.sup(X[i, i]) for i in 1:n]
end

"""
    ellipsoid_half_widths(level, inv_diag_up) -> Vector{Float64}

Rigorous half-widths of the box enclosing `{xi : xi' P xi <= level}`:
`max |xi_i| = sqrt(level (P^-1)_ii)`, rounded up.
"""
ellipsoid_half_widths(level::Real, inv_diag_up::AbstractVector{Float64}) =
    [IA.sup(sqrt(ival(level) * ival(d))) for d in inv_diag_up]
