# Branch primitives (roadmap 2.3).
#
# Component code never branches on a state-dependent value with a bare `if`. Every such
# decision goes through the comparisons below (`gt`, `ge`, `lt`, `le`) or the limiter
# patterns built on them (`clamp_mode`, `nonwindup`, `guard_min`, `outband_relax`,
# `select`), and each decision is passed to a recorder:
#
#   - `NoModes()` in `rhs!` and `outputs!`: records nothing and compiles away;
#   - `ModeLog()` in `modes(...)`: keeps the sequence of decisions, which the P4 parity gate
#     compares with the branch sites of PHPS's kernels.
#
# On `Real` numbers (Float64, ForwardDiff.Dual) the comparisons are the ordinary ones. On
# interval types the ROA pipeline (P11) adds methods that return the decided outcome, or
# throw `UndecidedBranch` when the interval straddles the switching surface.
#
# Each primitive reproduces the exact PHPS C++ pattern, including the order of its tests,
# so the same branches are taken and the arithmetic is identical.

abstract type ModeRecorder end

"""Recorder that records nothing (used by `rhs!` and `outputs!`)."""
struct NoModes <: ModeRecorder end

"""Recorder that keeps every state-dependent decision, in evaluation order."""
struct ModeLog <: ModeRecorder
    decisions::Vector{Bool}
end
ModeLog() = ModeLog(Bool[])

@inline _record!(::NoModes, b::Bool) = b
@inline _record!(r::ModeLog, b::Bool) = (push!(r.decisions, b); b)

"""
    UndecidedBranch

Thrown by a branch primitive when its operands do not decide the branch (an interval that
straddles the switching surface).
"""
struct UndecidedBranch <: Exception
    msg::String
end
Base.showerror(io::IO, e::UndecidedBranch) = print(io, "UndecidedBranch: ", e.msg)

# Decided comparisons. Extended for interval types at P11.
@inline _gt(a, b) = a > b
@inline _ge(a, b) = a >= b
@inline _lt(a, b) = a < b
@inline _le(a, b) = a <= b

"""
    decide(rec, b::Bool) -> Bool

Record the outcome of one branch site (a combined condition such as `a >= hi && r > 0`
is one site, as in PHPS).
"""
@inline decide(rec::ModeRecorder, b::Bool) = _record!(rec, b)

"""`a > b` as one recorded branch site."""
@inline gt(rec::ModeRecorder, a, b) = decide(rec, _gt(a, b))
"""`a >= b` as one recorded branch site."""
@inline ge(rec::ModeRecorder, a, b) = decide(rec, _ge(a, b))
"""`a < b` as one recorded branch site."""
@inline lt(rec::ModeRecorder, a, b) = decide(rec, _lt(a, b))
"""`a <= b` as one recorded branch site."""
@inline le(rec::ModeRecorder, a, b) = decide(rec, _le(a, b))

"""
    select(rec, cond::Bool, a, b)

`cond ? a : b` for a state-dependent `cond` (one recorded site). Both branches are
evaluated, so both must be finite.
"""
@inline select(rec::ModeRecorder, cond::Bool, a, b) = decide(rec, cond) ? a : b

"""
    clamp_mode(rec, v, lo, hi)

PHPS clamp: `if (v > hi) v = hi; if (v < lo) v = lo;` (two sites, in that order).
"""
@inline function clamp_mode(rec::ModeRecorder, v::T, lo, hi) where {T}
    v = gt(rec, v, hi) ? convert(T, hi) : v
    v = lt(rec, v, lo) ? convert(T, lo) : v
    return v
end

"""
    nonwindup(rec, s, rate, lo, hi)

PHPS non-windup limited integrator: the rate is frozen at a limit it would cross,
`if (s >= hi && rate > 0) rate = 0; if (s <= lo && rate < 0) rate = 0;` (two sites).
`s` is the (clamped) integrator position.
"""
@inline function nonwindup(rec::ModeRecorder, s, rate::T, lo, hi) where {T}
    rate = decide(rec, _ge(s, hi) && _gt(rate, 0.0)) ? zero(T) : rate
    rate = decide(rec, _le(s, lo) && _lt(rate, 0.0)) ? zero(T) : rate
    return rate
end

"""
    outband_relax(rec, rate, s, lo, hi, T)

PHPS leak-proof integrator: a state outside `[lo, hi]` is pulled back on time constant `T`,
`if (s > hi) rate += (hi - s)/T; if (s < lo) rate += (lo - s)/T;` (two sites). `s` is the
raw (unclamped) state.
"""
@inline function outband_relax(rec::ModeRecorder, rate, s, lo, hi, Tc)
    rate = gt(rec, s, hi) ? rate + (hi - s) / Tc : rate
    rate = lt(rec, s, lo) ? rate + (lo - s) / Tc : rate
    return rate
end

"""
    guard_min(rec, v, floor)

PHPS guard `if (v < floor) v = floor;` (one site).
"""
@inline guard_min(rec::ModeRecorder, v::T, floor) where {T} =
    lt(rec, v, floor) ? convert(T, floor) : v
