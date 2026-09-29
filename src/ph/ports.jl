# Port passivity of a single component (roadmap P10; PHPS `study/src/audit_controller_kyp.py`
# and `audit_governor_nonpassivity_exact.py`).
#
# A component's own equations are linearised at its equilibrium for one input and one
# output (for a governor: speed in, mechanical power out). Its one-way reservoir states are
# left out: they do not act on any output, so the transfer is unchanged. The port is
# passive for some storage only if H = sign * G is positive real, in particular
# Re H(jw) >= 0 at every w. A frequency with Re H(jw) < 0 is therefore a certificate that
# no storage, quadratic or not, makes the port passive, and that the positive-real (KYP)
# LMI  P > 0, A'P + PA <= 0, PB = C_h'  is infeasible.

"""
    PortModel

The small-signal model of one component at one port: `x' = A x + B u`, `y = C x + D u`,
over the states `states`, for input `input` and output `output`.
"""
struct PortModel
    component::String
    states::Vector{String}
    input::String
    output::String
    A::Matrix{Float64}
    B::Vector{Float64}
    C::Vector{Float64}
    D::Float64
end

Base.show(io::IO, m::PortModel) =
    print(io, "PortModel(", m.component, ": ", m.input, " -> ", m.output, ", ",
          length(m.states), " states)")

"""
    port_model(c, x, u; input, output, contracts) -> PortModel

Linearise component `c` at its state `x` and inputs `u` from `input` to `output` (names
from `input_names(c)` and `output_names(c)`; the output as the component's output kernel
computes it). The contract's reservoir state is left out.
"""
function port_model(c::AbstractComponent, x::AbstractVector, u::AbstractVector;
                    input::AbstractString, output::AbstractString,
                    contracts::ContractSet = default_contracts())
    snames = state_names(c)
    res = get(contract(contracts, model_type(c)).raw, :reservoir, nothing)
    keep = [j for (j, s) in enumerate(snames)
            if !(res isa JSON3.Object && s == String(res[:state]))]
    ki = findfirst(==(input), input_names(c))
    ko = findfirst(==(output), output_names(c))
    ki === nothing && throw(ArgumentError("$(name(c)) has no input $input"))
    ko === nothing && throw(ArgumentError("$(name(c)) has no output $output"))
    function F(z)
        n = length(keep)
        xx = collect(promote_type(eltype(z), eltype(x)), x)
        xx[keep] .= view(z, 1:n)
        uu = collect(promote_type(eltype(z), eltype(u)), u)
        uu[ki] = z[n + 1]
        dx = rhs!(similar(xx), c, xx, uu)
        y = outputs!(zeros(eltype(xx), noutputs(c)), c, xx, uu)
        return vcat(dx[keep], y[ko])
    end
    n = length(keep)
    J = ForwardDiff.jacobian(F, vcat(collect(float(x[keep])), float(u[ki])))
    return PortModel(name(c), snames[keep], String(input), String(output),
                     J[1:n, 1:n], J[1:n, n + 1], J[n + 1, 1:n], J[n + 1, n + 1])
end

"""`G(s) = C (sI - A)^{-1} B + D` of a port model."""
transfer(m::PortModel, s::Number) = dot(m.C, (s * I - m.A) \ complex.(m.B)) + m.D

"""
    real_part_crossings(m; sign = -1, wmin = 1e-3, wmax = 1e3, n = 60001) -> Vector

The frequencies (rad/s) where `Re(sign * G(jw))` changes sign: sign changes on a
logarithmic grid, each refined by bisection to the resolution of Float64.
"""
function real_part_crossings(m::PortModel; sign::Real = -1, wmin::Real = 1e-3,
                             wmax::Real = 1e3, n::Integer = 60001)
    f(w) = real(sign * transfer(m, im * w))
    ws = exp10.(range(log10(wmin), log10(wmax); length = n))
    vals = f.(ws)
    roots = Float64[]
    for i in 1:n - 1
        (vals[i] == 0 || Base.sign(vals[i]) == Base.sign(vals[i + 1])) && continue
        a, b, fa = ws[i], ws[i + 1], vals[i]
        while true
            mid = (a + b) / 2
            (mid <= a || mid >= b) && break
            fm = f(mid)
            if Base.sign(fm) == Base.sign(fa)
                a, fa = mid, fm
            else
                b = mid
            end
        end
        push!(roots, (a + b) / 2)
    end
    return roots
end

"""
    passivity_certificate(m; sign = -1, wmin = 1e-3, wmax = 1e3, n = 60001) -> NamedTuple

For `H = sign * G`: `H0 = H(0)`, the smallest `Re H(jw)` on the grid and its frequency, and
`nonpassive = true` when that value is negative by more than `1e-8 * max(1, |H|)` (then no
storage makes the port passive and the positive-real LMI is infeasible; the margin covers
the Float64 evaluation of a well-conditioned small system). `shortage` is `-min Re H`
(the passivity shortage index, 0 when passive on the grid).
"""
function passivity_certificate(m::PortModel; sign::Real = -1, wmin::Real = 1e-3,
                               wmax::Real = 1e3, n::Integer = 60001)
    ws = exp10.(range(log10(wmin), log10(wmax); length = n))
    Hs = [sign * transfer(m, im * w) for w in ws]
    k = argmin(real.(Hs))
    rmin = real(Hs[k])
    scale = max(1.0, maximum(abs, Hs))
    return (H0 = real(sign * transfer(m, 0.0)), min_re = rmin, at = ws[k],
            nonpassive = rmin < -1e-8 * scale, shortage = max(0.0, -rmin))
end

"""
    port_zeros(m) -> Vector{ComplexF64}

The finite zeros of `G(s)`: the finite generalized eigenvalues of the pencil
`([A B; C' D], [I 0; 0 0])`, sorted by real part. A zero in the right half-plane rules out
passivity at the port whatever the storage (PHPS roadmap section 7, the hydro governor).
"""
function port_zeros(m::PortModel)
    n = length(m.B)
    Mz = [m.A m.B; m.C' m.D]
    N = [Matrix{Float64}(I, n, n) zeros(n); zeros(1, n + 1)]
    return sort!(filter(isfinite, eigvals(Mz, N)); by = real)
end
