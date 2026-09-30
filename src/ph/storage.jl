# Assembled storage and the network-reduced field (roadmap P10).
#
# The storage is the sum of the component Hamiltonians, as PHPS's HamiltonianAssembler forms
# it: every component whose Hamiltonian is nonzero (the machines, and the PHTRUE exciters and
# governors, whose accounting reservoirs carry storage); the loads contribute nothing.
#
# The network-reduced field is f(x) = f(x, V(x)) with V(x) the solution of KCL g(x, V) = 0
# near a given voltage, which is what PHPS's Python solver evaluates (its network solve
# converges to round-off). Its Jacobian is exact here: A = f_x - f_V g_V^{-1} g_x by the
# implicit function theorem, with the partial derivatives from ForwardDiff.

"""
    storage_components(sys) -> Vector{Int}

The components that carry storage: those whose Hamiltonian is nonzero at the test state
`x = 0.01` (PHPS's `HamiltonianAssembler` rule).
"""
storage_components(sys::DAESystem) =
    [k for (k, c) in enumerate(sys.comps)
     if nstates(c) > 0 && hamiltonian(c, fill(0.01, nstates(c))) != 0]

_state_range(sys::DAESystem, k::Integer) = sys.offsets[k]:(sys.offsets[k] + nstates(sys.comps[k]) - 1)

"""
    total_hamiltonian(sys, x) -> Real

`H(x)`, the sum of the storage components' Hamiltonians.
"""
function total_hamiltonian(sys::DAESystem, x::AbstractVector)
    h = zero(eltype(x))
    for k in storage_components(sys)
        h += hamiltonian(sys.comps[k], view(x, _state_range(sys, k)))
    end
    return h
end

"""
    grad_total_hamiltonian(sys, x) -> Vector

`grad H(x)` on the full state (zero outside the storage components).
"""
function grad_total_hamiltonian(sys::DAESystem, x::AbstractVector)
    g = zeros(eltype(x), length(x))
    for k in storage_components(sys)
        r = _state_range(sys, k)
        g[r] .= grad_hamiltonian(sys.comps[k], x[r])
    end
    return g
end

"""
    hessian_total_hamiltonian(sys, x) -> Matrix

`Hess H(x)` (ForwardDiff of the gradient, symmetrised).
"""
function hessian_total_hamiltonian(sys::DAESystem, x::AbstractVector)
    S = ForwardDiff.jacobian(z -> grad_total_hamiltonian(sys, z), collect(float(x)))
    return 0.5 .* (S .+ S')
end

"""
    solve_network(sys, x, V0; tol = 1e-13, maxiter = 50, faults_on = false) -> V

The bus voltages with KCL `g(x, V) = 0` at the states `x` (with every fault shunt of the
scenario when `faults_on`), by Newton from `V0`.
As in PHPS's network solve, one more Newton step is taken once the residual is below
`tol`, so the result is at round-off and does not depend on the stopping rule.
"""
function solve_network(sys::DAESystem, x::AbstractVector, V0::AbstractVector;
                       tol::Real = 1e-13, maxiter::Integer = 50, faults_on::Bool = false)
    nd = sys.n_diff
    V = collect(Float64, V0)
    G(v) = dae_residual(sys, x, v; faults_on)[2]
    for _ in 1:maxiter
        g = G(V)
        r = maximum(abs, g)
        isfinite(r) || error("solve_network: non-finite residual")
        V .-= ForwardDiff.jacobian(G, V) \ g
        r <= tol && return V
    end
    error("solve_network: KCL did not converge in $maxiter iterations")
end

"""
    reduced_field(sys, x, V0) -> (f, V)

The network-reduced field `f(x) = f(x, V(x))` and the voltages `V(x)`, solved from `V0`.
"""
function reduced_field(sys::DAESystem, x::AbstractVector, V0::AbstractVector; kwargs...)
    V = solve_network(sys, x, V0; kwargs...)
    return dae_residual(sys, x, V)[1], V
end

"""
    reduced_jacobian(sys, x, V) -> Matrix

The exact Jacobian of the network-reduced field at `(x, V)` (with `g(x, V) = 0`):
`A = f_x - f_V g_V^{-1} g_x`.
"""
function reduced_jacobian(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    nd = sys.n_diff
    F(z) = vcat(dae_residual(sys, view(z, 1:nd), view(z, nd + 1:length(z)))...)
    J = ForwardDiff.jacobian(F, vcat(collect(float(x)), collect(float(V))))
    fx, fV = J[1:nd, 1:nd], J[1:nd, nd + 1:end]
    gx, gV = J[nd + 1:end, 1:nd], J[nd + 1:end, nd + 1:end]
    return fx - fV * (gV \ gx)
end
