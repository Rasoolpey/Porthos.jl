# Joint machine storage in current-corrected rotor coordinates: Bregman Hessians and the
# local exchange shortage (TODO.md step 5a, the reviewer's calculations 1 to 3).
#
# On the KCL branch the rotor energy is a function of the section coordinates through the
# stator currents, w(eta) = z(eta) - D I(eta). The gauge family
#
#   H_alpha = sum over linear rotors of  w'Qw/2 - alpha I'D'QD I/2
#
# differs only by the exact scalar F = I'D'QD I/2 (alpha = 0: the nonnegative rotor energy
# H_w; alpha = 1: z'Qz/2 - z'QD I). The equilibrium-shifted (Bregman) form of any storage S,
# S(eta) - S(eta*) - grad S(eta*)'(eta - eta*), has Hessian Hess S(eta*) at eta*, which here
# includes the stator-current curvature sum_j (Qw)_j Hess w_j: Q > 0 alone does not decide it.
#
# Locally, the incremental rotor balance is
#   d/dt of the Bregman form of w'Qw/2 = -(w - w*)' sym(-QA) (w - w*) + (Efd - Efd*) (i_fd^H - i_fd^H*)
#                                        - (w - w*)' QD di/dt,
# and at linear order around eta*, with W = dw/deta, J_I = dI/deta and A the section Jacobian,
#   rotor loss  ~ deta' W' L W deta,          L = -sym(QA) (block-diagonal)
#   exchange    ~ -deta' W' QD J_I A deta     (the shortage; sign-indefinite)
#   network     ~ deta' J_V' G J_V deta       (incremental conductance loss, J_V = dV/deta)
# `exchange_dissipation_forms` builds these quadratic forms so the shortage can be compared
# with the dissipation that could dominate it.

# the linear-rotor units: (rotor_structure, component index)
_linear_rotors(sys::DAESystem) =
    [(rs, k) for (k, c) in enumerate(sys.comps) for rs in (rotor_structure(c),) if rs !== nothing]

"""
    rotor_energy(sys, x, V0; alpha = 0.0) -> Real

`sum w'Qw/2 - alpha I'D'QD I/2` over the machines with a linear rotor, with the stator
currents on the KCL branch at `x` (`kcl_solve` from the converged `V0`), generic in the
number type of `x`.
"""
function rotor_energy(sys::DAESystem, x::AbstractVector, V0::AbstractVector; alpha::Real = 0.0)
    I = _stator_currents(sys, x, V0)
    e = zero(eltype(I))
    for (u, (rs, k)) in enumerate(_linear_rotors(sys))
        z = x[_state_range(sys, k)][rs.states]
        i = I[2u-1:2u]
        D = -(rs.A \ rs.Bs)
        w = z .- D * i
        e += dot(w, rs.Q * w) / 2 - alpha * dot(D * i, rs.Q * (D * i)) / 2
    end
    return e
end

"""
    shifted_kinetic_energy(sys, x) -> Real

`sum H (omega - 1)^2` over the synchronous machines (the model's shifted kinetic storage).
"""
function shifted_kinetic_energy(sys::DAESystem, x::AbstractVector)
    e = zero(eltype(x))
    for (k, c) in enumerate(sys.comps)
        component_role(c) === :generator && "omega" in state_names(c) && hasparam(c, "H") || continue
        j = findfirst(==("omega"), state_names(c))
        e += param_value(param_dict(c)["H"]) * (x[_state_range(sys, k)][j] - 1)^2
    end
    return e
end

"""
    quotient_hessian(m, S; eta0 = zeros, tol = 1e-9) -> NamedTuple

The Hessian of `eta -> S(lift(m, eta))` at `eta0` (the Hessian of the Bregman form of `S`
there), its eigenvalues and inertia (`positive`, `negative`, `zero`, counted against
`tol * max|eig|`). `S(x)` must accept dual numbers.
"""
function quotient_hessian(m::SectionModel, S; eta0::AbstractVector = zeros(neta(m)), tol::Real = 1e-9)
    Hs = ForwardDiff.hessian(e -> S(lift(m, e)), collect(float(eta0)))
    Hs = (Hs + Hs') / 2
    ev = eigvals(Symmetric(Hs))
    t = tol * max(maximum(abs, ev), floatmin())
    return (hessian = Hs, eigenvalues = ev, positive = count(>(t), ev), negative = count(<(-t), ev),
            zero = count(v -> abs(v) <= t, ev))
end

"""
    exchange_dissipation_forms(m; sys = m.sys, V0 = m.V0) -> NamedTuple

The quadratic forms at the reference point of `m` (see the comment at the top of
`src/ph/joint_storage.jl`): `rotor_loss`, `exchange` (symmetrised) and `network_loss`, all
on the section coordinates, and the comparison of the exchange with the rotor and network
dissipation: `null_exchange_max`, the largest eigenvalue of the exchange form on the null
space of the dissipation (positive: directions where no dissipation of these two kinds
dominates it, however small), and `ratio_max`, the largest generalized eigenvalue of the
exchange against the dissipation on its range (> 1: the exchange exceeds that dissipation
in some direction). Controllers carry no storage yet, so their dissipation is not part of
the comparison.
"""
function exchange_dissipation_forms(m::SectionModel; sys::DAESystem = m.sys,
                                    V0::AbstractVector = m.V0, tol::Real = 1e-9)
    n = neta(m)
    e0 = zeros(n)
    Vc = solve_network(sys, lift(m, e0), V0)
    units = _linear_rotors(sys)
    wfun(e) = (x = lift(m, e); I = _stator_currents(sys, x, Vc);
               vcat([x[_state_range(sys, k)][rs.states] .- (-(rs.A \ rs.Bs)) * I[2u-1:2u]
                     for (u, (rs, k)) in enumerate(units)]...))
    W = ForwardDiff.jacobian(wfun, e0)
    JI = ForwardDiff.jacobian(e -> _stator_currents(sys, lift(m, e), Vc), e0)
    JV = ForwardDiff.jacobian(e -> kcl_solve(sys, lift(m, e), Vc), e0)
    nu = length(units)
    L = zeros(4nu, 4nu)
    QD = zeros(4nu, 2nu)
    for (u, (rs, _)) in enumerate(units)
        L[4u-3:4u, 4u-3:4u] .= -(rs.Q * rs.A + rs.A' * rs.Q) / 2
        QD[4u-3:4u, 2u-1:2u] .= rs.Q * (-(rs.A \ rs.Bs))
    end
    A = section_jacobian(m, e0, Vc)
    nb = nbus(sys)
    Gv = zeros(2nb, 2nb)                      # V' G V with V = [Vd1, Vq1, ...]
    for i in 1:nb, j in 1:nb
        Gv[2i-1, 2j-1] = sys.G[i, j]
        Gv[2i, 2j] = sys.G[i, j]
    end
    loss = Symmetric(W' * L * W)
    ex = -(W' * QD * JI * A)
    ex = Symmetric((ex + ex') / 2)
    net = Symmetric(JV' * ((Gv + Gv') / 2) * JV)
    diss = Symmetric(Matrix(loss) + Matrix(net))
    F = eigen(diss)
    t = tol * max(maximum(abs, F.values), floatmin())
    rng = findall(>(t), F.values)
    nul = findall(<=(t), F.values)
    Nb = F.vectors[:, nul]
    null_ex = isempty(nul) ? -Inf : eigmax(Symmetric(Nb' * ex * Nb))
    Rb = F.vectors[:, rng] * Diagonal(1 ./ sqrt.(F.values[rng]))
    ratio = isempty(rng) ? NaN : eigmax(Symmetric(Rb' * ex * Rb))
    return (rotor_loss = loss, exchange = ex, network_loss = net,
            dissipation_rank = length(rng), null_exchange_max = null_ex, ratio_max = ratio,
            exchange_eig = extrema(eigvals(ex)), rotor_loss_eig_max = eigmax(loss),
            network_loss_eig_max = eigmax(net))
end

"""
    rotor_gradient_metric(rs) -> NamedTuple

The Nishino-Chakrabortty-Ishizaki structure for a linear rotor (`rotor_structure`): per axis
the symmetric `M` with `M B_s = -Psi''` (transposed) and `M A` symmetric, so that with an
internal EMF `E'' = (-psi_q'', psi_d'')` behind `j x''` on a lossless network, whose potential
has `grad_z U_net = Psi''' i`, the rotor dynamics are the gradient flow
`M z' = -grad(U_rot + U_net) + M B_f Efd` of the rotor energy `U_rot = -z' M A z / 2`. Returns
`M` (block-diagonal over the axes), `U = -sym(MA)` (the Hessian of `U_rot`), whether `M > 0`
and `U > 0` (the metric and the convexity), and the residual of `M B_s + Psi''`. The speed
voltage (`E'' = omega psi''` in the model) is outside this structure.
"""
function rotor_gradient_metric(rs)
    M = zeros(4, 4)
    for (idx, col, row) in ((1:2, 1, 1), (3:4, 2, 2))
        A, b, psi = rs.A[idx, idx], rs.Bs[idx, col], rs.Psi[row, idx]
        # unknowns (m11, m12, m22): M b = -psi, (MA)_12 = (MA)_21
        E = [b[1] b[2] 0.0; 0.0 b[1] b[2]; A[1, 2] (A[2, 2] - A[1, 1]) -A[2, 1]]
        m = E \ [-psi[1], -psi[2], 0.0]
        M[idx, idx] .= [m[1] m[2]; m[2] m[3]]
    end
    U = -(M * rs.A + (M * rs.A)') ./ 2
    return (M = M, U = U, metric_positive = isposdef(Symmetric(M)), convex = isposdef(Symmetric(U)),
            residual = maximum(abs, M * rs.Bs .+ rs.Psi'), asymmetry = maximum(abs, M * rs.A .- (M * rs.A)'))
end
