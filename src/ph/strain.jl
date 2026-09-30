# The joint strain-energy identity (TODO.md step 5b).
#
# Rotor. Per machine with a linear rotor in the stator-flux form (GENROU; GENSAL with its
# saturation read as a field input), z' = A z + B_f u_f + B_s i, (psi_d'', psi_q'') = Psi z, the
# strain metric M (`strain_metric`, per axis M B_s = -Psi', M A symmetric) makes
#   M z' = -grad U_rot + M B_f u_f - Psi' i,      U_rot = -z' M A z / 2.
# Network. With the internal EMF at nominal speed, E1 = (-psi_q'', psi_d'') (dq frame, rotated to
# the network frame by delta - pi/2) behind j x'' (ra = 0, xd'' = xq''),
#   U_net = U_B(V) + sum_k [|E1_k|^2/(2x'') - E1_k . V_k / x''] + sum_L U_L(V_L, z_L),
# U_B the susceptance potential (`network_potential`, with the Norton -1/x'' on the diagonal)
# and U_L the ComplexLoads' reactive potential, int (Q_act(v)/v^2 - Q0/V0^2) v dv (the
# constant-impedance base Q0/V0^2 is in B). A load with the voltage filter (state z,
# t1 z' = v - s, s = Vini + z) has Q_act = Q_act(z) and U_L = Q_act ln v - Q0 v^2/(2 V0^2) - Phi(z),
# Phi = Q_act (ln s - 1/kqv), dPhi/dz = ln s dQ_act/dz. Then grad over psi'' of the machine term is the
# nominal-speed current i1 = (E1 - v)/(j x''), and grad over delta_k is Te_k = Pe_k/omega_k.
#
# Along the flow (V' from the KCL branch) U_ext = U_rot + U_net has, exactly,
#   dU_ext/dt = -sum z'Mz' + sum Efd B_f'Mz' + sum Te delta'              (the paper's identity)
#             + sum (u_f - Efd) B_f'Mz'                                     (GENSAL saturation)
#             + sum (i1 - i) . Psi z' + sum_k (omega_k - 1) E1_k . V_k' / x''  (speed voltage)
#             + sum_i [-j (G V)_i] . V_i'                                   (conductances)
#             + sum_L [-j Gp_L V_L] . V_L'                                  (ComplexLoad active part)
#             + sum_i [(dQ_i - j dP_i) V_i] . V_i'                          (frequency-dependent loads)
#             + sum_L (ln v_L - ln s_L) dQ_act/dz z_L'                      (ComplexLoad lag)
# (complex numbers read as 2-vectors, `.` the real dot product). The lag term has the sign of
# Q0 (ln is increasing and t1 z' = v - s): it creates energy at every inductive load. The first line is the
# gradient structure of Nishino, Chakrabortty and Ishizaki; the others are what the detailed
# model adds. On `lossless_variant` the conductance and active-load terms are zero.
# `strain_balance` evaluates every term and the identity error.

"""
    strain_rotor(c) -> NamedTuple or nothing

The rotor of a machine in the form the strain-energy structure needs:
`z' = A z + B_f u_f + B_s [id, iq]` with `(psi_d'', psi_q'') = Psi z`, the axes' state blocks,
and the subtransient reactance. GENROU: its `rotor_structure` (`u_f = Efd`). GENSAL: the
linear part of its rotor; its additive saturation enters as `u_f = Efd - Sat (xd - xl)`.
`nothing` for other components.
"""
strain_rotor(::AbstractComponent) = nothing
function strain_rotor(c::GENROU_PHTRUE)
    rs = rotor_structure(c)
    return (states = rs.states, A = rs.A, Bf = rs.Bf, Bs = rs.Bs, Psi = rs.Psi,
            axes = (1:2, 3:4), ra = c.p.ra, xd_pp = c.p.xd_double_prime, xq_pp = c.p.xq_double_prime)
end
function strain_rotor(c::GENSAL_PHTRUE)
    p = c.p
    kd = (p.xd_double_prime - p.xl) / (p.xd_prime - p.xl)
    A = [-1/p.Td0_prime 0.0 0.0;
         1/p.Td0_double_prime -1/p.Td0_double_prime 0.0;
         0.0 0.0 -1/p.Tq0_double_prime]
    Bf = [1/p.Td0_prime, 0.0, 0.0]
    Bs = [-(p.xd - p.xd_prime)/p.Td0_prime 0.0;
          -(p.xd_prime - p.xl)/p.Td0_double_prime 0.0;
          0.0 -(p.xq - p.xq_double_prime)/p.Tq0_double_prime]
    Psi = [kd 1.0-kd 0.0; 0.0 0.0 1.0]
    return (states = 3:5, A = A, Bf = Bf, Bs = Bs, Psi = Psi, axes = (1:2, 3:3),
            ra = p.ra, xd_pp = p.xd_double_prime, xq_pp = p.xq_double_prime)
end

# per axis: the symmetric m with m b = -psi and (m A) symmetric (one- or two-state axis)
function _axis_metric(A::AbstractMatrix, b::AbstractVector, psi::AbstractVector)
    length(b) == 1 && return fill(-psi[1] / b[1], 1, 1)
    # unknowns (m11, m12, m22): M b = -psi, (MA)_12 = (MA)_21
    E = [b[1] b[2] 0.0; 0.0 b[1] b[2]; A[1, 2] (A[2, 2] - A[1, 1]) -A[2, 1]]
    m = E \ [-psi[1], -psi[2], 0.0]
    return [m[1] m[2]; m[2] m[3]]
end

"""
    strain_metric(sr) -> NamedTuple

For a `strain_rotor`: the block-diagonal metric `M` (per axis `M B_s = -Psi'`, `M A`
symmetric), the Hessian `U = -sym(MA)` of `U_rot = -z' M A z / 2`, whether `M > 0` and
`U > 0`, and the residuals of the two conditions. `rotor_gradient_metric` is the GENROU case.
"""
function strain_metric(sr)
    n = length(sr.states)
    M = zeros(n, n)
    for (col, idx) in enumerate(sr.axes)
        M[idx, idx] .= _axis_metric(sr.A[idx, idx], sr.Bs[idx, col], sr.Psi[col, idx])
    end
    U = -(M * sr.A + (M * sr.A)') ./ 2
    return (M = M, U = U, metric_positive = isposdef(Symmetric(M)), convex = isposdef(Symmetric(U)),
            residual = maximum(abs, M * sr.Bs .+ sr.Psi'), asymmetry = maximum(abs, M * sr.A .- (M * sr.A)'))
end

# the strain rotors of sys: (component index, strain_rotor, strain_metric, bus index)
function _strain_units(sys::DAESystem)
    out = []
    for (k, c) in enumerate(sys.comps)
        sr = strain_rotor(c)
        sr === nothing && continue
        (sr.ra == 0.0 && sr.xd_pp == sr.xq_pp) ||
            throw(ArgumentError("$(name(c)): the internal-EMF potential needs ra = 0 and xd'' = xq''"))
        push!(out, (k = k, sr = sr, g = strain_metric(sr), bus = sys.inj[k][1]))
    end
    return out
end

# dq (machine) -> network frame for delta: multiply by exp(j (delta - pi/2))
_dq_to_net(a, b, delta) = (a * sin(delta) + b * cos(delta), -a * cos(delta) + b * sin(delta))

# nominal-speed internal EMF E1 = (-psi_q'', psi_d'') of a strain unit, in the network frame
function _nominal_emf(sr, xk)
    ps = sr.Psi * xk[sr.states]
    return _dq_to_net(-ps[2], ps[1], xk[1])
end

# ComplexLoad reactive quantities at voltage magnitude v: (Q_act, dQ_act/dz, Gp, U_L, lag),
# lag = (ln v - ln s) dQ_act/dz with s = Vini + z the filtered voltage (0 without the filter);
# the potential needs the load law's middle band (k(V) = 1). With the filter, U_L includes
# -Phi(z), Phi = Q_act (ln s - 1/kqv), so its z-rate is lag z' instead of ln(v) dQ/dz z'.
function _complexload_reactive(c::COMPLEXLOAD, xk, v)
    p = c.p
    p.has_band && !(p.udmin_c < v < p.udmax_c) &&
        throw(DomainError(v, "$(name(c)): |V| outside (udmin, udmax); the reactive potential covers k(V) = 1 only"))
    if p.t1 > 0.0
        vr = 1.0 + xk[1] / p.Vini_c
        vr > 1.0e-6 || throw(DomainError(vr, "$(name(c)): voltage-ratio guard active"))
        Q = p.Q0_c * vr^p.kqv_c
        dQ = p.Q0_c * p.kqv_c * vr^(p.kqv_c - 1) / p.Vini_c
        lns = log(p.Vini_c * vr)
        Phi = p.kqv_c == 0.0 ? zero(Q) : Q * (lns - 1 / p.kqv_c)
        UL = Q * log(v) - p.Q_base_c * v^2 / 2 - Phi
        lag = (log(v) - lns) * dQ
        P = p.P0_c * vr^p.kpv_c
    else
        vr = v / p.Vini_c
        Q = p.Q0_c * vr^p.kqv_c
        UL = (p.kqv_c == 0.0 ? p.Q0_c * log(v) : p.Q0_c * vr^p.kqv_c / p.kqv_c) - p.Q_base_c * v^2 / 2
        lag = zero(Q)
        P = p.P0_c * vr^p.kpv_c
    end
    return Q, lag, P / v^2 - p.G_load_c, UL
end

"""
    strain_energy(sys, x, V) -> Real

`U_ext = U_rot + U_net` (see the comment at the top of `src/ph/strain.jl`) at the states `x`
and bus voltages `V`, generic in the number type.
"""
function strain_energy(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    u = network_potential(sys, V)
    for un in _strain_units(sys)
        xk = x[_state_range(sys, un.k)]
        z = xk[un.sr.states]
        u -= dot(z, un.g.M * (un.sr.A * z)) / 2
        e = _nominal_emf(un.sr, xk)
        b = un.bus
        u += (e[1]^2 + e[2]^2) / (2 * un.sr.xd_pp) - (e[1] * V[2b-1] + e[2] * V[2b]) / un.sr.xd_pp
    end
    for (k, c) in enumerate(sys.comps)
        b = sys.inj[k][1]
        if c isa COMPLEXLOAD
            u += _complexload_reactive(c, x[_state_range(sys, k)], sqrt(V[2b-1]^2 + V[2b]^2))[4]
        elseif c isa ConstantPowerSink               # integrable_loss_variant
            u += sink_potential(c, V[2b-1], V[2b])
        end
    end
    return u
end

"""
    strain_balance(sys, x, V) -> Dict

At `(x, V)` (KCL satisfied, healthy network), along the flow with `V'` from the KCL branch:
`rate = dU_ext/dt` (by differentiation) and the terms of its decomposition, `dissipation`
(`-sum z'Mz'`), `field` (`sum Efd B_f'Mz'`), `angle` (`sum Te delta'`), `saturation`,
`speed_rotor`, `speed_network`, `conductance`, `active_load`, `frequency_load`,
`load_lag`; `paper_residual = rate - (dissipation + field + angle)`, `identity_error =
rate - sum of all terms`, `torque_error` (`max |Te - Pe/omega|`, `Te = dU_net/d delta`),
`structure_error` (`max |z' - A z - B_s i - B_f u_f|`) and per machine `speed_rotor`,
`dissipation`, `field`, `saturation`.
"""
function strain_balance(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    f, _ = dae_residual(sys, x, V)
    Vdot = ForwardDiff.derivative(t -> kcl_solve(sys, x .+ t .* f, V), 0.0)
    rate = ForwardDiff.derivative(t -> (xt = x .+ t .* f; strain_energy(sys, xt, kcl_solve(sys, xt, V))), 0.0)
    ins, _ = component_io(sys, x, V)
    diss = field = angle = sat = spr = spn = 0.0
    terr = serr = 0.0
    machines = Dict{String,Any}()
    for un in _strain_units(sys)
        c, sr, M = sys.comps[un.k], un.sr, un.g.M
        r = _state_range(sys, un.k)
        xk, fk = x[r], f[r]
        om = xk[2]
        z, dz = xk[sr.states], fk[sr.states]
        ps, dps = sr.Psi * z, sr.Psi * dz
        b = un.bus
        Vb, dVb = [V[2b-1], V[2b]], [Vdot[2b-1], Vdot[2b]]
        sd, cd = sin(xk[1]), cos(xk[1])
        vd, vq = Vb[1] * sd - Vb[2] * cd, Vb[1] * cd + Vb[2] * sd
        x2 = sr.xd_pp
        i = [(om * ps[1] - vq) / x2, (vd + om * ps[2]) / x2]       # the model's stator currents
        i1 = [(ps[1] - vq) / x2, (vd + ps[2]) / x2]                # at nominal speed
        Efd = ins[un.k][findfirst(==("Efd"), input_names(c))]
        res = dz .- sr.A * z .- sr.Bs * i
        uf = dot(sr.Bf, res) / dot(sr.Bf, sr.Bf)
        serr = max(serr, maximum(abs, res .- sr.Bf .* uf))
        Mdz = M * dz
        d_k = -dot(dz, Mdz)
        fl_k = Efd * dot(sr.Bf, Mdz)
        sa_k = (uf - Efd) * dot(sr.Bf, Mdz)
        sp_k = dot(i1 .- i, dps)
        e = _nominal_emf(sr, xk)
        Te = (e[2] * Vb[1] - e[1] * Vb[2]) / x2                    # Im(E1 conj V) / x''
        terr = max(terr, abs(Te - (vd * i[1] + vq * i[2]) / om))
        diss += d_k; field += fl_k; sat += sa_k; spr += sp_k
        angle += Te * fk[1]
        spn += (om - 1) * dot(e, dVb) / x2
        machines[name(c)] = Dict("dissipation" => d_k, "field" => fl_k, "saturation" => sa_k,
                                 "speed_rotor" => sp_k, "Te" => Te)
    end
    # network terms: conductances, ComplexLoad active part and filter, frequency loads
    nb = nbus(sys)
    rot(a, b) = (b, -a)                                             # -j (a + j b)
    cond = act = freq = filt = 0.0
    for i in 1:nb
        sys.slack[i] && continue
        gd = gq = 0.0
        for j in 1:nb
            gd += sys.G[i, j] * V[2j-1]; gq += sys.G[i, j] * V[2j]
        end
        w = rot(gd, gq)
        cond += w[1] * Vdot[2i-1] + w[2] * Vdot[2i]
    end
    coi = _coi_omega(sys, x)
    la = sys.load
    for i in 1:nb
        (la.kpf[i] != 0.0 || la.kqf[i] != 0.0) && !sys.slack[i] || continue
        dP = la.kpf[i] * (coi - 1) * la.G[i]
        dQ = la.kqf[i] * (coi - 1) * la.B[i]
        vd, vq = V[2i-1], V[2i]
        # (dQ - j dP) V
        freq += (dQ * vd + dP * vq) * Vdot[2i-1] + (dQ * vq - dP * vd) * Vdot[2i]
    end
    for (k, c) in enumerate(sys.comps)
        c isa COMPLEXLOAD || continue
        b = sys.inj[k][1]
        xk = x[_state_range(sys, k)]
        v = sqrt(V[2b-1]^2 + V[2b]^2)
        _, lag, Gp, _ = _complexload_reactive(c, xk, v)
        w = rot(Gp * V[2b-1], Gp * V[2b])
        act += w[1] * Vdot[2b-1] + w[2] * Vdot[2b]
        c.p.t1 > 0.0 && (filt += lag * f[_state_range(sys, k)][1])
    end
    total = diss + field + angle + sat + spr + spn + cond + act + freq + filt
    return Dict{String,Any}("rate" => rate, "dissipation" => diss, "field" => field, "angle" => angle,
                            "saturation" => sat, "speed_rotor" => spr, "speed_network" => spn,
                            "conductance" => cond, "active_load" => act, "frequency_load" => freq,
                            "load_lag" => filt, "paper_residual" => rate - (diss + field + angle),
                            "identity_error" => rate - total, "torque_error" => terr,
                            "structure_error" => serr, "machines" => machines)
end

# the centre-of-inertia speed as the DAE residual computes it
function _coi_omega(sys::DAESystem, x::AbstractVector)
    m = sys.coi_members
    isempty(m) && return one(eltype(x))
    length(m) == 1 && return x[sys.offsets[m[1]] + 1]
    return sum(sys.coi_weights[q] * x[sys.offsets[m[q]] + 1] for q in eachindex(m)) / sys.coi_total
end

# -- the supply-shifted candidate and its local decay forms --------------------------------
#
# Around an equilibrium (omega* = 1, Tm* = Te*), with the kinetic storage
# K = sum omega_b H (omega - 1)^2, torque-consistent for D = 0:
#   d/dt K + sum Te delta' = sum omega_b (omega - 1) Tm/omega - omega_b (omega_coi - 1) sum Te.
# The candidate
#   S' = U_ext + K - sum u_f* B_f'Mz - sum Tm* delta - c*' V,
# c* = -j (G V* + Gp* V*) (the equilibrium conductance and ComplexLoad active currents,
# rotated), removes every exact first-order supply, so its rate is a sum of purely
# second-order terms:
#   dissipation  -z'Mz'                             damping      -omega_b sum (Tm*/omega + D) (omega - 1)^2
#   coi          -omega_b (omega_coi - 1) sum dTe   conductance  [-j G dV] . V'
#   active_load  [-j d(Gp V)] . V'                  speed        (i1 - i) . Psi z' + (omega - 1) E1 . V' / x''
#   load_lag     (ln v - ln s) dQ/dz z'             field        du_f B_f'Mz' (exciters; GENSAL saturation)
#   governor     omega_b dTm (omega - 1)/omega
# `strain_decay_forms` returns the quadratic form of each term at the equilibrium (section
# coordinates, symmetrised) and the Hessian of S'; their sum is sym(Hess S' A).

"""
    strain_reference(sys, x, V) -> NamedTuple

The equilibrium constants of `shifted_strain_energy`: per strain unit `u_f*` and `Tm*`, and
the rotated conductance and ComplexLoad active currents `c*` (stacked like `V`).
"""
function strain_reference(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    fac = _strain_factors(sys, x, V)
    nb = nbus(sys)
    c = zeros(2nb)
    for i in 1:nb
        gd = gq = 0.0
        for j in 1:nb
            gd += sys.G[i, j] * V[2j-1]; gq += sys.G[i, j] * V[2j]
        end
        c[2i-1] += gq; c[2i] -= gd                                  # -j (G V)_i
    end
    for (l, (k, b)) in enumerate(fac.loads)
        c[2b-1] += fac.GpV[2l]; c[2b] -= fac.GpV[2l-1]             # -j Gp V
    end
    return (uf = fac.uf, Tm = fac.Tm, c = c)
end

"""
    shifted_strain_energy(sys, x, V, ref) -> Real

`S' = U_ext + K - sum u_f* B_f'Mz - sum Tm* delta - c*' V` (see above), `ref` from
`strain_reference` at the equilibrium; generic in the number type.
"""
function shifted_strain_energy(sys::DAESystem, x::AbstractVector, V::AbstractVector, ref)
    s = strain_energy(sys, x, V) - dot(ref.c, V)
    for (u, un) in enumerate(_strain_units(sys))
        xk = x[_state_range(sys, un.k)]
        z = xk[un.sr.states]
        H = param_value(param_dict(sys.comps[un.k])["H"])
        s += sys.omega_b * H * (xk[2] - 1)^2 - ref.uf[u] * dot(un.sr.Bf, un.g.M * z) - ref.Tm[u] * xk[1]
    end
    return s
end

# generic per-unit and per-load factors of the second-order terms at (x, V)
function _strain_factors(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    T = promote_type(eltype(x), eltype(V))
    f, _ = dae_residual(sys, x, V)
    ins, _ = component_io(sys, x, V)
    Te, di, uf, Tm = T[], T[], T[], T[]
    for un in _strain_units(sys)
        c, sr = sys.comps[un.k], un.sr
        r = _state_range(sys, un.k)
        xk = x[r]
        om = xk[2]
        z = xk[sr.states]
        ps = sr.Psi * z
        b = un.bus
        sd, cd = sin(xk[1]), cos(xk[1])
        vd, vq = V[2b-1] * sd - V[2b] * cd, V[2b-1] * cd + V[2b] * sd
        x2 = sr.xd_pp
        i = [(om * ps[1] - vq) / x2, (vd + om * ps[2]) / x2]
        i1 = [(ps[1] - vq) / x2, (vd + ps[2]) / x2]
        e = _nominal_emf(sr, xk)
        push!(Te, (e[2] * V[2b-1] - e[1] * V[2b]) / x2)
        append!(di, i1 .- i)
        res = f[r][sr.states] .- sr.A * z .- sr.Bs * i
        push!(uf, dot(sr.Bf, res) / dot(sr.Bf, sr.Bf))
        push!(Tm, ins[un.k][findfirst(==("Tm"), input_names(c))])
    end
    loads = Tuple{Int,Int}[]
    GpV, lagarg = T[], T[]
    for (k, c) in enumerate(sys.comps)
        c isa COMPLEXLOAD || continue
        b = sys.inj[k][1]
        push!(loads, (k, b))
        xk = x[_state_range(sys, k)]
        v = sqrt(V[2b-1]^2 + V[2b]^2)
        _, _, Gp, _ = _complexload_reactive(c, xk, v)
        append!(GpV, (Gp * V[2b-1], Gp * V[2b]))
        push!(lagarg, c.p.t1 > 0.0 ? log(v) - log(c.p.Vini_c + xk[1]) : zero(T))
    end
    return (Te = Te, di = di, uf = uf, Tm = Tm, GpV = GpV, lagarg = lagarg, loads = loads)
end

"""
    strain_decay_forms(m; V0 = m.V0) -> NamedTuple

At the reference point of the section model `m` (an equilibrium with omega* = 1, no
frequency-dependent loads), the quadratic forms (section coordinates, symmetrised) of the
second-order terms of dS'/dt listed above (`forms`, a Dict), the Hessian `H` of `S'`, the
norm of its gradient (0 at the equilibrium), the section Jacobian `A`, and `closure`, the
largest entry of `sum(forms) - sym(H A)` relative to `sym(H A)`.
"""
function strain_decay_forms(m::SectionModel; V0::AbstractVector = m.V0)
    sys = m.sys
    (any(!=(0), sys.load.kpf) || any(!=(0), sys.load.kqf)) &&
        throw(ArgumentError("frequency-dependent loads are not covered by strain_decay_forms"))
    n = neta(m)
    e0 = zeros(n)
    x0 = lift(m, e0)
    Vc = solve_network(sys, x0, V0)
    units = _strain_units(sys)
    maximum(abs(x0[_state_range(sys, un.k)][2] - 1) for un in units) < 1e-9 ||
        throw(ArgumentError("strain_decay_forms needs omega* = 1"))
    ref = strain_reference(sys, x0, Vc)
    Sp(e) = (x = lift(m, e); shifted_strain_energy(sys, x, kcl_solve(sys, x, Vc), ref))
    grad = ForwardDiff.gradient(Sp, e0)
    H = ForwardDiff.hessian(Sp, e0)
    H = (H + H') / 2
    Lx = ForwardDiff.jacobian(e -> lift(m, e), e0)                  # state deviations
    Jf = ForwardDiff.jacobian(e -> (x = lift(m, e); dae_residual(sys, x, kcl_solve(sys, x, Vc))[1]), e0)
    Vx = ForwardDiff.jacobian(x -> kcl_solve(sys, x, Vc), x0)
    JV = Vx * Lx
    JVd = Vx * Jf                                                   # V' at linear order
    Jc = ForwardDiff.jacobian(e -> (x = lift(m, e); fa = _strain_factors(sys, x, kcl_solve(sys, x, Vc));
                                    vcat(fa.Te, fa.di, fa.uf, fa.Tm, fa.GpV, fa.lagarg)), e0)
    fac = _strain_factors(sys, x0, Vc)
    nu, nl = length(fac.Te), length(fac.loads)
    JTe = Jc[1:nu, :]; Jdi = Jc[nu+1:3nu, :]; Juf = Jc[3nu+1:4nu, :]; JTm = Jc[4nu+1:5nu, :]
    JGpV = Jc[5nu+1:5nu+2nl, :]; Jlag = Jc[5nu+2nl+1:end, :]
    sy(K) = (K + K') / 2
    F = Dict(k => zeros(n, n) for k in ("dissipation", "damping", "coi", "conductance", "active_load",
                                       "speed", "load_lag", "field", "governor"))
    wb = sys.omega_b
    lcoi = ForwardDiff.gradient(e -> _coi_omega(sys, lift(m, e)), e0)
    F["coi"] .= -wb .* sy(lcoi * vec(sum(JTe; dims = 1))')
    for (u, un) in enumerate(units)
        r = _state_range(sys, un.k)
        zr = r[un.sr.states]
        M, Bf = un.g.M, un.sr.Bf
        Jz = Jf[zr, :]
        lw = Lx[r[2], :]
        D = param_value(param_dict(sys.comps[un.k])["D"])
        F["dissipation"] .-= Jz' * M * Jz
        F["damping"] .-= wb * (fac.Tm[u] + D) .* (lw * lw')
        b = un.bus
        e = _nominal_emf(un.sr, x0[r])
        F["speed"] .+= sy(Jdi[2u-1:2u, :]' * un.sr.Psi * Jz) .+
                       sy(lw * (e[1] .* JVd[2b-1, :] .+ e[2] .* JVd[2b, :])') ./ un.sr.xd_pp
        F["field"] .+= sy(Juf[u, :] * (Bf' * M * Jz))
        F["governor"] .+= wb .* sy(JTm[u, :] * lw')
    end
    Gh = kron(sys.G, [0.0 1.0; -1.0 0.0])
    F["conductance"] .= sy(JV' * Gh' * JVd)
    Rj = [0.0 1.0; -1.0 0.0]
    for (l, (k, b)) in enumerate(fac.loads)
        F["active_load"] .+= sy((Rj * JGpV[2l-1:2l, :])' * JVd[2b-1:2b, :])
        c = sys.comps[k]
        c.p.t1 > 0.0 || continue
        p = c.p
        vr = 1.0 + x0[_state_range(sys, k)][1] / p.Vini_c
        dQ = p.Q0_c * p.kqv_c * vr^(p.kqv_c - 1) / p.Vini_c
        F["load_lag"] .+= dQ .* sy(Jlag[l, :] * Jf[_state_range(sys, k)[1], :]')
    end
    A = section_jacobian(m, e0, Vc)
    DHA = sy(H * A)
    return (forms = F, H = H, gradient_norm = norm(grad), A = A,
            closure = maximum(abs, sum(values(F)) .- DHA) / maximum(abs, DHA))
end
