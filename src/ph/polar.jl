# The polar network balance (TODO.md step 5, the reviewer's step 1).
#
# Network. For the susceptance part B of the Y-bus (symmetric), the scalar
#   U_B(V) = -1/2 sum_ij B_ij (Vd_i Vd_j + Vq_i Vq_j) = -1/2 sum_ij B_ij |V_i||V_j| cos(theta_i - theta_j)
# satisfies, for any voltage path,
#   dU_B/dt = sum_i (P_i^B theta_i' + Q_i^B d ln|V_i|/dt),   P^B + j Q^B = V conj(j B V),
# the power the susceptances absorb at each bus: the conjugate pairs follow from the identity,
# (P, theta) and (Q, ln|V|). The conductance part contributes
#   sum_i (P_i^G theta_i' + Q_i^G d ln|V_i|/dt) = -Im(V'^H G V),
# a one-form whose curl in Cartesian coordinates is 2 (G kron J), J = [0 1; -1 0]: nonzero for
# every transfer or shunt conductance, so no scalar network potential covers it. KCL makes
# the network's absorbed power equal what the components inject.
#
# Machines. A GENROU with ra = 0 and xd'' = xq'' is an internal EMF E'' = omega (-psi_q'', psi_d'')
# (dq frame) behind j x''; at that internal node the network delivers
#   S_int = E'' conj(I) = P_int + j Q_int,   I = id + j iq,
# and the polar supply it offers the rotor magnetics is (P_int phi' + Q_int d ln|E''|/dt)/omega_b,
# phi the angle of E'' in the dq frame; algebraically this is Im(conj(I) dE''/dt)/omega_b =
# omega i . dpsi''/dt/omega_b + omega' i . psi''/omega_b (the transformer power). A storage built
# from the network potential and the machine energies closes exactly when the rotor's own
# exchange -w'QD di/dt equals minus that supply; `polar_balance` measures both.

"""
    network_potential(sys, V) -> Real

`U_B = -1/2 sum B_ij (Vd_i Vd_j + Vq_i Vq_j)` for the susceptance part of the Y-bus.
"""
function network_potential(sys::DAESystem, V::AbstractVector)
    nb = nbus(sys)
    u = zero(eltype(V))
    for i in 1:nb, j in 1:nb
        u -= sys.B[i, j] * (V[2i-1] * V[2j-1] + V[2i] * V[2j]) / 2
    end
    return u
end

# power absorbed at each bus by the susceptances (B) and the conductances (G) of the Y-bus
function _bus_absorbed(sys::DAESystem, V::AbstractVector)
    nb = nbus(sys)
    T = eltype(V)
    PB, QB, PG, QG = zeros(T, nb), zeros(T, nb), zeros(T, nb), zeros(T, nb)
    for i in 1:nb
        bd = bq = gd = gq = zero(T)                       # (B V)_i and (G V)_i, complex parts
        for j in 1:nb
            bd += sys.B[i, j] * V[2j-1]; bq += sys.B[i, j] * V[2j]
            gd += sys.G[i, j] * V[2j-1]; gq += sys.G[i, j] * V[2j]
        end
        vd, vq = V[2i-1], V[2i]
        # I = j B V = (-bq) + j bd;  S = V conj(I)
        PB[i] = vd * (-bq) + vq * bd
        QB[i] = vq * (-bq) - vd * bd
        # I = G V = gd + j gq
        PG[i] = vd * gd + vq * gq
        QG[i] = vq * gd - vd * gq
    end
    return PB, QB, PG, QG
end

"""
    polar_balance(sys, x, V; faults_on = false) -> Dict

At `(x, V)` (KCL satisfied), along the flow `x' = f(x, V)` with `V'` from the KCL branch:
- `lossless_identity_error`: `dU_B/dt` against `sum (P^B theta' + Q^B d ln|V|/dt)`;
- `conductance_rate`: `sum (P^G theta' + Q^G d ln|V|/dt)` and its closed form `-Im(V'^H G V)`;
- per machine with a linear rotor: `P_int`, `Q_int`, the polar supply
  `(P_int phi' + Q_int d ln|E''|/dt)/omega_b`, the transformer power `i . dpsi''/dt / omega_b`,
  the rotor exchange `w'QD di/dt`, and `mismatch = exchange - polar supply`.
"""
function polar_balance(sys::DAESystem, x::AbstractVector, V::AbstractVector; faults_on::Bool = false)
    f, _ = dae_residual(sys, x, V; faults_on)
    Vdot = ForwardDiff.derivative(t -> kcl_solve(sys, x .+ t .* f, V; faults_on), 0.0)
    Udot = ForwardDiff.derivative(t -> network_potential(sys, V .+ t .* Vdot), 0.0)
    PB, QB, PG, QG = _bus_absorbed(sys, V)
    nb = nbus(sys)
    polarB = polarG = 0.0
    imG = 0.0
    for i in 1:nb
        vd, vq, dd, dq = V[2i-1], V[2i], Vdot[2i-1], Vdot[2i]
        m2 = vd^2 + vq^2
        thdot = (vd * dq - vq * dd) / m2
        dlnv = (vd * dd + vq * dq) / m2
        polarB += PB[i] * thdot + QB[i] * dlnv
        polarG += PG[i] * thdot + QG[i] * dlnv
        for j in 1:nb   # -Im(conj(Vdot_i) G_ij V_j)
            imG -= sys.G[i, j] * (dd * V[2j] - dq * V[2j-1])
        end
    end
    I = _stator_currents(sys, x, V; faults_on)
    Idot = ForwardDiff.derivative(t -> _stator_currents(sys, x .+ t .* f, V; faults_on), 0.0)
    machines = Dict{String,Any}()
    for (u, (rs, k)) in enumerate(_linear_rotors(sys))
        c = sys.comps[k]
        r = _state_range(sys, k)
        xk, fk = x[r], f[r]
        om, dom = xk[2], fk[2]
        z, dz = xk[rs.states], fk[rs.states]
        ps, dps = rs.Psi * z, rs.Psi * dz
        E = om .* [-ps[2], ps[1]]
        dE = dom .* [-ps[2], ps[1]] .+ om .* [-dps[2], dps[1]]
        i, di = I[2u-1:2u], Idot[2u-1:2u]
        P = E[1] * i[1] + E[2] * i[2]
        Q = E[2] * i[1] - E[1] * i[2]
        E2 = dot(E, E)
        phidot = (E[1] * dE[2] - E[2] * dE[1]) / E2
        dlnE = dot(E, dE) / E2
        wb = params(c).omega_b
        D = -(rs.A \ rs.Bs)
        w = z .- D * i
        supply = (P * phidot + Q * dlnE) / wb
        exch = dot(w, rs.Q * (D * di))
        machines[name(c)] = Dict("P_int" => P, "Q_int" => Q, "polar_supply" => supply,
                                 "transformer_power" => dot(i, dps) / wb,
                                 "supply_vs_transformer" => supply - (om * dot(i, dps) + dom * dot(i, ps)) / wb,
                                 "exchange" => exch, "mismatch" => exch - supply)
    end
    return Dict{String,Any}("Udot" => Udot, "lossless_polar_rate" => polarB,
                            "lossless_identity_error" => Udot - polarB,
                            "conductance_rate" => polarG, "conductance_closed_form_error" => polarG - imG,
                            "machines" => machines)
end

"""
    conductance_curl(sys) -> Matrix

The curl of the conductance one-form in the stacked Cartesian voltages `[Vd_1, Vq_1, ...]`:
`2 (G kron J)`, `J = [0 1; -1 0]`; it vanishes only without conductances.
"""
conductance_curl(sys::DAESystem) = 2 .* kron(sys.G, [0.0 1.0; -1.0 0.0])
