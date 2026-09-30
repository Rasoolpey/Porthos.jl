# The integrable-loss equivalent of a system (TODO.md, the controller-storage stage on a
# network whose every supply is exact).
#
# The plain `lossless_variant` has no equilibrium at a case's dispatch: without active
# consumption the machines' electrical power must sum to zero while the mechanical power does
# not. The integrable-loss equivalent keeps the equilibrium (x*, V*): it removes the Y-bus
# conductances and the ComplexLoads' active part, and at every bus draws the complex power
# they absorbed at (x*, V*), S_i* = V_i* conj((G V*)_i + Gp_i* V_i*), through a constant-power
# sink. A constant-power sink has the exact potential
#   U_S = P theta + Q ln|V|,   grad_V U_S = -j I_S,
# so on this network every supply is integrable; what stays non-exact is the speed voltage
# and the ComplexLoad lag. It is a certificate-only variant for staging the storage
# construction: the case's model and its losses are unchanged.

struct SinkParams
    P::Float64
    Q::Float64
    Vd0::Float64          # the equilibrium phasor, the zero of the sink's angle
    Vq0::Float64
end

"""
    ConstantPowerSink

A stateless certificate-only element drawing the constant complex power `P + jQ` at its
bus: injected current `-(P - jQ) V / |V|^2`. Built by `integrable_loss_variant`.
"""
struct ConstantPowerSink <: AbstractComponent
    name::String
    bus::Int
    p::SinkParams
    params::ParamDict
end

model_type(::ConstantPowerSink) = "CONSTANT_POWER_SINK"
component_role(::ConstantPowerSink) = :load
state_names(::ConstantPowerSink) = String[]
input_names(::ConstantPowerSink) = ["Vd", "Vq"]
output_names(::ConstantPowerSink) = ["Id", "Iq"]
bus(c::ConstantPowerSink) = c.bus
certificate_only(::ConstantPowerSink) = true
# the current turns with V; P and Q are invariant
rotation_action(::ConstantPowerSink) = (angles = String[], inputs = [("Vd", "Vq")], outputs = [("Id", "Iq")])

@inline function _sink_current(p::SinkParams, Vd, Vq)
    m2 = Vd * Vd + Vq * Vq
    return -(p.P * Vd + p.Q * Vq) / m2, -(p.P * Vq - p.Q * Vd) / m2
end

function _outputs!(y, c::ConstantPowerSink, x, u, p::SinkParams, rec::ModeRecorder)
    y[1], y[2] = _sink_current(p, u[1], u[2])
    return y
end

function _step!(dx, y, c::ConstantPowerSink, x, u, p::SinkParams, rec::ModeRecorder)
    y === nothing || _outputs!(y, c, x, u, p, rec)
    return dx
end

hamiltonian(c::ConstantPowerSink, x, p::SinkParams) = 0.0
grad_hamiltonian!(g, c::ConstantPowerSink, x, p::SinkParams) = g
injection(c::ConstantPowerSink, x, V, p::SinkParams = c.p) = _sink_current(p, V[1], V[2])

"""
    sink_potential(c, Vd, Vq) -> Real

`P theta + Q ln|V|`, `theta` the angle of `V` from the sink's equilibrium phasor.
"""
sink_potential(c::ConstantPowerSink, Vd, Vq) =
    c.p.P * atan(Vq * c.p.Vd0 - Vd * c.p.Vq0, Vd * c.p.Vd0 + Vq * c.p.Vq0) + c.p.Q * log(Vd * Vd + Vq * Vq) / 2

"""
    integrable_loss_variant(sys, x, V; tol = 1e-12) -> DAESystem

The integrable-loss equivalent of `sys` at its equilibrium `(x, V)` (see the comment at the
top of `src/ph/sink.jl`): `lossless_variant(sys)` plus one `ConstantPowerSink` per bus whose
absorbed power exceeds `tol`. `(x, V)` is an equilibrium of the result as well.
"""
function integrable_loss_variant(sys::DAESystem, x::AbstractVector, V::AbstractVector; tol::Real = 1e-12)
    L = lossless_variant(sys)
    nb = nbus(sys)
    Iab = zeros(2nb)                                   # current absorbed by the removed elements
    for i in 1:nb, j in 1:nb
        Iab[2i-1] += sys.G[i, j] * V[2j-1]
        Iab[2i] += sys.G[i, j] * V[2j]
    end
    for (k, c) in enumerate(sys.comps)
        c isa COMPLEXLOAD || continue
        b = sys.inj[k][1]
        v = sqrt(V[2b-1]^2 + V[2b]^2)
        _, _, Gp, _ = _complexload_reactive(c, x[_state_range(sys, k)], v)
        Iab[2b-1] += Gp * V[2b-1]
        Iab[2b] += Gp * V[2b]
    end
    busno = Dict(i => b for (b, i) in sys.net.index)
    comps = copy(L.comps)
    offsets = copy(L.offsets)
    sources = copy(L.sources)
    inj = copy(L.inj)
    for i in 1:nb
        vd, vq, id, iq = V[2i-1], V[2i], Iab[2i-1], Iab[2i]
        P = vd * id + vq * iq                          # V conj(I)
        Q = vq * id - vd * iq
        abs(P) + abs(Q) > tol || continue
        push!(comps, ConstantPowerSink("SINK_$(busno[i])", busno[i], SinkParams(P, Q, vd, vq), ParamDict()))
        push!(offsets, L.delta_coi)                    # stateless: an empty state range
        push!(sources, [InputSource(SRC_VD, 0.0, i, 0), InputSource(SRC_VQ, 0.0, i, 0)])
        push!(inj, (i, 1, 2))
    end
    return DAESystem(L.case, L.net, comps, offsets, L.n_diff, L.delta_coi, sources, inj, L.G, L.B,
                     L.load, L.slack, L.Vd_ref, L.Vq_ref, L.coi_members, L.coi_weights, L.coi_total,
                     L.omega_b, L.faults, L.pf)
end
