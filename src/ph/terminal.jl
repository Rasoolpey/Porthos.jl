# Route A: the terminal cut (TODO.md, "Plan: storage that keeps the dynamics").
#
# The system is cut at the machine terminals into
#   - units: a machine with its own controllers (the components whose outputs feed it:
#     governor, exciter), port = its terminal bus: input the bus voltage (Vd, Vq), output the
#     current the unit injects into the network, I = I_norton - y_norton V;
#   - the network: the Y-bus without the machines' Norton admittances, with the loads (their
#     states and voltage-dependent currents): input the bus voltages, output the current it
#     absorbs, which by KCL equals the units' injected currents.
# The linearisation is taken in the synchronous frame: the COI correction of the rotor-angle
# rows only fixes the reference frame and couples every machine to every other, so it is
# removed; the common rotation of all angles and phasors is then a neutral mode (eigenvalue 0).
# Supplies are incremental terminal powers dV'dI; the network's and the units' cancel by KCL.

"""
    TerminalModel

A unit (a machine and its controllers) at its terminal bus: `x' = A x + B dV`,
`dI = C x + D dV`, with `dV = (dVd, dVq)` the bus voltage and `dI` the current the unit
injects into the network.
"""
struct TerminalModel
    name::String
    components::Vector{String}
    bus::Int
    states::Vector{String}
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
end

Base.show(io::IO, m::TerminalModel) =
    print(io, "TerminalModel(", m.name, " at bus ", m.bus, ": ", join(m.components, "+"), ", ",
          length(m.states), " states)")

"""
    NetworkModel

The network with the loads, seen from the buses: `x' = A x + B dV`, `dI = C x + D dV`, with
`dV` all bus voltages `(Vd_1, Vq_1, ...)` and `dI` the current it absorbs at each bus.
"""
struct NetworkModel
    states::Vector{String}
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
end

"""
    sync_jacobian(sys, x, V) -> J

The Jacobian of `[f; g]` with respect to `[x; V]` in the synchronous frame: the COI
correction `-omega_b (omega_COI - 1)` is taken out of every rotor-angle row.
"""
function sync_jacobian(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    nd = sys.n_diff
    F(z) = vcat(dae_residual(sys, view(z, 1:nd), view(z, nd + 1:length(z)))...)
    J = ForwardDiff.jacobian(F, vcat(collect(float(x)), collect(float(V))))
    m = sys.coi_members
    if length(m) > 1
        for k in m, (q, w) in zip(m, sys.coi_weights)
            J[sys.offsets[k], sys.offsets[q] + 1] += sys.omega_b * w / sys.coi_total
        end
    end
    return J
end

"""
    terminal_models(sys, x, V; projection) -> (units, network)

Cut the system at the machine terminals (see the file header). Every coupling not through a
terminal voltage is checked to be absent.
"""
function terminal_models(sys::DAESystem, x::AbstractVector, V::AbstractVector;
                         projection::PhysicalProjection = physical_projection(sys, x, V))
    nd, nb = sys.n_diff, nbus(sys)
    J = sync_jacobian(sys, x, V)
    tol = 1e-9 * max(1.0, maximum(abs, J))
    keep = projection.keep
    names = projection.names
    vcol(b) = nd .+ (2b - 1:2b)                     # columns of bus b's voltage
    grow(b) = nd .+ (2b - 1:2b)                     # rows of bus b's KCL
    ins, _ = component_io(sys, x, V)
    bidx(bus_id) = findfirst(==(bus_id), sys.net.bus_ids)

    units = TerminalModel[]
    used = Set{Int}()
    for mk in sys.coi_members
        mach = sys.comps[mk]
        feeders = unique([s.index for s in sys.sources[mk] if s.kind == SRC_OUTPUT])
        members = vcat(mk, feeders)
        S = [i for k in members for i in _state_range(sys, k) if i in keep]
        b = bidx(bus(mach))
        others = setdiff(keep, S)
        maximum(abs, J[S, others]; init = 0.0) <= tol ||
            error("$(name(mach)): the unit depends on states outside it")
        vother = setdiff(nd + 1:nd + 2nb, vcol(b))
        maximum(abs, J[S, vother]; init = 0.0) <= tol ||
            error("$(name(mach)): the unit depends on another bus's voltage")
        # the machine's own V -> I (Norton current's voltage dependence, minus its admittance)
        jd, jq = sys.inj[mk][2], sys.inj[mk][3]
        u0 = ins[mk]
        Iout(v) = begin
            uu = collect(promote_type(eltype(v), eltype(u0)), u0)
            uu[1], uu[2] = v[1], v[2]
            y = outputs!(zeros(eltype(uu), noutputs(mach)), mach, x[_state_range(sys, mk)], uu)
            [y[jd], y[jq]]
        end
        dIdV = ForwardDiff.jacobian(Iout, collect(float(u0[1:2])))
        y = norton_admittance(mach)
        Ym = [real(y) -imag(y); imag(y) real(y)]
        push!(units, TerminalModel(name(mach), [name(sys.comps[k]) for k in members], bus(mach),
                                   names[S], J[S, S], J[S, vcol(b)], J[grow(b), S], dIdV - Ym))
        union!(used, S)
    end
    L = [i for i in keep if !(i in used)]
    maximum(abs, J[L, collect(used)]; init = 0.0) <= tol ||
        error("the loads depend on unit states")
    grows = nd + 1:nd + 2nb
    D = -J[grows, grows]
    for u in units
        b = bidx(u.bus)
        D[2b - 1:2b, 2b - 1:2b] .+= u.D
    end
    net = NetworkModel(names[L], J[L, L], J[L, grows], -J[grows, L], D)
    return units, net
end

"""`dI = G(jw) dV` of a unit (the current it injects)."""
transfer(m::TerminalModel, s::Number) = m.C * ((s * I - m.A) \ complex.(m.B)) + m.D

"""`dI = Y(jw) dV` of the network (the current it absorbs at every bus)."""
transfer(m::NetworkModel, s::Number) =
    isempty(m.states) ? complex.(m.D) : m.C * ((s * I - m.A) \ complex.(m.B)) + m.D

"""
    terminal_margins(units, net, bus_ids; ws) -> Dict

Frequency tests of the terminal cut (incremental supply dV'dI, power into each part):
- per unit, `min_w eigmin Herm(Y_u(jw))`, `Y_u = -G_u` its absorbing admittance: >= 0 means
  the unit alone is passive at its terminal; the most negative value is its shortage;
- the network, `min_w eigmin Herm(Y_net(jw))`;
- the cut, `min_w eigmin Herm(Y_net(jw) + sum_u E_u Y_u(jw) E_u')`: > 0 at every frequency
  means a storage split unit by unit exists with frequency-dependent terms; < 0 rules
  such a split out.
Each with the frequency where the minimum is reached.
"""
function terminal_margins(units::AbstractVector{TerminalModel}, net::NetworkModel,
                          bus_ids::AbstractVector{<:Integer}; ws = exp10.(range(-3, 3; length = 6001)))
    nb = length(bus_ids)
    pos = Dict(b => i for (i, b) in enumerate(bus_ids))
    Fu = [hessenberg(u.A) for u in units]
    Fn = isempty(net.states) ? nothing : hessenberg(net.A)
    umin = fill(Inf, length(units)); uat = fill(NaN, length(units))
    nmin, nat, tmin, tat = Inf, NaN, Inf, NaN
    for w in ws
        Yn = Fn === nothing ? complex.(net.D) :
             net.C * ((Fn - (im * w) * I) \ complex.(net.B)) + net.D
        l = eigmin(Hermitian((Yn + Yn') / 2))
        l < nmin && ((nmin, nat) = (l, w))
        T = copy(Yn)
        for (j, u) in enumerate(units)
            Yu = -(u.C * ((Fu[j] - (im * w) * I) \ complex.(u.B)) + u.D)
            lu = eigmin(Hermitian((Yu + Yu') / 2))
            lu < umin[j] && ((umin[j], uat[j]) = (lu, w))
            b = pos[u.bus]
            T[2b - 1:2b, 2b - 1:2b] .+= Yu
        end
        lt = eigmin(Hermitian((T + T') / 2))
        lt < tmin && ((tmin, tat) = (lt, w))
    end
    return Dict{String,Any}(
        "units" => Dict(u.name => Dict("min_eig_herm" => umin[j], "at_rad_s" => uat[j],
                                       "components" => u.components, "bus" => u.bus)
                        for (j, u) in enumerate(units)),
        "network" => Dict("min_eig_herm" => nmin, "at_rad_s" => nat),
        "cut" => Dict("min_eig_herm" => tmin, "at_rad_s" => tat, "passes" => tmin > 0))
end
