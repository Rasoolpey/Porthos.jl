# The GENROU stator exchange through the KCL branch (TODO.md step 5a).
#
# GENROU's rotor is a linear port system (`rotor_structure`): z' = A z + B_f Efd + B_s i,
# i = [id, iq]. In the current-corrected rotor coordinates w = z - D i, D = -A^-1 B_s,
#
#   z' = A w + B_f Efd,    d/dt (w'Qw/2) = Efd B_f'Qw - w'(-sym QA)w - w'QD di/dt,
#
# where Qw is the rotor current (coenergy variable), B_f'Qw = `i_fd_energy` the field flow
# of this storage (an audit quantity: the runtime `i_fd = B_f'Qz` output and the IEEET1
# reservoir it feeds are unchanged), and -w'QD di/dt the stator exchange as a
# transformer-type port. The stator currents are algebraic: on the KCL branch they are
# functions of the state, i(x) = i(x, V(x)), so di/dt = (di/dx) x' with
# di/dx = i_x - i_V g_V^-1 g_x. The complementary term the stator/network side must supply
# is therefore the one-form
#
#   omega = a(x) . dx,    a(x) = sum over GENROU units of (di_k/dx)' D_k' Q_k w_k,
#
# and a network storage U with dU/dt = +w'QD di/dt exists exactly when omega is exact (a
# gradient). `exchange_one_form` evaluates it in the section coordinates of a
# `SectionModel`, `one_form_exactness` tests reference-angle invariance and the symmetry of
# its Jacobian (zero curl), and `one_form_path_test` compares its integral along two paths.
# Nothing here assumes that U exists.

"""
    current_correction(c) -> Matrix or nothing

`D = -A^-1 B_s` for a machine with a linear `rotor_structure`: `z' = A (z - D i) + B_f Efd`.
"""
function current_correction(c::AbstractComponent)
    rs = rotor_structure(c)
    return rs === nothing ? nothing : -(rs.A \ rs.Bs)
end

"""
    kcl_solve(sys, x, V0; faults_on = false, iterations = 3) -> V

The KCL voltages at `x`, in the element type of `x`, from `V0`, a converged Float64 solution
at the value of `x`: iterations `V <- V - g_V^-1 g(x, V)` with `g_V` the Float64 Jacobian at
`(value(x), V0)`. The value part is already converged, and for the dual parts this map has
zero contraction, so each derivative order is exact after one more step: with duals the
result carries the derivatives of the KCL branch `V(x)` (the implicit function), up to
order `iterations - 1`.
"""
function kcl_solve(sys::DAESystem, x::AbstractVector{T}, V0::AbstractVector{Float64};
                   faults_on::Bool = false, iterations::Integer = 3) where {T}
    x0 = _plain.(x)
    Jinv = inv(ForwardDiff.jacobian(v -> dae_residual(sys, x0, v; faults_on)[2], collect(V0)))
    V = convert(Vector{T}, collect(V0))
    for _ in 1:iterations
        V = V .- Jinv * dae_residual(sys, x, V; faults_on)[2]
    end
    return V
end

"""
    rotor_current_coordinates(sys, x, V) -> Dict{String,NamedTuple}

For every machine with a linear rotor: `z`, the stator currents `i = [id, iq]`, `w = z - D i`,
the rotor currents `Qw`, `i_fd_energy = B_f'Qw` (against the runtime `i_fd = B_f'Qz`), the
field supply `Efd i_fd_energy` and the loss `-w' sym(QA) w`.
"""
function rotor_current_coordinates(sys::DAESystem, x::AbstractVector, V::AbstractVector)
    ins, outs = component_io(sys, x, V)
    out = Dict{String,Any}()
    for (k, c) in enumerate(sys.comps)
        rs = rotor_structure(c)
        rs === nothing && continue
        on, inn = output_names(c), input_names(c)
        z = x[_state_range(sys, k)][rs.states]
        i = [outs[k][findfirst(==(rs.currents[1]), on)], outs[k][findfirst(==(rs.currents[2]), on)]]
        D = -(rs.A \ rs.Bs)
        w = z .- D * i
        Qw = rs.Q * w
        Efd = ins[k][findfirst(==("Efd"), inn)]
        L = -(rs.Q * rs.A .+ rs.A' * rs.Q) ./ 2
        out[name(c)] = (z = z, i = i, w = w, Qw = Qw, i_fd_energy = dot(rs.Bf, Qw),
                        i_fd_runtime = outs[k][findfirst(==("i_fd"), on)],
                        field_supply = Efd * dot(rs.Bf, Qw), loss = dot(w, L * w))
    end
    return out
end

# the stacked stator currents of the linear-rotor machines at x (voltages on the KCL branch)
function _stator_currents(sys::DAESystem, x::AbstractVector, V0; faults_on::Bool = false)
    V = kcl_solve(sys, x, V0; faults_on)
    _, outs = component_io(sys, x, V)
    I = eltype(V)[]
    for (k, c) in enumerate(sys.comps)
        rs = rotor_structure(c)
        rs === nothing && continue
        on = output_names(c)
        push!(I, outs[k][findfirst(==(rs.currents[1]), on)], outs[k][findfirst(==(rs.currents[2]), on)])
    end
    return I
end

# the stacked coefficients D'Qw of the linear-rotor machines at x
function _exchange_weights(sys::DAESystem, x::AbstractVector, I::AbstractVector)
    c_ = eltype(I)[]
    j = 0
    for (k, c) in enumerate(sys.comps)
        rs = rotor_structure(c)
        rs === nothing && continue
        z = x[_state_range(sys, k)][rs.states]
        i = I[j+1:j+2]
        j += 2
        D = -(rs.A \ rs.Bs)
        append!(c_, D' * (rs.Q * (z .- D * i)))
    end
    return c_
end

"""
    exchange_one_form(m::SectionModel, eta; sys = m.sys, V0 = m.V0) -> Vector

The coefficients of the stator-exchange one-form in section coordinates,
`a(eta) = (dI/deta)' c`, with `I` the stacked stator currents of the linear-rotor machines
on the KCL branch of `sys` (by default the model's; a variant such as `lossless_variant`
can be passed) and `c` the stacked `D'Qw`. `a . eta'` is `w'QD di/dt` summed over the units.
`V0` must be a converged KCL solution near `lift(m, eta)`.
"""
function exchange_one_form(m::SectionModel, eta::AbstractVector; sys::DAESystem = m.sys,
                           V0::AbstractVector = m.V0)
    x = lift(m, eta)
    Vc = solve_network(sys, _plain.(x), _plain.(V0))
    I = _stator_currents(sys, x, Vc)
    J = ForwardDiff.jacobian(e -> _stator_currents(sys, lift(m, e), Vc), eta)
    return J' * _exchange_weights(sys, x, I)
end

_plain(v::Real) = Float64(v)
_plain(v::ForwardDiff.Dual) = _plain(ForwardDiff.value(v))

"""
    one_form_exactness(m, eta; sys = m.sys, V0 = m.V0, pairs = 6, seed = 3) -> NamedTuple

The curl test of `exchange_one_form` at `eta`: for random unit direction pairs `(u, v)`,
`v' Ja u` against `u' Ja v`, with `Ja u` the directional derivative of `a` (nested duals
through the KCL branch). An exact one-form has a symmetric Jacobian, so the two agree;
`curl` is the largest relative difference `|v'Ja u - u'Ja v| / (|v'Ja u| + |u'Ja v|)`.
Also returns `a` and the pairs' values.
"""
function one_form_exactness(m::SectionModel, eta::AbstractVector; sys::DAESystem = m.sys,
                            V0::AbstractVector = m.V0, pairs::Integer = 6, seed::Integer = 3)
    rng = Random.Xoshiro(seed)
    e0 = collect(float(eta))
    a = exchange_one_form(m, e0; sys, V0)
    dir(u) = ForwardDiff.derivative(t -> exchange_one_form(m, e0 .+ t .* u; sys, V0), 0.0)
    vals = NamedTuple[]
    for _ in 1:pairs
        u = normalize(randn(rng, length(e0)))
        v = normalize(randn(rng, length(e0)))
        vu, uv = dot(v, dir(u)), dot(u, dir(v))
        push!(vals, (vJu = vu, uJv = uv, relative = abs(vu - uv) / max(abs(vu) + abs(uv), floatmin())))
    end
    return (a = a, pairs = vals, curl = maximum(p -> p.relative, vals))
end

"""
    one_form_path_test(m, eta1; sys = m.sys, V0 = m.V0, order = 8) -> NamedTuple

The integral of the one-form from `eta = 0` to `eta1` along the straight segment and along a
two-leg path (first the coordinates of the first half of the entries, then the rest), by
Gauss-Legendre quadrature of order `order` per leg; equal for an exact form.
"""
function one_form_path_test(m::SectionModel, eta1::AbstractVector; sys::DAESystem = m.sys,
                            V0::AbstractVector = m.V0, order::Integer = 8)
    t, wq = _gauss_legendre(order)
    seg(a, b) = sum(wq[q] / 2 * dot(exchange_one_form(m, a .+ (1 + t[q]) / 2 .* (b .- a); sys, V0), b .- a)
                    for q in eachindex(t))
    n = length(eta1)
    mid = copy(eta1)
    mid[n÷2+1:end] .= 0
    straight = seg(zeros(n), eta1)
    twoleg = seg(zeros(n), mid) + seg(mid, eta1)
    return (straight = straight, two_leg = twoleg, difference = straight - twoleg)
end

function _gauss_legendre(n::Integer)
    # Golub-Welsch
    b = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    F = eigen(SymTridiagonal(zeros(n), b))
    return F.values, 2 .* F.vectors[1, :] .^ 2
end

"""
    lossless_variant(sys) -> DAESystem

The same system with a lossless network: every Y-bus conductance, the constant-impedance
load conductances, the frequency-dependent active loads and the ComplexLoads' active
power (`P0 = 0`) removed; susceptances, reactive loads and the machines unchanged. For
testing the network identities (step 5), not a model of the case; its KCL solution at a
state of the original system is generally near, not at, the original voltages.
"""
function lossless_variant(sys::DAESystem)
    comps = AbstractComponent[component_role(c) === :load && hasparam(c, "P0") ?
                              with_params(c, Dict("P0" => 0.0)) : c for c in sys.comps]
    la = sys.load
    load = LoadAdmittances(zero(la.G), copy(la.B), zero(la.P), copy(la.Q), zero(la.kpf),
                           copy(la.kqf), la.has_complex_loads)
    return DAESystem(sys.case, sys.net, comps, sys.offsets, sys.n_diff, sys.delta_coi, sys.sources,
                     sys.inj, zero(sys.G), copy(sys.B), load, sys.slack, sys.Vd_ref, sys.Vq_ref,
                     sys.coi_members, sys.coi_weights, sys.coi_total, sys.omega_b, sys.faults, sys.pf)
end

"""
    exchange_curl(m, eta; sys = m.sys, V0 = m.V0) -> Matrix

The curl of the stator-exchange one-form in section coordinates, `Ja - Ja'`, in closed
form: with `c = D'Q(z - D I)`, the part `-(D'QD I)' dI` is exact and `(D'Qz)' dI` is not,
so `Ja - Ja' = (dI/deta)' D'Q (dz/deta) - (dz/deta)' Q D (dI/deta)` (block-diagonal `D`, `Q`
over the units). It vanishes only where the rotor states and the stator currents do not vary
independently on the KCL branch, so no scalar storage has the one-form as its gradient.
"""
function exchange_curl(m::SectionModel, eta::AbstractVector; sys::DAESystem = m.sys,
                       V0::AbstractVector = m.V0)
    x = lift(m, eta)
    Vc = solve_network(sys, _plain.(x), _plain.(V0))
    dI = ForwardDiff.jacobian(e -> _stator_currents(sys, lift(m, e), Vc), collect(float(eta)))
    dz = ForwardDiff.jacobian(e -> _rotor_states(sys, lift(m, e)), collect(float(eta)))
    blocks = [(-(rs.A \ rs.Bs), rs.Q) for rs in (rotor_structure(c) for c in sys.comps) if rs !== nothing]
    nz, ni = 4 * length(blocks), 2 * length(blocks)
    DQ = zeros(ni, nz)                       # block-diagonal D'Q
    for (u, (D, Q)) in enumerate(blocks)
        DQ[2u-1:2u, 4u-3:4u] .= D' * Q
    end
    K = dI' * DQ * dz
    return K .- K'
end

function _rotor_states(sys::DAESystem, x::AbstractVector)
    out = eltype(x)[]
    for (k, c) in enumerate(sys.comps)
        rs = rotor_structure(c)
        rs === nothing && continue
        append!(out, x[_state_range(sys, k)][rs.states])
    end
    return out
end
