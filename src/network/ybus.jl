# Bus admittance matrix.
#
# Stamps and their order follow PHPS (ybus.py YBusBuilder, dirac/dae_compiler.py
# `_build_full_ybus`, `_parse_bus_faults`), and the arithmetic reproduces CPython's, so the
# matrices agree with PHPS to the last bit, not just within tolerance:
#
#   lines (file order) -> shunts -> constant-impedance loads (optional)
#   -> generator Norton admittances (optional, component order)
#
# A line is a pi section with an off-nominal complex tap a = tap * e^{j phi} on `bus1`:
#
#   Y11 += y/|a|^2 + j b/2     Y12 -= y/conj(a)
#   Y22 += y + j b/2           Y21 -= y/a            with y = 1/(r + jx).

"""
    Network

Bus ordering shared by every network matrix: bus ids sorted ascending (as PHPS does), and
the map from bus id to row.
"""
struct Network
    bus_ids::Vector{Int}
    index::Dict{Int,Int}
end

function Network(bus_ids::AbstractVector{<:Integer})
    ids = sort!(collect(Int, bus_ids))
    return Network(ids, Dict(b => i for (i, b) in enumerate(ids)))
end

Network(case::Case) = Network([b.idx for b in case.buses])

nbus(net::Network) = length(net.bus_ids)
bus_index(net::Network, bus::Integer) = net.index[bus]

# CPython complex division `_Py_c_quot` (Objects/complexobject.c), used by PHPS for
# 1.0 / complex(r, x). Julia's `inv(::ComplexF64)` rounds differently in the last bit.
function _py_cdiv(a::ComplexF64, b::ComplexF64)
    ar, ai, br, bi = reim(a)..., reim(b)...
    abs_br, abs_bi = abs(br), abs(bi)
    if abs_br >= abs_bi
        abs_br == 0.0 && throw(DivideError())
        ratio = bi / br
        denom = br + bi * ratio
        return complex((ar + ai * ratio) / denom, (ai - ar * ratio) / denom)
    elseif abs_bi >= abs_br
        ratio = br / bi
        denom = br * ratio + bi
        return complex((ar * ratio + ai) / denom, (ai * ratio - ar) / denom)
    end
    return complex(NaN, NaN)
end
_py_cdiv(a::Real, b::ComplexF64) = _py_cdiv(complex(Float64(a), 0.0), b)
# complex / float in CPython 3.13: promote the float to complex and use `_Py_c_quot`.
_py_cdiv(a::ComplexF64, b::Real) = _py_cdiv(a, complex(Float64(b), 0.0))

"""Accumulates stamps in order, starting each entry from zero (like a zeroed numpy array)."""
struct _Stamps
    n::Int
    vals::Dict{Tuple{Int,Int},ComplexF64}
    order::Vector{Tuple{Int,Int}}
end
_Stamps(n::Int) = _Stamps(n, Dict{Tuple{Int,Int},ComplexF64}(), Tuple{Int,Int}[])

function _add!(s::_Stamps, i::Int, j::Int, v::ComplexF64)
    k = (i, j)
    if haskey(s.vals, k)
        s.vals[k] += v
    else
        s.vals[k] = zero(ComplexF64) + v
        push!(s.order, k)
    end
    return s
end
_sub!(s::_Stamps, i::Int, j::Int, v::ComplexF64) =
    haskey(s.vals, (i, j)) ? (s.vals[(i, j)] -= v; s) : _add!(s, i, j, zero(ComplexF64) - v)

function SparseArrays.sparse(s::_Stamps)
    I = [k[1] for k in s.order]
    J = [k[2] for k in s.order]
    V = [s.vals[k] for k in s.order]
    return sparse(I, J, V, s.n, s.n)
end

function _stamp_line!(s::_Stamps, net::Network, l::LineData)
    (haskey(net.index, l.bus1) && haskey(net.index, l.bus2)) || return s
    i, j = net.index[l.bus1], net.index[l.bus2]
    y_series = _py_cdiv(1.0, complex(l.r, l.x))
    y_shunt = complex(0.0, l.b / 2.0)
    # a = tap * (cos(phi) + 1j*sin(phi)) with CPython float*complex semantics
    c, sn = cos(l.phi), sin(l.phi)
    a = complex(l.tap * c - 0.0 * sn, l.tap * sn + 0.0 * c)
    a_conj = conj(a)
    mag_a2 = hypot(real(a), imag(a))^2
    _add!(s, i, i, _py_cdiv(y_series, mag_a2) + y_shunt)
    _add!(s, j, j, y_series + y_shunt)
    _sub!(s, i, j, _py_cdiv(y_series, a_conj))
    _sub!(s, j, i, _py_cdiv(y_series, a))
    return s
end

function _stamp_shunts!(s::_Stamps, net::Network, case::Case)
    for sh in case.shunts
        haskey(net.index, sh.bus) || continue
        i = net.index[sh.bus]
        _add!(s, i, i, complex(sh.g, sh.b))
    end
    return s
end

"""
    load_power_by_bus(case, net) -> Vector{Tuple{Int,Float64,Float64}}

Total `PQ` consumption per bus as `(bus, P, Q)`, buses in order of first appearance in the
`PQ` table (PHPS accumulates in a defaultdict with that order).
"""
function load_power_by_bus(case::Case, net::Network)
    acc = Dict{Int,Tuple{Float64,Float64}}()
    order = Int[]
    for pq in case.pq
        haskey(net.index, pq.bus) || continue
        if !haskey(acc, pq.bus)
            acc[pq.bus] = (0.0, 0.0)
            push!(order, pq.bus)
        end
        P, Q = acc[pq.bus]
        acc[pq.bus] = (P + pq.p0, Q + pq.q0)
    end
    return [(b, acc[b]...) for b in order]
end

_bus_v0(case::Case, bus::Int) =
    (i = findfirst(b -> b.idx == bus, case.buses); i === nothing ? 1.0 : case.buses[i].v0)

function _stamp_loads!(s::_Stamps, net::Network, case::Case)
    for (bus, P, Q) in load_power_by_bus(case, net)
        (P == 0 && Q == 0) && continue
        i = net.index[bus]
        V0 = _bus_v0(case, bus)
        V02 = max(V0 * V0, 1e-12)
        # (P - 1j*Q) / V02, CPython: 1j*Q = complex(0*Q - 1*0, 0*0 + 1*Q)
        jQ = complex(0.0 * Q - 1.0 * 0.0, 0.0 * 0.0 + 1.0 * Q)
        _add!(s, i, i, _py_cdiv(complex(P, 0.0) - jQ, V02))
    end
    return s
end

"""
    NortonStamp

A generator's Norton admittance `1/(ra + j xd'')` at its bus.
"""
struct NortonStamp
    component::String
    bus::Int
    ra::Float64
    xd_pp::Float64
end

"""Component types that PHPS treats as generators with a Norton admittance."""
const NORTON_TYPES = Set(["GENROU", "GENROU_PHS", "GENROU_PHTRUE", "GENSAL", "GENSAL_PHTRUE",
                          "GENCLS", "GFM_VSM_PHTRUE", "GFM_DROOP_PHTRUE", "GFM_VOC_PHTRUE"])

"""
    norton_stamps(case) -> Vector{NortonStamp}

The Norton admittances PHPS adds to the DAE Y-bus, in component order: every generator
except the current-source converters (`GFL_PHTRUE`, `GFL_ZIF_PHTRUE`). `xd''` falls back to
`xd1`, then 0.2, as in PHPS. (At P4 this moves into each model's `injection`.)
"""
function norton_stamps(case::Case)
    out = NortonStamp[]
    for spec in case.components
        spec.type in NORTON_TYPES || continue
        p = component_params(case, spec)
        ra = haskey(p, "ra") ? param_value(p["ra"]) : 0.0
        xd_pp = haskey(p, "xd_double_prime") ? param_value(p["xd_double_prime"]) :
                haskey(p, "xd1") ? param_value(p["xd1"]) : 0.2
        push!(out, NortonStamp(spec.name, Int(p["bus"]), ra, xd_pp))
    end
    return out
end

function _stamp_norton!(s::_Stamps, net::Network, stamps)
    for st in stamps
        haskey(net.index, st.bus) || continue
        i = net.index[st.bus]
        z = complex(st.ra, st.xd_pp)
        abs(z) < 1e-6 && (z = complex(0.0, 0.0001))
        _add!(s, i, i, _py_cdiv(1.0, z))
    end
    return s
end

"""
    ybus(case; net = Network(case), lines = case.lines, loads = false, norton = ())
        -> SparseMatrixCSC{ComplexF64}

Bus admittance matrix in the bus order of `net`:

- `ybus(case)`: lines and shunts, the power-flow matrix (PHPS `build(include_loads=False)`);
- `loads = true`: adds each PQ load as the constant impedance `(P - jQ)/v0^2`;
- `norton`: generator Norton admittances to add (see [`norton_stamps`](@ref)).

See [`ybus_dae`](@ref) for the matrix the DAE uses.
"""
function ybus(case::Case; net::Network = Network(case),
              lines::AbstractVector{LineData} = case.lines, loads::Bool = false,
              norton = ())
    s = _Stamps(nbus(net))
    for l in lines
        _stamp_line!(s, net, l)
    end
    _stamp_shunts!(s, net, case)
    loads && _stamp_loads!(s, net, case)
    _stamp_norton!(s, net, norton)
    return sparse(s)
end

"""
    ybus_pf(case) -> SparseMatrixCSC{ComplexF64}

The power-flow Y-bus: lines and shunts only; loads enter the power flow as P/Q
specifications.
"""
ybus_pf(case::Case; kwargs...) = ybus(case; kwargs...)

"""
    ybus_dae(case) -> SparseMatrixCSC{ComplexF64}

The DAE network matrix of PHPS (`DiracCompiler.Y_full`): lines, shunts, constant-impedance
PQ loads at the case `v0`, and the generator Norton admittances. No Kron reduction; every
bus stays an algebraic variable.
"""
ybus_dae(case::Case; net::Network = Network(case), lines = case.lines) =
    ybus(case; net, lines, loads = true, norton = norton_stamps(case))

"""
    LoadAdmittances

Per-bus load data the DAE residual uses (bus order of the network): `G = P0/V0^2`,
`B = -Q0/V0^2`, the scheduled `P`, `Q`, and the frequency coefficients `kpf`, `kqf`.
`COMPLEXLOAD` components override `G`, `B`, `kpf`, `kqf` at their bus.
"""
struct LoadAdmittances
    G::Vector{Float64}
    B::Vector{Float64}
    P::Vector{Float64}
    Q::Vector{Float64}
    kpf::Vector{Float64}
    kqf::Vector{Float64}
    has_complex_loads::Bool
end

"""
    load_admittances(case; net = Network(case)) -> LoadAdmittances

Port of the load block of PHPS `_build_full_ybus` and the kpf/kqf fallback of
`DiracCompiler.build`.
"""
function load_admittances(case::Case; net::Network = Network(case))
    n = nbus(net)
    G, B, P, Q = zeros(n), zeros(n), zeros(n), zeros(n)
    for (bus, p, q) in load_power_by_bus(case, net)
        i = net.index[bus]
        V0 = _bus_v0(case, bus)
        V02 = max(V0 * V0, 1e-12)
        G[i] = p / V02
        B[i] = -q / V02
        P[i] = p
        Q[i] = q
    end
    kpf, kqf = zeros(n), zeros(n)
    has_complex = false
    for spec in case.components
        spec.type == "COMPLEXLOAD" || continue
        has_complex = true
        bus = Int(spec.params[:bus])
        haskey(net.index, bus) || continue
        i = net.index[bus]
        kpf[i] = param(spec, :kpf, 0.0)
        kqf[i] = param(spec, :kqf, 0.0)
        P0, Q0 = param(spec, :P0, 0.0), param(spec, :Q0, 0.0)
        V0 = param(spec, :V0, 1.0)
        V02 = max(V0 * V0, 1e-12)
        if P0 != 0.0 || Q0 != 0.0
            G[i] = P0 / V02
            B[i] = -Q0 / V02
        end
    end
    if !has_complex
        fill!(kpf, case.config.kpf)
        fill!(kqf, case.config.kqf)
    end
    return LoadAdmittances(G, B, P, Q, kpf, kqf, has_complex)
end
