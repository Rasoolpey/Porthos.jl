# Network changes caused by events.

"""
    FaultShunt

A fault admittance `g + jb` at bus `bus` (row `index` of the network) on `[t_start, t_end)`.
"""
struct FaultShunt
    bus::Int
    index::Int
    t_start::Float64
    t_end::Float64
    g::Float64
    b::Float64
end

"""
    fault_admittance(r, x) -> (g, b)

`1/(r + jx) = (r - jx)/(r^2 + x^2)`, computed as PHPS's DAE compiler does (`z^2` clipped at
`1e-20`).
"""
function fault_admittance(r::Float64, x::Float64)
    z_sq = r * r + x * x
    z_sq < 1e-20 && (z_sq = 1e-20)
    return r / z_sq, -x / z_sq
end

"""
    fault_shunts(net, events) -> Vector{FaultShunt}

The bus-fault shunts of a scenario, in event order. Faults at buses outside the network are
skipped, as in PHPS.
"""
function fault_shunts(net::Network, events)
    out = FaultShunt[]
    for ev in events
        ev isa BusFault || continue
        haskey(net.index, ev.bus) || continue
        g, b = fault_admittance(ev.r, ev.x)
        push!(out, FaultShunt(ev.bus, net.index[ev.bus], ev.t_start, ev.t_end, g, b))
    end
    return out
end

"""
    with_fault(Y, shunts) -> SparseMatrixCSC{ComplexF64}

`Y` with the fault shunts added on the diagonal (the fault-on network).
"""
function with_fault(Y::SparseMatrixCSC{ComplexF64}, shunts)
    Yf = copy(Y)
    for f in shunts
        Yf[f.index, f.index] += complex(f.g, f.b)
    end
    return Yf
end

"""
    split_line_for_fault(case, lf::LineFault) -> (lines, fault_bus, bus_fault)

Topology for a fault part-way along a line, as PHPS `SimulationRunner` builds it: the line
is replaced by two sections through a new bus (id 99, or one above the largest bus id if
99 is taken), whose series and shunt parameters are split by `distance`, and a
[`BusFault`](@ref) is placed at the new bus. The new bus starts at `v0 = 1`, `a0 = 0`.
"""
function split_line_for_fault(case::Case, lf::LineFault)
    k = findfirst(l -> l.idx == lf.line_idx, case.lines)
    k === nothing && throw(ArgumentError("Line $(lf.line_idx) not found for LineFault"))
    line = case.lines[k]
    max_bus = maximum(b.idx for b in case.buses)
    fault_bus = 99 <= max_bus ? max_bus + 1 : 99
    d = lf.distance
    part1 = LineData(line.idx * "_part1", line.bus1, fault_bus, line.r * d, line.x * d,
                     line.b * d, line.tap, line.phi)
    part2 = LineData(line.idx * "_part2", fault_bus, line.bus2, line.r * (1 - d),
                     line.x * (1 - d), line.b * (1 - d), line.tap, line.phi)
    lines = [l for l in case.lines if l.idx != lf.line_idx]
    append!(lines, (part1, part2))
    return lines, fault_bus, BusFault(fault_bus, lf.r, lf.x, lf.t_start, lf.t_end)
end
