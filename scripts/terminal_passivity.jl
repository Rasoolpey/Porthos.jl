# Route A (TODO.md): cut the system at the machine terminals and test passivity with the
# incremental terminal power dV'dI. One command:
#
#   julia --project=. scripts/terminal_passivity.jl [scenario.json]
#
# At Porthos's equilibrium (healthy network), in the synchronous frame:
#   - units (machine + governor + exciter) and the network with its loads, from the exact
#     Jacobian; closing them again reproduces the system's eigenvalues (checked);
#   - per unit: min over w of eigmin Herm(Y_unit(jw)) (Y_unit = its absorbing admittance),
#     also for the machine alone, with its governor, with its exciter (which layer costs
#     passivity);
#   - the network, and the whole cut: eigmin Herm(Y_net + sum Y_unit) > 0 at every frequency
#     would allow a split unit by unit.
# Written to outputs/ph_audit/<case>/terminal_passivity.json.

println("loading Porthos...")
flush(stdout)
using Porthos
using LinearAlgebra
using Printf

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
case = load_case(sc.system_path)
eq = solve_equilibrium(case, sc)
sys, x, V = eq.sys, eq.x, eq.V
proj = physical_projection(sys, x, V)
units, net = terminal_models(sys, x, V; projection = proj)
ws = exp10.(range(-3, 3; length = 6001))

m = terminal_margins(units, net, sys.net.bus_ids; ws)

function herm_min(A, B, C, D, ws)
    F = hessenberg(A)
    best = (Inf, NaN)
    for w in ws
        Y = -(C * ((F - (im * w) * I) \ complex.(B)) + D)
        l = eigmin(Hermitian((Y + Y') / 2))
        l < best[1] && (best = (l, w))
    end
    return best
end
layers = Dict{String,Any}()
println()
@printf("%-10s %4s  %-26s %-26s %-26s %-26s\n", "unit", "bus", "machine alone", "+ governor", "+ exciter", "full unit")
for u in units
    sel(p) = [i for (i, s) in enumerate(u.states) if any(startswith(s, q * ".") for q in p)]
    mach = sel([u.name])
    gov = sel([c for c in u.components if startswith(c, "IEEEG")])
    exc = sel([c for c in u.components if startswith(c, "IEEET")])
    row = Dict{String,Any}()
    cells = String[]
    for (label, idx) in (("machine", mach), ("machine+governor", vcat(mach, gov)),
                         ("machine+exciter", vcat(mach, exc)), ("unit", collect(eachindex(u.states))))
        l, w = herm_min(u.A[idx, idx], u.B[idx, :], u.C[:, idx], u.D, ws)
        row[label] = Dict("min_eig_herm" => l, "at_rad_s" => w)
        push!(cells, @sprintf("%+.3e @ %.3g", l, w))
    end
    layers[u.name] = row
    @printf("%-10s %4d  %-26s %-26s %-26s %-26s\n", u.name, u.bus, cells...)
end
@printf("\nnetwork with loads: min eig Herm(Y_net) = %+.4e at %.4g rad/s\n",
        m["network"]["min_eig_herm"], m["network"]["at_rad_s"])
@printf("whole cut: min eig Herm(Y_net + sum Y_unit) = %+.4e at %.4g rad/s  (%s)\n",
        m["cut"]["min_eig_herm"], m["cut"]["at_rad_s"], m["cut"]["passes"] ? "passes" : "fails")

out = joinpath(root, "outputs", "ph_audit", splitext(basename(sc.system_path))[1])
path = Porthos.write_json(joinpath(out, "terminal_passivity.json"),
                          Dict("case" => sc.system_path, "margins" => m, "layers" => layers))
println("report -> ", path)
