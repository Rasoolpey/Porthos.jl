# Run a fault scenario in PowerFactory and in Porthos and compare them. One command, from
# the Porthos.jl folder, with the PowerFactory window closed:
#
#   julia --project=. scripts/pf_compare_fault.jl [scenario.json]
#
# Default scenario: cases/IEEE39Bus_PF/bus_fault_bus16_150ms.json. PowerFactory runs the
# scenario's bus faults in a copy of the configured study case (pf/config.json); Porthos
# runs the scenario with IDA at its tolerances, as PHPS's validation did. Written to
# outputs/pf/<scenario name>/: pf_results.csv and run.json (PowerFactory), model.json (the
# PowerFactory model and load flow), compare.json (load flow and PHPS's trajectory metrics
# per machine) and figures.

println("loading Porthos...")
flush(stdout)
using Porthos
using Printf

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
name = splitext(basename(scenario))[1]
outdir = joinpath(root, "outputs", "pf", name)

# -- PowerFactory ----------------------------------------------------------------------
println("PowerFactory: model and load flow...")
flush(stdout)
model = pf_inspect(outdir = outdir)
println("PowerFactory: RMS simulation...")
flush(stdout)
pf = pf_simulate(scenario; outdir = outdir)
println("  ", pf)

# -- Porthos -----------------------------------------------------------------------------
println("Porthos: equilibrium and IDA...")
flush(stdout)
sc = load_scenario(scenario)
case = load_case(sc.system_path)
eq = solve_equilibrium(case, sc)
s = sc.solver
r = simulate_ida(eq.sys, vcat(eq.x, eq.V); dt = s.dt, duration = s.duration,
                 log_dt = s.log_dt === nothing ? s.dt : s.log_dt, rtol = s.rtol, atol = s.atol)
cols = csv_columns(eq.sys)
table = reduce(hcat, [csv_row(eq.sys, r.t[k], view(r.Y, :, k)) for k in eachindex(r.t)])
colidx = Dict{String,Int}()
for (j, c) in enumerate(cols)
    haskey(colidx, c) || (colidx[c] = j)
end
sim(n) = table[colidx[n], :]

# -- load flow ---------------------------------------------------------------------------
pfbus = Dict(b["name"] => b for b in model["load_flow"]["buses"])
ref = findfirst(==(31), eq.sys.net.bus_ids)            # reference machine's bus (G 02)
lf = Dict{String,Any}[]
for (i, b) in enumerate(eq.sys.net.bus_ids)
    p = pfbus[@sprintf("Bus %02d", b)]
    dth = rad2deg(eq.sys.pf.theta[i] - eq.sys.pf.theta[ref]) -
          (p["phiu"] - pfbus["Bus 31"]["phiu"])
    push!(lf, Dict("bus" => b, "dV_pu" => eq.sys.pf.V[i] - p["u"], "dtheta_deg" => dth))
end
dV = maximum(abs(x["dV_pu"]) for x in lf)
dth = maximum(abs(x["dtheta_deg"]) for x in lf)
@printf("load flow: max |V_Porthos - V_PF| = %.2e pu, max |angle difference| = %.2e deg (relative to bus 31)\n", dV, dth)

# -- trajectories ------------------------------------------------------------------------
faults = [e for e in sc.events if e isa BusFault]
t_fault = minimum(f.t_start for f in faults)
t_clear = maximum(f.t_end for f in faults)
m = pf_compare(r.t, sim, pf, case; t_fault, t_clear)
println()
println("Porthos (IDA) against PowerFactory, PHPS's metrics (rotor angle relative to G 01):")
@printf("  %-6s %-22s %12s %12s %12s %12s %10s\n", "PF", "Porthos", "angle rms", "angle max",
        "speed rms", "speed max", "P rms")
for (g, units) in pf_machine_map(case)
    x = m[g]
    @printf("  %-6s %-22s %9.4f deg %9.4f deg %9.2e pu %9.2e pu %7.4f pu\n", g, join(units, "+"),
            x["angle_rms_deg"], x["angle_max_deg"], x["speed_rms_pu"], x["speed_max_pu"], x["p_rms_pu"])
end
Porthos.write_json(joinpath(outdir, "compare.json"),
                   Dict("scenario" => scenario, "powerfactory" => pf.meta["powerfactory_dir"],
                        "project" => pf.meta["project"], "study_case" => pf.meta["base_study_case"],
                        "load_flow" => Dict("max_dV_pu" => dV, "max_dtheta_deg" => dth, "buses" => lf),
                        "machines" => m))
println("compare.json -> ", joinpath(outdir, "compare.json"))

# -- figures -----------------------------------------------------------------------------
if Base.find_package("CairoMakie") === nothing
    println("installing CairoMakie into Julia's default environment (first time only)...")
    run(`$(Base.julia_cmd()) --startup-file=no -e "using Pkg; Pkg.add(\"CairoMakie\")"`)
end
println("drawing the figures...")
flush(stdout)
@eval using CairoMakie
tpf, _ = Porthos._pf_series(pf.t, pf.t)
series(obj, var) = Porthos._pf_series(pf.t, pf_signal(pf, obj, var))[2]
mm = pf_machine_map(case)
refu = Dict(mm)["G 01"][1]
pre_pf = tpf .< t_fault - 0.05
pre_ph = r.t .< t_fault - 0.05
for (what, ylabel, file) in (("angle", "rotor angle relative to G 01 [deg]", "pf_angles.png"),
                             ("speed", "speed [pu]", "pf_speeds.png"))
    fig = Figure(size = (1200, 1150))
    for (k, (g, units)) in enumerate(mm)
        ax = Axis(fig[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1], title = g * " / " * join(units, "+"),
                  xlabel = "time [s]", ylabel = ylabel)
        vspan!(ax, t_fault, t_clear, color = (:red, 0.1))
        if what == "angle"
            a_pf = Porthos._unwrap_deg(series(g, "s:firot")) .- Porthos._unwrap_deg(series("G 01", "s:firot"))
            a_ph = rad2deg.(sim(units[1] * ".delta") .- sim(refu * ".delta"))
            y_pf = a_pf .- sum(a_pf[pre_pf]) / count(pre_pf)
            y_ph = a_ph .- sum(a_ph[pre_ph]) / count(pre_ph)
        else
            y_pf = series(g, "s:speed")
            y_ph = sim(units[1] * ".omega")
        end
        lines!(ax, tpf, y_pf, color = :black, label = "PowerFactory")
        lines!(ax, r.t, y_ph, color = :darkorange, linestyle = :dash, label = "Porthos")
        k == 1 && axislegend(ax, position = :rb, labelsize = 10)
    end
    Label(fig[0, :], "$(name): PowerFactory ($(basename(pf.meta["powerfactory_dir"]))) and Porthos (IDA)",
          fontsize = 16)
    save(joinpath(outdir, file), fig)
end
println("figures -> ", outdir)
