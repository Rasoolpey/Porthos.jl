# Simulate the IEEE-39 base case with a 3-phase fault at bus 16 (applied at 1.0 s, cleared
# after 150 ms) with BDF1 and IDA, and plot the results. One command, from the Porthos.jl
# folder:
#
#   julia --project=. scripts/run_base_fault.jl
#
# Results: outputs/IEEE39Bus_PF_base_bus16_150ms_bdf1/ and ..._ida/
#          (simulation_results.csv, simulation_results.jld2, run.json)
# Figures: outputs/figures/*.png
#
# Plotting uses CairoMakie. If it is missing, it is installed once into Julia's default
# environment (always available, nothing to activate).

println("loading Porthos (about 30 s)...")
flush(stdout)
using Porthos

root = normpath(joinpath(@__DIR__, ".."))
scenario = joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json")

# -- simulate ----------------------------------------------------------------------------
dirs = Dict{Symbol,String}()
for m in (:bdf1, :ida)
    out = joinpath(root, "outputs", "IEEE39Bus_PF_base_bus16_150ms_" * string(m))
    println("simulating 6 s with ", m, m === :bdf1 ? " (about 2 minutes)..." : " (under a minute)...")
    flush(stdout)
    t = @elapsed simulate(scenario; method = m, outdir = out)
    println(rpad(string(m), 5), " -> ", out, "  (", round(t; digits = 1), " s)")
    dirs[m] = out
end

# -- plot --------------------------------------------------------------------------------
if Base.find_package("CairoMakie") === nothing
    println("installing CairoMakie into Julia's default environment (first time only)...")
    run(`$(Base.julia_cmd()) --startup-file=no -e "using Pkg; Pkg.add(\"CairoMakie\")"`)
end
println("loading CairoMakie and drawing the figures (about a minute)...")
flush(stdout)
@eval using CairoMakie

function load(dir)
    lines_ = readlines(joinpath(dir, "simulation_results.csv"))
    header = split(lines_[1], ',')
    data = permutedims(reshape([parse(Float64, v) for l in lines_[2:end] for v in split(l, ',')],
                               length(header), :))
    return data, Dict(String(h) => j for (j, h) in enumerate(header))
end
col(d, name) = d[1][:, d[2][name]]
tt(d) = d[1][:, 1]

bdf1 = load(dirs[:bdf1])
ida = load(dirs[:ida])
figdir = joinpath(root, "outputs", "figures")
mkpath(figdir)
# eleven distinguishable colours (tab10 and black), so every machine has its own
colors = [CairoMakie.Makie.to_colormap(:tab10); CairoMakie.Makie.RGBf(0, 0, 0)]
machines = ["GENROU_$k" for k in 1:11]
fault!(ax) = vspan!(ax, 1.0, 1.15, color = (:red, 0.1))

# rotor angles and speeds (one legend, right of both axes)
fig = Figure(size = (1200, 820))
ax1 = Axis(fig[1, 1], title = "Rotor angles relative to GENROU_2 (slack machine, bus 31); fault at bus 16 from 1.00 to 1.15 s",
           xlabel = "time [s]", ylabel = "angle difference [deg]")
ax2 = Axis(fig[2, 1], title = "Rotor speeds", xlabel = "time [s]", ylabel = "omega [pu]")
fault!(ax1); fault!(ax2)
ref = col(bdf1, "GENROU_2.delta_deg")
for (k, m) in enumerate(machines)
    m == "GENROU_2" ||
        lines!(ax1, tt(bdf1), col(bdf1, m * ".delta_deg") .- ref, color = colors[k], label = m)
    lines!(ax2, tt(bdf1), col(bdf1, m * ".omega"), color = colors[k], label = m)
end
Legend(fig[1:2, 2], ax2, "machine")
save(joinpath(figdir, "base_bus16_rotor.png"), fig)

# voltages and controls
fig = Figure(size = (1200, 1000))
ax = Axis(fig[1, 1], title = "Bus voltage magnitudes", xlabel = "time [s]", ylabel = "|V| [pu]")
fault!(ax)
for (k, b) in enumerate((16, 15, 17, 24, 21, 39, 31, 30))
    lines!(ax, tt(bdf1), col(bdf1, "Vterm_Bus$b"), color = colors[k], label = "bus $b")
end
Legend(fig[1, 2], ax)
ax = Axis(fig[2, 1], title = "Field voltage of exciter IEEET1_4 (machine GENROU_4, bus 33)",
          xlabel = "time [s]", ylabel = "Efd [pu]")
fault!(ax)
lines!(ax, tt(bdf1), col(bdf1, "IEEET1_4.Efd"), color = colors[1], label = "IEEET1_4")
Legend(fig[2, 2], ax)
ax = Axis(fig[3, 1], title = "Mechanical torque: steam governor IEEEG1_4 (GENROU_4) and hydro governor IEEEG3_11 (GENROU_11, bus 30)",
          xlabel = "time [s]", ylabel = "Tm [pu]")
fault!(ax)
lines!(ax, tt(bdf1), col(bdf1, "IEEEG1_4.Tm"), color = colors[2], label = "IEEEG1_4")
lines!(ax, tt(bdf1), col(bdf1, "IEEEG3_11.Tm"), color = colors[3], label = "IEEEG3_11")
Legend(fig[3, 2], ax)
save(joinpath(figdir, "base_bus16_voltages_controls.png"), fig)

# BDF1 against IDA
fig = Figure(size = (1200, 480))
ax = Axis(fig[1, 1], title = "Same fault, two integrators: BDF1 (dt = 0.5 ms, solid) and IDA (adaptive, dashed)",
          xlabel = "time [s]", ylabel = "angle difference to GENROU_2 [deg]")
fault!(ax)
refi = col(ida, "GENROU_2.delta_deg")
for (k, m) in enumerate(("GENROU_4", "GENROU_5", "GENROU_9"))
    lines!(ax, tt(bdf1), col(bdf1, m * ".delta_deg") .- ref, color = colors[k], label = m * " BDF1")
    lines!(ax, tt(ida), col(ida, m * ".delta_deg") .- refi, color = colors[k], linestyle = :dash,
           label = m * " IDA")
end
Legend(fig[1, 2], ax)
save(joinpath(figdir, "base_bus16_bdf1_vs_ida.png"), fig)

println("figures -> ", figdir)
