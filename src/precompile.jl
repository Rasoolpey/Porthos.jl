# Precompile workload (roadmap 2.6): run the simulation path once while the package is
# precompiled, so that the compiled code is cached and a fresh `julia --project=.` process
# does not spend its first 20 s compiling. The workload is the base case with the
# bus-16 fault scenario: load, equilibrium, two BDF1 steps, a few milliseconds of IDA (fault
# segments included), the CSV and JLD2 writers and run.json. Nothing it computes is kept.

using PrecompileTools: @setup_workload, @compile_workload

# the system with its faults moved to [t0, t1)
function _with_fault_window(sys::DAESystem, t0, t1)
    faults = [FaultShunt(f.bus, f.index, t0, t1, f.g, f.b) for f in sys.faults]
    return DAESystem((k === :faults ? faults : getfield(sys, k) for k in fieldnames(DAESystem))...)
end

@setup_workload begin
    scenario_path = joinpath(PORTHOS_ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json")
    if isfile(scenario_path)
        @compile_workload begin
            sc = load_scenario(scenario_path)
            case = load_case(sc.system_path)
            eq = solve_equilibrium(case, sc)
            y0 = vcat(eq.x, eq.V)
            dt = sc.solver.dt
            r1 = simulate_bdf1(eq.sys, y0; dt, duration = 2dt)
            # a fault switching inside the window, so the segment restarts are compiled too
            fsys = _with_fault_window(eq.sys, 0.001, 0.002)
            r2 = simulate_ida(fsys, y0; dt, duration = 0.003)
            mktempdir() do dir
                write_results_csv(joinpath(dir, "simulation_results.csv"), r2)
                write_results_jld2(joinpath(dir, "simulation_results.jld2"), r1)
                write_json(joinpath(dir, "run.json"), run_metadata(r1, sc))
            end
        end
    end
end
