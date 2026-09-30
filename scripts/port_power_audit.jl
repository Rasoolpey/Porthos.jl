# Nonlinear port-power residual audit (TODO.md step 4; PHPS work package 1 item 3). One
# command:
#
#   julia --project=. scripts/port_power_audit.jl [scenario.json]
#
# About 15 s after compilation. At the equilibrium, at random KCL-consistent states near it, and
# along the scenario's fault trajectory (IDA, log grid of the scenario):
#   - per component, the storage rate grad H' f against the sum of the contract's port
#     powers (evaluated from the component's own signals); the residual is the power it
#     dissipates (> 0) or creates (< 0); for machines also against the physical kinetic
#     storage H omega^2 instead of the shifted H (omega - 1)^2;
#   - for GENROU (a linear rotor, `rotor_structure`): the rotor balance dH_mag/dt = field supply
#     + stator exchange [id, iq]' B_s' Q z - rotor loss, and a proof that the rotor loss matrix
#     -sym(QA) is positive definite;
#   - the network identities: injected power = V' G V (+ fault and frequency-load terms) + V . g,
#     machine terminal power, ComplexLoad declared against drawn power.
# Report: outputs/ph_audit/<case>/port_power.json (summary and every sample).

using Porthos
using Printf

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
say("equilibrium...")
eq = solve_equilibrium(load_case(sc.system_path), sc)
sys = eq.sys
s = sc.solver
log_dt = something(s.log_dt, s.dt)
say("IDA run of the scenario ($(s.duration) s)...")
traj = simulate_ida(sys, vcat(eq.x, eq.V); dt = s.dt, duration = s.duration, log_dt = max(log_dt, 0.01),
                    rtol = s.rtol, atol = s.atol)
window = isempty(sys.faults) ? nothing : (minimum(f.t_start for f in sys.faults), maximum(f.t_end for f in sys.faults))
samples = power_audit_samples(sys, eq.x, eq.V; n_random = 16, trajectory = traj, fault_window = window)
say(length(samples), " samples (equilibrium, 16 random, ", length(traj.t), " on the trajectory); auditing...")
audit = port_power_audit(sys, samples)

outdir = joinpath(root, "outputs", "ph_audit", basename(dirname(sc.system_path)))
mkpath(outdir)
path = joinpath(outdir, "port_power.json")
report = Dict{String,Any}("scenario" => scenario, "fault_window" => window, "audit" => audit)
write_certificate(path, report)          # JSON with non-finite numbers as strings
say("report: ", path)

println("\nIdentities (largest error over all samples):")
for (k, v) in sort(collect(audit["identities_max_error"]); by = first)
    @printf("  %-26s %.3g\n", k, v)
end
println("\nPer model type: residual = sum of contract port powers - grad H' f (> 0 dissipates)")
for (t, d) in sort(collect(audit["by_type"]); by = first)
    @printf("  %-15s %2d units, ports %s: residual [%.4g, %.4g]%s\n", t, d["components"],
            d["complete"] ? "complete" : "INCOMPLETE", d["residual_min"], d["residual_max"],
            isempty(d["components_creating_energy"]) ? "" :
            ", creates energy: " * join(sort(d["components_creating_energy"]), " "))
    if haskey(d, "residual_physical_kinetic_min")
        cr = get(d, "components_creating_energy_physical_kinetic", String[])
        @printf("  %-15s    with the physical kinetic storage H omega^2: [%.4g, %.4g]%s\n", "",
                d["residual_physical_kinetic_min"], d["residual_physical_kinetic_max"],
                isempty(cr) ? ", dissipative at every sample" : ", creates energy: " * join(sort(cr), " "))
    end
    for key in ("mechanical_residual", "magnetic_residual", "terminal_residual",
                "stator_exchange", "rotor_loss")
        haskey(d, key * "_min") && @printf("  %-15s      %-20s [%.4g, %.4g]\n", "", key,
                                           d[key * "_min"], d[key * "_max"])
    end
    haskey(d, "rotor_loss_matrix_min_eig_lower_bound") &&
        @printf("  %-15s      rotor loss matrix -sym(QA) positive definite: lambda_min >= %.3g (rigorous)\n",
                "", d["rotor_loss_matrix_min_eig_lower_bound"])
end
