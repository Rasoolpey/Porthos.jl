# Port-Hamiltonian audit of a case at Porthos's own equilibrium. One command, from the
# Porthos.jl folder:
#
#   julia --project=. scripts/ph_audit.jl [scenario.json]
#
# Default scenario: cases/IEEE39Bus_PF/bus_fault_bus16_150ms.json (only its case is used; the
# audit is on the healthy network). It reports:
#   - the shifted-storage audit of the assembled storage H (positivity: rank of Hess H on
#     the physical coordinates; decay: eigenvalues of sym(S A) on the common-angle section;
#     exact dH_s/dt along the worst direction, by storage component);
#   - every governor's port (speed -> Tm, H = -G): G(0), where Re H(jw) changes sign, the
#     passivity shortage (-min Re H), and the zeros of G (a right-half-plane zero rules out
#     passivity at the port for any storage).
# Written to outputs/ph_audit/<case name>/ph_audit.json. A diagnostic, not a proof.

println("loading Porthos...")
flush(stdout)
using Porthos
using Printf

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
case = load_case(sc.system_path)
println("equilibrium of ", basename(sc.system_path), "...")
flush(stdout)
eq = solve_equilibrium(case, sc)
sys, x, V = eq.sys, eq.x, eq.V

a = shifted_storage_audit(sys, x, V)
println()
@printf("storage: %d components; %d states -> %d physical (%d reservoirs, %d held, delta_COI), section %d\n",
        length(a["storage_components"]), a["n_states"], a["n_keep"], length(a["reservoir_states"]),
        length(a["held_states"]), a["n_section"])
@printf("positivity: rank of Hess H = %d of %d  (%s)\n", a["storage_rank"], a["n_keep"],
        a["storage_rank"] == a["n_keep"] ? "positive definite" : "rank deficient: H cannot bound every state")
@printf("decay: sym(S A) on the section has %d positive eigenvalues, range [%.6g, %.6g]\n",
        a["sym_SA_positive"], a["sym_SA_eig_min"], a["sym_SA_eig_max"])
worst = a["exact_along_top_eigenvector"][end]
@printf("       exact dH_s/dt = %+.3e at eps = %g along the top eigenvector (H_s = %.3e)\n",
        worst["dHs_dt"], worst["eps"], worst["Hs"])

ins, _ = component_io(sys, x, V)
ports = Dict{String,Any}()
println()
println("governor ports, speed -> Tm, H = -G:")
@printf("  %-10s %-14s %10s %22s %12s  %s\n", "governor", "type", "H(0)", "Re H(jw) sign change", "shortage", "zeros of G")
for (k, c) in enumerate(sys.comps)
    model_type(c) in ("IEEEG1_PHTRUE", "IEEEG3_PHTRUE") || continue
    r = sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
    m = port_model(c, x[r], ins[k]; input = "omega", output = "Tm")
    pc = passivity_certificate(m)
    cr = real_part_crossings(m)
    z = port_zeros(m)
    nm = Porthos.name(c)
    ports[nm] = Dict("type" => model_type(c), "H0" => pc.H0, "crossings_rad_s" => cr,
                     "min_re_H" => pc.min_re, "at_rad_s" => pc.at, "nonpassive" => pc.nonpassive,
                     "shortage" => pc.shortage,
                     "zeros" => [[real(v), imag(v)] for v in z], "states" => m.states)
    @printf("  %-10s %-14s %10.4f %22s %12.4f  %s\n", nm, model_type(c), pc.H0,
            isempty(cr) ? "none" : join((@sprintf("%.6f", w) for w in cr), ", "), pc.shortage,
            join((@sprintf("%.3f%+.3fi", real(v), imag(v)) for v in z), ", "))
end

out = joinpath(root, "outputs", "ph_audit", splitext(basename(sc.system_path))[1])
path = Porthos.write_json(joinpath(out, "ph_audit.json"),
                          Dict("case" => sc.system_path, "equilibrium_residual" => eq.residual,
                               "shifted_storage" => a, "governor_ports" => ports))
println()
println("report -> ", path)
