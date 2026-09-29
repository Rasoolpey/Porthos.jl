# Certify a local region of attraction with the quadratic candidate V_P (roadmap P11;
# TODO.md step 3). One command:
#
#   julia --project=. scripts/certify_roa.jl [scenario.json]
#
# About 8 minutes (the coordinate-by-coordinate centered hull takes about 25 s per level).
# At Porthos's equilibrium (healthy network):
#   1. the section model (common-angle section, bus voltages explicit) and the candidate
#      V_P = xi' P xi, P from A'P + PA = -I at the equilibrium;
#   2. the certificate: positivity of P, the equilibrium enclosure (Krawczyk), and a
#      log-bisection over levels c, each gated by the KCL branch (parametric Krawczyk), the
#      decay (mean-value matrix in the centered hull, interval definiteness) and the
#      single-mode containment audit with the contract domain clauses, on the same boxes;
#      the first-order hull and the single-direction centered hull run first, for
#      comparison, and the latter's level starts the main search;
#   3. ROACheck on the written record: digests recomputed from the files, the level
#      re-proved with interval Cholesky.
# Records: outputs/roa/<case>/certificate_quadratic.json, roa_check.json.

using Porthos
using Printf

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
case = load_case(sc.system_path)
say("equilibrium of ", basename(dirname(sc.system_path)), "...")
eq = solve_equilibrium(case, sc)
m = section_model(eq)
say(string(m), "; scale 2^", Int(log2(minimum(m.scale))), " to 2^", Int(log2(maximum(m.scale))))
c = quadratic_candidate(m)
say(@sprintf("candidate: P condition %.3g, spectral abscissa %.4g", c.construction["P_condition"],
             c.construction["spectral_abscissa"]))

say("first-order hull (comparison):")
fo = certify_roa(c; level_low = 1e-16, level_high = 1e-8, bisections = 20, time_limit = 300,
                 hull = :first_order, resolution = 0.01)
say("centered hull, one interval direction (comparison):")
cd = certify_roa(c; level_low = 1e-14, level_high = 1e-8, bisections = 20, time_limit = 300,
                 second_order = :direction, resolution = 0.01)
say("centered hull, coordinate by coordinate (about 25 s per level):")
lo = something(cd["verified_valid_level"], 1e-14)
rec = certify_roa(c; level_low = lo, level_high = 1e-8, bisections = 12, time_limit = 600,
                  resolution = 0.01)
rec["scenario_path"] = scenario
rec["first_order_verified_valid_level"] = fo["verified_valid_level"]
rec["direction_verified_valid_level"] = cd["verified_valid_level"]

outdir = joinpath(root, "outputs", "roa", basename(dirname(sc.system_path)))
path = write_certificate(joinpath(outdir, "certificate_quadratic.json"), rec)
say("record: ", path)

say("ROACheck (independent recheck from the written record)...")
check = roa_check(path)
write_certificate(joinpath(outdir, "roa_check.json"), check)
for (k, v) in sort(collect(check["checks"]); by = first)
    say(@sprintf("   %-32s %s", k, v))
end
say("ROACheck ", check["passed"] ? "PASSED" : "FAILED", get(check, "decay_method", "") == "" ? "" :
    " (decay: " * check["decay_method"] * ")")

lvl = rec["verified_valid_level"]
println()
if lvl === nothing
    println("No level certified.")
else
    r = rec["certified_row"]
    @printf("verified_valid_level = %.4g (centered hull, coordinate by coordinate)\n", lvl)
    println("comparison: centered hull along one interval direction ",
            something(cd["verified_valid_level"], "none"), "; first-order hull ",
            something(fo["verified_valid_level"], "none"))
    @printf("on the certified set: |xi_i| <= %.3g (section coordinates), decay lambda_max <= %.3g\n",
            r["xi_half_width_max"], r["decay_matrix_lambda_max_upper"])
    println("claim: ", rec["claim"])
end
