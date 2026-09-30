# Joint machine storage in current-corrected rotor coordinates (TODO.md step 5a, the
# reviewer's calculations 1 to 3). One command:
#
#   julia --project=. scripts/joint_storage_check.jl [scenario.json]
#
# About 4 minutes (each quotient Hessian differentiates twice through the KCL branch). At
# the equilibrium of the scenario's case, on the common-angle section:
#   1-2. the Hessians of the Bregman forms of the gauge family H_alpha = H_w - alpha F
#        (alpha = 0: sum w'Qw/2; alpha = 1: minus F = sum I'D'QD I/2), alone and with the
#        shifted kinetic storage sum H (omega - 1)^2: inertia and most negative eigenvalues;
#   3. the local exchange shortage -W'QD J_I A against the rotor loss W'LW and the
#      incremental network conductance loss J_V'G J_V: the exchange on the null space of the
#      dissipation, and its largest generalized eigenvalue against the dissipation.
# Controllers carry no storage yet, so their dissipation is not in the comparison.

using Porthos
using LinearAlgebra, Printf

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

function main(scenario)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    m = section_model(eq; scale = :none)
    sys, V0 = eq.sys, eq.V
    say("section: ", Porthos.neta(m), " coordinates; quotient Hessians (about 50 s each)...")
    for (label, S) in (("H_w (alpha = 0)", x -> rotor_energy(sys, x, V0; alpha = 0.0)),
                       ("H_w - F (alpha = 1)", x -> rotor_energy(sys, x, V0; alpha = 1.0)),
                       ("kinetic + H_w", x -> shifted_kinetic_energy(sys, x) + rotor_energy(sys, x, V0; alpha = 0.0)),
                       ("kinetic + H_w - F", x -> shifted_kinetic_energy(sys, x) + rotor_energy(sys, x, V0; alpha = 1.0)))
        r = quotient_hessian(m, S)
        say(@sprintf("  %-20s inertia +%d / -%d / 0: %d; most negative %s", label, r.positive, r.negative,
                     r.zero, join((@sprintf("%.3g", v) for v in r.eigenvalues[1:min(3, r.negative)]), ", ")))
    end
    f = exchange_dissipation_forms(m)
    say(@sprintf("3. exchange form eig [%.3g, %.3g]; rotor loss eig max %.3g; network loss eig max %.3g",
                 f.exchange_eig..., f.rotor_loss_eig_max, f.network_loss_eig_max))
    say(@sprintf("   on the null space of rotor + network dissipation the exchange is %.2g (vanishes);", f.null_exchange_max))
    say(@sprintf("   on its range the exchange exceeds that dissipation by up to %.4g times", f.ratio_max))
end

main(isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") : ARGS[1])
