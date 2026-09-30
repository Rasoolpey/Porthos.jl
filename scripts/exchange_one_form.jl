# The GENROU stator exchange through the KCL branch (TODO.md step 5a). One command:
#
#   julia --project=. scripts/exchange_one_form.jl [scenario.json]
#
# About 1 minute. At the equilibrium of the scenario's case:
#   1. the current-corrected rotor coordinates w = z - D i: damper currents (Qw) and the
#      balance of field power against rotor loss;
#   2. the stator-exchange one-form a . d eta (a = (dI/d eta)' D'Qw, section coordinates):
#      invariance of the stator currents under the common rotation; the curl test on random
#      direction pairs, on the case's network and on its lossless variant; the closed-form
#      curl (dI/d eta)' D'Q (dz/d eta) - transpose against the measured one; a path test.
# An exact one-form would be the gradient of a network storage U_net cancelling the
# exchange; a nonzero curl means it is a non-integrable exchange / passivity shortage, not
# a scalar cross-term or network energy.

using Porthos
using LinearAlgebra, Printf, Random

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

function main(scenario)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    m = section_model(eq; scale = :none)
    n = Porthos.neta(m)
    rc = rotor_current_coordinates(eq.sys, eq.x, eq.V)
    say(@sprintf("1. %d linear rotors; at the equilibrium: damper currents max %.2g, field power - loss max %.2g",
                 length(rc), maximum(maximum(abs, r.Qw[2:end]) for r in values(rc)),
                 maximum(abs(r.field_supply - r.loss) for r in values(rc))))
    for (nm, r) in sort(collect(rc); by = first)
        @printf("     %-9s i_fd runtime %.5g, i_fd_energy %.5g\n", nm, r.i_fd_runtime, r.i_fd_energy)
    end
    rot = zeros(length(eq.x))
    for (k, i) in enumerate(m.keep)
        m.angle[k] && (rot[i] = 1.0)
    end
    th = 0.1
    Vr = vcat([[cos(th) * eq.V[2i-1] - sin(th) * eq.V[2i], sin(th) * eq.V[2i-1] + cos(th) * eq.V[2i]]
               for i in 1:length(eq.V)÷2]...)
    Vr = solve_network(eq.sys, eq.x .+ th .* rot, Vr)
    dI = maximum(abs, Porthos._stator_currents(eq.sys, eq.x .+ th .* rot, Vr) .-
                      Porthos._stator_currents(eq.sys, eq.x, eq.V))
    say(@sprintf("2. stator currents under a common rotation by %.1f rad: change %.2g (the one-form lives on the quotient)", th, dI))
    sysL = lossless_variant(eq.sys)
    VL = solve_network(sysL, eq.x, eq.V)
    for (label, sys, V0) in (("case network", eq.sys, eq.V), ("lossless variant", sysL, VL))
        r = one_form_exactness(m, zeros(n); sys, V0, pairs = 4)
        W = exchange_curl(m, zeros(n); sys, V0)
        rng = Random.Xoshiro(3)                  # the pairs of one_form_exactness (same seed)
        closed = map(r.pairs) do p
            u = normalize(randn(rng, n))
            v = normalize(randn(rng, n))
            abs((p.vJu - p.uJv) - dot(v, W * u))
        end
        pt = one_form_path_test(m, 1e-2 .* randn(Random.Xoshiro(4), n); sys, V0)
        say(@sprintf("   %-16s curl %.3g (relative, 4 pairs); closed form matches to %.1e; path difference %.2g of %.3g",
                     label, r.curl, maximum(closed), pt.difference, pt.straight))
    end
    say("=> not exact: the exchange is a non-integrable passivity shortage, not the gradient of a scalar storage")
end

main(isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") : ARGS[1])
