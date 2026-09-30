# The storage completion H_ext = H_physical + H_controller + H_repair (TODO.md, method
# decision of 2026-09-30, steps 1 to 4). One command:
#
#   julia --project=. scripts/storage_completion.jl [scenario.json]
#
# About 15 minutes; every solve is capped and each stage prints. Solver: Hypatia (interior
# point), installed into Julia's default environment on first use, as JuMP.
#   1. Stage 1, on the integrable-loss equivalent (`integrable_loss_variant`: same
#      equilibrium, every network supply exact): the completion LMI with the physical Hessian
#      and the field cross-terms fixed, the own-unit controller blocks free, the swing block
#      and the recorded inter-unit blocks as repair, closed under `invariant_closure`, and the
#      per-load filter terms; a strict quadratic Lyapunov function is checked in interval
#      arithmetic (`verified_lyapunov`).
#   2. Stage 2, the same pattern on the real plant (conductances and ComplexLoad active part
#      restored), checked the same way.
#   3. The lift: the explicit nonlinear H_ext (`storage_completion`, `completion_energy`)
#      against the LMI matrix, and the size of each part.
# The inter-unit blocks below were found by a dual-guided greedy search on stage 1 (from the
# own-unit pattern and the swing block, three blocks per round ranked by the margin gradient
# Z_P - (A Z_R + Z_R A'), 15 rounds until t > 0; 2026-09-30).

for p in ("JuMP", "Hypatia")
    if Base.find_package(p) === nothing
        println("installing $p into Julia's default environment (first time only)...")
        run(`$(Base.julia_cmd()) --startup-file=no -e "using Pkg; Pkg.add(\"$p\")"`)
    end
end
using Porthos
using JuMP, Hypatia
using LinearAlgebra, Printf

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

const REPAIR_PAIRS = [
    ("GENROU_10", "GENROU_11"), ("IEEET1_10", "IEEET1_11"), ("GENROU_1", "IEEET1_10"), ("IEEET1_5", "IEEET1_6"),
    ("IEEET1_7", "IEEET1_8"), ("GENROU_11", "GENROU_8"), ("IEEET1_4", "IEEET1_6"), ("IEEET1_4", "IEEET1_5"),
    ("IEEET1_4", "IEEET1_8"), ("GENROU_11", "GENROU_9"), ("IEEET1_10", "IEEET1_9"), ("GENROU_10", "GENROU_9"),
    ("GENROU_1", "IEEET1_5"), ("GENROU_1", "IEEET1_6"), ("GENROU_11", "GENROU_2"), ("IEEET1_10", "IEEET1_7"),
    ("GENROU_10", "GENROU_8"), ("IEEET1_10", "IEEET1_4"), ("GENROU_3", "GENROU_8"), ("GENROU_11", "GENROU_3"),
    ("GENROU_10", "GENROU_3"), ("GENROU_10", "GENROU_4"), ("GENROU_4", "GENROU_8"), ("GENROU_8", "GENROU_8"),
    ("GENROU_11", "GENROU_7"), ("IEEET1_10", "IEEET1_2"), ("IEEET1_2", "IEEET1_5"), ("GENROU_10", "GENROU_10"),
    ("IEEET1_6", "IEEET1_9"), ("IEEET1_5", "IEEET1_9"), ("GENROU_4", "GENROU_7"), ("GENROU_3", "GENROU_4"),
    ("IEEET1_4", "IEEET1_7"), ("GENROU_1", "GENROU_11"), ("GENROU_11", "IEEET1_9"), ("GENROU_7", "GENROU_9"),
    ("GENROU_6", "GENROU_7"), ("GENROU_5", "GENROU_7"), ("GENROU_2", "GENROU_3")]

function stage(label, m, pairs, optimizer)
    say(label, ": the exact decomposition and the physical Hessian (about 75 s)...")
    prob = completion_problem(m)
    n = Porthos.neta(m)
    mask = completion_pattern(prob, pairs)
    A = prob.A
    T = pow2_scaling(diag(Porthos.lyap(Matrix(A'), Matrix(1.0I, n, n))))
    As = T \ A * T
    Hs = T * prob.Hfix * T
    Bs = [T * b.matrix * T for b in prob.basis]
    nn = [b.nonnegative for b in prob.basis]
    say(@sprintf("   free entries %d (own-unit controller %d, repair %d of %d), %d repair basis terms; solving (cap 600 s)...",
                 div(count(mask), 2), div(count(mask .& prob.own), 2), div(count(mask .& .!prob.own), 2),
                 div(n * (n + 1), 2), length(Bs)))
    r = structured_completion(As, Hs, mask, Bs; optimizer, nonnegative = nn)
    v = verified_lyapunov(r.P, As)
    say(@sprintf("   %s, margin t = %.4g; interval check: P_min >= %.3g, decay_min >= %.3g, Lyapunov %s; cond(P) %.3g",
                 r.record["termination_status"], r.t, v.P_min, v.decay_min, v.lyapunov, cond(Symmetric(r.P))))
    return prob, r, T, v
end

function main(scenario)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    optimizer = optimizer_with_attributes(Hypatia.Optimizer, "time_limit" => 600.0)
    sysI = integrable_loss_variant(eq.sys, eq.x, eq.V)
    say(@sprintf("integrable-loss equivalent: %d constant-power sinks (total P %.4g)",
                 count(c -> c isa ConstantPowerSink, sysI.comps), sum(c.p.P for c in sysI.comps if c isa ConstantPowerSink)))
    stage("1. stage 1 (integrable-loss equivalent)", section_model(sysI, eq.x, eq.V; scale = :none), REPAIR_PAIRS, optimizer)
    prob, r, T, v = stage("2. stage 2 (real plant)", section_model(eq; scale = :none), REPAIR_PAIRS, optimizer)
    v.lyapunov || (say("stage 2 did not certify; stopping"); return)
    say("3. the lift: explicit nonlinear H_ext (a Hessian, about 60 s)...")
    scm = storage_completion(prob, r; T = T)
    Ti = inv(T)
    hc = completion_hessian_check(scm, Ti * r.P * Ti)
    say(@sprintf("   Hessian of H_ext against the LMI matrix: relative %.2e; gradient %.2e", hc.relative, hc.gradient_norm))
    nf(X) = norm(T * X * T)
    Bp = sum(t.coefficient .* b.matrix for (t, b) in zip(scm.terms, prob.basis))
    say(@sprintf("   parts (scaled Frobenius): physical %.3g, controller %.3g, repair blocks %.3g, repair terms %.3g",
                 nf(prob.Hfix), nf(scm.Fc), nf(scm.Fr), nf(Bp)))
    say("done")
end
main(isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") : ARGS[1])
