# The joint strain-energy identity and its local decay (TODO.md step 5b). One command:
#
#   julia --project=. scripts/strain_energy_check.jl [scenario.json]
#
# About 4 minutes. Prints progress.
#   1. The strain metric M of every machine rotor (GENROU and GENSAL): M > 0, U_rot convex.
#   2. The identity dU_ext/dt = dissipation + field + Te delta' + residual terms, on the
#      lossless variant and on the case network, at the equilibrium, random KCL-consistent
#      states and the healthy part of the fault trajectory: identity error and the size of
#      each term (speed voltage, conductance, ComplexLoad active part and lag).
#   3. The joint candidate U_ext + sum omega_b H (omega - 1)^2 at the equilibrium: the inertia of
#      its Bregman Hessian on the quotient by block, and the local decay of the plant
#      (machines and loads, controller states frozen) with a quadratic gauge on the load
#      filters, attributed to the terms of the identity (`strain_decay_forms`).

using Porthos
using LinearAlgebra, Printf

const T_START = time()
say(s...) = (println(@sprintf("%6.1f s  ", time() - T_START), s...); flush(stdout))

function main(scenario)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys = eq.sys
    s = sc.solver
    say("1. strain metric per rotor")
    for c in sys.comps
        sr = strain_rotor(c)
        sr === nothing && continue
        g = strain_metric(sr)
        @printf("     %-10s %-14s M > 0 %-5s convex %-5s  min eig M %.3g, U %.3g; |M B_s + Psi'| %.1e\n",
                Porthos.name(c), Porthos.model_type(c), g.metric_positive, g.convex,
                eigmin(Symmetric(g.M)), eigmin(Symmetric(g.U)), g.residual)
    end

    say("2. identity along the flow (simulating the fault for the trajectory samples)")
    traj = simulate_ida(sys, vcat(eq.x, eq.V); dt = s.dt, duration = 3.0, log_dt = 0.05, rtol = s.rtol, atol = s.atol)
    window = isempty(sys.faults) ? nothing :
             (minimum(f.t_start for f in sys.faults), maximum(f.t_end for f in sys.faults))
    smp = filter(sm -> !sm.faults_on,
                 power_audit_samples(sys, eq.x, eq.V; n_random = 8, trajectory = traj, fault_window = window))
    sysL = lossless_variant(sys)
    terms = ("dissipation", "field", "angle", "saturation", "speed_rotor", "speed_network",
             "conductance", "active_load", "load_lag")
    for (label, S) in (("lossless variant", sysL), ("case network", sys))
        err = 0.0; skipped = 0; used = 0
        mx = Dict(k => 0.0 for k in (terms..., "rate", "paper_residual"))
        for sm in smp
            r = try
                V = S === sys ? sm.V : solve_network(S, sm.x, sm.V)
                strain_balance(S, sm.x, V)
            catch e
                e isa DomainError || rethrow()
                skipped += 1
                continue
            end
            used += 1
            err = max(err, abs(r["identity_error"]) / max(1.0, abs(r["rate"])))
            for k in keys(mx)
                mx[k] = max(mx[k], abs(r[k]))
            end
        end
        say(@sprintf("   %-16s %d samples (%d outside the ComplexLoad middle band); identity error %.1e (relative)",
                     label, used, skipped, err))
        println("     max |term|: ", join([@sprintf("%s %.3g", k, mx[k]) for k in ("rate", terms..., "paper_residual")], ", "))
    end

    say("3. joint candidate at the equilibrium (two quotient Hessians, about 2 minutes)")
    m = section_model(eq; scale = :none)
    r = strain_decay_forms(m)
    say(@sprintf("   S' (supply-shifted): decomposition closes to %.1e, gradient %.1e", r.closure, r.gradient_norm))
    kin(x) = sum(sys.omega_b * Porthos.param_value(Porthos.param_dict(c)["H"]) * (x[Porthos._state_range(sys, k)][2] - 1)^2
                 for (k, c) in enumerate(sys.comps) if strain_rotor(c) !== nothing)
    hb = quotient_hessian(m, x -> strain_energy(sys, x, kcl_solve(sys, x, eq.V)) + kin(x)).hessian
    C = hb - r.H                          # the exact scalar c*' V(eta): Bregman form minus S'
    A = r.A
    names = m.names
    isl(n) = occursin("CLOAD", n)
    ism(n) = occursin("GEN", n)
    P = findall(n -> isl(n) || ism(n), names)
    im = findall(ism, names)
    il = findall(isl, names)
    neg(X) = count(<(-1e-9), eigvals(Symmetric(X)))
    for (label, Hc) in (("Bregman form of U_ext + K", hb), ("S'", r.H))
        @printf("     %-26s negative directions: machines %d, loads %d, plant %d (of %d, %d, %d)\n", label,
                neg(Hc[im, im]), neg(Hc[il, il]), neg(Hc[P, P]), length(im), length(il), length(P))
    end
    DC = (C * A + A' * C) ./ 2
    cL = 5.0
    say(@sprintf("   plant decay, controller states frozen, load gauge %.0f dz^2/2 (rates in 1/s, against the candidate's Hessian):", cL))
    for alpha in (0.5, 1.0)
        Hg = r.H .+ alpha .* C
        Fg = Dict(k => copy(v) for (k, v) in r.forms)
        Fg["loss_curvature"] = alpha .* DC
        G = zeros(size(Hg))
        for i in il
            Hg[i, i] += cL
            ei = zeros(size(Hg, 1)); ei[i] = 1
            G .+= cL .* (ei * A[i, :]' .+ A[i, :] * ei') ./ 2
        end
        Fg["load_gauge"] = G
        Hp = Symmetric(Hg[P, P])
        isposdef(Hp) || (println("     alpha = $alpha: Hessian not positive on the plant"); continue)
        L = cholesky(Hp).L
        K = eigen(Symmetric(L \ sum(values(Fg))[P, P] / L'))
        @printf("     S' + %.1f c*'V: Hessian min %.3g; max rate %.3g, %d creating directions of %d\n",
                alpha, eigmin(Hp), K.values[end], count(>(1e-9), K.values), length(P))
        v = L' \ K.vectors[:, end]
        contrib = sort([(k, dot(v, X[P, P] * v)) for (k, X) in Fg]; by = p -> -abs(p[2]))
        println("       worst direction: ", join([@sprintf("%s %+.3g", k, c) for (k, c) in contrib[1:6]], ", "))
        println("       each term alone (max rate): ",
                join([@sprintf("%s %.3g", k, eigmax(Symmetric(L \ X[P, P] / L'))) for (k, X) in sort(collect(Fg); by = first)], ", "))
    end
    say("done")
end
main(isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") : ARGS[1])
