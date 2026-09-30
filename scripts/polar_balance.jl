# The polar network balance (TODO.md step 5, the reviewer's step 1). One command:
#
#   julia --project=. scripts/polar_balance.jl [scenario.json]
#
# About 1 minute. At the equilibrium, random KCL-consistent states and the healthy part of
# the fault trajectory (IDA, voltages re-solved on the KCL branch):
#   - the lossless identity dU_B/dt = sum (P^B theta' + Q^B d ln|V|/dt) with
#     U_B = -1/2 sum B_ij (Vd_i Vd_j + Vq_i Vq_j): the polar pairs (P, theta), (Q, ln|V|);
#   - the conductance part against its closed form -Im(V'^H G V) (curl 2 G kron J);
#   - per GENROU, the polar supply at its internal EMF node (P_int phi' + Q_int d ln|E''|)/omega_b
#     (= omega times the transformer power, plus an omega' term) against the rotor exchange
#     w'QD di/dt: whether the network's energy flow can cancel the exchange.

using Porthos
using LinearAlgebra, Printf

function main(scenario)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys = eq.sys; s = sc.solver
    traj = simulate_ida(sys, vcat(eq.x, eq.V); dt = s.dt, duration = 3.0, log_dt = 0.05, rtol = s.rtol, atol = s.atol)
    window = isempty(sys.faults) ? nothing :
             (minimum(f.t_start for f in sys.faults), maximum(f.t_end for f in sys.faults))
    smp = power_audit_samples(sys, eq.x, eq.V; n_random = 8, trajectory = traj, fault_window = window)
    idl = 0.0; gcf = 0.0; grate = 0.0; urate = 0.0; svt = 0.0
    mm = Dict{String,Vector{Tuple{Float64,Float64}}}()
    for sm in smp
        r = polar_balance(sys, sm.x, sm.V; faults_on = sm.faults_on)
        sm.faults_on && continue      # the fault shunt is outside B and G of the Y-bus
        idl = max(idl, abs(r["lossless_identity_error"])); gcf = max(gcf, abs(r["conductance_closed_form_error"]))
        grate = max(grate, abs(r["conductance_rate"])); urate = max(urate, abs(r["Udot"]))
        for (n, d) in r["machines"]
            svt = max(svt, abs(d["supply_vs_transformer"]))
            push!(get!(mm, n, Tuple{Float64,Float64}[]), (d["exchange"], d["polar_supply"]))
        end
    end
    @printf("samples %d; lossless identity error %.2g (|dU_B/dt| up to %.3g); conductance closed-form error %.2g (rate up to %.3g)\n",
            length(smp), idl, urate, gcf, grate)
    @printf("polar supply = omega x transformer power + omega' i.psi''/omega_b: max error %.2g\n", svt)
    for (n, v) in sort(collect(mm); by = first)
        ex = first.(v); po = last.(v)
        @printf("  %-9s exchange rms %.3g, polar supply rms %.3g, mismatch rms %.3g, correlation %.3f\n", n,
                norm(ex) / sqrt(length(ex)), norm(po) / sqrt(length(po)), norm(ex - po) / sqrt(length(ex)),
                dot(ex .- sum(ex)/length(ex), po .- sum(po)/length(po)) / (norm(ex .- sum(ex)/length(ex)) * norm(po .- sum(po)/length(po))))
    end
    C = conductance_curl(sys)
    @printf("conductance curl: max |2 G kron J| %.3g; |G| max %.3g, |B| max %.3g\n", maximum(abs, C), maximum(abs, sys.G), maximum(abs, sys.B))
end
main(isempty(ARGS) ? joinpath(@__DIR__, "..", "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") : ARGS[1])
