# Joint storage for the governors that pass the passivity test. One command:
#
#   julia --project=. scripts/joint_governor_storage.jl [scenario.json]
#
# At Porthos's equilibrium of the scenario's case (healthy network):
#  1. per governor g: its shortage nu_g and, with the loop opened at g alone, the rest's
#     excess rho_g (as scripts/governor_passivity.jl); the governors with rho_g > nu_g are
#     taken;
#  2. their loops are opened at once: the rest is an m-port (Tm in, speed out). The
#     constant-split test is min_w eig(Herm(Z(jw)) - diag(n)) > 0 with n_g = nu_g + 2 delta_g;
#  3. storages, from Riccati equations (no SDP): each governor, dV_g/dt <= -Tm w + d_g w^2
#     - eta Tm^2 - strict |x_g|^2 (d_g = nu_g + delta_g); the rest, dV_r/dt <= sum Tm w -
#     n_g w^2 + eta Tm^2 - strict |z|^2. The supplies cancel in the sum, leaving
#     -sum delta_g w_g^2 - strict |.|^2 < 0;
#  4. check, as the P10 audit judges a storage: P = blockdiag(P_r, P_g...) is positive
#     definite and sym(P A) is negative definite on the common-angle section (A the exact
#     reduced Jacobian). Then V = z'Pz is a local Lyapunov function of the unchanged model
#     whose governor blocks are the governors' own storages.
# Written to outputs/ph_audit/<case>/joint_governor_storage.json.

println("loading Porthos...")
flush(stdout)
using Porthos
using LinearAlgebra
using Printf

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
case = load_case(sc.system_path)
eq = solve_equilibrium(case, sc)
sys, x, V = eq.sys, eq.x, eq.V
proj = physical_projection(sys, x, V)
ins, _ = component_io(sys, x, V)
ws = exp10.(range(-3, 3; length = 20001))

# 1. single-loop indices
govs = [k for (k, c) in enumerate(sys.comps) if model_type(c) in ("IEEEG1_PHTRUE", "IEEEG3_PHTRUE")]
nu, rho = Dict{Int,Float64}(), Dict{Int,Float64}()
for k in govs
    c = sys.comps[k]
    r = sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
    mg = port_model(c, x[r], ins[k]; input = "omega", output = "Tm")
    mr = loop_port_model(sys, x, V, k; input = "omega", output = "Tm", projection = proj)
    nu[k] = max(0.0, -minimum(real, -frequency_response(mg, ws)))
    Z = 1 ./ frequency_response(mr, ws)
    rho[k] = min(minimum(real, Z), -dot(mr.C, mr.A * mr.B) / dot(mr.C, mr.B)^2)
end
ks = [k for k in govs if rho[k] > nu[k]]
names = [Porthos.name(sys.comps[k]) for k in ks]
println("governors taken (rho > nu): ", join(names, ", "))
println("left inside the rest:       ", join([Porthos.name(sys.comps[k]) for k in govs if !(k in ks)], ", "))

# 2. all their loops opened at once
rest, ports, T = open_loops_model(sys, x, V, ks; input = "omega", output = "Tm", projection = proj)
delta = [0.1 * (rho[k] - nu[k]) for k in ks]
d = [nu[k] for k in ks] .+ delta
n = [nu[k] for k in ks] .+ 2 .* delta
# eta: the smallest of a few values for which both sides pass (it lets the governors pay for
# the differential speed directions that the lossy network cannot dissipate at low frequency)
mm, eta = nothing, NaN
for e in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2)
    r = multiport_margin(rest, n; eta = e, ws)
    g = [port_margin(m, dk; eta = e, ws) for (m, dk) in zip(ports, d)]
    @printf("eta = %-6g rest margin %+.4e at %.4f rad/s; governors' worst margin %+.4e  %s
", e, r.margin, r.at,
            minimum(x.margin for x in g), r.passes && all(x.passes for x in g) ? "<- both pass" : "")
    if r.passes && all(x.passes for x in g) && isnan(eta)
        global mm, eta = r, e
    end
end
lm = loop_margin(rest, ports; ws)
@printf("frequency-wise loop margin (each governor pays its own Re H(jw)): %+.4e at %.4f rad/s  (%s)
",
        lm.margin, lm.at, lm.passes ? "passes" : "fails")
out = joinpath(root, "outputs", "ph_audit", splitext(basename(sc.system_path))[1])
if isnan(eta)
    println("
=> no joint storage split port by port at these governors' Tm ports",
            lm.passes ? " with a constant split (a frequency-dependent split exists)" :
            ", constant or frequency-dependent: the cut must be elsewhere")
    path = Porthos.write_json(joinpath(out, "joint_governor_storage.json"), Dict(
        "case" => sc.system_path, "governors" => names, "nu" => [nu[k] for k in ks],
        "rho_single_loop" => [rho[k] for k in ks], "d" => d, "n" => n,
        "constant_split" => "fails for every eta scanned",
        "loop_margin" => Dict("margin" => lm.margin, "at_rad_s" => lm.at, "passes" => lm.passes)))
    println("report -> ", path)
    exit(0)
end
println("using eta = ", eta)

# 3. storages
strict = 1e-6
Pr, resr = rest_storage(rest, n; eta, strict)
Pgs = [port_storage(m, dk; eta, strict) for (m, dk) in zip(ports, d)]
@printf("rest storage: %d x %d, min eig %.3e, KYP residual %.2e\n", size(Pr)..., eigmin(Symmetric(Pr)), resr)
for (nm, (Pg, res)) in zip(names, Pgs)
    @printf("  %-10s storage %d x %d, eig in [%.3e, %.3e], KYP residual %.2e\n", nm, size(Pg)...,
            eigmin(Symmetric(Pg)), eigmax(Symmetric(Pg)), res)
end

# 4. the joint storage against the unchanged dynamics
A = reduced_jacobian(sys, x, V)[proj.keep, proj.keep]
Aw = T' * A * T
Pw = cat(Pr, (P for (P, _) in Pgs)...; dims = (1, 2))
epP = eigvals(Symmetric(Pw))
M = Symmetric(Pw * Aw + Aw' * Pw)
eM = eigvals(M)
@printf("\njoint storage on the section (%d coordinates): eig(P) in [%.3e, %.3e]; eig(sym(P A)) in [%.3e, %.3e]\n",
        size(Pw, 1), first(epP), last(epP), first(eM), last(eM))
ok = first(epP) > 0 && last(eM) < 0
println(ok ? "=> a local Lyapunov function of the unchanged model, with the governors' own storages" :
             "=> not a Lyapunov function (see the eigenvalues above)")

path = Porthos.write_json(joinpath(out, "joint_governor_storage.json"), Dict(
    "case" => sc.system_path, "governors" => names,
    "nu" => [nu[k] for k in ks], "rho_single_loop" => [rho[k] for k in ks], "d" => d, "n" => n,
    "eta" => eta, "strict" => strict,
    "multiport" => Dict("margin" => mm.margin, "at_rad_s" => mm.at),
    "rest_storage" => Dict("n" => size(Pr, 1), "min_eig" => eigmin(Symmetric(Pr)), "kyp_residual" => resr),
    "governor_storages" => Dict(nm => Dict("P" => [P[i, :] for i in 1:size(P, 1)], "states" => m.states,
                                          "kyp_residual" => res)
                                for (nm, (P, res), m) in zip(names, Pgs, ports)),
    "joint" => Dict("n" => size(Pw, 1), "P_min_eig" => first(epP), "P_max_eig" => last(epP),
                    "symPA_min_eig" => first(eM), "symPA_max_eig" => last(eM), "lyapunov" => ok)))
println("report -> ", path)
