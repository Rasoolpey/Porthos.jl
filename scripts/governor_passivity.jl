# Can the rest of the system pay for each governor's passivity shortage? One command:
#
#   julia --project=. scripts/governor_passivity.jl [scenario.json]
#
# For each governor g, at Porthos's equilibrium of the scenario's case (healthy network):
#   - the governor at its port: omega -> Tm, H_gov = -G_gov; its shortage
#     nu = -min_w Re H_gov(jw) (P10 audit);
#   - the rest of the system seen from the governor (loop opened at g, exact): input Tm,
#     output omega of g's machine, G_rest; its excess rho = min_w Re(1/G_rest(jw)) (output
#     feedback passivity index), over the grid and the limit w -> infinity,
#     Re(1/G_rest) -> -c'Ab / (c'b)^2 (relative degree one: the rotor inertia);
#   - check: closing the two models again reproduces the full system's eigenvalues;
#   - index test: rho >= nu means a joint storage exists with the classical split (rest
#     dissipative with supply Tm*omega - rho*omega^2, governor with -Tm*omega + nu*omega^2);
#   - frequency-wise margin: min_w [Re(1/G_rest(jw)) + Re H_gov(jw)] > 0 means the loop at
#     this port is strictly positive real, so a joint storage exists, but its split between
#     the two sides then needs a frequency-dependent multiplier.
# Nothing in the model changes; this only measures how the existing interconnection shares
# its damping. Written to outputs/ph_audit/<case>/governor_passivity.json.

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

# the full system on the section, for the reconstruction check
A = reduced_jacobian(sys, x, V)[proj.keep, proj.keep]
Afull = proj.basis' * A * proj.basis
lam_full = eigvals(Afull)
@printf("full system: %d section states, max Re(eig) = %.4g\n", size(Afull, 1), maximum(real, lam_full))

ws = exp10.(range(-3, 3; length = 60001))
rows = Dict{String,Any}()
println()
@printf("%-10s %-14s %9s %10s %10s %9s  %-22s %10s %s\n", "governor", "type", "nu (gov)",
        "rho (rest)", "rho - nu", "index", "freq-wise margin", "at rad/s", "recon")
for (k, c) in enumerate(sys.comps)
    model_type(c) in ("IEEEG1_PHTRUE", "IEEEG3_PHTRUE") || continue
    nm = Porthos.name(c)
    r = sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
    mg = port_model(c, x[r], ins[k]; input = "omega", output = "Tm")
    mr = loop_port_model(sys, x, V, k; input = "omega", output = "Tm", projection = proj)

    # reconstruction: rest (input Tm, output omega) closed with the governor (omega -> Tm)
    Acl = [mr.A mr.B * mg.C'; mg.B * mr.C' mg.A]
    lam_cl = eigvals(Acl)
    recon = maximum(minimum(abs.(lam_cl .- l)) for l in lam_full) / max(1.0, maximum(abs, lam_full))
    stable_rest = maximum(real, eigvals(mr.A))

    Hg = -frequency_response(mg, ws)
    Gr = frequency_response(mr, ws)
    Zr = 1 ./ Gr
    cb = dot(mr.C, mr.B)
    zinf = -dot(mr.C, mr.A * mr.B) / cb^2              # Re(1/G_rest) as w -> infinity
    nu = max(0.0, -minimum(real, Hg))
    rho = min(minimum(real, Zr), zinf)
    margin = real.(Zr) .+ real.(Hg)
    km = argmin(margin)
    minmargin, atm = margin[km] < zinf ? (margin[km], ws[km]) : (zinf, Inf)
    nearzeros = [v for v in port_zeros(mr) if real(v) > -0.1]
    @printf("%-10s %-14s %9.4f %10.4f %10.4f %9s  %+22.4f %10.4f %.1e\n", nm, model_type(c), nu, rho,
            rho - nu, rho >= nu ? "passes" : "fails", minmargin, atm, recon)
    rows[nm] = Dict("type" => model_type(c), "nu_governor" => nu, "rho_rest" => rho,
                    "rho_rest_limit_infinity" => zinf, "index_test" => rho >= nu,
                    "freqwise_margin" => minmargin,
                    "freqwise_margin_at_rad_s" => isfinite(atm) ? atm : "infinity",
                    "rest_zeros_near_axis" => [[real(v), imag(v)] for v in nearzeros],
                    "rest_max_real_eig" => stable_rest, "G_rest_0" => real(transfer(mr, 0.0)),
                    "H_gov_0" => real(-transfer(mg, 0.0)), "reconstruction_error" => recon)
end
out = joinpath(root, "outputs", "ph_audit", splitext(basename(sc.system_path))[1])
path = Porthos.write_json(joinpath(out, "governor_passivity.json"),
                          Dict("case" => sc.system_path, "grid_rad_s" => [ws[1], ws[end], length(ws)],
                               "governors" => rows))
println()
println("report -> ", path)
