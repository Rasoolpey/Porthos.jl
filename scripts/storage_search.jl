# Route B (TODO.md, "Route B plan"): which couplings does a local quadratic storage need,
# keeping the dynamics? One command:
#
#   julia --project=. scripts/storage_search.jl [scenario.json]
#
# About 30 to 45 minutes; every solve is capped (`time_limit`) and each stage prints (and
# flushes) its label before it starts.
#
# At Porthos's equilibrium (healthy network), on the common-angle section in reference-angle
# coordinates (plain physical states; the reference angle is -l_z'z / l_ref). For a pattern S
# the decay margin is
#   gamma(S) = min { lambda_max(As'P + PAs) : P on S, P >= 0, tr P = 1 }
# (`decay_margin`): S carries a strict quadratic Lyapunov function exactly when gamma(S) < 0.
# Patterns are written on the physical coordinates and mapped to the section
# (`section_pattern`), so terms with the reference machine's own angle are kept.
#   0. references: the dense Lyapunov solution in the same normalisation, and the bound
#      gamma(full) >= 2 max Re eig(As);
#   1. the block-local patterns: "component" and "unit" (machine + governor + exciter), all
#      rotor angles coupled;
#   2. hypotheses, each the unit pattern plus one class of inter-unit terms;
#   3. greedy from the unit pattern: add the inter-component blocks the dual ranks highest
#      (`rank_couplings`), until the margin is substantial;
#   4. the fewest inter-unit terms (weighted L1) that keep half the margin of the first
#      substantial pattern (or, when none is, of the first verified Lyapunov pattern).
# Every margin comes with its dual certificate checked in interval arithmetic
# (`pattern_certificate`: a rigorous lower bound on gamma) and, for the primal P, a rigorous
# Lyapunov check (`verified_lyapunov`). Rigorous means: for the Float64 matrix As (the exact
# reduced Jacobian evaluated in floating point); the scaling T is by powers of two, so the
# scaled matrix is exact. The states are scaled (T = diag of the dense solution^(-1/2)); the
# normalisation tr P = 1, the L1 weights and the rankings act in these scaled coordinates,
# so magnitudes depend on the scaling (signs of gamma do not). Solver: Hypatia (interior
# point), installed into Julia's default environment on first use, as JuMP.
# Report: outputs/ph_audit/<case>/storage_search.json.

for p in ("JuMP", "Hypatia")
    if Base.find_package(p) === nothing
        println("installing $p into Julia's default environment (first time only)...")
        run(`$(Base.julia_cmd()) --startup-file=no -e "using Pkg; Pkg.add(\"$p\")"`)
    end
end
println("loading Porthos, JuMP, Hypatia...")
flush(stdout)
using Porthos
using JuMP, Hypatia
using LinearAlgebra
using Printf

const T_START = time()
elapsed() = @sprintf("%5.0f s", time() - T_START)
say(s...) = (println(elapsed(), "  ", s...); flush(stdout))

root = normpath(joinpath(@__DIR__, ".."))
scenario = isempty(ARGS) ? joinpath(root, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json") :
           abspath(ARGS[1])
sc = load_scenario(scenario)
case = load_case(sc.system_path)
eq = solve_equilibrium(case, sc)
sys, x, V = eq.sys, eq.x, eq.V
proj = physical_projection(sys, x, V)
groups = state_groups(sys, proj)
z, U, As, ref = reference_section(sys, x, V, proj)
zg = groups[z]
n = length(z)

time_limit = 300.0
optimizer = optimizer_with_attributes(Hypatia.Optimizer, "time_limit" => time_limit)

# ---- 0. references ----------------------------------------------------------------------
Q0 = lyap(Matrix(As'), Matrix(1.0I, n, n))
T = pow2_scaling(diag(Q0))
Ass = T \ As * T                                   # exact: T is a power of two per state
Qd = lyap(Matrix(Ass'), Matrix(1.0I, n, n))
Qd ./= tr(Qd)
gamma_dense = eigmax(Symmetric(Ass' * Qd + Qd * Ass))
abscissa = maximum(real, eigvals(As))
say(@sprintf("section: %d coordinates (reference angle %s.%s); spectral abscissa %.5g",
             n, groups[ref].component, groups[ref].state, abscissa))
say(@sprintf("scaling: condition of the dense solution %.2e -> %.2e; max |As| %.3g -> %.3g",
             cond(Symmetric(Q0)), cond(Symmetric(T * Q0 * T)), maximum(abs, As), maximum(abs, Ass)))
say(@sprintf("0. references: dense Lyapunov solution gamma = %.4g; gamma(full) >= 2 max Re eig = %.4g",
             gamma_dense, 2abscissa))
# a margin counts as substantial below 1% of the dense one
substantial = 0.01 * gamma_dense

finite(x) = x
finite(x::AbstractFloat) = isfinite(x) ? x : string(x)
finite(d::AbstractDict) = Dict{String,Any}(string(k) => finite(v) for (k, v) in d)
finite(v::AbstractVector) = map(finite, v)

report = Dict{String,Any}(
    "case" => sc.system_path, "scenario" => scenario, "n" => n,
    "reference_angle" => groups[ref].component * "." * groups[ref].state,
    "software" => Dict("julia" => string(VERSION), "jump" => string(pkgversion(JuMP)),
                       "hypatia" => string(pkgversion(Hypatia)),
                       "porthos" => string(pkgversion(Porthos)),
                       "note" => "JuMP and Hypatia are loaded from Julia's default environment, outside the project Manifest"),
    "time_limit_s" => time_limit,
    "normalisation" => "gamma(S) = min lambda_max(As'P + PAs) over P on S, P >= 0, tr P = 1, in the scaled coordinates (As -> T^-1 As T, T = diag(Q0)^(-1/2) rounded to powers of two, Q0 the dense Lyapunov solution); signs of gamma do not depend on the scaling, magnitudes, L1 weights and rankings do",
    "references" => Dict("gamma_dense" => gamma_dense, "gamma_full_lower_bound" => 2abscissa,
                         "spectral_abscissa" => abscissa, "substantial_below" => substantial),
    "stages" => Any[])

unit_phys = storage_pattern(groups, :unit)

function margin(label, mask_phys; stage)
    mz = section_pattern(mask_phys, U)
    free = (count(mz) + n) ÷ 2
    say(label, ": ", free, " free entries; solving (cap ", Int(time_limit), " s)...")
    r = decay_margin(Ass, mz; optimizer)
    cert = all(isfinite, r.Z) ? pattern_certificate(Ass, mz, r.Z, r.M) :
           (float_bound = NaN, verified_bound = -Inf, infeasible = false)
    lyap_ = all(isfinite, r.P) ? verified_lyapunov(r.P, Ass) : (P_min = NaN, decay_min = NaN, lyapunov = false)
    verdict = lyap_.lyapunov ? (r.gamma <= substantial ? "Lyapunov, substantial margin" : "Lyapunov, small margin") :
              cert.infeasible ? "no Lyapunov function (certified)" : "margin about 0 (boundary)"
    say(@sprintf("   gamma %.4g in [%.3g (rigorous), %.3g]; %s; %s, %.0f s, %s iterations",
                 r.gamma, cert.verified_bound, r.gamma, verdict, r.record["termination_status"],
                 something(r.record["solve_time_s"], NaN), something(r.record["iterations"], "?")))
    rec = Dict{String,Any}("label" => label, "stage" => stage, "free_entries" => free,
                           "gamma" => r.gamma, "gamma_relative_to_dense" => r.gamma / abs(gamma_dense),
                           "certificate" => Dict(pairs(cert)), "lyapunov" => Dict(pairs(lyap_)),
                           "verdict" => verdict, "solver" => r.record)
    push!(report["stages"], rec)
    return (; r, mz, rec)
end

# ---- 1. block-local patterns ------------------------------------------------------------
say("1. block-local patterns")
margin("component", storage_pattern(groups, :component); stage = 1)
unit = margin("unit", unit_phys; stage = 1)
first_lyap = nothing           # the first pattern with a verified Lyapunov function
function note_lyap!(label, mphys, res)
    if first_lyap === nothing && res.rec["lyapunov"][:lyapunov]
        global first_lyap = (label, copy(mphys), res)
    end
    return nothing
end

# ---- 2. hypotheses: the unit pattern plus one class of inter-unit terms -----------------
say("2. hypotheses")
is_machine(g) = g.role === :machine
hypotheses = [
    ("unit + speed-speed", g -> g.state == "omega", g -> g.state == "omega"),
    ("unit + speed-angle", g -> g.state == "omega", g -> g.angle),
    ("unit + machine-angle", is_machine, g -> g.angle),
    ("unit + machine-machine", is_machine, is_machine),
]
first_ok = nothing
for (label, f, g) in hypotheses
    res = margin(label, couple_states!(copy(unit_phys), groups, f, g); stage = 2)
    note_lyap!(label, couple_states!(copy(unit_phys), groups, f, g), res)
    if first_ok === nothing && res.rec["lyapunov"][:lyapunov] && res.r.gamma <= substantial
        global first_ok = (label, couple_states!(copy(unit_phys), groups, f, g), res)
    end
end

# ---- 3. greedy, ranked by the dual ------------------------------------------------------
say("3. greedy from the unit pattern (the dual ranks the missing inter-component blocks)")
mask = copy(unit_phys)
cur = unit
added = Tuple{String,String}[]
greedy = Any[]
for it in 1:10
    time() - T_START > 3600 && (say("   stopping: one hour spent"); break)
    G = Ass * cur.r.Z + cur.r.Z * Ass' - cur.r.M
    ranked = rank_couplings(G, cur.mz, zg; by = :component)
    top = ranked[1:min(4, length(ranked))]
    for t in top
        add_coupling!(mask, groups, t.a, t.b; by = :component)
        push!(added, (t.a, t.b))
    end
    say("   iteration $it adds: ", join(["$(t.a)-$(t.b) ($(@sprintf("%.3g", t.norm)))" for t in top], ", "))
    global cur = margin("greedy $it", mask; stage = 3)
    note_lyap!("greedy $it", mask, cur)
    push!(greedy, Dict("iteration" => it, "added" => [[t.a, t.b, t.norm] for t in top],
                       "gamma" => cur.r.gamma))
    if cur.rec["lyapunov"][:lyapunov] && cur.r.gamma <= substantial
        if first_ok === nothing
            global first_ok = ("greedy $it", copy(mask), cur)
        end
        break
    end
end
report["greedy"] = Dict("steps" => greedy, "added" => [[a, b] for (a, b) in added])

# ---- 4. fewest inter-unit terms -----------------------------------------------------------
target = first_ok === nothing ? first_lyap : first_ok
if target === nothing
    say("4. skipped: no pattern with a verified Lyapunov function was found")
    report["sparsest"] = "skipped: no pattern with a verified Lyapunov function"
else
    first_ok === nothing &&
        say("4. no substantial margin; using the first verified Lyapunov pattern")
    label, mphys, res = target
    rate = 0.5 * abs(res.r.gamma)
    mz = res.mz
    muz = section_pattern(unit_phys, U)
    W = [mz[i, j] && !muz[i, j] ? 1.0 : 0.0 for i in 1:n, j in 1:n]
    say("4. fewest inter-unit terms on '", label, "' keeping half its margin (rate ",
        @sprintf("%.3g", rate), "); solving (cap ", Int(time_limit), " s)...")
    s4 = decay_margin(Ass, mz; optimizer, weights = W, rate)
    P = s4.P
    ok = all(isfinite, P)
    v4 = ok ? verified_lyapunov(P, Ass) : (P_min = NaN, decay_min = NaN, lyapunov = false)
    blocks = Dict{String,Float64}()
    if ok
        scale = maximum(abs, P)
        for i in 1:n, j in (i + 1):n
            W[i, j] > 0 || continue
            a, b = sort([zg[i].component, zg[j].component])
            key = a * " - " * b
            blocks[key] = max(get(blocks, key, 0.0), abs(P[i, j]) / scale)
        end
    end
    needed = sort([(k, v) for (k, v) in blocks if v > 1e-6]; by = x -> -x[2])
    say(@sprintf("   %s; rigorous Lyapunov %s; %d of %d inter-unit component blocks used (> 1e-6 max|P|)",
                 s4.record["termination_status"], v4.lyapunov, length(needed), length(blocks)))
    for (k, v) in needed[1:min(25, end)]
        @printf("       %-30s %.3e\n", k, v)
    end
    flush(stdout)
    report["sparsest"] = Dict("pattern" => label, "rate" => rate, "solver" => s4.record,
                              "lyapunov" => Dict(pairs(v4)),
                              "inter_unit_blocks_used" => Dict(k => v for (k, v) in needed),
                              "inter_unit_blocks_available" => length(blocks))
end

out = joinpath(root, "outputs", "ph_audit", splitext(basename(sc.system_path))[1])
path = Porthos.write_json(joinpath(out, "storage_search.json"), finite(report))
say("report -> ", path)
