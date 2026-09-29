# The ROA certificate (roadmap P11): one level, the level search, and the record.
#
# Claim for a certified level c (certificate discipline, roadmap section 0; the exact wording
# is `certificate_claim`): on Omega_c = {(eta, V) : V_cand(eta - eta*) <= c, V = v(eta)} (the
# certified KCL branch) dV_cand/dt < 0 except at the equilibrium, in one smooth mode of every
# limiter and inside every contract domain. For a power system the statement is about the
# retained physical quotient (the section coordinates, modulo the common rotation): it is
# attracted to the enclosed equilibrium while the excluded one-way reservoirs and the
# monitor stay in the ranges over which they are proved not to feed back. The reservoirs
# themselves are not claimed to converge, so neither is the full DAE state.
#
# Proof of the decay (mean-value form). With xi = eta - eta*, h(eta*) = 0 and V_cand's
# gradient vanishing at eta*: h(eta* + xi) = M(xi) xi, M = int_0^1 J(eta* + t xi) dt in the
# Jacobian hull over the box E = [eta*] + [-w, w] (w from `sublevel_half_widths`), and
# grad V_cand = N(xi) xi, N in `gradient_matrix_hull`; so dV/dt = xi' sym(N'M) xi, negative
# for xi != 0 when every symmetric matrix in the interval enclosure of sym(N'M) is proved
# negative definite (`verified_max_eig` < 0).
#
# `verified_valid_level` is set only when, at the same level and on the same boxes (checked
# by digest): the candidate is positive, the equilibrium is enclosed, the KCL branch is
# unique, the decay is proved, and every containment category passes.

const ROA_SCHEMA_VERSION = "1.0"

"""
    certify_level(c::LyapunovCandidate, level, eq::EquilibriumEnclosure; contracts,
                  hull = :centered, second_order = :coordinates) -> Dict

All the gates at one level, on one set of boxes. `hull`: the enclosure of the mean-value
matrix, `:centered` (`centered_hull`, with `second_order` passed on) or `:first_order`
(`jacobian_hull`). Returns the row of the certificate record;
`row["certified"]` is the verdict.
"""
function certify_level(c::LyapunovCandidate, level::Real, eq::EquilibriumEnclosure;
                       contracts::ContractSet = default_contracts(), hull::Symbol = :centered,
                       second_order::Symbol = :coordinates)
    hull in (:centered, :first_order) || throw(ArgumentError("hull: :centered or :first_order"))
    m = candidate_model(c)
    t0 = time()
    row = Dict{String,Any}("level" => Float64(level), "decay" => false, "kcl_branch" => false,
                           "containment" => false, "certified" => false)
    w = sublevel_half_widths(c, level)
    Xi = [ival(-v, v) for v in w]
    E = eq.eta .+ Xi
    row["xi_half_width_max"] = maximum(w)
    row["xi_half_width_min"] = minimum(w)
    try
        br = enclose_kcl_branch(m, E, IA.mid.(eq.V))
        row["kcl_branch"] = true
        row["kcl_record"] = br.record
        row["V_box_width_max"] = maximum(IA.diam, br.V)
        M, hrec = hull === :centered ? centered_hull(m, eq, E, br.V, Xi; second_order, X = br.X) : jacobian_hull(m, E, br.V)
        row["hull"] = hrec
        N = gradient_matrix_hull(c, Xi)
        S = permutedims(N) * M
        S = (S .+ permutedims(S)) ./ 2
        wb = weyl_max_eig(S)
        lam = min(wb.bound, verified_max_eig(S))           # both rigorous
        row["decay_matrix_lambda_max_upper"] = lam
        row["weyl"] = Dict("midpoint_lambda_max_upper" => wb.midpoint, "radius_rho_upper" => wb.radius)
        row["decay_matrix_width_2norm_upper"] = sqrt(opnorm(IA.diam.(S), 1) * opnorm(IA.diam.(S), Inf))
        row["decay"] = lam < 0
        audit = containment_audit(m, E, br.V; contracts)
        row["containment"] = passed(audit)
        row["containment_categories"] = audit.passed
        row["containment_records"] = audit.records
        row["modes"] = audit.modes
        row["excluded_ranges"] = audit.excluded
        row["box_digest"] = box_digest(E, br.V)
        row["same_boxes"] = audit.box_digest == row["box_digest"]
        row["certified"] = row["decay"] && row["kcl_branch"] && row["containment"] && row["same_boxes"]
    catch e
        e isa ProofFailure || e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
        row["failure"] = sprint(showerror, e)
    end
    row["time_s"] = time() - t0
    return row
end

"""
    certify_roa(c::LyapunovCandidate; level_low = 1e-16, level_high = 1.0, bisections = 40,
                time_limit = 900, contracts, hull = :centered, second_order = :coordinates,
                resolution = 1e-3, io = stdout) -> Dict

The certificate record of `c`: the positivity proof, the equilibrium enclosure, then
`certify_level` at `level_high`, else at `level_low` and a log-bisection between (at most
`bisections` further levels, until the bracket is narrower than `resolution` in log10, and
no new level after `time_limit` seconds). Every tried level
is kept. `verified_valid_level` is the largest certified level found: a certified level that
depends on the search, not a proved maximum. Progress goes to `io` (flushed).
"""
function certify_roa(c::LyapunovCandidate; level_low::Real = 1e-16, level_high::Real = 1.0,
                     bisections::Integer = 40, time_limit::Real = 900.0,
                     contracts::ContractSet = default_contracts(), hull::Symbol = :centered,
                     second_order::Symbol = :coordinates, resolution::Real = 1e-3,
                     io::IO = stdout)
    0 < level_low <= level_high < Inf || throw(ArgumentError("need 0 < level_low <= level_high < Inf"))
    t0 = time()
    say(s) = (println(io, Printf.@sprintf("%6.1f s  ", time() - t0), s); flush(io))
    m = candidate_model(c)
    pos = positivity_proof(c)
    say("candidate positive: $(pos.proved) (lambda_min >= $(pos.lower_bound))")
    eq = enclose_equilibrium(m)
    say(Printf.@sprintf("equilibrium enclosed: eta width %.2e, V width %.2e",
                        eq.record["eta_width_max"], eq.record["V_width_max"]))
    rows = Dict{String,Any}[]
    function attempt(level)
        r = certify_level(c, level, eq; contracts, hull, second_order)
        push!(rows, r)
        say(Printf.@sprintf("level %.4e: %s (decay %s, KCL %s, containment %s; lambda_max <= %s) %.1f s",
                            level, r["certified"] ? "CERTIFIED" : "not certified", r["decay"],
                            r["kcl_branch"], r["containment"],
                            string(get(r, "decay_matrix_lambda_max_upper", "-")), r["time_s"]))
        haskey(r, "failure") && say("    " * r["failure"])
        return r["certified"]
    end
    best = nothing
    if pos.proved
        if attempt(level_high)
            best = level_high
        elseif level_low < level_high && attempt(level_low)
            best = level_low
            lo, hi = log10(level_low), log10(level_high)
            for _ in 1:bisections
                time() - t0 > time_limit && (say("time limit reached"); break)
                hi - lo < resolution && break
                mid = (lo + hi) / 2
                if attempt(10.0^mid)
                    best, lo = 10.0^mid, mid
                else
                    hi = mid
                end
            end
        end
    end
    best_row = best === nothing ? nothing : rows[findlast(r -> r["level"] == best && r["certified"], rows)]
    lvl(key) = best === nothing ? nothing : (best_row[key] ? best : nothing)
    return Dict{String,Any}(
        "schema_version" => ROA_SCHEMA_VERSION,
        "method" => "mean-value matrix enclosed by the $(hull == :centered ? "centered form (second-order jets) intersected with the first-order hull" : "first-order Jacobian hull") + interval definiteness (Weyl and Gershgorin bounds); KCL branch by the parametric Krawczyk test; single-mode containment audit on the same boxes",
        "hull" => string(hull),
        "second_order" => string(second_order),
        "claim" => best === nothing ? "no level certified" : certificate_claim(m),
        "assumptions" => model_assumptions(m),
        "excluded_ranges" => best === nothing ? nothing : get(best_row, "excluded_ranges", nothing),
        "candidate" => candidate_record(c),
        "model" => section_record(m),
        "positivity" => Dict("proved" => pos.proved, "lower_bound" => pos.lower_bound, "method" => pos.method),
        "equilibrium" => merge(Dict{String,Any}("eta_lo" => IA.inf.(eq.eta), "eta_hi" => IA.sup.(eq.eta),
                                                "V_lo" => IA.inf.(eq.V), "V_hi" => IA.sup.(eq.V)), eq.record),
        "search" => Dict("level_low" => level_low, "level_high" => level_high,
                         "bisections" => bisections, "time_limit_s" => time_limit,
                         "note" => "levels are nested: every certified level is a certificate; the search finds a large one, not the largest"),
        "derivative_verified_level" => lvl("decay"),
        "dae_limited_level" => lvl("kcl_branch"),
        "limiter_limited_level" => best === nothing ? nothing : (get(best_row["containment_categories"], "limiter_mode", true) ? best : nothing),
        "domain_limited_level" => best === nothing ? nothing : (get(best_row["containment_categories"], "domain", true) ? best : nothing),
        "verified_valid_level" => best,
        "certified_row" => best_row,
        "rows" => [Dict(k => v for (k, v) in r if !(k in ("containment_records", "modes"))) for r in rows],
        "software" => software_record(),
        "time_s" => time() - t0)
end

"""
    certificate_claim(m) -> String

The theorem a certified level proves, for the model class of `m`.
"""
certificate_claim(m::SectionModel) =
    "Retained physical quotient: the $(neta(m)) section coordinates eta (the physical states, " *
    "modulo the common rotation of every rotor angle and bus voltage phasor) with the bus " *
    "voltages V on the certified KCL branch v(eta). Let c = verified_valid_level and " *
    "Omega_c = {(eta, V): V_cand(eta - eta*) <= c, V = v(eta)}. For every initial state with " *
    "its retained part in Omega_c, every held state at its equilibrium value and every " *
    "excluded state (the one-way reservoirs and the delta_COI monitor) in its range in " *
    "excluded_ranges: as long as every excluded state stays in its range, the retained state " *
    "stays in Omega_c, which lies in one smooth mode of every limiter and inside every contract " *
    "domain, and dV_cand/dt < 0 there except at eta*; the retained dynamics do not depend on " *
    "the excluded states on their ranges. Hence the retained quotient converges to the enclosed " *
    "equilibrium eta* whenever the excluded states stay in their ranges for all time; an " *
    "excluded state drifts at most at its drift_rate_upper, so from r0 it stays in range for at " *
    "least min(r0 - lower, upper - r0) / drift_rate_upper. Not claimed: convergence of the " *
    "excluded states (the reservoirs are one-way accounts and in general settle elsewhere on an " *
    "equilibrium manifold), hence not attraction of the full DAE state; nor independence from " *
    "the excluded states outside their recorded ranges."

"""
    model_assumptions(m) -> Vector{String}

What the claim rests on beyond the proof steps, for the model class of `m`.
"""
model_assumptions(m::SectionModel) = [
    "equivariance of every model type under the common rotation, as declared by rotation_action (read from each model's equations); check_rotation_symmetry verifies the wiring, the injections and the absence of fixed-voltage buses against it",
    "held states at their equilibrium values (their derivative is proved exactly zero on the box)",
    "the model as evaluated with its Float64 parameters and constants, in exact real arithmetic"]

"""
    certify_roa(scenario::AbstractString; kwargs...) -> Dict

The quadratic candidate `V_P` at the equilibrium of a scenario's case (healthy network),
certified by `certify_roa(candidate; kwargs...)`; the record also names the scenario.
"""
function certify_roa(scenario::AbstractString; kwargs...)
    sc = load_scenario(scenario)
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    rec = certify_roa(quadratic_candidate(section_model(eq)); kwargs...)
    rec["scenario_path"] = abspath(scenario)
    return rec
end

"""
    section_record(m::SectionModel) -> Dict

The coordinates of a certificate: case and contract fingerprints, the reference point, the
kept states, the section and the scale.
"""
function section_record(m::SectionModel)
    sys = m.sys
    names = state_names(sys)
    p = m.projection
    return Dict{String,Any}(
        "fingerprint" => model_fingerprint(m),
        "system_digest" => system_digest(sys),
        "init_params" => m.init_params,
        "case_path" => sys.case.path,
        "case_sha256" => isfile(sys.case.path) ? bytes2hex(open(SHA.sha256, sys.case.path)) : nothing,
        "contracts_sha256" => bytes2hex(open(SHA.sha256, CONTRACTS_PATH)),
        "n_states" => sys.n_diff, "n_buses" => nbus(sys),
        "reservoir_states" => names[p.reservoir], "held_states" => names[p.held],
        "monitor_states" => names[p.monitor],
        "reference_angle" => names[m.keep[m.ref]],
        "section_coordinates" => m.names, "scale" => m.scale,
        "x0" => m.x0, "V0" => m.V0,
        "coordinates" => "x[keep[z_j]] = x0 + scale_j eta_j; reference angle = x0_ref - sum_j l_j scale_j eta_j / l_ref; other states at x0; field projected along the common rotation")
end

"""
    software_record() -> Dict

Julia, Porthos (version and git commit when the source is a git checkout) and the interval
library versions.
"""
function software_record()
    root = pkgdir(@__MODULE__)
    commit = try
        readchomp(Cmd(`git -C $root rev-parse HEAD`; ignorestatus = true))
    catch
        ""
    end
    dirty = try
        !isempty(readchomp(Cmd(`git -C $root status --porcelain -- src`; ignorestatus = true)))
    catch
        nothing
    end
    return Dict{String,Any}("julia" => string(VERSION), "porthos" => string(pkgversion(@__MODULE__)),
                            "porthos_commit" => isempty(commit) ? nothing : commit,
                            "porthos_src_modified" => dirty,
                            "IntervalArithmetic" => string(pkgversion(IA)),
                            "ForwardDiff" => string(pkgversion(ForwardDiff)))
end

"""
    write_certificate(path, record) -> path

Write a certificate record as JSON (non-finite numbers as strings).
"""
function write_certificate(path::AbstractString, record::AbstractDict)
    mkpath(dirname(path))
    write_json(path, _json_safe(record))
    return path
end

_json_safe(x) = x
_json_safe(x::AbstractFloat) = isfinite(x) ? x : string(x)
_json_safe(x::Symbol) = string(x)
_json_safe(d::AbstractDict) = Dict{String,Any}(string(k) => _json_safe(v) for (k, v) in d)
_json_safe(v::AbstractVector) = Any[_json_safe(x) for x in v]
_json_safe(t::Tuple) = Any[_json_safe(x) for x in t]
_json_safe(t::NamedTuple) = Dict{String,Any}(string(k) => _json_safe(v) for (k, v) in pairs(t))
