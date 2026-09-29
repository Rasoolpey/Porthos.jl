# Results: simulation_results.csv with PHPS's columns, a binary copy, and run.json
# (roadmap 2.4 and principle 6), and the `simulate` entry point.

"""
    csv_columns(sys) -> Vector{String}

PHPS's `simulation_results.csv` header: `t`; for each component its states then its
observables (`"COMP.name"`); then `Vd_Bus<id>`, `Vq_Bus<id>`, `Vterm_Bus<id>` per bus.
`delta_COI` is not written, as in PHPS.
"""
function csv_columns(sys::DAESystem)
    cols = ["t"]
    for c in sys.comps
        for s in state_names(c)
            push!(cols, name(c) * "." * s)
        end
        for o in observable_names(c)
            push!(cols, name(c) * "." * o)
        end
    end
    for b in sys.net.bus_ids
        push!(cols, "Vd_Bus$b", "Vq_Bus$b", "Vterm_Bus$b")
    end
    return cols
end

"""
    csv_row(sys, t, y) -> Vector{Float64}

One row of `simulation_results.csv` at time `t` and state `y = [x; V]`.
"""
function csv_row(sys::DAESystem, t::Real, y::AbstractVector)
    nd = sys.n_diff
    x = view(y, 1:nd)
    V = view(y, nd + 1:length(y))
    ins, outs = component_io(sys, x, V)
    row = Float64[t]
    for (k, c) in enumerate(sys.comps)
        r = sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
        append!(row, x[r])
        obs = zeros(length(observable_names(c)))
        observable_values!(obs, c, view(x, r), ins[k], outs[k])
        append!(row, obs)
    end
    for i in 1:nbus(sys)
        Vd, Vq = V[2i - 1], V[2i]
        push!(row, Vd, Vq, sqrt(Vd * Vd + Vq * Vq))
    end
    return row
end

"""
    write_results_csv(path, r::SimResult)

`simulation_results.csv` in PHPS's format (`%.15e`, comma separated, header line).
"""
function write_results_csv(path::AbstractString, r::SimResult)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, join(csv_columns(r.sys), ","))
        for k in eachindex(r.t)
            row = csv_row(r.sys, r.t[k], view(r.Y, :, k))
            for (j, v) in enumerate(row)
                j > 1 && write(io, ',')
                Printf.@printf(io, "%.15e", v)
            end
            write(io, '\n')
        end
    end
    return path
end

"""
    write_results_jld2(path, r::SimResult)

The logged states in binary: `t`, `Y` (one column per time) and the row `names` (the
states, then `Vd_Bus<id>`, `Vq_Bus<id>` per bus).
"""
function write_results_jld2(path::AbstractString, r::SimResult)
    mkpath(dirname(abspath(path)))
    JLD2.jldsave(path; t = r.t, Y = r.Y,
                 names = [state_names(r.sys); ["V$(q)_Bus$b" for b in r.sys.net.bus_ids
                                               for q in ("d", "q")]])
    return path
end

_sha256_hex(path::AbstractString) = isfile(path) ? _sha256_file(path) : nothing

function _git_state(root::AbstractString)
    commit = try
        strip(read(`git -C $root rev-parse HEAD`, String))
    catch
        nothing
    end
    dirty = try
        !isempty(strip(read(`git -C $root status --porcelain --untracked-files=no`, String)))
    catch
        nothing
    end
    return commit, dirty
end

const PORTHOS_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    run_metadata(r, scenario; extra...) -> Dict

What `run.json` records: hashes of the case, scenario, source (git commit and a dirty flag)
and Manifest, the Julia and Porthos versions, the solver settings the integrator used
(`r.settings`, then `settings`) and statistics, the wall time.
"""
function run_metadata(r::SimResult, sc::Scenario; settings::AbstractDict = Dict())
    commit, dirty = _git_state(PORTHOS_ROOT)
    return Dict{String,Any}(
        "porthos" => Dict("commit" => commit, "dirty" => dirty,
                          "manifest_sha256" => _sha256_hex(joinpath(PORTHOS_ROOT, "Manifest.toml"))),
        "julia" => string(VERSION),
        "case" => Dict("path" => r.sys.case.path, "sha256" => _sha256_hex(r.sys.case.path)),
        "scenario" => Dict("path" => sc.path, "sha256" => _sha256_hex(sc.path)),
        "solver" => merge(Dict{String,Any}("method" => string(r.method)), r.settings, settings),
        "result" => Dict("records" => length(r.t), "t_end" => isempty(r.t) ? 0.0 : r.t[end],
                         "nonconverged_steps" => r.nonconverged,
                         "first_nonconverged_t" => r.first_nonconverged_t,
                         "stopped_early" => r.stopped_early),
        "wall_time_s" => r.wall_time,
        "created_utc" => Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ"),
    )
end

"""
    simulate(scenario_path; outdir = nothing, method = nothing) -> String

Run a scenario: load the case, solve the equilibrium, integrate with the scenario's solver
(`bdf1` or `ida`; `method` overrides it), and write `simulation_results.csv`,
`simulation_results.jld2` and `run.json` into `outdir` (default: `outputs/<scenario
output directory name>`). Returns the path of `run.json`.
"""
function simulate(scenario_path::AbstractString; outdir = nothing, method = nothing)
    sc = load_scenario(scenario_path)
    case = load_case(sc.system_path)
    m = method === nothing ? sc.solver.method : Symbol(method)
    eq = solve_equilibrium(case, sc)
    y0 = vcat(eq.x, eq.V)
    s = sc.solver
    log_dt = s.log_dt === nothing ? s.dt : s.log_dt
    r = if m === :bdf1
        simulate_bdf1(eq.sys, y0; dt = s.dt, duration = s.duration, log_dt)
    elseif m === :ida
        simulate_ida(eq.sys, y0; dt = s.dt, duration = s.duration, log_dt, rtol = s.rtol,
                     atol = s.atol)
    else
        throw(ArgumentError("solver method $m is not supported (PHPS's DAE path has bdf1 and " *
                            "ida; rk4 is an explicit ODE method, see TODO.md)"))
    end
    if outdir === nothing
        sub = sc.output_dir === nothing ? splitext(basename(sc.path))[1] : basename(sc.output_dir)
        outdir = joinpath(PORTHOS_ROOT, "outputs", sub)
    end
    mkpath(outdir)
    write_results_csv(joinpath(outdir, "simulation_results.csv"), r)
    write_results_jld2(joinpath(outdir, "simulation_results.jld2"), r)
    meta = run_metadata(r, sc; settings = Dict{String,Any}(
        "equilibrium_residual" => eq.residual))
    path = joinpath(outdir, "run.json")
    write_json(path, meta)
    return path
end
