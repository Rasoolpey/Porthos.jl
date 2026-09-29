# PowerFactory interface: run the `pf/` driver (Python, PowerFactory's own interpreter) and
# read what it writes. The driver only operates PowerFactory; everything computed from its
# results is computed here. See pf/README.md.

"""
    PFResults

A PowerFactory RMS run written by `pf/run.py simulate`: times `t` (s), the result matrix
`data` (one row per time, one column per recorded variable), the column index by
`(object, variable)`, and the run record `meta` (`run.json`). PowerFactory writes two rows
at the time of each event (before and after it).
"""
struct PFResults
    dir::String
    meta::Dict{String,Any}
    t::Vector{Float64}
    data::Matrix{Float64}
    index::Dict{Tuple{String,String},Int}
end

Base.show(io::IO, r::PFResults) =
    print(io, "PFResults(", basename(r.dir), ": ", length(r.t), " rows, ", size(r.data, 2),
          " signals, t = ", isempty(r.t) ? 0.0 : r.t[1], " .. ", isempty(r.t) ? 0.0 : r.t[end], ")")

_plain(x::JSON3.Object) = Dict{String,Any}(String(k) => _plain(v) for (k, v) in x)
_plain(x::JSON3.Array) = Any[_plain(v) for v in x]
_plain(x) = x

"""
    read_pf_results(dir) -> PFResults

Read `pf_results.csv` and `run.json` from `dir`. The column map in `run.json` comes from the
result object itself; the object names in the CSV's first header row are checked against it.
"""
function read_pf_results(dir::AbstractString)
    meta = _plain(read_json(joinpath(dir, "run.json")))
    cols = meta["columns"]
    lines = readlines(joinpath(dir, meta["csv"]); keep = false)
    objs = split(lines[1], ',')
    length(objs) == length(cols) + 1 ||
        error("$(meta["csv"]) has $(length(objs) - 1) signal columns, run.json lists $(length(cols))")
    index = Dict{Tuple{String,String},Int}()
    for (k, c) in enumerate(cols)
        c["column"] == k - 1 || error("run.json column $(k - 1) is out of order")
        strip(objs[k + 1], '"') == c["object"] ||
            error("column $k of $(meta["csv"]) is $(objs[k + 1]), run.json says $(c["object"])")
        index[(c["object"], c["variable"])] = k
    end
    body = lines[3:end]
    data = Matrix{Float64}(undef, length(body), length(cols) + 1)
    for (i, l) in enumerate(body)
        for (j, v) in enumerate(eachsplit(l, ','))
            data[i, j] = parse(Float64, v)
        end
    end
    return PFResults(String(dir), meta, data[:, 1], data[:, 2:end], index)
end

"""The recorded values of `variable` of PowerFactory object `object`."""
function pf_signal(r::PFResults, object::AbstractString, variable::AbstractString)
    k = get(r.index, (String(object), String(variable)), 0)
    k == 0 && throw(KeyError((object, variable)))
    return r.data[:, k]
end

function _pf_config()
    cfg = _plain(read_json(joinpath(PORTHOS_ROOT, "pf", "config.json")))
    local_ = joinpath(PORTHOS_ROOT, "pf", "config.local.json")
    isfile(local_) && merge!(cfg, _plain(read_json(local_)))
    return cfg
end

"""
    pf_command(args...)

Run `pf/run.py args...` with PowerFactory's Python (`py -<python_version>`). PowerFactory
starts as an engine, so its window must be closed.
"""
function pf_command(args...)
    ver = _pf_config()["python_version"]
    cmd = `py -$ver $(joinpath(PORTHOS_ROOT, "pf", "run.py")) $(collect(String.(args)))`
    run(Cmd(cmd; dir = PORTHOS_ROOT))
    return nothing
end

"""
    pf_simulate(scenario_path; outdir = outputs/pf/<scenario name>, dt_ms = config) -> PFResults

Run the scenario's bus faults in PowerFactory (RMS) and read the results.
"""
function pf_simulate(scenario_path::AbstractString; outdir = nothing, dt_ms = nothing)
    name = splitext(basename(scenario_path))[1]
    out = outdir === nothing ? joinpath(PORTHOS_ROOT, "outputs", "pf", name) : String(outdir)
    args = ["simulate", abspath(scenario_path), "--out", out]
    dt_ms === nothing || push!(args, "--dt-ms", string(dt_ms))
    pf_command(args...)
    return read_pf_results(out)
end

"""
    pf_inspect(; outdir = outputs/pf/model) -> Dict

Dump the PowerFactory model (network, machines, controllers, load flow) and read it.
"""
function pf_inspect(; outdir = nothing)
    out = outdir === nothing ? joinpath(PORTHOS_ROOT, "outputs", "pf", "model") : String(outdir)
    pf_command("inspect", "--out", out)
    return _plain(read_json(joinpath(out, "model.json")))
end

"""
    pf_machine_map(case) -> Vector{Pair{String,Vector{String}}}

PowerFactory machine name => the Porthos machines that represent it, from the case's `_pf`
keys (`"G 05 (unit 1/2)"` is a unit of `G 05`), sorted by PowerFactory name.
"""
function pf_machine_map(case::Case)
    m = Dict{String,Vector{String}}()
    for s in case.components
        s.type in ("GENROU_PHTRUE", "GENSAL_PHTRUE") || continue
        tag = get(case.raw[:components][Symbol(s.name)], :_pf, nothing)
        tag === nothing && continue
        push!(get!(m, String(split(String(tag), " (")[1]), String[]), s.name)
    end
    return sort!(collect(m); by = first)
end

# PowerFactory series on increasing times: of two rows at the same time (an event), the later
function _pf_series(t, y)
    keep = [i == length(t) || t[i + 1] > t[i] for i in eachindex(t)]
    return t[keep], y[keep]
end

function _interp(tq, t, y)
    out = similar(tq, Float64)
    j = 1
    for (i, x) in enumerate(tq)
        (x < t[1] || x > t[end]) && (out[i] = NaN; continue)
        while j < length(t) - 1 && t[j + 1] < x
            j += 1
        end
        w = t[j + 1] == t[j] ? 0.0 : (x - t[j]) / (t[j + 1] - t[j])
        out[i] = (1 - w) * y[j] + w * y[j + 1]
    end
    return out
end

function _unwrap_deg(a)
    out = copy(a)
    for i in 2:length(a)
        d = a[i] - a[i - 1]
        out[i] = out[i - 1] + d - 360.0 * round(d / 360.0)
    end
    return out
end

_rms(x) = sqrt(sum(abs2, x) / length(x))

"""
    pf_compare(sim_t, sim, pf, case; t_fault, t_clear, ref = "G 01") -> Dict

PHPS's PowerFactory comparison metrics (`study/src/model_review_pf_compare.py`), per
PowerFactory machine, on the Porthos times `sim_t` inside PowerFactory's run:

- rotor angle relative to `ref`, offset removed over `t < t_fault - 0.05`, in degrees
  (PowerFactory `s:firot`, unwrapped; Porthos `delta` of the first unit);
- speed (`s:speed` against `omega`), in pu;
- terminal active power (`m:Psum:bus1` / MVA base against the sum of the units' `Pe`), in pu;

each as RMS and maximum over the run, and as RMS over three windows: fault to clearing +
0.5 s, then to 2.5 s, then to the end. `sim(name)` returns a Porthos signal
(`"GENROU_4.delta"`) on `sim_t`. PowerFactory is interpolated linearly onto `sim_t`.
"""
function pf_compare(sim_t::AbstractVector, sim, pf::PFResults, case::Case; t_fault::Real,
                    t_clear::Real, ref::AbstractString = "G 01")
    mva = case.config.mva_base
    mmap = pf_machine_map(case)
    units = Dict(mmap)
    haskey(units, ref) || throw(ArgumentError("no Porthos machine for the reference $ref"))
    tq = [t for t in sim_t if pf.t[1] <= t <= pf.t[end]]
    keep = [pf.t[1] <= t <= pf.t[end] for t in sim_t]
    pre = tq .< t_fault - 0.05
    pfs(obj, var) = _interp(tq, _pf_series(pf.t, pf_signal(pf, obj, var))...)
    rel(a, r) = (d = a .- r; d .- sum(d[pre]) / count(pre))
    ang_pf_ref = _unwrap_deg(pfs(ref, "s:firot"))
    ang_ph_ref = rad2deg.(sim(units[ref][1] * ".delta")[keep])
    wins = ["fault_to_clear+0.5s" => (tq .>= t_fault) .& (tq .< t_clear + 0.5),
            "clear+0.5s_to_2.5s" => (tq .>= t_clear + 0.5) .& (tq .< 2.5),
            "2.5s_to_end" => tq .>= 2.5]
    out = Dict{String,Any}()
    for (g, us) in mmap
        ea = rel(rad2deg.(sim(us[1] * ".delta")[keep]), ang_ph_ref) .-
             rel(_unwrap_deg(pfs(g, "s:firot")), ang_pf_ref)
        ew = sim(us[1] * ".omega")[keep] .- pfs(g, "s:speed")
        ep = sum(sim(u * ".Pe")[keep] for u in us) .- pfs(g, "m:Psum:bus1") ./ mva
        row = Dict{String,Any}("units" => us,
            "angle_rms_deg" => _rms(ea), "angle_max_deg" => maximum(abs, ea),
            "speed_rms_pu" => _rms(ew), "speed_max_pu" => maximum(abs, ew),
            "p_rms_pu" => _rms(ep), "p_max_pu" => maximum(abs, ep))
        row["windows"] = Dict(w => Dict("angle_rms_deg" => _rms(ea[m]),
                                        "speed_rms_pu" => _rms(ew[m]), "p_rms_pu" => _rms(ep[m]))
                              for (w, m) in wins if any(m))
        out[g] = row
    end
    return out
end
