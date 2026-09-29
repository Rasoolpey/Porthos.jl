# P7 gate: the bus-16 fault trajectories against PHPS's compiled runs (pack section `sim`),
# all started from PHPS's initialised state and parameters (pack section `dae`), so that the
# comparison is of the integrators alone.
#
# - BDF1 against PHPS BDF1 (dt = 5e-4, same Newton tolerance and switching times): every
#   state and bus voltage within 1e-9 at every 5 ms point of the run, and the same number
#   of steps where the Newton does not converge.
# - IDA against PHPS IDA, both at rtol = atol = 1e-10: within 1e-6 until the first limiter
#   enters a sliding mode, reported after it, over 15 s. The start of the sliding mode is
#   taken from PHPS's BDF1 record: the first run of at least P7_SLIDE_STEPS consecutive
#   non-converged steps on the same equation (a limiter crossing costs 1 to 3 steps; the
#   IEEEG3 pilot valve in the base case chatters for 239).
# - simulation_results.csv: the header equals PHPS's, and the complete rows the pack keeps
#   (states and observables) agree within 1e-9.
# - IDA at production tolerances: the difference to PHPS is reported, not gated (the step
#   sequences differ).

const P7_BDF1_TOL = 1e-9
const P7_IDA_TOL = 1e-6
const P7_CSV_TOL = 1e-9
const P7_SLIDE_STEPS = 10
const P7_GRID = 0.005

# Porthos's rows of the state vector, by name
p7_row_names(sys) = [state_names(sys);
                     ["V$(q)_Bus$b" for b in sys.net.bus_ids for q in ("d", "q")]]

"""
Largest difference between a Porthos run and a pack run at each of the pack's times:
(times, max abs difference per time, the row where it occurs).
"""
function p7_differences(r::SimResult, m, P)
    rows = p7_row_names(r.sys)
    idx = Dict(string(n) => j for (j, n) in enumerate(m[:columns]))
    sel = [(i, idx[n]) for (i, n) in enumerate(rows) if haskey(idx, n)]
    ours = Dict(round(Int, t / P7_GRID) => k for (k, t) in enumerate(r.t)
                if abs(t - P7_GRID * round(t / P7_GRID)) <= 1e-9)
    tp = Float64.(m[:times])
    err = zeros(length(tp))
    where_ = fill("", length(tp))
    for (k, t) in enumerate(tp)
        c = ours[round(Int, t / P7_GRID)]
        for (i, j) in sel
            e = abs(r.Y[i, c] - P[k, j])
            if !(e <= err[k])                  # NaN counts as the largest
                err[k] = e
                where_[k] = rows[i]
            end
        end
    end
    return tp, err, where_, length(sel)
end

"""Start of the first run of `n` or more consecutive non-converged BDF1 steps on one equation."""
function p7_sliding_start(nonconverged, dt; n = P7_SLIDE_STEPS)
    k = 1
    while k <= length(nonconverged)
        j = k
        while j < length(nonconverged) &&
              string(nonconverged[j + 1][:row]) == string(nonconverged[k][:row]) &&
              abs(Float64(nonconverged[j + 1][:t]) - Float64(nonconverged[j][:t]) - dt) <= 1e-9
            j += 1
        end
        j - k + 1 >= n && return Float64(nonconverged[k][:t])
        k = j + 1
    end
    return Inf
end

@testset "P7 simulation parity" begin
    @test Porthos.has_section(PACK, "sim")
    runs = Dict{String,Vector{String}}()
    for k in keys(PACK.manifest[:files])
        m = match(r"^sim/([^/]+)/([^/]+)\.json$", string(k))
        m === nothing || push!(get!(runs, m[1], String[]), m[2])
    end
    @test haskey(runs, "base")
    for cname in sort(collect(keys(runs)))
        @testset "$cname" begin
            @test issetequal(runs[cname], ["bdf1", "ida_prod", "ida_tight"])
            d = Porthos.pack_json(PACK, "dae/$cname.json")
            case = load_case(case_path(string(d[:system])))
            sc = load_scenario(case_path(string(d[:scenario])))
            sys = assemble(case, sc; init_params = phps_init_params(case, d))
            y0 = phps_initial_state(d)
            meta(run) = Porthos.pack_json(PACK, "sim/$cname/$run.json")
            data(run, m) = Porthos.pack_binary(PACK, "sim/$cname/$run.bin",
                                               Tuple(Int.(m[:shape])))

            @testset "BDF1" begin
                m = meta("bdf1")
                s = m[:settings]
                r = simulate_bdf1(sys, y0; dt = Float64(s[:dt]), duration = Float64(s[:duration]),
                                  log_dt = Float64(s[:log_dt]))
                @test !r.stopped_early
                tp, err, where_, nsel = p7_differences(r, m, data("bdf1", m))
                @test nsel == length(p7_row_names(sys)) - 1          # all but delta_COI
                @test tp[end] == Float64(s[:duration])
                @test maximum(err) <= P7_BDF1_TOL
                @test r.nonconverged == Int(m[:nonconverged_steps])
                @test r.first_nonconverged_t ≈ Float64(m[:nonconverged][1][:t])
                k = argmax(err)
                @info "P7 $cname BDF1: max |Porthos - PHPS| = $(err[k]) at t = $(tp[k]) " *
                      "($(where_[k])); non-converged steps $(r.nonconverged) (PHPS " *
                      "$(m[:nonconverged_steps]))"

                # simulation_results.csv: PHPS's header and its complete rows
                @test csv_columns(sys) == string.(m[:csv_header])
                logged = Dict(round(Int, t / 1e-3) => k for (k, t) in enumerate(r.t))
                worst = 0.0
                for row in m[:full_rows]
                    ref = Float64.(row)
                    ours = csv_row(sys, ref[1], view(r.Y, :, logged[round(Int, ref[1] / 1e-3)]))
                    @test close_to(ours, ref; atol = P7_CSV_TOL)
                    worst = max(worst, maximum(abs, ours .- ref))
                end
                @info "P7 $cname CSV: $(length(m[:full_rows])) complete rows, largest " *
                      "difference $worst"
            end

            @testset "IDA at rtol = atol = 1e-10" begin
                m = meta("ida_tight")
                s = m[:settings]
                r = simulate_ida(sys, y0; dt = Float64(s[:dt]), duration = Float64(s[:duration]),
                                 log_dt = Float64(s[:log_dt]), rtol = Float64(s[:rtol]),
                                 atol = Float64(s[:atol]))
                @test !r.stopped_early
                tp, err, where_, _ = p7_differences(r, m, data("ida_tight", m))
                @test tp[end] == Float64(s[:duration])
                t_slide = p7_sliding_start(meta("bdf1")[:nonconverged],
                                           Float64(meta("bdf1")[:settings][:dt]))
                before = tp .< t_slide
                @test any(before)
                @test maximum(err[before]) <= P7_IDA_TOL
                kb = findall(before)[argmax(err[before])]
                after = any(.!before) ? maximum(err[.!before]) : 0.0
                @info "P7 $cname IDA 1e-10: first sliding mode at t = $t_slide; before it " *
                      "max |Porthos - PHPS| = $(err[kb]) at t = $(tp[kb]) ($(where_[kb])); " *
                      "after it (not gated) $after"
            end

            @testset "IDA at production tolerances (reported)" begin
                m = meta("ida_prod")
                s = m[:settings]
                tol(x) = x === nothing ? nothing : Float64(x)
                r = simulate_ida(sys, y0; dt = Float64(s[:dt]), duration = Float64(s[:duration]),
                                 log_dt = Float64(s[:log_dt]), rtol = tol(s[:rtol]),
                                 atol = tol(s[:atol]))
                @test !r.stopped_early
                tp, err, where_, _ = p7_differences(r, m, data("ida_prod", m))
                k = argmax(err)
                @info "P7 $cname IDA production (rtol $(r.settings["rtol"]), atol " *
                      "$(r.settings["atol"])): max |Porthos - PHPS| = $(err[k]) at t = $(tp[k]) " *
                      "($(where_[k])), not gated"
            end
        end
    end
end
