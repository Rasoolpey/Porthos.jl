# The PowerFactory reader and comparison helpers, on a synthetic result (no PowerFactory
# needed).

@testset "PowerFactory interface" begin
    dir = mktempdir()
    cols = [Dict("column" => 0, "object" => "G 01", "class" => "ElmSym", "variable" => "s:firot",
                 "unit" => "deg", "description" => "Rotor angle"),
            Dict("column" => 1, "object" => "Bus 16", "class" => "ElmTerm", "variable" => "m:u",
                 "unit" => "p.u.", "description" => "Voltage, Magnitude")]
    Porthos.write_json(joinpath(dir, "run.json"),
                       Dict("csv" => "pf_results.csv", "columns" => cols, "project" => "x"))
    open(joinpath(dir, "pf_results.csv"), "w") do io
        println(io, "Porthos results,G 01,Bus 16")
        println(io, "\"Time in s\",\"firot in deg\",\"u, Magnitude in p.u.\"")
        for (t, a, u) in ((0.0, 170.0, 1.0), (1.0, 179.0, 1.0), (1.0, 179.0, 0.1), (2.0, -178.0, 0.9))
            println(io, join((t, a, u), ","))
        end
    end
    r = read_pf_results(dir)
    @test r.t == [0.0, 1.0, 1.0, 2.0]
    @test pf_signal(r, "Bus 16", "m:u") == [1.0, 1.0, 0.1, 0.9]
    @test_throws KeyError pf_signal(r, "Bus 16", "m:phiu")
    # of two rows at an event, the later one is kept
    t, u = Porthos._pf_series(r.t, pf_signal(r, "Bus 16", "m:u"))
    @test t == [0.0, 1.0, 2.0] && u == [1.0, 0.1, 0.9]
    @test Porthos._interp([0.5, 1.5, 3.0], t, u) ≈ [0.55, 0.5, NaN] nans = true
    # PowerFactory's rotor angle is wrapped to +-180 deg
    @test Porthos._unwrap_deg(pf_signal(r, "G 01", "s:firot")) ≈ [170.0, 179.0, 179.0, 182.0]

    # a mismatch between the CSV and the column map is refused
    Porthos.write_json(joinpath(dir, "run.json"),
                       Dict("csv" => "pf_results.csv", "columns" => reverse(cols)))
    @test_throws ErrorException read_pf_results(dir)

    # machine map from the case's _pf keys
    m = Dict(pf_machine_map(load_case(joinpath(ROOT, "cases", "IEEE39Bus_PF", "system_phtrue.json"))))
    @test sort(collect(keys(m))) == ["G " * lpad(k, 2, '0') for k in 1:10]
    @test m["G 05"] == ["GENROU_5", "GENROU_6"]
    @test m["G 10"] == ["GENROU_11"]
end
