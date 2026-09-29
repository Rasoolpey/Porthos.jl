import JSON3

@testset "case loader" begin
    case = load_case(case_path("IEEE39Bus_PF/system_phtrue.json"))
    @test length(case.buses) == 39
    @test length(case.lines) == 46
    @test length(case.pq) == 19 && length(case.pv) == 9 && length(case.slack) == 1
    @test case.config.mva_base == 100.0
    @test case.config.omega_b === 2.0 * π * 60.0
    @test case.buses[1] == BusData(1, "BUS1", 345.0, 1.047356, -0.1472830995880455)
    @test case.lines[1].idx == "Line_1" && case.lines[1].x == 0.0411
    @test component(case, "CLOAD_3").type == "COMPLEXLOAD"
    @test param(component(case, "CLOAD_3"), :P0) == 3.22
    @test Wire("BUS_3.Vd", "CLOAD_3.Vd") in case.connections
    @test_throws KeyError component(case, "nope")
end

@testset "schema validation" begin
    good = read(case_path("IEEE39Bus_PF/system_phtrue.json"), String)
    obj = JSON3.read(good, Dict{String,Any})
    # missing required table
    bad = copy(obj); delete!(bad, "Line")
    @test_throws Porthos.SchemaError validate_json(:system, bad)
    # bus idx must be an integer
    bad = deepcopy(obj); bad["Bus"][1]["idx"] = "1"
    @test_throws Porthos.SchemaError validate_json(:system, bad)
    # wire target must be COMP.port
    bad = deepcopy(obj); bad["connections"][1]["to"] = "CLOAD_3"
    @test_throws Porthos.SchemaError validate_json(:system, bad)
    @test validate_json(:system, obj) === nothing

    sc = JSON3.read(read(case_path("IEEE39Bus_PF/bus_fault_bus16_150ms.json"), String),
                    Dict{String,Any})
    @test validate_json(:scenario, sc) === nothing
    bad = deepcopy(sc); bad["solver"]["method"] = "euler"
    @test_throws Porthos.SchemaError validate_json(:scenario, bad)
    bad = deepcopy(sc); delete!(bad["events"][1], "bus")
    @test_throws Porthos.SchemaError validate_json(:scenario, bad)
end

@testset "scenario loader" begin
    sc = load_scenario(case_path("IEEE39Bus_PF/bus_fault_bus16_150ms.json"))
    @test sc.solver.method === :ida
    @test sc.solver.dt == 5e-4 && sc.solver.duration == 6.0 && sc.solver.log_dt == 1e-3
    @test sc.solver.rtol === nothing
    @test sc.events == [BusFault(16, 1e-4, 1e-4, 1.0, 1.0 + 0.15)]
    @test sc.system_path == normpath(case_path("IEEE39Bus_PF/system_phtrue.json"))
    @test load_json_input(sc.system_path) isa Case
    @test load_json_input(sc.path) isa Scenario
end

@testset "event defaults (PHPS DAE)" begin
    ev(d) = Porthos._event(JSON3.read(JSON3.write(d)))
    @test ev(Dict("type" => "BusFault", "bus" => 3, "t_start" => 1.0)) ==
          BusFault(3, 0.0, 1e-5, 1.0, 1.1)
    @test ev(Dict("type" => "BusFault", "bus" => 3, "t_start" => 1.0, "t_end" => 1.2)) ==
          BusFault(3, 0.0, 1e-5, 1.0, 1.2)
    lf = ev(Dict("type" => "LineFault", "line_idx" => "Line_1", "t_start" => 1.0,
                 "t_duration" => 0.1))
    @test lf == LineFault("Line_1", 0.5, 0.0, 0.001, 1.0, 1.1)
    @test ev(Dict("type" => "Toggler", "model" => "Line")) isa OtherEvent
end

@testset "contracts" begin
    cs = load_contracts()
    @test cs.schema_version == "2.1-reservoir-corrected"
    @test cs.power_convention == "positive_into_component"
    for t in ("GENROU_PHTRUE", "GENSAL_PHTRUE", "IEEET1_PHTRUE", "IEEEG1_PHTRUE",
              "IEEEG3_PHTRUE", "COMPLEXLOAD", "GFL_PHTRUE", "GFM_VSM_PHTRUE",
              "GFM_DROOP_PHTRUE", "GFM_VOC_PHTRUE")
        @test haskey(cs.entries, contract_key(t))
    end
    g = contract(cs, "GENROU_PHTRUE")
    @test g.domain_clauses[1].id == "omega_positive"
    @test g.domain_clauses[1].check === :state_bounds
    @test g.domain_clauses[1].state == "omega" && g.domain_clauses[1].lower == "0"
    @test :parameter_conditions in [c.check for c in g.domain_clauses]
    load = contract(cs, "COMPLEXLOAD")
    @test load.model == "ComplexLoad"
    @test only(load.domain_clauses).threshold == 1e-8
    @test sum(length(e.domain_clauses) for e in values(cs.entries)) == 31
end

@testset "machine base normalisation" begin
    case = load_case(case_path("IEEE39Bus_PF/system_phtrue.json"))
    spec = component(case, "GENROU_2")                      # Sn = 700 MVA
    p = component_params(case, spec)
    Z = 100.0 / 700.0
    @test p["xd_double_prime"] == param(spec, :xd_double_prime) * Z
    @test p["H"] == param(spec, :H) * (700.0 / 100.0)
    @test p["omega_b"] == "2.0 * M_PI * 60.0"
    @test p["_params_normalized"] === true
    # already on the system base: not scaled again
    q = Porthos.ParamDict("Sn" => 700.0, "xd1" => 0.3, "xd_prime" => 0.3 * Z, "xd" => 1.8 * Z,
                          "xl" => 0.15 * Z)
    r = Porthos.normalise_machine_params(q, 100.0, 60.0)
    @test r["xd_prime"] == 0.3 * Z && r["xd"] == 1.8 * Z
    # no D declared: 2 pu on the machine base
    @test Porthos.normalise_machine_params(Porthos.ParamDict("Sn" => 700.0), 100.0, 60.0)["D"] ==
          2.0 * 7.0
    # an explicit D = 0 is honoured
    @test Porthos.normalise_machine_params(Porthos.ParamDict("Sn" => 700.0, "D" => 0.0), 100.0,
                                           60.0)["D"] == 0.0
end
