import JSON3

# A small synthetic case built in memory: 3 buses, a line, a phase-shifting transformer,
# a shunt, a load and one machine.
function toy_case(; extra = Dict{String,Any}())
    d = Dict{String,Any}(
        "config" => Dict("mva_base" => 100.0, "fn" => 50.0),
        "Bus" => [Dict("idx" => 3, "v0" => 1.0), Dict("idx" => 1, "v0" => 1.02),
                  Dict("idx" => 2, "v0" => 0.98)],
        "PQ" => [Dict("idx" => "L2", "bus" => 2, "p0" => 0.5, "q0" => 0.2)],
        "PV" => [Dict("idx" => 3, "bus" => 3, "p0" => 0.3, "v0" => 1.01)],
        "Slack" => [Dict("idx" => 1, "bus" => 1, "v0" => 1.02, "a0" => 0.0)],
        "Line" => [Dict("idx" => "L12", "bus1" => 1, "bus2" => 2, "r" => 0.01, "x" => 0.1,
                        "b" => 0.02),
                   Dict("idx" => "T23", "bus1" => 2, "bus2" => 3, "r" => 0.0, "x" => 0.05,
                        "tap" => 1.05, "phi" => 0.1)],
        "Shunt" => [Dict("idx" => "S3", "bus" => 3, "g" => 0.0, "b" => 0.1)],
        "components" => Dict("G1" => Dict("type" => "GENROU_PHTRUE",
                                          "params" => Dict("bus" => 1, "Sn" => 200.0,
                                                           "ra" => 0.002,
                                                           "xd_double_prime" => 0.2,
                                                           "xl" => 0.1))),
        "connections" => Any[],
    )
    merge!(d, extra)
    raw = JSON3.read(JSON3.write(d))
    Porthos.validate_json(:system, raw)
    return Porthos.Case("<toy>", raw)
end

@testset "Y-bus stamps" begin
    case = toy_case()
    net = Network(case)
    @test net.bus_ids == [1, 2, 3]                       # sorted, like PHPS
    Y = Matrix(ybus_pf(case))
    y12 = 1 / complex(0.01, 0.1)
    y23 = 1 / complex(0.0, 0.05)
    a = 1.05 * cis(0.1)
    @test Y[1, 1] ≈ y12 + 0.01im
    @test Y[1, 2] ≈ -y12
    @test Y[2, 2] ≈ y12 + 0.01im + y23 / abs2(a)
    @test Y[2, 3] ≈ -y23 / conj(a)
    @test Y[3, 2] ≈ -y23 / a
    @test Y[3, 3] ≈ y23 + 0.1im
    @test Y[1, 3] == 0
    # with loads: (P - jQ)/v0^2 at bus 2 (v0 = 0.98)
    YL = Matrix(ybus(case; loads = true))
    @test YL[2, 2] - Y[2, 2] ≈ complex(0.5, -0.2) / 0.98^2
    # DAE: plus the machine Norton admittance on the system base (Sn = 200 -> x 0.5)
    st = only(norton_stamps(case))
    @test st.bus == 1 && st.ra ≈ 0.001 && st.xd_pp ≈ 0.1
    YD = Matrix(ybus_dae(case))
    @test YD[1, 1] - YL[1, 1] ≈ 1 / complex(0.001, 0.1)
    # no Norton for grid-following converters
    gfl = toy_case(extra = Dict("components" => Dict(
        "C1" => Dict("type" => "GFL_PHTRUE", "params" => Dict("bus" => 3)))))
    @test isempty(norton_stamps(gfl))
end

@testset "CPython complex division" begin
    # Agrees with the mathematically correct value to rounding, whichever branch runs.
    for z in (complex(0.01, 0.1), complex(0.3, 0.02), complex(0.0, 1e-4), complex(2.0, 2.0))
        @test Porthos._py_cdiv(1.0, z) ≈ 1 / z rtol = 4eps()
    end
    @test_throws DivideError Porthos._py_cdiv(1.0, complex(0.0, 0.0))
end

@testset "faults" begin
    @test fault_admittance(1e-4, 1e-4) == (1e-4 / 2e-8, -1e-4 / 2e-8)
    @test fault_admittance(0.0, 0.0) == (0.0, -0.0)           # z^2 clipped at 1e-20
    case = toy_case()
    net = Network(case)
    sh = fault_shunts(net, [BusFault(2, 0.0, 1e-5, 1.0, 1.1), BusFault(9, 0.0, 1e-5, 1.0, 1.1)])
    @test length(sh) == 1 && sh[1].index == 2
    Y = ybus_dae(case)
    Yf = with_fault(Y, sh)
    @test Yf[2, 2] == Y[2, 2] + complex(fault_admittance(0.0, 1e-5)...)
    @test Yf[1, 1] == Y[1, 1]

    lines, fb, bf = split_line_for_fault(case, LineFault("L12", 0.25, 0.0, 0.001, 1.0, 1.1))
    @test fb == 99 && bf == BusFault(99, 0.0, 0.001, 1.0, 1.1)
    @test [l.idx for l in lines] == ["T23", "L12_part1", "L12_part2"]
    @test lines[2].x ≈ 0.025 && lines[3].x ≈ 0.075 && lines[2].bus2 == 99 && lines[3].bus1 == 99
    @test_throws ArgumentError split_line_for_fault(case, LineFault("nope", 0.5, 0.0, 0.0, 0, 1))
end

@testset "load admittances" begin
    case = toy_case()
    la = load_admittances(case)
    @test la.G ≈ [0.0, 0.5 / 0.98^2, 0.0] && la.B ≈ [0.0, -0.2 / 0.98^2, 0.0]
    @test !la.has_complex_loads && all(iszero, la.kpf)
end

@testset "only PHTRUE models" begin
    for t in ("GENROU_PHS", "GENROU", "GENCLS", "IEEEX1_PHS", "TGOV1_PHS")
        old = toy_case(extra = Dict("components" => Dict(
            "G1" => Dict("type" => t, "params" => Dict("bus" => 1)))))
        @test_throws UnsupportedModelError norton_stamps(old)
        @test_throws UnsupportedModelError component_params(old, only(old.components))
    end
    @test "GENROU_PHTRUE" in MODEL_TYPES && !("GENROU_PHS" in MODEL_TYPES)
end
