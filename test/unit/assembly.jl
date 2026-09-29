using ForwardDiff

@testset "wiring semantics" begin
    case = toy_case(extra = Dict(
        "components" => Dict(
            "G1" => Dict("type" => "GENROU_PHTRUE",
                         "params" => Dict("bus" => 1, "Sn" => 100.0, "ra" => 0.0, "xd" => 1.8,
                                          "xq" => 1.7, "xd_prime" => 0.3, "xq_prime" => 0.55,
                                          "xd_double_prime" => 0.25, "xq_double_prime" => 0.25,
                                          "xl" => 0.15, "Td0_prime" => 8.0, "Tq0_prime" => 0.4,
                                          "Td0_double_prime" => 0.03,
                                          "Tq0_double_prime" => 0.05, "H" => 3.0, "D" => 0.0,
                                          "omega_b" => "2.0 * M_PI * 50.0")),
            "AVR" => Dict("type" => "IEEET1_PHTRUE",
                          "params" => Dict("TR" => 0.0, "KA" => 50.0, "TA" => 0.05, "KE" => 1.0,
                                           "TE" => 0.5, "KF" => 0.05, "TF" => 1.0, "VRMAX" => 5.0,
                                           "VRMIN" => -5.0, "SAT_A" => 0.0, "SAT_B" => 0.0))),
        "connections" => [
            Dict("from" => "BUS_1.Vd", "to" => "G1.Vd"),
            Dict("from" => "BUS_1.Vq", "to" => "G1.Vq"),
            Dict("from" => "AVR.Efd", "to" => "G1.Efd"),
            Dict("from" => "BUS_1.Vterm", "to" => "AVR.Vterm"),
            Dict("from" => "CONST:1.02", "to" => "AVR.Vref"),
            Dict("from" => "CONST:1.05", "to" => "AVR.Vref"),        # later wire wins
            Dict("from" => "G1.nope", "to" => "AVR.upss"),           # unknown port -> 0
            Dict("from" => "BUS_1.Freq", "to" => "AVR.i_fd"),        # unknown signal -> 0
            Dict("from" => "PARAM:G1.H", "to" => "G1.Tm"),
        ]))
    net = Network(case)
    comps = [build_component(case, s) for s in case.components]
    src = resolve_wiring(case, comps, net)
    k = Dict(Porthos.name(c) => i for (i, c) in enumerate(comps))
    g, a = src[k["G1"]], src[k["AVR"]]
    @test g[1] == Porthos.InputSource(Porthos.SRC_VD, 0.0, 1, 0)
    @test g[3] == Porthos.const_source(3.0)                  # PARAM:G1.H (system base)
    @test g[4] == Porthos.InputSource(Porthos.SRC_OUTPUT, 0.0, k["AVR"], 1)
    @test a[1] == Porthos.InputSource(Porthos.SRC_VTERM, 0.0, 1, 0)
    @test a[2] == Porthos.const_source(1.05)
    @test a[3] == Porthos.const_source(0.0)
    @test a[4] == Porthos.const_source(0.0)
    # refresh after initialisation: Vref and an undriven Tm take the initialised values
    comps2 = [with_params(comps[k["G1"]], Dict("Tm0" => 0.7, "Efd0" => 1.9)),
              with_params(comps[k["AVR"]], Dict("Vref" => 1.01))]
    comps2 = k["G1"] == 1 ? comps2 : reverse(comps2)
    src2 = resolve_wiring(case, comps2, net)
    @test src2[k["G1"]][3] == Porthos.const_source(0.7)     # Tm: constant, refreshed
    @test src2[k["G1"]][4].kind === Porthos.SRC_OUTPUT       # Efd: driven, kept
    @test src2[k["AVR"]][2] == Porthos.const_source(1.01)
end

@testset "DAE residual and sparsity" begin
    case = load_case(case_path("IEEE39Bus_PF/system_phtrue.json"))
    sc = load_scenario(case_path("IEEE39Bus_PF/bus_fault_bus16_150ms.json"))
    sys = assemble(case, sc)
    @test length(state_names(sys)) == sys.n_diff == 203
    @test state_names(sys)[end] == "delta_COI"
    x = zeros(sys.n_diff)
    V = zeros(nalg(sys))
    for (i, b) in enumerate(sys.net.bus_ids)
        V[2i - 1] = sys.pf.V[i] * cos(sys.pf.theta[i])
        V[2i] = sys.pf.V[i] * sin(sys.pf.theta[i])
    end
    # a plausible state: machines near synchronism, controllers and reservoirs mid-range
    for (k, c) in enumerate(sys.comps)
        o = sys.offsets[k]
        n = nstates(c)
        x[o:o + n - 1] .= 0.5
        component_role(c) === :generator && (x[o + 1] = 1.0)
        "x_field" in state_names(c) && (x[o + 4] = 1000.0)
        "x_steam" in state_names(c) && (x[o + 6] = 1000.0)
        "x_water" in state_names(c) && (x[o + 4] = 1000.0)
    end
    f, g = dae_residual(sys, x, V)
    @test all(isfinite, f) && all(isfinite, g)
    # generic in the number type, and the structural pattern covers the Jacobian
    P = jacobian_pattern(sys)
    y = vcat(x, V)
    nd = sys.n_diff
    F(y) = vcat(dae_residual(sys, y[1:nd], y[nd + 1:end]; faults_on = true)...)
    J = ForwardDiff.jacobian(F, y)
    @test all(P[i, j] for (i, j) in zip(findnz(sparse(J))[1:2]...))
    @test nnz(P) < 0.1 * length(P)                           # it is sparse
end
