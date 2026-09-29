# P3 gate: bus voltages match PHPS and the case v0/a0 within 1e-8.

@testset "P3 power-flow parity" begin
    for c in Porthos.pack_cases(PACK)
        @testset "$(c.name)" begin
            ref = Porthos.pack_json(PACK, "powerflow/$(c.name).json")
            case = load_case(case_path(c.system))
            pf = solve_powerflow(case; tol = Float64(ref[:tol]), max_iter = Int(ref[:max_iter]))
            @test pf.spec.net.bus_ids == Int.(ref[:bus_ids])
            @test Int.(pf.spec.types) == Int.(ref[:bus_types])
            @test pf.spec.P_spec == Porthos.pack_vector(ref[:P_spec])
            @test pf.spec.Q_spec == Porthos.pack_vector(ref[:Q_spec])
            @test pf.spec.V0 == Porthos.pack_vector(ref[:V_start])
            @test pf.spec.theta0 == Porthos.pack_vector(ref[:theta_start])
            @test pf.from_case == ref[:skip_pf_solve]
            @test pf.converged == ref[:converged]
            ref[:iterations] === nothing || @test pf.iterations == ref[:iterations]

            # Gate, part 1: PHPS voltages.
            @test close_to(pf.V, Porthos.pack_vector(ref[:V]); atol = 1e-8)
            @test close_to(pf.theta, Porthos.pack_vector(ref[:theta]); atol = 1e-8)

            # Gate, part 2: the case v0/a0. PHPS itself misses this: the cases store
            # PowerFactory voltages rounded to 6 decimals (|V - v0| up to 4.6e-7), and in the
            # converter cases a0 keeps PowerFactory's angle reference (bus 31) while the slack
            # is bus 39, a uniform 0.1755 rad shift. Recorded as broken, not loosened; see
            # TODO.md (open question on the P3 gate).
            v0 = Porthos.pack_vector(ref[:case_v0])
            a0 = Porthos.pack_vector(ref[:case_a0])
            @test_broken close_to(pf.V, v0; atol = 1e-8) && close_to(pf.theta, a0; atol = 1e-8)
            # Supplementary, not the gate: agreement with v0/a0 at the precision they are
            # stored with, angles compared relative to the slack bus (reference-invariant).
            s = findfirst(==(SLACK_BUS), pf.spec.types)
            @test close_to(pf.V, v0; atol = 1e-6)
            @test close_to(pf.theta .- pf.theta[s], a0 .- a0[s]; atol = 2e-6)
        end
    end
end
