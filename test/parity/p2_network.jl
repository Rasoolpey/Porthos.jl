# P2 gate: Y-bus equal to PHPS within 1e-12 on all parity cases (power-flow and DAE
# matrices, fault shunts, load admittances, Norton stamps).

@testset "P2 network parity" begin
    for c in Porthos.pack_cases(PACK)
        @testset "$(c.name)" begin
            ref = Porthos.pack_json(PACK, "network/$(c.name).json")
            case = load_case(case_path(c.system))
            sc = load_scenario(case_path(c.scenario))
            net = Network(case)
            @test net.bus_ids == Int.(ref[:bus_ids])

            Ypf = Matrix(ybus_pf(case))
            Ypf_ref = Porthos.pack_matrix(ref[:Y_pf])
            @test close_to(Ypf, Ypf_ref; atol = 1e-12)
            Ydae = Matrix(ybus_dae(case))
            Ydae_ref = Porthos.pack_matrix(ref[:Y_dae])
            @test close_to(Ydae, Ydae_ref; atol = 1e-12)
            # Stronger than the gate: the arithmetic reproduces PHPS bit for bit.
            @test bitdiff(Ypf, Ypf_ref) == 0
            @test bitdiff(Ydae, Ydae_ref) == 0
            @test isempty(ref[:excluded_dyn_lines])

            st = norton_stamps(case)
            @test [s.component for s in st] == [string(n[:name]) for n in ref[:norton]]
            @test [s.bus for s in st] == [Int(n[:bus]) for n in ref[:norton]]
            @test [s.ra for s in st] == [Float64(n[:ra]) for n in ref[:norton]]
            @test [s.xd_pp for s in st] == [Float64(n[:xd_pp]) for n in ref[:norton]]

            fs = fault_shunts(net, sc.events)
            @test length(fs) == length(ref[:faults])
            for (f, r) in zip(fs, ref[:faults])
                @test f.bus == Int(r[:bus]) && f.index == Int(r[:bus_index])
                @test f.t_start == Float64(r[:t_start]) && f.t_end == Float64(r[:t_end])
                @test f.g == Float64(r[:g]) && f.b == Float64(r[:b])
            end

            la = load_admittances(case; net)
            for k in (:G, :B, :P, :Q, :kpf, :kqf)
                refk = Porthos.pack_vector(ref[Symbol("load_", k)])
                @test close_to(getfield(la, k), refk; atol = 1e-12)
                @test bitdiff(getfield(la, k), refk) == 0
            end
        end
    end
end
