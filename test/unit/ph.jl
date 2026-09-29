# Port passivity tools on transfer functions with known answers.

@testset "port passivity tools" begin
    # G(s) = 1/(s + 1): positive real
    lag = PortModel("lag", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [1.0], 0.0)
    @test transfer(lag, 0.0) ≈ 1.0
    @test transfer(lag, 1.0im) ≈ 0.5 - 0.5im
    @test isempty(real_part_crossings(lag; sign = 1))
    pc = passivity_certificate(lag; sign = 1)
    @test pc.H0 ≈ 1.0 && !pc.nonpassive && pc.shortage == 0.0
    @test isempty(port_zeros(lag))

    # G(s) = (1 - s)/(1 + s) = 2/(s + 1) - 1: zero at +1, Re G(jw) = (1 - w^2)/(1 + w^2)
    allpass = PortModel("allpass", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [2.0], -1.0)
    @test only(port_zeros(allpass)) ≈ 1.0
    @test only(real_part_crossings(allpass; sign = 1)) ≈ 1.0 rtol = 1e-12
    pc = passivity_certificate(allpass; sign = 1)
    @test pc.nonpassive
    @test pc.min_re ≈ -1.0 rtol = 1e-5
    @test pc.shortage ≈ 1.0 rtol = 1e-5
end

@testset "loop opening at a governor" begin
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys, x, V = eq.sys, eq.x, eq.V
    proj = physical_projection(sys, x, V)
    ins, _ = component_io(sys, x, V)
    k = findfirst(c -> Porthos.name(c) == "IEEEG1_4", sys.comps)
    r = sys.offsets[k]:(sys.offsets[k] + nstates(sys.comps[k]) - 1)
    mg = port_model(sys.comps[k], x[r], ins[k]; input = "omega", output = "Tm")
    mr = loop_port_model(sys, x, V, k; input = "omega", output = "Tm", projection = proj)
    @test size(mr.A, 1) == size(proj.basis, 2) - length(mg.states)
    # closing the governor and the rest again gives the full system's eigenvalues
    lam = eigvals(proj.basis' * reduced_jacobian(sys, x, V)[proj.keep, proj.keep] * proj.basis)
    lcl = eigvals([mr.A mr.B * mg.C'; mg.B * mr.C' mg.A])
    @test maximum(minimum(abs.(lcl .- l)) for l in lam) <= 1e-10 * maximum(abs, lam)
    # the frequency response equals direct solves
    ws = [1e-3, 0.7, 1.55, 40.0]
    @test frequency_response(mr, ws) ≈ [transfer(mr, im * w) for w in ws] rtol = 1e-10
    @test frequency_response(mg, ws) ≈ [transfer(mg, im * w) for w in ws] rtol = 1e-12
    # a component with more than one channel to the rest is refused
    kg = findfirst(c -> Porthos.name(c) == "GENROU_4", sys.comps)
    @test_throws ErrorException loop_port_model(sys, x, V, kg; input = "Tm", output = "Pe", projection = proj)
end
