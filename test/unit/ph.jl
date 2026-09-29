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

@testset "storages by Riccati and multi-port loop opening" begin
    # H = -G = 1/(s + 1) is passive: a storage with dV/dt <= u y + 0.1 u^2 exists
    m = PortModel("g", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [-1.0], 0.0)
    P, res = port_storage(m, 0.1)
    @test P[1, 1] > 0 && res <= 1e-10
    @test port_margin(m, 0.1).passes
    # -G = 2/(s + 1): Re(-G) + 0.05 > 0 at every frequency
    @test port_margin(PortModel("b", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [-2.0], 0.0), 0.05).passes
    # -G = -2/(s + 1): Re(-G) + 0.5 < 0 at low frequency, so no storage exists
    @test_throws ErrorException port_storage(PortModel("n", ["x"], "u", "y", fill(-1.0, 1, 1), [1.0], [2.0], 0.0), 0.5)

    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys, x, V = eq.sys, eq.x, eq.V
    proj = physical_projection(sys, x, V)
    ks = [findfirst(c -> Porthos.name(c) == "IEEEG1_$i", sys.comps) for i in 2:4]
    rest, ports, T = open_loops_model(sys, x, V, ks; input = "omega", output = "Tm", projection = proj)
    @test size(rest.B, 2) == 3 && size(rest.C, 1) == 3
    # closing the three loops again gives the full system; so does T' A T
    nr = size(rest.A, 1)
    blocks = [zeros(length(p.states), 0) for p in ports]
    Acl = rest.A
    for (j, p) in enumerate(ports)
        nk = length(p.states)
        top = hcat(Acl, zeros(size(Acl, 1), nk))
        top[1:nr, end - nk + 1:end] .= rest.B[:, j] * p.C'
        bottom = hcat(zeros(nk, size(Acl, 2)), p.A)
        bottom[:, 1:nr] .= p.B * rest.C[j, :]'
        Acl = vcat(top, bottom)
    end
    A = reduced_jacobian(sys, x, V)[proj.keep, proj.keep]
    lam = eigvals(proj.basis' * A * proj.basis)
    for M in (Acl, T' * A * T)
        l = eigvals(M)
        @test maximum(minimum(abs.(l .- v)) for v in lam) <= 1e-10 * maximum(abs, lam)
    end
end

@testset "terminal cut (route A)" begin
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys, x, V = eq.sys, eq.x, eq.V
    proj = physical_projection(sys, x, V)
    units, net = terminal_models(sys, x, V; projection = proj)
    @test length(units) == 11 && length(net.states) == 19
    # close the units and the network through KCL: the section's eigenvalues plus the common
    # rotation (0) of the synchronous frame
    nb = nbus(sys)
    pos = Dict(b => i for (i, b) in enumerate(sys.net.bus_ids))
    nu = sum(length(u.states) for u in units)
    n = nu + length(net.states)
    Ax, Bv, Cx, Kv = zeros(n, n), zeros(n, 2nb), zeros(2nb, n), copy(net.D)
    off = 0
    for u in units
        k = length(u.states)
        r = off + 1:off + k
        vb = 2pos[u.bus] - 1:2pos[u.bus]
        Ax[r, r] .= u.A
        Bv[r, vb] .= u.B
        Cx[vb, r] .+= u.C
        Kv[vb, vb] .-= u.D
        off += k
    end
    rl = nu + 1:n
    Ax[rl, rl] .= net.A
    Bv[rl, :] .= net.B
    Cx[:, rl] .-= net.C
    lam = eigvals(Ax + Bv * (Kv \ Cx))
    lsec = eigvals(proj.basis' * reduced_jacobian(sys, x, V)[proj.keep, proj.keep] * proj.basis)
    @test length(lam) == length(lsec) + 1
    @test maximum(minimum(abs.(lam .- l)) for l in lsec) <= 1e-10 * maximum(abs, lsec)
    @test minimum(abs, lam) <= 1e-10 * maximum(abs, lsec)
end
