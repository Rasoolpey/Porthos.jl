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

@testset "structure search (route B)" begin
    IA = Porthos.IntervalArithmetic
    # interval eigenvalue bound
    @test 0.99 < verified_min_eig(IA.interval.(Matrix(Diagonal([1.0, 2.0, 3.0])))) <= 1.0
    @test verified_min_eig(IA.interval.([1.0 2.0; 2.0 1.0])) < 0
    # A is Hurwitz but has no diagonal Lyapunov function ((As'P + PAs)_11 = 0.2 p_1 >= 0).
    # A hand-made dual certificate: Z >= 0 with the diagonal of AZ + ZA' positive.
    A = [0.1 1.0; -1.0 -1.0]
    diagonal = BitMatrix([true false; false true])
    Z = [1.0 -0.05; -0.05 0.01]
    c = pattern_certificate(A, diagonal, Z, zeros(2, 2))
    @test c.infeasible
    @test c.float_bound ≈ 0.08 / 1.01 rtol = 1e-12
    @test 0 < c.verified_bound <= c.float_bound
    # the full pattern: the Lyapunov-equation solution passes the rigorous check, I does not
    Q = lyap(Matrix(A'), Matrix(1.0I, 2, 2))
    @test verified_lyapunov(Q, A).lyapunov
    @test !verified_lyapunov(Matrix(1.0I, 2, 2), A).lyapunov
    @test !pattern_certificate(A, trues(2, 2), Z, zeros(2, 2)).infeasible
    # powers of two: the scaling is exact
    T = pow2_scaling([3.0, 1e-5, 4.4e5])
    @test all(t -> t == exp2(round(log2(t))), diag(T))
    B = [1.3 -2.7 0.1; 5.5 0.3 -1.1; 2.2 7.9 -3.3]
    @test T \ B * T == [B[i, j] * T[j, j] / T[i, i] for i in 1:3, j in 1:3]

    # the base case: patterns, the section map, and the dense storage
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    sys, x, V = eq.sys, eq.x, eq.V
    proj = physical_projection(sys, x, V)
    groups = state_groups(sys, proj)
    z, U, As, ref = reference_section(sys, x, V, proj)
    @test length(groups) == length(z) + 1 == 172
    @test groups[ref].angle && groups[ref].component == "GENROU_1"
    @test maximum(real, eigvals(As)) < 0
    for kind in (:component, :unit)
        m = storage_pattern(groups, kind)
        @test m == m' && all(m[i, i] for i in axes(m, 1))
        mz = section_pattern(m, U)
        # terms with the reference angle become couplings with every rotor angle
        ang = findall(g -> g.angle, groups[z])
        own = findall(g -> g.unit == "GENROU_1", groups[z])
        @test all(mz[i, j] for i in own, j in ang)
        @test mz == mz' && all(mz[i, j] for i in 1:length(z), j in 1:length(z) if m[z[i], z[j]])
    end
    @test count(storage_pattern(groups, :unit)) > count(storage_pattern(groups, :component))
    @test all(section_pattern(storage_pattern(groups, :full), U))
    m = storage_pattern(groups, :component)
    add_coupling!(m, groups, "GENROU_2", "GENROU_3")
    i2 = findall(g -> g.component == "GENROU_2", groups)
    i3 = findall(g -> g.component == "GENROU_3", groups)
    @test all(m[i2, i3]) && all(m[i3, i2])
    m = couple_states!(storage_pattern(groups, :unit), groups, g -> g.state == "omega", g -> g.angle)
    iw = findall(g -> g.state == "omega", groups)
    ia = findall(g -> g.angle, groups)
    @test all(m[iw, ia]) && all(m[ia, iw])
    # the dense Lyapunov solution, scaled by powers of two, passes the rigorous check
    n = length(z)
    Q0 = lyap(Matrix(As'), Matrix(1.0I, n, n))
    T = pow2_scaling(diag(Q0))
    Ass = T \ As * T
    @test lyapunov_check(Q0, As).lyapunov
    @test verified_lyapunov(T * Q0 * T, Ass).lyapunov
end

@testset "port-power residual audit" begin
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    smp = power_audit_samples(eq.sys, eq.x, eq.V; n_random = 3)
    a = port_power_audit(eq.sys, smp)
    id = a["identities_max_error"]
    @test id["network_balance"] < 1e-11 && id["kcl_residual_power"] < 1e-10
    @test id["machine_terminal"] < 1e-11 && id["machine_Pe_vs_terminal"] < 1e-11
    @test id["load_declared_vs_drawn"] < 1e-9           # PHPS's 12-digit load constants
    bt = a["by_type"]
    # the one-way reservoirs are lossless accounts: supply equals the storage rate exactly
    for t in ("IEEEG1_PHTRUE", "IEEEG3_PHTRUE", "IEEET1_PHTRUE")
        @test bt[t]["complete"] && abs(bt[t]["residual_min"]) < 1e-9 && abs(bt[t]["residual_max"]) < 1e-9
    end
    # machines: the swing equation and the terminal balance close exactly; the residual
    # against the physical kinetic storage is all magnetic
    for t in ("GENROU_PHTRUE", "GENSAL_PHTRUE")
        @test max(abs(bt[t]["mechanical_residual_min"]), abs(bt[t]["mechanical_residual_max"])) < 1e-10
        @test max(abs(bt[t]["terminal_residual_min"]), abs(bt[t]["terminal_residual_max"])) < 1e-10
        @test bt[t]["residual_physical_kinetic_min"] ≈ bt[t]["magnetic_residual_min"] atol = 1e-10
    end
    # at the equilibrium the storage is at rest and every machine dissipates (field losses)
    at_eq = a["per_sample"][1]["components"]
    @test all(d -> abs(d["storage_rate"]) < 1e-9, at_eq)
    @test all(d -> d["magnetic_residual"] > 0, filter(d -> haskey(d, "magnetic_residual"), at_eq))
    @test !bt["COMPLEXLOAD"]["complete"]                 # its contract port is not an expression
    # GENROU's rotor as a linear port system: the form matches the model, the field port is
    # colocated, the balance closes with the stator-exchange port, and the rotor loss matrix
    # is proved positive definite (the magnetic residual is exchange, not creation)
    @test id["rotor_model"] < 1e-12 && id["field_colocation"] < 1e-12 && id["rotor_identity"] < 1e-12
    g = bt["GENROU_PHTRUE"]
    @test g["rotor_loss_matrix_min_eig_lower_bound"] > 0 && g["rotor_loss_min"] > 0
    rows = [d for smp in a["per_sample"] for d in smp["components"] if haskey(d, "rotor_loss")]
    @test length(rows) == 10 * length(smp)
    @test all(d -> abs(d["magnetic_residual"] - (d["rotor_loss"] - d["stator_exchange"])) < 1e-12, rows)
    @test rotor_structure(eq.sys.comps[findfirst(c -> model_type(c) == "GENSAL_PHTRUE", eq.sys.comps)]) === nothing
    # a fault-on KCL solve: the voltages differ from the healthy ones and satisfy the fault KCL
    Vf = solve_network(eq.sys, eq.x, eq.V; faults_on = true)
    @test maximum(abs, dae_residual(eq.sys, eq.x, Vf; faults_on = true)[2]) < 1e-10
    @test maximum(abs, Vf .- eq.V) > 1e-2
end

@testset "stator exchange through the KCL branch" begin
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    # current-corrected rotor coordinates: no damper current and field power = loss at rest
    rc = rotor_current_coordinates(eq.sys, eq.x, eq.V)
    @test length(rc) == 10
    @test all(r -> maximum(abs, r.Qw[2:4]) < 1e-12 && abs(r.field_supply - r.loss) < 1e-12, values(rc))
    @test all(r -> r.i_fd_runtime != r.i_fd_energy, values(rc))     # the runtime i_fd is unchanged
    # the KCL branch derivatives through kcl_solve match finite differences
    m = section_model(eq; scale = :none)
    n = Porthos.neta(m)
    u = normalize(randn(Random.Xoshiro(8), n))
    h = 1e-6
    fd = (exchange_one_form(m, h .* u) .- exchange_one_form(m, -h .* u)) ./ (2h)
    ad = ForwardDiff.derivative(t -> exchange_one_form(m, t .* u), 0.0)
    @test norm(ad - fd) <= 1e-6 * norm(ad)
    # the one-form is not exact, on the lossless variant too, and its curl has the closed form
    sysL = lossless_variant(eq.sys)
    @test iszero(sysL.G) && iszero(sysL.load.G)
    VL = solve_network(sysL, eq.x, eq.V)
    for (sys, V0) in ((eq.sys, eq.V), (sysL, VL))
        r = one_form_exactness(m, zeros(n); sys, V0, pairs = 1, seed = 5)
        @test r.curl > 1e-2
        W = exchange_curl(m, zeros(n); sys, V0)
        rng = Random.Xoshiro(5)
        uu, vv = normalize(randn(rng, n)), normalize(randn(rng, n))
        @test only(r.pairs).vJu - only(r.pairs).uJv ≈ dot(vv, W * uu) rtol = 1e-8
    end
end
