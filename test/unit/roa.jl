# ROA certificate pipeline (src/roa/): the interval core on small known problems, then the
# whole pipeline on IEEE-39 at a small level, with sampled sanity checks of the enclosures.

import IntervalArithmetic as IA
using Random: Xoshiro

@testset "interval branch primitives" begin
    a, b = IA.interval(1.0, 2.0), IA.interval(3.0, 4.0)
    @test Porthos.lt(NoModes(), a, b) === true
    @test Porthos.gt(NoModes(), a, b) === false
    @test Porthos.ge(NoModes(), b, 3.0) === true
    @test_throws UndecidedBranch Porthos.gt(NoModes(), IA.interval(1.0, 3.5), b)
    # a clamp decided on the box, and one that straddles a limit
    @test IA.isequal_interval(Porthos.clamp_mode(NoModes(), a, 0.0, 5.0), a)
    @test_throws UndecidedBranch Porthos.clamp_mode(NoModes(), a, 0.0, 1.5)
    # duals over intervals compare by their value
    d = ForwardDiff.Dual(IA.interval(2.0, 2.5), IA.interval(1.0))
    @test Porthos.gt(NoModes(), d, 1.0) === true
    @test is_interval_type(typeof(d)) && !is_interval_type(ForwardDiff.Dual{Nothing,Float64,1})
    # the recorder keeps the decided outcomes
    log = Porthos.ModeLog()
    Porthos.guard_min(log, IA.interval(0.5, 0.7), 1e-6)
    @test log.decisions == [false]
end

@testset "interval expressions" begin
    env = Dict("PMAX" => 0.1, "R_base" => 3.0)
    v = eval_param_expr(IA.Interval{Float64}, "PMAX*R_base", env)
    @test IA.in_interval(0.1 * 3.0, v) && IA.diam(v) <= 1e-15
    @test parse_param_expr("PMAX*R_base", env) === 0.1 * 3.0
end

@testset "interval linear algebra" begin
    rng = Xoshiro(3)
    A0 = randn(rng, 6, 6) + 8I
    A = IA.interval.(A0 .- 1e-3, A0 .+ 1e-3)
    B = IA.interval.(randn(rng, 6, 2))
    X, beta = interval_solve(A, B)
    @test beta < 1
    for _ in 1:20
        As = A0 .+ 1e-3 .* (2 .* rand(rng, 6, 6) .- 1)
        @test all(IA.in_interval.(As \ IA.mid.(B), X))
    end
    @test_throws ProofFailure interval_solve(IA.interval.([1.0 1.0; 1.0 1.0] .+ [0 0; 0 1e-3] .* IA.interval(-1, 1)),
                                             IA.interval.(ones(2, 1)))
    # definiteness: sampled members stay below the rigorous bounds
    S0 = Symmetric(randn(rng, 8, 8))
    S0 = Matrix(S0) - (eigmax(S0) + 0.5) * I
    S = IA.interval.(S0 .- 0.02, S0 .+ 0.02)
    wb = weyl_max_eig(S)
    ub = verified_max_eig(S)
    for _ in 1:20
        D = 0.02 .* (2 .* rand(rng, 8, 8) .- 1)
        @test eigmax(Symmetric(S0 + (D + D') / 2)) <= min(wb.bound, ub)
    end
    @test wb.bound < 0
    @test interval_cholesky(.-S).positive_definite
    @test !interval_cholesky(IA.interval.([1.0 2.0; 2.0 1.0])).positive_definite
    # an ill-conditioned positive definite matrix needs the preconditioned Cholesky
    Q = qr(randn(rng, 40, 40)).Q
    P = Matrix(Q * Diagonal(exp10.(range(-6, 0; length = 40))) * Q')
    P = (P + P') / 2
    @test cholesky_positive_definite(IA.interval.(P)).positive_definite
    @test all(verified_inverse_diagonal(P) .>= diag(inv(P)) .* (1 - 1e-8))
    # Krawczyk on x^2 = 2
    F(x) = [x[1]^2 - 2]
    K, X, rec = krawczyk(F, x -> ForwardDiff.jacobian(F, x), "sqrt 2", [1.4], [0.1]; C = fill(1 / 2.8, 1, 1))
    @test IA.in_interval(sqrt(2), K[1]) && IA.diam(K[1]) < 1e-12 + 0.03
end

@testset "analytic certificates" begin
    certified(m, c, level; kw...) = certify_level(c, level, enclose_equilibrium(m); kw...)["certified"]
    methods = ((hull = :first_order,), (hull = :centered, second_order = :direction),
               (hull = :centered, second_order = :coordinates))
    # 1-D, x' = -x + v^3 with 0 = v - x: x' = -x + x^3, attracted for |x| < 1. P = 1/2, so
    # V = x^2/2 <= c is |x| <= w = sqrt(2c). All three hulls bound the mean-value matrix by
    # -1 + 3 w^2 (hull: 3v^2 over the box; centered: (1/2) w max|6x|), so they certify
    # exactly the levels c < 1/6; the true answer, c < 1/2, is never exceeded.
    m1 = AnalyticModel("cubic", (x, v) -> [-x[1] + v[1]^3], (x, v) -> [v[1] - x[1]], 1, [0.0])
    c1 = quadratic_candidate(m1)
    @test c1.P ≈ fill(0.5, 1, 1)
    for kw in methods
        @test certified(m1, c1, 0.16; kw...)
        @test !certified(m1, c1, 0.17; kw...)
        @test !certified(m1, c1, 0.5; kw...)
    end
    r1 = certify_roa(c1; level_low = 1e-3, level_high = 1.0, io = devnull)
    @test 0.16 < r1["verified_valid_level"] < 1 / 6
    @test occursin("region of attraction", r1["claim"])

    # 2-D, x' = -x (1 - v) with 0 = v - |x|^2: attracted for |x| < 1, V = |x|^2/2 < 1/2
    m2 = AnalyticModel("radial", (x, v) -> -x .* (1 - v[1]), (x, v) -> [v[1] - (x[1]^2 + x[2]^2)],
                       2, [0.0])
    c2 = quadratic_candidate(m2)
    @test c2.P ≈ Matrix(0.5I, 2, 2)
    for kw in methods
        @test certified(m2, c2, 0.05; kw...)
        @test !certified(m2, c2, 0.5; kw...)
    end
    r2 = certify_roa(c2; level_low = 1e-3, level_high = 1.0, io = devnull)
    @test 0.05 <= r2["verified_valid_level"] < 0.5

    # a limiter: certified inside one mode, not across the switching surface
    m3 = AnalyticModel("clamped", (x, v) -> [-x[1] + Porthos.clamp_mode(NoModes(), v[1]^3, -1e-3, 1e-3)],
                       (x, v) -> [v[1] - x[1]], 1, [0.0])
    c3 = quadratic_candidate(m3)
    e3 = enclose_equilibrium(m3)
    @test certify_level(c3, 1e-3, e3)["certified"]
    r = certify_level(c3, 0.1, e3)
    @test !r["certified"] && occursin("UndecidedBranch", r["failure"])

    # an unstable equilibrium has no quadratic candidate
    @test_throws ArgumentError quadratic_candidate(AnalyticModel("unstable", (x, v) -> [x[1] + v[1]],
                                                                 (x, v) -> [v[1]], 1, [0.0]))
end

@testset "rotation action" begin
    # a model type without a declared rotation action is rejected (VOC: Cartesian states)
    scv = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF_gfm-voc", "no_fault_gfm_voc.json"))
    eqv = solve_equilibrium(load_case(scv.system_path), scv)
    err = try
        section_model(eqv); nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("rotation action", err.msg)
end

@testset "ROA pipeline on IEEE-39" begin
    sc = load_scenario(joinpath(ROOT, "cases", "IEEE39Bus_PF", "bus_fault_bus16_150ms.json"))
    eq = solve_equilibrium(load_case(sc.system_path), sc)
    m = section_model(eq)
    n = Porthos.neta(m)
    @test n == 171 && all(s -> exp2(round(log2(s))) == s, m.scale)
    # coordinates: lift stays on the section and inverts section_coordinates
    rng = Xoshiro(5)
    eta = 1e-3 .* randn(rng, n)
    x = lift(m, eta)
    @test maximum(abs, section_coordinates(m, x) .- eta) <= 1e-13 * maximum(abs, m.x0)
    @test dot(m.l, x[m.keep]) ≈ dot(m.l, m.x0[m.keep]) rtol = 1e-14
    # the projected field keeps the section invariant: l' (lifted h) = 0
    h, _ = section_field(m, eta, eq.V)
    hk = zeros(length(m.keep))
    hk[m.z] .= m.scale .* h
    hk[m.ref] = -dot(m.l[m.z], hk[m.z]) / m.l[m.ref]
    f, _ = state_field(m, x, eq.V)
    @test abs(dot(m.l, hk)) <= 1e-12 * maximum(abs, f)

    c = quadratic_candidate(m)
    @test positivity_proof(c).proved
    e = enclose_equilibrium(m)
    @test maximum(IA.diam, e.eta) < 1e-9 && maximum(IA.diam, e.V) < 1e-12

    # enclosures at a small level, checked against sampled points
    level = 1e-12
    w = sublevel_half_widths(c, level)
    Xi = [IA.interval(-v, v) for v in w]
    E = e.eta .+ Xi
    br = enclose_kcl_branch(m, E, IA.mid.(e.V))
    M1, _ = jacobian_hull(m, E, br.V)
    Mc, crec = centered_hull(m, e, E, br.V, Xi; X = br.X)
    Md, drec = centered_hull(m, e, E, br.V, Xi; second_order = :direction, X = br.X)
    @test crec["centered_radius_max"] < drec["centered_radius_max"]
    @test drec["hull_width_max"] < drec["first_order_hull_width_max"]
    etastar = IA.mid.(e.eta)
    solveV(et) = (V = copy(IA.mid.(e.V)); for _ in 1:4
                      G(v) = section_field(m, et, v)[2]
                      V .-= ForwardDiff.jacobian(G, V) \ G(V)
                  end; V)
    gl = ((-0.906179845938664, 0.236926885056189), (-0.538469310105683, 0.478628670499366),
          (0.0, 0.568888888888889), (0.538469310105683, 0.478628670499366),
          (0.906179845938664, 0.236926885056189))
    for _ in 1:3
        xi = w .* (2 .* rand(rng, n) .- 1)
        et = etastar .+ xi
        V = solveV(et)
        @test all(IA.in_interval.(V, br.V))
        J = section_jacobian(m, et, V)
        @test all(IA.in_interval.(J, M1))
        # the mean-value matrix (5-point Gauss-Legendre) lies in the centered hull
        Mv = sum(wq / 2 .* section_jacobian(m, etastar .+ (1 + tq) / 2 .* xi, solveV(etastar .+ (1 + tq) / 2 .* xi))
                 for (tq, wq) in gl)
        @test all(IA.inf.(Mc) .- 1e-9 .<= Mv .<= IA.sup.(Mc) .+ 1e-9)
        @test all(IA.inf.(Md) .- 1e-9 .<= Mv .<= IA.sup.(Md) .+ 1e-9)
    end

    # one level: every gate on the same boxes
    row = certify_level(c, level, e)
    @test row["certified"] && row["decay"] && row["kcl_branch"] && row["containment"] && row["same_boxes"]
    recs = row["containment_records"]
    @test any(r -> r["id"] == "omega_positive" && r["component"] == "GENROU_1" && r["passed"], recs)
    @test any(r -> r["id"] == "reservoir_conjugate_positive" && r["component"] == "IEEEG1_2" && r["passed"], recs)
    @test any(r -> r["id"] == "voltage_above_correction_clip" && r["passed"], recs)
    @test all(r -> r["passed"], recs)
    @test certify_level(c, level, e; hull = :first_order)["certified"]
    @test row["excluded_ranges"]["IEEEG1_2.x_steam"]["drift_rate_upper"] < Inf
    @test all(r -> r["lower"] < r["x0"] < r["upper"], values(row["excluded_ranges"]))
    @test length(check_rotation_symmetry(eq.sys)) == 11
    # a large level leaves the smooth mode: not certified, with the reason
    big = certify_level(c, 1e-6, e)
    @test !big["certified"] && occursin("UndecidedBranch", big["failure"])

    # the record and the independent check
    rec = certify_roa(c; level_low = level, level_high = level, second_order = :direction, io = devnull)
    @test rec["verified_valid_level"] == level
    @test occursin("Retained physical quotient", rec["claim"]) && occursin("Not claimed", rec["claim"])
    path = write_certificate(joinpath(mktempdir(), "certificate.json"), rec)
    ch = roa_check(path)
    @test ch["passed"]
    rec2 = Porthos.JSON3.read(read(path, String), Dict{String,Any})
    rec2["candidate"]["P_columns"][1][1] *= 1 + 1e-12
    @test !roa_check(rec2)["checks"]["candidate_fingerprint"]
end
