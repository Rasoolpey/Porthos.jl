using ForwardDiff

@testset "branch primitives" begin
    rec = ModeLog()
    @test Porthos.clamp_mode(rec, 5.0, 0.0, 1.0) === 1.0
    @test Porthos.clamp_mode(rec, -5.0, 0.0, 1.0) === 0.0
    @test Porthos.clamp_mode(rec, 0.5, 0.0, 1.0) === 0.5
    # PHPS order: upper test, then lower test (on the updated value)
    @test rec.decisions == [true, false, false, true, false, false]

    rec = ModeLog()
    @test Porthos.nonwindup(rec, 1.0, 2.0, 0.0, 1.0) === 0.0     # at the top, pushing up
    @test Porthos.nonwindup(rec, 1.0, -2.0, 0.0, 1.0) === -2.0   # at the top, leaving
    @test Porthos.nonwindup(rec, 0.0, -2.0, 0.0, 1.0) === 0.0    # at the bottom, pushing down
    @test rec.decisions == [true, false, false, false, false, true]

    rec = ModeLog()
    @test Porthos.outband_relax(rec, 0.0, 1.5, 0.0, 1.0, 2.0) === -0.25
    @test Porthos.outband_relax(rec, 0.0, -1.0, 0.0, 1.0, 2.0) === 0.5
    @test rec.decisions == [true, false, false, true]

    rec = ModeLog()
    @test Porthos.guard_min(rec, 1e-9, 1e-6) === 1e-6
    @test Porthos.guard_min(rec, 1.0, 1e-6) === 1.0
    @test rec.decisions == [true, false]

    # NoModes records nothing and keeps the type
    d = ForwardDiff.Dual(0.5, 1.0)
    @test Porthos.clamp_mode(NoModes(), d, 0.0, 1.0) === d
    @test Porthos.clamp_mode(NoModes(), ForwardDiff.Dual(2.0, 1.0), 0.0, 1.0) isa ForwardDiff.Dual
end

# A representative operating point per model type (near the base-case equilibria).
const UNIT_POINTS = Dict(
    "GENROU_PHTRUE" => ([0.6, 1.001, 1.0, 0.95, 0.2, -0.3], [0.9, 0.3, 5.0, 2.2]),
    "GENSAL_PHTRUE" => ([0.6, 0.999, 1.0, 0.95, -0.2], [0.9, 0.3, 2.5, 1.2]),
    "IEEET1_PHTRUE" => ([1.0, 0.5, 2.2, 2.2, 1000.0], [0.98, 0.94, 0.0, 0.003]),
    "IEEEG1_PHTRUE" => ([5.2, 5.2, 5.2, 5.2, 5.2, 5.2, 1000.0], [1.001, 2.04, 0.0]),
    "IEEEG3_PHTRUE" => ([0.0, 2.5, 2.5, 2.5, 1000.0], [0.999, 0.1, 0.0]),
    "COMPLEXLOAD" => ([-0.02], [0.95, -0.3]),
)

alloc_rhs(dx, c, x, u) = @allocated rhs!(dx, c, x, u)
alloc_out(y, c, x, u) = @allocated outputs!(y, c, x, u)

@testset "component models" begin
    case = load_case(case_path("IEEE39Bus_PF/system_phtrue.json"))
    for (type, (x, u)) in UNIT_POINTS
        spec = case.components[findfirst(s -> s.type == type, case.components)]
        c = build_component(case, spec)
        @testset "$type" begin
            n = nstates(c)
            @test length(x) == n && length(u) == ninputs(c)
            dx = zeros(n)
            y = zeros(noutputs(c))
            # type-stable and allocation-free
            @inferred rhs!(dx, c, x, u)
            @inferred outputs!(y, c, x, u)
            @test alloc_rhs(dx, c, x, u) == 0
            @test alloc_out(y, c, x, u) == 0
            @test all(isfinite, dx) && all(isfinite, y)
            # generic in the number type: ForwardDiff Jacobian = finite differences
            f(x) = rhs!(similar(x, n), c, x, u)
            J = ForwardDiff.jacobian(f, x)
            h = 1e-6
            Jfd = hcat([(f(x .+ h .* (1:n .== j)) .- f(x .- h .* (1:n .== j))) ./ (2h)
                        for j in 1:n]...)
            @test isapprox(J, Jfd; rtol = 1e-5, atol = 1e-6)
            @test modes(c, ForwardDiff.Dual.(x, 1.0), u) == modes(c, x, u)
            # storage
            g = grad_hamiltonian(c, x)
            @test isapprox(g, ForwardDiff.gradient(z -> hamiltonian(c, z), x); rtol = 1e-12,
                           atol = 1e-14)
            @test contract(c).model == contract_key(type)
        end
    end
end

@testset "Norton admittance matches the Y-bus stamps" begin
    case = load_case(case_path("IEEE39Bus_PF/system_phtrue.json"))
    for st in norton_stamps(case)
        c = build_component(case, component(case, st.component))
        @test norton_admittance(c) == Porthos._py_cdiv(1.0, complex(st.ra, st.xd_pp))
        @test bus(c) == st.bus
    end
    load = build_component(case, component(case, "CLOAD_3"))
    @test norton_admittance(load) === nothing
end

@testset "retired and unknown types are rejected" begin
    @test_throws UnsupportedModelError build_component("GENROU_PHS", "G", Porthos.ParamDict())
    @test_throws UnsupportedModelError build_component("NOPE", "G", Porthos.ParamDict())
end

@testset "grid-forming converter models" begin
    cases = ("GFM_VSM_PHTRUE" => "IEEE39Bus_PF_gfm-vsm/system_gfm_vsm.json",
             "GFM_DROOP_PHTRUE" => "IEEE39Bus_PF_gfm-droop/system_gfm_droop.json",
             "GFM_VOC_PHTRUE" => "IEEE39Bus_PF_gfm-voc/system_gfm_voc.json")
    # parameters that switch on the branches the cases keep off
    modes_on = Dict("GFM_VSM_PHTRUE" => Dict("pf_frame" => 1.0, "K_field" => 0.5, "Dq" => 5.0,
                                             "n_vi" => 2.0, "D_direct" => 2.0),
                    "GFM_DROOP_PHTRUE" => Dict("pf_frame" => 1.0, "adapt_droop" => 1.0,
                                               "adapt_vi" => 1.0, "vi_mode" => 1.0,
                                               "x_vi" => 0.01),
                    "GFM_VOC_PHTRUE" => Dict("pf_frame" => 1.0, "i_ref_max" => 5.0))
    for (type, file) in cases
        case = load_case(case_path(file))
        pf = solve_powerflow(case)
        x0, init, _ = first_pass(case; pf)
        spec = case.components[findfirst(s -> s.type == type, case.components)]
        k = findfirst(s -> s.name == spec.name, case.components)
        c = with_params(build_component(case, spec), init[spec.name])
        o = 1 + sum(nstates(build_component(case, s)) for s in case.components[1:k - 1]; init = 0)
        x = x0[o:o + nstates(c) - 1]
        i = pf.spec.net.index[bus(c)]
        V = pf.V[i] * cis(pf.theta[i])
        u = [real(V); imag(V); ones(ninputs(c) - 2)]
        @testset "$type" begin
            n = nstates(c)
            # the first pass is an equilibrium of the converter at its operating point
            dx = rhs!(zeros(n), c, x, u)
            @test maximum(abs, dx) <= 1e-9
            I = converter_current(c, x, V)
            _, pk, _ = converter_init(c, V, I)
            @test isapprox(pk["p_set"], init[spec.name]["p_set"]; rtol = 1e-12)
            @test norton_admittance(c) == Porthos._py_cdiv(1.0, complex(0.0, c.p.Zseries))
            for (label, cc) in (("case", c), ("modes on", with_params(c, modes_on[type])))
                # perturb away from the equilibrium, into the virtual-impedance regime too
                for xp in (x .* (1 .+ 1e-3 .* sin.(1:n)), x .+ 0.05 .* cos.(1:n))
                    y = zeros(noutputs(cc))
                    dxp = zeros(n)
                    @inferred rhs!(dxp, cc, xp, u)
                    @inferred outputs!(y, cc, xp, u)
                    @test alloc_rhs(dxp, cc, xp, u) == 0
                    @test alloc_out(y, cc, xp, u) == 0
                    @test all(isfinite, dxp) && all(isfinite, y)
                    f(z) = rhs!(similar(z, n), cc, z, u)
                    J = ForwardDiff.jacobian(f, xp)
                    h = 1e-7
                    Jfd = hcat([(f(xp .+ h .* (1:n .== j)) .- f(xp .- h .* (1:n .== j))) ./ (2h)
                                for j in 1:n]...)
                    @test isapprox(J, Jfd; rtol = 1e-4, atol = 1e-4 * max(1.0, maximum(abs, J)))
                    @test modes(cc, ForwardDiff.Dual.(xp, 1.0), u) == modes(cc, xp, u)
                    g = grad_hamiltonian(cc, xp)
                    @test isapprox(g, ForwardDiff.gradient(z -> hamiltonian(cc, z), xp);
                                   rtol = 1e-12, atol = 1e-14)
                end
            end
            @test contract(c).model == contract_key(type)
            # the passivity-based VOC has no parity reference and is rejected
            if type == "GFM_VOC_PHTRUE"
                @test_throws ArgumentError with_params(c, Dict("pvoc_mode" => 1.0))
            end
        end
    end
end
