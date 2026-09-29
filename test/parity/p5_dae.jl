# P5 gate: f and g equal the parity pack (PHPS's compiled DAE residual) at 50 random states
# per case within 1e-12, fault off and fault on; the state layout, the wiring and the
# network constants are PHPS's.
#
# Until P6, the parameters that PHPS's initialisation sets are taken from the pack. The test
# checks that no other parameter differs from Porthos's own processing.

const P5_TOL = 1e-12
const INIT_SET_KEYS = Set(["Efd0", "Tm0", "PFD_REF", "Vref", "PM_REF", "Pref", "V0", "Vini"])

# a PHPS C++ wiring expression as an InputSource
function phps_source(expr, names::Dict{String,Int}, outputs)
    expr === nothing && return Porthos.InputSource(Porthos.SRC_ZERO)
    e = string(expr)
    m = match(r"^(Vd|Vq|Vterm)_net\[(\d+)\]$", e)
    if m !== nothing
        kind = m[1] == "Vd" ? Porthos.SRC_VD : m[1] == "Vq" ? Porthos.SRC_VQ : Porthos.SRC_VTERM
        return Porthos.InputSource(kind, 0.0, parse(Int, m[2]) + 1, 0)
    end
    m = match(r"^outputs_(.+)\[(\d+)\]$", e)
    m !== nothing && return Porthos.InputSource(Porthos.SRC_OUTPUT, 0.0, names[m[1]],
                                                 parse(Int, m[2]) + 1)
    m = match(r"^v(d|q)_dq_(.+)$", e)
    m !== nothing && return Porthos.InputSource(m[1] == "d" ? Porthos.SRC_DQ_VD :
                                                 Porthos.SRC_DQ_VQ, 0.0, names[m[2]], 0)
    return Porthos.const_source(parse(Float64, e))
end

@testset "P5 DAE parity" begin
    @test Porthos.has_section(PACK, "dae")
    cases = Dict(c.name => c for c in Porthos.pack_cases(PACK))
    for file in sort([string(k) for k in keys(PACK.manifest[:files])
                      if startswith(string(k), "dae/")])
        d = Porthos.pack_json(PACK, file)
        cname = string(d[:case])
        @testset "$cname" begin
            case = load_case(case_path(string(d[:system])))
            sc = load_scenario(case_path(string(d[:scenario])))
            missing_types = unique(s.type for s in case.components
                                   if !haskey(Porthos.COMPONENT_CONSTRUCTORS, s.type))
            if !isempty(missing_types)
                @info "P5 $cname pending: model types not ported yet: $(join(missing_types, ", "))"
                @test_skip isempty(missing_types)
                continue
            end

            # initialisation parameters from PHPS, and nothing else
            init = Dict{String,Dict{String,Float64}}()
            for c in d[:components]
                nm = string(c[:name])
                ours = build_component(case, component(case, nm))
                diff = Dict{String,Float64}()
                for (k, v) in c[:params]
                    key = string(k)
                    pd = param_dict(ours)
                    if !haskey(pd, key) || !(pd[key] isa Real) || Float64(pd[key]) !== Float64(v)
                        diff[key] = Float64(v)
                    end
                end
                @test issubset(keys(diff), INIT_SET_KEYS)
                init[nm] = diff
            end
            sys = assemble(case, sc; init_params = init)

            @testset "layout and wiring" begin
                @test sys.n_diff == Int(d[:n_diff])
                @test sys.delta_coi == Int(d[:delta_coi_index])
                @test sys.net.bus_ids == Int.(d[:bus_ids])
                @test [Porthos.name(c) for c in sys.comps] == [string(c[:name]) for c in d[:components]]
                @test sys.offsets == [Int(c[:offset]) for c in d[:components]]
                for (c, r) in zip(sys.comps, d[:components])
                    @test state_names(c) == string.(r[:states])
                    @test input_names(c) == string.(r[:inputs])
                    @test output_names(c) == string.(r[:outputs])
                    @test component_role(c) === Symbol(r[:role])
                end
                names = Dict(Porthos.name(c) => k for (k, c) in enumerate(sys.comps))
                nbad = 0
                for (k, c) in enumerate(sys.comps), (j, p) in enumerate(input_names(c))
                    ref = phps_source(d[:wiring][Symbol(Porthos.name(c) * "." * p)], names, nothing)
                    nbad += sys.sources[k][j] != ref
                end
                @test nbad == 0
            end

            @testset "network constants" begin
                Y = Matrix(ybus_dae(case; v0 = sys.pf.V))
                @test close_to(Y, Porthos.pack_matrix(d[:Y_full]); atol = P5_TOL, rtol = P5_TOL)
                la = load_admittances(case; v0 = sys.pf.V,
                                      load_params = Dict(Porthos.name(c) => param_dict(c)
                                                         for c in sys.comps
                                                         if component_role(c) === :load))
                for k in (:G, :B, :P, :Q, :kpf, :kqf)
                    @test close_to(getfield(la, k), Float64.(d[Symbol("load_", k)]);
                                   atol = P5_TOL, rtol = P5_TOL)
                end
                @test [Porthos.name(sys.comps[k]) for k in sys.coi_members] ==
                      [string(c[:component]) for c in d[:coi]]
                @test sys.coi_weights == [Porthos._c6f(Float64(c[:weight])) for c in d[:coi]]
                @test sys.omega_b == parse_param_expr(string(d[:omega_b_sys]))
                @test count(sys.slack) == length(d[:slack])
                @test length(sys.faults) == length(d[:faults])
                for (f, r) in zip(sys.faults, d[:faults])
                    @test f.index == Int(r[:index])
                    @test f.g == Porthos._r10(Float64(r[:g]), true)
                    @test f.b == Porthos._r10(Float64(r[:b]), true)
                end
            end

            @testset "residual" begin
                vvec(e) = vec(permutedims(hcat(Float64.(e[:Vd]), Float64.(e[:Vq]))))
                e = d[:equilibrium]
                f, g = dae_residual(sys, Float64.(e[:x]), vvec(e))
                @test close_to(f, Float64.(e[:f]); atol = P5_TOL, rtol = P5_TOL)
                @test close_to(g, Float64.(e[:g]); atol = P5_TOL, rtol = P5_TOL)
                nbad = 0
                ndiff = 0
                nval = 0
                for s in d[:samples]
                    x, V = Float64.(s[:x]), Float64.(s[:V])
                    for (on, fk, gk) in ((false, :f, :g), (true, :f_fault, :g_fault))
                        haskey(s, fk) || continue
                        f, g = dae_residual(sys, x, V; faults_on = on)
                        rf, rg = Float64.(s[fk]), Float64.(s[gk])
                        ok = close_to(f, rf; atol = P5_TOL, rtol = P5_TOL) &&
                             close_to(g, rg; atol = P5_TOL, rtol = P5_TOL)
                        nbad += !ok
                        ndiff += bitdiff(f, rf) + bitdiff(g, rg)
                        nval += length(f) + length(g)
                    end
                end
                @test length(d[:samples]) == 50
                @test nbad == 0
                @info "P5 $cname: $ndiff of $nval residual values differ from PHPS's compiled C++ in the last bits (sin/cos of different math libraries); the rest are bit-identical"
            end
        end
    end
end
