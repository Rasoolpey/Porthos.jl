# P6 gate: on every parity case, the Porthos DAE residual at its equilibrium is at most
# 1e-12, and x*, V* match PHPS within 1e-10 or, where PHPS's residual is larger, within the
# distance that residual allows (||J^+|| r, computed here). The parameters initialisation
# sets are compared with PHPS's too.
#
# r is Porthos's residual evaluated at PHPS's point (x*_PHPS, V*_PHPS), on the rows the
# equilibrium solve uses: it contains PHPS's own residual (its multi-pass initialisation is
# not tight) and the effect of the small set-point differences. J^+ is the pseudo-inverse
# of the reduced Jacobian at Porthos's solution (same rows and unknowns as the solve).

const P6_RES_TOL = 1e-12
const P6_MATCH = 1e-10
const P6_PARAM_ATOL = 1e-8
const INIT_KEYS = ("Efd0", "Tm0", "PFD_REF", "Vref", "PM_REF", "Pref", "V0", "Vini")

@testset "P6 initialisation parity" begin
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
                @info "P6 $cname pending: model types not ported yet: $(join(missing_types, ", "))"
                @test_skip isempty(missing_types)
                continue
            end
            r = solve_equilibrium(case, sc)
            sys = r.sys
            nd = sys.n_diff

            # the first pass reproduces PHPS's initialisation before its refinement
            @test close_to(r.x_initial, Float64.(d[:equilibrium][:x_initial]); atol = 1e-12,
                           rtol = 1e-12)

            # residual at Porthos's equilibrium
            @test r.residual <= P6_RES_TOL

            # distance to PHPS's equilibrium
            e = d[:equilibrium]
            xp = Float64.(e[:x])
            Vp = vec(permutedims(hcat(Float64.(e[:Vd]), Float64.(e[:Vq]))))
            dist = max(maximum(abs, r.x .- xp), maximum(abs, r.V .- Vp))
            res_states = Porthos._reservoir_states(sys)
            fixed = Set(res_states)
            push!(fixed, sys.delta_coi)
            g0 = Porthos._gauge_state(sys)
            g0 > 0 && push!(fixed, g0)
            n = nd + Porthos.nalg(sys)
            cols = [j for j in 1:n if !(j in fixed)]
            rows = [i for i in 1:n if !(i in res_states)]
            F(z) = vcat(dae_residual(sys, view(z, 1:nd), view(z, nd + 1:n))...)
            J = ForwardDiff.jacobian(F, vcat(r.x, r.V))[rows, cols]
            Jpinv_norm = opnorm(pinv(J), Inf)
            r_phps = maximum(abs, F(vcat(xp, Vp))[rows])
            allowed = max(P6_MATCH, Jpinv_norm * r_phps)
            @test dist <= allowed
            @info "P6 $cname: residual $(r.residual); |(x*,V*) - PHPS| = $dist, allowed " *
                  "$allowed (||J^+|| = $Jpinv_norm, residual at PHPS's point = $r_phps; " *
                  "PHPS's own residual there = $(max(maximum(abs, Float64.(e[:f])), maximum(abs, Float64.(e[:g])))))"

            # the parameters initialisation sets
            worst = 0.0
            for c in d[:components], (k, v) in c[:params]
                string(k) in INIT_KEYS || continue
                ours = r.init_params[string(c[:name])][string(k)]
                worst = max(worst, abs(ours - Float64(v)))
                @test abs(ours - Float64(v)) <= P6_PARAM_ATOL
            end
            @info "P6 $cname: largest difference of an initialised parameter from PHPS: $worst"
        end
    end
end
