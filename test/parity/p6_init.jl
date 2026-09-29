# P6 gate: on every parity case, the Porthos DAE residual at its equilibrium is at most
# 1e-12, and x*, V* match PHPS within 1e-10 or, where PHPS's residual is larger, within the
# distance that residual allows (||J^+|| r, computed here). The parameters initialisation
# sets are compared with PHPS's too.
#
# r is Porthos's residual evaluated at PHPS's point (x*_PHPS, V*_PHPS), on the rows the
# equilibrium solve uses: it contains PHPS's own residual (its multi-pass initialisation is
# not tight) and the effect of the small set-point differences. J^+ is the pseudo-inverse
# of the reduced Jacobian at Porthos's solution (same rows and unknowns as the solve).
#
# One named exception to the 1e-12 residual (user decision, 2026-10-01): the measured-current
# lag rows Id_meas and Iq_meas of GFM_DROOP_PHTRUE may reach 2e-12. Their row
# (I - I_meas)/Tmeas_i carries the rounding of I = (u_out - V)/(j Zseries) amplified by
# 1/(Zseries Tmeas_i), and u_out depends on I_meas through r_vi = Zseries (loop gain 1), so no
# floating-point I_meas zeroes both rows (droop: 1.42e-12, one to two ulps of |I| = 7 pu).
# The test checks that every row above 1e-12 is one of these.

const P6_RES_TOL = 1e-12
const P6_RES_EXCEPTION_TOL = 2e-12
p6_excepted(c, state) = model_type(c) == "GFM_DROOP_PHTRUE" && state in ("Id_meas", "Iq_meas")
const P6_MATCH = 1e-10
const P6_PARAM_ATOL = 1e-8
const INIT_KEYS = ("Efd0", "Tm0", "PFD_REF", "Vref", "PM_REF", "Pref", "V0", "Vini",
                   "p_set", "u_set", "q_set", "v_set", "PSET_REF", "V_nom")

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

            # residual at Porthos's equilibrium: 1e-12 on every row but the named exception
            f, g = dae_residual(sys, r.x, r.V)
            @test maximum(abs, g) <= P6_RES_TOL
            excepted = falses(nd)
            for (k, c) in enumerate(sys.comps), (j, st) in enumerate(state_names(c))
                excepted[sys.offsets[k] + j - 1] = p6_excepted(c, st)
            end
            @test all(abs(f[i]) <= P6_RES_TOL for i in 1:nd if !excepted[i])
            @test all(abs(f[i]) <= P6_RES_EXCEPTION_TOL for i in 1:nd if excepted[i])
            above = [state_names(sys)[i] for i in 1:nd if abs(f[i]) > P6_RES_TOL]
            isempty(above) || @info "P6 $cname: rows above 1e-12 (all within the droop " *
                                    "current-lag exception): $(join(above, ", "))"

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
            comps = Dict(Porthos.name(c) => c for c in sys.comps)
            for c in d[:components], (k, v) in c[:params]
                string(k) in INIT_KEYS || continue
                # the value Porthos simulates with (set by its initialisation, or the case's)
                ours = Porthos.param_value(param_dict(comps[string(c[:name])])[string(k)])
                worst = max(worst, abs(ours - Float64(v)))
                @test abs(ours - Float64(v)) <= P6_PARAM_ATOL
            end
            @info "P6 $cname: largest difference of an initialised parameter from PHPS: $worst"
        end
    end
end
