# P10 gate: the port-Hamiltonian audits reproduce PHPS's records on the IEEE-39 base case,
# at PHPS's initialised state and parameters (pack section `dae`).
#
# - Shifted storage (`records/certificates/hs_decay__base_loadfix.json`, PHPS
#   `study/src/audit_hs_decay39.py`): 20 reservoirs removed, 171 section coordinates, rank of
#   Hess H 54, 41 positive eigenvalues of sym(S A) on the section, the extreme eigenvalues
#   within P10_EIG_RTOL (PHPS differenced its Jacobian and Hessian with steps of 1e-6;
#   Porthos's are exact), and the exact H_s and dH_s/dt along the top eigenvector within
#   P10_EXACT_RTOL (the eigenvector's sign is arbitrary, so +eps and -eps are matched as a
#   pair). The model is assembled without PHPS's C++ constant rounding, as PHPS's audit
#   evaluates its Python solver.
# - Governor ports (`records/model_review/audit_controller_kyp.json` and
#   `audit_governor_nonpassivity_exact.json`): every IEEEG1, linearised from its own
#   equations at speed -> Tm: G(0), the crossing of Re(-G(jw)), its minimum on the same grid,
#   -G(2j) and -G(j 15/4) within P10_PORT_RTOL; and each certified non-passive (Re(-G(jw)) < 0
#   by far more than the evaluation error), which makes the positive-real (KYP) LMI
#   infeasible: the PHPS records say "infeasible" (SCS) and "negative" (exact arithmetic).

const P10_EIG_RTOL = 1e-8
const P10_EXACT_RTOL = 1e-6
const P10_PORT_RTOL = 1e-12

@testset "P10 port-Hamiltonian audits" begin
    d = Porthos.pack_json(PACK, "dae/base.json")
    case = load_case(case_path(string(d[:system])))
    sc = load_scenario(case_path(string(d[:scenario])))
    sys = assemble(case, sc; init_params = phps_init_params(case, d), phps_rounding = false)
    y0 = phps_initial_state(d)
    x, V = y0[1:sys.n_diff], y0[sys.n_diff + 1:end]

    @testset "shifted storage" begin
        rec = Porthos.pack_json(PACK, "records/certificates/hs_decay__base_loadfix.json")
        a = shifted_storage_audit(sys, x, V)
        @test length(a["reservoir_states"]) == rec[:reservoirs_held_at_equilibrium]
        @test a["n_section"] == rec[:n_states_section]
        @test a["storage_components"] == string.(rec[:storage_components])
        @test a["storage_rank"] == rec[:rank_hessian_H]
        @test a["sym_SA_positive"] == rec[:n_positive_eigs]
        @test isapprox(a["sym_SA_eig_max"], rec[:M_eig_max]; rtol = P10_EIG_RTOL)
        @test isapprox(a["sym_SA_eig_min"], rec[:M_eig_min]; rtol = P10_EIG_RTOL)
        ours, theirs = a["exact_along_top_eigenvector"], rec[:exact_along_top_eigenvector]
        @test length(ours) == length(theirs)
        close(p, q) = isapprox(p["Hs"], q[:Hs]; rtol = P10_EXACT_RTOL) &&
                      isapprox(p["dHs_dt"], q[:dHs_dt]; rtol = P10_EXACT_RTOL)
        for k in 1:2:length(ours)
            @test ours[k]["eps"] == theirs[k][:eps]
            @test (close(ours[k], theirs[k]) && close(ours[k + 1], theirs[k + 1])) ||
                  (close(ours[k], theirs[k + 1]) && close(ours[k + 1], theirs[k]))
        end
        @test all(e -> e["dHs_dt"] > 0, ours)          # H_s increases: not a Lyapunov function
        @info "P10 shifted storage: section $(a["n_section"]), rank $(a["storage_rank"]), " *
              "$(a["sym_SA_positive"]) positive eigenvalues of sym(SA) in " *
              "[$(a["sym_SA_eig_min"]), $(a["sym_SA_eig_max"])] (PHPS $(rec[:M_eig_min]), " *
              "$(rec[:M_eig_max])); held: $(join(a["held_states"], ", "))"
    end

    @testset "governor ports" begin
        kyp = Porthos.pack_json(PACK, "records/model_review/audit_controller_kyp.json")
        ex = Porthos.pack_json(PACK, "records/model_review/audit_governor_nonpassivity_exact.json")
        ins, _ = component_io(sys, x, V)
        govs = [k for (k, c) in enumerate(sys.comps) if model_type(c) == "IEEEG1_PHTRUE"]
        @test sort([Porthos.name(sys.comps[k]) for k in govs]) == sort(string.(keys(kyp[:governors])))
        w0 = parse(Float64, "3.75")                     # 15/4 rad/s, as PHPS's exact audit
        @test string(ex[:w0_rad_s]) == "15/4"
        crossings = Float64[]
        for k in govs
            c = sys.comps[k]
            g = Symbol(Porthos.name(c))
            r = sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
            m = port_model(c, x[r], ins[k]; input = "omega", output = "Tm")
            ref, e = kyp[:governors][g], ex[:governors][g]
            @test isapprox(real(transfer(m, 0.0)), ref[:G0]; rtol = P10_PORT_RTOL)
            cr = real_part_crossings(m)
            append!(crossings, cr)
            @test length(cr) == length(ref[:re_minus_G_sign_changes_rad_s])
            @test all(isapprox.(cr, Float64.(ref[:re_minus_G_sign_changes_rad_s]); rtol = P10_PORT_RTOL))
            pc = passivity_certificate(m)
            @test isapprox(pc.min_re, ref[:min_re_minus_G]; rtol = P10_PORT_RTOL)
            @test isapprox(pc.at, ref[:at_rad_s]; rtol = P10_PORT_RTOL)
            h2 = -transfer(m, 2im)
            @test isapprox(real(h2), ref[:minus_G_at_2_rad_s][1]; rtol = P10_PORT_RTOL)
            @test isapprox(imag(h2), ref[:minus_G_at_2_rad_s][2]; rtol = P10_PORT_RTOL)
            h = -transfer(m, w0 * im)
            @test isapprox(real(h), e[:Re_H_jw0]; rtol = P10_PORT_RTOL)
            @test isapprox(imag(h), e[:Im_H_jw0]; rtol = P10_PORT_RTOL)
            @test pc.H0 > 0 && e[:H0_exact_is_positive] == true
            # non-passive, so the positive-real LMI is infeasible (as PHPS's solver and exact
            # arithmetic found)
            @test pc.nonpassive && real(h) < 0
            @test ref[:kyp_lmi][:status] == "infeasible" && e[:Re_H_jw0_sign] == "negative"
        end
        @info "P10 governor ports: $(length(govs)) IEEEG1 sets, each non-passive at speed -> Tm; " *
              "Re(-G(jw)) changes sign at $(extrema(crossings)) rad/s"
    end
end
