# P4 gate (per type): at the 200 random states of the parity pack, the limiter branches agree
# first (both sides exercised); then rhs, outputs, injection, H and grad H agree with
# rtol = 1e-12 and the per-quantity atol. The contract entry is identical.
#
# Parameters set by PHPS's initialisation (reservoir references, load Vini, Pref/Vref) are
# taken from the pack here; P6 checks that Porthos's own initialisation reproduces them.

const P4_RTOL = 1e-12

# per-quantity absolute tolerance: 1e-12 of the largest magnitude that quantity takes
function p4_atol(refs::Vector{Vector{Float64}})
    (isempty(refs) || isempty(refs[1])) && return Float64[]
    n = length(refs[1])
    return [1e-12 * max(maximum(abs(r[i]) for r in refs), 1e-300) for i in 1:n]
end

# the switching sites: state-dependent, not inside a loop (a root solve's comparisons), and
# reachable by the pack's samples
state_dependent(sites) = Set(Int(s[:id]) for s in sites
                             if !s[:param_only] && !get(s, :in_loop, false) &&
                                !haskey(s, :unreached))

phps_modes(branches, dep) = Bool[b[2] for b in branches if Int(b[1]) in dep]

numeric_params(d) = Dict(string(k) => Float64(v) for (k, v) in d
                         if v isa Real && !(v isa Bool))

@testset "P4 component parity" begin
    @test Porthos.has_section(PACK, "components")
    cases = Dict(c.name => c for c in Porthos.pack_cases(PACK))
    for file in sort(filter(f -> endswith(f, ".json"),
                            [string(k) for k in keys(PACK.manifest[:files])]))
        startswith(file, "components/") || continue
        rec = Porthos.pack_json(PACK, file)
        ctype = string(rec[:type])
        @testset "$ctype" begin
            case = load_case(case_path(cases[string(rec[:case])].system))
            comps = Dict{String,AbstractComponent}()
            for (nm, inst) in rec[:instances]
                # a parameter variant is its case instance with a few parameters overridden
                base = string(get(inst, :variant_of, nm))
                c = build_component(case, component(case, base))
                set = Dict(string(k) => Float64(v) for (k, v) in inst[:init_set])
                for (k, v) in get(inst, :overrides, Dict())
                    set[string(k)] = Float64(v)
                end
                comps[string(nm)] = with_params(c, set)
            end
            c1 = first(values(comps))
            @test model_type(c1) == ctype
            @test state_names(c1) == string.(rec[:states])
            @test input_names(c1) == string.(rec[:inputs])
            @test output_names(c1) == string.(rec[:outputs])

            # parameters: exactly the values PHPS computes with
            @testset "parameters" begin
                for (nm, inst) in rec[:instances]
                    ours = numeric_params(param_dict(comps[string(nm)]))
                    theirs = Dict(string(k) => Float64(v) for (k, v) in inst[:params_used])
                    @test ours == theirs
                end
            end

            samples = rec[:samples]
            @test length(samples) == 200
            dep_step = state_dependent(rec[:sites][:step])
            dep_out = state_dependent(rec[:sites][:out])

            # 1. branch agreement, with both sides of every site exercised
            @testset "branches" begin
                nbad = 0
                for s in samples
                    c = comps[string(s[:component])]
                    x, u = Float64.(s[:x]), Float64.(s[:u])
                    ok = modes(c, x, u; kernel = :step) == phps_modes(s[:branches_step], dep_step) &&
                         modes(c, x, u; kernel = :out) == phps_modes(s[:branches_out], dep_out)
                    nbad += !ok
                end
                @test nbad == 0
                for (kernel, dep) in ((:step, dep_step), (:out, dep_out)), k in dep
                    seen = Set(b[2] for s in samples for b in s[Symbol("branches_", kernel)]
                               if Int(b[1]) == k)
                    @test seen == Set([false, true])
                end
                # sites the pack lists as unreached (with the reason) were never reached
                for kernel in (:step, :out), st in rec[:sites][kernel]
                    haskey(st, :unreached) || continue
                    @test !any(Int(b[1]) == Int(st[:id]) for s in samples
                               for b in s[Symbol("branches_", kernel)])
                    @info "$ctype: $kernel site $(st[:id]) ($(st[:text])) not reached: $(st[:unreached])"
                end
            end

            # 2. values
            @testset "values" begin
                refs(key) = [Float64.(s[key]) for s in samples]
                atol_dx = p4_atol(refs(:dxdt))
                atol_y1 = p4_atol(refs(:outputs_out))
                atol_y2 = p4_atol(refs(:outputs_step))
                atol_g = p4_atol(refs(:grad_H))
                atol_H = 1e-12 * max(maximum(abs(Float64(s[:H])) for s in samples), 1e-300)
                bad = Dict(:rhs => 0, :out => 0, :step_out => 0, :inj => 0, :H => 0, :gradH => 0)
                exact = 0
                for s in samples
                    c = comps[string(s[:component])]
                    x, u = Float64.(s[:x]), Float64.(s[:u])
                    dx = rhs!(zeros(nstates(c)), c, x, u)
                    y1 = outputs!(zeros(noutputs(c)), c, x, u)
                    y2 = step_outputs!(zeros(noutputs(c)), c, x, u)
                    ref_dx, ref_y1 = Float64.(s[:dxdt]), Float64.(s[:outputs_out])
                    close_to(dx, ref_dx; atol = atol_dx, rtol = P4_RTOL) || (bad[:rhs] += 1)
                    close_to(y1, ref_y1; atol = atol_y1, rtol = P4_RTOL) || (bad[:out] += 1)
                    close_to(y2, Float64.(s[:outputs_step]); atol = atol_y2, rtol = P4_RTOL) ||
                        (bad[:step_out] += 1)
                    if bus(c) !== nothing
                        Ire, Iim = injection(c, x, (u[1], u[2]))
                        close_to([Ire, Iim], ref_y1[1:2]; atol = atol_y1[1:2], rtol = P4_RTOL) ||
                            (bad[:inj] += 1)
                    end
                    close_to(hamiltonian(c, x), Float64(s[:H]); atol = atol_H, rtol = P4_RTOL) ||
                        (bad[:H] += 1)
                    close_to(grad_hamiltonian(c, x), Float64.(s[:grad_H]); atol = atol_g,
                             rtol = P4_RTOL) || (bad[:gradH] += 1)
                    exact += bitdiff(dx, ref_dx) == 0
                end
                for (k, v) in bad
                    @test (k, v) == (k, 0)
                end
                @info "$ctype: rhs bit-identical to PHPS at $exact/200 samples"
            end

            @testset "contract" begin
                cs = default_contracts()
                @test Porthos.json_identical(contract(c1).raw,
                                             cs.raw[:contracts][Symbol(contract_key(ctype))])
            end
        end
    end
end
