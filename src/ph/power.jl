# Nonlinear port-power residual audit (roadmap P10; PHPS work package 1 item 3; TODO.md
# step 4). At any state (x, V) with KCL satisfied it reconstructs, independently of each
# other:
#
#   storage rate   dH_k/dt = grad H_k(x_k)' f_k(x, V), from the component's declared storage
#                  and its own right-hand side;
#   supply         the sum of the port powers the contract declares (`power_into_component`,
#                  power into the component), each evaluated from the component's runtime
#                  signals (inputs, outputs, states, parameters) by the expression reader;
#   residual       supply - rate: the power the component dissipates (> 0) or creates (< 0)
#                  given its declared storage and ports;
#
# and the network side:
#
#   injected       sum over buses of V . I_inj (every KCL injection: Norton currents and
#                  load corrections), from the components' outputs;
#   absorbed       V' G V split into the Norton admittances, the rest of the network (lines,
#                  shunts, constant-impedance loads) and the active fault shunts, plus the
#                  frequency-dependent load terms; injected - absorbed = V . g (KCL);
#   terminal       per machine, the power out of its terminal, V . I_t (its It output) against
#                  V . I_norton - Re(y_n) |V|^2;
#   loads          per ComplexLoad, its declared P against the power its bus draws,
#                  G_load |V|^2 - V . I_corr.
#
# Nothing here decides which storage or which ports are right; it measures how the declared
# ones balance, with signs and units, so that storage candidates (Route A', V_ext) start from
# reliable accounting.

# the value of every name a contract expression can use, for component k
function _signal_env(c::AbstractComponent, xk, ins, outs)
    env = Dict{String,Float64}()
    for (k, v) in param_dict(c)
        val = try
            param_value(v)
        catch
            nothing
        end
        val === nothing || (env[string(k)] = val)
    end
    for (n, v) in zip(state_names(c), xk)
        env[n] = v
    end
    for (n, v) in zip(input_names(c), ins)
        env[n] = v
    end
    for (n, v) in zip(output_names(c), outs)
        env[n] = v
    end
    return env
end

"""
    component_power(sys, x, V; contracts = default_contracts(), faults_on = false)
        -> Vector{Dict}

For every component at `(x, V)`: `storage_rate` (`grad H' f` over its own states),
`ports` (name, declared expression, value or `nothing` when the contract gives no
expression), `supply` (the evaluable ports), `complete` (every power port evaluable),
`residual = supply - storage_rate`. Machines also get `kinetic_rate_shift = 2 H omega'`
(the rate of the physical kinetic storage `H omega^2` minus that of the shifted
`H (omega - 1)^2` the model declares) and, against the physical kinetic storage, the split
of the residual into `mechanical_residual = Tm - Pe - d(H omega^2)/dt`,
`magnetic_residual = Efd i_fd - dH_mag/dt` and `terminal_residual = Pe - V . I_norton`.
Machines with a `rotor_structure` (GENROU) also get the rotor port split `field_supply`,
`stator_exchange = [id, iq]' B_s' Q z`, `rotor_loss = -z' sym(QA) z`, and the errors of the
field colocation, of the linear rotor form against the right-hand side, and of the identity
`dH_mag/dt = field_supply + stator_exchange - rotor_loss`.
"""
function component_power(sys::DAESystem, x::AbstractVector, V::AbstractVector;
                         contracts::ContractSet = default_contracts(), faults_on::Bool = false)
    f, _ = dae_residual(sys, x, V; faults_on)
    ins, outs = component_io(sys, x, V)
    out = Dict{String,Any}[]
    for (k, c) in enumerate(sys.comps)
        r = _state_range(sys, k)
        xk = collect(x[r])
        rate = nstates(c) > 0 ? dot(grad_hamiltonian(c, xk), f[r]) : 0.0
        env = _signal_env(c, xk, ins[k], outs[k])
        entry = get(contracts.entries, contract_key(model_type(c)), nothing)
        ports = Dict{String,Any}[]
        complete = entry !== nothing
        if entry !== nothing
            for p in get(entry.raw, :ports, ())
                haskey(p, :power_into_component) || continue      # a signal port
                expr = String(p[:power_into_component])
                val = try
                    eval_param_expr(Float64, expr, env)
                catch e
                    e isa ParamExprError || rethrow()
                    nothing
                end
                val === nothing && (complete = false)
                push!(ports, Dict{String,Any}("name" => String(p[:name]), "expression" => expr,
                                              "value" => val))
            end
        end
        supply = sum((p["value"] for p in ports if p["value"] !== nothing); init = 0.0)
        d = Dict{String,Any}("component" => name(c), "type" => model_type(c),
                             "storage_rate" => rate, "ports" => ports, "supply" => supply,
                             "complete" => complete, "residual" => supply - rate)
        if component_role(c) === :generator && "omega" in state_names(c) && hasparam(c, "H")
            j = findfirst(==("omega"), state_names(c))
            H, w, dw = param_value(param_dict(c)["H"]), xk[j], f[r[j]]
            d["kinetic_rate_shift"] = 2H * dw
            # split against the physical kinetic storage H omega^2: mechanical side
            # Tm - Pe - d(H omega^2)/dt (= D omega (omega - 1) by the swing equation), magnetic
            # side Efd i_fd - dH_mag/dt (H_mag = declared storage - H (omega - 1)^2), and the
            # terminal part Pe - V . I_norton (-Re(y_n)|V|^2)
            if all(k -> haskey(env, k), ("Tm", "Pe", "Efd", "i_fd", "Vd", "Vq", "Id", "Iq"))
                d["mechanical_residual"] = env["Tm"] - env["Pe"] - 2H * w * dw
                d["magnetic_residual"] = env["Efd"] * env["i_fd"] - (rate - 2H * (w - 1) * dw)
                d["terminal_residual"] = env["Pe"] - (env["Vd"] * env["Id"] + env["Vq"] * env["Iq"])
            end
            # the rotor as a linear port system (`rotor_structure`): dH_mag/dt =
            # field supply + stator exchange - rotor loss, with the stator exchange
            # [id, iq]' B_s' Q z and the loss -z' sym(QA) z
            rs = rotor_structure(c)
            if rs !== nothing
                z = xk[rs.states]
                i = [env[rs.currents[1]], env[rs.currents[2]]]
                Qz = rs.Q * z
                field = env["Efd"] * dot(rs.Bf, Qz)
                exchange = dot(i, rs.Bs' * Qz)
                loss = -dot(z, rs.Q * (rs.A * z))
                dHmag = rate - 2H * (w - 1) * dw
                d["field_supply"] = field
                d["stator_exchange"] = exchange
                d["rotor_loss"] = loss
                d["field_colocation_error"] = abs(dot(rs.Bf, Qz) - env["i_fd"])
                d["rotor_model_error"] = maximum(abs, f[r[rs.states]] .- (rs.A * z .+ rs.Bf .* env["Efd"] .+ rs.Bs * i))
                d["rotor_identity_error"] = abs(dHmag - (field + exchange - loss))
            end
        end
        push!(out, d)
    end
    return out
end

hasparam(c::AbstractComponent, k::AbstractString) = haskey(param_dict(c), k)

"""
    network_power(sys, x, V; faults_on = false) -> Dict

The network's power balance at `(x, V)` (see the comment at the top of `src/ph/power.jl`):
`injected`, `absorbed` (`norton`, `network`, `fault`, `frequency_loads`), `kcl_residual_power
= V . g` (their difference), and per machine and per ComplexLoad the terminal identities.
"""
function network_power(sys::DAESystem, x::AbstractVector, V::AbstractVector; faults_on::Bool = false)
    nb = nbus(sys)
    Vd = [V[2i - 1] for i in 1:nb]
    Vq = [V[2i] for i in 1:nb]
    V2 = Vd .^ 2 .+ Vq .^ 2
    _, g = dae_residual(sys, x, V; faults_on)
    ins, outs = component_io(sys, x, V)
    injected = 0.0
    inj_by = Dict{String,Float64}()
    for (k, c) in enumerate(sys.comps)
        b, jd, jq = sys.inj[k]
        b > 0 || continue
        p = Vd[b] * outs[k][jd] + Vq[b] * outs[k][jq]
        injected += p
        inj_by[name(c)] = p
    end
    # Y-bus conductance part (B cancels: the Y-bus is symmetric), split into Norton stamps
    Gn = zeros(nb)
    stamps = Dict(s.component => s for s in norton_stamps(sys.case))
    for c in sys.comps
        s = get(stamps, name(c), nothing)
        s === nothing && continue
        z = complex(s.ra, s.xd_pp)
        abs(z) < 1e-6 && (z = complex(0.0, 0.0001))
        Gn[sys.net.index[s.bus]] += real(1 / z)
    end
    total_G = sum(Vd[i] * sys.G[i, j] * Vd[j] + Vq[i] * sys.G[i, j] * Vq[j] for i in 1:nb, j in 1:nb)
    norton = sum(Gn .* V2)
    fault = faults_on ? sum(fs.g * V2[fs.index] for fs in sys.faults; init = 0.0) : 0.0
    la = sys.load
    m = sys.coi_members
    coi_omega = length(m) > 1 ?
        sum(sys.coi_weights[q] * x[sys.offsets[m[q]] + 1] for q in eachindex(m)) / sys.coi_total : 1.0
    freq = sum(la.kpf[i] * (coi_omega - 1.0) * la.G[i] * V2[i] for i in 1:nb)
    kcl = sum(Vd[i] * g[2i - 1] + Vq[i] * g[2i] for i in 1:nb if !sys.slack[i])
    machines = Dict{String,Any}()
    loads = Dict{String,Any}()
    for (k, c) in enumerate(sys.comps)
        b, jd, jq = sys.inj[k]
        b > 0 || continue
        on = output_names(c)
        if component_role(c) === :generator && "It_Re" in on
            It = outs[k][findfirst(==("It_Re"), on)], outs[k][findfirst(==("It_Im"), on)]
            s = stamps[name(c)]
            z = complex(s.ra, s.xd_pp)
            pt = Vd[b] * It[1] + Vq[b] * It[2]
            machines[name(c)] = Dict("terminal_out" => pt,
                                     "norton_minus_admittance" => inj_by[name(c)] - real(1 / z) * V2[b],
                                     "Pe_output" => "Pe" in on ? outs[k][findfirst(==("Pe"), on)] : nothing)
        elseif component_role(c) === :load && "Pload" in on
            drawn = la.G[b] * V2[b] - inj_by[name(c)]
            loads[name(c)] = Dict("declared_P" => outs[k][findfirst(==("Pload"), on)], "drawn" => drawn,
                                  "loads_at_bus" => count(q -> sys.inj[q][1] == b && component_role(sys.comps[q]) === :load,
                                                          eachindex(sys.comps)))
        end
    end
    return Dict{String,Any}("injected" => injected,
                            "absorbed" => Dict("total_G" => total_G, "norton" => norton,
                                               "network" => total_G - norton, "fault" => fault,
                                               "frequency_loads" => freq),
                            "kcl_residual_power" => kcl,
                            "balance_error" => injected - (total_G + fault + freq),
                            "machines" => machines, "loads" => loads)
end

"""
    power_audit_samples(sys, x, V; n_random = 8, scale = 1e-2, seed = 11, trajectory = nothing)
        -> Vector{NamedTuple}

The states the audit visits: the equilibrium, `n_random` random perturbations of it
(`scale` relative to each state's magnitude, at least `scale`, reservoirs, held states and
the monitor unchanged, voltages from `solve_network`), and, if given, the states of a
`SimResult` (`trajectory`) with their fault flag (`fault_window`; a record at an event time
is the left limit) and their voltages re-solved on that KCL branch, so every sample
satisfies KCL to round-off.
"""
function power_audit_samples(sys::DAESystem, x::AbstractVector, V::AbstractVector;
                             n_random::Integer = 8, scale::Real = 1e-2, seed::Integer = 11,
                             trajectory = nothing, fault_window = nothing,
                             projection::PhysicalProjection = physical_projection(sys, x, V))
    out = [(label = "equilibrium", x = collect(float(x)), V = collect(float(V)), faults_on = false)]
    rng = Random.Xoshiro(seed)
    for s in 1:n_random
        xs = collect(float(x))
        for i in projection.keep
            xs[i] += scale * max(abs(x[i]), 1.0) * randn(rng)
        end
        Vs = try
            solve_network(sys, xs, V)
        catch
            continue
        end
        push!(out, (label = "random $s", x = xs, V = Vs, faults_on = false))
    end
    if trajectory !== nothing
        nd = sys.n_diff
        for (k, t) in enumerate(trajectory.t)
            # a record at an event time is the left limit (the state just before the switch)
            on = fault_window !== nothing && fault_window[1] < t <= fault_window[2]
            y = trajectory.Y[:, k]
            xs = y[1:nd]
            # the integrator's voltages satisfy KCL only to its tolerance: re-solve them on the
            # healthy or the fault-on branch, from the recorded ones
            Vs = try
                solve_network(sys, xs, y[nd+1:end]; faults_on = on)
            catch
                continue
            end
            push!(out, (label = "t = $(round(t; digits = 4))", x = xs, V = Vs, faults_on = on))
        end
    end
    return out
end

"""
    port_power_audit(sys, samples; contracts = default_contracts()) -> Dict

`component_power` and `network_power` at every sample (from `power_audit_samples`), and a
summary per component and per model type: the range of the residual (`supply - rate`) and
of the storage rate, for machines also the residual against the physical kinetic storage
`H omega^2` (`residual - kinetic_rate_shift`), whether the contract's ports are complete,
and the largest errors of the network, terminal and load identities.
"""
function port_power_audit(sys::DAESystem, samples; contracts::ContractSet = default_contracts())
    per = Dict{String,Any}[]
    comp = Dict{String,Dict{String,Any}}()
    ident = Dict("network_balance" => 0.0, "kcl_residual_power" => 0.0,
                 "machine_terminal" => 0.0, "machine_Pe_vs_terminal" => 0.0,
                 "load_declared_vs_drawn" => 0.0, "field_colocation" => 0.0,
                 "rotor_model" => 0.0, "rotor_identity" => 0.0)
    for s in samples
        cp = component_power(sys, s.x, s.V; contracts, faults_on = s.faults_on)
        np = network_power(sys, s.x, s.V; faults_on = s.faults_on)
        push!(per, Dict{String,Any}("label" => s.label, "faults_on" => s.faults_on,
                                    "components" => cp, "network" => np))
        ident["network_balance"] = max(ident["network_balance"], abs(np["balance_error"] - np["kcl_residual_power"]))
        ident["kcl_residual_power"] = max(ident["kcl_residual_power"], abs(np["kcl_residual_power"]))
        for m in values(np["machines"])
            ident["machine_terminal"] = max(ident["machine_terminal"], abs(m["terminal_out"] - m["norton_minus_admittance"]))
            m["Pe_output"] === nothing ||
                (ident["machine_Pe_vs_terminal"] = max(ident["machine_Pe_vs_terminal"], abs(m["Pe_output"] - m["terminal_out"])))
        end
        for l in values(np["loads"])
            l["loads_at_bus"] == 1 &&
                (ident["load_declared_vs_drawn"] = max(ident["load_declared_vs_drawn"], abs(l["declared_P"] - l["drawn"])))
        end
        for d in cp
            for (key, id) in (("field_colocation_error", "field_colocation"), ("rotor_model_error", "rotor_model"),
                              ("rotor_identity_error", "rotor_identity"))
                haskey(d, key) && (ident[id] = max(ident[id], d[key]))
            end
            e = get!(comp, d["component"]) do
                Dict{String,Any}("type" => d["type"], "complete" => d["complete"],
                                 "ports" => [p["name"] * " = " * p["expression"] for p in d["ports"]],
                                 "residual_min" => Inf, "residual_max" => -Inf,
                                 "rate_min" => Inf, "rate_max" => -Inf, "negative_residual_samples" => 0)
            end
            for key in ("mechanical_residual", "magnetic_residual", "terminal_residual",
                        "stator_exchange", "rotor_loss")
                haskey(d, key) || continue
                e[key * "_min"] = min(get(e, key * "_min", Inf), d[key])
                e[key * "_max"] = max(get(e, key * "_max", -Inf), d[key])
            end
            if haskey(d, "kinetic_rate_shift")
                rp = d["residual"] - d["kinetic_rate_shift"]
                e["residual_physical_kinetic_min"] = min(get(e, "residual_physical_kinetic_min", Inf), rp)
                e["residual_physical_kinetic_max"] = max(get(e, "residual_physical_kinetic_max", -Inf), rp)
                e["negative_physical_kinetic_samples"] = get(e, "negative_physical_kinetic_samples", 0) + (rp < -1e-9)
            end
            e["residual_min"] = min(e["residual_min"], d["residual"])
            e["residual_max"] = max(e["residual_max"], d["residual"])
            e["rate_min"] = min(e["rate_min"], d["storage_rate"])
            e["rate_max"] = max(e["rate_max"], d["storage_rate"])
            d["residual"] < -1e-9 && (e["negative_residual_samples"] += 1)
        end
    end
    # the rotor loss matrix -sym(QA) of each linear rotor, proved positive definite
    # (a rigorous lower bound on its smallest eigenvalue, for the Float64 matrices)
    for c in sys.comps
        rs = rotor_structure(c)
        rs === nothing && continue
        L = -(rs.Q * rs.A + rs.A' * rs.Q) / 2
        comp[name(c)]["rotor_loss_matrix_min_eig_lower_bound"] =
            verified_min_eig(IntervalArithmetic.interval.((L + L') / 2))
    end
    bytype = Dict{String,Any}()
    for (n, e) in comp
        t = get!(bytype, e["type"]) do
            Dict{String,Any}("components" => 0, "complete" => e["complete"], "ports" => e["ports"],
                             "residual_min" => Inf, "residual_max" => -Inf,
                             "components_creating_energy" => String[])
        end
        t["components"] += 1
        t["residual_min"] = min(t["residual_min"], e["residual_min"])
        t["residual_max"] = max(t["residual_max"], e["residual_max"])
        e["negative_residual_samples"] > 0 && push!(t["components_creating_energy"], n)
        if haskey(e, "rotor_loss_matrix_min_eig_lower_bound")
            t["rotor_loss_matrix_min_eig_lower_bound"] = min(get(t, "rotor_loss_matrix_min_eig_lower_bound", Inf),
                                                            e["rotor_loss_matrix_min_eig_lower_bound"])
        end
        if haskey(e, "residual_physical_kinetic_min")
            t["residual_physical_kinetic_min"] = min(get(t, "residual_physical_kinetic_min", Inf), e["residual_physical_kinetic_min"])
            t["residual_physical_kinetic_max"] = max(get(t, "residual_physical_kinetic_max", -Inf), e["residual_physical_kinetic_max"])
            e["negative_physical_kinetic_samples"] > 0 &&
                push!(get!(t, "components_creating_energy_physical_kinetic", String[]), n)
            for key in ("mechanical_residual", "magnetic_residual", "terminal_residual",
                        "stator_exchange", "rotor_loss")
                haskey(e, key * "_min") || continue
                t[key * "_min"] = min(get(t, key * "_min", Inf), e[key * "_min"])
                t[key * "_max"] = max(get(t, key * "_max", -Inf), e[key * "_max"])
            end
        end
    end
    return Dict{String,Any}("samples" => length(samples), "identities_max_error" => ident,
                            "by_type" => bytype, "by_component" => comp, "per_sample" => per)
end
