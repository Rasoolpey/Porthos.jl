# Single-mode containment audit (roadmap P11; TODO.md decision "limiters: single mode
# first"). On the same boxes as the decay proof (section coordinates E, bus voltages V on
# the certified KCL branch) it proves, fail closed:
#
#   evaluation   the field and KCL residual are well defined on the box (every interval
#                decorated at least `dac`: no division by an interval containing zero, no
#                square root of a negative part), and bounded;
#   limiter_mode every branch site of every component is decided on the box (one smooth
#                mode; the decisions are recorded), so the field there is one smooth formula;
#   held         the identically held states have a derivative that is exactly zero on it;
#   nonbinding   the retained rows and KCL do not depend on the one-way reservoir states or
#                the delta_COI monitor: their interval partial derivatives are exactly zero,
#                with each such state anywhere in a recorded range (the widest of
#                x0 * [2^-k, 2^k] that passes); the drift rate of each over the box is
#                bounded too. Global independence (for every value) is not claimed;
#   domain       every contract `domain_clause` of every component, evaluated with interval
#                bounds (state bounds, parameter conditions, bus voltage above a threshold;
#                the single-mode and nonbinding clauses point to the categories above).
#
# A component type without a contract entry, or a clause check this audit does not know,
# fails the audit.

"""
    ContainmentAudit

The result of `containment_audit`: per-category verdicts, the clause records, the branch
decisions of every component on the box, the ranges over which the excluded states are
proved not to feed back (with a bound on their drift rate while the retained state is in the
box), and the SHA-256 of the boxes it ran on.
"""
struct ContainmentAudit
    passed::Dict{String,Bool}
    records::Vector{Dict{String,Any}}
    modes::Dict{String,Any}
    excluded::Dict{String,Any}           # per excluded state: range, x0, drift-rate bound
    box_digest::String
end

passed(a::ContainmentAudit) = all(values(a.passed))

"""
    box_digest(E, V) -> String

SHA-256 of the interval boxes (their bounds), which ties the gates of one certificate to
the same boxes.
"""
box_digest(E::AbstractVector{<:IA.Interval}, V::AbstractVector{<:IA.Interval}) =
    _digest("E_lo" => IA.inf.(E), "E_hi" => IA.sup.(E), "V_lo" => IA.inf.(V), "V_hi" => IA.sup.(V))

_component_env(c::AbstractComponent) =
    Dict{String,Float64}(k => v for (k, v) in
                         ((k, try param_value(v) catch; nothing end) for (k, v) in param_dict(c))
                         if v !== nothing)

const _CONDITION_RE = r"^(.*?)(>=|<=|!=|==|>|<)(.*)$"

# A parameter condition "lhs op rhs", both sides in interval arithmetic; decided or false.
function _condition(s::AbstractString, env)
    mt = match(_CONDITION_RE, s)
    mt === nothing && return (passed = false, detail = "cannot read the condition")
    a = eval_param_expr(Ival, strip(mt.captures[1]), env)
    b = eval_param_expr(Ival, strip(mt.captures[3]), env)
    op = mt.captures[2]
    ok = op == ">" ? IA.inf(a) > IA.sup(b) :
         op == ">=" ? IA.inf(a) >= IA.sup(b) :
         op == "<" ? IA.sup(a) < IA.inf(b) :
         op == "<=" ? IA.sup(a) <= IA.inf(b) :
         op == "!=" ? (IA.sup(a) < IA.inf(b) || IA.sup(b) < IA.inf(a)) :
         IA.isthin(a) && IA.isthin(b) && IA.inf(a) == IA.inf(b)
    return (passed = ok, detail = "$(strip(mt.captures[1])) = $a $op $(strip(mt.captures[3])) = $b")
end

"""
    containment_audit(m, E, V; contracts = default_contracts(), spans = (60, 40, 20, 10))
        -> ContainmentAudit

The audit on the section box `E` and the voltage box `V` (see the comment at the top of
`src/roa/containment.jl`). Each reservoir and monitor state gets the widest range
`x0 * [2^-k, 2^k]` (`[-2^k, 2^k]` when `x0 = 0`), `k` in `spans`, over which the retained
rows and KCL are proved independent of it.
"""
function containment_audit(m::SectionModel, E::AbstractVector{<:IA.Interval},
                           V::AbstractVector{<:IA.Interval};
                           contracts::ContractSet = default_contracts(),
                           spans = (60, 40, 20, 10))
    sys = m.sys
    p = m.projection
    names = state_names(sys)
    records = Dict{String,Any}[]
    verdict = Dict("evaluation" => true, "limiter_mode" => true, "held" => true,
                   "nonbinding" => true, "domain" => true)
    add!(cat, id, comp, ok, detail) = (push!(records, Dict{String,Any}(
        "category" => cat, "id" => id, "component" => comp, "passed" => ok, "detail" => detail));
        ok || (verdict[cat] = false))

    # evaluation and limiter modes
    x = lift(m, E)
    logs = (out = [ModeLog() for _ in sys.comps], step = [ModeLog() for _ in sys.comps],
            current = Ref((:none, 0)))
    f = g = nothing
    undecided = nothing
    try
        f, g = state_field(m, x, V; logs)
    catch e
        e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
        kernel, k = logs.current[]
        undecided = (component = k > 0 ? name(sys.comps[k]) : "?", kernel = kernel,
                     message = sprint(showerror, e))
    end
    modes = Dict{String,Any}()
    if undecided === nothing
        add!("limiter_mode", "all_branches_decided", "", true,
             "$(sum(length(l.decisions) for l in logs.out) + sum(length(l.decisions) for l in logs.step)) branch sites decided on the box")
        for (k, c) in enumerate(sys.comps)
            d = vcat(logs.out[k].decisions, logs.step[k].decisions)
            isempty(d) || (modes[name(c)] = Dict("out" => logs.out[k].decisions,
                                                 "step" => logs.step[k].decisions))
        end
        ok = well_defined(f) && well_defined(g)
        add!("evaluation", "field_well_defined", "", ok,
             ok ? "f and g decorated >= dac and bounded on the box" :
                  "an entry of f or g is not well defined on the box")
        held_ok = all(i -> IA.inf(f[i]) == 0 && IA.sup(f[i]) == 0, p.held)
        add!("held", "held_states_constant", "", held_ok,
             "$(length(p.held)) held states: " * (held_ok ? "derivative exactly 0 on the box" :
                  "nonzero derivative at " * join(names[filter(i -> !(IA.inf(f[i]) == 0 && IA.sup(f[i]) == 0), p.held)], ", ")))
    else
        add!("limiter_mode", "all_branches_decided", undecided.component, false,
             "undecided branch in the $(undecided.kernel) kernel: $(undecided.message)")
        add!("evaluation", "field_well_defined", "", false, "not evaluated (undecided branch)")
    end

    # nonbinding: reservoirs and the monitor do not feed the retained rows or KCL. For each
    # such state, the widest range x0 * [2^-k, 2^k] (k from `spans`, largest first) over which
    # the partials of the retained rows and KCL are exactly zero; then all ranges jointly, and
    # the rate at which each state can drift while the retained state is in the box.
    free = vcat(p.reservoir, p.monitor)
    nb = Dict{Int,Bool}(i => false for i in free)
    excluded = Dict{String,Any}()
    if undecided === nothing
        function zero_partials(idxs, R)
            try
                J = ForwardDiff.jacobian(R) do r
                    xr = Vector{eltype(r)}(x)
                    xr[idxs] .= r
                    fr, gr = state_field(m, xr, V)
                    vcat(fr[m.keep], gr)
                end
                return [all(v -> IA.inf(v) == 0 && IA.sup(v) == 0, view(J, :, j)) for j in eachindex(idxs)]
            catch e
                e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
                return fill(false, length(idxs))
            end
        end
        range_for(i, k) = m.x0[i] == 0 ? ival(-exp2(k), exp2(k)) :
            ival(min(m.x0[i] * exp2(-k), m.x0[i] * exp2(k)), max(m.x0[i] * exp2(-k), m.x0[i] * exp2(k)))
        ranges = Dict{Int,Ival}()
        for i in free
            for k in spans
                if only(zero_partials([i], [range_for(i, k)]))
                    ranges[i] = range_for(i, k)
                    break
                end
            end
        end
        if length(ranges) == length(free)
            joint = zero_partials(free, [ranges[i] for i in free])
            for (j, i) in enumerate(free)
                nb[i] = joint[j]
            end
            rates = try
                xr = copy(x)
                xr[free] .= [ranges[i] for i in free]
                fr, _ = state_field(m, xr, V)
                [IA.mag(fr[i]) for i in free]
            catch e
                e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
                fill(Inf, length(free))
            end
            for (j, i) in enumerate(free)
                excluded[names[i]] = Dict("lower" => IA.inf(ranges[i]), "upper" => IA.sup(ranges[i]),
                                          "x0" => m.x0[i], "drift_rate_upper" => rates[j])
            end
        end
        bad = [names[i] for i in free if !nb[i]]
        add!("nonbinding", "reservoirs_and_monitor_do_not_feed_back", "", isempty(bad),
             isempty(bad) ? "$(length(free)) excluded states: retained rows and KCL have exactly zero partials in each, jointly over its recorded range (excluded_ranges)" :
                            "no feedback-free range found for " * join(bad, ", "))
    else
        add!("nonbinding", "reservoirs_and_monitor_do_not_feed_back", "", false, "not evaluated")
    end

    # contract domain clauses
    for (k, c) in enumerate(sys.comps)
        entry = get(contracts.entries, contract_key(model_type(c)), nothing)
        if entry === nothing
            add!("domain", "contract", name(c), false, "no contract entry for $(model_type(c))")
            continue
        end
        env = _component_env(c)
        for cl in entry.domain_clauses
            if cl.check === :state_bounds
                j = findfirst(==(cl.state), state_names(c))
                if j === nothing
                    add!("domain", cl.id, name(c), false, "no state $(cl.state)")
                    continue
                end
                xi = x[sys.offsets[k] + j - 1]
                lo = cl.lower === nothing ? nothing : eval_param_expr(Ival, cl.lower, env)
                hi = cl.upper === nothing ? nothing : eval_param_expr(Ival, cl.upper, env)
                ok = (lo === nothing || IA.inf(xi) > IA.sup(lo)) && (hi === nothing || IA.sup(xi) < IA.inf(hi))
                add!("domain", cl.id, name(c), ok,
                     "$(cl.state) in $xi; bounds $(something(lo, "-")) .. $(something(hi, "-"))")
            elseif cl.check === :parameter_conditions
                res = [_condition(s, env) for s in cl.conditions]
                add!("domain", cl.id, name(c), all(r -> r.passed, res), join((r.detail for r in res), "; "))
            elseif cl.check === :bus_voltage_squared_above
                b = bus(c)
                i = b === nothing ? nothing : get(sys.net.index, b, nothing)
                if i === nothing
                    add!("domain", cl.id, name(c), false, "component has no bus")
                    continue
                end
                v2 = V[2i - 1]^2 + V[2i]^2
                add!("domain", cl.id, name(c), IA.inf(v2) > cl.threshold,
                     "Vd^2 + Vq^2 in $v2 > $(cl.threshold)")
            elseif cl.check === :nonbinding_accounting_state
                j = findfirst(==(cl.state), state_names(c))
                i = j === nothing ? 0 : sys.offsets[k] + j - 1
                ok = i > 0 && i in p.reservoir && get(nb, i, false)
                add!("domain", cl.id, name(c), ok,
                     ok ? "$(cl.state) excluded from the section and feeds no retained row" :
                          "$(cl.state) is not an excluded, nonbinding reservoir")
            elseif cl.check === :single_smooth_mode
                add!("domain", cl.id, name(c), undecided === nothing,
                     undecided === nothing ? "every branch decided on the box" : "undecided branch")
            else
                add!("domain", cl.id, name(c), false, "unknown clause check $(cl.check)")
            end
        end
    end
    return ContainmentAudit(verdict, records, modes, excluded, box_digest(E, V))
end
