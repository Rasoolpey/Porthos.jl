# Structural sparsity pattern of the DAE Jacobian d[f; g]/d[x; V].
#
# Derived from the structure, not detected at a point, so it holds in every limiter mode:
#
#   - a component's outputs depend on its own states and on what its inputs read;
#   - its derivatives depend on the same, plus (through its inputs) the outputs of the
#     components it is wired to, followed transitively;
#   - an input reading a bus voltage depends on (Vd, Vq) of that bus (Vterm on both, the
#     dq-frame voltage also on the machine angle);
#   - KCL at bus i depends on the Y-bus row i, the fault shunt at i, the injections of the
#     components at i, and (frequency-dependent loads) the swing-source speeds;
#   - every rotor angle and delta_COI depend on all swing-source speeds.

"""
    jacobian_pattern(sys) -> SparseMatrixCSC{Bool,Int}

Pattern of `d[f; g] / d[x; V]` (rows: `n_diff` derivatives then `2 nbus` KCL residuals;
columns: `n_diff` states then `[Vd_1, Vq_1, ...]`). Valid in every limiter mode and with
the faults on.
"""
function jacobian_pattern(sys::DAESystem)
    nd = sys.n_diff
    comps = sys.comps
    vcol(i, q) = nd + 2i - 2 + q                       # q = 1: Vd, 2: Vq

    # columns each component's outputs depend on (own states + inputs, transitively)
    own(k) = Set(sys.offsets[k]:(sys.offsets[k] + nstates(comps[k]) - 1))
    deps = [Set{Int}() for _ in comps]
    function input_cols!(acc::Set{Int}, k::Int, visiting::Set{Int})
        for s in sys.sources[k]
            if s.kind in (SRC_VD, SRC_VQ, SRC_VTERM)
                push!(acc, vcol(s.index, 1), vcol(s.index, 2))
            elseif s.kind in (SRC_DQ_VD, SRC_DQ_VQ)
                i = sys.net.index[bus(comps[s.index])]
                push!(acc, vcol(i, 1), vcol(i, 2), sys.offsets[s.index])
            elseif s.kind === SRC_OUTPUT
                j = s.index
                union!(acc, own(j))
                j in visiting || (push!(visiting, j); input_cols!(acc, j, visiting))
            end
        end
        return acc
    end
    for k in eachindex(comps)
        deps[k] = union(own(k), input_cols!(Set{Int}(), k, Set([k])))
    end

    I, J = Int[], Int[]
    add!(r, cols) = for c in cols
        push!(I, r); push!(J, c)
    end
    speeds = [sys.offsets[k] + 1 for k in sys.coi_members]
    for (k, c) in enumerate(comps)
        for r in sys.offsets[k]:(sys.offsets[k] + nstates(c) - 1)
            add!(r, deps[k])
        end
    end
    if length(sys.coi_members) > 1
        for k in sys.coi_members
            add!(sys.offsets[k], speeds)
        end
        add!(sys.delta_coi, speeds)
    end
    nb = nbus(sys)
    for i in 1:nb
        rows = (nd + 2i - 1, nd + 2i)
        if sys.slack[i]
            for r in rows
                add!(r, (vcol(i, 1), vcol(i, 2)))
            end
            continue
        end
        cols = Set{Int}()
        for j in 1:nb
            (sys.G[i, j] != 0 || sys.B[i, j] != 0) && push!(cols, vcol(j, 1), vcol(j, 2))
        end
        push!(cols, vcol(i, 1), vcol(i, 2))                 # fault shunt
        for (k, (b, _, _)) in enumerate(sys.inj)
            b == i && union!(cols, deps[k])
        end
        if sys.load.kpf[i] != 0 || sys.load.kqf[i] != 0
            union!(cols, speeds)
        end
        for r in rows
            add!(r, cols)
        end
    end
    n = nd + nalg(sys)
    return sparse(I, J, trues(length(I)), n, n, |)
end
