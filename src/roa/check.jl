# ROACheck: an independent second check of a certificate record (roadmap P11).
#
# It starts from the record only (the case file, the reference point, the scale and P),
# recomputes every digest and compares it with the record, then re-proves the certified
# level with a different definiteness method: interval Cholesky (the LDL' route) instead of
# the eigenvector-based Gershgorin and Weyl bounds of `certify_level`. If the interval
# Cholesky algorithm runs to the end on an interval matrix with positive pivots, every
# symmetric matrix in it is positive definite (Alefeld and Mayer). The enclosures of the
# field (equilibrium, KCL branch, hull, containment) are the same code as the certificate's:
# the check guards against a wrong record, a changed model or case, and a wrong
# definiteness argument, not against an error in the model code.

"""
    interval_cholesky(A) -> NamedTuple

The interval Cholesky algorithm on the symmetric interval matrix `A` (only the lower
triangle and diagonal are read). `positive_definite = true` proves every symmetric matrix
in `A` positive definite; `min_pivot` is the smallest lower bound of a squared pivot.
"""
function interval_cholesky(A::AbstractMatrix{<:IA.Interval})
    n = size(A, 1)
    L = fill(ival(0.0), n, n)
    minpiv = Inf
    for j in 1:n
        s = A[j, j]
        for k in 1:j-1
            s -= L[j, k]^2
        end
        minpiv = min(minpiv, IA.inf(s))
        IA.inf(s) > 0 || return (positive_definite = false, min_pivot = IA.inf(s), failed_at = j)
        L[j, j] = sqrt(s)
        for i in j+1:n
            t = A[i, j]
            for k in 1:j-1
                t -= L[i, k] * L[j, k]
            end
            L[i, j] = t / L[j, j]
        end
    end
    return (positive_definite = true, min_pivot = minpiv, failed_at = 0)
end

"""
    cholesky_positive_definite(A) -> NamedTuple

Proves every symmetric matrix in the interval matrix `A` positive definite by interval
Cholesky: directly, and, when the widening of the elimination defeats that, after the
congruence `B = X A X'` with `X` the inverse of the floating-point Cholesky factor of
mid(A) (`B` is near the identity; `B > 0` for all members implies `A > 0`, since `X` is a
fixed real matrix and `X A X' > 0` forces `X` nonsingular). `method` says which succeeded.
"""
function cholesky_positive_definite(A::AbstractMatrix{<:IA.Interval})
    d = interval_cholesky(A)
    d.positive_definite && return (positive_definite = true, method = "interval Cholesky", min_pivot = d.min_pivot)
    Am = IA.mid.(A)
    F = cholesky(Symmetric((Am + Am') / 2); check = false)
    issuccess(F) || return (positive_definite = false, method = "midpoint not positive definite in floating point", min_pivot = d.min_pivot)
    X = ival.(inv(F.L))
    B = X * A * permutedims(X)
    p = interval_cholesky(B)
    return (positive_definite = p.positive_definite,
            method = "interval Cholesky of X A X', X = inv(chol(mid A)) (direct Cholesky failed at pivot $(d.failed_at))",
            min_pivot = p.min_pivot)
end

"""
    roa_check(record; contracts = default_contracts()) -> Dict

Rebuild the candidate of a certificate record (from `certify_roa`, with the candidate data
it stores), check the case and contract hashes and the model and candidate fingerprints
against the record, and re-prove its `verified_valid_level`: positivity of `P` and the decay
by interval Cholesky (on `-sym(N'M)` directly, or, when its width defeats it, on the
midpoint shifted by the Collatz bound of the radius), on freshly computed enclosures with
the same hull method. `passed` requires every item.
"""
function roa_check(record::AbstractDict; contracts::ContractSet = default_contracts())
    t0 = time()
    out = Dict{String,Any}("checks" => Dict{String,Any}())
    chk = out["checks"]
    rm = record["model"]
    case = load_case(rm["case_path"])
    chk["case_sha256"] = bytes2hex(open(SHA.sha256, rm["case_path"])) == rm["case_sha256"]
    chk["contracts_sha256"] = bytes2hex(open(SHA.sha256, CONTRACTS_PATH)) == rm["contracts_sha256"]
    x0 = Float64.(rm["x0"])
    V0 = Float64.(rm["V0"])
    sys = assemble(case; init_params = _record_init_params(record))
    chk["system_digest"] = system_digest(sys) == rm["system_digest"]
    m = section_model(sys, x0, V0; scale = Float64.(rm["scale"]), init_params = _record_init_params(record))
    chk["model_fingerprint"] = model_fingerprint(m) == rm["fingerprint"]
    rc = record["candidate"]
    P = reduce(hcat, [Float64.(col) for col in rc["P_columns"]])
    c = QuadraticCandidate(m, P, Dict{String,Any}("method" => "from record"))
    chk["candidate_fingerprint"] = candidate_fingerprint(c) == rc["fingerprint"]
    chk["P_symmetric"] = P == P'
    pc = cholesky_positive_definite(ival.(P))
    chk["P_positive_definite_cholesky"] = pc.positive_definite
    out["P_method"] = pc.method
    level = record["verified_valid_level"]
    if level === nothing
        out["passed"] = false
        out["detail"] = "the record certifies no level"
        return out
    end
    hull = Symbol(get(record, "hull", "centered"))
    eq = enclose_equilibrium(m)
    w = sublevel_half_widths(c, level)
    Xi = [ival(-v, v) for v in w]
    E = eq.eta .+ Xi
    try
        br = enclose_kcl_branch(m, E, IA.mid.(eq.V))
        chk["kcl_branch"] = true
        so = Symbol(get(record, "second_order", "coordinates"))
        M, _ = hull === :centered ? centered_hull(m, eq, E, br.V, Xi; second_order = so, X = br.X) :
               jacobian_hull(m, E, br.V)
        N = gradient_matrix_hull(c, Xi)
        S = permutedims(N) * M
        S = (S .+ permutedims(S)) ./ 2
        direct = cholesky_positive_definite(.-S)
        if direct.positive_definite
            chk["decay_cholesky"] = true
            out["decay_method"] = direct.method * " on -sym(N'M)"
        else
            Sc = IA.mid.(S)
            Sc = (Sc + Sc') / 2
            Rad = [IA.mag(s - ival(v)) for (s, v) in zip(S, Sc)]
            Rad = max.(Rad, Rad')
            rho = _collatz_bound(Rad)
            shifted = .-ival.(Sc)
            for i in axes(shifted, 1)
                shifted[i, i] -= ival(rho)
            end
            sc = cholesky_positive_definite(shifted)
            chk["decay_cholesky"] = sc.positive_definite
            out["decay_method"] = sc.method * " on -mid(S) - rho(Rad) I (on -sym(N'M) itself it failed)"
            out["collatz_rho_upper"] = rho
        end
        audit = containment_audit(m, E, br.V; contracts)
        chk["containment"] = passed(audit)
        chk["box_digest_matches_record"] = record["certified_row"] !== nothing &&
                                           box_digest(E, br.V) == record["certified_row"]["box_digest"]
    catch e
        e isa ProofFailure || e isa UndecidedBranch || e isa IA.InconclusiveBooleanOperation || rethrow()
        out["failure"] = sprint(showerror, e)
        chk["decay_cholesky"] = false
    end
    out["level"] = level
    out["passed"] = all(v -> v === true, values(chk))
    out["time_s"] = time() - t0
    return out
end

"""
    roa_check(path::AbstractString; kwargs...) -> Dict

`roa_check` on a certificate record read from its JSON file.
"""
roa_check(path::AbstractString; kwargs...) =
    roa_check(JSON3.read(read(path, String), Dict{String,Any}); kwargs...)

function _collatz_bound(Rad::AbstractMatrix{Float64})
    v = ones(size(Rad, 1))
    for _ in 1:50
        w = Rad * v
        nw = maximum(w)
        nw > 0 || break
        v = w ./ nw
    end
    v .= max.(v, 1e-3 * maximum(v))
    Rv = ival.(Rad) * ival.(v)
    return maximum(IA.sup(Rv[i] / ival(v[i])) for i in eachindex(v))
end

_record_init_params(record::AbstractDict) =
    Dict{String,Dict{String,Float64}}(string(k) => Dict{String,Float64}(string(p) => Float64(v) for (p, v) in d)
                                      for (k, d) in record["model"]["init_params"])
