# Study p3-abstention (Amendment D95): the data of the figures F2 and F3 (E72). No figure is drawn (P110); each block carries the SHA-256 of its table written with 12 significant digits.

# The points of the rule family around theta*: q_C and q_R each moved by one step inside their ranges (V*, w, v fixed), the centre first.
function _p3t_family(theta)
    n = length(theta.V)
    pts = [(theta.qC, theta.qR)]
    for d in (-1, 1)
        1 <= theta.qC + d <= n && push!(pts, (theta.qC + d, theta.qR))
    end
    for d in (-1, 1)
        (n ÷ 2 + 1) <= theta.qR + d <= n && push!(pts, (theta.qC, theta.qR + d))
    end
    return pts
end

# The decided share and the error of a list of outputs against rho: chaotic for rho = R and regular for rho = C contradict; NaN where nothing is decided.
function _p3t_risk(outs, rows, chaotic, regular)
    dec = count(o -> o === chaotic || o === regular, outs)
    contra = count(i -> (outs[i] === chaotic && rows[i].rho === :R) || (outs[i] === regular && rows[i].rho === :C), eachindex(outs))
    n = length(outs)
    return (n == 0 ? NaN : Float64(P3T_R(dec, n)), dec == 0 ? NaN : Float64(P3T_R(contra, dec)))
end

"""
    p3t_f2(rows, qt, theta, M) -> NamedTuple

The data of F2 (E72): the points of error against decided share at b_P of the rule family around theta* on E_test, then each single indicator of I_7 at its matched budget (a chaotic verdict
as chaotic-supported, a regular one as regular-supported). `(columns, rows, digest)`.
"""
function p3t_f2(rows, qt, theta, M)
    kP = _p3_kP(qt)
    out = Vector{Any}[]
    for (qC, qR) in _p3t_family(theta)
        th = merge(theta, (qC = qC, qR = qR))
        s, e = _p3t_risk([p3_rule_a(r.V, kP, th) for r in rows], rows, :chaotic_supported, :regular_supported)
        push!(out, Any["rule", "q_C=$qC q_R=$qR", s, e])
    end
    for X in Symbol.(qt.indicators_i7)
        b = X === :LE ? M.le_read : M.X[X].b
        kx = b === nothing ? nothing : findfirst(==(b), qt.budgets)
        s, e = kx === nothing ? (NaN, NaN) : _p3t_risk([r.V[X][kx] for r in rows], rows, :chaotic, :regular)
        push!(out, Any["single", string(X), s, e])
    end
    t = _tbl(["curve", "point", "decided_share", "error"], out)
    return (columns = t.columns, rows = t.rows, digest = p2_digest(t))
end

"""
    p3t_f3(hc) -> NamedTuple

The data of F3 (E72): the ROC points of the scores s_I and s_P for the failures f on the HC set, one row per distinct score in decreasing order `(curve, threshold, fpr, tpr)`, and the AUC
row of each curve (threshold `AUC`, the AUC in the third column). `(columns, rows, digest)`.
"""
function p3t_f3(hc)
    f = Bool[r.f for r in hc.rows]
    out = Vector{Any}[]
    for (curve, sc, auc) in (("s_I", [r.sI for r in hc.rows], hc.auc_i), ("s_P", [r.sP for r in hc.rows], hc.auc_p))
        for r in p3t_roc(sc, f)
            push!(out, Any[curve, Float64(r.threshold), r.fpr === nothing ? NaN : Float64(r.fpr), r.tpr === nothing ? NaN : Float64(r.tpr)])
        end
        push!(out, Any[curve, "AUC", auc === nothing ? NaN : Float64(auc), "n.a."])
    end
    t = _tbl(["curve", "threshold", "fpr", "tpr"], out)
    return (columns = t.columns, rows = t.rows, digest = p2_digest(t))
end
