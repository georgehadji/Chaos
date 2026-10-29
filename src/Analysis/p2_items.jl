# Study p2-indicators (Amendment D89): the items T1 to T17 and the data tables of F1 to F6 (E37). Each item is `(columns, rows)`; a cell is a number, a text or `nothing`
# (undefined). T1 may carry `extra` text lines.

_tbl(cols, rows; extra = String[]) = (columns = String[cols...], rows = Vector{Any}[Any[r...] for r in rows], extra = extra)

_src(o) = o.source === :panel ? "panel" : "grid"

# Groups for the rows by source and by K index: (source label, K label, predicate on an orbit record).
function _groups(orbits)
    ks = sort!(unique(o.kindex for o in orbits if !o.control))
    return vcat([("all", "all", _ -> true)], [(s, "all", o -> _src(o) == s) for s in ("panel", "grid")], [("all", k, o -> o.kindex == k) for k in ks])
end

function p2_t1(orbits, S, reports)
    dropped = Dict(d.id => d.reason for d in S.dropped)
    keys_ = sort!(unique((_src(o), o.kindex) for o in orbits); by = k -> (k[1] == "grid", k[2]))
    rows = Vector{Any}[]
    for (src, k) in keys_
        os = [o for o in orbits if _src(o) == src && o.kindex == k]
        ids = Set(o.id for o in os)
        cnt(v) = count(in(ids), v)
        push!(rows, [src, k, length(os), count(o -> o.rho === :C, os), count(o -> o.rho === :R, os), count(o -> o.rho === :U, os), count(o -> o.control, os),
                     count(o -> get(dropped, o.id, "") in ("drift", "excluded"), os), cnt(S.E_all), cnt(S.E8), cnt(S.W2), cnt(S.S2), cnt(S.overlap), cnt(S.T_s),
                     cnt(S.T_s_design), cnt(S.regular), cnt(S.gate), src == "grid" ? length(os) : 0])
    end
    extra = String["", "Orbits omitted per spec (manifest not_applicable, recomputed from the store): " *
                   (isempty(reports.omissions) ? "none" : join(["$s: $(join(v, ", "))" for (s, v) in sort!(collect(reports.omissions); by = first)], "; ")::String),
                   "Orbits dropped from the sets: " * (isempty(S.dropped) ? "none" : join(["$(d.id) ($(d.reason))" for d in S.dropped], ", ")::String)]
    cols = ["source", "K index", "orbits", "C", "R", "U", "controls", "flagged", "E_all", "E8", "W2", "S2", "W2 and S2", "T_s", "T_s_design", "regular", "gate", "ref_weak"]
    return _tbl(cols, rows; extra = extra)
end

# T2: counts of the verdicts R1 of an indicator at a budget on a set, pooled, by source and by K index, with the step or window used.
function p2_t2(orbits, S, vt, p)
    sets = (("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("regular", S.regular))
    orb = Dict(o.id => o for o in orbits)
    groups = _groups(orbits)
    rows = Vector{Any}[]
    for X in P2_I8, k in 0:(length(p.budgets) - 1)
        anystep = nothing
        for o in orbits
            haskey(vt.steps, (o.id, X)) && (anystep = vt.steps[(o.id, X)][k + 1]; break)
        end
        for (name, ids) in sets, (src, kk, pred) in groups
            sub = [id for id in ids if pred(orb[id])]
            vs = [get(vt.verdicts, (id, X, :R1), nothing) for id in sub]
            v = [x[k + 1] for x in vs if x !== nothing]
            st = unique(vt.steps[(id, X)][k + 1] for id in sub if haskey(vt.steps, (id, X)))
            push!(rows, [String(X), k, p.budgets[k + 1], name, src, kk, count(==(:chaotic), v), count(==(:regular), v), count(==(:undecided), v), count(isnothing, vs),
                         isempty(st) ? anystep : length(st) == 1 ? st[1] : join(st, "/")])
        end
    end
    return _tbl(["indicator", "k", "budget", "set", "source", "K index", "chaotic", "regular", "undecided", "not applicable", "step"], rows)
end

function _verdict_vec(vt, ids, X, k)
    return Symbol[haskey(vt.verdicts, (id, X, :R1)) ? vt.verdicts[(id, X, :R1)][k + 1] : :na for id in ids]
end

# F1: the 2 x 2 counts, the agreement, kappa and the Spearman correlation of every pair of I_8 at b_P and the last budget on W2 and S2.
function p2_f1(S, vt, p)
    rows = Vector{Any}[]
    for (name, ids) in (("W2", S.W2), ("S2", S.S2)), k in (_kp(p) - 1, length(p.budgets) - 1), a in 1:8, b in (a + 1):8
        X, Y = P2_I8[a], P2_I8[b]
        t = p2_contingency(_verdict_vec(vt, ids, X, k), _verdict_vec(vt, ids, Y, k))
        both = [id for id in ids if haskey(vt.scores, (id, X)) && haskey(vt.scores, (id, Y)) && vt.scores[(id, X)][k + 1] !== nothing && vt.scores[(id, Y)][k + 1] !== nothing]
        sp = p2_spearman(Float64[vt.scores[(id, X)][k + 1] for id in both], Float64[vt.scores[(id, Y)][k + 1] for id in both])
        push!(rows, [name, p.budgets[k + 1], String(X), String(Y), t.cc, t.cr, t.rc, t.rr, t.n, t.A, t.kappa, sp])
    end
    return _tbl(["set", "budget", "X", "Y", "cc", "cr", "rc", "rr", "N", "A", "kappa", "spearman"], rows)
end

# F2: the distribution of the detection time t_det, of the censored time and of the outcomes (E28) of S2, W2 and T_s.
function p2_f2(S, vt, p)
    rows = Vector{Any}[]
    vals = vcat([b for b in p.budgets if b * p.persist_divisor <= p.N_c], p.censor_factor * p.N_c)
    for (name, ids) in (("S2", S.S2), ("W2", S.W2), ("T_s", S.T_s)), X in P2_I8
        d = [vt.detection[(id, X)] for id in ids if haskey(vt.detection, (id, X))]
        for v in vals
            push!(rows, [name, String(X), "t_censored", v, count(x -> x.censored == v, d)])
        end
        for o in (:correct, :wrong, :undetected)
            push!(rows, [name, String(X), "outcome", String(o), count(x -> x.outcome === o, d)])
        end
        for b in p.budgets
            push!(rows, [name, String(X), "t_det", b, count(x -> x.t_det == b, d)])
        end
    end
    return _tbl(["set", "indicator", "kind", "value", "count"], rows)
end

# F3: the fraction of chaotic verdicts among the decided ones, per indicator and budget, for W2, S2, T_s and the regular set.
function p2_f3(S, vt, p)
    rows = Vector{Any}[]
    for (name, ids) in (("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("regular", S.regular)), X in P2_I8, k in 0:(length(p.budgets) - 1)
        v = _verdict_vec(vt, ids, X, k)
        dec = count(_decided, v)
        push!(rows, [name, String(X), k, dec, dec == 0 ? nothing : count(==(:chaotic), v) / dec])
    end
    return _tbl(["set", "indicator", "k", "decided", "fraction"], rows)
end

# F4 (and T9): the AUC rows.
function p2_f4(auc)
    rows = [[String(X), k, kidx, src, a.n_C, a.n_R, a.auc] for ((X, k, kidx, src), a) in sort!(collect(auc); by = kv -> (findfirst(==(kv[1][1]), P2_ALL), kv[1][2], kv[1][3], kv[1][4]))]
    return _tbl(["indicator", "k", "K index", "source", "n_C", "n_R", "AUC"], rows)
end

# F5: the lattice of the grid: one row per grid orbit with its class and the verdicts of I_8 at b_P and the last budget.
function p2_f5(grid_rows, vt, p)
    kp, nb = _kp(p), length(p.budgets)
    cols = vcat(["id", "label", "K index", "i", "j", "rho"], [c for X in P2_I8 for c in ("$X b_P", "$X N_c")])
    rows = Vector{Any}[]
    for g in grid_rows
        vs = [get(vt.verdicts, (g.id, X, :R1), nothing) for X in P2_I8]
        push!(rows, vcat(Any[g.id, g.label, g.kindex, g.i, g.j, String(g.rho)], Any[x for v in vs for x in (v === nothing ? ("n.a.", "n.a.") : (String(v[kp]), String(v[nb])))]))
    end
    return _tbl(cols, rows)
end

# F6: the count and the 5th, 25th, 50th, 75th and 95th percentile of the stretching distance D at each recorded step for the classes C and R and for T_s.
function p2_f6(orbits, S, vt, ctx)
    p = ctx.p
    rho = Dict(o.id => o.rho for o in orbits)
    classes = (("C", [id for id in S.E_all if rho[id] === :C]), ("R", S.regular), ("T_s", S.T_s))
    first_it = ctx.idx[p2_key(orbits[1].id, "st", "Float64", p.N_c, 1, p.directions[1])]
    t = first_it.t
    rows = Vector{Any}[]
    for (name, ids) in classes
        its = [ctx.idx[p2_key(id, "st", "Float64", p.N_c, 1, p.directions[1])] for id in ids]
        col = _p2_col(first_it, first_it.columns, "stretching.D")
        for (i, n) in enumerate(t)
            d = sort!(Float64[it.series[i, col] for it in its])
            q = isempty(d) ? fill(nothing, 5) : [_p2_percentile(d, x) for x in (0.05, 0.25, 0.5, 0.75, 0.95)]
            push!(rows, Any[name, n, length(d), q...])
        end
    end
    return _tbl(["class", "step", "count", "p5", "p25", "p50", "p75", "p95"], rows)
end
