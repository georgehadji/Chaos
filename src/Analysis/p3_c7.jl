# Study p3-abstention (Amendment D93): the reported items C7 of E54 for theta* (fitted values, not used for the selection) and the votes of HC (E57).

# The verdict of indicator X of an orbit at its matched budget (LE at the budget it is read at).
_p3_kx(q, m, X) = findfirst(==(X === :LE ? m.le_read : m.X[X].b), q.budgets)

# One row of the HB analogue on Cal (C7): the comparator X against the rule on the set `set` of Cal orbits.
function _p3_hb_row(q, m, th, kP, set, X)
    info = m.X[X]
    bX = X === :LE ? m.le_read : info.b
    kx = _p3_kx(q, m, X)
    n = length(set)
    rfr = Bool[p3_rule_a(o.V, kP, th) === :regular_supported for o in set]
    xfr = Bool[o.V[X][kx] === :regular for o in set]
    n10, n01 = count(rfr .& .!xfr), count(xfr .& .!rfr)
    frx, frr = count(xfr) // n, count(rfr) // n
    return (X = X, b = bX, cost = bX * _p3_cost(q, _p3_spec(q, X)), capped = info.capped, cost_unmatched = X === :LE ? m.le_unmatched : false, fr_x = frx, fr_rule = frr,
            n10 = n10, n01 = n01, p = p3_mcnemar_p(n10, n01), margin = frx != 0 && frr <= p3_rat(q.rho_hb) * frx,
            coverage = count(o -> p3_rule_a(o.V, kP, th) !== :unresolved, set) // n)
end

# The risk-coverage points (C7): the rule family around theta* (q_C and q_R each varied by one step inside their ranges, V, w, v fixed) and each single indicator at its matched budget.
function _p3_risk_coverage(rows, q, sel, th, m)
    ncal = count(o -> o.in_cal, rows)
    fam = [i for (i, t) in enumerate(sel.thetas) if t.V == th.V && t.w == th.w && t.v == th.v && ((t.qC == th.qC && abs(t.qR - th.qR) <= 1) || (t.qR == th.qR && abs(t.qC - th.qC) <= 1))]
    family = [(qC = sel.thetas[i].qC, qR = sel.thetas[i].qR, coverage = sel.stats[i].decided // ncal, e = sel.stats[i].e) for i in fam]
    single = map(Symbol.(q.indicators_i7)) do X
        kx = _p3_kx(q, m, X)
        cal = [o for o in rows if o.in_cal]
        dec = [o for o in cal if o.V[X][kx] !== :undecided]
        contra = count(o -> (o.V[X][kx] === :chaotic && o.rho === :R) || (o.V[X][kx] === :regular && o.rho === :C), dec)
        (X = X, b = X === :LE ? m.le_read : m.X[X].b, coverage = length(dec) // ncal, e = isempty(dec) ? nothing : contra // length(dec))
    end
    return (family = family, single = single)
end

# The rule's false-regular and decided shares on T_s and W2 by source, and e and the unresolved share on the Cal orbits with K > 0 (P48): reported.
function _p3_by_source(rows, q, th, kP)
    out(o) = p3_rule_a(o.V, kP, th)
    share(os, f) = isempty(os) ? nothing : count(f, os) // length(os)
    stats(os) = (n = length(os), fr_ts = share([o for o in os if o.Ts], o -> out(o) === :regular_supported), fr_w2 = share([o for o in os if o.W2], o -> out(o) === :regular_supported),
                 dec_ts = share([o for o in os if o.Ts], o -> out(o) !== :unresolved), dec_w2 = share([o for o in os if o.W2], o -> out(o) !== :unresolved))
    cal = [o for o in rows if o.in_cal]
    pos = [o for o in cal if !o.k0]
    dec = [o for o in pos if out(o) !== :unresolved]
    contra = count(o -> (out(o) === :chaotic_supported) == (o.rho === :R), dec)
    return (panel = stats([o for o in cal if o.source === :panel]), grid = stats([o for o in cal if o.source === :grid]),
            k_positive = (n = length(pos), e = isempty(dec) ? nothing : contra // length(dec), u = share(pos, o -> out(o) === :unresolved)))
end

"""
    p3_c7(rows, q, sel) -> NamedTuple

The items C7 of E54 for theta* (`sel.star`): `top` (the first 20 elements of Theta_tol in the order of C6), `hb[:T_s]` (the six comparators of I_6) and `hb[:W2]` (the seven of I_7)
with `(X, b, cost, capped, cost_unmatched, fr_x, fr_rule, n10, n01, p, margin, coverage)`, the fitted values `e, u_cal, u_s2, g, J, cost`, the `risk_coverage` points and the shares
by source (reported only).
"""
function p3_c7(rows, q, sel)
    th = sel.thetas[sel.star]
    st = sel.stats[sel.star]
    m = p3_matched(q, th.V)
    kP = _p3_kP(q)
    cal = [o for o in rows if o.in_cal]
    key(i) = _p3_key(sel.stats[i].cost, length(sel.thetas[i].V), sel.thetas[i].w, sel.thetas[i].v, sel.thetas[i].mask, sel.thetas[i].qC, sel.thetas[i].qR)
    hb = Dict(:T_s => [_p3_hb_row(q, m, th, kP, [o for o in cal if o.Ts], X) for X in Symbol.(q.voter_bits)],
              :W2 => [_p3_hb_row(q, m, th, kP, [o for o in cal if o.W2], X) for X in Symbol.(q.indicators_i7)])
    return (top = sort(sel.tol; by = key)[1:min(20, length(sel.tol))], hb = hb, e = st.e, u_cal = st.u_cal, u_s2 = st.u_s2, g = st.g, J = sel.J[sel.star], cost = st.cost,
            risk_coverage = _p3_risk_coverage(rows, q, sel, th, m), by_source = _p3_by_source(rows, q, th, kP))
end

"""
    p3_hc_population(rows, q, th) -> NamedTuple

The Cal orbits of HC (E57): `rho` in `{C, R}` and K > 0 (the gate and all K = 0 orbits are out), the decided Float64 votes of V* at b_P (not filtered by stability, P64), their plain
label, the failure indicator `f` (a tie is a failure) and the score `s`. Returns `fails` and `nonfails` (the scores), their counts and `auc` (`Float64`, `nothing` if a class is empty).
"""
function p3_hc_population(rows, q, th)
    kP = _p3_kP(q)
    fails, nonf = Float64[], Float64[]
    for o in rows
        (o.in_cal && !o.k0) || continue
        s = p3_vote_summary([o.V[X][kP + 1] for X in th.V if o.V[X][kP + 1] !== :undecided], o.rho)
        push!(s.f ? fails : nonf, s.s)
    end
    a = p3_auc(fails, nonf)
    return (fails = fails, nonfails = nonf, n_fail = length(fails), n_nonfail = length(nonf), auc = a === nothing ? nothing : Float64(a))
end
