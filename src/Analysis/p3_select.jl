# Study p3-abstention (Amendment D93): calibration C3 to C8 (E54). The verdicts of an orbit at b_P are held as bit masks (one bit per voter) so that the output of R_a on an
# element of Theta is three population counts; all counts are integers and every share and comparison is an exact rational (I32).

# Per-orbit data of C3: for each depth w of `w_range` the voters that are chaotic and stable (`cs`) and regular and stable (`rs`), the voters that are chaotic at all (`ca`).
struct P3Compact
    n::Int
    incal::Vector{Bool}
    rhoR::Vector{Bool}
    S2::Vector{Bool}
    Ts::Vector{Bool}
    W2::Vector{Bool}
    gate::Vector{Bool}
    cs::Vector{Vector{UInt64}}
    rs::Vector{Vector{UInt64}}
    ca::Vector{UInt64}
end

_p3_kP(q) = findfirst(==(q.b_primary), q.budgets) - 1

function _p3_masks(V, voters, k, ws)
    ca = UInt64(0)
    cs = zeros(UInt64, length(ws)); rs = zeros(UInt64, length(ws))
    for (b, X) in enumerate(voters)
        v = V[X]
        bit = UInt64(1) << (b - 1)
        v[k + 1] === :chaotic && (ca |= bit)
        for (iw, w) in enumerate(ws)
            s = p3_stable(v, k, w)
            s && v[k + 1] === :chaotic && (cs[iw] |= bit)
            s && v[k + 1] === :regular && (rs[iw] |= bit)
        end
    end
    return ca, cs, rs
end

# The compact form of orbit rows `(in_cal, rho, S2, Ts, W2, gate, V)`; a row without verdicts is allowed only outside Cal and the gate set.
function _p3_compact(rows, q)
    voters, ws, k = Symbol.(q.voter_bits), q.w_range, _p3_kP(q)
    m = length(rows)
    ca = zeros(UInt64, m); cs = [zeros(UInt64, m) for _ in ws]; rs = [zeros(UInt64, m) for _ in ws]
    for (i, o) in enumerate(rows)
        (o.in_cal || o.gate) || continue
        a, c, r = _p3_masks(o.V, voters, k, ws)
        ca[i] = a
        for iw in eachindex(ws); cs[iw][i] = c[iw]; rs[iw][i] = r[iw]; end
    end
    f(s) = Bool[getproperty(o, s) for o in rows]
    return P3Compact(m, f(:in_cal), Bool[o.rho === :R for o in rows], f(:S2), f(:Ts), f(:W2), f(:gate), cs, rs, ca)
end

_p3_sub(c::P3Compact, idx) = P3Compact(length(idx), c.incal[idx], c.rhoR[idx], c.S2[idx], c.Ts[idx], c.W2[idx], c.gate[idx], [x[idx] for x in c.cs], [x[idx] for x in c.rs], c.ca[idx])

# Integer counts of one element on the orbits of Cal (C3): decided, contradicting, unresolved, unresolved in S2, regular-supported and decided in T_s and W2, and the gate count.
function _p3_count(c::P3Compact, V::UInt64, qC::Int, qR::Int, iw::Int, v::Int)
    dec = contra = unres = unresS2 = frTs = frW2 = decTs = decW2 = g = 0
    cs, rs, ca = c.cs[iw], c.rs[iw], c.ca
    @inbounds for i in 1:c.n
        out = count_ones(cs[i] & V) >= qC ? 1 : (count_ones(ca[i] & V) <= v && count_ones(rs[i] & V) >= qR) ? 2 : 0
        (out == 1 && c.gate[i]) && (g += 1)
        c.incal[i] || continue
        if out == 0
            unres += 1
            c.S2[i] && (unresS2 += 1)
            continue
        end
        dec += 1
        (out == 1) == c.rhoR[i] && (contra += 1)
        c.Ts[i] && (decTs += 1; out == 2 && (frTs += 1))
        c.W2[i] && (decW2 += 1; out == 2 && (frW2 += 1))
    end
    return (dec = dec, contra = contra, unres = unres, unresS2 = unresS2, frTs = frTs, frW2 = frW2, decTs = decTs, decW2 = decW2, g = g)
end

_p3_share(a, n) = n == 0 ? nothing : a // n

# The statistics of one element (E54) in the order of the contract; `sz` holds the sizes of Cal and of S2, T_s, W2 in Cal.
function _p3_stat(n, sz, cost)
    return (e = n.dec == 0 ? nothing : n.contra // n.dec, decided = n.dec, contradicting = n.contra, u_cal = _p3_share(n.unres, sz.cal), u_s2 = _p3_share(n.unresS2, sz.S2),
            fr_ts = _p3_share(n.frTs, sz.Ts), fr_w2 = _p3_share(n.frW2, sz.W2), dec_ts = _p3_share(n.decTs, sz.Ts), dec_w2 = _p3_share(n.decW2, sz.W2), g = n.g, cost = cost)
end

_p3_sizes(c::P3Compact) = (cal = count(c.incal), S2 = count(i -> c.incal[i] && c.S2[i], 1:c.n), Ts = count(i -> c.incal[i] && c.Ts[i], 1:c.n), W2 = count(i -> c.incal[i] && c.W2[i], 1:c.n))

# F1 to F6 of C4, each an exact comparison with the key of the table.
function _p3_flags(s, t, cap)
    le(a, x) = a !== nothing && a <= x
    ge(a, x) = a !== nothing && a >= x
    return (s.cost <= cap, le(s.e, t.e), le(s.u_s2, t.s), s.g == 0, le(s.u_cal, t.a), ge(s.dec_w2, t.c) && ge(s.dec_ts, t.c))
end

# The total order of C6: lowest C_R, smallest |V|, w, v, bitmask of V, q_C, q_R (the tuple identifies the element: no tie, I35).
_p3_key(cost, nV, w, v, mask, qC, qR) = (cost, nV, w, v, mask, qC, qR)

"""
    p3_choose(q, candidates) -> Int

The position of the first candidate `(cost, V, qC, qR, w, v)` in the order of C6 (E54): lowest `C_R`, then smallest `|V|`, `w`, `v`, bitmask of `V`, `q_C`, `q_R`.
"""
function p3_choose(q, candidates)
    bit = Dict(Symbol(X) => i - 1 for (i, X) in enumerate(q.voter_bits))
    key(c) = _p3_key(c.cost, length(c.V), c.w, c.v, sum(1 << bit[X] for X in c.V), c.qC, c.qR)
    return argmin([key(c) for c in candidates])
end

# C3 to C6 on a compact table. `check_sets` applies the floor `min_set` of C8 to the sets of Cal (false inside the optimism bootstrap, P66).
function _p3_run(c::P3Compact, q, th, cost; check_sets::Bool = true)
    sz = _p3_sizes(c)
    iw = Dict(w => i for (i, w) in enumerate(q.w_range))
    cap = p3_cost_cap(q)
    stats = [_p3_stat(_p3_count(c, UInt64(t.mask), t.qC, t.qR, iw[t.w], t.v), sz, cost[i]) for (i, t) in enumerate(th)]
    thr = (e = p3_rat(q.e_cal), s = p3_rat(q.s_cal), a = p3_rat(q.a_cal), c = p3_rat(q.c_min))
    flags = [_p3_flags(s, thr, cap) for s in stats]
    first_fail = Int[(j = findfirst(!, f); j === nothing ? 0 : j) for f in flags]
    F = [i for i in eachindex(th) if first_fail[i] == 0]
    base = (thetas = th, stats = stats, flags = flags, first_fail = first_fail, alone = ntuple(j -> count(f -> f[j], flags), 6), feasible = F, cost_cap = cap, sizes = sz)
    none = (J = Dict{Int, Rational{Int}}(), J_min = nothing, tol = Int[], star = nothing)
    if check_sets && any(<(q.min_set), (sz.Ts, sz.W2, sz.S2))
        return merge(base, (outcome = :sets_too_small,), none)
    elseif isempty(F)
        return merge(base, (outcome = :no_feasible_rule,), none)
    end
    J = Dict{Int, Rational{Int}}(i => (stats[i].fr_ts + stats[i].fr_w2) // 2 for i in F)
    Jmin = minimum(values(J))
    tol = [i for i in F if J[i] <= Jmin + p3_rat(q.delta_tol)]
    star = tol[argmin([_p3_key(stats[i].cost, length(th[i].V), th[i].w, th[i].v, th[i].mask, th[i].qC, th[i].qR) for i in tol])]
    return merge(base, (outcome = :selected, J = J, J_min = Jmin, tol = tol, star = star))
end

_p3_costs(q, th) = [p3_rule_cost(q, t.V, q.b_primary) for t in th]

"""
    p3_select(rows, q) -> NamedTuple

C3 to C8 (E54) on orbit rows `(id, source, kindex, in_cal, rho, S2, Ts, W2, gate, V)`: `outcome` (`:selected`, `:no_feasible_rule`, `:sets_too_small`), `thetas`, `stats[i]`
(`e, decided, contradicting, u_cal, u_s2, fr_ts, fr_w2, dec_ts, dec_w2, g, cost`, exact rationals), `flags[i]` (F1 to F6), `first_fail[i]`, `alone` (the count of each
item alone), `feasible`, `J`, `J_min`, `tol` (Theta_tol), `star` (the index of theta*), `cost_cap` and `sizes`.
"""
p3_select(rows, q) = (th = p3_enumerate(q); _p3_run(_p3_compact(rows, q), q, th, _p3_costs(q, th)))
