# Study p3-abstention (Amendment D93): PP for HB, the part that uses the Cal data only (E56 a to c): the optimism bootstrap, the corrected rates, the strata, the model of one orbit.

"""
    p3_optimism(rows, q; threads) -> (optimism, used, skipped, thetas)

The optimism of E56 a. For `b = 1` to `b_opt` the stream (optimism, b) draws as many orbits as Cal has, with replacement, from Cal listed by increasing id; C3 to C6 run on the
resample (the floor `min_set` is not applied, P66) and give theta_b, or `nothing` if F is empty (the resample is then skipped and counted). `optimism = (T_s, W2)` is the mean
over the valid resamples, summed in increasing `b`, of `FR(theta_b; Cal set) - FR(theta_b; resample set)` (binary64); `thetas[b]` is the index of theta_b.
"""
function p3_optimism(rows, q; threads::Integer = Threads.nthreads())
    th = p3_enumerate(q)
    cost = _p3_costs(q, th)
    c = _p3_compact(rows, q)
    cal = [i for i in sortperm([o.id for o in rows]) if c.incal[i]]
    iw = Dict(w => i for (i, w) in enumerate(q.w_range))
    function unit(b)
        g = p3_stream(q.test_seed, "optimism", b)
        idx = [cal[p3_randint!(g, length(cal))] for _ in cal]
        r = _p3_run(_p3_sub(c, idx), q, th, cost; check_sets = false)
        return r.star === nothing ? nothing : (star = r.star, fr_ts = r.stats[r.star].fr_ts, fr_w2 = r.stats[r.star].fr_w2)
    end
    res = _pmap(unit, 1:q.b_opt, threads)
    sz = _p3_sizes(c)
    full(i) = (n = _p3_count(c, UInt64(th[i].mask), th[i].qC, th[i].qR, iw[th[i].w], th[i].v); (_p3_share(n.frTs, sz.Ts), _p3_share(n.frW2, sz.W2)))
    fulls = Dict(i => full(i) for i in unique(r.star for r in res if r !== nothing))
    sT = sW = 0.0
    for r in res
        r === nothing && continue
        sT += Float64(fulls[r.star][1]) - Float64(r.fr_ts)
        sW += Float64(fulls[r.star][2]) - Float64(r.fr_w2)
    end
    used = count(!isnothing, res)
    return (optimism = used == 0 ? nothing : (T_s = sT / used, W2 = sW / used), used = used, skipped = q.b_opt - used,
            thetas = Union{Int, Nothing}[r === nothing ? nothing : r.star for r in res])
end

# Counts of the rule and of the comparators on the sets T_s and W2 of Cal at theta (E56 b): decided, false-regular, and for each comparator its false-regular count and the
# count shared with the rule.
function _p3_hb_counts(rows, q, theta)
    kP = _p3_kP(q)
    m = p3_matched(q, theta.V)
    cal = [o for o in rows if o.in_cal]
    comps = Symbol.(q.indicators_i7)
    kx(X) = findfirst(==(X === :LE ? m.le_read : m.X[X].b), q.budgets)
    cnt(set) = begin
        os = [o for o in cal if set(o)]
        out = [p3_rule_a(o.V, kP, theta) for o in os]
        rfr = out .=== :regular_supported
        xfr(X) = Bool[o.V[X][kx(X)] === :regular for o in os]
        (n = length(os), dec = count(!=(:unresolved), out), rfr = count(rfr), x = Dict(X => count(xfr(X)) for X in comps), shared = Dict(X => count(rfr .& xfr(X)) for X in comps))
    end
    return (T = cnt(o -> o.Ts), W = cnt(o -> o.W2))
end

# The rule rate of a stratum: the corrected rate if it is at least the floor, else the larger of factor times it and the floor.
_p3_rule_rate(r, q) = r >= q.rule_rate_floor ? r : max(q.rule_rate_factor * r, q.rule_rate_floor)

# R_W by counts by difference (the share of W2 minus that of T_s, over the size of R_W); an empty R_W gives `default`.
_p3_diff(aW, aT, nW, nT, default = 0.0) = nW - nT == 0 ? default : Float64((aW - aT) // (nW - nT))

"""
    p3_hb_rates(rows, q, theta, optimism) -> NamedTuple

The rates of E56 a and b for the element `theta` (`V, qC, qR, w, v`) from the Cal rows and the optimism `(T_s, W2)` (a number for each; `nothing` counts as zero): `r` (rule
false-regular share on T_s and W2), `r_corr` (also R_W, by difference, clamped), `p_rule` (T_s, R_W, the floor rule), `p_comp[(X, S)]` and `s_comp[(X, S)]` for each
indicator of I_7 and the strata `:T_s`, `:R_W` (comparator rate at its matched budget, shared share or `shared_default`), `d` (decided share), and the comparator lists.
"""
function p3_hb_rates(rows, q, theta, optimism)
    k = _p3_hb_counts(rows, q, theta)
    T, W = k.T, k.W
    r = Dict(:T_s => Float64(T.rfr // T.n), :W2 => Float64(W.rfr // W.n))
    opt(a) = optimism === nothing ? 0.0 : max(getfield(optimism, a), 0)
    rc = Dict(:T_s => r[:T_s] + opt(:T_s), :W2 => r[:W2] + opt(:W2))
    rc[:R_W] = W.n - T.n == 0 ? 0.0 : clamp((rc[:W2] * W.n - rc[:T_s] * T.n) / (W.n - T.n), 0.0, 1.0)
    pc, sc = Dict{Tuple{Symbol, Symbol}, Float64}(), Dict{Tuple{Symbol, Symbol}, Float64}()
    for X in Symbol.(q.indicators_i7)
        pc[(X, :T_s)] = Float64(T.x[X] // T.n); pc[(X, :R_W)] = _p3_diff(W.x[X], T.x[X], W.n, T.n)
        sc[(X, :T_s)] = T.rfr == 0 ? q.shared_default : Float64(T.shared[X] // T.rfr)
        sc[(X, :R_W)] = _p3_diff(W.shared[X], T.shared[X], W.rfr, T.rfr, q.shared_default)
    end
    return (r = r, r_corr = rc, p_rule = Dict(:T_s => _p3_rule_rate(rc[:T_s], q), :R_W => _p3_rule_rate(rc[:R_W], q)), p_comp = pc, s_comp = sc,
            d = Dict(:T_s => Float64(T.dec // T.n), :R_W => _p3_diff(W.dec, T.dec, W.n, T.n)), ts = Symbol.(q.voter_bits), w2 = Symbol.(q.indicators_i7))
end

"The model of one orbit (E56 c) from the rates: the rule rate and decided share per stratum, the comparator rate and shared share per stratum (zero for a comparator outside I_7 of the table)."
function p3_hb_model(h)
    z = (p = 0.0, s = 0.0)
    cm(X) = haskey(h.p_comp, (X, :T_s)) ? (T = (p = h.p_comp[(X, :T_s)], s = h.s_comp[(X, :T_s)]), RW = (p = h.p_comp[(X, :R_W)], s = h.s_comp[(X, :R_W)])) : (T = z, RW = z)
    return (rule = Dict(:T => (p = h.p_rule[:T_s], d = h.d[:T_s]), :RW => (p = h.p_rule[:R_W], d = h.d[:R_W])), comp = Dict(X => cm(X) for X in P3_COMPS), ts = h.ts, w2 = h.w2)
end
