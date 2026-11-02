# Study p3-abstention (Amendment D95): the exact statistics of the test side (E66 to E70). Pure functions of tables; counts, shares, p-values, bootstrap statistics, percentiles and
# every comparison with a threshold are exact rationals (I42). The resamples come from the named streams of E66 a (one stream per replicate, no unit depends on another).

const P3T_R = Rational{Int}

"The type 7 percentile of `xs` at the rational probability `p` in exact arithmetic (E66 e), or `nothing` for an empty list."
function p3t_percentile(xs, p)
    isempty(xs) && return nothing
    s = sort(collect(xs))
    pr = p isa Rational ? p : p // 1
    h = (length(s) - 1) * pr + 1
    j = Int(floor(h))
    h == j && return s[j] // 1
    return s[j] + (h - j) * (s[j + 1] - s[j])
end

"The ROC rows `(threshold, fpr, tpr)` of the scores `score` for the failures `f`: one row per distinct score in decreasing order, the rates of the orbits at or above it."
function p3t_roc(score, f)
    nf, nn = count(f), count(!, f)
    return [(threshold = t, fpr = nn == 0 ? nothing : count(i -> !f[i] && score[i] >= t, eachindex(f)) // nn,
             tpr = nf == 0 ? nothing : count(i -> f[i] && score[i] >= t, eachindex(f)) // nf) for t in sort(unique(score); rev = true)]
end

"One replicate: one index vector (1-based inside the stratum) per stratum, in the given order, all drawn from the stream (`purpose`, `unit`); an empty stratum draws nothing (E66 b)."
function p3t_resample(test_seed::Integer, purpose::AbstractString, unit, sizes)
    g = p3_stream(test_seed, purpose, unit)
    return [Int[p3_randint!(g, n) for _ in 1:n] for n in sizes]
end

"The stratum of HA (E64 d, P95): 1 T_s, 2 W2 minus T_s, 3 S2, 4 the other orbits with rho = C, 5 the orbits with rho = R; an orbit is in the first it qualifies for."
p3t_stratum(row) = row.rho === :R ? 5 : row.Ts ? 1 : row.W2 ? 2 : row.S2 ? 3 : 4

"The count rule of the stopping rule (E62 b): `close` iff both counts reach their target or the round is `r_max`."
p3t_decide(t_s, w2, target_t, target_w2, r, r_max) = (r >= r_max || (t_s >= target_t && w2 >= target_w2)) ? "close" : "continue"

# A bound of E66 d, e: the percentile of the valid replicate values, or `nothing` if more than 5 percent are dropped (dropped * 20 > replicates).
function _p3t_bound(vals, nrep, p)
    dropped = nrep - length(vals)
    return (value = dropped * 20 > nrep ? nothing : p3t_percentile(vals, p), dropped = dropped)
end

_frac(a, b) = b == 0 ? nothing : P3T_R(a, b)
_sub(a, b) = (a === nothing || b === nothing) ? nothing : a - b

# ---------------------------------------------------------------- HA (E67)

# The flags of the orbits of a stratum as bits: 1 decided, 2 contradicts rho, 4 unresolved, 8 W2, 16 T_s, 32 S2.
function _p3t_ha_bits(r)
    dec = r.out !== :unresolved
    contra = (r.out === :chaotic_supported && r.rho === :R) || (r.out === :regular_supported && r.rho === :C)
    return UInt8(dec) | UInt8(contra) << 1 | UInt8(!dec) << 2 | UInt8(r.W2) << 3 | UInt8(r.Ts) << 4 | UInt8(r.S2) << 5
end

# The integer counts of a list of bit words: n, decided, contradicting, unresolved, then (n, unresolved) of W2, T_s, S2.
function _p3t_ha_counts(bits)
    c = zeros(Int, 10)
    for b in bits
        c[1] += 1; c[2] += b & 1; c[3] += (b >> 1) & 1; c[4] += (b >> 2) & 1
        for (k, sh) in ((5, 3), (7, 4), (9, 5))
            (b >> sh) & 1 == 1 && (c[k] += 1; c[k + 1] += (b >> 2) & 1)
        end
    end
    return c
end

# The statistics of the counts: e, a(E_test), a(W2), a(T_s), a(S2), dw, dt (`nothing` where undefined).
function _p3t_ha_stats(c)
    aW, aT, aS = _frac(c[6], c[5]), _frac(c[8], c[7]), _frac(c[10], c[9])
    return (e = _frac(c[3], c[2]), aE = _frac(c[4], c[1]), aW = aW, aT = aT, aS = aS, dw = _sub(aW, aS), dt = _sub(aT, aS))
end

# One replicate of HA: the strata are resampled in the fixed order (E66 b) and the statistics taken on the resampled orbits.
function _p3t_ha_replicate(strata, q, r)
    d = p3t_resample(q.test_seed, "boot-ha", r, length.(strata))
    return _p3t_ha_stats(_p3t_ha_counts(UInt8[strata[s][i] for s in eachindex(strata) for i in d[s]]))
end

"""
    p3t_ha(orbits, gate, q; threads) -> NamedTuple

HA (E67) on rows `(id, rho, out, W2, Ts, S2)` and the outputs `gate` of the K = 0 orbits (`nothing` = no verdict): `point`, `bounds` (`e_up`, `aE_up`, `dw_lo`, `dt_lo`, each
`(value, dropped)`), `strata` (sizes), `decision` (`"supported"`, `"not supported"`, `"not decidable"`) and `failed` (the numbers of the failing items 1 to 4).
"""
function p3t_ha(orbits, gate, q; threads::Integer = Threads.nthreads())
    rows = sort(collect(orbits); by = r -> r.id)
    strata = [UInt8[_p3t_ha_bits(r) for r in rows if p3t_stratum(r) == s] for s in 1:5]
    pt = _p3t_ha_stats(_p3t_ha_counts(UInt8[b for s in strata for b in s]))
    reps = _pmap(r -> _p3t_ha_replicate(strata, q, r), 1:q.bootstrap_resamples, threads)
    alpha = p3_rat(q.alpha_ci)
    vals(f) = [getfield(x, f) for x in reps if getfield(x, f) !== nothing]
    nrep = q.bootstrap_resamples
    bounds = (e_up = _p3t_bound(vals(:e), nrep, 1 - alpha), aE_up = _p3t_bound(vals(:aE), nrep, 1 - alpha), dw_lo = _p3t_bound(vals(:dw), nrep, alpha), dt_lo = _p3t_bound(vals(:dt), nrep, alpha))
    g = count(==(:chaotic_supported), gate)
    point = (e = pt.e, aE = pt.aE, a = Dict(:E => pt.aE, :W2 => pt.aW, :Ts => pt.aT, :S2 => pt.aS), dw = pt.dw, dt = pt.dt, g = g)
    sizes = (E = length(rows), W2 = count(r -> r.W2, rows), Ts = count(r -> r.Ts, rows), S2 = count(r -> r.S2, rows))
    return merge((point = point, bounds = bounds, strata = length.(strata)), _p3t_ha_decide(point, bounds, sizes, gate, q))
end

# The decision of E67: not decidable first (P99, P100, P106), then the four items.
function _p3t_ha_decide(pt, b, sizes, gate, q)
    undec = pt.e === nothing || any(<(q.min_set), values(sizes)) || any(x -> x.value === nothing, values(b)) || any(isnothing, gate)
    undec && return (decision = "not decidable", failed = Int[])
    da, ea, aa = p3_rat(q.delta_a), p3_rat(q.e_max), p3_rat(q.a_max)
    ok = [b.e_up.value <= ea, b.aE_up.value <= aa, pt.dw >= da && pt.dt >= da && b.dw_lo.value > 0 && b.dt_lo.value > 0, pt.g <= q.g_max]
    failed = Int[i for i in 1:4 if !ok[i]]
    return (decision = isempty(failed) ? "supported" : "not supported", failed = failed)
end

# ---------------------------------------------------------------- HB (E68)

"""
    p3t_hb_compare(rule, x, q) -> NamedTuple

One paired comparison (E68) of the false-regular indicators of the rule and of a comparator on the same orbits: `n11, n10, n01, n00`, the exact McNemar value `p`, `fr_rule`, `fr_x`,
`margin` (`FR_rule <= rho_hb FR_X`, with `FR_X = 0` failing) and `pass` (`p <= alpha_hb` and the margin).
"""
function p3t_hb_compare(rule, x, q)
    n11, n10 = count(rule .& x), count(rule .& .!x)
    n01, n00 = count(.!rule .& x), count(.!rule .& .!x)
    N = n11 + n10 + n01 + n00
    fr, fx = N == 0 ? 0 // 1 : P3T_R(n11 + n10, N), N == 0 ? 0 // 1 : P3T_R(n11 + n01, N)
    p = p3_mcnemar_p(n10, n01)
    margin = fx != 0 && fr <= p3_rat(q.rho_hb) * fx
    return (n11 = n11, n10 = n10, n01 = n01, n00 = n00, p = p, fr_rule = fr, fr_x = fx, margin = margin, pass = p <= p3_rat(q.alpha_hb) && margin)
end

"""
    p3t_hb(cmp, cov, sizes, q; capped = Symbol[]) -> NamedTuple

HB (E68) from the paired indicators `cmp[(set, X)] = (rule, x)` (`set` in `:Ts`, `:W2`), the coverage `cov[set]` and the final `sizes = (Ts, W2)`: `rows`, `coverage_ok`, `failing`,
`defects` (the capped comparators other than LE), `decision` (floors and defects first, P106), `underpowered`.
"""
function p3t_hb(cmp, cov, sizes, q; capped = Symbol[])
    rows = Dict(k => p3t_hb_compare(v.rule, v.x, q) for (k, v) in cmp)
    cmin = p3_rat(q.c_min)
    cok = cov[:Ts] >= cmin && cov[:W2] >= cmin
    failing = sort!([k for (k, r) in rows if !r.pass]; by = k -> (k[1] === :Ts ? 0 : 1, string(k[2])))
    defects = Symbol[X for X in capped if X !== :LE]
    undec = sizes.Ts < q.n_t_min || sizes.W2 < q.n_w_min || !isempty(defects)
    dec = undec ? "not decidable" : (isempty(failing) && cok) ? "supported" : "not supported"
    return (rows = rows, coverage_ok = cok, failing = failing, defects = defects, decision = dec, underpowered = sizes.Ts < q.target_t_s || sizes.W2 < q.target_w2)
end

# ---------------------------------------------------------------- HC (E70)

# The rows of HC: the plain label from the votes, the failure indicator f, the scores s_I and s_P (exact).
function _p3t_hc_rows(rows)
    return map(sort(collect(rows); by = r -> r.id)) do r
        tie = r.nC == r.nR
        f = tie || (r.nC > r.nR) != (r.rho === :C)
        sI = r.nC + r.nR == 0 ? 0 // 1 : 1 - P3T_R(abs(r.nC - r.nR), r.nC + r.nR)
        (id = r.id, rho = r.rho, nC = r.nC, nR = r.nR, f = f, sI = sI, sP = P3T_R(r.sP))
    end
end

# Dense ranks (1 = smallest) of a score vector, and their number.
function _p3t_ranks(s)
    u = sort(unique(s))
    return [searchsortedfirst(u, x) for x in s], length(u)
end

# The doubled Mann-Whitney count (failure above non-failure twice, tied once) from the rank counts of the failures and the non-failures.
function _p3t_auc2(cf, cn)
    s, below = 0, 0
    for a in eachindex(cf)
        s += cf[a] * (2 * below + cn[a])
        below += cn[a]
    end
    return s
end

# One replicate of HC: the paired resample of the orbits (one stream, n draws) and Delta_AUC, or `nothing` if a class is empty.
function _p3t_hc_replicate(rI, kI, rP, kP, f, q, r)
    n = length(f)
    g = p3_stream(q.test_seed, "boot-hc", r)
    cfI, cnI, cfP, cnP = zeros(Int, kI), zeros(Int, kI), zeros(Int, kP), zeros(Int, kP)
    nf = 0
    for _ in 1:n
        i = p3_randint!(g, n)
        if f[i]
            cfI[rI[i]] += 1; cfP[rP[i]] += 1; nf += 1
        else
            cnI[rI[i]] += 1; cnP[rP[i]] += 1
        end
    end
    (nf == 0 || nf == n) && return nothing
    return P3T_R(_p3t_auc2(cfI, cnI) - _p3t_auc2(cfP, cnP), 2 * nf * (n - nf))
end

"""
    p3t_hc(rows, q; threads) -> NamedTuple

HC (E70) on rows `(id, rho, nC, nR, sP)`: `rows` with `f, sI, sP`, `auc_i`, `auc_p`, `delta`, the bounds `lower` and `upper` of the paired bootstrap, `dropped`, `failures`,
`nonfailures`, `decision` (`"supported"`, `"reversed"`, `"not supported"`, `"not decidable"`) and `underpowered` (failures below `target_f`).
"""
function p3t_hc(rows, q; threads::Integer = Threads.nthreads())
    R = _p3t_hc_rows(rows)
    f = Bool[r.f for r in R]
    sI, sP = [r.sI for r in R], [r.sP for r in R]
    aI, aP = p3_auc(sI[f], sI[.!f]), p3_auc(sP[f], sP[.!f])
    delta = _sub(aI, aP)
    (rI, kI), (rP, kP) = _p3t_ranks(sI), _p3t_ranks(sP)
    vals = Any[x for x in _pmap(r -> _p3t_hc_replicate(rI, kI, rP, kP, f, q, r), 1:q.bootstrap_resamples, threads) if x !== nothing]
    alpha = p3_rat(q.alpha_ci)
    lo, hi = _p3t_bound(vals, q.bootstrap_resamples, alpha), _p3t_bound(vals, q.bootstrap_resamples, 1 - alpha)
    return merge((rows = R, auc_i = aI, auc_p = aP, delta = delta, lower = lo.value, upper = hi.value, dropped = lo.dropped, failures = count(f), nonfailures = count(!, f)),
                 _p3t_hc_decide(delta, lo.value, hi.value, count(f), count(!, f), q))
end

function _p3t_hc_decide(delta, lo, hi, nf, nn, q)
    under = nf < q.target_f
    (nf < q.min_set || nn < q.min_set || lo === nothing || hi === nothing || delta === nothing) && return (decision = "not decidable", underpowered = under)
    dc = p3_rat(q.delta_c)
    dec = (delta >= dc && lo > 0) ? "supported" : (delta <= -dc && hi < 0) ? "reversed" : "not supported"
    return (decision = dec, underpowered = under)
end
