# Study p3-abstention (Amendment D93): PP for HB (E56 d to f): the draw path, the simulation of one point, the targets, the type-I check. Every run owns a named stream
# (pp-hb, (point, run)); the counts are integers, so the thread count never changes a number (I31).

const P3_COMPS = P2_I7

"""
    p3_draw_path(q) -> Vector{NamedTuple}

The points of the draw path of E56 d: `(j, rho, nT, nRW, nW, allowed)` for `rho = rho_path_start + j step` up to `r_max` included (the step is the difference of
`rho_path_second` and `rho_path_start`); sizes by round half away from zero; a point is allowed iff `nT >= n_t_min` and `nW >= n_w_min`.
"""
function p3_draw_path(q)
    start, step, rmax = p3_rat(q.rho_path_start), p3_rat(q.rho_path_second) - p3_rat(q.rho_path_start), p3_rat(q.r_max)
    out = NamedTuple[]
    j = 0
    while start + j * step <= rmax
        rho = start + j * step
        nT, nRW = p3_round_half_away(rho * q.y_t), p3_round_half_away(rho * (q.y_w - q.y_t))
        push!(out, (j = j, rho = rho, nT = nT, nRW = nRW, nW = nT + nRW, allowed = nT >= q.n_t_min && nT + nRW >= q.n_w_min))
        j += 1
    end
    return out
end

# The model of one stratum prepared for the nine uniforms of an orbit (E56 c): the rule rate p, the probability of being decided given that the rule is not false-regular, and for
# each comparator the probability of false-regular given that the rule is (s) and given that it is not.
function _p3_stratum(model, st::Symbol)
    R = model.rule[st]
    p = Float64(R.p)
    dp = p >= 1 ? 1.0 : clamp((R.d - p) / (1 - p), 0.0, 1.0)
    qa = Float64[model.comp[X][st].s for X in P3_COMPS]
    qb = Float64[p >= 1 ? 0.0 : max(model.comp[X][st].p - model.comp[X][st].s * p, 0.0) / (1 - p) for X in P3_COMPS]
    return (p = p, dp = dp, qa = qa, qb = qb)
end

# One simulated run: the counts n11, n10, n01, n00 (rule false-regular and X; rule only; X only; neither) per stratum and comparator, and the decided counts.
function _p3_orbits!(cnt, dec, g, S, nT, nRW)
    fill!(cnt, 0); fill!(dec, 0)
    for s in 1:2, _ in 1:(s == 1 ? nT : nRW)
        st = S[s]
        rfr = p3_uniform!(g) < st.p
        for x in 1:length(P3_COMPS)
            fr = p3_uniform!(g) < (rfr ? st.qa[x] : st.qb[x])
            cnt[s, x, rfr ? (fr ? 1 : 2) : (fr ? 3 : 4)] += 1
        end
        (p3_uniform!(g) < st.dp || rfr) && (dec[s] += 1)
    end
end

# The comparison of one set: m >= 1, the exact tail at most alpha (critical value), c >= 1 and a <= rho_hb c.
function _p3_cmp_ok(n11, n10, n01, crit, rho)
    m = n10 + n01
    (m >= 1 && n01 >= crit[m + 1] && n01 + n11 >= 1) || return false
    return (n10 + n11) * denominator(rho) <= numerator(rho) * (n01 + n11)
end

# The critical values c(m), m = 0 to nmax, of the exact McNemar tail at level alpha (c(m) = m + 1 if none): c(m + 1) is c(m) or c(m) + 1 (the tail is the mean of two tails of m), so the table
# follows from the tail count at c(m) by two recurrences, in O(nmax) exact integer operations; it equals `p3_mcnemar_critical` for every m.
function _p3_crit_table(alpha::Rational, nmax::Integer)
    crit = Vector{Int}(undef, nmax + 1)
    crit[1] = 1
    T, B, c = big(0), big(1), 1          # T = sum of C(m, j), j >= c;  B = C(m, c - 1)
    for m in 0:(nmax - 1)
        T = 2T + B                                  # the tail count of m + 1 at c
        B = B * (m + 1) ÷ (m + 2 - c)               # C(m + 1, c - 1)
        if T > alpha * big(2)^(m + 1)
            Cc = B * (m + 2 - c) ÷ c                # C(m + 1, c)
            T -= Cc; B = Cc; c += 1
        end
        crit[m + 2] = c
    end
    return crit
end

# The prepared decision data: positions of the comparators of T_s and of W2, the critical table, the margin and the coverage threshold.
function _p3_decision_data(model, q, nmax)
    pos(L) = Int[findfirst(==(X), P3_COMPS) for X in L]
    return (ts = pos(model.ts), w2 = pos(model.w2), crit = _p3_crit_table(p3_rat(q.alpha_hb), nmax), rho = p3_rat(q.rho_hb), cmin = p3_rat(q.c_min))
end

# Success of the full decision of E56 e, and the single comparison of each comparator (all its comparisons without coverage and the others).
function _p3_decide(cnt, dec, D, nT, nRW, singles)
    okT(x) = _p3_cmp_ok(cnt[1, x, 1], cnt[1, x, 2], cnt[1, x, 3], D.crit, D.rho)
    okW(x) = _p3_cmp_ok(cnt[1, x, 1] + cnt[2, x, 1], cnt[1, x, 2] + cnt[2, x, 2], cnt[1, x, 3] + cnt[2, x, 3], D.crit, D.rho)
    cov = dec[1] * denominator(D.cmin) >= numerator(D.cmin) * nT && (dec[1] + dec[2]) * denominator(D.cmin) >= numerator(D.cmin) * (nT + nRW)
    good = all(okT, D.ts) && all(okW, D.w2)
    if singles !== nothing
        for x in eachindex(P3_COMPS)
            singles[x] += (x in D.ts ? okT(x) : true) && (x in D.w2 ? okW(x) : true)
        end
    end
    return good && cov
end

# The runs i0:i1 of the point: the successes (and the singles if asked).
function _p3_hb_chunk(S, D, q, nT, nRW, point, runs_range, want_singles)
    cnt = zeros(Int, 2, length(P3_COMPS), 4)
    dec = zeros(Int, 2)
    singles = want_singles ? zeros(Int, length(P3_COMPS)) : nothing
    ok = 0
    for i in runs_range
        _p3_orbits!(cnt, dec, p3_stream(q.test_seed, "pp-hb", (point, i)), S, nT, nRW)
        ok += _p3_decide(cnt, dec, D, nT, nRW, singles)
    end
    return ok, singles
end

# The runs 1:runs split into contiguous blocks, one task per block; the integer sums do not depend on the split.
function _p3_blocks(runs, threads)
    nb = threads > 1 && Threads.nthreads() > 1 ? min(runs, 4 * min(threads, Threads.nthreads())) : 1
    edges = round.(Int, range(0, runs; length = nb + 1))
    return [(edges[b] + 1):edges[b + 1] for b in 1:nb if edges[b + 1] > edges[b]]
end

function _p3_hb_sim(model, nT, nRW, q, point, runs, threads, want_singles)
    S = (_p3_stratum(model, :T), _p3_stratum(model, :RW))
    D = _p3_decision_data(model, q, nT + nRW)
    blocks = _p3_blocks(runs, threads)
    res = Vector{Any}(undef, length(blocks))
    if length(blocks) > 1
        Threads.@threads for b in eachindex(blocks)
            res[b] = _p3_hb_chunk(S, D, q, nT, nRW, point, blocks[b], want_singles)
        end
    else
        for b in eachindex(blocks); res[b] = _p3_hb_chunk(S, D, q, nT, nRW, point, blocks[b], want_singles); end
    end
    ok = sum(r[1] for r in res)
    return ok, want_singles ? sum(r[2] for r in res) : nothing
end

_p3_power(ok, runs) = (power = ok / runs, se = sqrt((ok / runs) * (1 - ok / runs) / runs))

"""
    p3_simulate_hb(model, nT, nRW, q; point, runs, threads) -> (successes, runs, power, se)

`runs` runs of the simulation of E56 e at one point of the path: run `i` draws the `nT` orbits of T_s and the `nRW` of R_W from the stream (pp-hb, (point, i)); the full
decision succeeds iff the 13 comparisons (exact McNemar tail at most `alpha_hb`, margin `rho_hb`) and the coverage of both sets hold. `model` is
`(rule = Dict(:T => (p, d), :RW => (p, d)), comp = Dict(X => (T = (p, s), RW = (p, s))), ts, w2)`.
"""
function p3_simulate_hb(model, nT::Integer, nRW::Integer, q; point::Integer, runs::Integer, threads::Integer = 1)
    ok, _ = _p3_hb_sim(model, nT, nRW, q, point, runs, threads, false)
    return merge((successes = ok, runs = runs), _p3_power(ok, runs))
end

"""
    p3_hb_path(model, q; runs, threads) -> Vector{NamedTuple}

The rows of the draw path (`p3_draw_path`) with the simulation of every allowed point (`successes, runs, power, se`); a point that is not allowed has `nothing` there.
"""
function p3_hb_path(model, q; runs::Integer, threads::Integer = 1)
    return map(p3_draw_path(q)) do x
        x.allowed || return merge(x, (successes = nothing, runs = runs, power = nothing, se = nothing))
        merge(x, p3_simulate_hb(model, x.nT, x.nRW, q; point = x.j, runs = runs, threads = threads))
    end
end

_p3_intish(x::Rational) = denominator(x) == 1 ? Int(numerator(x)) : x

"""
    p3_hb_targets(path, q) -> (outcome, cap, j, rho_star, n_T, n_W, power_cap)

The targets of E56 e from the simulated path: rho* is the smallest allowed point whose power reaches `p_hb` (exact comparison), `n_T = ceil5(y_t rho*)`, `n_W = ceil5(y_w rho*)`.
If no allowed point reaches it the targets are the cap `y_t r_max` and `y_w r_max`, `cap` is true and `power_cap` is the power of the last point; with no allowed point the
outcome is `:no_allowed_point`.
"""
function p3_hb_targets(path, q)
    ok = [x for x in path if x.allowed]
    isempty(ok) && return (outcome = :no_allowed_point, cap = false, j = nothing, rho_star = nothing, n_T = nothing, n_W = nothing, power_cap = nothing)
    p = p3_rat(q.p_hb)
    i = findfirst(x -> x.successes * denominator(p) >= numerator(p) * x.runs, ok)
    if i !== nothing
        x = ok[i]
        return (outcome = :ok, cap = false, j = x.j, rho_star = x.rho, n_T = p3_ceil_step(x.rho * q.y_t, q.ceil_step), n_W = p3_ceil_step(x.rho * q.y_w, q.ceil_step), power_cap = nothing)
    end
    last = path[end]
    return (outcome = :cap, cap = true, j = last.j, rho_star = last.rho, n_T = _p3_intish(p3_rat(q.r_max) * q.y_t), n_W = _p3_intish(p3_rat(q.r_max) * q.y_w),
            power_cap = last.allowed ? last.successes / last.runs : nothing)
end

"""
    p3_hb_type1(model, nT, nRW, q; point, runs, threads) -> (rates, max_joint)

The type-I check of E56 f at one point: for each comparator X in turn the same simulation (the same streams) with the rule rate of both strata set to the rate of X in that
stratum. `rates[X] = (joint, single)`: the rejection rate of the joint decision and of the comparisons of X alone; `max_joint` the largest joint rate.
"""
function p3_hb_type1(model, nT::Integer, nRW::Integer, q; point::Integer, runs::Integer, threads::Integer = 1)
    rates = Dict{Symbol, NamedTuple}()
    for (x, X) in enumerate(P3_COMPS)
        m = merge(model, (rule = Dict(:T => (p = model.comp[X][:T].p, d = model.rule[:T].d), :RW => (p = model.comp[X][:RW].p, d = model.rule[:RW].d)),))
        ok, singles = _p3_hb_sim(m, nT, nRW, q, point, runs, threads, true)
        rates[X] = (joint = ok / runs, single = singles[x] / runs)
    end
    return (rates = rates, max_joint = maximum(r.joint for r in values(rates)))
end
