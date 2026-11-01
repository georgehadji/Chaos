# Amendment D95, stage 02. Independent reference implementation of the test side of study p3-abstention (spec E62 to E70, O21), written from the amendment text. It never calls
# Chaoslab.Analysis, so it works before the module has any D95 functionality. Names: p3tref_*. Exact rationals (BigInt) for every count-based quantity (I42). The generator, the
# streams and the exact McNemar tail are those of the D93 reference (an_p3_ref.jl: p3ref_stream, p3ref_randint!, p3ref_mcnemar, p3ref_auc, p3ref_out, p3ref_matched).
# Verdict symbols :chaotic :regular :undecided, rule outputs :chaotic_supported :regular_supported :unresolved.

const P3T_R = Rational{BigInt}
p3tref_rat(x) = P3T_R(p3ref_rat(x))

# ---------------------------------------------------------------- the round log (E62 a, b): the text of an entry, the decision rule
# counts: Dict name => Integer with the seven names drift excluded orbits t_s target_t_s target_w2 w2; the count lines in ascending byte order of their names
function p3tref_entry(; round, decision, counts, panel, prereg, digest, prev)
    io = IOBuffer()
    println(io, "round=", round)
    println(io, "decision=", decision)
    for k in sort(collect(keys(counts)))
        println(io, "count.", k, "=", counts[k])
    end
    println(io, "panel_sha256=", panel)
    println(io, "prereg_sha256=", prereg)
    println(io, "ref_digest=", digest)
    println(io, "prev_sha256=", prev)
    return String(take!(io))
end
p3tref_decide(t_s, w2, target_t, target_w, r, r_max) = (t_s >= target_t && w2 >= target_w) || r >= r_max ? "close" : "continue"

# ---------------------------------------------------------------- type 7 percentile in exact arithmetic (E66 e)
function p3tref_percentile(xs, p)
    x = sort(P3T_R.(xs)); n = length(x)
    n == 0 && return nothing
    h = (n - 1) * P3T_R(p) + 1
    j = Int(floor(h))
    (h == j || j >= n) && return x[min(j, n)]
    return x[j] + (h - j) * (x[j + 1] - x[j])
end

# Mann-Whitney AUC of a score for the failures f (Bool vector); a tie counts one half; nothing for an empty class
function p3tref_auc(score, f)
    pos = [s for (s, y) in zip(score, f) if y]; neg = [s for (s, y) in zip(score, f) if !y]
    (isempty(pos) || isempty(neg)) && return nothing
    s = sum(x -> sum(y -> x > y ? 2 : x == y ? 1 : 0, neg), pos)
    return P3T_R(s) / (2 * length(pos) * length(neg))
end
# ROC rows (threshold, fpr, tpr): for each distinct score in decreasing order the share of non-failures and of failures with score >= threshold
function p3tref_roc(score, f)
    nf = count(f); nn = length(f) - nf
    return [(threshold = t, fpr = P3T_R(count(i -> score[i] >= t && !f[i], eachindex(f))) / nn, tpr = P3T_R(count(i -> score[i] >= t && f[i], eachindex(f))) / nf)
            for t in sort(unique(score); rev = true)]
end

# ---------------------------------------------------------------- the stratified resample (E66 b): sizes = strata sizes in the fixed order; replicate r from the stream (purpose, r)
function p3tref_resample(test_seed, purpose, r, sizes)
    g = p3ref_stream(test_seed, purpose, r)
    return [[p3ref_randint!(g, n) for _ in 1:n] for n in sizes]
end
# the five strata of E64 d: first T_s, then W2 minus T_s, then S2, then the other C, then R
p3tref_stratum(o) = o.rho === :R ? 5 : o.Ts ? 1 : o.W2 ? 2 : o.S2 ? 3 : 4

# ---------------------------------------------------------------- HA (E67). orbits: rows (id, rho, out, W2, Ts, S2) sorted by id; gate: the outputs of the gate orbits (nothing = no verdict)
# returns the point values, the bounds (nothing if undefined), the numbers of dropped replicates per statistic, the decision and the failed items
function p3tref_ha(orbits, gate, q)
    o = sort(collect(orbits); by = x -> x.id)
    A = Dict(:E => o, :W2 => [x for x in o if x.W2], :Ts => [x for x in o if x.Ts], :S2 => [x for x in o if x.S2])
    contra(xs) = count(x -> (x.out === :chaotic_supported && x.rho === :R) || (x.out === :regular_supported && x.rho === :C), xs)
    dec(xs) = count(x -> x.out !== :unresolved, xs)
    e_of(xs) = dec(xs) == 0 ? nothing : P3T_R(contra(xs)) / dec(xs)
    a_of(xs) = isempty(xs) ? nothing : P3T_R(count(x -> x.out === :unresolved, xs)) / length(xs)
    dif(xs, ys) = (a = a_of(xs); b = a_of(ys); (a === nothing || b === nothing) ? nothing : a - b)
    function stats(os)
        W2 = [x for x in os if x.W2]; Ts = [x for x in os if x.Ts]; S2 = [x for x in os if x.S2]
        return (e = e_of(os), aE = a_of(os), dw = dif(W2, S2), dt = dif(Ts, S2))
    end
    pt = (e = e_of(o), aE = a_of(o), a = Dict(k => a_of(v) for (k, v) in A if k !== :E), dw = dif(A[:W2], A[:S2]), dt = dif(A[:Ts], A[:S2]), g = count(x -> x === :chaotic_supported, gate))
    strata = [[x for x in o if p3tref_stratum(x) == s] for s in 1:5]
    reps = [stats(reduce(vcat, [st[idx] for (st, idx) in zip(strata, p3tref_resample(q.test_seed, "boot-ha", r, length.(strata)))])) for r in 1:q.bootstrap_resamples]
    alpha = p3tref_rat(q.alpha_ci)
    bound(name, p) = (v = [getfield(s, name) for s in reps if getfield(s, name) !== nothing]; dropped = length(reps) - length(v);
                      (value = dropped * 20 > length(reps) ? nothing : p3tref_percentile(v, p), dropped = dropped))
    b = (e_up = bound(:e, 1 - alpha), aE_up = bound(:aE, 1 - alpha), dw_lo = bound(:dw, alpha), dt_lo = bound(:dt, alpha))
    nd = any(k -> length(A[k]) < q.min_set, keys(A)) || pt.e === nothing || any(x -> x === nothing, gate) || any(v -> v.value === nothing, values(b))
    failed = Int[]
    if !nd
        b.e_up.value <= p3tref_rat(q.e_max) || push!(failed, 1)
        b.aE_up.value <= p3tref_rat(q.a_max) || push!(failed, 2)
        (pt.dw >= p3tref_rat(q.delta_a) && pt.dt >= p3tref_rat(q.delta_a) && b.dw_lo.value > 0 && b.dt_lo.value > 0) || push!(failed, 3)
        pt.g <= q.g_max || push!(failed, 4)
    end
    return (point = pt, bounds = b, decision = nd ? "not decidable" : isempty(failed) ? "supported" : "not supported", failed = failed, strata = length.(strata))
end

# ---------------------------------------------------------------- HB (E68). cmp[(A, X)] = (rule = Bool vector, x = Bool vector) of false-regular flags on the same orbits
function p3tref_hb_compare(rule, x, q)
    n11 = count(i -> rule[i] && x[i], eachindex(rule)); n10 = count(i -> rule[i] && !x[i], eachindex(rule))
    n01 = count(i -> !rule[i] && x[i], eachindex(rule)); n00 = length(rule) - n11 - n10 - n01
    n = length(rule)
    frr = P3T_R(n11 + n10) // n; frx = P3T_R(n11 + n01) // n
    p = p3ref_mcnemar(n10, n01)
    margin = frx > 0 && frr <= p3tref_rat(q.rho_hb) * frx
    return (n11 = n11, n10 = n10, n01 = n01, n00 = n00, p = p, margin = margin, pass = p <= p3tref_rat(q.alpha_hb) && margin, fr_rule = frr, fr_x = frx)
end
function p3tref_hb(cmp, cov, sizes, q; capped = Symbol[])
    rows = Dict(k => p3tref_hb_compare(v.rule, v.x, q) for (k, v) in cmp)
    covok = cov[:Ts] >= p3tref_rat(q.c_min) && cov[:W2] >= p3tref_rat(q.c_min)
    nd = sizes.Ts < q.n_t_min || sizes.W2 < q.n_w_min || !isempty(capped)
    failing = sort([k for (k, r) in rows if !r.pass]; by = string)
    return (rows = rows, coverage_ok = covok, decision = nd ? "not decidable" : (isempty(failing) && covok) ? "supported" : "not supported", failing = failing,
            underpowered = (sizes.Ts < q.target_t_s || sizes.W2 < q.target_w2))
end

# ---------------------------------------------------------------- HC (E70). rows: (id, rho, nC, nR, sP) ; the plain label, the failure f, the scores s_I and s_P
function p3tref_hc_rows(rows)
    out = NamedTuple[]
    for r in sort(collect(rows); by = x -> x.id)
        lab = r.nC > r.nR ? :chaotic : r.nR > r.nC ? :regular : :tie
        f = lab === :tie || (lab === :chaotic && r.rho === :R) || (lab === :regular && r.rho === :C)
        sI = r.nC + r.nR == 0 ? P3T_R(0) : 1 - P3T_R(abs(r.nC - r.nR)) / (r.nC + r.nR)
        push!(out, (id = r.id, f = f, sI = sI, sP = P3T_R(r.sP)))
    end
    return out
end
function p3tref_hc(rows, q)
    R = p3tref_hc_rows(rows)
    f = [r.f for r in R]
    aI = p3tref_auc([r.sI for r in R], f); aP = p3tref_auc([r.sP for r in R], f)
    delta = (aI === nothing || aP === nothing) ? nothing : aI - aP
    n = length(R)
    reps = P3T_R[]
    for r in 1:q.bootstrap_resamples
        g = p3ref_stream(q.test_seed, "boot-hc", r)
        idx = [p3ref_randint!(g, n) for _ in 1:n]
        ff = f[idx]
        a = p3tref_auc([R[i].sI for i in idx], ff); b = p3tref_auc([R[i].sP for i in idx], ff)
        (a === nothing || b === nothing) || push!(reps, a - b)
    end
    dropped = q.bootstrap_resamples - length(reps)
    alpha = p3tref_rat(q.alpha_ci)
    lo = dropped * 20 > q.bootstrap_resamples ? nothing : p3tref_percentile(reps, alpha)
    up = dropped * 20 > q.bootstrap_resamples ? nothing : p3tref_percentile(reps, 1 - alpha)
    nf = count(f); nn = n - nf
    nd = nf < q.min_set || nn < q.min_set || lo === nothing || up === nothing || delta === nothing
    dc = p3tref_rat(q.delta_c)
    dec = nd ? "not decidable" : (delta >= dc && lo > 0) ? "supported" : (delta <= -dc && up < 0) ? "reversed" : "not supported"
    return (rows = R, auc_i = aI, auc_p = aP, delta = delta, lower = lo, upper = up, dropped = dropped, failures = nf, nonfailures = nn, decision = dec,
            underpowered = nf < q.target_f)
end

# ---------------------------------------------------------------- the power at the achieved sizes (E69 c): one simulated run of the full HB decision of D93 E56 c with the stream (pp-achieved-hb, run index)
function p3tref_hb_run(q, model, nT, nRW, i)
    g = p3ref_stream(q.test_seed, "pp-achieved-hb", i)
    sets = Dict(X => Dict(:T => zeros(Int, 4), :RW => zeros(Int, 4)) for X in P2_X7)     # n11, n10, n01, n00
    dec = Dict(:T => 0, :RW => 0)
    for (st, n) in ((:T, nT), (:RW, nRW)), _ in 1:n
        o = p3ref_hb_orbit!(g, st, model)
        o.decided && (dec[st] += 1)
        for X in P2_X7
            idx = o.rfr ? (o.comp[X] ? 1 : 2) : (o.comp[X] ? 3 : 4)
            sets[X][st][idx] += 1
        end
    end
    alpha, rho_hb, c_min = p3ref_rat(q.alpha_hb), p3ref_rat(q.rho_hb), p3ref_rat(q.c_min)
    function ok(cnt)
        n11, n10, n01, _ = cnt
        m = n10 + n01; a = n10 + n11; c = n01 + n11
        return m >= 1 && p3ref_mcnemar(n10, n01) <= alpha && c >= 1 && a <= rho_hb * c
    end
    good = all(X -> ok(sets[X][:T]), model.ts) && all(X -> ok(sets[X][:T] .+ sets[X][:RW]), model.w2)
    cov = dec[:T] >= c_min * nT && dec[:T] + dec[:RW] >= c_min * (nT + nRW)
    return good && cov
end
p3tref_power_hb(q, model, nT, nRW, runs) = count(i -> p3tref_hb_run(q, model, nT, nRW, i), 1:runs)
