# Study p3-abstention (Amendment D93): PP for HC (E57). The two predictors are independent binormal samples; every run owns a named stream (pp-hc, (n_f, run)); the decision uses the
# Hanley and McNeil variances at the empirical AUCs (P64). The platform logarithm enters the normal variates (E55): byte identity across machines is not claimed.

# One run: the doubled Mann-Whitney counts of the two predictors from the draws in the order of E57 (failures of I, non-failures of I, failures of P, non-failures of P), and their common denominator.
function _p3_hc_run(g, dI, dP, nf, nn)
    fI = [dI + p3_normal!(g) for _ in 1:nf]; nI = [p3_normal!(g) for _ in 1:nn]
    fP = [dP + p3_normal!(g) for _ in 1:nf]; nP = [p3_normal!(g) for _ in 1:nn]
    return _p3_auc_count(fI, nI), _p3_auc_count(fP, nP), 2 * nf * nn
end

function _p3_hc_chunk(q, dI, dP, nf, runs_range)
    nn = q.nonfail_per_fail * nf
    ok = 0
    dc = p3_rat(q.delta_c)
    for i in runs_range
        cI, cP, d = _p3_hc_run(p3_stream(q.test_seed, "pp-hc", (nf, i)), dI, dP, nf, nn)
        AI, AP = cI / d, cP / d
        Delta = AI - AP
        SE = sqrt(p3_hanley_mcneil(AI, nf, nn) + p3_hanley_mcneil(AP, nf, nn))
        # Delta >= delta_c exactly: both AUCs share the denominator d, so (cI - cP) / d >= delta_c is an integer against a rational.
        ((cI - cP) * denominator(dc) >= numerator(dc) * d && Delta - q.z_one_sided * SE > 0) && (ok += 1)
    end
    return ok
end

"""
    p3_simulate_hc(q, auc_i, auc_p; nf, runs, threads) -> (successes, runs, power, se)

`runs` runs of the simulation of E57 with `nf` failures and `nonfail_per_fail * nf` non-failures: binormal scores with the shifts `d(auc_i)` and `d(auc_p)`; the decision succeeds
iff `Delta >= delta_c` and `Delta - z_one_sided SE > 0`. An `ArgumentError` naming the AUC outside the open unit interval.
"""
function p3_simulate_hc(q, auc_i::Real, auc_p::Real; nf::Integer, runs::Integer, threads::Integer = 1)
    dI, dP = p3_binormal_shift(auc_i), p3_binormal_shift(auc_p)
    blocks = _p3_blocks(runs, threads)
    res = zeros(Int, length(blocks))
    if length(blocks) > 1
        Threads.@threads for b in eachindex(blocks)
            res[b] = _p3_hc_chunk(q, dI, dP, nf, blocks[b])
        end
    else
        for b in eachindex(blocks); res[b] = _p3_hc_chunk(q, dI, dP, nf, blocks[b]); end
    end
    ok = sum(res)
    return merge((successes = ok, runs = runs), _p3_power(ok, runs))
end

"""
    p3_pp_hc(q, auc_i_cal; runs, threads) -> NamedTuple

PP for HC (E57): `auc_i_used = 0.5 + auc_shrink (auc_i_cal - 0.5)` against `auc_p`; the `table` of `(nf, successes, runs, power, se)` for `nf_start`, `nf_second`, ... up to
the smallest listed value at least the worst-case deliverable failures; the deliverable numbers `floor(d_total / (1 + nonfail_per_fail))` for the expected and the worst case; the
power at the expected number; the `target` (the first listed `nf` with power at least `p_hc`, exact) or `nothing`; `underpowered` if the target is `nothing` or above the expected
deliverable number; `power_cap` (the power at the last listed value) when the target is not reached.
"""
function p3_pp_hc(q, auc_i_cal::Real; runs::Integer, threads::Integer = 1)
    used = 0.5 + q.auc_shrink * (auc_i_cal - 0.5)
    dI, dP = p3_binormal_shift(used), p3_binormal_shift(q.auc_p)
    step, per = q.nf_second - q.nf_start, 1 + q.nonfail_per_fail
    expected, worst = fld(q.d_expected_total, per), fld(q.d_worst_total, per)
    last = q.nf_start + step * max(0, cld(worst - q.nf_start, step))
    row(nf) = merge((nf = nf,), p3_simulate_hc(q, used, q.auc_p; nf = nf, runs = runs, threads = threads))
    table = [row(nf) for nf in q.nf_start:step:last]
    pe = expected >= 1 ? (i = findfirst(r -> r.nf == expected, table); i === nothing ? row(expected) : table[i]) : nothing
    p = p3_rat(q.p_hc)
    hit = findfirst(r -> r.successes * denominator(p) >= numerator(p) * r.runs, table)
    target = hit === nothing ? nothing : table[hit].nf
    return (auc_i_used = used, auc_p = q.auc_p, true_difference = used - q.auc_p, d_i = dI, d_p = dP, table = table, deliverable_expected = expected, deliverable_worst = worst,
            power_expected = pe, target = target, target_power = hit === nothing ? nothing : table[hit].power, underpowered = target === nothing || target > expected,
            power_cap = hit === nothing ? table[end].power : nothing)
end
