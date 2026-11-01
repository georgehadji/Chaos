# D95h (HA), D95i (HB), D95k (HC) as pure functions on prescribed tables, D95m (streams and bootstrap), and the hand cases of O21 (a) to (g). Spec: analysis-spec.md, Amendment D95,
# E66 to E70, O21, criteria D95h, D95i, D95k, D95m. Expected values: hand tables, the exact binomial table, and the independent reference an_p3t_ref.jl (exact rationals, the
# generator of an_p3_ref.jl). API proposed by the tests: analysis-test-plan.md, "Amendment D95 (stage 02)", "Contract proposed by the tests".

const P3T_Q = merge(p3_params(; r_max = 3.0, n_t_min = 3, n_w_min = 5, y_t = 3, y_w = 5, bootstrap_resamples = 100), p3t_values())

@testset "D95i O21 (a) the exact one-sided McNemar value against the binomial table" begin
    for m in 0:12, n01 in 0:m
        n10 = m - n01
        want = m == 0 ? 1 // 1 : sum(binomial(big(m), i) for i in n01:m) // big(2)^m
        @test AN.p3_mcnemar_p(n10, n01) == want
    end
    @test AN.p3_mcnemar_p(5, 0) == 1 && AN.p3_mcnemar_p(0, 5) == 1 // 32 && AN.p3_mcnemar_p(0, 0) == 1
end

@testset "D95h O21 (b) the type 7 percentile in exact arithmetic" begin
    xs = collect(1 // 1:20 // 1)
    @test AN.p3t_percentile(xs, 19 // 20) == 19 + 1 // 20 && AN.p3t_percentile(xs, 1 // 20) == 1 + 19 // 20          # h = 19.05 and h = 1.95
    @test AN.p3t_percentile(Rational{Int}[1, 2, 3, 4, 5], 1 // 2) == 3 && AN.p3t_percentile(Rational{Int}[1, 2, 3, 4, 5], 1 // 4) == 2 && AN.p3t_percentile(Rational{Int}[3, 1, 2], 1) == 3
    @test AN.p3t_percentile(Rational{Int}[7], 19 // 20) == 7 && AN.p3t_percentile(Rational{Int}[], 1 // 2) === nothing
    # a replicate list with one dropped replicate: the percentile is taken over the 19 valid values
    v = vcat(xs[1:9], xs[11:20])                                           # 19 values, x = 10 dropped
    @test AN.p3t_percentile(v, 19 // 20) == p3tref_percentile(v, 19 // 20) && AN.p3t_percentile(v, 1 // 20) == p3tref_percentile(v, 1 // 20)
    # exact rationals out: a binary64 input is never needed, the result of 1/3 inputs is exact
    @test AN.p3t_percentile([1 // 3, 2 // 3, 1 // 1], 1 // 2) == 2 // 3 && AN.p3t_percentile([1 // 3, 2 // 3], 1 // 2) == 1 // 2
end

const P3T_ROC_SCORES = [1, 1 // 2, 1 // 2, 0, 0, 1 // 4, 1, 0]
const P3T_ROC_F = [true, true, false, false, false, true, false, false]

@testset "D95k O21 (c) the AUC and the ROC rows" begin
    @test AN.p3_auc([3, 3], [1, 1]) == 1 && AN.p3_auc([1, 1], [1, 1]) == 1 // 2 && AN.p3_auc([2, 2], [2, 1]) == 3 // 4 && AN.p3_auc(Int[], [1]) === nothing
    roc = AN.p3t_roc([3, 2, 2, 1], [true, true, false, false])
    @test [(r.threshold, r.fpr, r.tpr) for r in roc] == [(3, 0, 1 // 2), (2, 1 // 2, 1), (1, 1, 1)]
    r2 = AN.p3t_roc(P3T_ROC_SCORES, P3T_ROC_F)
    @test [(r.threshold, r.fpr, r.tpr) for r in r2] == [(r.threshold, r.fpr, r.tpr) for r in p3tref_roc(P3T_ROC_SCORES, P3T_ROC_F)]
end

@testset "D95m O21 (d) the stratified resample, the streams, the labels" begin
    # two strata of 2 and 3 orbits: the integers of the streams (boot-ha, 1) and (boot-ha, 2) written by the independent Python generator (streams.py)
    @test AN.p3t_resample(20261030, "boot-ha", 1, [2, 3]) == [[2, 1], [3, 1, 3]] && AN.p3t_resample(20261030, "boot-ha", 2, [2, 3]) == [[1, 2], [1, 3, 2]]
    # the statistic by hand: share of unresolved orbits of the resample, a1 and b1 unresolved: replicate 1 draws a2 a1 | b3 b1 b3 (2 of 5), replicate 2 draws a1 a2 | b1 b3 b2 (2 of 5)
    unres = [[true, false], [true, false, false]]
    share(r) = (d = AN.p3t_resample(20261030, "boot-ha", r, [2, 3]); (sum(unres[1][i] for i in d[1]) + sum(unres[2][i] for i in d[2])) // 5)
    @test share(1) == 2 // 5 && share(2) == 2 // 5
    # a stream does not change when another unit or purpose is drawn, or when it is drawn again
    a = AN.p3t_resample(20261030, "boot-hc", 5, [7]); AN.p3t_resample(20261030, "boot-ha", 6, [7]); AN.p3t_resample(20261030, "boot-hc", 4, [7])
    @test AN.p3t_resample(20261030, "boot-hc", 5, [7]) == a == p3tref_resample(20261030, "boot-hc", 5, [7])
    # the first four 64-bit outputs of the four streams (unit 1) written by the independent generator; the seeds
    first4(p) = (g = AN.p3_stream(20261030, p, 1); [AN.p3_u64!(g) for _ in 1:4])
    @test first4("boot-ha") == UInt64[8392252990483039953, 10316736585100517932, 15353738929327571399, 11852175026992284072]
    @test first4("boot-hc") == UInt64[8318303221503667656, 5078820220052927605, 5646466915089531211, 4378077879993873153]
    @test first4("pp-achieved-hb") == UInt64[12428499709038559142, 891794636315231885, 11446961040307043218, 10850844614071461781]
    @test first4("pp-achieved-hc") == UInt64[9456558871754260945, 643522912293856473, 4585469879774012822, 17960362725938284023]
    @test AN.p3_stream_seed(20261030, "boot-ha", 1) == 6130294642297488880 && AN.p3_stream_seed(20261030, "boot-hc", 1) == 15727976254305454371 &&
          AN.p3_stream_seed(20261030, "pp-achieved-hb", 1) == 4500577164153853627 && AN.p3_stream_seed(20261030, "pp-achieved-hc", 1) == 4279710978371254201
    # the labels differ from those of D93 and D89 (optimism, pp-hb, pp-hc) and from each other
    seeds = [AN.p3_stream_seed(20261030, p, 1) for p in ("boot-ha", "boot-hc", "pp-achieved-hb", "pp-achieved-hc", "optimism", "pp-hb", "pp-hc")]
    @test allunique(seeds)
    g = AN.p3_stream(20261030, "boot-ha", 1)
    @test [AN.p3_randint!(g, 25) for _ in 1:8] == [4, 8, 25, 23, 10, 17, 23, 21]
end

@testset "D95h O21 (e) the strata and the count rule (g)" begin
    row(rho, Ts, W2, S2) = (id = 1, rho = rho, out = :unresolved, Ts = Ts, W2 = W2, S2 = S2)
    @test AN.p3t_stratum(row(:C, true, true, true)) == 1 && AN.p3t_stratum(row(:C, true, true, false)) == 1 && AN.p3t_stratum(row(:C, false, true, true)) == 2 &&
          AN.p3t_stratum(row(:C, false, true, false)) == 2 && AN.p3t_stratum(row(:C, false, false, true)) == 3 && AN.p3t_stratum(row(:C, false, false, false)) == 4 && AN.p3t_stratum(row(:R, false, false, false)) == 5
    # the count rule: targets reached, not reached, round r_max, one short by one in either count
    @test AN.p3t_decide(5, 11, 5, 11, 1, 3) == "close" && AN.p3t_decide(6, 12, 5, 11, 1, 3) == "close" && AN.p3t_decide(4, 11, 5, 11, 1, 3) == "continue" &&
          AN.p3t_decide(5, 10, 5, 11, 2, 3) == "continue" && AN.p3t_decide(0, 0, 5, 11, 3, 3) == "close" && AN.p3t_decide(4, 10, 5, 11, 3, 3) == "close" && AN.p3t_decide(4, 10, 5, 11, 2, 3) == "continue"
end

# ---------------------------------------------------------------- HA on prescribed tables
function p3t_ha_rows(; nT = 6, nT_unres = 6, nW = 6, nW_unres = 3, nS = 8, nS_unres = 0, nR = 60, nR_unres = 0, nR_bad = 0, overlap_S2 = 0)
    rows = NamedTuple[]; id = 0
    add(rho, out, Ts, W2, S2) = push!(rows, (id = (id += 1), rho = rho, out = out, W2 = W2, Ts = Ts, S2 = S2))
    for i in 1:nT; add(:C, i <= nT_unres ? :unresolved : :chaotic_supported, true, true, i <= overlap_S2); end
    for i in 1:nW; add(:C, i <= nW_unres ? :unresolved : :chaotic_supported, false, true, false); end
    for i in 1:nS; add(:C, i <= nS_unres ? :unresolved : :chaotic_supported, false, false, true); end
    for i in 1:nR; add(:R, i <= nR_unres ? :unresolved : i <= nR_unres + nR_bad ? :chaotic_supported : :regular_supported, false, false, false); end
    return rows
end
const P3T_GATE_OK = [:regular_supported, :regular_supported]

@testset "D95h HA: all items pass, each item fails alone, not decidable" begin
    q = P3T_Q
    function check(name, rows, gate, want_dec, want_failed; q = q)
        r = AN.p3t_ha(rows, gate, q); x = p3tref_ha(rows, gate, q)
        @test r.decision == want_dec == x.decision && r.failed == want_failed == x.failed
        @test r.point.e == x.point.e && r.point.aE == x.point.aE && r.point.dw == x.point.dw && r.point.dt == x.point.dt && r.point.g == x.point.g && r.point.a[:S2] == x.point.a[:S2] && r.point.a[:W2] == x.point.a[:W2] && r.point.a[:Ts] == x.point.a[:Ts]
        for k in (:e_up, :aE_up, :dw_lo, :dt_lo)
            @test getfield(r.bounds, k).value == getfield(x.bounds, k).value && getfield(r.bounds, k).dropped == getfield(x.bounds, k).dropped
        end
        @test r.strata == x.strata
        return r
    end
    r = check("all pass", p3t_ha_rows(), P3T_GATE_OK, "supported", Int[])
    @test r.strata == [6, 6, 8, 0, 60] && sum(r.strata) == 80 && r.point.g == 0 && r.point.e == 0 && r.point.aE == 9 // 80
    check("item 1: error rate", p3t_ha_rows(nR_bad = 30), P3T_GATE_OK, "not supported", [1])
    check("item 2: unresolved share", p3t_ha_rows(nR_unres = 50), P3T_GATE_OK, "not supported", [2])
    check("item 3: the T_s difference", p3t_ha_rows(nT_unres = 0), P3T_GATE_OK, "not supported", [3])
    check("item 4: the gate", p3t_ha_rows(), [:chaotic_supported, :chaotic_supported], "not supported", [4])
    check("g at g_max", p3t_ha_rows(), [:chaotic_supported, :regular_supported], "supported", Int[])                  # g = 1 = g_max passes
    # not decidable: S2 below min_set (4), no decided orbit (e undefined), a gate orbit without a verdict (P100), a stratum overlap that empties S2 in more than 5 percent of the replicates
    check("S2 below min_set", p3t_ha_rows(nS = 4), P3T_GATE_OK, "not decidable", Int[])
    check("S2 at min_set", p3t_ha_rows(nS = 5), P3T_GATE_OK, "supported", Int[])
    check("T_s below min_set", p3t_ha_rows(nT = 4, nT_unres = 4), P3T_GATE_OK, "not decidable", Int[])
    check("no decided orbit", p3t_ha_rows(nT_unres = 6, nW_unres = 6, nS_unres = 8, nR_unres = 60), P3T_GATE_OK, "not decidable", Int[])
    check("gate orbit without verdict", p3t_ha_rows(), Union{Symbol, Nothing}[:regular_supported, nothing], "not decidable", Int[])
    rows = p3t_ha_rows(nT = 20, nT_unres = 10, nW = 6, nS = 0, nR = 40, overlap_S2 = 1)               # S2 is the single orbit 1 inside the T_s stratum: a replicate that misses it has no S2
    nd = check("dropped replicates", rows, P3T_GATE_OK, "not decidable", Int[]; q = merge(q, (; min_set = 1)))
    @test nd.bounds.dw_lo.value === nothing && nd.bounds.dw_lo.dropped * 20 > 100
end

# ---------------------------------------------------------------- HB on prescribed tables
const P3T_HB_T = [:SALI, :FLI, :MEGNO, :E1, :E2, :WBA]
const P3T_HB_W = vcat(:LE, P3T_HB_T)
function p3t_hb_base(nT, nW)
    cmp = Dict{Tuple{Symbol, Symbol}, Any}()
    for X in P3T_HB_T; cmp[(:Ts, X)] = (rule = fill(false, nT), x = fill(true, nT)); end
    for X in P3T_HB_W; cmp[(:W2, X)] = (rule = fill(false, nW), x = fill(true, nW)); end
    return cmp
end
# a comparison with the counts n11 n10 n01 n00 in this order
function p3t_hb_counts(n11, n10, n01, n00)
    rule = vcat(fill(true, n11), fill(true, n10), fill(false, n01), fill(false, n00)); x = vcat(fill(true, n11), fill(false, n10), fill(true, n01), fill(false, n00))
    return (rule = rule, x = x)
end

@testset "D95i HB: the 13 comparisons, the margin, the p value, coverage, precedence" begin
    q = P3T_Q
    cov = Dict{Symbol, Any}(:Ts => 1 // 1, :W2 => 1 // 1)
    nT, nW = 12, 14
    base = p3t_hb_base(nT, nW)
    r = AN.p3t_hb(base, cov, (Ts = nT, W2 = nW), q); x = p3tref_hb(base, cov, (Ts = nT, W2 = nW), q)
    @test r.decision == "supported" == x.decision && isempty(r.failing) && r.coverage_ok && length(r.rows) == 13 && r.underpowered == (nW < q.target_w2)
    @test r.rows[(:Ts, :SALI)].n01 == nT && r.rows[(:Ts, :SALI)].p == 1 // 2^nT && r.rows[(:W2, :LE)].p == 1 // 2^nW
    # each of the 13 fails alone: by the p value (n01 = 3 of 12, n10 = 0: p = 1/8) and by the margin only (n11 = 7, n01 = 5: p = 1/32 passes, FR_rule 7/12 > half of FR_X = 1)
    for (A, X) in sort(collect(keys(base)); by = string)
        n = A === :Ts ? nT : nW
        for (kind, cnt) in (("p value", p3t_hb_counts(0, 0, 3, n - 3)), ("margin", p3t_hb_counts(7, 0, 5, n - 12)), ("FR_X = 0", p3t_hb_counts(0, 0, 0, n)), ("n10 above n01", p3t_hb_counts(0, 4, 2, n - 6)))
            cm = copy(base); cm[(A, X)] = cnt
            r = AN.p3t_hb(cm, cov, (Ts = nT, W2 = nW), q); x = p3tref_hb(cm, cov, (Ts = nT, W2 = nW), q)
            @test r.decision == "not supported" == x.decision && r.failing == [(A, X)] == x.failing
            @test r.rows[(A, X)].p == x.rows[(A, X)].p && r.rows[(A, X)].margin == x.rows[(A, X)].margin && r.rows[(A, X)].pass == false
        end
    end
    # the rows of a comparison: counts sum to |A| (I43), m = 0 gives p = 1, FR_X = 0 fails even when the p value is 1
    c0 = AN.p3t_hb_compare(fill(false, 5), fill(false, 5), q)
    @test c0.p == 1 && !c0.pass && c0.n00 == 5 && c0.n11 + c0.n10 + c0.n01 + c0.n00 == 5
    c1 = AN.p3t_hb_compare(p3t_hb_counts(2, 1, 4, 3).rule, p3t_hb_counts(2, 1, 4, 3).x, q)
    @test (c1.n11, c1.n10, c1.n01, c1.n00) == (2, 1, 4, 3) && c1.p == p3ref_mcnemar(1, 4) && c1.fr_rule == 3 // 10 && c1.fr_x == 6 // 10 && c1.margin
    # the coverage condition fails alone: not supported
    r = AN.p3t_hb(base, Dict{Symbol, Any}(:Ts => 1 // 2, :W2 => 49 // 100), (Ts = nT, W2 = nW), q)
    @test r.decision == "not supported" && !r.coverage_ok && isempty(r.failing)
    r = AN.p3t_hb(base, Dict{Symbol, Any}(:Ts => 1 // 2, :W2 => 1 // 2), (Ts = nT, W2 = nW), q)
    @test r.decision == "supported"                                                                                       # coverage at c_min exactly passes (>=)
    # not decidable: the floors (|T_s| < n_t_min, |W2| < n_w_min) and a capped non-LE comparator (a defect, named); the floors win over a failing comparison
    @test AN.p3t_hb(p3t_hb_base(2, nW), cov, (Ts = 2, W2 = nW), q).decision == "not decidable" && AN.p3t_hb(p3t_hb_base(nT, 4), cov, (Ts = nT, W2 = 4), q).decision == "not decidable"
    @test AN.p3t_hb(p3t_hb_base(3, 5), cov, (Ts = 3, W2 = 5), q).decision == "not supported"                             # exactly at the floors: decidable (p = 1/8 > alpha with 3 orbits)
    cm = copy(base); cm[(:Ts, :E1)] = p3t_hb_counts(0, 0, 0, nT)
    @test AN.p3t_hb(cm, cov, (Ts = 2, W2 = nW), q).decision == "not decidable"
    d = AN.p3t_hb(base, cov, (Ts = nT, W2 = nW), q; capped = [:E1])
    @test d.decision == "not decidable" && :E1 in d.defects
    # the underpowered flag follows the final sizes against the targets, whatever the outcome
    @test AN.p3t_hb(base, cov, (Ts = 4, W2 = 14), q).underpowered == true && AN.p3t_hb(p3t_hb_base(5, 11), cov, (Ts = 5, W2 = 11), q).underpowered == false &&
          AN.p3t_hb(p3t_hb_base(5, 10), cov, (Ts = 5, W2 = 10), q).underpowered == true
end

# ---------------------------------------------------------------- HC on prescribed tables
function p3t_hc_rows(; nf = 6, nn = 6, kind = :supported, off = 0)
    rows = NamedTuple[]; id = off
    for i in 1:nf          # failures: a tie of votes (kind :supported/:not_supported) or a wrong plain label (kind :reversed)
        kind === :reversed ? push!(rows, (id = (id += 1), rho = :R, nC = 1, nR = 0, sP = 1)) : push!(rows, (id = (id += 1), rho = :C, nC = 2, nR = 2, sP = kind === :not_supported ? 1 : 0))
    end
    for i in 1:nn          # non-failures: a plain label equal to rho
        kind === :reversed ? push!(rows, (id = (id += 1), rho = :C, nC = 1, nR = 0, sP = 0)) : push!(rows, (id = (id += 1), rho = :C, nC = 4, nR = 0, sP = 0))
    end
    return rows
end

@testset "D95k HC: f, s_I, s_P, the AUCs, the bounds, the outcomes" begin
    q = P3T_Q
    function check(rows, want; q = q)
        r = AN.p3t_hc(rows, q); x = p3tref_hc(rows, q)
        @test r.decision == want == x.decision && r.failures == x.failures && r.nonfailures == x.nonfailures && r.auc_i == x.auc_i && r.auc_p == x.auc_p && r.delta == x.delta &&
              r.lower == x.lower && r.upper == x.upper && r.dropped == x.dropped && r.underpowered == x.underpowered
        return r
    end
    r = check(p3t_hc_rows(), "supported")
    @test r.auc_i == 1 && r.auc_p == 1 // 2 && r.delta == 1 // 2 && r.lower > 0 && r.failures == 6 && r.nonfailures == 6 && !r.underpowered
    r = check(p3t_hc_rows(kind = :reversed), "reversed")
    @test r.auc_i == 1 // 2 && r.auc_p == 1 && r.delta == -1 // 2 && r.upper < 0
    r = check(p3t_hc_rows(kind = :not_supported), "not supported")
    @test r.delta == 0 && r.auc_i == r.auc_p == 1
    check(p3t_hc_rows(nf = 4), "not decidable"); check(p3t_hc_rows(nn = 4), "not decidable")                           # failures or non-failures at min_set - 1
    r = check(p3t_hc_rows(nf = 5), "supported"); @test r.underpowered && r.failures == 5                                    # at min_set: decided, underpowered (target_f = 6)
    check(p3t_hc_rows(nf = 5, nn = 5), "supported")
    # the plain label and the scores by hand: votes (nC, nR), rho
    t = [(id = 1, rho = :C, nC = 2, nR = 2, sP = 0), (id = 2, rho = :C, nC = 0, nR = 0, sP = 0), (id = 3, rho = :R, nC = 3, nR = 1, sP = 1 // 4), (id = 4, rho = :C, nC = 3, nR = 1, sP = 1 // 2),
         (id = 5, rho = :R, nC = 0, nR = 2, sP = 0), (id = 6, rho = :C, nC = 1, nR = 2, sP = 1 // 4)]
    rows = AN.p3t_hc(vcat(t, p3t_hc_rows(nf = 5, nn = 5, off = 100)), q).rows
    @test [rows[i].f for i in 1:6] == [true, true, true, false, false, true]                                             # tie, no vote (a tie), wrong label, right, right, wrong
    @test [rows[i].sI for i in 1:6] == [1, 0, 1 // 2, 1 // 2, 0, 2 // 3] && [rows[i].sP for i in 1:6] == [0, 0, 1 // 4, 1 // 2, 0, 1 // 4]
    # a two-orbit class that is empty after a resample is dropped, not counted: all failures in one orbit
    few = vcat([(id = 1, rho = :C, nC = 2, nR = 2, sP = 0)], [(id = 1 + i, rho = :C, nC = 4, nR = 0, sP = 0) for i in 1:19])
    r = AN.p3t_hc(few, merge(q, (; min_set = 1))); x = p3tref_hc(few, merge(q, (; min_set = 1)))
    @test r.dropped == x.dropped > 5 && r.decision == "not decidable" == x.decision && r.lower === nothing
end
