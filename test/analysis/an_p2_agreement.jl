# D89i agreement, kappa, U, A_min, Spearman (E30) and D89j AUC, regular range, leaving budget (E31). Spec: analysis-spec.md, Amendment D89, O7, criteria D89i, D89j.
# Oracles: the reference of an_p2_ref.jl (independent re-implementation) on 1000 random tables, and hand tables with the arithmetic in the comments.

const P2_VSYM = (:chaotic, :regular, :undecided, :na)
p2_close(a, b; rtol = 1e-12) = (a === nothing || b === nothing) ? a === b : isapprox(a, b; rtol = rtol, atol = 1e-15)

@testset "D89i contingency table, agreement and Cohen kappa: hand tables (O7 a)" begin
    tab(a, b, c, d) = AN.p2_contingency(vcat(fill(:chaotic, a + b), fill(:regular, c + d)), vcat(fill(:chaotic, a), fill(:regular, b), fill(:chaotic, c), fill(:regular, d)))
    @test (t = tab(10, 0, 0, 10); (t.cc, t.cr, t.rc, t.rr, t.n, t.A, t.kappa) == (10, 0, 0, 10, 20, 1.0, 1.0))             # perfect agreement: kappa = 1
    @test tab(25, 25, 25, 25).kappa == 0.0 && tab(25, 25, 25, 25).A == 0.5                                               # p_o = p_e = 0.5: kappa = 0
    @test tab(10, 0, 0, 0).kappa == 0.0 && tab(0, 0, 0, 10).kappa == 0.0                                               # one verdict only: p_e = 1, kappa = 0 by definition
    @test tab(0, 7, 7, 0).kappa == -1.0 && tab(0, 7, 7, 0).A == 0.0                                                      # a = d = 0, b = c: kappa = -1
    # [[20, 5], [10, 15]]: N = 50, p_o = 35 / 50 = 0.7, p_e = (25 x 30 + 25 x 20) / 2500 = 0.5, kappa = 0.2 / 0.5 = 0.4
    t = tab(20, 5, 10, 15)
    @test t.n == 50 && t.A == 0.7 && isapprox(t.kappa, 0.4; atol = 1e-15)
    # orbits undecided (or without the item) in either indicator are left out; N = 0 is undefined
    t = AN.p2_contingency([:chaotic, :undecided, :na, :regular, :chaotic], [:chaotic, :chaotic, :regular, :undecided, :regular])
    @test t.n == 2 && t.cc == 1 && t.cr == 1 && t.A == 0.5
    t = AN.p2_contingency(Symbol[:undecided, :na], Symbol[:chaotic, :chaotic])
    @test t.n == 0 && t.A === nothing && t.kappa === nothing
end

@testset "D89i unanimity failure U, A_min and the mean kappa: a hand table" begin
    # three indicators (LE, SALI, FLI), five orbits:        LE SALI FLI
    #   o1 c c c   o2 c c r   o3 r r r   o4 c u r (dropped: SALI undecided)   o5 c r r
    V = Dict(:LE => [:chaotic, :chaotic, :regular, :chaotic, :chaotic], :SALI => [:chaotic, :chaotic, :regular, :undecided, :regular],
             :FLI => [:chaotic, :regular, :regular, :regular, :regular])
    u = AN.p2_unanimity(V, [:LE, :SALI, :FLI])
    # used: o1, o2, o3, o5 = 4; not all equal: o2, o5 -> U = 2 / 4; one orbit dropped
    @test u.used == 4 && u.dropped == 1 && u.U == 0.5
    # (LE, SALI): o1 cc, o2 cc, o3 rr, o5 c-r (o4 out): N = 4, agree 3 -> 0.75, kappa 0.5 (p_e = (3 x 2 + 1 x 2) / 16 = 0.5)
    # (LE, FLI): o1 cc, o2 c-r, o3 rr, o4 c-r, o5 c-r: N = 5, agree 2 -> 0.4, p_e = (4 x 1 + 1 x 4) / 25 = 0.32, kappa = 0.08 / 0.68
    # (SALI, FLI): o1 cc, o2 c-r, o3 rr, o5 r-r: N = 4, agree 3 -> 0.75, p_e = (2 x 1 + 2 x 3) / 16 = 0.5, kappa = 0.5
    @test u.A[(:LE, :SALI)] == 0.75 && u.A[(:LE, :FLI)] == 0.4 && u.A[(:SALI, :FLI)] == 0.75 && length(u.A) == 3 && u.A_min == 0.4
    @test isapprox(u.kappa[(:LE, :SALI)], 0.5; atol = 1e-15) && isapprox(u.kappa[(:LE, :FLI)], 0.08 / 0.68; rtol = 1e-14) && isapprox(u.kappa[(:SALI, :FLI)], 0.5; atol = 1e-15)
    @test isapprox(u.kappa_mean, (0.5 + 0.08 / 0.68 + 0.5) / 3; rtol = 1e-14)
    # U is 0 if every orbit has all verdicts equal; an indicator without item (:na) drops the orbit from U
    @test AN.p2_unanimity(Dict(:LE => [:chaotic, :regular], :SALI => [:chaotic, :regular]), [:LE, :SALI]).U == 0.0
    w = AN.p2_unanimity(Dict(:LE => [:chaotic, :regular], :SALI => [:na, :regular]), [:LE, :SALI])
    @test w.used == 1 && w.dropped == 1 && w.U == 0.0
    # the number of pairs: 21 for the seven indicators, 28 for eight
    Z(n) = Dict(X => [:chaotic, :regular, :chaotic] for X in P2_X8[1:n])
    @test length(AN.p2_unanimity(Z(7), P2_X7).A) == 21 && length(AN.p2_unanimity(Z(8), P2_X8).A) == 28
end

@testset "D89i 1000 random verdict tables against the independent implementation" begin
    rng = StableRNG(p2_seed('i', 1))                    # testset D89i: sub 1; the seed is recorded in the plan
    pr = (0.45, 0.35, 0.12, 0.08)                      # chaotic, regular, undecided, no item
    draw() = (u = rand(rng); u < pr[1] ? :chaotic : u < pr[1] + pr[2] ? :regular : u < pr[1] + pr[2] + pr[3] ? :undecided : :na)
    bad = String[]
    for t in 1:1000
        I = t % 2 == 0 ? P2_X8 : P2_X7
        n = rand(rng, 4:40)
        V = Dict(X => [draw() for _ in 1:n] for X in I)
        a, b = AN.p2_unanimity(V, I), p2ref_unanimity(V, I)
        ok = a.used == b.used && a.dropped == b.dropped && p2_close(a.U, b.U) && p2_close(a.A_min, b.A_min) && p2_close(a.kappa_mean, b.kappa_mean) &&
             all(pr -> p2_close(a.A[pr], b.A[pr]) && p2_close(a.kappa[pr], b.kappa[pr]), keys(b.A)) && length(a.A) == length(b.A)
        ok || push!(bad, "table $t")
        x, y = V[I[1]], V[I[2]]
        c, d = AN.p2_contingency(x, y), p2ref_table(x, y)
        (c.n == d.n && p2_close(c.A, d.A) && p2_close(c.kappa, d.kappa)) || push!(bad, "pair $t")
    end
    @test isempty(bad)
end

@testset "D89i Spearman correlation: midranks, ties, undefined cases, 1000 random vectors" begin
    # ties get the mean of their ranks: x = (1, 2, 2, 3), y = (1, 3, 2, 4): ranks (1, 2.5, 2.5, 4) and (1, 3, 2, 4); Pearson = 4.5 / sqrt(4.5 x 5) = sqrt(0.9)
    @test isapprox(AN.p2_spearman([1.0, 2.0, 2.0, 3.0], [1.0, 3.0, 2.0, 4.0]), sqrt(0.9); rtol = 1e-14)
    @test AN.p2_spearman([1.0, 2, 3, 4], [10.0, 20, 30, 40]) == 1.0 && AN.p2_spearman([1.0, 2, 3, 4], [4.0, 3, 2, 1]) == -1.0
    @test AN.p2_spearman([1.0, 1, 1], [1.0, 2, 3]) === nothing && AN.p2_spearman([1.0, 2, 3], [5.0, 5, 5]) === nothing        # a constant rank vector: undefined
    @test AN.p2_spearman([1.0, 2], [2.0, 1]) === nothing && AN.p2_spearman(Float64[], Float64[]) === nothing                  # fewer than 3 orbits: undefined
    @test isapprox(AN.p2_spearman([-Inf, 1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0]), 1.0; atol = 1e-15)                            # minus infinity keeps its rank
    rng = StableRNG(p2_seed('i', 2))
    bad = 0
    for t in 1:1000
        n = rand(rng, 2:30)
        x = Float64.(rand(rng, 1:6, n)); y = Float64.(rand(rng, 1:6, n))
        t % 7 == 0 && (x[rand(rng, 1:n)] = -Inf)
        bad += !p2_close(AN.p2_spearman(x, y), p2ref_spearman(x, y); rtol = 1e-11)
    end
    @test bad == 0
end

@testset "D89j AUC: hand cases (O7 b), ties, minus infinity, 1000 random score vectors, the min_set rule" begin
    @test AN.p2_auc([5.0, 6.0], [1.0, 2.0]) == 1.0 && AN.p2_auc([1.0, 2.0], [3.0, 4.0]) == 0.0 && AN.p2_auc([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == 0.5
    @test AN.p2_auc([1.0, 2.0, 3.0], [2.0]) == 0.5                                      # (1 + 0.5 + 0) ... C = {1, 2, 3} against R = {2}: pairs 0, 1/2, 1 -> 0.5
    @test AN.p2_auc([2.0, 3.0], [2.0, 1.0]) == (0.5 + 1 + 1 + 1) / 4                    # one tied pair counts one half
    @test AN.p2_auc([-Inf, 1.0], [-Inf]) == 0.75                                         # (1/2 + 1) / 2
    rng = StableRNG(p2_seed('j', 1))
    bad = 0; neg = 0
    for t in 1:1000
        c = Float64.(rand(rng, 1:8, rand(rng, 1:25))); r = Float64.(rand(rng, 1:8, rand(rng, 1:25)))
        t % 5 == 0 && (c[1] = -Inf; r[1] = -Inf)
        a = AN.p2_auc(c, r)
        bad += !p2_close(a, p2ref_auc(c, r))
        neg += !isapprox(a + AN.p2_auc(-c, -r), 1.0; atol = 1e-15)                      # the AUC of a score and of its negation sum to 1
    end
    @test bad == 0 && neg == 0
    # computed only if both classes have at least min_set orbits: 19 and 20
    @test AN.p2_auc(fill(1.0, 19), fill(0.0, 20); min_set = 20) === nothing && AN.p2_auc(fill(1.0, 20), fill(0.0, 19); min_set = 20) === nothing
    @test AN.p2_auc(fill(1.0, 20), fill(0.0, 20); min_set = 20) == 1.0
end

@testset "D89j regular range, leaving budget and the share below the regular median (M10, P25)" begin
    p = p2_params()                                                                      # min_set 5, percentiles 5 and 95, budgets 40 x 2^k
    # type 7 on 1..20: h = 19 x 0.05 + 1 = 1.95 -> 1.95; h = 19 x 0.95 + 1 = 19.05 -> 19.05; median 10.5
    rr = AN.p2_regular_range(collect(1.0:20.0), p)
    @test isapprox(rr.lo, 1.95; atol = 1e-12) && isapprox(rr.hi, 19.05; atol = 1e-12) && rr.median == 10.5 && rr.n == 20
    # exactly min_set orbits: computed; fewer: no range
    @test AN.p2_regular_range([1.0, 2.0, 3.0, 4.0, 5.0], p) !== nothing && AN.p2_regular_range([1.0, 2.0, 3.0, 4.0], p) === nothing
    # leaving budget: the smallest k from which the score is strictly above the upper percentile at every later index up to 9
    up = fill(3.0, 10)
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[1, 2, 5, 6, 7, 8, 9, 10, 11, 12], up, p) == 160                      # from k = 2: b_2 = 160
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[1, 2, 5, 6, 7, 3, 9, 10, 11, 12], up, p) == 2560                     # 3.0 is not strictly above: from k = 6
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[1, 1, 1, 1, 1, 1, 1, 1, 1, 12], up, p) == 20480                      # only the last index
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[9, 9, 9, 9, 9, 9, 9, 9, 9, 2], up, p) === nothing                    # below at k = 9: undefined
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[9, 9, 9, 9, 9, 9, 9, 9, 9, nothing], up, p) === nothing              # no available score at k = 9
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[9, 9, nothing, 9, 9, 9, 9, 9, 9, 9], up, p) == 320                   # an unavailable score breaks the run: from k = 3, b_3 = 320
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[11, 12, 13, 14, 15, 16, 17, 18, 19, 20], collect(1.0:10.0), p) == 40                 # above at every index: from k = 0
    @test AN.p2_leaving_budget(Union{Nothing, Float64}[4, 4, 4, 4, 4, 4, 4, 4, 4, 4], collect(1.0:10.0), p) === nothing                    # 4 > 4 fails from k = 3 on
    # the share of the sticky orbits whose score at k = 9 is below the median of the regular class: scores 1, 2, 3, 10 against 5 -> 3 of 4
    @test AN.p2_share_below([1.0, 2.0, 3.0, 10.0], 5.0) == 0.75 && AN.p2_share_below([6.0, 7.0], 5.0) == 0.0 && AN.p2_share_below([5.0], 5.0) == 0.0
end
