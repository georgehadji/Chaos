# D93g rule evaluator, D93h search space, D93i cost and matched budgets. Spec: analysis-spec.md, Amendment D93, E51 to E53, criteria D93g, D93h, D93i.
# The oracle is the independent reference of an_p3_ref.jl (written from the amendment text) and hand tables. The functions are pure: no store is read.
# Seeds: p3_seed('g', 1) (random verdict tables), p3_seed('g', 2) (random elements), p3_seed('h', 1) (random keys), p3_seed('i', 1) (random cost tables).

const P3_VS = (:chaotic, :regular, :undecided)
p3_rand_V(rng, names) = Dict{Symbol, Vector{Symbol}}(X => Symbol[P3_VS[rand(rng, 1:3)] for _ in 1:10] for X in names)
p3_th_nt(t) = (V = t.V, qC = t.qC, qR = t.qR, w = t.w, v = t.v)

@testset "D93g evidence, stability and the first applicable line (hand cases)" begin
    # the three voters A, B, C of the plan: A = rrrrrrcccc, B = rrrrrrrrcc, C = uuuuuurrrr (index 0 to 9)
    V = Dict(:A => p2_sym("rrrrrrcccc"), :B => p2_sym("rrrrrrrrcc"), :C => p2_sym("uuuuuurrrr"))
    th(; qC = 1, qR = 2, w = 0, v = 1) = (V = [:A, :B, :C], qC = qC, qR = qR, w = w, v = v)
    # stability at depth w: decided and equal to the verdicts at b / 2, ..., b / 2^w; w = 0 means decided
    @test AN.p3_stable(V[:A], 8, 2) && !AN.p3_stable(V[:A], 7, 2) && AN.p3_stable(V[:A], 7, 1) && AN.p3_stable(V[:A], 6, 0) && !AN.p3_stable(V[:C], 5, 0) && !AN.p3_stable(V[:C], 6, 1)
    @test !AN.p3_stable(V[:B], 8, 1) && AN.p3_stable(V[:B], 9, 1) && AN.p3_stable(V[:B], 7, 2) && !AN.p3_stable(V[:B], 9, 2)
    # n_C, n_R, n_Cany (an unstable chaotic vote counts in n_Cany only)
    @test AN.p3_evidence(V, 6, th(w = 0)) == (nC = 1, nR = 2, nCany = 1) && AN.p3_evidence(V, 6, th(w = 1)) == (nC = 0, nR = 1, nCany = 1)
    @test AN.p3_evidence(V, 7, th(w = 1)) == (nC = 1, nR = 2, nCany = 1) && AN.p3_evidence(V, 8, th(w = 2)) == (nC = 1, nR = 1, nCany = 2)
    # the lines: chaotic-supported first (precedence), then regular-supported on n_Cany <= v and n_R >= q_R, else unresolved
    @test AN.p3_rule_a(V, 6, th()) === :chaotic_supported && AN.p3_rule_a(V, 6, th(qC = 2)) === :regular_supported && AN.p3_rule_a(V, 6, th(qC = 2, v = 0)) === :unresolved
    @test AN.p3_rule_a(V, 6, th(qC = 2, qR = 3)) === :unresolved && AN.p3_rule_a(V, 6, th(qC = 2, w = 1, v = 1)) === :unresolved && AN.p3_rule_a(V, 7, th(qC = 2, w = 1)) === :regular_supported
    # chaotic precedence: n_C >= q_C wins although the regular line would also hold
    W = Dict(:A => p2_sym("cccccccccc"), :B => p2_sym("rrrrrrrrrr"), :C => p2_sym("rrrrrrrrrr"))
    @test AN.p3_rule_a(W, 9, (V = [:A, :B, :C], qC = 1, qR = 2, w = 0, v = 1)) === :chaotic_supported && AN.p3_rule_a(W, 9, (V = [:A, :B, :C], qC = 2, qR = 2, w = 0, v = 1)) === :regular_supported
    @test AN.p3_rule_a(W, 9, (V = [:A, :B, :C], qC = 2, qR = 2, w = 0, v = 0)) === :unresolved
    # a budget b_k / 2^w that is not an element of B is an error of the caller
    @test an_rejects(() -> AN.p3_stable(V[:A], 1, 2), ["budget"]) && an_rejects(() -> AN.p3_rule_a(V, 1, th(w = 2)), ["budget"]) && an_rejects(() -> AN.p3_evidence(V, 0, th(w = 1)), ["budget"])
    @test AN.p3_stable(V[:A], 2, 2) && AN.p3_stable(V[:A], 0, 0)
end

@testset "D93g R_a on all 3^9 verdict histories of three indicators at depth 2 (exhaustive), all parameter values of V = SALI, FLI, MEGNO" begin
    names = [:SALI, :FLI, :MEGNO]
    els = [(V = names, qC = qC, qR = qR, w = w, v = v) for qC in 1:3, qR in 2:3, w in 0:2, v in 0:1][:]
    @test length(els) == 36
    syms = P3_VS
    bad = 0; total = 0
    for h in 0:(3^9 - 1)
        d = h; hist = Int[]
        for _ in 1:9; push!(hist, d % 3 + 1); d ÷= 3; end
        V = Dict{Symbol, Vector{Symbol}}(X => vcat(Symbol[syms[hist[3 * (j - 1) + i]] for i in 1:3], fill(:undecided, 7)) for (j, X) in enumerate(names))
        for t in els
            total += 1
            AN.p3_rule_a(V, 2, t) === p3ref_out(V, 2, t) || (bad += 1)
        end
    end
    @test total == 36 * 3^9 && bad == 0
end

@testset "D93g R_a on every element with |V| = 3 (symmetric in V: 200 random histories each) and on 1000 random elements with 10^5 random verdict tables" begin
    q = p3_params()
    th = p3ref_theta(q)
    rng = StableRNG(p3_seed('g', 1))
    three = [t for t in th if length(t.V) == 3]
    @test length(three) == 720
    bad = 0
    for t in three, _ in 1:200
        V = p3_rand_V(rng, t.V)
        k = 2 + rand(rng, 0:7)
        AN.p3_rule_a(V, k, p3_th_nt(t)) === p3ref_out(V, k, t) || (bad += 1)
    end
    @test bad == 0
    rng2 = StableRNG(p3_seed('g', 2))
    picks = [th[rand(rng2, 1:length(th))] for _ in 1:1000]
    @test any(t -> length(t.V) == 6, picks) && any(t -> length(t.V) == 4, picks) && any(t -> t.w == 2 && t.v == 1, picks)
    bad = 0; fired = Dict(:chaotic_supported => 0, :regular_supported => 0, :unresolved => 0)
    for t in picks, _ in 1:100
        V = Dict{Symbol, Vector{Symbol}}(X => (rand(rng2) < 0.5 ? Symbol[P3_VS[rand(rng2, 1:3)] for _ in 1:10] : fill(P3_VS[rand(rng2, 1:2)], 10)) for X in Symbol.(q.voter_bits))
        k = 2 + rand(rng2, 0:7)
        o = AN.p3_rule_a(V, k, p3_th_nt(t))
        fired[o] += 1
        o === p3ref_out(V, k, t) || (bad += 1)
    end
    @test bad == 0 && all(>(1000), values(fired))                      # all three outputs occur often: the comparison is not vacuous
end

@testset "D93g R_ab: the override fires iff a decided Float64 and Double64 pair differs; undefined without Double64 items" begin
    V = Dict(:SALI => p2_sym("rrrrrrcccc"), :FLI => p2_sym("rrrrrrcccc"), :MEGNO => p2_sym("rrrrrrcccc"))
    t = (V = [:SALI, :FLI, :MEGNO], qC = 2, qR = 2, w = 0, v = 0)
    same = Dict(X => copy(V[X]) for X in keys(V))
    @test AN.p3_rule_ab(V, same, 6, t) === :chaotic_supported === AN.p3_rule_a(V, 6, t)
    d1 = Dict(X => copy(V[X]) for X in keys(V)); d1[:FLI][7] = :regular                 # a decided pair that differs at b_6: unresolved whatever R_a says
    @test AN.p3_rule_ab(V, d1, 6, t) === :unresolved && AN.p3_rule_a(V, 6, t) === :chaotic_supported
    d2 = Dict(X => copy(V[X]) for X in keys(V)); d2[:FLI][7] = :undecided                # an undecided Double64 verdict is not a difference
    @test AN.p3_rule_ab(V, d2, 6, t) === :chaotic_supported
    d3 = Dict(X => copy(V[X]) for X in keys(V)); d3[:FLI][6] = :regular                  # a difference at another budget does not count
    @test AN.p3_rule_ab(V, d3, 6, t) === :chaotic_supported
    Vu = Dict(X => copy(V[X]) for X in keys(V)); Vu[:SALI][7] = :undecided
    @test AN.p3_rule_ab(Vu, d1, 6, t) === :unresolved && AN.p3_rule_ab(Vu, same, 6, t) === AN.p3_rule_a(Vu, 6, t)       # the override acts on R_a's result of the Float64 verdicts only
    # the stability of evidence (a) uses the Float64 verdicts only: a Double64 history that is unstable does not change it
    d4 = Dict(X => p2_sym("rcrcrcrrcc") for X in keys(V)); d4[:SALI][7] = :chaotic; d4[:FLI][7] = :chaotic; d4[:MEGNO][7] = :chaotic
    @test AN.p3_rule_ab(V, d4, 6, merge(t, (w = 2,))) === AN.p3_rule_a(V, 6, merge(t, (w = 2,)))
    # an orbit without Double64 items: R_ab is not defined, never silently R_a
    @test an_rejects(() -> AN.p3_rule_ab(V, nothing, 6, t), ["double64"])
    @test an_rejects(() -> AN.p3_rule_ab(V, Dict(:SALI => V[:SALI]), 6, t), ["double64"])
    # against the reference, on random pairs of histories
    rng = StableRNG(p3_seed('g', 3))
    q = p3_params(); th = p3ref_theta(q)
    bad = 0
    for _ in 1:3000
        t = th[rand(rng, 1:length(th))]
        Vf, Vd = p3_rand_V(rng, Symbol.(q.voter_bits)), p3_rand_V(rng, Symbol.(q.voter_bits))
        k = 2 + rand(rng, 0:7)
        AN.p3_rule_ab(Vf, Vd, k, p3_th_nt(t)) === p3ref_out_ab(Vf, Vd, k, t) || (bad += 1)
    end
    @test bad == 0
end

@testset "D93h the search space: 2088 distinct valid elements, the order of E52, the closed form for other keys" begin
    q = p3_params()
    th = AN.p3_enumerate(q)
    ref = p3ref_theta(q)
    @test length(th) == 2088 == AN.p3_space_size(q) && [t.index for t in th] == 1:2088
    @test [count(t -> length(t.V) == k, th) for k in 3:6] == [720, 720, 540, 108]
    @test all(i -> (th[i].index, th[i].V, th[i].mask, th[i].qC, th[i].qR, th[i].w, th[i].v) == (ref[i].index, ref[i].V, ref[i].mask, ref[i].qC, ref[i].qR, ref[i].w, ref[i].v), 1:2088)
    @test allunique((t.V, t.qC, t.qR, t.w, t.v) for t in th) && all(t -> 1 <= t.qC <= length(t.V) && length(t.V) ÷ 2 + 1 <= t.qR <= length(t.V) && t.w in 0:2 && t.v in 0:1 && length(t.V) >= 3, th)
    # the bit of an indicator of I_6 in the bitmask of V: 0 SALI, 1 FLI, 2 MEGNO, 3 E1, 4 E2, 5 WBA; V lists its indicators in the bit order
    @test all(t -> t.mask == sum(1 << (findfirst(==(X), Symbol.(q.voter_bits)) - 1) for X in t.V) && issorted([findfirst(==(X), Symbol.(q.voter_bits)) for X in t.V]), th)
    # the order: by |V| ascending, then the bitmask, then q_C, q_R, w, v, each ascending (the order of C6 is a different one)
    key(t) = (length(t.V), t.mask, t.qC, t.qR, t.w, t.v)
    @test issorted(th; by = key) && (th[1].V, th[1].qC, th[1].qR, th[1].w, th[1].v) == ([:SALI, :FLI, :MEGNO], 1, 2, 0, 0) && th[2].v == 1 && th[3].w == 1 && th[end].V == Symbol.(q.voter_bits) && th[end].qC == 6 && th[end].qR == 6
    # the closed form equals the enumeration for other key values: w_range of 2 and 4 elements, v_range of 1 element, min_voters 1 and 4, five voters
    rng = StableRNG(p3_seed('h', 1))
    five = q.voter_bits[[1, 2, 3, 5, 6]]
    cases = [merge(q, (; w_range = [0, 1])), merge(q, (; w_range = [0, 1, 2, 3])), merge(q, (; v_range = [0])), merge(q, (; min_voters = 1)), merge(q, (; min_voters = 4)),
             merge(q, (; voter_bits = five, indicators_i7 = vcat("LE", five))), merge(q, (; voter_bits = five, indicators_i7 = vcat("LE", five), min_voters = 1, w_range = [0, 2], v_range = [0, 1, 2])),
             merge(q, (; min_voters = 6)), merge(q, (; min_voters = 2, v_range = [1], w_range = [1]))]
    sizes = Int[]
    for c in cases
        a, r = AN.p3_enumerate(c), p3ref_theta(c)
        push!(sizes, length(a))
        @test length(a) == AN.p3_space_size(c) == p3ref_count(c) && all(i -> (a[i].V, a[i].qC, a[i].qR, a[i].w, a[i].v) == (r[i].V, r[i].qC, r[i].qR, r[i].w, r[i].v), eachindex(a)) &&
              allunique((t.V, t.qC, t.qR, t.w, t.v) for t in a)
    end
    @test sizes[1] == 1392 && sizes[2] == 2784 && sizes[3] == 1044 && sizes[8] == 18 * 6 && all(s -> s > 0, sizes)
    # for five voters the sizes by |V|: C(5, k) k (k - floor(k / 2)) x 3 x 2 for k = 3, 4, 5
    @test sizes[6] == (10 * 3 * 2 + 5 * 4 * 2 + 1 * 5 * 3) * 3 * 2 == 690
end

@testset "D93i cost of a rule, matched budgets (ceiling rule, exact ties), round-down, caps, Double64 costs" begin
    q = p3_params()
    # exact decimals: the rational of the decimal text of a key (2076 / 10, not the nearest binary fraction)
    @test AN.p3_rat(207.6) == 2076 // 10 && AN.p3_rat(0.1) == 1 // 10 && AN.p3_rat(1e-4) == 1 // 10000 && AN.p3_rat(5) == 5 && AN.p3_rat(16850.4) == 168504 // 10
    # hand table, V = SALI, FLI, MEGNO (the spec tg counted once): C_R = 2560 x 2658.0 = 6804480
    cr(V, b = 2560) = AN.p3_rule_cost(q, V, b)
    @test cr([:SALI, :FLI, :MEGNO]) == 6804480 && cr([:SALI]) == 6804480 && cr([:FLI, :MEGNO]) == 6804480 && cr([:SALI, :FLI, :MEGNO], 5120) == 2 * 6804480
    @test cr([:SALI, :E2]) == 2560 * (26580 + 18105) // 10 && cr([:SALI, :FLI, :E2, :WBA]) == 2560 * (26580 + 18105 + 30006) // 10 && cr([:E1]) == 2560 * 168504 // 10
    m = AN.p3_matched(q, [:SALI, :FLI, :MEGNO])
    r = p3ref_matched(q, [:SALI, :FLI, :MEGNO])
    @test m.CR == 6804480 && m.X[:SALI].b == 2560 && m.X[:FLI].b == 2560 && m.X[:MEGNO].b == 2560 && m.X[:E1].b == 640 && m.X[:E2].b == 5120 && m.X[:WBA].b == 2560
    @test m.X[:E1].round_down == 320 && m.X[:E2].round_down == 2560 && m.X[:WBA].round_down == 1280 && m.X[:SALI].round_down == 2560 && !m.X[:E1].round_down_none
    @test all(X -> !m.X[X].capped, (:SALI, :FLI, :MEGNO, :E1, :E2, :WBA)) && all(X -> (m.X[X].b, m.X[X].capped, m.X[X].round_down) == (r.X[X].b, r.X[X].capped, r.X[X].round_down), (:SALI, :E1, :E2, :WBA, :LE))
    # LE: never matched (capped at the ceiling rule), read at N_c, flagged cost-unmatched (read budget x cost below C_R)
    @test m.X[:LE].capped && m.X[:LE].b === nothing && m.le_read == q.N_c && m.le_unmatched && m.le_read * AN.p3_rat(q.cost_l1) < m.CR
    # a V that contains E1: the comparator E2 is capped (no element of B reaches C_R), F1 of the cost cap rejects it; the cost cap of F1 is min cost x N_c = 1810.5 x N_c
    me = AN.p3_matched(q, [:SALI, :FLI, :E1])
    @test me.X[:E2].capped && me.X[:E2].b === nothing && !me.X[:E1].capped && AN.p3_cost_cap(q) == 18105 // 10 * 20480 && me.CR > AN.p3_cost_cap(q) && m.CR <= AN.p3_cost_cap(q)
    th = AN.p3_enumerate(q)
    @test count(t -> AN.p3_rule_cost(q, t.V, q.b_primary) <= AN.p3_cost_cap(q), th) == 690 && all(t -> (:E1 in t.V) == (AN.p3_rule_cost(q, t.V, q.b_primary) > AN.p3_cost_cap(q)), th)
    # round-down: the largest element with b x cost <= C_R, or none (with the flag)
    qh = merge(q, (; cost_pa1 = 1.0e7))
    mh = AN.p3_matched(qh, [:SALI, :FLI, :MEGNO])
    @test mh.X[:E2].round_down === nothing && mh.X[:E2].round_down_none && mh.X[:E2].b == 40 && !mh.X[:E2].capped
    # the ceiling rule at an exact tie selects that b; binary64 would miss it (found by search: a pair of decimal costs whose float product is below the exact one)
    found = nothing
    for ctg in (1:99) ./ 10, cx in (1:99) ./ 100, kk in 0:9
        b = 10 * (kk + 1)
        (p3ref_rat(b) * p3ref_rat(cx) == p3ref_rat(10) * p3ref_rat(ctg)) && (b * cx < 10 * ctg) && (found = (ctg, cx, b); break)
        found === nothing || break
    end
    @test found !== nothing
    if found !== nothing
        ctg, cx, b = found
        qt = merge(q, (; budgets = collect(10:10:100), b_primary = 10, N_c = 100, cost_tg = ctg, cost_pa1 = cx))
        mt = AN.p3_matched(qt, [:SALI])
        @test mt.CR == 10 * AN.p3_rat(ctg) && mt.X[:E2].b == b && mt.X[:E2].round_down == b             # b x cost == C_R exactly: the ceiling rule and the round-down rule both select it
    end
    # the cost-unmatched flag is on for every feasible element (the amendment says none carries it), the capped flag of LE too
    @test all(t -> (mm = AN.p3_matched(q, t.V); mm.X[:LE].capped && mm.le_unmatched), [t for t in th if !(:E1 in t.V)])
    # Double64 costs (reported): pair cost b_P x sum(c + r c) over the distinct specs, Float64 vote cost b_P x sum c, ratio of the second to the first
    d = AN.p3_d64_costs(q, [:SALI, :FLI, :E2])
    @test d.float == cr([:SALI, :FLI, :E2]) && d.pair == 2560 * (AN.p3_rat(2658.0) * (1 + AN.p3_rat(10.5)) + AN.p3_rat(1810.5) * (1 + AN.p3_rat(7.28))) && d.ratio == d.float // d.pair
    d1 = AN.p3_d64_costs(q, [:SALI, :FLI, :MEGNO])
    @test d1.ratio == 1 // (1 + 105 // 10)                                                                    # one spec: 1 / (1 + r), the amendment's remark
    # random cost tables against the reference: the matched budgets and the round-down budgets of random elements
    rng = StableRNG(p3_seed('i', 1))
    for _ in 1:30
        qq = merge(q, (; cost_tg = round(rand(rng) * 5000; digits = 1), cost_pa1 = round(rand(rng) * 5000; digits = 1), cost_bk1 = round(rand(rng) * 5000; digits = 1),
                       cost_fq1 = round(rand(rng) * 20000; digits = 1)))
        t = th[rand(rng, 1:length(th))]
        a, b = AN.p3_matched(qq, t.V), p3ref_matched(qq, t.V)
        @test a.CR == b.CR && all(X -> (a.X[X].b, a.X[X].capped, a.X[X].round_down, a.X[X].round_down_none) == (b.X[X].b, b.X[X].capped, b.X[X].round_down, b.X[X].round_down_none), keys(b.X)) &&
              a.le_read == b.le_read && a.le_unmatched == b.le_unmatched
    end
end
