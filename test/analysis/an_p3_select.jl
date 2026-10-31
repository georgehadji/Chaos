# D93j calibration C3 to C8 on prescribed tables (O18 a to e): theta* by construction, F1 to F6 each excluding the prescribed elements, the order of C6 key by key, exact
# fractions and the exact tie at delta_tol, "no feasible rule", C3 to C6 on the designed 36 orbits against the independent reference. Spec: analysis-spec.md, Amendment D93,
# E54, I32, I35, O18, criterion D93j (S). Seeds: p3_seed('j', 1) (random candidates for C6).

@testset "D93j (a) theta* is the element built by construction; F1 to F6 each exclude exactly the prescribed elements" begin
    q = p3_toy_q()
    rows = p3_toy_rows(P3_TOY_A)
    sel = AN.p3_select(rows, q)
    @test sel.outcome === :selected && length(sel.thetas) == 55 == AN.p3_space_size(q) && sel.feasible == [11] && sel.star == 11 && sel.tol == [11]
    t = sel.thetas[11]
    @test (t.V, t.qC, t.qR, t.w, t.v) == ([:SALI, :FLI, :MEGNO, :E2], 2, 3, 0, 0) && t.mask == 1 + 2 + 4 + 16 && sel.J_min == sel.J[11]
    # F1: exactly the 8 elements of V = SALI, FLI, MEGNO, E2 satisfy the cost cap alone; every element with E1 fails first at F1
    @test sel.alone[1] == 8 && all(i -> sel.first_fail[i] == 1, [i for (i, t) in enumerate(sel.thetas) if :E1 in t.V]) && all(i -> sel.first_fail[i] != 1, [i for (i, t) in enumerate(sel.thetas) if !(:E1 in t.V)])
    # the 8 elements of V (indices 9 to 16: q_C 1 to 4 and q_R 3, 4): which items hold, by the construction above
    holds = Dict((1, 3) => (true, false, true, false, true, true), (1, 4) => (true, true, true, false, true, false), (2, 3) => (true, true, true, true, true, true),
                 (2, 4) => (true, true, true, true, false, false), (3, 3) => (true, true, false, true, true, true), (3, 4) => (true, true, false, true, false, false),
                 (4, 3) => (true, false, false, true, false, false), (4, 4) => (true, true, false, true, false, false))
    for (j, i) in enumerate(9:16)
        el = sel.thetas[i]
        @test el.V == [:SALI, :FLI, :MEGNO, :E2] && sel.flags[i] == holds[(el.qC, el.qR)]
    end
    # exact fractions (I32): e, shares and the gate count of E*; the stray votes are the contradictions of q_C = 1
    s = sel.stats[11]
    @test s.e isa Rational && s.u_cal isa Rational && s.fr_ts isa Rational && s.e == 1 // 3 && s.u_cal == 1 // 6 && s.u_s2 == 0 && s.fr_ts == 1 // 3 && s.fr_w2 == 5 // 7 && s.dec_ts == 1 && s.dec_w2 == 1 && s.g == 0
    @test s.decided == 15 && s.contradicting == 5 && sel.stats[9].e == 4 // 9 && sel.stats[9].g == 1 && sel.stats[9].contradicting == 8
    @test all(i -> sel.stats[i].cost == AN.p3_rule_cost(q, sel.thetas[i].V, q.b_primary), eachindex(sel.thetas)) && sel.stats[11].cost == 2560 * (26580 + 18105) // 10
    # the whole selection equals the independent reference, element by element
    @test same_sel(sel, p3ref_select(rows, q))
    # the outputs of R_a on the orbits of the table are what the construction says: the stray orbits with one chaotic vote
    @test AN.p3_rule_a(rows[16].V, 6, sel.thetas[9]) === :chaotic_supported && AN.p3_rule_a(rows[16].V, 6, sel.thetas[11]) === :unresolved && AN.p3_rule_a(rows[8].V, 6, sel.thetas[11]) === :regular_supported
end

@testset "D93j the ties of C6: the exact tie at delta_tol, binary64 would choose another element; J, Theta_tol and the order are exact" begin
    # T_s = W2 = 10 orbits: two chaotic, seven with four regular votes, one with three: q_R = 3 calls 8 of 10 regular (J = 8/10), q_R = 4 calls 7 of 10 (J = 7/10);
    # delta_tol = 0.1: J_min + delta_tol = 7/10 + 1/10 = 8/10 exactly, so both are in Theta_tol; binary64 gives 0.7 + 0.1 = 0.7999999999999999 < 0.8 and would drop q_R = 3
    q = p3_toy_q(; e_cal = 1.0, a_cal = 1.0, delta_tol = 0.1, c_min = 0.0)
    spec = vcat([(4, 0, (S2 = true,)) for _ in 1:4], [(4, 0, (Ts = true,)) for _ in 1:2], [(0, 4, (Ts = true,)) for _ in 1:7], [(0, 3, (Ts = true,))], [(0, 4, (rho = :R,)) for _ in 1:4])
    rows = p3_toy_rows(spec)
    sel = AN.p3_select(rows, q)
    @test 0.7 + 0.1 < 0.8 && sel.J_min == 7 // 10 && sel.J[10] == 7 // 10 && sel.J[9] == 8 // 10 && 9 in sel.tol && 10 in sel.tol
    @test all(i -> sel.J[i] <= sel.J_min + 1 // 10, sel.tol) && all(i -> sel.J[i] > sel.J_min + 1 // 10, setdiff(sel.feasible, sel.tol))
    # the first of the order: lowest C_R, |V|, w, v, bitmask (all equal here), then q_C = 1, then q_R = 3
    @test sel.star == 9 && (sel.thetas[9].qC, sel.thetas[9].qR) == (1, 3) && same_sel(sel, p3ref_select(rows, q))
    # one orbit more in R_W moves J by 1/40: W2 = T_s (10, all chaotic) + R_W (10); q_R = 4 calls 2 of the R_W orbits regular (J = 2/40 = J_min), q_R = 3 calls 4 (J = 4/40 = J_min + delta_tol
    # exactly with delta_tol = 0.05: in Theta_tol); with a fifth regular call (J = 5/40) it is out and q_R = 4 is chosen
    q2 = p3_toy_q(; e_cal = 1.0, a_cal = 1.0, delta_tol = 0.05, c_min = 0.0)
    tbl(nthree, nchaotic) = vcat([(4, 0, (S2 = true,)) for _ in 1:4], [(4, 0, (Ts = true,)) for _ in 1:10], [(0, 4, (W2 = true,)) for _ in 1:2], [(0, 3, (W2 = true,)) for _ in 1:nthree],
                                 [(4, 0, (W2 = true,)) for _ in 1:nchaotic], [(0, 4, (rho = :R,)) for _ in 1:3])
    sel2 = AN.p3_select(p3_toy_rows(tbl(2, 6)), q2)
    @test sel2.J[10] == 2 // 40 && sel2.J[9] == 4 // 40 && sel2.J_min == 1 // 20 && 9 in sel2.tol && sel2.star == 9 && same_sel(sel2, p3ref_select(p3_toy_rows(tbl(2, 6)), q2))
    sel3 = AN.p3_select(p3_toy_rows(tbl(3, 5)), q2)
    @test sel3.J[9] == 5 // 40 && sel3.J[10] == 2 // 40 && !(9 in sel3.tol) && 10 in sel3.tol && sel3.star == 10 && same_sel(sel3, p3ref_select(p3_toy_rows(tbl(3, 5)), q2))
end

@testset "D93j the order of C6 key by key: lowest C_R, |V|, w, v, bitmask of V, q_C, q_R (a fixture per key)" begin
    q = p3_params()
    mk(cost, V, qC, qR, w, v) = (cost = cost, V = V, qC = qC, qR = qR, w = w, v = v)
    S, F, M, E2, W = :SALI, :FLI, :MEGNO, :E2, :WBA
    # in each pair the first element wins on the key under test and is worse on every later key; equal on every earlier key
    pairs = [
        ("lowest C_R", mk(99 // 1, [S, F, M, E2, W], 5, 5, 2, 1), mk(100 // 1, [S, F, M], 1, 2, 0, 0)),
        ("smallest |V|", mk(100 // 1, [S, F, M], 3, 3, 2, 1), mk(100 // 1, [S, F, M, E2], 1, 3, 0, 0)),
        ("smallest w", mk(100 // 1, [F, M, E2], 3, 3, 0, 1), mk(100 // 1, [S, F, M], 1, 2, 1, 0)),
        ("smallest v", mk(100 // 1, [F, M, E2], 3, 3, 1, 0), mk(100 // 1, [S, F, M], 1, 2, 1, 1)),
        ("smallest bitmask of V", mk(100 // 1, [S, F, M], 3, 3, 1, 1), mk(100 // 1, [S, F, E2], 1, 2, 1, 1)),
        ("smallest q_C", mk(100 // 1, [S, F, M], 1, 3, 1, 1), mk(100 // 1, [S, F, M], 2, 2, 1, 1)),
        ("smallest q_R", mk(100 // 1, [S, F, M], 2, 2, 1, 1), mk(100 // 1, [S, F, M], 2, 3, 1, 1))]
    for (name, a, b) in pairs
        @test AN.p3_choose(q, [a, b]) == 1 && AN.p3_choose(q, [b, a]) == 2
    end
    # the bitmask of V in the bit order of the voters (SALI 0, FLI 1, MEGNO 2, E1 3, E2 4, WBA 5): {SALI, E2} = 17 < {FLI, MEGNO} = 6 is false; 6 < 17
    @test AN.p3_choose(q, [mk(1 // 1, [S, E2], 1, 2, 0, 0), mk(1 // 1, [F, M], 1, 2, 0, 0)]) == 2 && AN.p3_choose(q, [mk(1 // 1, [M, W], 1, 2, 0, 0), mk(1 // 1, [F, E2, W], 1, 2, 0, 0)]) == 1
    # exact costs: 1/3 + 1/3 < 2/3 + tiny is decided exactly, a tie in cost goes to the next key
    @test AN.p3_choose(q, [mk(2 // 3, [S, F, M], 1, 2, 0, 0), mk(1 // 3 + 1 // 3, [S, F, M], 1, 3, 0, 0)]) == 1
    # a total order: the same element from every permutation of 40 distinct random candidates, and the one the reference key gives
    rng = StableRNG(p3_seed('j', 1))
    th = p3ref_theta(q)
    cands = [(t = th[rand(rng, 1:length(th))], c = rand(rng, 1:4) // 2) for _ in 1:40]
    cands = unique(c -> (c.t.V, c.t.qC, c.t.qR, c.t.w, c.t.v), cands)
    el = [mk(c.c, c.t.V, c.t.qC, c.t.qR, c.t.w, c.t.v) for c in cands]
    best = el[argmin([p3ref_key(c.t, c.c) for c in cands])]
    for _ in 1:20
        perm = el[sortperm(rand(rng, length(el)))]
        @test perm[AN.p3_choose(q, perm)] == best
    end
end

@testset "D93j (c) F is empty: the outcome is 'no feasible rule', no theta*, no selection" begin
    q = p3_toy_q()
    # every orbit is called against its class: the error of every element is 1 (F2 excludes all)
    spec = vcat([(0, 4, (S2 = true,)) for _ in 1:4], [(0, 4, (Ts = true,)) for _ in 1:4], [(0, 4, (W2 = true,)) for _ in 1:3], [(4, 0, (rho = :R,)) for _ in 1:4])
    sel = AN.p3_select(p3_toy_rows(spec), q)
    @test sel.outcome === :no_feasible_rule && sel.star === nothing && isempty(sel.feasible) && isempty(sel.tol) && sel.J_min === nothing && sel.alone[1] == 8 && sel.alone[2] == 0
    @test all(i -> sel.first_fail[i] in 1:6 && sel.first_fail[i] != 0, eachindex(sel.thetas)) && same_sel(sel, p3ref_select(p3_toy_rows(spec), q))
    # an element with no decided orbit is not feasible (the error is undefined): every orbit votes undecided
    sel0 = AN.p3_select(p3_toy_rows(vcat([(0, 0, (S2 = true,)) for _ in 1:4], [(0, 0, (Ts = true,)) for _ in 1:4], [(0, 0, (W2 = true,)) for _ in 1:3], [(0, 0, (rho = :R,)) for _ in 1:4])), q)
    @test sel0.outcome === :no_feasible_rule && sel0.feasible == Int[] && all(s -> s.e === nothing && s.decided == 0, sel0.stats)
end

@testset "D93j the whole of C3 to C6 on the 36 orbits of the designed study equals the independent reference, element by element (2088 elements)" begin
    q = p3_params()
    rows = p3_design_rows(p3_design())
    @test length([r for r in rows if r.in_cal]) == 36
    sel = AN.p3_select(rows, q)
    ref = p3ref_select(rows, q)
    @test same_sel(sel, ref) && sel.outcome === :selected && sel.star == 147 && length(sel.feasible) == 62 && length(sel.tol) == 42 && sel.alone == (690, 1745, 1044, 1840, 681, 554)
    # P48: all Cal orbits enter C3 to C6, the K = 0 orbits (the 10 panel orbits of K index 1) too: 32 decided and 4 unresolved make the 36
    @test count(r -> r.in_cal && r.source === :panel && r.kindex == 1, rows) == 10 && sel.stats[147].decided + 4 == 36 && sel.stats[147].u_cal == 1 // 9
    # P49: shares are exact rationals; an undefined share (empty set) is not feasible and is not a number
    @test all(s -> s.u_cal isa Rational && s.fr_ts isa Rational && s.cost isa Rational, sel.stats)
    # the cost cap of F1 is the minimum over I_6 of the cost x N_c (E53, C4): 1810.5 x N_c
    @test sel.cost_cap == 18105 // 10 * 20480
end
