# D93k PP for HB, D93l PP for HC, D93m streams; hand and analytic cases of O17, statistical oracles of O16. Spec: analysis-spec.md, Amendment D93, E55 to E57, O16, O17, criteria
# D93k, D93l, D93m (S). Seeds: the test_seed of the table (20261030, the stream rule of E55) for the module streams; p3_seed('l', 1) for the independent simulation of the HC oracle.
# The stream vectors below were computed with Python 3 (hashlib, statistics.NormalDist, own SplitMix64 and xoshiro256++), never with this repository: see the test plan, "Stream vectors".

const P3_TS = 20261030

@testset "D93m streams: SplitMix64, xoshiro256++, the stream seed of E55, uniforms, integers, normals (stored constants from an independent implementation)" begin
    # the generators against their published first outputs
    @test AN.p3_splitmix_state(UInt64(0)) == (0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4, 0x06c45d188009454f, 0xf88bb8a8724c81ec)
    g = AN.p3_generator((UInt64(1), UInt64(2), UInt64(3), UInt64(4)))
    @test [AN.p3_u64!(g) for _ in 1:4] == UInt64[41943041, 58720359, 3588806011781223, 3591011842654386]
    # the seed: first 8 bytes, big-endian, of SHA-256 of "chaoslab-analysis-1|<test_seed>|p3|<purpose>|<unit>"; a pair is written "a,b" (a reading pinned by the tests)
    @test AN.p3_stream_seed(P3_TS, "optimism", 1) == UInt64(5692226309264838358) && AN.p3_stream_seed(P3_TS, "pp-hb", (0, 1)) == UInt64(3437253884877108471) &&
          AN.p3_stream_seed(P3_TS, "pp-hc", (20, 1)) == UInt64(15836479185517542896) && AN.p3_stream_seed(P3_TS, "pp-hb", (3, 17)) == UInt64(4982238177473391303) &&
          AN.p3_stream_seed(P3_TS, "pp-hc", (49, 1)) == UInt64(16471504700301138262) && AN.p3_stream_seed(P3_TS, "optimism", 6) == UInt64(982380449943276992)
    @test AN.p3_splitmix_state(AN.p3_stream_seed(P3_TS, "optimism", 1)) == (UInt64(15300154784399154531), UInt64(1595416676562361740), UInt64(9699484540660912394), UInt64(17671760824607871275))
    g = AN.p3_stream(P3_TS, "optimism", 1)
    @test [AN.p3_u64!(g) for _ in 1:4] == UInt64[1960354936839448390, 11115532021063746583, 17181476028355880092, 5824767452774032242]
    # a uniform is (k + 1/2) / 2^53 with k the top 53 bits of an output: in the open unit interval
    g = AN.p3_stream(P3_TS, "optimism", 1)
    @test [AN.p3_uniform!(g) for _ in 1:4] == [0.10627105406820064, 0.6025741982784751, 0.931409681822553, 0.3157612763260232]
    # normals: one inverse-function call per variate on one uniform (the generator advances by one output per normal), AS241 accuracy against the Python library to 1e-12
    g = AN.p3_stream(P3_TS, "pp-hb", (0, 1))
    z = [AN.p3_normal!(g) for _ in 1:4]
    @test isapprox(z, [0.13039229363167548, -0.36652077601596667, 0.17412024595138367, -0.38487600750808904]; rtol = 1e-12)
    h = AN.p3_stream(P3_TS, "pp-hb", (0, 1))
    for _ in 1:4; AN.p3_u64!(h); end
    @test AN.p3_u64!(g) == AN.p3_u64!(h)
    u = [1e-12, 1e-6, 0.01, 0.1, 0.3, 0.5, 0.7, 0.975, 0.9999, 1 - 1e-9]
    zq = [-7.034483825301132, -4.753424308822899, -2.3263478740408408, -1.2815515655446008, -0.5244005127080407, 0.0, 0.5244005127080407, 1.9599639845400536, 3.7190164854557084, 5.997807019601638]
    @test all(i -> isapprox(AN.p3_normal_quantile(u[i]), zq[i]; rtol = 1e-12, atol = 1e-15), eachindex(u)) && AN.p3_normal_quantile(0.5) == 0.0
    # an integer in 1..m by rejection, unbiased: modulo with rejection of the incomplete block (a reading pinned by the tests)
    g = AN.p3_stream(P3_TS, "optimism", 1)
    @test [AN.p3_randint!(g, 36) for _ in 1:8] == [23, 20, 5, 19, 31, 7, 13, 10]
    g = AN.p3_stream(P3_TS, "optimism", 1)
    @test [AN.p3_randint!(g, 3) for _ in 1:8] == [2, 2, 2, 1, 1, 1, 1, 1]
    g = AN.p3_stream(P3_TS, "optimism", 2)
    @test [AN.p3_randint!(g, 36) for _ in 1:8] == [35, 13, 30, 8, 5, 24, 12, 30]
    g = AN.p3_stream(P3_TS, "optimism", 9)
    cnt = zeros(Int, 3)
    for _ in 1:30_000; cnt[AN.p3_randint!(g, 3)] += 1; end
    @test all(c -> abs(c - 10_000) < 350, cnt)                                   # 4.3 standard errors
    big = Int(2)^62 + 1                                                         # rejection happens in about one draw of two
    @test all(x -> 1 <= x <= big, [AN.p3_randint!(g, big) for _ in 1:200])
    # the stream of a unit does not change when another unit is drawn from, or another purpose, or when the same unit is drawn again (a function of the unit only)
    a = AN.p3_stream(P3_TS, "optimism", 2)
    for _ in 1:100; AN.p3_u64!(AN.p3_stream(P3_TS, "optimism", 1)); AN.p3_u64!(AN.p3_stream(P3_TS, "pp-hb", (2, 1))); end
    @test [AN.p3_u64!(a) for _ in 1:2] == UInt64[0x02f85a4a319f736e, 0xe08e0a65f476bd94] && AN.p3_stream_seed(P3_TS, "optimism", 2) != AN.p3_stream_seed(P3_TS, "optimism", 1)
    @test length(unique(AN.p3_stream_seed(P3_TS, pu, un) for pu in ("optimism", "pp-hb", "pp-hc") for un in (1, 2, (1, 1), (1, 2), (2, 1)))) == 15
    # distinct from every p1 and p2 label of D89 E34 and the base E12: the text carries the prefix p3 that no other has
    others = ["H2/W", "H3/S", "F2/C/1", "T3/C/1", "p2|H1|W2|E_all|R1|all", "p2|H3|S2|E_all|R2|minus-E1", "optimism|1", "pp-hb|0,1", "p2|optimism|1"]
    @test all(l -> AN.p3_stream_seed(P3_TS, "optimism", 1) != AN.stream_seed(l, P3_TS), others) && all(l -> AN.p3_stream_seed(AN_SEED, "optimism", 1) != AN.stream_seed(l, AN_SEED), others)
    # the rule is the base rule of E12 (D89 E34) with the label p3|<purpose>|<unit>: the same seed as the existing stream function gives for that label
    @test AN.p3_stream_seed(P3_TS, "optimism", 1) == AN.stream_seed("p3|optimism|1", P3_TS) && AN.p3_stream_seed(P3_TS, "pp-hb", (3, 17)) == AN.stream_seed("p3|pp-hb|3,17", P3_TS)
    # the same draws as the reference generator written in the suite
    g = AN.p3_stream(P3_TS, "pp-hc", (49, 7)); r = p3ref_stream(P3_TS, "pp-hc", (49, 7))
    @test [AN.p3_u64!(g) for _ in 1:50] == [p3ref_next!(r) for _ in 1:50]
end

@testset "D93k, D93l, O17: exact McNemar tail, rounding, Hanley and McNeil, AUC, the binormal shift, the votes of HC (hand and analytic cases)" begin
    # the exact one-sided tail P(Bin(m, 1/2) >= n01), m = n10 + n01: every (n10, n01) up to m = 10 against the binomial table, and p = 1 for m = 0
    ok = true
    for m in 0:10, n01 in 0:m
        ok &= AN.p3_mcnemar_p(m - n01, n01) == p3ref_mcnemar(m - n01, n01)
    end
    @test ok && AN.p3_mcnemar_p(0, 0) == 1 && AN.p3_mcnemar_p(0, 1) == 1 // 2 && AN.p3_mcnemar_p(1, 1) == 3 // 4 && AN.p3_mcnemar_p(0, 5) == 1 // 32 && AN.p3_mcnemar_p(5, 0) == 1 && AN.p3_mcnemar_p(2, 3) == 1 // 2
    # the critical-value table of PP (decided by exact integer arithmetic): the smallest n01 with p <= alpha equals the direct sum for m up to 200
    alpha = 5 // 100
    bad = Int[]
    for m in 1:200
        direct = nothing
        for n01 in 0:m
            sum(binomial(big(m), j) for j in n01:m) * 100 <= big(2)^m * 5 && (direct = n01; break)
        end
        AN.p3_mcnemar_critical(m, alpha) == direct || push!(bad, m)
    end
    @test isempty(bad) && AN.p3_mcnemar_critical(1, alpha) === nothing && AN.p3_mcnemar_critical(4, alpha) === nothing && AN.p3_mcnemar_critical(5, alpha) == 5 && AN.p3_mcnemar_critical(10, alpha) == 9
    # rounding: half away from zero (71.25 to 71, 308.75 to 309, 237.5 to 238, -2.5 to -3), ceil5 (71 to 75, 75 to 75, 76 to 80); exact products of keys
    @test [AN.p3_round_half_away(x) for x in (7125 // 100, 30875 // 100, 2375 // 10, -5 // 2, 1 // 2, 0 // 1, 99 // 100)] == [71, 309, 238, -3, 1, 0, 1]
    @test [AN.p3_ceil_step(x // 1, 5) for x in (71, 75, 76, 1)] == [75, 75, 80, 5] && AN.p3_ceil_step(9975 // 100, 5) == 100 && AN.p3_ceil_step(43225 // 100, 5) == 435
    # Hanley and McNeil: A = 1/2 with n_a = n_n = 2 gives 5/48 (Q1 = Q2 = 1/3); A = 1 gives 0; hand cases computed with exact fractions in Python
    hm(A, na, nn) = AN.p3_hanley_mcneil(A, na, nn)
    @test isapprox(hm(0.5, 2, 2), 5 / 48; rtol = 1e-13) && hm(1.0, 2, 2) == 0 && isapprox(hm(0.8, 10, 20), 197 / 22500; rtol = 1e-12) && isapprox(hm(0.6, 5, 100), 6597 / 350000; rtol = 1e-12) &&
          isapprox(hm(0.7, 40, 800), 7119 / 3200000; rtol = 1e-12)
    @test all(t -> isapprox(hm(t[1], t[2], t[3]), Float64(p3ref_hm(Rational{BigInt}(t[1]), t[2], t[3])); rtol = 1e-12), ((0.55, 20, 400), (0.7, 62, 1240), (0.9, 49, 980), (0.5, 3, 7)))
    # AUC (D89 O7 b, Mann-Whitney, ties one half): separated 1, identical 0.5, one tied pair one half, reversed 0; empty class: no value
    @test AN.p3_auc([5, 6], [1, 2]) == 1 && AN.p3_auc([1, 2, 3], [1, 2, 3]) == 1 // 2 && AN.p3_auc([2, 3], [2, 1]) == 7 // 8 && AN.p3_auc([1, 2], [3, 4]) == 0 && AN.p3_auc([2.0], [2.0]) == 1 // 2
    @test AN.p3_auc(Float64[], [1.0]) === nothing && AN.p3_auc([1.0], Float64[]) === nothing
    # the binormal shift d(A) = sqrt(2) Phi^-1(A), d(1/2) = 0; A outside the open unit interval is refused
    @test AN.p3_binormal_shift(0.5) == 0 && isapprox(AN.p3_binormal_shift(0.75), 0.9538725524089398; rtol = 1e-12) && isapprox(AN.p3_binormal_shift(0.55), 0.1777119809885155; rtol = 1e-12) &&
          isapprox(AN.p3_binormal_shift(0.8), 1.1902321628999901; rtol = 1e-12) && isapprox(AN.p3_binormal_shift(0.9), 1.8123876048736471; rtol = 1e-12)
    @test an_rejects(() -> AN.p3_binormal_shift(0.0), ["auc"]) && an_rejects(() -> AN.p3_binormal_shift(1.0), ["auc"]) && an_rejects(() -> AN.p3_binormal_shift(1.2), ["auc"]) && an_rejects(() -> AN.p3_binormal_shift(-0.1), ["auc"])
    # the votes of an orbit of HC: plain label by the decided votes (n_C > n_R chaotic, n_R > n_C regular, a tie a failure), f = 1 iff the label differs from rho, s_I = 1 - |n_C - n_R| / (n_C + n_R)
    vs(v, rho) = AN.p3_vote_summary(v, rho)
    c, r, u = :chaotic, :regular, :undecided
    a = vs([c, c, r], :R); @test (a.nC, a.nR, a.label, a.f) == (2, 1, :chaotic, true) && a.s == 1 - 1 / 3
    a = vs([c, c, r], :C); @test (a.label, a.f) == (:chaotic, false) && a.s == 1 - 1 / 3
    a = vs([c, r], :C); @test (a.nC, a.nR, a.label, a.f, a.s) == (1, 1, :tie, true, 1.0)             # a tie is a failure, also for rho = C
    a = vs([c, r, u], :R); @test (a.label, a.f, a.s) == (:tie, true, 1.0)
    a = vs([u, u, u], :C); @test (a.nC, a.nR, a.label, a.f, a.s) == (0, 0, :tie, true, 0.0)           # no decided vote: s = 0 and a failure
    a = vs([r, r, u], :R); @test (a.label, a.f, a.s) == (:regular, false, 0.0)
    a = vs([c, c, c], :C); @test (a.label, a.f, a.s) == (:chaotic, false, 0.0)
    a = vs([c, r, r, r], :C); @test (a.label, a.f, a.s) == (:regular, true, 1 - 2 / 4)
end

@testset "D93k the draw path of E56 d: derived count, allowed points, sizes, rounding" begin
    q = p3_params()
    path = AN.p3_draw_path(q)
    ref = p3ref_path(q)
    @test length(path) == 17 == length(ref) && all(i -> (path[i].j, path[i].rho, path[i].nT, path[i].nRW, path[i].nW, path[i].allowed) == (ref[i].j, ref[i].rho, ref[i].nT, ref[i].nRW, ref[i].nW, ref[i].allowed), 1:17)
    @test path[4].rho == 7 // 4 && (path[4].nT, path[4].nRW, path[4].nW) == (100, 333, 433) && path[end].rho == 5 && (path[end].nT, path[end].nW) == (285, 1235) && path[1].rho == 1 && (path[1].nT, path[1].nW) == (57, 247)
    @test [x.j for x in path if x.allowed] == 3:16 && !path[3].allowed && path[3].nT == 86 && path[2].nT == 71 && path[2].nRW == 238 && path[2].nW == 309             # 71.25 -> 71, 237.5 -> 238, 308.75 -> 309 (O17 d)
    # the point count is derived from the keys: step 0.5 from 1 to 5 gives 9 points; a start above r_max gives none; the floors act (n_t_min, n_w_min)
    @test length(AN.p3_draw_path(merge(q, (; rho_path_second = 1.5)))) == 9 && isempty(AN.p3_draw_path(merge(q, (; rho_path_start = 6.0, rho_path_second = 6.25)))) &&
          [x.j for x in AN.p3_draw_path(merge(q, (; n_t_min = 200))) if x.allowed] == 10:16 && count(x -> x.allowed, AN.p3_draw_path(merge(q, (; n_w_min = 1000)))) == 4
    # the targets: the smallest allowed point whose power reaches p_hb; ceil5 of y x rho*; the cap and the no-allowed-point cases
    mk(powers; qq = q) = [merge((j = x.j, rho = x.rho, nT = x.nT, nRW = x.nRW, nW = x.nW, allowed = x.allowed), x.allowed ? (successes = Int(round(powers[x.j + 1] * 1000)), runs = 1000) : (successes = nothing, runs = 1000)) for x in AN.p3_draw_path(qq)]
    pw = [0.1 * min(8, i) for i in 1:17]; pw[6] = 0.8                      # j = 5 (rho = 2.25) is the first allowed point that reaches 0.8; j = 3, 4 are at 0.4, 0.5
    t = AN.p3_hb_targets(mk(pw), q)
    @test t.outcome === :ok && t.j == 5 && t.rho_star == 9 // 4 && t.n_T == 130 && t.n_W == 560 && t.cap === false                                  # 57 x 2.25 = 128.25 -> 130; 247 x 2.25 = 555.75 -> 560
    pw2 = fill(0.0, 17); pw2[4] = 0.799; pw2[5] = 0.8                       # exact comparison with p_hb: 799 of 1000 does not reach 0.8, 800 of 1000 does (j = 3 below, j = 4 reaches)
    @test AN.p3_hb_targets(mk(pw2), q).j == 4
    pw3 = fill(0.0, 17); pw3[1] = 0.95; pw3[2] = 0.9                       # a point below the floors that would reach p_hb is not allowed, hence not chosen
    t3 = AN.p3_hb_targets(mk(pw3), q)
    @test t3.outcome === :cap && t3.cap === true && t3.n_T == 285 && t3.n_W == 1235 && t3.power_cap == 0.0 && t3.rho_star == 5
    pw4 = [0.7 for _ in 1:17]
    t4 = AN.p3_hb_targets(mk(pw4), q)
    @test t4.outcome === :cap && t4.cap && (t4.n_T, t4.n_W) == (285, 1235) && t4.power_cap == 0.7
    qn = merge(q, (; n_t_min = 1000))
    @test AN.p3_hb_targets(mk(fill(0.9, 17); qq = qn), qn).outcome === :no_allowed_point
end

# ------------------------------------------------------------------ the simulation of one point (E56 c, e): the exact power oracle and the reference simulation
@testset "D93k, O16 a to c: the simulated power of one comparison equals the exact multinomial sum within 3 standard errors on 14 parameter sets; type-I below alpha" begin
    q = merge(p3_params(), (; alpha_hb = 0.05, rho_hb = 0.5, c_min = 0.5))
    sets = [(6, 1 // 20, 3 // 5, 1 // 2), (8, 1 // 20, 1 // 2, 1 // 2), (9, 1 // 20, 2 // 5, 1 // 2), (10, 1 // 10, 1 // 2, 1 // 2), (12, 1 // 10, 2 // 5, 1 // 2), (12, 1 // 10, 1 // 2, 1 // 4),
            (14, 1 // 10, 3 // 10, 1 // 2), (16, 1 // 5, 1 // 2, 1 // 2), (20, 3 // 10, 3 // 5, 1 // 2), (14, 1 // 20, 3 // 5, 1 // 2), (18, 1 // 10, 3 // 5, 1 // 4), (22, 1 // 10, 1 // 2, 0 // 1),
            (24, 1 // 10, 3 // 5, 1 // 2), (12, 0 // 1, 3 // 10, 1 // 2)]
    runs = 20_000
    zs = Float64[]; powers = Float64[]
    for (i, (n, p, pX, s)) in enumerate(sets)
        ex = Float64(p3ref_exact_power(n, p, pX, s, 5 // 100, 1 // 2))
        model = p3_model(rule_T = (p = Float64(p), d = 1.0), comp = Dict(:SALI => (T = (p = Float64(pX), s = Float64(s)), RW = (p = 0.0, s = 0.0))), ts = [:SALI])
        r = AN.p3_simulate_hb(model, n, 0, q; point = 0, runs = runs, threads = 2)
        se = sqrt(ex * (1 - ex) / runs)
        push!(zs, (r.power - ex) / se); push!(powers, ex)
        @test r.runs == runs && r.successes isa Integer && r.power == r.successes / runs && isapprox(r.se, sqrt(r.power * (1 - r.power) / runs); rtol = 1e-12)
    end
    @test all(z -> abs(z) <= 3, zs) && minimum(powers) < 0.2 && maximum(powers) > 0.9 && count(x -> 0.1 <= x <= 0.99, powers) == 14
    # the same draws as the reference simulation (the order of the nine uniforms, the event u < probability, the streams (pp-hb, (point, i)) with i from 1): equal counts, all 13 comparisons
    rate = (T = (p = 0.35, s = 0.6), RW = (p = 0.5, s = 0.7))
    rate2 = (T = (p = 0.5, s = 0.4), RW = (p = 0.6, s = 0.5))
    comp = Dict(X => (X in (:SALI, :FLI, :E2, :LE) ? rate : rate2) for X in P2_X7)
    model = p3_model(rule_T = (p = 0.12, d = 0.9), rule_RW = (p = 0.08, d = 0.8), comp = comp, ts = [:SALI, :FLI, :MEGNO, :E1, :E2, :WBA], w2 = P2_X7)
    for (point, nT, nRW) in ((3, 40, 70), (0, 20, 30), (9, 57, 95))
        a = AN.p3_simulate_hb(model, nT, nRW, q; point = point, runs = 400, threads = 1)
        @test a.successes == p3ref_hb_power(q.test_seed, point, model, nT, nRW, q, 400)
    end
    # one thread and several threads give identical counts (every run owns its stream, E55)
    b1 = AN.p3_simulate_hb(model, 40, 70, q; point = 3, runs = 400, threads = 1)
    b4 = AN.p3_simulate_hb(model, 40, 70, q; point = 3, runs = 400, threads = 4)
    @test b1.successes == b4.successes && b1.power == b4.power
    # the joint decision is at most the smallest exact single power plus 3 standard errors, and at least one minus the sum of one minus the single powers (Bonferroni) minus 3 standard errors
    n = 16
    cs = [(:SALI, 1 // 10, 7 // 10, 1 // 2), (:FLI, 1 // 10, 3 // 5, 1 // 2), (:E2, 1 // 10, 3 // 5, 1 // 4)]
    comp3 = Dict{Symbol, Any}(X => (T = (p = Float64(pX), s = Float64(s)), RW = (p = 0.0, s = 0.0)) for (X, _, pX, s) in cs)
    mj = p3_model(rule_T = (p = 0.1, d = 1.0), comp = comp3, ts = [:SALI, :FLI, :E2])
    singles = [Float64(p3ref_exact_power(n, 1 // 10, pX, s, 5 // 100, 1 // 2)) for (_, _, pX, s) in cs]
    rj = AN.p3_simulate_hb(mj, n, 0, q; point = 0, runs = runs, threads = 2)
    @test rj.power <= minimum(singles) + 3 * rj.se && rj.power >= 1 - sum(1 .- singles) - 3 * rj.se && 0 < rj.power < minimum(singles)
    # the coverage condition: with no orbit decided outside the rule's calls the decision cannot succeed when the decided share is below c_min
    low = p3_model(rule_T = (p = 0.1, d = 0.2), comp = Dict(:SALI => (T = (p = 0.7, s = 0.5), RW = (p = 0.0, s = 0.0))), ts = [:SALI])
    @test AN.p3_simulate_hb(low, 24, 0, q; point = 0, runs = 2000, threads = 1).successes == p3ref_hb_power(q.test_seed, 0, low, 24, 0, q, 2000) &&
          AN.p3_simulate_hb(low, 24, 0, q; point = 0, runs = 2000, threads = 1).power < AN.p3_simulate_hb(merge(low, (rule = Dict(:T => (p = 0.1, d = 1.0), :RW => low.rule[:RW]),)), 24, 0, q; point = 0, runs = 2000, threads = 1).power
    # type I (O16 c, E56 f): the rule rate set equal to the rate of each comparator in turn, everything else unchanged; the same streams as the main simulation at that point
    ti = AN.p3_hb_type1(model, 40, 70, q; point = 3, runs = 400, threads = 2)
    @test Set(keys(ti.rates)) == Set(P2_X7) && ti.max_joint == maximum(ti.rates[X].joint for X in P2_X7)
    for X in P2_X7
        m2 = merge(model, (rule = Dict(:T => (p = model.comp[X][:T].p, d = model.rule[:T].d), :RW => (p = model.comp[X][:RW].p, d = model.rule[:RW].d)),))
        @test ti.rates[X].joint == p3ref_hb_power(q.test_seed, 3, m2, 40, 70, q, 400) / 400 && ti.rates[X].single >= ti.rates[X].joint
    end
    one = p3_model(rule_T = (p = 0.2, d = 1.0), comp = Dict(:SALI => (T = (p = 0.2, s = 0.5), RW = (p = 0.0, s = 0.0))), ts = [:SALI])
    t1 = AN.p3_hb_type1(one, 20, 0, q; point = 0, runs = 20_000, threads = 2)
    @test t1.rates[:SALI].joint < q.alpha_hb && abs(t1.rates[:SALI].joint - Float64(p3ref_exact_power(20, 1 // 5, 1 // 5, 1 // 2, 5 // 100, 1 // 2))) <= 3 * sqrt(0.0085 / 20_000) &&
          t1.rates[:SALI].single < q.alpha_hb       # (the other indicators are not compared in this model: their type-I joint rates are not about a comparison of the model)
end

@testset "D93l PP for HC: shrinkage, the n_f list, the deliverable numbers, the target, the flag; O16 d" begin
    q = merge(p3_params(), (; pp_runs = 120))
    # AUC_I_used = 0.5 + auc_shrink (AUC_I_cal - 0.5); AUC_P = auc_p; the true difference is their difference; the shifts d = sqrt(2) Phi^-1(A)
    hc = AN.p3_pp_hc(q, 0.9; runs = 120, threads = 2)
    @test isapprox(hc.auc_i_used, 0.7; atol = 1e-15) && hc.auc_p == 0.55 && isapprox(hc.true_difference, 0.15; atol = 1e-15) && isapprox(hc.d_i, 0.7416143171871158; rtol = 1e-12) && isapprox(hc.d_p, 0.1777119809885155; rtol = 1e-12)
    # the n_f list: nf_start, nf_second, ... up to the smallest listed value >= the worst-case deliverable failures (floor(1309 / 21) = 62: 65); the search does not stop at the expected number (49)
    @test [x.nf for x in hc.table] == collect(20:5:65) && hc.deliverable_expected == 49 == fld(1048, 21) && hc.deliverable_worst == 62 == fld(1309, 21)
    @test all(x -> 0 <= x.power <= 1 && x.runs == 120 && x.successes isa Integer && x.power == x.successes / 120 && isapprox(x.se, sqrt(x.power * (1 - x.power) / 120); rtol = 1e-12), hc.table)
    # the power at the expected deliverable number (49 is not on the list) is a separate simulation with its own stream unit
    @test hc.power_expected.nf == 49 && 0 <= hc.power_expected.power <= 1 && hc.power_expected.runs == 120
    # the target: the first listed n_f whose power is at least p_hc (exact comparison), or "not reached"; the flag underpowered if not reached or above the expected deliverable number
    function target_ok(h)
        first_ok = findfirst(x -> x.successes * 10 >= 8 * x.runs, h.table)
        first_ok === nothing && return h.target === nothing && h.underpowered === true
        x = h.table[first_ok]
        return h.target == x.nf && h.target_power == x.power && h.underpowered == (x.nf > h.deliverable_expected)
    end
    @test target_ok(hc)
    hi = AN.p3_pp_hc(q, 1.0; runs = 120, threads = 1)                         # AUC_I_used = 0.75, a true difference of 0.2: the power is high from small n_f
    @test target_ok(hi) && hi.target !== nothing && hi.target <= 40 && hi.underpowered === false
    lo = AN.p3_pp_hc(q, 0.5; runs = 120, threads = 1)                         # AUC_I_used = 0.5 below AUC_P: the difference is negative, never reached
    @test target_ok(lo) && lo.target === nothing && lo.underpowered === true && all(x -> x.power < 0.8, lo.table) && lo.power_cap isa Float64
    # the same counts with one and several threads (E55)
    @test [x.successes for x in AN.p3_pp_hc(q, 0.9; runs = 120, threads = 1).table] == [x.successes for x in hc.table]
    # O16 d: the rate of the decision with a true difference of zero is below alpha_ci plus 3 standard errors; with a true difference of 0.3 (AUC 0.8 against 0.5) the power
    # agrees within 3 standard errors with an independent simulation that has its own generator (Box-Muller on a StableRNG) and the Hanley and McNeil variance of the reference
    runs = 3000
    h0 = AN.p3_simulate_hc(q, 0.55, 0.55; nf = 30, runs = runs, threads = 2)
    @test h0.power <= q.alpha_ci + 3 * sqrt(q.alpha_ci * (1 - q.alpha_ci) / runs) && h0.runs == runs
    nf = 10; nn = 20 * nf
    h3 = AN.p3_simulate_hc(q, 0.8, 0.5; nf = nf, runs = runs, threads = 2)
    rng = StableRNG(p3_seed('l', 1))
    gauss() = sqrt(-2 * log(1 - rand(rng))) * cos(2 * pi * rand(rng))
    di, dp = sqrt(2) * 0.8416212335729143, 0.0
    hit = 0
    for _ in 1:runs
        fI = [di + gauss() for _ in 1:nf]; nI = [gauss() for _ in 1:nn]; fP = [dp + gauss() for _ in 1:nf]; nP = [gauss() for _ in 1:nn]
        AI, AP = Float64(p3ref_auc(fI, nI)), Float64(p3ref_auc(fP, nP))
        Delta = AI - AP
        SE = sqrt(Float64(p3ref_hm(Rational{BigInt}(AI), nf, nn)) + Float64(p3ref_hm(Rational{BigInt}(AP), nf, nn)))
        (Delta >= q.delta_c && Delta - q.z_one_sided * SE > 0) && (hit += 1)
    end
    pi_, pj = h3.power, hit / runs
    @test abs(pi_ - pj) <= 3 * sqrt(pi_ * (1 - pi_) / runs + pj * (1 - pj) / runs) && 0.3 < pi_ < 0.99
    # an AUC outside the open unit interval is refused
    @test an_rejects(() -> AN.p3_simulate_hc(q, 1.0, 0.5; nf = 20, runs = 10), ["auc"]) && an_rejects(() -> AN.p3_pp_hc(merge(q, (; auc_p = 1.0)), 0.9; runs = 10), ["auc"])
end

@testset "D93k the optimism bootstrap (E56 a): resamples from the streams (optimism, b), C3 to C6 on each, skipped resamples counted, min_set not applied inside" begin
    q = merge(p3_toy_q(), (; b_opt = 12, test_seed = P3_TS))
    rows = p3_toy_rows(P3_TOY_A)
    ref = p3ref_optimism(rows, q)
    @test ref.used >= 2 && ref.skipped >= 1                                              # the scenario has both kinds of resample
    @test any(b -> ref.thetas[b] !== nothing && ref.sizes[b].Ts < q.min_set, 1:q.b_opt) # a resample with a set below min_set is used as is (P66)
    o = AN.p3_optimism(rows, q)
    @test o.used == ref.used && o.skipped == ref.skipped && o.used + o.skipped == q.b_opt && o.thetas == ref.thetas
    @test isapprox(o.optimism.T_s, ref.optimism.T_s; rtol = 1e-12, atol = 1e-15) && isapprox(o.optimism.W2, ref.optimism.W2; rtol = 1e-12, atol = 1e-15)
    # the draws come from the streams of the units b: the resample of unit b is the same whatever the other units and the thread count (a unit owns its stream)
    @test AN.p3_optimism(rows, merge(q, (; b_opt = 5))).thetas == ref.thetas[1:5] && AN.p3_optimism(rows, q; threads = 3).thetas == ref.thetas
    # the same on the designed study (36 orbits, 2088 elements, real thresholds scaled): the function agrees with the reference on its used count and optimism
    q2 = merge(p3_params(), (; b_opt = 6))
    rows2 = p3_design_rows(p3_design())
    o2 = AN.p3_optimism(rows2, q2)
    @test o2.used + o2.skipped == 6 && o2.used >= 1 && o2.optimism.T_s isa Float64 && o2.optimism.W2 isa Float64
end

@testset "D93k the rates of PP for HB from Cal at theta*: corrected rule rate, floor and clamp, strata by difference, comparator rates at matched budgets, shared share, decided share" begin
    q = p3_toy_q(; shared_default = 0.37)
    rows = p3_toy_rows(P3_TOY_A)
    # orbit 7 (a T_s orbit, three regular votes) is regular for the comparators at exactly their matched budgets of theta* (E* costs 11439360): E1 at b_5 = 1280 (k = 5), E2 at
    # b_8 = 10240, LE at N_c = b_9 (never matched: capped), and not at b_P; WBA is not an indicator of the toy table
    rows[7].V[:E1] = p2_sym("uuuuuruuuu"); rows[7].V[:E2] = p2_sym("uuuuuuuuru"); rows[7].V[:LE] = p2_sym("uuuuuuuuur")
    ref = p3ref_select(rows, q)
    th = ref.thetas[11]
    @test th == (index = 11, V = [:SALI, :FLI, :MEGNO, :E2], mask = 23, qC = 2, qR = 3, w = 0, v = 0)
    m = AN.p3_matched(q, th.V)
    @test m.CR == 11439360 && (m.X[:E1].b, m.X[:E2].b, m.X[:SALI].b, m.X[:MEGNO].b) == (1280, 10240, 5120, 5120) && m.X[:LE].capped && m.le_read == q.N_c
    th_q = (index = th.index, V = th.V, mask = th.mask, qC = th.qC, qR = th.qR, w = th.w, v = th.v)
    # optimism (T_s: 0.1, W2: -0.05): r_corr = r + max(optimism, 0); R_W by counts by difference, then clamped to [0, 1]; the rule rate is r_corr above the floor
    h = AN.p3_hb_rates(rows, q, th_q, (T_s = 0.1, W2 = -0.05))
    @test isapprox(h.r[:T_s], 1 / 3; rtol = 1e-14) && isapprox(h.r[:W2], 5 / 7; rtol = 1e-14) && isapprox(h.r_corr[:T_s], 1 / 3 + 0.1; rtol = 1e-14) && isapprox(h.r_corr[:W2], 5 / 7; rtol = 1e-14) &&
          isapprox(h.r_corr[:R_W], (5 / 7 * 7 - (1 / 3 + 0.1) * 3) / 4; rtol = 1e-12) && isapprox(h.r_corr[:R_W], 0.925; rtol = 1e-12)
    @test isapprox(h.p_rule[:T_s], 1 / 3 + 0.1; rtol = 1e-14) && isapprox(h.p_rule[:R_W], 0.925; rtol = 1e-12)
    # the comparator rates p_X(S) = FR_X(S) at the matched budgets (R_W by counts by difference): SALI, FLI, MEGNO regular on orbit 7 and on the four R_W orbits; E1 and E2 and LE only on orbit 7
    for X in (:SALI, :FLI, :MEGNO)
        @test isapprox(h.p_comp[(X, :T_s)], 1 / 3; rtol = 1e-14) && isapprox(h.p_comp[(X, :R_W)], 1.0; rtol = 1e-14)
    end
    for X in (:E1, :E2, :LE)
        @test isapprox(h.p_comp[(X, :T_s)], 1 / 3; rtol = 1e-14) && h.p_comp[(X, :R_W)] == 0.0
    end
    # the shared share: rule false-regular orbits that are also false-regular under X over the rule false-regular orbits (orbit 7 in T_s; the four orbits in R_W)
    for X in (:SALI, :FLI, :MEGNO)
        @test h.s_comp[(X, :T_s)] == 1.0 && h.s_comp[(X, :R_W)] == 1.0
    end
    for X in (:E1, :E2, :LE)
        @test h.s_comp[(X, :T_s)] == 1.0 && h.s_comp[(X, :R_W)] == 0.0
    end
    @test h.d[:T_s] == 1.0 && h.d[:R_W] == 1.0 && Set(k[1] for k in keys(h.p_comp)) == Set(Symbol.(q.indicators_i7))
    # the floor branch and the clamp: theta = (q_C 2, q_R 4) calls nothing regular; r_corr(T_s) = 0.004 < 0.01 gives max(2 x 0.004, 0.01); r_corr(R_W) = -0.003 clamps to 0 and gives the floor
    th4 = ref.thetas[12]
    th4q = (index = th4.index, V = th4.V, mask = th4.mask, qC = th4.qC, qR = th4.qR, w = th4.w, v = th4.v)
    @test (th4.qC, th4.qR) == (2, 4)
    h4 = AN.p3_hb_rates(rows, q, th4q, (T_s = 0.004, W2 = 0.0))
    @test h4.r[:T_s] == 0 && h4.r[:W2] == 0 && h4.r_corr[:T_s] == 0.004 && h4.r_corr[:R_W] == 0.0 && h4.p_rule[:T_s] == 0.01 && h4.p_rule[:R_W] == 0.01
    @test AN.p3_hb_rates(rows, q, th4q, (T_s = 0.007, W2 = 0.0)).p_rule[:T_s] == 0.014 && AN.p3_hb_rates(rows, q, th4q, (T_s = 0.02, W2 = 0.0)).p_rule[:T_s] == 0.02 &&
          AN.p3_hb_rates(rows, q, th4q, (T_s = 0.0, W2 = 0.0)).p_rule[:T_s] == 0.01
    # the rule has no false-regular orbit: the shared share is shared_default (here 0.37) for every comparator in both strata; the decided share of theta* by counts by difference
    @test all(v -> v == 0.37, values(h4.s_comp)) && isapprox(h4.d[:T_s], 2 / 3; rtol = 1e-14) && isapprox(h4.d[:R_W], 0.0; atol = 1e-14)
    # the model of one orbit (E56 c) built from the rates: the rule rate and decided share per stratum, comparator rate and shared share per stratum, the two comparator sets
    md = AN.p3_hb_model(h)
    @test md.rule[:T].p == h.p_rule[:T_s] && md.rule[:RW].p == h.p_rule[:R_W] && md.rule[:T].d == h.d[:T_s] && md.rule[:RW].d == h.d[:R_W] &&
          md.comp[:E2][:T].p == h.p_comp[(:E2, :T_s)] && md.comp[:E2][:RW].s == h.s_comp[(:E2, :R_W)] && md.ts == [:SALI, :FLI, :MEGNO, :E1, :E2] && md.w2 == Symbol.(q.indicators_i7)
end

@testset "D93k the whole path of E56: every allowed point simulated, equal to the reference counts; the targets from the powers" begin
    q = merge(p3_params(), (; pp_runs = 30))
    rate = (T = (p = 0.3, s = 0.5), RW = (p = 0.4, s = 0.6))
    comp = Dict(X => rate for X in P2_X7)
    model = p3_model(rule_T = (p = 0.1, d = 0.9), rule_RW = (p = 0.12, d = 0.9), comp = comp, ts = [:SALI, :FLI, :MEGNO, :E1, :E2, :WBA], w2 = P2_X7)
    path = AN.p3_hb_path(model, q; runs = 30, threads = 2)
    ref = p3ref_path(q)
    @test length(path) == 17 && all(i -> path[i].allowed == ref[i].allowed && (path[i].nT, path[i].nRW, path[i].nW) == (ref[i].nT, ref[i].nRW, ref[i].nW), 1:17)
    @test all(x -> x.successes === nothing, [x for x in path if !x.allowed]) && length([x for x in path if x.allowed]) == 14
    @test all(x -> x.successes == p3ref_hb_power(q.test_seed, x.j, model, x.nT, x.nRW, q, 30), [x for x in path if x.allowed])
    @test all(x -> x.power == x.successes / 30 && isapprox(x.se, sqrt(x.power * (1 - x.power) / 30); rtol = 1e-12), [x for x in path if x.allowed])
    t = AN.p3_hb_targets(path, q)
    ok = [x for x in path if x.allowed && x.successes * 10 >= 8 * x.runs]
    if isempty(ok)
        @test t.cap === true && (t.n_T, t.n_W) == (285, 1235)
    else
        x = first(ok)
        @test t.cap === false && t.j == x.j && t.n_T == p3ref_ceil_step(57 * x.rho, 5) && t.n_W == p3ref_ceil_step(247 * x.rho, 5)
    end
end
