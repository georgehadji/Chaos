# Amendment D94, function level: D94c (screening points), D94f (stream seeds), D94h (D sizes and indices, O23c), D94i (D_B allocation and
# selection, O23c). Spec: experiments-spec.md, "Amendment D94": E32, E33, E37, E38, O23b, O23c. API assumed: the test plan.

@testset "D94c screening points equal the independent generator bitwise" begin
    for seed in (D94_FA.seed, D94_FB.seed), i in 1:3
        pts = EX.random_screen_points(seed, i, 1000)
        ref = d94_points(seed, i, 1000)
        @test length(pts) == 1000
        @test all(j -> (pts[j][1], pts[j][2]) == ref[j], 1:1000)                         # tolerance 0
        @test all(p -> 0 <= p[1] < 2pi && 0 <= p[2] < 2pi, pts)
        @test allunique([(p[1], p[2]) for p in pts])
    end
    # a change of seed or of K index changes the points; the rule does not depend on K (K = 0 and K > 0 are index values only)
    a = [(p[1], p[2]) for p in EX.random_screen_points(D94_FA.seed, 1, 50)]
    @test a != [(p[1], p[2]) for p in EX.random_screen_points(D94_FA.seed + 1, 1, 50)]
    @test a != [(p[1], p[2]) for p in EX.random_screen_points(D94_FA.seed, 2, 50)]
    @test a == [(p[1], p[2]) for p in EX.random_screen_points(D94_FA.seed, 1, 1000)][1:50]               # a prefix: M only cuts the list
    # the generation is sequential on one stream: u_1, u_2, ... in order, x first then p (P78), RN of the 256-bit value
    u = rand(d94_stream(d94_screen_text(D94_FA.seed, 3)), 4)
    setprecision(BigFloat, 256) do
        p1 = EX.random_screen_points(D94_FA.seed, 3, 2)
        @test p1[1][1] == Float64(2 * BigFloat(pi) * BigFloat(u[1])) && p1[1][2] == Float64(2 * BigFloat(pi) * BigFloat(u[2]))
        @test p1[2][1] == Float64(2 * BigFloat(pi) * BigFloat(u[3])) && p1[2][2] == Float64(2 * BigFloat(pi) * BigFloat(u[4]))
    end
end

@testset "D94f stream seeds equal the independent SHA-256 for 1000 seeded tuples each (O23b)" begin
    rng = StableRNG(d94_seed('f'))
    S = [rand(rng, 0:10^15) for _ in 1:1000]
    I = [rand(rng, 1:12) for _ in 1:1000]
    C = [rand(rng, collect(D94_CLASSES)) for _ in 1:1000]
    R = [rand(rng, 1:9) for _ in 1:1000]
    d94_prefetch_many(vcat([d94_screen_text(S[k], I[k]) for k in 1:1000], [d94_round_text(S[k], I[k], C[k], R[k]) for k in 1:1000],
                           [d94_d_text(S[k], R[k], I[k], C[k]) for k in 1:1000], [d94_db_text(S[k], R[k], I[k]) for k in 1:1000]))
    @test all(k -> EX.screen_stream_seed(S[k], I[k]) == d94_sseed(d94_screen_text(S[k], I[k])), 1:1000)
    @test all(k -> EX.round_stream_seed(S[k], I[k], C[k], R[k]) == d94_sseed(d94_round_text(S[k], I[k], C[k], R[k])), 1:1000)
    @test all(k -> EX.d_stream_seed(S[k], R[k], I[k], C[k]) == d94_sseed(d94_d_text(S[k], R[k], I[k], C[k])), 1:1000)
    @test all(k -> EX.db_stream_seed(S[k], R[k], I[k]) == d94_sseed(d94_db_text(S[k], R[k], I[k])), 1:1000)
    # the round 1 text differs from the text of E15, so the seed differs from the base cell seed
    @test all(i -> all(c -> EX.round_stream_seed(EX_SEED, i, c, 1) != EX.cell_stream_seed(EX_SEED, i, c), D94_CLASSES), 1:6)
    # no two streams of the amendment coincide (I35): screening, round cells, D cells, D_B strata of one seed
    seeds = UInt64[EX.screen_stream_seed(EX_SEED, i) for i in 1:6]
    append!(seeds, vec([EX.round_stream_seed(EX_SEED, i, c, r) for i in 1:6, c in D94_CLASSES, r in 1:5]))
    append!(seeds, vec([EX.d_stream_seed(EX_SEED, r, i, c) for i in 1:6, c in D94_CLASSES, r in 1:5]))
    append!(seeds, vec([EX.db_stream_seed(EX_SEED, R, i) for i in 1:6, R in 1:5]))
    @test allunique(seeds) && length(seeds) == 6 + 150 + 150 + 30
    # classes may be given as symbols
    @test EX.round_stream_seed(EX_SEED, 2, :C, 3) == EX.round_stream_seed(EX_SEED, 2, "C", 3)
end

@testset "D94h D sizes and indices equal the exact rational reference for all n up to 300 and 40 fractions (O23c)" begin
    rng = StableRNG(d94_seed('h'))
    fr = Set{Rational{Int}}()
    while length(fr) < 40
        b = rand(rng, 1:16)
        push!(fr, rand(rng, 0:b) // b)
    end
    fracs = sort(collect(fr))
    bad = 0
    for f in fracs, n in 0:300
        EX.d_size(f, n) == d94_dsize(Rational{BigInt}(f), n) || (bad += 1)
    end
    @test bad == 0
    @test EX.d_size(1 // 2, 3) == 2 && EX.d_size(1 // 4, 2) == 1 && EX.d_size(1 // 5, 7) == 1 && EX.d_size(0 // 1, 9) == 0 && EX.d_size(1 // 1, 9) == 9
    # 20 fixtures with short cells and n up to 300: the selected draw indices equal the reference stream
    d94_prefetch_many([d94_d_text(EX_SEED, r, i, c) for r in 1:3, i in 1:4, c in ("T", "C")])
    for k in 1:20
        n, f = rand(rng, 0:300), rand(rng, fracs)
        r, i, c = rand(rng, 1:3), rand(rng, 1:4), rand(rng, ("T", "C"))
        got = EX.d_indices(EX_SEED, r, i, c, n, f)
        @test got == d94_dindices(EX_SEED, r, i, c, n, Rational{BigInt}(f))
        @test issorted(got) && allunique(got) && all(j -> 1 <= j <= n, got) && length(got) == d94_dsize(Rational{BigInt}(f), n)
    end
    @test isempty(EX.d_indices(EX_SEED, 1, 1, "T", 0, 1 // 2)) && isempty(EX.d_indices(EX_SEED, 1, 1, "T", 5, 0 // 1))
end

@testset "D94i D_B allocation and selection" begin
    rng = StableRNG(d94_seed('i'))
    # largest remainder allocation against the exact rational reference: 2000 seeded cases
    bad = 0
    for _ in 1:2000
        n = [rand(rng, 0:300) for _ in 1:rand(rng, 1:6)]
        b = rand(rng, 0:400)
        al = EX.db_allocation(n, b)
        (al == d94_alloc(n, b) && sum(al) == min(b, sum(n)) && all(i -> 0 <= al[i] <= n[i], eachindex(n))) || (bad += 1)
    end
    @test bad == 0
    @test EX.db_allocation([1, 6, 2], 5) == [1, 3, 1] && EX.db_allocation([1, 4, 1], 5) == [1, 3, 1] && EX.db_allocation([1, 2, 0], 5) == [1, 2, 0]
    @test EX.db_allocation([3, 3], 4) == [2, 2] && EX.db_allocation([5, 5, 5], 4) == [2, 1, 1]
    @test EX.db_allocation([4, 4], 0) == [0, 0]
    # selection: R = 1, 2, 3 on the hand D sets (hand tables of the plan) and on a larger seeded D
    for R in 1:3
        DR = d94_hand_D(R)
        got = EX.db_subset(EX_RUN_SEED, R, DR, 5)
        @test got == d94_db(EX_RUN_SEED, R, DR, 5)
        @test length(got) == min(5, length(DR)) && issubset(got, [x.id for x in DR])
        @test EX.db_subset(EX_RUN_SEED, R, reverse(DR), 5) == got                       # a function of the set, not of the order of the list
        @test isempty(EX.db_subset(EX_RUN_SEED, R, DR, 0))
    end
    @test EX.db_subset(EX_RUN_SEED, 1, d94_hand_D(1), 5) == [x.id for x in d94_hand_D(1)]      # |D_1| = 3 <= 5: D_B = D_R
    big = [(id = k, k_index = rand(rng, 1:5)) for k in 1:240]
    d94_prefetch_many([d94_db_text(EX_SEED, R, i) for R in 1:3, i in 1:5])
    for R in 1:3, b in (7, 50, 239, 240, 500)
        got = EX.db_subset(EX_SEED, R, big, b)
        @test got == d94_db(EX_SEED, R, big, b) && length(got) == min(b, 240)
        @test EX.db_allocation([count(x -> x.k_index == i, big) for i in 1:5], b) == [count(id -> big[id].k_index == i, got) for i in 1:5]
    end
    # R enters the stream text: the selection for R = 2 differs from R = 3 on the same D only through R
    a2, a3 = EX.db_subset(EX_SEED, 2, big, 50), EX.db_subset(EX_SEED, 3, big, 50)
    @test a2 == d94_db(EX_SEED, 2, big, 50) && a3 == d94_db(EX_SEED, 3, big, 50) && a2 != a3
end
