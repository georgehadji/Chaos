# D89o bootstrap of the P2 statistics (E34, E12 reused, I20, O7 c). Spec: analysis-spec.md, Amendment D89, E34, criterion D89o.
# The independent oracle is p2ref_bootstrap (an_p2_ref.jl: stream from the SHA-256 rule, draws in the order of E12, type 7 percentiles).

const P2_LAB1 = "p2|H1|W2|E_all|R1|all"
const P2_LAB2 = "p2|H3|S2|E_all|R2|minus-E1"

function p2_random_recs(rng, n, I; pT = 1.0)
    syms = (:chaotic, :regular, :undecided)
    return [(id = o, V = Dict(X => Symbol[syms[rand(rng, 1:3)] for _ in 1:10] for X in I), T = Dict(X => 40 * 2^rand(rng, 0:8) for X in I)) for o in 1:n]
end

@testset "D89o stream: label grammar, platform-independent constants, distinct from the P1 streams" begin
    @test AN.p2_label(:H1, "W2", :E_all, :R1, "all") == P2_LAB1 && AN.p2_label(:H3, "S2", :E8, :R2, "minus-E1") == "p2|H3|S2|E8|R2|minus-E1"
    @test AN.p2_label(:H4design, "T_s_design", :E_all, :R1, "all") == "p2|H4design|T_s_design|E_all|R1|all"
    # the seed of the stream is the first 8 bytes, big-endian, of SHA-256 of chaoslab-analysis-1|<seed>|<label> (E12); the first four outputs are stored constants
    @test AN.stream_seed(P2_LAB1, AN_SEED) == UInt64(3444927369723232618) == p2_stream_seed(P2_LAB1, AN_SEED)
    rng = AN.bootstrap_rng(P2_LAB1, AN_SEED)
    @test [rand(rng, UInt64) for _ in 1:4] == UInt64[0x986079ec7c98ba24, 0x1086a453d2baa0ac, 0x2c374cbb5b651d9f, 0x8a250f9f138a33db]
    @test AN.stream_seed(P2_LAB2, AN_SEED) == UInt64(17828611372690436285)
    rng = AN.bootstrap_rng(P2_LAB2, AN_SEED)
    @test [rand(rng, UInt64) for _ in 1:4] == UInt64[0x0943a72e539ace21, 0xe51db9067f4b005b, 0x8ccebb39745cc993, 0x346ec6584d67d5f2]
    # the prefix p2| separates every P2 stream from every P1 stream (P1 labels are H2/W, H3/S, F2/<class>/<k>, T3/<class>/<k>)
    @test all(l -> AN.stream_seed(l, AN_SEED) != AN.stream_seed("p2|" * l, AN_SEED), ("H2/W", "H3/W", "H3/S", "F2/C/1", "T3/C/1"))
    # replicate r draws m indices uniformly with replacement from the stream, in order: column r of the matrix
    ix = AN.p2_resample(P2_LAB1, 11, 5; seed = AN_SEED)
    @test size(ix) == (11, 5) && ix == p2_draws(P2_LAB1, 11, 5, AN_SEED) && all(1 .<= ix .<= 11)
    rng = StableRNG(p2_stream_seed(P2_LAB1, AN_SEED))
    @test [rand(rng, 1:11) for _ in 1:11] == ix[:, 1] && [rand(rng, 1:11) for _ in 1:11] == ix[:, 2]
    @test AN.p2_resample(P2_LAB1, 11, 5; seed = AN_SEED) == ix && AN.p2_resample(P2_LAB2, 11, 5; seed = AN_SEED) != ix && AN.p2_resample(P2_LAB1, 11, 5; seed = AN_SEED + 1) != ix
end

@testset "D89o intervals equal the independent implementation (H1 to H5); same seed and label give identical intervals" begin
    p = merge(p2_params(), (; bootstrap_resamples = 300))
    rng = StableRNG(p2_seed('o', 3))
    S = p2_random_recs(rng, 12, P2_X7); W = p2_random_recs(rng, 14, P2_X7)
    for (h, sets, labels) in ((:H1, (W,), P2_LAB1), (:H2, (S,), "p2|H2|S2|E_all|R1|all"), (:H3, (W, S), ("p2|H3|W2|E_all|R1|all", "p2|H3|S2|E_all|R1|all")),
                              (:H4, (W,), "p2|H4|T_s|E_all|R1|all"), (:H5, (W,), "p2|H5|W2|E_all|R1|all"))
        a = AN.p2_bootstrap(h, sets, labels, p)
        b = p2ref_bootstrap(h, sets, labels, p)
        @test a.dropped == b.dropped && a.defined == b.defined
        if b.interval === nothing
            @test a.interval === nothing
        else
            @test keys(a.interval) == keys(b.interval) && all(k -> isapprox(a.interval[k].lo, b.interval[k].lo; rtol = 1e-12) && isapprox(a.interval[k].hi, b.interval[k].hi; rtol = 1e-12), keys(b.interval))
        end
        @test a == AN.p2_bootstrap(h, sets, labels, p)                                                  # same seed and label: identical
    end
    # another label gives other intervals (the end points are discrete, so several labels are compared and not all may agree); another seed too
    ends(l; seed = AN_SEED) = (r = AN.p2_bootstrap(:H5, (W,), l, p; seed = seed); (r.interval[:F_med].lo, r.interval[:F_med].hi))
    @test length(unique([ends("p2|H5|W2|E_all|R1|v$k") for k in 1:6])) > 1 && length(unique([ends("p2|H5|W2|E_all|R1|all"; seed = s) for s in AN_SEED .+ (0:5)])) > 1
    # a smaller indicator set changes the statistic but not the draws (I20): the interval equals the independent one on the same draws
    I6 = P2_X7[1:6]
    a, b = AN.p2_bootstrap(:H2, (S,), "p2|H2|S2|E_all|R1|minus-WBA", p; I = I6), p2ref_bootstrap(:H2, (S,), "p2|H2|S2|E_all|R1|minus-WBA", p; I = I6)
    @test isapprox(a.interval[:Q].lo, b.interval[:Q].lo; rtol = 1e-12) && isapprox(a.interval[:Q].hi, b.interval[:Q].hi; rtol = 1e-12)
    # computing another label in between does not change the first (the draws of a replicate are a function of the label and the replicate only)
    first_ = AN.p2_bootstrap(:H2, (S,), "p2|H2|S2|E_all|R1|all", p)
    AN.p2_bootstrap(:H1, (W,), P2_LAB1, p); AN.p2_bootstrap(:H2, (S,), "p2|H2|S2|E_all|R2|all", p)
    @test AN.p2_bootstrap(:H2, (S,), "p2|H2|S2|E_all|R1|all", p) == first_
    # H3 resamples W and S independently: changing the S label changes only the S draws, hence the ratio intervals change with it, the W medians do not
    r1 = AN.p2_bootstrap(:H3, (W, S), ("p2|H3|W2|E_all|R1|all", "p2|H3|S2|E_all|R1|all"), p)
    r2 = AN.p2_bootstrap(:H3, (W, S), ("p2|H3|W2|E_all|R1|all", "p2|H3|S2|E_all|R1|other"), p)
    b2 = p2ref_bootstrap(:H3, (W, S), ("p2|H3|W2|E_all|R1|all", "p2|H3|S2|E_all|R1|other"), p)
    @test r1.interval != r2.interval && all(k -> isapprox(r2.interval[k].lo, b2.interval[k].lo; rtol = 1e-12), keys(b2.interval))
end

@testset "D89o exact enumeration of 6^6 resamples (O7 c): the percentiles of U and Q at B = 1e4 agree within 4 standard errors" begin
    p = merge(p2_params(), (; bootstrap_resamples = 10_000))
    B = 10_000; se = sqrt(0.025 * 0.975 / B)
    # U: six orbits, two indicators, verdicts at b_P (c,c) (r,r) (c,c) (c,r) (r,c) (c,r): three of six are split
    vv = [(c_, c_), (r_, r_), (c_, c_), (c_, r_), (r_, c_), (c_, r_)]
    recs = [p2_rec(o, Dict(:LE => (vv[o][1], c_), :SALI => (vv[o][2], c_))) for o in 1:6]
    split = [vv[o][1] !== vv[o][2] for o in 1:6]
    exactU = sort!(vec([count(i -> split[i], idx) / 6 for idx in Iterators.product(ntuple(_ -> 1:6, 6)...)]))
    # Q: six orbits, seven indicators, censored times fixed by a recorded seed
    rng = StableRNG(p2_seed('o', 4))
    Tm = Dict(X => [40 * 2^rand(rng, 0:7) for _ in 1:6] for X in P2_X7)
    recq = [p2_rec(o, Dict(X => (r_, c_) for X in P2_X7); T = Dict(X => Tm[X][o] for X in P2_X7)) for o in 1:6]
    med6(v) = (s = sort(v); (s[3] + s[4]) / 2)
    exactQ = sort!(vec([(m = [med6([Tm[X][i] for i in idx]) for X in P2_X7]; maximum(m) / minimum(m)) for idx in Iterators.product(ntuple(_ -> 1:6, 6)...)]))
    @test length(exactU) == 6^6 == length(exactQ)
    for (ex, iv) in ((exactU, AN.p2_bootstrap(:H1, (recs,), "p2|O7c|U", p; I = [:LE, :SALI]).interval[:U]),
                     (exactQ, AN.p2_bootstrap(:H2, (recq,), "p2|O7c|Q", p).interval[:Q]))
        cdf(v) = count(<=(v), ex) / length(ex)
        cdfm(v) = count(<(v), ex) / length(ex)
        @test cdfm(iv.lo) <= 0.025 + 4se && cdf(iv.lo) >= 0.025 - 4se
        @test cdfm(iv.hi) <= 0.975 + 4se && cdf(iv.hi) >= 0.975 - 4se
    end
end

@testset "D89o coverage: 1000 samples of 40 orbits, the 95 percent interval of U contains the true value in 900 to 1000 cases" begin
    p = merge(p2_params(), (; bootstrap_resamples = 200))
    rng = StableRNG(p2_seed('o', 2))
    hit = 0
    for k in 1:1000
        recs = [(u = rand(rng); u < 0.3 ? p2_rec(o, Dict(:LE => (c_, c_), :SALI => (r_, c_))) : u < 0.65 ? p2_rec(o, Dict(:LE => (c_, c_), :SALI => (c_, c_))) :
                 p2_rec(o, Dict(:LE => (r_, c_), :SALI => (r_, c_)))) for o in 1:40]
        iv = AN.p2_bootstrap(:H1, (recs,), "p2|cov|$k", p; I = [:LE, :SALI]).interval[:U]       # the population fraction of split orbits is 0.3
        hit += iv.lo <= 0.3 <= iv.hi
    end
    @test 900 <= hit <= 1000
end

@testset "D89o a statistic undefined on more than 5 percent of the replicates: undefined interval, not decidable" begin
    p = merge(p2_params(), (; bootstrap_resamples = 1000, min_set = 2))
    # eight orbits; E1 is decided at b_P and at k = 9 on only `d` of them: F_E1 is undefined when a resample misses all of them
    function recs8(d)
        return [p2_rec(o, Dict(X => (X === :E1 && o > d ? (u_, c_) : (r_, c_)) for X in P2_X7)) for o in 1:8]
    end
    ok = AN.p2_bootstrap(:H5, (recs8(3),), "p2|undef|3", p)           # (5/8)^8 = 2.3 percent of the replicates are dropped: the interval is defined
    @test ok.defined && ok.interval !== nothing && 0 < ok.dropped <= 50
    bad = AN.p2_bootstrap(:H5, (recs8(2),), "p2|undef|2", p)          # (6/8)^8 = 10 percent: more than 5 percent are dropped: undefined
    @test !bad.defined && bad.interval === nothing && bad.dropped > 50
    ref = p2ref_bootstrap(:H5, (recs8(3),), "p2|undef|3", p)
    @test ok.dropped == ref.dropped && p2ref_bootstrap(:H5, (recs8(2),), "p2|undef|2", p).dropped == bad.dropped
    # the decision of that hypothesis is NOT DECIDABLE
    @test AN.p2_decide(:H5, (n = 8, F_med = 0.3, F_med_lo = nothing), p).verdict === :not_decidable
end
