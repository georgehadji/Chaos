# AC9 bootstrap (E12), oracle O2. Spec: analysis-spec.md, Equations E12 and Oracle O2.

an_stream_seed(label, seed) = reinterpret(UInt64, reverse(hex2bytes(an_sha256("chaoslab-analysis-1|$(seed)|$(label)"))[1:8]))[1]

@testset "AC9 stream keyed by seed and label" begin
    for (label, seed) in (("H2/W", AN_SEED), ("H3/S", AN_SEED), ("F2/C/3", 1234))
        s = an_stream_seed(label, seed)                                  # independent SHA-256, first 8 bytes big-endian
        @test AN.stream_seed(label, seed) == s
        rng, ref = AN.bootstrap_rng(label, seed), StableRNG(s)
        @test [rand(rng, UInt64) for _ in 1:4] == [rand(ref, UInt64) for _ in 1:4]
    end
    @test AN.stream_seed("H2/W", AN_SEED) != AN.stream_seed("H3/W", AN_SEED)
    @test AN.stream_seed("H2/W", AN_SEED) != AN.stream_seed("H2/W", AN_SEED + 1)
end

@testset "AC9 determinism of the interval" begin
    rng = StableRNG(an_seed(9, 1))
    v = rand(rng, 40)
    kw = (; seed = AN_SEED, resamples = 2000, level = 0.95)
    a = AN.bootstrap_interval(v, "T/x"; kw...)
    @test a == AN.bootstrap_interval(v, "T/x"; kw...)
    # another label or seed means another stream: the end points of a median of 40 are discrete, so several streams are compared and not all may agree
    ends(lab, s) = (r = AN.bootstrap_interval(v, lab; seed = s, resamples = 2000, level = 0.95); (r.lo, r.hi))
    @test length(unique([ends(l, AN_SEED) for l in ("T/x", "T/y", "T/z", "T/w", "T/v")])) > 1
    @test AN.bootstrap_interval(v, "T/y"; kw...).point == a.point
    @test isapprox(a.point, an_median(v); atol = 1e-15) && a.lo <= a.point <= a.hi          # the median of 40
    @test length(unique([ends("T/x", AN_SEED + k) for k in 0:4])) > 1
    # a constant sample has a degenerate interval
    c = AN.bootstrap_interval(fill(3.0, 30), "T/c"; kw...)
    @test c.point == 3.0 && c.lo == 3.0 && c.hi == 3.0
end

@testset "AC9 exact enumeration O2: samples of size 6" begin
    x = [0.7, 1.9, 2.2, 3.1, 5.5, 8.0]
    meds = Float64[]
    for idx in Iterators.product(ntuple(_ -> 1:6, 6)...)
        s = sort([x[i] for i in idx])
        push!(meds, (s[3] + s[4]) / 2)
    end
    sort!(meds)
    @test length(meds) == 6^6
    B = 10_000
    se = sqrt(0.025 * 0.975 / B)                    # standard error of the probability level of a quantile from B draws
    cdf(v) = count(<=(v), meds) / length(meds)
    cdfm(v) = count(<(v), meds) / length(meds)
    for lab in ("O2/a", "O2/b", "O2/c")
        r = AN.bootstrap_interval(x, lab; seed = AN_SEED, resamples = B, level = 0.95)
        # the module's end points sit at probability levels 0.025 and 0.975 of the exact distribution, within 4 standard errors
        @test cdfm(r.lo) <= 0.025 + 4se && cdf(r.lo) >= 0.025 - 4se
        @test cdfm(r.hi) <= 0.975 + 4se && cdf(r.hi) >= 0.975 - 4se
    end
end

@testset "AC9 coverage of the 95 percent interval of a median" begin
    rng = StableRNG(an_seed(9, 2))
    hit = 0
    for k in 1:1000
        v = rand(rng, 40)                           # uniform on (0, 1): the true median is 0.5
        r = AN.bootstrap_interval(v, "cov/$k"; seed = AN_SEED, resamples = 1000, level = 0.95)
        hit += r.lo <= 0.5 <= r.hi
    end
    @test 900 <= hit <= 1000
end

@testset "AC9 level and ratio" begin
    v = collect(1.0:50.0)
    r90 = AN.bootstrap_interval(v, "lvl"; seed = AN_SEED, resamples = 4000, level = 0.90)
    r95 = AN.bootstrap_interval(v, "lvl"; seed = AN_SEED, resamples = 4000, level = 0.95)
    @test r95.lo <= r90.lo && r95.hi >= r90.hi
    w, s = collect(10.0:10.0:300.0), collect(1.0:30.0)
    q = AN.bootstrap_ratio(w, s, "H3/W", "H3/S"; seed = AN_SEED, resamples = 2000, level = 0.95)
    @test isapprox(q.point, an_median(w) / an_median(s); atol = 1e-12) && q.lo <= q.point <= q.hi
    @test q == AN.bootstrap_ratio(w, s, "H3/W", "H3/S"; seed = AN_SEED, resamples = 2000, level = 0.95)
end
