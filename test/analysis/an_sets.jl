# AC7 (sets E7), AC7/AC8 spread, convergence time, decisions H2 and H3 (E10, E11). Spec: analysis-spec.md.

an_so(id, K, class; control = false, exponent = 0.0, flags = String[], instability = "positive growth supported") =
    (; id, K, class, control, exponent, flags, instability)

@testset "AC7 orbit sets E7" begin
    p = an_params()
    Kc = p.K_c
    os = [an_so(1, 0.5, "C"), an_so(2, Kc, "C"), an_so(3, pi / 2, "C"), an_so(4, 0.5, "C"; instability = "unresolved"),
          an_so(5, 0.5, "S"), an_so(6, 2 * pi, "X"), an_so(7, 0.0, "S"), an_so(8, 0.5, "T"),
          an_so(9, 2 * pi, "C"), an_so(10, 10.0, "C"; instability = "none detected on the tested horizon"), an_so(11, 0.0, "C"),
          an_so(12, 0.5, "control"; control = true, exponent = 0.69), an_so(13, 0.5, "control"; control = true, exponent = 0.0),
          an_so(14, 0.5, "C"; flags = ["drift"]), an_so(15, 2 * pi, "C"; flags = ["excluded"]),
          an_so(16, 0.5, "S"; flags = ["diagnostics exceeded"])]
    s = AN.orbit_sets(os, p)
    @test sort(s.E_ch) == [1, 2, 3, 4, 9, 10, 12]
    # W: positive-growth orbits of class S or X at K > 0, and class C at K <= K_c (K_c itself counts); 16 is S at K > 0 and positive
    @test sort(s.W) == [1, 2, 5, 6, 16]
    # S: class C at K = 2 pi and K = 10, no label condition; the one not positive is counted
    @test sort(s.S) == [9, 10] && s.S_nonpositive == 1
    @test isempty(intersect(s.W, s.S))
    # drift and excluded orbits are removed from the sets and listed with the reason
    @test Set(d.id for d in s.dropped) == Set([14, 15]) && Dict(d.id => d.reason for d in s.dropped) == Dict(14 => "drift", 15 => "excluded")
    @test 14 ∉ s.E_ch && 14 ∉ s.W && 15 ∉ s.S && 15 ∉ s.E_ch
    # each condition switched: K_c above/below
    @test 2 ∈ AN.orbit_sets([an_so(2, Kc, "C")], p).W && 2 ∉ AN.orbit_sets([an_so(2, nextfloat(Kc), "C")], p).W
    @test 7 ∉ AN.orbit_sets([an_so(7, 0.0, "X")], p).W                 # K > 0 required for S and X
    # K_c comes from the parameter table
    @test 3 ∈ AN.orbit_sets([an_so(3, pi / 2, "C")], an_params(K_c = 2.0)).W
end

@testset "AC7 finite-time spread R(n) E10" begin
    steps = [1, 10, 100, 1000]
    lam = Dict{String, Vector{Float64}}(m => [0.5, 0.4, 0.31, 0.3] for m in AN_V)
    lam["L9"] = [0.6, 0.45, 0.32, 0.3001]
    lam["L2"] = [9.0, 9.0, 9.0, 9.0]                                       # L2 is not in V: it must not enter Delta_V
    lam["L4"] = [0.45, 0.39, 0.30, 0.3]
    R = AN.relative_spread(lam, 0.3)
    hi = [0.6, 0.45, 0.32, 0.3001]; lo = [0.45, 0.39, 0.30, 0.3]
    @test R ≈ (hi .- lo) ./ 0.3 atol = 1e-15
    @test AN.relative_spread(Dict(m => [0.3, 0.3] for m in AN_V), 0.3) == [0.0, 0.0]
end

@testset "AC7 H2 decision E10 and E12" begin
    p = an_params()
    seed = an_seed(7)
    r = AN.decide_h2(fill(0.2, 25), p; seed = seed)
    @test r.verdict === :supported && r.median == 0.2 && r.lo == 0.2 && r.hi == 0.2 && r.n == 25
    r = AN.decide_h2(fill(0.01, 25), p; seed = seed)
    @test r.verdict === :not_supported && r.median == 0.01
    # median above 0.10 but the interval reaches below 0.05: not supported (13 of 25 high, 12 low: the median flips in about 40 percent of resamples)
    r = AN.decide_h2(vcat(fill(0.2, 13), fill(0.01, 12)), p; seed = seed)
    @test r.median == 0.2 && r.lo == 0.01 && r.verdict === :not_supported
    # exactly at the limits passes, the next float fails (limits set to the computed median and interval end)
    v = vcat(fill(0.2, 20), fill(0.12, 10))
    r = AN.decide_h2(v, p; seed = seed)
    @test r.verdict === :supported
    @test AN.decide_h2(v, an_params(h2_median = r.median, h2_interval = r.lo); seed = seed).verdict === :supported
    @test AN.decide_h2(v, an_params(h2_median = nextfloat(r.median)); seed = seed).verdict === :not_supported
    @test AN.decide_h2(v, an_params(h2_interval = nextfloat(r.lo)); seed = seed).verdict === :not_supported
    # 19 orbits: not decidable with the count; 20: a verdict
    r = AN.decide_h2(fill(0.2, 19), p; seed = seed)
    @test r.verdict === :not_decidable && r.n == 19
    @test AN.decide_h2(fill(0.2, 20), p; seed = seed).verdict === :supported
end

@testset "AC8 convergence time E11" begin
    steps = an_grid(10_000)
    N = steps[end]
    lref = 0.3
    band = 0.05
    # lam = lref (1 + c / n): inside the band for n >= 20 |c|
    c = 20.0
    r = AN.convergence_time(steps, lref .* (1 .+ c ./ steps), lref, band)
    first_in = steps[findfirst(k -> abs(c / k) <= band, steps)]
    @test r.defined && r.n == first_in == 501
    @test AN.convergence_time(steps, lref .* (1 .+ 2.2 ./ steps), lref, band).n == 50      # 2.2 / k <= 0.05 from k = 44; the grid has 40 and 50
    # a late excursion: n_conv is the step after the last exit from the band
    lam = fill(lref, length(steps)); lam[end - 3] = 2 * lref
    r = AN.convergence_time(steps, lam, lref, band)
    @test r.defined && r.n == steps[end - 2]
    # always inside: the first recorded step
    @test AN.convergence_time(steps, fill(lref, length(steps)), lref, band) == (n = 1, defined = true)
    # outside the band at N itself: undefined, set to N
    lam = fill(lref, length(steps)); lam[end] = 2 * lref
    r = AN.convergence_time(steps, lam, lref, band)
    @test !r.defined && r.n == N
    # lambda_ref <= 0: the relative band is empty, undefined
    @test !AN.convergence_time(steps, fill(0.1, length(steps)), 0.0, band).defined
    @test !AN.convergence_time(steps, fill(0.1, length(steps)), -0.2, band).defined
    # the maximum over V; L2 is not in V
    lam = Dict{String, Vector{Float64}}(m => lref .* (1 .+ 2.0 ./ steps) for m in AN_V)
    lam["L9"] = lref .* (1 .+ 20.0 ./ steps)
    lam["L2"] = fill(5.0, length(steps))
    r = AN.orbit_convergence(steps, lam, lref, band)
    @test r.n == 501 && r.defined
    lam["L6"] = fill(5.0, length(steps))                                   # one method never converges: undefined, N, counted
    r = AN.orbit_convergence(steps, lam, lref, band)
    @test !r.defined && r.n == N
    # the band comes from the parameter table argument
    @test AN.convergence_time(steps, lref .* (1 .+ 2.2 ./ steps), lref, 0.2).n == 13       # 2.2 / k <= 0.2 from k = 11; the grid has 10 and 13
end

@testset "AC10 H3 ratio E11 and E12" begin
    p = an_params()
    seed = an_seed(10)
    r = AN.decide_h3(fill(100.0, 30), fill(5.0, 30), p; seed = seed)
    @test r.verdict === :supported && r.ratio == 20.0 && r.lo == 20.0 && r.n_w == 30 && r.n_s == 30
    # ratio 10 exactly (>=) passes; the interval end is strict (> 3)
    r = AN.decide_h3(fill(100.0, 30), fill(10.0, 30), p; seed = seed)
    @test r.ratio == 10.0 && r.verdict === :supported
    @test AN.decide_h3(fill(100.0, 30), fill(10.0, 30), an_params(h3_ratio = r.ratio); seed = seed).verdict === :supported
    @test AN.decide_h3(fill(100.0, 30), fill(10.0, 30), an_params(h3_ratio = nextfloat(r.ratio)); seed = seed).verdict === :not_supported
    @test AN.decide_h3(fill(100.0, 30), fill(10.0, 30), an_params(h3_interval = r.lo); seed = seed).verdict === :not_supported      # lo == limit: strict
    @test AN.decide_h3(fill(100.0, 30), fill(10.0, 30), an_params(h3_interval = prevfloat(r.lo)); seed = seed).verdict === :supported
    @test AN.decide_h3(fill(100.0, 30), fill(50.0, 30), p; seed = seed).verdict === :not_supported       # ratio 2
    # W and S are resampled independently: with S constant the ratio interval is the W interval over that constant
    w = collect(range(50.0, 400.0; length = 31))
    r = AN.decide_h3(w, fill(10.0, 25), p; seed = seed)
    bw = AN.bootstrap_interval(w, "H3/W"; seed = seed, resamples = p.bootstrap_resamples, level = p.bootstrap_level)
    @test isapprox(r.lo, bw.lo / 10; atol = 1e-12) && isapprox(r.hi, bw.hi / 10; atol = 1e-12) && isapprox(r.ratio, bw.point / 10; atol = 1e-12)
    # changing the S values does not change the W resamples (and the reverse)
    s2 = AN.decide_h3(w, fill(20.0, 25), p; seed = seed)
    @test s2.lo ≈ bw.lo / 20 atol = 1e-12
    # either set below the minimum: not decidable with both counts
    r = AN.decide_h3(fill(100.0, 19), fill(5.0, 30), p; seed = seed)
    @test r.verdict === :not_decidable && r.n_w == 19
    r = AN.decide_h3(fill(100.0, 30), fill(5.0, 19), p; seed = seed)
    @test r.verdict === :not_decidable && r.n_s == 19
end
