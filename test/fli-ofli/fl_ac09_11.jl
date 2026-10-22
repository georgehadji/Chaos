# AC9 to AC11: growth laws (ensemble level), tau dependence, precision ladder and representation. Ensemble-level and loose-interval criteria (K18): the
# intervals and allowances are the spec's proposals (P6), not source values.

@testset "AC9a regular growth law" begin
    _ = (FLI,)
    # 20 regular-class starts at K = 0.5 (S7), Float64, FLI_B, n up to 1e6, seeded frames: the least-squares slope of FLI_B against ln n over 21 log-spaced n in
    # [1e4, 1e6] lies in [0.7, 1.3] for at least 80 percent of the starts (Q8: slope 1). Transient starts are reported and at most 20 percent may be excluded (they
    # count as not passing). Deviation (test plan, item 8): the spec's detector "FLI_B below ln n at n = 1e4" is not a transient detector (8 of 20 starts are below ln n
    # with slope 1.00 to 1.03, FLI_B = ln(n c) with c below 1 for a weak twist or a small frame component); here a start is in the transient when its slope over the first
    # decade [1e4, 1e5] is outside [0.7, 1.3]. The number below ln n is reported.
    m = 20
    c = fl_case(:half, Float64, m)
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    run = fl_run((FLI(),), c.sys, c.x0, c.W0, 1_000_000, 1; steps = ns)
    first_decade = findall(<=(100_000), ns)
    transient, passed, below = 0, 0, 0
    slopes = Float64[]
    for r in run
        sn = r[1].snaps
        @test length(sn) == length(ns)
        slope = fl_slope(log.(Float64.(ns)), Float64[s.fli for s in sn])
        early = fl_slope(log.(Float64.(ns[first_decade])), Float64[sn[i].fli for i in first_decade])
        push!(slopes, slope)
        below += sn[1].fli < log(1e4)
        if !(0.7 <= early <= 1.3)
            transient += 1
        elseif 0.7 <= slope <= 1.3
            passed += 1
        end
    end
    @info "AC9a reported: slopes of FLI_B against ln n, starts in the transient, starts with FLI_B below ln n at n = 1e4" slopes transient below
    @test transient <= 0.2 * m
    @test passed >= ceil(Int, 0.8 * m)
end

@testset "AC9b chaotic rate against lambda_1" begin
    _ = (FLI, Benettin)
    # K in {2 pi, 5, 10}, Float64, 1000 starts: the least-squares slope of FLI_B against n over n in [20, 200] and the L1 series X(1e4) of a k = 1 run on the
    # identical orbit (first seeded column). Pass: |mean(slope) - mean(X)| <= 4 combined standard errors (L1 S5) + 4/180 (finite-window allowance of the spec,
    # larger than the 2/180 of sali because the O(1) part includes the drawdown and the basis maximum). (c) K = 5 and 10: the ensemble means of FLI_1(n)/n and
    # FLI_B(n)/n at n = 1e3 are within 4 combined standard errors + 4/n of X(1e4). Reported: R3 (sanity rates against the Chirikov law and the K = 2 pi row of
    # benchmarks.md, never pass/fail), R4 (first handoff at which FLI_B reaches ln 1e20, the sources' saturation cap).
    m = 1000
    sanity = Dict(:five => log(5 / 2), :ten => log(10 / 2), :two_pi => 1.172)      # benchmarks.md rows "large K" (ln(K/2)) and "single chaotic orbit" (K = 2 pi)
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = fl_seed(9, ci)
        sys = StandardMap(fl_K(Float64, label))
        x0 = [fl_cast(Float64, s) for s in fl_starts(StableRNG(seed), m)]
        W0 = [fl_frame(seed, j, Float64) for j in 1:m]
        w = collect(20:200)
        run = fl_run((FLI(),), sys, x0, W0, 1_000, 1; steps = vcat(w, [1_000]))
        l1 = fl_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        slopes = [fl_slope(Float64.(w), Float64[s.fli for s in r[1].snaps[1:181]]) for r in run]
        X = [(r[1].snaps[1].sum[1] + r[1].snaps[1].comp[1]) / 10_000 for r in l1]
        @test all(r -> length(r[1].snaps) == 182 && r[1].snaps[182].n == 1_000, run)
        @test fl_agree(slopes, X; allowance = 4 / 180)
        if label !== :two_pi
            @test fl_agree([r[1].snaps[182].fli1 / 1_000 for r in run], X; allowance = 4 / 1_000)
            @test fl_agree([r[1].snaps[182].fli / 1_000 for r in run], X; allowance = 4 / 1_000)
        end
        hit = [(h = findfirst(s -> s.fli >= log(1e20), r[1].snaps); h === nothing ? nothing : r[1].snaps[h].n) for r in run]
        @info "AC9b reported (R3, R4)" label mean_slope = fl_mean(slopes) mean_X = fl_mean(X) sanity = sanity[label] first_handoff_ln1e20_median = fl_median([n for n in hit if n !== nothing]) never_reached = count(isnothing, hit)
    end
end

@testset "AC10 tau dependence" begin
    _ = (FLI,)
    # Same start and frame, tau in {1, 2, 5, 10}, K in {0.5 (regular class), 2 pi}, n in {100, 1e3}, 20 starts, both forms. (a) the sums at a common
    # n: |S_i^tau(n) - S_i^1(n)| <= B_S^tau(n) + B_S^1(n), each 8 eps Psi_n^(i) (L1 criterion 6, Amendment D35); (b) FLI_tau(n) <= FLI_1(n) + B^tau + B^1 and
    # FLI_1(n) - FLI_tau(n) <= floor(tau/2) ln sigma(K) + B^tau + B^1, B the two-point B_FLI of the run and sigma(K) of L1 E14. Reported: the distribution of FLI_1 - FLI_tau.
    for T in FL_TYPES, (label, N, m) in fl_panel(T)
        label === :ten && continue
        epsT = fl_eps(T)
        lns = log(fl_sigma(Float64(fl_K(T, label))))
        O = fl_orbit_ref(label, T, m, N)
        base = Dict((r.j, r.n) => r for r in fl_rows(label, T, 1, m, N))
        for tau in (2, 5, 10)
            rows = [r for r in fl_rows(label, T, tau, m, N) if r.n in (100, 1_000)]
            @test length(rows) == 2m
            gap = Float64[]
            for r in rows
                r1 = base[(r.j, r.n)]
                o = O.per[r.j]
                @test abs(r.S1 - r1.S1) <= 16 * epsT * fl_Psi(o, 1, r.n)                          # (a)
                @test abs(r.S2 - r1.S2) <= 16 * epsT * fl_Psi(o, 2, r.n)
                for (a, a1, B, B1) in ((r.fli, r1.fli, r.BB, r1.BB), (r.fli1, r1.fli1, r.B1, r1.B1))
                    @test a <= a1 + (B + B1)                                                        # (b)
                    @test a1 - a <= (tau ÷ 2) * lns + (B + B1)
                end
                push!(gap, Float64(r1.fli1 - r.fli1))
            end
            @info "AC10 reported: FLI_1 - FLI_tau (median, max)" label T tau gap = (fl_median(gap), maximum(gap))
        end
    end
end

@testset "AC11 precision ladder and representation" begin
    _ = (FLI, fli_unlogged)
    # (a) K = 2 pi and 5, 1000 starts, n = 1e3, FLI_B/n, Float64 vs Double64 on the same initial points: ensemble means within 4 combined standard errors (L1 S5);
    # K in {5, 10}, 200 starts, n = 1e4: same rule. (b) K = 2 pi, 50 starts, n up to 60, BigFloat(256) on its own orbit against the BigFloat
    # reference at its realized points (computed at 512 bits, so the reference is not of the precision it checks): within B_FLI. (c) BigFloat(256) vs Double64,
    # n = 300, 100 starts, ensemble level, same rule. (d) K = 10, n = 1e4, Float64: FLI_B is finite, equals the reference within B_FLI, exceeds 709.78 (so the unlogged
    # exp overflows) and is not clipped; the display conversion reports the overflow instead of an infinity.
    for label in (:two_pi, :five)
        @test fl_agree(fl_ens(Float64, label, 1000, 1_000), fl_ens(Double64, label, 1000, 1_000))
    end
    for label in (:five, :ten)
        @test fl_agree(fl_ens(Float64, label, 200, 10_000), fl_ens(Double64, label, 200, 10_000))
    end
    setprecision(BigFloat, FL_BITS) do
        rows = fl_rows(:two_pi, BigFloat, 1, 50, 60; bits = 512)
        @test length(rows) == 50 * length(fl_checkpoints(60))
        @test maximum(r.ratio1 for r in rows) <= 1
        @test maximum(r.ratioB for r in rows) <= 1
        @info "AC11(b) realized ratio max |error| / B_FLI (BigFloat)" worst1 = maximum(r.ratio1 for r in rows) worstB = maximum(r.ratioB for r in rows)
        for label in (:two_pi, :five)
            @test fl_agree(fl_ens(BigFloat, label, 100, 300), fl_ens(Double64, label, 100, 300))
        end
    end
    # (d)
    O = fl_orbit_ref(:ten, Float64, 5, 10_000)
    c = O.case
    runs = fl_run((FLI(),), c.sys, c.x0, c.W0, 10_000, 1; steps = [10_000])
    for j in 1:5
        s = runs[j][1].snaps[1]
        r = only(fl_ref_at(O.per[j], 1, [10_000]))
        @test isfinite(s.fli) && s.fli > 709.78
        @test fl_err(s.fli, r.FB) <= fl_BFLI(fl_eps(Float64), O.per[j], s.colB, s.stepB, r.cB, r.jB)       # not clipped, equal to the reference
        @test fli_unlogged(s.fli) === nothing                                                              # the display says "overflow", never Inf
    end
    @test fli_unlogged(1.0) ≈ exp(1.0) && fli_unlogged(700.0) ≈ exp(700.0) && fli_unlogged(0.0) == 1.0
    @test fli_unlogged(709.0) isa Float64 && fli_unlogged(710.0) === nothing
end
