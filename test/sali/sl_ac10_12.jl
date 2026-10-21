# AC10 to AC12: chaotic rate against lambda_1 (ensemble), regular decay law, precision ladder. Ensemble-level and loose-interval criteria (K18):
# the intervals and allowances are the spec's proposals (P6), not source values.

@testset "AC10 chaotic rate against lambda_1" begin
    _ = (SALI,)
    # K in {2 pi, 5, 10}, Float64 (1000 starts): the least-squares slope m of -L(n) against n over n in [20, 200] (log form) and the L1 series
    # X(1e4) of a k = 1 run on the identical orbit (first seeded column). Pass: |mean(m/2) - mean(X)| <= 4 combined standard errors (L1 S5) + 2/180
    # (finite-window allowance, derived in the spec, to be replaced by the realized values at stage 04). The per-orbit values are reported.
    m = 1000
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = sl_seed(10, ci)
        sys = StandardMap(sl_K(Float64, label))
        x0 = [sl_cast(Float64, s) for s in sl_starts(StableRNG(seed), m)]
        W0 = [sl_frame(seed, j, Float64) for j in 1:m]
        run = sl_run(SALI(), sys, x0, W0, 200, 1; steps = 20:200)
        l1 = sl_run(Benettin(), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        half_slopes = [sl_slope(Float64[s.n for s in r.snaps], Float64[-s.L for s in r.snaps]) / 2 for r in run]
        X = [(r.snaps[1].sum[1] + r.snaps[1].comp[1]) / 10_000 for r in l1]
        @test all(r -> length(r.snaps) == 181, run)
        @test sl_agree(half_slopes, X; allowance = 2 / 180)
        @info "AC10 reported" label mean_half_slope = sl_mean(half_slopes) mean_X = sl_mean(X) sd_half_slope = sqrt(sl_var(half_slopes))
    end
end

@testset "AC11 regular decay law" begin
    _ = (SALI,)
    # 20 regular-class starts at K = 0.5 (S7), Float64, log form, n up to 1e6, seeded generic frames: the least-squares slope of L(n) against ln n over
    # 21 log-spaced n in [1e4, 1e6] lies in [-2.4, -1.6] for at least 80 percent of the starts. Starts still in the transient at n = 1e4 (SALI > 1e-2) are
    # reported and at most 20 percent may be excluded (they count as not passing). Loose proposal (K18): the interval is derived from A9.
    m = 20
    c = sl_case(:half, Float64, m)
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000_000, 1; steps = ns)
    transient, passed = 0, 0
    slopes = Float64[]
    for r in run
        @test length(r.snaps) == length(ns)
        slope = sl_slope(log.(Float64.(ns)), Float64[s.L for s in r.snaps])
        push!(slopes, slope)
        if r.snaps[1].sali > 1e-2
            transient += 1
        elseif -2.4 <= slope <= -1.6
            passed += 1
        end
    end
    @info "AC11 reported: slopes of L against ln n, and starts in the transient at n = 1e4" slopes transient
    @test transient <= 0.2 * m
    @test passed >= ceil(Int, 0.8 * m)
end

# Ensemble of -L(n)/(2n) in type T on m starts of a class (the same Float64 starts widened exactly), through the per-orbit body with the bare reducer.
function sl_ens(::Type{T}, label::Symbol, m::Int, n::Int) where {T}
    c = sl_case(label, T, m)
    return map(1:m) do j
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], n, 1, (SALI(),); orbit = j)
        -Float64(reducer_final(SALI(), st[1]).L) / (2n)
    end
end

@testset "AC12 precision ladder" begin
    _ = (SALI, sali_floor)
    # (a) floor scaling: K in {5, 2 pi}, 100 starts, n = 150, direct form Double64 vs Float64 on the same initial points: n_floor(T) = first handoff with
    # SALI < B_meas(T) (or SALI = 0); the median over starts of n_floor(Double64) / n_floor(Float64) lies in [1.5, 2.4] (A14: ln(1/(10 B_meas)) is 32.4 and 65.6).
    # (b) K = 2 pi, 50 starts (50 starts), n = 60, BigFloat(256) direct form on its own orbit against the BigFloat reference at its realized points
    # (computed at 512 bits, so that the reference is not of the precision it checks): within B_SALI(n). (c) log form: ensemble means of -L(n)/(2n) at
    # n = 1e3 for K = 2 pi and 5, Float64 vs Double64 on the same initial points (1000 starts), and BigFloat(256) vs Double64 at
    # n = 300 (100 starts), within 4 combined standard errors (L1 S5; ensemble level only, the orbits decorrelate).
    for label in (:five, :two_pi)
        nfl = Dict{DataType, Vector{Int}}()
        for T in SL_TYPES
            c = sl_case(label, T, 100)
            run = sl_run(SALI(), c.sys, c.x0, c.W0, 150, 1; all = true)
            nfl[T] = map(run) do r
                i = findfirst(s -> s.sali < sl_meas(T), r.snaps)
                i === nothing ? 0 : r.snaps[i].n
            end
        end
        @test all(>(0), nfl[Float64]) && all(>(0), nfl[Double64])
        ratios = nfl[Double64] ./ nfl[Float64]
        med = sl_median(ratios)
        @info "AC12(a) reported" label median_ratio = med median_nfloor_Float64 = sl_median(nfl[Float64]) median_nfloor_Double64 = sl_median(nfl[Double64])
        @test 1.5 <= med <= 2.4
    end
    setprecision(BigFloat, SL_BITS) do
        m = 50
        ns = [10, 20, 40, 60]
        c = sl_case(:two_pi, BigFloat, m)
        epsT = sl_eps(BigFloat)
        run = sl_run(SALI(), c.sys, c.x0, c.W0, 60, 1; steps = ns)
        nres = 0
        for j in 1:m
            pts, _ = lb_pts(c.sys, c.x0[j], 60)
            ref = sl_ref(sl_mats_std(c.K, pts), c.W0[j], 60; keep = ns, bits = 512)
            Ms = lb_tangents(c.sys, c.x0[j], 60)
            for s in run[j].snaps
                at = sl_at(ref, s.n)
                p1, p2 = lb_Psi(view(Ms, 1:s.n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:s.n), c.W0[j][:, 2])
                saliref = Float64(at.S)
                nres += saliref >= 10 * 4 * epsT
                @test sl_err(s.sali, at.S; bits = 512) <= sl_BSALI(epsT, sl_Lambda(ref.lnsali, s.n), saliref, p1, p2)
            end
        end
        @test nres > 0
    end
    for label in (:two_pi, :five)
        m = 1000
        @test sl_agree(sl_ens(Float64, label, m, 1_000), sl_ens(Double64, label, m, 1_000))
    end
    setprecision(BigFloat, SL_BITS) do
        for label in (:two_pi, :five)
            m = 100
            @test sl_agree(sl_ens(BigFloat, label, m, 300), sl_ens(Double64, label, m, 300))
        end
    end
end
