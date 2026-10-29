# D89e budget rule, D89f verdicts R1, D89g verdicts R2, D89y pilot rule check (O8). Pure functions of the amendment, no store.
# Spec: analysis-spec.md, Amendment D89, E26, E27, E31, O8, criteria D89e, D89f, D89g, D89y. API: test plan, "API proposed by the tests".

# workspaces/study/stages/03-run/output/p2-indicators-pilot-evidence.md, section 3 (n = 1e5; E1, E2, WBA at W = 51200). Columns: k, orbit, lambda_1, r_LE,
# r_SALI, r_FLI, r_MEGNO, SALI, RLI, D_E1, D_E2, dig_min, and the tabulated R1 verdicts LE SALI FLI MEGNO RLI E1 E2 WBA (C chaotic, r regular).
const P2_PILOT = [
    (1, 3188, 1.12e-04, 1.26, 1.26, 1.26, 1.00, 2.00e-10, 0.00e+00, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (1, 3189, 1.12e-04, 1.26, 1.26, 1.26, 1.00, 2.00e-10, 0.00e+00, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (1, 3190, 1.12e-04, 1.26, 1.26, 1.26, 1.00, 2.00e-10, 0.00e+00, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (1, 3191, 1.12e-04, 1.26, 1.26, 1.26, 1.00, 2.00e-10, 0.00e+00, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (2, 3192, 7.05e-05, 1.50, 1.52, 1.48, 1.00, 1.18e-06, 4.78e-15, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (2, 3193, 6.75e-05, 2.56, 2.82, 1.44, 1.14, 2.99e-06, 1.56e-13, -14.15, -14.15, 14.15, "rrrrrrrr"),
    (2, 3194, 1.13e-04, 1.23, 1.24, 1.25, 1.00, 2.64e-10, 1.20e-13, -14.15, -14.15, 13.76, "rrrrrrrr"),
    (2, 3195, 8.90e-05, 1.32, 1.29, 1.28, 1.00, 5.91e-09, 7.71e-14, -14.15, -14.15, 13.17, "rrrrrrrr"),
    (3, 3196, 1.34e-01, 13.36, 13.36, 13.36, 13.02, 0.00e+00, 1.55e-02, -1.08, -1.74, 1.03, "CCCCCCCC"),
    (3, 3197, 9.82e-05, 1.30, 1.31, 1.29, 1.00, 4.36e-09, 6.07e-14, -14.15, -14.15, 13.11, "rrrrrrrr"),
    (3, 3198, 1.10e-04, 1.25, 1.37, 1.25, 1.00, 8.36e-08, 1.23e-13, -14.15, -14.15, 13.32, "rrrrrrrr"),
    (3, 3199, 4.96e-03, 6.41, 6.56, 6.35, 6.21, 0.00e+00, 2.50e-03, -9.95, -9.86, 3.60, "CCCCCrrC"),
    (4, 3200, 2.65e-01, 9.05, 9.05, 9.04, 9.26, 0.00e+00, 6.52e-02, -0.48, -0.40, 0.62, "CCCCCCCC"),
    (4, 3201, 4.59e-03, 8.35, 8.34, 8.29, 8.79, 0.00e+00, 7.71e-04, -8.02, -7.75, 2.42, "CCCCCrrC"),
    (4, 3202, 3.04e-01, 9.76, 9.76, 9.76, 10.07, 0.00e+00, 6.24e-03, -3.79, -0.71, 1.74, "CCCCCCCC"),
    (4, 3203, 3.05e-01, 11.35, 11.35, 11.34, 12.64, 0.00e+00, 1.75e-02, -3.95, -1.45, 1.79, "CCCCCCCC"),
    (5, 3204, 1.17e+00, 10.09, 10.09, 10.09, 10.09, 0.00e+00, 4.54e-03, -0.66, -1.45, 1.81, "CCCCCCCC"),
    (5, 3205, 1.17e+00, 9.88, 9.88, 9.88, 9.82, 0.00e+00, 7.44e-03, -1.52, -1.05, 1.88, "CCCCCCCC"),
    (5, 3206, 1.18e+00, 9.87, 9.87, 9.87, 9.89, 0.00e+00, 1.97e-02, -0.42, -0.48, 1.77, "CCCCCCCC"),
    (5, 3207, 1.17e+00, 9.97, 9.97, 9.97, 9.94, 0.00e+00, 5.95e-03, -0.67, -0.31, 2.04, "CCCCCCCC"),
    (6, 3208, 1.63e+00, 9.99, 9.99, 9.99, 9.93, 0.00e+00, 7.12e-03, -1.93, -0.98, 2.11, "CCCCCCCC"),
    (6, 3209, 1.62e+00, 10.00, 10.01, 10.00, 9.95, 0.00e+00, 1.27e-02, -0.61, -0.49, 2.21, "CCCCCCCC"),
    (6, 3210, 1.61e+00, 10.04, 10.04, 10.04, 10.05, 0.00e+00, 4.31e-03, -0.68, -0.46, 2.06, "CCCCCCCC"),
    (6, 3211, 1.61e+00, 10.00, 10.00, 10.00, 9.97, 0.00e+00, 6.81e-03, -0.31, -0.63, 2.11, "CCCCCCCC")]

@testset "D89e budget rule: n(b), n-prime(n), ties, windows" begin
    bud = [200 * 2^k for k in 0:9]
    # the two real grids of E26, built by the grid builder of the runner (20 points per decade, every step to 100, plus N / 10 and N / 2)
    big = [200, 447, 891, 1778, 3548, 7079, 14125, 28184, 56234, 112202]
    for (N, nb, npr) in ((10^7, big, [20, 45, 89, 178, 355, 708, 1413, 2818, 5623, 11220]), (10^6, big, [20, 45, 89, 178, 355, 708, 1413, 2818, 5623, 11220]),
                         (102400, vcat(big[1:8], [51200, 102400]), [20, 45, 89, 178, 355, 708, 1413, 2818, 5012, 10240]))
        g = AX.record_grid(N, 1; points_per_decade = 20, first_steps = 100, extra_fractions = [0.1, 0.5])
        @test [AN.p2_budget_step(g, b) for b in bud] == nb
        @test [AN.p2_nearest_step(g, n / 10) for n in nb] == npr
        @test all(AN.p2_budget_step(g, b) <= 1.13 * b for b in bud)                    # P15
    end
    # the fixture grids (hand table of the plan): n(b_k) = 44 x 2^k on GA, n' = 4, 9, 18, 35, ...
    @test [AN.p2_budget_step(P2_GA, b) for b in P2_BUDGETS] == P2_NK && [AN.p2_nearest_step(P2_GA, n / 10) for n in P2_NK] == P2_NP
    @test [AN.p2_budget_step(P2_GB, b) for b in P2_BUDGETS] == P2_NKB
    # a budget exactly on the grid maps to itself; one more maps to the next step
    @test AN.p2_budget_step([1, 5, 10], 5) == 5 && AN.p2_budget_step([1, 5, 10], 6) == 10
    # a tie goes to the larger step
    @test AN.p2_nearest_step([1, 3], 2) == 3 && AN.p2_nearest_step([10, 20], 15) == 20 && AN.p2_nearest_step([2, 4, 8], 3.0) == 4 && AN.p2_nearest_step([2, 4, 8], 2.9) == 2
    @test AN.p2_nearest_step([5, 6, 7], 6.4) == 6 && AN.p2_nearest_step([5, 6, 7], 100) == 7
    # a budget above the horizon is refused
    @test an_rejects(() -> AN.p2_budget_step([1, 5, 10], 11), ["budget"])
    # the window row of b is W = b / 2, which must be an element of the windows
    @test [AN.p2_window_row(P2_WINDOWS, b) for b in P2_BUDGETS] == collect(1:10)
    @test an_rejects(() -> AN.p2_window_row(P2_WINDOWS, 30), ["window"]) && an_rejects(() -> AN.p2_window_row(P2_WINDOWS, 41), ["window"])
    @test an_rejects(() -> AN.p2_window_row(P2_WINDOWS, 40960), ["window"])
end

@testset "D89f growth ratio and verdict R1: boundaries and undecided causes" begin
    p = p2_params()
    vd(X, v) = AN.p2_verdict(X, :R1, v, p)
    # g_LE = n lambda, g_SALI = -L, g_FLI = fli, g_MEGNO = Ybar20_1; the ratio is undefined unless both are finite and g(n') is positive
    @test AN.p2_growth(:LE, 44, 0.5) == 22.0 && AN.p2_growth(:SALI, 44, -2.0) == 2.0 && AN.p2_growth(:FLI, 44, 3.0) == 3.0 && AN.p2_growth(:MEGNO, 44, 2.0) == 2.0
    @test AN.p2_ratio(10.0, 2.0) == 5.0 && AN.p2_ratio(-3.0, 2.0) == -1.5
    @test AN.p2_ratio(10.0, 0.0) === nothing && AN.p2_ratio(10.0, -1.0) === nothing && AN.p2_ratio(NaN, 1.0) === nothing && AN.p2_ratio(1.0, NaN) === nothing && AN.p2_ratio(Inf, 1.0) === nothing
    # a positive SALI.L makes g_SALI negative: the ratio is undefined and the verdict undecided
    gp = AN.p2_growth(:SALI, 4, 0.5)
    @test gp == -0.5 && AN.p2_ratio(AN.p2_growth(:SALI, 44, -3.0), gp) === nothing
    # LE: lambda_1 >= le_floor and r >= ratio_threshold (>= as written)
    le(lam, r) = vd(:LE, (n = 100, lambda = lam, r = r))
    @test le(1e-4, 5.0) === :chaotic && le(prevfloat(1e-4), 5.0) === :regular && le(1e-4, prevfloat(5.0)) === :regular && le(nextfloat(1e-4), nextfloat(5.0)) === :chaotic
    @test le(0.3, -3.0) === :regular && le(0.3, 0.0) === :regular && le(-0.1, 9.0) === :regular
    @test le(0.3, nothing) === :undecided && le(NaN, 9.0) === :undecided && le(Inf, 9.0) === :undecided
    # SALI, FLI, MEGNO: r >= ratio_threshold
    for (X, extra) in ((:SALI, (sali = 0.1,)), (:FLI, (n = 100, fli = 3.0)), (:MEGNO, (ybar = 2.0,)))
        v(r) = vd(X, merge((r = r,), extra))
        @test v(5.0) === :chaotic && v(prevfloat(5.0)) === :regular && v(nextfloat(5.0)) === :chaotic && v(-2.0) === :regular && v(nothing) === :undecided
    end
    # RLI >= rli_threshold, undecided if not finite
    rli(x) = vd(:RLI, (rli = x,))
    @test rli(1e-7) === :chaotic && rli(prevfloat(1e-7)) === :regular && rli(nextfloat(1e-7)) === :chaotic && rli(0.0) === :regular && rli(NaN) === :undecided && rli(Inf) === :undecided
    # E1 and E2: D >= -4; the second flag name (degenerate, undefined window: bit value 2) or a non-finite D is undecided; other bits change nothing
    for X in (:E1, :E2)
        d(D, w = 0) = vd(X, (D = D, flags = w))
        @test d(-4.0) === :chaotic && d(prevfloat(-4.0)) === :regular && d(nextfloat(-4.0)) === :chaotic && d(0.0) === :chaotic && d(-14.15) === :regular
        @test d(NaN) === :undecided && d(-Inf) === :undecided && d(-1.0, 2) === :undecided && d(-9.0, 2) === :undecided && d(-1.0, 3) === :undecided && d(-1.0, 2 + 512) === :undecided
        for w in (1, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1 + 512 + 1024)           # below_floor, failed/near_zero/mirror/tie/refined: data, no verdict change
            @test d(-1.0, w) === :chaotic && d(-9.0, w) === :regular
        end
    end
    # WBA: dig_min = minimum over the four rows < wba_digits_threshold; the flag zero (bit value 2) on any one row, or a non-finite dig, is undecided
    wba(dig, w = [0, 0, 0, 0]) = vd(:WBA, (dig = dig, flags = w))
    @test wba([6.0, 7.0, 8.0, 9.0]) === :regular && wba([prevfloat(6.0), 7.0, 8.0, 9.0]) === :chaotic && wba([9.0, 8.0, 7.0, nextfloat(6.0)]) === :regular && wba([9.0, 3.0, 8.0, 9.0]) === :chaotic
    for r in 1:4
        w = zeros(Int, 4); w[r] = 2
        @test wba([3.0, 9.0, 8.0, 7.0], w) === :undecided && wba([12.0, 9.0, 8.0, 7.0], w) === :undecided
        nan = [12.0, 9.0, 8.0, 7.0]; nan[r] = NaN
        @test wba(nan) === :undecided
    end
    for w in (1, 4, 1 + 4)
        @test wba([12.0, 9.0, 8.0, 7.0], [w, 0, w, 0]) === :regular && wba([3.0, 9.0, 8.0, 7.0], [0, w, 0, 0]) === :chaotic
    end
    # STR has no verdict: it enters only the scores (E31)
    @test_throws ArgumentError AN.p2_verdict(:STR, :R1, (D = 0.1,), p)
    @test_throws ArgumentError AN.p2_verdict(:LE, :R3, (n = 100, lambda = 0.3, r = 9.0), p)
end

@testset "D89g verdict R2: boundaries, non-finite, the other four rules equal R1" begin
    p = p2_params()
    vd(X, v) = AN.p2_verdict(X, :R2, v, p)
    # LE: n lambda_1(n) >= r2_le_factor ln n (n a power of two keeps n * lambda exact); FLI: fli(n) >= factor ln n
    for n in (64, 1024, 16384)
        thr = p.r2_le_factor * log(n)
        @test vd(:LE, (n = n, lambda = thr / n, r = 1.0)) === :chaotic && vd(:LE, (n = n, lambda = prevfloat(thr / n), r = 9.0)) === :regular
        @test vd(:LE, (n = n, lambda = nextfloat(thr / n), r = nothing)) === :chaotic                     # the ratio is not used by R2
        thf = p.r2_fli_factor * log(n)
        @test vd(:FLI, (n = n, r = nothing, fli = thf)) === :chaotic && vd(:FLI, (n = n, r = 9.0, fli = prevfloat(thf))) === :regular
    end
    # SALI(n) <= r2_sali_threshold, Ybar20_1(n) >= r2_megno_threshold
    @test vd(:SALI, (r = 9.0, sali = 1e-8)) === :chaotic && vd(:SALI, (r = 9.0, sali = nextfloat(1e-8))) === :regular && vd(:SALI, (r = nothing, sali = 0.0)) === :chaotic
    @test vd(:MEGNO, (r = 1.0, ybar = 2.5)) === :chaotic && vd(:MEGNO, (r = 9.0, ybar = prevfloat(2.5))) === :regular && vd(:MEGNO, (r = nothing, ybar = 3.0)) === :chaotic
    # a non-finite value is undecided
    @test vd(:LE, (n = 64, lambda = NaN, r = 9.0)) === :undecided && vd(:FLI, (n = 64, r = 9.0, fli = NaN)) === :undecided && vd(:SALI, (r = 9.0, sali = Inf)) === :undecided &&
          vd(:MEGNO, (r = 9.0, ybar = NaN)) === :undecided
    # E1, E2, WBA and RLI keep their R1 rules
    cases = [(:E1, (D = -4.0, flags = 0)), (:E1, (D = prevfloat(-4.0), flags = 0)), (:E1, (D = -1.0, flags = 2)), (:E2, (D = -4.0, flags = 0)), (:E2, (D = NaN, flags = 2)),
             (:WBA, (dig = [5.9, 7.0, 8.0, 9.0], flags = zeros(Int, 4))), (:WBA, (dig = [6.0, 7.0, 8.0, 9.0], flags = zeros(Int, 4))), (:WBA, (dig = [9.0, 7.0, 8.0, 9.0], flags = [0, 2, 0, 0])),
             (:RLI, (rli = 1e-7,)), (:RLI, (rli = prevfloat(1e-7),)), (:RLI, (rli = NaN,))]
    @test all(c -> AN.p2_verdict(c[1], :R2, c[2], p) === AN.p2_verdict(c[1], :R1, c[2], p), cases)
end

@testset "D89 chaos scores (E31): orientation and availability" begin
    sc(X, v) = AN.p2_score(X, v)
    @test sc(:LE, (n = 100, lambda = 0.3, r = 9.5)) == 9.5 && sc(:SALI, (r = 3.0, sali = 0.1)) == 3.0 && sc(:FLI, (n = 100, r = 7.0, fli = 3.0)) == 7.0 && sc(:MEGNO, (r = 2.0, ybar = 2.0)) == 2.0
    @test sc(:LE, (n = 100, lambda = 0.3, r = nothing)) === nothing
    @test sc(:RLI, (rli = 1e-3,)) == log10(1e-3) && sc(:RLI, (rli = 0.0,)) == -Inf && sc(:RLI, (rli = NaN,)) === nothing
    @test sc(:E1, (D = -3.5, flags = 0)) == -3.5 && sc(:E2, (D = -3.5, flags = 1)) == -3.5 && sc(:E1, (D = -3.5, flags = 2)) === nothing && sc(:E2, (D = NaN, flags = 2)) === nothing
    @test sc(:WBA, (dig = [9.0, 3.0, 8.0, 7.0], flags = [0, 0, 0, 0])) == -3.0 && sc(:WBA, (dig = [9.0, 3.0, 8.0, 7.0], flags = [0, 2, 0, 0])) === nothing
    @test sc(:STR, (D = 0.01,)) == -log10(0.01) && sc(:STR, (D = 0.0,)) === nothing && sc(:STR, (D = -1.0,)) === nothing && sc(:STR, (D = NaN,)) === nothing
end

@testset "D89y pilot rule check (O8): the 192 R1 verdicts of the 24 tabulated starts" begin
    p = p2_params()
    n = 100_000
    bad = 0; total = 0
    for (k, orbit, lam, rle, rsali, rfli, rmeg, sali, rli, de1, de2, dig, want) in P2_PILOT
        got = [AN.p2_verdict(:LE, :R1, (n = n, lambda = lam, r = rle), p), AN.p2_verdict(:SALI, :R1, (r = rsali, sali = sali), p),
               AN.p2_verdict(:FLI, :R1, (n = n, r = rfli, fli = 1.0), p), AN.p2_verdict(:MEGNO, :R1, (r = rmeg, ybar = 2.0), p), AN.p2_verdict(:RLI, :R1, (rli = rli,), p),
               AN.p2_verdict(:E1, :R1, (D = de1, flags = 0), p), AN.p2_verdict(:E2, :R1, (D = de2, flags = 0), p),
               AN.p2_verdict(:WBA, :R1, (dig = [dig, 14.15, 14.15, 14.15], flags = zeros(Int, 4)), p)]
        exp = [c == 'C' ? :chaotic : :regular for c in want]
        bad += count(got .!= exp); total += 8
    end
    @test total == 192 && bad == 0 && length(P2_PILOT) == 24
end
