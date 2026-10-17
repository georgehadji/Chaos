# AC10 predicates T, C, S on synthetic screening records (E11). The thresholds here are the fixture's own; every one is tested
# at its value and one least-significant step below and above. Expected outcomes are hand computed from the E11 comparison
# directions: T: lambda <= theta_T, S(N/10) > 0, g < g_T (strict), d >= r_T; C: lambda >= theta_C, d >= r_C; S: lambda >= theta_S,
# e_med > 0, some run of consecutive blocks with e_i < f_S e_med and length (blocks x w) >= L_S.

@testset "AC10 predicates" begin
    thT, gT, rT, thC, rC, thS, w, fS, LS = 0.01, 5.0, 0.05, 0.05, 0.3, 0.02, 100, 0.2, 500
    T(lam, st, sn, d) = EX.pred_T(lam, st, sn, d; theta_T = thT, g_T = gT, r_T = rT)
    C(lam, d) = EX.pred_C(lam, d; theta_C = thC, r_C = rC)
    S(lam, blocks) = EX.pred_S(lam, blocks; theta_S = thS, w = w, f_S = fS, L_S = LS)
    up(x) = nextfloat(x)
    dn(x) = prevfloat(x)

    @testset "T" begin
        @test T(thT, 1.0, 2.0, rT)                    # every threshold at its value qualifies where the rule is non-strict
        @test T(dn(thT), 1.0, 2.0, rT)
        @test !T(up(thT), 1.0, 2.0, rT)
        @test T(thT, 1.0, 2.0, up(rT))
        @test !T(thT, 1.0, 2.0, dn(rT))
        @test !T(thT, 1.0, gT, rT)                    # g = S(N) / S(N/10) = g_T: the block test is strict
        @test T(thT, 1.0, dn(gT), rT)
        @test !T(thT, 1.0, up(gT), rT)
        @test !T(thT, 0.0, 2.0, rT)                   # S(N/10) must be positive
        @test !T(thT, -1.0, -2.0, rT)
        @test T(thT, floatmin(Float64), 2 * floatmin(Float64), rT)
        @test T(-1.0, 1.0, 1.5, 1.0)                  # negative lambda (finite-time noise) is still <= theta_T
    end

    @testset "C" begin
        @test C(thC, rC)
        @test C(up(thC), rC)
        @test !C(dn(thC), rC)
        @test C(thC, up(rC))
        @test !C(thC, dn(rC))
        @test !C(thC, 0.0)
        @test C(1.0, pi)
    end

    # blocks with a run of nlow low values (0.0 < 0.2 e_med) in a background of ones
    blk(nlow; nhigh = 15, low = 0.0, gap = false) = begin
        v = fill(1.0, nhigh)
        run = fill(low, nlow)
        gap ? vcat(run[1:(nlow ÷ 2)], [1.0], run[(nlow ÷ 2 + 1):end], v) : vcat(v[1:5], run, v[6:end])
    end

    @testset "S" begin
        @test S(thS, blk(5))                          # run of exactly L_S / w = 5 blocks: length L_S qualifies
        @test !S(thS, blk(4))                         # L_S - w does not
        @test S(thS, blk(6))
        @test !S(thS, blk(6; gap = true))             # 3 + 1 high + 3: no contiguous stretch of 5
        @test S(up(thS), blk(5))
        @test S(thS, blk(5))
        @test !S(dn(thS), blk(5))                     # lambda >= theta_S, one step below fails
        @test !S(thS, blk(5; low = fS))               # a block exactly at f_S e_med (e_med = 1) is not below
        @test S(thS, blk(5; low = dn(fS)))            # one step below is
        @test !S(thS, blk(5; low = up(fS)))
        @test S(thS, vcat(fill(5.0, 15), fill(0.99, 5)))      # e_med = 5, stuck below 1.0: 0.99 < 0.2 * 5
        @test !S(thS, vcat(fill(5.0, 15), fill(1.0, 5)))      # exactly 0.2 * 5 is not below
        @test !S(thS, fill(0.0, 20))                  # e_med = 0: the comparison is meaningless, not S
        @test !S(thS, fill(-1.0, 20))                 # e_med < 0
        @test !S(thS, fill(1.0, 20))                  # no stretch at all
    end

    @testset "median and stuck stretch" begin
        @test EX.block_median([3.0, 1.0, 2.0]) == 2.0                  # odd count
        @test EX.block_median([4.0, 1.0, 3.0, 2.0]) == 2.5             # even count: mean of the two middle values
        @test EX.block_median([7.0]) == 7.0
        @test EX.block_median([2.0, 2.0, 2.0, 9.0]) == 2.0
        @test EX.stuck_length(blk(5), fS, w) == 5 * w
        @test EX.stuck_length(blk(6; gap = true), fS, w) == 3 * w
        @test EX.stuck_length(fill(1.0, 12), fS, w) == 0
        # two separate stretches: the longest counts
        v = vcat(fill(1.0, 6), zeros(2), fill(1.0, 6), zeros(4), fill(1.0, 6))
        @test EX.stuck_length(v, fS, w) == 4 * w
    end

    @testset "parabolic periodic points count as elliptic for C (distance to the elliptic subset)" begin
        # K = 2: the period-2 family A {(0, pi), (pi, pi)} has trace 2 - K^2 = -2, parabolic; (pi, 0) is elliptic, (0, 0) hyperbolic
        ps = ex_pset(2.0)
        z = (0.0, pi + 0.1)                                            # 0.1 from the parabolic point (0, pi), far from (pi, 0)
        d_ell = EX.periodic_distance(z, ps; elliptic_only = true)
        @test abs(d_ell - 0.1) <= 1e-9
        @test EX.periodic_distance(z, ps) <= d_ell                     # the full set is at least as close
        # the hyperbolic fixed point (0, 0) is not in the elliptic subset
        z0 = (0.01, 0.0)
        @test EX.periodic_distance(z0, ps) <= 0.0101
        @test EX.periodic_distance(z0, ps; elliptic_only = true) > 0.5
    end
end
