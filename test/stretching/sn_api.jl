# API tests: the constructor and its parameters, the reducer protocol on hand arithmetic, the pure functions (spectrum_cell, spectral_distance, spectrum_density, stretching_range), the
# workspace surface of the Tangent layer (additive: reducer_workspace, the six-argument reducer_step and reducer_update, the three-argument reducer_final, reducer_bytes, the `workspace`
# keyword of propagate_orbit). AC14 to AC16 are stage 04 (see the test plan). Every testset starts with a line naming the new symbols it uses (UndefVarError when missing).

@testset "API constructor and parameters" begin
    _ = (StretchingNumbers,)
    # StretchingNumbers(; range, bins = 1000, interval = 100, min_per_cell = 100, replace = true): the range is required (the reducer cannot see the map), the defaults are the keys proposed for
    # numerics-defaults.toml (spectrum_bins, spectrum_interval, spectrum_min_per_cell, P4); the keyword values are what reducer_provenance reports; an isbits reducer.
    r = StretchingNumbers(; range = (-1.0, 1.0))
    pr = reducer_provenance(r)
    @test pr.spectrum_bins == 1000 && pr.spectrum_interval == 100 && pr.spectrum_min_per_cell == 100 && pr.replace === true
    @test pr.spectrum_lo == -1.0 && pr.spectrum_hi == 1.0 && pr.spectrum_cell_width == 2.0 / 1000
    @test isbitstype(typeof(r))
    r2 = StretchingNumbers(; range = (-2.5, 3.5), bins = 60, interval = 30, min_per_cell = 7, replace = false)
    pr2 = reducer_provenance(r2)
    @test pr2.spectrum_bins == 60 && pr2.spectrum_interval == 30 && pr2.spectrum_min_per_cell == 7 && pr2.replace === false && pr2.spectrum_lo == -2.5 && pr2.spectrum_hi == 3.5
    @test pr2.spectrum_cell_width == 6.0 / 60 && pr2.indicator == "stretching" && pr2.convention == SN_CONVENTION
    @test sn_error(() -> StretchingNumbers()) isa Exception                                        # no range, no reducer
    # a Float64 range widens exactly into Double64; a Double64 range does not narrow into Float64
    rd = StretchingNumbers(; range = (Double64(-1), Double64(1)), bins = 10, interval = 10)
    @test reducer_init(r2, Double64, Val(2), Val(2)) isa Any
    @test sn_error(() -> reducer_init(rd, Float64, Val(2), Val(2))) isa ArgumentError
end

@testset "API reducer protocol" begin
    _ = (StretchingNumbers, reducer_workspace)
    # reducer_init is total except k != 2 (ArgumentError "StretchingNumbers needs k = 2, got k = $k"); the workspace length is 2 (B + 2) unsigned integers; the n = 0 final (everything zero, the
    # scalar part also from the two-argument reducer_final); the arithmetic of three handoffs at the map steps 7, 14, 21 (tau = 7 = n / K) on the columns with norms (alpha_1, alpha_2) = (5, 13),
    # (10, 10), (5/8, 1/2): a = ln alpha, z = a / tau, the compensated sums of a and of a^2, the cell by SN4, D by SN5, the mean and the standard deviation of z; the replacer returns the
    # column-normalised matrix, the non-replacing position returns W unchanged; a failure (zero or non-finite norm) is recorded without throwing in reducer_step, sums, counts, n and W unchanged,
    # thrown by reducer_update (orbit 0, backend :cpu), and a failed state stays failed and unchanged.
    B = 20
    lo, hi = -1.0, 1.0
    for k in (1, 3)
        e = sn_error(() -> reducer_init(StretchingNumbers(; range = (lo, hi), bins = B), Float64, Val(2), Val(k)))
        @test e isa ArgumentError && occursin("StretchingNumbers needs k = 2, got k = $k", e.msg)
        @test sn_error(() -> reducer_workspace(StretchingNumbers(; range = (lo, hi), bins = B), Float64, Val(2), Val(k))) isa ArgumentError
    end
    for T in SN_TYPES, replace in (true, false)
        epsT = sn_eps(T)
        r = StretchingNumbers(; range = (T(lo), T(hi)), bins = B, interval = 14, min_per_cell = 1, replace = replace)
        @test reducer_workspace(r, T, Val(2), Val(2)) == 2 * (B + 2)
        st0 = reducer_init(r, T, Val(2), Val(2))
        ws = zeros(UInt32, 2 * (B + 2))
        f0 = reducer_final(r, st0, ws)
        @test f0.n == 0 && f0.handoffs == 0 && f0.D == 0 && size(f0.counts) == (B + 2, 2) && all(iszero, f0.counts) && f0.below == [0, 0] && f0.above == [0, 0]
        @test f0.sum isa SVector{2, T} && all(iszero, f0.sum) && all(iszero, f0.comp) && all(iszero, f0.sum2) && all(iszero, f0.comp2) && all(iszero, f0.mean) && all(iszero, f0.std)
        g0 = reducer_final(r, st0)
        @test g0.n == 0 && g0.handoffs == 0 && all(iszero, g0.sum) && all(iszero, g0.mean) && !haskey(g0, :counts)
        x = SVector{2, T}(0, 0)
        cols = (((3.0, 4.0), (5.0, 12.0)), ((6.0, 8.0), (6.0, 8.0)), ((0.375, 0.5), (0.5, 0.0)))
        alpha = ((5.0, 13.0), (10.0, 10.0), (0.625, 0.5))
        Wm(ca, cb) = SMatrix{2, 2, T}(T(ca[1]), T(ca[2]), T(cb[1]), T(cb[2]))
        st = st0
        S = [zero(BigFloat), zero(BigFloat)]
        S2 = [zero(BigFloat), zero(BigFloat)]
        zs = [BigFloat[], BigFloat[]]
        tau = 7
        for (k, ((ca, cb), (a1, a2))) in enumerate(zip(cols, alpha))
            n = tau * k
            Wk = Wm(ca, cb)
            st, W2 = reducer_step(r, st, n, x, Wk, ws)
            f = reducer_final(r, st, ws)
            setprecision(BigFloat, SN_BITS) do
                for (i, a) in enumerate((a1, a2))
                    al = log(BigFloat(a))
                    S[i] += al
                    S2[i] += al^2
                    push!(zs[i], al / tau)
                end
                @test f.n == n && f.handoffs == k
                for i in 1:2
                    @test abs(sn_big(f.sum[i]) + sn_big(f.comp[i]) - S[i]) <= 4 * epsT * max(1, abs(S[i]))
                    @test abs(sn_big(f.sum2[i]) + sn_big(f.comp2[i]) - S2[i]) <= 4 * epsT * max(1, abs(S2[i]))
                    mref = S[i] / n
                    @test abs(sn_big(f.mean[i]) - mref) <= 4 * epsT * max(abs(mref), 1e-3)
                    sdref = sqrt(max(sum(z -> z^2, zs[i]) / k - (sum(zs[i]) / k)^2, zero(BigFloat)))
                    @test abs(sn_big(f.std[i]) - sdref) <= 3 * sqrt(epsT) * (abs(Float64(mref)) + Float64(sdref)) + 4 * epsT
                    ref_counts = sn_hist_ref(zs[i], T(lo), T(hi), B)
                    @test sn_cnt(f, i) == ref_counts
                end
                Dref = sn_D_ref(sn_hist_ref(zs[1], T(lo), T(hi), B), sn_hist_ref(zs[2], T(lo), T(hi), B), T(lo), T(hi), B)
                @test Dref == 0 ? f.D == 0 : abs(sn_big(f.D) - Dref) <= 4 * epsT * Dref
                @test f.undersampled == (k < 1 * B) && f.below == [0, 0] && f.above == [0, 0]
            end
            if replace
                Wref = SMatrix{2, 2, T}(T(ca[1]) / T(a1), T(ca[2]) / T(a1), T(cb[1]) / T(a2), T(cb[2]) / T(a2))
                @test all(abs.(W2 - Wref) .<= 4 * epsT)
            else
                @test isequal(W2, Wk)
            end
        end
        # failures: a zero column 1, then a non-finite column 2, both bad (column 1 reported), an infinite column 1
        good = Wm((3, 4), (5, 12))
        wsg = zeros(UInt32, 2 * (B + 2))
        stg, _ = reducer_step(r, st0, 7, x, good, wsg)
        for (Wbad, col, val, thr) in ((Wm((0, 0), (5, 12)), 1, 0.0, 0.0), (Wm((3, 4), (T(NaN), 1)), 2, NaN, Inf), (Wm((0, 0), (0, 0)), 1, 0.0, 0.0), (Wm((T(Inf), 1), (5, 12)), 1, Inf, Inf))
            wsf = copy(wsg)
            stf, Wf = reducer_step(r, stg, 14, x, Wbad, wsf)
            @test reducer_failed(r, stf) && !reducer_failed(r, stg)
            @test isequal(Wf, Wbad) && wsf == wsg                                        # W and the counts unchanged by the failure
            @test isequal(reducer_final(r, stf, wsf), reducer_final(r, stg, wsg))        # sums, counts, n and handoffs unchanged
            err = reducer_error(r, stf, 9, :cpu)
            @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 9 && err.step == 14 && err.backend === :cpu && err.column == col
            @test (isnan(val) ? isnan(err.value) : err.value == val) && err.threshold == thr
            e = sn_error(() -> reducer_update(r, stg, 14, x, Wbad, copy(wsg)))
            @test e isa NumericalHealthError && e.orbit == 0 && e.backend === :cpu && e.step == 14 && e.column == col
            wsf2 = copy(wsf)
            stf2, Wf2 = reducer_step(r, stf, 21, x, good, wsf2)                          # a failed state stays failed and unchanged
            @test reducer_failed(r, stf2) && isequal(stf2, stf) && isequal(Wf2, good) && wsf2 == wsf
        end
    end
end

@testset "API pure functions" begin
    _ = (spectrum_cell, spectral_distance, spectrum_density, stretching_range)
    # spectrum_cell(z, lo, hi, B): the cell rule of SN4 (0 below, B + 1 above, 1..B inside, the last cell closed at hi); nondecreasing in z. spectral_distance(N1, N2, h): SN5 on two integer
    # count vectors of equal totals, result in the type of h; 0 for a zero total. spectrum_density(counts, h): S_c = N_c / (k h) for c = 1..B (the overflow counts enter k only).
    # stretching_range(sys): the symmetric derived range (-hi, hi) of SN12 padded outward by 4 to 16 eps relative (so that the rounding of a stretching number that reaches the bound, as for
    # Y5, cannot leave the range), for StandardMap, WrappedStandardMap and LinearSymplecticMap; ArgumentError for a map without stretching (singular values at most 1).
    for T in SN_TYPES
        epsT = sn_eps(T)
        zs = [T(-1.5) + T(i) * T(0.01) for i in 0:300]
        cells = [spectrum_cell(z, T(-1), T(1), 7) for z in zs]
        @test issorted(cells) && minimum(cells) == 0 && maximum(cells) == 8 && all(c -> 0 <= c <= 8, cells)
        @test spectrum_cell(T(0), T(-1), T(1), 1) == 1 && spectrum_cell(T(-2), T(-1), T(1), 1) == 0 && spectrum_cell(T(2), T(-1), T(1), 1) == 2
        h = T(1) / 2
        D = spectral_distance([0, 3, 0], [0, 0, 3], h)
        @test D isa T && abs(D - 2) <= 4 * epsT * 2 && spectral_distance([0, 3, 0], [0, 3, 0], h) == 0 && spectral_distance([0, 0, 0], [0, 0, 0], h) == 0
        S = spectrum_density([1, 2, 3, 4], h)
        @test length(S) == 2 && all(s -> s isa T, S) && abs(S[1] - T(2) / 5) <= 2 * epsT * S[1] && abs(S[2] - T(3) / 5) <= 2 * epsT * S[2]
        # stretching_range: the printed bounds of F6 and the closed forms (BigFloat)
        for (label, Kt) in ((:zero, exact(T, 0)), (:half, exact(T, 1 // 2)), (:two_pi, exact(T, pi) * 2), (:five, exact(T, 5)), (:ten, exact(T, 10)))
            L = sn_Lsigma(Kt)
            for sys in (StandardMap(Kt), WrappedStandardMap(Kt))
                lo, hi = stretching_range(sys)
                @test lo isa T && hi isa T && lo == -hi
                @test sn_big(hi) >= L * (1 + 4 * BigFloat(epsT)) && sn_big(hi) <= L * (1 + 16 * BigFloat(epsT))
            end
        end
        y5 = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        lo, hi = stretching_range(y5)
        Lmu = log(sn_mu())
        @test lo == -hi && sn_big(hi) >= Lmu * (1 + 4 * BigFloat(epsT)) && sn_big(hi) <= Lmu * (1 + 16 * BigFloat(epsT))
        @test sn_error(() -> stretching_range(LinearSymplecticMap(rotation_block(T, 0.3)))) isa ArgumentError
    end
end

@testset "API workspace surface of the stream" begin
    _ = (StretchingNumbers, reducer_workspace, reducer_bytes)
    # The additive Tangent surface. (1) A reducer without a workspace is untouched: reducer_workspace is 0, the six-argument reducer_step and the three-argument reducer_final default to the
    # existing forms (MEGNO). (2) propagate_orbit(...; workspace = ws): a workspace of the wrong length (shorter or longer than the chain's total) is an ArgumentError naming the workspace before any
    # work; without the keyword the body allocates its own and gives the same orbit and states; a stale workspace is zeroed at the start (the n0 = 0 rule of the kernel). (3) At most one reducer of a
    # chain may request a workspace (ArgumentError before any work). (4) reducer_final returns a copy of the counts. (5) Prefix invariance: the final at the handoff j of a long run is bitwise the
    # final of a separate run of length j tau. K = 2 pi, tau = 3, N = 60, both types.
    seed = sn_seed(15, 1)
    for T in SN_TYPES
        m = MEGNO()
        st = reducer_init(m, T, Val(2), Val(2))
        W = SMatrix{2, 2, T}(3, 4, 5, 12)
        x = SVector{2, T}(0, 0)
        a, Wa = reducer_step(m, st, 5, x, W)
        b, Wb = reducer_step(m, st, 5, x, W, UInt32[])
        @test isequal(a, b) && isequal(Wa, Wb) && isequal(reducer_final(m, a), reducer_final(m, a, UInt32[]))
        c = sn_case(:two_pi, T, 3)
        r = sn_red(c.sys, 100, 3)
        L = reducer_workspace(r, T, Val(2), Val(2))
        @test L == 2 * 102
        x1, W1, s1 = propagate_orbit(c.sys, c.x0[1], c.W0[1], 60, 3, (r,); orbit = 1)          # allocates its own workspace
        ws = zeros(UInt32, L)
        x2, W2, s2 = propagate_orbit(c.sys, c.x0[1], c.W0[1], 60, 3, (r,); orbit = 1, workspace = ws)
        @test x1 === x2 && W1 === W2 && isequal(s1, s2)
        f = reducer_final(r, s2[1], ws)
        @test sum(sn_cnt(f, 1)) == 20 && sum(sn_cnt(f, 2)) == 20
        wsx = fill(0xffffffff, L)                                                              # stale content
        _, _, s3 = propagate_orbit(c.sys, c.x0[1], c.W0[1], 60, 3, (r,); orbit = 1, workspace = wsx)
        @test wsx == ws && isequal(s3, s2)
        for bad in (zeros(UInt32, L - 1), zeros(UInt32, L + 1), UInt32[])
            e = sn_error(() -> propagate_orbit(c.sys, c.x0[1], c.W0[1], 60, 3, (r,); orbit = 1, workspace = bad))
            @test e isa ArgumentError && occursin("workspace", e.msg)
        end
        e = sn_error(() -> propagate_orbit(c.sys, c.x0[1], c.W0[1], 60, 3, (sn_red(c.sys, 100, 3; replace = false), r); orbit = 1))
        @test e isa ArgumentError && occursin("workspace", e.msg)
        # the counts are a copy
        f1 = reducer_final(r, s2[1], ws)
        f1.counts[2, 1] = 0x0000ffff
        f2 = reducer_final(r, s2[1], ws)
        @test f2.counts[2, 1] != 0x0000ffff && wsx == ws
        # prefix invariance
        long = sn_run((r,), c.sys, [c.x0[1]], [c.W0[1]], 60, 3; steps = collect(3:3:60))[1][1].snaps
        for k in 1:20
            wsk = zeros(UInt32, L)
            _, _, st_k = propagate_orbit(c.sys, c.x0[1], c.W0[1], 3k, 3, (r,); orbit = 1, workspace = wsk)
            @test isequal(reducer_final(r, st_k[1], wsk), long[k])
        end
    end
end

@testset "API final fields against the sums on real orbits" begin
    _ = (StretchingNumbers, stretching_range)
    # On real orbits (K = 2 pi and the regular class, both types, 5 starts, tau in {1, 4}, N = 400, every handoff): n = k tau, handoffs = k, the counts of the cells 0..B + 1 sum to k and no count
    # decreases; the mean equals (sum + comp) / n and the standard deviation is the root of (sum2 + comp2) / (tau n) - mean^2 (clamped at 0), both from the module's own totals in BigFloat within
    # 2 eps relative (the std within 4 eps of the radicand's scale); D is the pure function of the counts (bitwise).
    for T in SN_TYPES, label in (:two_pi, :half), tau in (1, 4)
        epsT = sn_eps(T)
        c = sn_case(label, T, 5)
        lo, hi = stretching_range(c.sys)
        B = 100
        h = (hi - lo) / B
        res = sn_run((sn_red(c.sys, B, tau),), c.sys, c.x0, c.W0, 400, tau; steps = collect(tau:tau:400))
        nbad = 0
        for j in 1:5
            prev = [zeros(Int, B + 2), zeros(Int, B + 2)]
            for (k, f) in enumerate(res[j][1].snaps)
                n = k * tau
                ok = f.n == n && f.handoffs == k
                for i in 1:2
                    cnt = sn_cnt(f, i)
                    ok &= sum(cnt) == k && all(cnt .>= prev[i])
                    prev[i] = cnt
                    setprecision(BigFloat, SN_BITS) do
                        A = sn_big(f.sum[i]) + sn_big(f.comp[i])
                        A2 = sn_big(f.sum2[i]) + sn_big(f.comp2[i])
                        ok &= abs(sn_big(f.mean[i]) - A / n) <= 2 * epsT * abs(A / n)
                        rad = A2 / (tau * n) - (A / n)^2
                        ok &= abs(sn_big(f.std[i])^2 - max(rad, zero(BigFloat))) <= 8 * epsT * (A2 / (tau * n) + (A / n)^2) + 4 * epsT * abs(sn_big(f.std[i])^2)
                    end
                end
                ok &= f.D === spectral_distance(sn_cnt(f, 1), sn_cnt(f, 2), h)
                nbad += !ok
            end
        end
        @test nbad == 0
    end
end
