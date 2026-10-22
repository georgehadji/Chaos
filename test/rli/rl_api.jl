# API tests: the reducer protocol (init, the n = 0 final, the arithmetic of three handoffs, failures), memory of the mean (RL8), consistency of the final fields with the sums on
# real orbits, prefix invariance. AC14 and AC15 are stage 04 (see the test plan).

@testset "API reducer protocol" begin
    _ = (RLI,)
    # init is total except k != 2 (ArgumentError "RLI needs k = 2, got k = $k") and an odd dimension; the n = 0 final (all zero, handoffs 0); the arithmetic of three handoffs at the map steps
    # 7, 14, 21 (tau = 7) on the columns with norms (alpha_A, alpha_B) = (5, 13), (10, 10), (5/8, 1/2): a_k = ln alpha, s_A, s_B the compensated sums, DL_k = |s_B - s_A| / n, Q = sum DL,
    # RLI = Q / K (RL4); at the first handoff RLI = DL_1; the replacer returns the column-normalised matrix, the non-replacing position returns W unchanged; a failure (zero or non-finite
    # norm) is recorded without throwing in reducer_step, sums, series and W unchanged, thrown by reducer_update (orbit 0, backend :cpu), and a failed state stays failed and unchanged.
    @test_throws ArgumentError reducer_init(RLI(), Float64, Val(4), Val(1))
    @test_throws ArgumentError reducer_init(RLI(replace = false), Float64, Val(4), Val(3))
    @test_throws ArgumentError reducer_init(RLI(), Float64, Val(3), Val(2))
    e1 = try
        reducer_init(RLI(), Float64, Val(4), Val(1))
    catch ex
        ex
    end
    @test e1 isa ArgumentError && occursin("RLI needs k = 2, got k = 1", e1.msg)
    for T in RL_TYPES, replace in (true, false)
        epsT = rl_eps(T)
        r = RLI(replace = replace)
        st0 = reducer_init(r, T, Val(4), Val(2))
        f0 = reducer_final(r, st0)
        @test f0.n == 0 && f0.handoffs == 0 && f0.DL == 0 && f0.RLI == 0 && f0.Q == 0
        @test f0.sum isa SVector{2, T} && f0.comp isa SVector{2, T} && all(iszero, f0.sum) && all(iszero, f0.comp)
        x = zero(SVector{8, T})
        cols = (((3.0, 4.0), (5.0, 12.0)), ((6.0, 8.0), (6.0, 8.0)), ((0.375, 0.5), (0.5, 0.0)))
        W(ca, cb) = SMatrix{4, 2, T}(T(ca[1]), T(ca[2]), 0, 0, 0, 0, T(cb[1]), T(cb[2]))
        alpha = ((5.0, 13.0), (10.0, 10.0), (0.625, 0.5))
        st = st0
        sA = sB = zero(BigFloat)
        Q = zero(BigFloat)
        for (k, ((ca, cb), (aA, aB))) in enumerate(zip(cols, alpha))
            n = 7k
            Wk = W(ca, cb)
            st, W2 = reducer_step(r, st, n, x, Wk)
            f = reducer_final(r, st)
            setprecision(BigFloat, RL_BITS) do
                sA += log(BigFloat(aA)); sB += log(BigFloat(aB))
                DLref = abs(sB - sA) / n
                Q += DLref
                @test f.n == n && f.handoffs == k
                @test abs(rl_big(f.sum[1]) + rl_big(f.comp[1]) - sA) <= 4 * epsT * max(1, abs(sA)) && abs(rl_big(f.sum[2]) + rl_big(f.comp[2]) - sB) <= 4 * epsT * max(1, abs(sB))
                @test abs(rl_big(f.DL) - DLref) <= 8 * epsT * max(abs(sA), abs(sB), 1) / n + 4 * epsT * DLref
                @test abs(rl_big(f.RLI) - Q / k) <= 8 * epsT * max(abs(sA), abs(sB), 1) * rl_H(k) / n + 4 * epsT * Q / k
                @test abs(rl_big(f.Q) - Q) <= 8 * epsT * max(abs(sA), abs(sB), 1) * rl_H(k) + 4 * epsT * Q
                k == 1 && @test f.RLI == f.DL
            end
            if replace
                Wref = SMatrix{4, 2, T}(T(ca[1]) / T(aA), T(ca[2]) / T(aA), 0, 0, 0, 0, T(cb[1]) / T(aB), T(cb[2]) / T(aB))
                @test all(abs.(W2 - Wref) .<= 4 * epsT)
            else
                @test isequal(W2, Wk)
            end
        end
        # failures: a zero column 1, then a non-finite column 2, then both bad (column 1 reported)
        good = W((3, 4), (5, 12))
        stg, _ = reducer_step(r, st0, 7, x, good)
        for (Wbad, col, val, thr) in ((W((0, 0), (5, 12)), 1, 0.0, 0.0), (W((3, 4), (T(NaN), 1)), 2, NaN, Inf), (W((0, 0), (0, 0)), 1, 0.0, 0.0),
                                      (W((T(Inf), 1), (5, 12)), 1, Inf, Inf))
            stf, Wf = reducer_step(r, stg, 14, x, Wbad)
            @test reducer_failed(r, stf) && !reducer_failed(r, stg)
            @test isequal(Wf, Wbad)
            ff, fg = reducer_final(r, stf), reducer_final(r, stg)
            @test isequal(ff, fg)                                           # sums, series, n and handoffs unchanged by the failure
            err = reducer_error(r, stf, 9, :cpu)
            @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 9 && err.step == 14 && err.backend === :cpu && err.column == col
            @test (isnan(val) ? isnan(err.value) : err.value == val) && err.threshold == thr
            e = try
                reducer_update(r, stg, 14, x, Wbad)
                nothing
            catch ex
                ex
            end
            @test e isa NumericalHealthError && e.orbit == 0 && e.backend === :cpu && e.step == 14 && e.column == col
            stf2, Wf2 = reducer_step(r, stf, 21, x, good)                   # a failed state stays failed and unchanged
            @test reducer_failed(r, stf2) && isequal(stf2, stf) && isequal(Wf2, good)
        end
    end
end

@testset "API memory of the mean and the range" begin
    _ = (RLI,)
    # RL8 / I1: RLI is a running mean of nonnegative DL_k, not a function of DL at the final handoff. Two series with the same final DL_4: DL = (c, c, c, c) (the log ratio of the two
    # stretches grows by c tau per handoff) and DL = (0, 0, 0, c) (the ratio is 1 until the last handoff, then exp(4 c tau)): the final DL agree and RLI = c and c / 4. Also DL >= 0 and
    # RLI <= max DL on a real orbit.
    for T in RL_TYPES
        epsT = rl_eps(T)
        tau, c = 5, 0.03
        r = RLI()
        x = zero(SVector{8, T})
        run_series(ratios) = begin
            st = reducer_init(r, T, Val(4), Val(2))
            for (k, q) in enumerate(ratios)
                W = SMatrix{4, 2, T}(1, 0, 0, 0, 0, 0, T(q), 0)           # alpha_A = 1 (s_A stays 0), alpha_B = q
                st, _ = reducer_step(r, st, tau * k, x, W)
            end
            reducer_final(r, st)
        end
        q = exp(c * tau)
        fa = run_series([q, q, q, q])
        fb = run_series([1, 1, 1, exp(4 * c * tau)])
        @test abs(fa.DL - c) <= 1e-12 && abs(fb.DL - c) <= 1e-12 && abs(fa.DL - fb.DL) <= 1e-12
        @test abs(fa.RLI - c) <= 1e-12 && abs(fb.RLI - c / 4) <= 1e-12
        sys = rl_sys(rl_K(T, :two_pi), true)
        xi = rl_xi0(rl_seed(15, 3), 1, T)
        pr = rl_pair(sys, rl_cast(T, SVector(1.0, 2.0)), xi, 1e-9)
        snaps = rl_run((RLI(),), PairMap(sys), [pr.z], [rl_frame(xi)], 400, 1; all = true)[1][1].snaps
        dls = [s.DL for s in snaps]
        @test all(>=(0), dls) && all(s -> 0 <= s.RLI <= maximum(dls[1:s.n]) + 4 * epsT * maximum(dls), snaps) && maximum(dls) <= 2 * log(rl_sigma(Float64(rl_K(T, :two_pi))))
    end
end

@testset "API final fields against the sums on real orbits" begin
    _ = (RLI, PairMap)
    # On real pair orbits (K = 2 pi and the regular class, both types, 5 starts, tau in {1, 4}, N = 400, every handoff): n = k tau and handoffs = k; DL = |S_B - S_A| / n and the mean RLI
    # = Q / K from the module's own sums S_i = sum + comp, evaluated in BigFloat (the post-processing of a stored series): |DL - DL_post| <= 16 eps max(1, |S|) / n + 4 eps DL_post and
    # |RLI - RLI_post| <= 16 eps max(1, max |S|) H_K / n + 4 eps RLI_post; the internal identity RLI = Q / K within 2 eps.
    for T in RL_TYPES, label in (:two_pi, :half), tau in (1, 4)
        epsT = rl_eps(T)
        c = rl_case(label, T, 5)
        prs = rl_pairs(c, 1e-9)
        res = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), 400, tau; all = true)
        for j in 1:5
            sn = res[j][1].snaps
            @test length(sn) == 400 ÷ tau
            sAv, sBv = zeros(BigFloat, 400), zeros(BigFloat, 400)
            setprecision(BigFloat, RL_BITS) do
                for s in sn
                    sAv[s.n] = rl_big(s.sum[1]) + rl_big(s.comp[1])
                    sBv[s.n] = rl_big(s.sum[2]) + rl_big(s.comp[2])
                end
            end
            post = rl_ref_rli(sAv, sBv, tau)
            smax = maximum(max(abs(sAv[s.n]), abs(sBv[s.n])) for s in sn)
            nbad = 0
            for (k, s) in enumerate(sn)
                ok = s.n == k * tau && s.handoffs == k
                ok &= rl_err(s.DL, post.DL[k]) <= 16 * epsT * max(1.0, Float64(max(abs(sAv[s.n]), abs(sBv[s.n])))) / s.n + 4 * epsT * Float64(post.DL[k])
                ok &= rl_err(s.RLI, post.RLI[k]) <= 16 * epsT * max(1.0, Float64(smax)) * rl_H(k) / s.n + 4 * epsT * Float64(post.RLI[k])
                ok &= abs(s.RLI - s.Q / T(k)) <= 2 * epsT * s.RLI
                nbad += !ok
            end
            @test nbad == 0
        end
    end
end

@testset "API prefix invariance" begin
    _ = (RLI, PairMap)
    # The series read at every handoff of one run of N steps is bitwise the final of separate runs of length j tau (nothing depends on the run length), K = 2 pi, tau = 3, N = 60,
    # Float64 and Double64.
    for T in RL_TYPES
        c = rl_case(:two_pi, T, 3)
        prs = rl_pairs(c, 1e-9)
        long = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), 60, 3; all = true)
        for j in 1:3, k in 1:20
            _, _, st = propagate_orbit(PairMap(c.sys), prs[j].z, rl_frame(c.xi0[j]), 3k, 3, (RLI(),); orbit = j)
            @test isequal(reducer_final(RLI(), st[1]), long[j][1].snaps[k])
        end
    end
end
