# API tests: the reducer protocol (init, the n = 0 final, the arithmetic of three handoffs, failures), consistency of the final fields with the sums on real orbits,
# prefix invariance. AC14 and AC15 are stage 04 (see the test plan).

@testset "API reducer protocol" begin
    _ = (MEGNO,)
    # init is total except k != 2 (ArgumentError "MEGNO needs k = 2, got k = $k"); the n = 0 final (all zero, handoffs 0); the arithmetic of three handoffs at the map steps
    # 7, 14, 21 (tau = 7) on the columns with norms (5, 5), (10, 10), (5/8, 1/2): the log increments a_k = ln alpha_k, the sums s, R1 = sum t_k a_k, R2 = sum t_k^2 a_k,
    # Y = 2 R1/n, Q = sum of the Y, Ybar = Q/K, Q2 = sum of 3 R2, Ybar20 = tau Q2/n^3 (M5); at the first handoff Y = Ybar = 2 a_1 and Ybar20 = 3 a_1; the replacer returns the
    # column-normalised matrix, the non-replacing position returns W unchanged; a failure (zero or non-finite norm) is recorded without throwing in reducer_step, sums, series
    # and W unchanged, thrown by reducer_update (orbit 0, backend :cpu).
    @test_throws ArgumentError reducer_init(MEGNO(), Float64, Val(2), Val(1))
    @test_throws ArgumentError reducer_init(MEGNO(replace = false), Float64, Val(2), Val(3))
    e1 = try
        reducer_init(MEGNO(), Float64, Val(2), Val(1))
    catch ex
        ex
    end
    @test e1 isa ArgumentError && occursin("MEGNO needs k = 2, got k = 1", e1.msg)
    for T in MG_TYPES, replace in (true, false)
        epsT = mg_eps(T)
        lg(x) = setprecision(() -> log(BigFloat(x)), BigFloat, MG_BITS)
        r = MEGNO(replace = replace)
        st0 = reducer_init(r, T, Val(2), Val(2))
        f0 = reducer_final(r, st0)
        @test f0.n == 0 && f0.handoffs == 0
        @test all(iszero, f0.sum) && all(iszero, f0.comp) && all(iszero, f0.Y) && all(iszero, f0.Ybar) && all(iszero, f0.Ybar20)
        @test all(iszero, f0.R1) && all(iszero, f0.R2) && all(iszero, f0.Q) && all(iszero, f0.Q2)
        for fld in (:Y, :Ybar, :Ybar20, :sum, :comp, :R1, :R2, :Q, :Q2)
            @test getproperty(f0, fld) isa SVector{2, T}
        end
        @test !reducer_failed(r, st0)
        x = SVector{2, T}(0, 0)
        Ws = (SMatrix{2, 2, T}(3, 4, 0, 5), SMatrix{2, 2, T}(6, 8, 0, 10), SMatrix{2, 2, T}(0.375, 0.5, 0, 0.5))
        alphas = ((5, 5), (10, 10), (BigFloat(5) / 8, BigFloat(1) / 2))
        st1, W1 = reducer_step(r, st0, 7, x, Ws[1])
        if replace
            @test all(abs.(W1 - SMatrix{2, 2, T}(T(3) / 5, T(4) / 5, 0, 1)) .<= 2 * epsT)
        else
            @test W1 == Ws[1]
        end
        stu, Wu = reducer_update(r, st0, 7, x, Ws[1])
        @test isequal(stu, st1) && isequal(Wu, W1)
        f1 = reducer_final(r, st1)
        @test f1.n == 7 && f1.handoffs == 1
        for i in 1:2
            @test mg_err(f1.Y[i], 2 * lg(5)) <= 4 * epsT && mg_err(f1.Ybar[i], 2 * lg(5)) <= 4 * epsT && mg_err(f1.Ybar20[i], 3 * lg(5)) <= 6 * epsT
        end
        st2, _ = reducer_step(r, st1, 14, x, Ws[2])
        st3, _ = reducer_step(r, st2, 21, x, Ws[3])
        for (st, K) in ((st2, 2), (st3, 3))
            f = reducer_final(r, st)
            @test f.n == 7 * K && f.handoffs == K
            for i in 1:2
                a = setprecision(BigFloat, MG_BITS) do
                    [log(BigFloat(alphas[k][i])) for k in 1:K]
                end
                sh = setprecision(() -> cumsum(a), BigFloat, MG_BITS)
                ref = mg_M5(sh, 7)
                tol = 16 * epsT * max(1, maximum(abs, sh) * K)
                @test mg_err(f.sum[i] + f.comp[i], sh[K]) <= 4 * epsT * max(1, maximum(abs, sh))
                @test mg_err(f.Y[i], ref.Y[K]) <= tol && mg_err(f.Ybar[i], ref.Ybar[K]) <= tol && mg_err(f.Ybar20[i], ref.Ybar20[K]) <= tol
                R1, R2, Q, Q2 = setprecision(BigFloat, MG_BITS) do
                    Yk = [2 * sum(7 * k * a[k] for k in 1:kk) / (7 * kk) for kk in 1:K]
                    R2k = [sum((7 * k)^2 * a[k] for k in 1:kk) for kk in 1:K]
                    (sum(7 * k * a[k] for k in 1:K), R2k[K], sum(Yk), sum(3 * R2k))
                end
                @test mg_err(f.R1[i], R1) <= tol && mg_err(f.R2[i], R2) <= tol * 49 * K && mg_err(f.Q[i], Q) <= tol * K && mg_err(f.Q2[i], Q2) <= tol * 49 * K * K
            end
        end
        # failure: column 2 zero, column 1 zero, a non-finite column
        Wz = SMatrix{2, 2, T}(3, 4, 0, 0)
        stf, Wf = reducer_step(r, st2, 28, x, Wz)
        f2 = reducer_final(r, st2)
        @test reducer_failed(r, stf) && !reducer_failed(r, st2)
        @test Wf === Wz
        ff = reducer_final(r, stf)
        @test isequal(ff, f2)                                                                # nothing stored, sums and series unchanged
        e = reducer_error(r, stf, 9, :cpu)
        @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 9 && e.step == 28
        @test e.backend === :cpu && e.column == 2 && e.value == 0.0
        e2 = mg_health_error(() -> reducer_update(r, st2, 28, x, Wz))
        @test e2 isa NumericalHealthError && e2.what == "deviation matrix" && e2.orbit == 0 && e2.backend === :cpu && e2.step == 28 && e2.column == 2
        st1z, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(0, 0, 0, 5))
        @test reducer_failed(r, st1z) && reducer_error(r, st1z, 1, :cpu).column == 1
        stn, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(3, 4, T(NaN), 5))
        @test reducer_failed(r, stn)
        sti, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(3, 4, T(Inf), 5))
        @test reducer_failed(r, sti)
        # a failed state stays failed and unchanged
        stf2, _ = reducer_step(r, stf, 35, x, Ws[3])
        @test isequal(reducer_final(r, stf2), f2) && reducer_failed(r, stf2)
    end
    @test reducer_provenance(MEGNO()).indicator == "megno"
end

@testset "API final fields against the sums on real orbits" begin
    _ = (MEGNO,)
    # Y, Ybar and Ybar20 of the final are functions of the module's own sums s_i (M6): K in {0.5 (regular class), 2 pi}, N = 400, tau in {1, 4}, every handoff, both columns,
    # Float64 (10 starts) and Double64 (4 starts): the three series equal the BigFloat Abel form of the module's own totals s_i = sum + comp within 16 eps max(1, max |s|)
    # (a few ulp of the largest partial sum); the field `n` is the map step, `handoffs` = n/tau; the weighted sums are consistent with Y = 2 R1/n, Ybar = Q/K,
    # Ybar20 = tau Q2/n^3 to 8 ulp.
    for T in MG_TYPES, label in (:half, :two_pi), tau in (1, 4)
        epsT = mg_eps(T)
        m = T === Float64 ? 10 : 4
        c = mg_case(label, T, m)
        run = mg_run((MEGNO(),), c.sys, c.x0, c.W0, 400, tau; all = true)
        for r in run
            sn = r[1].snaps
            @test length(sn) == 400 ÷ tau
            for i in 1:2
                S = mg_S(sn, i)
                post = mg_M6(setprecision(() -> [mg_big(s) for s in S], BigFloat, MG_BITS), tau)
                tol = 16 * epsT * max(1, maximum(abs, S))
                @test mg_nbad(k -> mg_err(sn[k].Y[i], post.Y[k]) <= tol && mg_err(sn[k].Ybar[i], post.Ybar[k]) <= tol && mg_err(sn[k].Ybar20[i], post.Ybar20[k]) <= tol, eachindex(sn)) == 0
            end
            for (k, s) in enumerate(sn)
                @test s.n == k * tau && s.handoffs == k
                for i in 1:2
                    @test abs(s.Y[i] - 2 * s.R1[i] / s.n) <= 8 * epsT * max(1, abs(s.Y[i]))
                    @test abs(s.Ybar[i] - s.Q[i] / k) <= 8 * epsT * max(1, abs(s.Ybar[i]))
                    @test abs(s.Ybar20[i] - tau * s.Q2[i] / T(s.n)^3) <= 8 * epsT * max(1, abs(s.Ybar20[i]))
                end
            end
        end
    end
end

@testset "API prefix invariance" begin
    _ = (MEGNO,)
    # The state after handoff j is bitwise identical for every N >= j tau (same seed, orbit, start, tau): a run of length j tau returns exactly the series value at j tau.
    # All-handoff snapshots of one run of 400 steps against separate runs, tau in {1, 4}, both types; the whole final is compared.
    seed = mg_seed(14, 20)
    for T in MG_TYPES, tau in (1, 4)
        sys = StandardMap(mg_K(T, :two_pi))
        x0 = mg_cast(T, SVector(1.0, 2.0))
        W0 = mg_frame(seed, 3, T)
        probe = mg_Probe(MEGNO(); all = true)
        _, _, st = propagate_orbit(sys, x0, W0, 400, tau, (probe,); orbit = 3)
        snaps = reducer_final(probe, st[1]).snaps
        @test length(snaps) == 400 ÷ tau
        for j in (1, 5, 25, 100 ÷ tau)
            N = j * tau
            _, _, s2 = propagate_orbit(sys, x0, W0, N, tau, (MEGNO(),); orbit = 3)
            f = reducer_final(MEGNO(), s2[1])
            @test snaps[j].n == N && f.n == N
            @test isequal(snaps[j], f)
        end
    end
end
