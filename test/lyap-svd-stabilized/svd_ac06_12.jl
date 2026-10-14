# AC6 to AC9 and AC12 of lyap-svd-stabilized (AC10, AC11: stage 04). Seeds sv_seed(ac).

@testset "AC6 degenerate and regular cases" begin
    # V3: rotation by 0.7311 (exponents 0, 0, orthogonal tangent), n = 1e4, tau in {1, 10}: |lambda_i| <= theta_sum,
    # no health error, nothing non-finite (D is a sum of squares, safe at sigma_1 = sigma_2). K = 0 is in AC2.
    seed = sv_seed(6)
    for T in gs_types(), tau in (1, 10)
        A = rotation_block(T, 0.7311)
        sys = LinearSymplecticMap(A)
        res = sv_run(sys, [gs_zero(T, 2) for _ in 1:4], 10^4; tau, seed)
        th = gs_theta_lin(T, A, tau)
        for j in 1:4
            f = sv_fin(res, j)
            @test all(isfinite, f.lambda) && isfinite(f.c) && all(isfinite, res.W[j])
            @test all(i -> Float64(abs(f.lambda[i])) <= th, 1:2)
        end
    end
end

@testset "AC7 pairing and sum" begin
    # V10: |lambda_1 + lambda_2| <= theta_sum and equal to (S_1 + S_2)/n within 8 eps; eps_pair = eps_sum (diagnostics I3)
    seed = sv_seed(7)
    check(f, T, th) = begin
        s = f.lambda[1] + f.lambda[2]
        S = f.sum .+ f.comp
        Float64(abs(s)) <= th &&
            gs_absdiff(s, gs_big(S[1] + S[2]) / f.n) <= 8 * gs_eps(T) &&
            eps_pair(f.lambda) == eps_sum(f.lambda) && f.lambda[1] >= f.lambda[2]
    end
    for T in gs_types()
        sys = gs_cat(T)
        res = sv_run(sys, [gs_zero(T, 2) for _ in 1:4], 10^5; tau = 1, seed)
        @test all(j -> check(sv_fin(res, j), T, gs_theta(T, sv_fin(res, j).G_max, 1, 2.0)), 1:4)
    end
    for K in gs_KS()
        sys = gs_stdmap(Float64, K)
        x0s = gs_starts(Float64, seed + 1, 10)
        pts = [gs_points(sys, x, 10^4) for x in x0s]
        res = sv_run(sys, x0s, 10^4; tau = 1, seed)
        @test all(j -> check(sv_fin(res, j), Float64, gs_theta(Float64, sv_fin(res, j).G_max, 1, gs_Mmax(sys, pts[j]))), 1:10)
    end
end

@testset "AC8 wrapped state" begin
    # State length 4, d = 2: S_i and c bitwise equal to a plain per-orbit loop of the same map and tangent
    seed = sv_seed(8)
    for K in (2.0, 5.0), tau in (1, 3)
        sys = WrappedStandardMap(K)
        rng = StableRNG(seed + 1)
        x0s = [SVector{4,Float64}(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng), 0.0, 0.0) for _ in 1:20]
        res = propagate(sys, x0s, 1200; k = 2, tau, seed, reducers = (SVDSpectrum(),))   # the stream path
        r = SVDSpectrum()
        for j in 1:20
            W0 = initial_deviations(seed, j, Val(2), Val(2), Float64)
            xm, Wm, st = gs_manual_loop(r, sys, x0s[j], W0, 1200, tau)
            @test isequal(sv_fin(res, j), reducer_final(r, st))
            @test isequal(res.x[j], xm) && isequal(res.W[j], Wm) && length(res.x[j]) == 4
        end
    end
    # 10^6 steps, K = 5, 8 orbits, tau = 1: finite, |lambda_1 + lambda_2| <= theta_sum (M_max <= 1 + K = 6), no health
    # error, 0 bytes allocated after warm-up
    sys = WrappedStandardMap(5.0)
    rng = StableRNG(seed + 2)
    x0s = [SVector{4,Float64}(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng), 0.0, 0.0) for _ in 1:8]
    res = propagate(sys, x0s, 10^6; k = 2, tau = 1, seed, reducers = (SVDSpectrum(),))
    for j in 1:8
        f = sv_fin(res, j)
        @test all(isfinite, f.lambda) && isfinite(f.c)
        @test Float64(abs(f.lambda[1] + f.lambda[2])) <= gs_theta(Float64, f.G_max, 1, 6.0)
    end
    W0 = initial_deviations(seed, 1, Val(2), Val(2), Float64)
    @test gs_alloc_bytes(sys, x0s[1], W0, 10^4, 1, (SVDSpectrum(),)) == 0
end

@testset "AC9 thread and chunk invariance" begin
    # propagate (chunk 64, 1000, default) equals the serial map of propagate_orbit, bitwise: S_i, compensation, c, Q
    seed = sv_seed(9)
    n_orb, N = 4096, 1000
    sys = gs_stdmap(Float64, :twopi)
    x0s = gs_starts(Float64, seed + 1, n_orb)
    W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n_orb]
    reducers = (SVDSpectrum(),)
    serial = [propagate_orbit(sys, x0s[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
    xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
    fs = [map(reducer_final, reducers, s[3]) for s in serial]
    for chunk in (64, 1000, nothing)
        res = propagate(sys, x0s, N; k = 2, tau = 1, seed, reducers, chunk, W0)
        @test isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs)
    end
    # zero allocation of the per-orbit body after warm-up (isbits state)
    W = initial_deviations(seed, 1, Val(2), Val(2), Float64)
    @test gs_alloc_bytes(sys, x0s[1], W, 10^4, 1, (SVDSpectrum(),)) == 0
end

@testset "AC12 scope guard and provenance" begin
    # d != 2 or k != d is rejected before any work (the 4-dimensional test map throws if it is ever stepped);
    # exact message, exact record
    seed = sv_seed(12)
    msg(e) = e isa ArgumentError ? e.msg : ""
    for (sys, d, k) in ((sv_Boom4(), 4, 4), (sv_Boom4(), 4, 2), (sv_Boom4(), 4, 1), (gs_cat(Float64), 2, 1))
        x0s = [gs_zero(Float64, d)]
        err = try
            propagate(sys, x0s, 4; k, tau = 2, seed, reducers = (SVDSpectrum(),))
            nothing
        catch e
            e isa ArgumentError ? e : rethrow()
        end
        @test err isa ArgumentError && msg(err) == SV_GUARD_MESSAGE
        err = try
            reducer_init(SVDSpectrum(), Float64, Val(d), Val(k))
            nothing
        catch e
            e isa ArgumentError ? e : rethrow()
        end
        @test err isa ArgumentError && msg(err) == SV_GUARD_MESSAGE
    end
    @test reducer_init(SVDSpectrum(), Float64, Val(2), Val(2)) !== nothing
    # provenance
    @test SVDSpectrum().orthogonality_tol === 8.0
    for tol in (8.0, 4.0, Inf)
        @test reducer_provenance(SVDSpectrum(; orthogonality_tol = tol)) ==
              (estimator = "svd", orthonormalizer = "RGS", d = 2, orthogonality_tol = tol)
    end
    function check(T, tol, tau, extra)
        r = SVDSpectrum(; orthogonality_tol = tol)
        sys = gs_stdmap(T, 2.0)
        res = propagate(sys, gs_starts(T, seed + 1, 3), 12; k = 2, tau, seed, reducers = (r, extra...))
        p = res.provenance
        @test p.estimator == "svd" && p.orthonormalizer == "RGS" && p.d == 2 && p.orthogonality_tol === tol
        @test p.k == 2 && p.tau == tau && p.seed == seed && p.precision == precision_label(T)
        @test p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == ""
        @test !haskey(p, :iteration_path) && !haskey(p, :implementation)
    end
    for (tol, tau) in ((8.0, 1), (4.0, 3), (Inf, 3))
        check(Float64, tol, tau, ())
    end
    check(Float64, 8.0, 1, (ResidualMonitor(gs_stdmap(Float64, 2.0), 4.0, Inf),))   # a chained monitor adds no key
end
