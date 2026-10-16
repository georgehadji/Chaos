# AC6 to AC9 of lyap-gram-schmidt. Seeds gs_seed(ac).

@testset "AC6 orthogonality" begin
    seed = gs_seed(6)
    # (a) RGS: orth_max <= 8 k eps; (b) CGS and MGS with k = 2: orth_max <= 8 k G_w eps (G_w <= G + 1 at d = k = 2)
    for K in (0.5, 2.0, 5.0, 10.0, :twopi), tau in (1, 2, 4), alg in gs_algs()
        sys = gs_stdmap(Float64, K)
        res = gs_run(sys, gs_starts(Float64, seed + 1, 20), 1000; k = 2, tau, seed, alg)
        @test all(j -> Float64(gs_fin(res, j).orth_max) <= gs_orth_bound(alg, Float64, 2, gs_fin(res, j).G_max), 1:20)
    end
    # (c) discrimination on A_S, k = 4, 200 windows, Float64, 20 seeded frames (one orbit each)
    sys = LinearSymplecticMap(gs_AS(Float64))
    x0 = [gs_zero(Float64, 4)]
    for tau in (3, 1), j in 1:20
        W0 = [initial_deviations(seed + 2, j, Val(4), Val(4), Float64)]
        loss = Dict(gs_name(alg) => Float64(gs_fin(gs_run(sys, x0, 200 * tau; k = 4, tau, seed, alg, tol = Inf, W0), 1).orth_max)
                    for alg in gs_algs())
        if tau == 3
            @test loss["CGS"] / loss["MGS"] >= 100
            @test loss["MGS"] / loss["RGS"] >= 100
            # with the key value 8: CGS raises the orthogonality error, MGS and RGS do not
            err = gs_health_error(() -> gs_run(sys, x0, 600; k = 4, tau, seed, alg = CGS(), W0))
            @test err isa NumericalHealthError && err.what == "orthogonality" && err.step % tau == 0 && err.column == 0
            @test gs_health_error(() -> gs_run(sys, x0, 600; k = 4, tau, seed, alg = MGS(), W0)) === nothing
            @test gs_health_error(() -> gs_run(sys, x0, 600; k = 4, tau, seed, alg = RGS(), W0)) === nothing
        else
            @test loss["CGS"] > loss["MGS"] > loss["RGS"]
        end
    end
    # a non-finite orthogonality_tol still raises on a non-finite loss: see AC7 (NaN in W)
end

@testset "AC7 breakdown" begin
    # Integer one-step map A^p (exact): Float64 A^20 and Double64 A^40 (cond(W) eps_eff >= 1, derived 5e16 and mu^80) raise
    # "orthonormalisation breakdown" at handoff 1 for every variant and frame; A^16 and A^32 do not (eps G_w = 5e-3, 4e-4).
    # If A^40 does not break down in Double64 at GREEN, the spec says the power 44 is used (plan, "Conditional").
    seed = gs_seed(7)
    for (T, pbad, pok) in ((Float64, 20, 16), (Double64, 40, 32)), alg in gs_algsT(T)
        r = QRSpectrum(alg)
        sysb, syso = LinearSymplecticMap(gs_Apow(T, pbad)), LinearSymplecticMap(gs_Apow(T, pok))
        for j in 1:200
            W0 = initial_deviations(seed, j, Val(2), Val(2), T)
            err = gs_health_error(() -> propagate_orbit(sysb, gs_zero(T, 2), W0, 1, 1, (r,); orbit = j))
            @test err isa NumericalHealthError && err.what == "orthonormalisation breakdown"
            @test err isa NumericalHealthError && err.orbit == j && err.step == 1 && err.backend === :cpu
        end
        res = gs_serial(syso, [gs_zero(T, 2) for _ in 1:20], 200; k = 2, tau = 1, seed, reducers = (r,))   # no breakdown, n = 200
        G = maximum(j -> Float64(gs_fin(res, j).G_max), 1:20)
        @test all(j -> Float64(gs_fin(res, j).orth_max) <= gs_orth_bound(alg, T, 2, G), 1:20)
    end
    let syso = LinearSymplecticMap(gs_Apow(Float64, 16)), x0s = [gs_zero(Float64, 2) for _ in 1:10]
        rc, rm = (gs_run(syso, x0s, 200; k = 2, tau = 1, seed, alg) for alg in (CGS(), MGS()))
        same = isequal(rc.finals, rm.finals) && isequal(rc.W, rm.W)       # I6 at k = 2: reported, not gated (D43)
        @info "lyap-gram-schmidt AC7, I6 reported (D43): CGS and MGS bitwise equal at A^16: $same"
    end
    # Injected: a window with a zero column raises with the right column index; through propagate and through reducer_update
    for T in (Float64,), alg in gs_algsT(T)
        r = QRSpectrum(alg)
        sys4 = LinearSymplecticMap(one(SMatrix{4,4,T}))
        x4 = gs_zero(T, 4)
        for (zero_col, k) in ((3, 4), (1, 4))
            cols = ntuple(c -> c == zero_col ? zero(SVector{4,T}) : SVector{4,T}(ntuple(i -> i == c ? one(T) : zero(T), Val(4))), Val(k))
            W0 = SMatrix{4,k,T}(reduce(vcat, cols)...)
            err = gs_health_error(() -> gs_serial(sys4, [x4], 2; k, tau = 2, seed, reducers = (r,), W0 = [W0]))
            @test err isa NumericalHealthError && err.what == "orthonormalisation breakdown"
            @test err isa NumericalHealthError && err.column == zero_col && err.step == 2 && err.orbit == 1 && err.backend === :cpu
            st = reducer_init(r, T, Val(4), Val(k))
            err = gs_health_error(() -> reducer_update(r, st, 5, x4, W0))
            @test err isa NumericalHealthError && err.what == "orthonormalisation breakdown"
            @test err isa NumericalHealthError && err.column == zero_col && err.step == 5 && err.orbit == 0
        end
        # NaN in W: tangent-core's finiteness check fires first inside propagate ("deviation matrix"); the estimator
        # alone (reducer_update) raises the breakdown with an infinite threshold and the column of the NaN.
        Wn = SMatrix{4,4,T}(1, 0, 0, 0, 0, 1, 0, 0, NaN, 0, 0, 0, 0, 0, 0, 1)
        err = gs_health_error(() -> gs_serial(sys4, [x4], 2; k = 4, tau = 2, seed, reducers = (r,), W0 = [Wn]))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 2
        st = reducer_init(r, T, Val(4), Val(4))
        err = gs_health_error(() -> reducer_update(r, st, 2, x4, Wn))
        @test err isa NumericalHealthError && err.what == "orthonormalisation breakdown" && err.column == 3 && isinf(err.threshold)
    end
end

# Independent S_i loop: orthonormalize + compensated add, the contract of the step (api-batch-b.md, QRSpectrum step 4).
function gs_s_loop(alg, sys, x0::SVector{s,T}, W0::SMatrix{d,k,T}, N, tau) where {s,d,k,T}
    x, W = x0, W0
    acc = ntuple(_ -> CompensatedSum{T}(), Val(k))
    for _ in 1:(N ÷ tau)
        for _ in 1:tau
            W = tangent(sys, x) * W
            x = forward(sys, x)
        end
        Q, R = orthonormalize(alg, W)
        acc = ntuple(i -> add(acc[i], log(R[i, i])), Val(k))
        W = Q
    end
    return acc
end

@testset "AC8 wrapped state" begin
    # State length 4, d = 2: S_i bitwise equal to a plain per-orbit loop (reducer driven by hand, and orthonormalize + add)
    seed = gs_seed(8)
    for K in (2.0, 5.0), tau in (1, 3), alg in gs_algs()
        sys = WrappedStandardMap(K)
        rng = StableRNG(seed + 1)
        x0s = [SVector{4,Float64}(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng), 0.0, 0.0) for _ in 1:20]
        res = propagate(sys, x0s, 1200; k = 2, tau, seed, reducers = (QRSpectrum(alg),))   # the stream path
        r = QRSpectrum(alg)
        for j in 1:20
            W0 = initial_deviations(seed, j, Val(2), Val(2), Float64)
            xm, Wm, st = gs_manual_loop(r, sys, x0s[j], W0, 1200, tau)
            @test isequal(gs_fin(res, j), reducer_final(r, st))
            @test Float64(gs_fin(res, j).orth_max) <= gs_orth_bound(alg, Float64, 2, gs_fin(res, j).G_max)   # AC6(a), (b)
            @test isequal(res.x[j], xm) && isequal(res.W[j], Wm) && length(res.x[j]) == 4
            acc = gs_s_loop(alg, sys, x0s[j], W0, 1200, tau)
            @test all(i -> isequal(gs_fin(res, j).sum[i], acc[i].s) && isequal(gs_fin(res, j).comp[i], acc[i].c), 1:2)
        end
    end
    # 10^6 steps, K = 5, 8 orbits, tau = 1, RGS: finite, |lambda_1 + lambda_2| <= theta_sum, no health error, no allocation.
    # M_max <= 1 + K for this map (largest tangent entry 1 + |K cos x|), so theta_sum uses 1 + K.
    sys = WrappedStandardMap(5.0)
    rng = StableRNG(seed + 2)
    x0s = [SVector{4,Float64}(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng), 0.0, 0.0) for _ in 1:8]
    res = propagate(sys, x0s, 10^6; k = 2, tau = 1, seed, reducers = (QRSpectrum(RGS()),))
    for j in 1:8
        f = gs_fin(res, j)
        @test all(isfinite, f.lambda) && all(isfinite, res.W[j])
        @test Float64(abs(f.lambda[1] + f.lambda[2])) <= gs_theta(Float64, f.G_max, 1, 6.0)
    end
    W0 = initial_deviations(seed, 1, Val(2), Val(2), Float64)
    @test gs_alloc_bytes(sys, x0s[1], W0, 10^4, 1, (QRSpectrum(RGS()),)) == 0
end

@testset "AC9 thread and chunk invariance, zero allocation" begin
    # propagate (chunk 64, 1000, default) equals the serial map of propagate_orbit, bitwise; chain (estimator, monitor)
    seed = gs_seed(9)
    N = 1000
    sys = gs_stdmap(Float64, :twopi)
    for alg in gs_algs()
        n_orb = alg isa RGS ? 4096 : 512       # spec size for RGS; CGS and MGS reduced for the 3-minute budget (plan)
        x0s = gs_starts(Float64, seed + 1, n_orb)
        W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n_orb]
        reducers = (QRSpectrum(alg),)
        serial = [propagate_orbit(sys, x0s[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (alg isa RGS ? (64, 1000, nothing) : (64, nothing))
            res = propagate(sys, x0s, N; k = 2, tau = 1, seed, reducers, chunk, W0)
            @test isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs)
        end
    end
    # zero allocation after warm-up: n = 10^4 steps, k = d = 4 (Y5 pair), one orbit, all variants
    sys4 = LinearSymplecticMap(gs_AS(Float64))
    W4 = initial_deviations(seed, 1, Val(4), Val(4), Float64)
    for alg in gs_algs()
        @test gs_alloc_bytes(sys4, gs_zero(Float64, 4), W4, 10^4, 1, (QRSpectrum(alg),)) == 0
    end
end
