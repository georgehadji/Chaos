# AC5 to AC8: allocation, thread/chunk invariance, seeded initialization, numerical health.

@testset "AC5 zero allocations" begin
    # propagate_orbit allocates 0 bytes after warm-up: k = 2, chain (column normalizer, reader), n = 1e4.
    sys = tc_Std(2pi)
    x0 = SVector(1.0, 0.5)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    @test tc_alloc_bytes(sys, x0, W0, 10_000, 1, (tc_LogNorm(), tc_Reader(1))) == 0
    @test tc_alloc_bytes(sys, x0, W0, 10_000, 4, (tc_LogNorm(), tc_Reader(4))) == 0
    @test tc_alloc_bytes(sys, x0, W0, 100, 1, ()) == 0   # empty chain; short horizon, W is not renormalized
end

@testset "AC6 thread and chunk invariance" begin
    # propagate (chunk 64, 1000, default) equals the serial map of propagate_orbit, bitwise: x, W and accumulated logs.
    seed = TC_SEED + 6
    n_orb, N = 1024, 1000   # D38: 4096 orbits (and chunk 1000) are in verify/tangent-core.jl
    sys = tc_Std(2pi)
    rng = StableRNG(seed)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:n_orb]
    W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n_orb]
    reducers = (tc_LogNorm(), tc_Reader(1))
    serial = [propagate_orbit(sys, x0[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
    xs = [s[1] for s in serial]
    Ws = [s[2] for s in serial]
    fs = [map(reducer_final, reducers, s[3]) for s in serial]
    for chunk in (64, nothing)   # D38: chunk 1000 moved to verify/tangent-core.jl (its AC6 digests cover it)
        res = propagate(sys, x0, N; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
        @test isequal(res.x, xs)
        @test isequal(res.W, Ws)
        @test isequal(res.finals, fs)
    end
    # Default W0 (seeded initial_deviations) gives the same result as passing W0 explicitly.
    res = propagate(sys, x0, N; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = 64)
    @test isequal(res.W, Ws)
end

@testset "AC7 seeded initialization" begin
    # Q^T Q = I to 8 k eps; (s, j) -> matrix is a pure function; uniform angle for d = 2, k = 1 (KS, 1% level 1.63/sqrt(n)).
    seed = TC_SEED + 7
    n = 10_000
    # D38: Float64 d = k = 8 and Double64 d = 2 k = 1, d = k = 4, d = k = 8 moved to verify/tangent-core.jl (also its AC7a rows).
    for (d, k) in ((2, 1), (2, 2), (4, 4))
        worst = 0.0
        for j in 1:n
            Q = initial_deviations(seed, j, Val(d), Val(k), Float64)
            worst = max(worst, maximum(abs, Q' * Q - I))
        end
        @test worst <= 8 * k * eps(Float64)
        @test initial_deviations(seed, 1, Val(d), Val(k), Float64) isa SMatrix{d,k,Float64}
    end
    for (d, k) in ((2, 2),)
        worst = 0.0
        for j in 1:200
            Q = initial_deviations(seed, j, Val(d), Val(k), Double64)
            worst = max(worst, Float64(maximum(abs, Q' * Q - I)))
        end
        @test worst <= 8 * k * Float64(eps(Double64))
    end
    # (b) repeated calls, reversed order and worker threads give identical bits; distinct orbits and seeds differ.
    ref = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n]
    @test isequal(ref, [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n])
    @test isequal(ref, reverse([initial_deviations(seed, j, Val(2), Val(2), Float64) for j in n:-1:1]))
    thr = Vector{SMatrix{2,2,Float64,4}}(undef, n)
    Threads.@threads for j in 1:n
        thr[j] = initial_deviations(seed, j, Val(2), Val(2), Float64)
    end
    @test isequal(ref, thr)
    @test length(unique(Tuple.(ref))) == n
    @test initial_deviations(seed + 1, 1, Val(2), Val(2), Float64) != ref[1]
    # (b) chunk positions: through propagate with the identity map the default W0 equals initial_deviations(seed, j).
    ident = tc_Lin(SMatrix{2,2,Float64}(1, 0, 0, 1))
    for chunk in (64, 1000)
        res = propagate(ident, [SVector(0.0, 0.0) for _ in 1:300], 4; k = 2, tau = 2, seed = seed, chunk = chunk, reducers = (tc_Reader(2),))   # D38: AC1's kernel specialization (a read-only reducer leaves W unchanged)
        @test all(j -> isequal(res.W[j], ref[j]), 1:300)
    end
    # (c) angle of the single deviation vector is uniform on [0, 2 pi)
    u = [mod2pi(atan(Q[2, 1], Q[1, 1])) / (2pi) for Q in (initial_deviations(seed, j, Val(2), Val(1), Float64) for j in 1:n)]
    @test tc_ks_uniform(u) <= 1.63 / sqrt(n)
end

@testset "AC8 numerical health" begin
    # Non-finite x or W is caught at the first handoff after it appears, with orbit, step, backend and which quantity.
    seed = TC_SEED + 8
    x0s(flagged) = [SVector(0.0, 0.0, j == flagged ? 1.0 : 0.0) for j in 1:6]
    # control: no flagged orbit, no error
    ok = propagate(tc_Fault(:inf_state, 7), x0s(0), 20; k = 1, tau = 5, seed = seed, chunk = 2)
    @test all(x -> x[1] == 20 && isfinite(x[2]), ok.x)
    # (a) Inf in the state at step 7, tau = 5 -> error at handoff 10, nothing stored
    res = nothing
    err = tc_health_error(() -> (res = propagate(tc_Fault(:inf_state, 7), x0s(3), 20; k = 1, tau = 5, seed = seed, chunk = 2)))
    @test err isa NumericalHealthError
    @test res === nothing
    @test err.orbit == 3 && err.step == 10 && err.what == "state" && err.backend === :cpu
    @test err.threshold == Inf && !isfinite(err.value)
    # (b) NaN in W injected at step 3, tau = 5 -> error at handoff 5
    res = nothing
    err = tc_health_error(() -> (res = propagate(tc_Fault(:nan_W, 3), x0s(5), 20; k = 1, tau = 5, seed = seed, chunk = 2)))
    @test err isa NumericalHealthError
    @test res === nothing
    @test err.orbit == 5 && err.step == 5 && err.what == "deviation matrix" && err.backend === :cpu
    @test err.threshold == Inf && !isfinite(err.value)
    # (c) overflow (A = diag(e^3, e^-3), a k = 1 kernel specialization of its own): moved to verify/tangent-core.jl (D38)
    z0 = [SVector(0.0, 0.0)]
    # (d) contract: N must be a positive multiple of tau; ArgumentError before any work (tc_Boom throws if the map is called)
    @test_throws ArgumentError propagate(tc_Boom(), z0, 7; k = 1, tau = 5, seed = seed)
    @test_throws ArgumentError propagate(tc_Boom(), z0, 0; k = 1, tau = 1, seed = seed)
    @test_throws ArgumentError propagate(tc_Boom(), z0, -5; k = 1, tau = 5, seed = seed)
end
