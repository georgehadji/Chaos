# Amendment D28, criterion W9: state longer than the tangent dimension (s = 4, d = 2). Spec: tangent-core-spec.md, "Amendment D28".

@testset "W9 state longer than the tangent dimension" begin
    # s = 4, d = 2: memory formula first, then per-orbit body vs manual loop, chunked vs serial, health, allocation, reducer, k check
    seed = TC_SEED + 9
    # (f) memory: sizeof(T) (s + d k + 2k) + 16; default s = d
    @test bytes_per_orbit(Val(2), Val(2), Float64, 4) == 8 * (4 + 4 + 4) + 16
    @test bytes_per_orbit(Val(2), Val(2), Float64, 4) == 112
    @test bytes_per_orbit(Val(2), Val(2), Float64) == 96
    @test bytes_per_orbit(Val(2), Val(2), Double64, 4) == 16 * 12 + 16
    free = 2^30
    @test chunk_size(free, Val(2), Val(2), Float64; s = 4) == fld(2^29, 112)
    @test chunk_size(free, Val(2), Val(2), Float64) == 5592405
    sys = tc_Long(2.0)
    rng = StableRNG(seed)
    x0s = [SVector(2pi * rand(rng), 2pi * rand(rng), 0.0, 0.0) for _ in 1:257]
    # (a) per-orbit body equals the manual iteration bitwise; x has 4 components, W is 2 x k
    # D38: k = 1 of (a) and (b) moved to verify/tangent-core.jl; k = 1 stays in (c).
    for k in (2,), tau in (1, 3)
        W0 = initial_deviations(seed, 5, Val(2), Val(k), Float64)
        x, W, st = propagate_orbit(sys, x0s[5], W0, 30, tau, ())
        xm, Wm = tc_manual(sys, x0s[5], W0, 30)
        @test x isa SVector{4, Float64} && W isa SMatrix{2, k, Float64}
        @test isequal(x, xm) && isequal(W, Wm)
        @test x[3] == 30 && x[4] != 0
    end
    # (b) chunked propagate equals the serial map of propagate_orbit, bitwise, for chunk sizes 16 and default
    for k in (2,)
        reducers = (tc_StateLen(), tc_LogNorm(), tc_Reader(3))   # D38: StateLen (a read-only reducer) joins the chain so that (e) shares this kernel specialization
        W0 = [initial_deviations(seed, j, Val(2), Val(k), Float64) for j in 1:257]
        serial = [propagate_orbit(sys, x0s[j], W0[j], 30, 3, reducers; orbit = j) for j in 1:257]
        for chunk in (16, nothing)
            res = propagate(sys, x0s, 30; k = k, tau = 3, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, [s[1] for s in serial])
            @test isequal(res.W, [s[2] for s in serial])
            @test isequal(res.finals, [map(reducer_final, reducers, s[3]) for s in serial])
        end
        res = propagate(sys, x0s, 30; k = k, tau = 3, seed = seed, reducers = reducers, chunk = 16)   # seeded W0
        @test isequal(res.W, [propagate_orbit(sys, x0s[j], W0[j], 30, 3, reducers)[2] for j in 1:257])
    end
    # (c) non-finite value in component 3 or 4: error "state" at handoff 10 (poison at step 7, tau = 5), lowest orbit id first
    for bad in (3, 4)
        psys = tc_Long(2.0, bad, 7)
        # D38: k = 2 with the chain of (b), to share its kernel specialization; the original k = 1, no-reducer call moved to verify
        err = tc_health_error(() -> propagate(psys, x0s[1:5], 20; k = 2, tau = 5, seed = seed, chunk = 2, reducers = (tc_StateLen(), tc_LogNorm(), tc_Reader(3))))
        @test err isa NumericalHealthError
        @test err.what == "state" && err.step == 10 && err.orbit == 1 && err.backend === :cpu
        err = tc_health_error(() -> propagate_orbit(psys, x0s[1], initial_deviations(seed, 1, Val(2), Val(1), Float64), 20, 5, ()))
        @test err isa NumericalHealthError && err.what == "state" && err.step == 10
    end
    # (d) zero allocation after warm-up
    W0 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
    @test tc_alloc_bytes(sys, x0s[1], W0, 1000, 1, (tc_LogNorm(), tc_Reader(1))) == 0
    @test tc_alloc_bytes(sys, x0s[1], W0, 100, 1, ()) == 0
    # (e) a reducer sees the full 4-component state
    res = propagate(sys, x0s[1:3], 6; k = 2, tau = 3, seed = seed, reducers = (tc_StateLen(), tc_LogNorm(), tc_Reader(3)))
    @test all(f -> f[1].len == 4, res.finals)
    # (g) k larger than the tangent dimension is rejected before any work
    @test_throws ArgumentError propagate(sys, x0s, 6; k = 3, tau = 1, seed = seed)
end
