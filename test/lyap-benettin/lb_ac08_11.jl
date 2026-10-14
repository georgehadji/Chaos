# AC8 to AC11: norm evaluation, numerical health, orbit and execution invariance, wrapped versus lifted state.

@testset "AC8 norm evaluation" begin
    # safe_norm of c (3, 4) for c = 2^600, 2^-600, 1, 2^1000 (and 2^-1000, the lower end of the contract range): finite, relative
    # error <= 2 eps against the BigFloat value 5 c. 200 seeded random vectors per type (d = 2 and 4, common scale 2^e, |e| <= 900): same bound.
    # Zero gives zero; NaN or Inf entries give a non-finite norm.
    rng = StableRNG(LB_SEED + LB_OFF + 80)
    for T in LB_TYPES
        epsT = lb_eps(T)
        relerr(v, ref) = setprecision(BigFloat, LB_BITS) do
            Float64(abs(lb_big(safe_norm(v)) - ref) / ref)
        end
        for k in (600, -600, 0, 1000, -1000)
            c = ldexp(one(T), k)
            v = SVector{2, T}(3c, 4c)
            @test isfinite(safe_norm(v))
            @test relerr(v, 5 * BigFloat(2)^k) <= 2 * epsT
        end
        for _ in 1:200
            d = rand(rng, (2, 4))
            e = rand(rng, -900:900)
            v = d == 2 ? SVector{2, T}(ldexp(T(2 * rand(rng) - 1), e), ldexp(T(2 * rand(rng) - 1), e)) :
                SVector{4, T}(ntuple(_ -> ldexp(T(2 * rand(rng) - 1), e), 4))
            ref = setprecision(BigFloat, LB_BITS) do
                sqrt(sum(x -> lb_big(x)^2, v))
            end
            iszero(ref) && continue
            @test relerr(v, ref) <= 2 * epsT
        end
        @test safe_norm(SVector{2, T}(0, 0)) == 0
        @test !isfinite(safe_norm(SVector{2, T}(T(NaN), 1)))
        @test !isfinite(safe_norm(SVector{2, T}(T(Inf), 1)))
    end
end

@testset "AC9 numerical health" begin
    # (a) tangent zero for flagged orbits: error at the first handoff, kind "deviation matrix", orbit, step, backend, column 1,
    #     value 0, nothing stored. Threshold of the zero case is not asserted (Q11). (b) non-finite W is tangent-core AC8; here only
    #     that the Benettin chain does not mask it (column 0). The reducer protocol for the failure is tested in "API".
    #     Float64 through propagate (kernel path, chunk 2); Double64 through the per-orbit body.
    seed = LB_SEED + LB_OFF + 90
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, nan_map = lb_Fault(:zero, 0), lb_Fault(:nan, 3)
    ok = propagate(zero_map, x0s(()), 12; k = 1, tau = 3, seed = seed, reducers = (Benettin(),), chunk = 2)
    @test length(ok.finals) == 6
    res = nothing
    err = lb_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 1, tau = 3, seed = seed,
                                                reducers = (Benettin(),), chunk = 2)))
    @test err isa NumericalHealthError
    @test res === nothing
    @test err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu
    @test err.column == 1 && err.value == 0.0
    res = nothing
    err = lb_health_error(() -> (res = propagate(zero_map, x0s((5, 2)), 12; k = 1, tau = 3, seed = seed,
                                                reducers = (Benettin(),), chunk = 2)))
    @test err isa NumericalHealthError && err.orbit == 2 && err.step == 3 && res === nothing   # lowest orbit id first
    res = nothing
    err = lb_health_error(() -> (res = propagate(nan_map, x0s((3,)), 20; k = 1, tau = 5, seed = seed,
                                                reducers = (Benettin(),), chunk = 2)))
    @test err isa NumericalHealthError && res === nothing
    @test err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.backend === :cpu
    @test err.column == 0 && err.threshold == Inf && !isfinite(err.value)
    x = SVector{3, Double64}(0, 0, 1)
    W0 = initial_deviations(seed, 1, Val(2), Val(1), Double64)
    err = lb_health_error(() -> propagate_orbit(zero_map, x, W0, 12, 3, (Benettin(),); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3
    @test err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = lb_health_error(() -> propagate_orbit(nan_map, x, W0, 20, 5, (Benettin(),); orbit = 3))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 3 && err.step == 5
    @test err.column == 0 && err.threshold == Inf && !isfinite(err.value)
    # Benettin needs k = 1: ArgumentError before any work (lb_Boom throws if the map is called).
    e = try
        propagate(lb_Boom(), [SVector(0.0, 0.0)], 4; k = 2, tau = 2, seed = seed, reducers = (Benettin(),))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("Benettin needs k = 1, got k = 2", e.msg)
    @test_throws ArgumentError reducer_init(Benettin(), Float64, Val(2), Val(2))
end

@testset "AC10 orbit and execution invariance" begin
    # (a) the orbit at every handoff is bitwise equal for the chains {}, {L1}, {L1, det}; (b) propagate equals the serial map of
    # propagate_orbit bitwise (x, W, finals) for chunks 64 and 1000: 4096 orbits, K = 2 pi, n = 1e3, tau in {1, 4}, Float64
    # (thread count is the run's, never changed); (c) 0 bytes per orbit after warm-up, n = 1e4, chain {L1, det}, both types.
    seed = LB_SEED + LB_OFF + 100
    T = Float64
    for tau in (1, 4)
        sys = StandardMap(lb_K(T, :two_pi))
        mon = ResidualMonitor(sys, LB_DET_TOL, Inf)
        for (j, s) in enumerate(lb_starts(StableRNG(seed + 1), 5))
            x0 = lb_cast(T, s)
            W0 = initial_deviations(seed, j, Val(2), Val(1), T)
            x_a, _, st_a = propagate_orbit(sys, x0, W0, 400, tau, (lb_Trace(),); orbit = j)
            x_b, _, st_b = propagate_orbit(sys, x0, W0, 400, tau, (Benettin(), lb_Trace()); orbit = j)
            x_c, _, st_c = propagate_orbit(sys, x0, W0, 400, tau, (Benettin(), mon, lb_Trace()); orbit = j)
            tr_a = reducer_final(lb_Trace(), st_a[1]).xs
            tr_b = reducer_final(lb_Trace(), st_b[2]).xs
            tr_c = reducer_final(lb_Trace(), st_c[3]).xs
            @test length(tr_a) == 400 ÷ tau
            @test isequal(tr_a, tr_b) && isequal(tr_a, tr_c)
            @test x_a === x_b && x_a === x_c
        end
    end
    for tau in (1, 4)
        sys = StandardMap(lb_K(T, :two_pi))
        reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf))
        m, N = 4096, 1000
        x0 = [lb_cast(T, s) for s in lb_starts(StableRNG(seed + 2), m)]
        W0 = [initial_deviations(seed, j, Val(2), Val(1), T) for j in 1:m]
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs = [s[1] for s in serial]
        Ws = [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, xs)
            @test isequal(res.W, Ws)
            @test isequal(res.finals, fs)
            @test res.provenance.estimator == "benettin"
        end
    end
    for T in LB_TYPES, tau in (1, 4)
        sys = StandardMap(lb_K(T, :two_pi))
        reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf))
        x0 = SVector{2, T}(1, 2)
        W0 = initial_deviations(seed, 1, Val(2), Val(1), T)
        @test lb_alloc_bytes(sys, x0, W0, 10_000, tau, reducers) == 0
    end
end

@testset "AC11 wrapped versus lifted state" begin
    # (a) K = 0, dyadic starts: wrapped (state length 4, D28) and lifted series bitwise identical, tolerance 0.
    # (b) regular class K = 0.5 (S7), n in {10, 100}, 50 starts, Float64: |n X_wrapped - n X_lifted| <= 64 K M n^2 eps rho_n (M = M_state as in AC7, Amendment D35).
    # Float64 only (one more argument-type combination per element type costs seconds of compile time).
    seed = LB_SEED + LB_OFF + 110
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        resL = lb_run(StandardMap(exact(T, 0)), zs, 300; k = 1, tau = tau, seed = seed)
        resW = lb_run(WrappedStandardMap(exact(T, 0)), torus_state.(zs), 300; k = 1, tau = tau, seed = seed)
        @test all(j -> isequal(resL.finals[j][1], resW.finals[j][1]), 1:50)
    end
    ids, starts, rejected, rseed = lb_regular()
    @info "AC11: regular-class screening" replaced = rejected
    W0 = [initial_deviations(rseed, id, Val(2), Val(1), Float64) for id in ids]
    for n in (10, 100)
        resL = lb_run(StandardMap(0.5), starts, n; k = 1, tau = 1, seed = seed, W0 = W0)
        resW = lb_run(WrappedStandardMap(0.5), torus_state.(starts), n; k = 1, tau = 1, seed = seed, W0 = W0)
        excluded = 0
        for j in 1:50
            fL, fW = resL.finals[j][1], resW.finals[j][1]
            _, rho, M = lb_ref_full(0.5, starts[j], W0[j][:, 1], n, 0.0, 0.0)
            if rho > 50
                excluded += 1
                continue
            end
            diff = abs((fW.sum[1] + fW.comp[1]) - (fL.sum[1] + fL.comp[1]))
            @test diff <= 64 * 0.5 * M * n^2 * eps(Float64) * rho
        end
        @test excluded <= 5
    end
end
