# AC10 to AC12: health and recording, execution invariance and allocation, wrapped versus lifted state.

@testset "AC10 health and recording" begin
    seed = l2_seed(10)
    d0 = 1e-8
    # (a) d0 = 1e-20 at a start of magnitude ~1 in Float64 is not resolved: rejected before any work (the map throws if it is called)
    for call in (() -> propagate_shadow(l2_Boom(), [SVector(1.0, 2.0), SVector(1.0, 2.0)], 4; d0 = 1e-20, tau = 1, seed = seed),
                 () -> propagate_shadow_orbit(l2_Boom(), SVector(1.0, 2.0), SVector(1.0, 0.0), 4, 1, 1e-20; orbit = 5))
        err = l2_health_error(call)
        @test err isa NumericalHealthError && err.what == "separation" && err.step == 0 && err.backend === :cpu && err.value == 0.0
        @test err.orbit in (1, 5)
    end
    # (b) a map that collapses the pair to one point at the second handoff, orbit 3 only: error there, with orbit and step, nothing stored
    x0 = [SVector(0.5, i == 3 ? 1.0 : 0.0) for i in 1:5]
    U0 = fill(SVector(0.6, 0.8), 5)
    sysc = l2_Fault(:collapse, 2.25)
    err = l2_health_error(() -> propagate_shadow(sysc, x0, 4; d0 = d0, tau = 1, seed = seed, U0 = U0, chunk = 2))
    @test err isa NumericalHealthError && err.what == "separation" && err.orbit == 3 && err.step == 2 && err.backend === :cpu
    err = l2_health_error(() -> propagate_shadow_orbit(sysc, x0[3], U0[3], 4, 1, d0; orbit = 3))
    @test err isa NumericalHealthError && err.what == "separation" && err.orbit == 3 && err.step == 2
    ok = propagate_shadow(sysc, x0[[1, 2, 4, 5]], 4; d0 = d0, tau = 1, seed = seed, U0 = U0[1:4])
    @test all(f -> abs(f[1].X2) <= 1e-6 && f[1].n == 4, ok.finals)
    # (c) non-finite values in the reference, the shadow or both orbit: kind "state" as tangent-core AC8 (threshold Inf, non-finite value)
    for kind in (:inf_ref, :inf_shadow, :inf_both)
        sysi = l2_Fault(kind, 2.25)
        for call in (() -> propagate_shadow(sysi, x0, 4; d0 = d0, tau = 1, seed = seed, U0 = U0, chunk = 2),
                     () -> propagate_shadow_orbit(sysi, x0[3], U0[3], 4, 1, d0; orbit = 3))
            err = l2_health_error(call)
            @test err isa NumericalHealthError && err.what == "state" && err.orbit == 3 && err.step == 2 && err.backend === :cpu
            @test err.threshold == Inf && !isfinite(err.value)
        end
    end
    # (d) d0, tau and the estimator are part of the run record handed to the store (the grid index belongs to the Experiments record)
    for (T, mode) in ((Float64, :lifted), (Float64, :wrapped), (Double64, :wrapped))
        sys = l2_sys(mode, T(5))
        x0s = [l2_state(mode, T, z) for z in l2_starts(StableRNG(seed), 8)]
        dT = T === Float64 ? d0 : exact(Double64, "1e-8")
        p = propagate_shadow(sys, x0s, 20; d0 = dT, tau = 2, seed = seed).provenance
        @test p.estimator == "two-trajectory" && p.state_mode == string(mode) && p.k == 0
        @test p.d0 isa Float64 && p.d0 == Float64(dT) && p.tau == 2 && p.seed == seed && p.precision == string(T)
        @test p.backend in (:cpu, :gpu)
        @test all(k -> haskey(p, k), (:device, :backend_fallback_reason, :threads, :cpu, :chunk_size, :headroom, :gpu_chunks, :gpu_chunk_size))
    end
    @testset "AC10 contract checks" begin
        sys = StandardMap(2.0)
        x0s = [SVector(0.3, 0.4)]
        @test_throws ArgumentError propagate_shadow(sys, x0s, 7; d0 = d0, tau = 5, seed = seed)                 # N not a multiple of tau
        @test_throws ArgumentError propagate_shadow(sys, x0s, 0; d0 = d0, tau = 1, seed = seed)
        @test_throws ArgumentError propagate_shadow_orbit(sys, x0s[1], SVector(1.0, 0.0), 7, 5, d0)
        for bad in (0.0, -1e-6, NaN, Inf)
            @test_throws ArgumentError propagate_shadow(sys, x0s, 4; d0 = bad, tau = 1, seed = seed)
        end
        setprecision(BigFloat, 256) do
            @test_throws ArgumentError propagate_shadow(sys, x0s, 4; d0 = BigFloat("1e-8"), tau = 1, seed = seed)   # would narrow into Float64
        end
        # U0 is normalised and overrides the seeded direction; the seeded direction is the first column of the stream's frame
        a = propagate_shadow(sys, x0s, 100; d0 = 1e-6, tau = 1, seed = seed, U0 = [SVector(3.0, 4.0)]).finals[1][1].X2
        b = propagate_shadow(sys, x0s, 100; d0 = 1e-6, tau = 1, seed = seed, U0 = [SVector(0.6, 0.8)]).finals[1][1].X2
        @test abs(a - b) <= 1e-9
        D = l2_D2T(2.0, 1e-6, 1e3, Float64(effective_eps(Float64)))     # M_state bounded by 1e3 over 100 steps
        u = initial_deviations(seed, 1, Val(2), Val(1), Float64)[:, 1]
        c = propagate_shadow(sys, x0s, 100; d0 = 1e-6, tau = 1, seed = seed).finals[1][1].X2
        e = propagate_shadow(sys, x0s, 100; d0 = 1e-6, tau = 1, seed = seed, U0 = [u]).finals[1][1].X2
        @test abs(c - e) <= 2 * D
    end
end

function l2_alloc(sys, x0, u0, N, tau, d0)
    propagate_shadow_orbit(sys, x0, u0, N, tau, d0)
    return @allocated propagate_shadow_orbit(sys, x0, u0, N, tau, d0)
end

@testset "AC11 execution invariance" begin
    # K = 2 pi, 4096 seeded orbits, n = 1e3, d0 = 1e-8: propagate_shadow equals the serial map of propagate_shadow_orbit bitwise for chunk sizes 64, 1000
    # and the default; 0 bytes allocated per orbit after warm-up, n = 1e4, Float64 and Double64, wrapped and lifted
    seed = l2_seed(11)
    m, n, d0 = L2_FULL ? 4096 : 1024, 1000, 1e-8          # package suite 1024 orbits
    sys = WrappedStandardMap(L2_TWO_PI)
    x0 = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), m)]
    u0s = [SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
    serial = [propagate_shadow_orbit(sys, x0[j], u0s[j], n, 1, d0; orbit = j) for j in 1:m]
    for chunk in (L2_FULL ? (64, 1000, nothing) : (64, nothing))
        res = propagate_shadow(sys, x0, n; d0 = d0, tau = 1, seed = seed, chunk = chunk)
        @test isequal(res.x, [s[1] for s in serial])
        @test isequal(res.W, [s[2] for s in serial])
        @test isequal(res.finals, [(s[3],) for s in serial])
    end
    for T in (Float64, Double64), mode in (:wrapped, :lifted)
        dT = T === Float64 ? d0 : exact(Double64, "1e-8")
        s = l2_sys(mode, T(L2_TWO_PI))
        @test l2_alloc(s, l2_state(mode, T, SVector(1.0, 2.0)), SVector{2, T}(T(0.6), T(0.8)), 10_000, 1, dT) == 0
    end
end

@testset "AC12 wrapped versus lifted" begin
    # K = 5, n = 1e5, 100 starts, Float64, d0 = 1e-8, tau = 1: each mode agrees with L1 in the same mode within Delta_2T(d0) Psi_n/n + B_S^L1(n)/n (T5a, Amendment D35) with its own M_state
    # (wrapped 2 pi; lifted the largest coordinate magnitude of the run). Lifted points with E_ro > 0.5 are reported: they may run or end in a
    # "separation" health error (the double has no resolution left at |x| ~ 1e7), nothing else is allowed. Lifted points that are valid at n = 1e5
    # are tested; a lifted run at n = 1e3 (valid for every start) is added so that the lifted check cannot be vacuous.
    seed = l2_seed(12)
    K, d0, eps64 = 5.0, 1e-8, Float64(effective_eps(Float64))
    # package suite: the n = 1e5 run is wrapped only on 50 starts; verify/ runs both modes on 100 starts
    for (n, modes, nstart) in ((1000, (:lifted, :wrapped), 100), (100_000, L2_FULL ? (:wrapped, :lifted) : (:wrapped,), L2_FULL ? 100 : 50))
        starts = l2_starts(StableRNG(seed), nstart)
        for mode in modes
            sys = l2_sys(mode, K)
            x0s = [l2_state(mode, Float64, z) for z in starts]
            tested = reported = errored = 0
            if mode === :wrapped
                X = l2_X(sys, x0s, n, seed)
                X2 = l2_X2(sys, x0s, n, seed, d0)
                D = l2_D2T(K, d0, L2_TWO_PI, eps64)
                Psi = l2_Psis(sys, x0s, seed, Float64, n)
                @test all(j -> abs(X2[j] - X[j]) <= l2_pathbound(D, Psi[j], n, l2_BS(Float64, n, l2_Gfirst(sys, x0s[j], 1), K)), eachindex(x0s))
                continue
            end
            Xall = l2_X(sys, x0s, n, seed)          # one run for all starts: the seeded direction of start j is that of orbit id j (U0 below)
            for j in eachindex(x0s)
                M = l2_Mstate(sys, x0s[j], n)
                Ero = l2_Ero(K, d0, M, eps64)
                try
                    X = Xall[j]
                    X2 = propagate_shadow(sys, x0s[j:j], n; d0 = d0, tau = 1, seed = seed, U0 = [SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1])]).finals[1][1].X2
                    if Ero <= 0.5
                        tested += 1
                        Psi = lb_Psi(sys, x0s[j], initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1], n)
                        @test abs(X2 - X) <= l2_pathbound(l2_D2T(K, d0, M, eps64), Psi, n, l2_BS(Float64, n, l2_Gfirst(sys, x0s[j], 1), K))
                    else
                        reported += 1
                    end
                catch e
                    @test e isa NumericalHealthError && e.what == "separation"      # nothing but a health error
                    Ero > 0.5 ? (reported += 1) : (errored += 1)
                end
            end
            n == 1000 && @test tested == nstart && errored == 0
            @info "lyap-two-trajectory AC12: lifted, n = $n: tested $tested, reported (E_ro > 0.5) $reported, unexpected separation errors $errored"
        end
    end
end
