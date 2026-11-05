# AC17, AC18, AC20 (package part) and the API tests of the proposal. AC19 (ChaosTools), AC21 (GPU) and the timing part of AC20 are stage 04.
# Each testset starts by naming the new symbols it uses.

# ---- AC17: orbit and execution invariance (I7) ----
# Orbit trace (points at every handoff) and final point of a chain; `idx` is the position of the ga_Trace in the chain.
function ga_trace_run(sys, x0, W0, N, tau, chain, idx)
    x, _, st = propagate_orbit(sys, x0, W0, N, tau, chain)
    return Chaoslab.reducer_final(ga_Trace(), st[idx]).xs, x
end

@testset "AC17a orbit invariance across chains" begin
    _ = (GALI,)
    # The orbit points at every handoff are bitwise equal for the chains {empty}, {Benettin on k = 1}, {SALI on k = 2} and {the new reducer on K = 2, 3, 4} (a read-only trace
    # reducer reads the orbit; it precedes the replacer), tau in {1, 4}, 5 seeded chaotic starts of Y3 at GA_Y3B, n = 400, Float64 (L1 S2).
    rng = StableRNG(ga_seed(17))
    T = Float64
    sys = ga_sys(T, GA_Y3B)
    for tau in (1, 4), j in 1:5
        x0 = ga_state(T, GA_Y3B, ga_draw_chaotic(rng, GA_Y3B))
        W = Dict(K => ga_frame(ga_seed(17, 1), j, T, 4, K) for K in 1:4)
        base, xb = ga_trace_run(sys, x0, W[4], 400, tau, (ga_Trace(),), 1)
        @test length(base) == 400 ÷ tau
        runs = [ga_trace_run(sys, x0, W[1], 400, tau, (Benettin(), ga_Trace()), 2), ga_trace_run(sys, x0, W[2], 400, tau, (SALI(), ga_Trace()), 2)]
        for K in 2:4
            push!(runs, ga_trace_run(sys, x0, W[K], 400, tau, (ga_Trace(), GALI()), 1))
        end
        for (tr, x) in runs
            @test isequal(tr, base)
            @test x === xb
        end
    end
end

@testset "AC17b propagate equals the serial per-orbit body" begin
    _ = (GALI,)
    # propagate equals the serial map of propagate_orbit bitwise (x, W, finals) for chunk 64 and 1000, tau in {1, 4}, Y3 chaotic starts at GA_Y3B, K = 4, Float64 (D38: 512 orbits,
    # n = 200; spec 4096 orbits, n = 1e3); the thread count is the run's, never changed here (1 against N threads is stage 04); the run record names the indicator.
    seed = ga_seed(17, 2)
    rng = StableRNG(seed)
    T = Float64
    sys = ga_sys(T, GA_Y3B)
    m, N = ga_size(512, 4096), ga_size(200, 1000)
    x0 = [ga_state(T, GA_Y3B, ga_draw_chaotic(rng, GA_Y3B)) for _ in 1:m]
    W0 = [ga_frame(seed, j, T, 4, 4) for j in 1:m]
    reducers = (GALI(),)
    for tau in (1, 4)
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 4, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, xs)
            @test isequal(res.W, Ws)
            @test isequal(res.finals, fs)
            @test res.provenance.indicator == "gali"
        end
    end
end

@testset "AC17c no allocation per orbit" begin
    _ = (GALI,)
    # 0 bytes allocated by one propagate_orbit call after warm-up, Float64 and Double64, tau in {1, 4}: Y3 (K = 4, n = 1e4) and Y4 N = 100 (K = 4, fixed-size state of D105; D38: n = 300,
    # spec 1e4; Gaussian frames from a seeded StableRNG because a 200 x 200 static QR is out of scope).
    # D110: forced bounds checks (--check-bounds=yes, the CI default) move the MMatrix buffer to the heap; the claim is for default runs.
    forced = Base.JLOptions().check_bounds == 1
    rng = StableRNG(ga_seed(17, 3))
    for T in GA_TYPES, tau in (1, 4)
        if forced
            @test_skip ga_alloc_bytes(ga_sys(T, GA_Y3B), ga_state(T, GA_Y3B, (1.0, 2.0, 3.0, 4.0)), ga_frame(ga_seed(17, 4), 1, T, 4, 4), 10_000, tau, (GALI(),)) == 0
        else
            @test ga_alloc_bytes(ga_sys(T, GA_Y3B), ga_state(T, GA_Y3B, (1.0, 2.0, 3.0, 4.0)), ga_frame(ga_seed(17, 4), 1, T, 4, 4), 10_000, tau, (GALI(),)) == 0
        end
    end
    c = ga_y4(100, 3.0, 0.2)
    for T in GA_TYPES, tau in (1, 4)
        x0 = ga_state(T, c, ga_draw_chaotic(rng, c))
        if forced
            @test_skip ga_alloc_bytes(ga_sys(T, c), x0, ga_randset(rng, T, 200, 4; spread = false), ga_size(300, 10_000), tau, (GALI(),)) == 0
        else
            @test ga_alloc_bytes(ga_sys(T, c), x0, ga_randset(rng, T, 200, 4; spread = false), ga_size(300, 10_000), tau, (GALI(),)) == 0
        end
    end
end

@testset "AC17d pass-through reader first" begin
    _ = (GALI,)
    # The chain {pass-through reader, GALI} gives the same series bitwise as {GALI}: Y3 chaotic starts at GA_Y3B, K = 4, n = 300, tau in {1, 3}, Float64 and Double64.
    rng = StableRNG(ga_seed(17, 5))
    for T in GA_TYPES, tau in (1, 3), j in 1:3
        x0 = ga_state(T, GA_Y3B, ga_draw_chaotic(rng, GA_Y3B))
        W0 = ga_frame(ga_seed(17, 6), j, T, 4, 4)
        a = ga_run1(GALI(), ga_sys(T, GA_Y3B), x0, W0, 300, tau; all = true)
        probe = ga_Probe(GALI(); all = true)
        _, _, st = propagate_orbit(ga_sys(T, GA_Y3B), x0, W0, 300, tau, (ga_Trace(), probe))
        b = Chaoslab.reducer_final(probe, st[2])
        @test isequal(a.snaps, b.snaps)
    end
end

# ---- AC18: wrapped against lifted state (D28) ----
# Lifted run and wrapped run (the wrapped system with the torus state of the same start) of a K = d reducer, both probed at `steps`.
function ga_wl_runs(c::GACase, ::Type{T}, z, fid, N::Int, tau::Int, steps) where {T}
    d = ga_d(c)
    x0 = ga_state(T, c, z)
    W0 = ga_frame(fid[1], fid[2], T, d, d)
    wsys = ga_wsys(T, c)
    xw = torus_state(wsys, x0)
    return ga_run1(GALI(), ga_sys(T, c), x0, W0, N, tau; steps = steps, orbit = fid[2]), ga_run1(GALI(), wsys, xw, W0, N, tau; steps = steps, orbit = fid[2])
end

@testset "AC18a K = 0 shear, dyadic starts" begin
    _ = (GALI,)
    # Y3 (nu = kappa = mu = 0) and Y4 N = 2 (K = eps = 0), K = d: dyadic starts (multiples of 1/1024 in [0, 2]), G_k, L, sums and compensation bitwise identical in the lifted and the wrapped mode,
    # tau in {1, 3}, n in {10, 300}, Float64 and Double64 (D38: 10 starts per case, spec 50).
    rng = StableRNG(ga_seed(18))
    for T in GA_TYPES, c in (ga_y3(0.0, 0.0, 0.0), ga_y4(2, 0.0, 0.0)), tau in (1, 3)
        bad = 0
        for j in 1:ga_size(10, 50)
            z = [rand(rng, 0:2048) / 1024 for _ in 1:ga_d(c)]
            rl, rw = ga_wl_runs(c, T, z, (ga_seed(18, 1), j), 300, tau, [30, 300])
            for (a, b) in zip(rl.snaps, rw.snaps)
                bad += (isequal(a.gali, b.gali) && isequal(a.L, b.L) && isequal(a.sum, b.sum) && isequal(a.comp, b.comp)) ? 0 : 1
            end
        end
        @test bad == 0
    end
end

# K_sys (largest Lipschitz-type constant of the force), M_state (state magnitude along the lifted orbit) and sum_i rho_n^(i) (rho = sigma_max(Y_n) / norm(Y_n w_i)) for the
# bound of AC18(b) (the sali criterion 15 form, borrowed).
function ga_wl_bound(c::GACase, z, W0, n::Int)
    d = ga_d(c)
    sys, x0 = ga_sys(Float64, c), ga_state(Float64, c, z)
    Ksys = c.kind === :y3 ? max(abs(c.p[1]), abs(c.p[2])) + 2 * abs(c.p[3]) : abs(c.p[1]) + 4 * abs(c.p[2])
    M = max(1.0, abs(c.p[1]), abs(c.p[2]), abs(c.p[3]))
    Y = Matrix{Float64}(I, d, d)
    x = x0
    for _ in 1:n
        M = max(M, maximum(abs, x))
        Y = Matrix(tangent(sys, x)) * Y
        x = forward(sys, x)
    end
    M = max(M, maximum(abs, x))
    smax = opnorm(Y)
    rho = sum(smax / norm(Y * Float64.(W0[:, i])) for i in 1:d)
    return 64 * Ksys * M * n^2 * eps(Float64) * rho
end

@testset "AC18b wrapped against lifted on the orbits of AC9" begin
    _ = (GALI,)
    # K = d, Float64: |L_wrapped - L_lifted| <= 64 K_sys M_state n^2 eps sum_i rho_n^(i) (the sali criterion 15 form, borrowed, not derived), rho_n^(i) = sigma_max(Y_n) / norm(Y_n w_i).
    # The form bounds rounding differences along one orbit: wrapped and lifted orbits separate at the rate of the chaos (cm criteria 5 and 6), so the bound can only hold while
    # that separation is below the rounding scale. Orbits: Y3 orbit C and orbit D, 3 seeded regular starts at GA_Y3A and Y4 N = 2 at (-1, 0.3), n in {10, 100}; strongly chaotic starts
    # (GA_Y3B, Y4 N = 2 at (3, 0.2); D38: 4 and 2, spec 50) at n = 10 only (the sali criterion 15 uses the regular class too).
    rng = StableRNG(ga_seed(18, 2))
    T = Float64
    starts = Any[(GA_Y3A, GA_ORBIT_C, (10, 100)), (GA_Y3A, GA_ORBIT_D, (10, 100))]
    for _ in 1:3
        push!(starts, (GA_Y3A, ga_draw_regular(rng, GA_Y3A), (10, 100)))
        push!(starts, (ga_y4(2, -1.0, 0.3), ga_draw_regular(rng, ga_y4(2, -1.0, 0.3)), (10, 100)))
    end
    for _ in 1:ga_size(4, 25)
        push!(starts, (GA_Y3B, ga_draw_chaotic(rng, GA_Y3B), (10,)))
    end
    for _ in 1:ga_size(2, 25)
        push!(starts, (ga_y4(2, 3.0, 0.2), ga_draw_chaotic(rng, ga_y4(2, 3.0, 0.2)), (10,)))
    end
    worst = 0.0
    for (j, (c, z, ns)) in enumerate(starts)
        rl, rw = ga_wl_runs(c, T, z, (ga_seed(18, 3), j), maximum(ns), 1, collect(ns))
        W0 = ga_frame(ga_seed(18, 3), j, T, ga_d(c), ga_d(c))
        for (a, b, n) in zip(rl.snaps, rw.snaps, ns)
            worst = max(worst, abs(a.L - b.L) / ga_wl_bound(c, z, W0, n))
        end
    end
    @test worst <= 1
    @info "AC18b reported: worst ratio" worst
end

# ---- AC20: scale (package part) ----
@testset "AC20 scale: Y4 at N = 100 through the stream, and the chain at N = 10" begin
    _ = (GALI,)
    # Y4 at N = 100 (d = 200), open chain and ring, K = 4, through propagate with the fixed-size state (D105): G_k in [0, 1 + B_meas], nested (G_{k+1} <= G_k + 2 B_meas(k + 1)), finite
    # (D38: 4 orbits, n = 100; spec 200 orbits, n = 1e3; Gaussian frames from a seeded StableRNG, the 200 x 200 Haar QR is out of scope). The chain at N = 10 (d = 20), K = 4 and
    # K = d = 20, n <= 50, Float64 against the 256-bit product within the bounds of AC9 (B_G^k, B_L). The time per orbit-step against the L1 pass is stage 04.
    rng = StableRNG(ga_seed(19))
    T = Float64
    for ring in (false, true)
        c = ga_y4(100, 3.0, 0.2; ring = ring)
        m, N = ga_size(4, 200), ga_size(100, 1000)
        x0 = [ga_state(T, c, ga_draw_chaotic(rng, c)) for _ in 1:m]
        W0 = [ga_randset(rng, T, 200, 4; spread = false) for _ in 1:m]
        res = propagate(ga_sys(T, c), x0, N; k = 4, tau = 1, seed = ga_seed(19, 1), reducers = (GALI(),), W0 = W0)
        for f in res.finals
            g = f[1].gali
            @test f[1].n == N && length(g) == 3
            @test all(x -> isfinite(x) && 0 <= x <= 1 + ga_meas(T, 4), g)
            @test g[2] <= g[1] + 2 * ga_meas(T, 3) && g[3] <= g[2] + 2 * ga_meas(T, 4)
        end
    end
    c = ga_y4(10, 3.0, 0.2)
    z = ga_draw_chaotic(rng, c)
    b = ga_bundle(c, T, z, (ga_seed(19, 2), 1), 50)
    for K in (4, 20)
        r = ga_run1(GALI(), ga_sys(T, c), b.x0, ga_sm(b.W0[:, 1:K]), 50, 1; steps = b.cps)
        q = ga_check_orbit(b, r.snaps, K, T)
        @test q.G <= 1 && q.floorbad == 0
        K == 20 && @test q.L <= 1 && q.nonfinite == 0
    end
end

# ---- API tests of the proposal ----
const GA_HAND = (SMatrix{4, 4, Float64}(3, 4, 0, 0, 0, 5, 0, 0, 0, 0, 2, 0, 0, 0, 0, 7), [0.6, 0.6, 0.6])   # columns (3, 4, 0, 0), (0, 5, 0, 0), (0, 0, 2, 0), (0, 0, 0, 7); G_2 = G_3 = G_4 = 3/5

@testset "API init and the final at n = 0" begin
    _ = (GALI, gali_floor)
    # Isbits state for every 2 <= K <= d; reducer_final at n = 0: n 0, gali zero (SVector{K-1}), sum and comp zero (SVector{K}), L zero for K = d and nothing for K < d, floor = gali_floor(T, k)
    # for k = 2..K; not failed; K outside 2..d is an ArgumentError (AC7a).
    for T in GA_TYPES, (d, K) in ((2, 2), (4, 2), (4, 3), (4, 4), (6, 3))
        st = reducer_init(GALI(), T, Val(d), Val(K))
        @test isbitstype(typeof(st))
        f = reducer_final(GALI(), st)
        @test f.n == 0
        @test f.gali isa SVector{K - 1, T} && all(iszero, f.gali)
        @test f.sum isa SVector{K, T} && f.comp isa SVector{K, T} && all(iszero, f.sum) && all(iszero, f.comp)
        K == d ? (@test f.L isa T && iszero(f.L)) : (@test f.L === nothing)
        @test f.floor isa SVector{K - 1, T} && all(i -> f.floor[i] == gali_floor(T, i + 1), 1:(K - 1))
        @test !reducer_failed(GALI(), st)
    end
end

@testset "API one handoff against hand values" begin
    _ = (GALI, gali_values)
    # Columns (3, 4, 0, 0), (0, 5, 0, 0), (0, 0, 2, 0), (0, 0, 0, 7): norms 5, 5, 2, 7, unit columns (0.6, 0.8, 0, 0), e_2, e_3, e_4, G_2 = G_3 = G_4 = 0.6. One handoff at n = 1 with K = 4 and with K = 3:
    # the replaced W is the unit matrix within 4 eps, the sums S_i = ln 5, ln 5, ln 2, ln 7 (sum + comp) within 4 eps, L = -(sum S_i) (K = d only), G_k within 2 B_meas(T, k) of 0.6 (equal to
    # gali_values of the raw matrix); reducer_update returns the same as reducer_step.
    W64, expect = GA_HAND
    for T in GA_TYPES
        epsT = ga_eps(T)
        W = SMatrix{4, 4, T}(W64)
        U = SMatrix{4, 4, T}(T(3) / 5, T(4) / 5, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        tf = setprecision(() -> BigFloat(3) / 5, BigFloat, GA_BITS)
        lns = setprecision(() -> [log(BigFloat(5)), log(BigFloat(5)), log(BigFloat(2)), log(BigFloat(7))], BigFloat, GA_BITS)
        x = SVector{4, T}(0, 0, 0, 0)
        for K in (4, 3)
            st = reducer_init(GALI(), T, Val(4), Val(K))
            st1, W1 = reducer_step(GALI(), st, 1, x, W[:, 1:K] |> ga_sm)
            @test !reducer_failed(GALI(), st1)
            @test maximum(abs, Float64.(W1 - U[:, 1:K])) <= 4 * epsT
            f = reducer_final(GALI(), st1)
            @test f.n == 1 && f.gali isa SVector{K - 1, T}
            @test all(k -> ga_err(f.gali[k - 1], tf) <= 2 * ga_meas(T, k), 2:K)
            @test all(i -> ga_err(f.sum[i] + f.comp[i], lns[i]) <= 4 * epsT, 1:K)
            K == 4 ? (@test abs(Float64(f.L + sum(f.sum[i] + f.comp[i] for i in 1:4))) <= 8 * epsT * 4) : (@test f.L === nothing)
            su, Wu = reducer_update(GALI(), st, 1, x, W[:, 1:K] |> ga_sm)
            @test isequal(reducer_final(GALI(), su), f) && isequal(Wu, W1)
            g = gali_values(W[:, 1:K] |> ga_sm)
            @test all(k -> abs(g[k - 1] - f.gali[k - 1]) <= 2 * ga_meas(T, k), 2:K)
        end
    end
end

@testset "API failure in a column" begin
    _ = (GALI,)
    # A zero or non-finite column norm is recorded by reducer_step without throwing: reducer_failed true, W and the accumulators unchanged, reducer_final.n the last successful handoff (0 at a first
    # failure, 1 after a success); reducer_update throws reducer_error(r, st, 0, :cpu); the error is kind "deviation matrix", carries the first failing column, the value (0 or NaN), the step and the
    # orbit and backend given to reducer_error (threshold 0 for a zero norm, Inf for a non-finite one).
    T = Float64
    x = SVector{4, T}(0, 0, 0, 0)
    W = ga_frame(ga_seed(0, 10), 1, T, 4, 4)
    st = reducer_init(GALI(), T, Val(4), Val(4))
    Wz = SMatrix{4, 4, T}(W[:, 1]..., W[:, 2]..., 0, 0, 0, 0, W[:, 4]...)
    st1, W1 = reducer_step(GALI(), st, 5, x, Wz)
    @test reducer_failed(GALI(), st1) && W1 === Wz
    f = reducer_final(GALI(), st1)
    @test f.n == 0 && all(iszero, f.sum) && all(iszero, f.comp)
    e = ga_health_error(() -> reducer_update(GALI(), st, 5, x, Wz))
    @test e isa NumericalHealthError && e.what == "deviation matrix" && e.column == 3 && e.step == 5 && e.value == 0.0 && e.orbit == 0 && e.backend === :cpu
    e = reducer_error(GALI(), st1, 7, :gpu)
    @test e isa NumericalHealthError && e.orbit == 7 && e.backend === :gpu && e.column == 3 && e.threshold == 0.0
    Wn = SMatrix{4, 4, T}(W[:, 1]..., NaN, W[2, 2], W[3, 2], W[4, 2], W[:, 3]..., W[:, 4]...)
    e = reducer_error(GALI(), reducer_step(GALI(), st, 2, x, Wn)[1], 1, :cpu)
    @test e.column == 2 && isnan(e.value) && e.threshold == Inf && e.what == "deviation matrix"
    Wzz = SMatrix{4, 4, T}(0, 0, 0, 0, W[:, 2]..., 0, 0, 0, 0, W[:, 4]...)
    @test reducer_error(GALI(), reducer_step(GALI(), st, 2, x, Wzz)[1], 1, :cpu).column == 1        # the first failing column
    sa, Wa = reducer_step(GALI(), st, 1, x, W)                                                     # a failure after a success keeps the earlier state
    sb, Wb = reducer_step(GALI(), sa, 2, x, Wz)
    @test reducer_failed(GALI(), sb) && Wb === Wz
    fa, fb = reducer_final(GALI(), sa), reducer_final(GALI(), sb)
    @test fb.n == 1 && fb.sum == fa.sum && fb.comp == fa.comp && fb.gali == fa.gali
    @test reducer_step(GALI(), sb, 3, x, W)[2] === W                                               # a failed state is absorbing: it neither advances nor replaces W
end

@testset "API provenance and the run record" begin
    _ = (GALI, gali_floor)
    # reducer_provenance(GALI()).indicator == "gali"; the propagate record carries the indicator, tau, k and the precision; the reducer final carries the floors; the frame source and the chain
    # are recorded by Experiments, not by this module (as sali Amendment D67).
    @test reducer_provenance(GALI()).indicator == "gali"
    sys = ga_sys(Float64, GA_Y3B)
    res = propagate(sys, [ga_state(Float64, GA_Y3B, (1.0, 2.0, 0.5, 0.25))], 8; k = 3, tau = 2, seed = ga_seed(0, 11), reducers = (GALI(),))
    @test res.provenance.indicator == "gali"
    @test res.provenance.k == 3 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
    f = res.finals[1][1]
    @test f.n == 8 && f.sum isa SVector{3, Float64} && f.comp isa SVector{3, Float64} && f.gali isa SVector{2, Float64}
    @test f.floor == SVector(gali_floor(Float64, 2), gali_floor(Float64, 3))
end

@testset "API gali_floor" begin
    _ = (gali_floor,)
    # B_meas(T, k) = 2 k effective_eps(T) for every element type (Float64, Double64 with 2^-100, BigFloat at the precision in force), returned in T.
    for T in (Float64, Double64, BigFloat), k in 2:6
        v = gali_floor(T, k)
        @test v isa T
        @test v == 2 * k * effective_eps(T)
    end
    @test gali_floor(Float64, 2) == 4 * 2.0^-52 && gali_floor(Double64, 3) == 6 * Double64(2.0^-100)
end

@testset "API gali_values and gali_value" begin
    _ = (gali_values, gali_value)
    # gali_values: SVector{K-1, T} of G_2..G_K for an SMatrix (also BigFloat), the same numbers for the dense Matrix; gali_value(W, k) = G_k of the first k columns for k in 2..K, an ArgumentError
    # outside; hand values 0.6; an orthonormal set gives 1 for every order.
    W64, expect = GA_HAND
    for T in (Float64, Double64, BigFloat)
        W = SMatrix{4, 4, T}(W64)
        g = gali_values(W)
        @test g isa SVector{3, T}
        @test all(k -> ga_err(g[k - 1], setprecision(() -> BigFloat(3) / 5, BigFloat, GA_BITS)) <= 2 * ga_meas(T, k), 2:4)
        gm = gali_values(Matrix(W))
        @test length(gm) == 3 && all(k -> abs(gm[k - 1] - g[k - 1]) <= 2 * ga_meas(T, k), 2:4)
        @test all(k -> gali_value(W, k) == gali_values(ga_sm(W[:, 1:k]))[k - 1], 2:4)
        @test_throws ArgumentError gali_value(W, 1)
        @test_throws ArgumentError gali_value(W, 5)
        Q = ga_frame(ga_seed(0, 12), 1, T, 4, 4)
        @test all(k -> abs(gali_values(Q)[k - 1] - 1) <= 2 * ga_meas(T, k), 2:4)
    end
    @test gali_values(SMatrix{3, 2, Float64}(1, 0, 0, 0, 1, 0)) isa SVector{1, Float64}
end
