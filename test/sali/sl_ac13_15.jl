# AC13 to AC15 and the reducer-protocol tests (API). AC16 to AC18 are stage 04 (see the test plan).

@testset "AC13a zero tangent, column and orbit of the health error" begin
    _ = (SALI,)
    # A zero norm is not caught by the finiteness check of the stream: the reducer records it. Tangent zero for flagged orbits: error at the first handoff, kind
    # "deviation matrix", orbit id, step, backend, column 1, value 0, nothing stored. Column 2 collapsing alone (tangent diag(1, 0), frame I): column 2.
    # Float64 through propagate (kernel path, chunk 2); Double64 through the per-orbit body.
    seed = sl_seed(13, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map = sl_Fault(:zero, 0), sl_Fault(:col2, 0)
    ok = propagate(zero_map, x0s(()), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)
    @test length(ok.finals) == 6
    res = nothing
    err = sl_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)))
    @test err isa NumericalHealthError
    @test res === nothing
    @test err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu
    @test err.column == 1 && err.value == 0.0
    res = nothing
    err = sl_health_error(() -> (res = propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2)))
    @test err isa NumericalHealthError && err.orbit == 2 && err.step == 3 && res === nothing        # lowest orbit id first
    eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
    err = sl_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = (SALI(),), chunk = 2, W0 = eye))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 5 && err.step == 3
    @test err.column == 2 && err.value == 0.0
    x = SVector{3, Double64}(0, 0, 1)
    W0 = sl_frame(seed, 1, Double64)
    err = sl_health_error(() -> propagate_orbit(zero_map, x, W0, 12, 3, (SALI(),); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3
    @test err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = sl_health_error(() -> propagate_orbit(col2_map, x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (SALI(),); orbit = 4))
    @test err isa NumericalHealthError && err.column == 2 && err.step == 3
    # SALI needs k = 2: ArgumentError before any work (sl_Boom throws if the map is called).
    e = try
        propagate(sl_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (SALI(),))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("SALI needs k = 2, got k = 1", e.msg)
end

@testset "AC13b non-finite deviation matrix is not masked" begin
    _ = (SALI,)
    # Non-finite entries are tangent-core AC8; only that the SALI chain does not change that verdict (kind "deviation matrix", column 0, threshold Inf).
    seed = sl_seed(13, 2)
    nan_map = sl_Fault(:nan, 3)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    res = nothing
    err = sl_health_error(() -> (res = propagate(nan_map, x0s((3,)), 20; k = 2, tau = 5, seed = seed, reducers = (SALI(),), chunk = 2)))
    @test err isa NumericalHealthError && res === nothing
    @test err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.backend === :cpu
    @test err.column == 0 && err.threshold == Inf && !isfinite(err.value)
    x = SVector{3, Double64}(0, 0, 1)
    err = sl_health_error(() -> propagate_orbit(nan_map, x, sl_frame(seed, 1, Double64), 20, 5, (SALI(),); orbit = 3))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.column == 0
end

@testset "AC13c supplied frame check" begin
    _ = (SALI, sali_floor, check_frame)
    # A supplied frame is accepted only if SALI(0) >= 10 B_meas(T); otherwise a numerical health error of kind "initial frame" (step 0, the orbit id).
    # Dependent columns (SALI(0) = 0, identical or opposite) and SALI(0) = 5 B_meas are rejected; SALI(0) = 20 B_meas is accepted and check_frame returns SALI(0);
    # the seeded frame is accepted for every orbit id 1..1e4 (Double64: 1..1000).
    r = SALI()
    seed = sl_seed(13, 3)
    for T in SL_TYPES
        meas = sl_meas(T)
        frame(S) = (phi = 2 * asin(T(S) / 2); SMatrix{2, 2, T}(1, 0, cos(phi), sin(phi)))
        for W in (SMatrix{2, 2, T}(1, 0, 1, 0), SMatrix{2, 2, T}(1, 0, -1, 0), SMatrix{2, 2, T}(2, 3, 2, 3))
            e = sl_health_error(() -> check_frame(r, W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0
        end
        W5, W20 = frame(5 * meas), frame(20 * meas)
        S5, S20 = sl_sali_exact(W5[:, 1], W5[:, 2])[3], sl_sali_exact(W20[:, 1], W20[:, 2])[3]
        @test 4 * meas < S5 < 6 * meas && 19 * meas < S20 < 21 * meas           # the constructed frames have the intended SALI(0)
        e = sl_health_error(() -> check_frame(r, W5; orbit = 3))
        @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 3 && e.step == 0
        s0 = nothing
        @test sl_health_error(() -> (s0 = check_frame(r, W20; orbit = 3))) === nothing
        @test s0 !== nothing && sl_err(s0, S20) <= meas
        nid = T === Float64 ? 10_000 : 1_000
        rejected = count(j -> sl_health_error(() -> check_frame(r, sl_frame(seed, j, T); orbit = j)) !== nothing, 1:nid)
        @test rejected == 0
    end
end

@testset "AC13c rejection before any work" begin
    _ = (SALI, check_frame)
    # Amendment D67: the frame check runs before the stream, the tangent-core stream is unchanged. The caller checks every supplied frame first, so a
    # rejected frame never reaches the map (sl_Boom throws if forward or tangent is called); an accepted frame proceeds to the map.
    seed = sl_seed(13, 4)
    dep = SMatrix{2, 2, Float64}(1, 0, 1, 0)
    W0 = [sl_frame(seed, 1, Float64), dep]
    checked_run = () -> begin
        foreach(j -> check_frame(SALI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(sl_Boom(), SVector(0.0, 0.0), W0[2], 4, 2, (SALI(),))
    end
    err = sl_health_error(checked_run)
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 2 && err.step == 0
    err = sl_health_error(() -> check_frame(SALI(), dep; orbit = 5))
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 5 && err.step == 0
    @test sl_health_error(() -> check_frame(SALI(), W0[1]; orbit = 1)) === nothing
    @test_throws ErrorException propagate_orbit(sl_Boom(), SVector(0.0, 0.0), W0[1], 4, 2, (SALI(),))
end

@testset "AC13d overflow" begin
    _ = (SALI,)
    # A = diag(mu, 1/mu), mu = e^3, Float64 (as tangent-core AC8c): tau = 200 raises no error at n = 200 (the columns are normalised before they overflow);
    # tau = 300 raises the stream's error at handoff 300 (e^900 overflows). An exact SALI = 0 (both columns aligned with e_1) is not an error (A16).
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    W0 = sl_frame(sl_seed(13, 5), 1, Float64)
    x0 = SVector(0.0, 0.0)
    @test sl_health_error(() -> propagate_orbit(sys, x0, W0, 200, 200, (SALI(),))) === nothing
    err = sl_health_error(() -> propagate_orbit(sys, x0, W0, 300, 300, (SALI(),)))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 300 && err.backend === :cpu
end

@testset "AC13e run record" begin
    _ = (SALI, sali_floor)
    # The run record carries tau, precision and k (stream) and the indicator name (reducer_provenance); the floor B_meas(T) = 4 eps is a field of the reducer
    # final. The frame source and the chain order are recorded by Experiments, not by this module (Amendment D67).
    @test reducer_provenance(SALI()).indicator == "sali"
    sys = StandardMap(2.0)
    res = propagate(sys, [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = sl_seed(13, 6), reducers = (SALI(),))
    @test res.provenance.indicator == "sali"
    @test res.provenance.k == 2 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
    f = res.finals[1][1]
    @test f.n == 8 && f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64}
    for T in SL_TYPES
        _, _, st = propagate_orbit(StandardMap(exact(T, 2)), sl_cast(T, SVector(1.0, 2.0)), sl_frame(sl_seed(13, 6), 1, T), 8, 2, (SALI(),))
        @test reducer_final(SALI(), st[1]).floor == 4 * effective_eps(T)
        @test sali_floor(T) == 4 * effective_eps(T)
    end
end

@testset "AC14 orbit and execution invariance" begin
    _ = (SALI,)
    # (a) the orbit at every handoff is bitwise equal for the chains {}, {L1} on k = 1, {SALI} on k = 2, and for k = 1 vs k = 2 (L1 S2); (b) propagate equals the serial
    # map of propagate_orbit bitwise (x, W, finals) for chunk 64 and 1000: 4096 orbits, K = 2 pi, n = 1e3, tau in {1, 4}, Float64 (the thread count is the
    # run's, never changed here); (c) 0 bytes per orbit after warm-up, n = 1e4, both types; (d) absorbing zero (A16): once SALI is 0 it is 0 at every later
    # handoff, K = 2 pi, n = 1e3, Float64 (1000 starts).
    seed = sl_seed(14)
    T = Float64
    for tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        for (j, s) in enumerate(sl_starts(StableRNG(seed + 1), 5))
            x0 = sl_cast(T, s)
            W2 = sl_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (sl_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), sl_Trace()); orbit = j)
            xc, _, stc = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), sl_Trace()); orbit = j)
            xd, _, std = propagate_orbit(sys, x0, W2, 400, tau, (sl_Trace(),); orbit = j)
            tra, trb = reducer_final(sl_Trace(), sta[1]).xs, reducer_final(sl_Trace(), stb[2]).xs
            trc, trd = reducer_final(sl_Trace(), stc[2]).xs, reducer_final(sl_Trace(), std[1]).xs
            @test length(tra) == 400 ÷ tau
            @test isequal(tra, trb) && isequal(tra, trc) && isequal(tra, trd)
            @test xa === xb && xa === xc && xa === xd
        end
    end
    for tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        reducers = (SALI(),)
        m, N = 4096, 1000
        x0 = [sl_cast(T, s) for s in sl_starts(StableRNG(seed + 2), m)]
        W0 = [sl_frame(seed, j, T) for j in 1:m]
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, xs)
            @test isequal(res.W, Ws)
            @test isequal(res.finals, fs)
            @test res.provenance.indicator == "sali"
        end
    end
    for T in SL_TYPES, tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        @test sl_alloc_bytes(sys, sl_cast(T, SVector(1.0, 2.0)), sl_frame(seed, 1, T), 10_000, tau, (SALI(),)) == 0
    end
    m = 1000
    c = sl_case(:two_pi, Float64, m)
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000, 1; all = true)
    reached, violations = 0, 0
    for r in run
        sal = [s.sali for s in r.snaps]
        i = findfirst(==(0), sal)
        i === nothing && continue
        reached += 1
        violations += all(==(0), sal[i:end]) ? 0 : 1
    end
    @test violations == 0
    @test reached >= m ÷ 2                                       # the test is not vacuous: most chaotic starts reach the floor and stay there
end

@testset "AC15 wrapped versus lifted state" begin
    _ = (SALI,)
    # (a) K = 0, dyadic starts: wrapped (state length 4, D28) and lifted series bitwise identical (SALI, L, S_i, compensation), tolerance 0. (b) regular class K = 0.5
    # (S7), n in {10, 100}, Float64 (50 starts): |L_wrapped - L_lifted| <= 64 K M_state n^2 eps (rho_n^(1) + rho_n^(2)), rho_n^(i) as in L1
    # criterion 7 for vector i and M_state of L1 S1 (the L1 criterion 11(b) form, once per vector; borrowed, not derived for the pair).
    seed = sl_seed(15)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [sl_frame(seed, j, T) for j in 1:50]
        rl = sl_run(SALI(), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = sl_run(SALI(), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        @test all(j -> isequal(rl[j].snaps, rw[j].snaps), 1:50)
    end
    m = 50
    c = sl_case(:half, T, m)
    for n in (10, 100)
        rl = sl_run(SALI(), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = sl_run(SALI(), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0
        for j in 1:m
            _, rho1, M = lb_ref_full(0.5, c.x0[j], c.W0[j][:, 1], n, 0.0, 0.0)
            _, rho2, _ = lb_ref_full(0.5, c.x0[j], c.W0[j][:, 2], n, 0.0, 0.0)
            if max(rho1, rho2) > 50
                excluded += 1
                continue
            end
            @test abs(rw[j].snaps[1].L - rl[j].snaps[1].L) <= 64 * 0.5 * M * n^2 * eps(Float64) * (rho1 + rho2)
        end
        @test excluded <= m ÷ 10
    end
end

@testset "API pure functions" begin
    _ = (alignment_indices, sali_pair, sali_floor, sali_from_log)
    # sali_floor(T) = 4 eps (A14); alignment_indices returns (d-, d+) in T; sali_from_log is the A4 conversion with g = min(exp(L), 1) (no cancellation), clipped
    # and underflow-safe: relative error <= 8 eps (1 + |L|) (the rounding of L is magnified by exp), SALI of a log below -745 is 0 to Float64 resolution.
    for T in SL_TYPES
        epsT, meas = sl_eps(T), sl_meas(T)
        @test sali_floor(T) == 4 * effective_eps(T)
        a, b = SVector{2, T}(3, 4), SVector{2, T}(0, 5)
        dm, dp = alignment_indices(a, b)
        @test dm isa T && dp isa T
        @test sali_pair(a, b) == min(dm, dp)
        setprecision(BigFloat, SL_BITS) do
            @test abs(sl_big(dm) - sqrt(BigFloat(2) / 5)) <= meas && abs(sl_big(dp) - sqrt(BigFloat(18) / 5)) <= meas
        end
        for k in 1:14
            g = T(10)^(-k)
            S = sali_from_log(log(g))
            ex = setprecision(BigFloat, SL_BITS) do
                gb = sl_big(g)
                Float64(sqrt(BigFloat(2)) * gb / sqrt(1 + sqrt(1 - gb^2)))
            end
            @test S isa T
            @test abs(Float64(S) - ex) <= 8 * epsT * (1 + abs(Float64(log(g)))) * ex
        end
        @test abs(sali_from_log(T(0.1)) - sqrt(T(2))) <= 4 * epsT                  # g clipped to 1: the pair is orthogonal
        @test Float64(sali_from_log(T(-2000))) < 1e-300
    end
end

@testset "API reducer protocol" begin
    _ = (SALI,)
    # init is total except k != 2 (ArgumentError); the n = 0 final; step arithmetic on the columns (3, 4) and (0, 5): alpha = 5, 5, S_i = ln 5, L = -2 ln 5,
    # SALI = sqrt(2/5), W' = column-normalised matrix; a failure (zero or non-finite norm) is recorded without throwing in reducer_step, accumulators and W
    # unchanged, and thrown by reducer_update (orbit 0, backend :cpu); reducer_error carries orbit, backend, the failing column.
    @test_throws ArgumentError reducer_init(SALI(), Float64, Val(2), Val(1))
    for T in SL_TYPES
        epsT, meas = sl_eps(T), sl_meas(T)
        r = SALI()
        st0 = reducer_init(r, T, Val(2), Val(2))
        f0 = reducer_final(r, st0)
        @test f0.n == 0
        @test f0.sum isa SVector{2, T} && f0.comp isa SVector{2, T}
        @test all(iszero, f0.sum) && all(iszero, f0.comp) && f0.L == 0
        @test f0.floor == 4 * effective_eps(T)
        @test !reducer_failed(r, st0)
        x = SVector{2, T}(0, 0)
        W = SMatrix{2, 2, T}(3, 4, 0, 5)
        st1, W1 = reducer_step(r, st0, 7, x, W)
        f1 = reducer_final(r, st1)
        ln5 = setprecision(BigFloat, SL_BITS) do
            log(BigFloat(5))
        end
        @test f1.n == 7
        @test sl_err(f1.sum[1] + f1.comp[1], ln5) <= 2 * epsT && sl_err(f1.sum[2] + f1.comp[2], ln5) <= 2 * epsT
        @test sl_err(f1.L, -2 * ln5) <= 4 * epsT
        @test sl_err(f1.sali, setprecision(() -> sqrt(BigFloat(2) / 5), BigFloat, SL_BITS)) <= meas
        @test all(abs.(W1 - SMatrix{2, 2, T}(T(3) / 5, T(4) / 5, 0, 1)) .<= 2 * epsT)
        stu, Wu = reducer_update(r, st0, 7, x, W)
        @test isequal(stu, st1) && isequal(Wu, W1)
        st2, W2 = reducer_step(r, st1, 14, x, SMatrix{2, 2, T}(6, 8, 0, 10))
        f2 = reducer_final(r, st2)
        @test f2.n == 14
        @test sl_err(f2.sum[1] + f2.comp[1], ln5 + setprecision(() -> log(BigFloat(10)), BigFloat, SL_BITS)) <= 4 * epsT
        # failure: column 2 zero, then column 1 zero, then a non-finite column
        Wz = SMatrix{2, 2, T}(3, 4, 0, 0)
        stf, Wf = reducer_step(r, st2, 21, x, Wz)
        @test reducer_failed(r, stf) && !reducer_failed(r, st2)
        @test Wf === Wz
        ff = reducer_final(r, stf)
        @test ff.sum == f2.sum && ff.comp == f2.comp && ff.n == f2.n
        e = reducer_error(r, stf, 9, :cpu)
        @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 9 && e.step == 21
        @test e.backend === :cpu && e.column == 2 && e.value == 0.0
        e2 = sl_health_error(() -> reducer_update(r, st2, 21, x, Wz))
        @test e2 isa NumericalHealthError && e2.what == "deviation matrix" && e2.orbit == 0 && e2.backend === :cpu && e2.step == 21 && e2.column == 2
        st1z, _ = reducer_step(r, st2, 21, x, SMatrix{2, 2, T}(0, 0, 0, 5))
        @test reducer_failed(r, st1z) && reducer_error(r, st1z, 1, :cpu).column == 1
        stn, _ = reducer_step(r, st2, 21, x, SMatrix{2, 2, T}(3, 4, T(NaN), 5))
        @test reducer_failed(r, stn)
        sti, _ = reducer_step(r, st2, 21, x, SMatrix{2, 2, T}(3, 4, T(Inf), 5))
        @test reducer_failed(r, sti)
    end
    @test reducer_provenance(SALI()).indicator == "sali"
end

@testset "API prefix invariance" begin
    _ = (SALI,)
    # The state after handoff j is bitwise identical for every N >= j tau (same seed, orbit, start, tau): a run of length j tau returns exactly the series value
    # at j tau. All-handoff snapshots of one run of 400 steps against separate runs, tau in {1, 4}, both types.
    seed = sl_seed(15, 20)
    for T in SL_TYPES, tau in (1, 4)
        sys = StandardMap(sl_K(T, :two_pi))
        x0 = sl_cast(T, SVector(1.0, 2.0))
        W0 = sl_frame(seed, 3, T)
        probe = sl_Probe(SALI(); all = true)
        _, _, st = propagate_orbit(sys, x0, W0, 400, tau, (probe,); orbit = 3)
        snaps = reducer_final(probe, st[1]).snaps
        @test length(snaps) == 400 ÷ tau
        for j in (1, 5, 25, 100 ÷ tau)
            N = j * tau
            _, _, s2 = propagate_orbit(sys, x0, W0, N, tau, (SALI(),); orbit = 3)
            f = reducer_final(SALI(), s2[1])
            @test snaps[j].n == N && f.n == N
            @test isequal(snaps[j], f)
        end
    end
end
