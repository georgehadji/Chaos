# AC17 to AC19 and the reducer-protocol tests (API, prefix invariance). AC12 to AC16 are stage 04 (see the test plan).

@testset "AC17 regular-orbit decay" begin
    # K = 0.5 regular class (S7), 20 starts, n up to 1e6 (tau = 1): least-squares slope of ln X(n) against ln n over 21 log-spaced
    # n in [1e4, 1e6] lies in [-1.1, -0.8] for at least 80 percent of the starts (E5, Skokos2010 Eq. (62); loose sanity, K18).
    ids, starts, rejected, seed = lb_regular()
    ids, starts = ids[1:20], starts[1:20]
    @info "AC17: regular-class screening (first 50 regular starts shared with AC7 and AC11)" replaced = rejected
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    sys = StandardMap(0.5)
    good = 0
    for (id, s) in zip(ids, starts)
        probe = lb_Probe(Benettin(), ns; alphas = false)
        W0 = initial_deviations(seed, id, Val(2), Val(1), Float64)
        _, _, st = propagate_orbit(sys, SVector(s[1], s[2]), W0, 1_000_000, 1, (probe,); orbit = id)
        fin = reducer_final(probe, st[1])
        @test length(fin.snaps) == length(ns)
        X = [(S + C) / n for (n, S, C) in fin.snaps]
        all(>(0), X) || continue                       # a non-positive X has no logarithm: counts as outside the interval
        u = log.(Float64[n for (n, _, _) in fin.snaps])
        v = log.(X)
        slope = sum((u .- sum(u) / length(u)) .* (v .- sum(v) / length(v))) / sum((u .- sum(u) / length(u)) .^ 2)
        -1.1 <= slope <= -0.8 && (good += 1)
    end
    @test good >= 16
end

@testset "AC18 range invariant" begin
    # Every alpha_j of the runs of criteria 5 and 6 lies in [sigma^-tau, sigma^tau] within a relative 16 eps tau Phi_j (E14);
    # alpha_j is the norm the estimator sees (lb_Probe, before the inner step), Phi_j the squared Frobenius norm of the window product.
    rows5 = lb_ac5_data()
    rows6 = lb_ac6_data()
    @test !isempty(rows5) && !isempty(rows6)
    for r in rows5
        @test r.rv <= 1
    end
    for r in rows6
        @test r.rv <= 1
    end
end

@testset "AC19 symplectic dimension 2" begin
    # (a) chain (Benettin, ResidualMonitor): no health error in any run of criterion 5 (tau = 1 and 4), det_max under
    #     4 eps (1 + max|M|) <= 4 eps (2 + K) and nonzero over the ensemble. (b) injected det fault of 1 + delta with delta in [2 theta, 4 theta):
    #     "det" error with value delta, orbit 1, step 1; delta in (theta/4, theta/2]: no error over 100 steps.
    seed = LB_SEED + LB_OFF + 190
    for (ci, (label, N)) in enumerate(LB_AC5_CASES), T in LB_TYPES, tau in (1, 4)
        T === Double64 && tau == 4 && continue         # runtime: Double64 at tau = 1 and 8 starts
        K = lb_K(T, label)
        sys = StandardMap(K)
        x0 = [lb_cast(T, s) for s in lb_starts(StableRNG(seed + ci), T === Double64 ? 8 : 20)]
        res = lb_run(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf)), probe = false)
        dmax = maximum(f[2].det_max for f in res.finals)
        @test dmax <= LB_DET_TOL * effective_eps(T) * (2 + K)
        @test dmax > 0
    end
    for T in LB_TYPES
        A0 = SMatrix{2, 2, T}(2, 1, 1, 1)
        theta = Float64(det_threshold(LB_DET_TOL, T, A0))
        above = [exp2(j) for j in -130:20 if 2theta <= exp2(j) < 4theta]
        below = [exp2(j) for j in -130:20 if theta / 4 < exp2(j) <= theta / 2]
        @test !isempty(above) && !isempty(below)
        run(delta, N) = begin
            sys = LinearSymplecticMap(SMatrix{2, 2, T}(2 + T(delta), 1, 1, 1))
            mon = ResidualMonitor(sys, LB_DET_TOL, Inf)
            lb_health_error(() -> lb_run(sys, [SVector{2, T}(0, 0)], N; k = 1, tau = 1, seed = seed,
                                            reducers = (Benettin(), mon), W0 = [SMatrix{2, 1, T}(1, 0)], probe = false))
        end
        for delta in above
            e = run(delta, 1)
            @test e isa NumericalHealthError && e.what == "det" && e.value == delta && e.orbit == 1 && e.step == 1
        end
        for delta in below
            @test run(delta, 100) === nothing
        end
    end
end

@testset "API reducer protocol and provenance" begin
    # Batch B contract, Benettin: init total except k != 1; n = 0 final; step arithmetic (alpha = norm, W' = W / alpha, compensated
    # sum of ln alpha, n = handoff step); failure recorded without throwing in reducer_step, accumulators and W unchanged, thrown by
    # reducer_update with orbit 0 and backend :cpu; reducer_error carries orbit, backend, column 1; provenance.
    for T in LB_TYPES
        epsT = lb_eps(T)
        r = Benettin()
        st0 = reducer_init(r, T, Val(2), Val(1))
        f0 = reducer_final(r, st0)
        @test f0.n == 0
        @test f0.sum isa SVector{1, T} && f0.comp isa SVector{1, T} && f0.lambda isa SVector{1, T}
        @test f0.sum[1] == 0 && f0.comp[1] == 0 && f0.lambda[1] == 0
        @test !reducer_failed(r, st0)
        x = SVector{2, T}(0, 0)
        st1, W1 = reducer_step(r, st0, 7, x, SMatrix{2, 1, T}(3, 4))
        f1 = reducer_final(r, st1)
        ref = setprecision(BigFloat, LB_BITS) do
            log(BigFloat(5))
        end
        @test f1.n == 7
        @test lb_err(f1.sum[1], f1.comp[1], ref) <= 2 * epsT
        @test abs(W1[1] - T(3) / 5) <= 2 * epsT && abs(W1[2] - T(4) / 5) <= 2 * epsT
        @test abs(f1.lambda[1] - (f1.sum[1] + f1.comp[1]) / 7) <= 2 * epsT * abs(f1.lambda[1])
        stu, Wu = reducer_update(r, st0, 7, x, SMatrix{2, 1, T}(3, 4))
        @test isequal(stu, st1) && isequal(Wu, W1)
        st2, W2 = reducer_step(r, st1, 14, x, SMatrix{2, 1, T}(6, 8))
        ref2 = setprecision(BigFloat, LB_BITS) do
            log(BigFloat(5)) + log(BigFloat(10))
        end
        f2 = reducer_final(r, st2)
        @test f2.n == 14 && lb_err(f2.sum[1], f2.comp[1], ref2) <= 4 * epsT
        # failure: zero deviation vector
        Wz = SMatrix{2, 1, T}(0, 0)
        stf, Wf = reducer_step(r, st2, 21, x, Wz)
        @test reducer_failed(r, stf) && !reducer_failed(r, st2)
        @test Wf === Wz
        ff = reducer_final(r, stf)
        @test ff.sum == f2.sum && ff.comp == f2.comp
        e = reducer_error(r, stf, 9, :cpu)
        @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 9 && e.step == 21
        @test e.backend === :cpu && e.column == 1 && e.value == 0.0
        e2 = lb_health_error(() -> reducer_update(r, st2, 21, x, Wz))
        @test e2 isa NumericalHealthError && e2.what == "deviation matrix" && e2.orbit == 0 && e2.backend === :cpu && e2.step == 21
        # non-finite deviation vector is a failure too (not an exception inside reducer_step)
        stn, _ = reducer_step(r, st2, 21, x, SMatrix{2, 1, T}(T(NaN), 1))
        @test reducer_failed(r, stn)
    end
    @test reducer_provenance(Benettin()) == (estimator = "benettin",)
    sys = StandardMap(2.0)
    x0 = [SVector(1.0, 2.0)]
    for reducers in ((Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf)),)      # same argument types as AC10: one kernel compile
        res = propagate(sys, x0, 8; k = 1, tau = 2, seed = LB_SEED + LB_OFF + 200, reducers = reducers)
        @test res.provenance.estimator == "benettin"
        @test res.provenance.k == 1 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
        f = res.finals[1][1]
        @test f.n == 8 && f.lambda[1] == (f.sum[1] + f.comp[1]) / f.n
    end
end

@testset "API prefix invariance" begin
    # Q9: the state after handoff j is bitwise identical for every N >= j tau (same seed, orbit, start, tau), so a run of length
    # j tau returns exactly the series value at j tau. Series from the all-handoffs probe against separate runs of length j tau.
    seed = LB_SEED + LB_OFF + 210
    for T in LB_TYPES, tau in (1, 4)
        sys = StandardMap(lb_K(T, :two_pi))
        x0 = lb_cast(T, SVector(1.0, 2.0))
        W0 = initial_deviations(seed, 3, Val(2), Val(1), T)
        probe = lb_Probe(Benettin(); all = true)
        _, _, st = propagate_orbit(sys, x0, W0, 400, tau, (probe,); orbit = 3)
        snaps = reducer_final(probe, st[1]).snaps
        @test length(snaps) == 400 ÷ tau
        for j in (1, 5, 25, 100 ÷ tau)
            N = j * tau
            f = lb_final(sys, x0, W0, N, tau)
            @test snaps[j][1] == N && f.n == N
            @test isequal(snaps[j][2], f.sum[1]) && isequal(snaps[j][3], f.comp[1])
        end
    end
end
