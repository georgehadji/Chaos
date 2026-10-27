# AC10 to AC13: E2 on the stream (transparency, chain, memory, run length, against its reference), health / degenerate series / arguments, execution invariance, precision ladder.
# The stream reducer is WeightedPhaseAdvance (a reader: never replaces columns, never changes the orbit). A run needs 2 N + 1 handoffs (the first handoff is the sample 0 of the series; see the plan, Items for
# the owner 1). Chains whose other members are readers end with a test-side replacer (nf_Norm) so that a long run of an expanding map does not overflow.

@testset "AC10a transparency and chain" begin
    _ = (WeightedPhaseAdvance, LiftedAngleObservable, StretchingNumbers, stretching_range, FLI, SALI)
    # K_std = 2 pi, N = 1000, both types: orbit points bitwise equal at every handoff for the chains {Norm}, {E2, Norm}, {E2 then StretchingNumbers last} and {FLI(replace = false), E2, SALI}; the other
    # reducers' outputs bitwise their stand-alone values; the final columns of the chains with the same replacer are bitwise equal.
    N = 1000
    nsteps = 2N + 1
    for T in NF_TYPES
        c = nf_case(:two_pi, T, 1, 1)
        sys, x0 = c.sys, c.x0[1]
        W0 = nf_frame(nf_seed(10, 1), 1, T)
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        st = StretchingNumbers(; range = stretching_range(sys), bins = 100, interval = 100)
        xA, WA, fA = nf_run((nf_Trace(), nf_Norm()), sys, x0, W0, nsteps, 1)
        xB, WB, fB = nf_run((nf_Trace(), r, nf_Norm()), sys, x0, W0, nsteps, 1)
        xC, WC, fC = nf_run((nf_Trace(), r, st), sys, x0, W0, nsteps, 1)
        xD, WD, fD = nf_run((nf_Trace(), FLI(; replace = false), r, SALI()), sys, x0, W0, nsteps, 1)
        _, WS1, fS1 = nf_run((st,), sys, x0, W0, nsteps, 1)
        _, WS2, fS2 = nf_run((FLI(),), sys, x0, W0, nsteps, 1)
        _, WS3, fS3 = nf_run((SALI(),), sys, x0, W0, nsteps, 1)
        @test length(fA[1].xs) == nsteps
        @test all(f -> isequal(f[1].xs, fA[1].xs), (fB, fC, fD)) && xA == xB == xC == xD
        @test isequal(WA, WB) && isequal(WC, WS1) && isequal(WD, WS3)
        @test isequal(fC[3], fS1[1]) && isequal(fD[2], fS2[1]) && isequal(fD[4], fS3[1])
        # the two E2 results of the chains agree bitwise as well (the reader sees the same orbit)
        @test isequal(fB[2], fC[2]) && isequal(fB[2], fD[3])
        @test fB[2].complete && fB[2].handoffs == nsteps && fB[2].n == nsteps
    end
end

@testset "AC10b state size, isbits, no allocation" begin
    _ = (WeightedPhaseAdvance, LiftedAngleObservable, reducer_bytes)
    # The reducer bytes (state plus workspace, no workspace here) are at most the stated state size for one declared horizon (96 bytes Float64, 176 Double64) and grow by at most the stated increment per
    # extra horizon (64 and 128); the state is isbits; the CPU path allocates nothing after warm-up (tangent-core AC5 form).
    for (T, one_h, extra) in ((Float64, 96, 64), (Double64, 176, 128))
        r1 = WeightedPhaseAdvance(; window = 100, observable = LiftedAngleObservable(1, 3))
        r2 = WeightedPhaseAdvance(; window = (100, 200), observable = LiftedAngleObservable(1, 3))
        r3 = WeightedPhaseAdvance(; window = (100, 200, 400), observable = LiftedAngleObservable(1, 3))
        b1, b2, b3 = (reducer_bytes(r, T, Val(2), Val(2)) for r in (r1, r2, r3))
        @test b1 <= one_h && b2 - b1 <= extra && b3 - b2 <= extra
        @test isbitstype(typeof(reducer_init(r3, T, Val(2), Val(2)))) && reducer_workspace(r1, T, Val(2), Val(2)) == 0
        @info "AC10b reducer bytes (reported)" T one_horizon = b1 two = b2 three = b3
    end
    sys = WrappedStandardMap(nf_K(Float64, :two_pi))
    x0 = torus_state(SVector(1.0, 0.5))
    W0 = nf_frame(nf_seed(10, 2), 1, Float64)
    r = WeightedPhaseAdvance(; window = 100, observable = LiftedAngleObservable(1, 3))
    chain = (r, FLI())
    propagate_orbit(sys, x0, W0, 201, 1, chain)
    @test (@allocated propagate_orbit(sys, x0, W0, 201, 1, chain)) == 0
end

@testset "AC10c run length and tau" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, check_run, weighted_phase_frequency, phase_increments)
    # check_run(r, tau, N) is the caller's check before the stream (reducer_init receives neither tau nor N): an ArgumentError naming window unless the run has exactly 2 window + 1 handoffs (the largest
    # window of a tuple). With tau = 4 the windows count handoffs and the result equals E2 on the stored series sampled every 4 steps (planar observable: bitwise, one shared implementation).
    r = WeightedPhaseAdvance(; window = 100, observable = PlanarObservable(1, 2))
    @test check_run(r, 1, 201) === nothing && check_run(r, 4, 804) === nothing
    for (tau, n) in ((1, 200), (1, 202), (4, 800), (4, 805), (3, 603 + 1))
        e = nf_error(() -> check_run(r, tau, n))
        @test e isa ArgumentError && occursin("window", e.msg)
    end
    r2 = WeightedPhaseAdvance(; window = (50, 100), observable = PlanarObservable(1, 2))
    @test check_run(r2, 1, 201) === nothing && nf_error(() -> check_run(r2, 1, 101)) isa ArgumentError
    for T in NF_TYPES
        c = nf_case(:two_pi, T, 1, 1)
        sys, x0 = c.sys, c.x0[1]
        N, tau = 100, 4
        series = nf_obs(nf_obs_planar, nf_states_tau(sys, x0, 2N + 1, tau))
        nu1, nu2 = weighted_phase_frequency(phase_increments(series), N)
        _, _, fin = nf_run((WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2)), nf_Norm()), sys, x0, nf_frame(nf_seed(10, 3), 1, T), tau * (2N + 1), tau)
        f = fin[1]
        @test f.complete && f.handoffs == 2N + 1 && f.n == tau * (2N + 1)
        @test isequal(f.nu1[1], nu1) && isequal(f.nu2[1], nu2)
    end
end

@testset "AC10d the stream against NF6 at the realized samples" begin
    _ = (WeightedPhaseAdvance, LiftedAngleObservable)
    # K_std in {0.5 (regular class), 2 pi, 10}, N = 1000, both types: the stream frequencies of both windows are within B_E2 = 8 eps (1 + max_j |Delta_j| / (2 pi)) of the 512-bit evaluation of NF6 at the realized
    # samples (the lifted increments of the trace of the same run). D38: 4 starts per class; FULL: 20.
    N = 1000
    nstart = NF_FULL ? 20 : 4
    for T in NF_TYPES, (cl, sub) in ((:half, 0), (:two_pi, 1), (:ten, 3))
        epsT = nf_eps(T)
        c = nf_case(cl, T, nstart, sub)
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        worst = 0.0
        bad = String[]
        for j in 1:nstart
            _, _, fin = nf_run((nf_Trace(), r, nf_Norm()), c.sys, c.x0[j], nf_frame(nf_seed(10, 4), c.ids[j], T), 2N + 1, 1)
            xs = nf_vec(fin[1].xs)
            incr = nf_lifted_incr_big(xs)
            ref = nf_E2_ref(incr, N; bits = 512)
            B = nf_BE2(T, incr)
            f = fin[2]
            for (k, v) in enumerate((f.nu1[1], f.nu2[1]))
                e = Float64(abs(nf_big(v) - ref[k]))
                worst = max(worst, e / B)
                e <= B || push!(bad, "T=$T $cl start $j window $k error $(e / epsT) eps, bound $(B / epsT) eps")
            end
        end
        @test nf_report(bad, "AC10d $T $cl")
        @info "AC10d $T $cl: worst error / B_E2 (realized ratio for stage 04)" worst
    end
end

@testset "AC11a non-finite samples" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable, frequency_map, check_frequency_map)
    # Reducer: a non-finite orbit point at the handoff 4 is recorded in the state (no throw in reducer_step), the state is failed and unchanged by later handoffs, W is untouched, the host error carries
    # the orbit id, the map step and plane 1 ("series"); reducer_update throws at once with orbit 0. Series path: a non-finite sample is recorded with plane, window and sample index; the host check raises
    # with the orbit id. (A failed chunk storing nothing is the Experiments feed, stage 04.)
    for (obs, good, bad) in ((PlanarObservable(1, 2), SVector(0.3, 0.4), SVector(NaN, 0.4)), (LiftedAngleObservable(1, 3), SVector(0.3, 0.4, 1.0, 0.0), SVector(Inf, 0.4, 1.0, 0.0)))
        r = WeightedPhaseAdvance(; window = 16, observable = obs)
        W = SMatrix{2, 1, Float64}(1.0, 0.0)
        st = reducer_init(r, Float64, Val(2), Val(1))
        for n in 1:3
            st, W2 = reducer_step(r, st, n, good, W)
            @test W2 === W
        end
        st4, W4 = reducer_step(r, st, 4, bad, W)
        @test reducer_failed(r, st4) && !reducer_failed(r, st) && W4 === W
        e = reducer_error(r, st4, 9, :cpu)
        @test e isa NumericalHealthError && e.what == "series" && e.orbit == 9 && e.step == 4 && e.backend === :cpu && e.column == 1
        st5, _ = reducer_step(r, st4, 5, good, W)
        @test reducer_failed(r, st5) && isequal(st5, st4)
        ex = nf_health_error(() -> reducer_update(r, st, 4, bad, W))
        @test ex isa NumericalHealthError && ex.orbit == 0 && ex.step == 4 && ex.column == 1
    end
    N = 100
    z = nf_series(Float64, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), 2N))
    z[N + 5] = complex(NaN, 0.0)
    res = frequency_map(z, N)
    @test res.status[1] === :failed && res.failure.plane == 1 && res.failure.window == 2 && res.failure.sample == N + 5 && isnan(res.delta)
    e = nf_health_error(() -> check_frequency_map(res; orbit = 7))
    @test e isa NumericalHealthError && e.what == "series" && e.orbit == 7 && e.step == N + 5 && e.column == 1 && e.backend === :cpu
    zi = copy(z)
    zi[N + 5] = nf_series(Float64, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), 2N))[N + 5]
    zi[3] = complex(0.0, Inf)
    resi = frequency_map(zi, N)
    @test resi.failure.window == 1 && resi.failure.sample == 3 && nf_health_error(() -> check_frequency_map(resi; orbit = 2)).step == 3
    ok = frequency_map(nf_series(Float64, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), 2N)), N)
    @test check_frequency_map(ok; orbit = 1) === ok
    # Stream: a non-finite orbit among five, with the reducer in the chain, is raised by the host with the orbit id of the first bad orbit (the Tangent state check fires before the reducer sees the point; the
    # reducer's own record above is the second line of defence), the chunk stores nothing (propagate throws, no result).
    sysr = nf_rot_map(Float64, 0.7)
    xr = [SVector(0.1 * j, 0.2) for j in 1:5]
    xr[3] = SVector(NaN, 0.2)
    rr = WeightedPhaseAdvance(; window = 16, observable = PlanarObservable(1, 2))
    er = nf_health_error(() -> propagate(sysr, xr, 33; k = 1, tau = 1, seed = nf_seed(11, 2), reducers = (rr,), chunk = 2))
    @test er isa NumericalHealthError && er.orbit == 3 && er.backend === :cpu
end

@testset "AC11b degenerate series, fixed point" begin
    _ = (frequency_map, WeightedPhaseAdvance, LiftedAngleObservable)
    # An orbit started at the fixed point (0, 0) of the standard map: exp(i q) is the constant 1. With mean removal on the status is "degenerate" and the frequency undefined, recorded, not thrown; with mean
    # removal off the frequency is 0 (within 16 eps) with the "near zero" flag; E2 returns 0 exactly (the lifted increments are 0).
    N = 64
    for T in NF_TYPES
        epsT = nf_eps(T)
        sys = WrappedStandardMap(nf_K(T, :half))
        x0 = torus_state(SVector{2, T}(0, 0))
        st = nf_states(sys, x0, 2N + 1)
        @test all(s -> s[1] == 0 && s[2] == 0, st)
        z = nf_obs(nf_obs_cis, st)
        r = frequency_map(z, N; mean_removal = true)
        @test r.status[1] === :degenerate && r.w1[1].status === :degenerate && r.w2[1].status === :degenerate && isnan(r.w1[1].nu) && isnan(r.delta) && r.failure === nothing
        r0 = frequency_map(z, N)
        @test r0.status[1] === :ok && abs(r0.w1[1].nu) <= 16 * epsT && abs(r0.w2[1].nu) <= 16 * epsT && r0.w1[1].near_zero && r0.w2[1].near_zero
        red = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        _, _, fin = nf_run((red,), sys, x0, nf_frame(nf_seed(11, 1), 1, T, 2, 1), 2N + 1, 1)
        f = fin[1]
        @test f.nu1[1] == 0 && f.nu2[1] == 0 && f.delta[1] == 0 && f.below_floor[1]
    end
end

@testset "AC11c argument errors" begin
    _ = (window_frequency, frequency_map, weighted_phase_frequency, WeightedPhaseAdvance, PlanarObservable, check_run, window_weights)
    # N below the minimum window (8 (p + 1) samples for E1, 16 for E2), a series shorter than 2 N (E1) or 2 N + 1 (E2), a window length not matching the declared horizon, a non-positive pad, an order outside 0 to 5,
    # an unknown position convention, a mirror ratio or floor multiplier out of range: ArgumentError naming the parameter, before any work.
    named(f, word) = (e = nf_error(f); e isa ArgumentError && occursin(word, e.msg))
    z(n) = nf_series(Float64, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), n))
    @test named(() -> window_frequency(z(15)), "window")
    @test named(() -> window_frequency(z(23); order = 2), "window")
    @test nf_error(() -> window_frequency(z(16))) === nothing
    @test nf_error(() -> window_frequency(z(24); order = 2)) === nothing && nf_error(() -> window_frequency(z(8); order = 0)) === nothing
    @test named(() -> window_frequency(z(7); order = 0), "window")
    @test named(() -> frequency_map(z(200), 8), "window")
    @test named(() -> frequency_map(z(99), 50), "series") && nf_error(() -> frequency_map(z(100), 50)) === nothing && nf_error(() -> frequency_map(z(130), 50)) === nothing
    for pad in (0, -2)
        @test named(() -> window_frequency(z(64); pad = pad), "pad") && named(() -> frequency_map(z(128), 64; pad = pad), "pad")
    end
    for o in (-1, 6)
        @test named(() -> window_frequency(z(64); order = o), "order") && named(() -> window_weights(Float64, 64; order = o), "order")
    end
    @test named(() -> window_frequency(z(64); positions = :middle), "positions") && named(() -> window_weights(Float64, 64; positions = :middle), "positions")
    @test named(() -> frequency_map(z(128), 64; mirror_ratio = 0.0), "mirror_ratio") && named(() -> frequency_map(z(128), 64; mirror_ratio = 1.5), "mirror_ratio")
    @test named(() -> frequency_map(z(128), 64; floor_multiplier = 0), "floor_multiplier")
    # E2
    inc = fill(0.3, 100)
    @test named(() -> weighted_phase_frequency(inc, 8), "window") && named(() -> weighted_phase_frequency(inc, 64), "increments") && nf_error(() -> weighted_phase_frequency(inc, 50)) === nothing
    @test named(() -> WeightedPhaseAdvance(; window = 8, observable = PlanarObservable(1, 2)), "window") && named(() -> WeightedPhaseAdvance(; window = (100, 8), observable = PlanarObservable(1, 2)), "window")
    @test nf_error(() -> WeightedPhaseAdvance(; window = 100)) isa Exception && nf_error(() -> WeightedPhaseAdvance(; observable = PlanarObservable(1, 2))) isa Exception
    @test named(() -> WeightedPhaseAdvance(; window = 100, observable = PlanarObservable(1, 2), birkhoff_order = 0), "birkhoff_order")
    @test named(() -> check_run(WeightedPhaseAdvance(; window = 100, observable = PlanarObservable(1, 2)), 1, 150), "window")
end

@testset "AC11d tie of two equal lines" begin
    _ = (frequency_map,)
    # Two exactly equal lines (a real cosine, no mean removal): the coarse tie is flagged (values equal within the tie tolerance, 64 eps of Float64; the spec's 8 eps is below the rounding of a Float64 transform), the lowest index wins (the negative frequency), and the "mirror" flag is set.
    N = 1000
    nu = nf_dec("0.2")
    for T in NF_TYPES
        z = nf_series(T, nf_lines_big([(nu, BigFloat(1) / 2, BigFloat(0)), (-nu, BigFloat(1) / 2, BigFloat(0))], 2N))
        r = frequency_map(z, N)
        for w in (r.w1[1], r.w2[1])
            @test w.tie && w.mirror && w.nu < 0 && abs(abs(w.nu) - 0.2) < 1e-8
        end
    end
end

@testset "AC12 execution invariance" begin
    _ = (frequency_map, WeightedPhaseAdvance, LiftedAngleObservable, MEGNO)
    # E1 through a series (built serially) and E2 through propagate: outputs bitwise identical for 1 versus several threads and for chunk sizes 64 and 1000 on the CPU backend; E2 through propagate equals the
    # serial per-orbit body. The transform library's own thread count (1 when orbit threads exceed 1) is an implementation obligation checked in the stage 04 child-process run. D38: 128 orbits, N = 250
    # (spec: 4096 orbits, N = 1000, child processes at stage 04); K_std in {0.5, 2 pi} with uniform starts.
    norb = NF_FULL ? 4096 : 128
    N = NF_FULL ? 1000 : 250
    for (label, sub) in ((:half, 5), (:two_pi, 1))
        sys = WrappedStandardMap(nf_K(Float64, label))
        starts = nf_starts(StableRNG(nf_seed(12, sub)), norb)
        x0 = [torus_state(nf_cast(Float64, s)) for s in starts]
        # E1
        sers = [nf_obs(nf_obs_cis, nf_states(sys, x, 2N)) for x in x0]
        key(r) = (nf_wkey(r.w1[1]), nf_wkey(r.w2[1]), r.delta, r.D, r.below_floor, r.status)
        serial = [key(frequency_map(z, N)) for z in sers]
        thr = Vector{Any}(undef, norb)
        Threads.@threads :dynamic for j in 1:norb
            thr[j] = key(frequency_map(sers[j], N))
        end
        @test isequal(serial, thr)
        for chunk in (64, 1000)
            out = Vector{Any}(undef, norb)
            starts_ = collect(1:chunk:norb)
            Threads.@threads :dynamic for ci in reverse(eachindex(starts_))
                for j in starts_[ci]:min(starts_[ci] + chunk - 1, norb)
                    out[j] = key(frequency_map(sers[j], N))
                end
            end
            @test isequal(serial, out)
        end
        # E2
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        reducers = (r, MEGNO())
        seed = nf_seed(12, 10 + sub)
        W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:norb]
        ser = [propagate_orbit(sys, x0[j], W0[j], 2N + 1, 1, reducers; orbit = j) for j in 1:norb]
        fs = [map(reducer_final, reducers, s[3]) for s in ser]
        for chunk in (64, nothing)
            res = propagate(sys, x0, 2N + 1; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, [s[1] for s in ser]) && isequal(res.finals, fs)
        end
    end
end

@testset "AC13 precision ladder" begin
    _ = (frequency_map, WeightedPhaseAdvance, LiftedAngleObservable)
    # K_std = 0.5 regular class and K_std = 10, N = 1000, Float64 and Double64 each against its own BigFloat reference at its own realized samples: every E1 frequency (observable exp(i q), window means
    # removed) within the bound of criterion 9, every E2 frequency within B_E2 of criterion 10(d). The difference of the Float64 and Double64 frequencies on the same start is reported (R7), not gated. D38: 5
    # starts per class (FULL 20).
    N = 1000
    nstart = NF_FULL ? 20 : 5
    for T in NF_TYPES, (cl, sub) in ((:half, 0), (:ten, 3))
        epsT = nf_eps(T)
        c = nf_case(cl, T, nstart, sub)
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        bad = String[]
        counts = Dict(:ok => 0, :tie => 0, :nonunique => 0)
        for j in 1:nstart
            st = nf_states(c.sys, c.x0[j], 2N + 1)
            zb = Complex{BigFloat}[nf_big(x) for x in nf_obs(nf_obs_cis, st[1:(2N)])]
            z = nf_series(T, nf_demean_big(zb, N))
            res = frequency_map(z, N)
            for (k, (w, zw)) in enumerate(((res.w1[1], z[1:N]), (res.w2[1], z[(N + 1):(2N)])))
                s, _ = nf_compare_window(T, zw, w, 1, :half, 8)
                haskey(counts, s) ? (counts[s] += 1) : push!(bad, "T=$T $cl start $j E1 window $k $s")
            end
            _, _, fin = nf_run((nf_Trace(), r, nf_Norm()), c.sys, c.x0[j], nf_frame(nf_seed(13, 1), c.ids[j], T), 2N + 1, 1)
            incr = nf_lifted_incr_big(nf_vec(fin[1].xs))
            ref = nf_E2_ref(incr, N; bits = 512)
            B = nf_BE2(T, incr)
            for (k, v) in enumerate((fin[2].nu1[1], fin[2].nu2[1]))
                Float64(abs(nf_big(v) - ref[k])) <= B || push!(bad, "T=$T $cl start $j E2 window $k")
            end
        end
        @test nf_report(bad, "AC13 $T $cl")
        @info "AC13 $T $cl: E1 windows by status (listed, not failed)" counts
    end
    # R7 (reported): Float64 against Double64 first-window frequencies on the same regular starts (the orbits agree to ~1e-12 for the regular class; chaotic orbits are not compared)
    c64, c128 = nf_case(:half, Float64, nstart, 0), nf_case(:half, Double64, nstart, 0)
    d = Float64[]
    for j in 1:nstart
        z64 = nf_obs(nf_obs_cis, nf_states(c64.sys, c64.x0[j], N))
        z128 = nf_obs(nf_obs_cis, nf_states(c128.sys, c128.x0[j], N))
        push!(d, Float64(abs(nf_big(window_frequency(z64).nu) - nf_big(window_frequency(z128).nu))))
    end
    @info "R7 (reported, no threshold): |nu_Float64 - nu_Double64|, first window, regular class" maximum_difference = maximum(d)
    @test all(isfinite, d)
end
