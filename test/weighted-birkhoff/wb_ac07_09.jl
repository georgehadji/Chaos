# AC7 to AC9: consistency with E2, horizons / observables / state size (I4, P7), checks and edges.

@testset "AC7 consistency with E2" begin
    _ = (WeightedBirkhoff, check_run, weighted_birkhoff)
    # F3: the average of f(x) = p + K_std sin q (lifted p, a caller function) over each window equals the E2 sums A/B of the lifted increments of the same orbit, within c_1 eps S. K_std in {0.5 regular
    # class, 0.9}, tau = 1, N in {500, 1000, 4000}: E2 (WeightedPhaseAdvance, LiftedAngleObservable(1, 3)) has 2N + 1 handoffs, this reader 2N (it ignores the extra one). A/B = tp nu with tp = 2 pi_T.
    # D38: 4 starts, N in {500, 1000} (spec: 20 starts, N up to 4000).
    nstart = WB_FULL ? 20 : 4
    Ns = WB_FULL ? (500, 1000, 4000) : (500, 1000)
    for T in WB_TYPES
        epsT = wb_eps(T)
        bad = String[]
        worst = 0.0
        for (cl, sub) in ((:half, 0), (:nine, 2)), N in Ns
            c = wb_case(cl, T, nstart, sub)
            r = WeightedBirkhoff(; window = N, observables = (wb_StdKick(wb_Kf(cl)),))
            e2 = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
            check_run(r, 1, 2N)
            check_run(e2, 1, 2N + 1)
            tp = wb_tp(T)
            for j in 1:nstart
                _, _, fin = wb_run((r, e2, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(7, 1), c.ids[j], T), 2N + 1, 1)
                f, g = fin[1], fin[2]
                S = Float64(f.S[1])
                for (a, nu) in ((f.avg1[1], g.nu1[1]), (f.avg2[1], g.nu2[1]))
                    e = Float64(abs(wb_big(a) - wb_big(nu) * tp))
                    worst = max(worst, e / (epsT * S))
                    e <= wb_cp(1) * epsT * S || push!(bad, "T=$T $cl N=$N start $j: $(e / (epsT * S)) eps S")
                end
                (f.complete && g.complete && f.handoffs == 2N + 1) || push!(bad, "T=$T $cl N=$N start $j: runs not complete")
            end
        end
        @test wb_report(bad, "AC7 $T")
        @info "AC7 $T: worst |A/B - tp nu_E2| in eps S (spec dry run 0.5; gate 4)" worst
    end
end

@testset "AC8 horizons, observables, state" begin
    _ = (WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable, reducer_bytes, reducer_workspace)
    # H = 3 (N = 64, 256, 1000), K = 4: results bitwise equal to the 12 single runs and, with the observables and the horizons permuted, to the same entries permuted; a completed horizon never changes
    # (I4); the reducer bytes are at most 16 + 8 H (2 + 5 K) (Float64) and 16 + 16 H (2 + 5 K) (Double64) for (K, H) in {(1, 1), (4, 1), (4, 3), (4, 11)} (P7); the state is isbits, the CPU path
    # allocates nothing after warm-up (spec Compute).
    Ns = (64, 256, 1000)
    flds = (:avg1, :avg2, :S, :delta, :absdig, :reldig, :dig, :below_floor, :zero)
    perm, hperm = (4, 2, 1, 3), (3, 1, 2)
    for T in WB_TYPES
        bad = String[]
        obs = wb_std4()
        for (cl, sub) in ((:two_pi, 1), (:half, 0))
            c = wb_case(cl, T, 1, sub)
            sys, x0, W0 = c.sys, c.x0[1], wb_frame(wb_seed(8, 1), 1, T)
            full = WeightedBirkhoff(; window = Ns, observables = obs)
            shuf = WeightedBirkhoff(; window = map(i -> Ns[i], hperm), observables = map(i -> obs[i], perm))
            _, _, fin = wb_run((full, shuf, wb_Norm()), sys, x0, W0, 2000, 1)
            f, g = fin[1], fin[2]
            for fld in flds, k in 1:4, h in 1:3
                isequal(getfield(g, fld)[k, h], getfield(f, fld)[perm[k], hperm[h]]) || push!(bad, "T=$T $cl permuted: $fld obs $k horizon $h")
            end
            for (h, N) in enumerate(Ns), k in 1:4
                one = WeightedBirkhoff(; window = N, observables = (obs[k],))
                _, _, fo = wb_run((one, wb_Norm()), sys, x0, W0, 2000, 1)
                for fld in flds
                    isequal(getfield(fo[1], fld)[1, 1], getfield(f, fld)[k, h]) || push!(bad, "T=$T $cl single run: $fld obs $k horizon $h")
                end
                isequal(fo[1].B[1], f.B[h]) || push!(bad, "T=$T $cl single run: B horizon $h")
            end
            # a completed horizon never changes: the horizon 64 is complete after 128 handoffs and equals its value after 2000
            short = WeightedBirkhoff(; window = 64, observables = obs)
            _, _, fs = wb_run((short, wb_Norm()), sys, x0, W0, 128, 1)
            _, _, fl = wb_run((short, wb_Norm()), sys, x0, W0, 2000, 1)
            (fs[1].complete && fl[1].complete && all(fld -> isequal(getfield(fs[1], fld), getfield(fl[1], fld)), flds) && isequal(fs[1].B, fl[1].B)) || push!(bad, "T=$T $cl completed horizon changed")
            all(fld -> isequal(getfield(fs[1], fld)[:, 1], getfield(f, fld)[:, 1]), flds) || push!(bad, "T=$T $cl horizon 64 alone differs from the horizon 64 of three")
        end
        @test wb_report(bad, "AC8 $T")
        for (K, H) in ((1, 1), (4, 1), (4, 3), (4, 11))
            r = WeightedBirkhoff(; window = ntuple(i -> 16 * 2^(i - 1), H), observables = ntuple(i -> CosObservable(1, i), K))
            b = reducer_bytes(r, T, Val(2), Val(1))
            @test b <= 16 + sizeof(T) * H * (2 + 5K) && reducer_workspace(r, T, Val(2), Val(1)) == 0
            @test isbitstype(typeof(r)) && isbitstype(typeof(reducer_init(r, T, Val(2), Val(1))))
            @info "AC8 reducer bytes (reported)" T K H bytes = b bound = 16 + sizeof(T) * H * (2 + 5K)
        end
    end
    sys = WrappedStandardMap(wb_K(Float64, :two_pi))
    x0 = torus_state(SVector(1.0, 0.5))
    W0 = wb_frame(wb_seed(8, 2), 1, Float64, 2, 2)
    r = WeightedBirkhoff(; window = (50, 100), observables = wb_std4())
    chain = (r, FLI())
    propagate_orbit(sys, x0, W0, 200, 1, chain)
    @test (@allocated propagate_orbit(sys, x0, W0, 200, 1, chain)) == 0
end

@testset "AC9 checks and edges" begin
    _ = (WeightedBirkhoff, weighted_birkhoff, check_run, CosObservable, LiftedObservable, reducer_step, reducer_failed, reducer_error)
    # Argument errors name the parameter, before any work: a window below 16 (window), an order below 1 (birkhoff_order), K = 0 (observables), a run length other than 2 max(window) handoffs (window),
    # an observable non-finite at the initial point or not a real number (observable). The constructor and check_run are the host check; the errors may come from either.
    named(f, word) = (e = wb_error(f); e isa ArgumentError && occursin(word, e.msg))
    obs = (CosObservable(1),)
    x0 = SVector(0.3, 0.4, 0.0, 0.0)
    @test named(() -> WeightedBirkhoff(; window = 15, observables = obs), "window") && named(() -> WeightedBirkhoff(; window = (100, 8), observables = obs), "window")
    @test wb_error(() -> WeightedBirkhoff(; window = 16, observables = obs)) === nothing && wb_error(() -> WeightedBirkhoff(; window = (16, 17, 5000), observables = obs)) === nothing
    @test named(() -> WeightedBirkhoff(; window = 100, observables = obs, birkhoff_order = 0), "birkhoff_order") && named(() -> weighted_birkhoff(zeros(40), 16; order = 0), "birkhoff_order")
    @test named(() -> WeightedBirkhoff(; window = 100, observables = ()), "observables")
    @test named(() -> weighted_birkhoff(zeros(40), 15), "window") && named(() -> weighted_birkhoff(zeros(31), 16), "series") && wb_error(() -> weighted_birkhoff(zeros(32), 16)) === nothing
    @test wb_error(() -> WeightedBirkhoff(; window = 100)) isa Exception && wb_error(() -> WeightedBirkhoff(; observables = obs)) isa Exception
    r = WeightedBirkhoff(; window = 100, observables = obs)
    @test check_run(r, 1, 200) === nothing && check_run(r, 4, 800) === nothing && check_run(r, 1, 200; x0 = x0) === nothing
    for (tau, n) in ((1, 199), (1, 201), (1, 202), (4, 804), (4, 796), (3, 603), (3, 597))
        @test named(() -> check_run(r, tau, n), "window")
    end
    r2 = WeightedBirkhoff(; window = (50, 100), observables = obs)
    @test check_run(r2, 1, 200) === nothing && named(() -> check_run(r2, 1, 100), "window")
    for bad_obs in (wb_BadAt(0.1), wb_Const(Inf), wb_Const(NaN), wb_Tuple())
        rb = WeightedBirkhoff(; window = 100, observables = (CosObservable(1), bad_obs))
        @test named(() -> check_run(rb, 1, 200; x0 = x0), "observable")
    end
    # A non-finite observable at a handoff: recorded in the state with the observable index, never thrown in reducer_step; the state stays failed and unchanged; W is untouched; the host error carries
    # the orbit id, the step and the observable index; reducer_update (direct use) throws at once with orbit 0.
    rf = WeightedBirkhoff(; window = 16, observables = (CosObservable(1), wb_BadAt(0.5)))
    W = SMatrix{2, 1, Float64}(1.0, 0.0)
    good, bad = SVector(0.3, 0.4), SVector(0.9, 0.4)
    st = reducer_init(rf, Float64, Val(2), Val(1))
    for n in 1:3
        st, W2 = reducer_step(rf, st, n, good, W)
        @test W2 === W
    end
    st4, W4 = reducer_step(rf, st, 4, bad, W)
    @test reducer_failed(rf, st4) && !reducer_failed(rf, st) && W4 === W
    e = reducer_error(rf, st4, 9, :cpu)
    @test e isa NumericalHealthError && e.what == "observable" && e.orbit == 9 && e.step == 4 && e.backend === :cpu && e.column == 2
    st5, _ = reducer_step(rf, st4, 5, good, W)
    @test reducer_failed(rf, st5) && isequal(st5, st4)
    ex = wb_health_error(() -> reducer_update(rf, st, 4, bad, W))
    @test ex isa NumericalHealthError && ex.orbit == 0 && ex.step == 4 && ex.column == 2
    # Stream: five orbits of a rotation block, one of them (the first whose x exceeds 0.3) makes the observable non-finite; the host raises with that orbit id and the chunk stores nothing.
    sysr = wb_rot_map(Float64, 0.7)
    xr = [SVector(0.1 * j, 0.2) for j in 1:5]
    first_bad = findfirst(x -> any(s -> s[1] > 0.3, wb_states(sysr, x, 32)), xr)
    @test first_bad !== nothing && first_bad > 1
    rr = WeightedBirkhoff(; window = 16, observables = (wb_BadAt(0.3),))
    er = wb_health_error(() -> propagate(sysr, xr, 32; k = 1, tau = 1, seed = wb_seed(9, 2), reducers = (rr,), chunk = 2))
    @test er isa NumericalHealthError && er.what == "observable" && er.orbit == first_bad && er.column == 1 && er.backend === :cpu
    # A wrapped angle (a sawtooth: discontinuous) and a lifted angle (not periodic) are accepted and give finite output (R3): complete runs, finite averages, Delta and digits.
    for T in WB_TYPES
        c = wb_case(:half, T, 1, 0)
        ra = WeightedBirkhoff(; window = 100, observables = (wb_Coord(1), LiftedObservable(1, 3), wb_Coord(2)))
        _, _, fin = wb_run((ra, wb_Norm()), c.sys, c.x0[1], wb_frame(wb_seed(9, 3), 1, T), 200, 1)
        f = fin[1]
        @test f.complete && all(isfinite, f.avg1) && all(isfinite, f.avg2) && all(isfinite, f.delta) && all(isfinite, f.absdig) && all(isfinite, f.dig) && all(isfinite, f.S)
    end
end
