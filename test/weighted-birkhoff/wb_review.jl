# Tests of the stage 04 independent review of weighted-birkhoff (findings M1, M2, L1, L2, L3, L5; the module-build stage 04 report lists them with their dispositions). New file: the other files
# of this folder are byte-identical to stage 02. Each testset fails on the code before the fixes and passes now.

struct wb_Throws end                      # an observable that throws (a host-side check item; never used in a stream)
(o::wb_Throws)(x) = error("boom")
struct wb_BadLift                         # the lifted q of the wrapped standard map, NaN once it exceeds c (a monotone orbit: the bad sample is the first beyond c)
    c::Float64
end
(o::wb_BadLift)(x) = (v = x[1] + 2 * Chaoslab.Systems._pi(x[1]) * x[3]; v > o.c ? oftype(x[1], NaN) : v)

@testset "REV M1 birkhoff_order above 4" begin
    _ = (WeightedBirkhoff, weighted_birkhoff)
    # The centre weight exp(-4^p) underflows from p = 5 in Float64 and in Double64, so B = 0 and every average was NaN without a flag: now an ArgumentError naming birkhoff_order, in the
    # constructor and in the function; p = 4 still works in both types (B > 0, finite averages). A weight vector that sums to 0 is an ArgumentError naming weights.
    named(f, word) = (e = wb_error(f); e isa ArgumentError && occursin(word, e.msg))
    obs = (CosObservable(1),)
    for p in (5, 6, 10)
        @test named(() -> WeightedBirkhoff(; window = 100, observables = obs, birkhoff_order = p), "birkhoff_order")
        for T in WB_TYPES
            @test named(() -> weighted_birkhoff(ones(T, 200), 100; order = p), "birkhoff_order")
        end
    end
    @test wb_error(() -> WeightedBirkhoff(; window = 100, observables = obs, birkhoff_order = 4)) === nothing
    @test named(() -> WeightedBirkhoff(; window = 100, observables = obs, birkhoff_order = 0), "birkhoff_order")
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    for T in WB_TYPES, p in (1, 2, 3, 4)
        cs, _ = wb_rot_big(wb_dec("0.3"), alpha, 1, 2000)
        r = weighted_birkhoff(wb_round_vec(T, cs), 1000; order = p)
        @test r.B[1] > 0 && isfinite(r.avg1[1]) && isfinite(r.avg2[1]) && abs(r.avg1[1]) < 1e-3
    end
    for T in WB_TYPES
        @test named(() -> weighted_birkhoff(ones(T, 64), 32; weights = zeros(T, 32)), "weights")
        @test named(() -> weighted_birkhoff(ones(T, 64), 32; weights = -ones(T, 32)), "weights")
    end
    @test wb_error(() -> weighted_birkhoff(ones(64), 32; weights = ones(32))) === nothing
end

@testset "REV M2 nothing is evaluated after every horizon is frozen" begin
    _ = (WeightedBirkhoff, reducer_step, reducer_failed, reducer_final)
    # Hand protocol: after 2 max(window) good handoffs, handoffs beyond (the 2N + 1th of a chain with E2) with a non-finite observable do not set the failure code, change no field of the result
    # but the counters, and the horizons stay as they were (I4). Before the fix the extra handoff set the code and the host raised with the orbit id.
    rf = WeightedBirkhoff(; window = (16, 32), observables = (CosObservable(1), wb_BadAt(0.5)))
    W = SMatrix{2, 1, Float64}(1.0, 0.0)
    good, bad = SVector(0.3, 0.4), SVector(0.9, 0.4)
    st = reducer_init(rf, Float64, Val(2), Val(1))
    for n in 1:64
        st, _ = reducer_step(rf, st, n, n % 3 == 0 ? SVector(0.1 * n % 0.45, 0.4) : good, W)
    end
    @test !reducer_failed(rf, st) && reducer_final(rf, st).complete
    f0 = reducer_final(rf, st)
    st65, W65 = reducer_step(rf, st, 65, bad, W)
    st66, _ = reducer_step(rf, st65, 66, bad, W)
    @test !reducer_failed(rf, st65) && !reducer_failed(rf, st66) && W65 === W
    for st2 in (st65, st66)
        f = reducer_final(rf, st2)
        @test f.complete && all(fld -> isequal(getfield(f, fld), getfield(f0, fld)), (:avg1, :avg2, :S, :B, :delta, :absdig, :reldig, :dig, :below_floor, :zero))
    end
    @test reducer_final(rf, st65).handoffs == 65 && reducer_final(rf, st66).handoffs == 66 && reducer_final(rf, st66).n == 66
    # Before completion the failure is still recorded (handoff 64 is the last needed one).
    sb = reducer_init(rf, Float64, Val(2), Val(1))
    for n in 1:63
        sb, _ = reducer_step(rf, sb, n, good, W)
    end
    sb64, _ = reducer_step(rf, sb, 64, bad, W)
    @test reducer_failed(rf, sb64)
    # Stream: the K_std = 0 map (the lifted q grows by p = 1 per step) whose 33rd sample alone exceeds the bound of the observable, window 16, chain (reader, E2) of 2N + 1 = 33 handoffs: no error,
    # reader complete and equal to its 32-handoff run.
    sys0 = WrappedStandardMap(wb_K(Float64, :zero))
    x = torus_state(SVector(0.3, 1.0))
    c = wb_BadLift(Inf)(wb_states(sys0, x, 32)[32])
    rr = WeightedBirkhoff(; window = 16, observables = (wb_BadLift(c),))
    e2 = WeightedPhaseAdvance(; window = 16, observable = LiftedAngleObservable(1, 3))
    @test isnan(wb_BadLift(c)(wb_states(sys0, x, 33)[33])) && isfinite(wb_BadLift(c)(wb_states(sys0, x, 33)[32]))
    res33 = propagate(sys0, [x], 33; k = 1, tau = 1, seed = wb_seed(30, 1), reducers = (rr, e2))
    res32 = propagate(sys0, [x], 32; k = 1, tau = 1, seed = wb_seed(30, 1), reducers = (rr,))
    a, b = res33.finals[1][1], res32.finals[1][1]
    @test a.complete && b.complete && a.handoffs == 33 && b.handoffs == 32
    @test all(fld -> isequal(getfield(a, fld), getfield(b, fld)), (:avg1, :avg2, :S, :B, :delta, :absdig, :reldig, :dig, :below_floor, :zero))
end

@testset "REV L1 at most 31 observables" begin
    _ = (WeightedBirkhoff, check_run, reducer_bytes)
    # Beyond 31 the mapped tuple leaves the unrolled path (allocation, possible InvalidIRError on the device): the constructor refuses with an ArgumentError naming observables; 31 runs without
    # allocation.
    named(f, word) = (e = wb_error(f); e isa ArgumentError && occursin(word, e.msg))
    mk(K) = ntuple(i -> CosObservable(1, i), K)
    @test named(() -> WeightedBirkhoff(; window = 16, observables = mk(32)), "observables") && named(() -> WeightedBirkhoff(; window = 16, observables = mk(100)), "observables")
    r = WeightedBirkhoff(; window = 16, observables = mk(31))
    @test reducer_bytes(r, Float64, Val(2), Val(1)) == 16 + 8 * (2 + 5 * 31)
    sys = WrappedStandardMap(wb_K(Float64, :two_pi))
    x0 = torus_state(SVector(1.0, 0.5))
    W0 = wb_frame(wb_seed(31, 1), 1, Float64, 2, 2)
    chain = (r, FLI())
    propagate_orbit(sys, x0, W0, 32, 1, chain)
    @test (@allocated propagate_orbit(sys, x0, W0, 32, 1, chain)) == 0
end

@testset "REV L2 the floor of denormal and huge scales" begin
    _ = (birkhoff_index, weighted_birkhoff)
    # A denormal S made phi = 32 eps S underflow to 0 and the digits Inf without a flag: the floor is at least floatmin(T). Averages near floatmax overflowed m = (|a1| + |a2|) / 2 and gave
    # reldig = Inf: m is now formed from the halves. Delta overflowing (opposite averages near floatmax) gives undefined digits (NaN), unflagged, never an infinite number.
    for T in WB_TYPES
        tiny = T(1e-310)
        ix = birkhoff_index(tiny, tiny, tiny)
        @test ix.below_floor && ix.delta == 0 && isfinite(ix.absdig) && !ix.zero && ix.absdig <= -log10(floatmin(T)) + 1e-9
        ix2 = birkhoff_index(tiny, T(0.5) * tiny, tiny)                      # Delta = 5e-311, below the floor floatmin: flagged, finite digits
        @test ix2.below_floor && isfinite(ix2.absdig) && ix2.delta == tiny - T(0.5) * tiny
        big_ = T(1.7e308)
        ix3 = birkhoff_index(big_, big_, big_)
        @test ix3.below_floor && isfinite(ix3.absdig) && isfinite(ix3.reldig) && isfinite(ix3.dig) && ix3.reldig < -log10(32 * wb_eps(T)) + 1
        ix4 = birkhoff_index(big_, -big_, big_)
        @test isnan(ix4.absdig) && isnan(ix4.reldig) && isnan(ix4.dig) && !ix4.below_floor && !ix4.zero && !isfinite(ix4.delta)
        # ordinary scales are unchanged: S = 1, Delta = 1e-8, averages near 0.4 and 0.4 + 1e-8
        ix5 = birkhoff_index(T(0.4) + T(1e-8), T(0.4), T(1))
        @test abs(ix5.absdig - 8) < 1e-6 && !ix5.below_floor && abs(ix5.reldig - (8 + log10(0.4))) < 1e-6
    end
    # through the function: a series of denormal magnitude gives finite digits and the flag
    f = 1e-310 .* [cos(0.37 * j) for j in 0:63]
    r = weighted_birkhoff(f, 32)
    ix = birkhoff_index(r.avg1[1], r.avg2[1], r.S[1])
    @test isfinite(ix.absdig) && isfinite(ix.dig) && ix.below_floor
end

@testset "REV L3 the series of the pure function is checked" begin
    _ = (weighted_birkhoff,)
    # A NaN (or Inf) among the used samples gave NaN averages and NaN digits without any flag: an ArgumentError naming series. Samples beyond 2 max(window) are ignored, so a bad value there
    # is not an error.
    named(f, word) = (e = wb_error(f); e isa ArgumentError && occursin(word, e.msg))
    for T in (Float64, Double64, BigFloat), bad in (NaN, Inf, -Inf)
        f = T.(cos.(0.3 .* (0:63)))
        g = copy(f)
        g[20] = T(bad)
        @test named(() -> weighted_birkhoff(g, 32), "series") && named(() -> weighted_birkhoff(g, (16, 32)), "series")
        @test named(() -> weighted_birkhoff([f f.+1 g], 32), "series")
        g2 = vcat(f, T[bad])                                              # beyond the 2N = 64 used samples
        @test wb_error(() -> weighted_birkhoff(g2, 32)) === nothing
        g3 = copy(f); g3[40] = T(bad)                                     # inside horizon 32's second window but outside horizon 16's 2N = 32 samples
        @test named(() -> weighted_birkhoff(g3, (16, 32)), "series")
        @test wb_error(() -> weighted_birkhoff(g3, 16)) === nothing
    end
end

@testset "REV L5 check_run: tau and the names of the messages" begin
    _ = (WeightedBirkhoff, check_run)
    # tau = 0 raised a DivideError, a negative tau a wrong message: ArgumentError naming tau. The observable messages name observables (D84 item 2). Without x0 the observables are not evaluated,
    # so a throwing observable is not caught by the check (documented).
    named(f, word) = (e = wb_error(f); e isa ArgumentError && occursin(word, e.msg))
    x0 = SVector(0.3, 0.4, 0.0, 0.0)
    r = WeightedBirkhoff(; window = 100, observables = (CosObservable(1),))
    @test named(() -> check_run(r, 0, 200), "tau") && named(() -> check_run(r, -1, 200), "tau") && named(() -> check_run(r, -4, -800), "tau")
    @test check_run(r, 1, 200) === nothing
    for bad_obs in (wb_BadAt(0.1), wb_Const(Inf), wb_Const(NaN), wb_Tuple(), wb_Throws())
        rb = WeightedBirkhoff(; window = 100, observables = (CosObservable(1), bad_obs))
        e = wb_error(() -> check_run(rb, 1, 200; x0 = x0))
        @test e isa ArgumentError && occursin("observables", e.msg) && occursin("[2]", e.msg)
        @test check_run(rb, 1, 200) === nothing                        # no x0: no evaluation, a throwing observable is not seen
    end
    @test named(() -> WeightedBirkhoff(; window = 100, observables = ()), "observables")
end
