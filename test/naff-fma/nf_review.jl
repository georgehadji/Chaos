# Stage 04 review follow-up (module-build naff-fma, after the independent review): one testset per fixed item, each reproducing the reviewer's failing scenario.
# Added after the stage 02 tests, which are untouched. Items: MEDIUM 1 (increments near +-pi), 2 (tie at large N), 3 (step counts beyond Int32), LOW (magnitudes, zero weight,
# bracket, pad, observable indices). MEDIUM 4 (provenance key names) is not fixed: committed tests read the old names (impl notes).

@testset "Review M1 increments near pi are recorded" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable, frequency_map)
    # z_j = cis(2 pi 0.495 j) (1 + 0.05 cis(2 pi 0.1234 j)), N = 1000: E1 returns 0.495, the principal argument of E2 is biased (0.3205, status ok); the reducer must say so.
    N = 1000
    zs(nu0, n) = ComplexF64[cis(2 * pi * nu0 * j) * (1 + 0.05 * cis(2 * pi * 0.1234 * j)) for j in 0:(n - 1)]
    function runr(r, z, lifted = false)
        W = SMatrix{2, 1, Float64}(1, 0)
        st = reducer_init(r, Float64, Val(2), Val(1))
        for (n, x) in enumerate(z)
            st, _ = reducer_step(r, st, n, lifted ? SVector(real(x), 0.0, 0.0, 0.0) : SVector(real(x), imag(x)), W)
        end
        return reducer_final(r, st)
    end
    r = WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2))
    z = zs(0.495, 2N + 1)
    @test abs(frequency_map(z, N).w1[1].nu - 0.495) < 1e-6
    f = runr(r, z)
    @test f.complete && abs(f.nu1[1] - 0.495) > 0.1            # the estimator is biased: the data property the flag exists for
    @test f.near_pi > N                                        # most increments are near the branch cut
    g = runr(r, zs(0.1, 2N + 1))
    @test g.complete && g.near_pi == 0 && abs(g.nu1[1] - 0.1) < 1e-2
    # the lifted observable has exact increments: never flagged
    rl = WeightedPhaseAdvance(; window = 16, observable = LiftedAngleObservable(1, 3))
    @test runr(rl, ComplexF64[complex(0.1 * j, 0) for j in 1:33], true).near_pi == 0
    # increments that do not enter a window are not counted (the run continues past 2N + 1 handoffs)
    f2 = runr(r, zs(0.495, 2N + 11))
    @test f2.near_pi == f.near_pi
end

@testset "Review M2 tie of exactly symmetric spectra at large N" begin
    _ = (frequency_map,)
    # A real cosine has |phi(-nu)| = |phi(nu)| exactly; the direct sums of the coarse search must be accurate enough to flag the tie and take the lowest index (the negative frequency)
    # for N of 8000 and more (a plain running sum has rounding noise above the 64 eps tie tolerance there).
    for N in (8000, 20000, 100000)
        z = ComplexF64[complex(cos(2 * pi * 0.2 * j), 0.0) for j in 0:(2N - 1)]
        r = frequency_map(z, N)
        for w in (r.w1[1], r.w2[1])
            @test w.tie && w.mirror && w.nu < 0 && abs(abs(w.nu) - 0.2) < 1e-8
        end
    end
end

@testset "Review M3 step counts beyond Int32" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, reducer_error)
    # n is a map step as Int: no InexactError for n > 2^31 - 1 in a step, in a failure record, or in the final; the handoff counter saturates instead of wrapping or throwing.
    r = WeightedPhaseAdvance(; window = 16, observable = PlanarObservable(1, 2))
    W = SMatrix{2, 1, Float64}(1, 0)
    big = 2^31 + 5
    st = reducer_init(r, Float64, Val(2), Val(1))
    st, _ = reducer_step(r, st, big, SVector(0.3, 0.4), W)
    st, _ = reducer_step(r, st, big + 1, SVector(0.31, 0.41), W)
    f = reducer_final(r, st)
    @test f.n == big + 1 && f.handoffs == 2
    bad, _ = reducer_step(r, st, big + 2, SVector(NaN, 0.4), W)
    e = reducer_error(r, bad, 3, :cpu)
    @test reducer_failed(r, bad) && e.step == big + 2 && reducer_final(r, bad).n == big + 2
    stsat = Chaoslab.Indicators.WPAState{Float64, 1}(st.hor, st.prev1, st.prev2, big, typemax(Int32), Int16(0), Int8(0))
    st2, _ = reducer_step(r, stsat, big + 1, SVector(0.32, 0.42), W)
    @test reducer_final(r, st2).handoffs == typemax(Int32) && reducer_final(r, st2).complete
    # the size is unchanged by the Int64 n and the counter
    for (T, b) in ((Float64, 96), (Double64, 176))
        st1 = reducer_init(r, T, Val(2), Val(1))
        @test sizeof(st1) == b
        @test sizeof(reducer_init(WeightedPhaseAdvance(; window = (16, 32), observable = PlanarObservable(1, 2)), T, Val(2), Val(1))) == b + (T === Float64 ? 64 : 128)
    end
end

@testset "Review LOW magnitudes" begin
    _ = (window_frequency, phase_increments)
    # A window of samples of magnitude 1e300 or 1e-300 (the squares of the derivative sums overflowed or underflowed, and a transform of 1e307 samples was NaN): same frequency, amplitude scaled.
    N = 256
    z = ComplexF64[cis(2 * pi * 0.1234567891 * j) + 0.1 * cis(2 * pi * 0.3 * j) for j in 0:(N - 1)]
    w0 = window_frequency(z)
    for sc in (2.0^700, 2.0^-700, 1e300, 1e-300, 1e307)
        w = window_frequency(sc .* z)
        @test w.status === :ok && abs(w.nu - w0.nu) < 1e-12 && abs(w.amplitude / sc - w0.amplitude) < 1e-12 * w0.amplitude && w.bin == w0.bin
    end
    @test isequal(window_frequency(2.0^700 .* z).nu, w0.nu) && isequal(window_frequency(2.0^-700 .* z).nu, w0.nu)
    # the planar increment of samples of any magnitude (z' conj(z) overflowed or underflowed)
    zi = ComplexF64[cis(0.3 * j + 0.02 * j^2) * (1 + sin(j) / 2) for j in 0:40]
    inc0 = phase_increments(zi)
    for sc in (1e200, 1e-200, 1e300, 1e-300)
        @test maximum(abs.(phase_increments(sc .* zi) .- inc0)) <= 4 * eps()
    end
end

@testset "Review LOW zero weight, bracket, pad, indices" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable, weighted_phase_frequency, phase_increments, window_frequency, frequency_map)
    # A sample at the origin that only feeds the increment with the weight 0 (the first sample of window 1) does not spoil the window: the NaN increment times the weight 0 is skipped.
    N = 16
    z = ComplexF64[cis(0.3 * s + 0.02 * s^2) * (1 + sin(s) / 2) for s in 0:(2N)]
    z[1] = 0
    inc = phase_increments(z)
    @test isnan(inc[1])
    ref1, ref2 = weighted_phase_frequency(inc, N)
    @test isfinite(ref1) && isfinite(ref2)
    r = WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2))
    st = reducer_init(r, Float64, Val(2), Val(1))
    for (n, x) in enumerate(z)
        st, _ = reducer_step(r, st, n, SVector(real(x), imag(x)), SMatrix{2, 1, Float64}(1, 0))
    end
    f = reducer_final(r, st)
    @test f.complete && isequal(f.nu1[1], ref1) && isequal(f.nu2[1], ref2)
    # the converged stationary point lies in the bracket nu_c +- 1/M of the spec (NF2)
    M = 8 * 256
    for nu in (0.1234567891, -0.3, 0.45), a in (0.0, 0.3)
        zz = ComplexF64[cis(2 * pi * nu * j) + a * cis(2 * pi * 0.05 * j) for j in 0:255]
        w = window_frequency(zz)
        d = abs(w.nu - (-0.5 + w.bin / M))
        @test !w.refined || min(d, 1 - d) <= 1 / M * (1 + 1e-9)
    end
    # pad below 2 is an argument error naming pad (the grid must be finer than the bin); 2 is accepted
    zp = ComplexF64[cis(2 * pi * 0.2 * j) for j in 0:127]
    for pad in (1, 0)
        @test (e = try window_frequency(zp; pad = pad) catch err err end) isa ArgumentError && occursin("pad", e.msg)
        @test (e = try frequency_map(zp, 64; pad = pad) catch err err end) isa ArgumentError && occursin("pad", e.msg)
    end
    @test window_frequency(zp; pad = 2).status === :ok && frequency_map(zp, 64; pad = 2).status == [:ok]
    # observable indices outside the orbit point are an argument error naming the observable, from the caller-side check (reducer_init runs in the device kernel and must not build
    # messages); the point has d entries, 2 d for the wrapped maps
    named(f) = (e = try f() catch err err end; e isa ArgumentError && occursin("observable", e.msg))
    chk(o, d) = check_run(WeightedPhaseAdvance(; window = 16, observable = o), 1, 33; d = d)
    @test named(() -> chk(PlanarObservable(1, 9), 2))
    @test named(() -> chk(PlanarObservable(0, 2), 2))
    @test named(() -> chk(LiftedAngleObservable(1, 5), 2))
    @test chk(LiftedAngleObservable(1, 3), 2) === nothing && chk(PlanarObservable(1, 4), 4) === nothing && chk(PlanarObservable(1, 9), nothing) === nothing
end
