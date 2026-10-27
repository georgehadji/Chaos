# API tests: the pure functions (window_weights, circle_distance, frequency_distance, frequency_index, birkhoff_weight, birkhoff_weights, phase_increments, weighted_phase_frequency), the result types of
# window_frequency and frequency_map (fields, defaults, settings echo), and the WeightedPhaseAdvance reducer (constructor, observables, protocol, provenance, partial finals). Every testset starts with a line
# naming the new symbols it uses (UndefVarError when missing). Criteria 14 (independent orbit, ChaosTools, Python NAFF), 15 (CPU vs GPU) and 16 (series feed) are stage 04 or later (the plan maps them).

@testset "API pure functions" begin
    _ = (window_weights, circle_distance, frequency_distance, frequency_index, birkhoff_weight, birkhoff_weights, phase_increments, weighted_phase_frequency)
    for T in NF_TYPES
        epsT = nf_eps(T)
        # window_weights(T, N; order = 1, positions = :half): the Laskar normalisation (sum N); N = 4 and p = 1 gives 1 -+ cos(pi/4) (hand value), the default keywords are order 1 and :half
        w = window_weights(T, 4)
        @test w isa Vector{T} && length(w) == 4
        c = 1 - sqrt(BigFloat(2)) / 2
        @test Float64(abs(nf_big(w[1]) - c) / c) <= 4 * epsT && Float64(abs(nf_big(w[2]) - (2 - c)) / (2 - c)) <= 4 * epsT && w == window_weights(T, 4; order = 1, positions = :half)
        @test window_weights(T, 8; order = 0) == ones(T, 8)
        we = window_weights(T, 4; positions = :edge)
        @test we[1] == 0 && Float64(abs(nf_big(we[2]) - 1)) <= 4 * epsT && Float64(abs(nf_big(we[3]) - 2)) <= 8 * epsT && we[2] === we[4]
        # circle_distance(a, b) = |a - b - round(a - b)| in the type of the arguments; frequency_distance(nu1, nu2) the Euclidean norm of the plane distances (NF9)
        d = circle_distance(T(0.4999999), T(-0.4999999))
        dref = nf_cd(nf_big(T(0.4999999)), nf_big(T(-0.4999999)))
        @test d isa T && Float64(abs(nf_big(d) - dref)) <= 2 * epsT
        @test circle_distance(T(0.1), T(0.3)) == circle_distance(T(0.3), T(0.1)) && circle_distance(T(0.9), T(0.1)) isa T && circle_distance(T(0.3), T(0.3)) == 0
        @test nf_err(frequency_distance(T[0.1, 0.4], T[0.2, 0.1]), nf_vdist(nf_big.(T[0.1, 0.4]), nf_big.(T[0.2, 0.1]))) <= 4 * epsT && frequency_distance(T[0.1], T[0.1]) == 0
        # frequency_index(delta; floor_multiplier = 32): D = log10 max(delta, 32 eps), flag delta <= 32 eps (NF10); delta is not clipped by the function (it only reports)
        fi = frequency_index(T(1e-3))
        @test fi.D isa T && nf_err(fi.D, log10(nf_big(T(1e-3)))) <= 12 * epsT && !fi.below_floor
        fl = frequency_index(zero(T))
        @test fl.D == log10(T(32) * effective_eps(T)) && fl.below_floor
        @test frequency_index(T(32) * effective_eps(T)).below_floor && !frequency_index(T(33) * effective_eps(T)).below_floor
        @test frequency_index(T(1e-3); floor_multiplier = 1e30).below_floor
        # birkhoff_weight(t, p = 1) = exp(-1/(t^p (1-t)^p)) on (0, 1), 0 elsewhere (NF4)
        @test abs(birkhoff_weight(T(0.5)) - exp(-T(4))) <= 4 * epsT * exp(-T(4)) && birkhoff_weight(T(0)) == 0 && birkhoff_weight(T(1)) == 0 && birkhoff_weight(T(1.5)) == 0 && birkhoff_weight(T(-0.2)) == 0
        @test abs(birkhoff_weight(T(0.5), 2) - exp(-T(16))) <= 4 * epsT * exp(-T(16))
        bw = birkhoff_weights(T, 16)
        @test bw isa Vector{T} && length(bw) == 16 && bw[1] == 0
        ref = nf_wb_hat(16)
        @test maximum(abs.(Float64.([nf_big(bw[j]) - ref[j] for j in 1:16]))) <= 16 * epsT * Float64(maximum(ref))
        @test birkhoff_weights(T, 16; p = 2) isa Vector{T} && abs(sum(nf_big.(birkhoff_weights(T, 16; p = 2))) - 1) <= 4 * epsT
        # phase_increments: the principal argument of z_{j+1}/z_j in (-pi, pi]; hand values (phases 0.3, 1.3, 1.9 as T values: the differences 1.0 and 0.6; radii irrelevant) and the branch cut
        zi = Complex{T}[2 * nf_cis(T(0.3)), nf_cis(T(1.3)) / 2, 7 * nf_cis(T(1.9))]
        inc = phase_increments(zi)
        @test inc isa Vector{T} && length(inc) == 2 && abs(inc[1] - (T(1.3) - T(0.3))) <= 8 * epsT && abs(inc[2] - (T(1.9) - T(1.3))) <= 8 * epsT
        zc = Complex{T}[nf_cis(T(3.0)), nf_cis(T(-3.0))]                  # a step of 0.2830... through the branch cut at +-pi, not of -6
        @test abs(phase_increments(zc)[1] - (2 * T(pi) - 6)) <= 8 * epsT
        # weighted_phase_frequency(increments, N; weights = birkhoff_weights(T, N)): the pair (nu_1, nu_2); constant increments c give c/(2 pi) in both windows, a custom weight vector is used as given
        nu1, nu2 = weighted_phase_frequency(fill(T(0.7), 64), 32)
        @test nu1 isa T && nu2 isa T && abs(nu1 - T(0.7) / (2 * T(pi))) <= 4 * epsT * Float64(nu1) && nu1 == nu2
        eq = fill(inv(T(32)), 32)
        inc2 = T.(0:63)
        u1, u2 = weighted_phase_frequency(inc2, 32; weights = eq)
        @test nf_err(u1, nf_dec("15.5") / (2 * nf_pi())) <= 4 * epsT * Float64(u1) && nf_err(u2, nf_dec("47.5") / (2 * nf_pi())) <= 4 * epsT * Float64(u2)
    end
end

@testset "API result types and defaults (E1)" begin
    _ = (window_frequency, frequency_map)
    # window_frequency(z; order = 1, positions = :half, pad = 8, mean_removal = false, mirror_ratio = 0.9) -> a record with nu, amplitude, mirror_amplitude, second_ratio, bin (0-based index of the coarse grid
    # nu_m = -1/2 + m/(pad N)), tie, near_zero, mirror, refined, status. frequency_map(series, N; ..., floor_multiplier = 32) -> w1, w2 (one record per plane), dist (per plane, after the mirror rule),
    # delta, D, below_floor, status (per plane), failure, settings (the identifiers of P11).
    N = 256
    nu = nf_dec("0.1234567891")
    for T in NF_TYPES
        epsT = nf_eps(T)
        z = nf_series(T, nf_lines_big([(nu, BigFloat(1), BigFloat(0)), (nu + nf_dec("0.3"), nf_dec("0.1"), BigFloat(0.5))], 2N))
        w = window_frequency(z[1:N])
        @test w.nu isa T && w.amplitude isa T && w.mirror_amplitude isa T && w.second_ratio isa T
        @test w.bin isa Integer && 0 <= w.bin < 8 * N && w.status === :ok && w.refined && !w.near_zero && !w.mirror && !w.tie
        @test abs(w.bin - round(Int, (Float64(nu) + 0.5) * 8 * N)) <= 1
        @test Float64(abs(nf_big(w.amplitude) - 1)) < 1e-3                              # one line of amplitude 1 plus a weaker line 0.3 away (leakage far below)
        @test 0.095 < Float64(w.second_ratio) < 0.1005                                  # the weaker line: 0.1 of the peak (grid scalloping below 0.2 percent)
        @test Float64(w.mirror_amplitude) < 1e-4                                        # -nu holds no line
        wo = window_frequency(z[1:N]; order = 1, positions = :half, pad = 8, mean_removal = false, mirror_ratio = 0.9)
        @test isequal(nf_wkey(wo), nf_wkey(w))                                           # the defaults are these values
        r = frequency_map(z, N)
        @test length(r.w1) == 1 && length(r.w2) == 1 && length(r.dist) == 1 && length(r.status) == 1 && r.failure === nothing
        @test isequal(nf_wkey(r.w1[1]), nf_wkey(w)) && r.delta isa T && r.D isa T
        s = r.settings
        @test s.estimator == "windowed-fourier" && s.weight == "cosine" && s.window_order == 1 && s.positions === :half && s.pad == 8 && s.mean_removal === false
        @test s.mirror_ratio == 0.9 && s.floor_multiplier == 32 && s.window == N && s.eltype == precision_label(T) && s.planes == 1 && s.min_window == 16
        r2 = frequency_map(z, N; order = 3, positions = :edge, pad = 4, mean_removal = true, mirror_ratio = 0.5, floor_multiplier = 8)
        s2 = r2.settings
        @test s2.window_order == 3 && s2.positions === :edge && s2.pad == 4 && s2.mean_removal === true && s2.mirror_ratio == 0.5 && s2.floor_multiplier == 8 && s2.min_window == 32
        # mean_removal = true equals the series with the window means removed beforehand (within 16 eps)
        zb = Complex{BigFloat}[nf_big(x) for x in z]
        rm = frequency_map(z, N; mean_removal = true)
        rp = frequency_map(nf_series(T, nf_demean_big(zb, N)), N)
        @test nf_cderr(rm.w1[1].nu, nf_big(rp.w1[1].nu)) <= 16 * epsT && nf_cderr(rm.w2[1].nu, nf_big(rp.w2[1].nu)) <= 16 * epsT
        # a list of planes; the vector of one plane equals the plain series
        r3 = frequency_map([z, z], N)
        @test length(r3.w1) == 2 && length(r3.dist) == 2 && r3.status == [:ok, :ok] && r3.settings.planes == 2 && isequal(nf_wkey(r3.w2[2]), nf_wkey(r.w2[1]))
        @test isequal(frequency_map([z], N).delta, r.delta)
        # a series longer than 2N uses the first 2N samples
        @test isequal(frequency_map(vcat(z, z), N).delta, r.delta)
    end
end

@testset "API reducer constructor and observables" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable, reducer_provenance)
    # WeightedPhaseAdvance(; window, observable, birkhoff_order = 1): window is N, or a tuple of N for several horizons sharing one run (each horizon needs its own sums); an isbits reducer; provenance
    # keys merged into the propagate record.
    po, lo = PlanarObservable(1, 3), LiftedAngleObservable(1, 3)
    @test isbitstype(typeof(po)) && isbitstype(typeof(lo)) && isbitstype(typeof(LiftedAngleObservable(1)))
    r = WeightedPhaseAdvance(; window = 100, observable = po)
    @test isbitstype(typeof(r))
    pr = reducer_provenance(r)
    @test pr.indicator == "naff-fma" && pr.estimator == "weighted-phase-advance" && pr.weight == "wb_exp" && pr.birkhoff_order == 1 && pr.window == (100,)
    @test pr.observable == "planar" && pr.sample_position == "j/N" && pr.handoffs == 201 && pr.floor_multiplier == 32 && pr.min_window == 16
    pl = reducer_provenance(WeightedPhaseAdvance(; window = (50, 100), observable = lo, birkhoff_order = 2))
    @test pl.window == (50, 100) && pl.observable == "lifted-angle" && pl.birkhoff_order == 2 && pl.handoffs == 201
    @test reducer_workspace(r, Float64, Val(2), Val(1)) == 0 && reducer_workspace(r, Double64, Val(4), Val(3)) == 0
    # the reducer is a reader: any k >= 1 is accepted
    for k in (1, 2, 3)
        @test reducer_init(r, Float64, Val(4), Val(k)) isa Any
    end
    # the observables: the planar increment is the principal argument of z'/z with z = x[iq] + i x[ip]; the lifted increment is (q' - q) + 2 pi (m' - m); both are exposed through the reducer only
    # (see the protocol testset); a state too short for the indices is an error of the call, not of the constructor
    @test PlanarObservable(1, 2) == PlanarObservable(1, 2) && LiftedAngleObservable(1, 3) == LiftedAngleObservable(1, 3)
end

@testset "API reducer protocol on hand arithmetic" begin
    _ = (WeightedPhaseAdvance, PlanarObservable, LiftedAngleObservable)
    # N = 16 (the minimum window), planar observable on d = 2, phases phi_s = 0.3 s + 0.02 s^2 (s = 0..32), radii 1 + sin(s)/2: the increments Delta_j = 0.3 + 0.02 (2 j + 1) are linear in j, the weights
    # are symmetric (mean of j equal to N/2), so nu_1 = (0.3 + 0.02 (N + 1))/(2 pi) and nu_2 = (0.3 + 0.02 (2 (N + N/2) + 1))/(2 pi) (hand values), and they match the 512-bit evaluation of NF6.
    # The first handoff is the sample 0; nu_1 appears at the handoff N + 1, nu_2 and the pair distance at 2 N + 1; nothing is clipped; W is returned unchanged by every step.
    N = 16
    for T in NF_TYPES
        epsT = nf_eps(T)
        r = WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2))
        W = SMatrix{2, 1, T}(1, 0)
        xs = [SVector{2, T}(T(1 + sin(s) / 2) * T(cos(0.3 * s + 0.02 * s^2)), T(1 + sin(s) / 2) * T(sin(0.3 * s + 0.02 * s^2))) for s in 0:(2N)]
        st = reducer_init(r, T, Val(2), Val(1))
        f0 = reducer_final(r, st)
        @test f0.handoffs == 0 && f0.n == 0 && !f0.complete && all(isnan, f0.nu1) && all(isnan, f0.nu2) && all(isnan, f0.delta)
        zs = Complex{T}[complex(x[1], x[2]) for x in xs]
        ref = nf_E2_ref(nf_incr_planar(zs), N; bits = 512)
        B = nf_BE2(T, nf_incr_planar(zs))
        for n in 1:(2N + 1)
            st, W2 = reducer_step(r, st, n, xs[n], W)
            @test W2 === W
            f = reducer_final(r, st)
            @test f.handoffs == n && f.n == n && f.complete == (n == 2N + 1)
            n < N + 1 && @test isnan(f.nu1[1])
            if n >= N + 1
                @test !isnan(f.nu1[1]) && Float64(abs(nf_big(f.nu1[1]) - ref[1])) <= B
            end
            n < 2N + 1 && @test isnan(f.nu2[1]) && isnan(f.delta[1]) && isnan(f.D[1])
            reducer_failed(r, st) && break
        end
        f = reducer_final(r, st)
        @test !reducer_failed(r, st) && f.complete
        hand1 = (0.3 + 0.02 * (N + 1)) / (2 * pi)
        hand2 = (0.3 + 0.02 * (2 * (N + N / 2) + 1)) / (2 * pi)
        @test abs(Float64(f.nu1[1]) - hand1) < 1e-12 && abs(Float64(f.nu2[1]) - hand2) < 1e-12
        @test Float64(abs(nf_big(f.nu2[1]) - ref[2])) <= B
        @test Float64(abs(nf_big(f.delta[1]) - abs(ref[1] - ref[2]))) <= 2 * B && !f.below_floor[1]
        @test Float64(abs(nf_big(f.D[1]) - log10(abs(ref[1] - ref[2])))) <= 1e-6
        # one more handoff after the horizon is complete does not change the result
        st2, _ = reducer_step(r, st, 2N + 2, xs[2N + 1], W)
        @test isequal(reducer_final(r, st2).nu2, f.nu2) && isequal(reducer_final(r, st2).nu1, f.nu1)
        # two horizons share the run: window = (16, 32) needs 65 handoffs; the 16 horizon completes at 33 and its values equal the one-horizon reducer's
        r2 = WeightedPhaseAdvance(; window = (N, 2N), observable = PlanarObservable(1, 2))
        xl = vcat(xs, [SVector{2, T}(T(1 + sin(s) / 2) * T(cos(0.3 * s + 0.02 * s^2)), T(1 + sin(s) / 2) * T(sin(0.3 * s + 0.02 * s^2))) for s in (2N + 1):(4N)])
        s2 = reducer_init(r2, T, Val(2), Val(1))
        for n in 1:(4N + 1)
            s2, _ = reducer_step(r2, s2, n, xl[n], W)
        end
        g = reducer_final(r2, s2)
        @test length(g.nu1) == 2 && g.complete && isequal(g.nu1[1], f.nu1[1]) && isequal(g.nu2[1], f.nu2[1]) && !isnan(g.nu2[2])
    end
    # the lifted-angle observable: increments (q' - q) + 2 pi (m' - m) from the wrapped coordinate and its winding count (hand case: q advances 6.0 per step, wraps once every step)
    tp = 2 * pi
    r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
    W = SMatrix{2, 1, Float64}(1, 0)
    st = reducer_init(r, Float64, Val(2), Val(1))
    q, m = 0.5, 0
    for n in 1:(2N + 1)
        st, _ = reducer_step(r, st, n, SVector(q, 0.0, Float64(m), 0.0), W)
        q += 6.0
        if q >= tp
            q -= tp
            m += 1
        end
    end
    f = reducer_final(r, st)
    @test abs(f.nu1[1] - 6.0 / tp) < 4e-15 && abs(f.nu2[1] - 6.0 / tp) < 4e-15
end
