# Amendment D81: D81a (the printed bound of the equal-weight phase advance), D81b (the sampled Hanning spectrum, Bartolini1996 Eq. (24)) and the two clauses of criterion 11(b) that the amendment
# lists without a test (an identically zero series with mean removal off; the planar observable at the origin in the reducer).
include(joinpath(@__DIR__, "nf_d81_helpers.jl"))

@testset "D81 SELF Eq. (24) against the direct sum" begin
    # Green now (no new symbol): the reference Eq. (24) of nf_d81_helpers.jl equals the modulus of the direct sum (1/N) sum chi_j z_j exp(-2 pi i nu_j j), chi_j = 2 sin^2(pi j/N), at nu_j = -nu_0,
    # delta = 2 nu_0, for N in {16, 37, 100} (the spec's drafting check, here at 1e-60), and the two-sine generator has 2N + 1 samples.
    for N in (16, 37, 100), ns in ("0.1234567891", "-0.3711", "0.2")
        nu0 = nf_dec(ns)
        zb = nf_line_big(nu0, BigFloat(1), BigFloat(0), N)
        chi = nf_chi(1, N, :edge)
        direct = setprecision(BigFloat, NF_BITS) do
            abs(nf_S(zb, chi, -nu0)[1]) / N
        end
        @test abs(direct - nf_bart_eq24(N, 2 * nu0)) < BigFloat(10)^-60
    end
    @test length(nf_d81_two_sine_big(nf_dec("0.1"), nf_dec("0.5"), nf_dec("0.01"), 16)) == 33
end

@testset "D81a equal-weight phase advance, printed bound" begin
    _ = (weighted_phase_frequency, phase_increments)
    # z_n = exp(2 pi i nu_0 n) + b_0 exp(2 pi i nu_1 n), n = 0..2N, nu_0 = 0.1234567891, b_0 in {0.1, 0.5, 0.9}, nu_1 = nu_0 - Delta, Delta in {1e-4, 1e-3, 0.0113, 0.05, 0.3, -0.2}, N in {16, 100, 1000},
    # generated in BigFloat(256) and rounded once. E2 with equal weights on each window of N increments: |nu_est - nu_0| <= 1/(2N) + 16 eps in both windows (Bartolini1996 Eq. (7) and (A.5) with
    # N increments). The exponential weight on the same data is reported by the spec, not gated, and is not written here. The largest error times 2N is reported (the spec's drafting value: 0.63).
    nu0 = nf_dec("0.1234567891")
    for T in NF_TYPES
        epsT = nf_eps(T)
        bad = String[]
        worst = 0.0
        for N in (16, 100, 1000), b0s in ("0.1", "0.5", "0.9"), Ds in ("1e-4", "1e-3", "0.0113", "0.05", "0.3", "-0.2")
            z = nf_series(T, nf_d81_two_sine_big(nu0, nf_dec(b0s), nf_dec(Ds), N))
            nus = weighted_phase_frequency(phase_increments(z), N; weights = ones(T, N))
            for (k, nuk) in enumerate(nus)
                e = nf_err(nuk, nu0)
                worst = max(worst, e * 2N)
                e <= nf_d81a_bound(N, epsT) || push!(bad, "T=$T N=$N b0=$b0s Delta=$Ds window $k error*2N = $(e * 2N)")
            end
        end
        @test nf_report(bad, "D81a $T")
        @info "D81a (reported) largest error times 2N over the grid" T worst
    end
end

@testset "D81b sampled Hanning spectrum, Eq. (24)" begin
    _ = (frequency_map,)
    # Window order 1, sample positions j/N (:edge), mean removal off, z_j = exp(2 pi i nu_0 j), nu_0 in {0.1234567891, -0.3711, 0.2} generated in BigFloat(256) and rounded once, N in {64, 256, 1000}:
    # the mirror amplitude (|phi(-nu)| divided by the sum of the weights) of each window equals Bartolini1996 Eq. (24) at delta = 2 nu_0, evaluated in BigFloat(256) from the exact nu_0, within 64 eps
    # absolute. (The circle distance of nu_0 and -nu_0 is at least 0.24, above 4/N, so the mirror amplitude is defined.)
    for T in NF_TYPES
        epsT = nf_eps(T)
        bad = String[]
        for N in (64, 256, 1000), ns in ("0.1234567891", "-0.3711", "0.2")
            nu0 = nf_dec(ns)
            z = nf_series(T, nf_line_big(nu0, BigFloat(1), BigFloat(0), 2N))
            r = frequency_map(z, N; order = 1, positions = :edge, mean_removal = false)
            ref = nf_bart_eq24(N, 2 * nu0)
            for (k, w) in enumerate((r.w1[1], r.w2[1]))
                w.status === :ok || push!(bad, "T=$T N=$N nu0=$ns window $k status $(w.status)")
                nf_err(w.mirror_amplitude, ref) <= 64 * epsT || push!(bad, "T=$T N=$N nu0=$ns window $k mirror amplitude $(Float64(w.mirror_amplitude)) against $(Float64(ref))")
            end
        end
        @test nf_report(bad, "D81b $T")
    end
end

@testset "AC11b zero series, mean removal off" begin
    _ = (window_frequency, frequency_map)
    # Criterion 11(b) as amended (Item 8): a window is degenerate when its samples, after the optional removal of the window mean, are all equal and either mean removal is on or the common value is zero.
    # An identically zero series with mean removal OFF is therefore degenerate (the amendment lists it as a clause without a test): status :degenerate, nu, amplitude, mirror amplitude and second-line
    # ratio undefined (NaN), no failure recorded, no delta or D for the pair, the pair status naming the plane (6(d)). With mean removal on the same holds. A zero window next to a good one makes the
    # pair, and the plane, degenerate; the plane that is fine keeps its status in a two-plane series.
    N = 64
    for T in NF_TYPES, mr in (false, true)
        z0 = zeros(Complex{T}, 2N)
        w = window_frequency(z0[1:N]; mean_removal = mr)
        @test w.status === :degenerate && isnan(w.nu) && isnan(w.amplitude) && isnan(w.mirror_amplitude) && isnan(w.second_ratio)
        r = frequency_map(z0, N; mean_removal = mr)
        @test r.status[1] === :degenerate && r.w1[1].status === :degenerate && r.w2[1].status === :degenerate
        @test isnan(r.w1[1].nu) && isnan(r.w2[1].nu) && isnan(r.dist[1]) && isnan(r.delta) && isnan(r.D) && r.failure === nothing
        zl = nf_series(T, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), 2N))
        zh = copy(zl)
        zh[(N + 1):end] .= zero(Complex{T})
        rh = frequency_map(zh, N; mean_removal = mr)
        @test rh.w1[1].status === :ok && rh.w2[1].status === :degenerate && rh.status[1] === :degenerate && isnan(rh.delta) && isnan(rh.D) && rh.failure === nothing
        r2 = frequency_map([zl, z0], N; mean_removal = mr)
        @test r2.status == [:ok, :degenerate] && r2.w1[1].status === :ok && r2.w1[2].status === :degenerate && isnan(r2.delta) && isnan(r2.D) && r2.failure === nothing
    end
end

@testset "AC11b planar origin in the reducer" begin
    _ = (WeightedPhaseAdvance, PlanarObservable)
    # The planar observable at the origin (z = 0 has no phase increment) is recorded by the reducer as a degenerate window with the frequency undefined and no failure (the amendment lists it as a clause
    # without a test; stage 03 behaviour, a proposal). The reducer exposes no status field, so "degenerate" reads: reducer_failed false, the run complete, and nu of the window that holds the increments
    # to and from the origin sample NaN, with the pair distance and D NaN. Hand samples z_s = r_s exp(0.3 i s), s = 0..2N, N = 16 (the minimum window): (i) the origin at sample 5 (both increments in
    # window 1): nu1 NaN, nu2 bitwise that of the clean run; (ii) at sample 22 (both increments in window 2): nu2 NaN, nu1 bitwise that of the clean run; (iii) the whole orbit at the origin (the fixed point
    # of a rotation block, through the stream): both windows undefined, no failure, no throw.
    N = 16
    for T in NF_TYPES
        xs(s) = SVector{2, T}(T(1 + sin(s) / 2) * T(cos(0.3 * s)), T(1 + sin(s) / 2) * T(sin(0.3 * s)))
        samples = [xs(s) for s in 0:(2N)]
        r = WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2))
        W = SMatrix{2, 1, T}(1, 0)
        function runhand(pts)
            st = reducer_init(r, T, Val(2), Val(1))
            for n in eachindex(pts)
                st, _ = reducer_step(r, st, n, pts[n], W)
            end
            return (reducer_failed(r, st), reducer_final(r, st))
        end
        failed, clean = runhand(samples)
        @test !failed && clean.complete && abs(Float64(clean.nu1[1]) - 0.3 / (2 * pi)) < 1e-14 && abs(Float64(clean.nu2[1]) - 0.3 / (2 * pi)) < 1e-14
        for (s0, in_window_2) in ((5, false), (22, true))
            pts = copy(samples)
            pts[s0 + 1] = SVector{2, T}(0, 0)
            failed, f = runhand(pts)
            @test !failed && f.complete && f.handoffs == 2N + 1
            @test isnan(f.delta[1]) && isnan(f.D[1])
            if in_window_2
                @test isequal(f.nu1[1], clean.nu1[1]) && isnan(f.nu2[1])
            else
                @test isnan(f.nu1[1]) && isequal(f.nu2[1], clean.nu2[1])
            end
        end
        sys = nf_rot_map(T, 0.7)
        x0 = SVector{2, T}(0, 0)
        _, _, fin = nf_run((r,), sys, x0, nf_frame(nf_seed(23, 1), 1, T, 2, 1), 2N + 1, 1)
        f = fin[1]
        @test f.complete && f.handoffs == 2N + 1 && isnan(f.nu1[1]) && isnan(f.nu2[1]) && isnan(f.delta[1]) && isnan(f.D[1])
    end
end
