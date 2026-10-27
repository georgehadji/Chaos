# Self-test of the reference code of nf_helpers.jl against hand-computed values. It uses no name of the module under test (only verified modules: Systems, Tangent), so it must pass before and
# after stage 03 (GREEN in the RED run). It also records, as a derived and gated statement, what the series builders deliver (criterion 7(c) concerns the caller's series, not the module).

@testset "SELF reference code against hand values" begin
    # cosine windows NF20: p = 0 is all ones; p = 1, N = 4, positions (j + 1/2)/N: 1 - cos(pi (2 j + 1)/4); positions j/N: (0, 1, 2, 1); sums are N for every p, N > p (D80b(ii), derived)
    @test all(==(1), nf_chi(0, 8, :half))
    c1 = nf_chi(1, 4, :half)
    @test maximum(abs.(c1 .- [1 - sqrt(BigFloat(2)) / 2, 1 + sqrt(BigFloat(2)) / 2, 1 + sqrt(BigFloat(2)) / 2, 1 - sqrt(BigFloat(2)) / 2])) < 1e-70
    @test maximum(abs.(nf_chi(1, 4, :edge) .- [0, 1, 2, 1])) < 1e-70
    for p in 0:5, pos in (:half, :edge), N in (16, 17, 64)
        c = nf_chi(p, N, pos)
        @test abs(sum(c) - N) < 1e-60 && all(>=(0), c)
        # symmetric about the centre: (half) c_j = c_{N-1-j}; (edge) c_j = c_{N-j} for j = 1..N-1
        sym = pos === :half ? all(j -> abs(c[j] - c[N + 1 - j]) < 1e-70, 1:N) : all(j -> abs(c[j + 1] - c[N - j + 1]) < 1e-70, 1:(N - 1))
        @test sym
    end
    # the constant (2^p (p!)^2/(2p)!): p = 1 gives 1 (so chi_1 = 1 + cos), p = 2 gives 2/3 (chi_2 = (2/3)(1 + cos)^2 = (8/3) cos^4)
    @test abs(nf_chi(1, 4, :edge)[3] - 2) < 1e-70 && abs(nf_chi(2, 4, :edge)[3] - 8 * BigFloat(1) / 3) < 1e-70
    # the closed form Eq. (57): phi_p(0) = 1 (limit), phi_p''(0) = -(2/pi^2)(pi^2/6 - sum 1/k^2) (Eq. (62)), and the sampled transform of nf_chi at N = 4096 agrees with it within 1e-9 (D80b(iii))
    setprecision(BigFloat, NF_BITS) do
        pb = nf_pi()
        for p in 1:5
            h = BigFloat(10)^-20
            @test abs(nf_phi_closed(p, h) - 1) < 1e-30
            second = (nf_phi_closed(p, h) - 1) * 2 / h^2
            ref = -(2 / pb^2) * (pb^2 / 6 - sum(1 / BigFloat(k)^2 for k in 1:p))
            @test abs(second - ref) < 1e-10 * abs(ref)
        end
        for p in 1:5, x in (0.5, 1.0, 2.5, 7.0, 20.0), pos in (:half, :edge)
            xb = BigFloat(x)
            @test abs(nf_phi_sampled(nf_chi(p, 4096, pos), xb, pos) - nf_phi_closed(p, xb)) < 1e-9
        end
    end
    # circle helpers: reduction into [-1/2, 1/2) and the distance of NF9
    @test nf_red(0.5) == -0.5 && nf_red(-0.5) == -0.5 && nf_red(0.75) == -0.25 && nf_red(1.3) ≈ 0.3 && nf_red(-0.2) == -0.2
    @test isapprox(nf_cd(0.4999999, -0.4999999), 2e-7; atol = 1e-14) && nf_cd(0.1, 0.3) ≈ 0.2 && nf_cd(0.9, 0.1) ≈ 0.2 && nf_cd(0.3, 0.3) == 0
    @test nf_vdist([0.1, 0.4], [0.2, 0.1]) ≈ sqrt(0.1^2 + 0.3^2)
    # S_m(nu) and the derivatives f', f'' against central differences of f = |S_0|^2 in BigFloat, on a short series
    setprecision(BigFloat, NF_BITS) do
        rng = StableRNG(nf_seed(0, 1))
        zb = Complex{BigFloat}[complex(BigFloat(rand(rng) - 0.5), BigFloat(rand(rng) - 0.5)) for _ in 1:12]
        chi = nf_chi(1, 12, :half)
        nu = BigFloat(0.173)
        h = BigFloat(10)^-15
        f(x) = abs2(nf_S(zb, chi, x)[1])
        S = nf_S(zb, chi, nu)
        @test abs(nf_fp(S) - (f(nu + h) - f(nu - h)) / (2h)) < 1e-25 * abs(nf_fp(S))
        @test abs(nf_fpp(S) - (f(nu + h) - 2 * f(nu) + f(nu - h)) / h^2) < 1e-25 * abs(nf_fpp(S))
    end
    # the stationary point of a pure line is the line (N = 64, p = 1, both positions), with amplitude a; the coarse index is the grid node nearest to nu
    nu = nf_dec("0.1234567891")
    for pos in (:half, :edge), p in (0, 1, 2)
        z = nf_series(Float64, nf_line_big(nu, BigFloat(1), BigFloat(0), 64))
        r = nf_stat_ref(z, p, pos, 8)
        @test r.unique && abs(r.nu - nu) < 1e-14 && abs(r.amplitude - 1) < 1e-14 && r.refbin == round(Int, (Float64(nu) + 0.5) * 512)
    end
    zb = nf_line_big(nu, BigFloat(1), BigFloat(0), 64)
    r = setprecision(BigFloat, NF_BITS) do
        chi = nf_chi(1, 64, :half)
        x = nu + BigFloat(10)^-6
        for _ in 1:8
            S = nf_S(zb, chi, x)
            x -= nf_fp(S) / nf_fpp(S)
        end
        abs(x - nu)
    end
    @test r < 1e-60                                                  # the exact line: Newton in BigFloat converges to nu to 1e-60 (the reference is exact for the exact series)
    # two lines: the stationary point is near the stronger line (unique, negative curvature)
    z2 = nf_series(Float64, nf_lines_big([(nu, BigFloat(1), BigFloat(0)), (nu + nf_dec("0.05"), nf_dec("0.1"), BigFloat(1))], 256))
    r2 = nf_stat_ref(z2, 1, :half, 8)
    @test r2.unique && abs(r2.nu - nu) < 1e-5 && abs(r2.amplitude - 1) < 1e-2
    # the coarse search: the two-stage path (N > 1024) agrees with brute force over the whole grid (N = 1100, a chaotic window)
    sys = WrappedStandardMap(2 * pi)
    st = nf_states(sys, torus_state(SVector(1.0, 0.5)), 1100)
    zf = ComplexF64.(nf_obs(nf_obs_cis, st))
    chif = Float64.(nf_chi(1, 1100, :half))
    brute = argmax([nf_abs64(zf, chif, -0.5 + m / 8800) for m in 0:8799]) - 1
    @test nf_coarse(zf, chif, 8) == brute
    # signals: the recurrences agree with the direct exponentials; the chirp has the increments of NF11; the demeaned windows have zero mean
    ln = nf_line_big(nu, nf_dec("1e-3"), BigFloat(1.3), 20000)
    for j in (0, 1, 77, 19999)
        @test abs(ln[j + 1] - nf_dec("1e-3") * nf_cis(2 * nf_pi() * nu * j + BigFloat(1.3))) < 1e-70
    end
    ch = nf_chirp_big(nu, nf_dec("1e-6"), 5000)
    for j in (0, 1, 1000, 4999)
        @test abs(ch[j + 1] - nf_cis(2 * nf_pi() * (nu * j + nf_dec("1e-6") * j^2 / 2))) < 1e-65
    end
    dm = nf_demean_big(nf_lines_big([(nu, BigFloat(1), BigFloat(0))], 128), 64)
    @test abs(sum(dm[1:64])) < 1e-70 && abs(sum(dm[65:128])) < 1e-70
    # weighted Birkhoff weights and NF6: w(0) = 0, w(1/2) = exp(-4), symmetric, normalised; constant increments give c / (2 pi); NF15 tends to rho (1e-6 at N = 500, derived scale of the exponential weight)
    w = nf_wb_raw(16)
    @test w[1] == 0 && abs(w[9] - exp(-BigFloat(4))) < 1e-70 && all(j -> abs(w[j] - w[18 - j]) < 1e-70, 2:16)
    @test abs(sum(nf_wb_hat(16)) - 1) < 1e-70
    c = nf_dec("0.7")
    n1, n2 = nf_E2_ref(fill(c, 32), 16)
    @test abs(n1 - c / (2 * nf_pi())) < 1e-70 && abs(n2 - c / (2 * nf_pi())) < 1e-70
    incr, rho = nf_nf15_big(1000)
    e1, e2 = nf_E2_ref(incr, 500)
    @test abs(rho - (sqrt(BigFloat(2)) - 1)) < 1e-70 && abs(e1 - rho) < 1e-6 && abs(e2 - rho) < 1e-6
    # planar increments: the phases 0.3, 1.3, 1.9 give 1.0 and 0.6 (radii irrelevant)
    zi = [2.0 * nf_cis(0.3), 0.5 * nf_cis(1.3), 7.0 * nf_cis(1.9)]
    @test maximum(abs.(nf_incr_planar(zi) .- [1.0, 0.6])) < 1e-14
    # Y5 rotation block: z' = exp(-i theta) z (NF14), the reduced frequency -theta / (2 pi), and the direct sum of two blocks
    for th in (0.3, 1.0, 2.5, 3.1)
        m = nf_rot_map(Float64, th)
        z0 = SVector(0.4, -0.7)
        z1 = forward(m, z0)
        @test abs(complex(z1[1], z1[2]) - exp(-im * th) * complex(z0[1], z0[2])) < 4e-16
        @test abs(nf_rot_nu(th) + BigFloat(th) / (2 * nf_pi())) < 1e-70
    end
    m2 = nf_rot_map(Float64, 0.7, 2.3)
    z0 = SVector(0.1, 0.2, 0.3, 0.4)                                 # (q1, q2, p1, p2)
    z1 = forward(m2, z0)
    @test abs(complex(z1[1], z1[3]) - exp(-0.7im) * complex(z0[1], z0[3])) < 4e-16 && abs(complex(z1[2], z1[4]) - exp(-2.3im) * complex(z0[2], z0[4])) < 4e-16
    @test nf_rot_nu(3.3) > 0 && abs(nf_rot_nu(2pi - 0.1) - 0.1 / (2 * nf_pi())) < 1e-15           # a block angle above pi gives a positive frequency after the reduction
    # K = 0 map (NF13): p is constant exactly, the lifted increments are p0 within 4 eps (one rounding of the sum and of the reduction), and the series builders run
    for T in NF_TYPES
        epsT = nf_eps(T)
        sys0 = WrappedStandardMap(nf_K(T, :zero))
        s0 = torus_state(SVector{2, T}(T(1.7), T(3.3)))
        sts = nf_states(sys0, s0, 400)
        @test all(s -> s[2] == s0[2], sts)
        inc = nf_lifted_incr(sts)
        @test maximum(abs.(inc .- s0[2])) <= 4 * epsT * Float64(s0[2]) * 2 + 0.0
        incb = nf_lifted_incr_big(sts)
        @test maximum(abs.(incb .- nf_big(s0[2]))) <= 4 * epsT * 8
        @test length(nf_obs(nf_obs_cis, sts)) == 400 && length(nf_states_tau(sys0, s0, 10, 4)) == 10
    end
    # criterion 7(c) concerns the caller's series: the per-sample error of exp(i q) built from the wrapped coordinate of the K = 0 orbit against the exact orbit is NOT within 4 eps for long
    # windows (the wrapped orbit accumulates the rounding of 2 pi and of the sums); the derived bound below is what the builder delivers and is gated here (see the plan, Items for the owner).
    for T in NF_TYPES
        epsT = nf_eps(T)
        p0 = T(3.3)
        sys0 = WrappedStandardMap(nf_K(T, :zero))
        q0 = T(1.7)
        sts = nf_states(sys0, torus_state(SVector{2, T}(q0, p0)), 4096)
        worst, wraps = setprecision(BigFloat, NF_BITS) do
            tpT = nf_big(lift(SVector{4, T}(0, 0, 1, 0))[1])
            w = 0.0
            for (j, s) in enumerate(sts)
                qex = nf_big(q0) + j * nf_big(p0)                    # exact lifted orbit; the reduction modulo 2 pi is the exp
                zex = nf_cis(qex)
                w = max(w, Float64(abs(nf_big(nf_cis(s[1])) - zex)))
            end
            (w, Float64(nf_big(sts[end][3])))
        end
        tp_err = Float64(setprecision(BigFloat, NF_BITS) do
            abs(nf_big(lift(SVector{4, T}(0, 0, 1, 0))[1]) - 2 * nf_pi())
        end)
        # derived bound: 4 eps for the cosine and sine, |2 pi_T - 2 pi| per winding (the systematic part), 64 eps for the roundings of the additions
        @test worst <= 4 * epsT + wraps * tp_err + 64 * epsT
        @info "SELF 7(c): K = 0 builder, per-sample error of exp(i q), max over 4096 samples" T worst_in_eps = worst / epsT windings = wraps
    end
    # the stream helper: nf_run with an existing reducer and the trace reader gives the final point of the trace
    sysw = WrappedStandardMap(nf_K(Float64, :two_pi))
    x0 = torus_state(SVector(1.0, 0.5))
    x, W, fin = nf_run((nf_Trace(), MEGNO()), sysw, x0, nf_frame(nf_seed(0, 2), 1, Float64), 50, 1)
    @test length(fin[1].xs) == 50 && fin[1].xs[end] == x && fin[2].handoffs == 50
    @test nf_states(sysw, x0, 50)[end] == x
end
