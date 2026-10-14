# AC1 to AC5 of lyap-svd-stabilized, plus the V6 unit. theta_sum as GS G13 with eps_eff; seeds sv_seed(ac).

@testset "API V6 log_singular_values" begin
    # V6 against the singular values of P = [[a, b], [0, delta]], a = e^S1, delta = e^S2, b = c a, at 256 bits
    rng = StableRNG(sv_seed(0))
    for T in gs_types(), c in (T === Float64 ? (0.0, 1.27, -0.19, 3.0e3) : (1.27,)), _ in 1:(T === Float64 ? 20 : 4)
        S = SVector{2,T}(T(15 * rand(rng)), T(-15 * rand(rng)))
        out = log_singular_values(S, T(c))
        ref = setprecision(BigFloat, 256) do
            a, dl = exp(gs_big(S[1])), exp(gs_big(S[2]))
            b = BigFloat(c) * a
            F, Dt = a^2 + b^2 + dl^2, a * dl
            s1 = sqrt((F + sqrt(F^2 - 4 * Dt^2)) / 2)
            (log(s1), log(Dt / s1))
        end
        tol = 16 * gs_eps(T) * (1 + abs(Float64(S[1])) + abs(Float64(S[2])) + log1p(c^2))
        @test out isa SVector{2,T}
        @test gs_absdiff(out[1], ref[1]) <= tol && gs_absdiff(out[2], ref[2]) <= tol
    end
    # no overflow however large S_1: finite, decreasing order, product sigma_1 sigma_2 = e^(S1 + S2)
    out = log_singular_values(SVector(1.0e5, -1.0e5 + 3.0), 1.5)
    @test all(isfinite, out) && out[1] >= out[2] && out[1] >= 1.0e5
    @test abs(out[1] + out[2] - 3.0) <= 16 * eps(Float64) * 1.0e5
    @test all(isfinite, log_singular_values(SVector(0.0, 0.0), 0.0))   # sigma_1 = sigma_2 (V3)
end

@testset "AC1 symmetric map" begin
    # V9a: lambda_1^SVD = ln mu, lambda_2^SVD = -ln mu at every n, tau and frame; lambda^QR carries the offset of G12a.
    seed = sv_seed(1)
    lnmu = gs_lnmu()
    for T in gs_types()
        sys = gs_cat(T)
        x0s = [gs_zero(T, 2) for _ in 1:20]
        for n in (40, 400)
            l1 = Dict{Tuple{Int,Int},Tuple{T,Float64}}()   # (tau, frame) -> (lambda_1^SVD, theta)
            for tau in (1, 2, 4)
                res = sv_run(sys, x0s, n; tau, seed)
                for j in 1:20
                    f = sv_fin(res, j)
                    th = gs_theta(T, f.G_max, tau, 2.0)
                    q1 = initial_deviations(seed, j, Val(2), Val(2), T)[:, 1]
                    @test gs_absdiff(f.lambda[1], lnmu) <= th
                    @test gs_absdiff(f.lambda[2], -lnmu) <= th
                    @test gs_absdiff(f.lambda_qr[1], gs_cat_lambda1(q1, n)) <= th
                    @test f.lambda[1] >= f.lambda[2]
                    l1[(tau, j)] = (f.lambda[1], th)
                end
            end
            for a in keys(l1), b in keys(l1)    # any two frames and any two tau: pair bound of G13
                a < b && @test Float64(abs(l1[a][1] - l1[b][1])) <= l1[a][2] + l1[b][2]
            end
        end
    end
end

@testset "AC2 K = 0" begin
    # V9b: sigma_1 = (n + sqrt(n^2 + 4))/2, sigma_2 = 1/sigma_1; for the frame (e_2, e_1) the SVD/QR offset is
    # [ln sigma_1 - (1/2) ln(n^2 + 1)]/n (G12b). lambda_1(1e4) <= ln(1e4)/1e4 + 2e-12 (benchmarks.md, Skokos2010 Eq. (62)).
    seed = sv_seed(2)
    for T in gs_types()
        sys = StandardMap(zero(T))
        x0 = SVector{2,T}(T(0.3), T(0.7))
        frames = (SMatrix{2,2,T}(0, 1, 1, 0), SMatrix{2,2,T}(1, 0, 0, 1), initial_deviations(seed, 1, Val(2), Val(2), T))
        for (fi, W0) in enumerate(frames), n in (1, 10, 100, 1000, 10^4)
            f = sv_fin(sv_run(sys, [x0], n; tau = 1, seed, W0 = [W0]), 1)
            th = gs_theta(T, f.G_max, 1, 1.0)
            l1 = sv_V9b_lnsigma1(n)
            @test gs_absdiff(f.lambda[1], l1 / n) <= th
            @test Float64(abs(f.lambda[2] + f.lambda[1])) <= th
            if fi == 1
                lq = setprecision(BigFloat, 256) do
                    log(BigFloat(n)^2 + 1) / 2
                end
                off = (l1 - lq) / n
                @test gs_absdiff(f.lambda[1] - f.lambda_qr[1], off) <= th
            end
            n == 10^4 && @test Float64(f.lambda[1]) <= log(1e4) / 1e4 + 2e-12
        end
    end
end

@testset "AC3 identity with Gram-Schmidt and sandwich" begin
    seed = sv_seed(3)
    for K in gs_KS(), tau in (1, 2, 4)
        sys = gs_stdmap(Float64, K)
        x0s = gs_starts(Float64, seed + 1, 50)
        pts = [gs_points(sys, x0s[j], 1000) for j in 1:50]
        res = sv_run(sys, x0s, 1000; tau, seed)
        ref = gs_run(sys, x0s, 1000; k = 2, tau, seed, alg = RGS())
        res500 = sv_run(sys, x0s, 500; tau, seed)
        nshear = nslow = 0
        worstslow = 0.0
        for j in 1:50
            f, g = sv_fin(res, j), gs_fin(ref, j)
            # (a) the QR-family output is the output of lyap-gram-schmidt with RGS, bitwise
            @test isequal(f.sum, g.sum) && isequal(f.comp, g.comp) && isequal(f.lambda_qr, g.lambda)
            # (b) V7: 0 <= lambda_1^SVD - lambda_1^QR - ln(1 + c^2)/(2n) <= u/(2n), up to theta_sum
            th = gs_theta(Float64, f.G_max, tau, gs_Mmax(sys, pts[j]))
            n = 1000
            gap = Float64(f.lambda[1] - f.lambda_qr[1]) - log1p(Float64(f.c)^2) / (2n)
            @test -th <= gap <= Float64(f.u) / (2n) + th
            @test f.lambda[1] >= f.lambda[2] && f.lambda[1] + th >= f.lambda_qr[1]
            # (c) c is settled by the early windows only where S_1 grows (Amendment D35; classes by S_1(500), every orbit of the panel):
            # S_1(500) >= 25: the late increments are absorbed, 8 eps; 5 <= S_1(500) < 25 (shear-dominated): algebraic convergence, 0.06
            # relative (worst measured 0.0261); S_1(500) < 5 (near elliptic points and accelerator modes): c does not converge, reported only.
            h = sv_fin(res500, j)
            S1_500 = Float64(h.sum[1] + h.comp[1])
            dc = Float64(abs(f.c - h.c))
            if S1_500 >= 25
                @test dc <= 8 * eps(Float64) * max(1.0, Float64(abs(f.c)))
            elseif S1_500 >= 5
                nshear += 1
                @test dc <= 0.06 * max(1.0, Float64(abs(h.c)))
            else
                nslow += 1
                worstslow = max(worstslow, dc / max(1.0, Float64(abs(h.c))))
            end
        end
        K == 0.5 && @test nshear >= 4          # the scored shear class is not empty at K = 0.5
        @info "lyap-svd-stabilized AC3(c): K = $K, tau = $tau, shear class scored = $nshear, S_1(500) < 5 reported = $nslow (largest relative change $worstslow)"
    end
end

@testset "AC4 256-bit recomputation" begin
    # product of the 256-bit tangents at the stored orbit points, det Y = 1 used analytically (spec Oracle):
    # |lambda_i^SVD(n) - lambda_i^ref(n)| <= theta_sum, i = 1, 2
    seed = sv_seed(4)
    ns = (12, 20, 100, 1000)
    for (T, Ks, m, taus) in ((Float64, gs_KS(), 20, (1, 2, 4)), (Double64, (2.0, :twopi), 4, (1, 4)))
        for K in Ks
            sys = gs_stdmap(T, K)
            x0s = gs_starts(T, seed + 1, m)
            pts = [gs_points(sys, x0s[j], maximum(ns)) for j in 1:m]
            Ys = [gs_ref_std_Y(sys.K, pts[j], ns) for j in 1:m]
            for n in ns, tau in taus
                res = sv_run(sys, x0s, n; tau, seed)
                for j in 1:m
                    f = sv_fin(res, j)
                    th = gs_theta(T, f.G_max, tau, gs_Mmax(sys, pts[j][1:n]))
                    l1, l2 = sv_ref_logsigma(Ys[j][n])
                    @test gs_absdiff(f.lambda[1], l1 / n) <= th
                    @test gs_absdiff(f.lambda[2], l2 / n) <= th
                end
            end
        end
    end
end

@testset "AC5 negative control and overflow" begin
    # L10b: the plainly multiplied Float64 product loses sigma_2 at n = 20 (|ln sigma_2^naive - ln sigma_2^ref| > 10),
    # overflows before n = 700, while the estimator is finite and within theta_sum of the 256-bit value at n = 20 and 2000.
    # The orbit is the first of the seeded starts with ln sigma_1(20) >= 25 (Amendment D35: the gap exceeds 10 for every orbit with
    # ln sigma_1(20) >= 24.9, and fails for weakly stretching orbits, so the spec's 'one seeded orbit' needs this condition).
    seed = sv_seed(5)
    sys = gs_stdmap(Float64, :twopi)
    x0s = gs_starts(Float64, seed + 1, 60)
    lnsig1(x0) = Float64(sv_ref_logsigma(gs_ref_std_Y(sys.K, gs_points(sys, x0, 20), (20,))[20])[1])
    j = findfirst(x0 -> lnsig1(x0) >= 25, x0s)
    @test j !== nothing
    x0 = x0s[something(j, 1)]
    pts = gs_points(sys, x0, 2000)
    Ys = gs_ref_std_Y(sys.K, pts, (20, 2000))
    naive = svdvals(Matrix(sv_naive_product(sys, pts, 20)))
    l2ref = sv_ref_logsigma(Ys[20])[2]
    @test abs(log(naive[2]) - Float64(l2ref)) > 10
    @test any(!isfinite, sv_naive_product(sys, pts, 700))
    for n in (20, 2000)
        f = sv_fin(sv_run(sys, [x0], n; tau = 1, seed), 1)
        th = gs_theta(Float64, f.G_max, 1, gs_Mmax(sys, pts[1:n]))
        l1, l2 = sv_ref_logsigma(Ys[n])
        @test all(isfinite, f.lambda) && isfinite(f.c)
        @test gs_absdiff(f.lambda[1], l1 / n) <= th && gs_absdiff(f.lambda[2], l2 / n) <= th
    end
end
