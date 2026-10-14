# AC8: cross-estimator agreement on the identical orbit (L1 S2): same system, start, seed, initial basis and tau; Float64.
# Sibling names are taken from the Batch B API contract: Benettin, QRSpectrum(CGS/MGS/RGS), SVDSpectrum, propagate_shadow.
# Amendment D35: no sample-fraction tail rule; the bounds are the amended ones (theta_psi(n) = theta_sum + 8 eps Psi_n/n for the exponents, 16 eps Psi_n/n for k = 1).

@testset "AC8a Benettin, k = 1" begin
    # K in {0.5, 2, 5, 2 pi}, 50 orbits, n = 1e4, tau in {1, 4}: |lambda_1^QR - X| <= (B_S^Psi,QR + B_S^Psi,L1)/n = 16 eps Psi_n / n (S4(a) with the sharp form of
    # L1 S3 as amended, Amendment D35; Psi_n of the orbit and the seeded first column); every orbit, no neutral-orbit special case, no tail rule.
    seed = hh_seed(8, 1)
    n = 10_000
    eps64 = Float64(effective_eps(Float64))
    for (iK, K) in enumerate((0.5, 2.0, 5.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 4)
            X = propagate(sys, starts, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
            lam = [X.finals[j][1].lambda[1] for j in 1:50]
            bound = [16 * eps64 * Psi[j] / n for j in 1:50]
            for alg in hh_algs(Float64)
                Q = propagate(sys, starts, n; k = 1, tau = tau, seed = seed, reducers = (QRSpectrum(alg),))
                @test all(j -> abs(Q.finals[j][1].lambda[1] - lam[j]) <= bound[j], 1:50)
            end
        end
    end
end

@testset "AC8b Gram-Schmidt variants, k = 2" begin
    # K in {0.5, 2, 5, 10, 2 pi}, 50 orbits, n = 1e3, tau in {1, 2, 4}: |lambda_i^QR - lambda_i^GS| <= theta_psi^QR(n) + theta_psi^GS(n), i = 1, 2 (Amendment D35)
    seed = hh_seed(8, 2)
    n = 1000
    for (iK, K) in enumerate((0.5, 2.0, 5.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Mmax = [hh_Mmax(sys, z, n) for z in starts]
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 2, 4)
            runq(alg) = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals
            theta(f, j) = hh_theta_psi(Float64, f[j][1].G_max, tau, Mmax[j], Psi[j], n)       # theta_psi(n), Amendment D35
            gs = [(alg, runq(alg)) for alg in (CGS(), MGS(), RGS())]
            for alg in hh_algs(Float64)
                fq = runq(alg)
                for (_, fg) in gs
                    @test all(j -> maximum(abs(fq[j][1].lambda[i] - fg[j][1].lambda[i]) for i in 1:2) <= theta(fq, j) + theta(fg, j), 1:50)
                end
            end
        end
    end
end

@testset "AC8c SVD exponents, d = 2" begin
    # K in {0.5, 2, 5, 2 pi}, 50 orbits, tau in {1, 2, 4}, n = 1e3 (SVD V7, c_n and u_n from the SVD output):
    # -(th^QR + th^SVD) <= lambda_1^SVD - lambda_1^QR - ln(1 + c^2)/(2n) <= u/(2n) + th^QR + th^SVD with th = theta_psi(n) (Amendment D35),
    # and |lambda_1^SVD + lambda_2^SVD| <= theta_sum (strict)
    seed = hh_seed(8, 3)
    n = 1000
    for (iK, K) in enumerate((0.5, 2.0, 5.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Mmax = [hh_Mmax(sys, z, n) for z in starts]
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 2, 4)
            fs = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (SVDSpectrum(),)).finals
            ths = [hh_theta_psi(Float64, fs[j][1].G_max, tau, Mmax[j], Psi[j], n) for j in 1:50]
            @test all(j -> abs(fs[j][1].lambda[1] + fs[j][1].lambda[2]) <= hh_theta(Float64, fs[j][1].G_max, tau, Mmax[j]), 1:50)
            for alg in hh_algs(Float64)
                fq = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals
                function sandwich(j)
                    f = fs[j][1]
                    thq = hh_theta_psi(Float64, fq[j][1].G_max, tau, Mmax[j], Psi[j], n)
                    mid = f.lambda[1] - fq[j][1].lambda[1] - log(1 + f.c^2) / (2n)
                    lo, hi = -(thq + ths[j]), f.u / (2n) + thq + ths[j]
                    return mid < lo ? mid / lo : (mid > hi ? mid / hi : 0.0)     # violation ratio, 0 inside the sandwich
                end
                @test all(==(0.0), map(sandwich, 1:50))
            end
        end
    end
end

# Delta_2T at M_state = 2 pi (wrapped mode); statistical comparison S5
hh_mean(a) = sum(a) / length(a)
hh_var(a) = sum(abs2, a .- hh_mean(a)) / (length(a) - 1)
hh_se(a, b) = 4 * sqrt(hh_var(a) / length(a) + hh_var(b) / length(b))

@testset "AC8d two-trajectory" begin
    # ensemble: K = 2 pi, 1000 starts, n = 1e4, wrapped, d0 = d0* (T8): |mean lambda_1^QR - mean X2| <= 4 combined standard errors (S5);
    # pathwise: K = 0.5 (screened), 50 starts, n = 1e3, tau = 1, d0 = 1e-6: |lambda_1^QR - X2| <= Delta_2T(1e-6) Psi_n/n + (B_S^QR + 2 B_S^L1)/n
    seed = hh_seed(8, 4)
    eps64 = Float64(effective_eps(Float64))
    sysW = WrappedStandardMap(HH_TWO_PI)
    starts = hh_starts(StableRNG(seed), 1000)
    x0w = [torus_state(z) for z in starts]
    d0 = hh_d0star(HH_TWO_PI, HH_TWO_PI, eps64)
    q = propagate(sysW, x0w, 10_000; k = 2, tau = 1, seed = seed, reducers = (QRSpectrum(HouseholderQR()),))
    s = propagate_shadow(sysW, x0w, 10_000; d0 = d0, tau = 1, seed = seed)
    lq = [q.finals[j][1].lambda[1] for j in 1:1000]
    x2 = [s.finals[j][1].X2 for j in 1:1000]
    @test abs(hh_mean(lq) - hh_mean(x2)) <= hh_se(lq, x2)

    K = 0.5
    sysK = WrappedStandardMap(K)
    st, _ = hh_regular_starts(K, 50, StableRNG(seed + 1))
    x0k = [torus_state(z) for z in st]
    n = 1000
    q = propagate(sysK, x0k, n; k = 2, tau = 1, seed = seed, reducers = (QRSpectrum(HouseholderQR()),))
    s = propagate_shadow(sysK, x0k, n; d0 = 1e-6, tau = 1, seed = seed)
    D = hh_D2T(K, 1e-6, HH_TWO_PI, eps64)
    Psik = [lb_Psi(sysK, x0k[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], n) for j in 1:50]
    # Delta_2T Psi_n/n + (B_S^QR + 2 B_S^L1)/n (lyap-two-trajectory T5a; Amendment D35), B_S of L1 S3 as amended
    pathwise(j) = abs(q.finals[j][1].lambda[1] - s.finals[j][1].X2) <= D * Psik[j] / n + 3 * hh_BS(Float64, n, hh_Gfirst(sysK, x0k[j], 1), K) / n
    @test all(pathwise, 1:50)
end
