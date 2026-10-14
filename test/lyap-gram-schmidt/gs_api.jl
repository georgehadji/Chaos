# API units of the Gram-Schmidt contract (api-batch-b.md, "Shared helpers"): orthonormalize, window_G(w), thresholds,
# effective_eps, safe_norm, orth_loss. They carry no criterion number; the plan lists them as API1 to API3.

@testset "API1 orthonormalize contract" begin
    # W = Q R, Q'Q = I, R upper with exact zeros below the diagonal and R_ii > 0, status QR_OK, same R for all variants (G3)
    @test all(alg -> alg isa Orthonormalizer, gs_algs())
    rng = StableRNG(gs_seed(0))
    for (T, shapes) in ((Float64, ((2, 2), (4, 3))),), (d, k) in shapes, _ in 1:4
        # well conditioned: Gaussian plus a dominant diagonal
        W = SMatrix{d,k,T}(ntuple(l -> T(randn(rng)) + ((l - 1) % d == (l - 1) ÷ d ? T(4) : T(0)), Val(d * k)))
        fs = map(alg -> orthonormalize(alg, W), gs_algs())
        for (alg, f) in zip(gs_algs(), fs)
            Q, R, count, status, col = f
            @test status == QR_OK && col == 0
            @test Q isa SMatrix{d,k,T} && R isa SMatrix{k,k,T}
            @test all(R[i, j] === zero(T) for j in 1:k for i in (j + 1):k)   # exact zeros below the diagonal
            @test all(R[i, i] > 0 for i in 1:k)
            @test maximum(abs, Q' * Q - I) <= 8 * k * gs_eps(T) * (alg isa RGS ? 1 : 4)
            @test maximum(abs, Q * R - W) <= 16 * gs_eps(T) * maximum(abs, W) * k
            @test alg isa RGS ? (1 <= count <= 4) : count == 1
        end
        R1 = fs[1][2]                                   # uniqueness of R (G3): every variant gives the same R
        @test all(f -> maximum(abs, f[2] - R1) <= 1000 * gs_eps(T) * maximum(abs, R1), fs)
    end
    # zero column and non-finite column: QR_PIVOT with the column index
    for alg in gs_algs()
        f = orthonormalize(alg, SMatrix{4,3,Float64}(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0))             # columns e1, e2, 0
        @test f[4] == QR_PIVOT && f[5] == 3
        g = orthonormalize(alg, SMatrix{4,3,Float64}(1, 0, 0, 0, 0, 1, 0, 0, NaN, 0, 0, 0))
        @test g[4] == QR_PIVOT && g[5] == 3
    end
    # no allocation
    W = @SMatrix [4.0 1 0.5 0.2; 0.3 3 1 0.1; 0.2 0.4 5 1; 0.1 0.2 0.3 2]
    for alg in gs_algs()
        @test gs_alloc_orthonormalize(alg, W) == 0
    end
end

@testset "API2 window conditioning and thresholds" begin
    # G10: G_w = ||R||_F ||R^-1||_F; for d = k = 2 G_w = kappa + 1/kappa and window_G = kappa; k = 1: G = 1
    rng = StableRNG(gs_seed(0) + 1)
    for T in gs_types(), _ in 1:(T === Float64 ? 12 : 2)
        a, b, d = exp(T(4 * rand(rng))), T(randn(rng)), exp(-T(4 * rand(rng)))
        R = SMatrix{2,2,T}(a, 0, b, d)
        ref = setprecision(BigFloat, 256) do      # R = [a b; 0 d]: R^-1 = [1/a, -b/(a d); 0, 1/d]; sigma_1 +- sigma_2 = sqrt(F +- 2 a d)
            ab, bb, db = gs_big(a), gs_big(b), gs_big(d)
            F = ab^2 + bb^2 + db^2
            Fi = 1 / ab^2 + (bb / (ab * db))^2 + 1 / db^2
            s1, s2 = sqrt(F + 2 * ab * db), sqrt(F - 2 * ab * db)
            (sqrt(F) * sqrt(Fi), (s1 + s2) / (s1 - s2))
        end
        @test gs_absdiff(window_Gw(R), ref[1]) <= 16 * gs_eps(T) * Float64(ref[1])
        @test gs_absdiff(window_G(R), ref[2]) <= 64 * gs_eps(T) * Float64(ref[2])
        @test window_Gw(R) >= window_G(R)
        @test window_G(SMatrix{1,1,T}(T(3))) == 1
    end
    R3 = SMatrix{3,3,Float64}(2.0, 0, 0, 0.5, 1.0, 0, 0.1, 0.2, 0.05)   # k >= 3: G_w is an upper bound of the condition number
    @test window_G(R3) >= cond(Matrix(R3)) * (1 - 1e-12)
    # tol x eps_eff x k x g, g = G_w for CGS and MGS, 1 otherwise
    for T in gs_types()
        e = gs_eps(T)
        thr(alg, tol, k, Gw) = Float64(orthogonality_threshold(alg, tol, T, k, Gw))
        @test thr(RGS(), 8, 4, 123.0) ≈ 8 * e * 4 rtol = 1e-12
        @test thr(CGS(), 8, 4, 123.0) ≈ 8 * e * 4 * 123.0 rtol = 1e-12
        @test thr(MGS(), 8.0, 2, 50.0) ≈ 8 * e * 2 * 50.0 rtol = 1e-12
        @test thr(HouseholderQR(), 8, 3, 1e6) ≈ 8 * e * 3 rtol = 1e-12
    end
    @test effective_eps(Float64) == eps(Float64)
    @test Float64(effective_eps(Double64)) == 2.0^-100
end

@testset "API3 safe_norm and orth_loss" begin
    # L1 AC8: relative error <= 2 eps over the exponent range; zero and non-finite propagate; orth_loss equals eps_orth (D4)
    rng = StableRNG(gs_seed(0) + 2)
    for e in (-1000, -500, 0, 500, 1000)
        v = SVector{3,Float64}(ldexp(randn(rng), e), ldexp(randn(rng), e), ldexp(randn(rng), e))
        ref = setprecision(BigFloat, 256) do
            sqrt(sum(abs2, BigFloat.(v)))
        end
        @test gs_absdiff(safe_norm(v), ref) <= 2 * eps(Float64) * Float64(ref)
    end
    @test safe_norm(SVector(0.0, 0.0)) == 0
    @test isnan(safe_norm(SVector(NaN, 1.0)))
    @test !isfinite(safe_norm(SVector(Inf, 1.0)))
    Q = SMatrix{3,2,Float64}(ntuple(_ -> randn(rng), Val(6)))
    @test Chaoslab.Lyapunov.orth_loss(Q) == eps_orth(Q)
end
