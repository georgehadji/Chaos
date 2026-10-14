# AC4 to AC6: exact and analytic cases (K = 0, cat map, higher-dimensional Y5 maps).

@testset "AC4 exactness at K = 0" begin
    # (a) shear [[1,1],[0,1]], W_0 = I: Q = I and S_1 = S_2 = 0 bitwise at every handoff (n = 1e4); (b) standard map K = 0, k = 1 from e_2:
    # S_1(n) = (1/2) ln(n^2 + 1) to 8 n eps, S_1(n)/n <= ln(n)/n + 1e-12 (E12)
    seed = hh_seed(4)
    for T in (Float64, Double64), alg in hh_algs(T)
        eps_T = Float64(effective_eps(T))
        W0 = SMatrix{2, 2, T}(1, 0, 0, 1)
        shear = LinearSymplecticMap(SMatrix{2, 2, T}(1, 0, 1, 1))
        res = propagate(shear, [zero(SVector{2, T})], 10_000; k = 2, tau = 1, seed = seed, W0 = [W0], reducers = (QRSpectrum(alg), hh_Dev()))
        f, dv = res.finals[1]
        @test all(iszero, f.sum) && all(iszero, f.comp) && all(iszero, f.lambda)
        @test dv.maxdev == 0 && dv.count == 10_000
        @test res.W[1] == W0

        sys = StandardMap(T(0))
        x0 = SVector{2, T}(T(0.25), T(0.5))
        e2 = SMatrix{2, 1, T}(0, 1)
        for n in (10, 1000, 10_000)
            f = propagate(sys, [x0], n; k = 1, tau = 1, seed = seed, W0 = [e2], reducers = (QRSpectrum(alg),)).finals[1][1]
            S = f.sum[1] + f.comp[1]
            err = setprecision(BigFloat, 256) do
                Float64(abs(hh_big(S) - log(BigFloat(n)^2 + 1) / 2))
            end
            @test err <= 8 * n * eps_T
            n == 10_000 && @test Float64(S) / n <= log(n) / n + 1e-12
        end
    end
end

hh_fibbig(m) = (a = big(0); b = big(1); for _ in 1:m; a, b = b, a + b; end; a)

@testset "AC5 cat map finite-n identity" begin
    # A = [[2,1],[1,1]] at the fixed point 0, Q_0 = I: (a) |S_1 - ln||A^n e_1||| <= B_S; (b) |S_1 + S_2| <= n theta_sum (G_max = mu^2tau, M_max = 2);
    # (c) |lambda_1 - ln mu - ln(c)/n| <= B_S/n with c = (1 + ((sqrt 5 - 1)/2)^2)^(-1/2) (E13); (d) |lambda_1 + lambda_2| <= theta_sum
    seed = hh_seed(5)
    # E13 keeps the O(mu^-4n) term of E10 here: S_1(n) = n ln mu + ln c + (1/2) ln(1 + ((1 - c^2)/c^2) mu^-4n); at n = 10 it is 4e-18, far above B_S in Double64
    lnmu, lnc, r0, q = setprecision(BigFloat, 256) do
        mub = (3 + sqrt(BigFloat(5))) / 2
        cc = 1 / sqrt(1 + ((sqrt(BigFloat(5)) - 1) / 2)^2)
        (log(mub), log(cc), (1 - cc^2) / cc^2, mub^-4)
    end
    mu = Float64((3 + sqrt(5)) / 2)
    ref1 = Dict{Int, BigFloat}()
    for n in (10, 100, 1000, 10_000)
        ref1[n] = setprecision(BigFloat, 256) do
            log(sqrt(BigFloat(hh_fibbig(2n + 1))^2 + BigFloat(hh_fibbig(2n))^2))
        end
    end
    for T in (Float64, Double64)
        eps_T = Float64(effective_eps(T))
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        x0 = zero(SVector{2, T})
        W0 = SMatrix{2, 2, T}(1, 0, 0, 1)
        for tau in (1, 2, 5), n in (10, 100, 1000, 10_000), alg in hh_algs(T)
            f = propagate(sys, [x0], n; k = 2, tau = tau, seed = seed, W0 = [W0], reducers = (QRSpectrum(alg),)).finals[1][1]
            S1, S2 = f.sum[1] + f.comp[1], f.sum[2] + f.comp[2]
            Gfirst = Float64(hh_fibbig(2tau + 1)^2 + 2 * hh_fibbig(2tau)^2 + hh_fibbig(2tau - 1)^2)
            BS = hh_BS(T, n, Gfirst)
            theta = hh_theta(T, mu^(2tau), tau, 2.0)
            errs = setprecision(BigFloat, 256) do
                (Float64(abs(hh_big(S1) - ref1[n])), Float64(abs(hh_big(f.lambda[1]) - (lnmu + (lnc + log(1 + r0 * q^n) / 2) / n))))
            end
            @test errs[1] <= BS
            @test Float64(abs(S1 + S2)) <= n * theta
            @test errs[2] <= BS / n
            @test Float64(abs(f.lambda[1] + f.lambda[2])) <= theta
        end
    end
end

# 2N-dimensional hyperbolic Y5 case (function barrier on d = 2N); returns nothing, tests inside
function hh_y5(::Type{T}, mus, nonnormal::Bool, ::Val{d}, seed) where {T, d}
    N = d ÷ 2
    A0 = block_diagonal(map(m -> hyperbolic_block(T, m), mus)...)
    A = nonnormal ? conjugated(A0, N == 2 ? SMatrix{2, 2, T}(1, 2, 2, 1) : SMatrix{3, 3, T}(1, 2, 0, 2, 1, 1, 0, 1, 1)) : A0
    sys = LinearSymplecticMap(A)
    Q0 = initial_deviations(seed, 1, Val(d), Val(d), T)
    n = 500
    Af = Float64[Float64(A[i, j]) for i in 1:d, j in 1:d]
    sv = svdvals(Af)
    theta = hh_theta(T, sv[1] / sv[end], 1, maximum(abs, Af))
    lnvol, Lam, Bterm = setprecision(BigFloat, 2048) do
        Y = hh_matpow(hh_bigmat(A), n) * hh_bigmat(Q0)
        lv = [log(hh_det(Y[:, 1:i]' * Y[:, 1:i])) / 2 for i in 1:d]
        # exact exponents +-mu_j, eigenvectors (e_{q_j} +- e_{p_j})/sqrt 2, ordered by decreasing exponent (+0.45, +0.2, -0.2, -0.45)
        lst = Tuple{BigFloat, Int, Int}[]
        for (j, m) in enumerate(mus)
            push!(lst, (BigFloat(m), j, 1))
            push!(lst, (-BigFloat(m), j, -1))
        end
        sort!(lst; by = t -> -t[1])
        V = zeros(BigFloat, d, d)
        for (c, (_, j, sg)) in enumerate(lst)
            V[j, c] = 1 / sqrt(BigFloat(2))
            V[N + j, c] = sg / sqrt(BigFloat(2))
        end
        Q0b = hh_bigmat(Q0)
        Lam = cumsum([t[1] for t in lst])
        B = [log(abs(hh_det(V[:, 1:i]' * Q0b[:, 1:i]))) for i in 1:d]
        (lv, Lam, B)
    end
    for alg in hh_algs(T)
        f = propagate(sys, [zero(SVector{d, T})], n; k = d, tau = 1, seed = seed, W0 = [Q0], reducers = (QRSpectrum(alg),)).finals[1][1]
        setprecision(BigFloat, 2048) do
            C = cumsum([hh_big(f.sum[j]) + hh_big(f.comp[j]) for j in 1:d])
            for i in 1:d
                @test Float64(abs(C[i] - lnvol[i])) <= i * n * theta                                     # (a)
                nonnormal || @test Float64(abs(C[i] / n - Lam[i] - Bterm[i] / n)) <= i * theta           # (b)
            end
        end
        @test abs(sum(f.lambda)) <= theta                                                                # (c)
    end
    return nothing
end

@testset "AC6 Y5 higher dimension" begin
    # block-diagonal hyperbolic A, N = 2 mu = (0.45, 0.2) and N = 3 mu = (0.45, 0.3, 0.1), plus the non-normal A_S; k = d, tau = 1, seeded Q_0, n = 500, start 0:
    # (a) |C_i - ln vol_i| <= i n theta_sum (2048-bit Gram determinants); (b) A normal: |C_i/n - Lambda_i - B_i/n| <= i theta_sum; (c) |sum lambda| <= theta_sum
    seed = hh_seed(6)
    for (it, T) in enumerate((Float64, Double64)), (im, mus) in enumerate(((0.45, 0.2), (0.45, 0.3, 0.1))), (inn, nn) in enumerate((false, true))
        hh_y5(T, mus, nn, Val(2 * length(mus)), seed + 10 * it + 2 * im + inn)
    end
end
