# AC1 to AC4: analytic cases (Y5 hyperbolic map, rotation and elliptic point, K = 0 shear, hyperbolic fixed point).
# All through propagate with the chain (Benettin(),); n X(n) = sum + comp of the final; bound B_S(n) = 8 eps (n + G_first) (S3).

# Final (n, sum, comp, lambda) of the estimator for orbit j of a propagate result.
lb_fin(res, j) = res.finals[j][1]
lb_nX(f) = (f.sum[1], f.comp[1])

# 2 x 2 BigFloat matrices as (a, b, c, d) = [[a, b], [c, d]].
lb_mm(A, B) = (A[1] * B[1] + A[2] * B[3], A[1] * B[2] + A[2] * B[4], A[3] * B[1] + A[4] * B[3], A[3] * B[2] + A[4] * B[4])
function lb_matpow(A, n::Int)
    R = (one(A[1]), zero(A[1]), zero(A[1]), one(A[1]))
    B = A
    while n > 0
        isodd(n) && (R = lb_mm(R, B))
        B = lb_mm(B, B)
        n >>= 1
    end
    return R
end

@testset "AC1 Y5 hyperbolic map" begin
    # A = [[2,1],[1,1]], start at the origin (fixed point), ids 1..100, tau in {1,5,10}, n in {20,100,1e4}.
    # (a) |n X - ln(norm(A^n w0)/norm(w0))_ref| <= B_S(n); (b) |n X - closed form (E10)| <= B_S(n); G_first = ||A^tau||_F^2.
    seed = LB_SEED + LB_OFF + 10
    ns = (20, 100, 10_000)
    for T in LB_TYPES
        epsT = lb_eps(T)
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        nid = T === Double64 ? 50 : 100                  # Double64: ids 1..50 (runtime)
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        w0s = [lb_w0(seed, j, T) for j in 1:nid]
        setprecision(BigFloat, LB_BITS) do
            A = (BigFloat(2), BigFloat(1), BigFloat(1), BigFloat(1))
            mu = (3 + sqrt(BigFloat(5))) / 2
            phi = (1 + sqrt(BigFloat(5))) / 2
            nv = sqrt(phi^2 + 1)
            vp = (phi / nv, 1 / nv)           # eigenvector of mu
            vm = (-1 / nv, phi / nv)          # eigenvector of 1/mu
            powers = Dict(n => lb_matpow(A, n) for n in ns)
            excluded = 0
            for tau in (1, 5, 10)
                G = sum(abs2, Float64.(lb_matpow(A, tau)))
                for n in ns
                    res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                    P = powers[n]
                    for j in 1:nid
                        w = (lb_big(w0s[j][1]), lb_big(w0s[j][2]))
                        c1 = w[1] * vp[1] + w[2] * vp[2]
                        c2 = w[1] * vm[1] + w[2] * vm[2]
                        if abs(c1) < 0.05
                            tau == 1 && n == ns[1] && (excluded += 1)
                            continue
                        end
                        y = (P[1] * w[1] + P[2] * w[2], P[3] * w[1] + P[4] * w[2])
                        n0 = log(w[1]^2 + w[2]^2) / 2
                        refa = log(y[1]^2 + y[2]^2) / 2 - n0
                        refb = n * log(mu) + log(abs(c1)) + log(1 + (c2 / c1)^2 * mu^(-4n)) / 2 - n0
                        S, C = lb_nX(lb_fin(res, j))
                        B = lb_BS(epsT, n, G)
                        @test lb_fin(res, j).n == n
                        @test lb_err(S, C, refa) <= B
                        @test lb_err(S, C, refb) <= B
                    end
                end
            end
            @test excluded <= nid ÷ 10
        end
    end
end

@testset "AC2 rotation and elliptic point" begin
    # (a) rot4 [[0,1],[-1,0]], w0 = (1,0): every alpha = 1 exactly, so sum = comp = 0 (tolerance 0), tau in {1,3}, n up to 9996.
    # (b) K = -2 at (0,0), M = [[-1,1],[-2,1]], M^2 = -I (E13), tau in {1,2}, ids 1..100: even n |nX| <= B_S,
    #     odd n |nX - ln(norm(M w0)/norm(w0))| <= B_S and |nX| <= ln((3 + sqrt 5)/2) + B_S.
    seed = LB_SEED + LB_OFF + 20
    for T in LB_TYPES
        epsT = lb_eps(T)
        rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
        for tau in (1, 3), n in (12, 120, 1200, 9996)
            res = lb_run(rot, [SVector{2, T}(0, 0)], n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),),
                            W0 = [SMatrix{2, 1, T}(1, 0)])
            f = lb_fin(res, 1)
            @test f.sum[1] == 0 && f.comp[1] == 0 && f.lambda[1] == 0
        end
        sys = StandardMap(exact(T, -2))
        nid = T === Double64 ? 50 : 100                  # Double64: ids 1..50 (runtime)
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        w0s = [lb_w0(seed, j, T) for j in 1:nid]
        bound_odd = log((3 + sqrt(5.0)) / 2)
        Mpow = setprecision(BigFloat, LB_BITS) do
            Mb = (BigFloat(-1), BigFloat(1), BigFloat(-2), BigFloat(1))     # tangent at (0,0), K = -2: [[1 + c, 1], [c, 1]], c = -2
            Dict(n => lb_matpow(Mb, n) for n in (1, 2, 3, 4, 10, 11, 100, 101, 1000, 1001))
        end
        for (tau, ns) in ((1, (1, 2, 3, 4, 10, 11, 100, 101, 1000, 1001)), (2, (2, 4, 10, 100, 1000)))
            G = tau == 1 ? 7.0 : 2.0                      # G_first: ||M||_F^2 = 7 for tau = 1; M^2 = -I gives 2
            for n in ns
                res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                B = lb_BS(epsT, n, G)
                for j in 1:nid
                    S, C = lb_nX(lb_fin(res, j))
                    ref = setprecision(BigFloat, LB_BITS) do
                        P = Mpow[n]
                        w1, w2 = lb_big(w0s[j][1]), lb_big(w0s[j][2])
                        y1, y2 = P[1] * w1 + P[2] * w2, P[3] * w1 + P[4] * w2
                        (log(y1^2 + y2^2) - log(w1^2 + w2^2)) / 2
                    end
                    @test lb_err(S, C, ref) <= B
                    if iseven(n)
                        @test abs(Float64(S + C)) <= B
                    else
                        @test abs(Float64(S + C)) <= bound_odd + B
                    end
                end
            end
        end
    end
end

@testset "AC3 K = 0 shear" begin
    # w_n = (a + n b, b) (E11), dyadic start so the orbit is exact: |nX - ln(norm(w_n)/norm(w_0))| <= B_S(n), G_first = 2 + tau^2;
    # nX < ln(n + 1). 20 ids for n <= 1e4 (Double64 at n = 1e4: 8), 3 ids for n = 1e6 (runtime). w0 = (0, 1): |S_n - ln(sqrt(n^2 + 1))| <= 8 n eps at n = 1e4.
    seed = LB_SEED + LB_OFF + 30
    rng = StableRNG(seed)
    dy() = rand(rng, -2048:2048) / 1024        # dyadic in [-2, 2]
    starts = [(dy(), dy()) for _ in 1:20]
    for T in LB_TYPES
        epsT = lb_eps(T)
        sys = StandardMap(exact(T, 0))
        for tau in (1, 10), n in (10, 100, 10_000, 1_000_000)
            n == 1_000_000 && T !== Float64 && continue   # Double64 sin/cos at |x| ~ 1e6 cost 10 us per step: stage 04 list
            m = n == 1_000_000 ? 3 : (T === Double64 && n == 10_000 ? 8 : 20)
            x0 = [SVector{2, T}(T(s[1]), T(s[2])) for s in starts[1:m]]
            res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
            B = lb_BS(epsT, n, 2.0 + tau^2)
            for j in 1:m
                S, C = lb_nX(lb_fin(res, j))
                w = lb_w0(seed, j, T)
                ref = setprecision(BigFloat, LB_BITS) do
                    a, b = lb_big(w[1]), lb_big(w[2])
                    (log((a + n * b)^2 + b^2) - log(a^2 + b^2)) / 2
                end
                @test lb_err(S, C, ref) <= B
                @test Float64(S + C) < log(n + 1.0)
            end
        end
        for tau in (1, 10)
            res = lb_run(sys, [SVector{2, T}(T(0.25), T(0.5))], 10_000; k = 1, tau = tau, seed = seed,
                            reducers = (Benettin(),), W0 = [SMatrix{2, 1, T}(0, 1)])
            S, C = lb_nX(lb_fin(res, 1))
            ref = setprecision(BigFloat, LB_BITS) do
                log(BigFloat(10_000)^2 + 1) / 2
            end
            @test lb_err(S, C, ref) <= 8 * 10_000 * epsT
        end
    end
end

@testset "AC4 hyperbolic fixed point" begin
    # (0,0), K in {1,2,6} (E12): w_n = c1 mu^n v_u + c2 mu^-n v_s, so n X = (1/2) ln(norm(w_n)^2 / norm(w_0)^2) exactly (the
    # O(mu^-2n) cross term is kept); asymptote n ln mu + ln(abs(c1) norm(v_u)). tau in {1,5}, ids 1..100, n in {20,100,1e3}.
    seed = LB_SEED + LB_OFF + 40
    for T in LB_TYPES
        epsT = lb_eps(T)
        x0 = [SVector{2, T}(0, 0) for _ in 1:100]
        w0s = [lb_w0(seed, j, T) for j in 1:100]
        for K in (1, 2, 6)
            sys = StandardMap(exact(T, K))
            Mf = SMatrix{2, 2, Float64}(1 + K, K, 1, 1)             # tangent at the fixed point (c = K)
            for tau in (1, 5), n in (20, 100, 1000)
                G = sum(abs2, Mf^tau)
                B = lb_BS(epsT, n, G)
                res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                setprecision(BigFloat, LB_BITS) do
                    Kb = BigFloat(K)
                    s = sqrt(Kb * (Kb + 4))
                    mu = (2 + Kb + s) / 2
                    vu = (one(Kb), (s - Kb) / 2)
                    vs = (one(Kb), -(Kb + s) / 2)
                    cross(a, b) = a[1] * b[2] - a[2] * b[1]
                    for j in 1:100
                        w = (lb_big(w0s[j][1]), lb_big(w0s[j][2]))
                        c1 = cross(w, vs) / cross(vu, vs)
                        c2 = cross(vu, w) / cross(vu, vs)
                        nrm2 = c1^2 * mu^(2n) * (vu[1]^2 + vu[2]^2) + 2 * c1 * c2 * (vu[1] * vs[1] + vu[2] * vs[2]) +
                               c2^2 * mu^(-2n) * (vs[1]^2 + vs[2]^2)
                        ref = (log(nrm2) - log(w[1]^2 + w[2]^2)) / 2
                        S, C = lb_nX(lb_fin(res, j))
                        @test lb_err(S, C, ref) <= B
                    end
                end
            end
        end
    end
end
