# AC1 to AC4: exactness, renormalization algebra, K = 0 standard map, BigFloat tangent product.

@testset "AC1 Y5 exactness" begin
    # W_n = A^n (Fibonacci entries, E9) and x_n = A^n x0, bitwise, n = 1..35, k = 2, tau = 1.
    seed = TC_SEED + 1
    sys = tc_Lin(SMatrix{2,2,Float64}(2, 1, 1, 1))
    x0 = SVector(3.0, 5.0)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    for n in 1:35
        res = propagate(sys, [x0], n; k = 2, tau = 1, seed = seed, W0 = [W0], reducers = (tc_Reader(1),))
        An = tc_Apow(n)
        @test res.W[1] == An
        @test res.x[1] == An * x0
        @test res.finals[1][1] == (count = n, last = n, bad = 0, maxabs = maximum(abs, An))
    end
end

@testset "AC2 Y5 renormalization and exponent" begin
    # S_n + ln||w_n|| = ln||A^n w0|| (E7, I3) to 8 n eps; S_n / n -> ln((3 + sqrt 5)/2) (E9). x0 = 0: x_n = A^n x0 stays 0.
    seed = TC_SEED + 2
    N = 10_000
    sys = tc_Lin(SMatrix{2,2,Float64}(2, 1, 1, 1))
    x0 = [SVector(0.0, 0.0) for _ in 1:100]
    W0 = [initial_deviations(seed, j, Val(2), Val(1), Float64) for j in 1:100]
    ref = setprecision(BigFloat, 256) do
        AN = tc_matpow(BigFloat[2 1; 1 1], N)
        [Float64(log(sqrt(sum(abs2, AN * tc_bigmat(w))))) for w in W0]
    end
    lam = log((3 + sqrt(5)) / 2)
    for tau in (1, 5, 10)
        res = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, W0 = W0, reducers = (tc_LogNorm(), tc_Reader(tau)))
        worst_a = 0.0
        worst_b = 0.0
        for j in 1:100
            S = res.finals[j][1].logs[1]
            worst_a = max(worst_a, abs(S + log(norm(res.W[j])) - ref[j]))
            worst_b = max(worst_b, abs(S / N - lam))
        end
        @test worst_a <= 8 * N * eps(Float64)
        @test worst_b <= 2e-3
        # I2: handoffs exactly at multiples of tau, N/tau of them
        @test all(f -> f[2].count == N ÷ tau && f[2].last == N && f[2].bad == 0, res.finals)
    end
    # Chain order: a later reducer sees the matrix replaced by an earlier one (E7); reading first sees the unreplaced matrix.
    tau = 5
    after = propagate(sys, x0[1:10], N; k = 1, tau = tau, seed = seed, W0 = W0[1:10], reducers = (tc_LogNorm(), tc_Reader(tau)))
    before = propagate(sys, x0[1:10], N; k = 1, tau = tau, seed = seed, W0 = W0[1:10], reducers = (tc_Reader(tau), tc_LogNorm()))
    @test all(f -> f[2].maxabs <= 1 + 4 * eps(Float64), after.finals)
    @test all(f -> f[1].maxabs > 1.5, before.finals)
end

@testset "AC3 K = 0 standard map" begin
    # M = [[1,1],[0,1]] (E8, K = 0): W_n = [[1,n],[0,1]], x_n = x0 + n p0 bitwise; k = 1 from e_2: S_n = 0.5 ln(n^2 + 1).
    seed = TC_SEED + 3
    sys = tc_Std(0.0)
    x0 = SVector(0.25, 0.5)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    for n in (1, 10, 1000, 10^6)
        res = propagate(sys, [x0], n; k = 2, tau = 1, seed = seed, W0 = [W0])
        @test res.W[1] == SMatrix{2,2,Float64}(1, 0, n, 1)
        @test res.x[1] == SVector(0.25 + 0.5 * n, 0.5)   # lifted coordinates, never wrapped (I7)
    end
    n = 10_000
    res = propagate(sys, [x0], n; k = 1, tau = 1, seed = seed, W0 = [SMatrix{2,1,Float64}(0, 1)], reducers = (tc_LogNorm(),))
    S = res.finals[1][1].logs[1]
    @test abs(S - 0.5 * log(float(n)^2 + 1)) <= 8 * n * eps(Float64)
    @test S / n <= log(n) / n + 1e-12
end

@testset "AC4 BigFloat tangent product" begin
    # ||W_n - W_ref||_F <= 8 n eps(T) prod||M_i||_F ||W0||_F, W_ref at BigFloat(256) along the realized orbit points.
    seed = TC_SEED + 4
    for T in (Float64, Double64), (K, n) in ((0.5, 100), (2pi, 20))
        sys = tc_Std(T(K))
        rng = StableRNG(seed + (T === Float64 ? 0 : 1000) + n)
        x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:20]
        W0 = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:20]
        res = propagate(sys, x0, n; k = 2, tau = 1, seed = seed, W0 = W0)
        for j in 1:20
            pts = Vector{SVector{2,T}}(undef, n)
            z = x0[j]
            for i in 1:n
                pts[i] = z
                z = forward(sys, z)
            end
            @test isequal(res.x[j], z)   # the kernel realizes the orbit of forward
            err, bound = tc_std_residual(sys.K, pts, W0[j], res.W[j], eps(T))
            @test err <= bound
        end
    end
end
