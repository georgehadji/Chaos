# AC1 to AC5: analytic cases (linear map, shear, hyperbolic fixed point), minimal image, identical reference orbit.
# Spec: lyap-two-trajectory-spec.md, Acceptance criteria 1 to 5. l2_nX2(n) is read as sum + comp of the final record.

l2_nX2(f) = f.sum + f.comp

@testset "AC1 linear map at the origin" begin
    # Y5 map A = [[2,1],[1,1]], start 0, shadow at d0 u_0 with u_0 = (1,0), d0 = 2^-20: |n X2 - ln||A^n u_0||| <= B_S(n); G_first = ||A^tau||_F^2
    seed = l2_seed(1)
    for T in (Float64, Double64), tau in (1, 5, 10), n in (20, 100, 1000)
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        res = propagate_shadow(sys, [zero(SVector{2, T})], n; d0 = T(2.0^-20), tau = tau, seed = seed, U0 = [SVector{2, T}(1, 0)])
        f = res.finals[1][1]
        Gfirst = Float64(l2_fibbig(2tau + 1)^2 + 2 * l2_fibbig(2tau)^2 + l2_fibbig(2tau - 1)^2)
        err = setprecision(BigFloat, 256) do
            Float64(abs(l2_big(l2_nX2(f)) - l2_ln_cat(n)))
        end
        @test f.n == n
        @test err <= l2_BS(T, n, Gfirst)
    end
end

@testset "AC2 K = 0 shear at the origin" begin
    # StandardMap K = 0 at the origin (T9): |n X2 - ln sqrt((a + n b)^2 + b^2)| <= B_S(n), G_first = 2 + tau^2; seeded u_0 = (a, b) (unit to rounding; the reference divides by its exact norm)
    seed = l2_seed(2)
    rng = StableRNG(seed)
    angles = [2pi * rand(rng) for _ in 1:20]
    for T in (Float64, Double64), d0 in (T(2.0^-20), T === Float64 ? 1e-9 : exact(Double64, "1e-9")), tau in (1, 10), n in (10, 100, 10_000)
        sys = StandardMap(T(0))
        U0 = [SVector{2, T}(cos(th), sin(th)) for th in angles]
        res = propagate_shadow(sys, [zero(SVector{2, T}) for _ in angles], n; d0 = d0, tau = tau, seed = seed, U0 = U0)
        errs = setprecision(BigFloat, 256) do
            [Float64(abs(l2_big(l2_nX2(res.finals[j][1])) - log(sqrt((l2_big(U0[j][1]) + n * l2_big(U0[j][2]))^2 + l2_big(U0[j][2])^2) /
                                                              sqrt(l2_big(U0[j][1])^2 + l2_big(U0[j][2])^2)))) for j in eachindex(angles)]
        end
        @test all(errs .<= l2_BS(T, n, 2 + tau^2))
    end
end

@testset "AC3 hyperbolic fixed point" begin
    # (0,0), K in {1, 2, 6}, d0 in {1e-6, 1e-9}, tau = 1, 100 seeded u_0, n in {20, 100, 1e3}: |n X2 - [n ln mu_+ + ln(|c_1| ||v_u||)]| <= B_S + 2 n E_nl(d0) (E12)
    seed = l2_seed(3)
    rng = StableRNG(seed)
    angles = [2pi * rand(rng) for _ in 1:100]
    U0 = [SVector(cos(th), sin(th)) for th in angles]
    eps64 = Float64(effective_eps(Float64))
    for K in (1.0, 2.0, 6.0), d0 in (1e-6, 1e-9), n in (20, 100, 1000)
        sys = StandardMap(K)
        res = propagate_shadow(sys, [zeros(SVector{2, Float64}) for _ in angles], n; d0 = d0, tau = 1, seed = seed, U0 = U0)
        Gfirst = (1 + K)^2 + 1 + K^2 + 1
        bound = l2_BS(Float64, n, Gfirst) + 2 * n * l2_Enl(K, d0)
        errs = setprecision(BigFloat, 256) do
            Kb = BigFloat(K)
            s = sqrt(Kb * (Kb + 4))
            mu = (2 + Kb + s) / 2
            vu = sqrt(1 + ((s - Kb) / 2)^2)
            [Float64(abs(l2_big(l2_nX2(res.finals[j][1])) - (n * log(mu) + log(abs(BigFloat(U0[j][1]) * (Kb + s) / 2 + BigFloat(U0[j][2])) / s * vu /
                                                                          sqrt(BigFloat(U0[j][1])^2 + BigFloat(U0[j][2])^2))))) for j in eachindex(angles)]
        end
        @test all(errs .<= bound)
    end
end

@testset "AC4 minimal image" begin
    # wrapped state, K = 0, p = 0, d0 = 2^-20: a pair straddling x = 0 reads n X2 = 0 like the same pair at x = 1 (|X2| <= 2 E_ro, results within 2 E_ro);
    # negative control as amended (Amendment D35): the exponent of the plain-difference double is not observable (Algorithm item 1 rejects the pair), so (i) the double is
    # rejected before any work on the straddling pair, (ii) it agrees with the correct hook at x = 1, (iii) at hook level the plain difference of the straddling pair is
    # ~2 pi (ln(2 pi/d0) = 15.7 >= 10) while `separation` returns the pair within pi with norm d0;
    # mirror case across p = 0 with u_0 = (0,1): |n X2 - ln sqrt(n^2 + 1)| <= B_S + 2 n E_ro (E11 with a = 0, b = 1)
    seed = l2_seed(4)
    d0, eps64 = 2.0^-20, Float64(effective_eps(Float64))
    sys = WrappedStandardMap(0.0)
    Ero = l2_Ero(0.0, d0, L2_TWO_PI, eps64)
    u0x = SVector(1.0, 0.0)
    strad = torus_state(SVector(-d0 / 2, 0.0))        # x = 2 pi - d0/2 as the system computes it; the shadow lies across x = 0
    plain = torus_state(SVector(1.0, 0.0))
    @test strad[1] > 6.28 && plain[1] == 1.0
    for n in (1, 10, 100)
        res = propagate_shadow(sys, [strad, plain], n; d0 = d0, tau = 1, seed = seed, U0 = [u0x, u0x])
        X2 = [f[1].X2 for f in res.finals]
        @test all(abs.(X2) .<= 2 * Ero)
        @test abs(X2[1] - X2[2]) <= 2 * Ero
    end
    # hook-level checks on the straddling pair after one handoff
    res = propagate_shadow(sys, [strad], 1; d0 = d0, tau = 1, seed = seed, U0 = [u0x])
    xref, xsh = res.x[1], res.W[1]
    pl = SVector(xsh[1] - xref[1], xsh[2] - xref[2])
    sep = separation(sys, xref, xsh)
    @test log(sqrt(pl[1]^2 + pl[2]^2) / d0) >= 10                                              # (iii) the artifact of T3 itself
    @test all(abs.(sep) .<= pi) && abs(sqrt(sep[1]^2 + sep[2]^2) - d0) <= 4 * Ero * d0
    @test maximum(abs, separation(sys, xref, displaced(sys, xref, sep)) - sep) <= 4 * eps64 * L2_TWO_PI      # displaced and separation are inverse
    double = l2_Plain(sys)
    err = l2_health_error(() -> propagate_shadow(double, [strad], 1; d0 = d0, tau = 1, seed = seed, U0 = [u0x]))     # (i) rejected before any work
    @test err isa NumericalHealthError && err.what == "separation" && err.step == 0 && err.orbit == 1 && err.backend == :cpu
    ok = propagate_shadow(double, [plain], 10; d0 = d0, tau = 1, seed = seed, U0 = [u0x])     # (ii) the double is fine away from the boundary ...
    good = propagate_shadow(sys, [plain], 10; d0 = d0, tau = 1, seed = seed, U0 = [u0x])
    @test abs(ok.finals[1][1].X2 - good.finals[1][1].X2) <= 2 * Ero                           # ... and agrees with the correct hook
    # mirror case: p = 2 pi - d0/2, shadow across p = 0
    mirror = torus_state(SVector(1.0, -d0 / 2))
    u0p = SVector(0.0, 1.0)
    for n in (1, 10, 100)
        f = propagate_shadow(sys, [mirror], n; d0 = d0, tau = 1, seed = seed, U0 = [u0p]).finals[1][1]
        err = setprecision(BigFloat, 256) do
            Float64(abs(l2_big(l2_nX2(f)) - log(sqrt(BigFloat(n)^2 + 1))))
        end
        @test err <= l2_BS(Float64, n, 3) + 2 * n * Ero
    end
end

@testset "AC5 identical reference orbit" begin
    # reference orbit bitwise equal to the plain forward iteration and to the stream orbit of an L1 run (Benettin renormalises, so W cannot overflow), K in {0.5, 2 pi}, tau in {1, 4},
    # both modes, handoff counts N in {4, 8, 40, 400, 1000} (prefix invariance); Float64 100 starts, Double64 20 starts
    seed = l2_seed(5)
    for T in (Float64, Double64), mode in (:lifted, :wrapped), K in (0.5, L2_TWO_PI), tau in (1, 4)
        sys = l2_sys(mode, T(K))
        m = T === Float64 ? 100 : 20
        x0s = [l2_state(mode, T, z) for z in l2_starts(StableRNG(seed), m)]
        for N in (4, 8, 40, 400, 1000)
            plain = map(x -> (z = x; for _ in 1:N; z = forward(sys, z); end; z), x0s)
            shadow = propagate_shadow(sys, x0s, N; d0 = T(1e-8), tau = tau, seed = seed)
            l1 = propagate(sys, x0s, N; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
            @test isequal(shadow.x, plain)
            @test isequal(l1.x, plain)
        end
    end
end
