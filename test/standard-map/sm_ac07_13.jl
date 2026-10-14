# Known answers: criteria 7 to 13 (items 6 to 9 of the spec Equations).

# AC7 part 1: K = 0, grid data, x_n = x_0 + n p_0 and p_n = p_0 with zero difference for 1000 steps.
function sm_ac7_orbit_defects(::Type{T}, rng, m::Integer) where {T}
    sys = StandardMap(zero(T))
    bad = 0
    for _ in 1:m
        x0 = T(rand(rng, -8191:8191) / 1024)      # multiples of 2^-10, |value| < 8
        p0 = T(rand(rng, -8191:8191) / 1024)
        z = SVector(x0, p0)
        for k in 1:1000
            z = forward(sys, z)
            bad += (z[1] == x0 + k * p0) ? 0 : 1
            bad += (z[2] == p0) ? 0 : 1
        end
    end
    return bad
end

# AC7 part 2: the product of n one-step tangents equals [[1, n], [0, 1]] exactly, n in {1, 2, 10, 100, 1000, 10^4}.
function sm_ac7_product_defects(::Type{T}, rng, m::Integer) where {T}
    sys = StandardMap(zero(T))
    ns = (1, 2, 10, 100, 1000, 10_000)
    bad = 0
    for _ in 1:m
        z = SVector(T(rand(rng, -8191:8191) / 1024), T(rand(rng, -8191:8191) / 1024))
        P = SMatrix{2, 2, T}(one(T), zero(T), zero(T), one(T))
        for k in 1:last(ns)
            P = tangent(sys, z) * P
            z = forward(sys, z)
            if k in ns
                bad += (P == SMatrix{2, 2, T}(one(T), zero(T), T(k), one(T))) ? 0 : 1
            end
        end
    end
    return bad
end

# AC8: X(n) = log(|P w0| / |w0|) / n from the exact product P = [[1, n], [0, 1]], recorded at n in ns.
function sm_ac8_products(::Type{T}, ns) where {T}
    sys = StandardMap(zero(T))
    z = SVector(T(0.25), T(0.125))             # multiples of 2^-10, so the orbit stays exact
    P = SMatrix{2, 2, T}(one(T), zero(T), zero(T), one(T))
    out = Dict{Int, SMatrix{2, 2, T, 4}}()
    DT = tangent(sys, z)
    for k in 1:maximum(ns)
        if k <= 10_000                          # follow the orbit and take the tangent at each point
            DT = tangent(sys, z)
            z = forward(sys, z)
        end                                     # beyond 10^4 reuse DT: at K = 0 it is the same matrix at every point (AC7)
        P = DT * P
        k in ns && (out[k] = P)
    end
    return out
end

sm_expo(P, w0, n) = (w = P * w0; log(sqrt(w[1]^2 + w[2]^2) / sqrt(w0[1]^2 + w0[2]^2)) / n)

# AC9: fixed point class, trace and eigenvalue data at (0,0) and (pi,0) against the exact trace 2 + K cos x.
function sm_ac9_ratios(::Type{T}) where {T}
    worst_tr = 0.0
    worst_ev = 0.0
    wrong_class = 0
    piT = exact(T, pi)
    for K in (-6.0, -2.0, 1.0, 3.0, 4.5, 6.0), (x, sgn) in ((zero(T), 1), (piT, -1))
        sys = StandardMap(T(K))
        DT = tangent(sys, SVector(x, zero(T)))
        tr = DT[1, 1] + DT[2, 2]
        tol = 8 * sm_eps(T) * (1 + abs(K))
        setprecision(BigFloat, SM_BITS) do
            tr_exact = 2 + sgn * BigFloat(K)
            kind_exact, mag_exact = sm_class(tr_exact)
            kind, mag = sm_class(tr)
            wrong_class += (kind == kind_exact) ? 0 : 1
            worst_tr = max(worst_tr, Float64(abs(sm_big(tr) - tr_exact)) / tol)
            worst_ev = max(worst_ev, Float64(abs(sm_big(mag) - mag_exact)) / (4 * tol))
        end
    end
    return (worst_tr, worst_ev, wrong_class)
end

# AC9 Jordan blocks: max entry of DT^n over n is 4 (K = 4 at (pi,0)) and 1 (K = 0), within 1e-12.
function sm_ac9_jordan_defects(::Type{T}) where {T}
    piT = exact(T, pi)
    worst = 0.0
    for (K, x, ratio) in ((4.0, piT, 4.0), (0.0, zero(T), 1.0))
        DT = tangent(StandardMap(T(K)), SVector(x, zero(T)))
        P = DT
        for n in 2:1000
            P = DT * P
            n in (10, 100, 1000) && (worst = max(worst, abs(Float64(maximum(abs, P)) / n - ratio)))
        end
    end
    return worst
end

@testset "AC7 K=0 exactness" begin
    # x_n = x_0 + n p_0 exactly (0 ulp) over 1000 steps (100 orbits Float64, 20 Double64; D38: 1000 / 200 in verify); product of n tangents = [[1,n],[0,1]] exactly; grid data multiples of 2^-10
    @test sm_ac7_orbit_defects(Float64, StableRNG(SM_SEED + 7001), 100) == 0
    @test sm_ac7_orbit_defects(Double64, StableRNG(SM_SEED + 7002), 20) == 0
    # the tangent at K = 0 does not depend on the point, so 20 starting points suffice for the 10^4-step products
    @test sm_ac7_product_defects(Float64, StableRNG(SM_SEED + 7003), 20) == 0
    @test sm_ac7_product_defects(Double64, StableRNG(SM_SEED + 7004), 20) == 0
end

@testset "AC8 K=0 finite-time exponent" begin
    # X(n) = ln(sqrt(n^2+1))/n within 8 eps for w0=(0,1); n X/ln n - 1 <= 1e-9 at n=1e4; X = 0 exactly for w0=(1,0)
    ns = (1, 10, 100, 1000, 10_000, 1_000_000)
    for T in (Float64, Double64)
        prods = sm_ac8_products(T, ns)
        for n in ns
            P = prods[n]
            X = sm_expo(P, SVector(zero(T), one(T)), n)
            setprecision(BigFloat, SM_BITS) do
                ref = log(sqrt(BigFloat(n)^2 + 1)) / n
                @test Float64(abs(sm_big(X) - ref)) <= 8 * sm_eps(T)
            end
            @test sm_expo(P, SVector(one(T), zero(T)), n) == 0
        end
        X4 = sm_expo(prods[10_000], SVector(zero(T), one(T)), 10_000)
        @test abs(Float64(10_000 * X4 / log(T(10_000))) - 1) <= 1.0e-9
    end
end

@testset "AC9 fixed points and eigenvalues" begin
    # tr = 2 + K cos x to 8 eps (1+|K|); |lambda_+| or rotation angle to 32 eps (1+|K|); K=4 Jordan max entry 4n, K=0 gives n
    for T in (Float64, Double64)
        worst_tr, worst_ev, wrong_class = sm_ac9_ratios(T)
        @test wrong_class == 0
        @test worst_tr <= 1
        @test worst_ev <= 1
        @test sm_ac9_jordan_defects(T) <= 1.0e-12
    end
end

@testset "AC10 fixed point is exact" begin
    # (0,0) is mapped to itself: p' = 0 + K sin 0 = 0, x' = 0; 10^4 steps, every listed K, zero difference
    for T in (Float64, Double64), K in SM_KS
        sys = StandardMap(T(K))
        z = SVector(zero(T), zero(T))
        ok = true
        for _ in 1:10_000
            z = forward(sys, z)
            ok &= (z[1] == 0) && (z[2] == 0)
        end
        @test ok
    end
end

@testset "AC11 benchmark trace row" begin
    # unit-torus (1/2, 0) is (pi, 0): tr = 2 - K there and 2 + K at (0,0), to 8 eps (1+|K|); elliptic for K<4, hyperbolic for K>4
    # (benchmarks row "Stability of the central fixed point at (1/2, 0)"; conversion x = 2 pi x_u, spec item 6)
    for T in (Float64, Double64), Ks in ("1", "3", "3.9", "4.1", "5")
        K = exact(T, Ks)
        sys = StandardMap(K)
        tol = 8 * sm_eps(T) * (1 + abs(Float64(K)))
        xc = 2 * exact(T, pi) * exact(T, 1 // 2)
        DTc = tangent(sys, SVector(xc, zero(T)))
        DT0 = tangent(sys, SVector(zero(T), zero(T)))
        trc = DTc[1, 1] + DTc[2, 2]
        @test sm_absdiff((trc,), (2 - K,)) <= tol
        @test sm_absdiff((DT0[1, 1] + DT0[2, 2],), (2 + K,)) <= tol
        @test sm_class(trc)[1] == (Float64(K) < 4 ? :elliptic : :hyperbolic)
    end
end

@testset "AC12 Datseris eigenvectors at K=1" begin
    # DT v = lambda v gives v_p/v_x = (lambda - a)/b; derived (s-K)/2 = 0.618 and -(K+s)/2 = -1.618 to 1e-12, printed 0.62 within 0.005
    K = 1.0
    s = sqrt(K * (K + 4))
    DT = tangent(StandardMap(K), SVector(0.0, 0.0))
    tr = DT[1, 1] + DT[2, 2]
    lp = (tr + sqrt(tr^2 - 4)) / 2
    lm = (tr - sqrt(tr^2 - 4)) / 2
    ru = (lp - DT[1, 1]) / DT[1, 2]
    rs = (lm - DT[1, 1]) / DT[1, 2]
    @test abs(ru - (s - K) / 2) <= 1.0e-12
    @test abs(rs + (K + s) / 2) <= 1.0e-12
    @test maximum(abs, DT * SVector(1.0, ru) - lp * SVector(1.0, ru)) <= 1.0e-12
    @test maximum(abs, DT * SVector(1.0, rs) - lm * SVector(1.0, rs)) <= 1.0e-12
    @test abs(0.62 - abs(ru)) <= 0.005
    @test abs(0.62 - 1 / abs(rs)) <= 0.005
end

@testset "AC13 Contopoulos orbit data" begin
    # K=5, benchmarks row "Island centre and unstable 4/9 orbit point"; x = 2 pi x_u, p = 2 pi y (spec item 6); tolerance 2 pi 1e-4
    sys = StandardMap(5.0)
    tp = 2 * exact(Float64, pi)
    tol = tp * 1.0e-4
    O1 = SVector(tp * 0.68002, tp * 0.36002)
    P1 = SVector(tp * 0.642118, tp * 0.369344)
    zs = sm_orbit(sys, O1, 2)
    @test sm_tdist(zs[3], O1) <= tol
    zs = sm_orbit(sys, P1, 18)
    @test sm_tdist(zs[19], P1) <= tol
    @test sm_tdist(zs[10], P1) > tol
end
