# Known answers: criteria 7 to 12 (K = 0 shear, decoupled limit, symmetries, fixed points and exact origin spectra).

# [[I, nI], [0, I]] in element type T, dense
function cm_shear(::Type{T}, d::Integer, n::Integer) where {T}
    h = d ÷ 2
    S = Matrix{T}(I, d, d)
    for j in 1:h
        S[j, h + j] = T(n)
    end
    return S
end

# dyadic data: multiples of 2^-10 with magnitude < 8
cm_dyadic(rng, ::Type{T}, d::Integer) where {T} = [T(rand(rng, -8191:8191) / 1024) for _ in 1:d]

# ---- AC7: K = 0 exactness ----
function cm_ac7_orbit_defects(::Type{T}, c, rng, m) where {T}
    sys = cm_sys(T, c)
    d, h = cm_d(c), cm_h(c)
    bad = 0
    for _ in 1:m
        v = cm_dyadic(rng, T, d)
        z = cm_state(T, c, v)
        for k in 1:1000
            z = forward(sys, z)
            for j in 1:h
                bad += (z[j] == v[j] + k * v[h + j]) ? 0 : 1
                bad += (z[h + j] == v[h + j]) ? 0 : 1
            end
        end
    end
    return bad
end

function cm_ac7_product_defects(::Type{T}, c, rng) where {T}
    sys = cm_sys(T, c)
    d = cm_d(c)
    ns = (1, 2, 10, 100, 1000, 10_000)
    z = cm_state(T, c, cm_dyadic(rng, T, d))
    P = Matrix{T}(I, d, d)
    bad = 0
    for k in 1:last(ns)
        P = Matrix(tangent(sys, z)) * P
        z = forward(sys, z)
        k in ns && (bad += (P == cm_shear(T, d, k)) ? 0 : 1)
    end
    return bad
end

@testset "AC7 K=0 exactness" begin
    # q_j = q_j(0) + n p_j(0), p_j = p_j(0) with zero difference for 1000 steps (Y3; Y4 open N = 3 and ring N = 100); product of n tangents = [[I, nI], [0, I]] exactly;
    # orbits: 20 Float64, 4 Double64 per system (200 / 1000 in verify); products: one start per system (the K = 0 tangent is state independent)
    for c in (cm_y3(0, 0, 0), cm_y4(3, 0, 0), cm_y4(100, 0, 0; ring = true), cm_y4(2, 0, 0))
        @test cm_ac7_orbit_defects(Float64, c, StableRNG(CM_SEED + 7001), 20) == 0
        @test cm_ac7_orbit_defects(Double64, c, StableRNG(CM_SEED + 7002), 4) == 0
    end
    for c in (cm_y3(0, 0, 0), cm_y4(3, 0, 0), cm_y4(5, 0, 0; ring = true))
        @test cm_ac7_product_defects(Float64, c, StableRNG(CM_SEED + 7003)) == 0
    end
    for c in (cm_y3(0, 0, 0), cm_y4(3, 0, 0))
        @test cm_ac7_product_defects(Double64, c, StableRNG(CM_SEED + 7004)) == 0
    end
end

# ---- AC8: K = 0 finite-time exponent from the exact n-step product ----
cm_expo(P, w0, n) = (w = P * w0; log(sqrt(sum(abs2, w)) / sqrt(sum(abs2, w0))) / n)

@testset "AC8 K=0 finite-time exponent" begin
    # |X(n) - ln(sqrt(n^2 + 1)) / n| <= 8 eps for w0 = one momentum component; X(n) = 0 exactly for a position component;
    # the K = 0 tangent is state independent (checked on 100 orbit points), so the n-step product is M0^n (exact integers, repeated squaring)
    ns = (1, 10, 100, 1000, 10_000, 1_000_000)
    for (T, systems) in ((Float64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0), cm_y4(5, 0, 0; ring = true))), (Double64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0))))
        for c in systems
            sys = cm_sys(T, c)
            d, h = cm_d(c), cm_h(c)
            rng = StableRNG(CM_SEED + 8001)
            z = cm_state(T, c, cm_dyadic(rng, T, d))
            M0 = Matrix(tangent(sys, z))
            same = true
            for _ in 1:100
                z = forward(sys, z)
                same &= Matrix(tangent(sys, z)) == M0
            end
            @test same
            wp = zeros(T, d)
            wp[h + 1] = one(T)
            wq = zeros(T, d)
            wq[1] = one(T)
            for n in ns
                P = M0^n
                @test P == cm_shear(T, d, n)
                X = cm_expo(P, wp, n)
                setprecision(BigFloat, CM_BITS) do
                    ref = log(sqrt(BigFloat(n)^2 + 1)) / n
                    @test Float64(abs(cm_big(X) - ref)) <= 8 * cm_eps(T)
                end
                @test cm_expo(P, wq, n) == 0
            end
        end
    end
end

# ---- AC9: decoupled limit, bitwise ----
const CM_DECOUPLED_SETS = ((0.5, 0.1), (0.0, 0.0), (1.0, 1.0), (2.0, 0.5), (-1.5, 3.0))

# (a) Y3 with mu = 0 mapped by the shear C(x, p) = (x + p, p) equals the verified Y1 module with K = -nu (pair 1) and K = -kappa (pair 2)
function cm_ac9a_defects(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    bad = 0
    for (nu, ka) in CM_DECOUPLED_SETS, _ in 1:m
        v = cm_draw(rng, 4, 4 * Float64(pi))
        sys = FroeschleMap(T(nu), T(ka), zero(T))
        s1, s2 = StandardMap(T(-nu)), StandardMap(T(-ka))
        z = cm_state(T, cm_y3(nu, ka, 0), v)
        Y1 = SVector(z[1] + z[3], z[3])
        Y2 = SVector(z[2] + z[4], z[4])
        for _ in 1:1000
            z = forward(sys, z)
            Y1 = forward(s1, Y1)
            Y2 = forward(s2, Y2)
            ok = z[1] + z[3] == Y1[1] && z[3] == Y1[2] && z[2] + z[4] == Y2[1] && z[4] == Y2[2]
            bad += ok ? 0 : 1
        end
    end
    return bad
end

# (b) Y4 with eps = 0: every (q_j, p_j) equals the Y1 orbit with the same K
function cm_ac9b_defects(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    bad = 0
    for N in (2, 5, 10), K in (1.0, 3.0), ring in (false, true)
        (ring && N < 3) && continue
        c = cm_y4(N, K, 0.0; ring)
        sys = cm_sys(T, c)
        s1 = StandardMap(T(K))
        for _ in 1:m
            z = cm_state(T, c, cm_draw(rng, 2N, 4 * Float64(pi)))
            Y = [SVector(z[j], z[N + j]) for j in 1:N]
            for _ in 1:1000
                z = forward(sys, z)
                Y = [forward(s1, y) for y in Y]
                bad += all(j -> z[j] == Y[j][1] && z[N + j] == Y[j][2], 1:N) ? 0 : 1
            end
        end
    end
    return bad
end

# (c) independence: perturbing pair 2 (or site k) leaves the orbit of pair 1 (or all other sites) unchanged to the bit
function cm_ac9c_defects(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    bad = 0
    sys3 = FroeschleMap(T(2.0), T(0.5), zero(T))
    for _ in 1:m
        v = cm_state(T, cm_y3(2, 0.5, 0), cm_draw(rng, 4, 4 * Float64(pi)))
        w = SVector(v[1], v[2] + T(rand(rng) * 7 - 3.5), v[3], v[4] + T(rand(rng) * 7 - 3.5))
        a, b = v, w
        for _ in 1:1000
            a, b = forward(sys3, a), forward(sys3, b)
            bad += (a[1] == b[1] && a[3] == b[3]) ? 0 : 1
        end
    end
    for N in (3, 10), ring in (false, true)
        c = cm_y4(N, 3.0, 0.0; ring)
        sys = cm_sys(T, c)
        for _ in 1:m
            v = cm_state(T, c, cm_draw(rng, 2N, 4 * Float64(pi)))
            k = rand(rng, 1:N)
            w = copy(v)
            w[k] += T(rand(rng) * 7 - 3.5)
            w[N + k] += T(rand(rng) * 7 - 3.5)
            a, b = v, w
            for _ in 1:1000
                a, b = forward(sys, a), forward(sys, b)
                bad += all(j -> j == k || (a[j] == b[j] && a[N + j] == b[N + j]), 1:N) ? 0 : 1
            end
        end
    end
    return bad
end

@testset "AC9 decoupled limit bitwise" begin
    # (a) Y3 mu = 0 vs Y1 through the shear C, 1000 steps, zero difference; 3 orbits per (nu, kappa) Float64, 1 Double64 (200 in verify)
    @test cm_ac9a_defects(Float64, 3, CM_SEED + 9001) == 0
    @test cm_ac9a_defects(Double64, 1, CM_SEED + 9002) == 0
    # (b) Y4 eps = 0, N in {2, 5, 10}, K in {1, 3}, open and ring: every site equals Y1
    @test cm_ac9b_defects(Float64, 2, CM_SEED + 9003) == 0
    @test cm_ac9b_defects(Double64, 1, CM_SEED + 9004) == 0
    # (c) independence of the other pair or sites
    @test cm_ac9c_defects(Float64, 3, CM_SEED + 9005) == 0
    @test cm_ac9c_defects(Double64, 1, CM_SEED + 9006) == 0
    # control: with coupling the same perturbation does change the other pair or the neighbours (the independence test is not vacuous)
    sys3 = FroeschleMap(2.0, 0.5, 0.5)
    a, b = SVector(0.3, 0.4, 0.5, 0.6), SVector(0.3, 0.9, 0.5, 0.6)
    for _ in 1:3
        a, b = forward(sys3, a), forward(sys3, b)
    end
    @test a[1] != b[1]
    sys4 = StandardMapChain(3, 3.0, 0.6102)
    a, b = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8], [0.3, 0.9, 0.5, 0.6, 0.7, 0.8]
    for _ in 1:3
        a, b = forward(sys4, a), forward(sys4, b)
    end
    @test a[1] != b[1]
end

# ---- AC10: parameter symmetries ----
function cm_ac10_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    piT = exact(T, pi)
    for (c, z) in cases
        h = cm_h(c)
        cA = c
        cB = c.kind === :y3 ? cm_y3(-c.p[1], -c.p[2], c.p[3]) : cm_y4(c.N, -c.p[1], c.p[2]; ring = c.ring)
        sysA, sysB = cm_sys(T, cA), cm_sys(T, cB)
        nH = cm_nH(c)
        zA = cm_state(T, c, z)
        zB = copy(zA)
        zB = typeof(zB)([i <= h ? zB[i] + piT : zB[i] for i in 1:length(zB)])
        lifted = cm_orbit(sysA, zA, nH)
        for k in 1:nH
            zA = forward(sysA, zA)
            zB = forward(sysB, zB)
            back = [i <= h ? zB[i] - piT : zB[i] for i in 1:length(zB)]
            worst = max(worst, cm_tdist(back, zA) / cm_tau(cm_eps(T), c, lifted, k))
        end
    end
    return worst
end

@testset "AC10 parameter symmetries" begin
    # Y3: (nu, kappa) -> (-nu, -kappa), positions + pi, mu unchanged; Y4: K -> -K, positions + pi, eps unchanged; torus distance <= tau(k), k <= n_H;
    # 10 orbits per set Float64, 3 Double64 (200 in verify)
    @test cm_ac10_ratio(Float64, cm_cases(StableRNG(CM_SEED + 10001), 10)) <= 1
    @test cm_ac10_ratio(Double64, cm_cases(StableRNG(CM_SEED + 10002), 3; maxN = 5)) <= 1
end

# ---- AC11: Y3 fixed point ----
@testset "AC11 Y3 fixed point" begin
    # origin: 10^4 steps stay with zero difference for every parameter set (Float64 and Double64)
    for T in (Float64, Double64), p in (CM_Y3_SETS..., (-6.0, 3.0, 2.5), (4.0, 4.0, 5.0))
        sys = cm_sys(T, cm_y3(p...))
        z = zero(SVector{4, T})
        ok = true
        for _ in 1:10_000
            z = forward(sys, z)
            ok &= all(iszero, z)
        end
        @test ok
    end
    # blocks of M(0): trace 2 - nu and 2 - kappa within 8 eps (1 + |nu| + |kappa|); class, rotation angle or ln lambda_+ against items 7
    vals = (-6.0, -2.0, 1.0, 3.0, 4.5, 6.0)
    for T in (Float64, Double64)
        worst_tr, worst_ev, wrong_class, off = 0.0, 0.0, 0, 0.0
        for (i, nu) in enumerate(vals)
            ka = vals[mod1(i + 2, length(vals))]
            c = cm_y3(nu, ka, 1.0e-3)
            M = Matrix(tangent(cm_sys(T, c), zero(SVector{4, T})))
            off = max(off, Float64(maximum(abs, M[[1, 3], [2, 4]])), Float64(maximum(abs, M[[2, 4], [1, 3]])))   # pairs decouple at the origin (b = 0)
            for (pair, par) in (([1, 3], nu), ([2, 4], ka))
                tr = M[pair[1], pair[1]] + M[pair[2], pair[2]]
                setprecision(BigFloat, CM_BITS) do
                    trx = 2 - BigFloat(par)
                    kx, mx = cm_class(trx)
                    tol = 8 * cm_eps(T) * (1 + abs(nu) + abs(ka))
                    worst_tr = max(worst_tr, Float64(abs(cm_big(tr) - trx)) / tol)
                    kind, _ = cm_class(cm_big(tr))
                    wrong_class += kind == kx ? 0 : 1
                    _, mag = cm_class(cm_big(tr))
                    worst_ev = max(worst_ev, Float64(abs(mag - mx)) / (32 * cm_eps(T) * (1 + abs(par)) * max(1, 1 / Float64(sqrt(abs(trx^2 - 4))))))
                end
            end
        end
        @test off == 0
        @test worst_tr <= 1
        @test worst_ev <= 1
        @test wrong_class == 0
    end
    # nu = 4: Jordan block [[1, 1], [-4, -3]] in the pair (1, 3); max entry of M^n is 4n
    for T in (Float64, Double64)
        M = Matrix(tangent(cm_sys(T, cm_y3(4, 0.1, 1.0e-3)), zero(SVector{4, T})))
        B = M[[1, 3], [1, 3]]
        P = B
        worst = 0.0
        for n in 2:1000
            P = B * P
            n in (10, 100, 1000) && (worst = max(worst, abs(Float64(maximum(abs, P)) / n - 4)))
        end
        @test worst <= 1.0e-12
    end
end

# ---- AC12: Y4 fixed point and exact spectrum of the origin ----
function cm_ac12_data(N, ring, K, e)
    Tf = Float64
    c = cm_y4(N, K, e; ring)
    sys = cm_sys(Tf, c)
    M = Matrix(tangent(sys, zeros(Tf, 2N)))
    G = M[(N + 1):(2N), 1:N]
    ev = eigen(Symmetric(G))
    gk = cm_laplacian_eigs(N, ring)
    # closed-form eigenvalues K + eps l_k, ascending for eps > 0 and descending for eps < 0
    gx = sort!(setprecision(BigFloat, CM_BITS) do
        Float64.([BigFloat(K) + BigFloat(e) * BigFloat(l) for l in gk])
    end)
    return c, M, G, ev, gx
end

@testset "AC12 Y4 fixed point and exact spectrum" begin
    # origin stays for 10^4 steps (Float64 N = 3 ring, N = 100, Double64 N = 5); eigenvalues of G(0) vs K + eps l_k to 16 eps (|K| + 4|eps|);
    # off-block entries of (U + U)' M(0) (U + U) below 16 N eps (1 + a_max) (the factor N is reading R7: LAPACK's eigenvector backward error grows with N, raw ratio 44 at N = 100); class and exponent per mode (within a rounding band of the parabolic edge
    # |tr| = 2 the class is not asserted, away from it the angle or ln lambda_+ is checked to 32 eps (1 + a_max) times the condition factor)
    for (T, c) in ((Float64, cm_y4(3, 1, 0.6102; ring = true)), (Float64, cm_y4(100, 3, 0.2)), (Double64, cm_y4(5, 3, -0.5)))
        sys = cm_sys(T, c)
        z = zeros(T, cm_d(c))
        ok = true
        for _ in 1:10_000
            z = forward(sys, z)
            ok &= all(iszero, z)
        end
        @test ok
    end
    worst_ev, worst_off, worst_exp, wrong_class, tr_err = 0.0, 0.0, 0.0, 0, 0.0
    for N in (2, 3, 5, 10, 50, 100), ring in (false, true), K in (-2.0, 0.5, 1.0, 3.0), e in (0.2, 0.6102, -0.5)
        (ring && N < 3) && continue
        c, M, G, ev, gx = cm_ac12_data(N, ring, K, e)
        epsT = cm_eps(Float64)
        a = 1 + cm_amax(c)
        worst_ev = max(worst_ev, maximum(abs, ev.values .- gx) / (16 * epsT * (abs(K) + 4 * abs(e))))
        U = ev.vectors
        Z = zeros(N, N)
        UU = [U Z; Z U]
        Mt = UU' * M * UU
        offmask = trues(2N, 2N)
        for k in 1:N, (i, j) in ((k, k), (k, N + k), (N + k, k), (N + k, N + k))
            offmask[i, j] = false
        end
        worst_off = max(worst_off, maximum(abs.(Mt[offmask])) / (16 * N * epsT * a))
        for k in 1:N
            g = ev.values[k]
            trx = 2 + gx[k]
            band = abs(abs(trx) - 2)
            tr_err = max(tr_err, abs((2 + g) - trx) / (32 * epsT * a))
            if band > 64 * epsT * a
                kx, mx = cm_class(trx)
                kind, mag = cm_class(2 + g)
                wrong_class += kind == kx ? 0 : 1
                if band > 1.0e-6
                    cond = max(1.0, 1 / sqrt(abs(trx^2 - 4)))
                    worst_exp = max(worst_exp, abs(mag - mx) / (32 * epsT * a * cond))
                end
            end
        end
    end
    @test worst_ev <= 1
    @test worst_off <= 1
    @test tr_err <= 1
    @test worst_exp <= 1
    @test wrong_class == 0
end
