# Criteria 13 to 16: pairing and sum at finite n, total momentum, conversion of the printed form, element-type agreement.

# ---- AC13: pairing and sum of the log-singular-value exponents of the 256-bit tangent product ----
function cm_ac13_worst(cases, ns)
    wpair, wsum = 0.0, 0.0
    for (c, z) in cases, n in ns
        sys = cm_sys(Float64, c)
        d = cm_d(c)
        zs = cm_orbit(sys, cm_state(Float64, c, z), n)
        setprecision(BigFloat, CM_BITS) do
            Y = Matrix{BigFloat}(I, d, d)
            prodF = 1.0
            for k in 1:n
                M = Matrix(tangent(sys, zs[k]))
                prodF *= sqrt(sum(abs2, M))
                Y = BigFloat[cm_big(M[i, j]) for i in 1:d, j in 1:d] * Y
            end
            sv = cm_svdvals(Y)
            lam = log.(sv) ./ n
            B = 8 * n * cm_eps(Float64) * prodF * sqrt(Float64(d))      # tangent-core criterion 4 bound, ||Y0||_F = sqrt(d)
            for i in 1:(d ÷ 2)
                wpair = max(wpair, Float64(abs(lam[i] + lam[d + 1 - i])) / (B / n))
            end
            wsum = max(wsum, Float64(abs(sum(lam))) / (B / n))
        end
    end
    return wpair, wsum
end

@testset "AC13 pairing and sum at finite n" begin
    # |lambda_i(n) + lambda_{d+1-i}(n)| <= B(n)/n and |sum lambda_i(n)| <= B(n)/n; Y3 (five sets) and Y4 N in {2, 3, 5}; n in {5, n_H};
    # 2 orbits per case (the 256-bit product and a Jacobi singular-value solver are test helpers)
    cases = cm_cases(StableRNG(CM_SEED + 13001), 2; maxN = 5, uniform = false)
    wpair, wsum = cm_ac13_worst(cases, (5,))
    @test wpair <= 1
    @test wsum <= 1
    # n_H differs per system: run each case at its own horizon
    wp2, ws2 = 0.0, 0.0
    for (c, z) in cases
        a, b = cm_ac13_worst([(c, z)], (cm_nH(c),))
        wp2, ws2 = max(wp2, a), max(ws2, b)
    end
    @test wp2 <= 1
    @test ws2 <= 1
end

# ---- AC14: total momentum for K = 0 ----
function cm_ac14_ratio(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    worst = 0.0
    for N in (2, 3, 5, 10), e in (0.2, 0.6102, 1.0, -0.5), ring in (false, true)
        (ring && N < 3) && continue
        c = cm_y4(N, 0.0, e; ring)
        sys = cm_sys(T, c)
        nH = cm_nH(c)
        for _ in 1:m
            z0 = cm_state(T, c, cm_draw(rng, 2N, 4 * Float64(pi)))
            zs = cm_orbit(sys, z0, nH)
            for k in 1:nH
                dp = setprecision(BigFloat, CM_BITS) do
                    Float64(abs(sum(cm_big(v) for v in zs[k + 1][(N + 1):(2N)]) - sum(cm_big(v) for v in z0[(N + 1):(2N)])))
                end
                worst = max(worst, dp / (8 * N * cm_eps(T) * cm_Mn(c, zs, k) * k))
            end
        end
    end
    return worst
end

@testset "AC14 total momentum for K = 0" begin
    # |sum p_j(n) - sum p_j(0)| <= 8 N eps M_n n for every n <= n_H, open chain and ring, eps in {0.2, 0.6102, 1, -0.5}; 5 orbits Float64, 2 Double64
    @test cm_ac14_ratio(Float64, 5, CM_SEED + 14001) <= 1
    @test cm_ac14_ratio(Double64, 2, CM_SEED + 14002) <= 1
end

# ---- AC15: conversion of the printed unit-torus form (item 8), N = 2 ----
# printed lines in the symplectic reading (force at the pre-step positions, position update with the new momentum), labelling (x1, x2, x3, x4) = (q1, p1, q2, p2)
function cm_printed_orbit(K, beta, xu, n)
    x1, x2, x3, x4 = xu
    out = [(x1, x2, x3, x4)]
    for _ in 1:n
        x2n = x2 + (K / (2 * pi)) * sin(2 * pi * x1) - (beta / pi) * sin(2 * pi * (x3 - x1))
        x4n = x4 + (K / (2 * pi)) * sin(2 * pi * x3) - (beta / pi) * sin(2 * pi * (x1 - x3))
        x1, x2, x3, x4 = x1 + x2n, x2n, x3 + x4n, x4n
        push!(out, (x1, x2, x3, x4))
    end
    return out
end

@testset "AC15 conversion of the printed form" begin
    # q = 2 pi x_u, p = 2 pi y_u, eps = 2 beta: the printed orbit and the chain (N = 2, open) differ in torus distance by <= tau(n), n <= n_H; 20 orbits per (beta, K)
    rng = StableRNG(CM_SEED + 15001)
    worst = 0.0
    for beta in (0.1, 0.3051), K in (1.0, 3.0), _ in 1:20
        c = cm_y4(2, K, 2 * beta)
        z = cm_draw(rng, 4, 4 * Float64(pi))                      # library ordering (q1, q2, p1, p2)
        nH = cm_nH(c)
        zs = cm_orbit(cm_sys(Float64, c), cm_state(Float64, c, z), nH)
        pr = cm_printed_orbit(K, beta, (z[1], z[3], z[2], z[4]) ./ (2 * pi), nH)   # printed (q1, p1, q2, p2) = library (z1, z3, z2, z4)
        for k in 1:nH
            conv = [2 * pi * pr[k + 1][1], 2 * pi * pr[k + 1][3], 2 * pi * pr[k + 1][2], 2 * pi * pr[k + 1][4]]
            worst = max(worst, cm_tdist(conv, zs[k + 1]) / cm_tau(cm_eps(Float64), c, zs, k))
        end
    end
    @test worst <= 1
end

# ---- AC16: element-type agreement ----
# trig arguments seen by the module on an orbit: positions, neighbour differences (chain) or x1 + x2 (Y3)
function cm_trig_args(c::CMCase, z)
    h = cm_h(c)
    xs = Any[z[j] for j in 1:h]
    c.kind === :y3 ? push!(xs, z[1] + z[2]) : append!(xs, [z[j] - z[j + 1] for j in 1:(h - 1)])
    return xs
end

function cm_ac16_ratio(::Type{T}, cases, epsT) where {T}
    worst = 0.0
    for (c, z) in cases
        nH = cm_nH(c)
        ref = cm_orbit(cm_sys(BigFloat, c), cm_state(BigFloat, c, z), nH)
        low = cm_orbit(cm_sys(T, c), cm_state(T, c, z), nH)
        tau = cm_tau(epsT, c, ref, nH)
        worst = max(worst, cm_absdiff(low[end], ref[end]) / tau, cm_tdist(low[end], ref[end]) / tau)
    end
    return worst
end

@testset "AC16 element-type agreement" begin
    # Float64 and Double64 vs the 256-bit orbit from identical exactly representable data at n = n_H, bound tau(n_H) with eps of the lower precision;
    # 5 orbits per set Float64, 2 Double64 (200 in verify); Y4 N = 100: 2 orbits Float64 and 1 Double64 for K in {1, 5}, eps in {0.2, 1}
    rng = StableRNG(CM_SEED + 16001)
    big100 = Tuple{CMCase, Vector{Float64}}[]
    for K in (1.0, 5.0), e in (0.2, 1.0)
        c = cm_y4(100, K, e)
        push!(big100, (c, cm_draw(rng, 200, 4 * Float64(pi))))
        push!(big100, (c, cm_draw(rng, 200, 4 * Float64(pi))))
    end
    cases = vcat(cm_cases(StableRNG(CM_SEED + 16002), 5), big100)
    setprecision(BigFloat, CM_BITS) do
        @test cm_ac16_ratio(Float64, cases, cm_eps(Float64)) <= 1
        # Double64: effective eps = largest relative sin / cos error on the arguments of the sample (floor 2^-104); pass requires <= 2^-100
        dd = vcat(cm_cases(StableRNG(CM_SEED + 16003), 2; maxN = 5), big100[1:2:end])
        xs = Double64[]
        for (c, z) in dd
            for s in cm_orbit(cm_sys(Double64, c), cm_state(Double64, c, z), cm_nH(c))
                append!(xs, cm_trig_args(c, s))
            end
        end
        eff = max(CM_EPS_DD_FLOOR, cm_dd_trig_error(xs))
        @test cm_ac16_ratio(Double64, dd, eff) <= 1
        @test eff <= CM_EPS_DD
    end
end
