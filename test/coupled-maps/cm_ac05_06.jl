# Wrapped systems with winding counts: criteria 5 and 6 (policy D28, spec Algorithm item 3).
# State layout: (z_1..z_d, c_1..c_d) with z the reduced coordinates in [0, 2 pi_T) in the ordering (q, p) and c_j the integer
# count of coordinate j; the lifted point is z_j + 2 pi c_j. Length 8 for Y3, 4N for the chain. Helpers prefixed cmw_.

# 256-bit recomputation of the count increments of one wrapped step from the exact binary64 (or Double64) states s0 -> s1
# (teacher forced, as W3 of standard-map): returns true when every count of s1 equals the 256-bit one.
function cmw_forced_ok(c::CMCase, s0, s1)
    d = cm_d(c)
    h = d ÷ 2
    return setprecision(BigFloat, CM_BITS) do
        tp = 2 * BigFloat(pi)
        B(v) = cm_big(v)
        if c.kind === :y3
            nu, ka, mu = BigFloat.(c.p)
            b = [floor((B(s0[j]) + B(s0[j + 2])) / tp) for j in 1:2]                 # position reductions
            ok = all(j -> B(s1[d + j] - s0[d + j] - s0[d + j + 2]) == b[j], 1:2)          # m_j' = m_j + n_j + b_j (old n_j)
            x1, x2 = B(s1[1]), B(s1[2])                                              # forces see the reduced new positions
            f = mu * (1 - cos(x1 + x2))
            a1 = floor((B(s0[3]) - nu * sin(x1) - f) / tp)
            a2 = floor((B(s0[4]) - ka * sin(x2) - f) / tp)
            return ok && B(s1[d + 3] - s0[d + 3]) == a1 && B(s1[d + 4] - s0[d + 4]) == a2
        end
        N = c.N
        K, e = BigFloat(c.p[1]), BigFloat(c.p[2])
        F = cm_force(K, e, cm_nbrs(N, c.ring), [B(s0[j]) for j in 1:N])           # force at the pre-step reduced positions
        a = [floor((B(s0[N + j]) + F[j]) / tp) for j in 1:N]
        ok = all(j -> B(s1[d + N + j] - s0[d + N + j]) == a[j], 1:N)                  # n_j' = n_j + a_j
        for j in 1:N
            bj = floor((B(s0[j]) + B(s1[N + j])) / tp)                            # x~ = x + p' (new reduced momentum)
            ok &= B(s1[d + j] - s0[d + j] - s1[d + N + j]) == bj                     # m_j' = m_j + n_j' + b_j
        end
        return ok
    end
end

# AC5 helper: wrapped orbit vs lifted orbit (torus distance and lift at 256 bits), range and integrality of the state
function cm_ac5_ratios(::Type{T}, cases) where {T}
    wt, wl, ok = 0.0, 0.0, true
    tp = 2 * exact(T, pi)
    for (c, z) in cases
        d = cm_d(c)
        nH = cm_nH(c)
        zs = cm_orbit(cm_sys(T, c), cm_state(T, c, z), nH)
        wsys = cm_wsys(T, c)
        s = cm_wstate(T, c, z)
        for k in 1:nH
            s = forward(wsys, s)
            ok &= all(j -> zero(T) <= s[j] < tp, 1:d) && all(j -> s[d + j] == round(s[d + j]), 1:d)
            tau = cm_tau(cm_eps(T), c, zs, k)
            wt = max(wt, cm_tdist(s, zs[k + 1], d) / tau)
            L = setprecision(BigFloat, CM_BITS) do
                wb = cm_wsys(BigFloat, c)
                sb = c.kind === :y3 ? SVector{8, BigFloat}(ntuple(i -> cm_big(s[i]), 8)) : BigFloat[cm_big(v) for v in s]
                lift(wb, sb)
            end
            wl = max(wl, cm_absdiff(L, zs[k + 1]) / tau)
        end
    end
    return ok, wt, wl
end

# counts of the first 50 steps against the 256-bit increments, teacher forced; also returns the range check of the 50 states
function cm_ac5_counts(::Type{T}, cases) where {T}
    ok = true
    for (c, z) in cases
        wsys = cm_wsys(T, c)
        s = cm_wstate(T, c, z)
        for _ in 1:50
            s1 = forward(wsys, s)
            ok &= cmw_forced_ok(c, s, s1)
            s = s1
        end
    end
    return ok
end

@testset "AC5 lifted vs wrapped" begin
    # wrapped orbit vs lifted orbit, torus distance <= tau(k), k <= n_H; 256-bit lift of the state vs the lifted orbit <= tau(k);
    # counts equal to the 256-bit increments for 50 teacher-forced steps; 5 orbits per set Float64, 2 Double64 (200 in verify)
    ok, wt, wl = cm_ac5_ratios(Float64, cm_cases(StableRNG(CM_SEED + 5001), 5; maxN = 5, wrapped = true))
    @test ok
    @test wt <= 1
    @test wl <= 1
    ok, wt, wl = cm_ac5_ratios(Double64, cm_cases(StableRNG(CM_SEED + 5002), 2; maxN = 3, wrapped = true))
    @test ok
    @test wt <= 1
    @test wl <= 1
    @test cm_ac5_counts(Float64, cm_cases(StableRNG(CM_SEED + 5003), 2; maxN = 5, wrapped = true))
    @test cm_ac5_counts(Double64, cm_cases(StableRNG(CM_SEED + 5004), 1; maxN = 3, wrapped = true))
end

@testset "AC5 torus_state and lift" begin
    # torus_state(sys, z) reduces each coordinate with its count; lift(sys, s) returns z_j + 2 pi c_j in the element type
    rng = StableRNG(CM_SEED + 5011)
    for c in (cm_y3(0.5, 0.1, 1.0e-3), cm_y4(2, 3.0, 0.2), cm_y4(5, 1.0, 0.6102; ring = true))
        d = cm_d(c)
        wsys = cm_wsys(Float64, c)
        for _ in 1:20
            zr = cm_draw(rng, d, 100.0)
            s = torus_state(wsys, cm_state(Float64, c, zr))
            @test length(s) == 2d
            @test all(j -> 0 <= s[j] < 2 * Float64(pi), 1:d)
            @test all(j -> s[d + j] == round(s[d + j]), 1:d)
            @test all(abs.(lift(wsys, s) .- zr) .<= 16 * eps(100.0))
        end
    end
    # a known point: z_1 = 7 - 2 pi with count 1 and z_3 = -1 + 2 pi with count -1 (positions first for Y3)
    s = torus_state(WrappedFroeschleMap(0.5, 0.1, 1.0e-3), SVector(7.0, 0.5, -1.0, 0.25))
    @test s isa SVector{8, Float64}
    @test s[5] == 1 && s[7] == -1 && s[6] == 0 && s[8] == 0
    @test abs(s[1] - (7.0 - 2 * Float64(pi))) <= 4 * eps(7.0)
end

# ---- AC6: wrapped inverse ----
function cm_ac6_ratio(::Type{T}, cases) where {T}
    worst, counts_ok = 0.0, true
    for (c, z) in cases
        d = cm_d(c)
        lsys, wsys = cm_sys(T, c), cm_wsys(T, c)
        s0 = cm_wstate(T, c, z)
        for n in unique((1, cm_nH(c)))
            fw = cm_orbit(lsys, cm_state(T, c, z), n)
            iv = cm_orbit(lsys, cm_state(T, c, z), n; f = inverse)
            s = cm_orbit(wsys, s0, n)[end]
            s = cm_orbit(wsys, s, n; f = inverse)[end]
            worst = max(worst, cm_tdist(s, s0, d) / cm_tau(cm_eps(T), c, fw, n))
            counts_ok &= all(j -> s[d + j] == s0[d + j], 1:d)
            s = cm_orbit(wsys, s0, n; f = inverse)[end]
            s = cm_orbit(wsys, s, n)[end]
            worst = max(worst, cm_tdist(s, s0, d) / cm_tau(cm_eps(T), c, iv, n))
            counts_ok &= all(j -> s[d + j] == s0[d + j], 1:d)
        end
    end
    return counts_ok, worst
end

# overflow construction: the count that must increase sits at maxintfloat - 1 (NaN in all components) or maxintfloat - 2 (finite)
function cm_overflow_checks()
    lim = maxintfloat(Float64)
    res = Bool[]
    # Y3, nu = kappa = mu = 0: x1 + x3 = 9 >= 2 pi, so m_1' = m_1 + n_1 + 1
    w0 = WrappedFroeschleMap(0.0, 0.0, 0.0)
    up = SVector(6.0, 0.0, 3.0, 0.0, lim - 1, 0.0, 0.0, 0.0)
    push!(res, all(isnan, forward(w0, up)))
    ok = forward(w0, SVector(6.0, 0.0, 3.0, 0.0, lim - 2, 0.0, 0.0, 0.0))
    push!(res, all(isfinite, ok) && ok[5] == lim - 1)
    # Y3, nu = -6: x1' = pi/2, p~ = 1 + 6 = 7 >= 2 pi, so n_3' = n_3 + 1 (the position count stays below the limit: b = 0)
    w1 = WrappedFroeschleMap(-6.0, 0.0, 0.0)
    x1 = Float64(pi) / 2 - 1
    push!(res, all(isnan, forward(w1, SVector(x1, 0.0, 1.0, 0.0, 0.0, 0.0, lim - 1, 0.0))))
    ok = forward(w1, SVector(x1, 0.0, 1.0, 0.0, 0.0, 0.0, lim - 2, 0.0))
    push!(res, all(isfinite, ok) && ok[7] == lim - 1)
    # Y3 inverse: x3 = 3, x1' = 1 gives x1 = -2, m_1 = m_1' - n_1 - 1 = -2^53
    push!(res, all(isnan, inverse(w0, SVector(1.0, 0.0, 3.0, 0.0, -(lim - 1), 0.0, 0.0, 0.0))))
    # chain N = 2, eps = 0: K = 0, q_1 + p_1 = 9 >= 2 pi (position count), K = 6 with q_1 = pi/2, p_1 = 1 (momentum count)
    wc0 = WrappedStandardMapChain(2, 0.0, 0.0)
    up = [6.0, 0.0, 3.0, 0.0, lim - 1, 0.0, 0.0, 0.0]
    push!(res, all(isnan, forward(wc0, up)))
    ok = forward(wc0, [6.0, 0.0, 3.0, 0.0, lim - 2, 0.0, 0.0, 0.0])
    push!(res, all(isfinite, ok) && ok[5] == lim - 1)
    wc1 = WrappedStandardMapChain(2, 6.0, 0.0)
    push!(res, all(isnan, forward(wc1, [Float64(pi) / 2, 0.0, 1.0, 0.0, 0.0, 0.0, lim - 1, 0.0])))
    ok = forward(wc1, [Float64(pi) / 2, 0.0, 1.0, 0.0, 0.0, 0.0, lim - 2, 0.0])
    push!(res, all(isfinite, ok) && ok[7] == lim - 1)
    push!(res, all(isnan, inverse(wc0, [1.0, 0.0, 3.0, 0.0, -(lim - 1), 0.0, 0.0, 0.0])))
    return res
end

@testset "AC6 wrapped inverse" begin
    # n forward then n inverse steps (and inverse then forward) recover coordinates to tau(n) and counts exactly, n in {1, n_H};
    # 5 orbits per set Float64, 2 Double64 (200 in verify)
    counts_ok, worst = cm_ac6_ratio(Float64, cm_cases(StableRNG(CM_SEED + 6001), 5; maxN = 5, wrapped = true))
    @test counts_ok
    @test worst <= 1
    counts_ok, worst = cm_ac6_ratio(Double64, cm_cases(StableRNG(CM_SEED + 6002), 2; maxN = 3, wrapped = true))
    @test counts_ok
    @test worst <= 1
end

@testset "AC6 overflow guard" begin
    # a count one below 2^53 that must increase gives NaN in all components; two below stays finite and exact (W8)
    res = cm_overflow_checks()
    @test length(res) == 10
    @test all(res)
end
