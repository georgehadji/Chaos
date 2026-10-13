# Amendment D28: wrapped standard map with winding counts (criteria W1 to W8). State (x, p, m, n), prefix smw_ for helpers.
# Spec: workspaces/module-build/stages/01-spec/output/standard-map-spec.md, section "Amendment D28".

const SMW_KS = (0.5, 2.0, 5.0, 10.0)

smw_draw(rng) = SVector(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng))
smw_state(::Type{T}, z) where {T} = SVector{4, T}(T(z[1]), T(z[2]), zero(T), zero(T))
smw_xp(s) = SVector(s[1], s[2])
smw_isint(v) = v == round(v)

# W1: wrapped state vs lifted orbit, torus distance <= tau(k) for k <= n_H(K); also range and integrality of the state.
function smw_w1(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    worst, ok = 0.0, true
    tp = 2 * exact(T, pi)
    for K in SMW_KS, _ in 1:m
        z0 = smw_draw(rng)
        zs = sm_orbit(StandardMap(T(K)), sm_cast(T, z0), sm_nH(K))
        sys = WrappedStandardMap(T(K))
        s = smw_state(T, z0)
        for k in 1:sm_nH(K)
            s = forward(sys, s)
            ok &= (zero(T) <= s[1] < tp) && (zero(T) <= s[2] < tp) && smw_isint(s[3]) && smw_isint(s[4])
            worst = max(worst, sm_tdist(smw_xp(s), zs[k + 1]) / sm_tau(sm_eps(T), K, zs, k))
        end
    end
    return ok, worst
end

# W2: lift of the state at 256 bits (state widened exactly) vs the lifted orbit, plain max norm <= tau(k).
function smw_w2(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    worst = 0.0
    for K in SMW_KS, _ in 1:m
        z0 = smw_draw(rng)
        zs = sm_orbit(StandardMap(T(K)), sm_cast(T, z0), sm_nH(K))
        sys = WrappedStandardMap(T(K))
        s = smw_state(T, z0)
        for k in 1:sm_nH(K)
            s = forward(sys, s)
            L = setprecision(BigFloat, SM_BITS) do
                lift(SVector{4, BigFloat}(sm_big(s[1]), sm_big(s[2]), sm_big(s[3]), sm_big(s[4])))
            end
            worst = max(worst, sm_absdiff(L, zs[k + 1]) / sm_tau(sm_eps(T), K, zs, k))
        end
    end
    return worst
end

# Independent 256-bit wrapped step (formulas of the amendment, floor only; no Chaoslab call). Call inside setprecision.
function smw_big_step(K, s)
    tp = 2 * BigFloat(pi)
    pt = s[2] + K * sin(s[1])
    a = floor(pt / tp)
    p = pt - tp * a
    xt = s[1] + p
    b = floor(xt / tp)
    n = s[4] + a
    return (xt - tp * b, p, s[3] + n + b, n)
end

# W3 (a): 1e6 steps, returns (range ok, integer ok, counts below the exactness limit, first 51 states).
function smw_w3_run(K, s, N)
    first51 = Vector{typeof(s)}(undef, 51)
    sys = WrappedStandardMap(K)
    tp = 2 * Float64(pi)
    inrange, ints, bounded = true, true, true
    for k in 0:N
        k <= 50 && (first51[k + 1] = s)
        inrange &= (0.0 <= s[1] < tp) && (0.0 <= s[2] < tp)
        ints &= smw_isint(s[3]) && smw_isint(s[4])
        bounded &= abs(s[3]) < maxintfloat(Float64) && abs(s[4]) < maxintfloat(Float64)
        k < N && (s = forward(sys, s))
    end
    return inrange, ints, bounded, first51
end

# W3 (b): teacher-forced 256-bit increments from the exact binary64 state of each step.
function smw_w3_forced(K, states)
    return setprecision(BigFloat, SM_BITS) do
        tp, Kb = 2 * BigFloat(pi), BigFloat(K)
        all(1:(length(states) - 1)) do k
            s0, s1 = states[k], states[k + 1]
            a = floor((BigFloat(s0[2]) + Kb * sin(BigFloat(s0[1]))) / tp)
            b = floor((BigFloat(s0[1]) + BigFloat(s1[2])) / tp)
            a == s1[4] - s0[4] && b == s1[3] - s0[3] - s1[4]
        end
    end
end

# W3 (c): free-running 256-bit orbit vs the binary64 states for k <= n_H(K): torus ratio and exact counts.
function smw_w3_free(K, states, zs, nH)
    return setprecision(BigFloat, SM_BITS) do
        sb = (BigFloat(states[1][1]), BigFloat(states[1][2]), BigFloat(0), BigFloat(0))
        worst, same = 0.0, true
        for k in 1:nH
            sb = smw_big_step(BigFloat(K), sb)
            s = states[k + 1]
            worst = max(worst, sm_tdist(smw_xp(s), SVector(sb[1], sb[2])) / sm_tau(sm_eps(Float64), K, zs, k))
            same &= sb[3] == s[3] && sb[4] == s[4]
        end
        (worst, same)
    end
end

# W4: K = 0, m_k = floor((x0 + k p0) / 2 pi) at 256 bits, p and n unchanged, k <= 1000.
function smw_w4(seed)
    rng = StableRNG(seed)
    sys = WrappedStandardMap(0.0)
    ok = true
    for _ in 1:20
        x0, p0 = rand(rng, 0:6433) / 1024, rand(rng, 0:6433) / 1024
        s = SVector(x0, p0, 0.0, 0.0)
        ms = setprecision(BigFloat, SM_BITS) do
            tp = 2 * BigFloat(pi)
            [floor((BigFloat(x0) + k * BigFloat(p0)) / tp) for k in 1:1000]
        end
        for k in 1:1000
            s = forward(sys, s)
            ok &= s[2] == p0 && s[4] == 0 && s[3] == ms[k]
        end
    end
    return ok
end

# W5: tangent product along the wrapped orbit vs the 256-bit product at the same points (tangent-core AC4 bound), and vs the lifted product.
function smw_w5(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    worstA, worstB = 0.0, 0.0
    for K in SMW_KS, _ in 1:m
        z0 = smw_draw(rng)
        nH = sm_nH(K)
        lsys, wsys = StandardMap(T(K)), WrappedStandardMap(T(K))
        zs = sm_orbit(lsys, sm_cast(T, z0), nH)
        s = smw_state(T, z0)
        Y = Yl = SMatrix{2, 2, T}(one(T), zero(T), zero(T), one(T))
        setprecision(BigFloat, SM_BITS) do
            Yb, prodM = BigFloat[1 0; 0 1], one(BigFloat)
            for k in 1:nH
                c = BigFloat(K) * cos(sm_big(s[1]))
                Mb = [1 + c one(c); c one(c)]
                Yb, prodM = Mb * Yb, prodM * sqrt(sum(abs2, Mb))
                Y = tangent(wsys, s) * Y
                Yl = tangent(lsys, zs[k]) * Yl
                s = forward(wsys, s)
                boundA = 8 * k * BigFloat(sm_eps(T)) * prodM * sqrt(BigFloat(2))
                dA = sqrt(sum(abs2, [sm_big(Y[i, j]) - Yb[i, j] for i in 1:2, j in 1:2]))
                dB = sqrt(sum(abs2, [sm_big(Y[i, j]) - sm_big(Yl[i, j]) for i in 1:2, j in 1:2]))
                boundB = 2 * Float64(boundA) + k * abs(K) * sm_tau(sm_eps(T), K, zs, k) * sm_sigma(K)^k * sqrt(2.0)
                worstA = max(worstA, Float64(dA / boundA))
                worstB = max(worstB, Float64(dB) / boundB)
            end
        end
    end
    return worstA, worstB
end

# W6: n forward then n inverse (and inverse then forward): torus distance <= tau(n), counts exactly equal to the start (0, 0).
function smw_w6(::Type{T}, m, seed) where {T}
    rng = StableRNG(seed)
    worst, counts_ok = 0.0, true
    for K in SMW_KS, _ in 1:m, n in unique((1, sm_nH(K)))
        z0 = smw_draw(rng)
        lsys, sys = StandardMap(T(K)), WrappedStandardMap(T(K))
        s0 = smw_state(T, z0)
        fw = sm_orbit(lsys, sm_cast(T, z0), n)
        iv = sm_orbit(lsys, sm_cast(T, z0), n; f = inverse)
        s = s0
        for _ in 1:n
            s = forward(sys, s)
        end
        for _ in 1:n
            s = inverse(sys, s)
        end
        worst = max(worst, sm_tdist(smw_xp(s), smw_xp(s0)) / sm_tau(sm_eps(T), K, fw, n))
        counts_ok &= s[3] == s0[3] && s[4] == s0[4]
        s = s0
        for _ in 1:n
            s = inverse(sys, s)
        end
        for _ in 1:n
            s = forward(sys, s)
        end
        worst = max(worst, sm_tdist(smw_xp(s), smw_xp(s0)) / sm_tau(sm_eps(T), K, iv, n))
        counts_ok &= s[3] == s0[3] && s[4] == s0[4]
    end
    return counts_ok, worst
end

# W7: bytes allocated by one call after a warm-up call, inside a function.
function smw_allocs(sys, s, W0)
    forward(sys, s); inverse(sys, s); tangent(sys, s); inverse_tangent(sys, s); propagate_orbit(sys, s, W0, 100, 1, ())
    return (@allocated(forward(sys, s)), @allocated(inverse(sys, s)), @allocated(tangent(sys, s)),
            @allocated(inverse_tangent(sys, s)), @allocated(propagate_orbit(sys, s, W0, 100, 1, ())))
end

function smw_health(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

@testset "W1 torus agreement with the lifted map" begin
    # wrapped state vs lifted orbit, torus distance <= tau(k), k <= n_H(K), 200 orbits per K in {0.5, 2, 5, 10}
    ok, worst = smw_w1(Float64, 200, SM_SEED + 5001)
    @test ok
    @test worst <= 1
    ok, worst = smw_w1(Double64, 200, SM_SEED + 5002)
    @test ok
    @test worst <= 1
end

@testset "W2 lifted reconstruction" begin
    # lift(state) at 256 bits from the exact counts vs the lifted orbit, plain max norm <= tau(k); torus_state and lift round trip
    @test smw_w2(Float64, 200, SM_SEED + 5011) <= 1
    @test smw_w2(Double64, 200, SM_SEED + 5012) <= 1
    z = SVector(7.0, -1.0)   # x = 7 - 2 pi, m = 1; p = -1 + 2 pi, n = -1
    s = torus_state(z)
    @test s isa SVector{4, Float64}
    @test s[3] == 1 && s[4] == -1
    @test abs(s[1] - (7.0 - 2 * Float64(pi))) <= 4 * eps(7.0) && abs(s[2] - (-1.0 + 2 * Float64(pi))) <= 4 * eps(7.0)
    @test all(abs.(lift(s) - z) .<= 8 * eps(7.0))
    rng = StableRNG(SM_SEED + 5013)
    for _ in 1:100
        zr = SVector((2 * rand(rng) - 1) * 100, (2 * rand(rng) - 1) * 100)
        sr = torus_state(zr)
        @test 0 <= sr[1] < 2 * Float64(pi) && 0 <= sr[2] < 2 * Float64(pi)
        @test all(abs.(lift(sr) - zr) .<= 16 * eps(100.0))
    end
end

@testset "W3 counts exact and bounded" begin
    # 1e6 steps at K = 5: x, p in [0, 2 pi), m, n integers below 2^53; first 50 steps: teacher-forced 256-bit counts; free run to n_H
    z0 = SVector(1.0, 2.0)
    inrange, ints, bounded, first51 = smw_w3_run(5.0, SVector(z0[1], z0[2], 0.0, 0.0), 1_000_000)
    @test inrange
    @test ints
    @test bounded
    @test smw_w3_forced(5.0, first51)
    zs = sm_orbit(StandardMap(5.0), z0, sm_nH(5.0))
    worst, same = smw_w3_free(5.0, first51, zs, sm_nH(5.0))
    @test worst <= 1
    @test same
end

@testset "W4 K = 0 known answer" begin
    # p and n unchanged, m_k = floor((x0 + k p0) / 2 pi) at 256 bits, k <= 1000, dyadic x0, p0 = i / 1024 in [0, 2 pi)
    @test smw_w4(SM_SEED + 5041)
end

@testset "W5 tangent agreement" begin
    # product of wrapped tangents: 8 n eps prod||M||_F bound vs the 256-bit product at the same points; vs lifted product with orbit term
    a, b = smw_w5(Float64, 100, SM_SEED + 5051)
    @test a <= 1
    @test b <= 1
    a, b = smw_w5(Double64, 50, SM_SEED + 5052)
    @test a <= 1
    @test b <= 1
end

@testset "W6 inverse round trip" begin
    # n forward then n inverse (and reverse order): torus distance <= tau(n), counts exactly equal, n in {1, n_H(K)}, 200 orbits per K
    counts_ok, worst = smw_w6(Float64, 200, SM_SEED + 5061)
    @test counts_ok
    @test worst <= 1
    counts_ok, worst = smw_w6(Double64, 200, SM_SEED + 5062)
    @test counts_ok
    @test worst <= 1
end

@testset "W7 zero allocation" begin
    # forward, inverse, tangent, inverse_tangent and propagate_orbit allocate 0 bytes after warm-up
    for T in (Float64, Double64)
        sys = WrappedStandardMap(T(5))
        s = SVector{4, T}(T(1), T(2), zero(T), zero(T))
        W0 = SMatrix{2, 2, T}(one(T), zero(T), zero(T), one(T))
        @test smw_allocs(sys, s, W0) == (0, 0, 0, 0, 0)
    end
end

@testset "W8 overflow guard" begin
    # |m| or |n| reaching maxintfloat gives NaN in all components (forward and inverse); propagation raises "state"; two below is exact
    lim = maxintfloat(Float64)
    W0 = SMatrix{2, 2, Float64}(1.0, 0.0, 0.0, 1.0)
    sm_up = SVector(6.0, 3.0, lim - 1, 0.0)      # K = 0: x~ = 9 >= 2 pi, so m' = m + 1 = 2^53
    sn_up = SVector(Float64(pi) / 2, 5.0, 0.0, lim - 1)   # K = 3: p~ = 8 >= 2 pi, so n' = n + 1 = 2^53
    @test all(isnan, forward(WrappedStandardMap(0.0), sm_up))
    @test all(isnan, forward(WrappedStandardMap(3.0), sn_up))
    @test all(isnan, inverse(WrappedStandardMap(0.0), SVector(1.0, 3.0, -(lim - 1), 0.0)))   # x~ = -2: m = m' - 1 = -2^53
    err = smw_health(() -> propagate_orbit(WrappedStandardMap(0.0), sm_up, W0, 1, 1, ()))
    @test err isa NumericalHealthError && err.what == "state"
    err = smw_health(() -> propagate(WrappedStandardMap(3.0), [sn_up], 1; k = 2, seed = 1))
    @test err isa NumericalHealthError && err.what == "state" && err.orbit == 1 && err.step == 1
    ok = forward(WrappedStandardMap(0.0), SVector(6.0, 3.0, lim - 2, 0.0))
    @test all(isfinite, ok) && ok[3] == lim - 1
end
