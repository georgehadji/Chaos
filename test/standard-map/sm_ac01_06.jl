# Map invariants: criteria 1 to 6 of workspaces/module-build/stages/01-spec/output/standard-map-spec.md.
# Each helper returns the worst ratio (error / bound); a criterion passes when the ratio is <= 1.

# AC1: |DT - J| <= 8 eps (1 + |K|) entrywise, J from ForwardDiff (all T) and complex step h = 1e-30 (Float64 only).
function sm_ac1_ratios(::Type{T}, cases; complexstep::Bool) where {T}
    ad = 0.0
    cs = 0.0
    for (K, z) in cases
        sys = StandardMap(T(K))
        z0 = sm_cast(T, z)
        tol = 8 * sm_eps(T) * (1 + abs(K))
        DT = tangent(sys, z0)
        JA = ForwardDiff.jacobian(v -> forward(sys, v), z0)
        ad = max(ad, Float64(maximum(abs, DT - JA)) / tol)
        if complexstep
            h = 1.0e-30
            c1 = forward(sys, z0 + SVector(im * h, 0.0)) ./ h
            c2 = forward(sys, z0 + SVector(0.0, im * h)) ./ h
            JC = SMatrix{2, 2}(imag(c1[1]), imag(c1[2]), imag(c2[1]), imag(c2[2]))
            cs = max(cs, Float64(maximum(abs, DT - JC)) / tol)
        end
    end
    return (ad, cs)
end

# AC2: |det DT - 1| and max |DT' Omega DT - Omega| <= 8 eps (1 + |K|), via the module residuals and a plain recomputation.
function sm_ac2_ratios(::Type{T}, cases) where {T}
    worst = 0.0
    Om = sm_omega(T)
    for (K, z) in cases
        sys = StandardMap(T(K))
        DT = tangent(sys, sm_cast(T, z))
        tol = 8 * sm_eps(T) * (1 + abs(K))
        r = max(
            Float64(det_residual(DT)),
            Float64(symplectic_residual(DT)),
            Float64(abs(DT[1, 1] * DT[2, 2] - DT[1, 2] * DT[2, 1] - 1)),
            Float64(maximum(abs, DT' * Om * DT - Om)),
        )
        worst = max(worst, r / tol)
    end
    return worst
end

# AC3: n forward then n inverse (and inverse then forward) return to z0 within tau(n), n in {1, n_H(K)}.
function sm_ac3_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    for (K, z) in cases
        sys = StandardMap(T(K))
        z0 = sm_cast(T, z)
        for n in unique((1, sm_nH(K)))
            fw = sm_orbit(sys, z0, n)
            back = sm_orbit(sys, fw[end], n; f = inverse)
            worst = max(worst, sm_absdiff(back[end], z0) / sm_tau(sm_eps(T), K, fw, n))
            iv = sm_orbit(sys, z0, n; f = inverse)
            fwd = sm_orbit(sys, iv[end], n)
            worst = max(worst, sm_absdiff(fwd[end], z0) / sm_tau(sm_eps(T), K, iv, n))
        end
    end
    return worst
end

# AC4: max |DT(z) DT_inv(T z) - I| <= 8 eps (1 + |K|)^2.
function sm_ac4_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    for (K, z) in cases
        sys = StandardMap(T(K))
        z0 = sm_cast(T, z)
        E = tangent(sys, z0) * inverse_tangent(sys, forward(sys, z0)) - I
        worst = max(worst, Float64(maximum(abs, E)) / (8 * sm_eps(T) * (1 + abs(K))^2))
    end
    return worst
end

# AC5: wrap(lifted orbit) vs the wrapped-evolution orbit (wrap after every step), torus distance <= tau(k), k <= n_H(K).
function sm_ac5_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    for (K, z) in cases
        sys = StandardMap(T(K))
        z0 = sm_cast(T, z)
        nH = sm_nH(K)
        lifted = sm_orbit(sys, z0, nH)
        zw = wrap(z0)
        for k in 1:nH
            zw = wrap(forward(sys, zw))
            worst = max(worst, sm_tdist(zw, wrap(lifted[k + 1])) / sm_tau(sm_eps(T), K, lifted, k))
        end
    end
    return worst
end

# AC6: T_{-K}(x + pi, p) = (x' + pi, p'); shift back by pi built in T and compare on the torus, k <= n_H(K).
function sm_ac6_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    piT = exact(T, pi)
    for (K, z) in cases
        sysP = StandardMap(T(K))
        sysM = StandardMap(T(-K))
        z0 = sm_cast(T, z)
        nH = sm_nH(K)
        zA = z0
        zB = SVector(z0[1] + piT, z0[2])
        lifted = sm_orbit(sysP, z0, nH)
        for k in 1:nH
            zA = forward(sysP, zA)
            zB = forward(sysM, zB)
            back = SVector(zB[1] - piT, zB[2])
            worst = max(worst, sm_tdist(back, zA) / sm_tau(sm_eps(T), K, lifted, k))
        end
    end
    return worst
end

@testset "AC1 tangent vs derivative" begin
    # DT entries vs ForwardDiff and complex step h=1e-30, bound 8 eps (1+|K|); Float64: 200 points per listed K plus 200 with K uniform, Double64: 40 per K (D38: the 1000 / 200 sizes run in verify/standard-map.jl)
    cases = sm_cases(StableRNG(SM_SEED + 1001), 200)
    ad64, cs64 = sm_ac1_ratios(Float64, cases; complexstep = true)
    @test ad64 <= 1
    @test cs64 <= 1
    # Complex{Double64} is not assumed available, so Double64 uses ForwardDiff only
    ad_dd, _ = sm_ac1_ratios(Double64, sm_cases(StableRNG(SM_SEED + 1002), 40); complexstep = false)
    @test ad_dd <= 1
end

@testset "AC2 determinant and symplecticity" begin
    # |det DT - 1|, max|DT' Omega DT - Omega| <= 8 eps (1+|K|) in Float64 and Double64; 200 points per K (Double64: 40 per K; D38: 1000 / 200 in verify)
    @test sm_ac2_ratios(Float64, sm_cases(StableRNG(SM_SEED + 2001), 200)) <= 1
    @test sm_ac2_ratios(Double64, sm_cases(StableRNG(SM_SEED + 2002), 40)) <= 1
end

@testset "AC3 inverse round trip" begin
    # forward^n then inverse^n (and reverse order) within tau(n)=16 eps M_n n Gamma_n, n in {1, n_H}, 40 orbits per listed K (D38: 200 in verify)
    @test sm_ac3_ratio(Float64, sm_cases(StableRNG(SM_SEED + 3001), 40; uniform = false)) <= 1
    @test sm_ac3_ratio(Double64, sm_cases(StableRNG(SM_SEED + 3002), 40; uniform = false)) <= 1
end

@testset "AC4 inverse tangent" begin
    # max|DT(z) DT_inv(T z) - I| <= 8 eps (1+|K|)^2, 200 points per K (Double64: 40 per K; D38: 1000 / 200 in verify)
    @test sm_ac4_ratio(Float64, sm_cases(StableRNG(SM_SEED + 4001), 200)) <= 1
    @test sm_ac4_ratio(Double64, sm_cases(StableRNG(SM_SEED + 4002), 40)) <= 1
end

@testset "AC5 lifted vs wrapped" begin
    # wrapped(lifted orbit) vs wrap-after-each-step orbit, torus distance <= tau(k) for k <= n_H, 40 orbits per listed K (D38: 200 in verify)
    @test sm_ac5_ratio(Float64, sm_cases(StableRNG(SM_SEED + 5001), 40; uniform = false)) <= 1
    @test sm_ac5_ratio(Double64, sm_cases(StableRNG(SM_SEED + 5002), 40; uniform = false)) <= 1
end

@testset "AC6 K to -K symmetry" begin
    # T_{-K}(x+pi, p) = (x'+pi, p'); torus distance after shifting back <= tau(k), k <= n_H, 40 orbits per listed K (D38: 200 in verify)
    @test sm_ac6_ratio(Float64, sm_cases(StableRNG(SM_SEED + 6001), 40; uniform = false)) <= 1
    @test sm_ac6_ratio(Double64, sm_cases(StableRNG(SM_SEED + 6002), 40; uniform = false)) <= 1
end
