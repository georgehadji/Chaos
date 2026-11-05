# Map invariants: criteria 1 to 4 of workspaces/module-build/stages/01-spec/output/coupled-maps-spec.md, plus the API and
# scope checks (dimension, constructor validation, non-finite input). Each helper returns the worst ratio (error / bound).

# ---- API and scope ----
@testset "API dimension and constructors" begin
    @test dimension(FroeschleMap(0.5, 0.1, 1.0e-3)) == 4
    @test dimension(WrappedFroeschleMap(0.5, 0.1, 1.0e-3)) == 4
    for N in (2, 3, 10, 100)
        @test dimension(StandardMapChain(N, 3.0, 0.2)) == 2N
        @test dimension(WrappedStandardMapChain(N, 3.0, 0.2)) == 2N
    end
    @test dimension(StandardMapChain(3, 3.0, 0.2; boundary = :ring)) == 6
    @test FroeschleMap(0.5, 0.1, 1.0e-3) isa AbstractMap
    @test StandardMapChain(2, 3.0, 0.2) isa AbstractMap
    @test WrappedFroeschleMap(0.5, 0.1, 1.0e-3) isa AbstractMap
    @test WrappedStandardMapChain(2, 3.0, 0.2) isa AbstractMap
    # N = 1 is Y1, N above 100 is out of range, a ring needs N >= 3, the boundary is :open or :ring (reading R5 of the test plan)
    @test_throws ArgumentError StandardMapChain(1, 3.0, 0.2)
    @test_throws ArgumentError StandardMapChain(101, 3.0, 0.2)
    @test_throws ArgumentError StandardMapChain(2, 3.0, 0.2; boundary = :ring)
    @test_throws ArgumentError StandardMapChain(3, 3.0, 0.2; boundary = :periodic)
    @test_throws ArgumentError WrappedStandardMapChain(2, 3.0, 0.2; boundary = :ring)
    # accessible parameters (documented fields of the structs)
    @test FroeschleMap(0.5, 0.1, 1.0e-3).nu == 0.5
    @test FroeschleMap(0.5, 0.1, 1.0e-3).kappa == 0.1
    @test FroeschleMap(0.5, 0.1, 1.0e-3).mu == 1.0e-3
    @test StandardMapChain(4, 3.0, 0.2).N == 4
end

@testset "SC1 non-finite input" begin
    # scope: non-finite input gives non-finite output and is not an error (spec Scope, Domain)
    s3 = FroeschleMap(0.5, 0.1, 1.0e-3)
    s4 = StandardMapChain(3, 3.0, 0.2)
    for bad in (NaN, Inf, -Inf), j in (1, 3)
        z3 = SVector(0.3, 0.4, 0.5, 0.6)
        z3 = Base.setindex(z3, bad, j)
        @test !all(isfinite, forward(s3, z3))
        @test !all(isfinite, tangent(s3, z3))
        @test !all(isfinite, inverse(s3, z3))
        z4 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        z4[j] = bad
        @test !all(isfinite, forward(s4, z4))
        @test !all(isfinite, tangent(s4, z4))
        @test !all(isfinite, inverse(s4, z4))
        w4 = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        w4[j] = bad
        @test !all(isfinite, forward(WrappedStandardMapChain(3, 3.0, 0.2), w4))
    end
end

# ---- AC1: |M - J| <= 16 eps (1 + a_max) entrywise; J from ForwardDiff (all T) and complex step h = 1e-30 (Float64 only) ----
function cm_ac1_ratios(::Type{T}, cases; complexstep::Bool) where {T}
    ad, cs, rf = 0.0, 0.0, 0.0
    for (c, z) in cases
        sys = cm_sys(T, c)
        z0 = cm_state(T, c, z)
        d = cm_d(c)
        tol = 16 * cm_eps(T) * (1 + cm_amax(c))
        M = Matrix(tangent(sys, z0))
        JA = Matrix(ForwardDiff.jacobian(v -> forward(sys, v), z0))
        ad = max(ad, Float64(maximum(abs, M - JA)) / tol)
        MR = cm_ref_tangent(c, cm_params(T, c), z0)
        rf = max(rf, Float64(maximum(abs, M - MR)) / tol)
        if complexstep
            h = 1.0e-30
            JC = zeros(Float64, d, d)
            for k in 1:d
                JC[:, k] = imag.(forward(sys, cm_perturb(z0, k, im * h))) ./ h
            end
            cs = max(cs, Float64(maximum(abs, M - JC)) / tol)
        end
    end
    return (ad, cs, rf)
end

@testset "AC1 tangent vs derivative" begin
    # entries of M vs ForwardDiff, complex step (Float64) and the plain-formula tangent, bound 16 eps (1 + a_max);
    # Y3 sets and uniform parameters, Y4 panel N in {2, 3, 5, 10}; Float64 10 points per set, Double64 3 (D38-style: 1000 / 200 in verify)
    ad, cs, rf = cm_ac1_ratios(Float64, cm_cases(StableRNG(CM_SEED + 1001), 10); complexstep = true)
    @test ad <= 1
    @test cs <= 1
    @test rf <= 1
    ad, _, rf = cm_ac1_ratios(Double64, cm_cases(StableRNG(CM_SEED + 1002), 3; maxN = 5); complexstep = false)
    @test ad <= 1
    @test rf <= 1
end

# ---- AC2: inverse tangent ----
function cm_ac2_ratios(::Type{T}, cases) where {T}
    prod_r, der_r, ref_r = 0.0, 0.0, 0.0
    for (c, z) in cases
        sys = cm_sys(T, c)
        z0 = cm_state(T, c, z)
        zi = forward(sys, z0)
        a = 1 + cm_amax(c)
        M, Mi = Matrix(tangent(sys, z0)), Matrix(inverse_tangent(sys, zi))
        prod_r = max(prod_r, Float64(maximum(abs, M * Mi - I)) / (16 * cm_eps(T) * a^2))
        JI = Matrix(ForwardDiff.jacobian(v -> inverse(sys, v), zi))
        der_r = max(der_r, Float64(maximum(abs, Mi - JI)) / (16 * cm_eps(T) * a))
        MR = cm_ref_inverse_tangent(c, cm_params(T, c), zi)
        ref_r = max(ref_r, Float64(maximum(abs, Mi - MR)) / (16 * cm_eps(T) * a))
    end
    return (prod_r, der_r, ref_r)
end

@testset "AC2 inverse tangent" begin
    # max|M(z) M_inv(T z) - I| <= 16 eps (1 + a_max)^2; M_inv vs ForwardDiff of the inverse map and vs the plain formulas <= 16 eps (1 + a_max)
    for (T, m, maxN, seed) in ((Float64, 10, 10, 2001), (Double64, 3, 5, 2002))
        a, b, r = cm_ac2_ratios(T, cm_cases(StableRNG(CM_SEED + seed), m; maxN))
        @test a <= 1
        @test b <= 1
        @test r <= 1
    end
end

# ---- AC3: determinant and symplecticity, evaluated in 256 bits from the rounded matrix ----
function cm_ac3_ratios(::Type{T}, cases) where {T}
    wdet, wsym, wmod, wsymm = 0.0, 0.0, 0.0, 0.0
    for (c, z) in cases
        sys = cm_sys(T, c)
        z0 = cm_state(T, c, z)
        d = cm_d(c)
        M = Matrix(tangent(sys, z0))
        h = d ÷ 2
        wsymm = max(wsymm, M[(h + 1):d, 1:h] == M[(h + 1):d, 1:h]' ? 0.0 : 1.0)   # the lower-left block (A or G) is symmetric to the bit
        setprecision(BigFloat, CM_BITS) do
            Mb = BigFloat[cm_big(M[i, j]) for i in 1:d, j in 1:d]
            Mi = inv(Mb)
            dres = Float64(abs(det(Mb) - 1))
            fro = Float64(sqrt(sum(abs2, Mb)) * sqrt(sum(abs2, Mi)))
            wdet = max(wdet, dres / (8 * d * cm_eps(T) * fro))
            Om = cm_omega(BigFloat, d)
            sres = Float64(maximum(abs, Mb' * Om * Mb - Om))
            wsym = max(wsym, sres / (8 * d * cm_eps(T) * Float64(maximum(abs, Mb))^2))
            # the module residual functions (evaluated in T) stay within a generous multiple of the same scales
            wmod = max(wmod, Float64(det_residual(c.kind === :y3 ? SMatrix{4, 4}(M) : M)) / (64 * d * cm_eps(T) * fro))
        end
        wmod = max(wmod, Float64(symplectic_residual(c.kind === :y3 ? SMatrix{4, 4}(M) : M)) / (32 * d * cm_eps(T) * Float64(maximum(abs, M))^2))
    end
    return (wdet, wsym, wmod, wsymm)
end

@testset "AC3 determinant and symplecticity" begin
    # |det M - 1| <= 8 d eps ||M||_F ||M^-1||_F and |M' Omega M - Omega| <= 8 d eps ||M||_max^2, 256-bit evaluation, Y3 and Y4 N <= 10;
    # 10 points per set Float64, 3 Double64; module det_residual and symplectic_residual (also on a dense Matrix) within a multiple of the same scales
    for (T, m, seed) in ((Float64, 10, 3001), (Double64, 3, 3002))
        a, b, r, sy = cm_ac3_ratios(T, cm_cases(StableRNG(CM_SEED + seed), m))
        @test a <= 1
        @test b <= 1
        @test r <= 1
        @test sy == 0
    end
end

# ---- AC4: inverse round trip ----
function cm_ac4_ratio(::Type{T}, cases) where {T}
    worst = 0.0
    for (c, z) in cases
        sys = cm_sys(T, c)
        z0 = cm_state(T, c, z)
        for n in unique((1, cm_nH(c)))
            fw = cm_orbit(sys, z0, n)
            back = cm_orbit(sys, fw[end], n; f = inverse)
            worst = max(worst, cm_absdiff(back[end], z0) / cm_tau(cm_eps(T), c, fw, n))
            iv = cm_orbit(sys, z0, n; f = inverse)
            fwd = cm_orbit(sys, iv[end], n)
            worst = max(worst, cm_absdiff(fwd[end], z0) / cm_tau(cm_eps(T), c, iv, n))
        end
    end
    return worst
end

@testset "AC4 inverse round trip" begin
    # n forward then n inverse (and inverse then forward) return to the start within tau(n), n in {1, n_H}; 10 orbits per set Float64, 4 Double64 (200 in verify)
    @test cm_ac4_ratio(Float64, cm_cases(StableRNG(CM_SEED + 4001), 10)) <= 1
    @test cm_ac4_ratio(Double64, cm_cases(StableRNG(CM_SEED + 4002), 4; maxN = 5)) <= 1
    # the module orbit agrees with the independent plain-formula orbit and inverse step (256-bit, n_H steps, one orbit per set)
    setprecision(BigFloat, CM_BITS) do
        worst = 0.0
        for (c, z) in cm_cases(StableRNG(CM_SEED + 4003), 1; uniform = false)
            nH = cm_nH(c)
            sysb = cm_sys(BigFloat, c)
            z0 = cm_state(BigFloat, c, z)
            a = cm_orbit(sysb, z0, nH)
            b = cm_ref_orbit(c, cm_params(BigFloat, c), z0, nH)
            worst = max(worst, cm_absdiff(a[end], b[end]) / cm_tau(ldexp(1.0, 1 - CM_BITS), c, a, nH))
            ai = inverse(sysb, a[end])
            bi = cm_ref_inverse(c, cm_params(BigFloat, c), a[end])
            worst = max(worst, cm_absdiff(ai, bi) / cm_tau(ldexp(1.0, 1 - CM_BITS), c, a, nH))
        end
        @test worst <= 1
    end
end
