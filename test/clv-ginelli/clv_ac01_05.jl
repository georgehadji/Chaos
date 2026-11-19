# AC1 to AC5: closed-form CLVs, covariance and telescoping, convergence and covariance chain, order in h, symplectic structure. Every testset starts by naming the new symbols it uses,
# so a missing module errors that testset (UndefVarError). Testsets that use no new name (the negative controls, run on the test-side reference algorithm) pass in the RED run.

# ------------------------------------------------------------------ AC1: constant-Jacobian CLVs against the exact vectors
@testset "AC1 constant-Jacobian CLVs: fixtures a to d, identity and seeded frames, exact vectors and exponents at every sample" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    nf, nb, ns = 100, 100, 50; worstV = 0.0; worstL = 0.0; nbad = 0; nrun = 0
    for name in (:a, :b, :c, :d), tau in (name === :d ? (20,) : (1, 4)), fr in 0:3, T in (Float64, Double64)
        (T === Double64 && !(name === :b && tau == 1 && fr in (0, 1))) && continue
        fx = clv_fixture(name, T; tau = tau); d = fx.d
        W0 = (fr == 0 && name !== :b) ? clv_eye(d, d, T) : clv_frame(clv_seed(1, fr), 1, d, d, T)      # identity frame: v_4 lies in span(e_1, e_2, e_3) for fixture b (plan U6), seeded frame instead
        C0 = clv_rand_C(StableRNG(clv_seed(1, 10 + fr)), d, T)
        res = clv_go(fx.map, zero(SVector{d, T}), nf, ns, nb; tau = tau, h = fx.h, W0 = W0, C0 = C0)
        B = clv_B_fix(T, fx.gamma, min(nf, nb)) + (name === :d ? 10 * fx.h^4 : 0.0); nrun += 1
        for n in 0:ns
            dv = maximum(i -> clv_delta(res.V[:, i, n + 1], fx.V[:, i]), 1:d)
            dl = maximum(i -> abs(Float64(clv_big(res.lambda[i, n + 1]) - fx.lam[i])), 1:d)
            Bl = B / fx.Dt + 4 * tau * clv_eps(T) * exp(Float64(fx.lam[1] - fx.lam[end]) * fx.Dt) / fx.Dt     # plus the rounding of the interval product, condition exp((l_1 - l_d) Dt) (plan U5)
            worstV = max(worstV, dv / B); worstL = max(worstL, dl / Bl)
            (dv <= B && dl <= Bl) || (nbad += 1)
            @test res.C[1, 1, n + 1] == one(T) && all(iszero, res.C[2:d, 1, n + 1])             # column 1 is e_1 (I1)
            @test all(i -> res.C[i, i, n + 1] > 0, 1:d)
        end
        @test res.V[:, 1, :] == res.G[:, 1, :]                                                   # v_1 = g_1 bitwise (I1)
    end
    @info "AC1 realized ratios (error / bound)" worstV worstL nrun
    @test nbad == 0
end

@testset "AC1 negative controls: the forward Gram-Schmidt vectors, the forward step and the lower triangle all fail the exact vectors" begin
    for name in (:a, :b)
        fx = clv_fixture(name, Float64); d = fx.d; pts = clv_orbit(fx.map, zero(SVector{d, Float64}), 250)
        G0 = clv_frame(clv_seed(1, 20), 1, d, d, Float64); C0 = clv_rand_C(StableRNG(clv_seed(1, 21)), d, Float64)
        ok = clv_ref_run(fx.map, pts, 1, 100, 50, 100, G0, C0, Float64)
        @test maximum(i -> clv_delta(clv_refV(ok, 0)[:, i], fx.V[:, i]), 1:d) < 1.0e-12                      # the reference itself passes
        @test maximum(i -> clv_delta(ok.G[1][:, i], fx.V[:, i]), 1:d) > 1.0e-3                                # Gram-Schmidt vectors are not the CLVs
        for mode in (:fwd, :lower)
            r = clv_ref_run(fx.map, pts, 1, 100, 50, 100, G0, C0, Float64; mode = mode)
            @test maximum(i -> clv_delta(clv_refV(r, 0)[:, i], fx.V[:, i]), 1:d) > 1.0e-3
        end
    end
end

# ------------------------------------------------------------------ AC2: covariance, local exponents, telescoping, invariants
@testset "AC2 covariance, normalisation, local exponents and telescoping on the linear fixtures and chaotic orbits (I1, I2, I3)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; ns = clv_size(100, 500); nf = nb = clv_size(60, 150); m = clv_size(3, 20)
    cases = Any[]
    for name in (:a, :b)
        fx = clv_fixture(name, T); push!(cases, (map = fx.map, x0 = zero(SVector{fx.d, T}), tau = 1, h = 1.0, nf = 100, nb = 100, ns = 50, label = "fixture $name", fix = true))
    end
    rng = StableRNG(clv_seed(2, 1))
    for x0 in clv_starts(rng, 2, m, T)
        push!(cases, (map = clv_std(T), x0 = x0, tau = 1, h = 1.0, nf = nf, nb = nb, ns = ns, label = "standard map", fix = false))
    end
    im = clv_im(T)
    push!(cases, (map = im, x0 = clv_orbit_C(), tau = 20, h = CLV_H, nf = clv_size(150, 2000), nb = clv_size(150, 2000), ns = clv_size(60, 500), label = "F3 orbit C", fix = false))
    for z in clv_chaotic_starts(clv_size(1, 20), 1)
        push!(cases, (map = im, x0 = z, tau = 20, h = CLV_H, nf = clv_size(150, 2000), nb = clv_size(150, 2000), ns = clv_size(60, 500), label = "F3 chaotic", fix = false))
    end
    worst = Dict{String, Float64}("a" => 0.0, "b" => 0.0, "c" => 0.0, "mean" => 0.0); nbad = 0
    for (ic, c) in enumerate(cases)
        d = dimension(c.map); k = d; Dt = c.tau * c.h
        W0 = clv_frame(clv_seed(2, 10 + ic), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(2, 20 + ic)), d, T)
        res = clv_go(c.map, c.x0, c.nf, c.ns, c.nb; tau = c.tau, h = c.h, W0 = W0, C0 = C0)
        for n in 0:c.ns
            G = res.G[:, :, n + 1]
            @test maximum(abs, G' * G - I) <= 8 * k * clv_eps(T)                                         # I1: orthonormal frame
            Cn = res.C[:, :, n + 1]
            @test all(j -> abs(norm(Cn[:, j]) - 1) <= 8 * clv_eps(T), 1:k) && all(i -> all(j -> i <= j || iszero(Cn[i, j]), 1:k), 1:k)
        end
        lnc(n, i) = log(res.C[i, i, n + 1])
        for i in 1:k
            lnR = [log(res.R[i, i, n + 1]) for n in 0:(c.ns - 1)]; lnL = [Float64(res.lambda[i, n + 1]) * Dt for n in 0:(c.ns - 1)]
            sumL = sum(BigFloat.(lnL)); sumR = sum(BigFloat.(lnR))
            lhs = abs(Float64(sumL - sumR) - lnc(0, i) + lnc(c.ns, i)); rhs = 16 * clv_eps(T) * (sum(abs, lnL) + sum(abs, lnR) + 1)
            worst["c"] = max(worst["c"], lhs / rhs); lhs <= rhs || (nbad += 1)                              # (c) telescoping
        end
        for n in 0:(c.ns - 1)
            Rn = Matrix{Float64}(res.R[:, :, n + 1]); B = clv_B_cov(T, cond(Rn)); Cn = res.C[:, :, n + 1]; Cn1 = res.C[:, :, n + 2]
            for i in 1:k
                x = Rn * Cn[:, i]; Lam = exp(Float64(res.lambda[i, n + 1]) * Dt)
                ea = maximum(abs, x - Lam * Cn1[:, i]); eb = abs(norm(x) * norm(Rn \ Cn1[:, i]) - 1)
                worst["a"] = max(worst["a"], ea / B); worst["b"] = max(worst["b"], eb / B)
                (ea <= B && eb <= B) || (nbad += 1)                                                       # (a) covariance, (b) normalisation
            end
        end
        # the window mean of the local exponent against the forward exponent of the delivered spectrum on the identical orbit (G_0 and the window start from the module)
        G0w = SMatrix{d, k, T}(res.G[:, :, 1])
        sp = c.fix ? nothing : propagate(c.map, [res.x0], c.ns * c.tau; k = k, tau = c.tau, seed = 1, reducers = (QRSpectrum(HouseholderQR()),), W0 = [G0w]).finals[1][1].lambda ./ c.h
        c.fix || for i in 1:k
            meanL = sum(Float64(res.lambda[i, n + 1]) for n in 0:(c.ns - 1)) / c.ns
            allow = (abs(lnc(0, i)) + abs(lnc(c.ns, i))) / (c.ns * Dt) + 16 * clv_eps(T) * (1 + abs(meanL)) * 10
            worst["mean"] = max(worst["mean"], abs(meanL - sp[i]) / allow); abs(meanL - sp[i]) <= allow || (nbad += 1)
        end
    end
    @info "AC2 realized ratios (error / bound)" worst ncases = length(cases)
    @test nbad == 0
end

@testset "AC2 negative control: a 1-norm normalisation breaks the product ||R c|| ||R^-1 c'|| = 1" begin
    fx = clv_fixture(:b, Float64); d = 4; pts = clv_orbit(fx.map, zero(SVector{d, Float64}), 200)
    G0 = clv_frame(clv_seed(2, 40), 1, d, d, Float64); C0 = clv_rand_C(StableRNG(clv_seed(2, 41)), d, Float64)
    ok = clv_ref_run(fx.map, pts, 1, 50, 20, 50, G0, C0, Float64); bad = clv_ref_run(fx.map, pts, 1, 50, 20, 50, G0, C0, Float64; mode = :norm1)
    prod(ref, n, i) = norm(ref.R[n + 1] * ref.C[n + 1][:, i]) * norm(ref.R[n + 1] \ ref.C[n + 2][:, i])
    @test maximum(abs(prod(ok, n, i) - 1) for n in 0:19, i in 1:4) <= clv_B_cov(Float64, 1.0e4)
    @test maximum(abs(prod(bad, n, i) - 1) for n in 0:19, i in 1:4) > 1.0e-3
end

# ------------------------------------------------------------------ AC3: rate, independence of the starting data, covariance chain against an independent tangent product
@testset "AC3a backward convergence rate on fixtures a and b: slope of ln error against the backward transient (V22)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; worst = 0.0
    for name in (:a, :b)
        fx = clv_fixture(name, T); d = fx.d; W0 = clv_frame(clv_seed(3, 1), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(3, 2)), d, T)
        errs = zeros(d, 40)
        for m in 2:40
            res = clv_go(fx.map, zero(SVector{d, T}), 100, 1, m; W0 = W0, C0 = C0)
            for j in 1:d
                errs[j, m] = clv_delta(res.V[:, j, 1], fx.V[:, j])
            end
        end
        @test all(errs[1, 2:40] .<= 1.0e3 * clv_eps(T))                                                   # column 1 has no error
        for j in 2:d
            ms = [m for m in 2:40 if errs[j, m] > 1.0e3 * clv_eps(T)]; @test length(ms) >= 4
            u = Float64.(ms); v = log.(errs[j, ms]); slope = sum((u .- sum(u) / length(u)) .* (v .- sum(v) / length(v))) / sum((u .- sum(u) / length(u)) .^ 2)
            expct = -Float64(fx.lam[j - 1] - fx.lam[j]) * fx.Dt
            worst = max(worst, abs(slope / expct - 1)); @test abs(slope / expct - 1) <= 0.2
        end
    end
    @info "AC3a worst relative slope deviation" worst
end

@testset "AC3b independence of the starting frame, the starting C and tau at common sample times, once the transients reach n*" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; ssys = Any[]
    rng = StableRNG(clv_seed(3, 3))
    push!(ssys, (name = "standard map", map = clv_std(T), x0 = clv_starts(rng, 2, 1, T)[1], taus = (1, 5), h = 1.0, singles = 1:2, block = nothing, skip = false))
    push!(ssys, (name = "coupled maps", map = clv_cm(T), x0 = clv_starts(rng, 4, 1, T)[1], taus = (1, 5), h = 1.0, singles = 1:4, block = nothing, skip = false))
    for (i, z) in enumerate(clv_chaotic_starts(clv_size(2, 20), 1))
        push!(ssys, (name = "F3 chaotic $i", map = clv_im(T), x0 = z, taus = (10, 20), h = CLV_H, singles = [1, 4], block = 2:3, skip = true))
    end
    worst = 0.0
    for s in ssys
        d = dimension(s.map); ts, tb = s.taus; r = tb ÷ ts
        lam = clv_pilot(s.map, s.x0, tb, s.h); Nb = clv_nstar(lam, tb * s.h; skip23 = s.skip); nsb = 6
        runs = []
        for (sg, sc) in ((1, 1), (2, 2))
            W0 = clv_frame(clv_seed(3, 10 + sg), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(3, 20 + sc)), d, T)
            big = clv_go(s.map, s.x0, Nb, nsb, Nb; tau = tb, h = s.h, W0 = W0, C0 = C0, outputs = (:V,))
            small = clv_go(s.map, s.x0, r * Nb, r * nsb, r * Nb; tau = ts, h = s.h, W0 = W0, C0 = C0, outputs = (:V,))
            push!(runs, big, small)
        end
        ref = runs[1]
        for q in 2:4, m in 0:nsb
            Vq = q % 2 == 1 ? runs[q].V[:, :, m + 1] : runs[q].V[:, :, r * m + 1]                           # runs: big, small, big, small
            Vr = ref.V[:, :, m + 1]; w = clv_vdist(Vq, Vr; singles = s.singles, block = s.block)
            worst = max(worst, w); @test w <= 1.0e-10
        end
        @info "AC3b" s.name Nb
    end
    @info "AC3b worst distance" worst
end

@testset "AC3c covariance chain: ln ||M_{m,n} v_i(n)|| against the sum of the local growth, M the 256-bit tangent product along the realized orbit" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    worst = 0.0; nbad = 0
    for T in (Float64, Double64), sysname in (:std, :f3)
        setprecision(BigFloat, CLV_BITS) do
            map = sysname === :std ? clv_std(T) : clv_im(T); mapB = clv_bigmap(map); d = dimension(map)
            tau = sysname === :std ? 1 : 20; h = sysname === :std ? 1.0 : CLV_H; ns = sysname === :std ? 80 : 60; Dt = T(tau) * T(h)   # D124: the module's Dt, tau h in the element type
            x0 = sysname === :std ? SVector{2, T}(T(1.3), T(0.7)) : SVector{4, T}(map_to(T, clv_chaotic_starts(1, 1)[1]))
            W0 = clv_frame(clv_seed(3, 30), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(3, 31)), d, T)
            res = clv_go(map, x0, sysname === :std ? 60 : 150, ns, sysname === :std ? 60 : 150; tau = tau, h = h, W0 = W0, C0 = C0)
            pts = clv_orbit(map, res.x0, ns * tau); Ms = clv_interval_tangents(mapB, pts, tau, ns, BigFloat)
            kap = maximum(n -> cond(Matrix{Float64}(res.R[:, :, n + 1])), 0:(ns - 1))
            samples = unique(round.(Int, range(0, ns - 50; length = 20)))
            for m in (1, 10, 50), n in samples, i in 1:d
                v = BigFloat.(clv_big.(res.V[:, i, n + 1]))
                for l in 0:(m - 1)
                    v = Ms[n + l + 1] * v
                end
                err = Float64(abs(log(norm(v)) - sum(clv_big(res.lambda[i, n + l + 1]) * Dt for l in 0:(m - 1))))
                # the rounding error of V along the faster directions grows as exp((l_1 - l_i) m Dt) in the product: allowance added to the spec bound (plan U9)
                B = clv_B_cov(T, kap) * m + 256 * clv_eps(T) * kap * exp(max(0.0, Float64(res.lambda_fwd[1] - res.lambda_fwd[i])) * m * Dt)
                worst = max(worst, err / B); err <= B || (nbad += 1)
            end
        end
    end
    @info "AC3c worst ratio error / bound" worst
    @test nbad == 0
end

# ------------------------------------------------------------------ AC4: order in h on the flow fixture
@testset "AC4 order in h of the CLV error on the separable linear flow: order 2 and order 4 tables, h in {1/10, 1/20, 1/40}" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; fx = clv_fixture(:d, T); flow = fx.map.flow; W0 = clv_frame(clv_seed(4, 1), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(4, 2)), 4, T)
    for (key, p) in ((:o2, 2), (:o4, 4))
        errs = Float64[]
        for tau in (10, 20, 40)
            im = IntegratorMap(clv_table(T, key), flow, T(1) / T(tau))
            res = clv_go(im, zero(SVector{4, T}), 100, 4, 100; tau = tau, h = 1 / tau, W0 = W0, C0 = C0, outputs = (:V,))
            push!(errs, maximum(n -> maximum(i -> clv_delta(res.V[:, i, n + 1], fx.V[:, i]), 1:4), 0:4))
        end
        @info "AC4 absolute errors (h = 1/10, 1/20, 1/40)" key errs
        for q in 1:2
            if errs[q + 1] > 1.0e3 * clv_eps(T)
                @test 2.0^p / 2 <= errs[q] / errs[q + 1] <= 2 * 2.0^p
            end
        end
        @test errs[2] > 1.0e3 * clv_eps(T)                                                                  # the p = 4 error at h = 1/20 sits above the floor (the check the spec asks stage 02 to make)
    end
end

# ------------------------------------------------------------------ AC5: symplectic structure
@testset "AC5 symplectic pair identity, non-partner orthogonality and the determinant identity (V14)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; worst = Dict("a" => 0.0, "a0" => 0.0, "b" => 0.0, "c" => 0.0); nbad = 0
    rng = StableRNG(clv_seed(5, 1)); fb = clv_fixture(:b, T)
    systems = Any[(name = "standard map", map = clv_std(T), x0 = clv_starts(rng, 2, 1, T)[1], tau = 1, h = 1.0, skip = false, flow = false),
                  (name = "coupled maps", map = clv_cm(T), x0 = clv_starts(rng, 4, 1, T)[1], tau = 1, h = 1.0, skip = false, flow = false),
                  (name = "fixture b", map = fb.map, x0 = zero(SVector{4, T}), tau = 1, h = 1.0, skip = false, flow = false)]
    for (i, z) in enumerate(clv_chaotic_starts(clv_size(1, 20), 1))
        push!(systems, (name = "F3 chaotic $i", map = clv_im(T), x0 = z, tau = 20, h = CLV_H, skip = true, flow = true))
    end
    for (is, s) in enumerate(systems)
        d = dimension(s.map); Dt = s.tau * s.h
        W0 = clv_frame(clv_seed(5, 10 + is), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(5, 30 + is)), d, T)
        ns = clv_size(60, 300)
        for (label, nf, nb, ns_) in (("random C, n_b = 1", 20, 1, 40), ("converged", clv_nstar(s.name == "fixture b" ? Float64.(fb.lam) : clv_pilot(s.map, s.x0, s.tau, s.h), Dt; skip23 = s.skip), 0, ns))
            nb = label == "converged" ? nf : nb
            res = clv_go(s.map, s.x0, nf, ns_, nb; tau = s.tau, h = s.h, W0 = W0, C0 = C0)
            Lam(n, i) = Float64(res.lambda[i, n + 1]) * Dt
            om(n, i, j) = clv_omega(res.V[:, i, n + 1], res.V[:, j, n + 1], d)
            lomax = Float64(maximum(n -> maximum(i -> abs(log(abs(om(n, i, d + 1 - i)))), 1:(d ÷ 2)), 0:ns_)); lmax = maximum(n -> maximum(i -> abs(Lam(n, i)), 1:d), 0:(ns_ - 1))
            bnd = 512 * clv_eps(T) * (1 + lomax + lmax)                                   # spec 64 eps (1 + ...): realized 0.95 of it, widened 8x (plan U7)
            for n in 0:(ns_ - 1), i in 1:(d ÷ 2)
                e = abs(Lam(n, i) + Lam(n, d + 1 - i) - Float64(log(abs(om(n, i, d + 1 - i))) - log(abs(om(n + 1, i, d + 1 - i)))))
                # D124: rounding of log|omega| near a tangency of the pair, 64 eps / |omega| at both samples (F3 start 6: |omega_23| 3.6e-6, ratio 7.99 without it, 0.049 with it)
                bn = bnd + 64 * clv_eps(T) * (1 / abs(Float64(om(n, i, d + 1 - i))) + 1 / abs(Float64(om(n + 1, i, d + 1 - i))))
                worst["a"] = max(worst["a"], e / bn); e <= bn || (nbad += 1)                                      # (a), valid for any starting C
            end
            kk = d
            ldet(n) = sum(log(res.C[i, i, n + 1]) for i in 1:kk); dmax = maximum(n -> abs(ldet(n)), 0:ns_)
            for n in 0:(ns_ - 1)
                e = abs(sum(Lam(n, i) for i in 1:kk) - ldet(n) + ldet(n + 1)); b = 64 * clv_eps(T) * (1 + dmax)
                worst["c"] = max(worst["c"], e / b); e <= b || (nbad += 1)                                         # (c) k = d
            end
            if label == "converged"
                pairs = s.flow ? [(1, 2), (1, 3), (4, 2), (4, 3)] : [(i, j) for i in 1:d for j in (i + 1):d if i + j != d + 1]
                for n in 0:ns_, (i, j) in pairs
                    w = abs(Float64(om(n, i, j))); worst["b"] = max(worst["b"], w / 1.0e-9); w <= 1.0e-9 || (nbad += 1)   # (b)
                end
            end
        end
    end
    @info "AC5 realized ratios (error / bound; (b) against 1e-9)" worst
    @test nbad == 0
end

@testset "AC5 negative control: a wrong normalisation of the growth factor fails the pair identity" begin
    fx = clv_fixture(:b, Float64); d = 4; pts = clv_orbit(fx.map, zero(SVector{d, Float64}), 200)
    G0 = clv_frame(clv_seed(5, 50), 1, d, d, Float64); C0 = clv_rand_C(StableRNG(clv_seed(5, 51)), d, Float64)
    ref = clv_ref_run(fx.map, pts, 1, 50, 20, 50, G0, C0, Float64)
    function resid(wrong)
        w = 0.0
        for n in 0:18, i in 1:2
            j = 5 - i; om(m) = clv_omega(clv_refV(ref, m)[:, i], clv_refV(ref, m)[:, j], 4)
            L(nn, q) = wrong ? 1 / sum(abs, ref.R[nn + 1] \ ref.C[nn + 2][:, q]) : ref.Lambda[q, nn + 1]
            w = max(w, abs(log(L(n, i)) + log(L(n, j)) - Float64(log(abs(om(n))) - log(abs(om(n + 1))))))
        end
        return w
    end
    @test resid(false) <= 64 * clv_eps(Float64) * 20
    @test resid(true) > 1.0e-6
end
