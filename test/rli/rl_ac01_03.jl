# AC1 to AC3: the algebra of the definition on test map T1 (RL9(b)), the exact zero for constant tangents (RL5(d)), the same-orbit BigFloat reference (RL11).
# Every testset starts with a line naming the new symbols it uses (PairMap of the Systems layer, RLI of the Indicators layer), so a missing symbol errors that
# testset only (UndefVarError). The pair is one orbit of the product system on the existing stream with k = 2 and the supplied block frame (RL5, D75).

# Exact sums of test map T1 at the stored start (x0, p0): s(m) = c sum_{j<m} x_j, x_j = x0 + j p0 (BigFloat, exact for the dyadic starts used here), m = 1..N.
function rl_T1_sums(c, x0, p0, N::Int)
    return setprecision(BigFloat, RL_BITS) do
        cb, xb, pb = rl_big(c), rl_big(x0), rl_big(p0)
        s = Vector{BigFloat}(undef, N)
        acc = zero(BigFloat)
        for m in 1:N
            acc += cb * (xb + (m - 1) * pb)
            s[m] = acc
        end
        s
    end
end
# The closed form of RL9(b): RLI(n) = c (a + b (n + tau - 2)/4) and DL_k = c (a + b (t_k - 1)/2), BigFloat.
rl_T1_closed(c, a, b, n::Int, tau::Int) = setprecision(BigFloat, RL_BITS) do
    rl_big(c) * (BigFloat(a) + BigFloat(b) * (n + tau - 2) / 4)
end
# G_first of L1 S3 for T1 at the start (x0, p0): the squared Frobenius norm of the product of the first tau tangents, exp(2 C) + exp(-2 C), C = c sum_{j<tau} x_j.
rl_T1_G(c, x0, p0, tau::Int) = (C = Float64(c) * sum(Float64(x0) + j * Float64(p0) for j in 0:(tau - 1)); exp(2C) + exp(-2C))
# B_R^lin of RL11 (K_s = 0): 8 eps (n + G_first) per orbit, both orbits, times H_K / n, plus 4 eps RLI_ref.
rl_T1_B(epsT, GA, GB, n::Int, tau::Int, rli_ref) = 8 * epsT * (2n + GA + GB) * rl_H(n ÷ tau) / n + 4 * epsT * Float64(abs(rli_ref))

@testset "AC1 algebra of the definition" begin
    _ = (RLI, PairMap)
    # Test map T1 (x' = x + p, p' = p, tangent diag(exp(c x), exp(-c x)) at the pre-step point), c = 0.3, the supplied frames of RL5 with xi_0 = (1, 0) for both orbits,
    # the pair of T1 as the product system. The starts are dyadic (x_0 = 51/512 for 0.1, p_0 = 51/1024 for 0.05, the separations 5/512 for 0.01 and 5/256 for 0.02), so the orbit
    # is exact in the element type and the closed form is exact for the stored numbers; c is the Float64 0.3 widened exactly. (i) (a, b) in {(0.01, 0), (0, 0.02), (0.01, 0.02)}
    # (dyadic), tau in {1, 3, 5}, n in {1, 2, 5, 20, 60, 100} (the handoffs among them): |RLI - c (a + b (n + tau - 2)/4)| <= B_R^lin(n) of RL11 (K_s = 0), also DL_k = c (a + b (t_k - 1)/2);
    # (i') the same with the decimals of the spec (c = 3/10, x_0 = 1/10, p_0 = 1/20, a = 1/100, b = 1/50 rounded once into the element type), where the orbit carries its own rounding. (ii) mixed signs (0.01, -0.02) and (-0.01, 0.02):
    # DL_k and RLI against the BigFloat term-by-term sums of RL4. (iii) RLI(2a, 2b) = 2 RLI(a, b) within 2 B_R^lin. (iv) negative controls: a double with independent vectors
    # (the shadow started on (0, 1)) fails (i); a double that divides by n instead of K fails (i) for tau > 1; one that drops the absolute value fails (ii); one that forms the
    # mean over k = 0..K-1 differs in (i); the printed 1/t prefactor taken literally (mean / tau) differs for tau > 1. The BigFloat post-processor (rl_ref_rli) is checked
    # against the closed form to 1e-60 first.
    nlist = (1, 2, 5, 20, 60, 100)
    x0, p0 = 51 / 512, 51 / 1024
    seps_i = ((5 / 512, 0.0), (0.0, 5 / 256), (5 / 512, 5 / 256))
    seps_ii = ((5 / 512, -5 / 256), (-5 / 512, 5 / 256))
    for T in RL_TYPES
        epsT = rl_eps(T)
        c = T(0.3)
        P = PairMap(rl_T1(c))
        xi = SVector{2, T}(1, 0)
        W = rl_frame(xi)
        # run the module on the pair (x0, p0), (x0 + a, p0 + b): returns Dict n => final, for the handoffs among nlist
        run_pair(a, b, tau, W = W, xs = x0, ps = p0) = begin
            ns = filter(n -> n % tau == 0, collect(nlist))
            z = SVector{4, T}(T(xs), T(ps), T(xs + a), T(ps + b))
            r = rl_run((RLI(),), P, [z], [W], maximum(ns), tau; steps = ns)
            Dict(n => f for (n, f) in zip(ns, r[1][1].snaps))
        end
        ref_pair(a, b, tau, N) = (rl_ref_rli(rl_T1_sums(c, x0, p0, N), rl_T1_sums(c, x0 + a, p0 + b, N), tau),)
        # helper self-check: the post-processor on the exact sums equals the closed form (BigFloat, 1e-60)
        for (a, b) in seps_i, tau in (1, 3, 5)
            N = 60 - 60 % tau
            ref = ref_pair(a, b, tau, N)[1]
            for n in filter(n -> n % tau == 0 && n <= N, collect(nlist))
                cf = rl_T1_closed(c, a, b, n, tau)
                @test abs(ref.RLI[n ÷ tau] - cf) < big"1e-60"
                @test abs(ref.DL[n ÷ tau] - setprecision(() -> rl_big(c) * (BigFloat(a) + BigFloat(b) * (n - 1) / 2), BigFloat, RL_BITS)) < big"1e-60"
            end
        end
        # (i)
        worst = 0.0
        for (a, b) in seps_i, tau in (1, 3, 5)
            fin = run_pair(a, b, tau)
            for (n, f) in fin
                cf = rl_T1_closed(c, a, b, n, tau)
                GA, GB = rl_T1_G(c, x0, p0, tau), rl_T1_G(c, x0 + a, p0 + b, tau)
                B = rl_T1_B(epsT, GA, GB, n, tau, cf)
                e = rl_err(f.RLI, cf)
                @test e <= B
                worst = max(worst, e / B)
                dlcf = setprecision(() -> rl_big(c) * (BigFloat(a) + BigFloat(b) * (n - 1) / 2), BigFloat, RL_BITS)
                @test rl_err(f.DL, dlcf) <= 8 * epsT * (2n + GA + GB) / n + 4 * epsT * Float64(dlcf)
                @test f.n == n && f.handoffs == n ÷ tau
            end
        end
        @info "AC1 realized ratio max |error| / B_R^lin (R2)" T worst
        # (i') the decimals of the spec as rationals rounded once into T (exact(T, 1 // 10)), against the closed form of the real numbers
        for (ra, rb) in ((1 // 100, 0 // 1), (0 // 1, 1 // 50), (1 // 100, 1 // 50)), tau in (1, 3, 5)
            cT, xs, ps, a, b = exact(T, 3 // 10), exact(T, 1 // 10), exact(T, 1 // 20), exact(T, ra), exact(T, rb)
            ns = filter(n -> n % tau == 0, collect(nlist))
            z = SVector{4, T}(xs, ps, xs + a, ps + b)
            snaps = rl_run((RLI(),), PairMap(rl_T1(cT)), [z], [W], maximum(ns), tau; steps = ns)[1][1].snaps
            for (n, f) in zip(ns, snaps)
                cf = setprecision(() -> (BigFloat(3) / 10) * (BigFloat(numerator(ra)) / denominator(ra) + BigFloat(numerator(rb)) / denominator(rb) * (n + tau - 2) / 4), BigFloat, RL_BITS)
                B = rl_T1_B(epsT, rl_T1_G(cT, xs, ps, tau), rl_T1_G(cT, xs + a, ps + b, tau), n, tau, cf)
                @test rl_err(f.RLI, cf) <= B
            end
        end
        # (ii)
        for (a, b) in seps_ii, tau in (1, 3, 5)
            fin = run_pair(a, b, tau)
            N = maximum(keys(fin))
            ref = ref_pair(a, b, tau, N)[1]
            for (n, f) in fin
                cf = ref.RLI[n ÷ tau]
                GA, GB = rl_T1_G(c, x0, p0, tau), rl_T1_G(c, x0 + a, p0 + b, tau)
                @test rl_err(f.RLI, cf) <= rl_T1_B(epsT, GA, GB, n, tau, cf)
                @test rl_err(f.DL, ref.DL[n ÷ tau]) <= 8 * epsT * (2n + GA + GB) / n + 4 * epsT * Float64(ref.DL[n ÷ tau])
            end
        end
        # (iii) linearity in the separation (a, b) -> (2a, 2b), tau = 5
        for (a, b) in seps_i
            f1, f2 = run_pair(a, b, 5), run_pair(2a, 2b, 5)
            for n in (5, 60)
                cf = rl_T1_closed(c, a, b, n, 5)
                B = rl_T1_B(epsT, rl_T1_G(c, x0, p0, 5), rl_T1_G(c, x0 + 2a, p0 + 2b, 5), n, 5, 2cf)
                @test setprecision(() -> Float64(abs(rl_big(f2[n].RLI) - 2 * rl_big(f1[n].RLI))), BigFloat, RL_BITS) <= 2 * B
            end
        end
        # (iv) negative controls: each faulty double fails the predicate of (i) or (ii) (computed from the exact sums, the same predicate as above)
        for (a, b) in ((5 / 512, 5 / 256),), tau in (3, 5)
            N = 60
            sA, sB = rl_T1_sums(c, x0, p0, N), rl_T1_sums(c, x0 + a, p0 + b, N)
            cf(n) = rl_T1_closed(c, a, b, n, tau)
            B(n) = rl_T1_B(epsT, rl_T1_G(c, x0, p0, tau), rl_T1_G(c, x0 + a, p0 + b, tau), n, tau, cf(n))
            for kw in (:divn, :mean0, :literal)
                bad = rl_ref_rli(sA, sB, tau; NamedTuple{(kw,)}((true,))...)
                @test any(n -> Float64(abs(bad.RLI[n ÷ tau] - cf(n))) > B(n), (tau * 4, tau * 10))      # fails (i) for tau > 1
            end
            # independent vectors: the shadow started on (0, 1) (p-direction, whose stretch is exp(-c x), not exp(c x))
            Wind = SMatrix{4, 2, T}(1, 0, 0, 0, 0, 0, 0, 1)
            fin = run_pair(a, b, tau, Wind)
            @test any(n -> rl_err(fin[n].RLI, cf(n)) > B(n), intersect(keys(fin), (20, 60)))
        end
        for (a, b) in seps_ii[1:1], tau in (1,)
            N = 60
            sA, sB = rl_T1_sums(c, x0, p0, N), rl_T1_sums(c, x0 + a, p0 + b, N)
            good = rl_ref_rli(sA, sB, tau)
            bad = rl_ref_rli(sA, sB, tau; noabs = true)
            B(n) = rl_T1_B(epsT, rl_T1_G(c, x0, p0, tau), rl_T1_G(c, x0 + a, p0 + b, tau), n, tau, good.RLI[n])
            @test any(n -> Float64(abs(bad.RLI[n] - good.RLI[n])) > B(n), (20, 60))                 # fails (ii)
        end
    end
end

@testset "AC2 constant tangent gives exactly zero" begin
    _ = (RLI, PairMap)
    # RL5(d), I2. Y5 with A = [[2, 1], [1, 1]]: the pair at the origin and at Dx xi_0, Dx in {1e-12, 1e-6}, ids 1..100 (xi_0 seeded by the id), tau in {1, 5}, n in {1, 10, 200}
    # (the shadow reaches 1e-6 mu^200 below the overflow range); the K = 0 standard map (wrapped and lifted), panel starts, Dx in {1e-12, 1e-3}, n in {10, 1e3, 1e5}
    # (D38: Double64 on 4 starts to n = 1e3): every DL_k and every RLI is 0 with tolerance 0, and s_A = s_B bitwise (sum and compensation). Negative control: a double with
    # independent vectors (the shadow started on (0, 1) of the K = 0 shear) has DL_n = ln(sqrt(n^2 + 1)) / n.
    seed = rl_seed(2)
    for T in RL_TYPES
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        P = PairMap(sys)
        for Dx in (1e-12, 1e-6), tau in (1, 5)
            ns = filter(n -> n % tau == 0, (1, 10, 200))
            zs = Vector{SVector{4, T}}()
            Ws = Vector{SMatrix{4, 2, T, 8}}()
            for j in 1:100
                xi = rl_xi0(seed, j, T)
                push!(zs, rl_pair(sys, SVector{2, T}(0, 0), xi, Dx).z)
                push!(Ws, rl_frame(xi))
            end
            res = rl_run((RLI(),), P, zs, Ws, 200, tau; steps = ns)
            @test rl_nbad(r -> length(r[1].snaps) == length(ns) && all(f -> f.DL == 0 && f.RLI == 0 && f.sum[1] == f.sum[2] && f.comp[1] == f.comp[2], r[1].snaps), res) == 0
        end
    end
    for wrapped in (true, false), T in RL_TYPES, Dx in (1e-12, 1e-3)
        m, N = T === Float64 ? (RL_FULL ? 20 : 10, 100_000) : (RL_FULL ? 20 : 4, RL_FULL ? 100_000 : 1000)    # D38
        c = rl_case(:zero, T, m; wrapped = wrapped)
        prs = rl_pairs(c, Dx)
        ns = filter(<=(N), (10, 1000, 100_000))
        res = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), N, 1; steps = ns)
        @test rl_nbad(r -> length(r[1].snaps) == length(ns) && all(f -> f.DL == 0 && f.RLI == 0 && f.sum[1] == f.sum[2] && f.comp[1] == f.comp[2], r[1].snaps), res) == 0
    end
    # negative control: independent vectors on the K = 0 shear (A: (1, 0) stays; B: (0, 1) becomes (n, 1))
    sys = StandardMap(0.0)
    z = SVector(1.0, 2.0, 1.0 + 1e-3, 2.0)
    res = rl_run((RLI(),), PairMap(sys), [z], [SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 0, 0, 1)], 100, 1; steps = [10, 100])
    for (n, f) in zip((10, 100), res[1][1].snaps)
        @test f.DL > 0 && isapprox(f.DL, log(sqrt(n^2 + 1.0)) / n; rtol = 1e-12)
    end
end

@testset "AC3 same-orbit reference" begin
    _ = (RLI, PairMap)
    # Standard map, K in {0.5 (regular class), 2 pi, 10}, Dx in {1e-6, 1e-12}, tau in {1, 4}, wrapped state, both types, at the checkpoints (10 per decade; every multiple of
    # 10 with RL_FULL): |RLI - RLI_ref| <= B_R(n) = 8 eps (Psi_max^A + Psi_max^B) H_K / n + 4 eps RLI_ref (RL11), Psi_max the maximum of Psi_g over the checkpoints g <= n
    # for each orbit, the reference the BigFloat(256) exact sums at the realized points of both orbits. Also DL_k against the first display of RL11. Sizes D38: 6 starts per class
    # to n = 1e3 and 2 starts to n = 1e4 (the spec: 20 starts to 1e4; RL_FULL). Reported: the realized ratio error / B_R (R2). Negative control: a double with a Float32 accumulator
    # exceeds B_R at n = 1e3 for the chaotic starts.
    for T in RL_TYPES, (label, N, m) in rl_panel(), Dx in (1e-6, 1e-12), tau in (1, 4)
        rows = rl_rows(label, T, Dx, tau, m, N)
        @test length(rows) == m * count(n -> n % tau == 0, rl_checkpoints(N))
        @test rl_nbad(r -> r.eR <= r.BR, rows) == 0
        @test rl_nbad(r -> r.eD <= r.BD, rows) == 0
        worst, worstD = maximum(r.rR for r in rows), maximum(r.rD for r in rows)
        @info "AC3 realized ratio max |error| / B_R (R2)" label T Dx tau N worst worstD
    end
    for (label, N, m) in rl_panel()
        (label === :half || N != 1000) && continue
        rows = rl_rows(label, Float64, 1e-12, 1, m, N)
        refs = rl_pair_refs(label, Float64, 1e-12, m, N)
        for j in 1:min(m, 3)
            r = only(r for r in rows if r.j == j && r.n == 1000)
            sA, sB = refs[j].sA, refs[j].sB
            r32 = setprecision(RL_BITS) do
                a = b = q = Float32(0)
                pa = pb = zero(BigFloat)
                for k in 1:1000
                    a += Float32(sA[k] - pa); b += Float32(sB[k] - pb); pa, pb = sA[k], sB[k]
                    q += abs(b - a) / Float32(k)
                end
                Float64(q / Float32(1000))
            end
            @test abs(r32 - Float64(r.refRLI)) > r.BR                  # the Float32 accumulator is far outside the bound
        end
    end
end
