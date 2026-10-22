# AC1 to AC5: algebra of the definition, K = 0 shear, linear hyperbolic map, hyperbolic fixed point, rotation and elliptic cases.
# Each testset starts by naming the Indicators symbols it uses, so a missing module errors that testset (UndefVarError) and nothing else.
# Tolerances: "6 B_S(n)" with B_S(n) = 8 eps (n + G_first) of L1 S3 (the coefficient sum of the three series is below 6, M11), except for the supplied frames
# of AC3 with a tiny expanding coefficient (|c_1| = 1e-3, 1e-6), which use max(6 B_S(n), 48 eps Psi_n^(i)), the per-vector form of L1 Amendment D35
# (the plain B_S is exceeded there, L1 P6; fli-ofli Amendment D72 item 2).

# ------------------------------------------------------------------ AC1

# The criterion as a predicate on the three series (vectors indexed by handoff, Float64/T or BigFloat) against a reference with fields Y, Ybar, Ybar20 and
# the tolerances B[k], so that the negative controls can be run through it.
mg_alg_ok(Y, Yb, Y20, ref, B) = all(k -> mg_err(Y[k], ref.Y[k]) <= B[k] && mg_err(Yb[k], ref.Ybar[k]) <= B[k] && mg_err(Y20[k], ref.Ybar20[k]) <= B[k], eachindex(Y))

# Reference of the constant increments lam (M8(a)) at the handoffs 1..K: closed forms, as vectors.
function mg_const_ref(lam::BigFloat, tau::Int, K::Int)
    f = [mg_const_forms(lam, tau, k) for k in 1:K]
    return (Y = [x.Y for x in f], Ybar = [x.Ybar for x in f], Ybar20 = [x.Ybar20 for x in f])
end

@testset "AC1 algebra of the definition" begin
    _ = (MEGNO,)
    # Test maps with the prescribed diagonal tangent diag(e^(a_m), e^(-a_m)) and a supplied identity frame: column 1 carries +a_m, column 2 carries -a_m.
    # Reference helpers first: the weighted sum M5 and the Abel form M6 agree (the two algebraic forms of the spec), the harmonic form of Ybar too, and the
    # printed values of M7 and M8 come out (this checks the test code, not the module).
    setprecision(BigFloat, MG_BITS) do
        for tau in (1, 5)
            sh = [log(1 + BigFloat(k) * tau) for k in 1:200]
            a, b = mg_M5(sh, tau), mg_M6(sh, tau)
            @test all(k -> abs(a.Y[k] - b.Y[k]) < 1e-60 * max(1, abs(a.Y[k])) && abs(a.Ybar[k] - b.Ybar[k]) < 1e-60 * max(1, abs(a.Ybar[k])) &&
                           abs(a.Ybar20[k] - b.Ybar20[k]) < 1e-60 * max(1, abs(a.Ybar20[k])), 1:200)
            @test all(K -> abs(mg_Ybar_harmonic(sh, K) - a.Ybar[K]) < 1e-60 * max(1, abs(a.Ybar[K])), (1, 2, 7, 40, 200))
        end
        lam = BigFloat(0.05)
        for tau in (1, 5)
            sh = [lam * tau * k for k in 1:40]
            a = mg_M5(sh, tau)
            cf = mg_const_ref(lam, tau, 40)
            @test all(k -> abs(a.Y[k] - cf.Y[k]) < 1e-60 && abs(a.Ybar[k] - cf.Ybar[k]) < 1e-60 && abs(a.Ybar20[k] - cf.Ybar20[k]) < 1e-60, 1:40)
        end
        shA = [BigFloat(1) / 20 * k for k in 1:40]                                       # M7: s_A = 0.05 k
        shB = [k <= 20 ? BigFloat(1) / 10 * k : BigFloat(2) for k in 1:40]               # s_B = 0.1 k up to k = 20, then 2
        @test abs(mg_M5(shA, 1).Y[40] - BigFloat(41) / 20) < 1e-60 && abs(mg_M5(shB, 1).Y[40] - BigFloat(21) / 20) < 1e-60
        @test mg_M5(shA, 1).Y[1] == 2 * shA[1] && mg_M5(shB, 1).Y[1] == 2 * shB[1]       # first handoff: Y(tau) = 2 a_1
    end
    for T in MG_TYPES
        epsT = mg_eps(T)
        x0 = SVector{2, T}(0, 0)
        I2 = SMatrix{2, 2, T}(1, 0, 0, 1)
        snaps(kind, N, tau; W = I2) = mg_run((MEGNO(),), mg_Diag(kind), [x0], [W], N, tau; all = true)[1][1].snaps
        B6(kind, n, tau) = 6 * 8 * epsT * (n + mg_G_diag(kind, tau))
        # (i) a_m = 0.05 constant, j <= 40, tau in {1, 5}: Y, Ybar, Ybar20 equal the exact forms of M8(a) within 6 B_S(n), both columns (column 2: lambda = -0.05)
        lam0 = BigFloat(mg_c(:const, 1))
        for tau in (1, 5)
            sn = snaps(:const, 40, tau)
            K = 40 ÷ tau
            @test length(sn) == K && [s.n for s in sn] == collect(tau:tau:40)
            B = [B6(:const, k * tau, tau) for k in 1:K]
            sh = [lam0 * tau * k for k in 1:K]
            for (i, sg) in ((1, 1), (2, -1))
                ref = mg_const_ref(sg * lam0, tau, K)
                Y, Yb, Y20 = mg_col(sn, :Y, i), mg_col(sn, :Ybar, i), mg_col(sn, :Ybar20, i)
                @test mg_alg_ok(Y, Yb, Y20, ref, B)
                # negative controls: each double fails the predicate
                Ls = [sg * lam0 for _ in 1:K]                                       # s(j)/j, the L1 value: no memory
                @test !mg_alg_ok(Ls, Ls, Ls, ref, B)
                for kw in ((noweight = true,), (nofactor = true,), (mean0 = true,))
                    bad = mg_M5(sg .* sh, tau; kw...)
                    @test !mg_alg_ok(bad.Y, bad.Ybar, bad.Ybar20, ref, B)
                end
                @test mg_err(mg_M5(sg .* sh, tau; mean0 = true).Ybar[K], ref.Ybar[K]) > 0.9 * lam0 * tau     # the mean over 0..K-1 differs by about lambda tau
            end
            # exchange of the columns of the frame: the two columns exchange
            sw = snaps(:const, 40, tau; W = SMatrix{2, 2, T}(0, 1, 1, 0))
            @test isequal(mg_col(sw, :Y, 1), mg_col(sn, :Y, 2)) && isequal(mg_col(sw, :Ybar20, 2), mg_col(sn, :Ybar20, 1))
        end
        # (ii) a_m = ln((m + 1)/m), so s(k) = ln(1 + k), j <= 200: Y(j) = 2 (ln(j + 1) - ln(j!)/j) within 6 B_S(n); Ybar, Ybar20 against the M5 sums
        sn = snaps(:log1, 200, 1)
        B = [B6(:log1, j, 1) for j in 1:200]
        lnfact = setprecision(() -> accumulate(+, [log(BigFloat(k)) for k in 1:200]), BigFloat, MG_BITS)
        Ycf = setprecision(() -> [2 * (log(BigFloat(j + 1)) - lnfact[j] / j) for j in 1:200], BigFloat, MG_BITS)
        sh1 = setprecision(() -> [log(BigFloat(k + 1)) for k in 1:200], BigFloat, MG_BITS)
        ref1, ref2 = mg_M5(sh1, 1), mg_M5(-sh1, 1)
        @test all(j -> mg_err(sn[j].Y[1], Ycf[j]) <= B[j] && mg_err(sn[j].Y[2], -Ycf[j]) <= B[j], 1:200)
        @test mg_alg_ok(mg_col(sn, :Y, 1), mg_col(sn, :Ybar, 1), mg_col(sn, :Ybar20, 1), ref1, B)
        @test mg_alg_ok(mg_col(sn, :Y, 2), mg_col(sn, :Ybar, 2), mg_col(sn, :Ybar20, 2), ref2, B)
        @test all(j -> abs(ref1.Y[j] - Ycf[j]) < 1e-60, 1:200)                       # the M5 reference against the closed form of the spec
        # (iii) a_m = +0.2, -0.2 alternating: against the BigFloat sums, Ybar and Ybar20 tending to 0
        sn = snaps(:alt, 200, 1)
        c = setprecision(() -> cumsum([BigFloat(mg_c(:alt, m)) for m in 1:200]), BigFloat, MG_BITS)
        B = [B6(:alt, j, 1) for j in 1:200]
        @test mg_alg_ok(mg_col(sn, :Y, 1), mg_col(sn, :Ybar, 1), mg_col(sn, :Ybar20, 1), mg_M5(c, 1), B)
        @test mg_alg_ok(mg_col(sn, :Y, 2), mg_col(sn, :Ybar, 2), mg_col(sn, :Ybar20, 2), mg_M5(-c, 1), B)
        @test all(i -> abs(sn[200].Ybar[i]) <= 0.1 && abs(sn[200].Ybar20[i]) <= 0.1, 1:2)
        # (iv) memory: two paths with the same s(40) give Y(40) = 2.05 and 1.05 (M7), the L1 value s(40)/40 is the same for both
        sA, sB = snaps(:pathA, 40, 1), snaps(:pathB, 40, 1)
        B = [B6(:pathB, j, 1) for j in 1:40]
        cA = setprecision(() -> cumsum([BigFloat(mg_c(:pathA, m)) for m in 1:40]), BigFloat, MG_BITS)
        cB = setprecision(() -> cumsum([BigFloat(mg_c(:pathB, m)) for m in 1:40]), BigFloat, MG_BITS)
        refA, refB = mg_M5(cA, 1), mg_M5(cB, 1)
        @test mg_alg_ok(mg_col(sA, :Y, 1), mg_col(sA, :Ybar, 1), mg_col(sA, :Ybar20, 1), refA, B)
        @test mg_alg_ok(mg_col(sB, :Y, 1), mg_col(sB, :Ybar, 1), mg_col(sB, :Ybar20, 1), refB, B)
        @test abs(sA[40].Y[1] - 2.05) < 1e-12 && abs(sB[40].Y[1] - 1.05) < 1e-12           # the printed values of M7, loose (the stored c_m = 0.1 is not exactly 1/10)
        SA, SB = sA[40].sum[1] + sA[40].comp[1], sB[40].sum[1] + sB[40].comp[1]
        @test abs(SA - SB) <= 2 * 8 * epsT * (40 + mg_G_diag(:pathB, 1))                  # equal s(40)
        @test sA[40].Y[1] - sB[40].Y[1] > 0.99                                            # while Y differs by 1: a double returning s(j)/j cannot tell them apart
    end
    # (v) long run: constant a_m = 0.1, tau = 1, j in {1e5, 4e6}: the exact forms of M8(a) within 6 B_S(n), no integer overflow of the weights (K^3 at 4e6 is
    # above 2^63). D38: Float64 only (Double64 exp at 4e6 steps is slow; MG_FULL adds Double64 at 1e5). Negative control: Ybar20 formed with n^3 in 64-bit integers.
    for T in (MG_FULL ? MG_TYPES : (Float64,))
        epsT = mg_eps(T)
        ns = T === Float64 ? [100_000, 4_000_000] : [100_000]
        sn = mg_run((MEGNO(),), mg_Diag(:long), [SVector{2, T}(0, 0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], last(ns), 1; steps = ns)[1][1].snaps
        @test [s.n for s in sn] == ns
        lam = BigFloat(mg_c(:long, 1))
        for (s, n) in zip(sn, ns), (i, sg) in ((1, 1), (2, -1))
            f = mg_const_forms(sg * lam, 1, n)
            B = 6 * 8 * epsT * (n + mg_G_diag(:long, 1))
            @test mg_err(s.Y[i], f.Y) <= B && mg_err(s.Ybar[i], f.Ybar) <= B && mg_err(s.Ybar20[i], f.Ybar20) <= B
            @test s.handoffs == n
        end
        if T === Float64
            n = 4_000_000
            f = mg_const_forms(lam, 1, n)
            Q2 = setprecision(() -> Float64(f.Ybar20 * BigFloat(n)^3), BigFloat, MG_BITS)       # tau Q2 = Ybar20 n^3
            bad = Q2 / (n^3)                                                                    # n^3 in Int64: wraps at 2^63
            @test n^3 != BigInt(n)^3 && mg_err(bad, f.Ybar20) > 6 * 8 * epsT * (n + mg_G_diag(:long, 1))
        end
    end
end

# ------------------------------------------------------------------ AC2

# s_i(n) = (1/2) ln(q_i(n) / q_i(0)) at the steps `steps` for the shear: column (a, b) -> (a + n b, b) (L1 E11), BigFloat from the exact widened frame.
function mg_shear_s(W0, steps, i::Int)
    return setprecision(BigFloat, MG_BITS) do
        a, b = mg_big(W0[1, i]), mg_big(W0[2, i])
        q0 = a^2 + b^2
        [log(((a + n * b)^2 + b^2) / q0) / 2 for n in steps]
    end
end

@testset "AC2 K = 0 shear" begin
    _ = (MEGNO,)
    # The tangent is [[1, 1], [0, 1]] everywhere (tangent-core E8): the answer does not depend on the orbit. Seeded frames (20 starts) and W0 = I, tau in {1, 10},
    # n in {10, 100, 1e4, 1e6}: |M - M_cf| <= 6 B_S(n), M_cf the BigFloat evaluation of M5 on the exact norms (L1 E11), both columns, all three series. W0 = I:
    # column 1 exactly 0 (tolerance 0, s = ln 1); column 2 at tau = 1: Y within 6 B_S(n) of ln(n^2 + 1) - (1/n) ln prod (i^2 + 1) and, for n >= 100,
    # |Y(n) - 2 - (ln n - ln(2 sinh pi))/n| <= 3/n^2 + 6 B_S(n). D38: n = 1e6 at tau = 10 and Float64 only (3 seeded starts and the identity; a BigFloat reference
    # at tau = 1 costs 2e6 logarithms per frame); MG_FULL adds tau = 1 at 1e6.
    seed = mg_seed(2)
    rng = StableRNG(seed)
    M = SMatrix{2, 2, Float64}(1, 0, 1, 1)
    worst = 0.0                                                         # realized ratio error / (6 B_S) (R2)
    # the helpers against the closed form of the spec: M5 on the exact shear sums, W0 = I column 2, tau = 1
    setprecision(BigFloat, MG_BITS) do
        sh = [log(BigFloat(n)^2 + 1) / 2 for n in 1:300]
        lp = accumulate(+, [log(BigFloat(i)^2 + 1) for i in 0:299])          # lp[n] = sum_{i=0..n-1} ln(i^2 + 1)
        ref = mg_M5(sh, 1)
        @test all(n -> abs(ref.Y[n] - (log(BigFloat(n)^2 + 1) - lp[n] / n)) < 1e-60, (10, 100, 300))
        lns = log(2 * sinh(BigFloat(pi)))
        @test all(n -> abs(ref.Y[n] - 2 - (log(BigFloat(n)) - lns) / n) <= 2 / BigFloat(n)^2, (100, 300))   # the drafting-check remainder 1.8/n^2
    end
    for T in MG_TYPES
        epsT = mg_eps(T)
        sys = StandardMap(exact(T, 0))
        combos = [(1, 10_000), (10, 10_000)]
        T === Float64 && push!(combos, (10, 1_000_000))
        T === Float64 && MG_FULL && push!(combos, (1, 1_000_000))
        for (tau, N) in combos, (case, nid) in (("seeded", N == 1_000_000 ? 3 : 20), ("identity", 1))
            m = case == "identity" ? 1 : nid
            x0 = [mg_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
            W0 = case == "identity" ? [SMatrix{2, 2, T}(1, 0, 0, 1)] : [mg_frame(seed, j, T) for j in 1:m]
            ns = N == 1_000_000 ? [1_000_000] : [10, 100, 10_000]
            runs = mg_run((MEGNO(),), sys, x0, W0, N, tau; steps = ns)
            B(n) = 6 * 8 * epsT * (n + mg_Gfirst(M, tau))
            for j in 1:m
                snaps = runs[j][1].snaps
                @test [s.n for s in snaps] == ns
                steps = collect(tau:tau:N)
                s_ex = [mg_shear_s(W0[j], steps, 1), mg_shear_s(W0[j], steps, 2)]
                ref = [mg_M5(s_ex[1], tau), mg_M5(s_ex[2], tau)]
                for s in snaps
                    k = s.n ÷ tau
                    for i in 1:2
                        es = (mg_err(s.Y[i], ref[i].Y[k]), mg_err(s.Ybar[i], ref[i].Ybar[k]), mg_err(s.Ybar20[i], ref[i].Ybar20[k]))
                        @test all(<=(B(s.n)), es)
                        worst = max(worst, maximum(es) / B(s.n))
                    end
                    if case == "identity"
                        @test s.Y[1] == 0 && s.Ybar[1] == 0 && s.Ybar20[1] == 0                      # tolerance 0: s = ln 1
                        if tau == 1 && N <= 10_000
                            lp = setprecision(() -> sum(log(BigFloat(i)^2 + 1) for i in 0:(s.n - 1)), BigFloat, MG_BITS)
                            Ycf = setprecision(() -> log(BigFloat(s.n)^2 + 1) - lp / s.n, BigFloat, MG_BITS)
                            @test mg_err(s.Y[2], Ycf) <= B(s.n)
                            if s.n >= 100
                                asym = setprecision(() -> 2 + (log(BigFloat(s.n)) - log(2 * sinh(BigFloat(pi)))) / s.n, BigFloat, MG_BITS)
                                @test mg_err(s.Y[2], asym) <= 3 / s.n^2 + B(s.n)
                            end
                        end
                    end
                end
            end
        end
    end
    @info "AC2 realized ratio max |error| / (6 B_S) (R2)" worst
end

# ------------------------------------------------------------------ AC3 and AC4: constant tangents (Y5 and the fixed point of the standard map)

@testset "AC3 linear hyperbolic map" begin
    _ = (MEGNO,)
    # A = [[2, 1], [1, 1]] (Y5), start at the origin (a fixed point), seeded ids 1..100 (no exclusion, Amendment D74 item 2) and two supplied unit frames with |c_1| = 1e-3 and 1e-6 for column 1, tau in {1, 5}. (a) n in {1, 2, 4, 8, 20, 100, 200} (tau = 5: the handoffs
    # {5, 10, 20, 100, 200}): all three series, both columns, within 6 B_S(n) of the M5 sums of the BigFloat exact norms (the two supplied columns 1: max(6 B_S, 48 eps
    # Psi)); (b) the law of M10(b) at n in {20, 100, 200}: Y(n) = (n + tau) ln mu - (2/K) sum_{k<K} delta_k + 2 delta_K within the same bound, and
    # |Y(n) - (n + tau) ln mu| <= 2 |R|/K + bound + 2.01 (c_2/c_1)^2 mu^(-4n); (c) |c_1| = 1e-3 at n = 4, tau = 1: the column that starts nearly along the
    # contracting direction (column 1) has a negative Y, column 2 (the expanding one) Y > 0.5 (test plan, Deviations item 3).
    seed = mg_seed(3)
    worst = 0.0                                                         # realized ratio error / bound (R2)
    for T in MG_TYPES
        epsT = mg_eps(T)
        nid = 100
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        W0 = vcat([mg_frame(seed, j, T) for j in 1:nid], [mg_y5_frame(1e-3, T), mg_y5_frame(1e-6, T)])
        nfr = length(W0)
        x0 = [SVector{2, T}(0, 0) for _ in 1:nfr]
        M64 = SMatrix{2, 2, Float64}(2, 1, 1, 1)
        Ms = lb_tangents(sys, x0[1], 200)                       # the tangent is A everywhere
        sv = [mg_ref_s(mg_mats_const(2, 1, 1, 1), W0[j], 200) for j in 1:nfr]
        memo = [[Dict{Int, Float64}(), Dict{Int, Float64}()] for _ in 1:nfr]
        cs = [mg_y5_coeffs(W0[j]) for j in 1:nfr]
        # supplied frames as intended: |c_1| of column 1, and the second column carries the expanding part
        @test abs(abs(cs[nid + 1].c1[1]) - 1e-3) < 1e-12 && abs(abs(cs[nid + 2].c1[1]) - 1e-6) < 1e-12 && abs(cs[nid + 1].c1[2]) > 0.999
        # the iterated reference against the matrix power A^8 (independent path): the squared norms of Y_8 w_i(0) / norm(w_i(0))
        setprecision(BigFloat, MG_BITS) do
            Y8 = [BigFloat(2) 1; 1 1]^8
            for j in 1:3, i in 1:2
                w = (mg_big(W0[j][1, i]), mg_big(W0[j][2, i]))
                y = (Y8[1, 1] * w[1] + Y8[1, 2] * w[2], Y8[2, 1] * w[1] + Y8[2, 2] * w[2])
                s8 = log((y[1]^2 + y[2]^2) / (w[1]^2 + w[2]^2)) / 2
                @test abs(s8 - (i == 1 ? sv[j].s1[8] : sv[j].s2[8])) < 1e-60
            end
        end
        for tau in (1, 5)
            ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
            run = mg_run((MEGNO(),), sys, x0, W0, 200, tau; steps = ns)
            for j in 1:nfr, i in 1:2
                si = i == 1 ? sv[j].s1 : sv[j].s2
                ref = mg_M5(mg_handoffs(si, tau), tau)
                snaps = run[j][1].snaps
                @test [s.n for s in snaps] == ns
                c1, c2 = cs[j].c1[i], cs[j].c2[i]
                for s in snaps
                    n, k = s.n, s.n ÷ tau
                    B = mg_lin_B(epsT, M64, Ms, W0[j][:, i], n, tau, memo[j][i], j > nid && i == 1)
                    es = (mg_err(s.Y[i], ref.Y[k]), mg_err(s.Ybar[i], ref.Ybar[k]), mg_err(s.Ybar20[i], ref.Ybar20[k]))
                    @test all(<=(B), es)                                                              # (a)
                    worst = max(worst, maximum(es) / B)
                    if n in (20, 100, 200)                                                           # (b)
                        lam, dk, Ylaw, R, tail = setprecision(BigFloat, MG_BITS) do
                            mu = mg_y5_eigen().mu
                            x = (c2 / c1)^2
                            d(kk) = log(1 + x * mu^(-4 * kk * tau)) / 2
                            dsum = sum(d(kk) for kk in 0:(k - 1))
                            lam = log(mu)
                            Rinf = zero(BigFloat)                                                     # R = -2 sum_{k >= 0} delta_k, to a negligible tail
                            kk = 0
                            while true
                                t = d(kk)
                                Rinf -= 2 * t
                                (t < 1e-80 && kk > 3) && break
                                kk += 1
                            end
                            (lam, d(k), (n + tau) * lam - 2 * dsum / k + 2 * d(k), Rinf, Float64(2.01 * x * mu^(-4n)))
                        end
                        @test mg_err(s.Y[i], Ylaw) <= B
                        @test mg_err(s.Y[i], (n + tau) * lam) <= 2 * Float64(abs(R)) / k + B + tail
                    end
                    if j == nid + 1 && tau == 1 && n == 4                                             # (c)
                        i == 1 && @test ref.Y[4] < 0 && s.Y[1] < 0                                    # column 1 starts nearly along the contracting direction
                        i == 2 && @test ref.Y[4] > 0.5 && s.Y[2] > 0.5                               # column 2 carries the expanding part
                        i == 1 && @info "AC3(c) reported: Y, Ybar, Ybar20 of the two columns at n = 4 for |c_1| = 1e-3" T Y = Float64.(s.Y) Ybar = Float64.(s.Ybar) Ybar20 = Float64.(s.Ybar20)
                    end
                end
            end
        end
    end
    @info "AC3 realized ratio max |error| / bound (R2)" worst
end

@testset "AC4 hyperbolic fixed point" begin
    _ = (MEGNO,)
    # Fixed point (0, 0) of the standard map, K in {1, 2, 6}, seeded ids 1..100, tau in {1, 5}, n as in AC3 (L1 E12): all three series, both columns, within 6 B_S(n)
    # of the M5 sums of the BigFloat product norms; slope check at n = 200, tau = 1: |Y(200) - 201 ln mu_+| <= 2 |R|/K + 6 B_S + 0.1, R = -2 sum_k delta_k
    # from the exact sums (delta_k = s(k) - k lambda - b, b the limit of s(k) - k lambda).
    seed = mg_seed(4)
    worst = 0.0                                                         # realized ratio error / (6 B_S) (R2)
    for T in MG_TYPES, K in (1, 2, 6)
        epsT = mg_eps(T)
        nid = 100
        sys = StandardMap(exact(T, K))
        W0 = [mg_frame(seed, j, T) for j in 1:nid]
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        M64 = SMatrix{2, 2, Float64}(1 + K, K, 1, 1)
        Ms = lb_tangents(sys, x0[1], 200)
        sv = [mg_ref_s(mg_mats_const(1 + K, 1, K, 1), W0[j], 200) for j in 1:nid]
        memo = [[Dict{Int, Float64}(), Dict{Int, Float64}()] for _ in 1:nid]
        lam = setprecision(() -> log(((2 + K) + sqrt(BigFloat(K + 2)^2 - 4)) / 2), BigFloat, MG_BITS)
        for tau in (1, 5)
            ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
            run = mg_run((MEGNO(),), sys, x0, W0, 200, tau; steps = ns)
            for j in 1:nid, i in 1:2
                si = i == 1 ? sv[j].s1 : sv[j].s2
                ref = mg_M5(mg_handoffs(si, tau), tau)
                for s in run[j][1].snaps
                    k = s.n ÷ tau
                    B = mg_lin_B(epsT, M64, Ms, W0[j][:, i], s.n, tau, memo[j][i])
                    es = (mg_err(s.Y[i], ref.Y[k]), mg_err(s.Ybar[i], ref.Ybar[k]), mg_err(s.Ybar20[i], ref.Ybar20[k]))
                    @test all(<=(B), es)
                    worst = max(worst, maximum(es) / B)
                    if tau == 1 && s.n == 200
                        R = setprecision(BigFloat, MG_BITS) do
                            b = si[200] - 200 * lam
                            -2 * sum(si[kk] - kk * lam - b for kk in 1:200) - 2 * (0 - 0 - b)     # k = 0 term: delta_0 = s(0) - b = -b
                        end
                        @test mg_err(s.Y[i], 201 * lam) <= 2 * Float64(abs(R)) / 200 + B + 0.1
                    end
                end
            end
        end
    end
    @info "AC4 realized ratio max |error| / (6 B_S) (R2)" worst
end

# ------------------------------------------------------------------ AC5

@testset "AC5 rotation and elliptic cases" begin
    _ = (MEGNO,)
    # (a) the order-4 rotation [[0, 1], [-1, 0]], supplied identity frame (as L1 criterion 2(a)), n up to 1e4, tau in {1, 3}: all three series exactly 0, both
    # columns, at every handoff. (b) K = -2 at (0, 0) (M10(e), L1 E13), seeded ids 1..100, tau = 1: Y and Ybar against the exact forms (Y = -c for n even,
    # c (1 + 1/n) for n odd, Ybar = (c/n) ((n mod 2) + sum over odd m <= n of 1/m), c = ln norm(M w(0)) - ln norm(w(0))), |Y| <= ln((3 + sqrt 5)/2) (1 + 1/n) + bound,
    # Ybar20 against the M5 sums; n in 1..100 and every 10th to 1000 (D38; MG_FULL: every n). (c) the same at tau = 2: all three series 0 within the bound at every handoff.
    seed = mg_seed(5)
    worst = 0.0                                                         # realized ratio error / bound, tau = 1 (R2)
    for T in MG_TYPES
        epsT = mg_eps(T)
        rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
        for tau in (1, 3)
            N = tau == 1 ? 10_000 : 9_999
            snaps = mg_run((MEGNO(),), rot, [SVector{2, T}(0, 0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], N, tau; all = true)[1][1].snaps
            @test length(snaps) == N ÷ tau
            @test mg_nbad(s -> all(iszero, s.Y) && all(iszero, s.Ybar) && all(iszero, s.Ybar20), snaps) == 0
        end
        sys = StandardMap(exact(T, -2))
        nid = 100
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        W0 = [mg_frame(seed, j, T) for j in 1:nid]
        nsamp = MG_FULL ? collect(1:1000) : sort(unique(vcat(1:100, 110:10:1000)))
        M64 = SMatrix{2, 2, Float64}(-1, -2, 1, 1)
        Ms = lb_tangents(sys, x0[1], 1000)
        sv = [mg_ref_s(mg_mats_const(-1, 1, -2, 1), W0[j], 1000) for j in 1:nid]
        memo = [[Dict{Int, Float64}(), Dict{Int, Float64}()] for _ in 1:nid]
        ln0 = setprecision(() -> log((3 + sqrt(BigFloat(5))) / 2), BigFloat, MG_BITS)
        for tau in (1, 2)
            ns = filter(n -> n % tau == 0, nsamp)
            run = mg_run((MEGNO(),), sys, x0, W0, 1000, tau; steps = ns)
            for j in 1:nid, i in 1:2
                si = i == 1 ? sv[j].s1 : sv[j].s2
                ref = mg_M5(mg_handoffs(si, tau), tau)
                c = si[1]
                if tau == 1                                     # the helper against the closed forms of M10(e)
                    setprecision(BigFloat, MG_BITS) do
                        for n in (1, 2, 3, 10, 11, 100, 999, 1000)
                            Ycf = iseven(n) ? -c : c * (1 + inv(BigFloat(n)))
                            Ybcf = c / n * ((n % 2) + sum(inv(BigFloat(m)) for m in 1:2:n))
                            @test abs(ref.Y[n] - Ycf) < 1e-60 && abs(ref.Ybar[n] - Ybcf) < 1e-60
                        end
                    end
                end
                for s in run[j][1].snaps
                    n, k = s.n, s.n ÷ tau
                    B = mg_lin_B(epsT, M64, Ms, W0[j][:, i], n, tau, memo[j][i])
                    if tau == 1
                        Ycf = setprecision(() -> iseven(n) ? -c : c * (1 + inv(BigFloat(n))), BigFloat, MG_BITS)
                        Ybcf = setprecision(() -> c / n * ((n % 2) + sum(inv(BigFloat(m)) for m in 1:2:n)), BigFloat, MG_BITS)
                        es = (mg_err(s.Y[i], Ycf), mg_err(s.Ybar[i], Ybcf), mg_err(s.Ybar20[i], ref.Ybar20[k]))
                        @test all(<=(B), es)
                        worst = max(worst, maximum(es) / B)
                        @test abs(s.Y[i]) <= Float64(ln0) * (1 + 1 / n) + B
                    else
                        @test abs(s.Y[i]) <= B && abs(s.Ybar[i]) <= B && abs(s.Ybar20[i]) <= B       # (c)
                    end
                end
            end
        end
    end
    @info "AC5 realized ratio max |error| / bound, tau = 1 (R2)" worst
end
