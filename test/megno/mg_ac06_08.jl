# AC6 to AC8: same-orbit BigFloat reference with the bound B_M of M11, cross-check with L1, FLI and sali on the identical orbit (Mestre Eq. A2), growth laws.
# Shared computations (reference sums, Psi, runs, bounds) are memoized in MG_CACHE by mg_rows, so AC6 and AC9 reuse them.

# Y_1(n) at step n (tau = 1) of a double whose accumulator is a Float32: the increments of the exact sums, rounded to Float32, weighted by Float32(k) and added
# in Float32. Used only as the negative control of AC6 (it must exceed B_M).
function mg_single_Y(s::AbstractVector{BigFloat}, n::Int)
    return setprecision(BigFloat, MG_BITS) do
        R1 = Float32(0)
        prev = zero(BigFloat)
        for k in 1:n
            R1 += Float32(k) * Float32(s[k] - prev)
            prev = s[k]
        end
        Float64(2 * R1 / Float32(n))
    end
end

@testset "AC6 same-orbit reference" begin
    _ = (MEGNO,)
    # Standard map, K in {0.5 (regular class), 2 pi, 10}, 20 starts each, tau in {1, 4}, n up to 1e4, at the checkpoints (10 per decade; every multiple of 10
    # with MG_FULL), both columns, Float64 and Double64: |M - M_ref| <= B_M(n) = 48 eps Psi_max^(i)(n) for Y, Ybar and Ybar20 (M11), the reference the M5 sums of the
    # BigFloat(256) exact norms at the realized points, Psi_max the maximum of Psi_g over the checkpoints g <= n.
    # Reported: the realized ratio error / B_M (R2). Negative control: a double with a Float32 accumulator exceeds B_M at n = 1e3.
    for T in MG_TYPES, (label, N, m) in mg_panel(), tau in (1, 4)
        rows = mg_rows(label, T, tau, m, N)
        @test length(rows) == m * count(n -> n % tau == 0, mg_checkpoints(N))
        for f in (:Y, :Ybar, :Ybar20), i in 1:2
            @test mg_nbad(r -> r.e[f][i] <= r.B[i], rows) == 0
        end
        worst = maximum(r.ratio for r in rows)
        @test worst <= 1
        worstS = maximum(max(r.e.S[1] / (r.B[1] / 6), r.e.S[2] / (r.B[2] / 6)) for r in rows)
        @info "AC6 realized ratio max |error| / B_M (R2); the sums s_i against 8 eps Psi_max (reported)" label T tau worst worstS
    end
    for (label, N, m) in mg_panel()
        label === :half && continue
        rows = mg_rows(label, Float64, 1, m, N)
        c = mg_case(label, Float64, m)
        for j in 1:3
            r = only(r for r in rows if r.j == j && r.n == 1000)
            sv = mg_ref_s(mg_mats_std(c.K, mg_pts(c.sys, c.x0[j], 1000)), c.W0[j], 1000)
            @test abs(mg_single_Y(sv.s1, 1000) - Float64(r.r.Y[1])) > r.B[1]                  # the Float32 accumulator is far outside the bound
        end
    end
end

@testset "AC7 cross-check with L1, FLI and sali on the identical orbit" begin
    _ = (MEGNO, FLI, SALI, Benettin)
    # K in {0.5 (regular class), 2 pi, 10}, 20 starts, Float64, tau = 1, n up to 1e4, the chain {MEGNO (not the replacer), FLI (not the replacer), SALI} (the one-replacer
    # rule of D71, D72) at every multiple of 10; L1 (Benettin on k = 1 from the first and from the second seeded column) at every handoff. (a) the sums s_1, s_2 of the
    # MEGNO reducer are bitwise equal to those of the FLI and sali reducers (totals, and sum and compensation separately) and to the L1 sums per column; (b) column 1,
    # the three series against the BigFloat post-processing (M6) of the stored L1 series k X_L1(k): within 8 eps max(1, max |s_1|); (c) the Mestre relation (M13, Eq. A2)
    # with FLI(k) = k X_L1(k) and FLIbar of Eq. (3): |Y_1(n) - [2 (FLI(n) - FLIbar(n)) + (2/n) FLI(n)]| <= 8 eps max(1, max |s_1|), also with FLI(n) replaced by the final
    # sum S_1 of the FLI reducer, and the remainder equal to 2 X_L1(n); (d) link to the delivered running maximum: D = FLI_1 - k X_L1 >= -tol, Y_1 + 2 D >= -tol and
    # Y_1 = 2 (FLI_1 - D - FLIbar) + (2/n)(FLI_1 - D) within tol, tol = 8 eps max(1, max |s_1|, FLI_1(n)). Negative control: A2 with the running maximum FLI_1(n) in place of
    # FLI(n) fails (c) by 2 D (1 + 1/n) wherever D > 0. Reported (R1): the distribution of Y_1 + 2 D and the labels at n = 1e3.
    T = Float64
    epsT = mg_eps(T)
    N = 10_000
    cps = collect(10:10:N)
    ncontrol, nfail_control = 0, 0
    for label in (:half, :two_pi, :ten)
        m = 20
        c = mg_case(label, T, m)
        run = mg_run((MEGNO(replace = false), FLI(replace = false), SALI()), c.sys, c.x0, c.W0, N, 1; steps = cps)
        l1 = [mg_run((Benettin(),), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, i]) for W in c.W0], N, 1; all = true) for i in 1:2]
        yd = Dict(n => Float64[] for n in (100, 1_000, 10_000))
        lab = NamedTuple[]
        worst_b = worst_c = 0.0
        for j in 1:m
            msn, fsn, ssn = run[j][1].snaps, run[j][2].snaps, run[j][3].snaps
            X1, X2 = l1[1][j][1].snaps, l1[2][j][1].snaps
            @test length(msn) == length(cps) && length(fsn) == length(cps) && length(ssn) == length(cps) && length(X1) == N && length(X2) == N
            sL = [k * X1[k].lambda[1] for k in 1:N]                                  # the stored L1 series k X_L1(k) = s_1(k)
            smax = accumulate(max, abs.(sL))
            post, fbar_pre = setprecision(BigFloat, MG_BITS) do
                sh = BigFloat.(sL)
                (mg_M6(sh, 1), accumulate(+, sh))                                      # post-processing (M6); fbar_pre[n] = sum_{k=1..n} s(k)
            end
            for (ms, fs, ss) in zip(msn, fsn, ssn)
                n = ms.n
                @test fs.n == n && ss.n == n && X1[n].n == n && X2[n].n == n
                S1, S2 = ms.sum[1] + ms.comp[1], ms.sum[2] + ms.comp[2]
                # (a)
                @test ms.sum == fs.sum && ms.comp == fs.comp && ms.sum == ss.sum && ms.comp == ss.comp
                @test S1 == fs.sum[1] + fs.comp[1] && S2 == ss.sum[2] + ss.comp[2]
                @test S1 == X1[n].sum[1] + X1[n].comp[1] && S2 == X2[n].sum[1] + X2[n].comp[1]
                tol = 16 * epsT * max(1, smax[n])   # Amendment D74 item 3
                # (b)
                eb = (mg_err(ms.Y[1], post.Y[n]), mg_err(ms.Ybar[1], post.Ybar[n]), mg_err(ms.Ybar20[1], post.Ybar20[n]))
                @test all(<=(tol), eb)
                worst_b = max(worst_b, maximum(eb) / tol)
                # (c) FLIbar(n) = (1/n) sum_{k=0..n} FLI(k) (Eq. 3), FLI(0) = 0
                FLI_n = sL[n]
                a2, a2f, remd = setprecision(BigFloat, MG_BITS) do
                    fbar = fbar_pre[n] / n
                    fbarS = (fbar_pre[n] - BigFloat(FLI_n) + BigFloat(S1)) / n
                    (2 * (BigFloat(FLI_n) - fbar) + 2 * BigFloat(FLI_n) / n, 2 * (BigFloat(S1) - fbarS) + 2 * BigFloat(S1) / n,
                     BigFloat(ms.Y[1]) - 2 * (BigFloat(FLI_n) - fbar))
                end
                @test mg_err(ms.Y[1], a2) <= tol
                @test mg_err(ms.Y[1], a2f) <= tol
                worst_c = max(worst_c, mg_err(ms.Y[1], a2) / tol, mg_err(ms.Y[1], a2f) / tol)
                @test mg_err(remd, 2 * X1[n].lambda[1]) <= tol
                # (d) with FLI_1 the running maximum of the FLI reducer of the same chain
                D = fs.fli1 - FLI_n
                tolD = 16 * epsT * max(1, smax[n], fs.fli1)
                @test D >= -tolD
                @test ms.Y[1] + 2 * D >= -tolD
                a3 = setprecision(BigFloat, MG_BITS) do
                    fbar = fbar_pre[n] / n
                    f = BigFloat(fs.fli1) - BigFloat(D)
                    2 * (f - fbar) + 2 * f / n
                end
                @test mg_err(ms.Y[1], a3) <= tolD
                # negative control: the running maximum in place of FLI(n)
                if D > 10 * tolD
                    bad = setprecision(BigFloat, MG_BITS) do
                        fbar = fbar_pre[n] / n
                        2 * (BigFloat(fs.fli1) - fbar) + 2 * BigFloat(fs.fli1) / n
                    end
                    ncontrol += 1
                    nfail_control += mg_err(ms.Y[1], bad) > tol
                end
                if haskey(yd, n)
                    push!(yd[n], ms.Y[1] + 2 * D)
                    n == 1_000 && push!(lab, (Y20 = ms.Ybar20[1] > 0.5, Ybar = ms.Ybar[1] > 2, L1 = sL[n] > log(n), FLIB = fs.fli > log(n)))
                end
            end
        end
        @info "AC7 realized ratios error / (8 eps max(1, max |s_1|)) of (b) and (c) (R2)" label worst_b worst_c
        for n in (100, 1_000, 10_000)
            @info "AC7 reported (R1): Y_1 + 2 D_1 (min, median, max)" label n v = (minimum(yd[n]), mg_median(yd[n]), maximum(yd[n]))
        end
        @info "AC7 reported (R1): labels at n = 1e3 (Ybar20 > 0.5, Ybar > 2, n X_L1 > ln n, FLI_B > ln n), counts of true and of orbits where the four disagree" label counts = (count(l -> l.Y20, lab), count(l -> l.Ybar, lab), count(l -> l.L1, lab), count(l -> l.FLIB, lab)) disagree = count(l -> !(l.Y20 == l.Ybar == l.L1 == l.FLIB), lab)
    end
    @test ncontrol > 0 && nfail_control == ncontrol                                  # the control fails wherever D > 0
end

@testset "AC8 growth laws" begin
    _ = (MEGNO, Benettin)
    # (a) Regular class K = 0.5, 20 starts, Float64, seeded frames, n = 1e5: Ybar_1 in [1.5, 2.5] and Ybar20_1 in [0.4, 0.6] for at least 16 of the 20 starts; the
    # others are reported (transient near an elliptic point, F3). (b) Chaotic K in {2 pi, 5, 10}, 1000 starts, Float64: with m_Ybar, m_Y, m_20 the least-squares
    # slopes of Ybar_1, Y_1 and Ybar20_1 against n over n in [20, 200] and X(1e4) the L1 value of a k = 1 run on the identical orbit (first seeded column): the
    # ensemble means of 2 m_Ybar, m_Y and 4 m_20 each within 4 combined standard errors (L1 S5, sample spread of the scaled slope) of mean X, plus the allowance 4/180
    # (first two) or 8/180 (third) (proposals, K18). Reported: R3 (sanity rates against the Chirikov law and the K = 2 pi row of benchmarks.md, never pass/fail), R4 (the
    # first handoff at which Ybar_1 reaches 3).
    c = mg_case(:half, Float64, 20)
    run = mg_run((MEGNO(),), c.sys, c.x0, c.W0, 100_000, 1; steps = [100_000])
    ok = [(s = r[1].snaps[1]; 1.5 <= s.Ybar[1] <= 2.5 && 0.4 <= s.Ybar20[1] <= 0.6) for r in run]
    @info "AC8(a) reported: Ybar_1 and Ybar20_1 at n = 1e5, starts outside the intervals" Ybar = [r[1].snaps[1].Ybar[1] for r in run] Ybar20 = [r[1].snaps[1].Ybar20[1] for r in run] outside = findall(!, ok)
    @test mg_nbad(r -> r[1].snaps[1].n == 100_000, run) == 0
    @test count(ok) >= 16
    m = 1000
    sanity = Dict(:five => log(5 / 2), :ten => log(10 / 2), :two_pi => 1.172)        # benchmarks.md rows "large K" (ln(K/2)) and "single chaotic orbit" (K = 2 pi)
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = mg_seed(8, ci)
        sys = StandardMap(mg_K(Float64, label))
        x0 = [mg_cast(Float64, s) for s in mg_starts(StableRNG(seed), m)]
        W0 = [mg_frame(seed, j, Float64) for j in 1:m]
        w = collect(20:200)
        run = mg_run((MEGNO(),), sys, x0, W0, 200, 1; steps = collect(1:200))
        l1 = mg_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        X = [(r[1].snaps[1].sum[1] + r[1].snaps[1].comp[1]) / 10_000 for r in l1]
        @test mg_nbad(r -> length(r[1].snaps) == 200 && r[1].snaps[200].n == 200, run) == 0
        mYb = [mg_slope(Float64.(w), Float64[s.Ybar[1] for s in r[1].snaps[20:200]]) for r in run]
        mY = [mg_slope(Float64.(w), Float64[s.Y[1] for s in r[1].snaps[20:200]]) for r in run]
        m20 = [mg_slope(Float64.(w), Float64[s.Ybar20[1] for s in r[1].snaps[20:200]]) for r in run]
        @test mg_agree(2 .* mYb, X; allowance = 4 / 180)
        @test mg_agree(mY, X; allowance = 4 / 180)
        @test mg_agree(4 .* m20, X; allowance = 8 / 180)
        hit = [(h = findfirst(s -> s.Ybar[1] >= 3, r[1].snaps); h === nothing ? nothing : r[1].snaps[h].n) for r in run]
        @info "AC8 reported (R3, R4)" label mean_X = mg_mean(X) mean_2m_Ybar = mg_mean(2 .* mYb) mean_m_Y = mg_mean(mY) mean_4m_20 = mg_mean(4 .* m20) sanity = sanity[label] first_handoff_Ybar_3_median = mg_median([n for n in hit if n !== nothing]) never_reached = count(isnothing, hit)
    end
end

@testset "AC8 R5 initial vector sensitivity" begin
    _ = (MEGNO, check_frame)
    # Recorded line R5 (F3, M9): 20 regular-class starts at K = 0.5, each with 20 supplied unit vectors of uniformly random orientation (the frame is the vector and
    # its rotation by 90 degrees; frame check), the distribution of Ybar_1(1e4) and of the first n at which |Ybar_1 - 2| < 0.25, and the same for the two seeded
    # columns. Nothing is asserted beyond the frame check and the run length.
    c = mg_case(:half, Float64, 20)
    rng = StableRNG(mg_seed(8, 10))
    stats(sn) = (Ybar = sn[end].Ybar[1], first = (h = findfirst(s -> abs(s.Ybar[1] - 2) < 0.25, sn); h === nothing ? 0 : sn[h].n))
    supplied, seeded = NamedTuple[], NamedTuple[]
    for j in 1:20
        for _ in 1:20
            th = 2pi * rand(rng)
            W = SMatrix{2, 2, Float64}(cos(th), sin(th), -sin(th), cos(th))
            @test check_frame(MEGNO(), W; orbit = j) <= 2 * mg_eps(Float64)
            sn = mg_run((MEGNO(),), c.sys, [c.x0[j]], [W], 10_000, 1; steps = collect(10:10:10_000))[1][1].snaps
            push!(supplied, stats(sn))
        end
        sn = mg_run((MEGNO(),), c.sys, [c.x0[j]], [c.W0[j]], 10_000, 1; steps = collect(10:10:10_000))[1][1].snaps
        @test sn[end].n == 10_000
        push!(seeded, stats(sn))
    end
    q(v) = (minimum(v), mg_median(v), maximum(v))
    @info "AC8 R5 reported: supplied random unit vectors, Ybar_1(1e4) (min, median, max), first n with |Ybar_1 - 2| < 0.25 (0 = never)" Ybar = q([s.Ybar for s in supplied]) first = q([s.first for s in supplied])
    @info "AC8 R5 reported: the seeded frames, column 1" Ybar = q([s.Ybar for s in seeded]) first = q([s.first for s in seeded])
end

@testset "AC8 R6 Mestre orbits" begin
    _ = (MEGNO,)
    # Recorded line R6: Mestre2011 Sect. 5, standard map at K = 0.90, starts (x, p) = (0.628, pi/100) (chaotic) and (1.884, pi/100) (regular) in this repository's (x, p)
    # ordering, n up to 2000, tau = 1: Y_1/2, FLI - FLIbar (FLI = s_1) and the remainder (1/n) FLI at n in {10, 100, 1e3, 2e3}. Nothing is asserted beyond the run length
    # (the paper prints curves only, no numbers).
    sys = StandardMap(0.9)
    for (label, x) in (("chaotic", 0.628), ("regular", 1.884))
        sn = mg_run((MEGNO(),), sys, [SVector(x, pi / 100)], [SMatrix{2, 2, Float64}(1, 0, 0, 1)], 2000, 1; all = true)[1][1].snaps
        @test length(sn) == 2000
        S1 = [s.sum[1] + s.comp[1] for s in sn]
        for n in (10, 100, 1_000, 2_000)
            fbar = (sum(S1[1:n])) / n                                                      # (1/n) sum_{k=0..n} FLI(k), FLI(0) = 0
            @info "AC8 R6 reported" label n Y_half = sn[n].Y[1] / 2 FLI_minus_FLIbar = S1[n] - fbar remainder = S1[n] / n
        end
    end
end
