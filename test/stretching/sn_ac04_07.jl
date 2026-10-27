# AC4 to AC7: identities with L1 and the other column-wise reducers (SN2, I2, F3), conservation, normalisation and moments (I1, I3, SN4, SN7), the reference counts and the D envelope
# (SN14, SN15), range, overflow and parameters (SN12, SN4, P3, P4). Every testset starts with a line naming the new symbols it uses (UndefVarError when missing).
# Sizes follow D38 (sn_panel(): 6 starts to n = 1e3 and 2 starts to n = 1e4 for the reference criteria; SN_FULL: the 20 starts to n = 1e4 of the spec).

# Handoff indices k that are multiples of 10 up to N / tau, as map steps n = k tau (criterion 5: "every handoff that is a multiple of 10").
sn_every10(tau::Int, N::Int) = [10 * tau * j for j in 1:(N ÷ (10 * tau))]

@testset "AC4a column sums are the L1 sums and the mean of z is the L1 value" begin
    _ = (StretchingNumbers, stretching_range)
    # SN2, I2. K in {0.5 (regular class), 2 pi, 10}, 20 starts each, T in {Float64, Double64}, tau in {1, 4}, n in {100, 1e3, 1e4}, wrapped, the seeded pair as the two columns: (a) for each
    # column the sum A_j and its correction term are bitwise equal to the Benettin (L1) single-column run on the identical orbit with that column supplied as the frame; (b) the mean of
    # z equals the L1 value X(n) = total(A) / n within 2 eps relative. Negative control: a double that accumulates z instead of ln alpha has a sum different from the L1 sum (tau = 4:
    # the sum is divided by 4).
    ns = (100, 1000, 10_000)
    nbad_a = nbad_b = 0
    for T in SN_TYPES, label in (:half, :two_pi, :ten), tau in (1, 4)
        epsT = sn_eps(T)
        c = sn_case(label, T, 20)
        r = sn_red(c.sys, 1000, tau)
        for j in 1:20
            W0 = c.W0[j]
            sn = sn_run((r,), c.sys, [c.x0[j]], [W0], 10_000, tau; steps = collect(ns))[1][1].snaps
            b1 = sn_run((Benettin(),), c.sys, [c.x0[j]], [SMatrix{2, 1, T}(W0[1, 1], W0[2, 1])], 10_000, tau; steps = collect(ns))[1][1].snaps
            b2 = sn_run((Benettin(),), c.sys, [c.x0[j]], [SMatrix{2, 1, T}(W0[1, 2], W0[2, 2])], 10_000, tau; steps = collect(ns))[1][1].snaps
            for (i, n) in enumerate(ns)
                f = sn[i]
                for (col, b) in ((1, b1[i]), (2, b2[i]))
                    nbad_a += !(f.n == n && b.n == n && f.sum[col] == b.sum[1] && f.comp[col] == b.comp[1])
                    nbad_b += !(abs(f.mean[col] - b.lambda[1]) <= 2 * epsT * abs(b.lambda[1]))
                end
            end
        end
    end
    @test nbad_a == 0
    @test nbad_b == 0
    # the negative control: accumulating z instead of ln alpha (tau = 4) gives a different sum
    T = Float64
    c = sn_case(:two_pi, T, 2)
    r = sn_red(c.sys, 1000, 4)
    res = sn_run((r,), c.sys, [c.x0[1]], [c.W0[1]], 1000, 4; steps = [1000], light = true)[1][1]
    zacc = sum(sn_z_from_lights(res.lights, 1, 4))
    @test T(zacc) != res.snaps[1].sum[1] + res.snaps[1].comp[1]
end

@testset "AC4b chains with FLI and MEGNO" begin
    _ = (StretchingNumbers, stretching_range)
    # SN2, F3. In the chains {SN (reader), FLI}, {FLI (reader), SN} and {SN (reader), MEGNO} the column-1 sums of the stretching reducer and of the other reducer (sum and correction term)
    # are bitwise equal in the same run (n = 1e3, 20 starts per class, both types, tau in {1, 4}; D38: the structural check does not depend on n).
    nbad = 0
    for T in SN_TYPES, label in (:half, :two_pi, :ten), tau in (1, 4)
        c = sn_case(label, T, 20)
        for j in 1:20
            for (chain, ia, ib) in (((sn_red(c.sys, 1000, tau; replace = false), FLI()), 1, 2), ((FLI(replace = false), sn_red(c.sys, 1000, tau)), 2, 1),
                                    ((sn_red(c.sys, 1000, tau; replace = false), MEGNO()), 1, 2))
                res = sn_run(chain, c.sys, [c.x0[j]], [c.W0[j]], 1000, tau; steps = [1000])[1]
                fa, fb = res[ia].snaps[1], res[ib].snaps[1]
                nbad += !(fa.sum[1] == fb.sum[1] && fa.comp[1] == fb.comp[1])
            end
        end
    end
    @test nbad == 0
end

@testset "AC5 conservation, normalisation and moments" begin
    _ = (StretchingNumbers, stretching_range, spectrum_density)
    # I1, I3, SN4, SN7. The panel of criterion 4 (K in {0.5, 2 pi, 10}, 20 starts (D38: 10 in Float64, 4 in Double64), tau in {1, 4}, both types), B in {100, 1000}, N = 1000 (SN_FULL: 1e4), at every handoff k that is a multiple
    # of 10: (a) the counts of the cells 0..B + 1 of each column sum to k exactly and no count decreases; (b) with no overflow sample, |histogram mean - mean z| <= h/2 and |histogram
    # standard deviation - standard deviation of z| <= h/2, the mean and the standard deviation of z from the BigFloat post-processing of the stored series (the differences of the
    # module's compensated totals); (c) the normalised spectrum sums to 1 within 4 eps (sum of the module's S_c in BigFloat times h), each S_c equals N_c / (k h) within 2 eps relative;
    # (d) on T1 the standard deviation from the two compensated sums equals the closed form |c p_0 tau| sqrt((k^2 - 1) / 12) within 3 sqrt(eps) (|mean| + sd); on the panel the ratio
    # to the BigFloat value is reported. Negative control: a double that drops a sample outside the last cell fails (a).
    N = SN_FULL ? 10_000 : 1000
    nbad_a = nbad_b = nbad_c = ncheck = 0
    worst_sd = 0.0
    for T in SN_TYPES, label in (:half, :two_pi, :ten), tau in (1, 4), B in (100, 1000)
        epsT = sn_eps(T)
        m = SN_FULL ? 20 : (T === Float64 ? 10 : 4)               # D38: 10 starts in Float64, 4 in Double64
        c = sn_case(label, T, m)
        lo, hi = stretching_range(c.sys)
        h = setprecision(BigFloat, SN_BITS) do
            (sn_big(hi) - sn_big(lo)) / B
        end
        hF = Float64(h)
        hT = (hi - lo) / B
        ks = sn_every10(tau, N)
        r = sn_red(c.sys, B, tau)
        res = sn_run((r,), c.sys, c.x0, c.W0, N, tau; steps = ks, light = true)
        for j in 1:m
            snaps, lights = res[j][1].snaps, res[j][1].lights
            @test length(snaps) == length(ks)
            zs = [sn_z_from_lights(lights, col, tau) for col in 1:2]
            moments = [sn_prefix_moments(zs[col], [n ÷ tau for n in ks]) for col in 1:2]
            prev = [zeros(Int, B + 2), zeros(Int, B + 2)]
            for (f, n) in zip(snaps, ks)
                k = n ÷ tau
                ncheck += 1
                for col in 1:2
                    cnt = sn_cnt(f, col)
                    nbad_a += !(sum(cnt) == k && all(cnt .>= prev[col]) && f.handoffs == k)
                    prev[col] = cnt
                    if f.below[col] == 0 && f.above[col] == 0
                        m1, sd1 = sn_hist_moments(cnt, Float64(lo), hF)
                        mref, sdref = moments[col][k]
                        nbad_b += !(abs(m1 - Float64(mref)) <= hF / 2 * (1 + 1e-9) && abs(sd1 - Float64(sdref)) <= hF / 2 * (1 + 1e-9))
                    else
                        nbad_b += 1                       # the default range is padded: no overflow sample may appear (also AC7a)
                    end
                    S = spectrum_density(cnt, hT)
                    ok = length(S) == B
                    ssum, refs_c = setprecision(BigFloat, SN_BITS) do
                        (sum(sn_big(s) for s in S) * sn_big(hT), [BigFloat(cnt[c1 + 1]) / (k * sn_big(hT)) for c1 in 1:B])
                    end
                    ok &= abs(ssum - 1) <= 4 * epsT
                    ok &= all(c1 -> sn_err(S[c1], refs_c[c1]) <= 2 * epsT * Float64(refs_c[c1]), 1:B)
                    nbad_c += !ok
                end
            end
            if j <= 3                                       # the ratio to the BigFloat standard deviation of the stored series, reported
                for col in 1:2, (f, n) in zip(snaps, ks)
                    k = n ÷ tau
                    mref, sdref = moments[col][k]
                    scale = 3 * sqrt(epsT) * (abs(Float64(mref)) + Float64(sdref))
                    worst_sd = max(worst_sd, scale == 0 ? 0.0 : sn_err(f.std[col], sdref) / scale)
                end
            end
        end
    end
    @test nbad_a == 0
    @test nbad_b == 0
    @test nbad_c == 0
    @info "AC5(d) realized ratio |std - std_ref| / (3 sqrt(eps) (|mean| + sd)) on the panel" worst_sd checked = ncheck
    # (d) T1: z_{1,k} is an arithmetic progression of spacing delta = c p_0 tau, whose population standard deviation is |delta| sqrt((k^2 - 1) / 12); p_0 in {1/50, -1/20}, tau in {1, 3, 5},
    # k in {10, 100, 1000}, B = 1001 over [-1, 1] (the progression leaves the range; the standard deviation does not depend on the histogram)
    for T in SN_TYPES, p0r in (1 // 50, -1 // 20), tau in (1, 3, 5)
        epsT = sn_eps(T)
        c, x0, p0 = exact(T, 3 // 10), exact(T, 1 // 10), exact(T, p0r)
        r = StretchingNumbers(; range = (T(-1), T(1)), bins = 1001, interval = sn_interval(tau))
        ks = (10, 100, 1000)
        fs = sn_run((r,), sn_T1(c), [SVector{2, T}(x0, p0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], 1000 * tau, tau; steps = [k * tau for k in ks])[1][1].snaps
        for (f, k) in zip(fs, ks)
            n = k * tau
            sd_cf = setprecision(BigFloat, SN_BITS) do
                abs(sn_big(c) * sn_big(p0) * tau) * sqrt(BigFloat(k^2 - 1) / 12)
            end
            mean_cf = setprecision(BigFloat, SN_BITS) do
                sn_big(c) * (sn_big(x0) + sn_big(p0) * (n - 1) / 2)
            end
            tol = 3 * sqrt(epsT) * (abs(Float64(mean_cf)) + Float64(sd_cf))
            @test sn_err(f.std[1], sd_cf) <= tol && sn_err(f.std[2], sd_cf) <= tol
        end
    end
    # negative control on plain arrays: dropping the samples outside the range breaks the sum k (T1, p_0 = 1/50, tau = 1, range [-1, 1], k = 1000: z leaves the range after about 160 handoffs)
    zs = sn_T1_z1(exact(Float64, 3 // 10), exact(Float64, 1 // 10), exact(Float64, 1 // 50), 1, 1000)
    @test sum(sn_hist_ref(zs, -1.0, 1.0, 1001)) == 1000
    @test sum(sn_hist_fault(zs, -1.0, 1.0, 1001, :drop)) < 1000
end

@testset "AC6 reference counts and the D envelope" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # SN14, SN15. K in {0.5 (regular class), 2 pi, 10} (sn_panel(): D38 sizes; SN_FULL: 20 starts to n = 1e4), tau in {1, 4}, B in {100, 1000}, both seeded columns, the BigFloat(256) reference at
    # the realized points of the run: at every evaluation handoff (multiples of the interval 100, and the last), N-_c <= N_c <= N+_c in every cell of both columns, and |D - D_ref| <=
    # sqrt 2 (U_1 + U_2) / (k sqrt h) + 4 eps D_ref (SN15). Where max eps_k >= h/2 the line is counted as unresolved and not gated. The Psi_max of SN14 is the maximum of Psi_g over the grid
    # 1..9 and 10 per decade, g <= n (a test-side reading, as rli RL11). Reported: the unresolved fraction, the realized ratio |D - D_ref| / bound. Negative control: a double with a
    # single-precision accumulator of the log leaves the envelope at n = 1e3 for the chaotic starts.
    nbad_cnt = nbad_D = nlines = nunres = ndiff = maxU = 0
    worst_ratio = 0.0
    viol32 = 0
    for T in SN_TYPES, (label, N, m) in sn_panel(), tau in (1, 4), B in (100, 1000)
        epsT = sn_eps(T)
        c = sn_case(label, T, m)
        refs = sn_refs(label, T, m, N)
        lo, hi = stretching_range(c.sys)
        h = setprecision(BigFloat, SN_BITS) do
            (sn_big(hi) - sn_big(lo)) / B
        end
        interval = sn_interval(tau)
        evals = [n for n in sn_eval_steps(N, interval) if n % tau == 0]
        r = sn_red(c.sys, B, tau)
        res = sn_run((r,), c.sys, c.x0, c.W0, N, tau; steps = evals)
        for j in 1:m
            snaps = res[j][1].snaps
            K = N ÷ tau
            zs = [sn_z_from_s(refs[j].s1, tau), sn_z_from_s(refs[j].s2, tau)]
            eps = [sn_halfwidths(refs[j], i, epsT, tau, K, lo, hi) for i in 1:2]
            resolved = [maximum(eps[i]) < Float64(h) / 2 for i in 1:2]
            env = [sn_envelope(zs[i], eps[i], lo, hi, B, [n ÷ tau for n in evals]) for i in 1:2]
            for (f, n) in zip(snaps, evals)
                k = n ÷ tau
                nlines += 1
                if !(resolved[1] && resolved[2])
                    nunres += 1
                    continue
                end
                ok = true
                for i in 1:2
                    e = env[i][k]
                    cnt = sn_cnt(f, i)
                    ok &= all(e.Nm .<= cnt) && all(cnt .<= e.Np)
                end
                nbad_cnt += !ok
                ndiff += any(sn_cnt(f, i) != env[i][k].N for i in 1:2)
                maxU = max(maxU, env[1][k].U + env[2][k].U)
                Dref = sn_D_ref(env[1][k].N, env[2][k].N, lo, hi, B)
                bound = sqrt(BigFloat(2)) * (env[1][k].U + env[2][k].U) / (k * sqrt(h)) + 4 * epsT * Dref
                eD = abs(sn_big(f.D) - Dref)
                nbad_D += !(eD <= bound)
                worst_ratio = max(worst_ratio, bound == 0 ? (eD == 0 ? 0.0 : Inf) : Float64(eD / bound))
            end
            # negative control: Float32 accumulation of the log (chaotic classes, n = 1e3, tau = 1, B = 1000, Float64)
            if T === Float64 && label !== :half && N == minimum(p[2] for p in sn_panel()) && tau == 1 && B == 1000
                s1 = refs[j].s1
                S32 = cumsum(Float32.(s1 .- vcat(zero(BigFloat), s1[1:(end - 1)])))  # smallest panel N: 1e3 default, 1e4 under SN_FULL_RUN
                z32 = Float64.(S32 .- vcat(Float32(0), S32[1:(end - 1)]))
                e1 = env[1][N]
                cnt32 = sn_hist_ref(z32, lo, hi, B)
                viol32 += !(all(e1.Nm .<= cnt32) && all(cnt32 .<= e1.Np))
            end
        end
    end
    @test nbad_cnt == 0
    @test nbad_D == 0
    @test nlines > 0 && nunres <= nlines ÷ 4
    @test viol32 >= 1
    @info "AC6 realized" lines = nlines unresolved_fraction = nunres / max(nlines, 1) lines_with_counts_different_from_the_reference = ndiff largest_U1_plus_U2 = maxU worst_D_ratio = worst_ratio float32_double_violations = viol32
end

@testset "AC6b envelope stress at B = 10^4" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # Extra, beyond the spec sizes. At the spec sizes the edge certificate of SN14 is inactive: no sample lies within its half-width of a cell edge (U = 0), so criterion 6 reduces to equality with the
    # reference counts. With B = 10^4 (h = 1.4e-4, the largest cell number of the sensitivity line R2) on the regular class at n = 1e4 the half-widths reach 6e-8 and U > 0: the same predicates as
    # AC6 (counts inside N- .. N+, |D - D_ref| <= SN15) on 2 starts, tau in {1, 4}, Float64. Reported: U and the realized ratio.
    T = Float64
    epsT = sn_eps(T)
    N, B, m = 10_000, 10_000, 2
    c = sn_case(:half, T, m)
    refs = sn_refs(:half, T, m, N)
    lo, hi = stretching_range(c.sys)
    h = setprecision(BigFloat, SN_BITS) do
        (sn_big(hi) - sn_big(lo)) / B
    end
    nbad = nbadD = nlines = maxU = nunres = ndiff = 0
    worst = 0.0
    for tau in (1, 4)
        evals = [n for n in sn_eval_steps(N, sn_interval(tau)) if n % tau == 0]
        res = sn_run((sn_red(c.sys, B, tau),), c.sys, c.x0, c.W0, N, tau; steps = evals)
        for j in 1:m
            zs = [sn_z_from_s(refs[j].s1, tau), sn_z_from_s(refs[j].s2, tau)]
            eps = [sn_halfwidths(refs[j], i, epsT, tau, N ÷ tau, lo, hi) for i in 1:2]
            resolved = all(maximum(eps[i]) < Float64(h) / 2 for i in 1:2)
            env = [sn_envelope(zs[i], eps[i], lo, hi, B, [n ÷ tau for n in evals]) for i in 1:2]
            for (f, n) in zip(res[j][1].snaps, evals)
                k = n ÷ tau
                nlines += 1
                resolved || (nunres += 1; continue)
                ok = all(i -> all(env[i][k].Nm .<= sn_cnt(f, i)) && all(sn_cnt(f, i) .<= env[i][k].Np), 1:2)
                nbad += !ok
                maxU = max(maxU, env[1][k].U + env[2][k].U)
                ndiff += any(sn_cnt(f, i) != env[i][k].N for i in 1:2)
                Dref = sn_D_ref(env[1][k].N, env[2][k].N, lo, hi, B)
                bound = sqrt(BigFloat(2)) * (env[1][k].U + env[2][k].U) / (k * sqrt(h)) + 4 * epsT * Dref
                eD = abs(sn_big(f.D) - Dref)
                nbadD += !(eD <= bound)
                worst = max(worst, bound == 0 ? (eD == 0 ? 0.0 : Inf) : Float64(eD / bound))
            end
        end
    end
    @test nbad == 0
    @test nbadD == 0
    @test nunres == 0
    @info "AC6b realized (B = 10^4)" lines = nlines unresolved = nunres lines_with_counts_different_from_the_reference = ndiff largest_U1_plus_U2 = maxU worst_D_ratio = worst
end

@testset "AC7a default range leaves nothing outside" begin
    _ = (StretchingNumbers, stretching_range)
    # SN12, P3, I6. The default range of SN12 (stretching_range): both overflow counts are 0, tolerance 0, on the panel of criterion 4 (20 starts per class) and on Y5 hyperbolic (seeded ids
    # 1..20), tau in {1, 4, 16}, N = 1008 (a multiple of 16), B = 1000, both types.
    seed = sn_seed(7, 1)
    for T in SN_TYPES, tau in (1, 4, 16)
        for label in (:half, :two_pi, :ten)
            c = sn_case(label, T, 20)
            res = sn_run((sn_red(c.sys, 1000, tau),), c.sys, c.x0, c.W0, 1008, tau; steps = [1008])
            @test sn_nbad(r -> r[1].snaps[1].below == [0, 0] && r[1].snaps[1].above == [0, 0] && r[1].snaps[1].handoffs == 1008 ÷ tau, res) == 0
        end
        y5 = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        res = sn_run((sn_red(y5, 1000, tau),), y5, [SVector{2, T}(0, 0) for _ in 1:20], [sn_frame(seed, j, T) for j in 1:20], 1008, tau; steps = [1008])
        @test sn_nbad(r -> r[1].snaps[1].below == [0, 0] && r[1].snaps[1].above == [0, 0], res) == 0
    end
end

@testset "AC7b caller range, overflow cells" begin
    _ = (StretchingNumbers, spectral_distance)
    # SN4, P3. A caller range [-0.1, 0.1], K = 10 (the default range of SN12 is 2.7), B = 100, 20 starts, n = 1e3, tau = 1, both types (Double64: 6 starts, D38): the overflow counts and all
    # cell counts lie in the envelope of SN14; the counts of the cells 0..B + 1 still sum to k; the reported below / above equal the counts of cells 0 and B + 1; D equals the independent
    # plain-array D from the reference series within SN15. A double that clips the overflow values into the end cells fails the envelope.
    B, N = 100, 1000
    lo, hi = -0.1, 0.1
    viol_clip = 0
    for T in SN_TYPES
        epsT = sn_eps(T)
        m = T === Float64 ? 20 : 6
        c = sn_case(:ten, T, m)
        refs = sn_refs(:ten, T, m, N)
        r = StretchingNumbers(; range = (T(lo), T(hi)), bins = B, interval = 100)
        res = sn_run((r,), c.sys, c.x0, c.W0, N, 1; steps = [100, 500, 1000])
        h = setprecision(BigFloat, SN_BITS) do
            (BigFloat(hi) - BigFloat(lo)) / B
        end
        for j in 1:m
            zs = [sn_z_from_s(refs[j].s1, 1), sn_z_from_s(refs[j].s2, 1)]
            eps = [sn_halfwidths(refs[j], i, epsT, 1, N, T(lo), T(hi)) for i in 1:2]
            env = [sn_envelope(zs[i], eps[i], T(lo), T(hi), B, [100, 500, 1000]) for i in 1:2]
            for (f, k) in zip(res[j][1].snaps, (100, 500, 1000))
                for i in 1:2
                    e, cnt = env[i][k], sn_cnt(f, i)
                    @test all(e.Nm .<= cnt) && all(cnt .<= e.Np) && sum(cnt) == k
                    @test f.below[i] == cnt[1] && f.above[i] == cnt[B + 2]
                end
                Dref = sn_D_ref(env[1][k].N, env[2][k].N, T(lo), T(hi), B)
                bound = sqrt(BigFloat(2)) * (env[1][k].U + env[2][k].U) / (k * sqrt(h)) + 4 * epsT * Dref
                @test abs(sn_big(f.D) - Dref) <= bound
            end
            if T === Float64
                cnt_clip = sn_hist_fault(zs[1], lo, hi, B, :clip)
                e = env[1][N]
                viol_clip += !(all(e.Nm .<= cnt_clip) && all(cnt_clip .<= e.Np))
            end
        end
    end
    @test viol_clip == 20
end

@testset "AC7c cell rule at the edges" begin
    _ = (StretchingNumbers, spectrum_cell, reducer_workspace)
    # SN4. z set exactly to s_hi, exactly to s_lo and one unit in the last place above s_hi (below s_lo) falls in cell B, cell 1 and cell B + 1 (cell 0), tolerance 0. Both through the pure
    # function spectrum_cell and through the reducer: a hand-built deviation matrix whose first column has norm 2 gives z = log(2) (tau = 1, column 2 is a unit vector, z = 0); the range
    # is chosen so that z = log(2) is s_hi, s_lo, one ulp above s_hi, one ulp below s_lo. Float64 (nextfloat / prevfloat) and Double64 (an offset of 2^-104).
    B = 10
    for T in SN_TYPES
        z0 = log(T(2))
        tiny = T === Float64 ? eps(z0) : T(ldexp(1.0, -104))
        zq = T === Float64 ? nextfloat(z0) : z0 + tiny
        zm = T === Float64 ? prevfloat(z0) : z0 - tiny
        @test zq > z0 && zm < z0
        cases = ((T(1) / 4, z0, B), (z0, T(1), 1), (T(1) / 4, zm, B + 1), (zq, T(1), 0))      # (lo, hi, expected cell of z0)
        for (lo, hi, cell) in cases
            @test spectrum_cell(z0, lo, hi, B) == cell
            r = StretchingNumbers(; range = (lo, hi), bins = B, interval = 100)
            st = reducer_init(r, T, Val(2), Val(2))
            ws = zeros(UInt32, reducer_workspace(r, T, Val(2), Val(2)))
            W = SMatrix{2, 2, T}(2, 0, 0, 1)
            st, _ = reducer_step(r, st, 1, SVector{2, T}(0, 0), W, ws)
            f = reducer_final(r, st, ws)
            cnt = sn_cnt(f, 1)
            @test cnt[cell + 1] == 1 && sum(cnt) == 1
        end
        # interior points: the cell rule floor((z - lo) / h) + 1 on a grid of exact binary values (range [0, 1], B = 8, h = 1/8)
        for (z, cell) in ((T(0), 1), (T(1) / 16, 1), (T(1) / 8, 2), (T(7) / 8, 8), (T(1), 8), (T(-1) / 1024, 0), (T(1) + T(1) / 1024, 9))
            @test spectrum_cell(z, T(0), T(1), 8) == cell
        end
    end
end

@testset "AC7d argument errors before any work" begin
    _ = (StretchingNumbers, check_run)
    # P4. An argument error naming the parameter for s_lo >= s_hi, non-finite s_lo or s_hi, B = 0 (and negative), a non-positive interval, an interval that is not a multiple of tau, and a run
    # of more than 10^9 handoffs. The constructor checks what it can see; check_run(r, tau, N) is the caller's check before the stream (the reducer protocol does not give reducer_init the
    # step size or the run length); a valid reducer passes both.
    ok = StretchingNumbers(; range = (-1.0, 1.0), bins = 100, interval = 100)
    @test check_run(ok, 1, 10_000) === nothing && check_run(ok, 4, 10_000) === nothing && check_run(ok, 100, 10_000) === nothing
    for (kw, name) in (((range = (1.0, 1.0),), "range"), ((range = (2.0, 1.0),), "range"), ((range = (NaN, 1.0),), "range"), ((range = (-1.0, Inf),), "range"),
                       ((range = (-Inf, 1.0),), "range"), ((bins = 0,), "bins"), ((bins = -3,), "bins"), ((interval = 0,), "interval"), ((interval = -100,), "interval"))
        e = sn_error(() -> StretchingNumbers(; merge((range = (-1.0, 1.0), bins = 100, interval = 100), kw)...))
        @test e isa ArgumentError && occursin(name, e.msg)
    end
    e = sn_error(() -> check_run(ok, 3, 999))                      # interval 100 is not a multiple of tau = 3
    @test e isa ArgumentError && occursin("interval", e.msg)
    e = sn_error(() -> check_run(ok, 16, 160))
    @test e isa ArgumentError && occursin("interval", e.msg)
    e = sn_error(() -> check_run(ok, 1, 1_000_000_001))
    @test e isa ArgumentError && occursin("handoffs", e.msg)
    @test check_run(ok, 1, 1_000_000_000) === nothing
    e = sn_error(() -> check_run(ok, 4, 4_000_000_004))
    @test e isa ArgumentError && occursin("handoffs", e.msg)
end

@testset "AC7e under-sampling flag" begin
    _ = (StretchingNumbers, stretching_range)
    # P4. The flag is set exactly at the evaluation handoffs with k < min_per_cell * B and cleared after; it never changes a count or a D. B = 20, min_per_cell = 10 (threshold 200 handoffs),
    # tau = 1 and 4, N = 800, evaluation handoffs at the multiples of the interval 100; the same run with min_per_cell = 1 gives bitwise equal counts and D.
    for T in SN_TYPES, tau in (1, 4)
        c = sn_case(:two_pi, T, 3)
        evals = collect(100:100:800)
        ra = sn_red(c.sys, 20, tau; min_per_cell = 10)
        rb = sn_red(c.sys, 20, tau; min_per_cell = 1)
        a = sn_run((ra,), c.sys, c.x0, c.W0, 800, tau; steps = evals)
        b = sn_run((rb,), c.sys, c.x0, c.W0, 800, tau; steps = evals)
        for j in 1:3
            for (fa, fb, n) in zip(a[j][1].snaps, b[j][1].snaps, evals)
                k = n ÷ tau
                @test fa.undersampled == (k < 10 * 20)
                @test fb.undersampled == (k < 1 * 20)
                @test isequal(fa.counts, fb.counts) && fa.D === fb.D
            end
        end
    end
end
