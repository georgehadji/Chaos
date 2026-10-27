# AC8 to AC10: the algebra of D (SN5, SN6, SN10, SN13, I4, I10), the window length (SN16), the sample-size law of D (SN11).
# Every testset starts with a line naming the new symbols it uses (UndefVarError when missing).

@testset "AC8a identical columns, exchange of the columns, sign change" begin
    _ = (StretchingNumbers, stretching_range)
    # I4, I10. (a) The reducer with two identical supplied columns (the caller's frame check bypassed here) returns D = 0 at every evaluation handoff, bitwise, and the two count columns are
    # equal. (b) Exchange of the two columns: D bitwise equal, counts, sums and moments exchanged; a sign change of a start vector (column 1, column 2, both): every output bitwise equal.
    # K = 2 pi and the regular class, 5 starts, tau in {1, 4}, n = 1000, evaluation handoffs at the multiples of the interval, both types.
    for T in SN_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = sn_case(label, T, 5)
        r = sn_red(c.sys, 1000, tau)
        evals = [n for n in sn_eval_steps(1000, sn_interval(tau)) if n % tau == 0]
        run(Ws) = sn_run((r,), c.sys, c.x0, Ws, 1000, tau; steps = evals)
        base = run(c.W0)
        same = run([SMatrix{2, 2, T}(W[1, 1], W[2, 1], W[1, 1], W[2, 1]) for W in c.W0])
        swapped = run([SMatrix{2, 2, T}(W[1, 2], W[2, 2], W[1, 1], W[2, 1]) for W in c.W0])
        neg1 = run([SMatrix{2, 2, T}(-W[1, 1], -W[2, 1], W[1, 2], W[2, 2]) for W in c.W0])
        neg2 = run([SMatrix{2, 2, T}(W[1, 1], W[2, 1], -W[1, 2], -W[2, 2]) for W in c.W0])
        negb = run([-W for W in c.W0])
        nbad = 0
        for j in 1:5, i in eachindex(evals)
            fb, fs, fx = base[j][1].snaps[i], same[j][1].snaps[i], swapped[j][1].snaps[i]
            nbad += !(fs.D === zero(T) && fs.counts[:, 1] == fs.counts[:, 2] && fs.sum[1] == fs.sum[2] && fs.comp[1] == fs.comp[2])
            nbad += !(fx.D === fb.D && fx.counts[:, 1] == fb.counts[:, 2] && fx.counts[:, 2] == fb.counts[:, 1] && fx.sum == reverse(fb.sum) && fx.comp == reverse(fb.comp) &&
                      fx.sum2 == reverse(fb.sum2) && fx.mean == reverse(fb.mean) && fx.std == reverse(fb.std) && fx.below == reverse(fb.below) && fx.above == reverse(fb.above))
            for other in (neg1, neg2, negb)
                nbad += !isequal(other[j][1].snaps[i], fb)
            end
        end
        @test nbad == 0
    end
end

@testset "AC8c the pure distance function" begin
    _ = (spectral_distance,)
    # SN5, SN6, I4. 10^4 random pairs of integer count vectors with equal totals (1000 for Double64), B in {10, 100, 1000} (vectors of B + 2 cells, h = 2 / B): D equals the plain formula of SN5
    # in BigFloat within 4 eps relative, is bitwise symmetric, satisfies the triangle inequality within 4 eps (relative to the sum), lies in [0, sqrt(2 / h)] (within 4 eps), and attains
    # sqrt(2 / h) for two different single-cell vectors; D of equal vectors is exactly 0; unequal totals or lengths are an ArgumentError; the total 0 gives 0.
    rng = StableRNG(sn_seed(8, 3))
    function randcounts(k, ncell)
        v = zeros(Int, ncell)
        for _ in 1:k
            v[min(rand(rng, 1:ncell), rand(rng, 1:ncell))] += 1       # a skewed distribution over the cells
        end
        return v
    end
    for T in SN_TYPES, B in (10, 100, 1000)
        epsT = sn_eps(T)
        h = T(2) / B
        npairs = T === Float64 ? 10_000 : 1_000
        hmax = setprecision(BigFloat, SN_BITS) do
            sqrt(2 / sn_big(h))
        end
        nbad_ref = nbad_sym = nbad_tri = nbad_rng = 0
        for _ in 1:npairs
            k = rand(rng, 20:2000)
            N1, N2, N3 = randcounts(k, B + 2), randcounts(k, B + 2), randcounts(k, B + 2)
            D12, D13, D23 = spectral_distance(N1, N2, h), spectral_distance(N1, N3, h), spectral_distance(N2, N3, h)
            ref = sn_D_h(N1, N2, h)
            nbad_ref += !(typeof(D12) === T && (ref == 0 ? D12 == 0 : sn_err(D12, ref) <= 4 * epsT * Float64(ref)))
            nbad_sym += !(D12 === spectral_distance(N2, N1, h))
            nbad_tri += !(sn_big(D13) <= (sn_big(D12) + sn_big(D23)) * (1 + 4 * BigFloat(epsT)))
            nbad_rng += !(0 <= D12 && sn_big(D12) <= hmax * (1 + 4 * BigFloat(epsT)))
        end
        @test nbad_ref == 0
        @test nbad_sym == 0
        @test nbad_tri == 0
        @test nbad_rng == 0
        # extremes
        for k in (1, 7, 1000)
            a, b = zeros(Int, B + 2), zeros(Int, B + 2)
            a[3] = k
            b[B + 1] = k
            @test sn_err(spectral_distance(a, b, h), hmax) <= 4 * epsT * Float64(hmax)
            @test spectral_distance(a, a, h) == 0
        end
        # a count difference: one sample moved between two cells changes S by 2 and D to sqrt(2) / (k sqrt h)
        a = zeros(Int, B + 2); a[4] = 10
        b = copy(a); b[4] = 9; b[5] = 1
        @test sn_err(spectral_distance(a, b, h), sqrt(BigFloat(2)) / (10 * sqrt(sn_big(h)))) <= 4 * epsT * sqrt(2) / (10 * sqrt(Float64(h)))
        @test spectral_distance(zeros(Int, B + 2), zeros(Int, B + 2), h) == 0
        @test sn_error(() -> spectral_distance([1, 2, 3], [1, 2, 4], h)) isa ArgumentError            # unequal totals
        @test sn_error(() -> spectral_distance([1, 2, 3], [1, 2, 3, 0], h)) isa Exception              # unequal lengths
    end
end

@testset "AC8d alignment bound" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # SN10. D_k <= sqrt 2 M_k / (k sqrt h) (1 + 4 eps) at every evaluation handoff of the panel of criterion 4 (sn_panel(): D38 sizes), M_k the number of handoffs i <= k at which the cells of
    # z_{1,i} and z_{2,i} differ, from the BigFloat reference cells; the samples within their half-width of an edge (U_1 + U_2 of SN14) are added to M_k, so that the bound is the exact one for the
    # computed cells. tau in {1, 4}, B = 1000, Float64.
    T = Float64
    epsT = sn_eps(T)
    B = 1000
    nbad = nlines = 0
    for (label, N, m) in sn_panel(), tau in (1, 4)
        c = sn_case(label, T, m)
        refs = sn_refs(label, T, m, N)
        lo, hi = stretching_range(c.sys)
        h = setprecision(BigFloat, SN_BITS) do
            (sn_big(hi) - sn_big(lo)) / B
        end
        evals = [n for n in sn_eval_steps(N, sn_interval(tau)) if n % tau == 0]
        res = sn_run((sn_red(c.sys, B, tau),), c.sys, c.x0, c.W0, N, tau; steps = evals)
        K = N ÷ tau
        for j in 1:m
            zs = [sn_z_from_s(refs[j].s1, tau), sn_z_from_s(refs[j].s2, tau)]
            c1, c2 = sn_cells_ref(zs[1], lo, hi, B), sn_cells_ref(zs[2], lo, hi, B)
            Mk = cumsum(c1 .!= c2)
            eps = [sn_halfwidths(refs[j], i, epsT, tau, K, lo, hi) for i in 1:2]
            env = [sn_envelope(zs[i], eps[i], lo, hi, B, [n ÷ tau for n in evals]) for i in 1:2]
            for (f, n) in zip(res[j][1].snaps, evals)
                k = n ÷ tau
                nlines += 1
                bound = sqrt(BigFloat(2)) * (Mk[k] + env[1][k].U + env[2][k].U) / (k * sqrt(h)) * (1 + 4 * BigFloat(epsT))
                nbad += !(sn_big(f.D) <= bound)
            end
        end
    end
    @test nlines > 0
    @test nbad == 0
end

@testset "AC8e pooling" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # SN13. The pooled count vector of 8 orbits (K = 2 pi, N = 200, B = 100) equals the integer sum, independent of the order of the orbits (seeded shuffles); the pooled total is 8 k; the pooled D
    # equals the formula applied to the sum (the integer sum of squares exact, D within 4 eps).
    T = Float64
    epsT = sn_eps(T)
    c = sn_case(:two_pi, T, 8)
    lo, hi = stretching_range(c.sys)
    B = 100
    h = (hi - lo) / B
    res = sn_run((sn_red(c.sys, B, 1),), c.sys, c.x0, c.W0, 200, 1; steps = [200])
    cnts = [(sn_cnt(res[j][1].snaps[1], 1), sn_cnt(res[j][1].snaps[1], 2)) for j in 1:8]
    rng = StableRNG(sn_seed(8, 5))
    ref1 = sum(first, cnts)
    ref2 = sum(last, cnts)
    @test sum(ref1) == 1600 && sum(ref2) == 1600
    for _ in 1:5
        perm = sortperm(rand(rng, 8))
        p1 = foldl((a, j) -> a .+ cnts[j][1], perm; init = zeros(Int, B + 2))
        p2 = foldl((a, j) -> a .+ cnts[j][2], perm; init = zeros(Int, B + 2))
        @test p1 == ref1 && p2 == ref2
        D = spectral_distance(p1, p2, h)
        Dref = sn_D_h(ref1, ref2, h)
        @test sn_S(p1, p2) == sn_S(ref1, ref2)
        @test D === spectral_distance(ref1, ref2, h) && sn_err(D, Dref) <= 4 * epsT * Float64(Dref)
    end
end

@testset "AC9 window length" begin
    _ = (StretchingNumbers, stretching_range)
    # SN16. tau in {1, 2, 4, 8, 16}, K in {0.5 (regular class), 2 pi}, 20 starts (Float64; Double64 6, D38), identical orbit and start vectors, n = 1024: (a) |A_j^(tau)(n) - A_j^(1)(n)| <=
    # B_j^(tau)(n) + B_j^(1)(n) with B_j(n) = 8 eps Psi_n^(j) (L1 criterion 6 form); (b) for tau' = m tau (every pair of the set) and every k', |z^(tau')_{k'} - the mean of the m fine values| <=
    # eps^(tau')_{k'} + (1/m) sum of the m fine half-widths (SN14 half-widths from the Psi grid maximum); (c) the standard deviation of the coarse series is at most that of the fine series plus
    # the sum of the two largest half-widths of the fine series. Reported, not gated: D at n = 1024 against tau (R2). z series: the differences of the module's compensated totals.
    taus = (1, 2, 4, 8, 16)
    N = 1024
    nbad_a = nbad_b = nbad_c = 0
    worst_a = worst_b = 0.0
    for T in SN_TYPES, label in (:half, :two_pi)
        epsT = sn_eps(T)
        m = T === Float64 ? 20 : 6
        c = sn_case(label, T, m)
        refs = sn_refs(label, T, m, N)
        lo, hi = stretching_range(c.sys)
        runs = Dict(tau => sn_run((sn_red(c.sys, 1000, tau),), c.sys, c.x0, c.W0, N, tau; steps = [N], light = true) for tau in taus)
        Ds = Dict(tau => [Float64(runs[tau][j][1].snaps[1].D) for j in 1:m] for tau in taus)
        @info "AC9 D at n = 1024 against tau (median over starts, R2)" T label D_tau_1_2_4_8_16 = string([round(sort(Ds[tau])[(m + 1) ÷ 2]; sigdigits = 4) for tau in taus])
        for j in 1:m
            ref = refs[j]
            zs = Dict((tau, i) => sn_z_from_lights(runs[tau][j][1].lights, i, tau) for tau in taus, i in 1:2)
            hw = Dict((tau, i) => sn_halfwidths(ref, i, epsT, tau, N ÷ tau, lo, hi) for tau in taus, i in 1:2)
            for i in 1:2
                Psi = sn_Psi(ref, i, N)
                A(tau) = setprecision(BigFloat, SN_BITS) do
                    f = runs[tau][j][1].snaps[1]
                    sn_big(f.sum[i]) + sn_big(f.comp[i])
                end
                for tau in taus[2:end]
                    nbad_a += !(abs(A(tau) - A(1)) <= 2 * 8 * epsT * Psi)
                    worst_a = max(worst_a, Float64(abs(A(tau) - A(1))) / (2 * 8 * epsT * Psi))
                end
                for (a, tau) in enumerate(taus), tau2 in taus[(a + 1):end]
                    mm = tau2 ÷ tau
                    zf, zc = zs[(tau, i)], zs[(tau2, i)]
                    hf, hc = hw[(tau, i)], hw[(tau2, i)]
                    setprecision(BigFloat, SN_BITS) do
                        for k2 in eachindex(zc)
                            blk = ((k2 - 1) * mm + 1):(k2 * mm)
                            nbad_b += !(abs(zc[k2] - sum(zf[blk]) / mm) <= hc[k2] + sum(hf[blk]) / mm)
                            worst_b = max(worst_b, Float64(abs(zc[k2] - sum(zf[blk]) / mm) / (hc[k2] + sum(hf[blk]) / mm)))
                        end
                        sdf = sqrt(sum((z - sum(zf) / length(zf))^2 for z in zf) / length(zf))
                        sdc = sqrt(sum((z - sum(zc) / length(zc))^2 for z in zc) / length(zc))
                        two = sum(sort(hf; rev = true)[1:2])
                        nbad_c += !(sdc <= sdf + two)
                    end
                end
            end
        end
    end
    @test nbad_a == 0
    @test nbad_b == 0
    @test nbad_c == 0
    @info "AC9 realized ratios (a) |A_tau - A_1| / (16 eps Psi_n) and (b) |z_coarse - mean fine| / allowance" worst_a worst_b
end

@testset "AC10 sample-size law of D" begin
    _ = (spectral_distance,)
    # SN11. The pure distance function on seeded synthetic independent multinomial samples (exact inverse-cdf draws) with Dirichlet(1, ..., 1) cell probabilities, B in {20, 100}, Z in {200, 2000,
    # 20000}, h = 1/20, the cases p1 = p2 and p1 different from p2, 4000 repeats each (1000 at Z = 20000 in the package suite, D38): the mean of D^2 is within 4 standard errors of
    # D_true^2 + (2 - sum p1^2 - sum p2^2) / (Z h), D_true^2 = sum (p1 - p2)^2 / h. For p1 = p2 the fraction of repeats with D = 0 is reported (it is 0 for Z >= 200 and the cases tested).
    rng = StableRNG(sn_seed(10, 1))
    h = 0.05
    for B in (20, 100), Z in (200, 2000, 20000), case in (:same, :diff)
        reps = (Z == 20000 && !SN_FULL) ? 1000 : 4000
        p1 = sn_dirichlet(rng, B)
        p2 = case === :same ? p1 : sn_dirichlet(rng, B)
        c1, c2 = cumsum(p1), cumsum(p2)
        D2 = Vector{Float64}(undef, reps)
        zeros_D = 0
        for r in 1:reps
            D = spectral_distance(sn_multinomial(rng, Z, c1), sn_multinomial(rng, Z, c2), h)
            D2[r] = D^2
            zeros_D += D == 0
        end
        expected = sum((p1 .- p2) .^ 2) / h + (2 - sum(abs2, p1) - sum(abs2, p2)) / (Z * h)
        se = sqrt(sn_var(D2) / reps)
        @test abs(sn_mean(D2) - expected) <= 4 * se
        case === :same && @info "AC10 fraction of repeats with D = 0 (p1 = p2)" B Z fraction = zeros_D / reps
        case === :same && @test zeros_D == 0
    end
end
