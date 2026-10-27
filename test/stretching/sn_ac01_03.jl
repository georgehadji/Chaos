# AC1 to AC3: the algebra of the definition on test map T1 (SN8d), the constant and hyperbolic closed forms (SN8a, SN8b), the K_std = 0 shear (SN8c, SN10, F5).
# Every testset starts with a line naming the new symbols it uses, so a missing symbol errors that testset only (UndefVarError). The counts of the cells 0..B + 1 are compared
# with tolerance 0 against the plain-array histogram of the BigFloat closed form (sn_hist_ref), D within 4 eps relative of the exact integer formula (sn_D_ref).

# Check one snapshot against closed-form series: both columns' counts equal (tolerance 0) and D within 4 eps relative of the exact formula. Returns (counts_ok, S_ok, D_ok).
function sn_check_snapshot(f, zs1, zs2, k::Int, lo, hi, B::Int, epsT::Float64)
    N1, N2 = sn_hist_ref(view(zs1, 1:k), lo, hi, B), sn_hist_ref(view(zs2, 1:k), lo, hi, B)
    counts_ok = sn_cnt(f, 1) == N1 && sn_cnt(f, 2) == N2 && size(f.counts) == (B + 2, 2)
    S_ok = sn_S(sn_cnt(f, 1), sn_cnt(f, 2)) == sn_S(N1, N2)
    Dref = sn_D_ref(N1, N2, lo, hi, B)
    D_ok = Dref == 0 ? f.D == 0 : sn_err(f.D, Dref) <= 4 * epsT * Float64(Dref)
    return counts_ok, S_ok, D_ok
end

@testset "AC1 algebra of the definition on test map T1" begin
    _ = (StretchingNumbers, spectral_distance)
    # SN8(d): test map T1 (x' = x + p, p' = p, tangent diag(exp(c x), exp(-c x)) at the pre-step point), c = 3/10, x_0 = 1/10, p_0 in {0, 1/50, -1/20} (the decimals of the spec rounded
    # once into the element type), supplied columns (1, 0) and (0, 1), tau in {1, 3, 5}, k in {1, 2, 10, 100, 1000}, B0 in {101, 1000} over [-1, 1]. (i) the counts of both columns
    # equal the closed-form counts in every cell, tolerance 0 (the closed form is z_{1,k} = c [x_0 + p_0 (tau (2k - 1) - 1) / 2], z_{2,k} = -z_{1,k}; the overflow cells are filled
    # once z leaves [-1, 1]); (ii) the integer sum of the squared count differences equals the closed form's exactly and D within 4 eps relative (p_0 = 0: D = sqrt(2 / h), the cells
    # differ); (iii) the mean of z of column 1 equals the sum of c x_j over the realized orbit (the decimals p_0 = 1/50, -1/20 do not add exactly, so the closed form c (x_0 + p_0 (n - 1) / 2)
    # of the spec differs from the realized sum by the rounding of the orbit, up to 4500 eps at n = 5000: reported) within the same-orbit bound 8 eps (n + G_first) / n of L1 S3 (K_s = 0)
    # plus 2 eps |mean| for the division; with the closed form itself on dyadic parameters, for which the orbit is exact, below (Deviation: the printed "2 eps relative" cannot hold, z = 0.03
    # is formed as the log of a number near 1, an absolute error of about eps/2 is 17 eps relative).
    # The spec's own parameters put closed-form values exactly on cell edges (p_0 = 0, B = 1000: z = 0.03 is an edge; p_0 = 1/50, B = 1000: every z), so B is the smallest B >= B0
    # for which every closed-form value lies at least 1e-9 from every cell edge (the edge condition of the spec, verified from the BigFloat series; sn_pick_B), reported.
    ks = (1, 2, 10, 100, 1000)
    bumped = Tuple[]
    nbad_counts = nbad_S = nbad_D = nbad_mean = 0
    ncases = 0
    dev_cf = 0.0
    for T in SN_TYPES
        epsT = sn_eps(T)
        c, x0 = exact(T, 3 // 10), exact(T, 1 // 10)
        lo, hi = T(-1), T(1)
        sys = sn_T1(c)
        W0 = SMatrix{2, 2, T}(1, 0, 0, 1)
        for p0r in (0 // 1, 1 // 50, -1 // 20), tau in (1, 3, 5)
            p0 = exact(T, p0r)
            zs1 = sn_T1_z1(c, x0, p0, tau, 1000)
            zs2 = -zs1
            # the realized orbit of the stream (forward is bitwise the stream's): c x_j at the stored points, summed in BigFloat; the decimals p_0 = 1/50 and -1/20 do not add exactly in T, so the
            # closed form of the spec differs from the realized sums by the rounding of the orbit itself (up to 4500 eps at n = 5000)
            cumreal = setprecision(BigFloat, SN_BITS) do
                zq = SVector{2, T}(x0, p0)
                acc = zero(BigFloat)
                out = Vector{BigFloat}(undef, 1000 * tau)
                for j in 1:(1000 * tau)
                    acc += sn_big(c) * sn_big(zq[1])
                    out[j] = acc
                    zq = forward(sys, zq)
                end
                out
            end
            for B0 in (101, 1000)
                B, bump = sn_pick_B(B0, lo, hi, zs1, zs2)
                bump && push!(bumped, (T, p0r, tau, B0, B))
                r = StretchingNumbers(; range = (lo, hi), bins = B, interval = sn_interval(tau))
                res = sn_run((r,), sys, [SVector{2, T}(x0, p0)], [W0], 1000 * tau, tau; steps = [k * tau for k in ks])[1][1].snaps
                @test length(res) == length(ks)
                for (f, k) in zip(res, ks)
                    ncases += 1
                    n = k * tau
                    cok, Sok, Dok = sn_check_snapshot(f, zs1, zs2, k, lo, hi, B, epsT)
                    nbad_counts += !cok
                    nbad_S += !Sok
                    nbad_D += !Dok
                    @test f.n == n && f.handoffs == k
                    mean_real = cumreal[n] / n
                    mean_cf = setprecision(BigFloat, SN_BITS) do
                        sn_big(c) * (sn_big(x0) + sn_big(p0) * (n - 1) / 2)
                    end
                    dev_cf = max(dev_cf, Float64(abs(mean_real - mean_cf)) / epsT)
                    G = sn_T1_G(c, x0, p0, tau)
                    tol = 8 * epsT * (n + G) / n + 2 * epsT * Float64(abs(mean_real))
                    nbad_mean += !(sn_err(f.mean[1], mean_real) <= tol && sn_err(f.mean[2], -mean_real) <= tol)
                    if p0r == 0 && k == 1000
                        h = setprecision(BigFloat, SN_BITS) do
                            (sn_big(hi) - sn_big(lo)) / B
                        end
                        @test sn_err(f.D, sqrt(2 / h)) <= 4 * epsT * Float64(sqrt(2 / h))        # two single cells, D = sqrt(2 / h)
                        @test count(!iszero, sn_cnt(f, 1)) == 1 && count(!iszero, sn_cnt(f, 2)) == 1
                    end
                end
            end
        end
    end
    @test ncases == 2 * 3 * 3 * 2 * length(ks)
    @test nbad_counts == 0
    @test nbad_S == 0
    @test nbad_D == 0
    @test nbad_mean == 0
    @info "AC1 B bumped to satisfy the edge condition (type, p_0, tau, B0, B); largest deviation of the realized mean from the closed form, in eps" count = length(bumped) first = first(bumped, 4) dev_cf
    # (iii) with the closed form itself: dyadic parameters, for which the orbit x_j = x_0 + j p_0 is exact in T (x_0 = 51/512, p_0 in {0, 51/1024, -51/1024}), tau in {1, 3, 5}, k up to 1000:
    # |mean z - c (x_0 + p_0 (n - 1) / 2)| <= 8 eps (n + G_first) / n + 2 eps |mean|
    nbad_dy = 0
    for T in SN_TYPES, p0r in (0 // 1, 51 // 1024, -51 // 1024), tau in (1, 3, 5)
        epsT = sn_eps(T)
        c, x0, p0 = exact(T, 3 // 10), exact(T, 51 // 512), exact(T, p0r)
        r = StretchingNumbers(; range = (T(-1), T(1)), bins = 1001, interval = sn_interval(tau))
        res = sn_run((r,), sn_T1(c), [SVector{2, T}(x0, p0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], 1000 * tau, tau; steps = [k * tau for k in ks])[1][1].snaps
        for (f, k) in zip(res, ks)
            n = k * tau
            mean_cf = setprecision(BigFloat, SN_BITS) do
                sn_big(c) * (sn_big(x0) + sn_big(p0) * (n - 1) / 2)
            end
            tol = 8 * epsT * (n + sn_T1_G(c, x0, p0, tau)) / n + 2 * epsT * Float64(abs(mean_cf))
            nbad_dy += !(sn_err(f.mean[1], mean_cf) <= tol && sn_err(f.mean[2], -mean_cf) <= tol)
        end
    end
    @test nbad_dy == 0
    # (iv) negative controls on plain arrays (each fails the predicate it names): a double that divides the log by the step count n = k tau instead of tau, a double that builds the spectrum
    # from the last interval only, a double that stores |z|. Float64, p_0 = 1/50, tau = 3, B = 1001 over [-1, 1] (and p_0 = 0 for |z|).
    T = Float64
    c, x0 = exact(T, 3 // 10), exact(T, 1 // 10)
    tau = 3
    for (p0r, fault) in ((1 // 50, :divn), (1 // 50, :window), (0 // 1, :abs))
        p0 = exact(T, p0r)
        zs1 = sn_T1_z1(c, x0, p0, tau, 1000)
        B = sn_pick_B(1001, T(-1), T(1), zs1, -zs1)[1]
        for k in (10, 100, 1000)
            good1, good2 = sn_hist_ref(zs1[1:k], T(-1), T(1), B), sn_hist_ref(-zs1[1:k], T(-1), T(1), B)
            if fault === :divn
                bad1 = sn_hist_ref([zs1[i] / i for i in 1:k], T(-1), T(1), B)
                k > 1 && @test bad1 != good1                                   # fails (i) for k > 1
            elseif fault === :window
                bad1, bad2 = sn_hist_fault(zs1[1:k], T(-1), T(1), B, :window), sn_hist_fault(-zs1[1:k], T(-1), T(1), B, :window)
                k <= 100 ? (@test sn_S(bad1, bad2) == sn_S(good1, good2)) : (@test sn_S(bad1, bad2) != sn_S(good1, good2))   # fails (ii) for k > 100
            else
                bad1, bad2 = sn_hist_fault(zs1[1:k], T(-1), T(1), B, :abs), sn_hist_fault(-zs1[1:k], T(-1), T(1), B, :abs)
                @test sn_S(bad1, bad2) == 0 && sn_S(good1, good2) == 2 * k^2      # D = 0 instead of sqrt(2 / h)
            end
        end
    end
end

@testset "AC2a constant one-step matrices" begin
    _ = (StretchingNumbers, spectral_distance)
    # SN8(a). Constant one-step matrices diag(l1, l2), l in {exp(1/2), exp(-1/2), exp(3/2)} (formed in the element type), supplied eigenvector columns (1, 0) and (0, 1), tau in {1, 5},
    # k in {1, 10, 1000}, B = 1001 over [-2, 2] (the map is not a standard map, so the caller supplies the range; B odd puts the three values 0.5, -0.5, 1.5 at 0.125, 0.375 and 0.125
    # of a cell from the edges): every stretching number is in the cell of ln|l|, so each count vector has the single entry k, tolerance 0; D = sqrt(2 / h) when the two cells differ
    # and D = 0 when they coincide, within 4 eps.
    ls = (1 // 2, -1 // 2, 3 // 2)
    B, ks = 1001, (1, 10, 1000)
    for T in SN_TYPES
        epsT = sn_eps(T)
        lo, hi = T(-2), T(2)
        for e1 in ls, e2 in ls, tau in (1, 5)
            l1, l2 = exp(exact(T, e1)), exp(exact(T, e2))
            lns = setprecision(BigFloat, SN_BITS) do
                (log(sn_big(l1)), log(sn_big(l2)))
            end
            zs = (fill(lns[1], 1000), fill(lns[2], 1000))
            @test sn_edge_ok(zs[1], lo, hi, B) && sn_edge_ok(zs[2], lo, hi, B)
            r = StretchingNumbers(; range = (lo, hi), bins = B, interval = sn_interval(tau))
            res = sn_run((r,), sn_Diag(l1, l2), [SVector{2, T}(0, 0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], 1000 * tau, tau; steps = [k * tau for k in ks])[1][1].snaps
            for (f, k) in zip(res, ks)
                cok, Sok, Dok = sn_check_snapshot(f, zs[1], zs[2], k, lo, hi, B, epsT)
                @test cok && Sok && Dok
                @test count(!iszero, sn_cnt(f, 1)) == 1 && maximum(sn_cnt(f, 1)) == k && count(!iszero, sn_cnt(f, 2)) == 1 && maximum(sn_cnt(f, 2)) == k
                same = e1 == e2
                @test same ? f.D == 0 : sn_err(f.D, sqrt(2 / ((sn_big(hi) - sn_big(lo)) / B))) <= 4 * epsT * Float64(sqrt(2 / ((sn_big(hi) - sn_big(lo)) / B)))
            end
        end
    end
end

@testset "AC2b rotation block gives zero stretching numbers" begin
    _ = (StretchingNumbers, spectral_distance)
    # SN8(a) with l = 1. The Y5 rotation block (rotation_block(T, 0.3)), seeded pairs of columns (ids 1..10), tau = 1, N = 1000, B = 1001 over [-1, 1]: every z is within 2 eps of 0 (read from
    # the differences of the compensated totals at every handoff), all counts of both columns sit in the central cell (B odd), D = 0, tolerance 0.
    seed = sn_seed(2, 1)
    B = 1001
    for T in SN_TYPES
        epsT = sn_eps(T)
        sys = LinearSymplecticMap(rotation_block(T, 0.3))
        r = StretchingNumbers(; range = (T(-1), T(1)), bins = B, interval = 100)
        res = sn_run((r,), sys, [SVector{2, T}(0.3, 0.7) for _ in 1:10], [sn_frame(seed, j, T) for j in 1:10], 1000, 1; steps = [1000], light = true)
        for j in 1:10
            f = res[j][1].snaps[1]
            lights = res[j][1].lights
            @test length(lights) == 1000
            for col in 1:2
                zs = sn_z_from_lights(lights, col, 1)
                @test maximum(abs.(zs)) <= 2 * epsT * (1 + 1e-9)           # the decode (a difference of two compensated totals) adds about 1e-31
                cnt = sn_cnt(f, col)
                @test cnt[(B + 1) ÷ 2 + 1] == 1000 && sum(cnt) == 1000
            end
            @test f.D == 0
        end
    end
end

@testset "AC2c hyperbolic closed form (Y5)" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # SN8(b). Y5 A = [[2, 1], [1, 1]] (mu = (3 + sqrt 5) / 2), default range of SN12, supplied pairs with |c_1| in {0.05, 0.6, 0.99} (column 1 = c_1 e_1 + c_2 e_2, column 2 the
    # orthogonal vector) and the seeded pairs of ids 1..3, tau in {1, 2, 5}, k in {1, 10, 100, 1000}, B0 in {101, 1000}: the counts equal those of the closed form in every cell, tolerance 0
    # under the edge condition (at least 1e-9 from the interior edges; at the two end edges the closed-form values sit about 8 eps inside, so they must be at least 2 eps from them), D
    # within 4 eps relative of the closed form; the Float64 and Double64 results agree cell by cell.
    seed = sn_seed(2, 3)
    ks = (1, 10, 100, 1000)
    wsup(T, c1) = begin
        t = setprecision(BigFloat, SN_BITS) do
            mu = (3 + sqrt(BigFloat(5))) / 2
            tt = mu - 2
            nrm = sqrt(1 + tt^2)
            c2 = sqrt(1 - BigFloat(c1)^2)
            e1, e2 = (1 / nrm, tt / nrm), (-tt / nrm, 1 / nrm)
            w1 = (BigFloat(c1) * e1[1] + c2 * e2[1], BigFloat(c1) * e1[2] + c2 * e2[2])
            w2 = (-c2 * e1[1] + BigFloat(c1) * e2[1], -c2 * e1[2] + BigFloat(c1) * e2[2])
            (w1, w2)
        end
        SMatrix{2, 2, T}(T(t[1][1]), T(t[1][2]), T(t[2][1]), T(t[2][2]))
    end
    cases = Any[]
    for c1 in (0.05, 0.6, 0.99)
        push!(cases, (label = "supplied $c1", W0 = Dict(T => wsup(T, c1) for T in SN_TYPES)))
    end
    for j in 1:3
        push!(cases, (label = "seeded $j", W0 = Dict(T => sn_frame(seed, j, T) for T in SN_TYPES)))
    end
    for case in cases, tau in (1, 2, 5), B0 in (101, 1000)
        sys = Dict(T => LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1)) for T in SN_TYPES)
        rng = Dict(T => stretching_range(sys[T]) for T in SN_TYPES)
        zs = Dict(T => (sn_y5_z(case.W0[T][:, 1], tau, 1000), sn_y5_z(case.W0[T][:, 2], tau, 1000)) for T in SN_TYPES)
        # one B for both types, so that the cells can be compared across types
        B = nothing
        for Bc in B0:(B0 + 60)
            if all(T -> all(z -> sn_edge_ok(z, rng[T]..., Bc; endslack = 2 * sn_eps(T) * max(1.0, Float64(abs(rng[T][2])))), zs[T]), SN_TYPES)
                B = Bc
                break
            end
        end
        B === nothing && error("AC2c: no B satisfies the edge condition")
        cnts = Dict{Any, Any}()
        for T in SN_TYPES
            epsT = sn_eps(T)
            lo, hi = rng[T]
            r = StretchingNumbers(; range = (lo, hi), bins = B, interval = sn_interval(tau))
            res = sn_run((r,), sys[T], [SVector{2, T}(0, 0)], [case.W0[T]], 1000 * tau, tau; steps = [k * tau for k in ks])[1][1].snaps
            for (f, k) in zip(res, ks)
                cok, Sok, Dok = sn_check_snapshot(f, zs[T][1], zs[T][2], k, lo, hi, B, epsT)
                @test cok && Sok && Dok
                @test f.below == [0, 0] && f.above == [0, 0]
                cnts[(T, k)] = (sn_cnt(f, 1), sn_cnt(f, 2))
            end
        end
        @test all(k -> cnts[(Float64, k)] == cnts[(Double64, k)], ks)
    end
end

@testset "AC3 shear (K_std = 0)" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # SN8(c), SN10, F5. The K_std = 0 standard map (wrapped), supplied starts (a, b) in {(0.6, 0.8), (-0.6, 0.8), (1, 0), (0, 1)} for column 1 and the orthogonal vector (-b, a) for column 2,
    # tau in {1, 3}, k in {10, 1000, 10^4} (10^5 with SN_FULL, D38: the BigFloat closed form of 10^5 samples dominates the runtime), B = 1001 over the default range: (i) counts equal
    # the closed-form counts, tolerance 0 under the edge condition; for b = 0 every z is exactly 0 (sum, compensation, mean and second sum are exactly 0) and all counts sit in the
    # central cell; (ii) D_k <= sqrt 2 M_k / (k sqrt h) (1 + 4 eps) at every evaluation handoff, M_k the number of handoffs at which the two columns' cells differ; (iii) reported:
    # tau k z_k at the largest k and D k.
    kmax = SN_FULL ? 100_000 : 10_000
    ks = (10, 1000, kmax)
    starts = ((0.6, 0.8), (-0.6, 0.8), (1.0, 0.0), (0.0, 1.0))
    for T in SN_TYPES
        epsT = sn_eps(T)
        sys = WrappedStandardMap(exact(T, 0))
        lo, hi = stretching_range(sys)
        x0 = sn_astate(sys, SVector{2, T}(0.3, 0.7))
        for (a, b) in starts, tau in (1, 3)
            at, bt = T(a), T(b)
            W0 = SMatrix{2, 2, T}(at, bt, -bt, at)
            zs1 = sn_shear_z(at, bt, tau, kmax)
            zs2 = sn_shear_z(-bt, at, tau, kmax)
            B, bump = sn_pick_B(1001, lo, hi, zs1, zs2; endslack = 2 * epsT, tol = 1e-9)
            interval = sn_interval(tau)
            evals = sn_eval_steps(kmax * tau, interval)
            r = StretchingNumbers(; range = (lo, hi), bins = B, interval = interval)
            res = sn_run((r,), sys, [x0], [W0], kmax * tau, tau; steps = sort(unique(vcat(evals, [k * tau for k in ks]))))[1][1].snaps
            byn = Dict(f.n => f for f in res)
            for k in ks
                f = byn[k * tau]
                cok, Sok, Dok = sn_check_snapshot(f, zs1, zs2, k, lo, hi, B, epsT)
                @test cok && Sok && Dok
                if b == 0
                    @test iszero(f.sum[1]) && iszero(f.comp[1]) && iszero(f.sum2[1]) && iszero(f.mean[1])      # column 1 stays (1, 0): every z is exactly 0
                    @test sn_cnt(f, 1)[(B + 1) ÷ 2 + 1] == k
                end
                if a == 0
                    @test iszero(f.sum[2]) && iszero(f.comp[2]) && sn_cnt(f, 2)[(B + 1) ÷ 2 + 1] == k        # column 2 is (-1, 0)
                end
            end
            # (ii) the alignment bound at every evaluation handoff, with M_k from the reference cells
            c1, c2 = sn_cells_ref(zs1, lo, hi, B), sn_cells_ref(zs2, lo, hi, B)
            Mk = cumsum(c1 .!= c2)
            h = setprecision(BigFloat, SN_BITS) do
                (sn_big(hi) - sn_big(lo)) / B
            end
            nbad = 0
            for n in evals
                k = n ÷ tau
                f = byn[n]
                bound = sqrt(BigFloat(2)) * Mk[k] / (k * sqrt(h)) * (1 + 4 * BigFloat(epsT))
                nbad += !(sn_big(f.D) <= bound)
            end
            @test nbad == 0
            # (iii) reported
            f = byn[kmax * tau]
            zk = Float64(zs1[kmax])
            @info "AC3 shear (iii)" T a b tau B bumped = bump tau_k_zk = tau * kmax * zk D = Float64(f.D) D_times_k = Float64(f.D) * kmax
        end
    end
end
