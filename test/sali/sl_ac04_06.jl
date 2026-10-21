# AC4 to AC6: analytic cases (Y5 hyperbolic map, hyperbolic fixed point of the standard map, rotation and elliptic point).
# The reference is the BigFloat(256) product Y_n W0 (sl_ref), which is exact for these constant tangents; the closed forms of A11 and A12 are
# checked against it where the spec states them (asymptote or full expansion), so the helper itself is tested.

# Eigen data of A = [[2, 1], [1, 1]] in BigFloat (call inside setprecision): mu, orthonormal v_+ and v_- (L1 E10).
function sl_y5_eigen()
    mu = (3 + sqrt(BigFloat(5))) / 2
    phi = (1 + sqrt(BigFloat(5))) / 2
    nv = sqrt(phi^2 + 1)
    return (mu = mu, vp = (phi / nv, 1 / nv), vm = (-1 / nv, phi / nv))
end

@testset "AC4 linear hyperbolic map" begin
    _ = (SALI, sali_floor)
    # A = [[2, 1], [1, 1]] (Y5), start at the origin (a fixed point), ids 1..100 (Double64: 30, D38), tau in {1, 5}. Ids with
    # |c_1^(1) c_1^(2)| < 0.01 are excluded (at most 10 percent). (a) resolved n: |SALI - SALI_cf| <= B_SALI (n in {1, 2, 4, 8} for tau = 1; the handoffs
    # {5, 10} for tau = 5); (b) log form at n in {20, 100, 200}: |L - L_cf| <= B_L; (c) rate: |(-L) - (2 n ln mu - ln|D| + ln|c_1^(1) c_1^(2)|)| <= B_L + mu^(-4n)
    # at n in {20, 100}. L_cf is the A11 form |D| / (norm w1 norm w2) in BigFloat.
    seed = sl_seed(4)
    keep = [1, 2, 4, 5, 8, 10, 20, 100, 200]
    for T in SL_TYPES
        epsT, meas = sl_eps(T), sl_meas(T)
        nid = T === Double64 ? 30 : 100
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        W0 = [sl_frame(seed, j, T) for j in 1:nid]
        Ms = lb_tangents(sys, x0[1], 200)                       # the tangent is A everywhere
        refs = [sl_ref(sl_mats_const(2, 1, 1, 1), W0[j], 200; keep = keep) for j in 1:nid]
        psis = [sl_psis(Ms, W0[j], keep) for j in 1:nid]
        lam(j, n) = sl_Lambda(refs[j].lnsali, n)
        # c-coefficients and D of A11, per id, and the exclusion rule
        c = setprecision(BigFloat, SL_BITS) do
            e = sl_y5_eigen()
            map(1:nid) do j
                w1 = (sl_big(W0[j][1, 1]), sl_big(W0[j][2, 1]))
                w2 = (sl_big(W0[j][1, 2]), sl_big(W0[j][2, 2]))
                c1 = (w1[1] * e.vp[1] + w1[2] * e.vp[2], w2[1] * e.vp[1] + w2[2] * e.vp[2])
                c2 = (w1[1] * e.vm[1] + w1[2] * e.vm[2], w2[1] * e.vm[1] + w2[2] * e.vm[2])
                (c1 = c1, c2 = c2, D = c1[1] * c2[2] - c2[1] * c1[2])
            end
        end
        excluded = count(cj -> abs(cj.c1[1] * cj.c1[2]) < 0.01, c)
        @test excluded <= nid ÷ 10
        # the helper against the closed form of A11 (full expansion of the norms), n = 8
        setprecision(BigFloat, SL_BITS) do
            e = sl_y5_eigen()
            for j in 1:5
                n2(i) = e.mu^16 * c[j].c1[i]^2 + e.mu^(-16) * c[j].c2[i]^2
                Lcf = log(abs(c[j].D)) - log(n2(1)) / 2 - log(n2(2)) / 2
                @test abs(Lcf - sl_at(refs[j], 8).L) < 1e-60
            end
        end
        for tau in (1, 5)
            run = sl_run(SALI(), sys, x0, W0, 200, tau; steps = filter(n -> n % tau == 0, keep))
            nresolved = 0
            for j in 1:nid
                cj = c[j]
                abs(cj.c1[1] * cj.c1[2]) < 0.01 && continue
                for s in run[j].snaps
                    n = s.n
                    at = sl_at(refs[j], n)
                    (p1, p2) = psis[j][findfirst(==(n), keep)]
                    saliref = Float64(at.S)
                    if n in (tau == 1 ? (1, 2, 4, 8) : (5, 10)) && saliref >= 10 * meas              # (a) resolved only
                        nresolved += 1
                        @test sl_err(s.sali, at.S) <= sl_BSALI(epsT, lam(j, n), saliref, p1, p2)
                    end
                    if n in (20, 100, 200)                                                             # (b)
                        @test sl_err(s.L, at.L) <= sl_BL(epsT, p1, p2)
                    end
                    if n in (20, 100)                                                                  # (c)
                        asym = setprecision(BigFloat, SL_BITS) do
                            mu = sl_y5_eigen().mu
                            2 * n * log(mu) - log(abs(cj.D)) + log(abs(cj.c1[1] * cj.c1[2]))
                        end
                        mu4 = setprecision(BigFloat, SL_BITS) do
                            Float64(sl_y5_eigen().mu^(-4n))
                        end
                        @test sl_err(-s.L, asym) <= sl_BL(epsT, p1, p2) + mu4
                    end
                end
            end
            @test nresolved > 0
        end
    end
end

@testset "AC5 hyperbolic fixed point" begin
    _ = (SALI, sali_floor)
    # Fixed point (0, 0) of the standard map, K in {1, 2, 6}, ids 1..100 (Double64: 30, D38), tau in {1, 5} (A12, L1 E12). (a) resolved n
    # (>= 10 B_meas): SALI within B_SALI of the BigFloat product, n in {1, 2, 4, 8} (tau = 1) or the handoffs {5, 10} (tau = 5); (b) n in {20, 100, 200}:
    # |(-L) - [2 n ln mu_+ - ln|det W0| + sum_i ln(|a_i| norm v_u)]| <= B_L + mu_+^(-2n). The BigFloat product is first checked against the
    # eigen-expansion of A12 (cross term kept).
    seed = sl_seed(5)
    keep = [1, 2, 4, 5, 8, 10, 20, 100, 200]
    for T in SL_TYPES, K in (1, 2, 6)
        epsT, meas = sl_eps(T), sl_meas(T)
        nid = T === Double64 ? 30 : 100
        sys = StandardMap(exact(T, K))
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        W0 = [sl_frame(seed, j, T) for j in 1:nid]
        Ms = lb_tangents(sys, x0[1], 200)
        refs = [sl_ref(sl_mats_const(1 + K, 1, K, 1), W0[j], 200; keep = keep) for j in 1:nid]
        psis = [sl_psis(Ms, W0[j], keep) for j in 1:nid]
        # A12 data per id (BigFloat): mu_+, norm v_u, a_i, b_i, ln|det W0|
        d = setprecision(BigFloat, SL_BITS) do
            Kb = BigFloat(K)
            s = sqrt(Kb * (Kb + 4))
            mu = (2 + Kb + s) / 2
            vu, vs = (one(Kb), (s - Kb) / 2), (one(Kb), -(Kb + s) / 2)
            ai(w) = (w[2] + w[1] * (Kb + s) / 2) / s
            per = map(1:nid) do j
                w = ((sl_big(W0[j][1, 1]), sl_big(W0[j][2, 1])), (sl_big(W0[j][1, 2]), sl_big(W0[j][2, 2])))
                a = ntuple(i -> ai(w[i]), 2)
                b = ntuple(i -> w[i][1] - a[i], 2)
                (a = a, b = b, ldet = log(abs(w[1][1] * w[2][2] - w[1][2] * w[2][1])))
            end
            (mu = mu, vu = vu, vs = vs, per = per)
        end
        setprecision(BigFloat, SL_BITS) do
            vv = d.vu[1] * d.vs[1] + d.vu[2] * d.vs[2]
            for j in 1:nid, n in (8,)
                n2(i) = d.mu^(2n) * d.per[j].a[i]^2 * (d.vu[1]^2 + d.vu[2]^2) + 2 * d.per[j].a[i] * d.per[j].b[i] * vv +
                        d.mu^(-2n) * d.per[j].b[i]^2 * (d.vs[1]^2 + d.vs[2]^2)
                Lcf = d.per[j].ldet - log(n2(1)) / 2 - log(n2(2)) / 2
                @test abs(Lcf - sl_at(refs[j], n).L) < 1e-60
            end
        end
        for tau in (1, 5)
            run = sl_run(SALI(), sys, x0, W0, 200, tau; steps = filter(n -> n % tau == 0, keep))
            nresolved = 0
            for j in 1:nid, s in run[j].snaps
                n = s.n
                at = sl_at(refs[j], n)
                (p1, p2) = psis[j][findfirst(==(n), keep)]
                saliref = Float64(at.S)
                if n in (tau == 1 ? (1, 2, 4, 8) : (5, 10)) && saliref >= 10 * meas
                    nresolved += 1
                    @test sl_err(s.sali, at.S) <= sl_BSALI(epsT, sl_Lambda(refs[j].lnsali, n), saliref, p1, p2)
                end
                if n in (20, 100, 200)
                    asym, tail = setprecision(BigFloat, SL_BITS) do
                        nvu = sqrt(d.vu[1]^2 + d.vu[2]^2)
                        # tail: first-order coefficient of the neglected cross term, sum_i |b_i/a_i| |v_u . v_s| / norm(v_u)^2, times mu^-2n (the spec writes
                        # mu^-2n; the constant is made explicit, it exceeds 1 for some frames and then matters in Double64)
                        coef = sum(i -> abs(d.per[j].b[i] / d.per[j].a[i]), 1:2) * abs(d.vu[1] * d.vs[1] + d.vu[2] * d.vs[2]) / nvu^2
                        (2 * n * log(d.mu) - d.per[j].ldet + sum(i -> log(abs(d.per[j].a[i]) * nvu), 1:2), Float64(1.01 * (coef + 1) * d.mu^(-2n)))
                    end
                    @test sl_err(-s.L, asym) <= sl_BL(epsT, p1, p2) + tail
                end
            end
            @test nresolved > 0
        end
    end
end

@testset "AC6 elliptic and rotation cases" begin
    _ = (SALI, sali_floor)
    # (a) the order-4 rotation [[0, 1], [-1, 0]], seeded frame, n in {12, 120, 1200, 9996}, tau in {1, 3}: |SALI(n) - SALI(0)| <= B_SALI(n) (angles are
    # preserved). (b) K = -2 at (0, 0) (A13), ids 1..100 (Double64: 30, D38), tau in {1, 2}: even n |SALI(n) - SALI(0)| <= B_SALI, odd n within B_SALI of the
    # BigFloat SALI_cf (the exact product), and SALI(n) >= 0.146 - B_SALI(n). n in 1..100 and every 10th up to 1000 (D38; SL_FULL: every n).
    seed = sl_seed(6)
    for T in SL_TYPES
        epsT = sl_eps(T)
        rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
        Mrot = SMatrix{2, 2, Float64}(0, -1, 1, 0)
        W0 = [sl_frame(seed, 1, T)]
        S0 = sl_sali0(W0[1])
        ns = (12, 120, 1200, 9996)
        Ms = fill(Mrot, 9999)
        ref = sl_ref(sl_mats_const(0, 1, -1, 0), W0[1], 9999; keep = collect(ns))
        for tau in (1, 3)
            run = sl_run(SALI(), rot, [SVector{2, T}(0, 0)], W0, tau == 1 ? 10_000 : 9_999, tau; steps = collect(ns))
            @test length(run[1].snaps) == length(ns)
            for s in run[1].snaps
                n = s.n
                at = sl_at(ref, n)
                p1, p2 = lb_Psi(view(Ms, 1:n), W0[1][:, 1]), lb_Psi(view(Ms, 1:n), W0[1][:, 2])
                @test sl_err(s.sali, S0) <= sl_BSALI(epsT, sl_Lambda(ref.lnsali, n), Float64(S0), p1, p2)
                @test abs(at.S - S0) < 1e-70                                    # the reference itself keeps the angle
            end
        end
        sys = StandardMap(exact(T, -2))
        nid = T === Double64 ? 30 : 100
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        W0 = [sl_frame(seed, j, T) for j in 1:nid]
        nsamp = SL_FULL ? collect(1:1000) : sort(unique(vcat(1:100, 110:10:1000)))
        Ms = lb_tangents(sys, x0[1], 1000)
        refs = [sl_ref(sl_mats_const(-1, 1, -2, 1), W0[j], 1000; keep = nsamp) for j in 1:nid]
        psis = [sl_psis(Ms, W0[j], nsamp) for j in 1:nid]
        S0s = [sl_sali0(W0[j]) for j in 1:nid]
        for tau in (1, 2)
            run = sl_run(SALI(), sys, x0, W0, 1000, tau; steps = filter(n -> n % tau == 0, nsamp))
            worst_ref, worst_even, min_gap = 0.0, 0.0, Inf
            for j in 1:nid, s in run[j].snaps
                n = s.n
                at = sl_at(refs[j], n)
                (p1, p2) = psis[j][findfirst(==(n), nsamp)]
                saliref = Float64(at.S)
                B = sl_BSALI(epsT, sl_Lambda(refs[j].lnsali, n), saliref, p1, p2)
                worst_ref = max(worst_ref, sl_err(s.sali, at.S) / B)                                          # odd and even n against the exact product
                iseven(n) && (worst_even = max(worst_even, sl_err(s.sali, S0s[j]) / B))              # even n: SALI(0)
                min_gap = min(min_gap, Float64(s.sali) - (0.146 - B))
            end
            @test worst_ref <= 1
            @test worst_even <= 1
            @test min_gap >= 0
        end
    end
end
