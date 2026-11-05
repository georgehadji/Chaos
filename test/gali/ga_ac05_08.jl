# AC5 to AC8: determinant identity, floor and precision ladder of the evaluation, health and the frame check, reducer structure on the stream.
# Each testset starts by naming the new symbols it uses; the negative controls (test-side doubles) are separate testsets and pass now.

# ---- AC5: determinant identity (B5, I3) ----
# Worst ratio |G_d - abs(det W0) exp(L)| / (B_G^d(n) + 8 eps G_ref sum_i Psi_n^(i)) over the snapshots of a BigFloat run with K = d, against the same-orbit reference b of
# ga_bundle_core (the spec states the bound 2 B_meas / G_ref for the logarithm; that holds for the evaluation alone, not after the transport of n handoffs, whose error
# is the Lambda_n term of B10, so the bound of AC9(c) is used: it is the same identity, ln G_d = ln abs(det W0) - sum S_i, as an absolute bound on G_d).
function ga_identity_ratio(b, snaps, W0::AbstractMatrix)
    d, epsT = b.d, ga_eps(BigFloat)
    worst = 0.0
    for s in snaps
        ic = findfirst(==(s.n), b.cps)
        psi, lam, gd = b.psi[ic], b.lam[ic], Float64(b.ref.G[s.n][d - 1])
        B = ga_BG(epsT, d, lam[d], gd, psi) + 8 * epsT * gd * sum(psi)
        diff = setprecision(BigFloat, 2 * GA_BITS) do
            Float64(abs(s.gali[end] - abs(det(Matrix{BigFloat}(W0))) * exp(s.L)))
        end
        worst = max(worst, diff / B)
    end
    return worst
end

@testset "AC5a determinant identity, d = 2 standard map" begin
    _ = (GALI,)
    # d = 2, K = d = 2, standard map K in {0.5, 2 pi, 10}, BigFloat(256) module run (T = BigFloat), seeded and non-orthonormal unit-column frames, tau in {1, 5}:
    # G_2 = abs(det W0) exp(L), L = -(S_1 + S_2), within B_G^2(n) + 8 eps G_ref sum_i Psi (the ratio is reported). D38: 3 starts per K (spec 20), checkpoints
    # up to n = 40 (chaotic K = 2 pi, 10) and 100 (K = 0.5; spec 1000), tau in {1, 5} (the checkpoints are multiples of 10).
    rng = StableRNG(ga_seed(5))
    worst = 0.0
    setprecision(BigFloat, GA_BITS) do
        for (label, n) in ((:half, ga_size(100, 1000)), (:two_pi, 40), (:ten, 40)), j in 1:ga_size(3, 20), tau in (1, 5)
            Kb = label === :half ? BigFloat(0.5) : label === :ten ? BigFloat(10) : 2 * BigFloat(pi)
            sys = StandardMap(Kb)
            x0 = SVector{2, BigFloat}(BigFloat(2 * Float64(pi) * rand(rng)), BigFloat(2 * Float64(pi) * rand(rng)))
            ang = BigFloat(0.7)
            cps = ga_checkpoints(n)
            for W0 in (ga_frame(ga_seed(5, 1), j, BigFloat, 2, 2), SMatrix{2, 2, BigFloat}(1, 0, cos(ang), sin(ang)))
                r = ga_run1(GALI(), sys, x0, W0, n, tau; steps = cps, orbit = j)
                b = ga_bundle_core(sys, nothing, x0, W0, n, cps; bits = 2 * GA_BITS, sysB = sys)
                worst = max(worst, ga_identity_ratio(b, r.snaps, W0))
                @test length(r.snaps) == length(cps)
            end
        end
    end
    @test worst <= 1
    @info "AC5a reported: worst ratio" worst
end

@testset "AC5b determinant identity, products of random symplectic matrices" begin
    _ = (GALI,)
    # d in {4, 6}, K = d, products of 8 random symplectic BigFloat matrices from elementary factors (shears with symmetric blocks, block-diagonal [A, A^-T],
    # rotations), applied by a time-dependent test map; random orthonormal and non-orthonormal W0 (unit columns): the identity of AC5a at n = 8 within its bound.
    # D38: 6 products per d (spec: many).
    rng = StableRNG(ga_seed(5, 2))
    worst = 0.0
    setprecision(BigFloat, GA_BITS) do
        for d in (4, 6), j in 1:ga_size(6, 50)
            mats = [ga_symp_factor(rng, d ÷ 2) for _ in 1:8]
            sys = ga_Seq(d, mats)
            x0 = SVector{d + 1, BigFloat}(ntuple(_ -> zero(BigFloat), d + 1))
            Wr = ga_unitcols(BigFloat[BigFloat(randn(rng)) for _ in 1:d, _ in 1:d])
            for W0 in (ga_frame(ga_seed(5, 3), j, BigFloat, d, d), ga_sm(Wr))
                r = ga_run1(GALI(), sys, x0, W0, 8, 1; steps = [8], orbit = j)
                b = ga_bundle_core(sys, nothing, x0, W0, 8, [8]; bits = 2 * GA_BITS, sysB = sys)
                worst = max(worst, ga_identity_ratio(b, r.snaps, W0))
            end
        end
    end
    @test worst <= 1
end

@testset "AC5c determinant identity on Y3 and Y4 against the 256-bit product" begin
    _ = (GALI,)
    # Y3 orbits C and D (GA_Y3A) and Y4 chains N in {2, 3, 5} (K = 3, eps = 0.2, seeded chaotic starts), K = d, BigFloat(256) module run: the identity of AC5a at
    # n in {5, 18, 100} within its bound (reference: the 512-bit product of the tangents at the widened realised points), and G_d against that product within 1e-60.
    # D38: n = 100 only for the Y3 cases and N = 2, N = 3 and 5 only to n = 18 (spec: n in {5, 18, 100} for all).
    rng = StableRNG(ga_seed(5, 4))
    cases = [(GA_Y3A, GA_ORBIT_C, 100), (GA_Y3A, GA_ORBIT_D, 100), (ga_y4(2, 3.0, 0.2), ga_draw_chaotic(rng, ga_y4(2, 3.0, 0.2)), 100),
             (ga_y4(3, 3.0, 0.2), ga_draw_chaotic(rng, ga_y4(3, 3.0, 0.2)), ga_size(18, 100)), (ga_y4(5, 3.0, 0.2), ga_draw_chaotic(rng, ga_y4(5, 3.0, 0.2)), ga_size(18, 100))]
    worst, dref = 0.0, 0.0
    for (j, (c, z, N)) in enumerate(cases)
        steps = filter(n -> n <= N, [5, 18, 100])
        d = ga_d(c)
        b = ga_bundle(c, BigFloat, z, (ga_seed(5, 5), j), N; cps = steps)
        r = setprecision(BigFloat, GA_BITS) do
            ga_run1(GALI(), ga_sys(BigFloat, c), b.x0, b.W0, N, 1; steps = steps, orbit = j)
        end
        worst = max(worst, ga_identity_ratio(b, r.snaps, b.W0))
        for (s, n) in zip(r.snaps, steps)
            @test s.n == n
            dref = max(dref, setprecision(() -> Float64(abs(s.gali[end] - b.ref.G[n][end])), BigFloat, 2 * GA_BITS))
        end
    end
    @test worst <= 1
    @test dref <= 1e-60
    @info "AC5c reported" worst dref
end

@testset "AC5d negative control: the identity does not hold for k < d" begin
    # The same formula for k = 2 < d = 4, exp(ln vol(W0 columns) - S_1 - S_2), differs from G_2 by a factor of at least 10 for at least 90 percent of the
    # products of AC5b (test-side, exact BigFloat products with shear amplitude 6, so that the growth is large; no new symbol).
    rng = StableRNG(ga_seed(5, 6))
    far, total = 0, 60
    setprecision(BigFloat, GA_BITS) do
        for _ in 1:total
            Y = Matrix{BigFloat}(I, 4, 4)
            for _ in 1:8
                Y = ga_symp_factor(rng, 2; amp = 6.0) * Y
            end
            W0 = ga_unitcols(BigFloat[BigFloat(randn(rng)) for _ in 1:4, _ in 1:2])
            Z = Y * W0
            norms = [sqrt(sum(abs2, Z[:, j])) for j in 1:2]
            G2 = ga_vols(Z ./ reshape(norms, 1, 2))[2]
            vol0 = ga_vols(W0)[2]                           # vol of the two unit columns of W0
            f = vol0 / (norms[1] * norms[2])
            far += abs(log(f / G2)) >= log(10) ? 1 : 0
        end
    end
    @test far >= 0.9 * total
end

# ---- AC6: floor, pivot refinement and precision ladder (B9, F5) ----
# One near-dependent set of AC6 as an exact BigFloat(512) matrix rotated by Q, rounded to T; returns (W, U) with U the unit columns of the rounded W at 512 bits.
function ga_neardep_set(::Type{T}, Q, k::Int, m::Int) where {T}
    return setprecision(BigFloat, 512) do
        Ub = Q * ga_neardep(k, BigFloat(10)^(-m))
        W = ga_sm([ga_to(T, Ub[i, j]) for i in 1:4, j in 1:k])
        (W, ga_unitcols(BigFloat[ga_big(W[i, j]) for i in 1:4, j in 1:k]))
    end
end

@testset "AC6 floor, pivot refinement and precision ladder" begin
    _ = (gali_values, gali_floor)
    # Near-dependent unit sets in d = 4 rotated by a random orthogonal map (D38: 4 maps per type and order, spec 100): k = 2 with theta = 1e-m, k = 3 with two
    # columns at angle 1e-m and a third orthogonal, k = 4 a nearly coplanar quadruple, m = 1..40. (a) |G_T - exact| <= 2 B_meas(T, k) for all m, and a computed
    # value of at most 2 B_meas where exact < B_meas; (b) m_floor(T, k), the first m with relative error above 50 percent (minimum over the maps), is
    # at least log10(1/(2 B_meas)) - 1; (c) |error| <= 2 eps sum_i V_i + 8 eps^2 for at least 99 percent of the sets; (e) BigFloat(256) against a 512-bit exact
    # value, m = 5, 10, ..., 70 (D38: 4 maps per k, spec 200 sets): |error| <= 2 B_meas(BigFloat, k). The Gram-determinant control (d) is AC6d.
    rng = StableRNG(ga_seed(6))
    for T in GA_TYPES, k in 2:4
        epsT = ga_eps(T)
        bad_a, bad_floor, sharp_bad, total = 0, 0, 0, 0
        mfloor = Inf
        for _ in 1:ga_size(4, 100)
            Q = setprecision(() -> ga_randorth(rng, 4), BigFloat, 512)
            first_bad = Inf
            for m in 1:40
                W, U = ga_neardep_set(T, Q, k, m)
                g = gali_values(W)[end]
                ex, err, Vsum = setprecision(BigFloat, 512) do
                    ex = ga_vols(U)[end]
                    Vs = sum(ga_vols(U[:, setdiff(1:k, i)])[end] for i in 1:k)
                    (ex, ga_err(g, ex; bits = 512), Float64(Vs))
                end
                total += 1
                bad_a += err <= 2 * ga_meas(T, k) ? 0 : 1
                Float64(ex) < ga_meas(T, k) && (bad_floor += g <= 2 * ga_meas(T, k) ? 0 : 1)
                first_bad == Inf && err > 0.5 * Float64(ex) && (first_bad = m)
                sharp_bad += err <= 2 * epsT * Vsum + 8 * epsT^2 ? 0 : 1
            end
            mfloor = min(mfloor, first_bad)
        end
        @test bad_a == 0
        @test bad_floor == 0
        @test mfloor >= log10(1 / (2 * ga_meas(T, k))) - 1
        @test sharp_bad <= 0.01 * total
        @info "AC6 reported" T k mfloor sharp_bad total
    end
    for k in 2:4
        bad = 0
        for _ in 1:ga_size(4, 200), m in 5:5:70
            Q = setprecision(() -> ga_randorth(rng, 4), BigFloat, 512)
            W, U = ga_neardep_set(BigFloat, Q, k, m)
            g = gali_values(W)[end]
            ex = ga_exact(W; bits = 512)[end]
            bad += ga_err(g, ex; bits = 512) <= 2 * ga_meas(BigFloat, k) ? 0 : 1
        end
        @test bad == 0
    end
end

@testset "AC6d negative control: Gram-determinant double" begin
    # A double that evaluates through the Gram determinant fails the bound of AC6(a) in Float64 for the near-dependent k = 2 sets at m = 9..13 (and so the bound of
    # AC1(b) at the near-dependent sets); the unit-set value is sin(theta), the Gram value is 0 or about 1.5e-8. The spec says "at m >= 9": beyond m = 15 the
    # Gram value is 0 and the error 1e-m is below the bound, so the control is stated for m = 9..13.
    rng = StableRNG(ga_seed(6, 1))
    fail, total = 0, 0
    for _ in 1:30
        Q = setprecision(() -> ga_randorth(rng, 4), BigFloat, 512)
        for m in 9:13
            W, U = ga_neardep_set(Float64, Q, 2, m)
            ex = setprecision(() -> ga_vols(U)[end], BigFloat, 512)
            fail += ga_err(ga_ctl_gram(W, 2), ex; bits = 512) > 2 * ga_meas(Float64, 2) ? 1 : 0
            total += 1
        end
    end
    @test fail >= 0.5 * total
    @info "AC6d reported: failing fraction" fail total
end

# ---- AC7: health, rejection and the frame check ----
function ga_frame_with(::Type{T}, K::Int, g) where {T}
    # d = 4, columns e_1 .. e_{K-1} and a last column sqrt(1 - g^2) e_1 + g e_K, so that G_K(0) = g (g = 0: the last column equals column 1)
    M = zeros(T, 4, K)
    for j in 1:(K - 1)
        M[j, j] = one(T)
    end
    M[1, K] = sqrt(1 - T(g)^2)
    M[K, K] = T(g)
    return ga_sm(M)
end

@testset "AC7a rejection of K outside 2..d and of k outside 2..K" begin
    _ = (GALI, gali_value)
    # K < 2 and K > d are rejected before any work by an argument error (ga_Boom throws if the map is called); K = 2..d is accepted; gali_value rejects k outside 2..K.
    for T in GA_TYPES
        @test_throws ArgumentError reducer_init(GALI(), T, Val(4), Val(1))
        @test_throws ArgumentError reducer_init(GALI(), T, Val(4), Val(5))
        @test_throws ArgumentError reducer_init(GALI(), T, Val(2), Val(3))
        for K in 2:4
            @test reducer_init(GALI(), T, Val(4), Val(K)) !== nothing
        end
    end
    seed = ga_seed(7)
    @test_throws ArgumentError propagate(ga_Boom(4), [SVector(0.0, 0.0, 0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (GALI(),))
    W = ga_frame(seed, 1, Float64, 4, 3)
    for k in (0, 1, 4, 5)
        @test_throws ArgumentError gali_value(W, k)
    end
end

@testset "AC7b non-finite and zero columns never give a finite result" begin
    _ = (gali_values,)
    # A NaN entry, an Inf entry or a zero column gives an error or a non-finite value, never a finite number (documented behaviour); for every order.
    nonfinite(M) = (r = try gali_values(ga_sm(M)) catch; nothing end; r === nothing || any(!isfinite, r))
    for T in GA_TYPES
        W = Matrix(ga_frame(ga_seed(7, 1), 1, T, 4, 3))
        Wn, Wi, Wz = copy(W), copy(W), copy(W)
        Wn[2, 2] = T(NaN)
        Wi[3, 1] = T(Inf)
        Wz[:, 3] .= zero(T)
        @test nonfinite(Wn) && nonfinite(Wi) && nonfinite(Wz)
    end
end

@testset "AC7c frame check" begin
    _ = (GALI, check_frame, gali_floor)
    # G_K(0) = 0 and 5 B_meas(T, K) are rejected (NumericalHealthError "initial frame", step 0, the orbit id, backend :cpu, threshold 10 B_meas), 20 B_meas is
    # accepted and check_frame returns G_K(0), for K in {2, 3, 4}; the seeded frame is accepted for orbit ids 1..N (D38: 1000 Float64, 200 Double64; spec 1e4)
    # with the minimum of G_K(0) / (10 B_meas) reported; for K = 2 the decision agrees with the sali frame check on 200 frames of G_2(0) = 10 B_meas 10^u,
    # u in [-1, 1], except within 2 B_meas(T, 2) of the threshold (the spec's factor 1.0013 is untestable: the computed values carry an absolute error of
    # order eps, which is 2.5 percent of the threshold 40 eps).
    rng = StableRNG(ga_seed(7, 2))
    for T in GA_TYPES, K in 2:4
        meas = ga_meas(T, K)
        for W in (ga_frame_with(T, K, 0.0), ga_sm(hcat(ga_frame_with(T, K, 0.0)[:, 1:(K - 1)], -ga_frame_with(T, K, 0.0)[:, 1])), ga_frame_with(T, K, 5 * meas))
            e = ga_health_error(() -> check_frame(GALI(), W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0 && e.backend === :cpu
            @test e isa NumericalHealthError && abs(e.threshold - 10 * meas) <= 1e-9 * 10 * meas
        end
        W20 = ga_frame_with(T, K, 20 * meas)
        g0 = nothing
        @test ga_health_error(() -> (g0 = check_frame(GALI(), W20; orbit = 3))) === nothing
        @test g0 !== nothing && abs(Float64(g0) - 20 * meas) <= 2 * meas
    end
    for T in GA_TYPES
        nid = T === Float64 ? ga_size(1000, 10_000) : ga_size(200, 1000)
        rejected, ratio = 0, Inf
        for j in 1:nid
            W = ga_frame(ga_seed(7, 3), j, T, 4, 4)
            e = ga_health_error(() -> check_frame(GALI(), W; orbit = j))
            rejected += e === nothing ? 0 : 1
            e === nothing && (ratio = min(ratio, Float64(check_frame(GALI(), W; orbit = j)) / (10 * ga_meas(T, 4))))
        end
        @test rejected == 0
        @info "AC7c reported: min G_K(0) / (10 B_meas) of the seeded frames" T ratio
        dis = 0
        for _ in 1:200
            g = 10 * ga_meas(T, 2) * 10.0^(2 * rand(rng) - 1)
            W = ga_frame_with(T, 2, g)
            ga_ok = ga_health_error(() -> check_frame(GALI(), W)) === nothing
            sa_ok = ga_health_error(() -> check_frame(SALI(), W)) === nothing
            gex = ga_exact(W)[1]
            if ga_ok != sa_ok && abs(Float64(gex) - 10 * ga_meas(T, 2)) > 2 * ga_meas(T, 2)
                dis += 1
            end
        end
        @test dis == 0
    end
end

@testset "AC7c rejection before any work" begin
    _ = (GALI, check_frame)
    # The caller runs the frame check first, so a rejected frame never reaches the map (ga_Boom throws if forward or tangent is called); an accepted frame proceeds.
    dep = ga_frame_with(Float64, 4, 0.0)
    W0 = [ga_frame(ga_seed(7, 4), 1, Float64, 4, 4), dep]
    run = () -> begin
        foreach(j -> check_frame(GALI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(ga_Boom(4), SVector(0.0, 0.0, 0.0, 0.0), W0[2], 4, 2, (GALI(),))
    end
    e = ga_health_error(run)
    @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 2 && e.step == 0
    @test_throws ErrorException propagate_orbit(ga_Boom(4), SVector(0.0, 0.0, 0.0, 0.0), W0[1], 4, 2, (GALI(),))
end

# ---- AC8: reducer structure on the stream (B5, B12, I5) ----
ga_same2(a, b) = isequal(a[1], b[1]) && isequal(a[2], b[2])

@testset "AC8a nesting of the leading columns" begin
    _ = (GALI,)
    # Y3 orbits C and D and seeded starts (D38: 3 chaotic starts at GA_Y3B and 2 small-amplitude starts at GA_Y3A, spec 20 each), n in {10, 100, 1000}, Float64:
    # G_2 and G_3 of a K = 4 run against runs with K = 2 and K = 3 on the identical orbit (frames: the first K columns of one Haar frame, which nest bitwise):
    # bitwise equal if the tangent product is column independent (reported), else within 2 B_meas(Float64, k); the test gates the bound.
    rng = StableRNG(ga_seed(8))
    T = Float64
    cases = vcat([(GA_Y3A, GA_ORBIT_C), (GA_Y3A, GA_ORBIT_D)], [(GA_Y3B, ga_draw_chaotic(rng, GA_Y3B)) for _ in 1:ga_size(3, 20)],
                 [(GA_Y3A, ga_draw_regular(rng, GA_Y3A)) for _ in 1:ga_size(2, 20)])
    bitwise, total = 0, 0
    for (j, (c, z)) in enumerate(cases)
        sys, x0 = ga_sys(T, c), ga_state(T, c, z)
        W4, W3, W2 = (ga_frame(ga_seed(8, 1), j, T, 4, K) for K in (4, 3, 2))
        @test W4[:, 1:3] == W3 && W4[:, 1:2] == W2
        steps = [10, 100, 1000]
        r4, r3, r2 = (ga_run1(GALI(), sys, x0, W, 1000, 1; steps = steps) for W in (W4, W3, W2))
        for i in 1:3
            for (k, g, h) in ((2, r4.snaps[i].gali[1], r2.snaps[i].gali[1]), (2, r4.snaps[i].gali[1], r3.snaps[i].gali[1]), (3, r4.snaps[i].gali[2], r3.snaps[i].gali[2]))
                @test abs(g - h) <= 2 * ga_meas(T, k)
                total += 1
                bitwise += isequal(g, h) ? 1 : 0
            end
        end
    end
    @info "AC8a reported: bitwise equal comparisons" bitwise total
end

# Records of a GALI run (K columns) and of a SALI run on the first two columns of the same frame, at the handoffs `steps`.
function ga_pair_runs(c::GACase, ::Type{T}, z, fid, K::Int, N::Int, steps) where {T}
    d = ga_d(c)
    sys, x0 = ga_sys(T, c), ga_state(T, c, z)
    Wk, W2 = ga_frame(fid[1], fid[2], T, d, K), ga_frame(fid[1], fid[2], T, d, 2)
    return ga_run1(GALI(), sys, x0, Wk, N, 1; steps = steps, orbit = fid[2]), ga_run1(SALI(), sys, x0, W2, N, 1; steps = steps, orbit = fid[2])
end

# Worst of |G_2 - SALI sqrt(4 - SALI^2)/2| / (2 B_meas(T, 2) + sali_floor(T)) over the resolved handoffs, and the number of handoffs where the sums are not bitwise equal.
function ga_vs_sali(rg, rs, ::Type{T}) where {T}
    worst, unequal, bitwiseL = 0.0, 0, 0
    tol = 2 * ga_meas(T, 2) + Float64(sali_floor(T))
    for (g, s) in zip(rg.snaps, rs.snaps)
        unequal += (ga_same2(g.sum, s.sum) && ga_same2(g.comp, s.comp)) ? 0 : 1
        G2 = g.gali[1]
        G2 >= 10 * ga_meas(T, 2) && (worst = max(worst, Float64(abs(G2 - s.sali * sqrt(4 - s.sali^2) / 2)) / tol))
    end
    return worst, unequal
end

@testset "AC8b K = 4 run against the sali reducer" begin
    _ = (GALI,)
    # The sums S_1, S_2 (and their compensation terms) of the K = 4 run are bitwise equal to those of the sali reducer in a k = 2 pass on the identical orbit
    # (first two columns of the same frame); G_2 agrees with the B4 conversion of the sali SALI within 2 B_meas(T, 2) + sali_floor(T) at resolved handoffs.
    # Y3 orbits C and D and 3 seeded chaotic starts, n in {10, 100, 1000}, Float64 and Double64.
    rng = StableRNG(ga_seed(8, 2))
    for T in GA_TYPES
        cases = vcat([(GA_Y3A, GA_ORBIT_C), (GA_Y3A, GA_ORBIT_D)], [(GA_Y3B, ga_draw_chaotic(rng, GA_Y3B)) for _ in 1:ga_size(3, 20)])
        worst, unequal = 0.0, 0
        for (j, (c, z)) in enumerate(cases)
            w, u = ga_vs_sali(ga_pair_runs(c, T, z, (ga_seed(8, 3), j), 4, 1000, [10, 100, 1000])..., T)
            worst, unequal = max(worst, w), unequal + u
        end
        @test unequal == 0
        @test worst <= 1
        @info "AC8b reported" T worst unequal
    end
end

@testset "AC8c K = d = 2 equals the sali reducer" begin
    _ = (GALI,)
    # K = d = 2 on the standard map, K in {0.5, 2 pi, 10} (D38: 4 starts per K, spec 20), n in {100, 1000}: the log form L is bitwise equal to the sali L, G_2 agrees
    # with the sali conversion within 2 B_meas + sali_floor. Float64 and Double64.
    rng = StableRNG(ga_seed(8, 4))
    for T in GA_TYPES, (label, Kv) in ((:half, 0.5), (:two_pi, 2 * Float64(pi)), (:ten, 10.0))
        sys = StandardMap(label === :two_pi ? exact(T, pi) * 2 : exact(T, Kv))
        unequalL, worst = 0, 0.0
        for j in 1:ga_size(4, 20)
            x0 = SVector{2, T}(T(2 * Float64(pi) * rand(rng)), T(2 * Float64(pi) * rand(rng)))
            W0 = ga_frame(ga_seed(8, 5), j, T, 2, 2)
            rg = ga_run1(GALI(), sys, x0, W0, 1000, 1; steps = [100, 1000])
            rs = ga_run1(SALI(), sys, x0, W0, 1000, 1; steps = [100, 1000])
            for (g, s) in zip(rg.snaps, rs.snaps)
                unequalL += isequal(g.L, s.L) ? 0 : 1
                g.gali[1] >= 10 * ga_meas(T, 2) && (worst = max(worst, Float64(abs(g.gali[1] - s.sali * sqrt(4 - s.sali^2) / 2)) / (2 * ga_meas(T, 2) + Float64(sali_floor(T)))))
            end
        end
        @test unequalL == 0
        @test worst <= 1
    end
end

@testset "AC8d Y4 chains, first two columns" begin
    _ = (GALI,)
    # Y4 N in {2, 5} at (K, eps) = (3, 0.2) with K = d (d = 4, 10) and K = 4: the first two columns behave as AC8b (sums bitwise equal to a sali k = 2 pass; G_2
    # against the conversion). D38: 2 seeded chaotic starts per case, n in {10, 100} (spec: n up to the AC9 horizons).
    rng = StableRNG(ga_seed(8, 6))
    T = Float64
    for (N, Ks) in ((2, (4,)), (5, (10, 4)))
        c = ga_y4(N, 3.0, 0.2)
        for K in Ks, j in 1:ga_size(2, 20)
            z = ga_draw_chaotic(rng, c)
            w, u = ga_vs_sali(ga_pair_runs(c, T, z, (ga_seed(8, 7), j), K, 100, [10, 100])..., T)
            @test u == 0
            @test w <= 1
        end
    end
end
