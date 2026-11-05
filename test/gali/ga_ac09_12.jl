# AC9 to AC12: same-orbit bounds against the exact reference, normalisation interval and the K = 0 shear, chaotic law at hyperbolic fixed points, local dynamics
# near periodic orbits. Each testset starts by naming the new symbols it uses; the negative control (partial-sum log form) is a separate testset and passes now.

# Worst ratios (error / bound) of the series `snaps` (K-column run of type T) against the bundle b of the same start, over the checkpoints:
# G (direct form, B_G^k, k = 2..K), L (log form, B_L, K = d only), agree (the two forms above the floor, AC9c), floorbad (values above 2 B_meas where the
# reference is below B_meas, AC9d), nonfinite (L not finite).
function ga_check_orbit(b, snaps, K::Int, ::Type{T}) where {T}
    d, epsT = b.d, ga_eps(T)
    bits = T === BigFloat ? 2 * GA_BITS : GA_BITS            # the reference of an oracle-grade (BigFloat) run is at twice the precision
    r = (G = 0.0, L = 0.0, agree = 0.0, floorbad = 0, nonfinite = 0)
    for s in snaps
        ic = findfirst(==(s.n), b.cps)
        psi, lam, Gref = b.psi[ic], b.lam[ic], b.ref.G[s.n]
        for k in 2:K
            gr = Float64(Gref[k - 1])
            B = ga_BG(epsT, k, lam[k], gr, psi)
            r = merge(r, (G = max(r.G, ga_err(s.gali[k - 1], Gref[k - 1]; bits = bits) / B),))
            gr < ga_meas(T, k) && (r = merge(r, (floorbad = r.floorbad + (s.gali[k - 1] <= 2 * ga_meas(T, k) ? 0 : 1),)))
        end
        if K == d
            isfinite(s.L) || (r = merge(r, (nonfinite = r.nonfinite + 1,)))
            r = merge(r, (L = max(r.L, ga_err(s.L, ga_Lref(b.ref, s.n; bits = bits); bits = bits) / ga_BL(epsT, psi)),))
            gd = Float64(Gref[d - 1])
            if gd >= 10 * ga_meas(T, d)
                Bd = ga_BG(epsT, d, lam[d], gd, psi) + 8 * epsT * gd * sum(psi)
                r = merge(r, (agree = max(r.agree, Float64(abs(s.gali[d - 1] - exp(s.L))) / Bd),))
            end
        end
    end
    return r
end

# The cases of AC9 as (label, system, start, steps): orbits C, D, seeded starts of the chaotic and regular classes, Y4 chains.
# D38: n = 1000 (C) and 300 (D, spec 1e3), 150 for the chaotic starts, 2 starts per class (spec 20), Y4 N in {2, 3} (N = 5 only with GA_FULL), K in {1, 3} reduced to the
# chaotic panel point (3, 0.6102) / ring (1, 0.2) and a regular one (-1, 0.3).
function ga_ac9_cases()
    rng = StableRNG(ga_seed(9))
    m = ga_size(2, 20)
    cs = Any[("orbit C", GA_Y3A, GA_ORBIT_C, ga_size(1000, 10_000)), ("orbit D", GA_Y3A, GA_ORBIT_D, ga_size(300, 1000))]
    for _ in 1:m
        push!(cs, ("Y3 chaotic", GA_Y3B, ga_draw_chaotic(rng, GA_Y3B), ga_size(150, 1000)))
        push!(cs, ("Y4 N=2 chaotic", ga_y4(2, 3.0, 0.6102), ga_draw_chaotic(rng, ga_y4(2, 3.0, 0.6102)), ga_size(150, 1000)))
    end
    push!(cs, ("Y4 N=3 ring chaotic", ga_y4(3, 1.0, 0.2; ring = true), ga_draw_chaotic(rng, ga_y4(3, 1.0, 0.2; ring = true)), ga_size(150, 1000)))
    push!(cs, ("Y4 N=2 regular", ga_y4(2, -1.0, 0.3), ga_draw_regular(rng, ga_y4(2, -1.0, 0.3)), ga_size(1000, 10_000)))
    GA_FULL && push!(cs, ("Y4 N=5 chaotic", ga_y4(5, 3.0, 0.2), ga_draw_chaotic(rng, ga_y4(5, 3.0, 0.2)), 1000))
    return cs
end

@testset "AC9 same-orbit bounds against the exact reference" begin
    _ = (GALI,)
    # For each case, both element types: K-column runs with K = 4 (Y3, Y4) and K = d (Y3, Y4 N <= 3) on the first columns of the seeded frame, every checkpoint
    # (multiples of 10 at 10 per decade, D38). (a) |G_k,T - G_k,ref| <= B_G^k(n), k = 2..K; (b) |L - L_ref| <= B_L(n) and L finite everywhere (K = d);
    # (c) above the floor |G_d,direct - exp(L)| <= B_G^d(n) + 8 eps G_ref sum_i Psi_n^(i); (d) where G_k,ref < B_meas(T, k) the computed value is at most
    # 2 B_meas. Ratios error / bound are reported. The reference is the BigFloat(256) product at the realised orbit points (ga_bundle).
    worst = Dict{Any, Any}()
    for T in GA_TYPES, (j, (label, c, z, N)) in enumerate(ga_ac9_cases())
        d = ga_d(c)
        b = ga_bundle(c, T, z, (ga_seed(9, 1), j), N)
        for K in unique([4, d])
            r = ga_run1(GALI(), ga_sys(T, c), b.x0, b.W0[:, 1:K] |> ga_sm, N, 1; steps = b.cps, orbit = j)
            q = ga_check_orbit(b, r.snaps, K, T)
            @test length(r.snaps) == length(b.cps)
            @test q.G <= 1
            @test q.floorbad == 0
            if K == d
                @test q.L <= 1
                @test q.agree <= 1
                @test q.nonfinite == 0
            end
            w = get(worst, (T, K), (G = 0.0, L = 0.0))
            worst[(T, K)] = (G = max(w.G, q.G), L = max(w.L, q.L))
        end
    end
    @info "AC9 reported: worst ratios error / bound" worst
end

@testset "AC9e negative control: partial-sum log form" begin
    # A double that sums S_1..S_{d-1} only (it forgets the last column) differs from the exact L by S_d, of order l_1 n, and fails the bound of AC9(b) by more than
    # 10 B_L at n >= 50 for chaotic starts (test-side, from the 256-bit column sums of the reference; no new symbol).
    bad, total = 0, 0
    for (j, (label, c, z, N)) in enumerate(ga_ac9_cases())
        (startswith(label, "orbit C") || occursin("regular", label)) && continue
        b = ga_bundle(c, Float64, z, (ga_seed(9, 1), j), N)
        d = b.d
        for (ic, n) in enumerate(b.cps)
            n >= 50 || continue
            Lp, Lr = setprecision(BigFloat, GA_BITS) do
                (-sum(b.ref.S[n][1:(d - 1)]), ga_Lref(b.ref, n))
            end
            total += 1
            bad += ga_err(Lp, Lr) > 10 * ga_BL(ga_eps(Float64), b.psi[ic]) ? 1 : 0
        end
    end
    @test total > 0 && bad == total
end

@testset "AC10a normalisation interval" begin
    _ = (GALI,)
    # tau in {1, 2, 5, 10}, Y3: orbits C and D and 2 seeded chaotic starts (D38; spec 20 starts), n in {10, 100}, Float64: |G_k,tau - G_k,1| <= 2 B_G^k(n) for
    # k = 2, 3, 4 and |L_tau - L_1| <= 2 B_L(n) (K = d = 4); B from the tau = 1 reference of the bundle.
    rng = StableRNG(ga_seed(10))
    T = Float64
    cases = vcat([(GA_Y3A, GA_ORBIT_C), (GA_Y3A, GA_ORBIT_D)], [(GA_Y3B, ga_draw_chaotic(rng, GA_Y3B)) for _ in 1:ga_size(2, 20)])
    worstG, worstL = 0.0, 0.0
    for (j, (c, z)) in enumerate(cases)
        b = ga_bundle(c, T, z, (ga_seed(10, 1), j), 100)
        runs = Dict(tau => ga_run1(GALI(), ga_sys(T, c), b.x0, b.W0, 100, tau; steps = [10, 100], orbit = j) for tau in (1, 2, 5, 10))
        for (i, n) in enumerate((10, 100))
            ic = findfirst(==(n), b.cps)
            psi, lam, Gref = b.psi[ic], b.lam[ic], b.ref.G[n]
            for tau in (2, 5, 10), k in 2:4
                B = ga_BG(ga_eps(T), k, lam[k], Float64(Gref[k - 1]), psi)
                worstG = max(worstG, Float64(abs(runs[tau].snaps[i].gali[k - 1] - runs[1].snaps[i].gali[k - 1])) / (2 * B))
            end
            for tau in (2, 5, 10)
                worstL = max(worstL, Float64(abs(runs[tau].snaps[i].L - runs[1].snaps[i].L)) / (2 * ga_BL(ga_eps(T), psi)))
            end
        end
    end
    @test worstG <= 1
    @test worstL <= 1
    @info "AC10a reported" worstG worstL
end

# Closed forms and bounds of the K = 0 shear (B13): tangent [[I, I], [0, I]] each step, Y_n = [[I, nI], [0, I]], W0 = I. Columns 1, 2 keep their norm 1; columns 3, 4 have
# norm sqrt(j^2 + 1). G_2 = 1, G_3 = 1/sqrt(n^2 + 1), G_4 = 1/(n^2 + 1). Returns (closed forms, B_G^k for k = 2, 3, 4) in Float64 for epsT.
function ga_shear_bounds(epsT::Float64, n::Int)
    Gb = setprecision(BigFloat, GA_BITS) do
        (BigFloat(1), 1 / sqrt(BigFloat(n)^2 + 1), 1 / (BigFloat(n)^2 + 1))
    end
    G = Float64.(Gb)
    nw(i, j) = i <= 2 ? 1.0 : sqrt(j^2 + 1.0)
    psi = zeros(4)
    lam = (n * 1.0, sum(j -> sqrt((j^2 + 1.0) / (n^2 + 1.0)), 1:n), sum(j -> (j^2 + 1.0) / (n^2 + 1.0), 1:n))
    for j in 1:n
        a = n - j + 1.0
        op = (a + sqrt(a^2 + 4)) / 2
        for i in 1:4
            psi[i] += op * nw(i, j - 1) / nw(i, n)
        end
    end
    return Gb, ntuple(k -> ga_BG(epsT, k + 1, lam[k], G[k], psi), 3)
end

@testset "AC10b K = 0 shear and a seeded frame" begin
    _ = (GALI,)
    # Y3 with nu = kappa = mu = 0, W0 = I, tau in {1, 10}: |G_k - closed form| <= B_G^k(n) for k = 2, 3, 4 at n in {10, 100, 1e4, 1e6} (D38: tau = 1 to 1e5, tau = 10 to
    # 1e6, Double64 to 1e4); G_3 n and G_4 n^2 reported (they tend to 1). For a seeded frame the series is compared with the 256-bit product at n in {10, 100}
    # within B_G^k(n) (the check of AC9). Float64 and Double64.
    c = ga_y3(0.0, 0.0, 0.0)
    z = (0.3, 0.1, 0.25, -0.125)
    I4 = SMatrix{4, 4, Float64}(I)
    for T in GA_TYPES
        sys, x0 = ga_sys(T, c), ga_state(T, c, z)
        for tau in (1, 10)
            nmax = T === Float64 ? (tau == 1 ? ga_size(10^5, 10^6) : 10^6) : 10^4
            ns = filter(n -> n <= nmax, [10, 100, 10^4, 10^5, 10^6])
            r = ga_run1(GALI(), sys, x0, SMatrix{4, 4, T}(I4), nmax, tau; steps = ns)
            worst = 0.0
            lim = Tuple{Int, Float64, Float64}[]
            for (s, n) in zip(r.snaps, ns)
                G, B = ga_shear_bounds(ga_eps(T), n)
                for k in 1:3
                    worst = max(worst, ga_err(s.gali[k], G[k]) / B[k])
                end
                push!(lim, (n, Float64(s.gali[2]) * n, Float64(s.gali[3]) * n^2))
            end
            @info "AC10b reported: n, G_3 n, G_4 n^2 (tend to 1)" T tau lim
            @test worst <= 1
        end
        b = ga_bundle(c, T, z, (ga_seed(10, 2), 1), 100)
        r = ga_run1(GALI(), sys, x0, b.W0, 100, 1; steps = b.cps)
        q = ga_check_orbit(b, r.snaps, 4, T)
        @test q.G <= 1 && q.L <= 1
    end
end

# ---- AC11: chaotic law with exact exponents at hyperbolic fixed points ----
# Run at the origin (a fixed point of both systems) with K columns; returns the probe record at the handoffs ns. state zero, seeded frame (first K columns).
function ga_fixed_run(c::GACase, ::Type{T}, K::Int, fid, ns) where {T}
    d = ga_d(c)
    return ga_run1(GALI(), ga_sys(T, c), SVector{d, T}(ntuple(_ -> zero(T), d)), ga_frame(fid[1], fid[2], T, d, K), maximum(ns), 1; steps = ns, orbit = fid[2])
end

# Law check for the orders k = 2..K between two snapshots s1, s2 (steps n1, n2) against E_k: counts (ok, bad, unresolved) over the direct values resolved in both.
function ga_law_counts(s1, s2, E, K::Int, ::Type{T}) where {T}
    ok, bad, unres = 0, 0, 0
    for k in 2:K
        g1, g2 = s1.gali[k - 1], s2.gali[k - 1]
        if g1 >= 10 * ga_meas(T, k) && g2 >= 10 * ga_meas(T, k)
            dev = abs((log(Float64(g2)) - log(Float64(g1))) + (s2.n - s1.n) * E[k])
            dev <= 0.5 ? (ok += 1) : (bad += 1)
        else
            unres += 1
        end
    end
    return ok, bad, unres
end

@testset "AC11 chaotic law with exact exponents at hyperbolic fixed points" begin
    _ = (GALI,)
    # Y3 origin with (nu, kappa) in {(6, 0.1), (4.04, 0.1)} (mu = 0.3, b = 0 there) and Y4 origin, N in {2, 3, 5} (N = 10 with GA_FULL), open and ring (N >= 3), K = 0.05,
    # eps = 0.04, closed-form spectrum (ga_exponents: Laplacian eigenvalues, lambda + 1/lambda = |2 + g|); seeded frames (D38: 2 per system, spec 5), K = d for d <= 10 and
    # K = 6 otherwise. Pass: for the handoffs n_1 = 100, n_2 = 200 (spec: n_1 >= 40 and n_2 >= n_1 + 40; the windows [40, 80] are still in the transient for the tight spectra
    # of N = 5, where the log form deviates by 1.2, so the gated window is later and the [40, 80] deviations are only counted), |[ln G_k(n_2) - ln G_k(n_1)] + (n_2 - n_1) E_k| <= 0.5
    # for every order with both values resolved, in BigFloat(256) (D38: d <= 6) and in Double64; the log form (K = d) in Float64 additionally within B_L(n_1) + B_L(n_2) + 0.5.
    # Counts of unresolved orders are reported (BigFloat(256) resolves only about ln(1/(10 B_meas)) = 174, Double64 66, so the high orders are carried by the log form).
    cases = vcat([(ga_y3(6.0, 0.1, 0.3), "Y3 nu=6"), (ga_y3(4.04, 0.1, 0.3), "Y3 nu=4.04")],
                 [(ga_y4(N, 0.05, 0.04; ring = ring), "Y4 N=$N ring=$ring") for N in (GA_FULL ? (2, 3, 5, 10) : (2, 3, 5)) for ring in (false, true) if !(ring && N < 3)])
    ns = [40, 80, 100, 200]
    tot, early = (ok = 0, bad = 0, unres = 0), 0
    for (ci, (c, label)) in enumerate(cases), j in 1:ga_size(2, 5)
        d = ga_d(c)
        K = d <= 10 ? d : 6
        l = ga_exponents(c)
        E = Float64[Float64(ga_Ek(l, k)) for k in 1:K]
        for T in (BigFloat, Double64)
            T === BigFloat && !GA_FULL && d > 6 && continue
            r = ga_fixed_run(c, T, K, (ga_seed(11, 1), 10 * ci + j), ns)
            ok, bad, unres = ga_law_counts(r.snaps[3], r.snaps[4], E, K, T)
            tot = (ok = tot.ok + ok, bad = tot.bad + bad, unres = tot.unres + unres)
            @test bad == 0
            early += ga_law_counts(r.snaps[1], r.snaps[2], E, K, T)[2]
        end
        if K == d
            r = ga_fixed_run(c, Float64, K, (ga_seed(11, 1), 10 * ci + j), ns)
            Ms = fill(Float64.(Matrix(tangent(ga_sys(Float64, c), SVector{d, Float64}(ntuple(_ -> 0.0, d))))), 200)
            psi = ga_psis(Ms, ga_frame(ga_seed(11, 1), 10 * ci + j, Float64, d, d), [100, 200])
            dev = abs((r.snaps[4].L - r.snaps[3].L) + 100 * E[d])
            @test dev <= ga_BL(ga_eps(Float64), psi[1]) + ga_BL(ga_eps(Float64), psi[2]) + 0.5
        end
    end
    @test tot.ok > 0                                              # the law is not vacuous: some order is resolved
    @info "AC11 reported: orders checked, failed, unresolved in the gated window; direct orders failing in [40, 80]" tot early
end

# ---- AC12: local dynamics near periodic orbits ----
# ln G_k against ln n slopes over n in [n_lo, n_hi] from the snapshots of a run probed at every 10th handoff; NaN when a value in the window is below 10 B_meas.
function ga_slopes(snaps, n_lo, n_hi, K::Int, ::Type{T}) where {T}
    sel = filter(s -> n_lo <= s.n <= n_hi, snaps)
    return [(all(s -> s.gali[k - 1] >= 10 * ga_meas(T, k), sel) ?
             ga_slope(log.(Float64[s.n for s in sel]), Float64[log(Float64(s.gali[k - 1])) for s in sel]) : NaN) for k in 2:K]
end

@testset "AC12a stable fixed points" begin
    _ = (GALI,)
    # Y3 origin at (0.5, 0.1, 1e-3) and Y4 origin N in {2, 3} (N = 5 with GA_FULL), K = -1, eps = 0.3 (elliptic modes), seeded frames (D38: 2 per system), K = d, n <= 2000,
    # Float64: G_k in [0.05, 1 + B_meas] for all k = 2..K and every handoff, and the fitted slope of ln G_k against ln n over [100, 2000] has absolute value at most 0.2 (G_k
    # constant, B8 Eq. (38)); the slopes are reported.
    cases = vcat([(GA_Y3A, "Y3")], [(ga_y4(N, -1.0, 0.3), "Y4 N=$N") for N in (GA_FULL ? (2, 3, 5) : (2, 3))])
    for (ci, (c, label)) in enumerate(cases), j in 1:ga_size(2, 5)
        d = ga_d(c)
        r = ga_fixed_run(c, Float64, d, (ga_seed(12, 1), 10 * ci + j), collect(1:2000))
        @test all(s -> all(g -> 0.05 <= g <= 1 + ga_meas(Float64, d), s.gali), r.snaps)
        sl = ga_slopes(filter(s -> s.n % 10 == 0, r.snaps), 100, 2000, d, Float64)
        @test all(s -> abs(s) <= 0.2, sl)
        @info "AC12a reported: slopes" label j sl
    end
end

@testset "AC12b neighbouring regular orbit differs from the fixed point for k > N" begin
    _ = (GALI,)
    # Y3 orbit C (a neighbour of the stable origin), K = d = 4, n up to 2e4, Float64: the slopes over [1e3, 2e4] are within the band of AC14 (k = 2 in [-0.1, 0.1], k = 3 in
    # [-2.3, -1.7], k = 4 in [-4.3, -3.7]); the stable fixed point has slope 0 for every k (AC12a), so the difference for k = 3, 4 is 2 and 4 (reported). A Y4 N = 2 small
    # amplitude start is reported, not gated (it may be a lower-dimensional torus).
    N = 20_000
    c = GA_Y3A
    r = ga_run1(GALI(), ga_sys(Float64, c), ga_state(Float64, c, GA_ORBIT_C), ga_frame(ga_seed(12, 2), 1, Float64, 4, 4), N, 1; steps = collect(10:10:N))
    sl = ga_slopes(r.snaps, 1000, N, 4, Float64)
    @test -0.1 <= sl[1] <= 0.1 && -2.3 <= sl[2] <= -1.7 && -4.3 <= sl[3] <= -3.7
    c4 = ga_y4(2, -1.0, 0.3)
    z = ga_draw_regular(StableRNG(ga_seed(12, 3)), c4)
    r4 = ga_run1(GALI(), ga_sys(Float64, c4), ga_state(Float64, c4, z), ga_frame(ga_seed(12, 2), 2, Float64, 4, 4), N, 1; steps = collect(10:10:N))
    @info "AC12b reported: slopes of orbit C and of a Y4 N=2 neighbour" sl ga_slopes(r4.snaps, 1000, N, 4, Float64)
end
