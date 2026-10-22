# AC4 to AC7: the regular law (RL6), chaotic independence of Dx (RL7), the schedule dependence (RL10), the resolution and the realized separation (RL12, RL13), and
# the recorded lines R1, R5, R8 that need no threshold. Every testset starts with a line naming the new symbols it uses.

# Own-orbit runs of pairs in the element type T (BigFloat: inside setprecision(BigFloat, bits)), RLI at the snapshots `ns` for each start and each Dx. The direction of the
# shadow is xi_0 (u = :xi0, P3) or the second seeded column (u = :xi2); the vector of both orbits is always xi_0 (P2). Returns Dict (j, Dx) => (Dxr, rli::Vector, zA).
function rl_own(::Type{T}, label::Symbol, m::Int, Dxs, ns; u::Symbol = :xi0, wrapped::Bool = true, bits = RL_BITS) where {T}
    body() = begin
        c = rl_case(label, T, m; wrapped = wrapped)
        N = maximum(ns)
        out = Dict{Tuple{Int, Float64}, Any}()
        for Dx in Dxs, j in 1:m
            dir = u === :xi0 ? c.xi0[j] : rl_xi2(c.fseed, c.ids[j], T)
            pr = rl_pair(c.sys, c.x0[j], dir, Dx)
            r = rl_run((RLI(),), PairMap(c.sys), [pr.z], [rl_frame(c.xi0[j])], N, 1; steps = ns)
            out[(j, Float64(Dx))] = (Dxr = pr.Dxr, rli = [f.RLI for f in r[1][1].snaps], zA = pr.zA)
        end
        out
    end
    return T === BigFloat ? setprecision(body, BigFloat, bits) : body()
end
# log10 of a possibly-BigFloat positive value as Float64
rl_log(x) = Float64(log(x))

# The tests of the regular law for one start and one n: the slope of log RLI against log Dx_r in [0.99, 1.01], and the deviation of rho = RLI / Dx_r from its value at the smallest Dx,
# dev_i = |rho(Dx_i) / rho(Dx_1) - 1|, is O(Dx): dev_i <= 2 (Dx_i / Dx_{i+1}) dev_{i+1} for consecutive Dx (the proposal of the spec, |rho_2 / rho_1 - 1| <= 20 Dx_2, fails for 5 of
# 20 regular-class starts at n = 100 because the second-order coefficient is large near a resonance, Deviations 3; the realized coefficient max dev_i / Dx_i is reported).
function rl_regular_law(res, j::Int, ni::Int, Dxs)
    dxr = [res[(j, Float64(D))].Dxr for D in Dxs]
    rli = [res[(j, Float64(D))].rli[ni] for D in Dxs]
    rhoT = [r / d for (r, d) in zip(rli, dxr)]
    rho = [Float64(r) for r in rhoT]
    slope = rl_slope(rl_log.(dxr), rl_log.(rli))
    dev = [Float64(abs(rhoT[i] / rhoT[1] - 1)) for i in 2:length(Dxs)]
    law_ok = all(dev[i] <= 2 * (Dxs[i + 1] / Dxs[i + 2]) * dev[i + 1] for i in 1:(length(dev) - 1))
    return (slope = slope, law_ok = law_ok, rho = rho, coeff = maximum(dev[i] / Dxs[i + 1] for i in eachindex(dev)))
end

@testset "AC4 regular law" begin
    _ = (RLI, PairMap)
    # RL6, F3. K = 0.5, regular class (L1 S7), wrapped, tau = 1. (a) BigFloat(256) own orbits for both orbits, Dx in {1e-14, 1e-12, 1e-10, 1e-8, 1e-6}: the slope of log RLI
    # against log Dx_r (least squares over the five values) lies in [0.99, 1.01]; the deviation of rho = RLI / Dx_r from its value at the smallest Dx is O(Dx) (the spec's
    # |rho(Dx_2)/rho(Dx_1) - 1| <= 20 Dx_2 is replaced, Deviations 3). (b) Double64, Dx in {1e-20, 1e-16, 1e-12, 1e-8}: the slope, and rho within 1e-3 (relative) of the BigFloat rho
    # at the same Dx. (c) Direction: u = the second seeded column (P3): the slope holds in Double64 and BigFloat, the O(Dx) law in BigFloat, rho of Double64 within 1e-3 of BigFloat;
    # the ratio rho(u = xi_0) / rho(u = xi_2) is reported. Sizes D38: 10 starts, n in {100, 1e3} (the spec: 20 starts, n up to 1e4).
    m = RL_FULL ? 20 : 10
    ns = RL_FULL ? (100, 1000, 10_000) : (100, 1000)
    DxA = (1e-14, 1e-12, 1e-10, 1e-8, 1e-6)
    DxB = (1e-20, 1e-16, 1e-12, 1e-8)
    bfA = rl_own(BigFloat, :half, m, DxA, ns)
    worst_coeff = 0.0
    for j in 1:m, ni in eachindex(ns)
        r = rl_regular_law(bfA, j, ni, DxA)
        @test 0.99 <= r.slope <= 1.01
        @test r.law_ok
        worst_coeff = max(worst_coeff, r.coeff)
    end
    @info "AC4(a) realized second-order coefficient max |rho_i / rho_1 - 1| / Dx_i over the starts (the spec's proposal: 20)" worst_coeff
    bfB = rl_own(BigFloat, :half, m, DxB, ns)
    dB = rl_own(Double64, :half, m, DxB, ns)
    worst = 0.0
    for j in 1:m, ni in eachindex(ns)
        r = rl_regular_law(dB, j, ni, DxB)
        @test 0.99 <= r.slope <= 1.01            # the O(Dx) law cannot be seen below the Double64 floor (Deviations 3): rho is compared with the BigFloat rho instead
        rb = rl_regular_law(bfB, j, ni, DxB)
        rel = maximum(abs(a / b - 1) for (a, b) in zip(r.rho, rb.rho))
        worst = max(worst, rel)
        @test rel <= 1e-3
    end
    @info "AC4(b) worst relative gap of rho, Double64 against BigFloat" worst
    dX = rl_own(Double64, :half, m, DxB[2:4], ns; u = :xi2)
    bX = rl_own(BigFloat, :half, m, DxB[2:4], ns; u = :xi2)
    for j in 1:m, ni in eachindex(ns)
        r2, rb2 = rl_regular_law(dX, j, ni, DxB[2:4]), rl_regular_law(bX, j, ni, DxB[2:4])
        @test 0.99 <= r2.slope <= 1.01 && 0.99 <= rb2.slope <= 1.01 && rb2.law_ok
        @test maximum(abs(a / b - 1) for (a, b) in zip(r2.rho, rb2.rho)) <= 1e-3
    end
    ratios = [dB[(j, 1e-12)].rli[ni] / dB[(j, 1e-12)].Dxr / (dX[(j, 1e-12)].rli[ni] / dX[(j, 1e-12)].Dxr) for j in 1:m for ni in eachindex(ns)]
    @info "AC4(c) ratio rho(u = xi_0) / rho(u = xi_2) at Dx = 1e-12 (reported, P3)" extrema(Float64.(ratios))
end

@testset "AC5 chaotic independence of the separation" begin
    _ = (RLI, PairMap)
    # RL7, F3. K in {2 pi, 5, 10}, Float64 wrapped, tau = 1, the pairs (Dx_1, Dx_2) in {(1e-12, 1e-10), (1e-12, 1e-8)} on the same starts (paired design), d_j = RLI_j(Dx_1) -
    # RLI_j(Dx_2), s_d its sample standard deviation: |mean d| <= 4 s_d / sqrt(m) + mean RLI(Dx_1) |ln(Dx_1 / Dx_2)| / (X_bar n), X_bar the ensemble mean of the L1 value X(n)
    # of a k = 1 run on the identical orbits (orbit A, xi_0). Derived allowance, a proposal (K18). Sizes D38: 1000 starts at n = 1e3, 250 at n = 1e4 (the spec: 1000 for both).
    # Reported: the ensemble mean RLI at the three Dx and the ratio to the regular-class mean at Dx = 1e-12 (the printed 9 to 10 orders, F3).
    T = Float64
    for (n, m) in ((1000, 1000), (10_000, RL_FULL ? 1000 : 250))
        reg = rl_case(:half, T, 20)
        regp = rl_pairs(reg, 1e-12)
        regm = rl_mean([f[1].snaps[1].RLI for f in rl_run((RLI(),), PairMap(reg.sys), [p.z for p in regp], rl_frames(reg), n, 1; steps = [n])])
        for label in (:two_pi, :five, :ten)
            c = rl_case(label, T, m)
            rli = Dict(Dx => [f[1].snaps[1].RLI for f in rl_run((RLI(),), PairMap(c.sys), [p.z for p in rl_pairs(c, Dx)], rl_frames(c), n, 1; steps = [n])]
                       for Dx in (1e-12, 1e-10, 1e-8))
            Xs = [f[1].snaps[1].lambda[1] for f in rl_run((Benettin(),), c.sys, [rl_astate(c.sys, x) for x in c.x0], [rl_col(x) for x in c.xi0], n, 1; steps = [n])]
            Xbar = rl_mean(Xs)
            for Dx2 in (1e-10, 1e-8)
                d = rli[1e-12] .- rli[Dx2]
                allow = 4 * sqrt(rl_var(d) / m) + rl_mean(rli[1e-12]) * abs(log(1e-12 / Dx2)) / (Xbar * n)
                @test abs(rl_mean(d)) <= allow
                @info "AC5 paired difference of the means against the allowance" label n m Dx2 diff = rl_mean(d) allow
            end
            @info "AC5 ensemble mean RLI at Dx = 1e-12, 1e-10, 1e-8 and the ratio to the regular class (R3)" label n rl_mean(rli[1e-12]) rl_mean(rli[1e-10]) rl_mean(rli[1e-8]) ratio = rl_mean(rli[1e-12]) / regm
        end
    end
end

@testset "AC6 tau dependence" begin
    _ = (RLI, PairMap)
    # RL10, I7. The same pair and xi_0, tau in {1, 2, 5, 10}, K in {0.5 (regular class), 2 pi}, n in {100, 1e3}, Dx in {1e-8, 1e-12}, wrapped: (a) the sums at a common n:
    # |s_i^tau(n) - s_i^1(n)| <= 8 eps Psi_n^(i) + 8 eps Psi_n^(i) (the L1 criterion 6 form, B_S^Psi for each of the two runs); (b) |RLI_tau(n) - RLI_1(n)| <= 2 (tau - 1) L_sigma
    # H_K / n + B_R^tau + B_R^1 with B_R of RL11 and K = n / tau. Reported: the ratio RLI_tau / RLI_1 for the regular class, where (b) is vacuous. Sizes D38: 20 starts in Float64,
    # 6 in Double64 (the spec: 20).
    for T in RL_TYPES, label in (:half, :two_pi), Dx in (1e-8, 1e-12)
        m = T === Float64 ? 20 : (RL_FULL ? 20 : 6)
        epsT = rl_eps(T)
        c = rl_case(label, T, m)
        prs = rl_pairs(c, Dx)
        refs = rl_pair_refs(label, T, Dx, m, 1000)
        Lsig = log(rl_sigma(Float64(c.K)))
        runs = Dict(tau => rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), 1000, tau; steps = [100, 1000]) for tau in (1, 2, 5, 10))
        rows = Dict(tau => rl_rows(label, T, Dx, tau, m, 1000) for tau in (1, 2, 5, 10))
        ratios = Float64[]
        for tau in (2, 5, 10), (ni, n) in enumerate((100, 1000)), j in 1:m
            ft, f1 = runs[tau][j][1].snaps[ni], runs[1][j][1].snaps[ni]
            o = refs[j]
            for i in 1:2
                sdiff = setprecision(() -> Float64(abs((rl_big(ft.sum[i]) + rl_big(ft.comp[i])) - (rl_big(f1.sum[i]) + rl_big(f1.comp[i])))), BigFloat, RL_BITS)
                @test sdiff <= 16 * epsT * rl_Psi(o, i, n)
            end
            Bt = only(r.BR for r in rows[tau] if r.j == j && r.n == n)
            B1 = only(r.BR for r in rows[1] if r.j == j && r.n == n)
            rdiff = setprecision(() -> Float64(abs(rl_big(ft.RLI) - rl_big(f1.RLI))), BigFloat, RL_BITS)
            @test rdiff <= 2 * (tau - 1) * Lsig * rl_H(n ÷ tau) / n + Bt + B1
            label === :half && push!(ratios, Float64(ft.RLI / f1.RLI))
        end
        label === :half && @info "AC6 ratio RLI_tau / RLI_1, regular class (reported)" T Dx extrema(ratios)
    end
end

@testset "AC7a realized separation and the pair start" begin
    _ = (pair_start, pair_distance, PairMap)
    # RL13. Dx_r is computed from the stored states: pair_start returns the pair state, Dx, Dx_r and the mode; the pair state equals the construction with displaced() of
    # lyap-two-trajectory bitwise; the components of the minimal-image difference of the two stored initial points equal the BigFloat difference exactly (Sterbenz) and Dx_r equals
    # the norm of it within 2 ulp (Deviations 1: "tolerance 0" cannot apply to a square root); |Dx_r - Dx| <= sqrt(2) eps M_state; both modes; 1000 starts, Dx in {1e-14, 1e-12, 1e-8}
    # (Double64: 200 starts). In wrapped mode a pair that straddles 0 mod 2 pi has the minimal-image separation, not the plain difference of the stored coordinates.
    seed = rl_seed(7, 1)
    starts = rl_starts(StableRNG(seed), 1000)
    for T in RL_TYPES, wrapped in (true, false), Dx in (1e-14, 1e-12, 1e-8)
        epsT = rl_eps(T)
        sys = rl_sys(rl_K(T, :half), wrapped)
        nst = T === Float64 ? 1000 : 200
        bad = 0
        for j in 1:nst
            x0 = rl_cast(T, starts[j])
            xi = rl_xi0(seed, j, T)
            rec = pair_start(sys, x0, xi, Dx; orbit = j)
            ref = rl_pair(sys, x0, xi, Dx)
            d = rl_dist_ref(ref.zA, ref.zB)
            ok = rec.z == ref.z && rec.Dx == T(Dx) && rec.mode == (wrapped ? "wrapped" : "lifted")
            ok &= abs(rec.Dx_r - ref.Dxr) <= 4 * epsT * ref.Dxr           # the norm of the verified separation (safe_norm); Systems has no Foundation, so not bitwise
            ok &= setprecision(() -> abs(rl_big(rec.Dx_r) - d.d) <= 4 * rl_big(epsT) * rl_big(rec.Dx_r), BigFloat, RL_BITS)
            ok &= abs(rec.Dx_r - T(Dx)) <= sqrt(2) * epsT * 2pi
            ok &= pair_distance(sys, rec.z) == rec.Dx_r
            bad += !ok
        end
        @test bad == 0
    end
    # the components of the minimal-image difference are exact: Float64, 1000 starts, Dx = 1e-12 (the difference of two stored numbers within a factor two of each other)
    sys = WrappedStandardMap(0.5)
    xcomp_bad = 0
    for j in 1:1000
        ref = rl_pair(sys, starts[j], rl_xi0(seed, j, Float64), 1e-12)
        dsep = separation(sys, ref.zA, ref.zB)
        d = rl_dist_ref(ref.zA, ref.zB)
        xcomp_bad += !(setprecision(() -> rl_big(dsep[1]) == d.dx && rl_big(dsep[2]) == d.dp, BigFloat, RL_BITS))
    end
    @test xcomp_bad == 0
    # a pair across the seam 0 = 2 pi: x just below 2 pi, shadow just above; and across the seam in p
    for T in RL_TYPES, (x0, u) in ((SVector(2pi - 1e-13, 1.0), SVector(1.0, 0.0)), (SVector(0.5, 2pi - 1e-13), SVector(0.0, 1.0)),
                                   (SVector(1e-14, 0.5), SVector(-1.0, 0.0)), (SVector(0.5, 1e-14), SVector(0.0, -1.0)))
        sys = WrappedStandardMap(rl_K(T, :half))
        xT, uT = rl_cast(T, x0), rl_cast(T, u)
        rec = pair_start(sys, xT, uT, 1e-12; orbit = 3)
        zA, zB = rl_half(rec.z, 1), rl_half(rec.z, 2)
        naive = max(abs(zB[1] - zA[1]), abs(zB[2] - zA[2]))
        @test naive > 1                                                      # the stored coordinates differ by about 2 pi
        @test abs(rec.Dx_r - T(1e-12)) <= sqrt(2) * rl_eps(T) * 2pi
        @test abs(pair_distance(sys, rec.z) - rec.Dx_r) == 0
        lifted_sep = setprecision(() -> Float64(sqrt((rl_big(zB[1]) + 2 * BigFloat(pi) * zB[3] - rl_big(zA[1]) - 2 * BigFloat(pi) * zA[3])^2 +
                                                      (rl_big(zB[2]) + 2 * BigFloat(pi) * zB[4] - rl_big(zA[2]) - 2 * BigFloat(pi) * zA[4])^2)), BigFloat, RL_BITS)
        @test abs(lifted_sep - 1e-12) <= 4 * sqrt(2) * rl_eps(T) * 2pi      # the winding counts make the lifted separation the minimal-image one (L2 T3)
    end
    # pair_distance of final states, both modes: the minimal-image distance of the stored numbers (BigFloat reference), 4 eps scale of the coordinates
    for wrapped in (true, false), label in (:half, :two_pi)
        c = rl_case(label, Float64, 4; wrapped = wrapped)
        prs = rl_pairs(c, 1e-8)
        for j in 1:4
            x, _, _ = propagate_orbit(PairMap(c.sys), prs[j].z, rl_frame(c.xi0[j]), 400, 1, (RLI(),); orbit = j)
            d = rl_dist_ref(rl_half(x, 1), rl_half(x, 2))
            scale = max(1.0, maximum(abs, x))
            @test abs(pair_distance(c.sys, x) - Float64(d.d)) <= 8 * eps(Float64) * scale
        end
    end
end

@testset "AC7b rejection of a separation that is not realized" begin
    _ = (pair_start, PairMap)
    # RL13, criterion 7(b). Dx = 1e-18 at |x| near 3 in Float64 (Dx_r = 0), a non-finite Dx, and a Dx whose Dx_r falls outside [Dx/2, 2 Dx] raise the numerical health error of
    # kind "separation" with the pair id before any work and nothing is stored; Double64 accepts 1e-18. The spec's example for the third case (Dx near 5e-16 at |x| near 3) is
    # not reachable: a rounding to the grid of spacing h moves a coordinate by at most h/2, so with both coordinates on the grid Dx_r/Dx lies in [1/2, 2] (Deviations 2); the
    # case used here puts p on the coarse grid (p near 5, spacing 8.9e-16) and x on a fine one (x near 0.1): the p component of 4.4e-16 is lost and Dx_r = 9e-17 < Dx/2.
    for wrapped in (false, true)
        sys = rl_sys(0.5, wrapped)
        u = SVector(0.6, 0.8)
        e = rl_health_error(() -> pair_start(sys, SVector(3.0, 1.0), u, 1e-18; orbit = 17))
        @test e isa NumericalHealthError && e.what == "separation" && e.orbit == 17 && e.step == 0 && e.backend === :cpu && e.value == 0.0 && e.threshold == 0.0
        for bad in (NaN, Inf, -Inf, 0.0)
            e = rl_health_error(() -> pair_start(sys, SVector(3.0, 1.0), u, bad; orbit = 5))
            @test e isa NumericalHealthError && e.what == "separation" && e.orbit == 5 && e.step == 0
        end
        # outside [Dx/2, 2 Dx], nonzero: the precondition is computed from the stored states with the verified separation
        u2 = SVector(0.2, sqrt(1 - 0.04))
        x0 = SVector(0.1, 5.0)
        ref = rl_pair(sys, x0, u2, 4.5e-16)
        @test 0 < ref.Dxr < 4.5e-16 / 2
        e = rl_health_error(() -> pair_start(sys, x0, u2, 4.5e-16; orbit = 9))
        @test e isa NumericalHealthError && e.what == "separation" && e.orbit == 9 && e.step == 0 && isapprox(e.value, Float64(ref.Dxr); rtol = 1e-14) && e.threshold == 4.5e-16
        # Double64 accepts 1e-18 at the same point
        sysD = rl_sys(rl_K(Double64, :half), wrapped)
        rec = pair_start(sysD, rl_cast(Double64, SVector(3.0, 1.0)), rl_cast(Double64, u), 1e-18; orbit = 17)
        @test Double64(1e-18) / 2 <= rec.Dx_r <= 2 * Double64(1e-18) && rec.Dx == Double64(1e-18)
        # the Dx must widen exactly into the element type (the L2 rule)
        @test_throws ArgumentError pair_start(sys, SVector(3.0, 1.0), u, BigFloat(1e-12); orbit = 1)
        # a rejected start never reaches the map: the caller starts the pair, then the stream (rl_Boom throws if forward or tangent is called)
        run_pair = () -> begin
            rec = pair_start(rl_Boom(), SVector(3.0, 1.0), u, 1e-18; orbit = 4)
            propagate_orbit(PairMap(rl_Boom()), rec.z, rl_frame(u), 4, 2, (RLI(),))
        end
        e = rl_health_error(run_pair)
        @test e isa NumericalHealthError && e.what == "separation" && e.orbit == 4
        ok = pair_start(rl_Boom(), SVector(3.0, 1.0), u, 1e-6; orbit = 4)                          # accepted, and then the stream does reach the map
        @test ok.Dx_r > 0
        @test_throws ErrorException propagate_orbit(PairMap(rl_Boom()), ok.z, rl_frame(u), 4, 2, (RLI(),))
    end
end

@testset "AC7c Float64 against BigFloat own orbits" begin
    _ = (RLI, PairMap)
    # RL12, F5. Float64 wrapped against BigFloat(256) own orbits, K = 0.5 regular class, n in {100, 1e3}, Dx in {1e-8, 1e-10, 1e-12}: at least 90 percent of the starts (18 of 20)
    # have |RLI_F64 / RLI_BF - 1| <= 0.1 (drafting: 0.5 percent at Dx = 1e-12, n = 2000). 20 starts (the spec). Reported (not gated): the first-order resolution number
    # Theta_n = (|K| + 4) M_state eps Phi_n / Dx_r of each start and the realized ratio |RLI_F64 / RLI_BF - 1| / Theta_n, also for Dx = 1e-14.
    m = 20
    ns = (100, 1000)
    Dxs = (1e-14, 1e-12, 1e-10, 1e-8)
    bf = rl_own(BigFloat, :half, m, Dxs, ns)
    f64 = rl_own(Float64, :half, m, Dxs, ns)
    c = rl_case(:half, Float64, m)
    for Dx in (1e-12, 1e-10, 1e-8), (ni, n) in enumerate(ns)
        good = count(j -> abs(f64[(j, Dx)].rli[ni] / Float64(bf[(j, Dx)].rli[ni]) - 1) <= 0.1, 1:m)
        @test good >= ceil(Int, 0.9 * m)
    end
    for Dx in Dxs, (ni, n) in enumerate(ns)
        rel = [abs(f64[(j, Dx)].rli[ni] / Float64(bf[(j, Dx)].rli[ni]) - 1) for j in 1:m]
        theta = [(0.5 + 4) * 2pi * eps(Float64) * rl_Phi(lb_tangents(c.sys, f64[(j, Dx)].zA, n), c.xi0[j]) / f64[(j, Dx)].Dxr for j in 1:m]
        @info "AC7c realized ratio |RLI_F64 / RLI_BF - 1| (R2) and Theta_n (reported)" Dx n max_rel = maximum(rel) median_theta = sort(theta)[(m + 1) ÷ 2] max_ratio = maximum(rel ./ theta)
    end
end

@testset "AC7d lifted against wrapped (R8)" begin
    _ = (RLI, PairMap)
    # Recorded line R8 and criterion 7(d): Float64, K = 0.5, a rotational start (p = 1), Dx = 1e-12, n = 1e4, tau = 1: the fraction of handoffs at which the unit in the last place of the
    # lifted coordinate exceeds Dx_r / 8 (the lifted pair cannot hold the separation, F6), and RLI / Dx_r in the two modes. Asserted: both runs are finite and nonnegative.
    x0 = SVector(0.3, 1.0)
    xi = SVector(0.6, 0.8)
    out = map((true, false)) do wrapped
        sys = rl_sys(0.5, wrapped)
        pr = rl_pair(sys, x0, xi, 1e-12)
        _, _, st = propagate_orbit(PairMap(sys), pr.z, rl_frame(xi), 10_000, 1, (RLI(), rl_Trace()))
        xs = reducer_final(rl_Trace(), st[2]).xs
        f = reducer_final(RLI(), st[1])
        coarse = count(z -> maximum(eps(abs(z[i])) for i in (1, 2, 3, 4)) > pr.Dxr / 8, xs) / length(xs)      # lifted pair state (x_A, p_A, x_B, p_B)
        (rho = f.RLI / pr.Dxr, coarse = wrapped ? 0.0 : coarse, rli = f.RLI)
    end
    @test all(o -> isfinite(o.rli) && o.rli >= 0, out)
    @info "AC7d R8 rho = RLI / Dx_r (wrapped, lifted) and the fraction of handoffs with ulp > Dx_r / 8 in lifted mode" rho_wrapped = out[1].rho rho_lifted = out[2].rho coarse_fraction = out[2].coarse
end
