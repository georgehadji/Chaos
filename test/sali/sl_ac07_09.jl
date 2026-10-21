# AC7 to AC9: same-orbit BigFloat reference (direct form), log form and determinant identity, cross-check with the verified L1 estimator.
# Shared computations (reference products, Psi_n, runs) are memoized in SL_CACHE by sl_rows, so AC2, AC7 and AC8 reuse them.

# Panel of AC7 and AC8: (class, N, starts). D38: 8 Double64 starts instead of 20 (Float64 keeps 20); chaotic n up to 1e3
# (the direct form is at the floor after n_res, A14); checkpoints at 10 per decade (sl_checkpoints). SL_FULL restores 20 starts and every multiple of 10.
sl_panel(::Type{T}) where {T} = ((:half, 10_000, sl_size(T === Float64 ? 20 : 8, 20)),
                                 (:two_pi, 1_000, sl_size(T === Float64 ? 20 : 8, 20)),
                                 (:ten, 1_000, sl_size(T === Float64 ? 20 : 8, 20)))

@testset "AC7 same-orbit reference, direct form" begin
    _ = (SALI, sali_floor)
    # K in {0.5 (regular class, S7), 2 pi, 10}, tau = 1: |SALI_T(n) - SALI_ref(n)| <= B_SALI(n) at the checkpoints, with the BigFloat(256) product at the
    # realized points and B_SALI of A15 (Lambda_n from the reference, Psi_n from lb_Psi). Floor control for the chaotic cases: at every handoff with
    # SALI_ref < B_meas, SALI_T <= 2 B_meas (exact 0 included), and some handoffs are resolved (the test is not vacuous).
    for T in SL_TYPES, (label, N, m) in sl_panel(T)
        meas = sl_meas(T)
        rows = sl_rows(label, T, 1, m, N)
        @test length(rows) == m * length(sl_checkpoints(N))
        worst = maximum(r.errS / r.BS for r in rows)
        @test worst <= 1
        @info "AC7 realized ratio max |error| / B_SALI" label T worst
        if label !== :half
            below = [r for r in rows if r.saliref < meas]
            @test !isempty(below)
            @test all(r -> r.sali <= 2 * meas, below)
            @test any(r -> r.saliref >= 10 * meas, rows)
        end
    end
end

@testset "AC8 log form and determinant identity" begin
    _ = (SALI, sali_from_log)
    # Same orbits as AC7, T both, tau in {1, 4}: (a) |L(n) - L_ref(n)| <= B_L(n) at the checkpoints; (c) where the point is resolved (SALI_ref >= 10 B_meas)
    # the two forms agree: |SALI_direct - sali_from_log(L)| <= B_SALI + 8 eps (Psi1 + Psi2) SALI_ref; (d) negative control: a test double that sums only S_1
    # differs from L by S_2 (about lambda_1 n >= 40 for chaotic starts at n >= 50), so it fails (a) by more than 10 B_L. (b) is the next block.
    for T in SL_TYPES, (label, N, m) in sl_panel(T), tau in (1, 4)
        epsT, meas = sl_eps(T), sl_meas(T)
        rows = sl_rows(label, T, tau, m, N)
        @test maximum(r.errL / r.BL for r in rows) <= 1                                                   # (a)
        res = [r for r in rows if r.saliref >= 10 * meas]
        (tau == 1 || label === :half) && @test !isempty(res)              # chaotic Float64 at tau = 4 is below the floor at the first handoff n = 20
        @test all(r -> abs(Float64(r.sali - sali_from_log(r.L))) <= r.BS + 8 * epsT * (r.psi1 + r.psi2) * r.saliref, res)   # (c)
        if label !== :half
            late = [r for r in rows if r.n >= 50]
            @test !isempty(late)
            @test all(r -> abs(Float64(-r.S1) - r.Lref) > 10 * r.BL, late)                                # (d) the S_1-only double fails (a)
        end
    end
    # (b) L is finite at every handoff for K = 2 pi, n = 1e3, Float64, where exp(L) underflows
    c = sl_case(:two_pi, Float64, sl_size(8, 20))
    run = sl_run(SALI(), c.sys, c.x0, c.W0, 1_000, 1; all = true)
    Ls = [s.L for r in run for s in r.snaps]
    @test length(Ls) == length(run) * 1_000
    @test all(isfinite, Ls)
    @test minimum(Ls) < -745 && exp(minimum(Ls)) == 0
end

@testset "AC9 cross-check with L1 on the identical orbit" begin
    _ = (SALI, Benettin)
    # K in {0.5 (regular class), 2 pi, 10}, Float64, n in {100, 1e3, 1e4} (20 starts). The L1 series comes from Benettin on a k = 1 run on the
    # same orbit, started from the first (a) or the second (b) seeded column; B_S^Psi,i = 8 eps Psi_n^(i) (L1 Amendment D35), the bounds add (both sums are
    # computed bounds against the same reference): (a) |S_1 - n X^(1)| <= 2 B_S^Psi,1; (b) |S_2 - n X^(2)| <= 2 B_S^Psi,2; (c) |L + n (X^(1) + X^(2))| <= 2 B_L.
    T = Float64
    epsT = sl_eps(T)
    ns = [100, 1_000, 10_000]
    for label in (:half, :two_pi, :ten)
        m = 20
        c = sl_case(label, T, m)
        run = sl_run(SALI(), c.sys, c.x0, c.W0, 10_000, 1; steps = ns)
        l1a = sl_run(Benettin(), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 1]) for W in c.W0], 10_000, 1; steps = ns)
        l1b = sl_run(Benettin(), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 2]) for W in c.W0], 10_000, 1; steps = ns)
        gaps = Float64[]
        for j in 1:m
            Ms = lb_tangents(c.sys, c.x0[j], 10_000)
            for (i, n) in enumerate(ns)
                s = run[j].snaps[i]
                nX1 = l1a[j].snaps[i].sum[1] + l1a[j].snaps[i].comp[1]
                nX2 = l1b[j].snaps[i].sum[1] + l1b[j].snaps[i].comp[1]
                @test s.n == n && l1a[j].snaps[i].n == n && l1b[j].snaps[i].n == n
                p1, p2 = lb_Psi(view(Ms, 1:n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:n), c.W0[j][:, 2])
                @test abs((s.sum[1] + s.comp[1]) - nX1) <= 2 * 8 * epsT * p1                      # (a)
                @test abs((s.sum[2] + s.comp[2]) - nX2) <= 2 * 8 * epsT * p2                      # (b)
                @test abs(s.L + nX1 + nX2) <= 2 * sl_BL(epsT, p1, p2)                             # (c)
                n == 10_000 && push!(gaps, abs(nX2 - nX1) / 2)
            end
        end
        @info "AC9 reported: |ln(norm w2 / norm w1)| / 2 at n = 1e4 (A8, of order 1)" label median = sl_median(gaps)
    end
end
