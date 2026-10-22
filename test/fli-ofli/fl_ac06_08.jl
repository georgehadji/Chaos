# AC6 to AC8: same-orbit BigFloat reference with the two-point bound of Q11, basis-free sandwich of Q7, cross-check with L1 and sali on the identical orbit.
# Shared computations (reference products, Psi_n, runs, bounds) are memoized in FL_CACHE by fl_rows, so AC6, AC7 and AC10 reuse them.

# FLI_B at step n (tau = 1) of a double whose accumulator is a Float32: the increments of the exact sums, rounded to Float32 and added in Float32, with the
# running maximum taken in Float32. Used only as the negative control of AC6 (it must exceed B_FLI).
function fl_single_fli(o::fl_Orbit, n::Int)
    return setprecision(BigFloat, FL_BITS) do
        acc1 = acc2 = mx = Float32(0)
        p1 = p2 = BigFloat(0)
        for i in 1:n
            s1, s2 = log(o.q1[i]) / 2, log(o.q2[i]) / 2
            acc1 += Float32(s1 - p1)
            acc2 += Float32(s2 - p2)
            p1, p2 = s1, s2
            mx = max(mx, acc1, acc2)
        end
        Float64(mx)
    end
end

@testset "AC6 same-orbit reference" begin
    _ = (FLI,)
    # Standard map, K in {0.5 (regular class), 2 pi, 10}, 20 starts each, tau in {1, 4}, n up to 1e4, at the checkpoints (10 per decade; every
    # multiple of 10 with FL_FULL): |FLI_1 - FLI_ref,1| <= B_FLI(n) and |FLI_B - FLI_ref,B| <= B_FLI(n), the reference running maximum in BigFloat(256) at the realized
    # points and B_FLI = 8 eps max(Psi at the computed maximiser, Psi at the reference maximiser) (Q11). Reported: the realized ratio error / B_FLI and the number of
    # (start, checkpoint) pairs whose maximiser differs. Negative control: a double with a Float32 accumulator exceeds B_FLI at n = 1e3.
    for T in FL_TYPES, (label, N, m) in fl_panel(T), tau in (1, 4)
        rows = fl_rows(label, T, tau, m, N)
        @test length(rows) == m * count(n -> n % tau == 0, fl_checkpoints(N))
        worst1, worstB = maximum(r.ratio1 for r in rows), maximum(r.ratioB for r in rows)
        @test worst1 <= 1
        @test worstB <= 1
        last_n = maximum(r.n for r in rows)
        @info "AC6 realized ratio max |error| / B_FLI (R2); pairs and final-checkpoint orbits with a differing maximiser" label T tau worst1 worstB pairs = count(r -> r.maxdiff, rows) final_orbits = count(r -> r.maxdiff && r.n == last_n, rows)
    end
    for (label, N, m) in fl_panel(Float64)
        label === :half && continue
        O = fl_orbit_ref(label, Float64, m, N)
        rows = fl_rows(label, Float64, 1, m, N)
        for j in 1:3
            r = only(r for r in rows if r.j == j && r.n == 1000)
            @test abs(fl_single_fli(O.per[j], 1000) - Float64(r.FB)) > r.BB                  # the Float32 accumulator is far outside the bound
        end
    end
end

@testset "AC7 basis-free sandwich" begin
    _ = (FLI,)
    # Same orbits as AC6, tau = 1: F_*(n) - ln sqrt 2 - B_FLI(n) <= FLI_B(n) <= F_*(n) + B_FLI(n) and FLI_1(n) <= F_*(n) + B_FLI(n), F_* from the BigFloat prefix
    # products (the largest squared Frobenius norm over the handoffs; sigma_1^2 = (G + sqrt(G^2 - 4)) / 2). Reported: the distribution of (F_* - FLI_B) / ln sqrt 2
    # (in [0, 1]) and of F_* - FLI_1.
    lnsq2 = setprecision(() -> log(BigFloat(2)) / 2, BigFloat, FL_BITS)
    for T in FL_TYPES, (label, N, m) in fl_panel(T)
        rows = fl_rows(label, T, 1, m, N)
        @test all(r -> fl_diff(r.fli, r.Fstar) <= r.BB, rows)                       # FLI_B <= F_*
        @test all(r -> fl_diff(r.fli, r.Fstar - lnsq2) >= -r.BB, rows)              # FLI_B >= F_* - ln sqrt 2
        @test all(r -> fl_diff(r.fli1, r.Fstar) <= r.B1, rows)                      # FLI_1 <= F_*
        u = [-fl_diff(r.fli, r.Fstar) / Float64(lnsq2) for r in rows]
        v = [-fl_diff(r.fli1, r.Fstar) for r in rows]
        @info "AC7 reported: (F_* - FLI_B) / ln sqrt 2 (min, median, max) and F_* - FLI_1 (min, median, max)" label T u = (minimum(u), fl_median(u), maximum(u)) v = (minimum(v), fl_median(v), maximum(v))
    end
end

@testset "AC8 cross-check with L1 and sali on the identical orbit" begin
    _ = (FLI, SALI, Benettin)
    # K in {0.5 (regular class), 2 pi, 10}, 20 starts, Float64, tau = 1, n up to 1e4, the chain {FLI (not the replacer), SALI}. (a) the sums S_1, S_2 of the FLI
    # reducer are bitwise equal to those of the sali reducer (totals, and sum and compensation separately), and bitwise equal to the L1 sums of k = 1 runs from the
    # first and from the second seeded column (I6); (b) FLI_1(n) against the running maximum of n X_L1(i) over the stored series (every handoff, with the 0): within
    # 4 eps max(1, FLI_1(n)); (c) FLI_B(n) >= max over handoffs of (-L/2) - B_L(n)/2 and max(S_1, S_2) = -L/2 + |S_1 - S_2|/2 within B_L(n) at the checkpoints,
    # L the sali log form, B_L = 8 eps (Psi1 + Psi2). Reported (R1): the drawdown D(n) = FLI_1(n) - n X_L1(n) and FLI_B - FLI_1 at n in {1e2, 1e3, 1e4}.
    T = Float64
    epsT = fl_eps(T)
    N = 10_000
    cps = fl_checkpoints(N)
    for label in (:half, :two_pi, :ten)
        m = 20
        c = fl_case(label, T, m)
        run = fl_run((FLI(replace = false), SALI()), c.sys, c.x0, c.W0, N, 1; steps = cps, all = (false, true))
        l1 = [fl_run((Benettin(),), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, i]) for W in c.W0], N, 1; all = true) for i in 1:2]
        drawdown = Dict(n => Float64[] for n in (100, 1_000, 10_000))
        basis_gap = Dict(n => Float64[] for n in (100, 1_000, 10_000))
        for j in 1:m
            Ms = lb_tangents(c.sys, c.x0[j], N)
            fsn, ssn = run[j][1].snaps, run[j][2].snaps
            X1, X2 = l1[1][j][1].snaps, l1[2][j][1].snaps
            @test length(fsn) == length(cps) && length(ssn) == N && length(X1) == N && length(X2) == N
            rm = accumulate(max, [i * X1[i].lambda[1] for i in 1:N]; init = 0.0)     # running maximum of n X(n), with the s = 0 term
            runL = accumulate(max, [-s.L / 2 for s in ssn])                           # max over handoffs of -L/2
            for s in fsn
                n = s.n
                sa = ssn[n]
                S1, S2 = s.sum[1] + s.comp[1], s.sum[2] + s.comp[2]
                @test sa.n == n && X1[n].n == n && X2[n].n == n
                @test S1 == sa.sum[1] + sa.comp[1] && S2 == sa.sum[2] + sa.comp[2]                    # (a) sali
                @test s.sum == sa.sum && s.comp == sa.comp
                @test S1 == X1[n].sum[1] + X1[n].comp[1] && S2 == X2[n].sum[1] + X2[n].comp[1]        # (a) L1 per column
                @test abs(s.fli1 - rm[n]) <= 4 * epsT * max(1, s.fli1)                                # (b)
                p1, p2 = lb_Psi(view(Ms, 1:n), c.W0[j][:, 1]), lb_Psi(view(Ms, 1:n), c.W0[j][:, 2])
                BL = 8 * epsT * (p1 + p2)
                @test s.fli >= runL[n] - BL / 2                                                       # (c)
                @test abs(max(S1, S2) - (-sa.L / 2 + abs(S1 - S2) / 2)) <= BL
                if haskey(drawdown, n)
                    push!(drawdown[n], s.fli1 - n * X1[n].lambda[1])
                    push!(basis_gap[n], s.fli - s.fli1)
                end
            end
        end
        for n in (100, 1_000, 10_000)
            @info "AC8 reported (R1): FLI_1 - n X_L1 and FLI_B - FLI_1 (median, max)" label n D = (fl_median(drawdown[n]), maximum(drawdown[n])) FLIB_minus_FLI1 = (fl_median(basis_gap[n]), maximum(basis_gap[n]))
        end
    end
end
