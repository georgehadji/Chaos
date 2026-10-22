# AC9 to AC10: tau dependence, precision ladder and representation. Ensemble-level and loose-interval criteria (K18): the intervals and allowances are the
# spec's proposals (P7), not source values.

@testset "AC9 tau dependence" begin
    _ = (MEGNO,)
    # Same start and frame, tau in {1, 2, 5, 10}, K in {0.5 (regular class), 2 pi}, n in {100, 1e3}, 20 starts, both columns, Float64 and Double64. (a) the sums at a
    # common n: |s_i^tau(n) - s_i^1(n)| <= 16 eps Psi_n^(i) (L1 criterion 6 form, the two bounds 8 eps Psi_n^(i)); (b) |Y_tau(n) - Y_1(n)| <= (tau - 1) L + B_M^tau + B_M^1;
    # (c) |Ybar_tau(n) - Ybar_1(n)| <= 3 (tau - 1) L + B_M^tau + B_M^1, L = ln sigma(K_std) of L1 E14, B_M the bound of M11 of the run. Reported: the distributions of the
    # differences and of the Ybar20 differences (no bound).
    m, N = 20, 1000
    for T in MG_TYPES, label in (:half, :two_pi)
        epsT = mg_eps(T)
        lns = log(mg_sigma(Float64(mg_K(T, label))))
        O = mg_orbit_ref(label, T, m, N)
        base = Dict((r.j, r.n) => r for r in mg_rows(label, T, 1, m, N))
        for tau in (2, 5, 10)
            rows = [r for r in mg_rows(label, T, tau, m, N) if r.n in (100, 1_000)]
            @test length(rows) == 2m
            dY, dYb, dY20 = Float64[], Float64[], Float64[]
            for r in rows
                r1 = base[(r.j, r.n)]
                o = O.per[r.j]
                for i in 1:2
                    @test abs(r.S[i] - r1.S[i]) <= 16 * epsT * mg_Psi(o, i, r.n)                       # (a)
                    Bsum = r.B[i] + r1.B[i]
                    @test abs(r.Y[i] - r1.Y[i]) <= (tau - 1) * lns + Bsum                               # (b)
                    @test abs(r.Ybar[i] - r1.Ybar[i]) <= 3 * (tau - 1) * lns + Bsum                     # (c)
                    push!(dY, Float64(r.Y[i] - r1.Y[i])); push!(dYb, Float64(r.Ybar[i] - r1.Ybar[i])); push!(dY20, Float64(r.Ybar20[i] - r1.Ybar20[i]))
                end
            end
            q(v) = (minimum(v), mg_median(v), maximum(v))
            @info "AC9 reported: Y_tau - Y_1, Ybar_tau - Ybar_1, Ybar20_tau - Ybar20_1 (min, median, max); bound (tau - 1) L" label T tau dY = q(dY) dYb = q(dYb) dY20 = q(dY20) L = lns
        end
    end
end

@testset "AC10 precision ladder and representation" begin
    _ = (MEGNO,)
    # (a) K = 2 pi and 5, 1000 starts, n = 1e3, Y_1/n and Ybar_1/n, Float64 vs Double64 on the same initial points: ensemble means within 4 combined standard errors
    # (L1 S5); K in {5, 10}, 200 starts, n = 1e4: same rule. (b) K = 2 pi, 50 starts, n up to 60, BigFloat(256) on its own orbit against the BigFloat reference at its realized
    # points (computed at 512 bits, so the reference is not of the precision it checks): all three series within B_M, both columns. (c) BigFloat vs Double64, n = 300, 100
    # starts, K in {2 pi, 5}, ensemble level, same rule. (d) K = 10, n = 1e4, Float64 (the 20 starts of AC6): Y_1 is finite, of the order of lambda n (between 1.4e4 and
    # 1.9e4, spec: about 1.6e4), within B_M of the reference, and not clipped or replaced.
    for label in (:two_pi, :five)
        a, b = mg_ens(Float64, label, 1000, 1_000), mg_ens(Double64, label, 1000, 1_000)
        @test mg_agree(a.Y, b.Y) && mg_agree(a.Ybar, b.Ybar)
    end
    for label in (:five, :ten)
        a, b = mg_ens(Float64, label, 200, 10_000), mg_ens(Double64, label, 200, 10_000)
        @test mg_agree(a.Y, b.Y) && mg_agree(a.Ybar, b.Ybar)
    end
    setprecision(BigFloat, MG_BITS) do
        rows = mg_rows(:two_pi, BigFloat, 1, 50, 60; bits = 512)
        @test length(rows) == 50 * length(mg_checkpoints(60))
        for f in (:Y, :Ybar, :Ybar20), i in 1:2
            @test mg_nbad(r -> r.e[f][i] <= r.B[i], rows) == 0
        end
        @info "AC10(b) realized ratio max |error| / B_M (BigFloat)" worst = maximum(r.ratio for r in rows)
        for label in (:two_pi, :five)
            a, b = mg_ens(BigFloat, label, 100, 300), mg_ens(Double64, label, 100, 300)
            @test mg_agree(a.Y, b.Y) && mg_agree(a.Ybar, b.Ybar)
        end
    end
    rows = [r for r in mg_rows(:ten, Float64, 1, 20, 10_000) if r.n == 10_000]
    @test length(rows) == 20
    @test mg_nbad(r -> isfinite(r.Y[1]) && 1.4e4 < r.Y[1] < 1.9e4, rows) == 0
    @test mg_nbad(r -> r.e.Y[1] <= r.B[1] && r.e.Y[2] <= r.B[2], rows) == 0
    @test mg_nbad(r -> abs(r.Y[1] - Float64(r.r.Y[1])) <= r.B[1], rows) == 0                    # not clipped, equal to the reference
end
