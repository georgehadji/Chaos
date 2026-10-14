# AC5 to AC7 (and the data AC18 reuses): same-orbit tangent propagation, tau invariance, full BigFloat orbit.
# The series n X(n) at a handoff is read through lb_Probe (a wrapper reducer delegating to Benettin), because the stream has no
# emission channel (Q9); prefix invariance is tested separately in "API prefix invariance".

const LB_AC5_CASES = ((:half, 10_000), (:two_pi, 10_000), (:ten, 1_000))

# Runs of criterion 5: 20 seeded starts per case, tau in {1, 4}, T in {Float64, Double64}; one row per run.
function lb_ac5_data()
    return get!(LB_CACHE, :ac5) do
        seed = LB_SEED + LB_OFF + 50
        rows = NamedTuple[]
        for (ci, (label, N)) in enumerate(LB_AC5_CASES), T in LB_TYPES
            K = lb_K(T, label)
            sys = StandardMap(K)
            epsT = lb_eps(T)
            starts = lb_starts(StableRNG(seed + ci), T === Double64 ? (N > 1000 ? 6 : 10) : 20)   # Double64: 6 starts at n = 1e4, 10 at n = 1e3 (runtime)
            Kf = Float64(K)
            for (j, s) in enumerate(starts)
                x0 = lb_cast(T, s)
                W0 = initial_deviations(seed, j, Val(2), Val(1), T)
                pts, xN = lb_pts(sys, x0, N)
                Ms = [SMatrix{2, 2, Float64}(map(Float64, tangent(sys, z))) for z in pts]
                ref = lb_ref_std(K, pts, W0[:, 1], 100)                             # E1 at 256 bits at the realized points, every n (S3)
                for tau in (1, 4)
                    probe = lb_Probe(Benettin(), 100:100:N)
                    x, _, st = propagate_orbit(sys, x0, W0, N, tau, (probe,); orbit = j)
                    fin = reducer_final(probe, st[1])
                    G = lb_Gfirst(Ms, tau)
                    ratio = lb_worst_ratio(fin.snaps, ref, 100, epsT, G, Kf)
                    rv = lb_range_violation(fin.alphas, lb_Phis(Ms, tau), K, tau, epsT)
                    push!(rows, (; label, T, tau, j, N, ratio, G, nsnaps = length(fin.snaps), xeq = x == xN, rv))
                end
            end
        end
        rows
    end
end

# Runs of criterion 6: tau in {1, 2, 5, 10}, n in {100, 1000}, K in {0.5, 2 pi}, 20 starts; one row per (run, tau > 1).
function lb_ac6_data()
    return get!(LB_CACHE, :ac6) do
        seed = LB_SEED + LB_OFF + 60
        rows = NamedTuple[]
        for (ci, label) in enumerate((:half, :two_pi)), T in LB_TYPES
            K = lb_K(T, label)
            sys = StandardMap(K)
            epsT = lb_eps(T)
            starts = lb_starts(StableRNG(seed + ci), 20)
            for (j, s) in enumerate(starts)
                x0 = lb_cast(T, s)
                W0 = initial_deviations(seed, j, Val(2), Val(1), T)
                pts, _ = lb_pts(sys, x0, 1000)
                Ms = [SMatrix{2, 2, Float64}(map(Float64, tangent(sys, z))) for z in pts]
                for n in (100, 1000)
                    runs = map((1, 2, 5, 10)) do tau
                        probe = lb_Probe(Benettin())
                        _, _, st = propagate_orbit(sys, x0, W0, n, tau, (probe,); orbit = j)
                        fin = reducer_final(probe, st[1])
                        rv = lb_range_violation(fin.alphas, lb_Phis(Ms[1:n], tau), K, tau, epsT)
                        (; tau, S = fin.inner.sum[1], C = fin.inner.comp[1], G = lb_Gfirst(Ms, tau), rv)
                    end
                    for r in runs[2:end]
                        diff = setprecision(BigFloat, LB_BITS) do
                            Float64(abs(lb_big(r.S) + lb_big(r.C) - lb_big(runs[1].S) - lb_big(runs[1].C)))
                        end
                        push!(rows, (; label, T, n, tau = r.tau, j, diff,
                                     bound = lb_BS(epsT, n, r.G, Float64(K)) + lb_BS(epsT, n, runs[1].G, Float64(K)), rv = max(r.rv, runs[1].rv)))
                    end
                end
            end
        end
        rows
    end
end

@testset "AC5 same-orbit tangent propagation" begin
    # |n X(n) - ln(norm(Y_n w0)/norm(w0))_ref| <= B_S(n) at every handoff multiple of 100, Y_n the BigFloat(256) product of the E1
    # tangents at the realized points (all n), B_S of S3 as amended (Amendment D35: linear term plus 16 K n^2 eps). K = 0.5 (n = 1e4),
    # 2 pi (1e4), 10 (1e3); tau in {1, 4}; 20 starts (Double64: 6 at n = 1e4, 10 at n = 1e3, runtime). G_first > 1e4 is reported.
    rows = lb_ac5_data()
    @test length(rows) == 164                     # Float64 3 x 20 starts, Double64 6 + 6 + 10 starts, two tau each
    for r in rows
        @test r.ratio <= 1
        @test r.xeq                                   # the estimator leaves the orbit unchanged (I3)
        @test r.nsnaps == div(r.N, 100)
    end
    @info "AC5: runs with G_first > 1e4" count = count(r -> r.G > 1e4, rows) of = length(rows) worst_ratio = maximum(r.ratio for r in rows)
end

@testset "AC6 tau invariance" begin
    # |n X_tau(n) - n X_1(n)| <= B_S^tau(n) + B_S^1(n), same start and w0, tau in {2, 5, 10}, K in {0.5, 2 pi}, n in {100, 1e3}.
    rows = lb_ac6_data()
    @test length(rows) == 2 * 2 * 20 * 2 * 3
    for r in rows
        @test r.diff <= r.bound
    end
    @info "AC6: worst ratio" worst = maximum(r.diff / r.bound for r in rows)
end

@testset "AC7 full BigFloat orbit" begin
    # Regular class K = 0.5 (S7 screening), n in {10, 100}, 50 starts: |n X_T - n X_ref| <= 32 K M n^2 eps rho_n, rho_n =
    # sigma_max(Y_ref) / norm(Y_ref w0), M = M_state of S1 (max(1, |K|, largest lifted coordinate); Amendment D35: the printed bound
    # without it fails on 2 of 397 starts at n = 100); starts with rho_n > 50 are reported and at most 10 percent (5) may be excluded.
    ids, starts, rejected, seed = lb_regular()
    @info "AC7: regular-class screening" replaced = rejected
    for T in LB_TYPES, n in (10, 100)
        epsT = lb_eps(T)
        sys = StandardMap(exact(T, 0.5))
        excluded = 0
        for (id, s) in zip(ids, starts)
            W0 = initial_deviations(seed, id, Val(2), Val(1), T)
            f = lb_final(sys, lb_cast(T, s), W0, n, 1)
            err, rho, M = lb_ref_full(0.5, s, W0[:, 1], n, f.sum[1], f.comp[1])
            if rho > 50
                excluded += 1
                continue
            end
            @test err <= 32 * 0.5 * M * n^2 * epsT * rho
        end
        @test excluded <= 5
    end
end
