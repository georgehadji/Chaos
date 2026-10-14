# AC1 to AC5 of lyap-gram-schmidt. theta_sum as G13 with eps_eff; frames and starts as in the spec; seeds gs_seed(ac).

@testset "AC1 cat map" begin
    # G12a: |lambda_1 - [ln mu + (1/(2n)) ln(c_+^2 + c_-^2 mu^-4n)]| <= theta_sum, |lambda_1 + lambda_2| <= theta_sum;
    # eigenbasis frame: |lambda_1 - ln mu| <= theta_sum. 20 seeded frames (orbit ids 1..20), z_0 = 0 (fixed point).
    seed = gs_seed(1)
    lnmu = gs_lnmu()
    for T in gs_types(), alg in gs_algsT(T), tau in (1, 2, 4), n in (40, 400)
        sys = gs_cat(T)
        x0s = [gs_zero(T, 2) for _ in 1:20]
        res = gs_run(sys, x0s, n; k = 2, tau, seed, alg)
        for j in 1:20
            f = gs_fin(res, j)
            th = gs_theta(T, f.G_max, tau, 2.0)
            q1 = initial_deviations(seed, j, Val(2), Val(2), T)[:, 1]
            @test gs_absdiff(f.lambda[1], gs_cat_lambda1(q1, n)) <= th
            @test Float64(abs(f.lambda[1] + f.lambda[2])) <= th
            @test Float64(f.orth_max) <= gs_orth_bound(alg, T, 2, f.G_max)
        end
        W0 = gs_cat_eigframe(T)
        res = gs_run(sys, x0s[1:1], n; k = 2, tau, seed, alg, W0 = [W0])
        f = gs_fin(res, 1)
        @test gs_absdiff(f.lambda[1], lnmu) <= gs_theta(T, f.G_max, tau, 2.0)
    end
end

@testset "AC2 K = 0" begin
    # G12b: Q_0 = (e_2, e_1), Y_n = [[1, n], [0, 1]]: S_1 = (1/2) ln(n^2 + 1) to 8 n eps, S_2 = -S_1 to 8 n eps,
    # lambda_1(1e4) <= ln(1e4)/1e4 + 1e-12 (Skokos2010 Eq. (62) row of benchmarks.md, regular decay).
    seed = gs_seed(2)
    for T in (Float64,), alg in gs_algsT(T)
        sys = StandardMap(zero(T))
        x0 = SVector{2,T}(T(0.3), T(0.7))
        W0 = SMatrix{2,2,T}(0, 1, 1, 0)
        for n in (1, 10, 1000, 10^4)
            res = gs_run(sys, [x0], n; k = 2, tau = 1, seed, alg, W0 = [W0])
            f = gs_fin(res, 1)
            S = gs_S(f)
            ref = setprecision(BigFloat, 256) do
                log(BigFloat(n)^2 + 1) / 2
            end
            @test gs_absdiff(S[1], ref) <= 8 * n * gs_eps(T)
            @test Float64(abs(S[2] + S[1])) <= 8 * n * gs_eps(T)
            n == 10^4 && @test Float64(f.lambda[1]) <= log(1e4) / 1e4 + 1e-12
        end
    end
end

@testset "AC3 bitwise structure" begin
    seed = gs_seed(3)
    # (a) k = 2: CGS and MGS perform identical operations, so S_i, the final Q (and every window's Q, R) agree bitwise
    for K in (0.5, 2.0, 5.0, 10.0), tau in (1, 3)
        sys = gs_stdmap(Float64, K)
        x0s = gs_starts(Float64, seed + 1, 100)
        rc = gs_run(sys, x0s, 1200; k = 2, tau, seed, alg = CGS())
        rm = gs_run(sys, x0s, 1200; k = 2, tau, seed, alg = MGS())
        @test all(j -> isequal(gs_fin(rc, j).sum, gs_fin(rm, j).sum) && isequal(gs_fin(rc, j).comp, gs_fin(rm, j).comp), 1:100)
        @test isequal(rc.W, rm.W)
    end
    rng = StableRNG(seed + 2)
    for (T, d) in ((Float64, 2), (Float64, 4), (Double64, 2)), _ in 1:200
        W = SMatrix{d,2,T}(ntuple(_ -> T(randn(rng)), Val(2d)))
        fc, fm = orthonormalize(CGS(), W), orthonormalize(MGS(), W)
        @test isequal(fc[1], fm[1]) && isequal(fc[2], fm[2])
    end
    # (b) A_S, k = 4, tau = 3, 200 windows (window map A_S^3): CGS and MGS differ in Q for k >= 3; tolerance 0
    A3 = gs_AS(Float64)^3
    Qc = Qm = initial_deviations(seed + 4, 1, Val(4), Val(4), Float64)
    differs = false
    for _ in 1:200
        Qc = orthonormalize(CGS(), A3 * Qc)[1]
        Qm = orthonormalize(MGS(), A3 * Qm)[1]
        differs |= Qc != Qm
    end
    @test differs
    sys = LinearSymplecticMap(gs_AS(Float64))
    rc = gs_run(sys, [gs_zero(Float64, 4)], 600; k = 4, tau = 3, seed, alg = CGS(), tol = Inf)   # Inf: measure, do not raise
    rm = gs_run(sys, [gs_zero(Float64, 4)], 600; k = 4, tau = 3, seed, alg = MGS(), tol = Inf)
    @test rc.W != rm.W
end

@testset "AC4 256-bit recomputation" begin
    # tangent at 256 bits at the stored orbit points, same widened Q_0, 256-bit RGS chain (never a separate 256-bit orbit):
    # |lambda_i(n) - lambda_i^ref(n)| <= theta_sum, i = 1, 2, all variants. The reference does not depend on tau.
    seed = gs_seed(4)
    n = 1000
    for (T, Ks, m, taus) in ((Float64, gs_KS(), 20, (1, 2, 4)), (Double64, (2.0, :twopi), 4, (1, 4)))
        for K in Ks
            sys = gs_stdmap(T, K)
            x0s = gs_starts(T, seed + 1, m)
            pts = [gs_points(sys, x0s[j], n) for j in 1:m]
            Mmax = [gs_Mmax(sys, pts[j]) for j in 1:m]
            ref = [gs_ref_std_S(sys.K, pts[j], initial_deviations(seed, j, Val(2), Val(2), T), n) for j in 1:m]
            for alg in gs_algsT(T), tau in taus
                res = gs_run(sys, x0s, n; k = 2, tau, seed, alg)
                for j in 1:m
                    f = gs_fin(res, j)
                    th = gs_theta(T, f.G_max, tau, Mmax[j])
                    @test gs_absdiff(f.lambda[1], ref[j][1] / n) <= th
                    @test gs_absdiff(f.lambda[2], ref[j][2] / n) <= th
                    @test Float64(f.orth_max) <= gs_orth_bound(alg, T, 2, f.G_max)
                end
            end
        end
    end
end

# Fixed-point runs of the Y5 pair; a CGS run that raises the orthogonality error is reported, not scored (AC6c).
function gs_y5_runs(seed, taus, n, m)
    A = gs_AS(Float64)
    sys = LinearSymplecticMap(A)
    out = Dict{Tuple{String,Int},Any}()
    trips = Tuple{String,Int}[]
    for alg in gs_algs(), tau in taus
        res, err = gs_try_run(sys, [gs_zero(Float64, 4) for _ in 1:m], n; k = 4, tau, seed, alg)
        if err === nothing
            out[(gs_name(alg), tau)] = res
        else
            @test alg isa CGS && err.what == "orthogonality"
            push!(trips, (gs_name(alg), tau))
        end
    end
    return A, out, trips
end

@testset "AC5 cross variant and cross tau" begin
    # I1: every (variant, tau) pair agrees: |lambda_i^A - lambda_i^B| <= theta_sum^A + theta_sum^B (G13 pair bound)
    seed = gs_seed(5)
    n = 1000
    for K in (0.5, 2.0, 5.0, 10.0, :twopi)
        sys = gs_stdmap(Float64, K)
        x0s = gs_starts(Float64, seed + 1, 50)
        pts = [gs_points(sys, x0s[j], n) for j in 1:50]
        Mmax = [gs_Mmax(sys, p) for p in pts]
        runs = Dict((gs_name(alg), tau) => gs_run(sys, x0s, n; k = 2, tau, seed, alg) for alg in gs_algs(), tau in (1, 2, 4))
        algof = Dict(gs_name(a) => a for a in gs_algs())
        for ((nm, _), r) in runs    # AC6(a), (b) on the runs of this criterion
            @test all(j -> Float64(gs_fin(r, j).orth_max) <= gs_orth_bound(algof[nm], Float64, 2, gs_fin(r, j).G_max), 1:50)
        end
        keys_ = collect(keys(runs))
        for a in 1:length(keys_), b in (a + 1):length(keys_), j in 1:50
            fa, fb = gs_fin(runs[keys_[a]], j), gs_fin(runs[keys_[b]], j)
            th = gs_theta(Float64, fa.G_max, keys_[a][2], Mmax[j]) + gs_theta(Float64, fb.G_max, keys_[b][2], Mmax[j])
            @test all(i -> abs(fa.lambda[i] - fb.lambda[i]) <= th, 1:2)
        end
    end
    # Y5 pair A_S, k = 4, n = 600, tau in {1, 2}, 10 seeded frames (orbit ids)
    A, runs, trips = gs_y5_runs(seed + 2, (1, 2), 600, 10)
    @test ("CGS", 1) ∉ trips && ("MGS", 2) ∉ trips && ("RGS", 2) ∉ trips   # only CGS at tau = 2 may trip the check
    keys_ = collect(keys(runs))
    for a in 1:length(keys_), b in (a + 1):length(keys_), j in 1:10
        fa, fb = gs_fin(runs[keys_[a]], j), gs_fin(runs[keys_[b]], j)
        th = gs_theta_lin(Float64, A, keys_[a][2]) + gs_theta_lin(Float64, A, keys_[b][2])
        @test all(i -> abs(fa.lambda[i] - fb.lambda[i]) <= th, 1:4)
    end
end
