# AC10, AC11 and AC14 of lyap-gram-schmidt (AC12, AC13: stage 04). Seeds gs_seed(ac).

@testset "AC10 agreement with Benettin" begin
    # k = 1 (G6, Skokos2010 p.43): |lambda_1^GS - X| <= (B_S^GS + B_S^L1)/n, B_S of lyap-benettin S3 as amended (Amendment D35), G_first = ||P_1||_F^2;
    # the contract (shared safe_norm and compensated sum) makes the k = 1 run bitwise equal, which is asserted too.
    # k = 2 with the same first frame column obeys the same bound (lyap-benettin S4(a)).
    seed = gs_seed(10)
    function gfirst(sys, x0, tau)
        P = tangent(sys, x0)
        x = forward(sys, x0)
        for _ in 2:tau
            P = tangent(sys, x) * P
            x = forward(sys, x)
        end
        return Float64(sum(abs2, P))
    end
    for (T, Ks, m, n, taus) in ((Float64, gs_KS(), 50, 10^4, (1, 4)),), K in Ks, tau in taus
        sys = gs_stdmap(T, K)
        x0s = gs_starts(T, seed + 1, m)
        rb = propagate(sys, x0s, n; k = 1, tau, seed, reducers = (Benettin(),))
        for alg in gs_algsT(T)
            r1 = gs_run(sys, x0s, n; k = 1, tau, seed, alg)
            r2 = gs_run(sys, x0s[1:min(m, 10)], n; k = 2, tau, seed, alg)
            for j in 1:m
                X = rb.finals[j][1]
                bound = 2 * gs_BS(T, n, gfirst(sys, x0s[j], tau), K) / n
                f1 = gs_fin(r1, j)
                @test Float64(abs(f1.lambda[1] - X.lambda[1])) <= bound
                @test isequal(f1.sum[1], X.sum[1]) && isequal(f1.comp[1], X.comp[1])
                j <= 10 && @test Float64(abs(gs_fin(r2, j).lambda[1] - X.lambda[1])) <= bound
            end
        end
    end
end

@testset "AC11 volume and sum identities" begin
    # G7: sum_{i<=p} lambda_i(n) = (1/n) ln vol_p(Y_n Q_{0,p}); reference: A_S^n in BigInt, the frame as exact rationals,
    # the Gram determinant exact, the logarithm at 256 bits. Pass: <= p theta_sum (theta from cond_2(A_S^tau), API Q5) and
    # |sum_i lambda_i| <= theta_sum. The sandwich G8 is stage 04.
    seed = gs_seed(11)
    n = 60
    for (T, m) in ((Float64, 10),)
        A = gs_AS(T)
        Aint = gs_intmat(A)
        sys = LinearSymplecticMap(A)
        Q0s = [initial_deviations(seed, j, Val(4), Val(4), T) for j in 1:m]
        refs = [[gs_exact_logvol(Aint, n, Q0s[j], p) / n for p in 1:4] for j in 1:m]   # independent of variant and tau
        for alg in gs_algsT(T), tau in (1, 2, 3)
            res, err = gs_try_run(sys, [gs_zero(T, 4) for _ in 1:m], n; k = 4, tau, seed, alg)
            if err !== nothing    # the documented CGS failure at large G_w (AC6c) is reported, not scored
                @test alg isa CGS && tau >= 2 && err.what == "orthogonality"
                continue
            end
            th = gs_theta_lin(T, A, tau)
            for j in 1:m
                f = gs_fin(res, j)
                for p in 1:4
                    @test gs_absdiff(sum(f.lambda[1:p]), refs[j][p]) <= p * th
                end
                @test Float64(abs(sum(f.lambda))) <= th
            end
        end
    end
end

@testset "AC14 provenance" begin
    # exact equality of the record: estimator, orthonormalizer, orthogonality_tol (value in force), k, tau, precision,
    # backend, device, fallback reason, seed; no Householder keys for the Gram-Schmidt variants
    seed = gs_seed(14)
    for alg in gs_algs(), tol in (8.0, 4.0, Inf)
        @test reducer_provenance(QRSpectrum(alg; orthogonality_tol = tol)) ==
              (estimator = "qr", orthonormalizer = gs_name(alg), orthogonality_tol = tol)
    end
    @test QRSpectrum(RGS()).orthogonality_tol === 8.0                    # the key value, passed in by Experiments
    function check(T, alg, tol, tau, extra)
        r = QRSpectrum(alg; orthogonality_tol = tol)
        sys = gs_stdmap(T, 2.0)
        res = propagate(sys, gs_starts(T, seed + 1, 3), 12; k = 2, tau, seed, reducers = (r, extra...))
        p = res.provenance
        @test p.estimator == "qr" && p.orthonormalizer == gs_name(alg) && p.orthogonality_tol === tol
        @test p.k == 2 && p.tau == tau && p.seed == seed && p.precision == precision_label(T)
        @test p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == ""
        @test !haskey(p, :iteration_path) && !haskey(p, :implementation) && !haskey(p, :d)
    end
    for alg in gs_algs(), (tol, tau) in ((8.0, 1), (4.0, 3), (Inf, 3))
        check(Float64, alg, tol, tau, ())
    end
    check(Float64, RGS(), 8.0, 1, (ResidualMonitor(gs_stdmap(Float64, 2.0), 4.0, Inf),))   # a chained monitor adds no key
end
