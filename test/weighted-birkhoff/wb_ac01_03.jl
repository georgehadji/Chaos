# AC1 to AC3: constants, covariances, normalisation (I1 to I3), the exact rotation (WB16) and the K_std = 0 standard map (WB16).
# The pure function weighted_birkhoff(series, window; order, weights) and the reducer WeightedBirkhoff are both exercised; birkhoff_index forms Delta and the digits from an average pair.
# Every testset begins by naming the new symbols it uses, so a missing Chaoslab symbol errors that testset only (UndefVarError).

@testset "AC1a constants" begin
    _ = (weighted_birkhoff, birkhoff_index, WeightedBirkhoff, check_run)
    # f = c, c in {0.3, 1e3, 1e-3, pi}, N in {16, 1000, 10^5} (Double64: N <= 10^4), K = 2 (c and 2c), H = 2 (windows N and 2N): both averages equal c within 2 eps |c|; Delta = 0 bitwise; the flag
    # is set (and "zero" is not); absdig = -log10(32 eps S) within 2 eps relative (S = |c|). The pure function on the constant series and the reducer on a constant observable of the K_std = 0 orbit.
    for T in WB_TYPES
        epsT = wb_eps(T)
        bad = String[]
        for c in (0.3, 1e3, 1e-3, pi), N in (16, 1000, 100_000)
            T === Double64 && N > 10_000 && continue
            cT = T(c)
            c2T = 2 * cT
            ser = [fill(cT, 4N) fill(c2T, 4N)]
            res = weighted_birkhoff(ser, (N, 2N))
            sys = WrappedStandardMap(wb_K(T, :zero))
            x0 = torus_state(SVector{2, T}(T(1.0), T(2.0)))
            red = WeightedBirkhoff(; window = (N, 2N), observables = (wb_Const(cT), wb_Const(c2T)))
            check_run(red, 1, 4N; x0 = x0)
            _, _, fin = wb_run((red, wb_Norm()), sys, x0, wb_frame(wb_seed(1, 1), 1, T), 4N, 1)
            fr = fin[1]
            for (name, a1, a2, S, ix) in (("function", res.avg1, res.avg2, res.S, wb_idx(res)),
                                          ("reducer", fr.avg1, fr.avg2, fr.S, (delta = fr.delta, absdig = fr.absdig, below_floor = fr.below_floor, zero = fr.zero)))
                for (k, cc) in enumerate((cT, c2T)), h in 1:2
                    tag = "T=$T c=$c N=$N $name obs $k horizon $h"
                    cb = wb_big(cc)
                    (wb_err(a1[k, h], cb) <= 2 * epsT * abs(Float64(cc)) && wb_err(a2[k, h], cb) <= 2 * epsT * abs(Float64(cc))) || push!(bad, tag * " average")
                    (ix.delta[k, h] == 0 && ix.below_floor[k, h] && !ix.zero[k, h]) || push!(bad, tag * " delta/flag")
                    S[k, h] == abs(cc) || push!(bad, tag * " S")
                    ref = wb_ref_index(cb, cb, abs(cb), epsT).absdig
                    setprecision(BigFloat, WB_BITS) do
                        abs(wb_big(ix.absdig[k, h]) - ref) <= 2 * epsT * ref || push!(bad, tag * " absdig")
                    end
                end
            end
        end
        @test wb_report(bad, "AC1a $T")
    end
end

@testset "AC1b covariances" begin
    _ = (weighted_birkhoff, birkhoff_index, CosObservable)
    # K_std = 0.5 regular class, cos q, N = 1000, w^[1] (c_1 = 4): scaling by 2^k, k in {-20, 7}, scales each average, Delta and S bitwise; scaling by 3.7 changes each average by at most
    # 2 c_p eps (3.7 S); adding b = 1e3 changes Delta by at most 2 c_p eps (S + |b|). D38: 4 starts (spec: 20); the series are the samples of the module's own CosObservable.
    N = 1000
    nstart = WB_FULL ? 20 : 4
    cp = wb_cp(1)
    for T in WB_TYPES
        epsT = wb_eps(T)
        c = wb_case(:half, T, nstart, 0)
        bad = String[]
        for j in 1:nstart
            sts = wb_states(c.sys, c.x0[j], 2N)
            f = wb_series((CosObservable(1),), sts)[:, 1]
            base = weighted_birkhoff(f, N)
            ib = birkhoff_index(base.avg1[1], base.avg2[1], base.S[1])
            S = Float64(base.S[1])
            for k in (-20, 7)
                a = T(ldexp(1.0, k))
                r = weighted_birkhoff(f .* a, N)
                ig = birkhoff_index(r.avg1[1], r.avg2[1], r.S[1])
                (isequal(r.avg1[1], base.avg1[1] * a) && isequal(r.avg2[1], base.avg2[1] * a) && isequal(r.S[1], base.S[1] * a) && isequal(ig.delta, ib.delta * a)) || push!(bad, "T=$T start $j scale 2^$k not bitwise")
            end
            a = T(3.7)
            r = weighted_birkhoff(f .* a, N)
            for (new, old) in ((r.avg1[1], base.avg1[1]), (r.avg2[1], base.avg2[1]))
                setprecision(BigFloat, WB_BITS) do
                    Float64(abs(wb_big(new) - wb_big(a) * wb_big(old))) <= 2 * cp * epsT * Float64(a) * S || push!(bad, "T=$T start $j scale 3.7: average changed too much")
                end
            end
            b = T(1e3)
            r = weighted_birkhoff(f .+ b, N)
            ig = birkhoff_index(r.avg1[1], r.avg2[1], r.S[1])
            setprecision(BigFloat, WB_BITS) do
                Float64(abs(wb_big(ig.delta) - wb_big(ib.delta))) <= 2 * cp * epsT * (S + Float64(b)) || push!(bad, "T=$T start $j offset 1e3: Delta changed too much")
            end
        end
        @test wb_report(bad, "AC1b $T")
    end
end

@testset "AC1c exchanging the windows" begin
    _ = (weighted_birkhoff, birkhoff_index, CosObservable, SinObservable)
    # exchanging the two windows of a series negates Delta bitwise and leaves absdig, reldig, dig (and S) bitwise; also for the index function applied to a swapped pair of averages.
    N = 1000
    for T in WB_TYPES
        c = wb_case(:two_pi, T, 3, 2)
        bad = String[]
        for j in 1:3
            sts = wb_states(c.sys, c.x0[j], 2N)
            f = wb_series((CosObservable(1), SinObservable(1)), sts)
            swapped = vcat(f[(N + 1):(2N), :], f[1:N, :])
            r1, r2 = weighted_birkhoff(f, N), weighted_birkhoff(swapped, N)
            (isequal(r1.avg1, r2.avg2) && isequal(r1.avg2, r2.avg1) && isequal(r1.S, r2.S) && isequal(r1.B, r2.B)) || push!(bad, "T=$T start $j averages")
            i1, i2 = wb_idx(r1), wb_idx(r2)
            (isequal(i2.delta, -i1.delta) && isequal(i1.absdig, i2.absdig) && isequal(i1.reldig, i2.reldig) && isequal(i1.dig, i2.dig)) || push!(bad, "T=$T start $j index")
            for k in 1:2
                a, b, S = r1.avg1[k], r1.avg2[k], r1.S[k]
                x, y = birkhoff_index(a, b, S), birkhoff_index(b, a, S)
                (isequal(y.delta, -x.delta) && isequal(x.absdig, y.absdig) && isequal(x.reldig, y.reldig) && isequal(x.dig, y.dig)) || push!(bad, "T=$T start $j index function $k")
            end
        end
        @test wb_report(bad, "AC1c $T")
    end
end

@testset "AC1d normalisation" begin
    _ = (weighted_birkhoff,)
    # I3: the weights sum to 1 within 4 eps and are bitwise symmetric (the verified normalised weights birkhoff_weights, NF I5, D81 item 10); the module's weight total B of a horizon is the sum of the
    # sampled weights, within c_p eps relative of the independent 256-bit sum of the definition (WB18), one B for both windows (the constant series gives bitwise equal averages).
    for T in WB_TYPES
        epsT = wb_eps(T)
        for p in (1, 2), N in (16, 1000)
            w = birkhoff_weights(T, N; p = p)
            s = setprecision(BigFloat, WB_BITS) do
                sum(wb_big(x) for x in w)
            end
            @test Float64(abs(s - 1)) <= 4 * epsT && all(j -> w[j + 1] == w[N - j + 1], 1:(N - 1))
            ser = [fill(T(0.7), 2N) T[sin(T(j)) for j in 1:(2N)]]
            r = weighted_birkhoff(ser, N; order = p)
            Bref = sum(wb_wref(N, p))
            @test wb_err(r.B[1], Bref) <= wb_cp(p) * epsT * Float64(Bref) && r.avg1[1, 1] == r.avg2[1, 1]
        end
    end
end

@testset "AC2 exact rotation" begin
    _ = (weighted_birkhoff, birkhoff_index)
    # f = cos 2 pi theta, sin 2 pi theta, cos 4 pi theta; theta_j = frac(theta_0 + j alpha), alpha = sqrt 2 - 1, theta_0 in {1/7, 0.3}; signals generated in BigFloat(256) and rounded once. Float64 at
    # N in {500, 1000, 4000} for the k = 1 observables and N in {1000, 4000} for cos 4 pi theta; Double64 at N in {1000, 4000, 16000} (k = 1) and {4000, 16000} (k = 2); w^[1] and w^[2]: both window
    # averages within c_p eps (the limit is 0), |Delta| <= 2 c_p eps. The exact average is below a tenth of the bound at these N (WB6, F6 a: the premise, checked against the independent reference; Double64 N = 1000 with w^[1] has 2e-32 = 0.03 eps_T).
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    th0s = (setprecision(BigFloat, WB_BITS) do
                BigFloat(1) / 7
            end, wb_dec("0.3"))
    for T in WB_TYPES
        epsT = wb_eps(T)
        plan = T === Float64 ? ((1, (500, 1000, 4000)), (2, (1000, 4000))) : ((1, (1000, 4000, 16000)), (2, (4000, 16000)))
        bad = String[]
        worst = 0.0
        for th0 in th0s, p in (1, 2), (k, Ns) in plan, N in Ns
            cs, sn = wb_rot_big(th0, alpha, k, 2N)
            for (name, fb) in (("cos", cs), ("sin", sn))
                k == 2 && name == "sin" && continue
                r = weighted_birkhoff(wb_round_vec(T, fb), N; order = p)
                ix = birkhoff_index(r.avg1[1], r.avg2[1], r.S[1])
                cp = wb_cp(p)
                tag = "T=$T theta0=$(Float64(th0)) p=$p k=$k $name N=$N"
                (abs(r.avg1[1]) <= cp * epsT && abs(r.avg2[1]) <= cp * epsT) || push!(bad, tag * " average $(Float64(abs(r.avg1[1])) / epsT) $(Float64(abs(r.avg2[1])) / epsT) eps")
                abs(ix.delta) <= 2 * cp * epsT || push!(bad, tag * " Delta")
                ref = wb_ref_avgs(fb, N, p)
                (abs(ref.avg1) < cp * epsT * 0.1 && abs(ref.avg2) < cp * epsT * 0.1) || push!(bad, tag * " the exact average is not far below the bound (test premise)")
                worst = max(worst, Float64(abs(r.avg1[1])) / epsT / cp, Float64(abs(r.avg2[1])) / epsT / cp)
            end
        end
        @test wb_report(bad, "AC2 $T")
        @info "AC2 $T: worst |average| / (c_p eps) (realized ratio for stage 04)" worst
    end
end

@testset "AC3 K_std = 0 standard map" begin
    _ = (WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable, check_run)
    # WB16: q_j = q_0 + j p_0 mod 2 pi, p constant. Float64, p_0 in {1.0, 2.0, 3.3, 5.0}, random q_0, N in {1000, 4000}: the cos q and sin q averages within 4 eps and |Delta| <= 8 eps; the lifted p average
    # equals p_0 within 2 eps p_0. p_0 = pi: the average of cos 2q equals cos 2 q_0 within 16 eps and |Delta| <= 8 eps. Double64: p_0 in {1.0, 2.0, 3.3, 5.0}, N in {4000, 16000}: cos q within 4 eps_T.
    # D38: 4 random q_0 in Float64, 2 in Double64 (spec: 20), Double64 N = 16000 for 2 values of p_0 only; the full set is WB_FULL_RUN.
    for T in WB_TYPES
        epsT = wb_eps(T)
        nq = WB_FULL ? 20 : (T === Float64 ? 4 : 2)
        rng = StableRNG(wb_seed(3, T === Float64 ? 1 : 2))
        q0s = [T(2pi * rand(rng)) for _ in 1:nq]
        plan = T === Float64 ? ((1000, 4000), (1.0, 2.0, 3.3, 5.0)) : ((4000, 16000), (1.0, 2.0, 3.3, 5.0))
        sys = WrappedStandardMap(wb_K(T, :zero))
        bad = String[]
        worst = 0.0
        for N in plan[1], p0 in plan[2], q0 in q0s
            (T === Double64 && N == 16000 && !WB_FULL && p0 in (2.0, 5.0)) && continue          # D38
            pT = T(p0)
            x0 = torus_state(SVector{2, T}(q0, pT))
            red = WeightedBirkhoff(; window = N, observables = (CosObservable(1), SinObservable(1), LiftedObservable(2, 4)))
            check_run(red, 1, 2N; x0 = x0)
            _, _, fin = wb_run((red, wb_Norm()), sys, x0, wb_frame(wb_seed(3, 5), 1, T), 2N, 1)
            f = fin[1]
            tag = "T=$T p0=$p0 N=$N q0=$(Float64(q0))"
            for k in 1:2
                (abs(f.avg1[k]) <= 4 * epsT && abs(f.avg2[k]) <= 4 * epsT && abs(f.delta[k]) <= 8 * epsT) || push!(bad, tag * " trig $k: $(Float64(f.avg1[k]) / epsT), $(Float64(f.avg2[k]) / epsT) eps")
                worst = max(worst, Float64(abs(f.avg1[k])) / epsT, Float64(abs(f.avg2[k])) / epsT)
            end
            if T === Float64
                (abs(f.avg1[3] - x0[2]) <= 2 * epsT * abs(x0[2]) && abs(f.avg2[3] - x0[2]) <= 2 * epsT * abs(x0[2])) || push!(bad, tag * " lifted p")
            end
        end
        if T === Float64
            for N in plan[1], q0 in q0s                                    # p_0 = pi: cos 2q is the constant cos 2 q_0
                x0 = torus_state(SVector{2, Float64}(q0, pi))
                red = WeightedBirkhoff(; window = N, observables = (CosObservable(1, 2),))
                _, _, fin = wb_run((red, wb_Norm()), sys, x0, wb_frame(wb_seed(3, 6), 1, Float64), 2N, 1)
                f = fin[1]
                ref = setprecision(BigFloat, WB_BITS) do
                    cos(2 * wb_big(q0))
                end
                tag = "p0=pi N=$N q0=$q0"
                (wb_err(f.avg1[1], ref) <= 16 * epsT && wb_err(f.avg2[1], ref) <= 16 * epsT && abs(f.delta[1]) <= 8 * epsT) || push!(bad, tag * " cos 2q: $(wb_err(f.avg1[1], ref) / epsT) eps")
            end
        end
        @test wb_report(bad, "AC3 $T")
        @info "AC3 $T: worst |average| of cos q, sin q in eps (realized, spec dry run 0.27 for Delta)" worst
    end
end
