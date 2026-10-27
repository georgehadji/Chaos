# AC4 to AC6: the algebra against the independent 256-bit reference, the index (WB15, I5), the reducer and the pure function (I6).

@testset "AC4 algebra against the reference" begin
    _ = (WeightedBirkhoff, birkhoff_index, CosObservable, SinObservable, LiftedObservable)
    # K_std in {0.5 regular class, 2 pi, 10}, the five observables, p in {1, 2}, N in {256, 1000, 4096} (three horizons of one pass), both types: each window average within c_p eps S of the independent
    # BigFloat(256) sums (WB14) at the same realized samples (weights from the definition, observables evaluated in BigFloat at the widened wrapped states), Delta within 2 c_p eps S; where
    # |Delta_ref| >= 1e3 c_p eps S, absdig within 1e-2 of the reference absdig; every flag equals its definition from the module's own Delta and S. D38: 6 starts per class (spec: 20).
    Ns = (256, 1000, 4096)
    nstart = WB_FULL ? 20 : 6
    for T in WB_TYPES
        epsT = wb_eps(T)
        bad = String[]
        worst = Dict((p, what) => 0.0 for p in (1, 2), what in (:avg, :delta, :S))
        ndig = 0
        for (cl, sub) in ((:half, 0), (:two_pi, 1), (:ten, 3))
            c = wb_case(cl, T, nstart, sub)
            rs = (WeightedBirkhoff(; window = Ns, observables = wb_std5(), birkhoff_order = 1), WeightedBirkhoff(; window = Ns, observables = wb_std5(), birkhoff_order = 2))
            nmax = 2 * maximum(Ns)
            for j in 1:nstart
                x0 = c.x0[j]
                _, _, fin = wb_run((rs..., wb_Norm()), c.sys, x0, wb_frame(wb_seed(4, 1), c.ids[j], T), nmax, 1)
                fb = wb_fref_std(wb_states(c.sys, x0, nmax))
                cols = [fb[:, k] for k in 1:5]
                for (ir, p) in ((1, 1), (2, 2)), (h, N) in enumerate(Ns), k in 1:5
                    f = fin[ir]
                    ref = wb_ref_avgs(cols[k], N, p)
                    S = Float64(ref.S)
                    cp = wb_cp(p)
                    tag = "T=$T $cl start $j p=$p N=$N obs $k"
                    ea1, ea2 = wb_err(f.avg1[k, h], ref.avg1), wb_err(f.avg2[k, h], ref.avg2)
                    ed = wb_err(f.delta[k, h], ref.avg1 - ref.avg2)
                    (ea1 <= cp * epsT * S && ea2 <= cp * epsT * S) || push!(bad, tag * " average $(max(ea1, ea2) / (epsT * S)) eps S")
                    ed <= 2 * cp * epsT * S || push!(bad, tag * " Delta $(ed / (epsT * S)) eps S")
                    wb_err(f.S[k, h], ref.S) <= 4 * epsT * S || push!(bad, tag * " S")
                    worst[(p, :avg)] = max(worst[(p, :avg)], max(ea1, ea2) / (cp * epsT * S) * cp)
                    worst[(p, :delta)] = max(worst[(p, :delta)], ed / (epsT * S))
                    ix = wb_ref_index(ref.avg1, ref.avg2, ref.S, epsT)
                    if Float64(abs(ref.avg1 - ref.avg2)) >= 1e3 * cp * epsT * S
                        ndig += 1
                        wb_err(f.absdig[k, h], ix.absdig) <= 1e-2 || push!(bad, tag * " absdig")
                    end
                    (f.below_floor[k, h] == (abs(f.delta[k, h]) <= T(WB_MULT) * effective_eps(T) * f.S[k, h]) && f.zero[k, h] == (f.S[k, h] == 0)) || push!(bad, tag * " flags")
                end
            end
        end
        @test wb_report(bad, "AC4 $T")
        @test ndig > 0
        @info "AC4 $T: worst error in eps S (average, Delta) for p = 1 and p = 2 (spec dry run, Float64: 0.6 and 7.2, 0.7 and 8.2); entries with the digit comparison" avg_p1 = worst[(1, :avg)] avg_p2 = worst[(2, :avg)] delta_p1 = worst[(1, :delta)] delta_p2 = worst[(2, :delta)] ndig
    end
end

# The digits expected from the returned Delta, S, avg_1, avg_2 by WB15, in BigFloat. `got` is the NamedTuple (absdig, reldig, dig, below_floor, zero) of one entry of the module, `d` its Delta.
# Digit values below 1 cannot be matched to a relative 2 eps in T (the logarithm of an argument near 1 loses it), so the gate is 2 eps max(|value|, 1) (plan, Items for the owner).
function wb_check_index!(bad, tag, ::Type{T}, a1, a2, S, d, got) where {T}
    epsT = wb_eps(T)
    setprecision(BigFloat, WB_BITS) do
        Sb, a1b, a2b, db = wb_big(S), wb_big(a1), wb_big(a2), wb_big(d)
        phi = WB_MULT * BigFloat(epsT) * Sb
        if S == 0
            (isnan(got.absdig) && isnan(got.reldig) && isnan(got.dig) && got.zero) || push!(bad, tag * " zero case")
            return
        end
        got.zero && push!(bad, tag * " zero flag with S > 0")
        mx = max(abs(db), phi)
        m = (abs(a1b) + abs(a2b)) / 2
        ea = -log10(mx)
        near(x, e) = abs(wb_big(x) - e) <= 2 * epsT * max(abs(e), 1)
        near(got.absdig, ea) || push!(bad, tag * " absdig")
        if m > phi
            er = -log10(mx / m)
            (!isnan(got.reldig) && near(got.reldig, er) && near(got.dig, max(ea, er))) || push!(bad, tag * " reldig/dig (defined)")
        else
            (isnan(got.reldig) && near(got.dig, ea)) || push!(bad, tag * " reldig must be undefined where m <= phi")
        end
        got.below_floor == (abs(db) <= phi) || push!(bad, tag * " flag")
    end
    return nothing
end

@testset "AC5 index" begin
    _ = (WeightedBirkhoff, birkhoff_index, weighted_birkhoff, CosObservable, SinObservable, LiftedObservable)
    # (a) absdig, reldig, dig equal the values formed from the returned Delta, S, avg_1, avg_2 by WB15; (b) reldig is undefined exactly where m <= phi; (c) S = 0: digits undefined, flag "zero",
    # Delta = 0; (d) Delta is returned unclipped in every case, flagged ones included; (e) before the pair completes avg_2, Delta and the digits are NaN and `complete` is false.
    # Cases: reducer finals on K_std = 0.5 regular class, 2 pi, 10 and on the K_std = 0 orbit (digits at the floor) with the five observables; the pure function on exact rotation signals
    # (m far below phi: reldig undefined) and on constants.
    N = 500
    for T in WB_TYPES
        epsT = wb_eps(T)
        bad = String[]
        ncase = Dict(:reldig_undefined => 0, :reldig_defined => 0, :flagged_nonzero => 0, :flagged => 0, :unflagged => 0)
        function consume!(tag, a1, a2, S, d, ad, rd, dg, bf, zf)
            wb_check_index!(bad, tag, T, a1, a2, S, d, (absdig = ad, reldig = rd, dig = dg, below_floor = bf, zero = zf))
            isequal(d, a1 - a2) || push!(bad, tag * " Delta is not avg_1 - avg_2 (clipped?)")        # (d)
            S == 0 && return
            ncase[isnan(rd) ? :reldig_undefined : :reldig_defined] += 1
            ncase[bf ? :flagged : :unflagged] += 1
            bf && d != 0 && (ncase[:flagged_nonzero] += 1)
        end
        for (cl, sub) in ((:half, 0), (:two_pi, 1), (:ten, 3), (:zero, 4))
            c = wb_case(cl, T, 3, sub)
            red = WeightedBirkhoff(; window = (N ÷ 2, N), observables = wb_std5())
            for j in 1:3
                _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(5, 1), j, T), 2N, 1)
                f = fin[1]
                for k in 1:5, h in 1:2
                    consume!("T=$T $cl start $j obs $k horizon $h", f.avg1[k, h], f.avg2[k, h], f.S[k, h], f.delta[k, h], f.absdig[k, h], f.reldig[k, h], f.dig[k, h], f.below_floor[k, h], f.zero[k, h])
                end
            end
        end
        alpha = setprecision(BigFloat, WB_BITS) do
            sqrt(BigFloat(2)) - 1
        end
        cs, sn = wb_rot_big(wb_dec("0.3"), alpha, 1, 2N)
        # a spike at the sample 0 (weight zero, so it sets S = 1 only) over noise of 1/8 eps: the averages are far below the floor, Delta is nonzero and flagged, m <= phi (reldig undefined)
        noise = T[(j == 0 ? one(T) : zero(T)) + T(0.125) * effective_eps(T) * cos(T(j)) for j in 0:(2N - 1)]
        ser = hcat(wb_round_vec(T, cs), wb_round_vec(T, sn), fill(T(0.3), 2N), wb_round_vec(T, cs) .+ T(1e-9), noise)
        r = weighted_birkhoff(ser, N)
        ix = wb_idx(r)
        for k in 1:5
            consume!("T=$T function column $k", r.avg1[k], r.avg2[k], r.S[k], ix.delta[k], ix.absdig[k], ix.reldig[k], ix.dig[k], ix.below_floor[k], ix.zero[k])
        end
        @test wb_report(bad, "AC5 (a), (b), (d) $T")
        @test ncase[:reldig_undefined] > 0 && ncase[:reldig_defined] > 0 && ncase[:flagged_nonzero] > 0 && ncase[:unflagged] > 0
        @info "AC5 $T: cases by kind (premises of the test)" ncase
        # (c) S = 0: an identically zero observable, and sin q on the orbit started at the fixed point (0, 0) (the standard map keeps it there exactly)
        sys = WrappedStandardMap(wb_K(T, :half))
        x0 = torus_state(SVector{2, T}(0, 0))
        red = WeightedBirkhoff(; window = 16, observables = (wb_Const(0.0), SinObservable(1), CosObservable(1)))
        _, _, fin = wb_run((red, wb_Norm()), sys, x0, wb_frame(wb_seed(5, 2), 1, T), 32, 1)
        f = fin[1]
        for k in 1:2
            @test f.S[k] == 0 && f.zero[k] && f.delta[k] == 0 && isnan(f.absdig[k]) && isnan(f.reldig[k]) && isnan(f.dig[k]) && f.avg1[k] == 0 && f.avg2[k] == 0
        end
        @test f.S[3] == 1 && !f.zero[3] && f.below_floor[3] && f.delta[3] == 0                 # cos 0 = 1: a constant, at the floor, not zero
        z = weighted_birkhoff(zeros(T, 32), 16)
        iz = birkhoff_index(z.avg1[1], z.avg2[1], z.S[1])
        @test iz.zero && isnan(iz.absdig) && isnan(iz.reldig) && isnan(iz.dig) && iz.delta == 0
        # (e) before the pair completes: window 100, K_std = 0.5 orbit. At 60 handoffs avg_1 is not yet defined; at 150 avg_1 is, avg_2, Delta, digits are NaN and complete is false; at 200 all defined.
        c = wb_case(:half, T, 1, 0)
        red = WeightedBirkhoff(; window = 100, observables = wb_std5())
        for (n, a1def, a2def) in ((60, false, false), (150, true, false), (200, true, true))
            _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[1], wb_frame(wb_seed(5, 3), 1, T), n, 1)
            f = fin[1]
            @test f.handoffs == n && f.complete == a2def
            @test all(isnan, f.avg1) == !a1def && all(!isnan, f.avg1) == a1def
            for fld in (f.avg2, f.delta, f.absdig, f.dig)
                @test all(isnan, fld) == !a2def && all(!isnan, fld) == a2def
            end
        end
        # two horizons: at 150 handoffs the horizon 50 is complete and the horizon 100 is not; `complete` waits for the largest
        red2 = WeightedBirkhoff(; window = (50, 100), observables = wb_std5())
        _, _, fin = wb_run((red2, wb_Norm()), c.sys, c.x0[1], wb_frame(wb_seed(5, 3), 1, T), 150, 1)
        f = fin[1]
        @test !f.complete && all(!isnan, f.delta[:, 1]) && all(isnan, f.delta[:, 2]) && all(!isnan, f.avg1[:, 2])
    end
end

@testset "AC6 reducer and function" begin
    _ = (WeightedBirkhoff, weighted_birkhoff, CosObservable, SinObservable, LiftedObservable, check_run)
    # K = 4 observables, H = 3 horizons (64, 256, 1000), p in {1, 2}, tau in {1, 4}, both types, K_std in {0.5 regular class, 2 pi, 10}: the averages of both windows (and S, B) of the reducer are bitwise
    # those of the pure function on the series sampled at the handoffs; and the same for a weight vector given to the function equal to the module's weights (the unnormalised sampled weights of WB14).
    # D38: 3 starts per class (spec: 20).
    Ns = (64, 256, 1000)
    nstart = WB_FULL ? 20 : 3
    for T in WB_TYPES
        bad = String[]
        for (cl, sub) in ((:half, 0), (:two_pi, 1), (:ten, 3)), p in (1, 2), tau in (1, 4)
            c = wb_case(cl, T, nstart, sub)
            red = WeightedBirkhoff(; window = Ns, observables = wb_std4(), birkhoff_order = p)
            check_run(red, tau, 2000 * tau)
            for j in 1:nstart
                x0 = c.x0[j]
                _, _, fin = wb_run((red, wb_Norm()), c.sys, x0, wb_frame(wb_seed(6, 1), c.ids[j], T), 2000 * tau, tau)
                f = fin[1]
                ser = wb_series(wb_std4(), wb_states_tau(c.sys, x0, 2000, tau))
                g = weighted_birkhoff(ser, Ns; order = p)
                tag = "T=$T $cl p=$p tau=$tau start $j"
                (isequal(Matrix(f.avg1), g.avg1) && isequal(Matrix(f.avg2), g.avg2) && isequal(Matrix(f.S), g.S) && isequal(Vector(f.B), g.B)) || push!(bad, tag * " reducer vs function")
                for (h, N) in enumerate(Ns), k in 1:4
                    w = wb_wvec(T, N, p)
                    gw = weighted_birkhoff(ser[:, k], N; order = p, weights = w)
                    (isequal(gw.avg1[1], g.avg1[k, h]) && isequal(gw.avg2[1], g.avg2[k, h]) && isequal(gw.S[1], g.S[k, h]) && isequal(gw.B[1], g.B[h])) || push!(bad, tag * " weight vector N=$N obs $k")
                end
            end
        end
        @test wb_report(bad, "AC6 $T")
    end
end
