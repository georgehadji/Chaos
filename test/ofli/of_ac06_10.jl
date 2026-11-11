# AC6 to AC10: growth laws on Henon-Heiles orbits (quasi-periodic, chaotic), initial-basis independence, the basis-free sandwich, identities and chains.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).
# Sizes: D38 package sizes with the full sizes of the spec behind OF_FULL (the stage 04 script sets OF_FULL_RUN = true).

# A series of the OFLI reducer on the one-step map im, k columns, handoff every `tau` steps; returns the snapshots (handoff j is step j tau).
of_hh_snaps(im, z0, W0, N, tau) = of_snaps(OFLI(HenonHeilesFlow()), im, z0, W0, N, tau)
# The full running maximum (s = 0 term included) of the basis form at the snapshot s.
of_fullB(g0, s) = max(maximum(g0), s.ofli)

@testset "AC6 quasi-periodic orbit R and the regular class: log OFLI_B grows like ln t" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); im = of_im(Float64); Tmax = of_size(1.0e4, 1.0e5); W0 = SMatrix{4, 4, Float64}(I)
    function slope(z0, W)
        snaps = of_hh_snaps(im, z0, W, round(Int, Tmax * 20), 20)
        ts = collect(1.0e3:1.0:Tmax)
        return of_slope(log.(ts), [snaps[round(Int, t)].ofli for t in ts])                   # t = handoff number (20 steps of 1/20)
    end
    sR = slope(of_orbit_R(), W0)
    @test 0.7 <= sR <= 1.3                                                                  # the interval is wide because the constant still moves (O7)
    reg = of_classes(OF_NCLASS, 1).regular[1:of_size(6, 20)]
    slopes = [slope(z, of_frame(of_seed(6, 1), j, 4, 4, Float64)) for (j, z) in enumerate(reg)]
    inside = count(s -> 0.7 <= s <= 1.3, slopes)
    @test inside >= 0.8 * length(slopes)
    @info "AC6 reported" slope_R = sR slopes_regular = slopes sali_regular = [of_sali(z, 1000.0) for z in reg]
end

@testset "AC7 chaotic class: the slope of OFLI_B against t equals the mean finite-time exponent" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); im = of_im(Float64); m = OF_NCLASS; Tw = of_size(1000.0, 2000.0); n = round(Int, Tw * 20); h = OF_H
    starts = of_classes(m, 1).chaotic; seed = of_seed(7, 1)
    slopes = Float64[]; lams = Float64[]; nviol = 0; dfo = Float64[]
    ts = collect(200.0:1.0:Tw)
    for (j, z0) in enumerate(starts)
        W = of_frame(seed, j, 4, 4, Float64)
        snaps = of_hh_snaps(im, z0, W, n, 20)
        push!(slopes, of_slope(ts, [snaps[round(Int, t)].ofli for t in ts]))
        _, _, st = propagate_orbit(im, z0, of_frame(seed, j, 4, 1, Float64), n, 20, (Benettin(),))
        push!(lams, reducer_final(Benettin(), st[1]).lambda[1] / h)                      # exponents of the estimator are per step
        # OFLI_B <= FLI_B and OFLI_1 <= FLI_1 at every handoff, with the same S_i (the slack is the rounding of one sum and one logarithm)
        for s in snaps
            (s.ofli <= s.fli + 4 * of_eps(Float64) * max(1.0, abs(s.fli)) && s.ofli1 <= s.fli1 + 4 * of_eps(Float64) * max(1.0, abs(s.fli1))) || (nviol += 1)
        end
        push!(dfo, last(snaps).fli - last(snaps).ofli)
    end
    @test nviol == 0
    @test all(>(0), lams) && all(>(0), slopes)
    @test of_agree(slopes, lams; allowance = 4 / (Tw - 200))                               # 4 combined standard errors plus the finite-window allowance (L1 S5, proposal)
    @info "AC7 reported (R2)" median_FLI_minus_OFLI = of_median(dfo) mean_slope = of_mean(slopes) mean_lambda = of_mean(lams)
end

@testset "AC8 initial-basis independence: any two orthonormal frames agree within ln sqrt(4)" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); im = of_im(Float64); times = of_size((100.0, 1000.0), (1000.0, 10000.0)); Nmax = round(Int, last(times) * 20)
    for (name, z0) in ((:R, of_orbit_R()), (:C, of_orbit_C()))
        pts = of_pts(im, z0, Nmax); Ms = of_tangents(im, pts, Nmax); cache = Dict{Int, Vector{Float64}}()
        rng = StableRNG(of_seed(8, name === :R ? 1 : 2))
        frames = [SMatrix{4, 4, Float64}(I), of_rotated(rng, 4), of_rotated(rng, 4), of_rotated(rng, 4)]
        for W in frames, i in 1:4
            @test check_frame(OFLI(flow), W[:, i:i], z0) <= 8 * of_eps(Float64)           # unit columns (the single-vector form also needs a non-parallel start)
        end
        runs = [(W = W, snaps = of_hh_snaps(im, z0, W, Nmax, 20), g0 = ofli_initial(flow, z0, W)) for W in frames]
        for T in times
            idx = round(Int, T); vals = Float64[]; bnds = Float64[]
            for r in runs
                s = r.snaps[idx]; i0 = argmax(r.g0)
                push!(vals, of_fullB(r.g0, s))
                if r.g0[i0] > s.ofli
                    push!(bnds, of_bound_cheap(Ms, r.W[:, i0], 0, exp(r.g0[i0]), Float64, cache))
                else
                    push!(bnds, of_bound_cheap(Ms, r.W[:, s.colB], s.stepB, s.sinB, Float64, cache))
                end
            end
            spread = maximum(vals) - minimum(vals)
            @test spread <= log(sqrt(4)) + 2 * maximum(bnds)
            @info "AC8 reported" orbit = name t = T spread = spread bound = log(2.0)
        end
    end
end

# The largest singular value of P_x Y, ln, in Float64 from the BigFloat product (the scale is divided out in BigFloat first).
function of_lnsigma1(flow, Y, x)
    return setprecision(BigFloat, OF_BITS) do
        d = size(Y, 1); f = vector_field(flow, SVector{d, BigFloat}(map(of_big, x))); P = one(SMatrix{d, d, BigFloat}) - (f * f') / dot(f, f)
        A = P * Y; s = norm(A)
        Float64(log(BigFloat(maximum(svdvals(SMatrix{d, d, Float64}(map(Float64, A / s)))))) + log(s))
    end
end

@testset "AC9 basis-free sandwich: F* - ln sqrt(4) <= OFLI_B <= F*" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); N = of_size(1000, 10_000); m = of_size(3, 20); tau = 10
    keys_ = vcat([:R, :C], [(:reg, i) for i in 1:m], [(:cha, i) for i in 1:m])
    nbad = 0; fr = Float64[]
    for key in keys_, canonical in (true, false)
        c = of_case(key, Float64, :o8, N; frame = canonical ? Matrix{Float64}(I, 4, 4) : nothing)
        snaps = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau); g0 = ofli_initial(flow, c.x0, c.W0)
        Fstar = 0.0                                                                    # the s = 0 term: sigma_1(P_{x_0}) = 1
        for (kk, s) in enumerate(snaps)
            n = kk * tau
            Fstar = max(Fstar, of_lnsigma1(flow, c.Ys[n + 1], c.pts[n + 1]))
            n % 100 == 0 || n == N || continue
            G = of_fullB(g0, s); i0 = argmax(g0)
            B = g0[i0] > s.ofli ? of_bound(c, Float64, i0, 0) : of_bound(c, Float64, s.colB, s.stepB)
            (Fstar - log(sqrt(4.0)) - B <= G <= Fstar + B) || (nbad += 1)
            n == N && push!(fr, (Fstar - G) / log(sqrt(4.0)))
        end
    end
    @test nbad == 0
    @test all(x -> -1.0e-9 <= x <= 1 + 1.0e-9, fr)
    @info "AC9 reported: (F* - OFLI_B) / ln sqrt(4) at the last step" quantiles = (minimum(fr), of_median(fr), maximum(fr))
end

@testset "AC10a chains: the sums are bitwise those of FLI and SALI, the companions are the delivered FLI" begin
    _ = (OFLI,)
    flow = HenonHeilesFlow(); N = of_size(1000, 10_000); tau = 10; m = of_size(2, 20)
    keys_ = vcat([:R, :C], [(:reg, i) for i in 1:m], [(:cha, i) for i in 1:m])
    for key in keys_
        c = of_case(key, Float64, :o8, 10; frame = nothing)
        W2 = SMatrix{4, 2, Float64}(c.W0[:, 1:2]); im = c.im; x0 = c.x0
        fli = of_chain((FLI(),), im, x0, W2, N, tau)[1]
        sali = of_chain((SALI(),), im, x0, W2, N, tau)[1]
        o_alone = of_chain((OFLI(flow),), im, x0, W2, N, tau)[1]
        o_fli, f_last = of_chain((OFLI(flow; replace = false), FLI()), im, x0, W2, N, tau)
        o_sali, s_last = of_chain((OFLI(flow; replace = false), SALI()), im, x0, W2, N, tau)
        f_nr, o_last = of_chain((FLI(replace = false), OFLI(flow)), im, x0, W2, N, tau)
        @test all(i -> o_alone[i].sum == fli[i].sum && o_alone[i].comp == fli[i].comp, eachindex(fli))
        @test all(i -> o_alone[i].sum == sali[i].sum && o_alone[i].comp == sali[i].comp, eachindex(sali))   # SALI and FLI sums: the same operations
        @test all(i -> o_fli[i].sum == fli[i].sum && o_fli[i].comp == fli[i].comp && f_last[i].fli == fli[i].fli, eachindex(fli))
        @test all(i -> o_sali[i].sum == sali[i].sum && o_sali[i].comp == sali[i].comp && s_last[i].sali == sali[i].sali, eachindex(sali))
        @test all(i -> o_alone[i].fli == fli[i].fli && o_alone[i].fli1 == fli[i].fli1, eachindex(fli))           # companions, k = 2
        @test all(i -> o_last[i].sum == fli[i].sum && f_nr[i].fli == fli[i].fli && o_last[i].ofli == o_alone[i].ofli && o_last[i].ofli1 == o_alone[i].ofli1, eachindex(fli))
        @test all(i -> o_fli[i].ofli == o_alone[i].ofli && o_sali[i].ofli == o_alone[i].ofli, eachindex(fli))   # the position of a non-replacing OFLI changes nothing
    end
end

@testset "AC10b identities: instantaneous value, sin theta at the maximiser, OFLI <= FLI, the s = 0 floor" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); N = of_size(1000, 10_000); tau = 10; m = of_size(2, 20); eT = of_eps(Float64)
    keys_ = vcat([:R, :C], [(:reg, i) for i in 1:m], [(:cha, i) for i in 1:m])
    for key in keys_
        c = of_case(key, Float64, :o8, 10; frame = nothing)
        for W in (c.W0, SMatrix{4, 4, Float64}(I))
            _, _, st = propagate_orbit(c.im, c.x0, W, N, tau, (of_Trace(), of_Probe(OFLI(flow))))
            tr = Chaoslab.reducer_final(of_Trace(), st[1]); snaps = Chaoslab.reducer_final(of_Probe(OFLI(flow)), st[2]).snaps
            worst = 0.0; Sold = zeros(4); nbad_le = 0; sinB_bad = 0
            sinth(n, i) = (f = vector_field(flow, tr.xs[n]); w = tr.Ws[n][:, i]; norm(w - (dot(f, w) / dot(f, f)) * f) / norm(w))
            for (n, s) in enumerate(snaps)
                gs = [Sold[i] + log(sinth(n, i) * norm(tr.Ws[n][:, i])) for i in 1:4]
                worst = max(worst, abs(s.inst - maximum(gs)) / (4 * eT * max(1.0, abs(s.inst))))     # (b) g_i = S_i(old) + ln ||P w_arr||
                Sold = [s.sum[i] + s.comp[i] for i in 1:4]
                (s.ofli <= s.fli + 4 * eT * max(1.0, abs(s.fli)) && s.ofli1 <= s.fli1 + 4 * eT * max(1.0, abs(s.fli1))) || (nbad_le += 1)   # (c)
            end
            @test worst <= 1
            @test nbad_le == 0
            s = last(snaps)                                                                  # sin theta recorded at the maximiser
            @test abs(s.sinB - sinth(s.stepB ÷ tau, s.colB)) <= 4 * eT && abs(s.sin1 - sinth(s.step1 ÷ tau, 1)) <= 4 * eT
            @test s.sinmin <= s.sinB && s.nminf == 0
            # (d) the s = 0 floor: the squared orthogonal norms of an orthonormal full frame sum to d - 1 = 3, so the largest term is at least 3/4
            g0 = ofli_initial(flow, c.x0, W)
            @test abs(sum(exp.(2 .* g0)) - 3) <= 8 * eT && maximum(g0) >= 0.5 * log(3 / 4) - 8 * eT
            @test all(s -> max(maximum(g0), s.ofli) >= 0.5 * log(3 / 4) - 8 * eT, snaps)
        end
    end
end
