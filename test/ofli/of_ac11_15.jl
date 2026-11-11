# AC11 to AC15: the direction candidates, the same-orbit reference (rounding bound B_G), the schedule, health and frames, the precision ladder.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).

# ------------------------------------------------------------------ AC11 helpers: candidate B (the direction transported by the numerical tangent, one extra column)
# Run of the module reducer (candidate A, the exact field at the numerical state) with a trace of the arriving columns, and the cumulative transported field
# column Y_n f(x_0) (candidate B, read test-side from a run without replacement). Returns per-handoff vectors: G_A (the module's handoff maxima), G_B (the same
# maximum with the direction removed replaced by u_B), the angle sine s(n) between u_A and u_B, the running bound max_m s(m) ||w(m)||, the cumulative norm.
function of_candidates(flow, im, z0, W0, N, tau)
    d = size(W0, 1)
    _, _, st = propagate_orbit(im, z0, W0, N, tau, (of_Trace(), of_Probe(OFLI(flow))))
    trA = Chaoslab.reducer_final(of_Trace(), st[1]); snaps = Chaoslab.reducer_final(of_Probe(OFLI(flow)), st[2]).snaps
    f0 = vector_field(flow, z0); Wf = SMatrix{d, 1, Float64}(f0 / norm(f0))
    _, _, stB = propagate_orbit(im, z0, Wf, N, tau, (of_Trace(),))
    trB = Chaoslab.reducer_final(of_Trace(), stB[1])
    nh = length(snaps); GA = zeros(nh); GB = zeros(nh); sn = zeros(nh); bound = zeros(nh); gA = -Inf; gB = -Inf; bnd = 0.0; Sold = zeros(size(W0, 2))
    for n in 1:nh
        f = vector_field(flow, trA.xs[n]); uA = f / norm(f); uB = trB.Ws[n][:, 1] / norm(trB.Ws[n][:, 1])
        sn[n] = norm(uB - dot(uA, uB) * uA)
        for i in axes(W0, 2)
            w = trA.Ws[n][:, i]; gB = max(gB, Sold[i] + log(norm(w - dot(uB, w) * uB)))
        end
        Sold = [snaps[n].sum[i] + snaps[n].comp[i] for i in axes(W0, 2)]
        bnd = max(bnd, sn[n] * exp(maximum(Sold)))                                     # s(m) times the largest cumulative column norm
        GA[n] = snaps[n].ofli; GB[n] = gB; bound[n] = bnd
    end
    return (GA = GA, GB = GB, sin = sn, bound = bound)
end

@testset "AC11a direction candidates: the rigorous bound |OFLI_A - OFLI_B| <= sin(angle) ||w||" begin
    _ = (OFLI,)
    cases = Any[]
    ow = HarmonicOscillator(SVector(1.3)); push!(cases, (:oscillator, ow, IntegratorMap(ForestRuth(Float64), ow, 0.05), SVector(1.0, 0.3), 2000))
    pe = Pendulum(1.0); push!(cases, (:pendulum, pe, IntegratorMap(ForestRuth(Float64), pe, 0.05), SVector(2.0, 0.0), 2000))
    hh = HenonHeilesFlow(); im = of_im(Float64)
    push!(cases, (:R, hh, im, of_orbit_R(), of_size(20_000, 40_000))); push!(cases, (:C, hh, im, of_orbit_C(), of_size(20_000, 40_000)))
    for (name, flow, sys, z0, N) in cases
        d = dimension(flow); W0 = d == 2 ? SMatrix{2, 2, Float64}(I) : of_frame(of_seed(11, 1), 1, 4, 4, Float64)
        r = of_candidates(flow, sys, z0, W0, N, 20)
        a, b = exp.(r.GA), exp.(r.GB)
        # relative slack 1e-9: the two sides carry the rounding of two different recursions (spec: within B_G, about 1e-13 here)
        @test all(i -> abs(a[i] - b[i]) <= r.bound[i] * (1 + 1.0e-9) + 1.0e-9 * max(a[i], b[i]), eachindex(a))
        t1 = findfirst(>(0.1), r.sin)
        @info "AC11b reported (R1)" system = name first_time_angle_above_0_1 = (t1 === nothing ? Inf : Float64(t1)) ratio_candB_over_candA_at_end = exp(r.GB[end] - r.GA[end]) angle_end = r.sin[end]
    end
end

@testset "AC11c direction candidates: the angle at fixed t for the order-4 and order-8 tables (reported)" begin
    _ = (OFLI,)
    flow = HenonHeilesFlow(); z0 = of_orbit_R(); W0 = of_frame(of_seed(11, 2), 1, 4, 4, Float64); T = 50.0
    for (key, p) in ((:o2, 2), (:o4, 4), (:o8, 8))
        s = map((1 / 20, 1 / 40)) do h
            im = of_im(Float64, key, h); n = round(Int, T / h)
            of_candidates(flow, im, z0, W0, n, n).sin[1]
        end
        @info "AC11c reported (A1 prediction 2^p)" table = p angle_h = s[1] angle_half = s[2] ratio = s[1] / s[2] predicted = 2^p
        @test all(isfinite, s)
    end
end

# ------------------------------------------------------------------ AC12: same-orbit reference
@testset "AC12 same-orbit reference: |G - G_ref| <= B_G for both forms, Float64 and Double64" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); Ns = of_size((100, 1000), (100, 1000, 10_000)); N = maximum(Ns); m = of_size(2, 20)
    keys_ = vcat([:R, :C], [(:reg, i) for i in 1:m], [(:cha, i) for i in 1:m])
    ratios = Float64[]; nmax = 0; nrows = 0
    for T in (Float64, Double64), table in (:o8, :o4), key in keys_, tau in (1, 20)
        c = of_case(key, T, table, N)
        snaps = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau); g0 = ofli_initial(flow, c.x0, c.W0)
        cps = sort(unique(vcat(of_checkpoints(N, tau), filter(n -> n % tau == 0, collect(Ns)))))
        rows = of_rows(c, T, tau, snaps, g0, cps)
        @test all(r -> r.err <= r.bound, rows)
        append!(ratios, [r.ratio for r in rows]); nmax += count(r -> r.maxdiff, rows); nrows += length(rows)
    end
    @info "AC12 reported (R3)" max_error_over_bound = maximum(ratios) median_ratio = of_median(ratios) rows = nrows maximiser_differs = nmax
end

@testset "AC12 negative control: a single-precision projection exceeds B_G at n = 1e3 (test-side double)" begin
    flow = HenonHeilesFlow(); N = 1000; tau = 20; c = of_case(:C, Float64, :o8, N); mx = of_refmax(c.ref.g, tau, N); h = N ÷ tau
    for (mode, bad) in ((:ok, false), (:f32, true))
        snaps = of_snaps(of_Double(flow, mode), c.im, c.x0, c.W0, N, tau); g0 = of_initial(:ok, flow, c.x0, c.W0)
        G = max(maximum(g0), snaps[h].ofli)
        err = setprecision(BigFloat, OF_BITS) do
            Float64(abs(of_big(G) - mx.GB[h + 1]))
        end
        B = of_bound(c, Float64, mx.colB[h + 1], mx.stepB[h + 1])
        @test bad ? err > B : err <= 1.0e-10                       # the plain-sum double is not held to B_G; it must only be far below the single-precision error
    end
end

@testset "AC13 schedule: OFLI_tau(n) <= OFLI_1(n) + 2 B_G" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); N = 1000; nbad = 0
    for key in (:R, :C)
        c = of_case(key, Float64, :o8, N); g0 = ofli_initial(flow, c.x0, c.W0); full = Dict{Int, Any}()
        for tau in (1, 2, 5, 20)
            full[tau] = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau)
        end
        function gb(tau, n)                                                         # (value, bound) of the full basis series at step n
            s = full[tau][n ÷ tau]; i0 = argmax(g0)
            g0[i0] > s.ofli ? (g0[i0], of_bound(c, Float64, i0, 0)) : (s.ofli, of_bound(c, Float64, s.colB, s.stepB))
        end
        for n in (20, 100, 200, 1000), tau in (2, 5, 20)
            n % tau == 0 || continue
            v1, b1 = gb(1, n); vt, bt = gb(tau, n)
            vt <= v1 + b1 + bt || (nbad += 1)
        end
        @info "AC13 reported: OFLI_1 - OFLI_tau at n = 1000" orbit = key diffs = [gb(1, 1000)[1] - gb(tau, 1000)[1] for tau in (2, 5, 20)]
    end
    @test nbad == 0
end

# ------------------------------------------------------------------ AC14: health, frames, chain, run record
@testset "AC14a vector field: an equilibrium start is a health error at the first handoff, with orbit id and step" begin
    _ = (OFLI,)
    for (flow, im, good) in ((HarmonicOscillator(SVector(1.3)), IntegratorMap(Leapfrog(Float64), HarmonicOscillator(SVector(1.3)), 0.05), SVector(1.0, 0.3)),
                             (HenonHeilesFlow(), of_im(Float64), of_orbit_R()))
        d = dimension(flow); zero_start = zero(good)
        e = of_health_error(() -> propagate(im, [good, zero_start], 40; k = d, tau = 20, seed = 1, reducers = (OFLI(flow),)))
        @test e isa NumericalHealthError && e.what == "vector field" && e.orbit == 2 && e.step == 20
        e1 = of_health_error(() -> propagate_orbit(im, zero_start, SMatrix{d, d, Float64}(I), 40, 20, (OFLI(flow),); orbit = 9))
        @test e1 isa NumericalHealthError && e1.what == "vector field" && e1.orbit == 9 && e1.step == 20
    end
end

@testset "AC14b deviation matrix: a zero tangent is a health error naming the column and the step" begin
    _ = (OFLI,)
    flow = HarmonicOscillator(SVector(1.0)); z0 = SVector(1.0, 0.3); I2 = SMatrix{2, 2, Float64}(I)
    e = of_health_error(() -> propagate_orbit(of_Fault(:zero), z0, I2, 4, 2, (OFLI(flow),); orbit = 5))
    @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 5 && e.step == 2 && e.column == 1
    e = of_health_error(() -> propagate_orbit(of_Fault(:col2), z0, I2, 4, 2, (OFLI(flow),); orbit = 6))
    @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 6 && e.step == 2 && e.column == 2
    # non-finite state or deviation: the tangent-core rules apply and are not repeated; one line to show that the reducer does not mask them
    e = of_health_error(() -> propagate_orbit(of_Fault(:zero), SVector(NaN, 0.3), I2, 4, 2, (OFLI(flow),)))
    @test e isa NumericalHealthError && e.what == "state"
end

@testset "AC14c frames: a supplied single vector parallel to the flow is refused before any work, the basis is accepted" begin
    _ = (OFLI, check_frame)
    flow = HenonHeilesFlow(); z0 = of_orbit_R(); f0 = vector_field(flow, z0); u = f0 / norm(f0); eT = of_eps(Float64)
    e = of_health_error(() -> check_frame(OFLI(flow), hcat(u), z0; orbit = 3))
    @test e isa NumericalHealthError && e.what == "initial frame" && e.step == 0 && e.orbit == 3 && e.column == 1
    perp = normalize(of_perp(u))                                                       # a unit vector perpendicular to u, to move the start off the flow direction
    near = normalize(u + 4 * eT * perp)                                               # within 8 eps of parallel: still refused
    @test of_health_error(() -> check_frame(OFLI(flow), hcat(near), z0)) isa NumericalHealthError
    away = normalize(u + 1.0e-6 * perp)                                               # outside the threshold: accepted
    @test check_frame(OFLI(flow), hcat(away), z0) <= 8 * eT
    @test check_frame(OFLI(flow), hcat(u, perp), z0) <= 8 * eT                        # the basis form accepts a column parallel to the flow
    bad = of_health_error(() -> check_frame(OFLI(flow), hcat(perp, 1.0000001 * u), z0))   # a non-unit column is refused in either form
    @test bad isa NumericalHealthError && bad.what == "initial frame" && bad.column == 2
end

@testset "AC14d chains: the orbit is bitwise unchanged, a non-replacing OFLI changes nothing for the others" begin
    _ = (OFLI,)
    flow = HenonHeilesFlow(); im = of_im(Float64); z0 = of_orbit_C(); N = 4000; tau = 20
    W4 = of_frame(of_seed(14, 1), 1, 4, 4, Float64); W2 = SMatrix{4, 2, Float64}(W4[:, 1:2])
    function pts(rs, W)
        _, _, st = propagate_orbit(im, z0, W, N, tau, (of_Trace(), rs...))
        return Chaoslab.reducer_final(of_Trace(), st[1]).xs
    end
    base = pts((), W2)
    @test pts((OFLI(flow),), W2) == base
    @test pts((OFLI(flow; replace = false), SALI()), W2) == base
    @test pts((SALI(), OFLI(flow)), W2) == base                                      # a replaced matrix is the last reducer's business: the orbit does not change (L1 S2)
    @test pts((OFLI(flow),), W4) == pts((OFLI(flow),), W2) == base                   # k = 4 against k = 2
    x_empty, _, _ = propagate_orbit(im, z0, W2, N, tau, ())
    x_ofli, _, _ = propagate_orbit(im, z0, W4, N, tau, (OFLI(flow),))
    @test x_empty == x_ofli
    # the outputs of OFLI in {OFLI, SALI} are bitwise the stand-alone ones (every field, every handoff)
    alone = of_chain((OFLI(flow),), im, z0, W2, N, tau)[1]; chain = of_chain((OFLI(flow; replace = false), SALI()), im, z0, W2, N, tau)[1]
    @test all(i -> isequal(alone[i], chain[i]), eachindex(alone))
end

@testset "AC14e run record and dimension guard" begin
    _ = (OFLI, reducer_provenance)
    flow = HenonHeilesFlow(); im = of_im(Float64); z0 = of_orbit_R()
    res = propagate(im, [z0], 40; k = 2, tau = 20, seed = 1, reducers = (OFLI(flow),))
    p = res.provenance
    res0 = propagate(im, [z0], 40; k = 2, tau = 20, seed = 1)
    @test p.tau == 20 && p.k == 2 && p.precision == res0.provenance.precision
    @test p.indicator == "ofli" && p.ofli_flow == "HenonHeilesFlow" && p.ofli_replace === true
    @test reducer_provenance(OFLI(flow; replace = false)).ofli_replace === false
    # h and the flow of the map are held by the map (im.tau, im.flow): the record of the experiment layer compares them with the reducer's identifier (P10)
    @test im.tau == OF_H && string(nameof(typeof(im.flow))) == p.ofli_flow
    # the reducer refuses a state dimension different from the flow's (P11, extended phase space) and a mismatching flow
    @test_throws ArgumentError Chaoslab.reducer_init(OFLI(flow), Float64, Val(6), Val(2))
    @test_throws ArgumentError propagate(im, [z0], 40; k = 2, tau = 20, seed = 1, reducers = (OFLI(HarmonicOscillator(SVector(1.0))),))
end

# ------------------------------------------------------------------ AC15: precision ladder
@testset "AC15a precision ladder: Float64, Double64 and BigFloat(256) agree, BigFloat against its own reference" begin
    _ = (OFLI, ofli_initial)
    flow = HenonHeilesFlow(); N = 1200; tau = 20; m = min(OF_NCLASS, of_size(3, 50))   # t = 60: one Lyapunov-scale horizon
    keys_ = vcat([:C], [(:cha, i) for i in 1:m])
    vals = Dict(Float64 => Float64[], Double64 => Float64[], BigFloat => Float64[]); nbad = 0; worstB = 0.0
    for key in keys_, T in (Float64, Double64)
        c = of_case(key, T, :o8, N)
        snaps = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau); g0 = ofli_initial(flow, c.x0, c.W0)
        push!(vals[T], Float64(max(maximum(g0), last(snaps).ofli)))
        rows = of_rows(c, T, tau, snaps, g0, [N]); all(r -> r.err <= r.bound, rows) || (nbad += 1)
    end
    setprecision(BigFloat, OF_BITS) do
        for key in keys_[1:min(end, of_size(2, 6))]
            c = of_case(key, BigFloat, :o8, N)
            snaps = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau); g0 = ofli_initial(flow, c.x0, c.W0)
            push!(vals[BigFloat], Float64(max(maximum(g0), last(snaps).ofli)))
            rows = of_rows(c, BigFloat, tau, snaps, g0, [N]); all(r -> r.err <= r.bound, rows) || (nbad += 1)
            worstB = max(worstB, maximum(r.ratio for r in rows))
        end
    end
    @test nbad == 0
    @test of_agree(vals[Float64], vals[Double64]) && of_agree(vals[Float64][1:length(vals[BigFloat])], vals[BigFloat]) &&
          of_agree(vals[Double64][1:length(vals[BigFloat])], vals[BigFloat])
    @info "AC15a reported (R3)" bigfloat_error_over_bound = worstB
end

# a chaotic start of E = 1/6 with a large exponent, memoized: the draw with the largest finite-time Benettin exponent of 24 seeded draws
function of_overflow_start()
    return get!(OF_CACHE, :overflow_start) do
        flow = HenonHeilesFlow(); a, b = section_interval(flow, 1 / 6); w = sqrt(2 / 6); rng = StableRNG(of_seed(15, 1)); best = (-Inf, nothing)
        im = of_im(Float64, :o8, 0.25)
        for j in 1:24
            z, st = random_energy_start(flow, 1 / 6, (a, b), (-w, w), rng); st === :ok || continue
            _, _, s = propagate_orbit(im, z, of_frame(1, j, 4, 1, Float64), 2000, 4, (Benettin(),))
            lam = reducer_final(Benettin(), s[1]).lambda[1]
            lam > best[1] && (best = (lam, z))
        end
        best[2]
    end
end

function of_case_overflow(z0, N)
    return get!(OF_CACHE, (:overflow_case, N)) do
        T = Float64; im = of_im(T, :o8, 0.25); W0 = of_frame(of_seed(15, 2), 1, 4, 4, T)
        pts = of_pts(im, z0, N); Ys = of_refY(im, pts, N)
        (im = im, x0 = z0, W0 = W0, pts = pts, Ys = Ys, ref = of_refseries(HenonHeilesFlow(), Ys, pts, W0), Ms = of_tangents(im, pts, N), N = N)
    end
end

@testset "AC15b representation: OFLI_B above 709.78 is finite, not clipped, and the display conversion reports the overflow" begin
    _ = (OFLI, ofli_initial, fli_unlogged)
    flow = HenonHeilesFlow(); z0 = of_overflow_start(); h = 0.25; im = of_im(Float64, :o8, h); tau = 4
    W0 = of_frame(of_seed(15, 2), 1, 4, 4, Float64); N = of_size(26_000, 120_000)
    snaps = of_snaps(OFLI(flow), im, z0, W0, N, tau)
    G = last(snaps).ofli
    @test isfinite(G) && G > 709.78                                                  # the unlogged value exp(G) overflows in Float64; the log form does not
    @test fli_unlogged(G) === nothing && fli_unlogged(700.0) isa Float64
    @test issorted([s.ofli for s in snaps])                                          # nondecreasing, no clipping at 1e20 or anywhere
    # the same orbit against the BigFloat reference within B_G at the last handoff (cost: N BigFloat steps, package size chosen at stage 02 calibration)
    c = of_case_overflow(z0, N)
    rows = of_rows(c, Float64, tau, snaps, ofli_initial(flow, c.x0, c.W0), [N])
    @test all(r -> r.err <= r.bound, rows)
    @info "AC15b reported" G = G error_over_bound = maximum(r.ratio for r in rows)
end
