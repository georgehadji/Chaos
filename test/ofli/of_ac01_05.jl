# AC1 to AC5: the closed forms (free particle d = 4 and d = 2, harmonic oscillator, pendulum) and the periodic orbit of Henon-Heiles.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError). The "controls" testsets use no new name:
# they run the check of the criterion against test-side doubles with one defect each and pass on the code as it is today.

# ------------------------------------------------------------------ AC1 and AC2 helpers: the free particle in closed form (BigFloat)
# Closed form of g(n) = ln ||P w(n)|| for the column e_i of the canonical frame of the free particle with N degrees of freedom: w = (w_q + t w_p, w_p),
# field f = (p0, 0), t = n h. Call inside setprecision.
function of_free_g(i::Integer, p0b, t, ndof::Integer)
    d = 2 * ndof; w = zeros(BigFloat, d); w[i] = 1
    for j in 1:ndof
        w[j] += t * w[ndof + j]
    end
    f = zeros(BigFloat, d); f[1:ndof] .= p0b
    pw = w - (dot(f, w) / dot(f, f)) * f
    b = norm(pw)
    return b > 0 ? log(b) : -BigFloat(Inf)
end
# Worst ratio |G_comp - G_cf| / (16 eps max(1, |G_cf|)) over the checkpoints, for the single columns `cols` of the canonical frame (k = 1 runs read G_i). `mk(flow)`
# builds the reducer, `init(flow, x0, W0)` the s = 0 terms of the computed series. The closed form of the running maximum is g(n) itself (g is nondecreasing, checked
# by the closed form: the columns of the canonical frame give constant or increasing g).
function of_ac1_ratio(mk, init, ::Type{T}, p0, cols, N, tau; ndof = 2, h = 1 / 16) where {T}
    flow = of_FreeFlow{ndof}(); sys = of_FreeMap(flow, T(h)); d = 2 * ndof
    z0 = SVector{d, T}(ntuple(j -> j <= ndof ? zero(T) : T(p0[j - ndof]), Val(d)))
    worst = 0.0
    for i in cols
        W0 = SMatrix{d, 1, T}(ntuple(j -> j == i ? one(T) : zero(T), Val(d)))
        snaps = of_snaps(mk(flow), sys, z0, W0, N, tau)
        g0 = init(flow, z0, W0)[1]
        cps = of_checkpoints(N, tau)
        setprecision(BigFloat, OF_BITS) do
            p0b = [of_big(T(x)) for x in p0]; hb = of_big(T(h)); gcf0 = of_free_g(i, p0b, zero(BigFloat), ndof)
            for (kk, s) in enumerate(snaps)
                n = kk * tau; n in cps || continue
                gc = max(g0, s.ofli1); gf = max(gcf0, of_free_g(i, p0b, n * hb, ndof))
                if isinf(gf)
                    worst = max(worst, gc == gf ? 0.0 : Inf)
                else
                    worst = max(worst, Float64(abs(of_big(gc) - gf)) / (16 * of_eps(T) * max(1.0, abs(Float64(gf)))))
                end
            end
        end
    end
    return worst
end

@testset "AC1 algebra on the free particle, d = 4 (exact step, closed form)" begin
    _ = (OFLI, ofli_initial)
    # h = 2^-4, so the tangent [[I, h I], [0, I]] acts exactly on the canonical columns: no discretization error and the only rounding is the module's (spec point U1).
    N = 10_000
    for T in (Float64, Double64), tau in (1, 10)
        p0 = (T(6) / 10, T(8) / 10)                                       # the spec's p0 = (0.6, 0.8)
        @test of_ac1_ratio(f -> OFLI(f), ofli_initial, T, p0, 1:4, N, tau) <= 1
        q0 = (T(1) / 2, zero(T))                                           # p0 = (1/2, 0): column 3 has w_p parallel to p0 (constant), column 4 perpendicular (ln t)
        @test of_ac1_ratio(f -> OFLI(f), ofli_initial, T, q0, 3:4, N, tau) <= 1
    end
    # the shapes: constant for w_p parallel to p0, ln t growth for w_p perpendicular (p0 = (1/2, 0), tau = 10, Float64)
    flow = of_FreeFlow{2}(); sys = of_FreeMap(flow, 1 / 16); z0 = SVector(0.0, 0.0, 0.5, 0.0)
    W3 = SMatrix{4, 1, Float64}(0, 0, 1, 0); W4 = SMatrix{4, 1, Float64}(0, 0, 0, 1)
    s3 = of_snaps(OFLI(flow), sys, z0, W3, N, 10); s4 = of_snaps(OFLI(flow), sys, z0, W4, N, 10)
    @test all(s -> abs(s.ofli1) <= 16 * of_eps(Float64), s3)               # ln 1 at every handoff (the renormalised columns carry rounding)
    ts = [s.n / 16 for s in s4]; keep = findall(t -> 10 <= t <= 600, ts)
    @test 0.98 <= of_slope(log.(ts[keep]), [s4[i].ofli1 for i in keep]) <= 1.02
    # k = 4 canonical frame: OFLI_B is the maximum of the four closed forms
    p0 = (0.6, 0.8); W = SMatrix{4, 4, Float64}(I)
    sB = of_snaps(OFLI(flow), sys, SVector(0.0, 0.0, p0...), W, 1000, 10)
    gmax = setprecision(BigFloat, OF_BITS) do
        [maximum(of_free_g(i, of_big.(p0), s.n * of_big(1 / 16), 2) for i in 1:4) for s in sB]
    end
    gB0 = maximum(ofli_initial(flow, SVector(0.0, 0.0, p0...), W))
    @test all(i -> abs(max(gB0, sB[i].ofli) - Float64(gmax[i])) <= 16 * of_eps(Float64) * max(1.0, abs(Float64(gmax[i]))), eachindex(sB))
end

@testset "AC1 negative controls (test-side doubles, no module name)" begin
    # the check of AC1 passes for the correct double and fails for each defective one (N = 1000, tau = 1: plain sums in the doubles)
    p0 = (0.6, 0.8)
    ok = of_ac1_ratio(f -> of_Double(f, :ok), (f, x, W) -> of_initial(:ok, f, x, W), Float64, p0, 1:4, 1000, 1)
    @test ok <= 1
    for mode in (:unprojected, :parallel)
        @test of_ac1_ratio(f -> of_Double(f, mode), (f, x, W) -> of_initial(:ok, f, x, W), Float64, p0, 1:4, 1000, 1) > 1
    end
    # spec point U2: on the free particle g is nondecreasing, so a double without the s = 0 term and one with the field at the previous point cannot fail here
    # (the field is constant along the orbit); both controls are in AC3c, on the oscillator.
end

@testset "AC2 free particle d = 2" begin
    _ = (OFLI, ofli_initial, check_frame)
    N = 10_000; flow = of_FreeFlow{1}()
    for T in (Float64, Double64), tau in (1, 10)
        h = T(1) / 16; sys = of_FreeMap(flow, h); z0 = SVector{2, T}(0, T(1) / 2)
        W = SMatrix{2, 1, T}(0, 1)                                                  # w = (0, 1): w_p = 1, OFLI_1 = ln |w_p| = 0 for all n
        snaps = of_snaps(OFLI(flow), sys, z0, W, N, tau); g0 = ofli_initial(flow, z0, W)[1]
        @test g0 == 0
        @test all(s -> abs(max(g0, s.ofli1)) <= 16 * of_eps(T), snaps)
        # FLI_1 = ln sqrt(t^2 + 1) within the L1 sum bound 8 eps Psi_n (Psi from the product norms of the constant tangent)
        Ms = fill(SMatrix{2, 2, Float64}(1, 0, 1 / 16, 1), N)
        for n in (10, 100, 1000, N)
            n % tau == 0 || continue
            s = snaps[n ÷ tau]; lw = [0.5 * log((j / 16)^2 + 1) for j in 0:n]
            Fcf = 0.5 * log((n / 16)^2 + 1)
            bound = 8 * of_eps(T) * of_Psi(of_prodnorms(Ms, n), lw, n)
            @test abs(Float64(s.fli1) - Fcf) <= max(bound, 16 * of_eps(T) * Fcf)
        end
    end
    # a general unit start w = (0.6, 0.8): ||w_perp|| = |w_p| = 0.8 constant; within 8 eps (Psi + 1 / sin theta)
    sys = of_FreeMap(flow, 1 / 16); z0 = SVector(0.0, 0.5); W = SMatrix{2, 1, Float64}(0.6, 0.8)
    Ms = fill(SMatrix{2, 2, Float64}(1, 0, 1 / 16, 1), 1000)
    snaps = of_snaps(OFLI(flow), sys, z0, W, 1000, 1)
    for n in (10, 100, 1000)
        lw = zeros(n + 1); w = SVector(0.6, 0.8)
        for j in 1:n
            w = Ms[j] * w; lw[j + 1] = log(norm(w)); w = w / norm(w); lw[j + 1] += lw[j]
        end
        bound = of_BG(of_eps(Float64), of_Psi(of_prodnorms(Ms, n), lw, n), 0.8)
        @test abs(max(ofli_initial(flow, z0, W)[1], snaps[n].ofli1) - log(0.8)) <= bound
    end
    # a start parallel to the flow (w = (1, 0), f0 = (p0, 0)): the frame check refuses the single vector, accepts the basis (O13)
    x0 = SVector(0.0, 0.5); Wpar = SMatrix{2, 1, Float64}(1, 0)
    e = of_health_error(() -> check_frame(OFLI(flow), Wpar, x0; orbit = 7))
    @test e isa NumericalHealthError && e.what == "initial frame" && e.step == 0 && e.orbit == 7 && e.column == 1
    @test check_frame(OFLI(flow), SMatrix{2, 2, Float64}(1, 0, 0, 1), x0) <= 8 * of_eps(Float64)
    @test check_frame(OFLI(flow), SMatrix{2, 1, Float64}(0, 1), x0) <= 8 * of_eps(Float64)
    # the stream does not refuse it: every handoff is a minus infinity term, counted
    snaps = of_snaps(OFLI(flow), of_FreeMap(flow, 1 / 16), x0, Wpar, 100, 10)
    @test last(snaps).ofli1 == -Inf && last(snaps).nminf == 10
end

# ------------------------------------------------------------------ AC3: the harmonic oscillator
# Exact orbit and tangent of F1 in BigFloat from the Float64 frequency widened exactly: z(t), Y(t) w0 and the field norm; call inside setprecision.
function of_osc_exact(w::Float64, z0, w0, t)
    flow = HarmonicOscillator(SVector(of_big(w))); zb = SVector{2, BigFloat}(of_big.(z0)); wb = SVector{2, BigFloat}(of_big.(w0))
    zt = exact_flow(flow, zb, t); wt = exact_flow_matrix(flow, t) * wb
    return flow, zt, wt
end
# ln ||f(z(t_m))|| at the handoff times t_m = m * dt (m = 0 .. M), BigFloat
function of_osc_lnf(w::Float64, z0, M::Integer, dt)
    flow = HarmonicOscillator(SVector(of_big(w))); zb = SVector{2, BigFloat}(of_big.(z0))
    return [log(norm(vector_field(flow, exact_flow(flow, zb, m * of_big(dt))))) for m in 0:M]
end
# OFLI closed form of the columns of the frame W0 (d = 2) at the handoff steps n = 0, tau, 2 tau, ..: G_i = ln |gradH0 . w_i| - min ln ||f|| (running), the basis G_B = max_i.
function of_osc_G(w::Float64, z0, W0, tau::Integer, h, M::Integer)
    return setprecision(BigFloat, OF_BITS) do
        lnf = of_osc_lnf(w, z0, M, tau * of_big(h)); minf = accumulate(min, lnf)
        gH = [of_big(w)^2 * of_big(z0[1]), of_big(z0[2])]
        c = [log(abs(gH[1] * of_big(W0[1, i]) + gH[2] * of_big(W0[2, i]))) for i in 1:size(W0, 2)]
        (G = [[c[i] - minf[m + 1] for m in 0:M] for i in eachindex(c)], lnf = lnf)
    end
end
function of_ac3_errors(m, w, z0, W0, h, tmax)
    flow = HarmonicOscillator(SVector(w)); sys = IntegratorMap(m, flow, h)
    tau = round(Int, 1 / h); N = round(Int, tmax / h); M = N ÷ tau
    cf = of_osc_G(w, z0, W0, tau, h, M)
    x0 = SVector{2, Float64}(z0); g0 = ofli_initial(flow, x0, W0)
    snaps = of_snaps(OFLI(flow), sys, x0, W0, N, tau)
    GBcf = [maximum(cf.G[i][m + 1] for i in axes(W0, 2)) for m in 0:M]; G1cf = cf.G[1]
    errB = [abs(max(maximum(g0), snaps[m].ofli) - Float64(GBcf[m + 1])) for m in 1:M]
    err1 = [abs(max(g0[1], snaps[m].ofli1) - Float64(G1cf[m + 1])) for m in 1:M]
    return (errB = errB, err1 = err1)
end

@testset "AC3a oscillator: the exact-flow identity (the oracle itself, 256 bits)" begin
    # ||w_perp(t)|| = |gradH(z0) . w(0)| / ||f(z(t))|| to 1e-30 at 100 sample times (spec O5; closed form of the exact flow, not the module)
    w = 1.3; z0 = (1.0, 0.3); w0 = (0.6, 0.8); worst = 0.0
    setprecision(BigFloat, OF_BITS) do
        gH = [of_big(w)^2 * of_big(z0[1]), of_big(z0[2])]; c = abs(gH[1] * of_big(w0[1]) + gH[2] * of_big(w0[2]))
        for t in range(0.37, 99.0, length = 100)
            flow, zt, wt = of_osc_exact(w, z0, w0, of_big(t)); f = vector_field(flow, zt)
            pw = wt - (dot(f, wt) / dot(f, f)) * f
            worst = max(worst, Float64(abs(norm(pw) - c / norm(f))))
        end
    end
    @test worst <= 1.0e-30
end

@testset "AC3b oscillator: order of the discretization error of the OFLI (order-2 and order-4 tables)" begin
    _ = (OFLI, ofli_initial)
    w = 1.3; z0 = (1.0, 0.3); W1 = SMatrix{2, 1, Float64}(0.6, 0.8); W2 = SMatrix{2, 2, Float64}(0.6, 0.8, -0.8, 0.6)
    for (m, p) in ((Leapfrog(Float64), 2), (ForestRuth(Float64), 4)), W in (W1, W2)
        e1 = of_ac3_errors(m, w, z0, W, 0.1, 1000.0); e2 = of_ac3_errors(m, w, z0, W, 0.05, 1000.0)
        k = size(W, 2) == 1 ? :err1 : :errB
        a, b = getfield(e1, k)[100], getfield(e2, k)[100]                           # t = 100 (handoff m = t, one per time unit)
        @test b > 0 && 2^p / 2 <= a / b <= 2 * 2^p
        @test maximum(getfield(e2, k)) <= (p == 2 ? 1.0e-2 : 1.0e-3)               # sanity: the error stays small to t = 1000 (grows with t, O11 iii; calibrated)
        @info "AC3b reported" table = p form = k error_h01 = a error_h005 = b max_to_t1000 = maximum(getfield(e2, k))
    end
end

@testset "AC3c controls: doubles with one defect fail the closed form (no module name)" begin
    # the lagged double (field at the previous point, O(h)) and the double without the s = 0 term fail the closed form of the oscillator; the correct double does not
    # (h = 0.05, tau = 1, t = 20, order-4 table). The start (0, 0.3) is on q = 0, where ||f|| is smallest, so g(0) is the largest value and the first handoff is below it.
    w = 1.3; h = 0.05; N = 400; W = SMatrix{2, 1, Float64}(0.6, 0.8)
    flow = HarmonicOscillator(SVector(w)); sys = IntegratorMap(ForestRuth(Float64), flow, h)
    for (mode, z0, bad) in ((:ok, SVector(1.0, 0.3), false), (:ok, SVector(0.0, 0.3), false), (:lagged, SVector(1.0, 0.3), true), (:no_s0, SVector(0.0, 0.3), true))
        cf = of_osc_G(w, z0, W, 1, h, N)
        snaps = of_snaps(of_Double(flow, mode === :no_s0 ? :ok : mode), sys, z0, W, N, 1); g0 = of_initial(mode, flow, z0, W)[1]
        err = [abs(max(g0, snaps[m].ofli1) - Float64(cf.G[1][m + 1])) for m in 1:N]
        @test bad ? maximum(err) > 1.0e-3 : maximum(err) < 1.0e-5
        mode === :no_s0 && @test err[1] > 1.0e-3                                      # fails at the first handoff
    end
end

# ------------------------------------------------------------------ AC4: the pendulum (reference: Double64 order-10 table, the Taylor BigFloat reference is stage 04)
# Period of the libration of amplitude A of H = p^2/2 - cos q: P = 4 K(sin(A/2)) = 2 pi / agm(1, cos(A/2)).
function of_pend_period(A)
    a, b = 1.0, cos(A / 2)
    for _ in 1:40
        a, b = (a + b) / 2, sqrt(a * b)
    end
    return 2pi / a
end
# Reference run: Double64 order-10 table at h/40 = 1/40, the cumulative tangent of the unit start w0, handoffs every `per` steps. Returns the reference states,
# columns and ||f|| at every step (for the continuous minimum) and at every handoff.
function of_pend_ref(A::Float64, w0, T::Float64; per = 40)
    D = Double64; m = SofroniouSpaletta10(D); flow = Pendulum(D(1)); h = D(1) / per
    z = SVector{2, D}(D(A), D(0)); W = SMatrix{2, 1, D}(D(w0[1]), D(w0[2]))
    N = round(Int, T * per); zs = Vector{typeof(z)}(undef, N + 1); Ws = Vector{typeof(W)}(undef, N + 1); zs[1] = z; Ws[1] = W
    for n in 1:N
        z, W = tangent_advance(m, flow, z, W, h); zs[n + 1] = z; Ws[n + 1] = W
    end
    return (m = m, flow = flow, h = h, zs = zs, Ws = Ws, per = per, N = N)
end
of_pend_f(flow, z) = norm(vector_field(flow, z))
# The continuous minimum of ||f(z(s))|| over s in [lo, hi] (time units): scan the stored steps, then golden-section refine by sub-steps of the method.
function of_pend_minf(R, lo, hi)
    ilo, ihi = round(Int, lo * R.per) + 1, min(round(Int, hi * R.per) + 1, R.N + 1)
    k = argmin([of_pend_f(R.flow, R.zs[i]) for i in ilo:ihi]) + ilo - 1
    k = clamp(k, 2, R.N)
    z0 = R.zs[k - 1]; a, b = Double64(0), 2 * R.h; gr = (sqrt(Double64(5)) - 1) / 2
    ff(s) = of_pend_f(R.flow, advance(R.m, R.flow, z0, s))
    c, d = b - gr * (b - a), a + gr * (b - a); fc, fd = ff(c), ff(d)
    for _ in 1:70
        if fc < fd
            b, d, fd = d, c, fc; c = b - gr * (b - a); fc = ff(c)
        else
            a, c, fc = c, d, fd; d = a + gr * (b - a); fd = ff(d)
        end
    end
    return Float64(min(fc, fd))
end

@testset "AC4a pendulum: the exact-flow identity on the reference (Double64 order-10 stand-in; the 1e-30 Taylor form is stage 04)" begin
    for A in (2.0, 3.0)
        w0 = (0.6, 0.8); R = of_pend_ref(A, w0, 30.0)
        c0 = abs(dot(hamiltonian_gradient(R.flow, R.zs[1]), SVector(R.Ws[1][1], R.Ws[1][2])))
        worst = 0.0
        for n in 1:R.per:(R.N + 1)
            z, W = R.zs[n], R.Ws[n]; f = vector_field(R.flow, z); w = SVector(W[1], W[2])
            pw = w - (dot(f, w) / dot(f, f)) * f
            worst = max(worst, Float64(abs(norm(pw) - c0 / norm(f)) / (c0 / norm(f))))
        end
        @test worst <= 1.0e-18                                                      # calibrated at stage 02 (order-10 table, 1/40)
    end
end

@testset "AC4b pendulum: the exact OFLI_1 is periodic (continuous minimum), the exact FLI_1 grows as ln t" begin
    for A in (2.0, 3.0)
        P = of_pend_period(A); w0 = (0.6, 0.8); K = of_size(60, 600)             # periods covered by the reference (full: 1000 P in the spec, D38 here)
        R = of_pend_ref(A, w0, K * P)
        # the minimum of ||f|| in every one of the first 10 periods is the same number: OFLI_1 constant after the first period within 1e-12
        mins = [of_pend_minf(R, (k - 1) * P, k * P) for k in 1:10]
        @test maximum(abs.(log.(mins) .- log(mins[1]))) <= 1.0e-12
        # FLI_1 = ln ||w|| at the unit-time handoffs: slope against ln t over [10 P, K P] in [0.8, 1.2]
        ts = collect(Float64, (10P):1.0:(K * P)); idx = round.(Int, ts .* R.per) .+ 1
        S = [Float64(log(norm(R.Ws[i]))) for i in idx]
        @test 0.8 <= of_slope(log.(ts), S) <= 1.2
    end
end

@testset "AC4c pendulum: the numerical OFLI against the exact one, order-4 table at h = 1/20" begin
    _ = (OFLI, ofli_initial)
    for A in (2.0, 3.0)
        P = of_pend_period(A); w0 = (0.6, 0.8); W = SMatrix{2, 1, Float64}(w0); tmax = 10P
        R = of_pend_ref(A, w0, tmax + 1)
        flow = Pendulum(1.0); h = 1 / 20; sys = IntegratorMap(ForestRuth(Float64), flow, h); x0 = SVector(A, 0.0)
        M = floor(Int, tmax); snaps = of_snaps(OFLI(flow), sys, x0, W, 20M, 20); g0 = ofli_initial(flow, x0, W)[1]
        # exact discrete series (handoffs every unit time): G(m) = max over m' <= m of ln ||P w(m')|| from the reference columns
        gex = [(z = R.zs[1 + R.per * m]; w = SVector(R.Ws[1 + R.per * m][1], R.Ws[1 + R.per * m][2]); f = vector_field(R.flow, z);
                Float64(log(norm(w - (dot(f, w) / dot(f, f)) * f)))) for m in 0:M]
        Gex = accumulate(max, gex)
        err = [abs(max(g0, snaps[m].ofli1) - Gex[m + 1]) for m in 1:M]
        @test maximum(err) <= 1.0e-5                                                 # delta: calibrated at stage 02 (spec: order 1e-6)
        # R1: the plateau-loss time t* (first handoff where the computed OFLI exceeds the plateau by 0.05), at h and h/2, order-4 table: reported, not scored
        tstar(hh) = (s = IntegratorMap(ForestRuth(Float64), flow, hh); tau = round(Int, 1 / hh); n = of_size(5000, 20000);
                     sn = of_snaps(OFLI(flow), s, x0, W, n * tau, tau); pl = max(g0, sn[M].ofli1);
                     k = findfirst(x -> x.ofli1 > pl + 0.05, sn); k === nothing ? Inf : Float64(k))
        @info "AC4c reported (R1)" amplitude = A tstar_h_1_4 = tstar(1 / 4) tstar_h_1_8 = tstar(1 / 8) tstar_h_1_20 = tstar(1 / 20)
    end
end

# ------------------------------------------------------------------ AC5: a periodic orbit of Henon-Heiles (Algorithm item 8)
# Newton iteration for the fixed point of the first return map of the section q1 = 0, p1 > 0 at energy E, by the module's section_map (Jacobian of the leg).
# Returns (z0, T, A) with A the 2 x 2 return map Jacobian (the reduced monodromy) and the closure error; `nothing` if not converged.
function of_periodic(E, q2, p2; h = 1 / 80)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8)
    ret(q2, p2) = (z = energy_start(flow, E, q2, p2)[1]; z === nothing ? nothing : section_map(m, flow, z, h, 4000; jacobian = true, max_crossings = 2))
    x = SVector(q2, p2)
    for _ in 1:30
        r = ret(x[1], x[2]); (r === nothing || length(r.z) < 2) && return nothing
        y = SVector(r.z[2][2], r.z[2][4]); F = y - x
        if maximum(abs, F) <= 1.0e-13
            return (z0 = energy_start(flow, E, x[1], x[2])[1], T = r.t[2], A = r.J[2], closure = maximum(abs, F), x = x)
        end
        x = x - (r.J[2] - one(r.J[2])) \ F
    end
    return nothing
end

@testset "AC5 periodic orbit of F3: bounded OFLI against linearly growing FLI" begin
    _ = (OFLI, ofli_initial, section_map, energy_start)
    po = of_periodic(OF_E, 0.1, 0.0)
    @test po !== nothing && po.closure <= 1.0e-13
    A = po.A; T = po.T
    @test abs(tr(A)) < 2 && abs(det(A) - 1) < 1.0e-6                                 # A2: elliptic reduced monodromy
    dE = 1.0e-3; plo = of_periodic(OF_E - dE, po.x[1], po.x[2]); phi = of_periodic(OF_E + dE, po.x[1], po.x[2])
    Tprime = (phi.T - plo.T) / (2dE)
    @test abs(Tprime) > 1.0e-3                                                       # T'(E) not zero (the FLI of the same orbit grows)
    @info "AC5 reported" period = T trace = tr(A) Tprime = Tprime
    flow = HenonHeilesFlow(); n10, n100 = round(Int, 10T), round(Int, 100T)
    frames = (SMatrix{4, 4, Float64}(I), of_frame(of_seed(5, 1), 1, 4, 4, Float64))
    for W in frames
        # (a) the reference series: order-10 table at h = 1/40 (Float64; the Taylor BigFloat series is stage 04), handoffs every 40 steps = 1 time unit
        for (m, h, tau, tag) in ((SofroniouSpaletta10(Float64), 1 / 40, 40, :ref), (KahanLi(Float64, 8), 1 / 20, 20, :o8))
            sys = IntegratorMap(m, flow, h); snaps = of_snaps(OFLI(flow), sys, po.z0, W, n100 * tau, tau)
            G(n) = max(maximum(ofli_initial(flow, po.z0, W)), snaps[n].ofli)
            @test G(n100) - G(n10) <= 0.1                                            # (a), (b): ln OFLI_B(100 T) - ln OFLI_B(10 T)
            @test snaps[n100].fli - snaps[n10].fli >= 1.5                            # the FLI_B of the same columns grows like ln t
        end
        # (c) the order-4 table: reported with the plateau-loss time
        sys = IntegratorMap(ForestRuth(Float64), flow, 1 / 20); snaps = of_snaps(OFLI(flow), sys, po.z0, W, n100 * 20, 20)
        @info "AC5c reported (R1)" dG = snaps[n100].ofli - snaps[n10].ofli dFLI = snaps[n100].fli - snaps[n10].fli
    end
end
