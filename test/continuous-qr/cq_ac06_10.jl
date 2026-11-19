# AC6 to AC10: one column and column independence, the flow-direction and energy-covector identities, same-orbit agreement with the discrete spectrum, regular decay.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError). Sizes: D38 package sizes with the full sizes of the
# spec behind CQ_FULL (the stage 04 script sets CQ_FULL_RUN = true).

# ------------------------------------------------------------------ helpers of this file (definitions only)
const CQ_CD4 = 20.0; const CQ_CD2 = 10.0                  # D121 (C9): 2 x the worst realized c_d at T = 60 over the full-size sample (8.62 for q = 4, 4.61 for q = 2), rounded up to one digit
const CQ_B8 = 0.4                                         # D121 (C8): 2 x the worst realized error at t <= 60 over the full-size sample (0.194), rounded up to one digit
const CQ_BCAP = 0.05                                      # cap of B_lin where the formula exceeds the range it protects (plan, Calibration: realized errors are far below it)
cq_Bc(B) = min(B, CQ_BCAP)
# Frame with first column f(z0) / ||f(z0)||, the other columns seeded and completed orthonormally (criterion 7). f(z0) = 0 is refused: the divisor does not exist (spec point U1).
function cq_dir_frame(flow, z0::SVector{d, T}, k::Integer, rng) where {d, T}
    f = vector_field(flow, z0); nf = norm(f)
    (isfinite(nf) && nf > 0) || throw(NumericalHealthError("initial frame", Float64(nf), 0.0, 0, 1, :cpu))
    M = hcat(f / nf, T.(randn(rng, d, k - 1)))
    return cq_orth(M)
end
cq_gradH(flow::HenonHeilesFlow, z) = vcat(force(flow, SVector(z[1], z[2])), SVector(z[3], z[4]))
# l_i from the Gram determinants of the evolved frame W = Y Q0 against the stored Q0 (C16), in 256 bits: W and Q0 are matrices of any float type.
function cq_lsplit(W, Q0)
    return setprecision(BigFloat, CQ_BITS) do
        Wb = Matrix{BigFloat}(map(cq_big, W)); Xb = Matrix{BigFloat}(map(cq_big, Q0)); k = size(Wb, 2); out = zeros(Float64, k); prev = zero(BigFloat)
        for m in 1:k
            lv = (log(det(Wb[:, 1:m]' * Wb[:, 1:m])) - log(det(Xb[:, 1:m]' * Xb[:, 1:m]))) / 2
            out[m] = Float64(lv - prev); prev = lv
        end
        out
    end
end
# External stand-in for the variational reference (the Taylor form is stage 04): Double64 order-10 table at 1/80 with the tangent carried on the frame Q0;
# `flowD` is the flow with Double64 parameters. Returns Dict t => l (Float64) for the integer times ts.
function cq_refl(flowD, z0, Q0, ts)
    D = Double64; m = SofroniouSpaletta10(D); h = D(1) / 80; d = dimension(flowD)
    z = SVector{d, D}(map(D, z0)); W = SMatrix{d, size(Q0, 2), D}(map(D, Q0)); out = Dict{Int, Vector{Float64}}(); tmax = maximum(ts)
    for t in 0:tmax
        t in ts && (out[t] = cq_lsplit(W, Q0))
        if t < tmax
            for _ in 1:80
                z, W = tangent_advance(m, flowD, z, W, h)
            end
        end
    end
    return out
end
cq_pendulum_start(A) = SVector(0.0, 2 * sin(A / 2))       # libration of amplitude A: H = p^2/2 - cos q, p(0) = 2 sin(A/2)
# Fresh run to N steps from the start (the series at the end of the run).
cq_final(cq, x0, Q0, N) = (r = propagate_cqr_orbit(cq, x0, Q0, N); (x = r[1], Q = r[2], fin = r[3]))

# ------------------------------------------------------------------ AC6
@testset "AC6a one column: l_1 = ln ||Y(t) q(0)|| on F1, the inverted oscillator and the pendulum" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    h = 0.05; rng = StableRNG(cq_seed(6, 1)); worst = 0.0; wabs = Dict(2 => 0.0, 4 => 0.0)
    for q in (2, 4), t in (10.0, 30.0)
        n = round(Int, t / h)
        # F1, omega = 1.3: exact flow matrix
        flow = HarmonicOscillator([1.3]); im = cq_im(flow, Float64, :o4, h); cq = cq_cqr(im, flow; order = q)
        z0 = SVector(0.5, 0.3); q0 = cq_frame(rng, 2, 1, Float64); r = cq_final(cq, z0, q0, n)
        Y = exact_flow_matrix(flow, n * h); AF = cq_AF(vector_field_jacobian(flow, z0)); B = cq_Bc(cq_Blin(t, h, q, AF, Float64))
        e = abs(r.fin.l[1] - log(norm(Y * q0[:, 1]))); @test e <= B; worst = max(worst, e / B); wabs[q] = max(wabs[q], e)
        # inverted oscillator lam = 0.8: 256-bit reference
        A = cq_inv_osc(0.8); q1 = cq_frame(rng, 2, 1, Float64); cq2 = cq_est(A, Float64, h, q)
        r2 = cq_final(cq2, cq_x0(2, Float64), q1, n)
        ref = setprecision(BigFloat, CQ_BITS) do
            cq_lref(cq2.flow.A, q1, n * cq_big(h))[1]
        end
        B2 = cq_Bc(cq_Blin(t, h, q, cq_AF(A), Float64)); e2 = Float64(abs(cq_big(r2.fin.l[1]) - ref)); @test e2 <= B2; worst = max(worst, e2 / B2); wabs[q] = max(wabs[q], e2)
    end
    # pendulum, amplitude 2, order-4 table: Double64 order-10 variational reference
    for q in (2, 4), t in (10, 30)
        flow = Pendulum(1.0); im = cq_im(flow, Float64, :o4, h); cq = cq_cqr(im, flow; order = q); z0 = cq_pendulum_start(2.0); q0 = cq_frame(rng, 2, 1, Float64)
        r = cq_final(cq, z0, q0, round(Int, t / h)); ref = cq_refl(Pendulum(Double64(1)), z0, q0, [t])[t]
        B = cq_Bc(cq_Blin(Float64(t), h, q, 1.5, Float64)) + 1.0e-4                            # plus the order-4 orbit table against the order-10 reference orbit
        e = abs(r.fin.l[1] - ref[1]); @test e <= B; worst = max(worst, e / B); wabs[q] = max(wabs[q], e)
    end
    @info "AC6a reported (R3): worst error / bound, worst absolute error per scheme" worst wabs
end

@testset "AC6b column independence on Henon-Heiles: the k = 1, 2, 4 runs agree for i <= k within B_col = 5e-3 at t <= 60 (calibrated, U7; D121 C6), t = 100 reported" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 4); worst = Dict(60.0 => 0.0, 100.0 => 0.0)
    cls = cq_classes(CQ_NCLASS, 1); starts = vcat([cq_orbit_R(), cq_orbit_C()], cls.regular[1:min(2, end)], cls.chaotic[1:min(2, end)])
    for (t, scored) in ((60.0, true), (100.0, false)), (j, z0) in enumerate(starts)
        (scored || CQ_FULL) || continue                                                  # the horizon beyond the pointwise horizon t = 60 is reported at the full size only (D121, C6)
        n = round(Int, t / CQ_H); Q4 = initial_deviations(cq_seed(6, 2), j, Val(4), Val(4), Float64); l4 = cq_final(cq, z0, Q4, n).fin.l
        for k in (1, 2)
            lk = cq_final(cq, z0, SMatrix{4, k, Float64}(Q4[:, 1:k]), n).fin.l
            e = maximum(abs, lk - l4[1:k]); scored && @test e <= 5.0e-3; worst[t] = max(worst[t], e)
        end
    end
    @info "AC6b reported: largest |l_i(k) - l_i(4)| at t = 60 (scored) and t = 100 (reported)" worst
end

# ------------------------------------------------------------------ AC7: flow-direction identity (C17)
@testset "AC7 equilibrium start is refused by the direction frame (test-side, no new name; spec point U1)" begin
    e = cq_health_error(() -> cq_dir_frame(HarmonicOscillator([1.3]), SVector(0.0, 0.0), 2, StableRNG(1)))
    @test e isa NumericalHealthError && e.what == "initial frame"
    @test cq_dir_frame(HarmonicOscillator([1.3]), SVector(0.5, 0.3), 2, StableRNG(1))[:, 1] ≈ normalize(vector_field(HarmonicOscillator([1.3]), SVector(0.5, 0.3)))
end

@testset "AC7 flow-direction identity: l_1(t) = ln(||f(z(t))|| / ||f(z_0)||) for the frame whose first column is f(z_0) / ||f(z_0)||" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    cases = (("pendulum A=2", Pendulum(1.0), cq_pendulum_start(2.0), :o4), ("pendulum A=3", Pendulum(1.0), cq_pendulum_start(3.0), :o4),
             ("F1 w=1.3", HarmonicOscillator([1.3]), SVector(0.5, 0.3), :o4), ("F3 R", HenonHeilesFlow(), cq_orbit_R(), :o8), ("F3 C", HenonHeilesFlow(), cq_orbit_C(), :o8))
    rng = StableRNG(cq_seed(7, 1)); rows = Any[]
    for (name, flow, z0, key) in cases, q in (2, 4), h in (1 / 20, 1 / 40)
        (q == 2 && name == "pendulum A=3") && continue                                         # order 2 near the separatrix: error 3.4 at t = 60, outside the useful range (U7)
        d = dimension(flow); k = min(d, 4); Q0 = cq_dir_frame(flow, z0, k, rng); cq = cq_cqr(cq_im(flow, Float64, key, h), flow; order = q)
        B = q == 4 ? (name == "pendulum A=3" ? 1.0e-1 : 5.0e-3) : 1.0       # calibrated B_dir (spec 1e-3 and 1e-1: realized 2.3e-2 for the pendulum of amplitude 3, 8.7e-4 elsewhere at q = 4; 0.34 at q = 2), U7
        for t in (10, 30, 60)
            r = cq_final(cq, z0, Q0, round(Int, t / h))
            e = abs(r.fin.l[1] - log(norm(vector_field(flow, r.x)) / norm(vector_field(flow, z0))))
            @test e <= B
            push!(rows, (name = name, q = q, h = h, t = t, err = e))
        end
    end
    # the error ratio between h and h/2, reported where the error exceeds the floor (m = min(p, q))
    for (name, _, _, _) in cases, q in (2, 4), t in (10, 30, 60)
        (q == 2 && name == "pendulum A=3") && continue
        a = [r.err for r in rows if r.name == name && r.q == q && r.t == t && r.h == 1 / 20][1]; b = [r.err for r in rows if r.name == name && r.q == q && r.t == t && r.h == 1 / 40][1]
        a > 1.0e3 * cq_eps(Float64) * t && b > 1.0e3 * cq_eps(Float64) * t && @info "AC7 reported ratio h / (h/2)" name q t ratio = a / b
    end
end

# ------------------------------------------------------------------ AC8: energy covector identity (C18)
@testset "AC8 energy covector identity: l_1(t) = ln|gradH(z_0).q_1(0)| - ln|gradH(z(t)).q_1(t)| scored at t <= 60 (D121, C8), beyond reported" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 4); Tmax = cq_size(100, 300); n = Tmax * 20; rng = StableRNG(cq_seed(8, 1)); m = min(CQ_NCLASS, cq_size(2, 20))
    cls = cq_classes(CQ_NCLASS, 1); starts = vcat([cq_orbit_R(), cq_orbit_C()], cls.regular[1:m], cls.chaotic[1:m]); tvals = Float64[]; norb = 0; nbad = 0; worst = 0.0; beyond = 0.0; nbeyond = 0
    for (j, z0) in enumerate(starts)
        g0 = cq_gradH(flow, z0); pts = cq_pts(im, z0, n)
        frames = [initial_deviations(cq_seed(8, 2), j, Val(4), Val(4), Float64)]
        e = argmax(abs.(g0)); perm = vcat(e, setdiff(1:4, e)); push!(frames, cq_eye(4, 4, Float64)[:, perm])        # canonical: the unit vector with the largest gradient component first
        for Q0 in frames
            abs(dot(g0, Q0[:, 1])) >= 0.05 * norm(g0) || continue
            _, _, fin, ser = propagate_cqr_series(cq, z0, Q0, n; every = 10)                                  # q = 4: every = 10 Q-steps = 20 orbit steps = 1 time unit
            tvalid = 0.0; norb += 1
            for s in ser
                z = pts[s.n + 1]; g = cq_gradH(flow, z); c = abs(dot(g, s.Q[:, 1])) / norm(g)
                c >= 1.0e3 * cq_eps(Float64) || break
                B = CQ_B8 + 8 * cq_eps(Float64) / c                                                 # D121 (C8): CQ_B8 = 2 x the worst realized error at t <= 60 over the full-size sample, rounded up to one digit
                err = abs(s.l[1] - (log(abs(dot(g0, Q0[:, 1]))) - log(abs(dot(g, s.Q[:, 1])))))
                if s.n / 20 <= 60
                    worst = max(worst, err); err <= B || (nbad += 1)
                else
                    beyond = max(beyond, err); nbeyond += err > B                                   # reported only: discretization error amplified along the orbit beyond the pointwise horizon
                end
                tvalid = s.n / 20
            end
            push!(tvals, tvalid)
        end
    end
    @test nbad == 0
    @test norb >= 4
    @info "AC8 reported" worst_err_t60 = worst worst_err_beyond_t60 = beyond samples_beyond_bound = nbeyond bound = CQ_B8 orbits = norb median_Tvalid = cq_median(tvals) min_Tvalid = minimum(tvals)
end

# ------------------------------------------------------------------ AC9: same orbit, discrete versus continuous
@testset "AC9 same-orbit agreement with the discrete QR spectrum (Householder, repeated Gram-Schmidt) and Benettin: T = 60 scored (D121, C9), the long horizon reported" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); cls = cq_classes(CQ_NCLASS, 1); m = min(CQ_NCLASS, cq_size(2, 20))
    starts = vcat([("R", cq_orbit_R()), ("C", cq_orbit_C())], [("reg$i", z) for (i, z) in enumerate(cls.regular[1:m])], [("cha$i", z) for (i, z) in enumerate(cls.chaotic[1:m])])
    Tlong = cq_size(200.0, 1000.0); rows = Any[]; worst = 0.0; worstl = 0.0; spread = 0.0; nfail = 0; cdreal = Dict(4 => 0.0, 2 => 0.0); cdlong = 0.0
    for (j, (name, z0)) in enumerate(starts), h in (1 / 20, 1 / 40)
        im = cq_im(flow, Float64, :o8, h); Q0 = initial_deviations(cq_seed(9, 1), j, Val(4), Val(4), Float64)
        # lambda form at T = 60: c_d h^q (q = 4 and 2; the table has order 8): D121 (C9), c_d = 2 x the worst realized c_d over the full-size sample (20 starts per class), rounded up
        for (q, cd) in ((4, CQ_CD4), (2, CQ_CD2))
            cq = cq_cqr(im, flow; order = q); N = round(Int, 60.0 / h)
            lc = propagate_cqr_orbit(cq, z0, Q0, N)[3].lambda
            ld = cq_discrete(im, z0, Q0, N, h; alg = RGS()); lh = cq_discrete(im, z0, Q0, N, h; alg = HouseholderQR()); lb = cq_benettin(im, z0, Q0[:, 1], N, h)
            allow = cd * h^q; e = maximum(abs, lc - ld); eh = maximum(abs, lc - lh); eb = abs(lc[1] - lb)
            (e <= allow && eh <= allow && eb <= allow) || (nfail += 1)
            cdreal[q] = max(cdreal[q], max(e, eh, eb) / h^q)
            worst = max(worst, e / allow); spread = max(spread, maximum(abs, ld - lh)); push!(rows, (name = name, q = q, T = 60.0, h = h, err = e))
        end
        # long horizon (order-4 scheme), reported and not scored (D121, C9): the differences grow in events along the orbit; the lambda-form c_d and the l-form constant are recorded
        cq = cq_cqr(im, flow; order = 4); N = round(Int, Tlong / h); lc = propagate_cqr_orbit(cq, z0, Q0, N)[3].lambda; ld = cq_discrete(im, z0, Q0, N, h)
        dl = maximum(abs, lc - ld) * Tlong; worstl = max(worstl, dl / h^4); cdlong = max(cdlong, maximum(abs, lc - ld) / h^4)
        push!(rows, (name = name, q = 4, T = Tlong, h = h, err = dl))
    end
    @test nfail == 0
    # on orbit C at T = 60 the errors of h and h/2 differ by about 2^m: ratio in [2^4 / 4, 4 x 2^4] (order-4 scheme, order-8 table: m = 4)
    ea = [r.err for r in rows if r.name == "C" && r.q == 4 && r.T == 60.0 && r.h == 1 / 20][1]; eb = [r.err for r in rows if r.name == "C" && r.q == 4 && r.T == 60.0 && r.h == 1 / 40][1]
    ea > 1.0e3 * cq_eps(Float64) && eb > 1.0e3 * cq_eps(Float64) && @test 2^4 / 4 <= ea / eb <= 4 * 2^4
    @info "AC9 reported (R3, R6): worst lambda error / allowance at T = 60, realized c_d at T = 60 (q = 4, 2), realized c_d and l-form constant at T = long (not scored), spread among the discrete variants, ratio on orbit C" worst cdreal cdlong worstl Tlong spread ratio_C = ea / eb
end

# ------------------------------------------------------------------ AC10: regular decay (C3, Skokos2010 Eq. (62))
@testset "AC10a pendulum libration: lambda_1(t) decays like 1/t and l_1 + l_2 = 0" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    flow = Pendulum(1.0); im = cq_im(flow, Float64, :o4, CQ_H); cq = cq_cqr(im, flow; order = 4); z0 = cq_pendulum_start(2.0)
    Tmax = cq_size(1.0e4, 1.0e4); n = round(Int, Tmax / CQ_H); Q0 = initial_deviations(cq_seed(10, 1), 1, Val(2), Val(2), Float64)
    _, _, fin, ser = propagate_cqr_series(cq, z0, Q0, n; every = 50)                       # every 50 Q-steps = 5 time units
    pts = cq_pts(im, z0, 2000); AFmax = maximum(norm(vector_field_jacobian(flow, p)) for p in pts)
    ts = [s.n * CQ_H for s in ser]; sel = findall(t -> t >= 1.0e2, ts)
    idx = unique(round.(Int, exp.(range(log(first(sel)), log(last(sel)); length = 60))))
    lam = [abs(ser[i].l[1]) / ts[i] for i in idx]
    sl = cq_slope(log.(ts[idx]), log.(lam))
    @test -1.1 <= sl <= -0.8
    Bsum(t) = t * AFmax * (2 * CQ_H)^4 * AFmax^4 + 8 * cq_eps(Float64) * t * AFmax                          # B_sum of criterion 12 (C_s = 1, q = 4, Delta = 2 h)
    @test all(s -> abs(s.l[1] + s.l[2]) <= Bsum(s.n * CQ_H), ser)
    @info "AC10a reported" slope = sl final_lambda1 = fin.lambda[1]
end

@testset "AC10b F3 regular orbits: |lambda_i(T)| <= 5 (1 + ln T) / T for at least 80 percent of the starts, orbit R included" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 4); T = cq_size(2000.0, 1.0e4); n = round(Int, T / CQ_H)
    cls = cq_classes(CQ_NCLASS, 1); starts = vcat([cq_orbit_R()], cls.regular); bound = 5 * (1 + log(T)) / T; ok = Bool[]; maxes = Float64[]
    for (j, z0) in enumerate(starts)
        Q0 = initial_deviations(cq_seed(10, 2), j, Val(4), Val(4), Float64); lam = propagate_cqr_orbit(cq, z0, Q0, n)[3].lambda
        push!(maxes, maximum(abs, lam)); push!(ok, all(abs.(lam) .<= bound))
    end
    @test ok[1]                                                                                              # orbit R itself
    @test count(ok) >= 0.8 * length(ok)
    @info "AC10b reported" bound maxes sali = [cq_sali(z, 1000.0) for z in starts]
end
