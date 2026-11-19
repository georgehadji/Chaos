# AC11 to AC15: the chaotic class and the pairing, the sum identity, the orthogonality control, same-orbit rounding against a 256-bit reference, the precision ladder.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError). The two negative controls of AC12 and AC14 use no new name.

# ------------------------------------------------------------------ helpers (definitions only)
cq_AFmax(flow, pts) = maximum(norm(vector_field_jacobian(flow, p)) for p in pts)
cq_A2max(flow, pts) = maximum(opnorm(Matrix(vector_field_jacobian(flow, p))) for p in pts)
# The scheme of the spec in the type of `Q0` at the realized orbit points `pts` (any float type, widened exactly by the caller), l and Q at the end and at `checkpoints`.
function cq_ref_run(flow, pts, Q0b, hb, q::Integer, N::Integer; checkpoints = (), cnu = size(Q0b, 2) == 4 ? 1.0 : 2.0, form = size(Q0b, 2) == 4 ? :f84 : :f92, mode = :ok)
    Bt = eltype(Q0b)
    Afun(n) = vector_field_jacobian(flow, SVector{4, Bt}(map(cq_big, pts[n + 1])))
    return cq_scheme(Afun, Q0b, hb, q, N; cnu = cnu, form = form, mode = mode, checkpoints = checkpoints)
end
# Rounding bound B_C(t) of criterion 14 (derived worst case, proposal; L1 Amendment D35 form).
cq_BC(::Type{T}, AFmax, t, nQ) where {T} = 8 * cq_eps(T) * AFmax * t * (1 + nQ)

# ------------------------------------------------------------------ AC11
@testset "AC11 chaotic class: positive lambda_1, agreement with the discrete value, the zero pair and the pairing, independent ensembles" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 4); T = cq_size(500.0, 2000.0); n = round(Int, T / CQ_H); m = CQ_NCLASS
    A = cq_classes(m, 1).chaotic; B = cq_classes(m, 2).chaotic; allow = 5.0e4 * CQ_H^4 / T                   # l form of AC9 at the horizon T: |l^C - l^D| <= c_l h^4, c_l = 5e4 (U9)
    lc = Vector{Vector{Float64}}(); ld = Vector{Vector{Float64}}(); nzero = 0; npair = 0; npos = 0; ndisc = 0; dratio = Float64[]
    for (j, z0) in enumerate(A)
        Q0 = initial_deviations(cq_seed(11, 1), j, Val(4), Val(4), Float64)
        lam = propagate_cqr_orbit(cq, z0, Q0, n)[3].lambda; dis = cq_discrete(im, z0, Q0, n, CQ_H)
        push!(lc, lam); push!(ld, dis)
        lam[1] > 0 && (npos += 1)
        abs(lam[1] - dis[1]) <= allow && (ndisc += 1); push!(dratio, abs(lam[1] - dis[1]) / allow)
        (abs(lam[2]) <= 2.0e-2 && abs(lam[3]) <= 2.0e-2) && (nzero += 1)
        (abs(lam[1] + lam[4]) <= 2.0e-2 && abs(lam[2] + lam[3]) <= 2.0e-2) && (npair += 1)
    end
    @test npos == length(A) && ndisc >= 0.95 * length(A)                                                # D121 (C11): at least 95 percent of the chaotic starts within the allowance; the maximum ratio is reported
    @test nzero >= 0.95 * length(A) && npair >= 0.95 * length(A)
    # independent orbit sets: the module on the starts A against the discrete spectrum on the further starts B (L1 S5)
    lB = [cq_discrete(im, z0, initial_deviations(cq_seed(11, 2), j, Val(4), Val(4), Float64), n, CQ_H)[1] for (j, z0) in enumerate(B)]
    @test cq_agree([l[1] for l in lc], lB)
    @info "AC11 reported" max_ratio_to_allowance = maximum(dratio) within_allowance = ndisc mean_lambda1_cqr = cq_mean([l[1] for l in lc]) mean_lambda1_discrete_B = cq_mean(lB) max_zero_pair = maximum(max(abs(l[2]), abs(l[3])) for l in lc) max_pairing = maximum(max(abs(l[1] + l[4]), abs(l[2] + l[3])) for l in lc) discrete_zero = maximum(max(abs(l[2]), abs(l[3])) for l in ld)
end

# ------------------------------------------------------------------ AC12: sum identity
@testset "AC12 sum identity: |sum l_i| <= B_sum(T) and its order in Delta" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); T = cq_size(60.0, 100.0); worst = 0.0
    orders = Dict()
    for (name, z0) in (("R", cq_orbit_R()), ("C", cq_orbit_C())), q in (2, 4)
        sums = Float64[]
        for h in (1 / 20, 1 / 40)
            im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = q); N = round(Int, T / h)
            Q0 = initial_deviations(cq_seed(12, 1), 1, Val(4), Val(4), Float64)
            fin = propagate_cqr_orbit(cq, z0, Q0, N)[3]; AF = cq_AFmax(flow, cq_pts(im, z0, N)); Delta = (q / 2) * h
            Bs = T * AF * Delta^q * AF^q + 8 * cq_eps(Float64) * T * AF; s = abs(sum(fin.l))
            @test s <= Bs; worst = max(worst, s / Bs); push!(sums, s)
        end
        orders[(name, q)] = sums[1] / sums[2]                         # D121 (C12, owner "Keep M2, amend"): the order of the sum in Delta is reported, not scored
    end
    @info "AC12 reported (R3): worst |sum| / B_sum; sum ratio h = 1/20 to 1/40 per (orbit, q)" worst orders
end

@testset "AC12 negative control: a double that restores Q to its polar factor each step passes the sum bound without any penalty (test-side, no new name; spec point U5)" begin
    flow, im = cq_hh(); z0 = cq_orbit_C(); N = 400; pts = cq_pts(im, z0, N); Q0 = initial_deviations(cq_seed(12, 2), 1, Val(4), Val(4), Float64); t = N * CQ_H
    Afun(n) = vector_field_jacobian(flow, pts[n + 1]); AF = cq_AFmax(flow, pts)
    pol = cq_scheme(Afun, Q0, CQ_H, 4, N; mode = :polar, cnu = 0.0); ok = cq_scheme(Afun, Q0, CQ_H, 4, N)
    Bs = t * AF * (2 * CQ_H)^4 * AF^4 + 8 * cq_eps(Float64) * t * AF
    @test abs(sum(pol.l)) <= Bs && abs(sum(ok.l)) <= Bs         # the sum alone does not test the penalty (the stage values carry a defect of order Delta^2, C15: not within rounding)
    @info "AC12 control reported" polar_sum = sum(pol.l) scheme_sum = sum(ok.l) bound = Bs
end

# ------------------------------------------------------------------ AC13: orthogonality control
cq_sym_delta(rng, d, size_) = (x = randn(rng, d, d); x = (x + x') / 2; x * (size_ / maximum(abs, x)))
@testset "AC13a injected defect decays: coupled hyperbolic fixture (form 84, k = d) and two columns of orbit C (form 92, k = 2)" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    rng = StableRNG(cq_seed(13, 1)); worst = 0.0
    for q in (2, 4)
        # form 84 on the coupled fixture
        A = cq_coupled([0.8, 0.35], CQ_B); h = 0.05; T = Float64; cq = cq_est(A, T, h, q; orthogonality_tol = 0.1); AF = cq_AF(A)
        Qc = cq_frame(rng, 4, 4, T); Qi = SMatrix{4, 4, T}(Qc * (I + SMatrix{4, 4, Float64}(cq_sym_delta(rng, 4, 1.0e-3)))); E0 = cq_defect(Qi)
        nu = 1.0 * AF; Eeq = propagate_cqr_orbit(cq, cq_x0(4, T), Qc, 100)[3].Emax
        ser = propagate_cqr_series(cq, cq_x0(4, T), Qi, 100; every = 1)[4]
        for s in ser
            t = s.n * h; bound = 2E0 * exp(-2 * nu * (1 - 2E0) * t) + Eeq
            @test s.E <= bound; worst = max(worst, s.E / bound)
        end
        # form 92 on the first two columns of orbit C
        flow, im = cq_hh(); cqc = cq_cqr(im, flow; order = q, orthogonality_tol = 0.1); z0 = cq_orbit_C(); pts = cq_pts(im, z0, 100)
        Q2 = SMatrix{4, 2, Float64}(initial_deviations(cq_seed(13, 2), 1, Val(4), Val(4), Float64)[:, 1:2]); d2 = cq_sym_delta(rng, 2, 1.0e-3)
        Q2i = SMatrix{4, 2, Float64}(Q2 * (I + d2)); E02 = cq_defect(Q2i)
        numin = 2.0 * minimum(cq_AF(vector_field_jacobian(flow, p)) for p in pts); A2 = cq_A2max(flow, pts)
        Eeq2 = propagate_cqr_orbit(cqc, z0, Q2, 100)[3].Emax; ser2 = propagate_cqr_series(cqc, z0, Q2i, 100; every = 1)[4]
        for s in ser2
            t = s.n * CQ_H; bound = 2E02 * exp(-2 * (numin * (1 - 2E02) - (1 + 2E02) * A2) * t) + Eeq2
            @test s.E <= bound; worst = max(worst, s.E / bound)
        end
    end
    @info "AC13a reported: worst E / bound" worst
end

@testset "AC13b equilibrium defect: E_max seen is below the tolerance and scales like Delta^q" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); z0 = cq_orbit_C(); T = cq_size(60.0, 100.0); rec = Any[]
    for q in (2, 4)
        Es = Float64[]
        for h in (1 / 20, 1 / 40)
            im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = q); N = round(Int, T / h)
            Q0 = initial_deviations(cq_seed(13, 3), 1, Val(4), Val(4), Float64)
            E = propagate_cqr_orbit(cq, z0, Q0, N)[3].Emax; @test E <= (q == 4 ? CQ_EMAX : 5 * CQ_EMAX); push!(Es, E)
        end
        Es[1] > 1.0e3 * cq_eps(Float64) && Es[2] > 1.0e3 * cq_eps(Float64) && @test 2^q / 4 <= Es[1] / Es[2] <= 4 * 2^q
        push!(rec, (q = q, E_eq = Es, ratio = Es[1] / Es[2]))
    end
    @info "AC13b reported (R2)" rec
end

@testset "AC13c injected fault: Q scaled by 1.05 (E = 0.1) at step 500 raises the orthogonality error at the next Q-step, nothing stored" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr)
    flow, im = cq_hh(); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(13, 4), 1, Val(4), Val(4), Float64)
    for q in (2, 4)
        cq = cq_cqr(im, flow; order = q)
        x, Q, _ = propagate_cqr_orbit(cq, z0, Q0, 500)                                       # the run is healthy up to step 500
        e = cq_health_error(() -> propagate_cqr_orbit(cq, x, 1.05 * Q, 20; orbit = 7))      # the restart at the faulty frame counts steps from its own start
        @test e isa NumericalHealthError && e.what == "orthogonality" && e.orbit == 7 && e.step == q ÷ 2 && e.value > CQ_EMAX
        e2 = cq_health_error(() -> propagate_cqr(cq, [z0, x, z0], 20; k = 4, seed = 1, Q0 = [Q0, 1.05 * Q, Q0]))
        @test e2 isa NumericalHealthError && e2.what == "orthogonality" && e2.orbit == 2 && e2.step == q ÷ 2        # nothing stored: the call throws
    end
end

@testset "AC13d explicit stability limits of C14: E diverges above nu Delta = 1 (order 2) and 1.39 (order 4), not at 0.9 of it" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    A = cq_diag2(); h = 0.05; T = Float64; AF = cq_AF(A); rng = StableRNG(cq_seed(13, 5))
    for (q, lim) in ((2, 1.0), (4, 1.39))
        Delta = (q / 2) * h; cnu_lim = lim / Delta / AF
        Q0 = SMatrix{2, 2, T}(cq_frame(rng, 2, 2, T) * (I + SMatrix{2, 2, Float64}(cq_sym_delta(rng, 2, 1.0e-4))))
        above = cq_health_error(() -> propagate_cqr_orbit(cq_est(A, T, h, q; penalty_factor = 1.1 * cnu_lim, orthogonality_tol = 0.1), cq_x0(2, T), Q0, 600))
        @test above isa NumericalHealthError && above.what == "orthogonality"
        below = propagate_cqr_orbit(cq_est(A, T, h, q; penalty_factor = 0.9 * cnu_lim, orthogonality_tol = 0.1), cq_x0(2, T), Q0, 600)[3]
        @test below.E <= cq_defect(Q0)                                           # the injected defect decays (it rises transiently to about 2e-2 x tol, reported)
    end
end

# ------------------------------------------------------------------ AC14: same-orbit rounding against the scheme in 256 bits
@testset "AC14 same-orbit rounding: |l_i - l_i,ref| <= B_C(t) at every 50th Q-step, Float64 and Double64, reference in BigFloat(256) at the realized points" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    flow = HenonHeilesFlow(); cls = cq_classes(CQ_NCLASS, 1); m = min(CQ_NCLASS, cq_size(1, 20)); t_end = cq_size(60.0, 1.0e3)
    starts = vcat([("R", cq_orbit_R()), ("C", cq_orbit_C())], [("reg$i", z) for (i, z) in enumerate(cls.regular[1:m])], [("cha$i", z) for (i, z) in enumerate(cls.chaotic[1:m])])
    worst = Dict(Float64 => 0.0, Double64 => 0.0)
    for (j, (name, z0)) in enumerate(starts), T in (Float64, Double64)
        im = cq_im(flow, T); cq = cq_cqr(im, flow; order = 4); x0 = SVector{4, T}(map(T, z0)); N = round(Int, t_end / CQ_H)
        Q0 = initial_deviations(cq_seed(14, 1), j, Val(4), Val(4), T); pts = cq_pts(im, x0, N)
        ser = propagate_cqr_series(cq, x0, Q0, N; every = 50)[4]; AF = cq_AFmax(flow, pts)
        setprecision(BigFloat, CQ_BITS) do
            Q0b = SMatrix{4, 4, BigFloat}(map(cq_big, Q0)); hb = cq_big(im.tau)
            ref = cq_ref_run(flow, pts, Q0b, hb, 4, N; checkpoints = Set(s.n for s in ser))
            for s in ser, i in 1:4
                err = Float64(abs(cq_big(s.l[i]) - ref.cps[s.n].l[i])); B = cq_BC(T, AF, s.n * Float64(im.tau), s.n ÷ 2)
                @test err <= B
                worst[T] = max(worst[T], err / B)
            end
        end
    end
    @info "AC14 reported (R3): worst error / B_C" worst
end

@testset "AC14 negative control: a double with the product Q^T A Q in single precision exceeds B_C (test-side, no new name)" begin
    flow, im = cq_hh(); z0 = cq_orbit_C(); N = cq_size(1200, 20000); pts = cq_pts(im, z0, N); T = Float64
    Q0 = initial_deviations(cq_seed(14, 2), 1, Val(4), Val(4), T); Afun(n) = vector_field_jacobian(flow, pts[n + 1]); AF = cq_AFmax(flow, pts)
    f32 = cq_scheme(Afun, Q0, T(CQ_H), 4, N; mode = :f32); ok = cq_scheme(Afun, Q0, T(CQ_H), 4, N)
    B = cq_BC(T, AF, N * CQ_H, N ÷ 2)
    @test maximum(abs, f32.l - ok.l) > B
    @info "AC14 control reported" f32_error = maximum(abs, f32.l - ok.l) B
end

# ------------------------------------------------------------------ AC15: precision ladder
@testset "AC15 precision ladder: Float64, Double64 and BigFloat(256) agree; BigFloat against its own 512-bit reference" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); t_end = cq_size(20.0, 60.0); N = round(Int, t_end / CQ_H); m = min(CQ_NCLASS, cq_size(5, 50))
    starts = vcat([cq_orbit_C()], cq_classes(CQ_NCLASS, 1).chaotic[1:min(m, CQ_NCLASS)])
    vals = Dict(Float64 => Float64[], Double64 => Float64[], BigFloat => Float64[]); gaps = 0; nagree = 0; Wb = 0.0
    fin = Dict{Any, Any}()
    for (j, z0) in enumerate(starts), T in (Float64, Double64)
        im = cq_im(flow, T); cq = cq_cqr(im, flow; order = 4); x0 = SVector{4, T}(map(T, z0)); Q0 = initial_deviations(cq_seed(15, 1), j, Val(4), Val(4), T)
        r = propagate_cqr_orbit(cq, x0, Q0, N); fin[(j, T)] = r; push!(vals[T], Float64(r[3].l[1]))
    end
    for j in eachindex(starts)
        xa, xb = fin[(j, Float64)][1], fin[(j, Double64)][1]
        maximum(abs.(Float64.(xb) .- Float64.(xa))) <= 1.0e-9 || continue
        nagree += 1
        im = cq_im(flow, Float64); AF = cq_AFmax(flow, cq_pts(im, SVector{4, Float64}(starts[j]), N))
        gap = maximum(abs, fin[(j, Float64)][3].l - Float64.(fin[(j, Double64)][3].l)); gap <= cq_BC(Float64, AF, t_end, N ÷ 2) + cq_BC(Double64, AF, t_end, N ÷ 2) + 1.0e-9 || (gaps += 1)
    end
    @test gaps == 0 && nagree >= 1
    setprecision(BigFloat, CQ_BITS) do
        for (j, z0) in enumerate(starts[1:min(end, cq_size(2, 6))])
            im = cq_im(flow, BigFloat); cq = cq_cqr(im, flow; order = 4); x0 = SVector{4, BigFloat}(map(BigFloat, z0)); Q0 = initial_deviations(cq_seed(15, 1), j, Val(4), Val(4), BigFloat)
            r = propagate_cqr_orbit(cq, x0, Q0, N); push!(vals[BigFloat], Float64(r[3].l[1])); pts = cq_pts(im, x0, N)
            AF = cq_AFmax(flow, pts); l512 = setprecision(BigFloat, 2 * CQ_BITS) do
                cq_ref_run(flow, pts, SMatrix{4, 4, BigFloat}(map(cq_big, Q0)), cq_big(im.tau), 4, N).l
            end
            err = maximum(abs, [Float64(abs(r[3].l[i] - l512[i])) for i in 1:4]); B = cq_BC(BigFloat, AF, t_end, N ÷ 2)
            @test err <= B
            Wb = max(Wb, err / B)
        end
    end
    @test cq_agree(vals[Float64], vals[Double64])
    nb = length(vals[BigFloat]); @test cq_agree(vals[Float64][1:nb], vals[BigFloat])
    @info "AC15 reported" orbits_in_agreement = nagree bigfloat_error_over_bound = Wb
end
