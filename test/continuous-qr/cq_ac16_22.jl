# AC16 to AC22 and the recorded lines R1, R2, R5: orbit and execution invariance, cost counts and memory, health and the run record, the external-reference stand-in,
# the CPU part of the CPU-GPU comparison. Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).
# Criterion 21 is the verification environment only (ChaosTools); the device part of 22 is local GPU; the package carries a Double64 stand-in for 20 and the CPU part of 22.

# ------------------------------------------------------------------ AC16: orbit and chain invariance
@testset "AC16 orbit invariance: the orbit is bitwise that of the integrator alone, of the discrete spectrum and of Benettin, for every k and both schemes" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr)
    flow, im = cq_hh(); z0 = cq_orbit_C(); N = cq_size(400, 2000); Q4 = initial_deviations(cq_seed(16, 1), 1, Val(4), Val(4), Float64)
    xref = advance_n(im.method, flow, z0, im.tau, N)
    for q in (2, 4), k in (1, 2, 4)
        cq = cq_cqr(im, flow; order = q); Qk = SMatrix{4, k, Float64}(Q4[:, 1:k])
        x = propagate_cqr_orbit(cq, z0, Qk, N)[1]
        xd = propagate(im, [z0], N; k = k, tau = 1, seed = 1, reducers = (QRSpectrum(RGS()),), W0 = [Qk]).x[1]
        @test x == xref && x == xd && x == last(cq_pts(im, z0, N))
        k == 1 && @test x == propagate(im, [z0], N; k = 1, tau = 1, seed = 1, reducers = (Benettin(),), W0 = [Qk]).x[1]
    end
    # the k = 1 series does not depend on other estimators living in the same program
    cq = cq_cqr(im, flow; order = 4); q1 = SMatrix{4, 1, Float64}(Q4[:, 1]); a = propagate_cqr_orbit(cq, z0, q1, N)[3]
    cq_discrete(im, z0, Q4, N, CQ_H); cq_benettin(im, z0, Q4[:, 1], N, CQ_H); propagate_cqr_orbit(cq, z0, Q4, N)
    b = propagate_cqr_orbit(cq, z0, q1, N)[3]; c = propagate_cqr(cq, [z0], N; k = 1, seed = 1, Q0 = [q1]).finals[1]
    @test a.l == b.l && a.l == c.l
end

# ------------------------------------------------------------------ AC17: execution invariance, allocation
@testset "AC17a execution invariance: chunk sizes and the serial body give bitwise identical series, sums, final Q and orbit points" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr)
    flow, im = cq_hh(); n = cq_size(40, 4096); N = cq_size(1200, 20_000); seed = cq_seed(17, 1)
    a, b = section_interval(flow, CQ_E); w = sqrt(2 * CQ_E); rng = StableRNG(cq_seed(17, 2))
    xs = [random_energy_start(flow, CQ_E, (a, b), (-w, w), rng)[1] for _ in 1:n]
    @info "AC17 reported" threads = Threads.nthreads() starts = n
    for q in (2, 4)
        cq = cq_cqr(im, flow; order = q)
        serial = [propagate_cqr_orbit(cq, xs[j], initial_deviations(seed, j, Val(4), Val(4), Float64), N; orbit = j) for j in 1:n]
        for chunk in (cq_size(7, 64), cq_size(32, 1000), nothing)
            res = propagate_cqr(cq, xs, N; k = 4, seed = seed, chunk = chunk)
            @test all(j -> isequal(res.finals[j], serial[j][3]) && res.x[j] == serial[j][1] && res.W[j] == serial[j][2], 1:n)
        end
    end
end

@testset "AC17b the hot loop allocates nothing after warm-up" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); z0 = cq_orbit_C()
    for T in (Float64, Double64), q in (2, 4), k in (1, 2, 4)
        im = cq_im(flow, T); cq = cq_cqr(im, flow; order = q); x0 = SVector{4, T}(map(T, z0))
        Q0 = SMatrix{4, k, T}(initial_deviations(cq_seed(17, 3), 1, Val(4), Val(4), T)[:, 1:k])
        @test cq_alloc_bytes(propagate_cqr_orbit, cq, x0, Q0, 2000) == 0
    end
end

@testset "AC17c the series is a bitwise prefix of the run: the record at step n equals a fresh run to n" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    flow, im = cq_hh(); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(17, 4), 1, Val(4), Val(4), Float64)
    for q in (2, 4)
        cq = cq_cqr(im, flow; order = q); ser = propagate_cqr_series(cq, z0, Q0, 200; every = 10)[4]
        @test length(ser) == 200 ÷ (q ÷ 2) ÷ 10
        for s in ser[[1, end ÷ 2, end]]
            r = propagate_cqr_orbit(cq, z0, Q0, s.n)
            @test r[3].l == s.l && r[2] == s.Q
        end
    end
end

# ------------------------------------------------------------------ AC18: cost counts and memory
@testset "AC18 cost counts: one Jacobian and one orbit step per orbit step, two right-hand sides per orbit step, memory formula of C21" begin
    _ = (ContinuousQR, propagate_cqr_orbit, cqr_bytes)
    flow, im = cq_hh(); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(18, 1), 1, Val(4), Val(4), Float64)
    for q in (2, 4)
        cf = cq_CountFlow(flow); cm = cq_CountMap(im); cq = cq_cqr(cm, cf; order = q)
        r1 = propagate_cqr_orbit(cq, z0, Q0, 500)[3]; j1 = cf.n[]; m1 = cm.n[]
        r2 = propagate_cqr_orbit(cq, z0, Q0, 1500)[3]; j2 = cf.n[] - j1; m2 = cm.n[] - m1
        @test j1 == 501 && j2 == 1501 && m1 == 500 && m2 == 1500           # the initial Jacobian plus exactly one per orbit step; exactly one orbit step per orbit step
        @test r2.njac - r1.njac == 1000 && r2.nrhs - r1.nrhs == 2 * 1000
        @test r1.nrhs == 2 * 500
    end
    # C21: sizeof(T) (d + 4 d k + 2 k + 3 N^2) + 16 at d = k = 4, N = 2 (order 4): 720 bytes in Float64
    cq4 = cq_cqr(im, flow; order = 4)
    @test cqr_bytes(cq4, Val(4), Val(4), Float64) == 8 * (4 + 4 * 4 * 4 + 2 * 4 + 3 * 2^2) + 16 == 720
    @test cqr_bytes(cq4, Val(4), Val(4), Double64) == 16 * (4 + 4 * 4 * 4 + 2 * 4 + 3 * 2^2) + 16
    @test cqr_bytes(cq4, Val(4), Val(2), Float64) == 8 * (4 + 4 * 4 * 2 + 2 * 2 + 3 * 2^2) + 16
end

# ------------------------------------------------------------------ AC19: health and run record
@testset "AC19a an equilibrium start is accepted (the module does not divide by the field)" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HarmonicOscillator([1.0]); im = cq_im(flow, Float64, :o4, 0.05); cq = cq_cqr(im, flow; order = 4); z0 = SVector(0.0, 0.0)
    x, Q, fin = propagate_cqr_orbit(cq, z0, cq_eye(2, 2, Float64), 100)
    @test x == z0 && all(isfinite, fin.l) && all(i -> abs(fin.l[i]) <= 8 * cq_eps(Float64) * 5, 1:2)
end

@testset "AC19b non-finite state, Jacobian or frame: kind, orbit id and step, nothing stored" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr)
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 2); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(19, 1), 1, Val(4), Val(4), Float64)
    zs = [z0, z0, SVector(NaN, z0[2], z0[3], z0[4])]
    e = cq_health_error(() -> propagate_cqr(cq, zs, 10; k = 4, seed = 1)); @test e isa NumericalHealthError && e.what == "state" && e.orbit == 3 && e.step == 0
    e = cq_health_error(() -> propagate_cqr(cq, [z0, z0], 10; k = 4, seed = 1, Q0 = [Q0, SMatrix{4, 4, Float64}(fill(NaN, 16))]))
    @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 2 && e.step == 0
    # a Jacobian that turns non-finite at its 10th call (call 1 is the initial one, call j + 1 is the orbit step j, order 2): spec point U4
    cqn = cq_cqr(im, cq_NaNJac(flow, 10); order = 2)
    e = cq_health_error(() -> propagate_cqr_orbit(cqn, z0, Q0, 40; orbit = 5)); @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 5 && e.step == 9
end

@testset "AC19c, AC19d argument errors come before any work" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr)
    flow, im = cq_hh(); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(19, 2), 1, Val(4), Val(4), Float64)
    cf = cq_CountFlow(flow); cq4 = cq_cqr(im, cf; order = 4); cq2 = cq_cqr(im, cf; order = 2)
    @test_throws ArgumentError propagate_cqr_orbit(cq4, z0, Q0, 7)                       # odd number of orbit steps with the order-4 scheme
    @test_throws ArgumentError propagate_cqr_orbit(cq4, z0, Q0, 0)
    @test_throws ArgumentError propagate_cqr_orbit(cq2, z0, Q0, -3)
    @test_throws ArgumentError propagate_cqr(cq4, [z0], 7; k = 4, seed = 1)
    @test_throws ArgumentError propagate_cqr(cq4, [z0], 8; k = 0, seed = 1)
    @test_throws ArgumentError propagate_cqr(cq4, [z0], 8; k = 5, seed = 1)
    @test_throws ArgumentError propagate_cqr_orbit(cq4, z0, SMatrix{4, 5, Float64}(ones(20)), 8)
    @test_throws ArgumentError ContinuousQR(im, flow; order = 3)
    @test cf.n[] == 0                                                                    # no Jacobian was evaluated
end

@testset "AC19e the run record carries h, the integrator and its order, q, the form, c_nu, E_max, the element type, k, the seed, backend, device, fallback reason and the flow" begin
    _ = (ContinuousQR, propagate_cqr)
    flow, im = cq_hh(); z0 = cq_orbit_C()
    r4 = propagate_cqr(ContinuousQR(im, flow; order = 4), [z0], 40; k = 4, seed = 77); p = r4.provenance
    @test p.cqr_h == CQ_H && p.cqr_integrator == method_name(im.method) && p.cqr_integrator_order == 8 && p.cqr_order == 4 && p.cqr_form === :f84
    @test p.cqr_penalty_factor == 1.0 && p.cqr_orthogonality_tol == CQ_EMAX && p.precision == "Float64" && p.k == 4 && p.seed == 77
    @test p.backend === :cpu && p.device isa AbstractString && p.backend_fallback_reason isa AbstractString && p.cqr_flow == sprint(show, flow)
    p2 = propagate_cqr(ContinuousQR(im, flow; order = 2, penalty_factor = 3.0, orthogonality_tol = 0.01), [z0], 40; k = 2, seed = 5).provenance
    @test p2.cqr_order == 2 && p2.cqr_form === :f92 && p2.cqr_penalty_factor == 3.0 && p2.cqr_orthogonality_tol == 0.01 && p2.k == 2
    @test propagate_cqr(ContinuousQR(cq_im(flow, Double64), flow), [SVector{4, Double64}(map(Double64, z0))], 40; k = 4, seed = 1).provenance.precision == "Double64"
    @test_throws ArgumentError propagate_cqr(ContinuousQR(im, flow; form = :f84), [z0], 40; k = 2, seed = 1)       # form 84 needs k = d (P1)
end

# ------------------------------------------------------------------ AC20: external reference stand-in
@testset "AC20 external variational reference, package stand-in (Double64 order-10 table at 1/80; the Taylor form is stage 04)" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    h = CQ_H; rng = StableRNG(cq_seed(20, 1)); worst = Dict(2 => 0.0, 4 => 0.0); rows = Any[]
    cases = (("F3 R", HenonHeilesFlow(), HenonHeilesFlow(), cq_orbit_R(), :o8, 4, (10, 30, 60)), ("F3 C", HenonHeilesFlow(), HenonHeilesFlow(), cq_orbit_C(), :o8, 4, (10, 30, 60)),
             ("pend 2", Pendulum(1.0), Pendulum(Double64(1)), cq_pendulum_start(2.0), :o4, 2, (10, 30, 60, 100)), ("pend 3", Pendulum(1.0), Pendulum(Double64(1)), cq_pendulum_start(3.0), :o4, 2, (10, 30, 60, 100)))
    for (name, flow, flowD, z0, key, k, ts) in cases
        d = dimension(flow); Q0 = cq_frame(rng, d, k, Float64); ref = cq_refl(flowD, z0, Q0, collect(ts)); im = cq_im(flow, Float64, key, h)
        for q in (2, 4)
            cq = cq_cqr(im, flow; order = q)
            for t in ts
                r = cq_final(cq, z0, Q0, round(Int, t / h)); e = maximum(abs, r.fin.l - ref[t]); B = q == 4 ? (name == "pend 3" ? 1.0e-1 : 5.0e-3) : 1.0
                (q == 2 && name == "pend 3") || @test e <= B; worst[q] = max(worst[q], e / B)
                ld = cq_discrete(im, z0, Q0, round(Int, t / h), h) * t; push!(rows, (name = name, q = q, t = t, err = e, discrete_err = maximum(abs, ld - ref[t])))
            end
        end
    end
    @info "AC20 reported (R3, R4): worst error / delta_ext per scheme, and the error of the discrete spectrum against the same reference" worst rows
end

# ------------------------------------------------------------------ AC22: CPU part of the CPU-GPU comparison
@testset "AC22 CPU part: each backend against its own 256-bit reference, the gap where the orbits agree; the :gpu request falls back to the CPU with its reason" begin
    _ = (ContinuousQR, propagate_cqr)
    flow = HenonHeilesFlow(); m = CQ_NCLASS; t_end = cq_size(60.0, 400.0); N = round(Int, t_end / CQ_H); cls = cq_classes(m, 1); nbad = 0; ngap = 0; nagree = 0; ndiv = 0; divs = Float64[]
    for (cls_name, zs) in (("regular", cls.regular), ("chaotic", cls.chaotic))
        fin = Dict{Any, Any}(); bound = Dict{Any, Any}()
        for T in (Float64, Double64)
            im = cq_im(flow, T); cq = cq_cqr(im, flow; order = 4); xs = [SVector{4, T}(map(T, z)) for z in zs]
            Q0 = [initial_deviations(cq_seed(22, 1), j, Val(4), Val(4), T) for j in eachindex(zs)]
            res = propagate_cqr(cq, xs, N; k = 4, seed = 1, Q0 = Q0, backend = :cpu); @test res.provenance.backend === :cpu; fin[T] = res
            for j in eachindex(zs)
                pts = cq_pts(im, xs[j], N); AF = cq_AFmax(flow, pts); B = cq_BC(T, AF, t_end, N ÷ 2); bound[(T, j)] = B
                err = setprecision(BigFloat, CQ_BITS) do
                    ref = cq_ref_run(flow, pts, SMatrix{4, 4, BigFloat}(map(cq_big, Q0[j])), cq_big(im.tau), 4, N).l
                    maximum(i -> Float64(abs(cq_big(res.finals[j].l[i]) - ref[i])), 1:4)
                end
                err <= B || (nbad += 1)
            end
        end
        for j in eachindex(zs)
            dx = maximum(abs.(Float64.(fin[Double64].x[j]) .- Float64.(fin[Float64].x[j]))); push!(divs, dx)
            gap = maximum(abs, fin[Float64].finals[j].l - Float64.(fin[Double64].finals[j].l))
            if dx <= 1.0e-9
                nagree += 1
                cls_name == "regular" && !(gap <= bound[(Float64, j)] + bound[(Double64, j)] + 1.0e-9) && (ngap += 1)       # chaotic class: reported, not gated (P10)
            end
        end
    end
    @test nbad == 0 && ngap == 0
    @info "AC22 reported: orbits in agreement, orbit divergence" nagree divs
    if Chaoslab.Foundation.DEVICE[] isa Chaoslab.Foundation.NoDevice
        im = cq_im(flow, Float64); cq = cq_cqr(im, flow; order = 4); xs = [SVector{4, Float64}(z) for z in cls.regular[1:1]]
        r1 = propagate_cqr(cq, xs, 40; k = 4, seed = 1, backend = :gpu); r0 = propagate_cqr(cq, xs, 40; k = 4, seed = 1, backend = :cpu)
        @test r1.provenance.backend === :cpu && !isempty(string(r1.provenance.backend_fallback_reason)) && isequal(r1.finals, r0.finals)
    end
end

# ------------------------------------------------------------------ recorded lines R1, R2, R5 (reported, not scored)
@testset "R1 cost ratio against the discrete QR, reported (the 40 of Geist1990 Table III p.886 is a reproduction target, never a pass)" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow, im = cq_hh(); z0 = cq_orbit_C(); N = cq_size(1200, 20_000); Q0 = initial_deviations(cq_seed(21, 1), 1, Val(4), Val(4), Float64)
    tq = Dict{Any, Float64}()
    for q in (2, 4)
        cq = cq_cqr(im, flow; order = q); propagate_cqr_orbit(cq, z0, Q0, 40); tq[q] = @elapsed propagate_cqr_orbit(cq, z0, Q0, N)
    end
    cq_discrete(im, z0, Q0, 40, CQ_H); td = @elapsed cq_discrete(im, z0, Q0, N, CQ_H); th = @elapsed cq_discrete(im, z0, Q0, N, CQ_H; alg = HouseholderQR())
    @test all(isfinite, values(tq)) && isfinite(td) && isfinite(th)
    @info "R1 reported: seconds and ratios to the repeated Gram-Schmidt and Householder spectra" tq discrete_rgs = td discrete_householder = th ratio_q4_rgs = tq[4] / td ratio_q2_rgs = tq[2] / td
end

@testset "R2 equilibrium defect against h, c_nu and the form, and the exponent bias against it, reported" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    flow = HenonHeilesFlow(); z0 = cq_orbit_C(); Q0 = initial_deviations(cq_seed(21, 2), 1, Val(4), Val(4), Float64); rows = Any[]; T = 60.0
    for (k, form) in ((4, :f84), (2, :f92)), h in (1 / 20, 1 / 40), cnu in (1.0, 2.0, 4.0)
        im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = 4, penalty_factor = cnu, form = form, orthogonality_tol = 0.1)
        r = propagate_cqr_orbit(cq, z0, SMatrix{4, k, Float64}(Q0[:, 1:k]), round(Int, T / h))[3]; push!(rows, (k = k, form = form, h = h, cnu = cnu, E = r.Emax, l1 = r.l[1]))
    end
    @test all(r -> isfinite(r.E) && r.E <= 0.1, rows)
    @info "R2 reported" rows
end

@testset "R5 form 92 at nu = 0 and k = 2 on orbit C: the time at which E reaches 1e-2, with and without the penalty (the claim of Skokos2010 p.53 is tested, not assumed)" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    flow, im = cq_hh(); z0 = cq_orbit_C(); n = cq_size(4000, 20_000); Q0 = SMatrix{4, 2, Float64}(initial_deviations(cq_seed(21, 3), 1, Val(4), Val(4), Float64)[:, 1:2]); tcross = Dict{Float64, Float64}()
    for cnu in (0.0, 2.0)
        cq = cq_cqr(im, flow; order = 4, form = :f92, penalty_factor = cnu, orthogonality_tol = 1.0e6)
        ser = propagate_cqr_series(cq, z0, Q0, n; every = 5)[4]; i = findfirst(s -> s.E >= 1.0e-2, ser); tcross[cnu] = i === nothing ? Inf : ser[i].n * CQ_H
    end
    @test tcross[2.0] == Inf                                                              # with the penalty the defect never reaches 1e-2
    @info "R5 reported (if the no-penalty time is Inf, form 92 did not grow at nu = 0 and the finding goes to the owner)" tcross
end
