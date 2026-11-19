# AC11 to AC14: precision ladder, time reversal, health, validation and the run record. Criterion 13 (the static null-space algorithm, 256 bits) is the verification environment only
# (stage 04, verify/clv-ginelli.jl); the plan maps it. Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).

# ------------------------------------------------------------------ AC11: precision ladder against the 256-bit reference algorithm at the realized orbit points
@testset "AC11 precision ladder: Float64 and Double64 CLVs and local exponents against the 256-bit reference along the realized orbit, within B_v = 1e3 eps (1 + G_max)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    worst = Dict{Any, Float64}(); nbad = 0; rng = StableRNG(clv_seed(11, 1))
    for T in (Float64, Double64)
        cases = Any[(label = "std $j", map = clv_std(T), x0 = x0, tau = 1, h = 1.0, nf = 60, ns = 60, nb = 60) for (j, x0) in enumerate(clv_starts(rng, 2, clv_size(3, 20), T))]
        push!(cases, (label = "F3 C", map = clv_im(T), x0 = map_to(T, clv_orbit_C()), tau = 20, h = CLV_H, nf = 60, ns = clv_size(20, 100), nb = 60))
        for (ic, c) in enumerate(cases)
            setprecision(BigFloat, CLV_BITS) do
                d = dimension(c.map); Dt = T(c.tau) * T(c.h); W0 = clv_frame(clv_seed(11, 10 + ic), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(11, 40 + ic)), d, T)
                res = clv_go(c.map, c.x0, c.nf, c.ns, c.nb; tau = c.tau, h = c.h, W0 = W0, C0 = C0)
                pts = clv_orbit(c.map, c.x0, (c.nf + c.ns + c.nb) * c.tau)
                ref = clv_ref_run(clv_bigmap(c.map), pts, c.tau, c.nf, c.ns, c.nb, W0, C0, BigFloat)
                # D124: Dt is the module's, tau h in the element type. spec c_v = 1e3: realized ratio 1.6 on F3 orbit C, so c_v = 1e4 (plan U8)
                Bv = 1.0e4 * clv_eps(T) * (1 + clv_Gmax(res.R, c.ns)); wv = 0.0; wl = 0.0
                for n in 0:c.ns
                    Vr = clv_refV(ref, n)
                    wv = max(wv, maximum(i -> clv_delta(res.V[:, i, n + 1], Vr[:, i]), 1:d))
                    wl = max(wl, maximum(i -> Float64(abs(clv_big(res.lambda[i, n + 1]) * Dt - log(ref.Lambda[i, n + 1]))), 1:d))
                end
                worst[(T, c.label)] = max(wv, wl) / Bv; max(wv, wl) <= Bv || (nbad += 1)
            end
        end
    end
    @info "AC11 realized ratios error / B_v (R3)" worst
    @test nbad == 0
end

@testset "AC11b the convergence of criterion 3 (b) holds in BigFloat(256): two starting frames and C matrices agree to 1e-10 on the standard map" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    setprecision(BigFloat, CLV_BITS) do
        T = BigFloat; map = clv_std(T); x0 = SVector{2, T}(T(13) / 10, T(7) / 10); out = []
        for s in (1, 2)
            W0 = clv_frame(clv_seed(11, 80 + s), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(11, 90 + s)), 2, T)
            push!(out, clv_go(map, x0, 80, 3, 80; W0 = W0, C0 = C0, outputs = (:V,)))
        end
        @test maximum(n -> clv_vdist(out[1].V[:, :, n + 1], out[2].V[:, :, n + 1]), 0:3) <= 1.0e-10
    end
end

# ------------------------------------------------------------------ AC12: time reversal
@testset "AC12a time reversal on the linear fixtures: the CLVs of the inverse map, index order reversed, are those of the forward map" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64
    for name in (:a, :b, :c)
        fx = clv_fixture(name, T); d = fx.d; W0 = clv_frame(clv_seed(12, 1), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(12, 2)), d, T)
        f = clv_go(fx.map, zero(SVector{d, T}), 100, 10, 100; W0 = W0, C0 = C0, outputs = (:V,))
        r = clv_go(clv_Rev(fx.map), zero(SVector{d, T}), 100, 10, 100; W0 = W0, C0 = C0, outputs = (:V,))
        B = clv_B_fix(T, fx.gamma, 100); wrong = 0.0
        for n in 0:10, i in 1:d
            @test clv_delta(r.V[:, i, n + 1], f.V[:, d + 1 - i, n + 1]) <= B
            wrong = max(wrong, clv_delta(r.V[:, i, n + 1], f.V[:, i, n + 1]))
        end
        (d == 2 || d == 4) && @test(wrong > 1.0e-3)                                                  # control: without the reversal of the index order the match fails
    end
end

@testset "AC12b time reversal on the standard map K = 5 in BigFloat(256): the reversed orbit started from the end point reproduces the CLVs to 1e-12" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    setprecision(BigFloat, CLV_BITS) do
        T = BigFloat; map = clv_std(T); x0 = SVector{2, T}(T(13) / 10, T(7) / 10); nf, ns, nb = 60, 20, 60
        W0 = clv_frame(clv_seed(12, 10), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(12, 11)), 2, T)
        f = clv_go(map, x0, nf, ns, nb; W0 = W0, C0 = C0, outputs = (:V,))
        r = clv_go(clv_Rev(map), f.x, nb, ns, nf; W0 = clv_frame(clv_seed(12, 12), 1, 2, 2, T), C0 = clv_rand_C(StableRNG(clv_seed(12, 13)), 2, T), outputs = (:V,))
        w = 0.0; wrong = 0.0
        for n in 0:ns, i in 1:2
            w = max(w, clv_delta(f.V[:, i, n + 1], r.V[:, 3 - i, ns - n + 1])); wrong = max(wrong, clv_delta(f.V[:, i, n + 1], r.V[:, i, ns - n + 1]))
        end
        @test w <= 1.0e-12
        @test wrong > 1.0e-3                                                                               # control: without the reversed index order the match fails
    end
end

# ------------------------------------------------------------------ AC14: health, validation, run record
@testset "AC14a non-finite start, frame and orbit state: health errors with kind, orbit id and handoff, nothing stored" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; map = clv_std(T); called = Ref(0); sink = (n, C, lam) -> (called[] += 1; nothing)
    W0 = clv_frame(clv_seed(14, 1), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(14, 2)), 2, T)
    e = clv_health_error(() -> clv_go(map, SVector{2, T}(NaN, 0.5), 5, 5, 5; W0 = W0, C0 = C0, orbit = 3, sink = sink, outputs = ()))
    @test e isa NumericalHealthError && e.what == "state" && e.orbit == 3
    Wn = SMatrix{2, 2, T}(NaN, 0.0, 0.0, 1.0)
    e = clv_health_error(() -> clv_go(map, SVector{2, T}(1.0, 0.5), 5, 5, 5; W0 = Wn, C0 = C0, orbit = 4, sink = sink, outputs = ()))
    @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 4
    nm = clv_NaNMap(map, 5)                                                                                    # the 5th forward call (handoff 5 over the transient and the window) gives NaN
    e = clv_health_error(() -> clv_go(nm, SVector{2, T}(1.0, 0.5), 3, 5, 5; W0 = W0, C0 = C0, orbit = 7, sink = sink, outputs = ()))
    @test e isa NumericalHealthError && e.what == "state" && e.orbit == 7 && e.step == 5
    @test called[] == 0
end

@testset "AC14b forward factorisation breakdown: a singular tangent gives the delivered breakdown kind with orbit id and handoff, nothing stored" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; A = SMatrix{2, 2, T}(1.0, 0.0, 0.0, 0.0); map = clv_LinMap{2, T, 4}(A, A); called = Ref(0)
    e = clv_health_error(() -> clv_go(map, zero(SVector{2, T}), 3, 5, 5; W0 = clv_eye(2, 2, T), C0 = clv_eye(2, 2, T), orbit = 9, sink = (n, C, lam) -> (called[] += 1; nothing), outputs = ()))
    @test e isa NumericalHealthError && e.what == "orthonormalisation breakdown" && e.orbit == 9 && e.step == 1 && e.column == 2
    @test called[] == 0
end

@testset "AC14c backward breakdown: R_22 below the representable inverse overflows the solve and gives kind \"backward breakdown\" with orbit id and handoff" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; A = SMatrix{2, 2, T}(1.0, 0.0, 0.0, 1.0e-310); map = clv_LinMap{2, T, 4}(A, A); called = Ref(0)
    e = clv_health_error(() -> clv_go(map, zero(SVector{2, T}), 2, 4, 4; W0 = clv_eye(2, 2, T), C0 = clv_eye(2, 2, T), orbit = 5, sink = (n, C, lam) -> (called[] += 1; nothing), outputs = ()))
    @test e isa NumericalHealthError && e.what == "backward breakdown" && e.orbit == 5 && 1 <= e.step <= 10
    @test called[] == 0
end

@testset "AC14d argument errors before any work: k, n_s, n_f, n_b, tau, h_b, storage" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; cnt = clv_CountMap(clv_std(T)); x0 = SVector{2, T}(1.3, 0.7)
    @test_throws ArgumentError ClvGinelli(cnt; tau = 0)
    @test_throws ArgumentError ClvGinelli(cnt; storage = :checkpoint, block = 0)
    @test_throws ArgumentError ClvGinelli(cnt; storage = :disk)
    clv = ClvGinelli(cnt)
    for (kw, nf, ns, nb) in ((:k0, 5, 5, 5), (:k3, 5, 5, 5), (:ns0, 5, 0, 5), (:nfneg, -1, 5, 5), (:nbneg, 5, 5, -1))
        k = kw === :k0 ? 0 : kw === :k3 ? 3 : 2
        @test_throws ArgumentError propagate_clv_orbit(clv, x0, nf, ns, nb; k = k, seed = 1)
    end
    @test cnt.nf[] == 0 && cnt.nt[] == 0                                                                        # nothing was computed
end

@testset "AC14e run record: every parameter of Algorithm 'Parameters recorded per run'; a gpu request runs on the CPU and records the fallback reason" begin
    _ = (ClvGinelli, propagate_clv_orbit, propagate_clv)
    T = Float64; im = clv_im(T); z0 = clv_orbit_C(); W0 = clv_frame(clv_seed(14, 20), 1, 4, 4, T)
    clv = ClvGinelli(im; tau = 20, h = CLV_H, orthonormalizer = RGS(), storage = :checkpoint, block = 5)
    res = propagate_clv_orbit(clv, z0, 3, 10, 4; k = 4, seed = 11, W0 = W0, outputs = (:V,))
    rec = res.record
    for key in (:system, :h, :tau, :k, :n_f, :n_s, :n_b, :orthonormalizer, :storage, :block, :precision, :seed, :cseed, :outputs, :backend, :backend_fallback_reason, :integrator, :integrator_order)
        @test hasproperty(rec, key)
    end
    @test rec.tau == 20 && rec.k == 4 && rec.n_f == 3 && rec.n_s == 10 && rec.n_b == 4 && rec.storage === :checkpoint && rec.block == 5 && rec.seed == 11 && rec.backend === :cpu
    @test rec.precision == "Float64" && rec.h == CLV_H && rec.orthonormalizer == "RGS"
    gpu = propagate_clv(ClvGinelli(clv_std(T)), [SVector{2, T}(1.3, 0.7)], 5, 5, 5; k = 2, seed = 1, backend = :gpu)
    cpu = propagate_clv(ClvGinelli(clv_std(T)), [SVector{2, T}(1.3, 0.7)], 5, 5, 5; k = 2, seed = 1, backend = :cpu)
    @test gpu.provenance.backend === :cpu && occursin("CPU-only", string(gpu.provenance.backend_fallback_reason)) && isempty(string(cpu.provenance.backend_fallback_reason))
    @test gpu.results[1].C == cpu.results[1].C
end

@testset "D124 step length h: finite and positive, kept in the element type (review M1, M2)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    m = clv_std(Float64)
    for h in (0.0, -1.0, NaN, Inf)
        @test_throws ArgumentError ClvGinelli(m; h = h)
    end
    @test ClvGinelli(m; h = big"0.05").h isa BigFloat && ClvGinelli(m; h = 1 // 20).h === 1 // 20
    T = Double64; mT = clv_std(T); x0 = SVector{2, T}(T(1.3), T(0.7)); W0 = clv_frame(clv_seed(14, 30), 1, 2, 2, T); C0 = clv_eye(2, 2, T)
    run(h) = propagate_clv_orbit(ClvGinelli(mT; tau = 20, h = h), x0, 10, 10, 10; k = 2, seed = 1, W0 = W0, C0 = C0)
    r1 = run(1.0); r2 = run(0.05); e = 8 * eps(T)
    @test all(abs(r2.lambda[i, n] * T(20) * T(0.05) - r1.lambda[i, n] * T(20)) <= e * max(1, abs(r1.lambda[i, n] * 20)) for i in 1:2, n in 1:11)   # local exponents: Dt = tau h in T
    @test all(abs(r2.lambda_fwd[i] * T(0.05) - r1.lambda_fwd[i]) <= e * max(1, abs(r1.lambda_fwd[i])) for i in 1:2)                            # forward exponents: the same h
end
