# AC10 to AC13: chains and transparency, execution invariance, the independent oracles (Poisson form, BigFloat orbit), CPU against GPU.
# The reducer is a reader (it never replaces columns and never changes the orbit). Chains whose other members are readers end with a test-side replacer (wb_Norm) so that a long run of an expanding
# map does not overflow: the spec's chains {empty} and {this reader alone} become {Trace, Norm} and {Trace, reader, Norm}. The chains with E2 have 2N + 1 handoffs (Items for the owner).

@testset "AC10 chains and transparency" begin
    _ = (WeightedBirkhoff, FLI, SALI, MEGNO, check_run)
    # K_std = 2 pi, N = 1000, both types: orbit points bitwise equal at every handoff for the chains {Norm}, {reader, Norm}, {FLI(replace = false), reader, SALI}, {reader, E2, MEGNO} and {E2, reader,
    # MEGNO} (both orders of the reader and E2); the deviation matrix of a chain bitwise that of the chain without the reader; the other reducers' outputs bitwise their stand-alone values and the
    # reader's output bitwise its stand-alone value in every chain; the reader's averages after 2N + 1 handoffs bitwise those after 2N.
    N = 1000
    nsteps = 2N + 1
    for T in WB_TYPES
        c = wb_case(:two_pi, T, 1, 1)
        sys, x0 = c.sys, c.x0[1]
        W0 = wb_frame(wb_seed(10, 1), 1, T, 2, 2)
        r = WeightedBirkhoff(; window = N, observables = wb_std5())
        e2 = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        xA, WA, fA = wb_run((wb_Trace(), wb_Norm()), sys, x0, W0, nsteps, 1)
        xB, WB, fB = wb_run((wb_Trace(), r, wb_Norm()), sys, x0, W0, nsteps, 1)
        xC, WC, fC = wb_run((wb_Trace(), FLI(; replace = false), r, SALI()), sys, x0, W0, nsteps, 1)
        xD, WD, fD = wb_run((wb_Trace(), r, e2, MEGNO()), sys, x0, W0, nsteps, 1)
        xE, WE, fE = wb_run((wb_Trace(), e2, r, MEGNO()), sys, x0, W0, nsteps, 1)
        _, WS1, fS1 = wb_run((FLI(),), sys, x0, W0, nsteps, 1)
        _, WS2, fS2 = wb_run((SALI(),), sys, x0, W0, nsteps, 1)
        _, WS3, fS3 = wb_run((MEGNO(),), sys, x0, W0, nsteps, 1)
        _, _, fS4 = wb_run((e2, wb_Norm()), sys, x0, W0, nsteps, 1)
        _, _, fS5 = wb_run((r, wb_Norm()), sys, x0, W0, 2N, 1)
        @test length(fA[1].xs) == nsteps
        @test all(f -> isequal(f[1].xs, fA[1].xs), (fB, fC, fD, fE)) && xA == xB == xC == xD == xE
        @test isequal(WA, WB) && isequal(WC, WS2) && isequal(WD, WS3) && isequal(WE, WS3)                 # the reader returns the deviation matrix unchanged
        @test isequal(fC[2], fS1[1]) && isequal(fC[4], fS2[1]) && isequal(fD[4], fS3[1]) && isequal(fE[4], fS3[1]) && isequal(fD[3], fS4[1]) && isequal(fE[2], fS4[1])   # the other reducers
        stand = fB[2]
        @test stand.complete && stand.handoffs == nsteps && isequal(fC[3], stand) && isequal(fD[2], stand) && isequal(fE[3], stand)                   # the reader in every chain
        flds = (:avg1, :avg2, :S, :B, :delta, :absdig, :reldig, :dig, :below_floor, :zero)
        @test fS5[1].complete && fS5[1].handoffs == 2N && all(fld -> isequal(getfield(fS5[1], fld), getfield(stand, fld)), flds)
    end
end

@testset "AC11 execution invariance" begin
    _ = (WeightedBirkhoff, MEGNO)
    # The reducer through propagate (chunk sizes 64 and 1000, several threads) equals the serial per-orbit body bitwise; K = 4 observables, H = 2 horizons, K_std in {0.5, 2 pi} with uniform
    # starts. D38: 128 orbits, N = 250 (spec: 4096 orbits, N = 1000, 1 versus N threads in child processes at stage 04; the thread count of this process is whatever Julia was started with).
    norb = WB_FULL ? 4096 : 128
    N = WB_FULL ? 1000 : 250
    for (label, sub) in ((:half, 5), (:two_pi, 1))
        sys = WrappedStandardMap(wb_K(Float64, label))
        starts = wb_starts(StableRNG(wb_seed(11, sub)), norb)
        x0 = [torus_state(wb_cast(Float64, s)) for s in starts]
        r = WeightedBirkhoff(; window = (N ÷ 2, N), observables = wb_std4())
        reducers = (r, MEGNO())
        seed = wb_seed(11, 10 + sub)
        W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:norb]
        ser = [propagate_orbit(sys, x0[j], W0[j], 2N, 1, reducers; orbit = j) for j in 1:norb]
        fs = [map(reducer_final, reducers, s[3]) for s in ser]
        thr = Vector{Any}(undef, norb)
        Threads.@threads :dynamic for j in 1:norb
            s = propagate_orbit(sys, x0[j], W0[j], 2N, 1, reducers; orbit = j)
            thr[j] = (s[1], map(reducer_final, reducers, s[3]))
        end
        @test isequal([t[1] for t in thr], [s[1] for s in ser]) && isequal([t[2] for t in thr], fs)
        for chunk in (64, 1000)
            res = propagate(sys, x0, 2N; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            @test isequal(res.x, [s[1] for s in ser]) && isequal(res.finals, fs)
        end
    end
end

@testset "AC12a Poisson form" begin
    _ = (weighted_birkhoff,)
    # WB17: f = cos 2 pi (theta_0 + j alpha), theta_0 = 1/7, alpha = sqrt 2 - 1, N in {50, 100, 200}: the BigFloat(256) average from the module's function equals the Poisson form (trapezoid quadrature of
    # W(s), a different algorithm from every sum of the module, wb_poisson_avg) within 1e-18 absolute; and the independent plain-array sum of the reference within 1e-60 (the function accepts BigFloat).
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    th0 = setprecision(BigFloat, WB_BITS) do
        BigFloat(1) / 7
    end
    for N in (50, 100, 200)
        ser, _ = wb_rot_big(th0, alpha, 1, 2N)
        r = setprecision(BigFloat, WB_BITS) do
            weighted_birkhoff(ser, N)
        end
        pois = wb_poisson_avg(N, th0, alpha)
        ref = wb_ref_avgs(ser, N, 1)
        @test r.avg1 isa Matrix{BigFloat} && abs(pois) > 1e-30
        @test abs(r.avg1[1] - pois) <= 1e-18
        @test abs(r.avg1[1] - ref.avg1) <= 1e-60 && abs(r.avg2[1] - ref.avg2) <= 1e-60 && abs(r.S[1] - ref.S) <= 1e-70
        @info "AC12a N=$N: |function - Poisson|, |function - plain-array sum| (spec dry run: 1.9e-21, 2.1e-30, 3.8e-43 against the direct sum)" poisson = Float64(abs(r.avg1[1] - pois)) plain = Float64(abs(r.avg1[1] - ref.avg1))
    end
end

@testset "AC12b independent orbit" begin
    _ = (WeightedBirkhoff, CosObservable, SinObservable)
    # K_std = 0.5 regular class, Float64, four trigonometric observables, N in {1000, 4000}: |avg(module orbit) - avg(BigFloat orbit)| <= 1e-12 for both windows (order-of-magnitude bound, derived: a rounding
    # perturbation of a regular orbit is a nearby regular orbit; proposal). The BigFloat orbit is a plain lifted iteration of the standard map from the exactly widened start (wb_orbit_big).
    # D38: 4 starts, N = 1000 (spec: 20 starts, both N).
    nstart = WB_FULL ? 20 : 4
    Ns = WB_FULL ? (1000, 4000) : (1000,)
    c = wb_case(:half, Float64, nstart, 0)
    obs = (CosObservable(1), SinObservable(1), CosObservable(2), SinObservable(2))
    bad = String[]
    worst = 0.0
    for N in Ns, j in 1:nstart
        red = WeightedBirkhoff(; window = N, observables = obs)
        _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(12, 1), c.ids[j], Float64), 2N, 1)
        f = fin[1]
        orb = wb_orbit_big(0.5, c.starts[j], 2N)
        cols = setprecision(BigFloat, WB_BITS) do
            [[cos(o[1]) for o in orb], [sin(o[1]) for o in orb], [cos(o[2]) for o in orb], [sin(o[2]) for o in orb]]
        end
        for k in 1:4
            ref = wb_ref_avgs(cols[k], N, 1)
            for (a, rr) in ((f.avg1[k], ref.avg1), (f.avg2[k], ref.avg2))
                e = wb_err(a, rr)
                worst = max(worst, e)
                e <= 1e-12 || push!(bad, "N=$N start $j obs $k: $e")
            end
        end
    end
    @test wb_report(bad, "AC12b")
    @info "AC12b: worst |avg(module orbit) - avg(BigFloat orbit)| (spec dry run 1.6e-14, 3.4e-14)" worst
    # AC12(c): no external weighted-Birkhoff routine exists in the environment (spec, Oracle): there is nothing to test.
end

# A functional device exists in this process when the Foundation backend resolution (the CUDA extension is loaded by the caller: verify/weighted-birkhoff.jl) says so.
wb_device_ok() = try
    Chaoslab.resolve_backend(:gpu, Float64).backend === :gpu
catch
    false
end

@testset "AC13 CPU against GPU" begin
    _ = (WeightedBirkhoff, weighted_birkhoff, CosObservable, SinObservable, LiftedObservable)
    # Tag (G): needs a functional device, so it runs only when the caller has loaded the CUDA extension (verify/weighted-birkhoff.jl at stage 04); otherwise it is skipped here with a note.
    # (a) K_std = 0.5 regular class and 10, 1024 orbits each, Float64, N = 1000, K = 4 trigonometric observables, H = 1, p = 1, chain (reader, Benettin()): the kernel compiles and runs (provenance.backend
    # == :gpu), every orbit is complete and finite; (b) on IDENTICAL samples: the K_std = 0 orbit has bitwise equal samples on both backends for observables built from additions, products and
    # the wrap only (the coordinate, the lifted q), so the device reducer against the CPU pure function on the host samples differs only by the device exponential: |avg difference| <= 4 eps S and
    # |Delta difference| <= 8 eps S (the sample equality of the K = 0 orbit is itself asserted by the final states); (c) the orbit-level CPU-GPU gap on the chaotic class is reported, not gated;
    # (d) Double64: the kernel compiles, or the recorded fallback is the pass (reported).
    if !wb_device_ok()
        @info "AC13 skipped: no functional device in this process (stage 04 runs it through verify/weighted-birkhoff.jl)"
        @test_skip false
    else
        N = 1000
        m = WB_FULL ? 1024 : 256
        epsT = wb_eps(Float64)
        red(obs) = WeightedBirkhoff(; window = N, observables = obs)
        for (label, sub) in ((:half, 0), (:ten, 3))
            sys = WrappedStandardMap(wb_K(Float64, label))
            x0 = label === :half ? wb_case(:half, Float64, min(m, 50), 0).x0 : [torus_state(wb_cast(Float64, s)) for s in wb_starts(StableRNG(wb_seed(13, sub)), m)]
            gpu = propagate(sys, x0, 2N; k = 1, tau = 1, seed = wb_seed(13, 9), reducers = (red(wb_std4()), Benettin()), backend = :gpu)
            cpu = propagate(sys, x0, 2N; k = 1, tau = 1, seed = wb_seed(13, 9), reducers = (red(wb_std4()), Benettin()), backend = :cpu)
            @test gpu.provenance.backend === :gpu && all(f -> f[1].complete && all(isfinite, f[1].avg1) && all(isfinite, f[1].avg2), gpu.finals)
            gap = [maximum(abs.(Matrix(g[1].avg1) .- Matrix(c[1].avg1))) for (g, c) in zip(gpu.finals, cpu.finals)]
            @info "AC13(c) $label: orbit-level CPU-GPU gap of the averages (reported, not gated)" median = wb_quant(gap, 0.5) maximum = maximum(gap)
        end
        sys0 = WrappedStandardMap(wb_K(Float64, :zero))
        x0 = [torus_state(wb_cast(Float64, s)) for s in wb_starts(StableRNG(wb_seed(13, 20)), m)]
        obs0 = (wb_Coord(1), LiftedObservable(1, 3))
        gpu = propagate(sys0, x0, 2N; k = 1, tau = 1, seed = wb_seed(13, 9), reducers = (red(obs0), Benettin()), backend = :gpu)
        bad = String[]
        for j in eachindex(x0)
            sts = wb_states(sys0, x0[j], 2N)
            sts[end] == gpu.x[j] || push!(bad, "orbit $j: final state differs between backends")
            g = weighted_birkhoff(wb_series(obs0, sts), N)
            f = gpu.finals[j][1]
            for k in 1:2
                S = Float64(g.S[k])
                abs(f.avg1[k] - g.avg1[k]) <= 4 * epsT * S && abs(f.avg2[k] - g.avg2[k]) <= 4 * epsT * S || push!(bad, "orbit $j obs $k: average")
                abs(f.delta[k] - (g.avg1[k] - g.avg2[k])) <= 8 * epsT * S || push!(bad, "orbit $j obs $k: Delta")
            end
        end
        @test wb_report(bad, "AC13(b)")
        sysd = WrappedStandardMap(wb_K(Double64, :half))
        xd = wb_case(:half, Double64, 16, 0).x0
        r64 = try
            propagate(sysd, xd, 2N; k = 1, tau = 1, seed = wb_seed(13, 9), reducers = (red(wb_std4()), Benettin()), backend = :gpu)
        catch e
            e
        end
        @test !(r64 isa Exception) && all(f -> f[1].complete && all(isfinite, f[1].avg1), r64.finals)
        r64 isa Exception || @info "AC13(d) Double64: provenance (the pass is a compiled kernel or a recorded fallback reason)" r64.provenance
    end
end
