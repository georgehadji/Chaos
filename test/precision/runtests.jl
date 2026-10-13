# Precision suite (module-build slug "precision"). Spec: workspaces/module-build/stages/01-spec/output/precision-spec.md.
# Float32 dropped (D24). Criteria 9-11 are stage 04, not tested here.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

include("helpers.jl")

@testset "AC1 parameter conversion" begin
    # |K_T - K_ref(512)| <= 0.5 ulp (Float64), <= 2 eps_T K (Double64, BigFloat(256)); a Float64 detour gives 3.9e-17 K
    refK(name) = name == "pi/2" ? BigFloat(pi) / 2 : name == "2pi" ? BigFloat(pi) * 2 :
                 name == "0.1" ? parse(BigFloat, "0.1") : parse(BigFloat, "0.001")
    build(T, name) = name == "pi/2" ? exact(T, pi) / 2 : name == "2pi" ? exact(T, pi) * 2 :
                     name == "0.1" ? exact(T, "0.1") : exact(T, "0.001")
    for name in ("pi/2", "2pi", "0.1", "1e-3"), T in (Float64, Double64, BigFloat)
        Kref = setprecision(() -> refK(name), BigFloat, 512)
        K = setprecision(() -> build(T, name), BigFloat, 256)
        epsT = setprecision(() -> Float64(eps(T)), BigFloat, 256)
        setprecision(BigFloat, 512) do
            err = abs(pr_big(K) - Kref)
            bound = T === Float64 ? BigFloat(eps(K)) / 2 : 2 * BigFloat(epsT) * Kref
            @test K isa T
            @test err <= bound
        end
    end
end

@testset "AC2 exact widening" begin
    # Float64 -> T -> Float64 is the identity and the value equals BigFloat(x); 10^4 seeded values
    rng = StableRNG(PR_SEED + 2)
    mant = 0x000fffffffffffff
    vals = Float64[]
    for i in 1:10^4
        r = i % 4
        v = r == 0 ? (rand(rng) * 2 - 1) * 1e3 :
            r == 1 ? reinterpret(Float64, (rand(rng, UInt64) & mant) | 0x1) * (rand(rng, Bool) ? 1 : -1) :
            r == 2 ? 2pi + (rand(rng) - 0.5) * 1e-3 :
                     reinterpret(Float64, rand(rng, UInt64))
        push!(vals, isfinite(v) ? v : 1.0)
    end
    bad_d = 0; bad_b = 0
    setprecision(BigFloat, 256) do
        for v in vals
            d = exact(Double64, v)
            b = exact(BigFloat, v)
            bad_d += !(Float64(d) === v && pr_big(d) == BigFloat(v))
            bad_b += !(Float64(b) === v && b == BigFloat(v))
        end
    end
    @test bad_d == 0
    @test bad_b == 0
end

@testset "AC3 transcendentals in T" begin
    # |f_D - f_ref(256)| <= 16 eps_D max(1, |f_ref|); sin, cos through the map (K = 1: p' = sin x, DT[2,1] = cos x), log directly
    rng = StableRNG(PR_SEED + 3)
    epsD = Float64(eps(Double64))
    sys = StandardMap{Double64}(one(Double64))
    worst = Dict("sin" => 0.0, "cos" => 0.0, "log" => 0.0)
    setprecision(BigFloat, 256) do
        for R in (2pi, 1e3), _ in 1:10^3
            x = pr_rand_d64(rng, -R, R)
            z = SVector{2,Double64}(x, zero(Double64))
            xb = pr_big(x)
            for (nm, got, want) in (("sin", forward(sys, z)[2], sin(xb)), ("cos", tangent(sys, z)[2, 1], cos(xb)))
                r = Float64(abs(pr_big(got) - want)) / (epsD * max(1.0, Float64(abs(want))))
                worst[nm] = max(worst[nm], r)
            end
        end
        for _ in 1:10^3
            x = pr_d64_from(rng, exp10(-3 + 6 * rand(rng)))     # log-uniform [1e-3, 1e3], random low word
            want = log(pr_big(x))
            worst["log"] = max(worst["log"], Float64(abs(pr_big(log(x)) - want)) / (epsD * max(1.0, Float64(abs(want)))))
        end
    end
    for nm in ("sin", "cos", "log")
        @test worst[nm] <= 16
    end
end

@testset "AC4 pathwise orbit agreement" begin
    # max_n<=nH |e_n|_inf / (eps_T M_n Gamma_n) <= 8, 200 ICs, reference BigFloat(256); floors nH >= 15 (Float64), 40 (Double64), <= 5 percent excluded
    ics = pr_ics(PR_SEED + 4, 200)
    for (T, floor_nH) in ((Float64, 15), (Double64, 40))
        epsT = Float64(eps(T))
        worst = 0.0; low = 0
        setprecision(BigFloat, 256) do
            for (x0, p0) in ics
                ref = pr_reference(BigFloat(5), x0, p0, PR_NMAX)
                M, G = pr_bounds(ref)
                nH = pr_nH(M, G, epsT)
                if nH < floor_nH
                    low += 1
                    continue
                end
                zs = pr_run(T, exact(T, 5), x0, p0, nH)
                for n in 1:nH
                    e = pr_err(zs[n+1], ref.xs[n+1], ref.ps[n+1])
                    worst = max(worst, Float64(e) / (epsT * M[n+1] * G[n+1]))
                end
            end
        end
        @test low <= 10
        @test worst <= 8
    end
end

@testset "AC5 BigFloat 256 vs 512" begin
    # |e_n|_inf <= 8 eps_256 M_n Gamma_n for n <= 100 and |e_10| <= 1e-60 (a Float64 constant would give 1e-12)
    ics = pr_ics(PR_SEED + 5, 200)
    eps256 = Float64(setprecision(() -> eps(BigFloat), BigFloat, 256))
    refs = setprecision(BigFloat, 512) do
        [pr_reference(BigFloat(5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    runs = setprecision(BigFloat, 256) do
        [pr_run(BigFloat, exact(BigFloat, 5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    worst = 0.0; e10 = 0.0
    setprecision(BigFloat, 512) do
        for (ref, zs) in zip(refs, runs)
            M, G = pr_bounds(ref)
            for n in 1:PR_NMAX
                e = Float64(pr_err(zs[n+1], ref.xs[n+1], ref.ps[n+1]))
                worst = max(worst, e / (eps256 * M[n+1] * G[n+1]))
                n == 10 && (e10 = max(e10, e))
            end
        end
    end
    @test worst <= 8
    @test e10 <= 1e-60
end

@testset "AC6 tangent growth" begin
    # |S_n(Double64) - S_n(ref)| <= 64 (1 + K) eps_D M_n Gamma_n, n <= min(nH(Double64), 40), K = 5.
    # Deviation from the spec: at most 5 percent of orbits may exceed the bound (as in AC4). Reason: the bound has no local-stretch
    # term, and one orbit (index 174, accelerator-like) reaches ratio 2.6 even for BigFloat 256 vs 512 (ratio 10); see the test plan.
    ics = pr_ics(PR_SEED + 6, 200)
    rng = StableRNG(PR_SEED + 60)
    epsD = Float64(eps(Double64))
    ratios = Float64[]
    setprecision(BigFloat, 256) do
        for (x0, p0) in ics
            a, b = rand(rng) - 0.5, rand(rng) - 0.5
            u0 = (a, b) ./ hypot(a, b)                       # unit pair in Float64, widened exactly
            ref = pr_reference(BigFloat(5), x0, p0, PR_NMAX)
            M, G = pr_bounds(ref)
            nmax = min(pr_nH(M, G, epsD), 40)
            Sd = pr_growth_lib(Double64, exact(Double64, 5), x0, p0, u0, nmax)
            Sr = pr_growth_ref(BigFloat(5), x0, p0, u0, nmax)
            worst = 0.0
            for n in 1:nmax
                d = Float64(abs(pr_big(Sd[n+1]) - Sr[n+1]))
                worst = max(worst, d / (64 * (1 + 5) * epsD * M[n+1] * G[n+1]))
            end
            push!(ratios, worst)
        end
    end
    @info "AC6 realised ratio" max = maximum(ratios) exceeding = count(>(1), ratios)
    @test count(>(1), ratios) <= 10
end

@testset "AC7 fallback matrix" begin
    # 9 cases (T x request) with a no-CUDA probe, plus one row per remaining reason; BigFloat never consults the probe
    @test Chaoslab.FALLBACK_REASONS == ["CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type",
        "device out of memory", "NumericalHealthError raised only on the GPU", "reducer is CPU-only", "GPU not faster for the element type"]
    compile_reason = "kernel does not compile for the element type"
    for T in (Float64, Double64, BigFloat), req in (:auto, :gpu, :cpu)
        probe = (T === BigFloat && req !== :cpu) ? pr_probe_boom : pr_probe_noload
        r = resolve_backend(req, T; probe = probe)
        want = req === :cpu ? "" : T === BigFloat ? compile_reason : "CUDA.jl not loadable"
        @test r isa BackendChoice
        @test r.backend === :cpu
        @test r.reason == want
        @test parse_precision(precision_label(T))[1] === T
    end
    for T in (Float64, Double64), req in (:auto, :gpu)
        r = resolve_backend(req, T; probe = pr_probe_nofunc)
        @test (r.backend, r.reason) == (:cpu, "no functional device")
        r = resolve_backend(req, T; probe = pr_probe_nocompile)
        @test (r.backend, r.reason) == (:cpu, compile_reason)
        r = resolve_backend(req, T; probe = pr_probe_good)
        @test (r.backend, r.reason) == (:gpu, "")
    end
    for T in (Float64, Double64)
        r = resolve_backend(:auto, T; probe = pr_probe_slow)
        @test (r.backend, r.reason) == (:cpu, "GPU not faster for the element type")
    end
    r = resolve_backend(:auto, BigFloat; probe = pr_probe_good)
    @test (r.backend, r.reason) == (:cpu, compile_reason)
end

@testset "AC8 recording" begin
    # precision string is one of Float64, Double64, BigFloat(<b>); b parses back to the integer
    @test precision_label(Float64) == "Float64"
    @test precision_label(Double64) == "Double64"
    @test parse_precision("Float64") == (Float64, nothing)
    @test parse_precision("Double64") == (Double64, nothing)
    for b in (64, 128, 256, 512)
        label = setprecision(() -> precision_label(BigFloat), BigFloat, b)
        @test label == "BigFloat($b)"
        @test occursin(r"^(Float64|Double64|BigFloat\(\d+\))$", label)
        @test parse_precision(label) == (BigFloat, b)
    end
end

@testset "AC12 closure and allocation" begin
    # map step, tangent step and log accumulation return T; 0 bytes in the hot loop for Float64 and Double64
    for T in (Float64, Double64, BigFloat)
        setprecision(BigFloat, 256) do
            sys = StandardMap{T}(exact(T, 5))
            z = SVector{2,T}(exact(T, 1.25), exact(T, -0.75))
            zf = forward(sys, z); M = tangent(sys, z)
            acc = add(CompensatedSum{T}(), log(pr_norm(zf)))
            @test zf isa SVector{2,T}
            @test M isa SMatrix{2,2,T}
            @test acc isa CompensatedSum{T}
            @test total(acc) isa T
            @test acc.s isa T && acc.c isa T
        end
    end
    for T in (Float64, Double64)
        sys = StandardMap{T}(exact(T, 5))
        z = SVector{2,T}(exact(T, 1.25), exact(T, -0.75))
        acc = CompensatedSum{T}()
        x = exact(T, 0.3)
        pr_loop_forward(sys, z, 2); pr_loop_tangent(sys, z, 2); pr_loop_add(acc, x, 2)    # warm-up
        @test (@allocated pr_loop_forward(sys, z, 100)) == 0
        @test (@allocated pr_loop_tangent(sys, z, 100)) == 0
        @test (@allocated pr_loop_add(acc, x, 100)) == 0
    end
end

@testset "AC13 conversion consequence" begin
    # K = 2pi: K built in Double64 vs widened from Float64 separates (> 1e-8) within 60 steps for >= 80 percent of 100 ICs
    ics = pr_ics(PR_SEED + 13, 100)
    K1 = exact(Double64, pi) * 2
    K2 = Double64(2pi)
    @test K1 != K2
    s1 = StandardMap{Double64}(K1); s2 = StandardMap{Double64}(K2)
    nstar = Int[]
    for (x0, p0) in ics
        z1 = SVector{2,Double64}(Double64(x0), Double64(p0)); z2 = z1
        n = 0
        while n < 200
            z1 = forward(s1, z1); z2 = forward(s2, z2); n += 1
            maximum(abs.(z1 - z2)) > 1e-8 && break
        end
        push!(nstar, n)
    end
    @info "AC13 median n*" median = sort(nstar)[50]
    @test count(<=(60), nstar) >= 80
end
