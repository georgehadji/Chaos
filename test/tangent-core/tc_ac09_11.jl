# AC9 to AC11: compensated summation, backend resolution and provenance, chunk size.

@testset "AC9 compensated summation" begin
    # Exact = 1e8 fl(0.1) (rational); compensated error <= eps; plain error ~1.9e-9 (spec, measured) must exceed it.
    x = 0.1
    n = 10^8
    exact = n * Rational{BigInt}(x)
    relerr(v) = Float64(abs(Rational{BigInt}(v) - exact) / exact)
    @test relerr(tc_csum(n, x)) <= eps(Float64)
    # (c) Double64, 1e7 copies of fl(0.1)
    n7 = 10^7
    exact7 = n7 * Rational{BigInt}(x)
    c7 = tc_csum(n7, Double64(x))
    @test Float64(abs((Rational{BigInt}(c7.hi) + Rational{BigInt}(c7.lo)) - exact7) / exact7) <= Float64(eps(Double64))
    # Neumaier branches on the larger magnitude: [1, 1e100, 1, -1e100] sums to 2 (plain and Kahan give 0)
    acc = CompensatedSum{Float64}()
    for v in (1.0, 1e100, 1.0, -1e100)
        acc = add(acc, v)
    end
    @test total(acc) == 2.0
    # (b) negative control, last so that it only runs once the compensated pieces exist: the plain sum fails (a)
    @test relerr(tc_plainsum(n, x)) > 100 * relerr(tc_csum(n, x)) + 1e-10
end

@testset "AC10 backend resolution and provenance" begin
    # Check order: element type, cpu_only_reducer, loadable, functional, compiles, gpu_faster. Reasons are exact strings.
    seed = TC_SEED + 10
    good = () -> (loadable = true, functional = true, device = "TestGPU", compiles = true, gpu_faster = true)
    withp(; kw...) = () -> merge(good(), NamedTuple(kw))
    nocall = () -> error("probe must not be called")
    @test FALLBACK_REASONS == [
        "CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type",
        "device out of memory", "NumericalHealthError raised only on the GPU", "reducer is CPU-only",
        "GPU not faster for the element type",
    ]
    # (a) CI without CUDA: default probe
    c = resolve_backend(:auto, Float64)
    @test c isa BackendChoice && c.backend === :cpu && c.device == "cpu"
    @test c.reason in ("CUDA.jl not loadable", "no functional device")
    # (b) cpu request never falls back
    c = resolve_backend(:cpu, Float64; probe = good)
    @test c.backend === :cpu && c.device == "cpu" && c.reason == ""
    # a working device is used for auto and gpu, no reason
    for req in (:auto, :gpu), T in (Float64, Double64)
        c = resolve_backend(req, T; probe = good)
        @test c.backend === :gpu && c.device == "TestGPU" && c.reason == ""
    end
    # (c) gpu request without a device falls back
    c = resolve_backend(:gpu, Float64; probe = withp(functional = false))
    @test c.backend === :cpu && c.device == "cpu" && c.reason == "no functional device"
    # (d) BigFloat never reaches the probe
    for req in (:auto, :gpu)
        c = resolve_backend(req, BigFloat; probe = nocall)
        @test c.backend === :cpu && c.reason == "kernel does not compile for the element type"
    end
    # (e) each probe-forced reason, exact strings (reasons 4 and 5 are per-chunk: stage 04)
    @test resolve_backend(:auto, Float64; probe = withp(loadable = false)).reason == "CUDA.jl not loadable"
    @test resolve_backend(:gpu, Float64; probe = withp(functional = false)).reason == "no functional device"
    @test resolve_backend(:gpu, Float64; probe = withp(compiles = false)).reason == "kernel does not compile for the element type"
    @test resolve_backend(:auto, Float64; probe = good, cpu_only_reducer = true).reason == "reducer is CPU-only"
    @test resolve_backend(:auto, Float64; probe = withp(gpu_faster = false)).reason == "GPU not faster for the element type"
    for req in (:auto, :gpu)
        @test resolve_backend(req, Float64; probe = withp(loadable = false)).backend === :cpu
        @test resolve_backend(req, Float64; probe = withp(compiles = false)).backend === :cpu
    end
    # order of checks
    allbad = () -> (loadable = false, functional = false, device = "", compiles = false, gpu_faster = false)
    @test resolve_backend(:auto, Float64; probe = allbad).reason == "CUDA.jl not loadable"
    @test resolve_backend(:auto, Float64; probe = allbad, cpu_only_reducer = true).reason == "reducer is CPU-only"
    @test resolve_backend(:auto, BigFloat; probe = nocall, cpu_only_reducer = true).reason == "kernel does not compile for the element type"
    @test resolve_backend(:auto, Float64; probe = withp(functional = false, compiles = false)).reason == "no functional device"
    @test resolve_backend(:auto, Float64; probe = withp(compiles = false, gpu_faster = false)).reason == "kernel does not compile for the element type"
    # (f) provenance record of a run
    ident = tc_Lin(SMatrix{2,2,Float64}(1, 0, 0, 1))
    x0 = [SVector(0.0, 0.0) for _ in 1:200]
    for req in (:auto, :gpu, :cpu)
        p = propagate(ident, x0, 4; k = 2, tau = 2, seed = seed, backend = req, chunk = 64).provenance
        for name in (:backend, :device, :backend_fallback_reason, :precision, :threads, :cpu, :seed, :tau, :k, :chunk_size)
            @test hasproperty(p, name)
        end
        @test string(p.backend) in ("cpu", "gpu")   # never "auto"
        @test p.device == "cpu" && p.precision == "Float64"
        @test p.seed == seed && p.tau == 2 && p.k == 2 && p.chunk_size == 64
        @test p.threads isa Integer && p.threads >= 1
        @test p.cpu isa AbstractString && !isempty(p.cpu)
        if req === :cpu
            @test p.backend_fallback_reason == ""
        else
            @test p.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device")
        end
    end
end

@testset "AC11 chunk size" begin
    # chunk = floor(f free / (M + reducer bytes)), M = sizeof(T) (d + d k + 2k) + 16; free = 2^30, f = 0.5.
    free = 2^30
    @test bytes_per_orbit(Val(2), Val(2), Float64) == 96
    @test bytes_per_orbit(Val(4), Val(4), Float64) == 240
    @test bytes_per_orbit(Val(2), Val(2), Double64) == 16 * 10 + 16
    @test chunk_size(free, Val(2), Val(2), Float64) == 5592405
    @test chunk_size(free, Val(4), Val(4), Float64) == 2236962
    @test chunk_size(free ÷ 2, Val(2), Val(2), Float64) == fld(2^28, 96)   # halved: 2796202
    @test chunk_size(free ÷ 2, Val(4), Val(4), Float64) == fld(2^28, 240)   # halved: 1118481
    @test chunk_size(free, Val(2), Val(2), Float64; reducer_bytes = 48) == fld(2^29, 96 + 48)
    @test chunk_size(free, Val(2), Val(2), Float64; f = 0.25) == fld(2^28, 96)
    @test chunk_size(free, Val(2), Val(2), Double64) == fld(2^29, 176)
end
