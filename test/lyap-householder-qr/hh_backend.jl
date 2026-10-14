# AC15: backend resolution and provenance on CI (no GPU: fake devices drive the resolution), exact strings.

const HH_PROV = Dict("HouseholderQR" => ("V-H", "kernel"), "AugmentedQR" => ("V-A", "kernel"),
                     "GivensQR" => ("V-G", "kernel"), "LapackQR" => ("V-H", "library"))

@testset "AC15 backend resolution and provenance" begin
    seed = hh_seed(15)
    sys = StandardMap(2.0)
    x0 = hh_starts(StableRNG(seed), 8)
    run15(alg, req; extra = ()) = propagate(sys, x0, 20; k = 2, tau = 2, seed = seed, reducers = (QRSpectrum(alg), extra...), backend = req)
    none, present = hh_F.NoDevice(), hh_FakeDevice(true, true)

    # (a) library path: "reducer is CPU-only" for auto and gpu, with and without a device present
    for dev in (none, present), req in (:auto, :gpu)
        p = hh_with_device(() -> run15(LapackQR(), req), dev).provenance
        @test p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == "reducer is CPU-only"
    end
    # (b) kernel path without a GPU runtime: the two exact reasons for a runtime that cannot be loaded and for no functional device
    p = hh_with_device(() -> run15(HouseholderQR(), :gpu), hh_FakeDevice(false, false)).provenance
    @test p.backend === :cpu && p.backend_fallback_reason == "CUDA.jl not loadable"
    p = hh_with_device(() -> run15(HouseholderQR(), :gpu), hh_FakeDevice(true, false)).provenance
    @test p.backend === :cpu && p.backend_fallback_reason == "no functional device"
    p = hh_with_device(() -> run15(hh_pick(GivensQR(), HouseholderQR()), :auto), none).provenance   # CI default: either string, depending on whether CUDA.jl is findable
    @test p.backend === :cpu && p.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device")
    # (c) BigFloat with a gpu request: element-type reason in both paths (P9 precedence)
    setprecision(BigFloat, 128) do
        sysB = StandardMap(BigFloat(2))
        xB = [SVector{2, BigFloat}(BigFloat(0.5) + j, BigFloat(0.25)) for j in 1:3]
        for dev in (none, present)
            p = hh_with_device(() -> propagate(sysB, xB, 4; k = 2, tau = 2, seed = seed, reducers = (QRSpectrum(HouseholderQR()),), backend = :gpu), dev).provenance
            @test p.backend === :cpu && p.backend_fallback_reason == "kernel does not compile for the element type"
            @test resolve_backend(:gpu, BigFloat; cpu_only_reducer = true).reason == "kernel does not compile for the element type"
        end
    end
    @test_throws ArgumentError reducer_init(QRSpectrum(LapackQR()), Double64, Val(2), Val(2))   # library path is Float64 only (Q8)
    # (d) a cpu request records an empty reason for every variant
    for alg in hh_pick(hh_algs(Float64), (HouseholderQR(), LapackQR()))     # package suite: V-H kernel and library path
        @test run15(alg, :cpu).provenance.backend_fallback_reason == ""
    end
    # (e) provenance: estimator keys (exact values, order free) merged with tangent-core's fields; backend is never auto
    for alg in hh_pick(hh_algs(Float64), (HouseholderQR(), LapackQR()))
        name, impl = HH_PROV[hh_name(alg)]
        expected = (estimator = "qr", orthonormalizer = name, iteration_path = "diagonal", implementation = impl, orthogonality_tol = 8.0)
        @test hh_nt_eq(reducer_provenance(QRSpectrum(alg)), expected)
        for res in (run15(alg, :cpu), run15(alg, :auto; extra = (ResidualMonitor(sys, 4.0, Inf),)))
            p = res.provenance
            @test all(k -> getproperty(p, k) == getproperty(expected, k), keys(expected))
            @test p.backend in (:cpu, :gpu) && p.k == 2 && p.tau == 2 && p.seed == seed && p.precision == "Float64"
            @test all(k -> haskey(p, k), (:device, :backend_fallback_reason, :threads, :cpu, :chunk_size, :headroom, :gpu_chunks, :gpu_chunk_size))
        end
    end
    @test reducer_provenance(QRSpectrum(HouseholderQR(); orthogonality_tol = 4)).orthogonality_tol === 4.0
end
