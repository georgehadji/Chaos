# Amendment D29: device hook and per-chunk fallback. No GPU in CI: a FakeDevice runs the device path on the CPU backend.
# Spec: tangent-core-spec.md, "Amendment D29".
const tc_CPU = tc_Tg.KernelAbstractions.CPU

# `ns` is the recorded GPU time per orbit-step (nothing: no record), `free` the free device memory in bytes.
struct tc_FakeDevice <: tc_F.AbstractDevice
    faster::Bool
    ns::Union{Nothing, Float64}
    free::Int
end
tc_FakeDevice(faster::Bool) = tc_FakeDevice(faster, nothing, 2^28)
tc_F.probe(d::tc_FakeDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = d.faster)
tc_F.free_memory(d::tc_FakeDevice) = d.free
tc_F.gpu_step_ns(d::tc_FakeDevice, ::Type) = d.ns
tc_Tg.device_backend(::tc_FakeDevice) = tc_CPU()

# Devices that fail before the first chunk: free_memory throws an out-of-memory error / the buffer allocation does.
struct tc_FreeOOMDevice <: tc_F.AbstractDevice end
tc_F.probe(::tc_FreeOOMDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = true)
tc_F.free_memory(::tc_FreeOOMDevice) = throw(tc_FakeOOM())
struct tc_AllocOOMDevice <: tc_F.AbstractDevice end
struct tc_BadBackend <: tc_Tg.KernelAbstractions.Backend end
tc_F.probe(::tc_AllocOOMDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = true)
tc_F.free_memory(::tc_AllocOOMDevice) = 2^28
tc_Tg.device_backend(::tc_AllocOOMDevice) = tc_BadBackend()
tc_Tg.KernelAbstractions.allocate(::tc_BadBackend, ::Type, ::Vararg{Integer}) = throw(tc_FakeOOM())

struct tc_FakeOOM <: Exception end
struct tc_FakeCompileError <: Exception end
struct tc_FakeOther <: Exception end
struct tc_FakeLost <: Exception end
tc_F.is_oom(::tc_FakeOOM) = true
tc_F.is_compile_error(::tc_FakeCompileError) = true
tc_F.is_device_lost(::tc_FakeLost) = true

# W0 vector that throws `err` on its `at`-th element access (a fault inside the device chunk, on the host side).
mutable struct tc_FailAt{M} <: AbstractVector{M}
    W::Vector{M}
    err::Union{Exception, Nothing}
    at::Int
    calls::Int
end
Base.size(v::tc_FailAt) = size(v.W)
function Base.getindex(v::tc_FailAt, i::Int)
    v.calls += 1
    v.calls == v.at && v.err !== nothing && throw(v.err)
    return v.W[i]
end

# Standard map whose first `forward` call anywhere returns NaN (the flag is shared), so only the first run of a chunk fails.
struct tc_FlakyStd <: tc_AbstractMap
    K::Float64
    armed::Threads.Atomic{Bool}
end
Chaoslab.dimension(::tc_FlakyStd) = 2
function Chaoslab.forward(m::tc_FlakyStd, z)
    Threads.atomic_xchg!(m.armed, false) && return SVector(NaN, NaN)
    p2 = z[2] + m.K * sin(z[1])
    return SVector(z[1] + p2, p2)
end
function Chaoslab.tangent(m::tc_FlakyStd, z)
    c = m.K * cos(z[1])
    return SMatrix{2, 2}(1 + c, c, one(c), one(c))
end

# Reducer whose state is not isbits (BigFloat): cannot run in a GPU kernel.
struct tc_HeapSum end
Chaoslab.reducer_init(::tc_HeapSum, ::Type{T}, ::Val, ::Val) where {T} = (big(0.0),)
Chaoslab.reducer_update(::tc_HeapSum, st, n::Int, x, W) = ((st[1] + n,), W)
Chaoslab.reducer_final(::tc_HeapSum, st) = (sum = st[1],)

@testset "D29 device hook and fallbacks" begin
    dev0 = tc_F.DEVICE[]
    seed = TC_SEED + 29
    n, N, tau, chunk = 40, 100, 4, 16   # 3 chunks of 16, 16, 8
    sys = tc_Std(5.0)
    x0 = [SVector(0.1 * j, 0.3) for j in 1:n]
    red = (tc_LogNorm(), tc_Reader(tau))   # D38: the chain of AC6, so the kernel specialization is shared
    run(s = sys; kw...) = propagate(s, x0, N; k = 2, tau = tau, seed = seed, reducers = red, chunk = chunk, kw...)
    same(a, b) = a.x == b.x && a.W == b.W && a.finals == b.finals
    ref = run(backend = :cpu)
    ref0 = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :cpu)   # default chunk
    W0s = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n]
    failing(err, at) = tc_FailAt(copy(W0s), err, at, 0)

    @testset "defaults and hooks" begin
        @test tc_with_device(() -> tc_F.DEVICE[], tc_F.NoDevice()) isa tc_F.NoDevice
        @test tc_Tg.device_backend(tc_F.NoDevice()) isa tc_CPU
        @test tc_Tg.workgroup(tc_CPU()) == 16
        @test tc_F.free_memory(tc_F.NoDevice()) > 0
        @test !tc_F.is_oom(ErrorException("x")) && !tc_F.is_compile_error(ErrorException("x")) && !tc_F.is_device_lost(ErrorException("x"))
        @test tc_F.gpu_step_ns(tc_F.NoDevice(), Float64) === nothing
        @test length(FALLBACK_REASONS) == 8 && FALLBACK_REASONS[8] == "launch exceeds the display watchdog"
        p = tc_F.probe(tc_F.NoDevice(), Float64)
        @test keys(p) == (:loadable, :functional, :device, :compiles, :gpu_faster) && !p.functional && !p.gpu_faster
        for name in (:AbstractDevice, :NoDevice, :DEVICE, :probe, :free_memory, :is_oom, :is_compile_error, :is_device_lost, :gpu_step_ns)
            @test !(name in names(tc_F))
        end
    end

    @testset "device path equals the CPU path bitwise" begin
        tc_with_device(tc_FakeDevice(true)) do
            @test resolve_backend(:gpu, Float64).device == "FakeGPU"
            @test resolve_backend(:auto, Float64).backend === :gpu
            for req in (:gpu, :auto)
                r = run(backend = req)
                @test same(r, ref)
                @test r.provenance.backend === :gpu && r.provenance.device == "FakeGPU" && r.provenance.backend_fallback_reason == ""
            end
            @test same(run(backend = :cpu), ref)
        end
        # chunk size comes from the device's free memory (a fast recorded device, so the launch budget does not bind)
        tc_with_device(tc_FakeDevice(true, 1.0, 2^28)) do
            r = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :gpu)
            s0 = tc_Tg._init_states(red, Float64, Val(2), Val(2))
            @test r.provenance.chunk_size == chunk_size(2^28, Val(2), Val(2), Float64; reducer_bytes = sizeof(s0), s = 2)
        end
        tc_with_device(tc_FakeDevice(false)) do
            r = run(backend = :auto)
            @test same(r, ref) && r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "GPU not faster for the element type"
            @test run(backend = :gpu).provenance.backend === :gpu
        end
    end

    @testset "reason 6: reducer state is not isbits" begin
        h = (tc_LogNorm(), tc_HeapSum())
        base = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = h, chunk = chunk, backend = :cpu)
        tc_with_device(tc_FakeDevice(true)) do
            r = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = h, chunk = chunk, backend = :gpu)
            @test r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "reducer is CPU-only"
            @test r.x == base.x && r.finals == base.finals
        end
    end

    @testset "reasons 4 and 3: a chunk fails on the device" begin
        for (err, reason) in ((tc_FakeOOM(), "device out of memory"), (tc_FakeCompileError(), "kernel does not compile for the element type"))
            tc_with_device(tc_FakeDevice(true)) do
                # first chunk, and a later chunk: the failed chunk is rerun on the CPU, the rest of the run stays there
                for at in (1, 2chunk + 3, 4chunk + 3)   # W0 is read twice per orbit
                    r = run(backend = :gpu, W0 = failing(err, at))
                    @test same(r, ref)
                    @test r.provenance.backend === :cpu && r.provenance.device == "cpu" && r.provenance.backend_fallback_reason == reason
                    @test r.provenance.gpu_chunks == (at == 1 ? 0 : at == 2chunk + 3 ? 1 : 2)   # chunks that ran on the device before the switch
                end
            end
        end
    end

    @testset "reason 5: health error on the device chunk only" begin
        tc_with_device(tc_FakeDevice(true)) do
            flaky = tc_FlakyStd(5.0, Threads.Atomic{Bool}(true))
            r = run(flaky; backend = :gpu)
            @test same(r, ref)
            @test r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "NumericalHealthError raised only on the GPU"
            # an error that the CPU raises too propagates, from the CPU rerun
            x1 = [SVector(0.0, 0.0, j == 3 ? 1.0 : 0.0) for j in 1:6]
            err = tc_health_error(() -> propagate(tc_Fault(:inf_state, 7), x1, 20; k = 1, tau = 5, seed = seed, chunk = 2, backend = :gpu))
            @test err isa NumericalHealthError && err.orbit == 3 && err.step == 10 && err.backend === :cpu
        end
    end

    @testset "launch budget (display watchdog)" begin
        budget_ns(orbits) = 1e9 / (orbits * N)   # gpu_ns_per_step at which one launch of `orbits` orbits takes the 1 s budget
        # the budget caps the device chunk
        tc_with_device(tc_FakeDevice(true, budget_ns(1100.5), 2^28)) do
            r = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :gpu)
            @test same(r, ref0) && r.provenance.backend === :gpu && r.provenance.chunk_size == 1100
            r = run(backend = :gpu)   # an explicit chunk below the budget is kept
            @test r.provenance.backend === :gpu && r.provenance.chunk_size == chunk
        end
        # N alone exceeds the budget: even one wave of 1024 orbits does not fit
        tc_with_device(tc_FakeDevice(true, budget_ns(500), 2^28)) do
            r = run(backend = :gpu)
            @test same(r, ref) && r.provenance.backend === :cpu && r.provenance.device == "cpu"
            @test r.provenance.backend_fallback_reason == "launch exceeds the display watchdog" && r.provenance.gpu_chunks == 0
            @test run(backend = :cpu).provenance.backend_fallback_reason == ""
        end
        # no record: the documented default of 100 ns per orbit-step applies
        @test tc_Tg.DEFAULT_GPU_NS == 100.0
        tc_with_device(tc_FakeDevice(true)) do
            @test propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :gpu).provenance.backend === :gpu
            r = propagate(sys, x0, 10_000; k = 2, tau = 10, seed = seed, reducers = red, backend = :gpu)   # 1 s / (100 ns x 1e4) = 1000 orbits < 1024
            @test r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "launch exceeds the display watchdog"
        end
    end

    @testset "lost device" begin
        tc_with_device(tc_FakeDevice(true)) do
            r = run(backend = :gpu, W0 = failing(tc_FakeLost(), 2chunk + 3))
            @test same(r, ref) && r.provenance.backend === :cpu && r.provenance.gpu_chunks == 1
            @test r.provenance.backend_fallback_reason == "launch exceeds the display watchdog"
            @test tc_F.DEVICE[] isa tc_F.NoDevice   # later calls do not touch the dead context
            @test run(backend = :gpu).provenance.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device")
        end
    end

    @testset "out of memory before the first chunk" begin
        for dev in (tc_FreeOOMDevice(), tc_AllocOOMDevice())
            tc_with_device(dev) do
                r = run(backend = :gpu)
                @test same(r, ref) && r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "device out of memory"
                @test r.provenance.gpu_chunks == 0
            end
        end
    end

    @testset "provenance after a fallback" begin
        tc_with_device(tc_FakeDevice(true, nothing, 2^20)) do
            dchunk = chunk_size(2^20, Val(2), Val(2), Float64; reducer_bytes = sizeof(tc_Tg._init_states(red, Float64, Val(2), Val(2))), s = 2)
            r = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :gpu, W0 = failing(tc_FakeOOM(), 1))
            @test r.provenance.backend === :cpu && r.provenance.gpu_chunks == 0 && r.provenance.gpu_chunk_size == dchunk
            @test r.provenance.chunk_size > dchunk   # recomputed from host memory after the switch
            r = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = red, backend = :gpu)
            @test r.provenance.backend === :gpu && r.provenance.gpu_chunks == 1 && r.provenance.chunk_size == dchunk
        end
        @test ref.provenance.gpu_chunks == 0 && ref.provenance.gpu_chunk_size == 0
    end

    @testset "other failures propagate" begin
        tc_with_device(tc_FakeDevice(true)) do
            @test_throws tc_FakeOther run(backend = :gpu, W0 = failing(tc_FakeOther(), 5))
        end
    end

    @test tc_F.DEVICE[] === dev0
end
