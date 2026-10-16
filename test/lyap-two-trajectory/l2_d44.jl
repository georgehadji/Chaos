# Amendment D44, item 15 and criterion D44k of lyap-two-trajectory: step-chunked relaunch for the two-trajectory estimator.
# Spec: tangent-core-spec.md "Amendment D44" (items 1 to 15, D44k); estimator spec: lyap-two-trajectory-spec.md.
# Interface names fixed by the orchestrator (not present before stage 03): Chaoslab.Tangent._LAUNCH_STEPS (forced L, CPU and device),
# Chaoslab.Tangent._before_launch(device, chunk_index, launch, carried) with carried = (ends = Vector{ShadowEnd},) for propagate_shadow,
# provenance launch_steps and launches_per_chunk, FALLBACK_REASONS[9] == "device lost", Chaoslab.Tangent.DEFAULT_GPU_NS shared.
# Stage 04 (verify/gpu.jl, local GTX 1050): D44k(c) hardware part. Seeds: l2_seed(18, sub) (L2_BASE + 900 + sub, inside this slug's range).
# This file defines its own fakes (prefix l2_): the tangent-core helpers (tc_) are not part of this suite.

const l2_F = Chaoslab.Foundation
const l2_Tg = Chaoslab.Tangent
const l2_CPU = l2_Tg.KernelAbstractions.CPU

# Device on the CPU backend; hook(chunk_index, launch, carried) runs on the host before every device launch.
struct l2_HookDevice <: l2_F.AbstractDevice
    ns::Union{Nothing, Float64}
    hook::Function
end
l2_F.probe(::l2_HookDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = true)
l2_F.free_memory(::l2_HookDevice) = 2^28
l2_F.gpu_step_ns(d::l2_HookDevice, ::Type) = d.ns
l2_Tg.device_backend(::l2_HookDevice) = l2_CPU()
@static if isdefined(l2_Tg, :_before_launch)
    l2_Tg._before_launch(d::l2_HookDevice, ci::Int, li::Int, carried) = d.hook(ci, li, carried)
end

struct l2_FakeLost <: Exception end
struct l2_FakeOOM <: Exception end
l2_F.is_device_lost(::l2_FakeLost) = true
l2_F.is_oom(::l2_FakeOOM) = true

l2_nohook(ci, li, carried) = nothing
l2_fake(ns = nothing) = l2_HookDevice(ns, l2_nohook)

function l2_with_device(f, dev)
    old = l2_F.DEVICE[]
    l2_F.DEVICE[] = dev
    try
        return f()
    finally
        l2_F.DEVICE[] = old
    end
end
function l2_with_launch_steps(f, L)
    old = l2_Tg._LAUNCH_STEPS[]
    l2_Tg._LAUNCH_STEPS[] = L
    try
        return f()
    finally
        l2_Tg._LAUNCH_STEPS[] = old
    end
end

# Bitwise equality of final reference points, shadow points and finals (n, sum, comp, X2, d_prev); isequal tells -0.0 from 0.0.
l2_same(a, b) = isequal(a.x, b.x) && isequal(a.W, b.W) && isequal(a.finals, b.finals)

# Hook that records the launches reached, may corrupt the carried state at one launch and throw at another.
function l2_hook(seen; fail = nothing, corrupt = nothing)
    return function (ci, li, carried)
        push!(seen, (ci, li))
        if corrupt !== nothing && (ci, li) == corrupt
            e = carried.ends[1]   # finite but wrong: the previous separation doubles, so the next stretch ratio is wrong
            carried.ends[1] = Chaoslab.Lyapunov.ShadowEnd(e.x, e.y, e.sum, 2 * e.dprev, e.code, e.step, e.value)
        end
        fail !== nothing && (ci, li) == fail[1] && throw(fail[2])
        return nothing
    end
end

# 40 wrapped standard-map orbits (K = 5), chunks of 16, 16, 8; N = 100, tau = 4, d0 = 1e-8, forced L = 20 (5 launches per chunk).
const L2_44_N, L2_44_TAU, L2_44_CHUNK, L2_44_L, L2_44_D0 = 100, 4, 16, 20, 1e-8
l2_44_x0(seed) = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), 40)]
function l2_44_run(seed; backend = :gpu, L = L2_44_L, N = L2_44_N, tau = L2_44_TAU, x0 = l2_44_x0(seed), sys = l2_sys(:wrapped, 5.0))
    go() = propagate_shadow(sys, x0, N; d0 = L2_44_D0, tau = tau, seed = seed, chunk = L2_44_CHUNK, backend = backend)
    return L === nothing ? go() : l2_with_launch_steps(go, L)   # L = nothing: the run does not touch the hook
end

@testset "D44k(a) CPU bitwise, carried state" begin
    # K = 0.5 and 2 pi, 256 seeded orbits, d0 = 1e-8, tau in {1, 2, 10}, N = 2000, forced L in {tau, 20 tau, 130 tau}:
    # reference points, shadow points, running sums, compensation terms, previous separations and estimates equal the single-launch run
    seed = l2_seed(18, 1)
    d0, N = 1e-8, 2000
    for K in (0.5, L2_TWO_PI)
        sys = l2_sys(:wrapped, K)
        x0 = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), 256)]
        for tau in (1, 2, 10)
            ref = propagate_shadow(sys, x0, N; d0 = d0, tau = tau, seed = seed)
            @test ref.provenance.backend === :cpu
            for mult in (1, 20, 130)
                L = mult * tau
                r = l2_with_launch_steps(L) do
                    propagate_shadow(sys, x0, N; d0 = d0, tau = tau, seed = seed)
                end
                @test l2_same(r, ref)
                @test all(i -> r.finals[i][1].sum == ref.finals[i][1].sum && r.finals[i][1].comp == ref.finals[i][1].comp &&
                               r.finals[i][1].d_prev == ref.finals[i][1].d_prev && r.finals[i][1].X2 == ref.finals[i][1].X2, eachindex(ref.finals))
            end
        end
    end
    @test l2_Tg._LAUNCH_STEPS[] === nothing
    # L must be a positive multiple of tau (CPU backend and device path)
    sys, x0 = l2_sys(:wrapped, 2.0), [l2_state(:wrapped, Float64, SVector(1.0, 2.0))]
    bad(L, backend) = l2_with_launch_steps(L) do
        propagate_shadow(sys, x0, 20; d0 = d0, tau = 4, seed = seed, backend = backend)
    end
    for L in (0, -4, 3, 6)
        @test_throws ArgumentError bad(L, :cpu)
        l2_with_device(l2_fake()) do
            @test_throws ArgumentError bad(L, :gpu)
        end
    end
end

@testset "D44k(b) error step is the global step" begin
    # a pair that collapses to one point at global step 12 (tau = 5, handoff 15, launch 3 for L = 5), and a state that becomes Inf there
    seed = l2_seed(18, 2)
    d0 = 1e-8
    x0 = [SVector(0.5, i == 3 ? 1.0 : 0.0) for i in 1:5]
    U0 = fill(SVector(0.6, 0.8), 5)
    for (sys, what) in ((l2_Fault(:collapse, 12.25), "separation"), (l2_Fault(:inf_both, 12.25), "state"), (l2_Fault(:inf_ref, 12.25), "state"))
        call() = propagate_shadow(sys, x0, 20; d0 = d0, tau = 5, seed = seed, U0 = U0, chunk = 2)
        single = l2_health_error(call)
        @test single isa NumericalHealthError && single.what == what && single.step == 15 && single.orbit == 3 && single.backend === :cpu
        for L in (5, 10, 15, 20)
            err = l2_with_launch_steps(L) do
                l2_health_error(call)
            end
            @test err isa NumericalHealthError && err.step == single.step && err.orbit == 3 && err.what == what && err.backend === :cpu
        end
        # the device path: the device-side error reruns on the CPU, which raises the same error at the same global step
        l2_with_device(l2_fake()) do
            err = l2_with_launch_steps(5) do
                l2_health_error(() -> propagate_shadow(sys, x0, 20; d0 = d0, tau = 5, seed = seed, U0 = U0, chunk = 2, backend = :gpu))
            end
            @test err isa NumericalHealthError && err.step == 15 && err.orbit == 3 && err.what == what && err.backend === :cpu
        end
    end
end

@testset "D44k device path on the CPU backend, bitwise" begin
    # device on the CPU backend: forced L < N equals forced L = N (and the unforced run) bitwise; also equals the CPU run
    seed = l2_seed(18, 3)
    N, d0 = 2000, 1e-8
    for K in (0.5, L2_TWO_PI), tau in (2, 10)
        sys = l2_sys(:wrapped, K)
        x0 = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), 64)]
        run() = propagate_shadow(sys, x0, N; d0 = d0, tau = tau, seed = seed, backend = :gpu)
        l2_with_device(l2_fake()) do
            full = l2_with_launch_steps(run, N)
            @test full.provenance.backend === :gpu && full.provenance.launch_steps == N && full.provenance.launches_per_chunk == 1
            @test l2_same(run(), full)
            @test l2_same(full, propagate_shadow(sys, x0, N; d0 = d0, tau = tau, seed = seed, backend = :cpu))
            for L in (tau, 20 * tau, 130 * tau)
                r = l2_with_launch_steps(run, L)
                @test l2_same(r, full)
                @test r.provenance.backend === :gpu && r.provenance.launch_steps == L && r.provenance.launches_per_chunk == cld(N, L)
            end
        end
    end
end

@testset "D44k(d) budget and the shared default ns value" begin
    seed = l2_seed(18, 4)
    sys = l2_sys(:wrapped, L2_TWO_PI)
    x0 = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), 4)]
    estimate(ns, p) = ns * 1e-9 * p.launch_steps * max(p.chunk_size, 1024)
    @test l2_Tg.DEFAULT_GPU_NS == 100.0
    # no record: the estimator takes the tangent stream's default (L = 9765 for tau = 1, 9760 for tau = 10); N = 20000 splits
    l2_with_device(l2_fake(nothing)) do
        for (tau, L) in ((1, 9765), (10, 9760))
            r = propagate_shadow(sys, x0, 20_000; d0 = 1e-8, tau = tau, seed = seed, backend = :gpu)
            p = r.provenance
            @test p.backend === :gpu && p.backend_fallback_reason == "" && p.gpu_chunks == 1
            @test p.launch_steps == L && p.launches_per_chunk == cld(20_000, L) == 3
            @test estimate(l2_Tg.DEFAULT_GPU_NS, p) <= 1.0 + 1e-9
            @test l2_same(r, propagate_shadow(sys, x0, 20_000; d0 = 1e-8, tau = tau, seed = seed, backend = :cpu))
        end
    end
    # a record of 25 ns: L = 39062 (tau = 1) and 39060 (tau = 10); N = 1e5 gives 3 launches
    l2_with_device(l2_fake(25.0)) do
        for (tau, L) in ((1, 39062), (10, 39060))
            p = propagate_shadow(sys, x0, 100_000; d0 = 1e-8, tau = tau, seed = seed, backend = :gpu).provenance
            @test p.backend === :gpu && p.launch_steps == L && p.launches_per_chunk == cld(100_000, L) == 3
            @test estimate(25.0, p) <= 1.0 + 1e-9
        end
    end
    # even L = tau exceeds the budget for one wave (tau = 10: ns above 97656.25): CPU with the budget reason, provenance 0 and 0, CPU result
    cpu = propagate_shadow(sys, x0, 100; d0 = 1e-8, tau = 10, seed = seed, backend = :cpu)
    l2_with_device(l2_fake(97657.0)) do
        r = propagate_shadow(sys, x0, 100; d0 = 1e-8, tau = 10, seed = seed, backend = :gpu)
        @test l2_same(r, cpu) && r.provenance.backend === :cpu
        @test r.provenance.backend_fallback_reason == "launch exceeds the display watchdog"
        @test r.provenance.gpu_chunks == 0 && r.provenance.launch_steps == 0 && r.provenance.launches_per_chunk == 0
    end
    # a CPU request, unforced: 0 and 0; forced CPU step chunk: its L and cld(N, L) with backend :cpu
    p = propagate_shadow(sys, x0, 100; d0 = 1e-8, tau = 4, seed = seed).provenance
    @test p.launch_steps == 0 && p.launches_per_chunk == 0
    forced = l2_with_launch_steps(36) do
        propagate_shadow(sys, x0, 100; d0 = 1e-8, tau = 4, seed = seed)
    end
    p = forced.provenance
    @test p.backend === :cpu && p.launch_steps == 36 && p.launches_per_chunk == cld(100, 36)
end

@testset "D44k lost device and classified faults in launch j > 1" begin
    seed = l2_seed(18, 5)
    ref = l2_44_run(seed; backend = :cpu, L = nothing)
    for c in 1:3, corrupt in (false, true)   # first, middle and last chunk (16, 16, 8 orbits)
        seen = Tuple{Int, Int}[]
        hook = l2_hook(seen; fail = ((c, 3), l2_FakeLost()), corrupt = corrupt ? (c, 2) : nothing)
        l2_with_device(l2_HookDevice(nothing, hook)) do
            r = l2_44_run(seed)
            @test l2_same(r, ref)
            p = r.provenance
            @test p.backend === :cpu && p.backend_fallback_reason == "device lost" && p.gpu_chunks == c - 1
            @test p.launch_steps == L2_44_L && p.launches_per_chunk == cld(L2_44_N, L2_44_L)   # the device values stay after the fallback
            @test (c, 3) in seen && !((c, 4) in seen) && !any(s -> s[1] > c, seen)
            @test l2_F.DEVICE[] isa l2_F.NoDevice
        end
        # an out-of-memory error in launch 3 reruns the chunk on the CPU with reason 4; the device is kept
        seen = Tuple{Int, Int}[]
        l2_with_device(l2_HookDevice(nothing, l2_hook(seen; fail = ((c, 3), l2_FakeOOM()), corrupt = (c, 2)))) do
            r = l2_44_run(seed)
            @test l2_same(r, ref) && r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "device out of memory"
            @test r.provenance.gpu_chunks == c - 1 && l2_F.DEVICE[] isa l2_HookDevice
        end
    end
    # device, split: the values of items 13 and 14
    l2_with_device(l2_fake()) do
        r = l2_44_run(seed)
        @test l2_same(r, ref) && r.provenance.backend === :gpu
        @test r.provenance.launch_steps == L2_44_L && r.provenance.launches_per_chunk == 5 && r.provenance.gpu_chunks == 3
        p = l2_44_run(seed; L = nothing).provenance
        @test p.launch_steps == L2_44_N && p.launches_per_chunk == 1
    end
    @test FALLBACK_REASONS[9] == "device lost" && length(FALLBACK_REASONS) == 9
end
