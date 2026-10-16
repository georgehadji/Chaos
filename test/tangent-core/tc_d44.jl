# Amendment D44 (step-chunked relaunch): criteria D44a, D44b, D44c (CPU-backend part), D44d, D44e, D44f, D44g, D44h (package part).
# Spec: tangent-core-spec.md, "Amendment D44". Interface names fixed by the orchestrator (they do not exist before stage 03):
#   Chaoslab.Tangent._LAUNCH_STEPS   Ref{Union{Nothing, Int}}, forced launch length L (CPU backend and device path)
#   Chaoslab.Tangent._before_launch(device, chunk_index, launch, carried)   host hook before every device launch
#   Chaoslab.Tangent._budget_steps(ns, tau)   largest multiple of tau within the 1 s budget at one wave, 0 if none
#   provenance fields launch_steps, launches_per_chunk; FALLBACK_REASONS[9] == "device lost"
# Stage 04 (verify/gpu.jl, local GTX 1050): D44c hardware part, D44i (inspection), D44j (payoff). All seeds: TC_SEED + 4400 + sub.
# RED: UndefVarError for _LAUNCH_STEPS / _budget_steps, missing provenance keys, 8 reasons, :cpu where a relaunch is now expected.

const TC44_SEED = TC_SEED + 4400
const TC44_ORBITS = 256   # D44a orbits (spec: 1024, reduced for the package suite runtime; the test plan records it)

# Set the forced launch length for the duration of f, restore it afterwards.
function tc_with_launch_steps(f, L)
    old = tc_Tg._LAUNCH_STEPS[]
    tc_Tg._LAUNCH_STEPS[] = L
    try
        return f()
    finally
        tc_Tg._LAUNCH_STEPS[] = old
    end
end

# Bitwise equality of points, deviation matrices and reducer finals (isequal: NaN equals NaN, -0.0 differs from 0.0).
tc44_same(a, b) = isequal(a.x, b.x) && isequal(a.W, b.W) && isequal(a.finals, b.finals)

# Orbits whose single-launch run raises no health error (K = 2 pi with tau = 10 makes the Gram-Schmidt reducer break down on a few orbits):
# returns the kept starts and their seeded initial deviation matrices, passed as W0 so that the orbit ids keep their seeded frames.
function tc44_clean(sys, x0, N, k, tau, red, seed)
    keep = Int[]
    for j in eachindex(x0)
        try
            propagate_orbit(sys, x0[j], initial_deviations(seed, j, Val(2), Val(k), Float64), N, tau, red; orbit = j)
            push!(keep, j)
        catch e
            e isa NumericalHealthError || rethrow()
        end
    end
    return x0[keep], [initial_deviations(seed, j, Val(2), Val(k), Float64) for j in keep]
end

# Standard-map variant with a step counter (state length 3) whose tangent is the zero matrix when the old point has counter `at`:
# the QR spectrum reducer then fails at the next handoff (a reducer failure with a known global step).
struct tc44_ZeroTan <: tc_AbstractMap
    K::Float64
    at::Int
end
Chaoslab.dimension(::tc44_ZeroTan) = 2
function Chaoslab.forward(m::tc44_ZeroTan, z)
    p2 = z[2] + m.K * sin(z[1])
    return SVector(z[1] + p2, p2, z[3] + 1)
end
function Chaoslab.tangent(m::tc44_ZeroTan, z)
    c = m.K * cos(z[1])
    return z[3] == m.at ? zero(SMatrix{2, 2, Float64, 4}) : SMatrix{2, 2}(1 + c, c, one(c), one(c))
end

# A device whose launch hook runs a test closure (hook(chunk_index, launch, carried)); the backend is the KernelAbstractions CPU backend.
struct tc44_HookDevice <: tc_F.AbstractDevice
    ns::Union{Nothing, Float64}
    hook::Function
end
tc_F.probe(::tc44_HookDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = true)
tc_F.free_memory(::tc44_HookDevice) = 2^28
tc_F.gpu_step_ns(d::tc44_HookDevice, ::Type) = d.ns
tc_Tg.device_backend(::tc44_HookDevice) = tc_CPU()
@static if isdefined(tc_Tg, :_before_launch)
    tc_Tg._before_launch(d::tc44_HookDevice, ci::Int, li::Int, carried) = d.hook(ci, li, carried)
end

# A device whose free_memory call reports a lost device (setup failure).
struct tc44_FreeLostDevice <: tc_F.AbstractDevice end
tc_F.probe(::tc44_FreeLostDevice, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = true, gpu_faster = true)
tc_F.free_memory(::tc44_FreeLostDevice) = throw(tc_FakeLost())

# Hook builders. `seen` records every (chunk, launch) the run reached.
function tc44_hook(seen; fail = nothing, corrupt = nothing)
    return function (ci, li, carried)
        push!(seen, (ci, li))
        corrupt !== nothing && (ci, li) == corrupt && tc44_corrupt!(carried)
        fail !== nothing && (ci, li) == fail[1] && throw(fail[2])
        return nothing
    end
end
# Finite but wrong carried state: moves the first point and doubles the first deviation matrix.
function tc44_corrupt!(carried)
    carried.x[1] = carried.x[1] .+ 0.25
    carried.W[1] = carried.W[1] * 2
    return nothing
end

@testset "D44 interface" begin
    @test tc_Tg._LAUNCH_STEPS isa Base.RefValue{Union{Nothing, Int}}
    @test tc_Tg._LAUNCH_STEPS[] === nothing
    @test isdefined(tc_Tg, :_before_launch) && isdefined(tc_Tg, :_budget_steps)
    sys = tc_Std(2.0)
    x0 = [SVector(0.1, 0.2)]
    # L must be a positive multiple of tau, on the CPU backend and the device path
    bad(L, backend) = tc_with_launch_steps(L) do
        propagate(sys, x0, 20; k = 1, tau = 4, seed = TC44_SEED, backend = backend)
    end
    for L in (0, -4, 3, 6, 10)
        @test_throws ArgumentError bad(L, :cpu)
        tc_with_device(tc_FakeDevice(true)) do
            @test_throws ArgumentError bad(L, :gpu)
        end
    end
    @test tc_Tg._LAUNCH_STEPS[] === nothing   # restored after the throwing calls
end

@testset "D44a CPU bitwise, carried state" begin
    # 256 seeded orbits (spec 1024; orbits whose single-launch run raises a health error are dropped, see tc44_clean), K = 2 pi, N = 2000; k = 1 Benettin, k = 2 QRSpectrum(RGS()), k = 2 SVDSpectrum;
    # tau in {1, 2, 10}, forced L in {tau, 20 tau, 130 tau}: x, W and finals equal the single-launch run bitwise (tolerance 0).
    seed = TC44_SEED + 1
    N = 2000
    sys = StandardMap(2 * pi)
    rng = StableRNG(seed)
    x0 = [SVector(2 * pi * rand(rng), 2 * pi * rand(rng)) for _ in 1:TC44_ORBITS]
    cases = (("Benettin", 1, (Benettin(),)), ("QRSpectrum(RGS())", 2, (QRSpectrum(RGS()),)), ("SVDSpectrum", 2, (SVDSpectrum(),)))
    for (name, k, red) in cases, tau in (1, 2, 10)
        xc, W0 = tc44_clean(sys, x0, N, k, tau, red, seed)
        @test length(xc) >= 64
        ref = propagate(sys, xc, N; k = k, tau = tau, seed = seed, reducers = red, W0 = W0)
        @test ref.provenance.backend === :cpu
        for mult in (1, 20, 130)
            L = mult * tau
            r = tc_with_launch_steps(L) do
                propagate(sys, xc, N; k = k, tau = tau, seed = seed, reducers = red, W0 = W0)
            end
            @test tc44_same(r, ref)
            @test r.provenance.backend === :cpu
        end
    end
    @test tc_Tg._LAUNCH_STEPS[] === nothing
end

@testset "D44b health error step" begin
    seed = TC44_SEED + 2
    # state: Inf at global step 12, tau = 5: the error is raised at handoff 15, whatever L (step 12 lies in launch 3 for L = 5)
    x1 = [SVector(0.0, 0.0, j == 3 ? 1.0 : 0.0) for j in 1:6]
    for (sys, what) in ((tc_Fault(:inf_state, 12), "state"), (tc_Fault(:nan_W, 12), "deviation matrix"))
        single = tc_health_error(() -> propagate(sys, x1, 20; k = 1, tau = 5, seed = seed, chunk = 2))
        @test single isa NumericalHealthError && single.step == 15 && single.orbit == 3 && single.what == what && single.backend === :cpu
        for L in (5, 10, 15, 20)
            err = tc_with_launch_steps(L) do
                tc_health_error(() -> propagate(sys, x1, 20; k = 1, tau = 5, seed = seed, chunk = 2))
            end
            @test err isa NumericalHealthError && err.step == single.step == 15 && err.orbit == 3 && err.what == what && err.backend === :cpu
        end
        # the device path: the device-side error reruns on the CPU, which raises the same error at the same global step
        tc_with_device(tc_FakeDevice(true)) do
            err = tc_with_launch_steps(5) do
                tc_health_error(() -> propagate(sys, x1, 20; k = 1, tau = 5, seed = seed, chunk = 2, backend = :gpu))
            end
            @test err isa NumericalHealthError && err.step == 15 && err.orbit == 3 && err.what == what && err.backend === :cpu
        end
    end
    # reducer failure in launch 2: the tangent is zero at the old point with counter 7 (step 8, handoff 10); the QR spectrum fails there
    zsys = tc44_ZeroTan(2.0, 7)
    xz = [SVector(0.3 + 0.1 * j, 0.2, 0.0) for j in 1:5]
    red = (QRSpectrum(RGS()),)
    single = tc_health_error(() -> propagate(zsys, xz, 20; k = 2, tau = 5, seed = seed, reducers = red, chunk = 2))
    @test single isa NumericalHealthError && single.step == 10 && single.orbit == 1
    for L in (5, 10, 15)
        err = tc_with_launch_steps(L) do
            tc_health_error(() -> propagate(zsys, xz, 20; k = 2, tau = 5, seed = seed, reducers = red, chunk = 2))
        end
        @test err isa NumericalHealthError && err.step == single.step && err.orbit == single.orbit
        @test err.what == single.what && err.column == single.column && err.backend === :cpu
    end
    tc_with_device(tc_FakeDevice(true)) do
        err = tc_with_launch_steps(5) do
            tc_health_error(() -> propagate(zsys, xz, 20; k = 2, tau = 5, seed = seed, reducers = red, chunk = 2, backend = :gpu))
        end
        @test err isa NumericalHealthError && err.step == single.step && err.orbit == single.orbit && err.backend === :cpu
    end
end

@testset "D44c device path on the CPU backend, bitwise" begin
    # With a test device whose backend is the CPU backend: forced L < N equals forced L = N (and the unforced device run) bitwise.
    seed = TC44_SEED + 3
    N, m = 2000, 128
    sys = StandardMap(2 * pi)
    rng = StableRNG(seed)
    x0 = [SVector(2 * pi * rand(rng), 2 * pi * rand(rng)) for _ in 1:m]
    cases = (("Benettin", 1, (Benettin(),)), ("QRSpectrum(RGS())", 2, (QRSpectrum(RGS()),)), ("SVDSpectrum", 2, (SVDSpectrum(),)))
    tc_with_device(tc_FakeDevice(true)) do
        for (name, k, red) in cases, tau in (2, 10)
            xc, W0 = tc44_clean(sys, x0, N, k, tau, red, seed)
            run() = propagate(sys, xc, N; k = k, tau = tau, seed = seed, reducers = red, W0 = W0, backend = :gpu)
            full = tc_with_launch_steps(run, N)
            @test full.provenance.backend === :gpu && full.provenance.launch_steps == N && full.provenance.launches_per_chunk == 1
            @test tc44_same(run(), full)
            for L in (tau, 20 * tau, 130 * tau)
                r = tc_with_launch_steps(run, L)
                @test tc44_same(r, full)
                @test r.provenance.backend === :gpu && r.provenance.launch_steps == L && r.provenance.launches_per_chunk == cld(N, L)
            end
        end
    end
end

@testset "D44d budget" begin
    seed = TC44_SEED + 4
    # the budget function, exact values (D44d)
    @test tc_Tg._budget_steps(25.0, 1) == 39062 && tc_Tg._budget_steps(25.0, 10) == 39060
    @test tc_Tg._budget_steps(100.0, 1) == 9765 && tc_Tg._budget_steps(100.0, 10) == 9760
    @test tc_Tg._budget_steps(97657.0, 10) == 0   # just above 1 / (10 x 1024 x 1e-9) = 97656.25
    @test tc_Tg._budget_steps(25.0, 1) isa Int
    for ns in (3.0, 20.2, 25.0, 100.0, 5000.0), tau in (1, 2, 10)
        @test tc_Tg._budget_steps(ns, tau) == tc44_budget_L(ns, tau)
    end
    @test tc_Tg.DEFAULT_GPU_NS == 100.0
    # runs on a fake device: the launch length comes from the record (or the default), the estimate stays within the budget
    sys = tc_Std(5.0)
    x0 = [SVector(0.1 * j, 0.3) for j in 1:4]   # tc_LogNorm renormalizes W: a chaotic K = 5 orbit overflows an unrenormalized W within ~1000 steps
    estimate(ns, p) = ns * 1e-9 * p.launch_steps * max(p.chunk_size, 1024)
    for (rec, ns) in ((25.0, 25.0), (nothing, 100.0)), tau in (1, 10)
        tc_with_device(tc_FakeDevice(true, rec, 2^28)) do
            L = rec === nothing ? (tau == 1 ? 9765 : 9760) : (tau == 1 ? 39062 : 39060)
            N = 4 * L + 3 * tau   # not a multiple of L: the last launch is shorter
            r = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = (tc_LogNorm(),), backend = :gpu)
            p = r.provenance
            @test p.backend === :gpu && p.backend_fallback_reason == "" && p.gpu_chunks > 0
            @test p.launch_steps == L && p.launches_per_chunk == 5
            @test estimate(ns, p) <= 1.0 + 1e-9
            @test tc44_same(r, propagate(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = (tc_LogNorm(),), backend = :cpu))
        end
    end
    # N = 1e6 with the record of 25 ns: 26 launches per chunk for tau = 10 (and for tau = 1), device chunk count above 0
    tc_with_device(tc_FakeDevice(true, 25.0, 2^28)) do
        for tau in (1, 10)
            p = propagate(sys, x0, 1_000_000; k = 1, tau = tau, seed = seed, reducers = (tc_LogNorm(),), backend = :gpu).provenance
            @test p.backend === :gpu && p.gpu_chunks > 0 && p.launches_per_chunk == 26
            @test p.launch_steps == (tau == 1 ? 39062 : 39060) && estimate(25.0, p) <= 1.0 + 1e-9
        end
    end
    # even L = tau exceeds the budget for one wave: CPU with the budget reason, result equals the CPU run, provenance 0 and 0
    red = (tc_LogNorm(), tc_Reader(10))
    cpu = propagate(sys, x0, 100; k = 2, tau = 10, seed = seed, reducers = red, backend = :cpu)
    tc_with_device(tc_FakeDevice(true, 97657.0, 2^28)) do
        r = propagate(sys, x0, 100; k = 2, tau = 10, seed = seed, reducers = red, backend = :gpu)
        @test tc44_same(r, cpu) && r.provenance.backend === :cpu && r.provenance.device == "cpu"
        @test r.provenance.backend_fallback_reason == "launch exceeds the display watchdog"
        @test r.provenance.gpu_chunks == 0 && r.provenance.launch_steps == 0 && r.provenance.launches_per_chunk == 0
    end
end

# Shared fixture of D44e, D44f, D44g: 40 orbits in chunks of 16, 16, 8; N = 100, tau = 4, forced L = 20 (5 launches per chunk).
const TC44_N, TC44_TAU, TC44_CHUNK, TC44_L = 100, 4, 16, 20
tc44_x0() = [SVector(0.1 * j, 0.3) for j in 1:40]
tc44_red() = (tc_LogNorm(), tc_Reader(TC44_TAU))
function tc44_run(sys, seed; backend = :gpu, L = TC44_L)
    go() = propagate(sys, tc44_x0(), TC44_N; k = 2, tau = TC44_TAU, seed = seed, reducers = tc44_red(), chunk = TC44_CHUNK, backend = backend)
    return L === nothing ? go() : tc_with_launch_steps(go, L)   # L = nothing: the run does not touch the hook
end

@testset "D44e lost device in launch j > 1" begin
    seed = TC44_SEED + 5
    sys = tc_Std(5.0)
    ref = propagate(sys, tc44_x0(), TC44_N; k = 2, tau = TC44_TAU, seed = seed, reducers = tc44_red(), chunk = TC44_CHUNK, backend = :cpu)
    for c in 1:3, corrupt in (false, true)   # first, middle, last chunk; with and without a corrupted carried state at launch 2
        seen = Tuple{Int, Int}[]
        hook = tc44_hook(seen; fail = ((c, 3), tc_FakeLost()), corrupt = corrupt ? (c, 2) : nothing)
        tc_with_device(tc44_HookDevice(nothing, hook)) do
            r = tc44_run(sys, seed)
            @test tc44_same(r, ref)
            p = r.provenance
            @test p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == "device lost"
            @test p.gpu_chunks == c - 1
            @test (c, 3) in seen && !((c, 4) in seen) && !any(s -> s[1] > c, seen)   # the failed chunk is not relaunched, the rest stays on the CPU
            @test tc_F.DEVICE[] isa tc_F.NoDevice
            @test tc44_run(sys, seed).provenance.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device")
        end
    end
    # a lost device during setup is also "device lost"
    tc_with_device(tc44_FreeLostDevice()) do
        r = tc44_run(sys, seed)
        @test tc44_same(r, ref) && r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == "device lost"
        @test r.provenance.gpu_chunks == 0 && tc_F.DEVICE[] isa tc_F.NoDevice
    end
end

@testset "D44f other classified faults in launch j > 1" begin
    seed = TC44_SEED + 6
    sys = tc_Std(5.0)
    ref = propagate(sys, tc44_x0(), TC44_N; k = 2, tau = TC44_TAU, seed = seed, reducers = tc44_red(), chunk = TC44_CHUNK, backend = :cpu)
    for c in 1:3
        for (err, reason) in ((tc_FakeOOM(), "device out of memory"), (tc_FakeCompileError(), "kernel does not compile for the element type"))
            seen = Tuple{Int, Int}[]
            tc_with_device(tc44_HookDevice(nothing, tc44_hook(seen; fail = ((c, 3), err), corrupt = (c, 2)))) do
                r = tc44_run(sys, seed)
                @test tc44_same(r, ref)
                @test r.provenance.backend === :cpu && r.provenance.backend_fallback_reason == reason && r.provenance.gpu_chunks == c - 1
                @test tc_F.DEVICE[] isa tc44_HookDevice   # only a lost device resets the run-wide device
            end
        end
        # a health error that only the device run raises: the map fails once, armed by the hook at launch 3
        flaky = tc_FlakyStd(5.0, Threads.Atomic{Bool}(false))
        fref = propagate(flaky, tc44_x0(), TC44_N; k = 2, tau = TC44_TAU, seed = seed, reducers = tc44_red(), chunk = TC44_CHUNK, backend = :cpu)
        arm = (ci, li, carried) -> (ci == c && li == 3 && (flaky.armed[] = true); nothing)
        tc_with_device(tc44_HookDevice(nothing, arm)) do
            r = tc44_run(flaky, seed)
            @test tc44_same(r, fref)
            @test r.provenance.backend === :cpu && r.provenance.gpu_chunks == c - 1
            @test r.provenance.backend_fallback_reason == "NumericalHealthError raised only on the GPU"
        end
    end
end

@testset "D44g provenance: launch_steps and launches_per_chunk" begin
    seed = TC44_SEED + 7
    sys = tc_Std(5.0)
    N, tau, chunk, L = TC44_N, TC44_TAU, TC44_CHUNK, TC44_L
    # a CPU request, unforced: 0 and 0
    p = tc44_run(sys, seed; backend = :cpu, L = nothing).provenance
    @test p.launch_steps == 0 && p.launches_per_chunk == 0 && p.backend === :cpu
    # device, single launch: N and 1
    tc_with_device(tc_FakeDevice(true)) do
        p = tc44_run(sys, seed; L = nothing).provenance
        @test p.backend === :gpu && p.launch_steps == N && p.launches_per_chunk == 1
        # device, split: L and cld(N, L), the last launch shorter when L does not divide N
        for L in (tau, L, 36, 2tau * 3)
            p = tc44_run(sys, seed; L = L).provenance
            @test p.backend === :gpu && p.launch_steps == L && p.launches_per_chunk == cld(N, L)
        end
    end
    # forced CPU step chunk: its L and cld(N, L), backend stays :cpu
    for L in (tau, 20, 36)
        p = tc44_run(sys, seed; backend = :cpu, L = L).provenance
        @test p.backend === :cpu && p.launch_steps == L && p.launches_per_chunk == cld(N, L)
    end
    # a fallback before the first launch (setup failure, budget): 0 and 0
    for dev in (tc_FreeOOMDevice(), tc_AllocOOMDevice(), tc44_FreeLostDevice())
        tc_with_device(dev) do
            p = tc44_run(sys, seed; L = nothing).provenance
            @test p.backend === :cpu && p.gpu_chunks == 0 && p.launch_steps == 0 && p.launches_per_chunk == 0
        end
    end
    tc_with_device(tc_FakeDevice(true, 97657.0, 2^28)) do
        p = propagate(sys, tc44_x0(), 100; k = 1, tau = 10, seed = seed, backend = :gpu).provenance
        @test p.backend === :cpu && p.launch_steps == 0 && p.launches_per_chunk == 0
    end
    # after a mid-run fallback the device values stay (chunk 2 fails in launch 3 after chunk 1 ran on the device)
    for err in (tc_FakeLost(), tc_FakeOOM())
        tc_with_device(tc44_HookDevice(nothing, tc44_hook(Tuple{Int, Int}[]; fail = ((2, 3), err)))) do
            p = tc44_run(sys, seed).provenance
            @test p.backend === :cpu && p.gpu_chunks == 1
            @test p.launch_steps == L && p.launches_per_chunk == cld(N, L)
        end
    end
end

@testset "D44h ninth reason" begin
    @test length(FALLBACK_REASONS) == 9 && FALLBACK_REASONS[9] == "device lost"
    @test FALLBACK_REASONS[8] == "launch exceeds the display watchdog"
    @test length(unique(FALLBACK_REASONS)) == 9
    # a lost device (tc_FakeLost has is_device_lost true) gives the ninth reason and never the budget string; the budget case only the budget string
    seed = TC44_SEED + 8
    sys = tc_Std(5.0)
    tc_with_device(tc44_HookDevice(nothing, tc44_hook(Tuple{Int, Int}[]; fail = ((1, 2), tc_FakeLost())))) do
        @test tc44_run(sys, seed).provenance.backend_fallback_reason == FALLBACK_REASONS[9]
        @test tc_F.DEVICE[] isa tc_F.NoDevice
    end
    tc_with_device(tc_FakeDevice(true, 97657.0, 2^28)) do
        r = propagate(sys, tc44_x0(), 100; k = 1, tau = 10, seed = seed, backend = :gpu)
        @test r.provenance.backend_fallback_reason == FALLBACK_REASONS[8]
    end
    # ChaoslabCUDAExt table (codes 700, 702, ...): stage 04, verify/gpu.jl
end
