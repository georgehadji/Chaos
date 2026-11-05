# Stage 04 GPU verification (decision D29): every STAGE04-GPU line of the tangent-core, precision, standard-map and
# diagnostics verify logs. Run from the repo root on a machine with a CUDA device:
#   julia --project=verify -t 4 verify/gpu.jl > workspaces/module-build/stages/04-verify/output/gpu-verify-log.txt
# One line per criterion:  AC<n> | quantity | value | bound | PASS/FAIL/FALLBACK(<reason>)
# FALLBACK(r): the element type did not run on the GPU (probe or per-chunk fallback, reason r); that is a valid result,
# not a kernel pass (precision spec criterion 9). Every GPU-side line checks provenance.backend == :gpu, so a silent
# fallback can never make a CPU-vs-CPU comparison pass.
# It also writes data/gpu-benchmark.toml (per machine, git-ignored) from the AC13 measurement before AC12a, the 1e5-step run.
using Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, Printf, Statistics, TOML, Serialization

const CHILD = length(ARGS) >= 2 && ARGS[1] == "child"
CHILD || @eval using CUDA

const SEED = 20261012
const ROOT = normpath(joinpath(@__DIR__, ".."))
# GPU_ONLY="svd" (comma-separated tags record, tc, precision, diagnostics, standard-map, svd, gs, hh, hh-perf, l2, d44, d44c, sali, fli, megno, rli, stretching, naff, weighted-birkhoff, coupled-maps) runs only those sections; the
# default runs everything. Skipping "record" needs an existing data/gpu-benchmark.toml for the launch budget.
const ONLY = get(ENV, "GPU_ONLY", "")
want(tag) = isempty(ONLY) || tag in strip.(split(ONLY, ','))
const ROWS = Tuple{String,String}[]

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const MODTAG = Ref("")   # module whose log line this row answers
function row(id, qty, val, bound, verdict)
    push!(ROWS, (string(id), verdict))
    println(id, " | ", MODTAG[], qty, " | ", fmt(val), " | ", fmt(bound), " | ", verdict)
    flush(stdout)
end
pf(b::Bool) = b ? "PASS" : "FAIL"
fb(reason) = "FALLBACK($reason)"
epsf(::Type{T}) where {T} = Float64(eps(T))

# Column-normalising reducer with compensated log sums (same as verify/tangent-core.jl).
struct LogNorm end
function Chaoslab.reducer_init(::LogNorm, ::Type{T}, ::Val{d}, ::Val{k}) where {T,d,k}
    return ntuple(_ -> CompensatedSum{T}(), Val(k))
end
function Chaoslab.reducer_update(::LogNorm, st, n::Int, x, W::SMatrix{d,k,T}) where {d,k,T}
    s = ntuple(j -> sqrt(sum(abs2, ntuple(i -> W[(j - 1) * d + i], Val(d)))), Val(k))
    st2 = ntuple(j -> add(st[j], log(s[j])), Val(k))
    Wn = SMatrix{d,k,T}(ntuple(i -> W[i] / s[(i - 1) ÷ d + 1], Val(d * k)))
    return st2, Wn
end
Chaoslab.reducer_final(::LogNorm, st) = (logs = map(total, st),)

# AC13 ensemble: d = 2, k = 2, K = 2 pi, LogNorm, tau = 1 (as verify/tangent-core.jl). `backend` is :cpu or :gpu.
function ac13_run(T, norb, N, backend; chunk = nothing)
    seed = SEED + 13
    sys = StandardMap(T(2pi))
    rng = StableRNG(seed)
    x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:norb]
    kw = (k = 2, tau = 1, seed = seed, reducers = (LogNorm(),), backend = backend, chunk = chunk)
    propagate(sys, x0[1:64], 10; kw..., chunk = nothing)   # warm-up (compilation)
    t = @elapsed res = propagate(sys, x0, N; kw...)
    return t, res.provenance
end

# CPU sizes: 4 threads at the size of verify/tangent-core.jl; the 1-thread child runs smaller ensembles (ns per orbit-step is
# size independent to a few percent) so that the whole script stays within minutes.
const AC13_SIZES = Dict(Float64 => (10^5, 10^3), Double64 => (2000, 10^3))
const AC13_SIZES_1T = Dict(Float64 => (20_000, 10^3), Double64 => (500, 10^3))
const AC13_GPU_CHUNK = 8_000

# lyap-householder-qr AC14(c), K15: wall time of the kernel path chain at N = 1e4, K = 2 pi, tau = 2; warm-up excluded.
function hh_perf_run(backend, no)
    rng = StableRNG(SEED + 50000 + 1015)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:no]
    kw = (k = 2, tau = 2, seed = SEED, reducers = (QRSpectrum(HouseholderQR()),), backend = backend)
    propagate(StandardMap(2pi), x0[1:64], 10; kw...)
    t = @elapsed res = propagate(StandardMap(2pi), x0, 10^4; kw...)
    return t, res.provenance
end
if length(ARGS) >= 2 && ARGS[1] == "child-hh"
    serialize(ARGS[2], first(hh_perf_run(:cpu, 10^6)))
    exit(0)
end

if CHILD
    res = Dict{Any,Any}(:nthreads => Threads.nthreads())
    for (T, (no, N)) in AC13_SIZES_1T
        res[T] = first(ac13_run(T, no, N, :cpu)) / (no * N) * 1e9
    end
    serialize(ARGS[2], res)
    exit(0)
end

include(joinpath(ROOT, "test", "precision", "helpers.jl"))      # pr_*: BigFloat reference orbit, M_n, Gamma_n, n_H
include(joinpath(ROOT, "test", "standard-map", "sm_helpers.jl")) # sm_*: sigma, n_H(K), tau(n), torus distance
include(joinpath(ROOT, "test", "diagnostics", "faultmaps.jl"))   # dg_FaultMap, dg_Renorm

const F = Chaoslab.Foundation
const DEV = F.DEVICE[]
const DEVNAME = DEV isa F.NoDevice ? "" : CUDA.name(CUDA.device())
const RECORD = joinpath(pkgdir(Chaoslab), "data", "gpu-benchmark.toml")

function tdr_values()
    try
        out = read(`reg query HKLM\\SYSTEM\\CurrentControlSet\\Control\\GraphicsDrivers`, String)
        return join((strip(l) for l in split(out, '\n') if occursin("Tdr", l) && occursin("REG_DWORD", l)), "; ")
    catch
        return "unreadable"
    end
end

println("# gpu verification, julia ", VERSION, ", threads ", Threads.nthreads(), ", cpu ", Sys.CPU_NAME, ", base seed ", SEED)
println("# device type ", typeof(DEV), ", name '", DEVNAME, "'")
if !(DEV isa F.NoDevice)
    println("# CUDA.jl ", pkgversion(CUDA), ", runtime ", CUDA.runtime_version(), ", driver ", CUDA.driver_version(),
            ", capability ", CUDA.capability(CUDA.device()), ", free memory ", round(CUDA.free_memory() / 2^20), " MiB",
            ", stack limit ", CUDA.limit(CUDA.LIMIT_STACK_SIZE), " B (extension minimum D40)")
end
println("# Windows TDR registry values (GraphicsDrivers): ", tdr_values())
println("# record file: ", RECORD, " (", isfile(RECORD) ? "exists" : "absent", " at start)")
for T in (Float64, Double64)
    println("# probe ", T, ": ", F.probe(DEV, T))
end
println("# format: AC<n> | quantity | value | bound | verdict")
if DEV isa F.NoDevice
    println("# no functional CUDA device: nothing to verify")
    exit(1)
end

# Reason a GPU request for T ends on the CPU ("" when it runs on the device), from the probe.
gpu_reason(T) = F.resolve_backend(:gpu, T).reason

# propagate on the GPU; returns (result, reason) with reason == nothing when the run stayed on the device.
function on_gpu(args...; kw...)
    res = propagate(args...; backend = :gpu, kw...)
    p = res.provenance
    return res, p.backend === :gpu ? nothing : p.backend_fallback_reason
end

catchhe(f) = try; f(); nothing; catch e; e isa NumericalHealthError ? e : rethrow(); end
sdiff(a, b) = maximum(abs, a .- b)                      # max-norm difference of two states
slogs(res, N) = [f[1].logs[1] / N for f in res.finals]   # S/n of the first deviation vector

# AC13 and the record come first: the launch budget of propagate reads gpu_ns_per_step from the record, so a machine without one
# must measure it before the long-N lines.
# ------------------------------------------------------------------------------------------ AC13 (tangent-core) and record
MODTAG[] = "tangent-core: "
gpu_ns = Dict{Any,Float64}()
cpu_ns = Dict{Any,Float64}()
cpu1_ns = Dict{Any,Float64}()
want("record") && let
    tmp = tempname()
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "verify")) -t 1 $(@__FILE__) child $tmp`
    println("# starting 1-thread CPU child run (AC13 sizes: ", AC13_SIZES_1T, ")")
    run(cmd)
    child = deserialize(tmp)
    rm(tmp; force = true)
    BLAS.set_num_threads(1)
    today = Libc.strftime("%Y-%m-%d", time())
    for T in (Float64, Double64)
        no, N = AC13_SIZES[T]
        steps = no * N
        scale = (10^6 * 10^4) / steps
        tN, _ = ac13_run(T, no, N, :cpu)
        cpu_ns[T] = tN / steps * 1e9
        cpu1_ns[T] = child[T]
        reason = gpu_reason(T)
        entry = Dict{String,Any}("cpu_ns_per_step" => cpu_ns[T], "threads" => Threads.nthreads(), "date" => today)
        if isempty(reason)
            tg = Inf; prov = nothing
            for _ in 1:2
                t, prov = ac13_run(T, no, N, :gpu; chunk = AC13_GPU_CHUNK)
                tg = min(tg, t)
            end
            prov.backend === :gpu ? (gpu_ns[T] = tg / steps * 1e9) : (reason = prov.backend_fallback_reason)
        end
        if isempty(reason)
            # D44 item 9: the launch budget's value at one wave (1024 orbits x 1e4 steps, one launch), the slower of two runs
            tw = 0.0
            old = Chaoslab.Tangent._LAUNCH_STEPS[]
            Chaoslab.Tangent._LAUNCH_STEPS[] = 10^4
            try
                for _ in 1:2
                    t, pw = ac13_run(T, 1024, 10^4, :gpu)
                    pw.backend === :gpu || error("one-wave record run fell back: $(pw.backend_fallback_reason)")
                    (pw.gpu_chunks == 1 && pw.gpu_chunk_size >= 1024 && pw.launches_per_chunk == 1) ||
                        error("one-wave record run is not one chunk of one launch: $pw")
                    tw = max(tw, t)
                end
            finally
                Chaoslab.Tangent._LAUNCH_STEPS[] = old
            end
            entry["gpu_ns_per_step_wave"] = tw / (1024 * 10^4) * 1e9
            row("D44", "record: GPU ns/orbit-step at one wave, T=$T (1024 orbits x 1e4 steps, one launch, slower of 2)", fmt(entry["gpu_ns_per_step_wave"]), "recorded", "INFO")
            entry["gpu_ns_per_step"] = gpu_ns[T]
            entry["gpu_faster"] = gpu_ns[T] <= cpu_ns[T]
            for (nt, base) in ((1, cpu1_ns[T]), (Threads.nthreads(), cpu_ns[T]))
                row("AC13", "GPU wall s, T=$T, $no orbits x $N steps (d=2, k=2, K=2pi, LogNorm reducer, tau=1, chunk $AC13_GPU_CHUNK): $(fmt(gpu_ns[T])) ns/orbit-step vs CPU $nt thread(s) $(fmt(base)) ns/orbit-step; speedup $(fmt(base / gpu_ns[T])); extrapolated to 1e6 x 1e4: $(fmt(gpu_ns[T] * 1e-9 * 1e10)) s",
                    gpu_ns[T] * 1e-9 * steps, "GPU<=CPU: " * fmt(base * 1e-9 * steps) * " s", pf(gpu_ns[T] <= base))
            end
        else
            entry["gpu_faster"] = false
            entry["reason"] = reason
            for nt in (1, Threads.nthreads())
                row("AC13", "GPU wall s, T=$T, $no orbits x $N steps, vs CPU $nt thread(s): $(fmt(nt == 1 ? cpu1_ns[T] : cpu_ns[T])) ns/orbit-step",
                    "-", "GPU<=CPU", fb(reason))
            end
        end
        # the record (per machine, git-ignored): [<device name>."<precision label>"]
        rec = isfile(RECORD) ? TOML.parsefile(RECORD) : Dict{String,Any}()
        get!(rec, DEVNAME, Dict{String,Any}())[precision_label(T)] = entry
        mkpath(dirname(RECORD))
        tmp = RECORD * ".tmp"   # temp file then move, so a reader or a crash never sees half a record
        open(io -> TOML.print(io, rec), tmp, "w")
        mv(tmp, RECORD; force = true)
        println("# record ", precision_label(T), " -> ", entry)
    end
end

# ------------------------------------------------------------------------------------------------------ AC12b (tangent-core)
MODTAG[] = "tangent-core: "
want("tc") && for K in (2pi, 10.0)
    G = sqrt((1 + K)^2 + K^2 + 2)
    bound = 4 * 10 * epsf(Float64) * G^10
    rng = StableRNG(SEED + 1202)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:1024]
    sys = StandardMap(K)
    cpu = propagate(sys, x0, 10; k = 1, seed = SEED, backend = :cpu)
    gpu, why = on_gpu(sys, x0, 10; k = 1, seed = SEED)
    d = maximum(sdiff(a, b) for (a, b) in zip(cpu.x, gpu.x))
    row("AC12b", "K=$(K == 2pi ? "2pi" : "10") n=10 same-orbit |x_cpu - x_gpu|, max over 1024 orbits (bound 4 n eps G^n, G=$(fmt(G)))",
        d, bound, why === nothing ? pf(d <= bound) : fb(why))
end

# ------------------------------------------------------------------------------------------------------ AC12c (tangent-core)
want("tc") && let
    K, no, N, tau = 10.0, 1024, 10^4, 10
    rng = StableRNG(SEED + 1203)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:no]
    sys = StandardMap(K)
    kw = (k = 1, tau = tau, seed = SEED, reducers = (LogNorm(),))
    cpu = propagate(sys, x0, N; kw..., backend = :cpu)
    gpu, why = on_gpu(sys, x0, N; kw...)
    sc, sg = slogs(cpu, N), slogs(gpu, N)
    sig = std(sc)
    maxd = maximum(abs.(sc .- sg)) / sig
    se = sqrt((var(sc) + var(sg)) / no)
    dm = abs(mean(sc) - mean(sg)) / se
    row("AC12c", "chaotic ensemble K=10, $no orbits, n=1e4, tau=$tau: max per-orbit |S_cpu/n - S_gpu/n| / sigma_FT (sigma_FT=$(fmt(sig))); |difference of means| / SE; mean S/n cpu $(fmt(mean(sc))) gpu $(fmt(mean(sg)))",
        "$(fmt(maxd)); $(fmt(dm))", "6; 4", why === nothing ? pf(maxd <= 6 && dm <= 4) : fb(why))
end

# --------------------------------------------------------------------------------------------------------- AC9 (precision)
MODTAG[] = "precision: "
want("precision") && let
    ics = pr_ics(PR_SEED + 4, 200)
    refs = setprecision(BigFloat, 256) do
        [pr_reference(BigFloat(5), x0, p0, PR_NMAX) for (x0, p0) in ics]
    end
    MG = [pr_bounds(r) for r in refs]
    for T in (Float64, Double64)
        epsT = epsf(T)
        nH = [pr_nH(M, G, epsT) for (M, G) in MG]
        sys = StandardMap{T}(exact(T, 5))
        reason = gpu_reason(T)
        if !isempty(reason)
            row("AC9", "GPU vs CPU orbit agreement, $T, 200 seeded orbits K=5 (n_H min/median/max $(minimum(nH))/$(median(nH))/$(maximum(nH)))",
                "-", "2 c eps M_n Gamma_n", fb(reason))
            continue
        end
        worst, why = 0.0, nothing
        for cap in (10, typemax(Int)), n in unique(min.(nH, cap))
            idx = findall(==(n), min.(nH, cap))
            x0 = [SVector{2,T}(T(ics[j][1]), T(ics[j][2])) for j in idx]
            cpu = propagate(sys, x0, n; k = 1, seed = SEED, backend = :cpu)
            gpu, w = on_gpu(sys, x0, n; k = 1, seed = SEED)
            w === nothing || (why = w; break)
            for (i, j) in enumerate(idx)
                M, G = MG[j]
                bound = 2 * PR_C * epsT * M[n + 1] * G[n + 1]
                worst = max(worst, Float64(maximum(abs, cpu.x[i] .- gpu.x[i])) / bound)
            end
        end
        row("AC9", "GPU vs CPU orbit agreement, $T, 200 seeded orbits K=5, at n=min(n_H,10) and n=n_H (n_H min/median/max $(minimum(nH))/$(median(nH))/$(maximum(nH))): max |x_cpu - x_gpu| / (2 c eps M_n Gamma_n)",
            worst, "1 (c = 8)", why === nothing ? pf(worst <= 1) : fb(why))
    end
end

# ------------------------------------------------------------------------------------------------------- AC11 (diagnostics)
MODTAG[] = "diagnostics: "
want("diagnostics") && let
    DET_TOL, ORTH_TOL = 4.0, 1.0e-6
    above(theta) = [exp2(j) for j in -130:20 if 2theta <= exp2(j) < 4theta]
    below(theta) = [exp2(j) for j in -130:20 if theta / 4 < exp2(j) <= theta / 2]
    # item 1: integer symplectic maps, residual maxima exactly zero on both backends
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        if !isempty(reason)
            row("AC11", "item 1 integer maps x1000 steps, $T: residual maxima GPU vs CPU", "-", "0 (exact)", fb(reason))
            continue
        end
        mxc, mxg, why = zero(T), zero(T), nothing
        for ent in ((1, 0, 1, 1), (2, 1, 1, 1), (0, -1, 1, 0), (1, 1, -1, 0))
            sys = LinearSymplecticMap(SMatrix{2,2,T}(ent...))
            mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
            kw = (k = 1, tau = 1, seed = SEED + 11, reducers = (dg_Renorm(), mon), W0 = [SMatrix{2,1,T}(1, 0)])
            cpu = propagate(sys, [SVector{2,T}(0, 0)], 1000; kw..., backend = :cpu)
            gpu, w = on_gpu(sys, [SVector{2,T}(0, 0)], 1000; kw...)
            w === nothing || (why = w; break)
            f(r) = (g = r.finals[1][2]; max(g.det_max, g.symp_max, g.orth_max))
            mxc, mxg = max(mxc, f(cpu)), max(mxg, f(gpu))
        end
        row("AC11", "item 1 integer maps x1000 steps, $T: max residual GPU / CPU", "$(Float64(mxg)) / $(Float64(mxc))", 0.0,
            why === nothing ? pf(mxg == 0 && mxc == 0) : fb(why))
    end
    # item 5: det fault on both sides of the threshold. A raise on the GPU is a health error there, so the chunk reruns on
    # the CPU (reason 5) and the CPU raises the error that reaches the caller (backend cpu); a quiet run must stay on the GPU.
    A0 = SMatrix{2,2,Float64}(2, 1, 1, 1)
    theta = Float64(det_threshold(DET_TOL, Float64, A0))
    function run5(delta, N, backend)
        sys = dg_FaultMap(SMatrix{2,2,Float64}(2 + delta, 1, 1, 1))
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        res = nothing
        e = catchhe(() -> (res = propagate(sys, [SVector(0.5, 0.25)], N; k = 1, tau = 1, seed = SEED + 11,
                                            reducers = (dg_Renorm(), mon), backend = backend, W0 = [SMatrix{2,1,Float64}(1, 0)])))
        return e, res
    end
    ups = [(run5(d, 1, :gpu), run5(d, 1, :cpu)) for d in above(theta)]
    dns = [(run5(d, 1000, :gpu), run5(d, 1000, :cpu)) for d in below(theta)]
    okU = count(u -> u[1][1] isa NumericalHealthError && u[1][1].what == "det" && u[2][1] isa NumericalHealthError && u[2][1].what == "det", ups)
    okD = count(u -> u[1][1] === nothing && u[1][2].provenance.backend === :gpu && u[2][1] === nothing &&
                     abs(u[1][2].finals[1][2].det_max - u[2][2].finals[1][2].det_max) <= theta, dns)
    row("AC11", "item 5 det fault (theta=$(fmt(theta))): raise on both backends (GPU raises via CPU rerun) / quiet on the GPU with det_max within theta of the CPU", "$okU/$(length(ups)), $okD/$(length(dns))", "all",
        pf(okU == length(ups) > 0 && okD == length(dns) > 0))
    # item 6: standard-map panel, residual maxima GPU vs CPU within the threshold in force
    ratio, errs, why = 0.0, 0, nothing
    for (iK, K) in enumerate((0.5, 0.9716354, 2.0, 5.0, 10.0))
        sys = StandardMap(K)
        mon = ResidualMonitor(sys, DET_TOL, ORTH_TOL)
        rng = StableRNG(SEED + 1100 + iK)
        x0s = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:50]
        W0s = [(a = 2pi * rand(rng); SMatrix{2,1,Float64}(cos(a), sin(a))) for _ in 1:50]
        kw = (k = 1, tau = 1, seed = SEED + 11, reducers = (dg_Renorm(), mon), W0 = W0s)
        cpu = nothing; gpu = nothing; w = nothing
        e = catchhe(() -> (cpu = propagate(sys, x0s, 10_000; kw..., backend = :cpu); (gpu, w) = on_gpu(sys, x0s, 10_000; kw...)))
        e === nothing || (errs += 1; continue)
        w === nothing || (why = w; break)
        thr = Float64(det_threshold(DET_TOL, Float64, SMatrix{2,2,Float64}(1 + K, K, 1, 1)))
        ratio = max(ratio, maximum(abs(a[2].det_max - b[2].det_max) for (a, b) in zip(cpu.finals, gpu.finals)) / thr)
    end
    row("AC11", "item 6 panel (50 starts x 1e4 steps per K): health errors; max |det_max GPU - CPU| / det threshold in force", "$errs; $(fmt(ratio))", "0; 1",
        why === nothing ? pf(errs == 0 && ratio <= 1) : fb(why))
end

# --------------------------------------------------------------------------------------------------- AC18 (standard-map)
MODTAG[] = "standard-map: "
want("standard-map") && let
    rng = StableRNG(SM_SEED + 40181)
    lim = 4 * Float64(pi)
    Z = [SVector((2 * rand(rng) - 1) * lim, (2 * rand(rng) - 1) * lim) for _ in 1:100_000]
    I2 = SMatrix{2,2,Float64}(1.0, 0.0, 0.0, 1.0)
    W0 = fill(I2, length(Z))
    for K in (1.0, 5.0)
        sys = StandardMap(K); n = sm_nH(K)
        kw = (k = 2, tau = 1, seed = SM_SEED, W0 = W0)
        cpu = propagate(sys, Z, n; kw..., backend = :cpu)
        gpu, why = on_gpu(sys, Z, n; kw...)
        ratio = maximum(sm_tdist(cpu.x[j], gpu.x[j]) / sm_tau(sm_eps(Float64), K, sm_orbit(sys, Z[j], n), n) for j in eachindex(Z))
        dmax = maximum(sm_tdist(cpu.x[j], gpu.x[j]) for j in 1:2000)
        row("AC18", "GPU vs CPU, K=$K n=n_H=$n, $(length(Z)) orbits: max over orbits of torus distance / tau(n) (max distance over the first 2000: $(fmt(dmax)))",
            ratio, "1", why === nothing ? pf(ratio <= 1) : fb(why))
    end
    reason = gpu_reason(Double64)
    if isempty(reason)
        sys = StandardMap(Double64(5)); n = sm_nH(5.0)
        ZD = [SVector{2,Double64}(z) for z in Z[1:20_000]]
        kw = (k = 2, tau = 1, seed = SM_SEED, W0 = fill(SMatrix{2,2,Double64}(1, 0, 0, 1), length(ZD)))
        cpu = propagate(sys, ZD, n; kw..., backend = :cpu)
        gpu, why = on_gpu(sys, ZD, n; kw...)
        ratio = maximum(sm_tdist(cpu.x[j], gpu.x[j]) / sm_tau(sm_eps(Double64), 5.0, sm_orbit(sys, ZD[j], n), n) for j in eachindex(ZD))
        row("AC18", "Double64 kernel gate: GPU vs CPU, K=5 n=$n, 20000 orbits: max torus distance / tau(n)", ratio, "1", why === nothing ? pf(ratio <= 1) : fb(why))
    else
        probe_run = on_gpu(StandardMap(Double64(5)), [SVector{2,Double64}(1, 2)], 2; k = 1, seed = SM_SEED)[2]
        row("AC18", "Double64 kernel gate: compile or reason (explicit :gpu run ends on the CPU with the reason recorded in provenance)",
            probe_run === nothing ? "ran on GPU" : probe_run, "compile or reason", fb(reason))
    end
end

# ------------------------------------------------------------------------------------------------- AC11 (precision), from AC13
MODTAG[] = "precision: "
if !want("record")
elseif haskey(gpu_ns, Double64) && haskey(gpu_ns, Float64)
    row("AC11", "GPU Double64/Float64 time ratio, same ensemble as AC13 (per orbit-step)", gpu_ns[Double64] / gpu_ns[Float64], "reported, no threshold", "PASS")
else
    row("AC11", "GPU Double64/Float64 time ratio, same ensemble as AC13: Float64 GPU $(fmt(get(gpu_ns, Float64, NaN))) ns/orbit-step; CPU ratio $(fmt(cpu_ns[Double64] / cpu_ns[Float64])) (4 threads, this run)",
        "-", "reported, no threshold", fb(gpu_reason(Double64)))
end
println("# after the record: resolve_backend(:auto, Float64) = ", F.resolve_backend(:auto, Float64), "; (:auto, Double64) = ", F.resolve_backend(:auto, Double64))

# ---------------------------------------------------- Batch B: tangent-core AC12 tolerances on the Lyapunov estimators
# One protocol for lyap-svd-stabilized AC10, lyap-gram-schmidt AC12, lyap-householder-qr AC14(a) and lyap-two-trajectory AC17:
# (a) regular class K = 0.5 (lyap-benettin S7 screening: Benettin X(1e5) <= 5e-4 on the CPU, failing starts replaced by the next
# id), D33 symmetric form against a Double64 CPU reference, at n = 4e4 (the largest one-wave launch within the budget, as AC12a),
# tau = 10; (c) chaotic ensemble K = 10, 1024 orbits, n = 1e4, tau = 2 (RGS window), Float64 and Double64 where the kernel
# compiles: max per-orbit |cpu - gpu| <= 6 sigma_FT and |difference of means| <= 4 SE. `est(sys, x0, N, tau, backend)` returns a
# PropagationResult, `lam(final)` the exponents of one orbit.
const REGULAR = Ref{Any}(nothing)
function regular_starts(no)
    REGULAR[] === nothing || return REGULAR[]
    rng = StableRNG(SEED + 50000 + 1001)
    cand = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:3 * no]
    X = [f[1].lambda[1] for f in propagate(StandardMap(0.5), cand, 10^5; k = 1, tau = 10, seed = SEED, reducers = (Benettin(),), backend = :cpu).finals]
    reg = findall(<=(5e-4), X)
    m = min(no, length(reg))
    return REGULAR[] = (cand[reg[1:m]], length(cand), reg[m] - m)   # starts, candidates, non-regular candidates skipped
end
lams(res, lam) = [lam(f[1]) for f in res.finals]
gpu_why(res) = res.provenance.backend === :gpu ? nothing : res.provenance.backend_fallback_reason

function ac12_regular(ac, name, est, lam, nl; extra = nothing)
    x0, ncand, nrep = regular_starts(1024)
    N, tau = 4 * 10^4, 10
    gpu = est(StandardMap(0.5), x0, N, tau, :gpu)
    why = gpu_why(gpu)
    why === nothing || return row(ac, "(a) K=0.5 regular class, $(length(x0)) orbits, n=$(fmt(N)): $name, D33 symmetric form", "-", "count diff <= 3; gap <= 2e-6", fb(why))
    cpu = est(StandardMap(0.5), x0, N, tau, :cpu)
    ref = est(StandardMap(Double64(0.5)), [SVector{2,Double64}(z) for z in x0], N, tau, :cpu)
    for i in 1:nl
        lc, lg, lr = [l[i] for l in lams(cpu, lam)], [l[i] for l in lams(gpu, lam)], [Float64(l[i]) for l in lams(ref, lam)]
        dcr, dgr, dcg = abs.(lc .- lr), abs.(lg .- lr), abs.(lc .- lg)
        nc, ng = count(>(1e-6), dcr), count(>(1e-6), dgr)
        both = (dcr .<= 1e-6) .& (dgr .<= 1e-6)
        gap = any(both) ? maximum(dcg[both]) : Inf
        row(ac, "(a) K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced), n=$(fmt(N)), tau=$tau, $name lambda_$i: orbits above 1e-6 against the Double64 reference CPU / GPU; max CPU-GPU gap over the $(count(both)) orbits where both agree (overall max |cpu-gpu| $(fmt(maximum(dcg))))",
            "$nc / $ng (diff $(abs(nc - ng))); $(fmt(gap))", "count diff <= 3; gap <= 2e-6", pf(abs(nc - ng) <= 3 && gap <= 2e-6))
    end
    extra === nothing || extra(StandardMap(0.5), x0, N, tau, cpu, gpu)
end

# Returns (sys, x0, N, cpu, gpu) when both ran, else nothing (after the FALLBACK or health row).
function ac12_chaotic(ac, name, est, lam, nl, ::Type{T}) where {T}
    no, N, tau = 1024, 10^4, 2
    label = "(c) K=10 chaotic ensemble, $T: $name, 6 sigma_FT per orbit and 4 SE for the means"
    reason = gpu_reason(T)
    isempty(reason) || (row(ac, label, "-", "6 sigma; 4 SE", fb(reason)); return nothing)
    rng = StableRNG(SEED + 50000 + 1003)
    x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:no]
    sys = StandardMap(exact(T, 10))
    cpu = nothing; gpu = nothing
    e = catchhe(() -> (cpu = est(sys, x0, N, tau, :cpu); gpu = est(sys, x0, N, tau, :gpu)))
    e === nothing || (row(ac, "(c) K=10, $T, $name: health error", sprint(showerror, e), "none", "FAIL"); return nothing)
    why = gpu_why(gpu)
    why === nothing || (row(ac, label, "-", "6 sigma; 4 SE", fb(why)); return nothing)
    for i in 1:nl
        sc, sg = [Float64(l[i]) for l in lams(cpu, lam)], [Float64(l[i]) for l in lams(gpu, lam)]
        sig = std(sc)
        maxd = maximum(abs.(sc .- sg)) / sig
        dm = abs(mean(sc) - mean(sg)) / sqrt((var(sc) + var(sg)) / no)
        row(ac, "(c) K=10, $T, $no orbits, n=1e4, tau=$tau, $name lambda_$i: max per-orbit |cpu - gpu| / sigma_FT (sigma_FT=$(fmt(sig))); |difference of means| / SE; means cpu $(fmt(mean(sc))) gpu $(fmt(mean(sg)))",
            "$(fmt(maxd)); $(fmt(dm))", "6; 4", pf(maxd <= 6 && dm <= 4))
    end
    return (sys, x0, N, cpu, gpu)
end

qr_run(r; k = 2) = (sys, x0, N, tau, b) -> propagate(sys, x0, N; k = k, tau = tau, seed = SEED, reducers = (r,), backend = b)
lam_of(f) = f.lambda

# AC10 (lyap-svd-stabilized): the protocol above on lambda_i^SVD, plus the c check (spec Amendment D41): on the orbits with
# S_1(500) >= 25, max |c_cpu - c_gpu| / max(1, |c_cpu|) <= the max over the same orbits of the relative CPU change of c under a
# 1-ulp shift of x0[1] (c is ill-conditioned on some orbits, so a fixed bound does not hold).
MODTAG[] = "lyap-svd-stabilized: "
svdf(res, f) = [getproperty(r[1], f) for r in res.finals]
want("svd") && let
    est = qr_run(SVDSpectrum())
    ac12_regular("AC10", "SVDSpectrum", est, lam_of, 2)
    for T in (Float64, Double64)
        out = ac12_chaotic("AC10", "SVDSpectrum", est, lam_of, 2, T)
        out === nothing && continue
        sys, x0, N, cpu, gpu = out
        no = length(x0)
        S500 = [Float64(f[1].sum[1] + f[1].comp[1]) for f in est(sys, x0, 500, 2, :cpu).finals]
        sel = S500 .>= 25
        cc = Float64.(svdf(cpu, :c))
        dc = abs.(cc .- Float64.(svdf(gpu, :c))) ./ max.(1, abs.(cc))
        cp = Float64.(svdf(est(sys, [SVector(nextfloat(x[1]), x[2]) for x in x0], N, 2, :cpu), :c))
        dp = abs.(cc .- cp) ./ max.(1, abs.(cc))
        bnd = any(sel) ? maximum(dp[sel]) : NaN
        row("AC10", "(c) K=10, $T: max |c_cpu - c_gpu| / max(1, |c_cpu|) at n=1e4 on the $(count(sel)) of $no orbits with S_1(500) >= 25, against the max over the same orbits of the CPU change under a 1-ulp shift of x0[1] (D41; others reported: max $(any(.!sel) ? fmt(maximum(dc[.!sel])) : "-"); max |c| $(fmt(maximum(abs, cc))))",
            any(sel) ? maximum(dc[sel]) : NaN, bnd, pf(any(sel) && maximum(dc[sel]) <= bnd))
        w = argmax(dc)
        row("AC10", "(c) K=10, $T: conditioning of c, CPU change under a 1-ulp shift of x0[1] (relative form): max; median; at the orbit of the largest CPU-GPU difference (#$w, |c| $(fmt(abs(cc[w]))), CPU-GPU $(fmt(dc[w]))); CPU-GPU median $(fmt(sort(dc)[no ÷ 2]))",
            "$(fmt(maximum(dp))); $(fmt(sort(dp)[no ÷ 2])); $(fmt(dp[w]))", "-", "INFO")
    end
end

# AC12 (lyap-gram-schmidt): RGS and MGS, S_i / n for i = 1, 2 (the lambda of QRSpectrum).
MODTAG[] = "lyap-gram-schmidt: "
want("gs") && for (nm, alg) in (("RGS", RGS()), ("MGS", MGS()))
    est = qr_run(QRSpectrum(alg))
    ac12_regular("AC12", "QRSpectrum($nm)", est, lam_of, 2)
    for T in (Float64, Double64)
        ac12_chaotic("AC12", "QRSpectrum($nm)", est, lam_of, 2, T)
    end
end

# AC14 (lyap-householder-qr): (a) kernel path HouseholderQR, GivensQR, AugmentedQR, k = 2, plus |lambda_1 + lambda_2| <= theta_sum
# (GS G13, `hh_theta`: 16 eps (G_max / tau + M_max), M_max the largest tangent entry of the orbit) on both backends; (b) library
# path (LapackQR) with a gpu or auto request resolves to cpu, reason "reducer is CPU-only"; (c) K15 performance record (hh-perf).
MODTAG[] = "lyap-householder-qr: "
(want("hh") || want("hh-perf")) && include(joinpath(ROOT, "test", "lyap-householder-qr", "hh_helpers.jl"))   # hh_theta, hh_Mmax
function hh_sum_check(lbl, sys, x0, N, tau, cpu, gpu)
    T = eltype(eltype(x0))
    th = [hh_theta(T, f[1].G_max, tau, hh_Mmax(sys, z, N)) for (f, z) in zip(cpu.finals, x0)]
    r(res) = maximum(Float64(abs(f[1].lambda[1] + f[1].lambda[2])) / t for (f, t) in zip(res.finals, th))
    rc, rg = r(cpu), r(gpu)
    row("AC14", "$lbl, $T: max over orbits of |lambda_1 + lambda_2| / theta_sum, CPU; GPU", "$(fmt(rc)); $(fmt(rg))", "<= 1", pf(rc <= 1 && rg <= 1))
end
want("hh") && let
    for (nm, alg) in (("HouseholderQR", HouseholderQR()), ("GivensQR", GivensQR()), ("AugmentedQR", AugmentedQR()))
        est = qr_run(QRSpectrum(alg))
        ac12_regular("AC14", "QRSpectrum($nm)", est, lam_of, 2; extra = (a...) -> hh_sum_check("(a) $nm, K=0.5", a...))
        for T in (Float64, Double64)
            out = ac12_chaotic("AC14", "QRSpectrum($nm)", est, lam_of, 2, T)
            out === nothing || hh_sum_check("(c) $nm, K=10", out[1], out[2], out[3], 2, out[4], out[5])
        end
    end
    x0 = [SVector(0.1, 0.2), SVector(0.3, 0.4)]
    for b in (:gpu, :auto)
        p = propagate(StandardMap(2pi), x0, 4; k = 2, tau = 2, seed = SEED, reducers = (QRSpectrum(LapackQR()),), backend = b).provenance
        ok = p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == "reducer is CPU-only"
        row("AC14", "(b) library path QRSpectrum(LapackQR), $b request on a machine with a GPU: backend; device; reason", "$(p.backend); $(p.device); \"$(p.backend_fallback_reason)\"", "cpu; cpu; \"reducer is CPU-only\"", pf(ok))
    end
end

# AC14(c) (lyap-householder-qr, K15): kernel path chain QRSpectrum(HouseholderQR()), 10^6 orbits, d = k = 2, N = 10^4, K = 2 pi,
# tau = 2: wall time on the GPU, on the CPU with all threads, and on the CPU with 1 thread (child process). The gate (auto may pick
# the GPU only if GPU time <= CPU time) is reported, not changed here: the record file stays the AC13 one.
MODTAG[] = "lyap-householder-qr: "
want("hh-perf") && let
    no = 10^6
    tg, pg = hh_perf_run(:gpu, no)
    tc, pc = hh_perf_run(:cpu, no)
    f = tempname()
    Base.run(`$(Base.julia_cmd()) -t 1 --project=$(dirname(Base.active_project())) $(@__FILE__) child-hh $f`)
    t1 = deserialize(f)
    ok = pg.backend === :gpu && pc.backend === :cpu
    row("AC14", "(c) K15 record, QRSpectrum(HouseholderQR), $no orbits, N=1e4, tau=2, K=2pi: wall s GPU ($DEVNAME, CUDA runtime $(CUDA.runtime_version()), $(pg.gpu_chunks) chunks); CPU $(Threads.nthreads()) threads; CPU 1 thread ($(Sys.CPU_NAME)); gate GPU <= CPU $(Threads.nthreads()) threads",
        "$(fmt(tg)); $(fmt(tc)); $(fmt(t1)); $(tg <= tc)", "recorded", ok ? "INFO" : fb(something(pg.backend_fallback_reason, "?")))
end

# AC17 (lyap-two-trajectory): X2 of propagate_shadow, d0 = 1e-8 (as AC11). (a) regular class in L2 round-off units (spec Amendment
# D42), not the tangent-core 1e-6 count: same S7 starts, n = 4e4, tau = 1, lifted mode; per orbit |X2_cpu - X2_gpu| <= 4 E_ro(d0)
# Psi_n / n, the round-off part of T5a (2 E_ro Psi_n / n) once per backend, E_ro with the orbit's measured M_state; orbits with
# E_ro > 0.5 (outside the T5 validity) are reported, not gated. (c) chaotic ensemble: the common protocol.
MODTAG[] = "lyap-two-trajectory: "
want("l2") && include(joinpath(ROOT, "test", "lyap-two-trajectory", "l2_helpers.jl"))   # l2_Ero, l2_Mstate, l2_Psis
want("l2") && let
    est = (sys, x0, N, tau, b) -> propagate_shadow(sys, x0, N; d0 = eltype(eltype(x0))(1e-8), tau = tau, seed = SEED, backend = b)
    lam = f -> SVector(f.X2)
    x0, ncand, nrep = regular_starts(1024)
    sys, N, d0, K = StandardMap(0.5), 4 * 10^4, 1e-8, 0.5
    gpu = est(sys, x0, N, 1, :gpu)
    why = gpu_why(gpu)
    if why !== nothing
        row("AC17", "(a) K=0.5 regular class, n=$(fmt(N)), tau=1: X2 in L2 round-off units (D42)", "-", "<= 1", fb(why))
    else
        cpu = est(sys, x0, N, 1, :cpu)
        ero = [l2_Ero(K, d0, l2_Mstate(sys, z, N), eps(Float64)) for z in x0]
        psi = l2_Psis(sys, x0, SEED, Float64, N)
        ratio = [abs(c[1].X2 - g[1].X2) / (4 * e * p / N) for (c, g, e, p) in zip(cpu.finals, gpu.finals, ero, psi)]
        ok = ero .<= 0.5
        row("AC17", "(a) K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced), lifted, n=$(fmt(N)), tau=1, d0=1e-8: max over the $(count(ok)) valid orbits (E_ro <= 0.5) of |X2_cpu - X2_gpu| / (4 E_ro Psi_n / n) (D42); outside validity $(count(.!ok)) orbits, their max ratio $(any(.!ok) ? fmt(maximum(ratio[.!ok])) : "-"); max |X2_cpu - X2_gpu| $(fmt(maximum(abs(c[1].X2 - g[1].X2) for (c, g) in zip(cpu.finals, gpu.finals)))); E_ro max $(fmt(maximum(ero))); Psi_n/n median $(fmt(sort(psi)[length(psi) ÷ 2] / N))",
            any(ok) ? maximum(ratio[ok]) : NaN, "<= 1", pf(any(ok) && maximum(ratio[ok]) <= 1))
    end
    for T in (Float64, Double64)
        ac12_chaotic("AC17", "X2", est, lam, 1, T)
    end
end

# ------------------------------------------------------------------------------------- D44 (tangent-core Amendment D44, local)
# Step-chunked relaunch on the local device (GTX 1050 at D44, RTX 5060 from D46), after the record (the budget reads its one-wave value). Tag "d44".
MODTAG[] = "tangent-core: "
const Tg = Chaoslab.Tangent
function with_L(f, L)
    old = Tg._LAUNCH_STEPS[]
    Tg._LAUNCH_STEPS[] = L
    try
        return f()
    finally
        Tg._LAUNCH_STEPS[] = old
    end
end
same_run(a, b) = isequal(a.x, b.x) && isequal(a.W, b.W) && isequal(a.finals, b.finals)
on_dev(p, L, N) = p.backend === :gpu && p.backend_fallback_reason == "" && p.gpu_chunks > 0 && p.launch_steps == L && p.launches_per_chunk == cld(N, L)

# D44c [local]: N = 2e3 (one launch fits); forced L < N equals forced L = N bitwise on the device. Orbits whose single-launch run
# raises a health error are dropped (as the package test D44a; RGS breaks down on a few K = 2 pi, tau = 10 orbits): first on the CPU,
# then on the device (a device-only breakdown would send the whole run to the CPU with reason 5, making the comparison CPU vs CPU).
(want("d44") || want("d44c")) && let
    N = 2000
    rng = StableRNG(SEED + 4403)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:1024]
    for (name, k, red) in (("Benettin", 1, (Benettin(),)), ("QRSpectrum(RGS)", 2, (QRSpectrum(RGS()),)), ("SVDSpectrum", 2, (SVDSpectrum(),))), tau in (2, 10)
        Kc = tau == 10 ? 0.5 : 2pi   # D46: K = 2 pi with tau = 10 is outside the RGS window (58 of 1024 orbits survived)
        sys = StandardMap(Kc)
        w(j) = initial_deviations(SEED, j, Val(2), Val(k), Float64)
        keep = [j for j in eachindex(x0) if catchhe(() -> propagate_orbit(sys, x0[j], w(j), N, tau, red; orbit = j)) === nothing]
        ncpu = length(keep)
        filter!(j -> propagate(sys, [x0[j]], N; k = k, tau = tau, seed = SEED, reducers = red, W0 = [w(j)], backend = :gpu).provenance.backend === :gpu, keep)
        xs, W0 = x0[keep], [w(j) for j in keep]
        println("# D44c $name tau=$tau: kept $(length(keep)) of 1024 orbits ($(1024 - ncpu) fail on the CPU, $(ncpu - length(keep)) more only on the device)")
        run() = propagate(sys, xs, N; k = k, tau = tau, seed = SEED, reducers = red, W0 = W0, backend = :gpu)
        full = with_L(run, N)
        ok = on_dev(full.provenance, N, N) && same_run(run(), full)
        for L in (tau, 20tau, 130tau)
            r = with_L(run, L)
            ok &= on_dev(r.provenance, L, N) && same_run(r, full)
        end
        row("D44c", "device bitwise, $name, K=$(tau == 10 ? "0.5" : "2pi"), $(length(xs)) orbits, N=2e3, tau=$tau: forced L in (tau, 20 tau, 130 tau) equals L = N (and the unforced run), all on the device",
            ok ? "equal" : "differs", "bitwise", pf(ok))
    end
end

# D44k(c) [local]: the same for propagate_shadow, lifted and wrapped, K = 0.5 and 2 pi, d0 = 1e-8, 256 orbits.
want("d44") && let
    N = 2000
    rng = StableRNG(SEED + 4411)
    z0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:256]
    for mode in (:lifted, :wrapped), K in (0.5, 2pi), tau in (2, 10)
        sys = mode === :wrapped ? WrappedStandardMap(K) : StandardMap(K)
        x0 = mode === :wrapped ? [torus_state(z) for z in z0] : z0
        run() = propagate_shadow(sys, x0, N; d0 = 1e-8, tau = tau, seed = SEED, backend = :gpu)
        full = with_L(run, N)
        ok = on_dev(full.provenance, N, N) && same_run(run(), full)
        for L in (tau, 20tau, 130tau)
            r = with_L(run, L)
            ok &= on_dev(r.provenance, L, N) && same_run(r, full)
        end
        row("D44k", "(c) device bitwise, propagate_shadow, $mode, K=$(K == 2pi ? "2pi" : K), 256 orbits, N=2e3, tau=$tau: forced L in (tau, 20 tau, 130 tau) equals L = N, all on the device",
            ok ? "equal" : "differs", "bitwise", pf(ok))
    end
end

# D44h [local part]: every lost-device code of the extension's table yields "device lost" and never the budget string; codes that
# leave the context usable (out of memory, invalid value) are not lost-device errors.
want("d44") && let
    ext = Base.get_extension(Chaoslab, :ChaoslabCUDAExt)
    lost = ext.LOST_CODES
    got = [FALLBACK_REASONS[Tg._fallback_reason(CUDA.CuError(c))] for c in lost]
    other = [CUDA.ERROR_OUT_OF_MEMORY, CUDA.ERROR_INVALID_VALUE]
    ok = all(==("device lost"), got) && CUDA.ERROR_LAUNCH_TIMEOUT in lost && CUDA.ERROR_ILLEGAL_ADDRESS in lost &&
         all(c -> !F.is_device_lost(CUDA.CuError(c)), other)
    row("D44h", "lost-device codes $(join(Int.(lost), ", ")): reason of each; codes $(join(Int.(other), ", ")) not lost",
        join(unique(got), "; "), "\"device lost\" for all, 700 and 702 included", pf(ok))
end

# D44i: one stepping body (inspection; the counts of the tangent-core AC14 check).
want("d44") && let
    tdir = joinpath(ROOT, "src", "Tangent")
    code(l) = first(split(l, '#'))
    hits(pat) = [(basename(f), i) for f in filter(f -> endswith(f, ".jl"), readdir(tdir; join = true)) for (i, l) in enumerate(eachline(f)) if occursin(pat, code(l))]
    fwd, tan = hits("forward(sys"), hits("tangent(sys")
    calls = filter(h -> !occursin("@inline function _orbit_body", code(readlines(joinpath(tdir, h[1]))[h[2]])), hits("_orbit_body("))
    ok = length(fwd) == 1 && length(tan) == 1 && length(calls) == 2
    row("D44i", "single stepping body: forward(sys / tangent(sys sites in src/Tangent; callers of _orbit_body (kernel, propagate_orbit), $(join(("$(c[1]):$(c[2])" for c in calls), ", "))",
        "$(length(fwd)) / $(length(tan)); $(length(calls))", "1 / 1; 2", pf(ok))
end

# D44j [local]: the rows cut down or fallen back for the launch budget, rerun at n = 1e5 on the device, original tolerances. Per row:
# provenance on the device with at least one launch per chunk (D46: one launch when the device's budget fits N, as on the RTX 5060;
# the relaunch path is then covered by D44c and D44k), the per-launch estimate from provenance (<= 1 s), and the measured
# time of one full launch (the same row at N = L with one chunk of the run's device chunk size, <= 4 s).
function d44j(ac, name, est, lam, nl, x0; extra = nothing, N = 10^5, tau = 10)
    sys = StandardMap(0.5)
    gpu = est(sys, x0, N, tau, :gpu)
    p = gpu.provenance
    launched = p.backend === :gpu && p.backend_fallback_reason == "" && p.gpu_chunks > 0 && p.launches_per_chunk >= 1
    launched || return row(ac, "D44j $name, n=$(fmt(N)): on the device in launches", "backend $(p.backend), reason \"$(p.backend_fallback_reason)\", launches $(p.launches_per_chunk)", "gpu, at least 1 launch", "FAIL")
    ns = something(F.gpu_step_ns(DEV, Float64), Tg.DEFAULT_GPU_NS)
    estm = ns * 1e-9 * p.launch_steps * max(p.chunk_size, 1024)
    xs = x0[1:min(length(x0), p.gpu_chunk_size)]
    t1 = @elapsed g1 = with_L(() -> est(sys, xs, p.launch_steps, tau, :gpu), p.launch_steps)
    one = g1.provenance.backend === :gpu && g1.provenance.gpu_chunks == 1 && g1.provenance.launches_per_chunk == 1
    cpu = est(sys, x0, N, tau, :cpu)
    ref = est(StandardMap(Double64(0.5)), [SVector{2,Double64}(z) for z in x0], N, tau, :cpu)
    for i in 1:nl
        lc, lg, lr = [Float64(l[i]) for l in lams(cpu, lam)], [Float64(l[i]) for l in lams(gpu, lam)], [Float64(l[i]) for l in lams(ref, lam)]
        dcr, dgr, dcg = abs.(lc .- lr), abs.(lg .- lr), abs.(lc .- lg)
        nc, ng = count(>(1e-6), dcr), count(>(1e-6), dgr)
        both = (dcr .<= 1e-6) .& (dgr .<= 1e-6)
        gap = any(both) ? maximum(dcg[both]) : Inf
        ok = abs(nc - ng) <= 3 && gap <= 2e-6 && estm <= 1 && one && t1 <= 4
        row(ac, "D44j $name lambda_$i, K=0.5, $(length(x0)) orbits, n=$(fmt(N)), tau=$tau, L=$(p.launch_steps), $(p.launches_per_chunk) launches per chunk, $(p.gpu_chunks) device chunk(s) of $(p.gpu_chunk_size): orbits above 1e-6 vs the Double64 reference CPU / GPU; gap over the $(count(both)) agreeing orbits; per-launch estimate s; measured full launch s ($(length(xs)) orbits x L, one launch: $one); ratio measured/estimate (max |cpu-gpu| $(fmt(maximum(dcg))), max |gpu-ref| $(fmt(maximum(dgr))))",
            "$nc / $ng (diff $(abs(nc - ng))); $(fmt(gap)); $(fmt(estm)); $(fmt(t1)); $(fmt(t1 / estm))", "count diff <= 3; gap <= 2e-6; <= 1; <= 4", pf(ok))
    end
    extra === nothing || extra(sys, x0, N, tau, cpu, gpu)
end
if want("d44")
    let rng = StableRNG(SEED + 1201), N = 10^5   # tangent-core AC12a starts and reducer
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:1024]
        est = (sys, x, n, tau, b) -> propagate(sys, x, n; k = 1, tau = tau, seed = SEED, reducers = (LogNorm(),), backend = b)
        d44j("AC12a", "tangent-core LogNorm (S/n)", est, f -> f.logs ./ N, 1, x0)
    end
    let rng = StableRNG(SEED + 41000 + 50 * 16 + 1)   # lyap-benettin AC16 starts (vs(16, 1)) and seed (vs(16))
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:1024]
        est = (sys, x, n, tau, b) -> propagate(sys, x, n; k = 1, tau = tau, seed = SEED + 41800, reducers = (Benettin(),), backend = b)
        MODTAG[] = "lyap-benettin: "
        d44j("AC16", "Benettin", est, lam_of, 1, x0)
    end
    xr = first(regular_starts(1024))
    MODTAG[] = "lyap-gram-schmidt: "
    for (nm, alg) in (("RGS", RGS()), ("MGS", MGS()))
        d44j("AC12", "QRSpectrum($nm)", qr_run(QRSpectrum(alg)), lam_of, 2, xr)
    end
    MODTAG[] = "lyap-householder-qr: "
    isdefined(Main, :hh_theta) || include(joinpath(ROOT, "test", "lyap-householder-qr", "hh_helpers.jl"))
    for (nm, alg) in (("HouseholderQR", HouseholderQR()), ("GivensQR", GivensQR()), ("AugmentedQR", AugmentedQR()))
        d44j("AC14", "QRSpectrum($nm)", qr_run(QRSpectrum(alg)), lam_of, 2, xr; extra = (a...) -> hh_sum_check("(a) D44j $nm, K=0.5, n=1e5", a...))
    end
    MODTAG[] = "lyap-svd-stabilized: "
    d44j("AC10", "SVDSpectrum", qr_run(SVDSpectrum()), lam_of, 2, xr)
    MODTAG[] = "tangent-core: "
end

# ------------------------------------------------------------------------------------------------ AC12a (tangent-core), last
# D33: symmetric form. Each backend is compared with a Double64 CPU reference on all 1024 seeded orbits; the numbers of orbits above
# 1e-6 may differ by at most 3, and on the orbits where both backends agree with the reference the CPU-GPU gap is within 2e-6.
# Before D44 the spec's n = 1e5 (one launch of about 2 s) fell back to the CPU for the 1 s launch budget and the criterion was
# evaluated at n = 4e4; since D44 n = 1e5 runs on the device in launches (the D44j rows above add the launch checks).
MODTAG[] = "tangent-core: "
function ac12a(N)
    K, no, tau = 0.5, 1024, 10
    rng = StableRNG(SEED + 1201)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:no]
    sys = StandardMap(K)
    kw = (k = 1, tau = tau, seed = SEED, reducers = (LogNorm(),))
    tg = @elapsed ((gpu, why) = on_gpu(sys, x0, N; kw...))
    why === nothing || return row("AC12a", "K=0.5, $no orbits, n=$(fmt(N)), tau=$tau (D33 symmetric form)", "-", "count diff <= 3; gap <= 2e-6", fb(why))
    cpu = propagate(sys, x0, N; kw..., backend = :cpu)
    ref = propagate(StandardMap(Double64(K)), [SVector{2,Double64}(z) for z in x0], N;
                    k = 1, tau = tau, seed = SEED, reducers = (LogNorm(),), backend = :cpu)
    sc, sg, sr = slogs(cpu, N), slogs(gpu, N), Float64.(slogs(ref, N))
    dcr, dgr, dcg = abs.(sc .- sr), abs.(sg .- sr), abs.(sc .- sg)
    ncpu, ngpu = count(>(1e-6), dcr), count(>(1e-6), dgr)
    both = (dcr .<= 1e-6) .& (dgr .<= 1e-6)
    gap = maximum(dcg[both])
    ok = abs(ncpu - ngpu) <= 3 && gap <= 2e-6
    row("AC12a", "K=0.5, $no orbits, n=$(fmt(N)), tau=$tau, GPU run $(fmt(tg)) s ($(gpu.provenance.launches_per_chunk) launch(es) of $(gpu.provenance.launch_steps) steps): orbits above 1e-6 against the Double64 reference, CPU / GPU; max CPU-GPU gap over the $(count(both)) orbits where both agree with the reference",
        "$ncpu / $ngpu (diff $(abs(ncpu - ngpu))); $(fmt(gap))", "count diff <= 3; gap <= 2e-6", pf(ok))
    row("AC12a", "same, maxima over all $no orbits (chaotic layer included): |S_cpu - S_ref| $(fmt(maximum(dcr))), |S_gpu - S_ref| $(fmt(maximum(dgr))), |S_cpu - S_gpu| $(fmt(maximum(dcg)))",
        maximum(dcg), "-", "INFO")
end
if want("tc")
    ac12a(10^5)
    ac12a(4 * 10^4)
end

# ------------------------------------------------------------------------------------------------ AC18 (sali), tag "sali"
# Spec criterion 18 as amended by D68 (the same-orbit bound applies per backend, not across backends: the GPU and the CPU do not produce the
# same orbit). SALI reducer (k = 2), lifted, tau = 1; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64, where the kernel
# does not compile). Classes: regular K = 0.5 (S7 screening, regular_starts), 1024 orbits, n = 4e4; chaotic K = 10, 1024 orbits, n = 8 (inside
# the resolved window; direct SALI is gated only where SALI_ref >= 10 B_meas). (a) each backend against a BigFloat(256) reference formed at its
# own realized orbit points (L1 S3): max |L_b - L_ref,b| / (8 eps (Psi1 + Psi2)) <= 1 and max |SALI_b - SALI_ref,b| / B_SALI <= 1. (b) CPU-GPU
# gap only on the orbits whose final states agree to 1e-9 in every coordinate: max |L_cpu - L_gpu| / (16 eps (Psi1 + Psi2)) <= 1 and
# max |SALI_cpu - SALI_gpu| / (2 B_SALI) <= 1; their number is reported. (c) K = 10, n = 1e4, tau = 2: tangent-core AC12(c) tolerances on -L/(2n).
# Reported, not gated: orbits with final states apart by more than 1e-9 (CPU vs GPU, CPU Float64 vs Double64 CPU).
# Chaotic class of (b), owner addition 2026-10-22 (sali-spec.md, Amendment D68, "Addition to part (b) for the chaotic class"): the orbits do not
# agree there, so the gap gets the orbit-divergence term D_orb(n) (d68_orb below) and the 1e-9 final-state selection is replaced by the item 4
# horizon 16 eps M_n Gamma_n <= 1e-4 (precision-spec Equations 3, 4). Gates: |L_cpu - L_gpu| <= 16 eps (Psi1 + Psi2) + D_orb and, in the resolved
# window, |SALI_cpu - SALI_gpu| <= 2 B_SALI + SALI_ref D_orb; the item 4 bound itself is checked on the device orbit (a FAIL there voids the term).
# All 1024 orbits get a reference for both backends (no subsample). The GPU's own orbit is read with a trace kernel that runs the same
# `forward` on the device (the run itself keeps no orbit); its final state is checked bitwise against the propagate result.
MODTAG[] = "sali: "
want("sali") && include(joinpath(ROOT, "test", "sali", "sl_helpers.jl"))   # sl_*: B_SALI, B_L, BigFloat reference, Psi_n
sali_est(sys, x0, N, tau, b) = propagate(sys, x0, N; k = 2, tau = tau, seed = SEED, reducers = (SALI(),), backend = b)

# Device orbit trace: D[j, i] = point of orbit j before step i (the point where the step-i tangent is evaluated), X[j] advanced.
function gpu_trace_kernel!(D, X, sys, steps)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if j <= length(X)
        x = X[j]
        for i in 1:steps
            D[j, i] = x
            x = forward(sys, x)
        end
        X[j] = x
    end
    return nothing
end
function gpu_trace(sys, x0::Vector{SVector{2,Float64}}, N; chunk = 5000)
    m = length(x0)
    X = CuArray(x0)
    P = Matrix{SVector{2,Float64}}(undef, m, N)   # P[j, i]: orbit j before step i, on the host
    for s0 in 0:chunk:N-1
        s1 = min(s0 + chunk, N)
        D = CuArray{SVector{2,Float64}}(undef, m, s1 - s0)
        @cuda threads = 256 blocks = cld(m, 256) gpu_trace_kernel!(D, X, sys, s1 - s0)
        P[:, s0+1:s1] = Array(D)
        CUDA.unsafe_free!(D)
    end
    return P, Array(X)
end

# BigFloat(256) reference of orbit j at the points `pts` of backend result `res`: ratios to the single-backend bounds and the scales.
function d68_one(res, pts, Kt, sysT, j, N, T, epsT)
    W0 = sl_frame(SEED, j, T)
    ref = sl_ref(sl_mats_std(Kt, pts), W0, N; keep = [N])
    Ms = [SMatrix{2,2,Float64}(map(Float64, tangent(sysT, pt))) for pt in pts]
    p1, p2 = sl_psis(Ms, W0, [N])[1]
    at = sl_at(ref, N)
    sref = Float64(at.S)
    BS = sl_BSALI(epsT, sl_Lambda(ref.lnsali, N), sref, p1, p2)
    f = res.finals[j][1]
    return (rL = sl_err(f.L, at.L) / sl_BL(epsT, p1, p2), rS = sl_err(f.sali, at.S) / BS, resolved = sref >= 10 * 4 * epsT,
            BS = BS, bnd16 = 16 * epsT * (p1 + p2), Lref = Float64(at.L), Sref = sref)
end

# Largest singular value of a 2 x 2 matrix: sqrt of the larger root of s^2 - tr(A'A) s + det(A)^2.
function d68_sig2(A)
    t, d = sum(abs2, A), det(A)^2
    return sqrt((t + sqrt(max(t^2 - 4d, zero(t)))) / 2)
end

# D_orb(n) of sali-spec.md, Amendment D68, "Addition to part (b) for the chaotic class" (the one place the formula lives), in BigFloat(256) on the
# CPU reference orbit: pts[i] = z_{i-1} (i = 1..n), zN = z_n, W0 the seeded frame (w_1, w_2). With A_k = DT(z_k), P_{<k} = A_{k-1} ... A_0,
# P_{>k} = A_{n-1} ... A_{k+1}, bnd_k = 16 eps M_k Gamma_k (precision-spec Equations 3, 4: M_k = max(1, max_{j<=k} |x_j|, |p_j|),
# Gamma_k = sum_{j<k} norm(A_{k-1} ... A_{j+1})_2, empty product = identity, Gamma_0 = 0):
#   D_orb = sqrt 2 |K| sum_{i=1,2} sum_{k<n} norm(P_{>k})_2 bnd_k norm(P_{<k} w_i) / norm(P w_i).
# Returns D_orb, the item 4 bounds bnd_k for k = 0..n (the horizon test is bnd_n <= 1e-4), all as Float64.
function d68_orb(Kt, pts, zN, W0, epsT)
    return setprecision(BigFloat, 256) do
        n = length(pts)
        mats = sl_mats_std(Kt, pts)
        A = [(m = mats(k + 1); SMatrix{2,2,BigFloat}(m[1], m[3], m[2], m[4])) for k in 0:n-1]   # A[k+1] = DT(z_k)
        zs = vcat(pts, [zN])                                                                    # zs[k+1] = z_k
        Mk = accumulate(max, [max(one(BigFloat), abs(BigFloat(z[1])), abs(BigFloat(z[2]))) for z in zs])
        Gam = map(0:n) do k
            Q, g = SMatrix{2,2,BigFloat}(1, 0, 0, 1), zero(BigFloat)
            for j in k-1:-1:0            # Q = A_{k-1} ... A_{j+1}
                g += d68_sig2(Q)
                Q = Q * A[j+1]
            end
            g
        end
        bnd = [16 * BigFloat(epsT) * Mk[k+1] * Gam[k+1] for k in 0:n]
        Pgt = Vector{SMatrix{2,2,BigFloat,4}}(undef, n)   # Pgt[k+1] = P_{>k}
        Pgt[n] = SMatrix{2,2,BigFloat}(1, 0, 0, 1)
        for k in n-1:-1:1
            Pgt[k] = Pgt[k+1] * A[k+1]
        end
        D = zero(BigFloat)
        for i in 1:2
            v = [SVector{2,BigFloat}(W0[1, i], W0[2, i])]
            for k in 0:n-1
                push!(v, A[k+1] * v[end])
            end
            nv = norm.(v)
            D += sum(d68_sig2(Pgt[k+1]) * bnd[k+1] * nv[k+1] / nv[n+1] for k in 0:n-1)
        end
        (D = Float64(sqrt(BigFloat(2)) * abs(BigFloat(Kt)) * D), bnd = Float64.(bnd))
    end
end

# Chaotic class of part (b) (Amendment D68, addition of 2026-10-22). Entering orbits: 16 eps M_n Gamma_n <= 1e-4 (item 4 horizon). Gates:
# |L_cpu - L_gpu| <= 16 eps (Psi1 + Psi2) + D_orb and |SALI_cpu - SALI_gpu| <= 2 B_SALI + SALI_ref D_orb (resolved orbits); item 4 check on the
# device: max_k |e_k|_inf / (16 eps M_k Gamma_k) <= 1 (what D_orb uses) and max_k |e_k|_inf / (16 eps M_n Gamma_n) <= 1 (the literal n bound).
function d68_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    enter = findall(j -> orb[j].hor <= 1e-4, 1:m)
    rest = setdiff(1:m, enter)
    row("AC18b", "$lab: orbits that enter part (b): 16 eps M_n Gamma_n <= 1e-4 (horizon max over all orbits $(fmt(maximum(o.hor for o in orb))); reported only: $(length(rest)))", length(enter), "reported", "INFO")
    gap(j) = abs(Float64(cpu.finals[j][1].L) - Float64(gpu.finals[j][1].L))
    gapS(j) = abs(Float64(cpu.finals[j][1].sali) - Float64(gpu.finals[j][1].sali))
    rL(j) = gap(j) / (rc[j].bnd16 + orb[j].D)
    rS(j) = gapS(j) / (2 * rc[j].BS + rc[j].Sref * orb[j].D)
    res(js) = [j for j in js if rc[j].resolved && rg[j].resolved]
    entS = res(enter)
    gL, gS = [rL(j) for j in enter], [rS(j) for j in entS]
    row("AC18b", "$lab: entering orbits: max |L_cpu - L_gpu| / (16 eps (Psi1 + Psi2) + D_orb(n)) over $(length(enter)) orbits", isempty(gL) ? NaN : maximum(gL), "<= 1", pf(!isempty(gL) && maximum(gL) <= 1))
    row("AC18b", "$lab: entering orbits: max |SALI_cpu - SALI_gpu| / (2 B_SALI(n) + SALI_ref(n) D_orb(n)) over the $(length(entS)) resolved ones", isempty(gS) ? NaN : maximum(gS), "<= 1", pf(!isempty(gS) && maximum(gS) <= 1))
    med(v) = isempty(v) ? NaN : sort(v)[max(end ÷ 2, 1)]
    q = [gap(j) / orb[j].D for j in enter]
    nz = filter(>(0), q)
    row("AC18b", "$lab: entering orbits: gap / D_orb(n) = |L_cpu - L_gpu| / D_orb(n), median / max (D_orb median $(fmt(med([orb[j].D for j in enter]))), max $(fmt(maximum(orb[j].D for j in enter)))); nonzero gaps: $(length(nz)) of $(length(q)), median among them $(fmt(med(nz)))",
        "$(fmt(med(q))) / $(fmt(maximum(q)))", "reported (how loose the term is)", "INFO")
    row("AC18b", "$lab: entering orbits: max |L_cpu - L_gpu| / (16 eps (Psi1 + Psi2)) without the D_orb term (what the first D68 run gated), $(count(j -> gap(j) > rc[j].bnd16, enter)) of $(length(enter)) above 1",
        maximum(gap(j) / rc[j].bnd16 for j in enter), "reported", "INFO")
    pk, lt = maximum(orb[j].perk for j in enter), maximum(orb[j].lit for j in enter)
    row("AC18b", "$lab: item 4 check on the device, entering orbits: max over orbits and k = 1..n of |z_k^cpu - z_k^gpu|_inf / (16 eps M_k Gamma_k); same against the literal 16 eps M_n Gamma_n (orbits above 1: $(count(j -> orb[j].perk > 1, enter)) / $(count(j -> orb[j].lit > 1, enter)))",
        "$(fmt(pk)); $(fmt(lt))", "<= 1; <= 1", pf(pk <= 1 && lt <= 1))
    if !isempty(rest)
        row("AC18b", "$lab: reported-only orbits ($(length(rest))): max gap ratios (L; SALI over $(length(res(rest))) resolved); item 4 per-k ratio max", "$(fmt(maximum(rL, rest))); $(isempty(res(rest)) ? NaN : fmt(maximum(rS, res(rest)))); $(fmt(maximum(orb[j].perk for j in rest)))", "reported", "INFO")
    end
end

function sali_d68(tag, T, Kt, xT, N; chaotic = false)
    sysT = StandardMap(Kt)
    m = length(xT)
    lab = "$tag, $T, $m orbits, n=$(fmt(N)), tau=1"
    gpu = sali_est(sysT, xT, N, 1, :gpu)
    why = gpu_why(gpu)
    why === nothing || return row("AC18a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    cpu = sali_est(sysT, xT, N, 1, :cpu)
    epsT = Float64(effective_eps(T))
    P, XN = gpu_trace(sysT, xT, N)
    rc = Vector{NamedTuple}(undef, m); rg = Vector{NamedTuple}(undef, m)
    same_c = zeros(Bool, m)
    stepdiff = zeros(m)   # largest CPU-GPU difference of the orbit points over the steps (diagnostic)
    orb = Vector{Any}(undef, m)   # chaotic class: d68_orb result and the realized orbit differences e_k, k = 1..N (inf norm)
    Threads.@threads :dynamic for j in 1:m
        ptsc, zc = lb_pts(sysT, xT[j], N)
        same_c[j] = zc == cpu.x[j]
        stepdiff[j] = maximum(maximum(abs.(ptsc[i] - P[j, i])) for i in 1:N)
        rc[j] = d68_one(cpu, ptsc, Kt, sysT, j, N, T, epsT)
        rg[j] = d68_one(gpu, P[j, :], Kt, sysT, j, N, T, epsT)
        if chaotic
            od = d68_orb(Kt, ptsc, zc, sl_frame(SEED, j, T), epsT)
            e = [i < N ? maximum(abs.(ptsc[i+1] - P[j, i+1])) : maximum(abs.(zc - XN[j])) for i in 1:N]   # e_k = |z_k^cpu - z_k^gpu|_inf
            orb[j] = (D = od.D, hor = od.bnd[end], perk = maximum(e ./ od.bnd[2:end]), lit = maximum(e) / od.bnd[end])
        end
    end
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    row("AC18a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        rL = [r.rL for r in R]
        rS = [r.rS for r in R if r.resolved]
        nres = length(rS)
        notgated = count(i -> R[i].bnd16 / 2 > abs(R[i].Lref), 1:m)
        row("AC18a", "$lab: backend $nm (provenance $(res.provenance.backend)): max |L - L_ref| / (8 eps (Psi1 + Psi2)) over $m orbits (median $(fmt(sort(rL)[end ÷ 2])); orbits with the bound above |L|: $notgated, reported)",
            maximum(rL), "<= 1", pf(maximum(rL) <= 1))
        row("AC18a", "$lab: backend $nm: max |SALI - SALI_ref| / B_SALI(n) over the $nres of $m resolved orbits (median $(fmt(sort(rS)[max(end ÷ 2, 1)])))",
            maximum(rS), "<= 1", pf(nres > 0 && maximum(rS) <= 1))
    end
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = findall(<=(1e-9), dx)
    row("AC18b", "$lab: orbits with final states apart by more than 1e-9, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(sort(dx)[end ÷ 2])), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    Kd = exact(Double64, Kt)
    ref64 = sali_est(StandardMap(Kd), [SVector{2,Double64}(z) for z in xT], N, 1, :cpu)
    dd = [Float64(maximum(abs.(SVector{2,Double64}(cpu.x[j]) - ref64.x[j]))) for j in 1:m]
    row("AC18b", "$lab: orbits with final states apart by more than 1e-9, CPU Float64 vs Double64 CPU (max $(fmt(maximum(dd))))", count(>(1e-9), dd), "reported", "INFO")
    if chaotic
        row("AC18b", "$lab: orbits with final states within 1e-9 in every coordinate (the selection of the first D68 run; no longer the entering rule for this class)", length(agree), "reported", "INFO")
        return d68_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    end
    row("AC18b", "$lab: orbits that enter part (b) (final states within 1e-9 in every coordinate)", length(agree), "reported", "INFO")
    gL = [abs(Float64(cpu.finals[j][1].L) - Float64(gpu.finals[j][1].L)) / rc[j].bnd16 for j in agree]
    gS = [abs(Float64(cpu.finals[j][1].sali) - Float64(gpu.finals[j][1].sali)) / (2 * rc[j].BS) for j in agree if rc[j].resolved && rg[j].resolved]
    row("AC18b", "$lab: orbits that agree: max |L_cpu - L_gpu| / (16 eps (Psi1 + Psi2)) over $(length(agree)) orbits", isempty(gL) ? NaN : maximum(gL), "<= 1", pf(!isempty(gL) && maximum(gL) <= 1))
    row("AC18b", "$lab: orbits that agree: max |SALI_cpu - SALI_gpu| / (2 B_SALI(n)) over the $(length(gS)) resolved ones", isempty(gS) ? NaN : maximum(gS), "<= 1", pf(!isempty(gS) && maximum(gS) <= 1))
    sd = stepdiff[agree]
    row("AC18b", "$lab: diagnostic on the orbits that agree: largest CPU-GPU orbit-point difference over all steps, median / max (the final-state threshold 1e-9 does not bound it); orbits with the L ratio above 1: $(count(>(1), gL)), with the SALI ratio above 1: $(count(>(1), gS))", isempty(sd) ? "-" : "$(fmt(sort(sd)[max(end ÷ 2, 1)])) / $(fmt(maximum(sd)))", "reported", "INFO")
end

want("sali") && let
    x0, ncand, nrep = regular_starts(1024)
    N = 4 * 10^4
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC18a", "K=0.5 regular class, $T, n=$(fmt(N)): each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        sali_d68("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, exact(T, 0.5), xT, N)
    end
    for T in (Float64, Double64)
        ac12_chaotic("AC18c", "SALI, -L/(2n)", sali_est, f -> SVector(-f.L / (2 * 10^4)), 1, T)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC18a", "K=10 chaotic class, $T, n=8: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
        sali_d68("K=10 chaotic class", T, exact(T, 10), xT, 8; chaotic = true)
    end
end

# ------------------------------------------------------------------------------------------------ AC16 (fli-ofli), tag "fli"
# Spec criterion 16 of fli-ofli-spec.md, the D68 form from the start (the same-orbit bound B_FLI applies per backend; the GPU and the CPU do not produce
# the same orbit). FLI reducer (k = 2, chain {FLI()}), lifted, tau = 1; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64 if the kernel
# does not compile). Classes: regular K = 0.5 (S7 screening, regular_starts), 1024 orbits, n = 4e4; chaotic K = 10, 1024 orbits, n = 8 and n = 10 (FLI has no
# floor, so n is not limited by a resolved window; the perturbative horizon 16 eps M_n Gamma_n <= 1e-4 of the D68 term decides which orbits are gated).
# (a) each backend against a BigFloat(256) reference formed at its own realized orbit points: |FLI_b - FLI_ref,b| <= B_FLI,b (two-point rule of Q11, Psi from
# the backend's own tangents) for FLI_B and FLI_1. (b) CPU-GPU gap: regular class, on the orbits whose final states differ by at most 1e-9 in every
# coordinate, |FLI_cpu - FLI_gpu| <= B_FLI,cpu + B_FLI,gpu; chaotic class, + max(D_orb^(i_a)(j_a), D_orb^(i_b)(j_b)), D_orb^(i)(j) the first line of the D68
# term per column (the sum over i dropped) at the computed maximiser pairs of the two backends (the finals stepB, colB, step1 come back from the device). The
# gap in units of the bound is reported. (c) K = 10, n = 1e4, tau = 2: tangent-core AC12(c) tolerances on FLI_B/n and FLI_1/n.
# The GPU orbit is read with the trace kernel of the sali section (the run keeps no orbit); its final state is checked bitwise against the propagate result.
MODTAG[] = "fli: "
want("fli") && include(joinpath(ROOT, "test", "fli-ofli", "fl_helpers.jl"))   # fl_*: B_FLI, BigFloat reference, Psi_n, fl_Orbit
fli_est(sys, x0, N, tau, b) = propagate(sys, x0, N; k = 2, tau = tau, seed = SEED, reducers = (FLI(),), backend = b)
d16_rat(x, B) = x <= 0 ? 0.0 : (B == 0 ? Inf : x / B)

# BigFloat(256) reference of orbit j at the points `pts` of backend result `res`: errors, two-point bounds and the reference values for both forms.
function d16_one(res, pts, Kt, sysT, j, N, T, epsT)
    W0 = fl_frame(SEED, j, T)
    R = fl_ref_steps(fl_mats_std(Kt, pts), W0, N)
    Ms = [SMatrix{2,2,Float64}(map(Float64, tangent(sysT, pt))) for pt in pts]
    o = fl_Orbit(R.q1, R.q2, R.G, Ms, W0, Dict{Tuple{Int,Int},Float64}(), Dict{Int,Vector{NamedTuple}}())
    r = only(fl_ref_at(R, 1, [N]))
    f = res.finals[j][1]
    B1 = fl_BFLI(epsT, o, 1, f.step1, 1, r.j1)
    BB = fl_BFLI(epsT, o, f.colB, f.stepB, r.cB, r.jB)
    e1, eB = fl_err(f.fli1, r.F1), fl_err(f.fli, r.FB)
    return (r1 = d16_rat(e1, B1), rB = d16_rat(eB, BB), B1 = B1, BB = BB, e1 = e1, eB = eB, F1 = Float64(r.F1), FB = Float64(r.FB), j1 = r.j1, jB = r.jB, cB = r.cB)
end

# D_orb^(i)(j) of fli-ofli-spec.md criterion 16 (first line of the D68 term of sali-spec.md, Amendment D68, per column) in BigFloat(256) on the CPU reference
# orbit, for every j = 0..n and i = 1, 2: with A_k = DT(z_k), P_{<k} = A_{k-1} ... A_0, P^(j)_{>k} = A_{j-1} ... A_{k+1}, bnd_k = 16 eps M_k Gamma_k
# (precision-spec Equations 3, 4; the same as d68_orb of the sali section):
#   D_orb^(i)(j) = sqrt 2 |K| sum_{k<j} norm(P^(j)_{>k})_2 bnd_k norm(P_{<k} w_i) / norm(P_j w_i),  D_orb^(i)(0) = 0 (the s = 0 term is exact).
# Returns Dij (2 x (n+1), Float64), the item 4 bounds bnd_k for k = 0..n.
function d16_orb(Kt, pts, zN, W0, epsT; matsf = (K, p) -> fl_mats_std(K, p))   # matsf: the BigFloat step tangent of the module's helpers (megno passes mg_mats_std)
    return setprecision(BigFloat, 256) do
        n = length(pts)
        mats = matsf(Kt, pts)
        A = [(m = mats(k + 1); SMatrix{2,2,BigFloat}(m[1], m[3], m[2], m[4])) for k in 0:n-1]   # A[k+1] = DT(z_k)
        zs = vcat(pts, [zN])
        Mk = accumulate(max, [max(one(BigFloat), abs(BigFloat(z[1])), abs(BigFloat(z[2]))) for z in zs])
        Gam = map(0:n) do k
            Q, g = SMatrix{2,2,BigFloat}(1, 0, 0, 1), zero(BigFloat)
            for j in k-1:-1:0            # Q = A_{k-1} ... A_{j+1}
                g += d68_sig2(Q)
                Q = Q * A[j+1]
            end
            g
        end
        bnd = [16 * BigFloat(epsT) * Mk[k+1] * Gam[k+1] for k in 0:n]
        v = [[SVector{2,BigFloat}(W0[1, i], W0[2, i])] for i in 1:2]
        for i in 1:2, k in 0:n-1
            push!(v[i], A[k+1] * v[i][end])           # v[i][k+1] = P_{<k} w_i
        end
        nv = [norm.(v[i]) for i in 1:2]
        Dij = zeros(2, n + 1)
        for j in 1:n
            Pg = Vector{SMatrix{2,2,BigFloat,4}}(undef, j)     # Pg[k+1] = P^(j)_{>k}, k = 0..j-1
            Pg[j] = SMatrix{2,2,BigFloat}(1, 0, 0, 1)
            for k in j-1:-1:1
                Pg[k] = Pg[k+1] * A[k+1]
            end
            for i in 1:2
                Dij[i, j+1] = Float64(sqrt(BigFloat(2)) * abs(BigFloat(Kt)) * sum(d68_sig2(Pg[k+1]) * bnd[k+1] * nv[i][k+1] / nv[i][j+1] for k in 0:j-1))
            end
        end
        (Dij = Dij, bnd = Float64.(bnd))
    end
end

# Chaotic class of part (b). Entering orbits: 16 eps M_n Gamma_n <= 1e-4. Gate: |FLI_cpu - FLI_gpu| <= B_FLI,cpu + B_FLI,gpu + max(D_orb^(i_a)(j_a), D_orb^(i_b)(j_b)).
function d16_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    enter = findall(j -> orb[j].hor <= 1e-4, 1:m)
    rest = setdiff(1:m, enter)
    row("AC16b", "$lab: orbits that enter part (b): 16 eps M_n Gamma_n <= 1e-4 (horizon max over all orbits $(fmt(maximum(o.hor for o in orb))); reported only: $(length(rest)))", length(enter), "reported", "INFO")
    isempty(enter) && return
    Dat(j, i, s) = s == 0 ? 0.0 : orb[j].Dij[i, s+1]
    fc(j) = cpu.finals[j][1]
    fg(j) = gpu.finals[j][1]
    termB(j) = max(Dat(j, max(fc(j).colB, 1), fc(j).stepB), Dat(j, max(fg(j).colB, 1), fg(j).stepB))
    term1(j) = max(Dat(j, 1, fc(j).step1), Dat(j, 1, fg(j).step1))
    gapB(j) = abs(Float64(fc(j).fli) - Float64(fg(j).fli))
    gap1(j) = abs(Float64(fc(j).fli1) - Float64(fg(j).fli1))
    ratB(j) = gapB(j) / (rc[j].BB + rg[j].BB + termB(j))
    rat1(j) = gap1(j) / (rc[j].B1 + rg[j].B1 + term1(j))
    row("AC16b", "$lab: entering orbits: max |FLI_B,cpu - FLI_B,gpu| / (B_FLI,cpu + B_FLI,gpu + max(D_orb^(i_a)(j_a), D_orb^(i_b)(j_b))) over $(length(enter)) orbits", maximum(ratB(j) for j in enter), "<= 1", pf(maximum(ratB(j) for j in enter) <= 1))
    row("AC16b", "$lab: entering orbits: same for FLI_1 (column 1 at step1 of each backend)", maximum(rat1(j) for j in enter), "<= 1", pf(maximum(rat1(j) for j in enter) <= 1))
    med(v) = isempty(v) ? NaN : sort(v)[max(end ÷ 2, 1)]
    q = [gapB(j) / max(termB(j), eps()) for j in enter]
    nz = [gapB(j) for j in enter if gapB(j) > 0]
    row("AC16b", "$lab: entering orbits: FLI_B gap / D_orb term, median / max (D_orb term median $(fmt(med([termB(j) for j in enter]))), max $(fmt(maximum(termB(j) for j in enter)))); nonzero gaps: $(length(nz)) of $(length(enter)), median among them $(fmt(med(nz)))",
        "$(fmt(med(q))) / $(fmt(maximum(q)))", "reported (how loose the term is)", "INFO")
    row("AC16b", "$lab: entering orbits: FLI_B gap over (B_FLI,cpu + B_FLI,gpu) without the D_orb term, max ($(count(j -> gapB(j) > rc[j].BB + rg[j].BB, enter)) of $(length(enter)) above 1); realized |FLI_ref,cpu - FLI_ref,gpu| (BigFloat series on the two orbits) over the full bound, max",
        "$(fmt(maximum(gapB(j) / (rc[j].BB + rg[j].BB) for j in enter))); $(fmt(maximum(abs(rc[j].FB - rg[j].FB) / (rc[j].BB + rg[j].BB + termB(j)) for j in enter)))", "reported", "INFO")
    rref(j) = gapB(j) / (rc[j].BB + rg[j].BB + max(Dat(j, max(rc[j].cB, 1), rc[j].jB), Dat(j, max(rg[j].cB, 1), rg[j].jB)))
    row("AC16b", "$lab: entering orbits: same ratio with the maximisers of the BigFloat series (jB, cB of each backend's reference) in the D_orb term, max", maximum(rref(j) for j in enter), "reported", "INFO")
    nmd = count(j -> (fc(j).stepB, fc(j).colB) != (rc[j].jB, rc[j].cB) || (fg(j).stepB, fg(j).colB) != (rg[j].jB, rg[j].cB), enter)
    row("AC16b", "$lab: entering orbits: orbits whose computed maximiser (stepB, colB) differs from the reference maximiser on either backend; orbits whose CPU and GPU maximisers differ", "$nmd; $(count(j -> (fc(j).stepB, fc(j).colB) != (fg(j).stepB, fg(j).colB), enter))", "reported", "INFO")
    pk, lt = maximum(orb[j].perk for j in enter), maximum(orb[j].lit for j in enter)
    row("AC16b", "$lab: item 4 check on the device, entering orbits: max over orbits and k = 1..n of |z_k^cpu - z_k^gpu|_inf / (16 eps M_k Gamma_k); same against the literal 16 eps M_n Gamma_n (orbits above 1: $(count(j -> orb[j].perk > 1, enter)) / $(count(j -> orb[j].lit > 1, enter)))",
        "$(fmt(pk)); $(fmt(lt))", "<= 1; <= 1", pf(pk <= 1 && lt <= 1))
    if !isempty(rest)
        row("AC16b", "$lab: reported-only orbits ($(length(rest))): max gap ratios (FLI_B; FLI_1) over the full bound; item 4 per-k ratio max", "$(fmt(maximum(ratB, rest))); $(fmt(maximum(rat1, rest))); $(fmt(maximum(orb[j].perk for j in rest)))", "reported", "INFO")
    end
end

function fli_d16(tag, T, Kt, xT, N; chaotic = false)
    sysT = StandardMap(Kt)
    m = length(xT)
    lab = "$tag, $T, $m orbits, n=$(fmt(N)), tau=1"
    gpu = fli_est(sysT, xT, N, 1, :gpu)
    why = gpu_why(gpu)
    why === nothing || return row("AC16a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    cpu = fli_est(sysT, xT, N, 1, :cpu)
    epsT = Float64(effective_eps(T))
    P, XN = gpu_trace(sysT, xT, N)
    rc = Vector{NamedTuple}(undef, m); rg = Vector{NamedTuple}(undef, m)
    same_c = zeros(Bool, m)
    stepdiff = zeros(m)
    orb = Vector{Any}(undef, m)
    Threads.@threads :dynamic for j in 1:m
        ptsc = fl_pts(sysT, xT[j], N)
        zc = forward(sysT, ptsc[end])
        same_c[j] = zc == cpu.x[j]
        stepdiff[j] = maximum(maximum(abs.(ptsc[i] - P[j, i])) for i in 1:N)
        rc[j] = d16_one(cpu, ptsc, Kt, sysT, j, N, T, epsT)
        rg[j] = d16_one(gpu, P[j, :], Kt, sysT, j, N, T, epsT)
        if chaotic
            od = d16_orb(Kt, ptsc, zc, fl_frame(SEED, j, T), epsT)
            e = [i < N ? maximum(abs.(ptsc[i+1] - P[j, i+1])) : maximum(abs.(zc - XN[j])) for i in 1:N]   # e_k = |z_k^cpu - z_k^gpu|_inf
            orb[j] = (Dij = od.Dij, hor = od.bnd[end], perk = maximum(e ./ od.bnd[2:end]), lit = maximum(e) / od.bnd[end])
        end
    end
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    row("AC16a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        for (form, key) in (("FLI_B", :rB), ("FLI_1", :r1))
            r = [getproperty(x, key) for x in R]
            row("AC16a", "$lab: backend $nm (provenance $(res.provenance.backend)): max |$form - $form,ref| / B_FLI over $m orbits (median $(fmt(sort(r)[end ÷ 2])))", maximum(r), "<= 1", pf(maximum(r) <= 1))
        end
        nmd = count(j -> (res.finals[j][1].stepB, res.finals[j][1].colB) != (R[j].jB, R[j].cB), 1:m)
        row("AC16a", "$lab: backend $nm: orbits whose computed FLI_B maximiser (stepB, colB) differs from the BigFloat series' (reported); largest FLI_B $(fmt(maximum(R[j].FB for j in 1:m)))", nmd, "reported", "INFO")
    end
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = findall(<=(1e-9), dx)
    row("AC16b", "$lab: orbits with final states apart by more than 1e-9, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(sort(dx)[end ÷ 2])), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    Kd = exact(Double64, Kt)
    ref64 = fli_est(StandardMap(Kd), [SVector{2,Double64}(z) for z in xT], N, 1, :cpu)
    dd = [Float64(maximum(abs.(SVector{2,Double64}(cpu.x[j]) - ref64.x[j]))) for j in 1:m]
    row("AC16b", "$lab: orbits with final states apart by more than 1e-9, CPU Float64 vs Double64 CPU (max $(fmt(maximum(dd))))", count(>(1e-9), dd), "reported", "INFO")
    if chaotic
        row("AC16b", "$lab: orbits with final states within 1e-9 in every coordinate (the regular-class selection; not the entering rule for this class)", length(agree), "reported", "INFO")
        return d16_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    end
    row("AC16b", "$lab: orbits that enter part (b) (final states within 1e-9 in every coordinate)", length(agree), "reported", "INFO")
    gB = [abs(Float64(cpu.finals[j][1].fli) - Float64(gpu.finals[j][1].fli)) / (rc[j].BB + rg[j].BB) for j in agree]
    g1 = [abs(Float64(cpu.finals[j][1].fli1) - Float64(gpu.finals[j][1].fli1)) / (rc[j].B1 + rg[j].B1) for j in agree]
    row("AC16b", "$lab: orbits that agree: max |FLI_B,cpu - FLI_B,gpu| / (B_FLI,cpu + B_FLI,gpu) over $(length(agree)) orbits", isempty(gB) ? NaN : maximum(gB), "<= 1", pf(!isempty(gB) && maximum(gB) <= 1))
    row("AC16b", "$lab: orbits that agree: max |FLI_1,cpu - FLI_1,gpu| / (B_FLI,cpu + B_FLI,gpu)", isempty(g1) ? NaN : maximum(g1), "<= 1", pf(!isempty(g1) && maximum(g1) <= 1))
    sd = stepdiff[agree]
    neq = count(j -> (cpu.finals[j][1].stepB, cpu.finals[j][1].colB) == (gpu.finals[j][1].stepB, gpu.finals[j][1].colB), agree)
    row("AC16b", "$lab: diagnostic on the orbits that agree: largest CPU-GPU orbit-point difference over all steps, median / max (the final-state threshold 1e-9 does not bound it); orbits with the FLI_B ratio above 1: $(count(>(1), gB)), FLI_1: $(count(>(1), g1)); maximiser (stepB, colB) equal on CPU and GPU for $neq of $(length(agree))",
        isempty(sd) ? "-" : "$(fmt(sort(sd)[max(end ÷ 2, 1)])) / $(fmt(maximum(sd)))", "reported", "INFO")
end

want("fli") && let
    x0, ncand, nrep = regular_starts(1024)
    N = 4 * 10^4
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC16a", "K=0.5 regular class, $T, n=$(fmt(N)): each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        fli_d16("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, exact(T, 0.5), xT, N)
    end
    for T in (Float64, Double64)
        ac12_chaotic("AC16c", "FLI_B/n, FLI_1/n", fli_est, f -> SVector(f.fli / 10^4, f.fli1 / 10^4), 2, T)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC16a", "K=10 chaotic class, $T: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
        for nch in (8, 10)
            fli_d16("K=10 chaotic class", T, exact(T, 10), xT, nch; chaotic = true)
        end
    end
end


# ------------------------------------------------------------------------------------------------ AC15 (megno), tag "megno"
# Spec criterion 15 of megno-spec.md, the D68 form (the bound B_M applies per backend; the GPU and the CPU do not produce the same orbit). MEGNO reducer (k = 2, chain
# {MEGNO()}, the replacer), lifted, tau = 1; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64 if the kernel does not compile). Classes: regular
# K = 0.5 (S7 screening, regular_starts), 1024 orbits, n = 4e4; chaotic K = 10, 1024 orbits, n = 8 and n = 10 (MEGNO has no floor; as for FLI the perturbative horizon
# 16 eps M_n Gamma_n <= 1e-4 of the D68 term decides which orbits are gated; at these n the "grid" of M11 is every handoff).
# (a) each backend against a BigFloat(256) reference formed at its own realized orbit points: |M_b - M_ref,b| <= B_M,b(n) = 48 eps Psi_max^(i) (Psi_max over the 65-handoff
# grid g_k = round(n k / 64), k = 1..64, and g = 1, from the backend's own tangents), the three series Y, Ybar, Ybar20, both columns. (b) CPU-GPU gap: regular class, on the
# orbits whose final states differ by at most 1e-9 in every coordinate, |M_cpu - M_gpu| <= B_M,cpu + B_M,gpu; chaotic class, + 6 max over the grid of D_orb^(i)(j) (first
# line of the D68 term per column, d16_orb with the megno tangent), the gap in units of the bound reported. (c) K = 10, n = 1e4, tau = 2: tangent-core AC12(c) tolerances
# on Y_1/n and Ybar_1/n. The GPU orbit is read with the trace kernel of the sali section; its final state is checked bitwise against the propagate result.
MODTAG[] = "megno: "
want("megno") && include(joinpath(ROOT, "test", "megno", "mg_helpers.jl"))   # mg_*: BigFloat step tangent mg_mats_std, M5, B_M ingredients, Psi_n
megno_est(sys, x0, N, tau, b) = propagate(sys, x0, N; k = 2, tau = tau, seed = SEED, reducers = (MEGNO(),), backend = b)
m15_grid(N) = unique(vcat(1, [max(1, round(Int, N * k / 64)) for k in 1:64]))

# Streaming BigFloat(256) reference of both columns at the final step (tau = 1): the M5 sums Y, Ybar, Ybar20 and the log sum s, from the BigFloat step tangents.
function m15_ref(mats, W0, N)
    return setprecision(BigFloat, 256) do
        w = [(mg_big(W0[1, i]), mg_big(W0[2, i])) for i in 1:2]
        l0 = [log(w[i][1]^2 + w[i][2]^2) / 2 for i in 1:2]
        prev = [zero(BigFloat) for _ in 1:2]
        R1, R2, Q, Q2 = ([zero(BigFloat) for _ in 1:2] for _ in 1:4)
        for k in 1:N
            ma, mb, mc, md = mats(k)
            t = BigFloat(k)
            for i in 1:2
                w[i] = (ma * w[i][1] + mb * w[i][2], mc * w[i][1] + md * w[i][2])
                sk = log(w[i][1]^2 + w[i][2]^2) / 2 - l0[i]
                a = sk - prev[i]
                prev[i] = sk
                R1[i] += t * a
                R2[i] += t^2 * a
                Q[i] += 2 * R1[i] / t
                Q2[i] += 3 * R2[i]
            end
        end
        (Y = Tuple(2 .* R1 ./ N), Ybar = Tuple(Q ./ N), Ybar20 = Tuple(Q2 ./ BigFloat(N)^3), s = Tuple(prev))
    end
end

# Reference, bound and ratios of orbit j of backend result `res` at the points `pts`.
function m15_one(res, pts, Kt, sysT, j, N, T, epsT, grid)
    W0 = mg_frame(SEED, j, T)
    ref = m15_ref(mg_mats_std(Kt, pts), W0, N)
    Ms = [SMatrix{2,2,Float64}(map(Float64, tangent(sysT, pt))) for pt in pts]
    B = ntuple(i -> 48 * epsT * maximum(lb_Psi(view(Ms, 1:g), W0[:, i]) for g in grid if g <= N), 2)
    f = res.finals[j][1]
    e = (Y = ntuple(i -> mg_err(f.Y[i], ref.Y[i]), 2), Ybar = ntuple(i -> mg_err(f.Ybar[i], ref.Ybar[i]), 2), Ybar20 = ntuple(i -> mg_err(f.Ybar20[i], ref.Ybar20[i]), 2),
         S = ntuple(i -> mg_err(f.sum[i] + f.comp[i], ref.s[i]), 2))
    r(x, b) = b == 0 ? (x == 0 ? 0.0 : Inf) : x / b
    return (e = e, B = B, rY = maximum(r(e.Y[i], B[i]) for i in 1:2), rYb = maximum(r(e.Ybar[i], B[i]) for i in 1:2), rY20 = maximum(r(e.Ybar20[i], B[i]) for i in 1:2),
            rS = maximum(r(e.S[i], B[i] / 6) for i in 1:2), ref = ref)
end

function m15_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    enter = findall(j -> orb[j].hor <= 1e-4, 1:m)
    rest = setdiff(1:m, enter)
    row("AC15b", "$lab: orbits that enter part (b): 16 eps M_n Gamma_n <= 1e-4 (horizon max over all orbits $(fmt(maximum(o.hor for o in orb))); reported only: $(length(rest)))", length(enter), "reported", "INFO")
    isempty(enter) && return
    fc(j) = cpu.finals[j][1]
    fg(j) = gpu.finals[j][1]
    Dg(j, i) = maximum(orb[j].Dij[i, 2:end])                    # max over the grid (every handoff 1..n) of D_orb^(i)(g)
    gapr(j, f, withD) = maximum(abs(Float64(getproperty(fc(j), f)[i]) - Float64(getproperty(fg(j), f)[i])) / (rc[j].B[i] + rg[j].B[i] + (withD ? 6 * Dg(j, i) : 0.0)) for i in 1:2)
    for (f, nm) in ((:Y, "Y"), (:Ybar, "Ybar"), (:Ybar20, "Ybar20"))
        r = [gapr(j, f, true) for j in enter]
        row("AC15b", "$lab: entering orbits: max |$nm,cpu - $nm,gpu| / (B_M,cpu + B_M,gpu + 6 max_g D_orb^(i)(g)) over $(length(enter)) orbits, both columns", maximum(r), "<= 1", pf(maximum(r) <= 1))
    end
    med(v) = isempty(v) ? NaN : sort(v)[max(end ÷ 2, 1)]
    plain = [gapr(j, :Y, false) for j in enter]
    row("AC15b", "$lab: entering orbits: Y gap over (B_M,cpu + B_M,gpu) without the D_orb term, max ($(count(>(1), plain)) of $(length(enter)) above 1); max for Ybar; Ybar20; D_orb term (6 max_g D_orb, column 1) median / max",
        "$(fmt(maximum(plain))); $(fmt(maximum(gapr(j, :Ybar, false) for j in enter))); $(fmt(maximum(gapr(j, :Ybar20, false) for j in enter))); $(fmt(med([6 * Dg(j, 1) for j in enter]))) / $(fmt(maximum(6 * Dg(j, 1) for j in enter)))", "reported (how loose the term is)", "INFO")
    pk, lt = maximum(orb[j].perk for j in enter), maximum(orb[j].lit for j in enter)
    row("AC15b", "$lab: item 4 check on the device, entering orbits: max over orbits and k = 1..n of |z_k^cpu - z_k^gpu|_inf / (16 eps M_k Gamma_k); same against the literal 16 eps M_n Gamma_n (orbits above 1: $(count(j -> orb[j].perk > 1, enter)) / $(count(j -> orb[j].lit > 1, enter)))",
        "$(fmt(pk)); $(fmt(lt))", "<= 1; <= 1", pf(pk <= 1 && lt <= 1))
    if !isempty(rest)
        row("AC15b", "$lab: reported-only orbits ($(length(rest))): max gap ratios (Y; Ybar; Ybar20) over the full bound; item 4 per-k ratio max", "$(fmt(maximum(gapr(j, :Y, true) for j in rest))); $(fmt(maximum(gapr(j, :Ybar, true) for j in rest))); $(fmt(maximum(gapr(j, :Ybar20, true) for j in rest))); $(fmt(maximum(orb[j].perk for j in rest)))", "reported", "INFO")
    end
end

function megno_d15(tag, T, Kt, xT, N; chaotic = false)
    sysT = StandardMap(Kt)
    m = length(xT)
    lab = "$tag, $T, $m orbits, n=$(fmt(N)), tau=1"
    gpu = megno_est(sysT, xT, N, 1, :gpu)
    why = gpu_why(gpu)
    why === nothing || return row("AC15a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    cpu = megno_est(sysT, xT, N, 1, :cpu)
    epsT = Float64(effective_eps(T))
    grid = m15_grid(N)
    P, XN = gpu_trace(sysT, xT, N)
    rc = Vector{NamedTuple}(undef, m); rg = Vector{NamedTuple}(undef, m)
    same_c = zeros(Bool, m)
    stepdiff = zeros(m)
    orb = Vector{Any}(undef, m)
    Threads.@threads :dynamic for j in 1:m
        ptsc = mg_pts(sysT, xT[j], N)
        zc = forward(sysT, ptsc[end])
        same_c[j] = zc == cpu.x[j]
        stepdiff[j] = maximum(maximum(abs.(ptsc[i] - P[j, i])) for i in 1:N)
        rc[j] = m15_one(cpu, ptsc, Kt, sysT, j, N, T, epsT, grid)
        rg[j] = m15_one(gpu, P[j, :], Kt, sysT, j, N, T, epsT, grid)
        if chaotic
            od = d16_orb(Kt, ptsc, zc, mg_frame(SEED, j, T), epsT; matsf = mg_mats_std)
            e = [i < N ? maximum(abs.(ptsc[i+1] - P[j, i+1])) : maximum(abs.(zc - XN[j])) for i in 1:N]   # e_k = |z_k^cpu - z_k^gpu|_inf
            orb[j] = (Dij = od.Dij, hor = od.bnd[end], perk = maximum(e ./ od.bnd[2:end]), lit = maximum(e) / od.bnd[end])
        end
    end
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    row("AC15a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        for (form, key) in (("Y", :rY), ("Ybar", :rYb), ("Ybar20", :rY20))
            r = [getproperty(x, key) for x in R]
            row("AC15a", "$lab: backend $nm (provenance $(res.provenance.backend)): max |$form - $form,ref| / B_M over $m orbits, both columns (median $(fmt(sort(r)[end ÷ 2])))", maximum(r), "<= 1", pf(maximum(r) <= 1))
        end
        row("AC15a", "$lab: backend $nm: max |s_i - s_ref| / (8 eps Psi_max) over the grid, both columns (reported, the borrowed constant of D74 item 7)", maximum(x.rS for x in R), "reported", "INFO")
    end
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = findall(<=(1e-9), dx)
    row("AC15b", "$lab: orbits with final states apart by more than 1e-9, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(sort(dx)[end ÷ 2])), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    Kd = exact(Double64, Kt)
    ref64 = megno_est(StandardMap(Kd), [SVector{2,Double64}(z) for z in xT], N, 1, :cpu)
    dd = [Float64(maximum(abs.(SVector{2,Double64}(cpu.x[j]) - ref64.x[j]))) for j in 1:m]
    row("AC15b", "$lab: orbits with final states apart by more than 1e-9, CPU Float64 vs Double64 CPU (max $(fmt(maximum(dd))))", count(>(1e-9), dd), "reported", "INFO")
    if chaotic
        row("AC15b", "$lab: orbits with final states within 1e-9 in every coordinate (the regular-class selection; not the entering rule for this class)", length(agree), "reported", "INFO")
        return m15_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    end
    row("AC15b", "$lab: orbits that enter part (b) (final states within 1e-9 in every coordinate)", length(agree), "reported", "INFO")
    fc(j) = cpu.finals[j][1]
    fg(j) = gpu.finals[j][1]
    for (f, nm) in ((:Y, "Y"), (:Ybar, "Ybar"), (:Ybar20, "Ybar20"))
        g = [maximum(abs(Float64(getproperty(fc(j), f)[i]) - Float64(getproperty(fg(j), f)[i])) / (rc[j].B[i] + rg[j].B[i]) for i in 1:2) for j in agree]
        row("AC15b", "$lab: orbits that agree: max |$nm,cpu - $nm,gpu| / (B_M,cpu + B_M,gpu) over $(length(agree)) orbits, both columns (median $(fmt(sort(g)[max(end ÷ 2, 1)])))", isempty(g) ? NaN : maximum(g), "<= 1", pf(!isempty(g) && maximum(g) <= 1))
    end
    sd = stepdiff[agree]
    row("AC15b", "$lab: diagnostic on the orbits that agree: largest CPU-GPU orbit-point difference over all steps, median / max (the final-state threshold 1e-9 does not bound it)", isempty(sd) ? "-" : "$(fmt(sort(sd)[max(end ÷ 2, 1)])) / $(fmt(maximum(sd)))", "reported", "INFO")
end

# Kernel compile record and time per orbit-step against FLI and SALI (K = 2 pi, tau = 1, chain of one reducer), GPU and CPU.
function megno_timing()
    stsz(r, T) = sizeof(typeof(Chaoslab.reducer_init(r, T, Val(2), Val(2))))
    row("AC15", "reducer state size in bytes: MEGNO Float64 / Double64; FLI Float64; SALI Float64", "$(stsz(MEGNO(), Float64)) / $(stsz(MEGNO(), Double64)); $(stsz(FLI(), Float64)); $(stsz(SALI(), Float64))", "reported", "INFO")
    norb, N = 20_000, 1_000
    rng = StableRNG(SEED + 15)
    x0 = [SVector{2,Float64}(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:norb]
    sys = StandardMap(2pi)
    for (nm, r) in (("MEGNO", MEGNO()), ("FLI", FLI()), ("SALI", SALI()))
        kw = (k = 2, tau = 1, seed = SEED, reducers = (r,))
        cpu_t = gpu_t = Inf
        prov = nothing
        try
            propagate(sys, x0[1:64], 10; kw..., backend = :gpu, chunk = AC13_GPU_CHUNK)       # warm-up (compilation of the kernel)
            for _ in 1:3
                t = @elapsed res = propagate(sys, x0, N; kw..., backend = :gpu, chunk = AC13_GPU_CHUNK)
                prov = res.provenance
                gpu_t = min(gpu_t, t)
            end
        catch e
            row("AC15", "$nm kernel on the device: compile / run error", sprint(showerror, e)[1:min(end, 200)], "compiles", "FAIL")
            continue
        end
        propagate(sys, x0[1:64], 10; kw..., backend = :cpu)
        for _ in 1:2
            cpu_t = min(cpu_t, @elapsed propagate(sys, x0, N; kw..., backend = :cpu))
        end
        steps = norb * N
        row("AC15", "$nm kernel (provenance backend $(prov.backend)); $norb orbits x $N steps, K = 2 pi, tau = 1: GPU $(fmt(gpu_t / steps * 1e9)) ns/orbit-step, CPU ($(Threads.nthreads()) threads) $(fmt(cpu_t / steps * 1e9)) ns/orbit-step, speedup $(fmt(cpu_t / gpu_t))",
            prov.backend === :gpu ? gpu_t / steps * 1e9 : NaN, "recorded", prov.backend === :gpu ? "INFO" : fb(prov.backend_fallback_reason))
    end
end

want("megno") && let
    megno_timing()
    x0, ncand, nrep = regular_starts(1024)
    N = 4 * 10^4
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=0.5 regular class, $T, n=$(fmt(N)): each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        megno_d15("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, exact(T, 0.5), xT, N)
    end
    for T in (Float64, Double64)
        ac12_chaotic("AC15c", "Y_1/n, Ybar_1/n", megno_est, f -> SVector(f.Y[1] / 10^4, f.Ybar[1] / 10^4), 2, T)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=10 chaotic class, $T: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
        for nch in (8, 10)
            megno_d15("K=10 chaotic class", T, exact(T, 10), xT, nch; chaotic = true)
        end
    end
end

# ------------------------------------------------------------------------------------------------ AC15 (rli), tag "rli"
# Spec criterion 15 of rli-spec.md, the D68 form (the bound B_R applies per backend; the GPU and the CPU do not produce the same orbit, so each backend is compared with a BigFloat(256)
# reference formed at its own realized points of BOTH orbits of the pair). RLI reducer (k = 2, d = 4, chain {RLI()}, the replacer) on PairMap(WrappedStandardMap), supplied block
# frame, Dx = 1e-12 along xi_0, tau = 1; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64 if the kernel does not compile). Classes: regular K = 0.5 (S7
# screening, regular_starts), 1024 pairs, n = 4e4; chaotic K = 10, 1024 pairs, n = 8 and n = 10 (part (b) with the orbit-divergence term) and n = 100 (part (a) only).
# (a) each backend against its own reference: |RLI_b - RLI_ref,b| <= B_R,b(g) = 8 eps (Psi_max^A + Psi_max^B) H_g / g + 4 eps RLI_ref (RL11, Psi_max over the 65-handoff grid
# g_k = round(n k / 64), k = 1..64, and g = 1, from the backend's own tangents) at EVERY grid handoff g (one propagate call per g, the orbit is a prefix), and
# |DL_b - DL_ref,b| <= 8 eps (Psi_g^A + Psi_g^B) / g + 4 eps DL_ref. (b) CPU-GPU gap: regular class, on the pairs whose final states differ by at most 1e-9 in every coordinate (all 8
# components), |RLI_cpu - RLI_gpu| <= B_R,cpu + B_R,gpu; chaotic class (D68 term of the sali Amendment D68, first line, per orbit via d16_orb) + (1/n) sum_{k<=n} (D_orb^A(k) + D_orb^B(k)) / k,
# the mean of the s-gap over the orbit divergence of both orbits, on the pairs that enter by the item 4 horizon 16 eps M_n Gamma_n <= 1e-4 for both orbits. Spec 15(b) as written is
# statistical (reported, not gated): the distribution of RLI_cpu / RLI_gpu, the count of pairs with final states apart by more than 1e-9, and for K = 10, n = 1e4,
# tau = 2 the tangent-core AC12(c) tolerances (6 sigma per pair, 4 SE for the means). The GPU orbit is read with a trace kernel that runs the same `forward` on the device (the run keeps
# no orbit); its final state is checked bitwise against the propagate result. Each trace is of one orbit of the pair (the base map's wrapped forward on the device).
MODTAG[] = "rli: "
want("rli") && include(joinpath(ROOT, "test", "rli", "rl_helpers.jl"))   # rl_*: BigFloat step tangent rl_mats_std, RL11 bounds rl_BR, rl_BDL, rl_frame, Psi_n (lb_Psi)
r15_x(z) = z[1]
r15_id(z) = z
function r15_trace_kernel!(D, X, sys, steps, proj)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if j <= length(X)
        x = X[j]
        for i in 1:steps
            D[j, i] = proj(x)
            x = forward(sys, x)
        end
        X[j] = x
    end
    return nothing
end
# D[j, i] = proj(state of orbit j before step i), X advanced; chunked along the steps as gpu_trace.
function r15_trace(sys, x0::Vector{S}, N, proj; chunk = 5000) where {S}
    R = typeof(proj(x0[1]))
    m = length(x0)
    X = CuArray(x0)
    P = Matrix{R}(undef, m, N)
    for s0 in 0:chunk:N-1
        s1 = min(s0 + chunk, N)
        D = CuArray{R}(undef, m, s1 - s0)
        @cuda threads = 256 blocks = cld(m, 256) r15_trace_kernel!(D, X, sys, s1 - s0, proj)
        P[:, s0+1:s1] = Array(D)
        CUDA.unsafe_free!(D)
    end
    return P, Array(X)
end

# The pair ensemble of the class: base map, pair states, block frames, shared xi_0.
function r15_setup(T, Kt, x0, Dx)
    base = WrappedStandardMap(Kt)
    xis = [SVector{2,T}(initial_deviations(SEED, j, Val(2), Val(1), T)[:, 1]) for j in eachindex(x0)]
    ps = [pair_start(base, x0[j], xis[j], T(Dx); orbit = j) for j in eachindex(x0)]
    return (base = base, P = PairMap(base), xis = xis, z0 = [p.z for p in ps], W0 = [rl_frame(x) for x in xis], Dxr = [p.Dx_r for p in ps])
end
rli_est(S, N, tau, b; chunk = nothing) = propagate(S.P, S.z0, N; k = 2, tau = tau, seed = SEED, reducers = (RLI(),), W0 = S.W0, backend = b, chunk = chunk)

# Streaming BigFloat(256) reference of a pair: DL_k and RLI_k (RL4, tau = 1) at the handoffs in `grid`, from the step tangents of both orbits and the shared xi_0.
function r15_ref(matsA, matsB, xi0, N, grid)
    return setprecision(BigFloat, 256) do
        w = [(rl_big(xi0[1]), rl_big(xi0[2])) for _ in 1:2]
        l0 = log(w[1][1]^2 + w[1][2]^2) / 2
        s = [zero(BigFloat), zero(BigFloat)]
        Q = zero(BigFloat)
        gset = Set(grid)
        out = Dict{Int,NamedTuple}()
        for k in 1:N
            for (i, mats) in ((1, matsA), (2, matsB))
                ma, mb, mc, md = mats(k)
                w[i] = (ma * w[i][1] + mb * w[i][2], mc * w[i][1] + md * w[i][2])
                s[i] = log(w[i][1]^2 + w[i][2]^2) / 2 - l0
            end
            DL = abs(s[2] - s[1]) / k
            Q += DL
            k in gset && (out[k] = (RLI = Q / k, DL = DL))
        end
        out
    end
end
r15_Ms(Kt, xs) = [(c = Kt * cos(x); SMatrix{2,2,Float64}(1 + c, c, 1, 1)) for x in xs]    # Float64 step tangents (for Psi) from the x coordinates
mi2pi(u) = u - 2pi * round(u / 2pi)
rmed(v) = isempty(v) ? NaN : sort(v)[max((length(v) + 1) ÷ 2, 1)]

# Per backend: reference, bounds and ratios of pair j. xsA, xsB: the x coordinates of the realized orbits (steps 1..N), res_g: g -> result of the prefix run.
function r15_one(res_g, grid, j, xsA, xsB, xi0, Kt, epsT)
    ref = r15_ref(rl_mats_std(Kt, xsA), rl_mats_std(Kt, xsB), xi0, length(xsA), grid)
    MsA, MsB = r15_Ms(Kt, xsA), r15_Ms(Kt, xsB)
    psiA = [lb_Psi(view(MsA, 1:g), xi0) for g in grid]
    psiB = [lb_Psi(view(MsB, 1:g), xi0) for g in grid]
    pA, pB = accumulate(max, psiA), accumulate(max, psiB)
    rR = rD = 0.0
    BR = 0.0
    rat(e, B) = B == 0 ? (e == 0 ? 0.0 : Inf) : e / B
    for (q, g) in enumerate(grid)
        f = res_g[g].finals[j][1]
        rr = ref[g]
        BR = rl_BR(epsT, pA[q], pB[q], g, g, Float64(rr.RLI))
        BD = rl_BDL(epsT, psiA[q], psiB[q], g, Float64(rr.DL))
        rR = max(rR, rat(rl_err(f.RLI, rr.RLI), BR))
        rD = max(rD, rat(rl_err(f.DL, rr.DL), BD))
    end
    return (rR = rR, rD = rD, BR = BR, RLIref = Float64(ref[last(grid)].RLI))
end

function rli_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    enter = findall(j -> orb[j].hor <= 1e-4, 1:m)
    rest = setdiff(1:m, enter)
    row("AC15b", "$lab: pairs that enter part (b): 16 eps M_n Gamma_n <= 1e-4 for both orbits (horizon max over all pairs $(fmt(maximum(o.hor for o in orb))); reported only: $(length(rest)))", length(enter), "reported", "INFO")
    isempty(enter) && return
    fc(j) = cpu.finals[j][1]
    fg(j) = gpu.finals[j][1]
    gapr(j, withD) = abs(Float64(fc(j).RLI) - Float64(fg(j).RLI)) / (rc[j].BR + rg[j].BR + (withD ? orb[j].Dterm : 0.0))
    r = [gapr(j, true) for j in enter]
    row("AC15b", "$lab: entering pairs: max |RLI,cpu - RLI,gpu| / (B_R,cpu + B_R,gpu + (1/n) sum_k (D_orb^A(k) + D_orb^B(k))/k) over $(length(enter)) pairs", maximum(r), "<= 1", pf(maximum(r) <= 1))
    plain = [gapr(j, false) for j in enter]
    row("AC15b", "$lab: entering pairs: gap over (B_R,cpu + B_R,gpu) without the D_orb term, max ($(count(>(1), plain)) of $(length(enter)) above 1); D_orb term median / max",
        "$(fmt(maximum(plain))); $(fmt(rmed([orb[j].Dterm for j in enter]))) / $(fmt(maximum(orb[j].Dterm for j in enter)))", "reported (how loose the term is)", "INFO")
    pk, lt = maximum(orb[j].perk for j in enter), maximum(orb[j].lit for j in enter)
    row("AC15b", "$lab: item 4 check on the device, entering pairs: max over pairs, both orbits and k = 1..n of |z_k^cpu - z_k^gpu|_inf / (16 eps M_k Gamma_k); same against the literal 16 eps M_n Gamma_n (pairs above 1: $(count(j -> orb[j].perk > 1, enter)) / $(count(j -> orb[j].lit > 1, enter)))",
        "$(fmt(pk)); $(fmt(lt))", "<= 1; <= 1", pf(pk <= 1 && lt <= 1))
    if !isempty(rest)
        row("AC15b", "$lab: reported-only pairs ($(length(rest))): max gap ratio over the full bound; item 4 per-k ratio max", "$(fmt(maximum(gapr(j, true) for j in rest))); $(fmt(maximum(orb[j].perk for j in rest)))", "reported", "INFO")
    end
end

function rli_d15(tag, T, Kt, xT, N; chaotic = false, partb = true)
    m = length(xT)
    lab = "$tag, $T, $m pairs, n=$(fmt(N)), tau=1, Dx=1e-12"
    S = r15_setup(T, Kt, xT, 1e-12)
    grid = m15_grid(N)
    gres = Dict{Int,Any}()
    gres[N] = rli_est(S, N, 1, :gpu; chunk = AC13_GPU_CHUNK)
    why = gpu_why(gres[N])
    why === nothing || return row("AC15a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    gpu = gres[N]
    cpu = rli_est(S, N, 1, :cpu)
    cres = Dict{Int,Any}(N => cpu)
    for g in grid
        g == N && continue
        gres[g] = rli_est(S, g, 1, :gpu; chunk = AC13_GPU_CHUNK)
        cres[g] = rli_est(S, g, 1, :cpu)
    end
    epsT = Float64(effective_eps(T))
    base = S.base
    zA(z) = rl_half(z, 1); zB(z) = rl_half(z, 2)
    proj = chaotic ? r15_id : r15_x
    PA, XNA = r15_trace(base, [zA(z) for z in S.z0], N, proj)
    PB, XNB = r15_trace(base, [zB(z) for z in S.z0], N, proj)
    xof(v) = chaotic ? v[1] : v
    rc = Vector{NamedTuple}(undef, m); rg = Vector{NamedTuple}(undef, m)
    same_c = zeros(Bool, m)
    orb = Vector{Any}(undef, m)
    Threads.@threads :dynamic for j in 1:m
        ptsA, ptsB = rl_pts(base, zA(S.z0[j]), N), rl_pts(base, zB(S.z0[j]), N)
        zcA, zcB = forward(base, ptsA[end]), forward(base, ptsB[end])
        same_c[j] = vcat(zcA, zcB) == cpu.x[j]
        rc[j] = r15_one(cres, grid, j, [p[1] for p in ptsA], [p[1] for p in ptsB], S.xis[j], Kt, epsT)
        rg[j] = r15_one(gres, grid, j, [xof(P) for P in PA[j, :]], [xof(P) for P in PB[j, :]], S.xis[j], Kt, epsT)
        if chaotic
            W2 = SMatrix{2,2,T}(S.xis[j][1], S.xis[j][2], S.xis[j][1], S.xis[j][2])
            oA = d16_orb(Kt, ptsA, zcA, W2, epsT; matsf = rl_mats_std)
            oB = d16_orb(Kt, ptsB, zcB, W2, epsT; matsf = rl_mats_std)
            ek(P, pts, zc, XN) = [i < N ? maximum(abs.(mi2pi.(pts[i+1][1:2] - P[j, i+1][1:2]))) : maximum(abs.(mi2pi.(zc[1:2] - XN[j][1:2]))) for i in 1:N]
            eA, eB = ek(PA, ptsA, zcA, XNA), ek(PB, ptsB, zcB, XNB)
            Dterm = sum((oA.Dij[1, k+1] + oB.Dij[1, k+1]) / k for k in 1:N) / N
            orb[j] = (Dterm = Dterm, hor = max(oA.bnd[end], oB.bnd[end]), perk = maximum(max.(eA ./ oA.bnd[2:end], eB ./ oB.bnd[2:end])), lit = max(maximum(eA) / oA.bnd[end], maximum(eB) / oB.bnd[end]))
        end
    end
    same_g = [vcat(XNA[j], XNB[j]) == gpu.x[j] for j in 1:m]
    row("AC15a", "$lab: the reference orbits (both of the pair) reproduce the backend's final state bitwise, CPU / GPU (GPU orbits from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        row("AC15a", "$lab: backend $nm (provenance $(res.provenance.backend)): max |RLI - RLI,ref| / B_R over $m pairs and the $(length(grid)) grid handoffs (median over pairs $(fmt(rmed([x.rR for x in R]))))", maximum(x.rR for x in R), "<= 1", pf(maximum(x.rR for x in R) <= 1))
        row("AC15a", "$lab: backend $nm: max |DL - DL,ref| / (8 eps (Psi_g^A + Psi_g^B)/g + 4 eps DL,ref) over $m pairs and the grid handoffs", maximum(x.rD for x in R), "<= 1", pf(maximum(x.rD for x in R) <= 1))
        row("AC15a", "$lab: backend $nm: median RLI,ref(n); median B_R(n) (how far the worst-case bound is from the value)", "$(fmt(rmed([x.RLIref for x in R]))); $(fmt(rmed([x.BR for x in R])))", "reported", "INFO")
    end
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = findall(<=(1e-9), dx)
    row("AC15b", "$lab: pairs with final states apart by more than 1e-9 in some coordinate, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(rmed(dx))), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    fc(j) = cpu.finals[j][1]
    fg(j) = gpu.finals[j][1]
    qr = sort([Float64(fc(j).RLI) / Float64(fg(j).RLI) for j in 1:m])
    row("AC15b", "$lab: distribution of RLI_cpu / RLI_gpu over all $m pairs: min / 1% / median / 99% / max (spec 15(b), reported, not gated)", join(fmt.((qr[1], qr[max(1, round(Int, 0.01m))], rmed(qr), qr[round(Int, 0.99m)], qr[end])), " / "), "reported", "INFO")
    if chaotic
        partb || return
        row("AC15b", "$lab: pairs with final states within 1e-9 in every coordinate (the regular-class selection; not the entering rule for this class)", length(agree), "reported", "INFO")
        return rli_chaotic_b(lab, cpu, gpu, rc, rg, orb, m)
    end
    row("AC15b", "$lab: pairs that enter part (b) (final states within 1e-9 in every coordinate)", length(agree), "reported", "INFO")
    g = [abs(Float64(fc(j).RLI) - Float64(fg(j).RLI)) / (rc[j].BR + rg[j].BR) for j in agree]
    row("AC15b", "$lab: pairs that agree: max |RLI,cpu - RLI,gpu| / (B_R,cpu + B_R,gpu) over $(length(agree)) pairs (median $(fmt(rmed(g))))", isempty(g) ? NaN : maximum(g), "<= 1", pf(!isempty(g) && maximum(g) <= 1))
    rel = [abs(Float64(fc(j).RLI) - Float64(fg(j).RLI)) / Float64(fc(j).RLI) for j in agree]
    row("AC15b", "$lab: pairs that agree: relative gap |RLI,cpu - RLI,gpu| / RLI,cpu median / max (reported)", isempty(rel) ? "-" : "$(fmt(rmed(rel))) / $(fmt(maximum(rel)))", "reported", "INFO")
end

# 15(b) statistical form for the chaotic class: K = 10, n = 1e4, tau = 2, 1024 pairs, tangent-core AC12(c) tolerances on RLI(n).
function rli_ac12c(T)
    label = "(c) K=10 chaotic ensemble, $T: RLI(n), 6 sigma_FT per pair and 4 SE for the means"
    reason = gpu_reason(T)
    isempty(reason) || return row("AC15c", label, "-", "6 sigma; 4 SE", fb(reason))
    rng = StableRNG(SEED + 50000 + 1003)
    x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
    S = r15_setup(T, exact(T, 10), x0, 1e-12)
    cpu = gpu = nothing
    e = catchhe(() -> (cpu = rli_est(S, 10^4, 2, :cpu); gpu = rli_est(S, 10^4, 2, :gpu; chunk = AC13_GPU_CHUNK)))
    e === nothing || return row("AC15c", "(c) K=10, $T, RLI: health error", sprint(showerror, e), "none", "FAIL")
    why = gpu_why(gpu)
    why === nothing || return row("AC15c", label, "-", "6 sigma; 4 SE", fb(why))
    sc, sg = [Float64(f[1].RLI) for f in cpu.finals], [Float64(f[1].RLI) for f in gpu.finals]
    sig = std(sc)
    maxd = maximum(abs.(sc .- sg)) / sig
    dm = abs(mean(sc) - mean(sg)) / sqrt((var(sc) + var(sg)) / length(sc))
    row("AC15c", "(c) K=10, $T, 1024 pairs, n=1e4, tau=2, Dx=1e-12, RLI: max per-pair |cpu - gpu| / sigma_FT (sigma_FT=$(fmt(sig))); |difference of means| / SE; means cpu $(fmt(mean(sc))) gpu $(fmt(mean(sg)))", "$(fmt(maxd)); $(fmt(dm))", "6; 4", pf(maxd <= 6 && dm <= 4))
end

# Kernel compile record and R7 (GPU): time per orbit-step of the pair run against a single-orbit L1 (k = 1) and MEGNO (k = 2) run, K = 2 pi, wrapped, tau = 1.
function rli_timing()
    stsz(r, T, d) = sizeof(typeof(Chaoslab.reducer_init(r, T, Val(d), Val(2))))
    row("AC15", "reducer state size in bytes: RLI Float64 / Double64 (d = 4, k = 2)", "$(stsz(RLI(), Float64, 4)) / $(stsz(RLI(), Double64, 4))", "reported", "INFO")
    norb, N = 20_000, 1_000
    rng = StableRNG(SEED + 15)
    x0 = [SVector{2,Float64}(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:norb]
    S = r15_setup(Float64, 2pi, x0, 1e-9)
    xw = [torus_state(x) for x in x0]
    cases = (("RLI pair (k = 2, d = 4)", (b, c) -> propagate(S.P, S.z0, N; k = 2, tau = 1, seed = SEED, reducers = (RLI(),), W0 = S.W0, backend = b, chunk = c)),
             ("L1 single (k = 1)", (b, c) -> propagate(S.base, xw, N; k = 1, tau = 1, seed = SEED, reducers = (Benettin(),), backend = b, chunk = c)),
             ("MEGNO single (k = 2)", (b, c) -> propagate(S.base, xw, N; k = 2, tau = 1, seed = SEED, reducers = (MEGNO(),), backend = b, chunk = c)))
    for (nm, run) in cases
        cpu_t = gpu_t = Inf
        prov = nothing
        try
            run(:gpu, AC13_GPU_CHUNK)       # warm-up (compilation of the kernel on the full chunk shape)
            for _ in 1:3
                t = @elapsed res = run(:gpu, AC13_GPU_CHUNK)
                prov = res.provenance
                gpu_t = min(gpu_t, t)
            end
        catch e
            row("AC15", "$nm kernel on the device: compile / run error", sprint(showerror, e)[1:min(end, 300)], "compiles", "FAIL")
            continue
        end
        run(:cpu, nothing)
        for _ in 1:2
            cpu_t = min(cpu_t, @elapsed run(:cpu, nothing))
        end
        steps = norb * N
        row("AC15", "$nm kernel (provenance backend $(prov.backend)); $norb orbits x $N steps, K = 2 pi, wrapped, tau = 1: GPU $(fmt(gpu_t / steps * 1e9)) ns/orbit-step, CPU ($(Threads.nthreads()) threads) $(fmt(cpu_t / steps * 1e9)) ns/orbit-step, speedup $(fmt(cpu_t / gpu_t))",
            prov.backend === :gpu ? gpu_t / steps * 1e9 : NaN, "recorded", prov.backend === :gpu ? "INFO" : fb(prov.backend_fallback_reason))
    end
end

want("rli") && let
    rli_timing()
    x0, ncand, nrep = regular_starts(1024)
    N = 4 * 10^4
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=0.5 regular class, $T, n=$(fmt(N)): each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        rli_d15("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, exact(T, 0.5), xT, N)
    end
    for T in (Float64, Double64)
        rli_ac12c(T)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=10 chaotic class, $T: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
        for nch in (8, 10)
            rli_d15("K=10 chaotic class", T, exact(T, 10), xT, nch; chaotic = true)
        end
        rli_d15("K=10 chaotic class", T, exact(T, 10), xT, 100; chaotic = true, partb = false)
    end
end

# ------------------------------------------------------------------------------------------------ AC16 (stretching-numbers), tag "stretching"
# Spec criterion 16 of stretching-numbers-spec.md, the D68 form (the bounds of SN14 and SN15 apply per backend; the GPU and the CPU do not produce the same orbit, so each backend is compared
# with a BigFloat(256) reference formed at its own realized points). StretchingNumbers reducer (k = 2, chain {StretchingNumbers()}, the replacer) on WrappedStandardMap, seeded frames passed
# as W0, tau = 1, B = 1000 over the derived range, interval 100; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64 if the kernel does not compile). Classes: regular
# K = 0.5 (S7 screening, regular_starts), 1024 orbits, n = 4e4; chaotic K = 10, 1024 orbits, n = 1e4.
# (a) each backend against its own reference: at 65 handoffs g_k = round(n k / 64), k = 1..64, and g = 1 (one propagate call per g, the orbit is a prefix) every cell count lies in its
# envelope N-_c <= N_c <= N+_c of SN14 (Psi_max over the grid 1..9 and 10 per decade, from the backend's own tangents; orbits with max half-width >= h/2 are counted unresolved and not gated)
# and |D - D_ref| <= sqrt 2 (U_1 + U_2) / (g sqrt h) + 4 eps D_ref (SN15). (b) CPU-GPU gap: regular class, on the orbits whose final states differ by at most 1e-9 in every coordinate, at
# every handoff of the grid |N_cpu,c - N_gpu,c| <= (N+ - N-)_cpu,c + (N+ - N-)_gpu,c in every cell of both columns and |D_cpu - D_gpu| <= the SN15 bounds of both backends added; the number of
# orbits with any differing count at n and the largest per-cell difference, and the distance of the two pooled spectra (SN13) against the same-backend split-half distance, are reported, not
# gated. The chaotic class is statistical (spec 16(b)): K = 10, n = 1e4, tau = 2, 1024 orbits, the mean of z (column 1 and 2) with the tangent-core AC12(c) tolerances (6 sigma per orbit,
# 4 SE for the means). The GPU orbit is read with the trace kernel of the rli section (x coordinate; the tangent of the wrapped map depends on x only); its final state is checked bitwise
# against the propagate result. GPU_SN_M, GPU_SN_N, GPU_SNC_N shrink the orbit count and the run lengths (debugging aid; the default is the spec size).
MODTAG[] = "stretching-numbers: "
want("stretching") && include(joinpath(ROOT, "test", "stretching", "sn_helpers.jl"))   # sn_*: BigFloat reference sn_Ref, edge certificate sn_envelope, sn_halfwidths, sn_D_ref
sn_gest(sys, x0, W0, N, tau, b, B; chunk = nothing) = propagate(sys, x0, N; k = 2, tau = tau, seed = SEED, reducers = (sn_red(sys, B, tau),), W0 = W0, backend = b, chunk = chunk)

# The reference orbit of one backend from its realized x coordinates (steps 1..N, the state before step i) and the seeded columns.
function sn_g_ref(Kt, xs, W)
    N = length(xs)
    w1, w2 = SVector{2,Float64}(W[1, 1], W[2, 1]), SVector{2,Float64}(W[1, 2], W[2, 2])
    mats = sn_mats_std(Kt, xs)
    return sn_Ref(sn_ref_s(mats, w1, N), sn_ref_s(mats, w2, N), r15_Ms(Kt, xs), w1, w2, Dict{Tuple{Int,Int},Float64}(), Dict{Int,Vector{Float64}}())
end

# Envelope, bound and ratios of orbit j of one backend (res_g: g -> propagate result of the prefix run). Returns the envelopes for the pairwise comparison.
function sn_g_assess(res_g, grid, j, ref, lo, hi, B, h, epsT)
    N = length(ref.s1)
    zs = [sn_z_from_s(ref.s1, 1), sn_z_from_s(ref.s2, 1)]
    eps = [sn_halfwidths(ref, i, epsT, 1, N, lo, hi) for i in 1:2]
    resolved = all(maximum(eps[i]) < Float64(h) / 2 for i in 1:2)
    env = [sn_envelope(zs[i], eps[i], lo, hi, B, grid) for i in 1:2]
    rC = rD = 0.0
    ndiff = maxU = 0
    bnd = Dict{Int,Float64}()
    for g in grid
        f = res_g[g].finals[j][1]
        U = env[1][g].U + env[2][g].U
        maxU = max(maxU, U)
        Dref = sn_D_ref(env[1][g].N, env[2][g].N, lo, hi, B)
        bound = setprecision(BigFloat, 256) do
            sqrt(BigFloat(2)) * U / (g * sqrt(h)) + 4 * epsT * Dref
        end
        bnd[g] = Float64(bound)
        resolved || continue
        for i in 1:2
            cnt = sn_cnt(f, i)
            for c in eachindex(cnt)
                rC = max(rC, cnt[c] >= env[i][g].Nm[c] && cnt[c] <= env[i][g].Np[c] ? 0.0 : Inf)
            end
        end
        ndiff += any(sn_cnt(f, i) != env[i][g].N for i in 1:2)
        eD = setprecision(BigFloat, 256) do
            abs(sn_big(f.D) - Dref)
        end
        rD = max(rD, bound == 0 ? (eD == 0 ? 0.0 : Inf) : Float64(eD / bound))
    end
    return (env = env, resolved = resolved, rC = rC, rD = rD, ndiff = ndiff, maxU = maxU, bnd = bnd)
end

function sn_d16(tag, T, Kt, xT, N; chaotic = false)
    m = length(xT)
    lab = "$tag, $T, $m orbits, n=$(fmt(N)), tau=1, B=1000"
    base = WrappedStandardMap(exact(T, Kt))
    xW = [torus_state(x) for x in xT]
    Ws = [sn_frame(SEED, j, T) for j in 1:m]
    B = 1000
    lo, hi = stretching_range(base)
    h = setprecision(BigFloat, 256) do
        (sn_big(hi) - sn_big(lo)) / B
    end
    grid = m15_grid(N)
    gres = Dict{Int,Any}()
    gres[N] = sn_gest(base, xW, Ws, N, 1, :gpu, B; chunk = AC13_GPU_CHUNK)
    why = gpu_why(gres[N])
    why === nothing || return row("AC16a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    gpu = gres[N]
    cpu = sn_gest(base, xW, Ws, N, 1, :cpu, B)
    cres = Dict{Int,Any}(N => cpu)
    for g in grid
        g == N && continue
        gres[g] = sn_gest(base, xW, Ws, g, 1, :gpu, B; chunk = AC13_GPU_CHUNK)
        cres[g] = sn_gest(base, xW, Ws, g, 1, :cpu, B)
    end
    epsT = Float64(effective_eps(T))
    PA, XN = r15_trace(base, xW, N, r15_x)
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = Set(findall(<=(1e-9), dx))
    rc = Vector{Any}(undef, m); rg = Vector{Any}(undef, m); gap = Vector{Any}(undef, m)
    same_c = zeros(Bool, m)
    Threads.@threads :dynamic for j in 1:m
        pts = sn_pts(base, xW[j], N)
        same_c[j] = forward(base, pts[end]) == cpu.x[j]
        refc = sn_g_ref(Kt, [p[1] for p in pts], Ws[j])
        refg = sn_g_ref(Kt, PA[j, :], Ws[j])
        ac = sn_g_assess(cres, grid, j, refc, lo, hi, B, h, epsT)
        ag = sn_g_assess(gres, grid, j, refg, lo, hi, B, h, epsT)
        gap[j] = nothing
        if !chaotic && j in agree && ac.resolved && ag.resolved
            rN = rDg = 0.0
            for g in grid
                fc, fg = cres[g].finals[j][1], gres[g].finals[j][1]
                for i in 1:2
                    ec, eg = ac.env[i][g], ag.env[i][g]
                    cc, cg = sn_cnt(fc, i), sn_cnt(fg, i)
                    for c in eachindex(ec.N)
                        d = abs(cc[c] - cg[c])
                        w = (ec.Np[c] - ec.Nm[c]) + (eg.Np[c] - eg.Nm[c])
                        rN = max(rN, d == 0 ? 0.0 : (w == 0 ? Inf : d / w))
                    end
                end
                dD = abs(Float64(fc.D) - Float64(fg.D))
                bb = ac.bnd[g] + ag.bnd[g]
                rDg = max(rDg, dD == 0 ? 0.0 : (bb == 0 ? Inf : dD / bb))
            end
            gap[j] = (rN = rN, rD = rDg)
        end
        rc[j] = (resolved = ac.resolved, rC = ac.rC, rD = ac.rD, ndiff = ac.ndiff, maxU = ac.maxU)
        rg[j] = (resolved = ag.resolved, rC = ag.rC, rD = ag.rD, ndiff = ag.ndiff, maxU = ag.maxU)
    end
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    row("AC16a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        bad = count(x -> x.rC > 0, R)
        row("AC16a", "$lab: backend $nm (provenance $(res.provenance.backend)): orbits with a cell count outside its envelope N- .. N+ at some of the $(length(grid)) handoffs (of $m; resolved orbits only)", bad, "0", pf(bad == 0))
        row("AC16a", "$lab: backend $nm: max |D - D_ref| / (sqrt 2 (U1 + U2) / (g sqrt h) + 4 eps D_ref) over $m orbits and the grid handoffs (median over orbits $(fmt(rmed([x.rD for x in R]))))", maximum(x.rD for x in R), "<= 1", pf(maximum(x.rD for x in R) <= 1))
        row("AC16a", "$lab: backend $nm: unresolved orbits (max half-width >= h/2); orbit-handoffs with counts different from the reference; largest U1 + U2", "$(count(x -> !x.resolved, R)); $(sum(x.ndiff for x in R)); $(maximum(x.maxU for x in R))", "reported", "INFO")
    end
    row("AC16b", "$lab: orbits with final states apart by more than 1e-9 in some coordinate, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(rmed(dx))), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    chaotic && return
    G = [g for g in gap if g !== nothing]
    row("AC16b", "$lab: orbits that enter part (b) (final states within 1e-9 in every coordinate, both backends resolved)", length(G), "reported", "INFO")
    if !isempty(G)
        row("AC16b", "$lab: orbits that agree: max over cells, both columns and the grid handoffs of |N_cpu,c - N_gpu,c| / ((N+ - N-)_cpu,c + (N+ - N-)_gpu,c) ($(length(G)) orbits)", maximum(g.rN for g in G), "<= 1", pf(maximum(g.rN for g in G) <= 1))
        row("AC16b", "$lab: orbits that agree: max |D_cpu - D_gpu| / (SN15 bound of cpu + that of gpu)", maximum(g.rD for g in G), "<= 1", pf(maximum(g.rD for g in G) <= 1))
    end
    # reported: count differences at n, pooled spectra (SN13) and the split-half noise reference
    cn(res, j, i) = Int.(res.finals[j][1].counts[:, i])
    dmax = [maximum(abs.(cn(cpu, j, 1) .- cn(gpu, j, 1))) for j in 1:m]
    dmax2 = [maximum(abs.(cn(cpu, j, 2) .- cn(gpu, j, 2))) for j in 1:m]
    row("AC16b", "$lab: reported, not gated: orbits with any differing count at n (column 1; column 2); largest per-cell count difference (column 1; column 2)", "$(count(>(0), dmax)); $(count(>(0), dmax2)); $(maximum(dmax)); $(maximum(dmax2))", "reported", "INFO")
    pool(res, idx, i) = sum(cn(res, j, i) for j in idx)
    hT = (Float64(hi) - Float64(lo)) / B
    Dxb = spectral_distance(pool(cpu, 1:m, 1), pool(gpu, 1:m, 1), hT)
    Dcc = spectral_distance(pool(cpu, 1:2:m, 1), pool(cpu, 2:2:m, 1), hT)
    Dgg = spectral_distance(pool(gpu, 1:2:m, 1), pool(gpu, 2:2:m, 1), hT)
    row("AC16b", "$lab: reported, not gated: distance D of the two pooled spectra (all $m orbits, column 1) CPU against GPU; same-backend split-half distance (odd against even orbits) CPU, GPU", "$(fmt(Dxb)); $(fmt(Dcc)); $(fmt(Dgg))", "reported", "INFO")
end

# 16(b) statistical form for the chaotic class: K = 10, n = 1e4, tau = 2, 1024 orbits, tangent-core AC12(c) tolerances on the mean of z (both columns).
function sn_stat(T)
    label = "(b) K=10 chaotic ensemble, $T: mean of z, 6 sigma_FT per orbit and 4 SE for the means"
    reason = gpu_reason(T)
    isempty(reason) || return row("AC16c", label, "-", "6 sigma; 4 SE", fb(reason))
    rng = StableRNG(SEED + 50000 + 1003)
    x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:1024]
    base = WrappedStandardMap(exact(T, 10))
    xW = [torus_state(x) for x in x0]
    Ws = [sn_frame(SEED, j, T) for j in 1:1024]
    cpu = gpu = nothing
    e = catchhe(() -> (cpu = sn_gest(base, xW, Ws, 10^4, 2, :cpu, 1000); gpu = sn_gest(base, xW, Ws, 10^4, 2, :gpu, 1000; chunk = AC13_GPU_CHUNK)))
    e === nothing || return row("AC16c", "(b) K=10, $T, stretching: health error", sprint(showerror, e), "none", "FAIL")
    why = gpu_why(gpu)
    why === nothing || return row("AC16c", label, "-", "6 sigma; 4 SE", fb(why))
    for col in 1:2
        sc, sg = [Float64(f[1].mean[col]) for f in cpu.finals], [Float64(f[1].mean[col]) for f in gpu.finals]
        sig = std(sc)
        maxd = maximum(abs.(sc .- sg)) / sig
        dm = abs(mean(sc) - mean(sg)) / sqrt((var(sc) + var(sg)) / length(sc))
        row("AC16c", "(b) K=10, $T, 1024 orbits, n=1e4, tau=2, mean of z column $col: max per-orbit |cpu - gpu| / sigma_FT (sigma_FT=$(fmt(sig))); |difference of means| / SE; means cpu $(fmt(mean(sc))) gpu $(fmt(mean(sg)))", "$(fmt(maxd)); $(fmt(dm))", "6; 4", pf(maxd <= 6 && dm <= 4))
    end
    Dc = Float64[Float64(f[1].D) for f in cpu.finals]
    Dg = Float64[Float64(f[1].D) for f in gpu.finals]
    row("AC16c", "(b) K=10, $T, D at n = 1e4 (tau = 2, 5000 handoffs, under-sampled): median cpu / gpu; reported", "$(fmt(rmed(Dc))) / $(fmt(rmed(Dg)))", "reported", "INFO")
end

# Kernel compile record and R6 (GPU): time per orbit-step of the two-column stretching run against a single-orbit L1 (k = 1) and MEGNO (k = 2) run, K = 2 pi, wrapped, tau = 1, B = 100 and 1000.
function sn_timing()
    stsz(r, T) = sizeof(typeof(Chaoslab.reducer_init(r, T, Val(2), Val(2))))
    sys0 = WrappedStandardMap(2pi)
    mk(B) = StretchingNumbers(; range = stretching_range(sys0), bins = B, interval = 100)
    row("AC16", "reducer state size in bytes (without the workspace): StretchingNumbers Float64 / Double64; workspace bytes per orbit at B = 100 / 1000 (4 (B + 2) per column, two columns)", "$(stsz(mk(100), Float64)) / $(stsz(mk(100), Double64)); $(8 * 102) / $(8 * 1002)", "reported", "INFO")
    norb, N = 20_000, 1_000
    rng = StableRNG(SEED + 16)
    x0 = [SVector{2,Float64}(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:norb]
    xw = [torus_state(x) for x in x0]
    cases = (("StretchingNumbers B = 100 (k = 2)", (b, c) -> propagate(sys0, xw, N; k = 2, tau = 1, seed = SEED, reducers = (mk(100),), backend = b, chunk = c)),
             ("StretchingNumbers B = 1000 (k = 2)", (b, c) -> propagate(sys0, xw, N; k = 2, tau = 1, seed = SEED, reducers = (mk(1000),), backend = b, chunk = c)),
             ("L1 single (k = 1)", (b, c) -> propagate(sys0, xw, N; k = 1, tau = 1, seed = SEED, reducers = (Benettin(),), backend = b, chunk = c)),
             ("MEGNO single (k = 2)", (b, c) -> propagate(sys0, xw, N; k = 2, tau = 1, seed = SEED, reducers = (MEGNO(),), backend = b, chunk = c)))
    for (nm, run) in cases
        cpu_t = gpu_t = Inf
        prov = nothing
        try
            run(:gpu, AC13_GPU_CHUNK)       # warm-up (compilation of the kernel on the full chunk shape)
            for _ in 1:3
                t = @elapsed res = run(:gpu, AC13_GPU_CHUNK)
                prov = res.provenance
                gpu_t = min(gpu_t, t)
            end
        catch e
            row("AC16", "$nm kernel on the device: compile / run error", sprint(showerror, e)[1:min(end, 1500)], "compiles", "FAIL")
            continue
        end
        run(:cpu, nothing)
        for _ in 1:2
            cpu_t = min(cpu_t, @elapsed run(:cpu, nothing))
        end
        steps = norb * N
        row("AC16", "$nm kernel (provenance backend $(prov.backend)); $norb orbits x $N steps, K = 2 pi, wrapped, tau = 1: GPU $(fmt(gpu_t / steps * 1e9)) ns/orbit-step, CPU ($(Threads.nthreads()) threads) $(fmt(cpu_t / steps * 1e9)) ns/orbit-step, speedup $(fmt(cpu_t / gpu_t))",
            prov.backend === :gpu ? gpu_t / steps * 1e9 : NaN, "recorded", prov.backend === :gpu ? "INFO" : fb(prov.backend_fallback_reason))
    end
end

want("stretching") && let
    sn_timing()
    nreg = parse(Int, get(ENV, "GPU_SN_M", "1024"))
    x0, ncand, nrep = regular_starts(nreg)
    N = parse(Int, get(ENV, "GPU_SN_N", string(4 * 10^4)))
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC16a", "K=0.5 regular class, $T, n=$(fmt(N)): each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        sn_d16("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, 0.5, xT, N)
    end
    for T in (Float64, Double64)
        sn_stat(T)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC16a", "K=10 chaotic class, $T: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:nreg]
        sn_d16("K=10 chaotic class", T, 10.0, xT, parse(Int, get(ENV, "GPU_SNC_N", string(10^4))); chaotic = true)
    end
end

# ------------------------------------------------------------------------------------------------ AC15 (naff-fma), tag "naff"
# Spec criterion 15 of naff-fma-spec.md, the D68 form (the bound B_E2 of criterion 10(d) applies per backend; the GPU and the CPU do not produce the same orbit, so each backend is compared with a BigFloat
# reference formed at its own realized samples). WeightedPhaseAdvance reader (E2: window N = 1000, LiftedAngleObservable(1, 3), exponential weight) chained before Benettin() (the replacer that renormalises
# the one column; k = 1) on WrappedStandardMap, tau = 1, 2N + 1 = 2001 handoffs; Float64 on the GPU (a FALLBACK with its reason is the pass for Double64 if the kernel does not compile; a compiler error of
# the Float64 kernel is a FAIL and a finding). Classes: regular K = 0.5 (S7 screening, regular_starts), 1024 orbits; chaotic K = 10, 1024 orbits. (a) each backend against NF6 evaluated at 512 bits at that
# backend's own realized lifted increments (the GPU orbit is read with the trace kernel of the rli section, full state, its final state checked bitwise against the propagate result): every orbit's two
# frequencies within B_E2 = 8 eps (1 + max |Delta_j| / (2 pi)). (b) CPU-GPU gap only on the orbits whose final states differ by at most 1e-9 in every coordinate: |nu_cpu - nu_gpu| <= B_cpu + B_gpu; their
# number and the number of other orbits are reported. (c) The series path: the E1 transform (frequency_map, window order 1, exp(i q) with the window means removed in BigFloat and rounded once) applied to the
# series recorded on the device (1024 regular-class orbits, N = 1000, Float64) against the BigFloat(256) stationary point of that same series (nf_compare_window: 16 eps; coarse bin or tie); the device series
# against the CPU series is reported (bitwise equal orbits, orbits within 1e-9 at the last sample). GPU_NF_M, GPU_NF_N shrink the orbit count and N (debugging aid; the default is the spec size).
MODTAG[] = "naff-fma: "
want("naff") && include(joinpath(ROOT, "test", "naff-fma", "nf_helpers.jl"))   # nf_*: BigFloat reference nf_E2_ref / nf_stat_ref / nf_compare_window, nf_lifted_incr_big, nf_BE2
naff_red(N) = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
naff_est(sys, x0, N, b; chunk = nothing) = propagate(sys, x0, 2N + 1; k = 1, tau = 1, seed = SEED, reducers = (naff_red(N), Benettin()), backend = b, chunk = chunk)
naff_cis(z) = complex(cos(z[1]), sin(z[1]))
# Device states at the handoffs 1..n (after steps 1..n): the trace gives the state before step i, so handoff h is column h + 1 and the last one is the final state.
function naff_dev_states(sys, x0, n)
    P, XN = r15_trace(sys, x0, n, r15_id)
    return [(v = P[j, 2:n]; push!(v, XN[j]); v) for j in eachindex(x0)], XN
end

function naff_d15(tag, T, Kt, xT, N)
    m = length(xT)
    n = 2N + 1
    lab = "$tag, $T, $m orbits, N=$N ($n handoffs), tau=1"
    sys = WrappedStandardMap(exact(T, Kt))
    xW = [torus_state(x) for x in xT]
    gpu = try
        naff_est(sys, xW, N, :gpu; chunk = AC13_GPU_CHUNK)
    catch e
        return row("AC15a", "$lab: GPU kernel of the E2 reader chain compiles and runs", sprint(showerror, e)[1:min(end, 3000)], "compiles", "FAIL")
    end
    why = gpu_why(gpu)
    why === nothing || return row("AC15a", "$lab: each backend against its own BigFloat reference", "-", "<= B_E2", fb(why))
    cpu = naff_est(sys, xW, N, :cpu)
    Sg, XN = naff_dev_states(sys, xW, n)
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    wh = setprecision(BigFloat, 512) do
        w = [nf_wb(BigFloat(j) / N) for j in 0:(N - 1)]
        w ./ sum(w)
    end
    function ref(states)
        incr = nf_lifted_incr_big(states)
        nus = setprecision(BigFloat, 512) do
            tp = 2 * BigFloat(pi)
            (sum(wh[k] * incr[k] for k in 1:N) / tp, sum(wh[k] * incr[N + k] for k in 1:N) / tp)
        end
        return nus, nf_BE2(T, incr)
    end
    rc = Vector{Any}(undef, m); rg = Vector{Any}(undef, m)
    same_c = zeros(Bool, m)
    Threads.@threads :dynamic for j in 1:m
        sc = nf_states(sys, xW[j], n)
        same_c[j] = sc[end] == cpu.x[j]
        for (R, st, res) in ((rc, sc, cpu), (rg, Sg[j], gpu))
            (rf, B) = ref(st)
            f = res.finals[j][1]
            nu = (nf_big(f.nu1[1]), nf_big(f.nu2[1]))
            e = [Float64(abs(nu[1] - rf[1])), Float64(abs(nu[2] - rf[2]))]
            R[j] = (e = e, B = B, nu = nu, ratio = maximum(e) / B)
        end
    end
    row("AC15a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        bad = count(x -> x.ratio > 1, R)
        row("AC15a", "$lab: backend $nm (provenance $(res.provenance.backend)): orbits with a frequency (either window) outside B_E2 of its own BigFloat evaluation of NF6 at the backend's realized samples (of $m)", bad, "0", pf(bad == 0))
        row("AC15a", "$lab: backend $nm: max |nu - nu_ref| / B_E2 over $m orbits and both windows (median over orbits $(fmt(rmed([x.ratio for x in R]))))", maximum(x.ratio for x in R), "<= 1", pf(maximum(x.ratio for x in R) <= 1))
    end
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    agree = findall(<=(1e-9), dx)
    row("AC15b", "$lab: orbits with final states apart by more than 1e-9 in some coordinate, CPU vs GPU (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(rmed(dx))), max $(fmt(maximum(dx))))", count(>(1e-9), dx), "reported", "INFO")
    row("AC15b", "$lab: orbits that enter part (b) (final states within 1e-9 in every coordinate)", length(agree), "reported", "INFO")
    if !isempty(agree)
        g = [maximum(Float64(abs(rc[j].nu[k] - rg[j].nu[k])) for k in 1:2) / (rc[j].B + rg[j].B) for j in agree]
        row("AC15b", "$lab: orbits that agree: max |nu_cpu - nu_gpu| / (B_E2,cpu + B_E2,gpu) over both windows ($(length(agree)) orbits; median $(fmt(rmed(g))))", maximum(g), "<= 1", pf(maximum(g) <= 1))
    end
    allgap = [maximum(Float64(abs(rc[j].nu[k] - rg[j].nu[k])) for k in 1:2) for j in 1:m]
    # diagnosis of the gap (reported): the two backends read different orbits (final states apart by dx), so the gap contains the frequency change between nearby orbits, which B_E2 does not bound;
    # the algebra alone is gated on IDENTICAL samples: the device reducer against the pure function weighted_phase_frequency (CPU) on the device samples, |difference| <= B_E2,gpu + B_E2,pure.
    if !isempty(agree)
        gd = [allgap[j] / dx[j] for j in agree if dx[j] > 0]
        tight = [j for j in agree if dx[j] <= 1e-11]
        row("AC15b", "$lab: reported: max gap / final-state distance dx over the $(length(agree)) agreeing orbits (the frequency sensitivity to the orbit difference); orbits with dx <= 1e-11 and the max gap ratio to B_cpu + B_gpu among them", "$(fmt(maximum(gd))); $(length(tight)); $(isempty(tight) ? NaN : fmt(maximum(allgap[j] / (rc[j].B + rg[j].B) for j in tight)))", "reported", "INFO")
    end
    ss = Vector{Float64}(undef, m)
    Threads.@threads :dynamic for j in 1:m
        pure = weighted_phase_frequency(nf_lifted_incr(Sg[j]), N)
        ss[j] = maximum(Float64(abs(nf_big(pure[k]) - rg[j].nu[k])) for k in 1:2) / (2 * rg[j].B)
    end
    row("AC15b", "$lab: same samples, algebra only: max over $m orbits and both windows of |nu_gpu reducer - nu_cpu pure function on the device samples| / (2 B_E2) (median $(fmt(rmed(ss))))", maximum(ss), "<= 1", pf(maximum(ss) <= 1))
    row("AC15b", "$lab: reported, not gated: |nu_cpu - nu_gpu| over all orbits and both windows: median / max (cycles per sample)", "$(fmt(rmed(allgap))) / $(fmt(maximum(allgap)))", "reported", "INFO")
end

# (c) the series path on the device: E1 on the exp(i q) series recorded by the trace kernel (Float64, regular class).
const NAFF_CNT_LOCK = ReentrantLock()
function naff_e1_device(tag, xT, N)
    m = length(xT)
    lab = "$tag, Float64, $m orbits, N=$N, window order 1, exp(i q) with the window means removed"
    sys = WrappedStandardMap(exact(Float64, 0.5))
    xW = [torus_state(x) for x in xT]
    Sg, XN = try
        naff_dev_states(sys, xW, 2N)
    catch e
        return row("AC15c", "$lab: device trace kernel compiles and runs", sprint(showerror, e)[1:min(end, 3000)], "compiles", "FAIL")
    end
    cnt = Dict(:ok => 0, :tie => 0, :nonunique => 0, :binfail => 0, :nufail => 0)
    errs = fill(NaN, m, 2)
    same = zeros(Bool, m); near = zeros(Bool, m)
    Threads.@threads :dynamic for j in 1:m
        sc = nf_states(sys, xW[j], 2N)
        same[j] = all(sc .== Sg[j]); near[j] = maximum(abs.(sc[end] - Sg[j][end])) <= 1e-9
        zb = Complex{BigFloat}[nf_big(naff_cis(s)) for s in Sg[j]]
        z = nf_series(Float64, nf_demean_big(zb, N))
        fm = frequency_map(z, N)
        for (k, (w, zw)) in enumerate(((fm.w1[1], z[1:N]), (fm.w2[1], z[(N + 1):(2N)])))
            s, e = nf_compare_window(Float64, zw, w, 1, :half, 8)
            errs[j, k] = e
            lock(NAFF_CNT_LOCK) do
                cnt[s] += 1
            end
        end
    end
    bad = cnt[:binfail] + cnt[:nufail]
    row("AC15c", "$lab: windows by status against the BigFloat(256) stationary point of the same series (ok / tie within 64 eps of Float64 / non-unique reference, listed / bin fail / frequency fail)", "$(cnt[:ok]) / $(cnt[:tie]) / $(cnt[:nonunique]) / $(cnt[:binfail]) / $(cnt[:nufail]) of $(2m)", "0 fail", pf(bad == 0))
    fin = filter(isfinite, vec(errs))
    row("AC15c", "$lab: max frequency error against the stationary point / 16 eps (windows with a unique reference, $(length(fin)))", isempty(fin) ? NaN : maximum(fin) / (16 * epsf(Float64)), "<= 1", pf(isempty(fin) || maximum(fin) <= 16 * epsf(Float64)))
    row("AC15c", "$lab: device series against the CPU series: orbits bitwise equal over all $(2N) samples / equal in the last sample within 1e-9 (reported, recorder agreement D54 not restated)", "$(count(same)) / $(count(near)) of $m", "reported", "INFO")
end

# R5 on the device: the E2 reader against the stream alone, and the orbit pass of the series path (device trace of exp(i q), read back to the host) per orbit, at N in {1e3, 1e4, 1e5}, K = 2 pi, wrapped.
function naff_timing()
    sys0 = WrappedStandardMap(2pi)
    for (N, norb, ntr) in ((1000, 20_000, 1024), (10_000, 4_000, 512), (100_000, 400, 64))
        n = 2N + 1
        rng = StableRNG(SEED + 15)
        xw = [torus_state(SVector{2,Float64}(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:norb]
        cases = (("stream alone (Benettin, k = 1)", (b, c) -> propagate(sys0, xw, n; k = 1, tau = 1, seed = SEED, reducers = (Benettin(),), backend = b, chunk = c)),
                 ("E2 reader + Benettin (k = 1)", (b, c) -> propagate(sys0, xw, n; k = 1, tau = 1, seed = SEED, reducers = (naff_red(N), Benettin()), backend = b, chunk = c)))
        for (nm, run) in cases
            cpu_t = gpu_t = Inf
            prov = nothing
            try
                run(:gpu, AC13_GPU_CHUNK)
                for _ in 1:3
                    t = @elapsed res = run(:gpu, AC13_GPU_CHUNK)
                    prov = res.provenance
                    gpu_t = min(gpu_t, t)
                end
            catch e
                row("AC15", "$nm, N=$N: kernel on the device compile / run error", sprint(showerror, e)[1:min(end, 1500)], "compiles", "FAIL")
                continue
            end
            run(:cpu, nothing)
            for _ in 1:2
                cpu_t = min(cpu_t, @elapsed run(:cpu, nothing))
            end
            steps = norb * n
            row("AC15", "R5 $nm (provenance backend $(prov.backend)); N=$N, $norb orbits x $n steps, K = 2 pi, wrapped, tau = 1: GPU $(fmt(gpu_t / steps * 1e9)) ns/orbit-step ($(fmt(gpu_t / norb)) s per orbit), CPU ($(Threads.nthreads()) threads) $(fmt(cpu_t / steps * 1e9)) ns/orbit-step, speedup $(fmt(cpu_t / gpu_t))",
                prov.backend === :gpu ? gpu_t / steps * 1e9 : NaN, "recorded", prov.backend === :gpu ? "INFO" : fb(prov.backend_fallback_reason))
        end
        xt = xw[1:ntr]
        try
            r15_trace(sys0, xt, 2N, naff_cis; chunk = 2000)
            t = minimum(@elapsed(r15_trace(sys0, xt, 2N, naff_cis; chunk = 2000)) for _ in 1:2)
            row("AC15", "R5 E1 orbit pass on the device (trace kernel writing exp(i q) at 2N = $(2N) handoffs, read back to the host): N=$N, $ntr orbits: seconds per orbit; ns per orbit-step", "$(fmt(t / ntr)); $(fmt(t / (ntr * 2N) * 1e9))", "recorded", "INFO")
        catch e
            row("AC15", "R5 E1 orbit pass on the device, N=$N: trace kernel error", sprint(showerror, e)[1:min(end, 1500)], "runs", "FAIL")
        end
    end
end

want("naff") && let
    naff_timing()
    nreg = parse(Int, get(ENV, "GPU_NF_M", "1024"))
    x0, ncand, nrep = regular_starts(nreg)
    N = parse(Int, get(ENV, "GPU_NF_N", "1000"))
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=0.5 regular class, $T, N=$N: each backend against its own BigFloat reference", "-", "<= B_E2", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        naff_d15("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", T, 0.5, xT, N)
    end
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC15a", "K=10 chaotic class, $T, N=$N: each backend against its own BigFloat reference", "-", "<= B_E2", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:nreg]
        naff_d15("K=10 chaotic class", T, 10.0, xT, N)
    end
    naff_e1_device("K=0.5 regular class (S7)", [SVector{2,Float64}(z[1], z[2]) for z in x0], N)
end

# ------------------------------------------------------------------------------------------ AC13 (weighted-birkhoff), tag "weighted-birkhoff"
# Spec criterion 13 of weighted-birkhoff-spec.md, the D68 form as amended by naff-fma D82 (each backend against a BigFloat(256) reference formed at its own realized samples; the CPU-GPU comparison gated
# only on IDENTICAL samples; the orbit-level gap reported). WeightedBirkhoff (window N = 1000, one horizon, p = 1) chained before Benettin() (k = 1) on WrappedStandardMap, tau = 1, 2N = 2000 handoffs,
# Float64 on the device (Double64: FALLBACK with the reason if the kernel does not compile; a compiler error of the Float64 kernel is a FAIL and a finding). Classes: K_std = 0.5 regular (S7 screening,
# regular_starts), K_std = 10 chaotic, 1024 orbits each. Observable sets: the four trigonometric observables cos q, sin q, cos p, sin p of the spec (set "trig4") and the five of the criteria with the lifted p
# (set "std5"). (a) each backend's windows within c_1 eps S = 4 eps S and Delta within 8 eps S of the BigFloat(256) evaluation of WB14 at that backend's own samples (the observables evaluated in BigFloat at the
# exactly widened wrapped states; the device states from the trace kernel, checked bitwise against the propagate result). (b) the device reducer against the CPU pure function on the SAME samples,
# two readings: (b1) the samples the device itself computed (a trace kernel writing the K observable values at each handoff, read back: the exact input of the device reducer) and (b2) the CPU evaluation of the
# observables at the device states; gate |avg difference| <= 4 eps S, |Delta difference| <= 8 eps S for every orbit. (c) the orbit-level CPU-GPU gap of the averages is reported. (d) Double64. GPU_WB_M and
# GPU_WB_N shrink the orbit count and N (debugging aid; the default is the spec size).
MODTAG[] = "weighted-birkhoff: "
want("weighted-birkhoff") && include(joinpath(ROOT, "test", "weighted-birkhoff", "wb_helpers.jl"))   # wb_*: BigFloat reference wb_ref_avgs, observable columns wb_fref_std, wb_std5, wb_series
wbg_trig4() = (CosObservable(1), SinObservable(1), CosObservable(2), SinObservable(2))
struct WbProj{O}     # the K observable values of a state, evaluated by the device arithmetic (the exact input of the reducer)
    obs::O
end
(p::WbProj)(z) = map(o -> Float64(o(z)), p.obs)
wbg_red(N, obs; order = 1) = WeightedBirkhoff(; window = N, observables = obs, birkhoff_order = order)
wbg_est(sys, x0, N, b, obs; chunk = nothing) = propagate(sys, x0, 2N; k = 1, tau = 1, seed = SEED, reducers = (wbg_red(N, obs), Benettin()), backend = b, chunk = chunk)

function wbg_d13(tag, name, T, Kt, xT, N, obs, ncol)
    m = length(xT)
    K = length(obs)
    lab = "$tag, $name (K = $K), $T, $m orbits, N=$N ($(2N) handoffs), tau=1, p=1"
    sys = WrappedStandardMap(exact(T, Kt))
    xW = [torus_state(x) for x in xT]
    gpu = try
        wbg_est(sys, xW, N, :gpu, obs; chunk = AC13_GPU_CHUNK)
    catch e
        return row("AC13a", "$lab: GPU kernel of the reducer chain compiles and runs", sprint(showerror, e)[1:min(end, 3000)], "compiles", "FAIL")
    end
    why = gpu_why(gpu)
    why === nothing || return row("AC13a", "$lab: each backend against its own BigFloat reference", "-", "<= 1", fb(why))
    cpu = wbg_est(sys, xW, N, :cpu, obs)
    Sg, XN = naff_dev_states(sys, xW, 2N)                      # device states at the handoffs 1..2N (id trace) and the final state
    same_g = [XN[j] == gpu.x[j] for j in 1:m]
    # the observable values computed by the device at the handoffs 1..2N: trace of 2N + 1 steps, state before step i = handoff i - 1, columns 2..2N+1
    Pobs, _ = r15_trace(sys, xW, 2N + 1, WbProj(obs))
    epsT = epsf(T)
    wref = wb_wref(N, 1)                                       # fills the cache before the threaded loop
    rc = Vector{Any}(undef, m); rg = Vector{Any}(undef, m); b1 = Vector{Any}(undef, m); b2 = Vector{Any}(undef, m)
    same_c = zeros(Bool, m)
    Threads.@threads :dynamic for j in 1:m
        sc = wb_states(sys, xW[j], 2N)
        same_c[j] = sc[end] == cpu.x[j]
        for (R, st, res) in ((rc, sc, cpu), (rg, Sg[j], gpu))
            fb_ = wb_fref_std(st)[:, 1:ncol]
            f = res.finals[j][1]
            e1 = zeros(K); e2 = zeros(K); ed = zeros(K)
            for k in 1:K
                ref = wb_ref_avgs(fb_[:, k], N, 1)
                S = Float64(ref.S)
                e1[k] = wb_err(f.avg1[k], ref.avg1) / (4 * epsT * S)
                e2[k] = wb_err(f.avg2[k], ref.avg2) / (4 * epsT * S)
                ed[k] = setprecision(BigFloat, 256) do
                    Float64(abs(wb_big(f.delta[k]) - (ref.avg1 - ref.avg2))) / (8 * epsT * S)
                end
            end
            R[j] = (a = max(maximum(e1), maximum(e2)), d = maximum(ed))
        end
        # (b1): device-computed samples; (b2): CPU observables at the device states
        samp1 = Matrix{T}(undef, 2N, K)
        for i in 1:2N, k in 1:K
            samp1[i, k] = Pobs[j, i + 1][k]
        end
        samp2 = wb_series(obs, Sg[j])
        f = gpu.finals[j][1]
        for (R, samp) in ((b1, samp1), (b2, samp2))
            g = weighted_birkhoff(samp, N)
            ra = 0.0; rd = 0.0
            for k in 1:K
                S = Float64(g.S[k])
                ra = max(ra, abs(f.avg1[k] - g.avg1[k]) / (4 * epsT * S), abs(f.avg2[k] - g.avg2[k]) / (4 * epsT * S))
                rd = max(rd, abs(f.delta[k] - (g.avg1[k] - g.avg2[k])) / (8 * epsT * S))
            end
            R[j] = (a = ra, d = rd)
        end
    end
    row("AC13a", "$lab: the reference orbit reproduces the backend's final state bitwise, CPU / GPU (GPU orbit from the device trace kernel)", "$(count(same_c)) / $(count(same_g)) of $m", "$m / $m", "INFO")
    for (nm, R, res) in (("cpu", rc, cpu), ("gpu", rg, gpu))
        ma, md = maximum(x.a for x in R), maximum(x.d for x in R)
        row("AC13a", "$lab: backend $nm (provenance $(res.provenance.backend)): max over $m orbits, observables and both windows of |avg - avg_ref| / (4 eps S) (median over orbits $(fmt(rmed([x.a for x in R]))))", ma, "<= 1", pf(ma <= 1))
        row("AC13a", "$lab: backend $nm: max over orbits and observables of |Delta - Delta_ref| / (8 eps S) (median $(fmt(rmed([x.d for x in R]))))", md, "<= 1", pf(md <= 1))
        row("AC13a", "$lab: backend $nm: realized max |avg - avg_ref| in eps S (the spec's dry run: device 0.52)", ma * 4, "reported", "INFO")
    end
    for (nm, R) in (("b1: the samples the device computed (trace of the observable values), read back", b1), ("b2: the CPU evaluation of the observables at the device states", b2))
        ma, md = maximum(x.a for x in R), maximum(x.d for x in R)
        row("AC13b", "$lab: same samples $nm: device reducer against the CPU function weighted_birkhoff, max |avg difference| / (4 eps S) over $m orbits, observables and windows (median over orbits $(fmt(rmed([x.a for x in R]))))", ma, "<= 1", pf(ma <= 1))
        row("AC13b", "$lab: same samples $nm: max |Delta difference| / (8 eps S) (median $(fmt(rmed([x.d for x in R]))))", md, "<= 1", pf(md <= 1))
        row("AC13b", "$lab: same samples $nm: realized max |avg difference| and |Delta difference| in eps S (spec dry run: 1.0 and 1.0 regular, 0.13 chaotic)", "$(fmt(ma * 4)) / $(fmt(md * 8))", "reported", "INFO")
    end
    gapa = [maximum(abs.(Matrix(g[1].avg1) .- Matrix(c[1].avg1))) for (g, c) in zip(gpu.finals, cpu.finals)]
    gapd = [maximum(abs.(Matrix(g[1].delta) .- Matrix(c[1].delta))) for (g, c) in zip(gpu.finals, cpu.finals)]
    dx = [maximum(abs.(cpu.x[j] - gpu.x[j])) for j in 1:m]
    row("AC13c", "$lab: reported, not gated: orbit-level CPU-GPU gap over all orbits: |avg1_cpu - avg1_gpu| median / max; |Delta| median / max; final states apart (bitwise equal: $(count(==(0), dx)); median |dx| $(fmt(rmed(dx))), max $(fmt(maximum(dx))))",
        "$(fmt(rmed(gapa))) / $(fmt(maximum(gapa))); $(fmt(rmed(gapd))) / $(fmt(maximum(gapd)))", "reported", "INFO")
    cm = count(f -> f[1].complete && all(isfinite, f[1].avg1) && all(isfinite, f[1].avg2) && all(isfinite, f[1].delta), gpu.finals)
    row("AC13a", "$lab: orbits with complete and finite device finals", cm, "$m", pf(cm == m && gpu.provenance.backend === :gpu))
end

# Double64 on the device: the kernel compiles, or the recorded fallback reason is the pass.
function wbg_d13d(Kt, N, x0)
    sys = WrappedStandardMap(exact(Double64, Kt))
    xW = [torus_state(SVector{2,Double64}(Double64(z[1]), Double64(z[2]))) for z in x0[1:min(64, length(x0))]]
    reason = gpu_reason(Double64)
    res = try
        propagate(sys, xW, 2N; k = 1, tau = 1, seed = SEED, reducers = (wbg_red(N, wbg_trig4()), Benettin()), backend = :gpu)
    catch e
        return row("AC13d", "K=0.5 regular class, Double64, 64 orbits (or fewer): reducer chain on the device", sprint(showerror, e)[1:min(end, 1500)], "compiles or fallback", "FAIL")
    end
    if res.provenance.backend === :gpu
        row("AC13d", "K=0.5 regular class, Double64, 64 orbits (or fewer): the Double64 reducer kernel compiled and ran on the device (probe said: '$reason')", "compiled", "compiles or fallback", "PASS")
    else
        cpu = propagate(sys, xW, 2N; k = 1, tau = 1, seed = SEED, reducers = (wbg_red(N, wbg_trig4()), Benettin()), backend = :cpu)
        row("AC13d", "K=0.5 regular class, Double64, 64 orbits (or fewer): the result is the CPU result (finals bitwise equal to the :cpu run: $(isequal(res.finals, cpu.finals))), recorded fallback reason", res.provenance.backend_fallback_reason, "compiles or fallback", fb(res.provenance.backend_fallback_reason))
    end
end

# R5 on the device: ns per orbit-step of the reader for (K, H) in {(1, 1), (4, 1), (4, 11)} against the stream alone and the E2 reader, CPU (all threads) and device, K_std = 2 pi wrapped, 2000 steps.
function wbg_timing()
    sys0 = WrappedStandardMap(2pi)
    n = 2000
    norb = parse(Int, get(ENV, "GPU_WB_TIMING_M", "20000"))
    rng = StableRNG(SEED + 16)
    xw = [torus_state(SVector{2,Float64}(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:norb]
    wins(H) = H == 1 ? (1000,) : (16, 32, 64, 128, 256, 400, 500, 600, 700, 800, 1000)
    cases = Pair{String,Any}["stream alone (Benettin, k = 1)" => (Benettin(),), "E2 reader (window 1000) + Benettin" => (naff_red(1000), Benettin())]
    for (K, H) in ((1, 1), (4, 1), (4, 11))
        push!(cases, "WeightedBirkhoff (K, H) = ($K, $H) + Benettin" => (WeightedBirkhoff(; window = wins(H), observables = ntuple(i -> CosObservable(1, i), K)), Benettin()))
    end
    for (nm, rs) in cases
        run(b, c) = propagate(sys0, xw, n; k = 1, tau = 1, seed = SEED, reducers = rs, backend = b, chunk = c)
        cpu_t = gpu_t = Inf
        prov = nothing
        try
            run(:gpu, AC13_GPU_CHUNK)
            for _ in 1:3
                t = @elapsed res = run(:gpu, AC13_GPU_CHUNK)
                prov = res.provenance
                gpu_t = min(gpu_t, t)
            end
        catch e
            row("AC13", "R5 $nm: kernel on the device compile / run error", sprint(showerror, e)[1:min(end, 1500)], "compiles", "FAIL")
            continue
        end
        run(:cpu, nothing)
        for _ in 1:2
            cpu_t = min(cpu_t, @elapsed run(:cpu, nothing))
        end
        steps = norb * n
        row("AC13", "R5 $nm (provenance backend $(prov.backend)); $norb orbits x $n steps, K = 2 pi, wrapped, tau = 1: GPU $(fmt(gpu_t / steps * 1e9)) ns/orbit-step, CPU ($(Threads.nthreads()) threads) $(fmt(cpu_t / steps * 1e9)) ns/orbit-step, speedup $(fmt(cpu_t / gpu_t))",
            prov.backend === :gpu ? gpu_t / steps * 1e9 : NaN, "recorded", prov.backend === :gpu ? "INFO" : fb(prov.backend_fallback_reason))
    end
end

want("weighted-birkhoff") && let
    nreg = parse(Int, get(ENV, "GPU_WB_M", "1024"))
    N = parse(Int, get(ENV, "GPU_WB_N", "1000"))
    x0, ncand, nrep = regular_starts(nreg)
    sets = (("trig4", wbg_trig4(), 4), ("std5", wb_std5(), 5))
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC13a", "K=0.5 regular class, $T, N=$N: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        xT = [SVector{2,T}(T(z[1]), T(z[2])) for z in x0]
        for (nm, obs, ncol) in sets
            wbg_d13("K=0.5 regular class (S7: $(length(x0)) of $ncand candidates regular, $nrep replaced)", nm, T, 0.5, xT, N, obs, ncol)
        end
    end
    for T in (Float64, Double64)
        reason = gpu_reason(T)
        isempty(reason) || (row("AC13a", "K=10 chaotic class, $T, N=$N: each backend against its own BigFloat reference", "-", "<= 1", fb(reason)); continue)
        rng = StableRNG(SEED + 50000 + 1003)
        xT = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:nreg]
        for (nm, obs, ncol) in sets
            wbg_d13("K=10 chaotic class", nm, T, 10.0, xT, N, obs, ncol)
        end
    end
    wbg_d13d(0.5, N, x0)
    wbg_timing()
end

want("coupled-maps") && include(joinpath(ROOT, "test", "coupled-maps", "cm_helpers.jl"))   # cm_*: CMCase, n_H, tau, torus distance (definitions only)

# ---------------------------------------------------------------------------- AC18 GPU parts and D105d (coupled-maps)
# Y3 at the published parameters and at (2, 0.5, 0.1): 1e5 seeded orbits; Y4 (chain, N in {2, 5}, K = 3, eps = 0.2): 1e4 orbits through the tangent stream with
# fixed-size (SVector) states, n = n_H, k = 4, tau = 1, W0 = the first four columns of the identity; wrapped Y3 (1e4) and wrapped chain N = 2 (1e4) as extra rows.
# Bound tau(n) of the spec (Float64: eps 2^-52); the Double64 rows are the compile gate (run when the probe lets Double64 on the device, else the recorded reason).
MODTAG[] = "coupled-maps: "
cm_frame(::Type{T}, d) where {T} = SMatrix{d,4,T}(ntuple(i -> T(((i - 1) % d) == ((i - 1) ÷ d)), d * 4))
function cmg_row(ac, label, c, T, Z, epsT; wrapped = false)
    sys = wrapped ? cm_wsys(T, c) : cm_sys(T, c)
    n = cm_nH(c); d = cm_d(c)
    kw = (k = 4, tau = 1, seed = CM_SEED, W0 = fill(cm_frame(T, d), length(Z)))
    cpu = propagate(sys, Z, n; kw..., backend = :cpu)
    gpu, why = on_gpu(sys, Z, n; kw...)
    why === nothing || return row(ac, label, "-", "1", fb(why))
    lsys = cm_sys(T, c)
    ratio = 0.0; dmax = 0.0
    for j in eachindex(Z)
        if wrapped
            zl = cm_orbit(lsys, T.(lift(sys, Z[j])), n)
            diff = cm_absdiff(lift(sys, cpu.x[j]), lift(sys, gpu.x[j]))
            tau = cm_tau(epsT, c, zl, n)
        else
            diff = cm_tdist(cpu.x[j], gpu.x[j])
            tau = cm_tau(epsT, c, cm_orbit(sys, Z[j], n), n)
        end
        ratio = max(ratio, diff / tau); dmax = max(dmax, diff)
    end
    wd = maximum(maximum(abs, cpu.W[j] .- gpu.W[j]) / max(1.0, maximum(abs, cpu.W[j])) for j in 1:min(2000, length(Z)))
    row(ac, label * " (max $(wrapped ? "lifted" : "torus") distance over the ensemble $(fmt(dmax)); max relative frame difference over the first 2000 orbits $(fmt(wd)))", ratio, "1", pf(ratio <= 1))
end
function cmg_try(ac, label, f)
    try
        f()
    catch e
        row(ac, label * ": kernel on the device compile / run error", sprint(showerror, e)[1:min(end, 1500)], "compiles", "FAIL")
    end
end
want("coupled-maps") && let
    rng = StableRNG(CM_SEED + 41181)
    lim = 4 * Float64(pi)
    for p in ((0.5, 0.1, 1.0e-3), (2.0, 0.5, 0.1))
        c = cm_y3(p...); n = cm_nH(c)
        Z = [SVector{4,Float64}(ntuple(_ -> (2 * rand(rng) - 1) * lim, 4)) for _ in 1:100_000]
        cmg_try("AC18", "Y3 $p", () -> cmg_row("AC18", "Y3 $p GPU vs CPU, n=n_H=$n, 100000 orbits through the stream (k = 4): max over orbits of torus distance / tau(n)", c, Float64, Z, cm_eps(Float64)))
        Zw = [SVector{8,Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:10_000]
        cmg_try("AC18", "wrapped Y3 $p", () -> cmg_row("AC18", "wrapped Y3 $p GPU vs CPU, n=n_H=$n, 10000 orbits: max over orbits of lifted difference (256 bits) / tau(n)", c, Float64, Zw, cm_eps(Float64); wrapped = true))
    end
    for N in (2, 5)
        c = cm_y4(N, 3.0, 0.2); n = cm_nH(c); d = 2N
        Z = [SVector{d,Float64}(cm_draw(rng, d, lim)) for _ in 1:10_000]
        cmg_try("D105d", "Y4 N=$N", () -> cmg_row("D105d", "Y4 N=$N (K=3, eps=0.2) GPU vs CPU, n=n_H=$n, 10000 orbits through the stream with fixed-size states (k = 4): max over orbits of torus distance / tau(n)", c, Float64, Z, cm_eps(Float64)))
    end
    c = cm_y4(2, 3.0, 0.2); n = cm_nH(c)
    Zw = [SVector{8,Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:10_000]
    cmg_try("D105d", "wrapped Y4 N=2", () -> cmg_row("D105d", "wrapped Y4 N=2 GPU vs CPU, n=n_H=$n, 10000 orbits (state length 4N = 8): max over orbits of lifted difference (256 bits) / tau(n)", c, Float64, Zw, cm_eps(Float64); wrapped = true))
    # Double64 compile gate
    reason = gpu_reason(Double64)
    for (nm, c, d) in (("Y3 (0.5, 0.1, 1e-3)", cm_y3(0.5, 0.1, 1.0e-3), 4), ("Y4 N=2", cm_y4(2, 3.0, 0.2), 4), ("Y4 N=5", cm_y4(5, 3.0, 0.2), 10))
        if isempty(reason)
            ZD = [SVector{d,Double64}(Double64.(cm_draw(rng, d, lim))) for _ in 1:5_000]
            cmg_try("AC18", "Double64 $nm", () -> cmg_row("AC18", "Double64 kernel gate, $nm: GPU vs CPU, 5000 orbits: max torus distance / tau(n_H) (eps 2^-100)", c, Double64, ZD, cm_eps(Double64)))
        else
            sys = cm_sys(Double64, c)
            probe_run = on_gpu(sys, [SVector{d,Double64}(ntuple(i -> Double64(i) / 7, d))], 2; k = 4, seed = CM_SEED)[2]
            row("AC18", "Double64 kernel gate, $nm: compile or reason (explicit :gpu run ends on the CPU with the reason recorded in provenance)",
                probe_run === nothing ? "ran on GPU" : probe_run, "compile or reason", fb(reason))
        end
    end
end

counts = Dict{String,Int}()
for (_, v) in ROWS
    k = startswith(v, "FALLBACK") ? "FALLBACK" : v
    counts[k] = get(counts, k, 0) + 1
end
println("# SUMMARY ", join(("$k=$v" for (k, v) in sort!(collect(counts))), ", "), " (", length(ROWS), " lines)")
