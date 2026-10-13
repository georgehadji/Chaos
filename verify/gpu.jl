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
            ", capability ", CUDA.capability(CUDA.device()), ", free memory ", round(CUDA.free_memory() / 2^20), " MiB")
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
let
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
for K in (2pi, 10.0)
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
let
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
let
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
let
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
let
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
if haskey(gpu_ns, Double64) && haskey(gpu_ns, Float64)
    row("AC11", "GPU Double64/Float64 time ratio, same ensemble as AC13 (per orbit-step)", gpu_ns[Double64] / gpu_ns[Float64], "reported, no threshold", "PASS")
else
    row("AC11", "GPU Double64/Float64 time ratio, same ensemble as AC13: Float64 GPU $(fmt(get(gpu_ns, Float64, NaN))) ns/orbit-step; CPU ratio $(fmt(cpu_ns[Double64] / cpu_ns[Float64])) (4 threads, this run)",
        "-", "reported, no threshold", fb(gpu_reason(Double64)))
end
println("# after the record: resolve_backend(:auto, Float64) = ", F.resolve_backend(:auto, Float64), "; (:auto, Double64) = ", F.resolve_backend(:auto, Double64))

# ------------------------------------------------------------------------------------------------ AC12a (tangent-core), last
# D33: symmetric form. Each backend is compared with a Double64 CPU reference on all 1024 seeded orbits; the numbers of orbits above
# 1e-6 may differ by at most 3, and on the orbits where both backends agree with the reference the CPU-GPU gap is within 2e-6.
# The spec's n = 1e5 is one launch of about 2 s, over the 1 s launch budget of D33, so propagate runs it on the CPU with the
# reason "launch exceeds the display watchdog" (FALLBACK line); the same criterion is then evaluated at the largest n whose
# launch fits the budget with one wave of 1024 orbits (n = 4e4). It runs last, after the record is written.
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
    row("AC12a", "K=0.5, $no orbits, n=$(fmt(N)), tau=$tau, one launch ($(fmt(tg)) s): orbits above 1e-6 against the Double64 reference, CPU / GPU; max CPU-GPU gap over the $(count(both)) orbits where both agree with the reference",
        "$ncpu / $ngpu (diff $(abs(ncpu - ngpu))); $(fmt(gap))", "count diff <= 3; gap <= 2e-6", pf(ok))
    row("AC12a", "same, maxima over all $no orbits (chaotic layer included): |S_cpu - S_ref| $(fmt(maximum(dcr))), |S_gpu - S_ref| $(fmt(maximum(dgr))), |S_cpu - S_gpu| $(fmt(maximum(dcg)))",
        maximum(dcg), "-", "INFO")
end
ac12a(10^5)
ac12a(4 * 10^4)

counts = Dict{String,Int}()
for (_, v) in ROWS
    k = startswith(v, "FALLBACK") ? "FALLBACK" : v
    counts[k] = get(counts, k, 0) + 1
end
println("# SUMMARY ", join(("$k=$v" for (k, v) in sort!(collect(counts))), ", "), " (", length(ROWS), " lines)")
