# Stage 04 device check of the recording of the chaos indicators (Amendment D86, criteria D86p (G), the GPU part of D86y, oracle O21), local GPU only. Run from the repo root:
#   julia --project=verify -t 23 verify/experiments-indicators-device.jl
# Output (raw): workspaces/module-build/stages/04-verify/output/gpu-rtx5060-d86-verify-log.txt; figure: experiments-d86-device.png.
# The record file data/gpu-benchmark.toml is only read: the per-spec records of D54 are taken here into a file of the scratch folder (the GPU record is not retaken for D86).
# Sections (EXD_ONLY=a,b runs some of them): TABLE, RECORD, PARITY, O21, DOUBLE, CAMPAIGN, PLOT.
#   TABLE    backend table with the record in force; RECORD  per-spec records (ns per orbit-step device and CPU threads) into the scratch file, then the table following them;
#   PARITY   device recording against the CPU recording of the same items, every kind (orbit-level gap reported, not gated);
#   O21      the algebra gates: window kinds on the samples of the device orbit (trace kernel), tangent and pair chains bitwise against the verified device path (propagate with the
#            same frames), SALI against the BigFloat(256) reference at the device orbit (sali criterion 18 form);
#   DOUBLE   Double64 on the device records the fallback;  CAMPAIGN  the reduced campaign of D86y on the device (auto with the per-spec records, and a gpu request).
include(joinpath(@__DIR__, "experiments-indicators-prep.jl"))
using CUDA, HDF5
import CairoMakie as CM
const F = Chaoslab.Foundation
const ONLY = filter(!isempty, split(get(ENV, "EXD_ONLY", ""), ","))
mkpath(OUT)
const LOG = open(get(ENV, "EXD_LOG", joinpath(OUT, "gpu-rtx5060-d86-verify-log.txt")), isempty(ONLY) ? "w" : "a")
fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const COUNTS = Dict("PASS" => 0, "FAIL" => 0, "INFO" => 0, "FALLBACK" => 0)
function say(args...)
    s = string(args...)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(id, quantity, value, bound, verdict)
    v = startswith(verdict, "FAIL") ? "FAIL" : startswith(verdict, "PASS") ? "PASS" : startswith(verdict, "FALLBACK") ? "FALLBACK" : "INFO"
    COUNTS[v] += 1
    say(string(id, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
info(id, msg) = (COUNTS["INFO"] += 1; say("$id | INFO | $msg"))
function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 2000)], "-", "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end
bestof(f, reps) = (best = Inf; out = nothing; for _ in 1:reps; t = @elapsed out = f(); best = min(best, t); end; (best, out))

say("# D86 device check, julia ", VERSION, ", threads ", NT, ", cpu ", Sys.CPU_NAME)
say("# device ", typeof(F.DEVICE[]), ", ", CUDA.name(CUDA.device()), ", CUDA.jl ", pkgversion(CUDA), ", runtime ", CUDA.runtime_version(), ", free ", round(CUDA.free_memory() / 2^20), " MiB, stack limit ", CUDA.limit(CUDA.LIMIT_STACK_SIZE), " B")
F.DEVICE[] isa F.NoDevice && (say("# no functional device"); exit(1))
const DEVNAME = F.probe(F.DEVICE[], Float64).device
say("# probe(Float64): ", F.probe(F.DEVICE[], Float64))
const GIT_HEAD = strip(read(`git -C $ROOT rev-parse HEAD`, String))
say("# git head $GIT_HEAD; workaround F1 active (EXI_PATCH): ", get(ENV, "EXI_PATCH", "1") == "1")
const REAL_RECORD = joinpath(ROOT, "data", "gpu-benchmark.toml")
const REAL_RECORD_SHA = bytes2hex(open(sha256, REAL_RECORD))
say("# data/gpu-benchmark.toml sha256 before: $REAL_RECORD_SHA (read only)")
const PP = prepare()
const PLAN = EX.plan_run(PP.prereg, PP.panel)
const SPECS = Dict(s["id"] => s for s in PLAN.block.raw["indicator"])
const K19 = PLAN.block.raw["k19"]["threads"]
const ITEM_OF = Dict(it.hash => it for it in PLAN.items)
const SCRD = joinpath(SCR, "device"); mkpath(SCRD)
const REC_SCRATCH = joinpath(SCRD, "gpu-benchmark-d86.toml")
const SEED = 20261012
const FALLBACKS = F.FALLBACK_REASONS
say("# plan: ", PLAN.counts, "; FALLBACK_REASONS (", length(FALLBACKS), "): ", join(("$i: $r" for (i, r) in enumerate(FALLBACKS)), " | "))

# ---------------------------------------------------------------------------------------------------
# helpers shared with the CPU script (the same stand-alone constructions)
function with_prec(f, label::AbstractString)
    T, nb = Chaoslab.parse_precision(label)
    return nb === nothing ? f(T) : setprecision(() -> f(T), BigFloat, nb)
end
function sys_x(item, T)
    sys = if item.system == "linear"
        Chaoslab.y5_control(parse(Float64, item.param)).sys
    elseif item.coordinate_mode == "wrapped"
        Chaoslab.WrappedStandardMap(Chaoslab.exact(T, item.param))
    else
        Chaoslab.StandardMap(Chaoslab.exact(T, item.param))
    end
    x = SVector(Chaoslab.exact(T, item.x0), Chaoslab.exact(T, item.p0))
    return sys, (item.coordinate_mode == "wrapped" && item.system != "linear" ? Chaoslab.torus_state(x) : x)
end
frame(item, T) = SMatrix{2, 2, T}(EX.deviation_frame(T, item.direction, 2))
colvec(t, name) = t.series[:, findfirst(==(name), t.columns)]
isbitseq(a, b) = length(a) == length(b) && all(isequal.(a, b))
function member_cols(m, fin)
    m == "SALI" && return Float64[fin.sali, fin.L]
    m == "FLI" && return Float64[fin.fli, fin.fli1]
    m == "MEGNO" && return Float64[v for c in 1:2 for v in (fin.Y[c], fin.Ybar[c], fin.Ybar20[c])]
    m == "RLI" && return Float64[fin.DL, fin.RLI]
    error(m)
end
const MEMBER_NAMES = Dict("SALI" => ["SALI.sali", "SALI.L"], "FLI" => ["FLI.fli", "FLI.fli1"],
                          "MEGNO" => ["MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"], "RLI" => ["RLI.DL", "RLI.RLI"])
function ensemble(id, n; prec = "Float64", mode = "wrapped", standard = true)
    base = [it for it in PLAN.items if it.method == id && it.precision == prec && it.coordinate_mode == mode && (!standard || it.system == "standard")]
    return NamedTuple[base[mod1(i, length(base))] for i in 1:n]
end
runsteps(id) = (spec = SPECS[id]; spec["kind"] in ("frequency", "birkhoff") ? 2 * maximum(spec["windows"]) : spec["kind"] == "phase_advance" ? 2 * maximum(spec["windows"]) + 1 : spec["kind"] == "pair" ? 300 : 4000)

# device orbit trace (the kernel of verify/gpu.jl, section rli): proj of the state before step i, i = 1 .. steps; the final state
function dev_trace_kernel!(D, X, sys, steps, proj)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if j <= length(X)
        x = X[j]
        for i in 1:steps
            D[j, i] = proj(x)
            x = Chaoslab.forward(sys, x)
        end
        X[j] = x
    end
    return nothing
end
function dev_trace(sys, x0::Vector{S}, N, proj; chunk = 5000) where {S}
    R = typeof(proj(x0[1]))
    m = length(x0)
    X = CuArray(x0)
    P = Matrix{R}(undef, m, N)
    for s0 in 0:chunk:(N - 1)
        s1 = min(s0 + chunk, N)
        D = CuArray{R}(undef, m, s1 - s0)
        @cuda threads = 256 blocks = cld(m, 256) dev_trace_kernel!(D, X, sys, s1 - s0, proj)
        P[:, (s0 + 1):s1] = Array(D)
        CUDA.unsafe_free!(D)
    end
    return P, Array(X)
end
idproj(z) = z

# ---------------------------------------------------------------------------------------------------
section("TABLE") do
    rec = EX.gpu_record()
    info("D86y", "record in force (data/gpu-benchmark.toml, Float64 table of $DEVNAME): " * (rec === nothing ? "none" : "threads $(rec.threads), gpu_faster $(rec.gpu_faster), per-method records for " * join(sort(collect(keys(rec.recorder))), ", ")))
    tbl = EX.backend_table(PLAN; threads = NT, record = rec)
    for ((p, m), r) in sort!(collect(tbl); by = first)
        haskey(SPECS, m) || continue
        say("TABLE | ", p, " | ", m, " (", SPECS[m]["kind"], ") | ", r.request, " | ", r.backend, " | ", r.device, " | ", r.reason, " | ", r.note)
    end
    rows = [(k, r) for (k, r) in tbl if haskey(SPECS, k[2]) && k[1] == "Float64"]
    ok = all(x -> x[2].backend === :cpu, rows)
    emit("D86o/D86y", "backend table with the record in force at $NT threads (no per-spec record exists yet): every Float64 spec row resolves to the CPU; CPU-only kinds with the sixth reason, the others with the seventh and no new reason",
         join(unique(x[2].reason for x in rows), " | "), "reasons among the nine of FALLBACK_REASONS", pf(ok && all(x -> x[2].reason in FALLBACKS, rows)))
    for id in ("st", "fq")
        rr = [r for ((p, m), r) in tbl if m == id && p == "Float64"]
        emit("D86o/R16", "spec $id ($(SPECS[id]["kind"])) is CPU only: Float64 auto resolves to the CPU with the sixth reason", join(unique(r.reason for r in rr), ","), FALLBACKS[6], pf(all(r -> r.backend === :cpu && r.reason == FALLBACKS[6], rr)))
    end
end

# ---------------------------------------------------------------------------------------------------
const ENS_ORB = parse(Int, get(ENV, "EXD_ENS_ORB", "1024"))
const COST = Dict{String, NamedTuple}()
section("RECORD") do
    real = TOML.parsefile(REAL_RECORD)
    haskey(real, DEVNAME) || error("no table for $DEVNAME in the record file")
    new = Dict{String, Any}(DEVNAME => deepcopy(real[DEVNAME]))
    rt = new[DEVNAME]["Float64"]
    recs = get!(rt, "recorder", Dict{String, Any}())
    for id in ["tg", "pr", "st", "fq", "pa", "bk"]
        items = ensemble(id, ENS_ORB)
        steps = ENS_ORB * runsteps(id)
        ctx = EX.make_ctx(PLAN, items)
        EX.device_records(ctx, items; hook = false)
        tg, g = bestof(() -> EX.device_records(ctx, items), 2)
        EX.compute_records(ctx, items, NT)
        tc, _ = bestof(() -> EX.compute_records(ctx, items, NT), 2)
        if g.recs === nothing
            say("RECORD | $id | $(SPECS[id]["kind"]) | FALLBACK(", g.prov.reason, ") | cpu $(fmt(tc / steps * 1e9)) ns per orbit-step on $NT threads ($(fmt(tc / steps * 1e9 * NT)) thread-ns)")
            COST[id] = (kind = SPECS[id]["kind"], cpu = tc / steps * 1e9, gpu = NaN, fallback = g.prov.reason)
            continue
        end
        cpu, gpu = tc / steps * 1e9, tg / steps * 1e9
        COST[id] = (kind = SPECS[id]["kind"], cpu = cpu, gpu = gpu, fallback = "")
        recs[id] = Dict{String, Any}("threads" => NT, "cpu_ns_per_step" => cpu, "gpu_ns_per_step" => gpu, "gpu_faster" => gpu <= cpu, "date" => Libc.strftime("%Y-%m-%d", time()),
                                     "launch_steps" => g.prov.launch_steps, "orbits" => ENS_ORB, "steps" => runsteps(id), "spec_sha256" => EX.spec_hash(SPECS[id]))
        say("RECORD | $id | $(SPECS[id]["kind"]) | $ENS_ORB orbits x $(runsteps(id)) steps | device $(fmt(gpu)) ns per orbit-step | cpu $NT threads $(fmt(cpu)) ns per orbit-step ($(fmt(cpu * NT)) thread-ns) | device/cpu $(fmt(gpu / cpu)) | gpu_faster $(gpu <= cpu) | launches $(g.prov.launches) of $(g.prov.launch_steps) steps, batches $(g.prov.gpu_chunks) of $(g.prov.chunk_size)")
    end
    open(io -> TOML.print(io, new), REC_SCRATCH, "w")
    rec = EX.read_gpu_record(REC_SCRATCH, DEVNAME)
    info("D86y", "per-spec records written to $REC_SCRATCH (scratch folder; data/gpu-benchmark.toml untouched)")
    tbl = EX.backend_table(PLAN; threads = NT, record = rec)
    for ((p, m), r) in sort!(collect(tbl); by = first)
        haskey(SPECS, m) || continue
        say("TABLE-AFTER | ", p, " | ", m, " (", SPECS[m]["kind"], ") | ", r.request, " | ", r.backend, " | ", r.device, " | ", r.reason, " | ", r.note)
    end
    ok = true
    for id in ("tg", "pr", "pa", "bk")
        r = tbl[("Float64", id)]
        want = get(recs, id, nothing)
        faster = want !== nothing && want["gpu_faster"]
        ok &= faster ? (r.backend === :gpu && r.reason == "") : (r.backend === :cpu && r.reason in FALLBACKS)
    end
    for id in ("st", "fq")
        r = tbl[("Float64", id)]
        ok &= r.backend === :cpu && r.reason == FALLBACKS[6]
    end
    emit("D86o/D86y", "backend table at $NT threads following the per-spec records (gpu where the record says gpu_faster, else the cpu with the seventh reason); st and fq with the sixth reason", join(("$id:$(tbl[("Float64", id)].backend)" for id in sort(collect(keys(SPECS)))), " "), "follows the records", pf(ok))
    wrong = EX.read_gpu_record(REC_SCRATCH, DEVNAME)
    for id in ("tg",)
        haskey(wrong.recorder, id) || continue
        r2 = merge(wrong.recorder[id], (spec_sha256 = "0"^64,))
        tbl2 = EX.backend_table(PLAN; threads = NT, record = merge(wrong, (recorder = merge(wrong.recorder, Dict(id => r2)),)))
        emit("D86o/R16", "a record whose spec_sha256 differs from the spec hash resolves tg to the CPU with the note `record taken on another chain` and the seventh reason", "$(tbl2[("Float64", id)].backend), $(tbl2[("Float64", id)].reason), note: $(tbl2[("Float64", id)].note)", "cpu, seventh reason, note",
             pf(tbl2[("Float64", id)].backend === :cpu && occursin("another chain", tbl2[("Float64", id)].note) && tbl2[("Float64", id)].reason == FALLBACKS[7]))
    end
end

# ---------------------------------------------------------------------------------------------------
# PARITY: the device recording against the CPU recording of the same items, every kind. Reported, not gated (no cross-backend tolerance on chaotic orbits, R20 item 2).
const PAR = Dict{String, NamedTuple}()
function gap_stats(id, items, g, cpu)
    t = [(a.indicators[id], b.indicators[id]) for (a, b) in zip(g.recs, cpu)]
    sera = [x for (x, y) in t if x.series !== nothing]
    maxd = 0.0; neq = 0
    gaps = Float64[]
    for ((x, y), (a, b)) in zip(t, zip(g.recs, cpu))
        tabs = [(x.series, y.series), (x.final, y.final), (x.scalars, y.scalars)]
        d = 0.0
        for (u, v) in tabs
            u === nothing && continue
            d = max(d, maximum((abs(p - q) for (p, q) in zip(u, v) if isfinite(p) && isfinite(q)); init = 0.0))
        end
        same = all(u === nothing ? true : isbitseq(u, v) for (u, v) in tabs) && isbitseq(a.state, b.state) && isbitseq(a.t, b.t)
        neq += same
        push!(gaps, d)
        maxd = max(maxd, d)
    end
    dx = [maximum(abs.(a.state[end, :] .- b.state[end, :])) for (a, b) in zip(g.recs, cpu)]
    return (neq = neq, maxd = maxd, gaps = gaps, dx = dx, flags = count(i -> t[i][1].flags !== nothing && !isbitseq(t[i][1].flags, t[i][2].flags), eachindex(t)))
end
section("PARITY") do
    for id in ["tg", "pr", "pa", "bk", "fq", "st"]
        items = NamedTuple[it for it in PLAN.items if it.method == id && it.precision == "Float64" && it.coordinate_mode == "wrapped"]
        ctx = EX.make_ctx(PLAN, items)
        g = EX.device_records(ctx, items)
        if g.recs === nothing
            emit("D86p", "$id ($(SPECS[id]["kind"])): device recording of $(length(items)) items", "FALLBACK($(g.prov.reason))", id in ("fq", "st") ? "the sixth reason (CPU only)" : "device", id in ("fq", "st") ? (g.prov.reason == FALLBACKS[6] ? "PASS" : "FAIL") : "FAIL")
            continue
        end
        cpu = EX.compute_records(ctx, items, NT)
        s = gap_stats(id, items, g, cpu)
        PAR[id] = s
        nok = count(r -> all(isfinite, r.state), g.recs)
        emit("D86p", "$id ($(SPECS[id]["kind"])): $(length(items)) items on the device: provenance backend / launches / batches", "$(g.prov.backend) / $(g.prov.launches) / $(g.prov.gpu_chunks) of $(g.prov.chunk_size)", "backend gpu", pf(g.prov.backend === :gpu && g.prov.reason == ""))
        info("D86p", "$id: device records bitwise equal to the CPU records (series, tables, scalars, state, t): $(s.neq) of $(length(items)); flag words differing: $(s.flags); largest |difference| of a stored table value $(fmt(s.maxd)); final states apart by more than 1e-9: $(count(>(1e-9), s.dx)); median / max final-state gap $(fmt(sort(s.dx)[div(end, 2)])) / $(fmt(maximum(s.dx))) (orbit-level gap: reported, not gated)")
        fin = all(o -> all(isfinite, o.state), g.recs)
        emit("D86p/R18", "$id: every device record has a finite orbit point and passes the health rule of its kind (no non-finite value outside its flag)", fin && all(i -> isempty(EX.ind_nonfinite(id, g.recs[i].indicators[id])), eachindex(g.recs)), "true", pf(fin && all(i -> isempty(EX.ind_nonfinite(id, g.recs[i].indicators[id])), eachindex(g.recs))))
    end
end

# ---------------------------------------------------------------------------------------------------
# O21. The window kinds: the device reducer against the CPU pure function on the samples of the device orbit (the D82 form: algebra only, never trajectories).
include(joinpath(ROOT, "test", "sali", "sl_helpers.jl"))    # sl_*: BigFloat reference of SALI at given orbit points, Psi_n, B_SALI (as verify/gpu.jl, section sali)
epsT = eps(Float64)
section("O21") do
    # ---- pa (E2): same samples, algebra only: |nu_device reducer - nu CPU pure function on the device samples| <= 2 B_E2 (naff-fma criterion 15(b) form); B_E2 = 8 eps (1 + max|Delta_j| / (2 pi))
    PA = SPECS["pa"]
    items = NamedTuple[it for it in PLAN.items if it.method == "pa" && it.precision == "Float64" && it.coordinate_mode == "wrapped" && it.system == "standard"]
    ctx = EX.make_ctx(PLAN, items)
    g = EX.device_records(ctx, items)
    if g.recs === nothing
        emit("O21", "pa on the device", "FALLBACK($(g.prov.reason))", "device", "FAIL")
    else
        byK = Dict{String, Vector{Int}}()
        for (i, it) in enumerate(items); push!(get!(byK, it.param, Int[]), i); end
        Wm = maximum(PA["windows"]); nsamp = 2Wm + 1
        worst = 0.0; nrow = 0; nstate = 0; nbad = 0; nnan = 0
        for (param, idx) in byK
            sys = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
            x0 = [last(sys_x(items[i], Float64)) for i in idx]
            P, XN = dev_trace(sys, x0, nsamp, idproj)
            for (k, i) in enumerate(idx)
                # the samples: states after steps 1 .. nsamp = column 2 .. nsamp of the trace and the final state
                S = vcat([P[k, j] for j in 2:nsamp], [XN[k]])
                z = ComplexF64[Complex(s[1], s[2]) for s in S]
                inc = Chaoslab.phase_increments(z)
                t = g.recs[i].indicators["pa"]
                nstate += (g.recs[i].state[end, :] == XN[k][1:2] || true) ? 0 : 1
                BE2 = 8 * epsT * (1 + maximum(abs, inc) / (2pi))
                for (h, W) in enumerate(PA["windows"])
                    nu = Chaoslab.weighted_phase_frequency(inc, W; p = PA["birkhoff_order"])
                    for (c, col) in enumerate((1, 2))
                        a, b = t.final[h, col], nu[col]
                        if isnan(a) && isnan(b)         # a planar sample at the origin (the `undefined window` flag): both undefined
                            nnan += 1
                            continue
                        end
                        r = abs(a - b) / (2 * BE2)
                        worst = max(worst, r); nrow += 1; nbad += r <= 1 ? 0 : 1
                    end
                end
            end
        end
        emit("O21/D86p", "pa, Float64, $(length(items)) orbits x 2 windows x 2 frequencies ($nrow compared, $nnan undefined on both sides), same samples (trace of the device orbit): max |nu_device reducer - nu CPU pure function| / (2 B_E2), B_E2 = 8 eps (1 + max|Delta|/2 pi)", fmt(worst), "<= 1", pf(nbad == 0 && nrow > 100))
        # the trace's final state is the recorder's own orbit: recorded state at the last grid step equals the trace state at the last handoff
        eqs = 0; tot = 0
        for (param, idx) in byK
            sys = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
            x0 = [last(sys_x(items[i], Float64)) for i in idx]
            P, XN = dev_trace(sys, x0, nsamp, idproj)
            for (k, i) in enumerate(idx)
                tot += 1
                eqs += (g.recs[i].state[end, 1] == XN[k][1] && g.recs[i].state[end, 2] == XN[k][2]) ? 1 : 0
            end
        end
        emit("O21", "pa: the trace kernel's final state equals the recorder's stored orbit point at the last grid step (the samples are the recorder's own device orbit), bitwise", "$eqs of $tot", "all", pf(eqs == tot))
    end
    # ---- bk: (b1) the observable values computed by the device at each handoff, read back; (b2) the CPU evaluation of the observables at the device states. Gates 4 eps S and 8 eps S (weighted-birkhoff criterion 13(b))
    BK = SPECS["bk"]
    items = NamedTuple[it for it in PLAN.items if it.method == "bk" && it.precision == "Float64" && it.coordinate_mode == "wrapped" && it.system == "standard"]
    ctx = EX.make_ctx(PLAN, items)
    g = EX.device_records(ctx, items)
    if g.recs === nothing
        emit("O21", "bk on the device", "FALLBACK($(g.prov.reason))", "device", "FAIL")
    else
        obs = Tuple(map(BK["observables"]) do ob
            c = ob["coordinate"]
            ob["kind"] == "cos" ? Chaoslab.CosObservable(c, ob["frequency"]) : ob["kind"] == "sin" ? Chaoslab.SinObservable(c, ob["frequency"]) : Chaoslab.LiftedObservable(c, c + 2)
        end)
        K = length(obs); Wm = maximum(BK["windows"]); nsamp = 2Wm
        proj = z -> map(o -> Float64(o(z)), obs)
        byK = Dict{String, Vector{Int}}()
        for (i, it) in enumerate(items); push!(get!(byK, it.param, Int[]), i); end
        wa = Dict(:b1 => 0.0, :b2 => 0.0); wd = Dict(:b1 => 0.0, :b2 => 0.0); nbad = 0; nrow = 0
        for (param, idx) in byK
            sys = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
            x0 = [last(sys_x(items[i], Float64)) for i in idx]
            Pobs, _ = dev_trace(sys, x0, nsamp + 1, proj)      # observable values at the states before steps 1 .. nsamp + 1: samples 0 .. nsamp are columns 2 .. nsamp + 1
            Pst, XN = dev_trace(sys, x0, nsamp, idproj)
            for (k, i) in enumerate(idx)
                t = g.recs[i].indicators["bk"]
                S1 = Matrix{Float64}(undef, nsamp, K)
                for j in 1:nsamp, c in 1:K
                    S1[j, c] = Pobs[k, j + 1][c]
                end
                st = vcat([Pst[k, j] for j in 2:nsamp], [XN[k]])
                S2 = Matrix{Float64}(undef, nsamp, K)
                for j in 1:nsamp, c in 1:K
                    S2[j, c] = Float64(obs[c](st[j]))
                end
                for (nm, S) in ((:b1, S1), (:b2, S2))
                    r = Chaoslab.weighted_birkhoff(S, Tuple(BK["windows"]); order = BK["birkhoff_order"])
                    for (h, W) in enumerate(BK["windows"]), c in 1:K
                        row = (h - 1) * K + c
                        Sv = Float64(r.S[c, h])
                        ra = max(abs(t.final[row, 1] - r.avg1[c, h]), abs(t.final[row, 2] - r.avg2[c, h])) / (4 * epsT * Sv)
                        rd = abs(t.final[row, 5] - (r.avg1[c, h] - r.avg2[c, h])) / (8 * epsT * Sv)
                        isfinite(ra) && (wa[nm] = max(wa[nm], ra)); isfinite(rd) && (wd[nm] = max(wd[nm], rd))
                        nrow += 1
                        nbad += (ra <= 1 || !isfinite(ra)) && (rd <= 1 || !isfinite(rd)) ? 0 : 1
                    end
                end
            end
        end
        emit("O21/D86p", "bk, Float64, $(length(items)) orbits, 4 observables x 2 windows: (b1) the samples the device computed (trace of the observable values), device reducer against weighted_birkhoff: max |avg difference| / (4 eps S), max |Delta difference| / (8 eps S)", "$(fmt(wa[:b1])), $(fmt(wd[:b1]))", "<= 1, <= 1", pf(wa[:b1] <= 1 && wd[:b1] <= 1))
        emit("O21/D86p", "bk: (b2) the CPU evaluation of the observables at the device states", "$(fmt(wa[:b2])), $(fmt(wd[:b2]))", "<= 1, <= 1", pf(wa[:b2] <= 1 && wd[:b2] <= 1))
    end
    # ---- tg and pr: the recorder's device finals against the verified device path (propagate with the same frames and chain), bitwise
    items = NamedTuple[it for it in PLAN.items if it.method == "tg" && it.precision == "Float64" && it.coordinate_mode == "wrapped" && it.system == "standard"]
    ctx = EX.make_ctx(PLAN, items)
    g = EX.device_records(ctx, items)
    if g.recs === nothing
        emit("O21", "tg on the device", "FALLBACK($(g.prov.reason))", "device", "FAIL")
    else
        members = SPECS["tg"]["members"]
        red = (Chaoslab.FLI(replace = false), Chaoslab.MEGNO(replace = false), Chaoslab.SALI())
        byK = Dict{String, Vector{Int}}()
        for (i, it) in enumerate(items); push!(get!(byK, it.param, Int[]), i); end
        neq = ntot = 0; maxd = 0.0; ng = 0
        for (param, idx) in byK
            sys = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
            x0 = [last(sys_x(items[i], Float64)) for i in idx]
            W0 = [frame(items[i], Float64) for i in idx]
            res = Chaoslab.Tangent.propagate(sys, x0, 4000; k = 2, tau = 1, seed = SEED, reducers = red, backend = :gpu, W0 = W0)
            res.provenance.backend === :gpu && (ng += length(idx))
            for (k, i) in enumerate(idx)
                got = Float64[colvec(g.recs[i].indicators["tg"], nm)[end] for m in members for nm in MEMBER_NAMES[m]]
                ref = Float64[v for (m, f) in zip(members, res.finals[k]) for v in member_cols(m, f)]
                ntot += 1; neq += isbitseq(got, ref) ? 1 : 0
                maxd = max(maxd, maximum(abs.(got .- ref)))
            end
        end
        emit("O21/D86p", "tg (FLI, MEGNO, SALI) Float64: recorder device columns at N = 4000 against the verified device path (propagate, same chain, same frames, backend gpu; provenance gpu for $ng of $(length(items))): orbits bitwise equal / largest difference", "$neq of $ntot / $(fmt(maxd))", "all bitwise equal", pf(neq == ntot && ng == length(items)))
    end
    items = NamedTuple[it for it in PLAN.items if it.method == "pr" && it.precision == "Float64" && it.coordinate_mode == "wrapped" && it.system == "standard"]
    ctx = EX.make_ctx(PLAN, items)
    g = EX.device_records(ctx, items)
    if g.recs === nothing
        emit("O21", "pr on the device", "FALLBACK($(g.prov.reason))", "device", "FAIL")
    else
        sp = SPECS["pr"]
        red = (Chaoslab.RLI(replace = false), Chaoslab.MEGNO(replace = true))
        byK = Dict{String, Vector{Int}}()
        for (i, it) in enumerate(items); push!(get!(byK, it.param, Int[]), i); end
        neq = ntot = 0; maxd = 0.0; ng = 0
        for (param, idx) in byK
            sys = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
            starts = map(idx) do i
                it = items[i]
                x = SVector(Chaoslab.exact(Float64, it.x0), Chaoslab.exact(Float64, it.p0))
                u = SVector{2, Float64}(Float64.(g.recs[i].attrs["xi0"])...)
                ps = Chaoslab.pair_start(sys, x, u, sp["separation"]; orbit = it.orbit)
                (ps.z, SMatrix{4, 2, Float64}(u[1], u[2], 0.0, 0.0, 0.0, 0.0, u[1], u[2]))
            end
            res = Chaoslab.Tangent.propagate(Chaoslab.PairMap(sys), first.(starts), 300; k = 2, tau = 1, seed = SEED, reducers = red, backend = :gpu, W0 = last.(starts))
            res.provenance.backend === :gpu && (ng += length(idx))
            for (k, i) in enumerate(idx)
                got = Float64[colvec(g.recs[i].indicators["pr"], nm)[end] for m in sp["members"] for nm in MEMBER_NAMES[m]]
                ref = Float64[v for (m, f) in zip(sp["members"], res.finals[k]) for v in member_cols(m, f)]
                ntot += 1; neq += isbitseq(got, ref) ? 1 : 0
                maxd = max(maxd, maximum(abs.(got .- ref)))
            end
        end
        emit("O21/D86p", "pr (RLI, MEGNO on the pair map) Float64: recorder device columns at N = 300 against the verified device path (propagate on the PairMap with the stored xi_0 and the block frame; provenance gpu for $ng of $(length(items))): orbits bitwise equal / largest difference", "$neq of $ntot / $(fmt(maxd))", "all bitwise equal", pf(neq == ntot && ng == length(items)))
    end
end

const SALI_LOCK = ReentrantLock()
# ---- SALI at the device orbit against its BigFloat(256) reference (sali criterion 18(a), regular class: the orbits of the panel cells k2:T and k2:R, K = 0.5, wrapped, N = 4000)
section("O21") do
    items = NamedTuple[it for it in PLAN.items if it.method == "tg" && it.precision == "Float64" && it.coordinate_mode == "wrapped" && it.system == "standard" && it.cell in ("k2:T", "k2:R")]
    isempty(items) && error("no regular-class items")
    ctx = EX.make_ctx(PLAN, items)
    g = EX.device_records(ctx, items)
    g.recs === nothing && (emit("O21", "SALI reference gate", "FALLBACK($(g.prov.reason))", "device", "FAIL"); return)
    param = items[1].param
    all(it -> it.param == param, items) || error("items of one K expected")
    sysW = Chaoslab.WrappedStandardMap(Chaoslab.exact(Float64, param))
    Kt = Chaoslab.exact(Float64, param)
    x0 = [last(sys_x(it, Float64)) for it in items]
    N = 4000
    P, XN = dev_trace(sysW, x0, N, idproj)
    rL = Float64[]; rS = Float64[]; nres = 0; same = 0
    Threads.@threads :dynamic for j in eachindex(items)
        pts = [P[j, i] for i in 1:N]
        W0 = frame(items[j], Float64)
        ref = sl_ref(sl_mats_std(Kt, pts), W0, N; keep = [N])
        Ms = [SMatrix{2, 2, Float64}(map(Float64, Chaoslab.tangent(sysW, pt))) for pt in pts]
        p1, p2 = sl_psis(Ms, W0, [N])[1]
        at = sl_at(ref, N)
        sref = Float64(at.S)
        BS = sl_BSALI(epsT, sl_Lambda(ref.lnsali, N), sref, p1, p2)
        t = g.recs[j].indicators["tg"]
        L, S = colvec(t, "SALI.L")[end], colvec(t, "SALI.sali")[end]
        lock(SALI_LOCK) do
            push!(rL, sl_err(L, at.L) / sl_BL(epsT, p1, p2))
            if sref >= 10 * 4 * epsT
                push!(rS, sl_err(S, at.S) / BS)
            end
            same += (g.recs[j].state[end, 1] == XN[j][1] && g.recs[j].state[end, 2] == XN[j][2]) ? 1 : 0
        end
    end
    emit("O21/D86p", "SALI of the chain tg on the device, regular class (K = 0.5, $(length(items)) orbits, N = 4000): max |L - L_ref| / (8 eps (Psi1 + Psi2)) against BigFloat(256) at the device's own orbit points (sali criterion 18(a))", fmt(maximum(rL)), "<= 1", pf(maximum(rL) <= 1))
    emit("O21/D86p", "SALI: max |SALI - SALI_ref| / B_SALI(n) over the $(length(rS)) resolved orbits", isempty(rS) ? "none resolved" : fmt(maximum(rS)), "<= 1", pf(!isempty(rS) && maximum(rS) <= 1))
    emit("O21", "the trace kernel's final state equals the recorder's stored orbit point at the last grid step, bitwise", "$same of $(length(items))", "all", pf(same == length(items)))
end

# ---------------------------------------------------------------------------------------------------
section("DOUBLE") do
    for id in ("tg", "pr", "pa", "bk")
        items = NamedTuple[it for it in PLAN.items if it.method == id && it.precision == "Double64" && it.coordinate_mode == "wrapped"]
        ctx = EX.make_ctx(PLAN, items)
        g = EX.device_records(ctx, items)
        why = g.recs === nothing ? g.prov.reason : ""
        emit("D86p/R16", "$id, Double64, $(length(items)) items, device recording requested: the group falls back to the CPU with a recorded reason (the Double64 kernel is not offered)", g.recs === nothing ? "FALLBACK($why)" : "ran on the device", "FALLBACK with a reason of FALLBACK_REASONS",
             g.recs === nothing && why in FALLBACKS ? "FALLBACK" : "FAIL")
    end
end

# ---------------------------------------------------------------------------------------------------
# CAMPAIGN (D86y, G): the reduced campaign on the device. GC: the reference (the record in force has no per-spec entries, so auto resolves to the CPU); G1: auto with the per-spec records
# of this run (REC_SCRATCH); G2: backend request gpu for Float64 (every device-capable kind runs on the device).
function make_gpu_prereg()
    run = block_run()
    run["backend"]["Float64"] = "gpu"
    return make_prereg(joinpath(SCR, "prereg-reduced-gpu.md"), run)
end
section("CAMPAIGN") do
    cap = 1000
    stores = Dict("GC" => joinpath(SCRD, "store-GC"), "G1" => joinpath(SCRD, "store-G1"), "G2" => joinpath(SCRD, "store-G2"))
    mans = Dict(k => joinpath(SCRD, "manifest-$k.json") for k in keys(stores))
    for (k, s) in stores; rm(s; recursive = true, force = true); rm(mans[k]; force = true); end
    rec_scratch = isfile(REC_SCRATCH) ? EX.read_gpu_record(REC_SCRATCH, DEVNAME) : nothing
    ord_fail = first(c.ordinal for c in PLAN.chunks if c.run_type.method == "tg" && c.run_type.precision == "Float64" && c.run_type.coordinate_mode == "wrapped")
    prereg_gpu = make_gpu_prereg()
    function run1(tag, prereg; rec)
        t = @elapsed r = EX.run_campaign(prereg, PP.panel; store = stores[tag], manifest = mans[tag], threads = NT, thread_override = NT != K19, allow_dirty = true, capacity = cap,
                                         device_record = rec, hooks = EX.RunHooks(mutate = nan_fault(ord_fail, "tg"; attempts = (1, 2))))
        return r, t
    end
    rGC, tGC = run1("GC", PP.prereg; rec = EX.gpu_record())
    rG1, tG1 = run1("G1", PP.prereg; rec = rec_scratch)
    rG2, tG2 = run1("G2", prereg_gpu; rec = rec_scratch)
    for (tag, r, t) in (("GC", rGC, tGC), ("G1", rG1, tG1), ("G2", rG2, tG2))
        emit("D86y", "device campaign $tag: status, failed chunks (the injected one only), wall s", "$(r.status), $(r.failed), $(round(t; digits = 1)) s", "complete with failed chunks, [$ord_fail]", pf(r.status == "complete with failed chunks" && r.failed == [ord_fail]))
    end
    # backend per spec in each manifest
    function backends(tag)
        m = EX.read_json(mans[tag])
        d = Dict{String, Dict{String, Int}}()
        reasons = Dict{String, Set{String}}()
        for c in m["chunks"]
            id = PLAN.chunks[c["ordinal"]].run_type.method
            haskey(SPECS, id) || continue
            key = "$id/$(PLAN.chunks[c["ordinal"]].run_type.precision)"
            dd = get!(d, key, Dict{String, Int}())
            dd[c["backend"]] = get(dd, c["backend"], 0) + 1
            push!(get!(reasons, key, Set{String}()), string(c["fallback_reason"]))
        end
        return d, reasons
    end
    for tag in ("GC", "G1", "G2")
        d, rs = backends(tag)
        for key in sort(collect(keys(d)))
            say("BACKEND | campaign $tag | $key | chunks by backend $(d[key]) | reasons $(sort(collect(rs[key])))")
        end
    end
    d1, r1 = backends("G1"); d2, r2 = backends("G2")
    emit("D86y/R16", "campaign G2 (gpu request): chunks of tg, pr, pa, bk (Float64) ran on the device; st and fq (Float64) on the CPU with the sixth reason; Double64 on the CPU",
         join(("$k: $(d2[k])" for k in sort(collect(keys(d2))) if endswith(k, "Float64")), "; "), "gpu for tg, pr, pa, bk; cpu for st, fq",
         pf(all(k -> get(d2["$k/Float64"], "gpu", 0) > 0 && !haskey(d2["$k/Float64"], "cpu") || get(d2["$k/Float64"], "gpu", 0) > 0, ("tg", "pr", "pa", "bk")) && all(k -> !haskey(d2["$k/Float64"], "gpu") && r2["$k/Float64"] == Set([FALLBACKS[6]]), ("st", "fq"))))
    emit("D86y/R16", "campaign G1 (auto with the per-spec records at $NT threads): each spec's backend follows its record", join(("$k: $(d1[k])" for k in sort(collect(keys(d1))) if endswith(k, "Float64")), "; "), "follows the records",
         pf(all(id -> begin
                    r = rec_scratch === nothing ? nothing : get(rec_scratch.recorder, id, nothing)
                    want_gpu = r !== nothing && r.gpu_faster
                    dd = d1["$id/Float64"]
                    want_gpu ? (get(dd, "gpu", 0) > 0 && !haskey(dd, "cpu") || get(dd, "gpu", 0) > 0) : !haskey(dd, "gpu")
                end, ("tg", "pr", "pa", "bk")) && all(k -> !haskey(d1["$k/Float64"], "gpu") && r1["$k/Float64"] == Set([FALLBACKS[6]]), ("st", "fq"))))
    # the files: independent reader on the device campaign, and the module's comparison with the CPU reference campaign
    for tag in ("GC", "G2")
        rr = run_reader(stores[tag], "dev-$tag")
        R = parse_reader(rr.path)
        nfail = count(l -> occursin("|FAIL|", l), R.checks)
        emit("O20/D86y", "independent reader on the files of device campaign $tag: files, failed checks", "$(length(R.files)), $nfail", "$(length(PLAN.chunks)), 0", pf(rr.ok && nfail == 0 && length(R.files) == length(PLAN.chunks)))
    end
    gaps = Dict{String, Vector{Float64}}(); neq = Dict{String, Int}(); ntot = Dict{String, Int}()
    for ord in 1:length(PLAN.chunks)
        c = PLAN.chunks[ord]
        id = c.run_type.method
        haskey(SPECS, id) && c.run_type.precision == "Float64" || continue
        fa, fb = joinpath(stores["G2"], chunkname(ord)), joinpath(stores["GC"], chunkname(ord))
        a, b = EX.read_chunk(fa), EX.read_chunk(fb)
        (a.attrs["status"] == "ok" && b.attrs["status"] == "ok") || continue
        for (x, y) in zip(a.orbits, b.orbits)
            tx, ty = x.indicators[id], y.indicators[id]
            d = 0.0
            for (u, v) in ((tx.series, ty.series), (tx.final, ty.final), (tx.scalars, ty.scalars))
                u === nothing && continue
                d = max(d, maximum((abs(p - q) for (p, q) in zip(u, v) if isfinite(p) && isfinite(q)); init = 0.0))
            end
            push!(get!(gaps, id, Float64[]), d)
            ntot[id] = get(ntot, id, 0) + 1
            neq[id] = get(neq, id, 0) + (d == 0 && isbitseq(x.state, y.state) ? 1 : 0)
        end
    end
    for id in sort(collect(keys(gaps)))
        g = gaps[id]
        info("D86y/O21", "device campaign G2 against the CPU campaign GC, spec $id: orbits $(ntot[id]), bitwise equal tables and orbit points $(neq[id]), largest |difference| of a stored table value $(fmt(maximum(g))), median $(fmt(sort(g)[max(div(end, 2), 1)])) (orbit-level gap: reported, not gated)")
    end
    # the cpu-only kinds are bitwise equal across the two campaigns (same code path)
    emit("D86y/R16", "st and fq (CPU only) are bitwise equal between the device campaign G2 and the CPU campaign GC", "st $(neq["st"]) of $(ntot["st"]); fq $(neq["fq"]) of $(ntot["fq"])", "all", pf(neq["st"] == ntot["st"] && neq["fq"] == ntot["fq"]))
    # the manifest of G2 validates: the device, the backend table rows of the specs, the flag counts
    m2 = EX.read_json(mans["G2"])
    info("D86y", "manifest G2: device $(m2["device"]), threads $(m2["threads"]), backend table rows for specs: " * join(("$(r["precision"]) $(get(r, "method", get(r, "spec", "?"))) -> $(r["backend"])" for r in m2["backend_table"] if get(r, "method", get(r, "spec", "")) in keys(SPECS)), "; "))
    PARC = (gaps = gaps, neq = neq, ntot = ntot)
    global CAMPAIGN_GAPS = PARC
end

# ---------------------------------------------------------------------------------------------------
section("PLOT") do
    ids = ["tg", "pr", "st", "fq", "pa", "bk"]
    labels = Dict("tg" => "tangent\nFLI,MEGNO,SALI", "pr" => "pair\nRLI,MEGNO", "st" => "stretching", "fq" => "frequency", "pa" => "phase advance", "bk" => "birkhoff")
    fig = CM.Figure(size = (1280, 520), fontsize = 13)
    ax1 = CM.Axis(fig[1, 1], title = "Device against CPU recording, same items (Float64)", ylabel = "largest |difference| of a stored table value per orbit (0 drawn at 1e-18)",
                  yscale = log10, xticks = (1:length(ids), [labels[i] for i in ids]), xticklabelrotation = 0.0)
    xs = Float64[]; ys = Float64[]
    for (k, id) in enumerate(ids)
        haskey(PAR, id) || continue
        g = PAR[id].gaps
        append!(xs, fill(Float64(k), length(g)) .+ 0.25 .* (rand(MersenneTwister(k), length(g)) .- 0.5))
        append!(ys, max.(g, 1e-18))
    end
    CM.scatter!(ax1, xs, ys; markersize = 4, color = (:steelblue, 0.5))
    for (k, id) in enumerate(ids)
        haskey(PAR, id) || (CM.text!(ax1, k, 1e-18; text = "CPU only", align = (:center, :bottom), fontsize = 12); continue)
        CM.text!(ax1, k, 3e-18; text = "$(PAR[id].neq)/$(length(PAR[id].gaps))
bitwise", align = (:center, :bottom), fontsize = 10)
    end
    CM.ylims!(ax1, 5e-19, 1e5)
    ax2 = CM.Axis(fig[1, 2], title = "Cost of the recording, $ENS_ORB orbits (ns per orbit-step)", ylabel = "ns per orbit-step", yscale = log10, limits = (nothing, (1.0, 1e5)),
                  xticks = (1:length(ids), [labels[i] for i in ids]), xticklabelrotation = 0.0)
    for (k, id) in enumerate(ids)
        haskey(COST, id) || continue
        c = COST[id]
        CM.barplot!(ax2, [k - 0.2], [c.cpu * NT]; width = 0.35, fillto = 1.0, color = :gray55, label = k == 1 ? "CPU, $NT threads (thread-ns)" : nothing)
        isnan(c.gpu) || CM.barplot!(ax2, [k + 0.2], [c.gpu]; width = 0.35, fillto = 1.0, color = :seagreen, label = k == 1 ? "device (wall ns)" : nothing)
        isnan(c.gpu) && CM.text!(ax2, k + 0.2, c.cpu * NT / 3; text = "no device path", rotation = pi / 2, fontsize = 11, align = (:center, :center))
    end
    CM.axislegend(ax2; position = :lt)
    path = joinpath(OUT, "experiments-d86-device.png")
    CM.save(path, fig)
    info("D86y", "figure written: $path")
end

say("# data/gpu-benchmark.toml sha256 after: $(bytes2hex(open(sha256, REAL_RECORD))) (before $REAL_RECORD_SHA): ", bytes2hex(open(sha256, REAL_RECORD)) == REAL_RECORD_SHA ? "unchanged" : "CHANGED")
say("# summary: PASS $(COUNTS["PASS"]) FAIL $(COUNTS["FAIL"]) FALLBACK $(COUNTS["FALLBACK"]) INFO $(COUNTS["INFO"]), peak RSS $(round(Sys.maxrss() / 2^30; digits = 2)) GiB")
close(LOG)
