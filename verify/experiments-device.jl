# Stage 04 device check of the D53 runner (device path of the recorder), local GPU only. Run from the repo root:
#   julia --project=verify -t 23 verify/experiments-device.jl > workspaces/module-build/stages/04-verify/output/gpu-rtx5060-d53-verify-log.txt
# EXD_ONLY=table prints only the backend table and the record in force (used before and after the AC13 record is retaken).
# For each of the eight kernel-path methods (L1, L2, L3, L4, L4b, L5, L5b, L9) at Float64, lifted (all panel orbits, run type of arm
# "main") and wrapped (arm "wrapped"), horizon scaled down to EXD_N steps: the recording on the device against the CPU recording,
# chaotic cells (classes C and X) by the form of tangent-core AC12(c) with the values of [run.verify].backend_agreement, the launch
# provenance, a relaunch at a launch length that is not a grid step, and the cost in ns per orbit-step on the device and on the CPU threads.
using Chaoslab, CUDA, TOML, Statistics, Printf
const EX = Chaoslab.Experiments
const F = Chaoslab.Foundation
const ROOT = normpath(joinpath(@__DIR__, ".."))
const PRE = joinpath(ROOT, "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")
const PAN = joinpath(ROOT, "workspaces", "study", "stages", "02-panel", "output", "p1-le-benchmark-panel.toml")
const N = parse(Int, get(ENV, "EXD_N", "20000"))
const NT = Threads.nthreads()

say(args...) = println(args...)
fmt(x) = @sprintf("%.4g", x)

say("# D53 device check, julia ", VERSION, ", threads ", NT, ", cpu ", Sys.CPU_NAME, ", horizon ", N)
say("# device ", typeof(F.DEVICE[]), ", ", CUDA.name(CUDA.device()), ", CUDA.jl ", pkgversion(CUDA), ", runtime ", CUDA.runtime_version(),
    ", free ", round(CUDA.free_memory() / 2^20), " MiB, stack limit ", CUDA.limit(CUDA.LIMIT_STACK_SIZE), " B")
F.DEVICE[] isa F.NoDevice && (say("# no functional device"); exit(1))
plan = EX.plan_run(PRE, PAN)
say("# plan ", plan.counts, ", [run.k19].threads = ", plan.block.raw["k19"]["threads"], ", [run.verify] = ", plan.block.raw["verify"],
    ", [run.budget].Float64 = ", plan.block.raw["budget"]["Float64"])
rec = EX.gpu_record()
say("# record in force (data/gpu-benchmark.toml, Float64): ", rec)
tbl = EX.backend_table(plan; threads = NT, record = rec)
say("# backend table at ", NT, " threads (precision, method: request, backend, device, reason, note)")
for ((p, m), r) in sort!(collect(tbl); by = first)
    say("TABLE | ", p, " | ", m, " | ", r.request, " | ", r.backend, " | ", r.device, " | ", r.reason, " | ", r.note)
end
get(ENV, "EXD_ONLY", "") == "table" && exit(0)

const AGREE = plan.block.raw["verify"]["backend_agreement"]
const METHODS = ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L9"]

# the items of the run type (method, Float64, lifted or wrapped, tau 1) of the main arm, horizon scaled to N
function group_items(m, mode)
    arm = mode == "lifted" ? "main" : "wrapped"
    its = [plan.items[i] for c in plan.chunks for i in c.items if c.run_type.method == m && c.run_type.precision == "Float64" &&
           c.run_type.coordinate_mode == mode && c.run_type.tau == 1 && c.run_type.horizon == plan.block.raw["arm"][findfirst(a -> a["name"] == arm, plan.block.raw["arm"])]["horizons"][1] && arm in c.arms]
    return NamedTuple[merge(it, (horizon = N,)) for it in its]
end

same(a, b) = isequal(a.lyap, b.lyap) && isequal(a.state, b.state) && isequal(a.diagnostics, b.diagnostics) && isequal(a.winding, b.winding)

# AC12(c) form per ensemble (cell): max over grid steps of the per-orbit |dlambda_1| / sigma_FT (sigma over the CPU ensemble) and of
# |difference of ensemble means| / combined standard error
function ac12c(cpu, dev, items, m)
    worst_orbit, worst_mean, nens = 0.0, 0.0, 0
    for cell in unique(it.cell for it in items)
        idx = [i for (i, it) in enumerate(items) if it.cell == cell && it.design_class in ("C", "X")]
        length(idx) < 8 && continue
        nens += 1
        Lc = reduce(hcat, [cpu[i].lyap[m][:, 1] for i in idx])'
        Ld = reduce(hcat, [dev[i].lyap[m][:, 1] for i in idx])'
        for j in axes(Lc, 2)
            sig = std(Lc[:, j])
            sig > 0 || continue
            worst_orbit = max(worst_orbit, maximum(abs.(Lc[:, j] .- Ld[:, j])) / sig)
            se = sqrt((var(Lc[:, j]) + var(Ld[:, j])) / length(idx))
            se > 0 && (worst_mean = max(worst_mean, abs(mean(Lc[:, j]) - mean(Ld[:, j])) / se))
        end
    end
    return worst_orbit, worst_mean, nens
end

function bestof(f, reps)
    best, out = Inf, nothing
    for _ in 1:reps
        t = @elapsed out = f()
        best = min(best, t)
    end
    return best, out
end

say("# format: ROW | mode | method | orbits | steps | ... ")
for mode in ("lifted", "wrapped"), m in METHODS
    items = group_items(m, mode)
    isempty(items) && (say("ROW | $mode | $m | no items"); continue)
    ctx = EX.make_ctx(plan, items)
    steps = length(items) * N
    EX.device_records(ctx, items; hook = false)                      # compile
    tg, g = bestof(() -> EX.device_records(ctx, items), 3)
    if g.recs === nothing
        say("ROW | $mode | $m | $(length(items)) | $N | FALLBACK(", g.prov.reason, ")")
        continue
    end
    EX.compute_records(ctx, items, NT)                               # compile the CPU path
    tc, cpu = bestof(() -> EX.compute_records(ctx, items, NT), 2)
    eq = count(same.(g.recs, cpu))
    maxd = maximum(maximum(abs.(g.recs[i].lyap[m] .- cpu[i].lyap[m])) for i in eachindex(items))
    wo, wm, ne = ac12c(cpu, g.recs, items, m)
    ok = wo <= AGREE["per_orbit_sigma"] && wm <= AGREE["mean_standard_errors"]
    # relaunch: launches of 777 steps (not a grid step): the device record must equal the device record of one launch
    old = Chaoslab.Tangent._LAUNCH_STEPS[]
    Chaoslab.Tangent._LAUNCH_STEPS[] = 777
    gr = try
        EX.device_records(ctx, items)
    finally
        Chaoslab.Tangent._LAUNCH_STEPS[] = old
    end
    releq = gr.recs !== nothing && all(same.(gr.recs, g.recs))
    say("ROW | $mode | $m | $(length(items)) | $N | identical series $eq of $(length(items)) | max |dlambda| ", fmt(maxd), " | AC12c orbit/sigma ", fmt(wo), " mean/SE ", fmt(wm),
        " over $ne ensembles | ", ok ? "PASS" : "FAIL", " | relaunch 777 equal ", releq, " launches ", gr.prov.launches, " | device ns/orbit-step ", fmt(tg / steps * 1e9),
        " | cpu ", NT, " threads ns/orbit-step ", fmt(tc / steps * 1e9), " | thread ns ", fmt(tc / steps * 1e9 * NT))
    say("PROV | $mode | $m | ", g.prov)
end

# Cost at the ensemble size of the campaign: ENS_ORB orbits (the main-arm items repeated) x ENS_N steps, lifted, tau 1, ns per orbit-step end to end.
# (a) the launch plan of tangent-core (L from the one-wave value of the record, batches capped by the same budget), (b) launches of 10^4 steps
# (the batch then holds the whole ensemble), (c) the CPU threads on the same items; against [run.budget].Float64 (ns_cpu_thread, ns_gpu).
const ENS_ORB = parse(Int, get(ENV, "EXD_ENS_ORB", "2048"))
const ENS_N = parse(Int, get(ENV, "EXD_ENS_N", "200000"))
bg = plan.block.raw["budget"]["Float64"]
say("# ensemble cost, $ENS_ORB orbits x $ENS_N steps; [run.budget].Float64 = $bg")
for m in METHODS
    base = group_items(m, "lifted")
    items = NamedTuple[merge(base[mod1(i, length(base))], (horizon = ENS_N,)) for i in 1:ENS_ORB]
    ctx = EX.make_ctx(plan, items)
    steps = ENS_ORB * ENS_N
    EX.device_records(ctx, items; hook = false)
    ta, ga = bestof(() -> EX.device_records(ctx, items), 2)
    old = Chaoslab.Tangent._LAUNCH_STEPS[]
    Chaoslab.Tangent._LAUNCH_STEPS[] = 10^4
    tb, gb = try
        bestof(() -> EX.device_records(ctx, items), 2)
    finally
        Chaoslab.Tangent._LAUNCH_STEPS[] = old
    end
    tc, _ = bestof(() -> EX.compute_records(ctx, items, NT), 1)
    say("ENS | $m | (a) tangent-core plan ", ga.recs === nothing ? "FALLBACK($(ga.prov.reason))" : "$(fmt(ta / steps * 1e9)) ns per orbit-step, launches $(ga.prov.launches) of $(ga.prov.launch_steps) steps, batches $(ga.prov.gpu_chunks) of $(ga.prov.chunk_size)",
        " | (b) launches of 1e4 ", gb.recs === nothing ? "FALLBACK($(gb.prov.reason))" : "$(fmt(tb / steps * 1e9)) ns per orbit-step, batches $(gb.prov.gpu_chunks) of $(gb.prov.chunk_size)",
        " | (c) cpu $NT threads ", fmt(tc / steps * 1e9), " ns per orbit-step (", fmt(tc / steps * 1e9 * NT), " thread-ns)")
end
