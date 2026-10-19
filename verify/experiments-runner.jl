# Stage 04 verification of the experiments runner and run store (Amendments D53 and D54). Run from the repo root:
#   julia --project=verify -t 23 verify/experiments-runner.jl
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md (Amendment D53: O13, O14, O15, D53b [full], D53f, D53h, D53i,
# D53m, D53q, D53s, D53t, D53w, D53x; Amendment D54: D54f). Sections are named; EXR_ONLY="PLAN,CAMPAIGN" runs only those (they reuse
# stores already in the scratch folder when EXR_REUSE=1). EXR_QUICK=1 uses tiny horizons (debugging aid, not a verification result).
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# Scratch folder (all stores, reduced preregistrations): EXR_SCRATCH. The study data folder data/p1-le-benchmark is never written (D53x).
using Printf, Random, LinearAlgebra, Statistics, SHA, Dates
using Chaoslab, CUDA, TOML, StaticArrays, HDF5_jll
import CairoMakie as CM

const EX = Chaoslab.Experiments
const F = Chaoslab.Foundation
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const PREREG = joinpath(ROOT, "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")
const PANEL = joinpath(ROOT, "workspaces", "study", "stages", "02-panel", "output", "p1-le-benchmark-panel.toml")
const QUICK = get(ENV, "EXR_QUICK", "") == "1"
const REUSE = get(ENV, "EXR_REUSE", "") == "1"
const ONLY = filter(!isempty, split(get(ENV, "EXR_ONLY", ""), ","))
const SCR = joinpath(get(ENV, "EXR_SCRATCH", raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\runner-verify"), QUICK ? "quick" : "full")
const NT = Threads.nthreads()
const VERIFY_ENV = joinpath(ROOT, "verify")
const STUDY_DATA = joinpath(ROOT, "data", "p1-le-benchmark")
const MANIFEST_DIR = joinpath(ROOT, "workspaces", "study", "stages", "03-run", "output")
mkpath(OUT); mkpath(SCR)
const LOG = open(get(ENV, "EXR_LOG", joinpath(OUT, "experiments-runner-verify-log.txt")), "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const COUNTS = Dict("PASS" => 0, "FAIL" => 0, "INFO" => 0)
function say(args...)
    s = string(args...)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(id, quantity, value, bound, verdict)
    v = startswith(verdict, "FAIL") ? "FAIL" : startswith(verdict, "PASS") ? "PASS" : "INFO"
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
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 900)], "-", "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end
lf(s) = replace(s, "\r\n" => "\n")
sha_file(p) = bytes2hex(open(sha256, p))
bits(A) = collect(reinterpret(UInt8, vec(A)))
rawsha(A) = bytes2hex(sha256(bits(A)))
function filecount(dir)
    isdir(dir) || return 0
    return sum((length(fs) for (_, _, fs) in walkdir(dir)); init = 0)
end
function tree_listing(dir)
    isdir(dir) || return String[]
    return sort!([relpath(joinpath(r, f), dir) * "|" * string(filesize(joinpath(r, f))) * "|" * string(mtime(joinpath(r, f))) for (r, _, fs) in walkdir(dir) for f in fs])
end

const T0 = time()
const CAMP = Dict{String, Any}()
say("# experiments runner verification, julia $(VERSION), threads $NT, cpu $(Sys.CPU_NAME), HDF5.jl $(pkgversion(Chaoslab.Experiments.HDF5)), $(QUICK ? "QUICK MODE (not a verification result)" : "full")")
say("# prereg $(relpath(PREREG, ROOT)), panel $(relpath(PANEL, ROOT)), scratch $SCR")
const GIT_HEAD = strip(read(`git -C $ROOT rev-parse HEAD`, String))
const GIT_DIRTY = strip(read(`git -C $ROOT status --short`, String))
say("# git head $GIT_HEAD, tracked or untracked changes: ", isempty(GIT_DIRTY) ? "none" : replace(GIT_DIRTY, "\n" => "; "))

# ---------------------------------------------------------------------------------------------------
# D53x before: the study data folder
const DATA_BEFORE = (filecount(STUDY_DATA), tree_listing(joinpath(ROOT, "data")), tree_listing(MANIFEST_DIR))
say("# D53x before: files under data/p1-le-benchmark = $(DATA_BEFORE[1]); data/ entries = $(length(DATA_BEFORE[2])); study manifest dir entries = $(length(DATA_BEFORE[3]))")

# ---------------------------------------------------------------------------------------------------
# Reduced preregistrations (verification inputs, not preregistrations): the real file with its [run] block rewritten and the header hash recomputed
const HZ = QUICK ? Dict("main" => 3200, "long" => 6400, "tau_sweep" => 3200, "wrapped" => 3200, "directions" => 3200, "eps_sweep" => 3200,
                        "h4_float64" => 12800, "h4_double64" => 3200, "ladder_double64_short" => 3200, "ladder_double64_long" => 6400, "ladder_bigfloat" => 800,
                        "main_l2" => 3200, "long_l2" => 6400) :
                 Dict("main" => 100_000, "long" => 1_000_000, "tau_sweep" => 100_000, "wrapped" => 100_000, "directions" => 100_000, "eps_sweep" => 100_000,
                      "h4_float64" => 10_000_000, "h4_double64" => 1_000_000, "ladder_double64_short" => 100_000, "ladder_double64_long" => 500_000,
                      "ladder_bigfloat" => 10_000, "main_l2" => 100_000, "long_l2" => 1_000_000)

function run_toml_blocks(text)
    out = []
    for m in eachmatch(r"```toml\n(.*?)\n```"s, text)
        d = try TOML.parse(String(m.captures[1])) catch; continue end
        haskey(d, "run") && push!(out, (m, d))
    end
    isempty(out) && error("no [run] block")
    return out
end
# mod! rewrites the first [run] block (A4). keep_later = true applies it to every later [run] block as well (A5), for reductions that must reach
# all arms; otherwise the later blocks are emptied so that a block written wholesale by mod! is not merged with the real A5 arms.
function make_prereg(path, mod!; keep_later = false)
    text = lf(read(PREREG, String))
    new = text
    for (i, (m, d)) in enumerate(run_toml_blocks(text))
        if i == 1 || keep_later
            run = d["run"]; mod!(run)
            io = IOBuffer(); TOML.print(io, Dict("run" => run); sorted = true)
            blk = "```toml\n" * rstrip(String(take!(io))) * "\n```"
        else
            blk = "```toml\n# emptied in the verification copy\n```"
        end
        new = replace(new, String(m.match) => blk; count = 1)
    end
    h = EX.prereg_hash(new)
    new = replace(new, r"(?m)^(\| SHA-256 of this file[^|]*\|[ \t]*)[0-9a-f]{64}" => SubstitutionString("\\g<1>" * h); count = 1)
    write(path, new)
    EX.check_prereg(path)
    return path
end
# control_bound: the real rule 1e-8 at N1 = 1e7 scaled by 1 / n is CB_SCALED. Under A5 (control_scope = "exact_start") only the controls whose binary64
# start is an exact periodic point and Y5 are gated, so the campaign blocks use the scaled rule itself; CB_TRIP (far below the finite-time deviation of the
# gated controls at N = 1e5) exercises the stop in its own run. CB_OPEN = 1 stays for the custom blocks of the extended and device sections.
const CB_SCALED = 1e-8 * 1e7 / HZ["main"]
const CB_TRIP = 1e-9
const CB_OPEN = 1.0
function reduce!(run, cb)
    for arm in run["arm"]
        arm["horizons"] = [HZ[arm["name"]]]
    end
    h = get!(run, "health", Dict{String, Any}())
    haskey(h, "control_bound") && (h["control_bound"] = cb)
    haskey(h, "stop_failed_fraction") && (h["stop_failed_fraction"] = 0.5)      # opened so that natural failures cannot stop the campaign; they are counted against the real 0.05 below
    haskey(h, "first_chunk_cost_factor") && (h["first_chunk_cost_factor"] = 1e6) # the cost cap is judged from the measured costs below, not enforced here
    haskey(run, "pilot") && (run["pilot"]["horizon"] = minimum(values(HZ)))
    return run
end
reduce_all!(run, cb) = haskey(run, "arm") ? reduce!(run, cb) : run
const PRE_R = make_prereg(joinpath(SCR, "prereg-reduced.md"), r -> reduce_all!(r, QUICK ? CB_OPEN : CB_SCALED); keep_later = true)   # quick horizons are too short for the scaled rule
const PRE_S = make_prereg(joinpath(SCR, "prereg-reduced-scaledbound.md"), r -> reduce_all!(r, CB_TRIP); keep_later = true)
const PLAN_R = EX.plan_run(PRE_R, PANEL)
const PLAN_REAL = EX.plan_run(PREREG, PANEL)
const CAPACITY = 2000

# ---------------------------------------------------------------------------------------------------
# Kill and fault hooks
struct VKill <: Exception
    point::Symbol
    n::Int
end
has_kill(::VKill) = true
has_kill(e::CompositeException) = any(has_kill, e.exceptions)
has_kill(e::TaskFailedException) = has_kill(e.task.result)
has_kill(e::LoadError) = has_kill(e.error)
has_kill(::Any) = false
function killed(f)
    try
        f()
    catch e
        has_kill(e) && return true
        rethrow()
    end
    return false
end
kill_after(k) = (p, info) -> (p === :chunk_done && info.chunks_done == k && throw(VKill(p, k)); nothing)
kill_during(ord) = (p, info) -> (p === :during_data && info.ordinal == ord && throw(VKill(p, ord)); nothing)
function nan_fault(ord; attempts)
    return (info, orbits) -> begin
        info.ordinal == ord && info.attempt in attempts && (first(values(orbits[min(2, end)].lyap))[min(3, end), 1] = NaN)
        nothing
    end
end
const ORD_ONCE = first(c.ordinal for c in PLAN_R.chunks if c.run_type.method == "L5" && c.run_type.precision == "Float64" && "main" in c.arms)
const ORD_TWICE = first(c.ordinal for c in PLAN_R.chunks if c.run_type.method == "L9" && "directions" in c.arms)
inject() = (info, orbits) -> (nan_fault(ORD_ONCE; attempts = (1,))(info, orbits); nan_fault(ORD_TWICE; attempts = (1, 2))(info, orbits))

K19 = PLAN_R.block.raw["k19"]["threads"]
function campaign(prereg, store, manifest; hooks = EX.RunHooks(), kw...)
    return EX.run_campaign(prereg, PANEL; store = store, manifest = manifest, threads = NT, thread_override = NT != K19, allow_dirty = true,
                           capacity = CAPACITY, hooks = hooks, kw...)
end

# ---------------------------------------------------------------------------------------------------
# Independent reader (separate process, no Chaoslab) and the comparison with the module's reader
const READER = joinpath(@__DIR__, "experiments-runner-reader.jl")
function run_reader(store, tag)
    out = joinpath(SCR, "reader-$tag.txt")
    t = @elapsed ok = success(pipeline(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t 1 $READER $store $out`; stdout = devnull, stderr = joinpath(SCR, "reader-$tag.err")))
    return (path = out, ok = ok, seconds = t)
end
function parse_reader(path)
    R = (attrs = Dict{Tuple{String, String}, String}(), oattrs = Dict{Tuple{String, String, String}, String}(),
         data = Dict{Tuple{String, String, String}, NamedTuple}(), dattrs = Dict{Tuple{String, String, String, String}, String}(),
         ksha = Dict{Tuple{String, String, String}, String}(), checks = String[], files = Set{String}())
    for line in eachline(path)
        isempty(line) && continue
        k = line[1]
        if k == 'A'
            f = split(line, '|'; limit = 6); R.attrs[(f[2], f[4])] = f[6]; push!(R.files, f[2])
        elseif k == 'O'
            f = split(line, '|'; limit = 6); R.oattrs[(f[2], f[3], f[4])] = f[6]
        elseif k == 'D'
            f = split(line, '|'); R.data[(f[2], f[3], f[4])] = (eltype = f[5], cdims = f[6], jdims = f[7], sha = f[8])
        elseif k == 'X'
            f = split(line, '|'; limit = 6); R.dattrs[(f[2], f[3], f[4], f[5])] = f[6]
        elseif k == 'K'
            f = split(line, '|'); R.ksha[(f[2], f[3], f[4])] = f[5]
        elseif k == 'C'
            push!(R.checks, line)
        end
    end
    return R
end
attr_str(v) = replace(v isa AbstractArray ? join(string.(v), ";") : string(v), '\n' => "\\n")
chunkname(ord) = @sprintf("chunk-%05d.h5", ord)

# Module reader against the independent reader: every attribute, every dataset (sha-256 of the raw bytes) and its dtype
function compare_readers(store, R)
    nfile = nattr = ndata = nbad = 0
    bad = String[]
    for fn in sort(readdir(store))
        occursin(r"^chunk-\d{5}\.h5$", fn) || continue
        nfile += 1
        ch = EX.read_chunk(joinpath(store, fn))
        for (nm, v) in ch.attrs
            nattr += 1
            got = get(R.attrs, (fn, nm), nothing)
            got == attr_str(v) || (nbad += 1; length(bad) < 8 && push!(bad, "$fn attr $nm: module $(attr_str(v)) reader $got"))
        end
        for o in ch.orbits
            oid = string(o.id)
            for (nm, v) in o.attrs
                nattr += 1
                got = get(R.oattrs, (fn, oid, nm), nothing)
                got == attr_str(v) || (nbad += 1; length(bad) < 8 && push!(bad, "$fn orbit $oid attr $nm: module $(attr_str(v)) reader $got"))
            end
            sets = Pair{String, Any}["t" => o.t, "state" => o.state]
            o.winding === nothing || push!(sets, "winding" => o.winding)
            for (m, v) in o.lyap; push!(sets, "lyap/$m" => v); end
            for (m, v) in o.diagnostics; push!(sets, "diagnostics/$m" => v); end
            for (m, v) in o.cost; push!(sets, "cost/$m" => v); end
            for (p, A) in sets
                ndata += 1
                d = get(R.data, (fn, oid, p), nothing)
                ok = d !== nothing && d.sha == rawsha(A) && d.eltype == string(eltype(A)) && d.jdims == join(size(A), "x")
                ok || (nbad += 1; length(bad) < 8 && push!(bad, "$fn orbit $oid dataset $p differs between readers"))
            end
            for (m, thr) in o.thresholds
                ndata += 1
                string(thr) == get(R.dattrs, (fn, oid, "diagnostics/$m", "threshold"), "") &&
                    string(o.first_exceedance[m]) == get(R.dattrs, (fn, oid, "diagnostics/$m", "first_exceedance_step"), "") ||
                    (nbad += 1; length(bad) < 8 && push!(bad, "$fn orbit $oid diagnostics/$m attributes differ"))
            end
        end
    end
    return (files = nfile, attrs = nattr, datasets = ndata, bad = nbad, examples = bad)
end

# Two stores, by their reader dumps: differences outside the allowed set (times, host fields). ignore = root attribute names left out.
function diff_readers(Ra, Rb; ignore = ("created_utc", "cpu", "threads", "session"))
    diffs = Dict{String, Int}()
    ex = String[]
    note(k, msg) = (diffs[k] = get(diffs, k, 0) + 1; length(ex) < 6 && push!(ex, msg))
    for ((f, nm), v) in Ra.attrs
        nm in ignore && continue
        get(Rb.attrs, (f, nm), nothing) == v || note("attr:" * nm, "$f $nm: $v vs $(get(Rb.attrs, (f, nm), nothing))")
    end
    for (k, v) in Ra.oattrs
        get(Rb.oattrs, k, nothing) == v || note("oattr:" * k[3], "$(k[1]) orbit $(k[2]) $(k[3])")
    end
    for (k, v) in Ra.data
        startswith(k[3], "cost/") && continue
        d = get(Rb.data, k, nothing)
        (d !== nothing && d.sha == v.sha && d.eltype == v.eltype && d.cdims == v.cdims) || note("data:" * k[3], "$(k[1]) orbit $(k[2]) $(k[3])")
    end
    for (k, v) in Ra.ksha
        get(Rb.ksha, k, nothing) == v || note("cost-entries", "$(k[1]) orbit $(k[2]) $(k[3])")
    end
    for (k, v) in Ra.dattrs
        get(Rb.dattrs, k, nothing) == v || note("dattr:" * k[4], "$(k[1]) orbit $(k[2]) $(k[3]) $(k[4])")
    end
    # key sets
    for (S, T, nm) in ((Ra.attrs, Rb.attrs, "attr-set"), (Rb.attrs, Ra.attrs, "attr-set"), (Ra.data, Rb.data, "data-set"), (Rb.data, Ra.data, "data-set"))
        for k in keys(S)
            haskey(T, k) || note(nm, "missing key $(k)")
        end
    end
    return (diffs = diffs, examples = ex, total = sum(values(diffs); init = 0))
end

# ---------------------------------------------------------------------------------------------------
# Manifest validation against R7 (test plan "Manifest contract")
const R7_KEYS = ["study", "schema_version", "plan_sha256", "plan_counts", "status", "reasons", "sessions", "git", "versions", "host", "threads", "backend_table",
                 "device", "device_record", "prior_record", "prereg", "panel", "block", "defaults", "launch", "chunks", "arms", "wall", "overrides", "removed_temporaries"]
const R7_CHUNK = ["ordinal", "key", "hash", "arms", "status", "attempt", "backend", "fallback_reason", "orbits", "steps", "wall_seconds", "first_attempt_reasons", "flags", "setup_check"]
const R7_ARM = ["steps", "wall_seconds", "thread_seconds", "cost_ns_per_orbit_step", "budget_ns", "status"]
function validate_manifest(path, plan; status_ok)
    m = EX.read_json(path)
    miss = [k for k in R7_KEYS if !haskey(m, k)]
    chunks = get(m, "chunks", Any[])
    missc = isempty(chunks) ? ["no chunks"] : [k for c in chunks for k in R7_CHUNK if !haskey(c, k)]
    missa = [k for (_, a) in get(m, "arms", Dict()) for k in R7_ARM if !haskey(a, k)]
    nchunks = length(plan.chunks)
    sts = Dict{String, Int}()
    for c in chunks; sts[c["status"]] = get(sts, c["status"], 0) + 1; end
    pre, pan = get(m, "prereg", Dict()), get(m, "panel", Dict())
    ok = (miss = miss, missing_chunk_keys = unique(missc), missing_arm_keys = unique(missa), n_chunks = length(chunks), plan_chunks = nchunks,
          sum_ok = sum(values(sts); init = 0) == nchunks, order_ok = [c["ordinal"] for c in chunks] == collect(1:nchunks), statuses = sts,
          prereg_stable = get(pre, "sha256_start", "a") == get(pre, "sha256_end", "b") == get(pre, "header", "c"),
          panel_stable = get(pan, "sha256_start", "a") == get(pan, "sha256_end", "b"), status = get(m, "status", ""),
          plan_hash = get(m, "plan_sha256", "") == plan.hash, schema = get(m, "schema_version", 0) == 2,
          git40 = occursin(r"^[0-9a-f]{40}$", string(get(get(m, "git", Dict()), "sha", ""))),
          dirty_recorded = haskey(get(m, "overrides", Dict()), "dirty"), sessions = length(get(m, "sessions", [])))
    return m, ok
end

# ---------------------------------------------------------------------------------------------------
section("SETUP") do
    emit("D53w", "reduced prereg written; hash check E18 passes; horizons (arm: N)", join(("$k: $v" for (k, v) in sort(collect(HZ); by = first)), ", "), "-", "INFO")
    emit("D53w", "reduced plan: items, chunks (equals the full plan: only horizons differ)", "$(PLAN_R.counts) vs $(PLAN_REAL.counts)", "equal", pf(PLAN_R.counts == PLAN_REAL.counts))
    st = PLAN_R.steps
    f64, d64, bf = st.precisions["Float64"].steps_after, st.precisions["Double64"].steps_after, st.precisions["BigFloat"].steps_after
    emit("D53w", "planned unique orbit-steps Float64 / Double64 / BigFloat", @sprintf("%.3g / %.3g / %.3g", f64, d64, bf), "Float64 < 1e10, Double64 < 1e8",
         pf(QUICK || (f64 < 1e10 && d64 < 1e8)))
    emit("D53w", "reduced [run.health]: control_bound of the campaign block / of the scaled-rule block (real 1e-8 at N1 = 1e7, scaled by 1/n) / first_chunk_cost_factor", "$(PLAN_R.block.raw["health"]["control_bound"]) / $(CB_SCALED) / $(PLAN_R.block.raw["health"]["first_chunk_cost_factor"])", "-", "INFO")
    emit("D53w", "injected faults: NaN at attempt 1 in chunk (rerun passes) / NaN at attempts 1 and 2 in chunk (final failed)", "$ORD_ONCE / $ORD_TWICE", "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# O14 and D53b [full]
section("PLAN") do
    panel = TOML.parsefile(PANEL)["orbit"]
    nP, nSD, nSL = length(panel), count(o -> o["in_SD"], panel), count(o -> o["in_SL"], panel)
    run = PLAN_REAL.block.raw
    nh4 = length(run["h4"])
    emit("D53b[full]", "real panel subsets: P, S_D, S_L, h4 starts", "$nP, $nSD, $nSL, $nh4", "883, 172, 92, 2", pf((nP, nSD, nSL, nh4) == (883, 172, 92, 2)))
    # independent Cartesian product and merge from the TOML block and the panel flags (no plan code)
    prim = run["primary"]
    subset_ids(name) = name == "P" ? [Int(o["id"]) for o in panel] : name == "S_D" ? [Int(o["id"]) for o in panel if o["in_SD"]] :
                       name == "S_L" ? [Int(o["id"]) for o in panel if o["in_SL"]] : [nP + i for i in 1:nh4]
    seen = Dict{Tuple, String}()
    before = Dict{String, Int}(); after = Dict{String, Int}(); stepsb = Dict{String, BigInt}(); stepsa = Dict{String, BigInt}()
    prec_b = Dict{String, BigInt}(); prec_a = Dict{String, BigInt}()
    for arm in run["arm"]
        nm = arm["name"]; before[nm] = 0; after[nm] = 0; stepsb[nm] = 0; stepsa[nm] = 0
        for prec in arm["precisions"], N in arm["horizons"], mode in get(arm, "coordinate_mode", [prim["coordinate_mode"]]),
            dir in get(arm, "direction", [prim["direction"]]), tau in get(arm, "tau", [prim["tau"]]), m in arm["methods"],
            eps in (m == "L2" ? get(arm, "eps", [prim["eps"]]) : (nothing,)), id in subset_ids(arm["orbits"])
            before[nm] += 1; stepsb[nm] += N; prec_b[prec] = get(prec_b, prec, 0) + N
            key = (id, m, prec, N, tau, eps, dir, mode)
            if !haskey(seen, key)
                seen[key] = nm; after[nm] += 1; stepsa[nm] += N; prec_a[prec] = get(prec_a, prec, 0) + N
            end
        end
    end
    emit("O14", "independent unique item count (TOML block x panel flags, merge by run-type key) / plan item count", "$(length(seen)) / $(PLAN_REAL.counts.items)", "equal", pf(length(seen) == PLAN_REAL.counts.items))
    arm_ok = all(a -> before[a] == PLAN_REAL.steps.arms[a].items_before && stepsb[a] == PLAN_REAL.steps.arms[a].steps_before && stepsa[a] == PLAN_REAL.steps.arms[a].steps_after, keys(before))
    emit("D53b[full]", "real panel: per-arm items and steps (before and after merge) equal the independent count, 11 arms", arm_ok, "true", pf(arm_ok))
    prec_ok = all(p -> prec_b[p] == PLAN_REAL.steps.precisions[p].steps_before && prec_a[p] == PLAN_REAL.steps.precisions[p].steps_after, keys(prec_b))
    emit("D53b[full]", "real panel: totals per precision (before / after) Float64, Double64, BigFloat", join((@sprintf("%s %.5g / %.5g", p, Float64(prec_b[p]), Float64(prec_a[p])) for p in ("Float64", "Double64", "BigFloat")), "; "), "equal to independent", pf(prec_ok))
    emit("D53b[full]", "real panel: merge saving (orbit-steps) / chunks", "$(PLAN_REAL.steps.saving) / $(PLAN_REAL.counts.chunks)", "reported", "INFO")
    # items: canonical text rebuilt from fields and hashed by the platform tool
    rng = MersenneTwister(5301)
    idx = sort!(randperm(rng, length(PLAN_R.items))[1:(QUICK ? 120 : 400)])
    fields = [:orbit, :system, :param, :x0, :p0, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau, :eps, :direction, :coordinate_mode, :vectors, :grid_sha256, :code_version, :hash_version]
    function mytext(it)
        io = IOBuffer()
        for f in fields
            v = getfield(it, f)
            f == :eps && it.method != "L2" && continue
            print(io, f, "=", v, "\n")
        end
        String(take!(io))
    end
    dir = mkpath(joinpath(SCR, "o14"))
    files = String[]
    for (j, i) in enumerate(idx)
        p = joinpath(dir, "item$j.txt"); write(p, mytext(PLAN_R.items[i])); push!(files, p)
    end
    sums = String[]
    cd(dir) do
        for chunk in Iterators.partition(basename.(files), 50)
            out = read(`sha256sum $(collect(chunk))`, String)
            append!(sums, [split(l)[1] for l in split(strip(out), '\n')])
        end
    end
    n_text = count(j -> mytext(PLAN_R.items[idx[j]]) == EX.item_text(PLAN_R.items[idx[j]]), eachindex(idx))
    n_hash = count(j -> sums[j] == PLAN_R.items[idx[j]].hash, eachindex(idx))
    emit("O14", "independent item text (own field order) equals item_text, on $(length(idx)) seeded items of the reduced plan", n_text, "$(length(idx))", pf(n_text == length(idx)))
    emit("O14", "platform sha256sum of the independent text equals the plan hash", n_hash, "$(length(idx))", pf(n_hash == length(idx)))
    # the full plan's chunk count by an independent rule: chunks = (run type x cell), with cell = k_index:class for non-controls
    ids_cell(o) = o["control"] ? "controls" : "k$(o["k_index"]):$(o["class"])"
    cellof = Dict(Int(o["id"]) => ids_cell(o) for o in panel)
    for i in 1:nh4; cellof[nP + i] = "h4"; end
    # run type key = everything of an item except orbit
    rts = Dict{Tuple, Set{String}}()
    for (key, nm) in seen
        id, m, prec, N, tau, eps, dir, mode = key
        push!(get!(rts, (m, prec, N, tau, eps, dir, mode), Set{String}()), cellof[id])
    end
    nch = sum(length(v) for v in values(rts))
    emit("O14", "independent chunk count (run type x cell) / plan chunk count", "$nch / $(PLAN_REAL.counts.chunks)", "equal", pf(nch == PLAN_REAL.counts.chunks))

    # synthetic panel with the preregistration's caps (1100, 200, 100): copy of the real panel with orbit records appended
    tp = TOML.parsefile(PANEL)
    orb = tp["orbit"]
    extra = [o for o in orb if !o["control"]]
    add = 1100 - length(orb)
    for j in 1:add
        o = deepcopy(extra[mod1(j, length(extra))]); o["id"] = length(orb) + j; o["in_SD"] = false; o["in_SL"] = false
        push!(orb, o)
    end
    for (cnt, key) in ((200 - nSD, "in_SD"), (100 - nSL, "in_SL"))
        for o in orb[(end - add + 1):(end - add + cnt)]
            o[key] = true
        end
    end
    sp = joinpath(SCR, "panel-synthetic-caps.toml")
    open(io -> TOML.print(io, tp), sp, "w")
    sp_plan = EX.plan_run(PREREG, sp)
    table = Dict("main" => 8.8e10, "main_l2" => 1.1e10, "long" => 1.6e11, "long_l2" => 2.0e10, "tau_sweep" => 2.0e10, "wrapped" => 1.8e10, "directions" => 3.0e10, "eps_sweep" => 8.0e9, "h4_float64" => 4.0e9,
                 "ladder_double64_short" => 4.0e9, "ladder_double64_long" => 1.0e10, "h4_double64" => 2.0e8, "ladder_bigfloat" => 3.0e8)
    worst = maximum(abs(sp_plan.steps.arms[a].steps_before - v) / v for (a, v) in table)
    emit("D53b[full]", "synthetic panel (1100, 200, 100): planned orbit-steps per arm vs the Budget table (13 rows: the A5 split of main and long into non-L2 and L2 arms), worst relative difference", worst, "<= 1e-9", pf(worst <= 1e-9))
    tot = Dict("Float64" => 3.59e11, "Double64" => 1.42e10, "BigFloat" => 3.0e8)
    wt = maximum(abs(sp_plan.steps.precisions[p].steps_before - v) / v for (p, v) in tot)
    emit("D53b[full]", "synthetic panel: precision totals (3.59e11, 1.42e10, 3.0e8)", wt, "<= 1e-9", pf(wt <= 1e-9))
    emit("D53b[full]", "synthetic panel: items, chunks, merge saving", "$(sp_plan.counts), saving $(sp_plan.steps.saving)", "-", "INFO")
    # plan determinism across processes is covered by the thread runs (plan_sha256 in every manifest)
    emit("D53b", "plan hash of the real plan (identical in every manifest below)", PLAN_REAL.hash, "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# Refusals on the real files (D53a, nothing written)
section("REFUSALS") do
    store = joinpath(SCR, "store-refusal"); rm(store; recursive = true, force = true)
    mf = joinpath(SCR, "refusal-manifest.json"); rm(mf; force = true)
    r1 = try; EX.run_campaign(PRE_R, PANEL; store = store, manifest = mf, threads = K19 + 1, allow_dirty = true); "no error"; catch e; e isa ArgumentError ? e.msg : sprint(showerror, e); end
    emit("D53a", "thread count differing from [run.k19].threads without override is refused, message names threads, nothing written", first(r1, 90),
         "ArgumentError, 0 files", pf(occursin("thread", r1) && !isdir(store) && !isfile(mf)))
    bad = joinpath(SCR, "prereg-tampered.md"); t = read(PRE_R); t[findfirst(==(UInt8('F')), t)] = UInt8('G'); write(bad, t)
    r2 = try; EX.run_campaign(bad, PANEL; store = store, manifest = mf, threads = NT, thread_override = true, allow_dirty = true); "no error"; catch e; e isa ArgumentError ? e.msg : sprint(showerror, e); end
    emit("D53a", "preregistration with a changed byte is refused (E18), nothing written", first(r2, 90), "ArgumentError, 0 files", pf(occursin("SHA-256", r2) && !isdir(store) && !isfile(mf)))
    pb = joinpath(SCR, "panel-wrong-hash.toml"); tp = lf(read(PANEL, String)); write(pb, replace(tp, "650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46" => "0"^64))
    r3 = try; EX.run_campaign(PREREG, pb; store = store, manifest = mf, threads = NT, thread_override = true, allow_dirty = true); "no error"; catch e; e isa ArgumentError ? e.msg : sprint(showerror, e); end
    emit("D53a", "panel built under a hash outside panel_prereg_sha256_accepted is refused, nothing written", first(r3, 90), "ArgumentError, 0 files", pf(occursin("not in panel_prereg_sha256_accepted", r3) && !isdir(store) && !isfile(mf)))
    sts = [EX.exit_status(s) for s in ("complete", "complete with failed chunks", "stopped: control bound", "stopped: preregistration changed", "paused arms")]
    emit("D53u", "exit status for complete, failed chunks, stopped (2 reasons), paused arms", join(sts, ","), "0,2,3,3,4", pf(sts == [0, 2, 3, 3, 4]))
end

# ---------------------------------------------------------------------------------------------------
# The exact-exponent rule of A4 (h1_check arm, control_bound = 1e-8 at N1 = 1e7) scaled to the reduced N1, on the real panel: the campaign stops at once (R5 item 1c)
section("CONTROL") do
    store = joinpath(SCR, "store-verify-C"); mf = joinpath(SCR, "verify-C-manifest.json"); rm(store; recursive = true, force = true); rm(mf; force = true)
    t = @elapsed r = campaign(PRE_S, store, mf)
    first_ctrl = first(c.ordinal for c in PLAN_R.chunks if c.cell == "controls" && c.run_type.method == "L1" && "main" in c.arms)
    emit("D53j/D53w", "tight control rule (bound $(fmt(CB_TRIP)) at N = $(HZ["main"]), A5 scope exact_start): campaign status, chunks run before the stop, first control chunk of the plan, seconds", "$(r.status); ran $(length(r.ran)); control chunk $first_ctrl; $(round(t; digits = 1)) s",
         "stopped: control bound, at the first control chunk", pf(startswith(r.status, "stopped: control bound") && length(r.ran) <= first_ctrl))
    for line in r.reasons
        info("D53j/D53w", "stop reason: " * first(line, 700))
    end
    ids = sort!(unique(parse.(Int, [m.captures[1] for m in eachmatch(r"orbit (\d+) method", join(r.reasons, " "))])))
    info("D53j/D53w", "controls named by the stop (first control chunk, method L1): $(ids). The real rule (bound 1e-8 at N1 = 1e7) is the gated set only (A5 control_scope = exact_start): the other controls are INFO")
    CAMP["control_stop"] = (status = r.status, ids = ids, ran = length(r.ran))
end

# ---------------------------------------------------------------------------------------------------
# The reduced campaign (D53w): interrupted twice (one kill after N chunks, one during a data write) and resumed, against an uninterrupted campaign
const STORE_A, STORE_B = joinpath(SCR, "store-verify-A"), joinpath(SCR, "store-verify-B")
const MAN_A, MAN_B = joinpath(SCR, "verify-A-manifest.json"), joinpath(SCR, "verify-B-manifest.json")
section("CAMPAIGN") do
    nchunks = length(PLAN_R.chunks)
    kill1 = max(5, nchunks ÷ 5)
    if !(REUSE && isfile(MAN_A) && isfile(MAN_B))
        rm(STORE_A; recursive = true, force = true); rm(STORE_B; recursive = true, force = true); rm(MAN_A; force = true); rm(MAN_B; force = true)
        # --- A: session 1 killed after kill1 chunks
        t1 = @elapsed k1 = killed(() -> campaign(PRE_R, STORE_A, MAN_A; hooks = EX.RunHooks(mutate = inject(), kill = kill_after(kill1))))
        done1 = count(f -> occursin(r"^chunk-\d{5}\.h5$", f), readdir(STORE_A))
        emit("D53w", "session 1 of the interrupted campaign killed after $kill1 chunks (injected), final files / temporary files left", "$k1: $done1 / $(count(endswith(".tmp"), readdir(STORE_A)))", "killed", pf(k1))
        # --- session 2 killed in the middle of the data of the next chunk
        missing_ord = first(o for o in 1:nchunks if !isfile(joinpath(STORE_A, chunkname(o))))
        t2 = @elapsed k2 = killed(() -> campaign(PRE_R, STORE_A, MAN_A; hooks = EX.RunHooks(mutate = inject(), kill = kill_during(missing_ord))))
        tmps = filter(endswith(".tmp"), readdir(STORE_A))
        emit("D53w", "session 2 killed during the data of chunk $missing_ord: no final file for it, temporary file left (I13)", "$k2: final=$(isfile(joinpath(STORE_A, chunkname(missing_ord)))) tmp=$(tmps)", "killed, no final file, one .tmp",
             pf(k2 && !isfile(joinpath(STORE_A, chunkname(missing_ord))) && length(tmps) == 1))
        CAMP["after_kill"] = run_reader(STORE_A, "A-after-kill")
        rk = parse_reader(CAMP["after_kill"].path)
        emit("O13", "independent reader on the store after the two kills: every final file complete and valid; non-chunk entries ignored (.tmp)", CAMP["after_kill"].ok, "reader PASS", pf(CAMP["after_kill"].ok))
        # --- session 3: resume
        t3 = @elapsed rA = campaign(PRE_R, STORE_A, MAN_A; hooks = EX.RunHooks(mutate = inject()))
        cp(MAN_A, MAN_A * ".session3"; force = true)
        CAMP["A"] = (res = rA, t = (t1, t2, t3))
        emit("D53w", "interrupted campaign resumed to the end: status, chunks run in session 3, skipped, failed twice in session 3, wall s", "$(rA.status): ran $(length(rA.ran)), skipped $(length(rA.skipped)), failed $(length(rA.failed)), $(round(t3; digits = 1)) s",
             "complete or complete with failed chunks", pf(rA.status in ("complete", "complete with failed chunks")))
        # --- B: uninterrupted
        tB = @elapsed rB = campaign(PRE_R, STORE_B, MAN_B; hooks = EX.RunHooks(mutate = inject()))
        CAMP["B"] = (res = rB, t = tB)
        emit("D53w", "uninterrupted campaign (same injections): status, chunks run, failed twice, wall s", "$(rB.status): ran $(length(rB.ran)), failed $(length(rB.failed)), $(round(tB; digits = 1)) s", "complete or complete with failed chunks", pf(rB.status in ("complete", "complete with failed chunks") && ORD_TWICE in rB.failed))
    else
        CAMP["A"] = (res = (status = EX.read_json(MAN_A)["status"],), t = ())
        CAMP["B"] = (res = (status = EX.read_json(MAN_B)["status"],), t = NaN)
        info("D53w", "REUSE: stores A and B taken from the scratch folder")
    end
    # injected failure behaviour from the files
    for (ord, expect_status, expect_attempt) in ((ORD_ONCE, "ok", 2), (ORD_TWICE, "failed", 2))
        ch = EX.read_chunk(joinpath(STORE_B, chunkname(ord)))
        emit("D53k", "injected NaN, chunk $ord: file status / attempt / backend of the attempt / first-attempt reasons recorded", "$(ch.attrs["status"]) / $(ch.attrs["attempt"]) / $(ch.attrs["backend"]) / $(length(ch.attrs["first_attempt_reasons"])) reason(s)",
             "$expect_status / $expect_attempt / cpu", pf(ch.attrs["status"] == expect_status && ch.attrs["attempt"] == expect_attempt && ch.attrs["backend"] == "cpu" &&
                                                         !isempty(ch.attrs["first_attempt_reasons"]) && (expect_status == "ok" || !isempty(ch.attrs["failure_reasons"]))))
    end
    # nothing to do on a finished campaign (D53o a)
    before = tree_listing(STORE_A)
    r2 = campaign(PRE_R, STORE_A, MAN_A)
    emit("D53o", "second run on the finished campaign: chunks run, files modified", "$(length(r2.ran)), $(before == tree_listing(STORE_A) ? 0 : 1)", "0, 0", pf(isempty(r2.ran) && before == tree_listing(STORE_A)))
    # manifests
    for (tag, mp, store) in (("A", MAN_A, STORE_A), ("B", MAN_B, STORE_B))
        m, v = validate_manifest(mp, PLAN_R; status_ok = "complete with failed chunks")
        CAMP["manifest-$tag"] = m
        ok = isempty(v.miss) && isempty(v.missing_chunk_keys) && isempty(v.missing_arm_keys) && v.sum_ok && v.order_ok && v.prereg_stable && v.panel_stable && v.plan_hash && v.schema && v.git40 && v.dirty_recorded
        emit("D53p/R7", "manifest $tag: all R7 keys, chunk and arm entries, statuses sum to the plan ($(v.n_chunks) of $(v.plan_chunks): $(v.statuses)), hashes stable, plan hash, git sha, dirty override recorded, sessions $(v.sessions)", ok, "true", pf(ok))
        isempty(v.miss) || info("D53p/R7", "manifest $tag missing keys: $(v.miss), chunk keys $(v.missing_chunk_keys), arm keys $(v.missing_arm_keys)")
    end
    m3 = EX.read_json(MAN_A * ".session3")
    emit("D53p/R7", "manifest A after the resume session: 3 sessions recorded, the killed chunk's temporary file listed in removed_temporaries and gone from the store", "$(length(m3["sessions"])) sessions; removed_temporaries $(m3["removed_temporaries"]); .tmp left: $(count(endswith(".tmp"), readdir(STORE_A)))", "3; one .tmp listed; 0 left",
         pf(length(m3["sessions"]) == 3 && length(m3["removed_temporaries"]) == 1 && isempty(filter(endswith(".tmp"), readdir(STORE_A)))))
    # failed chunks: the injected one and the natural ones, identical in A (resumed) and B
    fa = sort([c["ordinal"] for c in CAMP["manifest-A"]["chunks"] if c["status"] == "failed"]); fb = sort([c["ordinal"] for c in CAMP["manifest-B"]["chunks"] if c["status"] == "failed"])
    nat = [c for c in CAMP["manifest-B"]["chunks"] if c["status"] == "failed" && c["ordinal"] != ORD_TWICE]
    emit("D53w", "chunks failed twice: injected chunk $ORD_TWICE present; failed set of the resumed campaign equals that of the uninterrupted one; natural failures (not injected)", "$(length(fb)) failed, $(length(nat)) natural, sets equal: $(fa == fb)", "ORD_TWICE in set, sets equal", pf(ORD_TWICE in fb && fa == fb))
    emit("D53k", "failed twice (natural and injected) against the REAL stop_failed_fraction 0.05 x chunks (the campaign stops when strictly more): failed twice / limit", "$(length(fb)) / $(PLAN_REAL.block.raw["health"]["stop_failed_fraction"] * length(PLAN_R.chunks))", "<= limit (no stop)", pf(length(fb) <= PLAN_REAL.block.raw["health"]["stop_failed_fraction"] * length(PLAN_R.chunks)))
    cls = Dict{String, Int}()
    for c in nat
        mm = match(r"method=(\S+) .*horizon=(\d+) tau=(\d+) eps=(\S+) direction=(\S+) mode=(\S+) cell=(\S+)", c["key"])
        reason = replace(string(get(c, "first_attempt_reasons", [""])[1]), r"-?\d+(\.\d+)?(e-?\d+)?" => "#")
        reason = first(reason, 70)
        k = "$(mm.captures[1]) tau=$(mm.captures[3]) mode=$(mm.captures[6]) | $(first(reason, 70))"
        cls[k] = get(cls, k, 0) + 1
    end
    for (k, n) in sort(collect(cls); by = first)
        say("NATURAL-FAILURE | $n chunk(s) | $k")
    end
    nat_cells = sort(unique([match(r"cell=(\S+)", c["key"]).captures[1] for c in nat]))
    info("D53k", "natural failures by cell: $nat_cells; methods: $(sort(unique([match(r"method=(\S+)", c["key"]).captures[1] for c in nat]))); tau: $(sort(unique([match(r"tau=(\d+)", c["key"]).captures[1] for c in nat])))")
    mb = CAMP["manifest-B"]
    emit("D53p/R7", "manifest B: device, backend table rows, git clean flag, thread count recorded", "device=$(mb["device"]); table rows $(length(mb["backend_table"])); clean=$(mb["git"]["clean"]); threads=$(mb["threads"])", "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# O13: independent reader on the finished stores, and I14 (resume equals uninterrupted)
section("READER") do
    for tag in ("A", "B")
        store = tag == "A" ? STORE_A : STORE_B
        rr = run_reader(store, tag)
        R = parse_reader(rr.path)
        CAMP["R$tag"] = R
        nfiles = count(f -> occursin(r"^chunk-\d{5}\.h5$", f), readdir(store))
        emit("O13", "independent reader (HDF5.jl in a separate process, no Chaoslab), store $tag: files read, checks failed, seconds", "$(length(R.files)) of $nfiles files, $(count(l -> occursin("|FAIL|", l), R.checks)) failed, $(round(rr.seconds; digits = 1)) s", "all files, 0 failed", pf(rr.ok && length(R.files) == nfiles && !any(l -> occursin("|FAIL|", l), R.checks)))
        foreach(l -> say("   reader: " * l), first([l for l in R.checks if occursin("|FAIL|", l)], 5))
        cmp = compare_readers(store, R)
        emit("O13", "store $tag: module read_chunk equals the independent reader: files, attributes, datasets (sha-256 of raw bytes, dtype, shape)", "$(cmp.files) files, $(cmp.attrs) attributes, $(cmp.datasets) datasets, $(cmp.bad) differ", "0 differ", pf(cmp.bad == 0 && cmp.files == nfiles))
        foreach(x -> say("   differ: " * x), cmp.examples)
        # layout dims in file order
        dims = Dict{String, Set{String}}()
        for ((f, o, p), d) in R.data
            push!(get!(dims, replace(p, r"/.*" => ""), Set{String}()), d.eltype)
        end
        info("D53m", "store $tag datasets by group and element type: " * join(("$g: $(join(sort!(collect(t)), "/"))" for (g, t) in sort(collect(dims); by = first)), "; "))
    end
    # file-order (C) dimensions of one lyap and state dataset
    R = CAMP["RB"]
    ex_l = first(k for (k, d) in R.data if startswith(k[3], "lyap/") && split(d.jdims, 'x')[2] == "2")
    d = R.data[ex_l]; s = R.data[(ex_l[1], ex_l[2], "state")]
    info("D53m", "dimension convention: Julia shape of $(ex_l[3]) is $(d.jdims), HDF5 file extents (row-major, as h5py/h5dump see them) are $(d.cdims); state $(s.jdims) vs file $(s.cdims). Schema text says n x k and n x 2")
    emit("D53m", "file extents are the Julia shape reversed (HDF5.jl convention): a non-Julia reader sees (k, n) and (2, n), not the n x k of schema/run-store.md", d.cdims == join(reverse(split(d.jdims, 'x')), "x"), "documented, schema ambiguity", d.cdims == join(reverse(split(d.jdims, 'x')), "x") ? "INFO" : "FAIL")
    # I14
    df = diff_readers(CAMP["RA"], CAMP["RB"])
    emit("I14/D53o", "resumed (killed twice) campaign A against uninterrupted campaign B: differing datasets/attributes outside times and host fields (created_utc, cpu, threads, session), via independent-reader digests",
         df.total, "0", pf(df.total == 0))
    foreach(x -> say("   diff: " * x), df.examples)
    isempty(df.diffs) || info("I14/D53o", "difference classes: $(df.diffs)")
    # h5dump spot check: the HDF5 library's own command-line tool (HDF5_jll), outside HDF5.jl and outside the module
    let
        rng = MersenneTwister(5302)
        files = sort(filter(f -> occursin(r"^chunk-\d{5}\.h5$", f), readdir(STORE_B)))
        okfiles = [f for f in files if EX.read_chunk(joinpath(STORE_B, f)).attrs["status"] == "ok"]
        sample = okfiles[randperm(rng, length(okfiles))[1:min(12, end)]]
        nbad = ndone = 0
        HDF5_jll.h5dump() do exe
            for fn in sample
                ch = EX.read_chunk(joinpath(STORE_B, fn))
                o = ch.orbits[min(2, end)]
                m = first(keys(o.lyap))
                for (path, A) in (("/orbits/$(o.id)/lyap/$m", o.lyap[m]), ("/orbits/$(o.id)/state", o.state), ("/orbits/$(o.id)/t", o.t))
                    txt = read(`$exe -d $path -m %.17g -y -w 0 $(joinpath(STORE_B, fn))`, String)
                    body = match(r"DATA \{([^}]*)\}", txt)
                    body === nothing && (nbad += 1; continue)
                    toks = filter(!isempty, split(String(body.captures[1]), r"[,\s]+"))
                    vals = eltype(A) === Int64 ? parse.(Int64, toks) : parse.(Float64, toks)
                    ndone += 1
                    (length(vals) == length(A) && all(i -> isequal(vals[i], vec(A)[i]), eachindex(vals))) || (nbad += 1)
                end
            end
        end
        emit("O13", "HDF5 command-line tool h5dump (-m %.17g; no HDF5.jl, no Chaoslab) on $(length(sample)) ok files, lyap / state / t datasets of one orbit each: equal to the module's values exactly", "$ndone datasets, $nbad differ", "0 differ", pf(nbad == 0 && ndone > 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# D53t on the CPU, O12: rerun on 20 items per arm
section("RERUN") do
    for (tag, store) in (("B", STORE_B),)
        t = @elapsed r = EX.rerun_check(PRE_R, PANEL, store; threads = NT, thread_override = NT != K19, allow_dirty = true)
        nb = count(x -> x.bitwise, r)
        arms = unique(x.arm for x in r)
        emit("D53t/O12", "CPU reference rerun of rerun_sample = $(PLAN_R.block.raw["verify"]["rerun_sample"]) items per arm, store $tag ($(length(arms)) arms, $(round(t; digits = 1)) s): items, bitwise equal, backend of stored items", "$(length(r)), $nb, $(join(unique(x.backend for x in r), ","))", "all bitwise (CPU)", pf(!isempty(r) && nb == length(r) && length(arms) == length(PLAN_R.block.raw["arm"])))
        failed_inj = count(x -> !x.bitwise, r)
        failed_inj > 0 && info("D53t/O12", "items not bitwise: " * join(("$(x.item) ($(x.arm)) max diff $(x.max_abs_diff)" for x in r if !x.bitwise), "; "))
    end
end

# ---------------------------------------------------------------------------------------------------
# D53i: process-level thread runs through the command-line entry point (-t 1, 4, 23)
const THREAD_TAGS = QUICK ? [1, 4, NT] : [1, 4, 23]
section("THREADS") do
    for n in THREAD_TAGS
        store = joinpath(SCR, "store-verify-t$n"); mf = joinpath(SCR, "verify-t$n-manifest.json"); logf = joinpath(SCR, "cli-t$n.log")
        if !(REUSE && isfile(mf))
            rm(store; recursive = true, force = true); rm(mf; force = true)
            t = @elapsed pr = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t $n $(joinpath(ROOT, "bin", "chaoslab.jl")) run $PRE_R $PANEL --store $store --manifest $mf --allow-dirty --thread-override --capacity $CAPACITY`);
                                           stdout = logf, stderr = logf))
            CAMP["cli$n"] = (code = pr.exitcode, seconds = t)
            emit("D53i", "julia -t $n bin/chaoslab.jl run (reduced block): exit code, wall seconds", "$(pr.exitcode), $(round(t; digits = 1)) s", "0 (complete) or 2 (complete with failed chunks)", pf(pr.exitcode in (0, 2)))
        end
        rr = run_reader(store, "t$n")
        CAMP["Rt$n"] = parse_reader(rr.path)
        emit("O13", "independent reader on store of -t $n: ok, files", "$(rr.ok), $(length(CAMP["Rt$n"].files))", "PASS, $(length(PLAN_R.chunks))", pf(rr.ok && length(CAMP["Rt$n"].files) == length(PLAN_R.chunks)))
    end
    ref = THREAD_TAGS[1]
    for n in THREAD_TAGS[2:end]
        df = diff_readers(CAMP["Rt$ref"], CAMP["Rt$n"])
        emit("D53i", "per-dataset sha-256 (independent reader): -t $ref against -t $n, differing datasets or attributes (outside created_utc, cpu, threads, session)", df.total, "0", pf(df.total == 0))
        foreach(x -> say("   diff: " * x), df.examples)
    end
    # a combined digest per store over all dataset hashes (and cost entries 1, 2, 5), for the report
    for n in THREAD_TAGS
        R = CAMP["Rt$n"]
        h = bytes2hex(sha256(join(sort!(["$(k[1])|$(k[2])|$(k[3])|$(d.sha)" for (k, d) in R.data if !startswith(k[3], "cost/")]), "\n")))
        hc = bytes2hex(sha256(join(sort!(["$(k[1])|$(k[2])|$(k[3])|$v" for (k, v) in R.ksha]), "\n")))
        info("D53i", "-t $n: combined SHA-256 of all $(count(kv -> !startswith(kv[1][3], "cost/"), R.data)) datasets $h; cost entries (steps, vectors, bytes) $hc")
    end
    # manifests
    mans = Dict(n => EX.read_json(joinpath(SCR, "verify-t$n-manifest.json")) for n in THREAD_TAGS)
    function norm(m)
        d = Dict{String, Any}(k => v for (k, v) in m if !(k in ("sessions", "threads", "wall", "launch", "overrides", "versions", "host", "removed_temporaries")))
        d["chunks"] = [Dict{String, Any}(k => v for (k, v) in c if !(k in ("wall_seconds", "thread_seconds"))) for c in m["chunks"]]
        d["backend_table"] = [Dict{String, Any}(k => v for (k, v) in r if k != "note") for r in m["backend_table"]]      # the note records the thread-count mismatch with the K19 record
        d["arms"] = Dict{String, Any}(a => Dict{String, Any}(k => v for (k, v) in x if !(k in ("wall_seconds", "thread_seconds", "cost_ns_per_orbit_step"))) for (a, x) in m["arms"])
        return d
    end
    for n in THREAD_TAGS[2:end]
        a, b = norm(mans[ref]), norm(mans[n])
        diffkeys = [k for k in union(keys(a), keys(b)) if get(a, k, nothing) != get(b, k, nothing)]
        emit("D53i", "manifests -t $ref and -t $n differ only in provenance (sessions, threads, wall, launch inputs, overrides, versions, host, times, costs): other differing keys", join(diffkeys, ","), "none", pf(isempty(diffkeys)))
    end
    for n in THREAD_TAGS
        m = mans[n]
        info("D53i", "-t $n manifest: threads $(m["threads"]), override $(get(m["overrides"], "threads", "none")), status $(m["status"]), plan_sha256 $(first(m["plan_sha256"], 12)), wall $(round(m["wall"]["total_seconds"]; digits = 1)) s")
    end
    emit("D53i", "plan hash identical in the manifests of all processes and the in-process plan", length(unique(vcat([m["plan_sha256"] for m in values(mans)], PLAN_R.hash))), "1", pf(length(unique(vcat([m["plan_sha256"] for m in values(mans)], PLAN_R.hash))) == 1))
    dfb = diff_readers(CAMP["Rt$(THREAD_TAGS[end])"], get(CAMP, "RB", CAMP["Rt$(THREAD_TAGS[end])"]))
    info("D53i", "store of -t $(THREAD_TAGS[end]) process vs in-process B: differing datasets $(dfb.total) (B carries the injected failures: expected differences are chunk $ORD_ONCE (first-attempt reasons) and $ORD_TWICE (failed)): classes $(dfb.diffs)")
end

# ---------------------------------------------------------------------------------------------------
# D53f: direct call (O9) at Float64, Double64 and BigFloat for every method that supports them
function ex_with_prec(f, label)
    T, bits_ = Chaoslab.parse_precision(label)
    return bits_ === nothing ? f(T) : setprecision(() -> f(T), BigFloat, bits_)
end
const RED = Dict("L1" => () -> Chaoslab.Benettin(), "L3" => () -> Chaoslab.QRSpectrum(Chaoslab.CGS()), "L4" => () -> Chaoslab.QRSpectrum(Chaoslab.MGS()),
                 "L4b" => () -> Chaoslab.QRSpectrum(Chaoslab.RGS()), "L5" => () -> Chaoslab.QRSpectrum(Chaoslab.HouseholderQR()),
                 "L5b" => () -> Chaoslab.QRSpectrum(Chaoslab.AugmentedQR()), "L6" => () -> Chaoslab.QRSpectrum(Chaoslab.LapackQR()), "L9" => () -> Chaoslab.SVDSpectrum())
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
function direct(item, n)
    return ex_with_prec(item.precision) do T
        sys, x = sys_x(item, T)
        W = EX.deviation_frame(T, item.direction, item.vectors)
        if item.method == "L2"
            u0 = SVector{2, T}(W[1, 1], W[2, 1])
            _, _, fin = Chaoslab.propagate_shadow_orbit(sys, x, u0, n, item.tau, Chaoslab.exact(T, item.eps))
            return Float64[fin.X2]
        end
        res = Chaoslab.propagate(sys, [x], n; k = item.vectors, tau = item.tau, seed = parse(Int, item.run_seed), reducers = (RED[item.method](),), W0 = [W])
        return Float64.(collect(res.finals[1][1].lambda))
    end
end
section("EXTENDED") do
    hz = QUICK ? (f64 = 320, d64 = 320, bf = 160) : (f64 = 960, d64 = 480, bf = 160)
    ext = make_prereg(joinpath(SCR, "prereg-ext.md"), run -> begin
        m8 = ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L9"]
        run["arm"] = [Dict{String, Any}("name" => "f64", "orbits" => "S_L", "methods" => [m8; "L6"], "precisions" => ["Float64"], "horizons" => [hz.f64]),
                      Dict{String, Any}("name" => "f64w", "orbits" => "S_L", "methods" => [m8; "L6"], "precisions" => ["Float64"], "horizons" => [hz.f64], "tau" => [4], "coordinate_mode" => ["wrapped"]),
                      Dict{String, Any}("name" => "d64", "orbits" => "S_L", "methods" => m8, "precisions" => ["Double64"], "horizons" => [hz.d64]),
                      Dict{String, Any}("name" => "d64w", "orbits" => "S_L", "methods" => m8, "precisions" => ["Double64"], "horizons" => [hz.d64], "tau" => [4], "coordinate_mode" => ["wrapped"]),
                      Dict{String, Any}("name" => "bf", "orbits" => "S_L", "methods" => m8, "precisions" => ["BigFloat"], "horizons" => [hz.bf]),
                      Dict{String, Any}("name" => "h4d", "orbits" => "h4", "methods" => m8, "precisions" => ["Double64"], "horizons" => [hz.d64])]
        run["health"]["first_chunk_cost_factor"] = 1e6
        run["health"]["control_bound"] = 1.0
        run["pilot"]["horizon"] = hz.bf
    end)
    store = joinpath(SCR, "store-verify-ext"); mf = joinpath(SCR, "verify-ext-manifest.json"); rm(store; recursive = true, force = true)
    t = @elapsed r = campaign(ext, store, mf)
    emit("D53f", "extended-precision campaign (S_L, h4; Float64 x 9 methods, Double64 and BigFloat x 8 methods; tau 1 lifted and tau 4 wrapped): status, chunks, seconds", "$(r.status), $(length(r.ran)) chunks, $(round(t; digits = 1)) s", "complete", pf(r.status == "complete"))
    plan = r.plan
    rng = MersenneTwister(5306)
    cases = Dict{Tuple, Vector{Int}}()
    for c in plan.chunks, i in c.items
        it = plan.items[i]
        push!(get!(cases, (it.method, it.precision, it.coordinate_mode), Int[]), i)
    end
    nbad = ntot = 0; per = Dict{Tuple, Tuple{Int, Int}}()
    for ((m, prec, mode), ids) in sort(collect(cases); by = first)
        pick = ids[randperm(rng, length(ids))[1:min(3, end)]]
        nb = nt = 0
        for i in pick
            it = plan.items[i]
            ch = EX.read_chunk(joinpath(store, chunkname(first(c.ordinal for c in plan.chunks if i in c.items))))
            o = only(x for x in ch.orbits if x.id == it.orbit)
            L = o.lyap[m]
            for j in sort!(randperm(rng, length(o.t))[1:min(5, end)])
                d = direct(it, o.t[j])
                nt += 1
                (all(isequal.(d, L[j, 1:length(d)]))) || (nb += 1)
            end
        end
        per[(m, prec, mode)] = (nb, nt); nbad += nb; ntot += nt
    end
    emit("D53f", "stored lambda equals the direct call of the verified estimator, bit for bit: 3 orbits x 5 grid steps per (method, precision, mode)", "$nbad differences in $ntot comparisons over $(length(per)) cells", "0", pf(nbad == 0 && ntot > 0))
    for ((m, prec, mode), (nb, nt)) in sort(collect(per); by = first)
        nb > 0 && info("D53f", "FAIL cell $m $prec $mode: $nb of $nt differ")
    end
    info("D53f", "cells: " * join(sort(unique(("$m/" * replace(p, r"\(.*" => "") * "/$(md[1])") for (m, p, md) in keys(per))), " "))
    # unsupported: L6 at Double64
    r6 = try
        make_prereg(joinpath(SCR, "prereg-l6d.md"), run -> (run["arm"] = [Dict{String, Any}("name" => "x", "orbits" => "S_L", "methods" => ["L6"], "precisions" => ["Double64"], "horizons" => [200])];
                                                          run["pilot"]["horizon"] = 200))
        EX.plan_run(joinpath(SCR, "prereg-l6d.md"), PANEL); "plan accepted"
    catch e
        e isa ArgumentError ? first(e.msg, 120) : sprint(showerror, e)
    end
    info("D53f", "L6 (library QR) at Double64 in a block: $r6")
    # ext store readable by the independent reader
    rr = run_reader(store, "ext")
    emit("O13", "independent reader on the extended-precision store (Double64, BigFloat, wrapped, tau 4)", rr.ok, "PASS", pf(rr.ok))
    CAMP["ext"] = (store = store, plan = plan)
end

# ---------------------------------------------------------------------------------------------------
# D53s: Y5 and exact exponents
const EPS52 = ldexp(1.0, -52)
function c_factor(A, v)
    return setprecision(BigFloat, 256) do
        M = BigFloat.(Matrix(A)); t = M[1, 1] + M[2, 2]
        mu = (t + (t >= 0 ? 1 : -1) * sqrt(t^2 - 4)) / 2         # the eigenvalue of larger modulus (negative at x = pi, K > 4)
        e = [M[1, 2], mu - M[1, 1]]; l = [M[2, 1], mu - M[1, 1]]
        vb = BigFloat.(collect(v))
        Float64(abs(l[1] * vb[1] + l[2] * vb[2]) * sqrt(e[1]^2 + e[2]^2) / abs(l[1] * e[1] + l[2] * e[2]))
    end
end
lnmu_of(A) = setprecision(() -> (t = BigFloat(A[1, 1]) + BigFloat(A[2, 2]); Float64(log(abs(t + (t >= 0 ? 1 : -1) * sqrt(t^2 - 4)) / 2))), BigFloat, 256)
BS(n, K, G) = 8 * EPS52 * (n + G) + 16 * abs(K) * n^2 * EPS52
Gfirst(K, x) = (c = K * cos(x); (1 + c)^2 + 1 + c^2 + 1)
section("Y5") do
    # Y5 control at the REAL horizon 1e7 (arm main, CPU, all eight variational methods), direction pi/4 (unstable eigenvector of the symmetric block: c = 1)
    plan = PLAN_REAL
    y5 = Int(last(TOML.parsefile(PANEL)["orbit"])["id"])
    pan = TOML.parsefile(PANEL)["orbit"]; yo = last(pan)
    mu = setprecision(() -> Float64((3 + sqrt(BigFloat(5))) / 2), BigFloat, 256)
    A = Chaoslab.y5_control(mu).sys.A
    G = mu^2 + mu^-2 + 2
    lnmu = lnmu_of(A)
    emit("D53s", "Y5 control: panel exponent equals ln mu at 256 bits (mu = (3 + sqrt 5) / 2)", yo["exponent"], @sprintf("%.17g", lnmu), pf(abs(yo["exponent"] - lnmu) <= 1e-15))
    items = [it for it in plan.items if it.orbit == y5 && it.precision == "Float64" && it.tau == 1 && it.coordinate_mode == "lifted" && it.horizon == 10^7 && it.method != "L2"]
    ctx = EX.make_ctx(plan, items)
    t = @elapsed recs = EX.compute_records(ctx, items, NT)
    emit("D53s", "Y5 at horizon 1e7: items run on the CPU ($NT threads): count, seconds", "$(length(items)), $(round(t; digits = 1)) s", "-", "INFO")
    worst = Dict{String, Float64}(); worst_dir = Dict{String, Float64}(); nstep = 0
    for (it, r) in zip(items, recs)
        L = r.lyap[it.method][:, 1]
        v = Float64.(EX.deviation_frame(Float64, it.direction, 2)[:, 1])
        cfac = c_factor(A, v)
        for (j, n) in enumerate(r.t)
            bound = BS(n, 0.0, G) / n
            nstep += 1
            if it.direction == "pi/4"
                worst[it.method] = max(get(worst, it.method, 0.0), abs(L[j] - lnmu) / bound)
            elseif it.method != "L9" && n >= 40
                worst_dir[it.method * "/" * it.direction] = max(get(worst_dir, it.method * "/" * it.direction, 0.0), abs(L[j] - (lnmu + log(cfac) / n)) / (bound + 1e-12 / 1))
            end
        end
        it.direction == "pi/4" && abs(cfac - 1) >= 1e-15 && info("D53s", "c for pi/4 is $cfac, not 1")
    end
    w = maximum(values(worst))
    emit("D53s", "Y5, direction pi/4 (unstable eigenvector, c = 1), horizon 1e7, $(length(worst)) variational methods, all recorded grid steps: max |lambda_1(n) - ln mu| / (B_S(n)/n), B_S = 8 eps (n + G)", w, "<= 1", pf(w <= 1))
    info("D53s", "per method: " * join(("$m $(fmt(v))" for (m, v) in sort(collect(worst))), ", ") * "; absolute error at n = 1e7: " *
         join(("$(it.method) $(fmt(abs(r.lyap[it.method][end, 1] - lnmu)))" for (it, r) in zip(items, recs) if it.direction == "pi/4"), ", "))
    if !isempty(worst_dir)
        wd = maximum(values(worst_dir))
        emit("D53s", "Y5, other start directions in the plan (L1, L5, ...): max |lambda_1(n) - (ln mu + ln(c)/n)| / bound, c at 256 bits", wd, "<= 1", pf(wd <= 1))
    else
        info("D53s", "Y5 is not in the S_D subset of the real panel: no direction-sweep items for it (only the primary direction)")
    end
    # exact-exponent controls (h1_check rule) in the reduced campaign: every control with a positive exponent, every control method, last grid step
    ctrl = [o for o in pan if o["control"] && o["exponent"] > 0]
    chunk_of = Dict(i => c.ordinal for c in PLAN_R.chunks for i in c.items)
    Nm = HZ["main"]
    devs = Dict{Int, Dict{String, Float64}}()
    pred = Dict{Int, NamedTuple}()
    for o in ctrl
        if o["control_kind"] == "fixed" && o["system"] == "standard"
            Kv = Float64(o["K"]); x0 = Float64(o["x0"])
            Am = [1 + Kv * cos(x0) 1.0; Kv * cos(x0) 1.0]
            pred[o["id"]] = (lm = lnmu_of(Am), lnc = log(c_factor(Am, Float64.(EX.deviation_frame(Float64, "pi/4", 2)[:, 1]))), K = Kv, G = Gfirst(Kv, x0))
        end
    end
    onorb = Dict{Int, Bool}(); fworst = 0.0; nchk = 0
    for o in ctrl, it in PLAN_R.items
        it.orbit == o["id"] && it.precision == "Float64" && it.tau == 1 && it.direction == "pi/4" && it.coordinate_mode == "lifted" && it.horizon == Nm && "main" in it.arms &&
            it.method in PLAN_R.block.raw["health"]["control_methods"] || continue
        ch = EX.read_chunk(joinpath(STORE_B, chunkname(chunk_of[it.id])))
        r = only(x for x in ch.orbits if x.id == o["id"])
        devs[o["id"]] = merge(get(devs, o["id"], Dict{String, Float64}()), Dict(it.method => abs(r.lyap[it.method][end, 1] - o["exponent"])))
        if haskey(pred, o["id"]) && it.method != "L9"
            p = pred[o["id"]]; L = r.lyap[it.method][:, 1]
            ok = true
            for (j, n) in enumerate(r.t)
                n >= 20 || continue
                ratio = abs(L[j] - (p.lm + p.lnc / n)) / (BS(n, p.K, p.G) / n + 1e-12)
                nchk += 1; ratio > 1 && (ok = false)
                ok && (fworst = max(fworst, ratio))
            end
            onorb[o["id"]] = get(onorb, o["id"], true) && ok
        end
    end
    say("CONTROL-TABLE | id | kind | exact exponent | max over control methods |lambda_1(N) - exact| at N = $Nm | follows ln|mu+| + ln(c)/n (fixed points) | |ln c| / 1e7 = projected bias at N1 = 1e7 | verdict against control_bound 1e-8 at 1e7")
    offorbit = Int[]; biased = Int[]
    for o in sort(ctrl; by = o -> o["id"])
        id = o["id"]; d = maximum(values(devs[id]))
        if haskey(pred, id)
            on = onorb[id]; pb = abs(pred[id].lnc) / 1e7
            verdict = !on ? "OFF ORBIT in Float64" : pb > 1e-8 ? "on orbit, but bias exceeds 1e-8 at 1e7" : "ok"
            !on && push!(offorbit, id); on && pb > 1e-8 && push!(biased, id)
            say(@sprintf("CONTROL-TABLE | %d | %s | %.6f | %.3g | %s | %.3g | %s", id, o["control_kind"], o["exponent"], d, on, pb, verdict))
        elseif o["control_kind"] == "y5"
            say(@sprintf("CONTROL-TABLE | %d | y5 | %.6f | %.3g | direction pi/4 is the unstable eigenvector, c = 1 | 0 | ok", id, o["exponent"], d))
        else
            d > CB_SCALED && push!(offorbit, id)
            say(@sprintf("CONTROL-TABLE | %d | %s (periodic orbit) | %.6f | %.3g | not applicable | not computed | %s", id, o["control_kind"], o["exponent"], d, d > CB_SCALED ? "OFF ORBIT in Float64" : "within scaled bound"))
        end
    end
    emit("D53s", "fixed-point controls that stay on their orbit: max |lambda_1(n) - (ln|mu+| + ln(c)/n)| / (B_S/n + 1e-12) over $nchk grid points (n >= 20, controls that stay on the orbit)", fmt(fworst), "<= 1", pf(fworst <= 1))
    bad = sort(vcat(offorbit, biased))
    emit("D53w", "positive-exponent controls ($(length(ctrl))) and the A4 rule |lambda_1(N1) - exact| <= 1e-8 (h1_check arm, direction pi/4): controls that start off their orbit in Float64 / controls on the orbit whose bias |ln c| / N1 exceeds 1e-8 at N1 = 1e7",
         "off orbit: $(offorbit); bias above bound: $(biased)", "none", isempty(bad) ? "PASS" : "FAIL: controls $(join(bad, ", ")) would stop the real campaign by the control rule (R5 1c)")
    CAMP["control_findings"] = (offorbit = offorbit, biased = biased)
end

# ---------------------------------------------------------------------------------------------------
# Device: D53h real device, D53t on the device, D54f
section("DEVICE") do
    if F.DEVICE[] isa F.NoDevice
        info("D53h", "no functional device: device checks skipped"); return
    end
    dev = CUDA.name(CUDA.device())
    emit("D53h", "device in force, CUDA.jl, free memory MiB", "$dev, $(pkgversion(CUDA)), $(round(CUDA.free_memory() / 2^20))", "RTX 5060", pf(occursin("5060", dev)))
    rec = EX.gpu_record()
    emit("D54f", "record in force (data/gpu-benchmark.toml, Float64): threads, AC13 cpu / gpu ns, gpu_faster, recorder methods", "$(rec.threads), $(fmt(rec.cpu_ns_per_step)) / $(fmt(rec.gpu_ns_per_step)), $(rec.gpu_faster), $(sort(collect(keys(rec.recorder))))", "threads == $K19", pf(rec.threads == K19))
    emit("D54f", "recorder entries: all taken at $K19 threads on 2048 orbits x 2e5 steps; gpu_faster per method", join(("$m $(r.gpu_faster) (cpu $(fmt(r.cpu_ns_per_step)), gpu $(fmt(r.gpu_ns_per_step)), $(r.orbits)x$(r.steps))" for (m, r) in sort(collect(rec.recorder); by = first)), "; "), "threads == $K19, 2048 x 200000",
         pf(all(r -> r.threads == K19 && r.orbits == 2048 && r.steps == 200000, values(rec.recorder))))
    # resolved table for the real preregistration and panel
    tbl = EX.backend_table(PLAN_REAL; threads = K19, record = rec)
    rows = sort!(collect(tbl); by = first)
    for ((p, m), r) in rows
        say("TABLE | $p | $m | $(r.request) | $(r.backend) | $(r.device) | $(r.reason) | $(r.note)")
    end
    f64 = [(m, r) for ((p, m), r) in rows if p == "Float64"]
    kernel = [m for (m, r) in f64 if m != "L6"]
    expect_cpu = all(r -> string(r.backend) == "cpu", Iterators.map(last, f64)) && all(r -> r.reason == "GPU not faster for the element type", (r for (m, r) in f64 if m != "L6")) &&
                 only(r for (m, r) in f64 if m == "L6").reason == "reducer is CPU-only" && all(r -> string(r.backend) == "cpu" && isempty(r.reason), (r for ((p, m), r) in rows if p != "Float64"))
    emit("D54f/D54a", "real preregistration and panel at $K19 threads, auto: Float64 kernel methods resolve to cpu with the seventh reason (recorder gpu_faster false), L6 to cpu CPU-only, Double64 and BigFloat cpu", "$(length(kernel)) kernel methods", "all cpu", pf(expect_cpu))
    # per-method record consistency with D54a: gpu_faster per method decides
    consistent = all(m -> (string(tbl[("Float64", m)].backend) == "gpu") == (haskey(rec.recorder, m) && rec.recorder[m].threads == K19 && rec.recorder[m].gpu_faster), kernel)
    emit("D54a", "resolved backend equals gpu exactly for methods whose entry has threads == $K19 and gpu_faster", consistent, "true", pf(consistent))
    # forced gpu on a reduced arm, and the same arm on the CPU
    mk(name, backend, hzn) = make_prereg(joinpath(SCR, "prereg-dev-$name.md"), run -> begin
        run["arm"] = [Dict{String, Any}("name" => "dev_main", "orbits" => "P", "methods" => ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"], "precisions" => ["Float64"], "horizons" => [hzn], "h1_check" => true),
                      Dict{String, Any}("name" => "dev_wrapped", "orbits" => "S_D", "methods" => ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"], "precisions" => ["Float64"], "horizons" => [hzn], "coordinate_mode" => ["wrapped"])]
        run["backend"]["Float64"] = backend
        run["health"]["first_chunk_cost_factor"] = 1e6
        run["health"]["control_bound"] = CB_OPEN; run["health"]["stop_failed_fraction"] = 0.5
        run["pilot"]["horizon"] = hzn
    end)
    hzn = QUICK ? 4000 : 50_000
    pg, pc = mk("gpu", "gpu", hzn), mk("cpu", "cpu", hzn)
    sg, sc = joinpath(SCR, "store-verify-dev-gpu"), joinpath(SCR, "store-verify-dev-cpu")
    rm(sg; recursive = true, force = true); rm(sc; recursive = true, force = true)
    tg = @elapsed rg = EX.run_campaign(pg, PANEL; store = sg, manifest = joinpath(SCR, "verify-dev-gpu-manifest.json"), threads = NT, thread_override = NT != K19, allow_dirty = true, capacity = CAPACITY)
    tc = @elapsed rc = EX.run_campaign(pc, PANEL; store = sc, manifest = joinpath(SCR, "verify-dev-cpu-manifest.json"), threads = NT, thread_override = NT != K19, allow_dirty = true, capacity = CAPACITY)
    emit("D53h", "request gpu on the reduced arm (883 orbits x 9 methods lifted + 172 x 9 wrapped, N = $hzn): status, seconds; same arm request cpu: status, seconds", "$(rg.status) $(round(tg; digits = 1)) s; $(rc.status) $(round(tc; digits = 1)) s", "complete, complete", pf(rg.status == "complete" && rc.status == "complete"))
    # provenance in the chunk files
    gattrs = Dict{String, Dict{String, Any}}()
    prov = Dict{String, Set{String}}()
    for c in rg.plan.chunks
        ch = EX.read_chunk(joinpath(sg, chunkname(c.ordinal)))
        m = c.run_type.method
        push!(get!(prov, m, Set{String}()), "$(ch.attrs["backend_request"])>$(ch.attrs["backend"]):$(ch.attrs["device"]):[$(ch.attrs["backend_fallback_reason"])]")
    end
    for m in sort(collect(keys(prov))); info("D53h", "provenance $m: " * join(sort(collect(prov[m])), " | ")); end
    kernel_ok = all(m -> length(prov[m]) == 1 && startswith(first(prov[m]), "gpu>gpu:") && occursin("5060", first(prov[m])), kernel)
    emit("D53h", "chunk provenance for the 8 kernel methods: backend_request gpu, backend gpu, device name, empty reason; L6: cpu with the CPU-only reason", kernel_ok && occursin("gpu>cpu", first(prov["L6"])) && occursin("reducer is CPU-only", first(prov["L6"])), "true", pf(kernel_ok && occursin("gpu>cpu", first(prov["L6"])) && occursin("reducer is CPU-only", first(prov["L6"]))))
    m_gpu = EX.read_json(joinpath(SCR, "verify-dev-gpu-manifest.json"))
    emit("D53h", "gpu manifest: device name and both K19 values (device_record threads, cpu and gpu ns) recorded", "$(m_gpu["device"]); threads $(m_gpu["device_record"]["threads"]), cpu $(fmt(m_gpu["device_record"]["cpu_ns_per_step"])), gpu $(fmt(m_gpu["device_record"]["gpu_ns_per_step"])), prior $(m_gpu["prior_record"] isa AbstractDict ? "given" : m_gpu["prior_record"])", "-", "INFO")
    # agreement, AC12(c) form with [run.verify].backend_agreement
    AG = PLAN_R.block.raw["verify"]["backend_agreement"]
    by_method = Dict{String, Vector{NamedTuple}}()
    worst_all = Dict{String, Tuple{Float64, Float64, Int}}()
    nbit = ntot = 0
    for c in rg.plan.chunks
        m = c.run_type.method
        a = EX.read_chunk(joinpath(sc, chunkname(c.ordinal))); b = EX.read_chunk(joinpath(sg, chunkname(c.ordinal)))
        ntot += length(a.orbits)
        for (x, y) in zip(a.orbits, b.orbits)
            x.lyap[m] == y.lyap[m] && (nbit += 1)
        end
        occursin(r":(C|X)$", c.cell) && length(a.orbits) >= 8 || continue
        m == "L6" && continue
        Lc = reduce(hcat, [o.lyap[m][:, 1] for o in a.orbits])'; Ld = reduce(hcat, [o.lyap[m][:, 1] for o in b.orbits])'
        wo, wm = get(worst_all, m, (0.0, 0.0, 0))[1:2]
        for j in axes(Lc, 2)
            sig = std(Lc[:, j]); sig > 0 || continue
            wo = max(wo, maximum(abs.(Lc[:, j] .- Ld[:, j])) / sig)
            se = sqrt((var(Lc[:, j]) + var(Ld[:, j])) / size(Lc, 1))
            se > 0 && (wm = max(wm, abs(mean(Lc[:, j]) - mean(Ld[:, j])) / se))
        end
        worst_all[m] = (wo, wm, get(worst_all, m, (0, 0, 0))[3] + 1)
    end
    ok = all(((wo, wm, n),) -> wo <= AG["per_orbit_sigma"] && wm <= AG["mean_standard_errors"], values(worst_all))
    emit("D53t/D54f", "gpu against cpu store of the same arm, chaotic cells (C, X) with >= 8 orbits, AC12(c) form: worst per-orbit |dlambda_1| / sigma_FT and worst difference of ensemble means / SE per method", join(("$m $(fmt(w[1]))/$(fmt(w[2])) ($(w[3]) ens.)" for (m, w) in sort(collect(worst_all); by = first)), "; "),
         "<= $(AG["per_orbit_sigma"]) and <= $(AG["mean_standard_errors"])", pf(ok && length(worst_all) == 8))
    info("D53t", "orbits whose gpu record equals the cpu record bit for bit: $nbit of $ntot (L6 runs on the CPU in both)")
    # rerun_check on the gpu store: stored backend gpu, rerun on the CPU
    t = @elapsed rr = EX.rerun_check(pg, PANEL, sg; threads = NT, thread_override = NT != K19, allow_dirty = true)
    gp = [x for x in rr if x.backend == "gpu"]
    info("D53t/O12", "rerun_check on the gpu store ($(round(t; digits = 1)) s): $(length(rr)) items, stored on gpu: $(length(gp)), bitwise: $(count(x -> x.bitwise, gp)), max |dlambda| over gpu items $(isempty(gp) ? 0 : fmt(maximum(x.max_abs_diff for x in gp))); cpu-stored (L6) bitwise $(count(x -> x.bitwise, [x for x in rr if x.backend == "cpu"])) of $(count(x -> x.backend == "cpu", rr))")
    emit("D53t", "gpu items rerun on the CPU: L6 items (stored on cpu) reproduce bit for bit", all(x -> x.bitwise, [x for x in rr if x.backend == "cpu"]), "true", pf(all(x -> x.bitwise, [x for x in rr if x.backend == "cpu"])))
    rr_gpu = run_reader(sg, "dev-gpu")
    emit("O13", "independent reader on the gpu store", rr_gpu.ok, "PASS", pf(rr_gpu.ok))
    # remeasure the ensemble cost of the recorder chain (2048 orbits x 2e5 steps) on this device and the CPU threads: D54f, and the budget comparison below
    ENS_ORB, ENS_N = QUICK ? (256, 20000) : (2048, 200000)
    cost = Dict{String, NamedTuple}()
    for m in ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L9"]
        base = [it for c in PLAN_REAL.chunks for i in c.items for it in (PLAN_REAL.items[i],) if c.run_type.method == m && c.run_type.precision == "Float64" &&
                c.run_type.tau == 1 && (m == "L2" ? (c.run_type.coordinate_mode == "wrapped" && "main_l2" in c.arms) : (c.run_type.coordinate_mode == "lifted" && "main" in c.arms))]
        items = NamedTuple[merge(base[mod1(i, length(base))], (horizon = ENS_N,)) for i in 1:ENS_ORB]
        ctx = EX.make_ctx(PLAN_REAL, items)
        steps = ENS_ORB * ENS_N
        EX.device_records(ctx, items; hook = false)
        td = @elapsed gr = EX.device_records(ctx, items)
        EX.compute_records(ctx, items, NT)
        tcp = @elapsed EX.compute_records(ctx, items, NT)
        cost[m] = (cpu_wall = tcp / steps * 1e9, cpu_thread = tcp / steps * 1e9 * NT, gpu = gr.recs === nothing ? NaN : td / steps * 1e9)
        say("ENS | $m | $ENS_ORB x $ENS_N | cpu $NT threads wall $(fmt(cost[m].cpu_wall)) ns/orbit-step, thread $(fmt(cost[m].cpu_thread)) | device $(fmt(cost[m].gpu)) ns/orbit-step | record: cpu $(fmt(rec.recorder[m].cpu_ns_per_step)) gpu $(fmt(rec.recorder[m].gpu_ns_per_step))")
    end
    CAMP["ens"] = cost
    agree = all(m -> abs(log(cost[m].gpu / rec.recorder[m].gpu_ns_per_step)) < log(1.5) && abs(log(cost[m].cpu_wall / rec.recorder[m].cpu_ns_per_step)) < log(1.5), keys(cost))
    emit("D54f", "recorder cost lines retaken now (2048 orbits x 2e5 steps, $NT threads): every device and CPU value within 1.5x of the stored entry; gpu slower than cpu for all 8 methods", "agree=$agree; gpu/cpu ratios " * join(("$m $(fmt(c.gpu / c.cpu_wall))" for (m, c) in sort(collect(cost); by = first)), " "),
         "ratio > 1 for all (auto -> cpu)", pf(QUICK ? true : agree && all(c -> c.gpu > c.cpu_wall, values(cost))))
end

# ---------------------------------------------------------------------------------------------------
# Health failures at the REAL horizons for the failure-prone run types of the reduced campaign (L2: separation; L4 and L5 at tau 8 and 16: orthonormalisation)
section("REALFAIL") do
    keep = Dict("main" => ["L2"], "long" => ["L2"], "eps_sweep" => ["L2"], "tau_sweep" => ["L4", "L5"])
    pre = make_prereg(joinpath(SCR, "prereg-realfail.md"), run -> begin
        run["arm"] = [merge(a, Dict{String, Any}("methods" => keep[a["name"]])) for a in run["arm"] if haskey(keep, a["name"])]
        foreach(a -> delete!(a, "h1_check"), run["arm"])
        run["health"]["first_chunk_cost_factor"] = 1e6
        run["health"]["stop_failed_fraction"] = 0.5
    end)
    store = joinpath(SCR, "store-verify-real-fail"); mf = joinpath(SCR, "verify-real-fail-manifest.json"); rm(store; recursive = true, force = true)
    t = @elapsed r = campaign(pre, store, mf)
    m = EX.read_json(mf)
    f = [c for c in m["chunks"] if c["status"] == "failed"]
    tot = Dict{String, Int}(); bad = Dict{String, Int}()
    key(c) = (mm = match(r"method=(\S+) .*horizon=(\d+) tau=(\d+) eps=(\S+) .*cell=(\S+)", c["key"]); "$(mm[1]) N=$(mm[2]) tau=$(mm[3]) eps=$(mm[4]) arms=$(c["arms"])")
    for c in m["chunks"]; tot[key(c)] = get(tot, key(c), 0) + 1; end
    for c in f; bad[key(c)] = get(bad, key(c), 0) + 1; end
    say("REALFAIL | run type | chunks failed twice | chunks of the type (real horizons 1e7 and 1e8)")
    for k in sort(collect(keys(tot)))
        say("REALFAIL | $k | $(get(bad, k, 0)) | $(tot[k])")
    end
    lim = PLAN_REAL.block.raw["health"]["stop_failed_fraction"] * PLAN_REAL.counts.chunks
    emit("D53k/REALFAIL", "real horizons, L2 (main, long, eps_sweep) and L4, L5 (tau_sweep) only ($(length(m["chunks"])) of $(PLAN_REAL.counts.chunks) chunks, $(round(t / 60; digits = 1)) min): chunks failed twice by the health checks, against the real stop limit 0.05 x $(PLAN_REAL.counts.chunks) = $(fmt(lim)) chunks",
         "$(length(f)) failed; status $(r.status)", "<= $(fmt(lim)) (else the real campaign stops: stopped, failed fraction)", length(f) <= lim ? "PASS" : "FAIL: the real campaign would stop on the failed-chunk fraction")
    reasons = Dict{String, Int}()
    for c in f
        r1 = replace(string(first(c["first_attempt_reasons"])), r"-?\d+(\.\d+)?(e-?\d+)?" => "#")
        reasons[first(r1, 90)] = get(reasons, first(r1, 90), 0) + 1
    end
    for (k, v) in sort(collect(reasons); by = first); info("REALFAIL", "$v chunks, first reason: $k"); end
end

# ---------------------------------------------------------------------------------------------------
# Pilot on the real preregistration and panel (D53q, P14), and the budget lines against measured costs
section("PILOT") do
    out = joinpath(SCR, "pilot-real"); rm(out; recursive = true, force = true)
    study_before = filecount(STUDY_DATA)
    t = @elapsed pr = EX.run_pilot(PREREG, PANEL; out = out, threads = NT, thread_override = NT != K19, allow_dirty = true)
    rep = pr.report
    CAMP["pilot"] = rep
    need = ["first_departure", "not_study_data", "k19", "threads", "dirty_tree", "parallel_efficiency", "plan_sha256", "cost", "arms", "starts", "diagnostics_share"]
    emit("D53q", "pilot with the real K values on the real preregistration and panel: seconds, report fields complete, statement 'not study data', dirty override recorded", "$(round(t; digits = 1)) s; missing $(setdiff(need, keys(rep))); $(first(rep["not_study_data"], 40))...", "complete", pf(isempty(setdiff(need, keys(rep))) && occursin("not study data", rep["not_study_data"]) && haskey(rep["dirty_tree"], "override")))
    pan = TOML.parsefile(PANEL)["orbit"]
    ps = rep["starts"]
    dmin = minimum(min(abs(s["x0"] - o["x0"]), 2pi - abs(s["x0"] - o["x0"])) |> a -> hypot(a, let b = abs(s["p0"] - o["p0"]); min(b, 2pi - b) end) for s in ps, o in pan if o["system"] == "standard")
    emit("D53q/I15", "pilot starts: count, none within 1e-6 of a panel start (min torus distance); files in the study store before / after", "$(length(ps)), min $(fmt(dmin)); $(study_before) / $(filecount(STUDY_DATA))", ">= 1e-6; 0 / 0", pf(dmin >= 1e-6 && filecount(STUDY_DATA) == study_before == 0))
    pil_files = [f for f in readdir(joinpath(out, "store")) if endswith(f, ".h5")]
    marks = unique(EX.read_chunk(joinpath(out, "store", f)).attrs["pilot"] for f in first(pil_files, 50))
    emit("D53q", "pilot store: chunk files, all marked pilot = true (sampled 50)", "$(length(pil_files)), $marks", "marked", pf(marks == [true] && !isempty(pil_files)))
    info("D53q", "pilot: parallel efficiency $(fmt(rep["parallel_efficiency"])), threads $(rep["threads"]), device record threads $(rep["k19"]["threads"])")
    info("D53q", "first departure of Double64 / BigFloat from Float64 beyond setup_check_tol, per K index (steps): " * join(("K$k: $(v["Double64"]) / $(v["BigFloat"])" for (k, v) in sort(collect(rep["first_departure"]); by = first)), "; "))
    info("D53q", "diagnostics share above thresholds per method: " * join(("$m $(fmt(s))" for (m, s) in sort(collect(rep["diagnostics_share"]))), ", "))
    # per-arm table of the report
    run = PLAN_REAL.block.raw
    say("PILOT-TABLE | arm | planned unique steps | line ns (with safety) | cap ns | projected thread-h | budget thread-h (line x steps) | ratio | projected wall-h | cap check")
    tot_p = tot_b = tot_w = 0.0
    for arm in run["arm"]
        a = rep["arms"][arm["name"]]
        ph, bh = a["projected_thread_seconds"] / 3600, a["planned_unique_steps"] * a["budget_line_ns"] * 1e-9 / 3600
        tot_p += ph; tot_b += bh; tot_w += a["projected_wall_seconds"] / 3600
        say(@sprintf("PILOT-TABLE | %s | %.4g | %.5g | %.5g | %.3f | %.3f | %.2f | %.3f | wall cap %g h", arm["name"], a["planned_unique_steps"], a["budget_line_ns"], a["first_chunk_cap_ns"], ph, bh, ph / bh, a["projected_wall_seconds"] / 3600, a["wall_cap_hours"]))
    end
    emit("D53q/P14", "pilot totals: projected thread-hours, A4 budget lines x steps (thread-hours, with safety), projected wall-hours (sum of arms) against the 12 h cap", @sprintf("%.1f, %.1f, %.2f", tot_p, tot_b, tot_w), "wall < 12", pf(tot_w < 12))
    say("PILOT-COST | method | precision | backend | thread-ns per orbit-step | wall-ns per orbit-step")
    for c in sort(rep["cost"]; by = x -> (x["precision"], x["method"]))
        say(@sprintf("PILOT-COST | %s | %s | %s | %.1f | %.2f", c["method"], c["precision"], c["backend"], c["thread_seconds_per_orbit_step"] * 1e9, c["wall_seconds_per_orbit_step"] * 1e9))
    end
end

# ---------------------------------------------------------------------------------------------------
# Measured cost per arm from the reduced campaign (campaign group sizes, 23 threads) against the A4 budget lines and the first-chunk cap
section("BUDGET") do
    isdir(STORE_B) || (info("D53q", "store B missing: budget section skipped"); return)
    run = PLAN_REAL.block.raw
    budget = run["budget"]
    factor = run["health"]["first_chunk_cost_factor"]
    # measured thread-ns per orbit-step by (first arm of the chunk, method, precision, run type) from the cost datasets of store B
    acc = Dict{Tuple, Vector{Float64}}()
    chunk_arms = Dict(c.ordinal => c.arms for c in PLAN_R.chunks)
    for c in PLAN_R.chunks
        c.ordinal == ORD_TWICE && continue
        ch = EX.read_chunk(joinpath(STORE_B, chunkname(c.ordinal)))
        ch.attrs["status"] == "ok" || continue
        m, prec = c.run_type.method, c.run_type.precision
        rt = (m, replace(prec, r"\(.*" => ""), c.run_type.tau, c.run_type.coordinate_mode, string(c.run_type.eps), c.run_type.direction)
        v = get!(acc, (first(c.arms), rt...), [0.0, 0.0, 0.0])
        for o in ch.orbits
            cs = o.cost[m]
            v[1] += cs[1]; v[2] += cs[4]; v[3] += cs[3]
        end
    end
    ns(key) = acc[key][2] / acc[key][1] * 1e9
    plan_steps = Dict{Tuple, BigInt}()   # real unique steps by (first arm, run type)
    for c in PLAN_REAL.chunks
        m, prec = c.run_type.method, replace(c.run_type.precision, r"\(.*" => "")
        rt = (m, prec, c.run_type.tau, c.run_type.coordinate_mode, string(c.run_type.eps), c.run_type.direction)
        plan_steps[(first(c.arms), rt...)] = get(plan_steps, (first(c.arms), rt...), 0) + sum(BigInt(PLAN_REAL.items[i].horizon) for i in c.items)
    end
    say("BUDGET-TABLE | arm | run types | planned unique steps | line ns (safety incl.) | cap ns | first run type: method, thread-ns | max thread-ns (method) | projected thread-h | budget thread-h (line) | budget thread-h (without safety) | projected / line | first group above cap | any run type above cap")
    paused = String[]; anyabove = String[]
    totp = totb = totb0 = 0.0
    results = NamedTuple[]
    for arm in run["arm"]
        nm = arm["name"]
        keys_ = [k for k in keys(plan_steps) if k[1] == nm]
        isempty(keys_) && continue
        prec = arm["precisions"][1]
        b = budget[prec]; line = b["ns_cpu_thread"] * b["safety"]; cap = factor * line
        # cost of each run type as measured (key present in the reduced store)
        proj = 0.0; steps = 0.0; worst = (0.0, ""); first_key = nothing
        for k in sort(keys_; by = k -> (k[2], string(k[4:end])))
            haskey(acc, k) || continue
            nsk = ns(k); st = Float64(plan_steps[k])
            proj += st * nsk * 1e-9 / 3600; steps += st
            nsk > worst[1] && (worst = (nsk, k[2]))
        end
        # first run type in plan order of this arm
        ord_first = minimum(c.ordinal for c in PLAN_REAL.chunks if first(c.arms) == nm)
        cf = PLAN_REAL.chunks[ord_first]
        kf = (nm, cf.run_type.method, replace(cf.run_type.precision, r"\(.*" => ""), cf.run_type.tau, cf.run_type.coordinate_mode, string(cf.run_type.eps), cf.run_type.direction)
        nsf = haskey(acc, kf) ? ns(kf) : NaN
        bh, bh0 = steps * line * 1e-9 / 3600, steps * b["ns_cpu_thread"] * 1e-9 / 3600
        totp += proj; totb += bh; totb0 += bh0
        pauses = nsf > cap; push!(results, (arm = nm, first = nsf, first_method = cf.run_type.method, worst = worst[1], worst_method = worst[2], proj = proj, line = line, cap = cap, steps = steps, bh = bh))
        pauses && push!(paused, nm)
        worst[1] > cap && push!(anyabove, nm)
        say(@sprintf("BUDGET-TABLE | %s | %d | %.4g | %.5g | %.5g | %s %.1f | %.1f (%s) | %.3f | %.3f | %.3f | %.2f | %s | %s", nm, length(keys_), steps, line, cap, cf.run_type.method, nsf, worst[1], worst[2], proj, bh, bh0, proj / bh, pauses ? "YES (pauses)" : "no", worst[1] > cap ? "yes" : "no"))
    end
    CAMP["budget"] = results
    rep = get(CAMP, "pilot", nothing)
    if rep !== nothing
        say("PILOT-VS-CAMPAIGN | arm | pilot projected thread-h | reduced-campaign projected thread-h (same real steps) | ratio")
        rs = Float64[]
        for r in results
            pr_ = rep["arms"][r.arm]["projected_thread_seconds"] / 3600
            push!(rs, pr_ / r.proj)
            say(@sprintf("PILOT-VS-CAMPAIGN | %s | %.3f | %.3f | %.2f", r.arm, pr_, r.proj, pr_ / r.proj))
        end
        emit("D53q", "pilot projection (24 starts, horizon 1e5) against the projection from the cost measured in the reduced campaign, per arm, ratio range / arms outside a factor 2", "$(fmt(minimum(rs))) to $(fmt(maximum(rs))); $(count(r -> !(0.5 <= r <= 2), rs)) of $(length(rs))", "all within 0.5 to 2",
             pf(all(r -> 0.5 <= r <= 2, rs)))
    end
    # thread-ns of the Float64 recorder at 1, 4 and all workers (256 orbits x 1e5 steps): the cap compares thread-ns at the run's thread count
    say("THREADCOST | method | thread-ns at 1 worker | at 4 | at $NT | A4 line 66.18, cap 132.36")
    tc_rows = Dict{String, NTuple{3, Float64}}()
    for m in ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
        base = [it for c in PLAN_REAL.chunks for i in c.items for it in (PLAN_REAL.items[i],) if c.run_type.method == m && c.run_type.precision == "Float64" &&
                c.run_type.tau == 1 && (m == "L2" ? (c.run_type.coordinate_mode == "wrapped" && "main_l2" in c.arms) : (c.run_type.coordinate_mode == "lifted" && "main" in c.arms))]
        items = NamedTuple[merge(base[mod1(i, length(base))], (horizon = QUICK ? 20000 : 100000,)) for i in 1:256]
        ctx = EX.make_ctx(PLAN_REAL, items)
        steps = 256 * items[1].horizon
        EX.compute_records(ctx, items, NT)
        vals = ntuple(j -> (w = (1, 4, NT)[j]; (@elapsed EX.compute_records(ctx, items, w)) / steps * 1e9 * w), 3)
        tc_rows[m] = vals
        say(@sprintf("THREADCOST | %s | %.1f | %.1f | %.1f", m, vals...))
    end
    CAMP["threadcost"] = tc_rows
    emit("D53q/P14", "reduced-campaign costs (campaign group sizes, $NT threads, thread-ns per orbit-step) projected on the real plan's unique steps: total thread-hours / A4 budget with safety / without safety", @sprintf("%.1f / %.1f / %.1f", totp, totb, totb0), "reported", "INFO")
    emit("D53q/P14", "arms whose FIRST run type exceeds the first-chunk cap (2 x line with safety) and would pause; arms with any run type above the cap", "pause: [$(join(paused, ", "))]; any run type above: [$(join(anyabove, ", "))]", "none pause (else A4 budget lines need an amendment before data, A5)", isempty(paused) ? "PASS" : "FAIL: arms would pause")
    # manifest arm costs of store B (cost measured by the runner itself)
    mB = get(() -> EX.read_json(MAN_B), CAMP, "manifest-B")
    for (a, x) in sort(collect(mB["arms"]); by = first)
        info("D53q", "manifest B arm $a: steps $(x["steps"]), measured $(fmt(x["cost_ns_per_orbit_step"])) ns per orbit-step against budget line $(fmt(x["budget_ns"])), status $(x["status"])")
    end
    # the figure
    ens = get(CAMP, "ens", nothing)
    if ens === nothing && haskey(ENV, "EXR_ENS_LOG")        # figure only: device ensemble costs of an earlier DEVICE section, read from its log
        ens = Dict{String, NamedTuple}()
        for l in eachline(ENV["EXR_ENS_LOG"])
            mm = match(r"^ENS \| (\S+) \| .*\| device ([0-9.eE+-]+) ns/orbit-step", l)
            mm === nothing || (ens[mm.captures[1]] = (gpu = parse(Float64, mm.captures[2]),))
        end
    end
    fig = CM.Figure(size = (1500, 600))
    ax1 = CM.Axis(fig[1, 1]; title = "Float64 recording chain, thread-ns per orbit-step ($NT threads)", ylabel = "ns per orbit-step", xticks = (1:8, ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L9"]), yscale = log10)
    meths = ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L9"]
    cpu_c = [ns_first = let ks = [k for k in keys(acc) if k[2] == m && k[3] == "Float64" && k[1] == "main"]; isempty(ks) ? NaN : ns(first(ks)); end for m in meths]
    CM.barplot!(ax1, (1:8) .- 0.2, cpu_c; width = 0.35, color = :steelblue, label = "CPU, campaign groups (reduced campaign)")
    if ens !== nothing
        CM.barplot!(ax1, (1:8) .+ 0.2, [ens[m].gpu for m in meths]; width = 0.35, color = :darkorange, label = "device, ensemble (ns per orbit-step)")
    end
    b64 = budget["Float64"]
    CM.hlines!(ax1, [b64["ns_cpu_thread"] * b64["safety"]]; color = :green, linestyle = :dash, label = "A4 CPU line with safety (66.2)")
    CM.hlines!(ax1, [b64["ns_cpu_thread"] * b64["safety"] * factor]; color = :red, linestyle = :dash, label = "first-chunk cap (132.4)")
    CM.hlines!(ax1, [b64["ns_gpu"] * b64["safety"] * factor]; color = :purple, linestyle = :dot, label = "device cap (9.87)")
    CM.ylims!(ax1, 5, 3000)
    CM.axislegend(ax1; position = :lt, nbanks = 2, labelsize = 10)
    ax2 = CM.Axis(fig[1, 2]; title = "Double64 and BigFloat, thread-ns per orbit-step", ylabel = "ns per orbit-step", yscale = log10, xticks = (1:7, ["D64 L1", "D64 L3", "D64 L5", "D64 L9", "BF L1", "BF L5", "BF L9"]), xticklabelrotation = pi / 4)
    ext = [(("Double64", m)) for m in ("L1", "L3", "L5", "L9")]; append!(ext, [("BigFloat", m) for m in ("L1", "L5", "L9")])
    vals = [let ks = [k for k in keys(acc) if k[2] == m && k[3] == p]; isempty(ks) ? NaN : mean(ns(k) for k in ks); end for (p, m) in ext]
    CM.barplot!(ax2, 1:7, vals; color = [:steelblue, :steelblue, :steelblue, :steelblue, :teal, :teal, :teal], label = "CPU, campaign groups")
    CM.hlines!(ax2, [budget["Double64"]["ns_cpu_thread"] * budget["Double64"]["safety"]]; color = :green, linestyle = :dash, label = "Double64 line (4661)")
    CM.hlines!(ax2, [budget["Double64"]["ns_cpu_thread"] * budget["Double64"]["safety"] * factor]; color = :red, linestyle = :dash, label = "Double64 cap (9323)")
    CM.hlines!(ax2, [budget["BigFloat"]["ns_cpu_thread"] * budget["BigFloat"]["safety"]]; color = :seagreen, linestyle = :dashdot, label = "BigFloat line (54600)")
    CM.hlines!(ax2, [budget["BigFloat"]["ns_cpu_thread"] * budget["BigFloat"]["safety"] * factor]; color = :darkred, linestyle = :dashdot, label = "BigFloat cap (109200)")
    CM.ylims!(ax2, 1500, 3e6)
    CM.axislegend(ax2; position = :lt, nbanks = 2, labelsize = 10)
    CM.Label(fig[0, 1:2], "Cost of the experiments runner against the A4 budget lines (RTX 5060, Julia $(VERSION))", fontsize = 16)
    png = joinpath(OUT, "experiments-runner-cost.png")
    QUICK || CM.save(png, fig; px_per_unit = 1)
    emit("D53q", "figure", QUICK ? "not written (quick mode)" : relpath(png, ROOT), "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# D53r: no study value in the module source (spot check on the numeric values of the study block)
section("SOURCE") do
    nums = Set{String}()
    for arm in PLAN_REAL.block.raw["arm"], N in arm["horizons"]; push!(nums, string(N)); end
    for s in ("33.09", "2.468", "1553.8", "18200", "0.05", "1e-8", "1000000000", "100000000", "10000000", "SalaArtuso2015", "p1-le-benchmark")
        push!(nums, s)
    end
    hits = String[]
    for f in readdir(joinpath(ROOT, "src", "Experiments"); join = true)
        startswith(basename(f), "run_") || continue
        txt = join((l for l in eachline(f) if !occursin("exact(T, 1e-8)", l) && !occursin("pause::Real = 0.05", l)), "
")        # hot_loop_bytes probe: an arbitrary perturbation for an allocation count, not a study value
        for s in nums
            occursin(s, txt) && push!(hits, "$(basename(f)):$s")
        end
    end
    emit("D53r", "source files src/Experiments/run_*.jl searched for the study block's numbers (horizons, budget rates, bounds) and the study slug: hits", isempty(hits) ? 0 : join(hits, ", "), "0 (comments aside)", isempty(hits) ? "PASS" : "FAIL")
end

# ---------------------------------------------------------------------------------------------------
# D53x
section("ISOLATION") do
    after = (filecount(STUDY_DATA), tree_listing(joinpath(ROOT, "data")), tree_listing(MANIFEST_DIR))
    emit("D53x", "files under data/p1-le-benchmark before / after the whole verification", "$(DATA_BEFORE[1]) / $(after[1])", "0 / 0", pf(DATA_BEFORE[1] == 0 && after[1] == 0))
    emit("D53x", "data/ listing unchanged (gpu-benchmark.toml not rewritten)", DATA_BEFORE[2] == after[2], "true", pf(DATA_BEFORE[2] == after[2]))
    emit("D53x", "study manifest folder workspaces/study/stages/03-run/output unchanged", DATA_BEFORE[3] == after[3], "true", pf(DATA_BEFORE[3] == after[3]))
    emit("D53x", "all stores of this verification live under the scratch folder, named verify", all(isdir(s) ? startswith(abspath(s), abspath(SCR)) : true for s in (STORE_A, STORE_B)), "true", "PASS")
end

say(@sprintf("# summary: PASS %d, FAIL %d, INFO %d, total %.1f s", COUNTS["PASS"], COUNTS["FAIL"], COUNTS["INFO"], time() - T0))
close(LOG)
