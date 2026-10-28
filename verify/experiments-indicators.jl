# Stage 04 verification of the recording of the chaos indicators by the experiments runner (Amendment D86, module-build experiments). CPU part. Run from the repo root:
#   julia --project=verify -t 23 verify/experiments-indicators.jl
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md (Amendment D86: oracles O17 to O20, O22; criteria D86b [full], D86s (V), D86w (V), D86y,
# and the (V) checks of D86e to D86v on real panel orbits). The device part is verify/experiments-indicators-device.jl.
# EXI_ONLY="SETUP,CAMPAIGN" runs only the named sections (they reuse the stores of the scratch folder when EXI_REUSE=1). One printed line per check:
#   <id> | quantity | value | bound | PASS/FAIL/INFO
# Scratch folder (all stores, reduced preregistration and panel): EXI_SCRATCH. The study data folder data/ is only read (D86y).
include(joinpath(@__DIR__, "experiments-indicators-prep.jl"))
using HDF5, HDF5_jll

const REUSE = get(ENV, "EXI_REUSE", "") == "1"
const ONLY = filter(!isempty, split(get(ENV, "EXI_ONLY", ""), ","))
mkpath(OUT)
const LOG = open(get(ENV, "EXI_LOG", joinpath(OUT, "experiments-d86-verify-log.txt")), isempty(ONLY) ? "w" : "a")
const Tg = Chaoslab.Tangent

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
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 1500)], "-", "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

# ---------------------------------------------------------------------------------------------------
const T0 = time()
say("# experiments D86 verification (CPU part), julia $(VERSION), threads $NT, cpu $(Sys.CPU_NAME), HDF5.jl $(pkgversion(Chaoslab.Experiments.HDF5)), $(Sys.total_memory() / 2^30) GiB host memory")
say("# prereg $(relpath(PREREG, ROOT)), panel $(relpath(PANEL, ROOT)), block $(relpath(BLOCK, ROOT)), scratch $SCR")
const GIT_HEAD = strip(read(`git -C $ROOT rev-parse HEAD`, String))
const GIT_DIRTY = strip(read(`git -C $ROOT status --short`, String))
say("# git head $GIT_HEAD, tracked or untracked changes: ", isempty(GIT_DIRTY) ? "none" : replace(GIT_DIRTY, "\n" => "; "))
say("# workaround for finding F1 active (EXI_PATCH): ", get(ENV, "EXI_PATCH", "1") == "1")

const STUDY_DATA = joinpath(ROOT, "data", "p1-le-benchmark")
const DATA_BEFORE = (filecount(STUDY_DATA), tree_listing(joinpath(ROOT, "data")))
say("# D86y before: files under data/p1-le-benchmark = $(DATA_BEFORE[1]) (the finished P1 campaign); entries under data/ = $(length(DATA_BEFORE[2]))")

const PP = prepare()
const PLAN = EX.plan_run(PP.prereg, PP.panel)
const NCH = length(PLAN.chunks)
const K19 = PLAN.block.raw["k19"]["threads"]
const ITEM_OF = Dict(it.hash => it for it in PLAN.items)
const SPECS = Dict(s["id"] => s for s in PLAN.block.raw["indicator"])
const ORD_FAIL = first(c.ordinal for c in PLAN.chunks if c.run_type.method == "tg" && c.run_type.precision == "Float64" && c.run_type.coordinate_mode == "wrapped")
const STORE_A, MAN_A = joinpath(SCR, "store-A"), joinpath(SCR, "manifest-A.json")
const STORE_B, MAN_B = joinpath(SCR, "store-B"), joinpath(SCR, "manifest-B.json")
const STORE_C, MAN_C = joinpath(SCR, "store-C"), joinpath(SCR, "manifest-C.json")
const STORE_T, MAN_T = joinpath(SCR, "store-T4"), joinpath(SCR, "manifest-T4.json")
const CAP_FIXED = 1000
inject() = nan_fault(ORD_FAIL, "tg"; attempts = (1, 2))
function campaign(store, manifest; hooks = EX.RunHooks(), kw...)
    return EX.run_campaign(PP.prereg, PP.panel; store = store, manifest = manifest, threads = NT, thread_override = NT != K19, allow_dirty = true, hooks = hooks, kw...)
end
const CAMP = Dict{String, Any}()

# ---------------------------------------------------------------------------------------------------
# helpers: element types, systems, stand-alone references (the verified reducers run directly on the identical orbit)
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
reducer(m; replace = true) = m == "SALI" ? Chaoslab.SALI() : getfield(Chaoslab, Symbol(m))(; replace = replace)
function member_cols(m, fin)
    m == "SALI" && return Float64[fin.sali, fin.L]
    m == "FLI" && return Float64[fin.fli, fin.fli1]
    m == "MEGNO" && return Float64[v for c in 1:2 for v in (fin.Y[c], fin.Ybar[c], fin.Ybar20[c])]
    m == "RLI" && return Float64[fin.DL, fin.RLI]
    error(m)
end
function stand_alone(it, m, n)
    return with_prec(it.precision) do T
        sys, x = sys_x(it, T)
        r = reducer(m)
        _, _, st = Tg.propagate_orbit(sys, x, frame(it, T), n, it.tau, (r,))
        member_cols(m, Tg.reducer_final(r, st[1]))
    end
end
function stretch_final(it, spec, n)
    return with_prec(it.precision) do T
        sys, x = sys_x(it, T)
        rng = spec["range"] == "derived" ? Chaoslab.stretching_range(sys) : Tuple(spec["range"])
        r = Chaoslab.StretchingNumbers(; range = rng, bins = spec["bins"], interval = spec["interval"], min_per_cell = spec["min_per_cell"])
        ws = zeros(UInt32, Tg.reducer_workspace(r, T, Val(2), Val(2)))
        _, _, st = Tg.propagate_orbit(sys, x, frame(it, T), n, it.tau, (r,); workspace = ws)
        (fin = Tg.reducer_final(r, st[1], ws), lo = r.lo, hi = r.hi)
    end
end
function pair_ref(it, spec, xi0, n, member)
    return with_prec(it.precision) do T
        sys, _ = sys_x(it, T)
        x = SVector(Chaoslab.exact(T, it.x0), Chaoslab.exact(T, it.p0))
        u0 = Chaoslab.initial_deviations(parse(Int, it.run_seed), it.orbit, Val(2), Val(2), T)       # the seeded xi_0 in the element type (the file stores it as binary64)
        u = SVector{2, T}(u0[1, 1], u0[2, 1])
        Float64.(u) == SVector{2, Float64}(xi0...) || error("the seeded xi_0 of orbit $(it.orbit) is not the stored xi0")
        ps = Chaoslab.pair_start(sys, x, u, spec["separation"]; orbit = it.orbit)
        W = SMatrix{4, 2, T}(u[1], u[2], zero(T), zero(T), zero(T), zero(T), u[1], u[2])
        r = reducer(member)
        _, _, st = Tg.propagate_orbit(Chaoslab.PairMap(sys), ps.z, W, n, it.tau, (r,))
        (cols = member_cols(member, Tg.reducer_final(r, st[1])), start = ps)
    end
end
colvec(t, name) = t.series[:, findfirst(==(name), t.columns)]
function grid_of(it)
    gr = PLAN.block.raw["grid"]
    return EX.record_grid(it.horizon, it.tau; points_per_decade = gr["points_per_decade"], first_steps = gr["first_steps"], extra_fractions = gr["extra_fractions"])
end
chunk_files(store) = sort(filter(f -> occursin(r"^chunk-\d{5}\.h5$", f), readdir(store)))
function load_store(store)
    return Dict(parse(Int, f[7:11]) => EX.read_chunk(joinpath(store, f)) for f in chunk_files(store))
end
isbitseq(a, b) = length(a) == length(b) && all(isequal.(a, b))
include(joinpath(@__DIR__, "experiments-indicators-series.jl"))      # plain_series: the plain iteration of the map (O17), at top level for the world age

# ---------------------------------------------------------------------------------------------------
section("SETUP") do
    emit("D86y", "reduced panel: drawn orbits (60 per K for 6 K), controls, orbits in S_L (one per cell)", "$(PP.info.selected), $(PP.info.controls), $(PP.info.sl)", "360, 13, one per populated cell",
         pf(PP.info.selected == 360 && PP.info.controls == 13))
    emit("D86y", "reduced block: spec kinds present (one spec of each kind)", join(sort([s["kind"] for s in values(SPECS)]), ","), "all six kinds", pf(sort([s["kind"] for s in values(SPECS)]) == sort(["tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff"])))
    emit("D86y", "plan of the reduced block: items, chunks", "$(PLAN.counts.items), $(PLAN.counts.chunks)", "-", "INFO")
    for (id, v) in sort(collect(PLAN.steps.specs); by = first)
        info("D86y", "plan, spec $id ($(SPECS[id]["kind"])): items $(v.items), planned orbit-steps $(v.steps)")
    end
    emit("D86y", "windows of the window kinds and the horizons of their arms", "fq/bk 250,1000 -> 2000; pa 250,1000 -> 2001", "horizon = tau x run length (R11 item 3)",
         pf(all(a["horizons"] == [2000] for a in PLAN.block.raw["arm"] if a["name"] in ("win", "win_d")) && all(a["horizons"] == [2001] for a in PLAN.block.raw["arm"] if a["name"] in ("pha", "pha_d"))))
    emit("D86c", "orbits without an item for a spec (the Y5 control has no torus: exp_i_q, cos, sin undefined), by spec", string(PLAN.not_applicable), "fq: [883], bk: [883]", pf(PLAN.not_applicable == Dict("fq" => [883], "bk" => [883])))
    p2 = EX.plan_run(PP.prereg, PP.panel)
    emit("D86c", "plan hash identical on a second planning", PLAN.hash == p2.hash ? "equal" : "different", "equal", pf(PLAN.hash == p2.hash))
    info("D86y", "injected failure: chunk $ORD_FAIL = $(PLAN.chunks[ORD_FAIL].key)")
    info("D86y", "chunks by spec: " * join(("$id: $(count(c -> c.run_type.method == id, PLAN.chunks))" for id in ["L1"; sort(collect(keys(SPECS)))]), ", "))
end

# ---------------------------------------------------------------------------------------------------
# Finding F1: unpatched Double64 stretching chunk
section("F1") do
    store, mf = joinpath(SCR, "store-F1"), joinpath(SCR, "manifest-F1.json")
    logf = joinpath(SCR, "f1.log")
    t = @elapsed ok = success(pipeline(addenv(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t $NT $(joinpath(@__DIR__, "experiments-indicators-campaign.jl")) $store $mf f1`, "EXI_PATCH" => "0");
                                       stdout = logf, stderr = logf))
    txt = read(logf, String)
    m = match(r"ERROR: LoadError: (MethodError: no method matching unsafe_convert\(::Type\{Ptr\{Nothing\}\}, ::Double64\))", txt)
    nfiles = length(chunk_files(store))
    emit("F1", "unpatched campaign: process exit ok / the first error / chunk files written before it ($t s)", "$ok / $(m === nothing ? "other: " * first(txt, 200) : m.captures[1]) / $nfiles",
         "FINDING: the campaign aborts with an exception at the first Double64 stretching chunk (run_indicator_store.jl:86, _h5val passes a Double64 to HDF5)", m !== nothing && !ok ? "INFO" : "FAIL")
    line = match(r"@ Chaoslab.Experiments .*run_indicator_store.jl:\d+", txt)
    info("F1", "stack: " * (line === nothing ? "?" : line.match))
    rm(store; recursive = true, force = true)
end

# ---------------------------------------------------------------------------------------------------
# D86b [full]: the real P1 preregistration and panel give the plan and chunks of the code before the amendment
section("D86B") do
    h5 = h5open(joinpath(STUDY_DATA, "chunk-00001.h5"), "r")
    cv = try read_attribute(h5, "code_version") finally close(h5) end
    man = EX.read_json(joinpath(ROOT, "workspaces", "study", "stages", "03-run", "output", "p1-le-benchmark-manifest.json"))
    defaults = joinpath(ROOT, "_config", "numerics-defaults.toml")
    info("D86b", "code_version of the running P1 campaign (attribute of data/p1-le-benchmark/chunk-00001.h5, read only): $cv; manifest plan_sha256 $(man["plan_sha256"]), chunks $(length(man["chunks"]))")
    wt = joinpath(SCR, "old-code")
    rm(wt; recursive = true, force = true)
    run(`git -C $ROOT worktree prune`)
    run(`git -C $ROOT worktree add --detach $wt ba21bcc`)
    try
        names = split(read(`git -C $ROOT diff --name-only 1c3de19 ba21bcc -- src ext Project.toml Manifest.toml Manifest-v1.13.toml`, String), '
'; keepempty = false)
        code = [n for n in names if endswith(n, ".jl") || n in ("Project.toml", "Manifest.toml", "Manifest-v1.13.toml")]     # what code_version hashes: *.jl under src and ext, Project.toml, the lock
        emit("D86b", "the code of commit ba21bcc (worktree: *.jl under src and ext, Project.toml, the lock) against 1c3de19, the commit before Amendment D86 touched src: differing code files (documentation files in src: $(length(names) - length(code)), not code)", length(code), "0", pf(isempty(code)))
        info("D86b", "differing files of src between ba21bcc and HEAD: " * strip(last(split(read(`git -C $ROOT diff --stat ba21bcc HEAD -- src`, String), '\n'; keepempty = false))))
        dd = read(`git -C $ROOT diff --stat ba21bcc HEAD -- _config/numerics-defaults.toml`, String)
        info("D86b", "_config/numerics-defaults.toml changed between ba21bcc and HEAD: $(isempty(strip(dd)) ? "no" : "yes")")
        childplan = joinpath(@__DIR__, "experiments-indicators-plan.jl")
        out_old, out_new = joinpath(SCR, "plan-old.txt"), joinpath(SCR, "plan-new.txt")
        t1 = @elapsed ok1 = success(pipeline(`$(Base.julia_cmd()) --project=$(joinpath(wt, "verify")) -t 4 $childplan $PREREG $PANEL $cv $defaults $out_old`; stdout = joinpath(SCR, "plan-old.log"), stderr = joinpath(SCR, "plan-old.log")))
        ok1 || (say(read(joinpath(SCR, "plan-old.log"), String)); error("child plan of the old code failed"))
        t2 = @elapsed ok2 = success(pipeline(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t 4 $childplan $PREREG $PANEL $cv $defaults $out_new`; stdout = joinpath(SCR, "plan-new.log"), stderr = joinpath(SCR, "plan-new.log")))
        ok2 || (say(read(joinpath(SCR, "plan-new.log"), String)); error("child plan of the new code failed"))
        old, new = readlines(out_old), readlines(out_new)
        info("D86b", "plans written by child processes: old code (worktree ba21bcc) $(round(t1; digits = 1)) s, current code $(round(t2; digits = 1)) s; $(length(old)) lines each")
        head(l) = Dict(split(x, ' ')[1] => join(split(x, ' ')[2:end], ' ') for x in l if !startswith(x, "item") && !startswith(x, "chunk"))
        ho, hn = head(old), head(new)
        emit("D86b [full]", "plan_sha256 of the real P1 preregistration and panel: old code / current code / running campaign's manifest", "$(ho["plan_sha256"][1:16]) / $(hn["plan_sha256"][1:16]) / $(man["plan_sha256"][1:16])", "all equal",
             pf(ho["plan_sha256"] == hn["plan_sha256"] == man["plan_sha256"]))
        emit("D86b [full]", "plan text hash, block hash, counts (items, chunks): old equals current", "$(ho["plan_text_sha256"] == hn["plan_text_sha256"]), $(ho["block_sha256"] == hn["block_sha256"]), $(hn["counts"])", "true, true, 15569 1363",
             pf(ho["plan_text_sha256"] == hn["plan_text_sha256"] && ho["block_sha256"] == hn["block_sha256"] && ho == hn && hn["counts"] == "15569 1363"))
        iold, inew = [l for l in old if startswith(l, "item")], [l for l in new if startswith(l, "item")]
        emit("D86b [full]", "item hashes (all work items): old code against current code, differences", count(i -> iold[i] != inew[i], eachindex(iold)) + abs(length(iold) - length(inew)), "0 of $(length(iold))", pf(iold == inew))
        cold, cnew = [l for l in old if startswith(l, "chunk")], [l for l in new if startswith(l, "chunk")]
        mchunks = ["chunk\t$(c["ordinal"])\t$(c["key"])\t$(c["hash"])" for c in man["chunks"]]
        nd_oc = count(i -> cold[i] != cnew[i], eachindex(cold)) + abs(length(cold) - length(cnew))
        nd_m = count(i -> mchunks[i] != cnew[i], eachindex(mchunks)) + abs(length(mchunks) - length(cnew))
        emit("D86b [full]", "ordinals, keys and hashes of all $(length(cnew)) chunks: old code against current code, differences", nd_oc, "0", pf(nd_oc == 0 && length(cnew) == 1363))
        emit("D86b [full]", "ordinals, keys and hashes of all chunks: current code against the manifest of the running (finished) P1 campaign, differences", nd_m, "0", pf(nd_m == 0))
        info("D86b", "no chunk key of an old method has a `kind=` or `spec=` field: $(!any(l -> occursin(" kind=", l) || occursin(" spec=", l), cnew))")
    finally
        run(`git -C $ROOT worktree remove --force $wt`)
        run(`git -C $ROOT worktree prune`)
        info("D86b", "worktree removed: $(!isdir(wt)); git worktree list: $(replace(strip(read(`git -C $ROOT worktree list`, String)), '\n' => "; "))")
    end
end

# ---------------------------------------------------------------------------------------------------
# D86y: the reduced campaign. C (child, wall time and peak memory), A (killed twice, resumed), B (uninterrupted), T4 (child, 4 threads and capacity 7)
function child_campaign(store, mf, tag, nthreads; cap = nothing)
    logf = joinpath(SCR, "child-$tag.log")
    args = cap === nothing ? String[] : [string(cap)]
    t = @elapsed ok = success(pipeline(`$(Base.julia_cmd()) --project=$VERIFY_ENV -t $nthreads $(joinpath(@__DIR__, "experiments-indicators-campaign.jl")) $store $mf $tag $args`; stdout = logf, stderr = logf))
    txt = read(logf, String)
    m = match(r"CAMPAIGN (\S+) status=(.*?) chunks=(\d+) ran=(\d+) failed=(\d+) threads=(\d+) wall_s=([\d.]+) maxrss_before_bytes=(\d+) maxrss_bytes=(\d+)", txt)
    return (ok = ok, seconds = t, log = txt, m = m)
end
section("CAMPAIGN") do
    if !(REUSE && isfile(MAN_B) && isfile(MAN_A) && isfile(MAN_C))
        # --- C: uninterrupted, own process, 23 threads, capacity read from the free memory: wall time and peak memory of the campaign
        rc = child_campaign(STORE_C, MAN_C, "C", K19)
        m = rc.m
        m === nothing && (say(rc.log); error("child campaign C gave no result line"))
        CAMP["C"] = (wall = parse(Float64, m.captures[7]), rss = parse(Int, m.captures[9]), status = m.captures[2], failed = parse(Int, m.captures[5]), ran = parse(Int, m.captures[4]), total = rc.seconds)
        emit("D86y", "uninterrupted reduced campaign in its own process ($(m.captures[6]) threads, cold start, capacity from the free memory): status, chunks run, chunks failed twice", "$(m.captures[2]), $(m.captures[4]) of $(m.captures[3]), $(m.captures[5])",
             "complete with failed chunks, all chunks, exactly the injected one", pf(m.captures[2] == "complete with failed chunks" && parse(Int, m.captures[4]) == NCH && parse(Int, m.captures[5]) == 1))
        emit("D86y", "wall time of the campaign (run_campaign only) / of the whole process including Julia start, package load and compilation", "$(round(CAMP["C"].wall; digits = 1)) s / $(round(rc.seconds; digits = 1)) s", "< 900 s on 23 threads (proposed 15 min)", pf(rc.seconds < 900))
        emit("D86y", "peak host memory of that process (Sys.maxrss, peak working set) / before the campaign", "$(round(CAMP["C"].rss / 2^30; digits = 2)) GiB / $(round(parse(Int, m.captures[8]) / 2^30; digits = 2)) GiB", "< 8 GiB (proposed)", pf(CAMP["C"].rss < 8 * 2^30))
        # --- A: session 1 killed after kill1 chunks, session 2 killed in the data of the next chunk, session 3 resumes
        rm(STORE_A; recursive = true, force = true); rm(STORE_B; recursive = true, force = true); rm(MAN_A; force = true); rm(MAN_B; force = true)
        kill1 = max(5, div(NCH, 5))
        t1 = @elapsed k1 = killed(() -> campaign(STORE_A, MAN_A; capacity = CAP_FIXED, hooks = EX.RunHooks(mutate = inject(), kill = kill_after(kill1))))
        done1 = length(chunk_files(STORE_A))
        emit("D86y", "session 1 killed after $kill1 chunks (injected kill): killed, final files / temporary files left", "$k1: $done1 / $(count(endswith(".tmp"), readdir(STORE_A)))", "killed", pf(k1 && done1 >= kill1))
        missing_ord = first(o for o in 1:NCH if !isfile(joinpath(STORE_A, chunkname(o))))
        kind_of_kill = PLAN.chunks[missing_ord].run_type.method
        t2 = @elapsed k2 = killed(() -> campaign(STORE_A, MAN_A; capacity = CAP_FIXED, hooks = EX.RunHooks(mutate = inject(), kill = kill_during(missing_ord))))
        tmps = filter(endswith(".tmp"), readdir(STORE_A))
        emit("D86y", "session 2 killed during the data of chunk $missing_ord (spec $kind_of_kill): no final file for it, one temporary file left (I13)", "$k2: final=$(isfile(joinpath(STORE_A, chunkname(missing_ord)))) tmp=$tmps", "killed, no final file, one .tmp",
             pf(k2 && !isfile(joinpath(STORE_A, chunkname(missing_ord))) && length(tmps) == 1))
        CAMP["A_kill"] = run_reader(STORE_A, "A-after-kill")
        emit("O20", "independent reader on the store after the two kills: every final file complete and valid, the .tmp ignored", CAMP["A_kill"].ok, "reader PASS", pf(CAMP["A_kill"].ok))
        t3 = @elapsed rA = campaign(STORE_A, MAN_A; capacity = CAP_FIXED, hooks = EX.RunHooks(mutate = inject()))
        emit("D86y", "interrupted campaign resumed to the end (session 3): status, chunks run in session 3, skipped, failed", "$(rA.status): ran $(length(rA.ran)), skipped $(length(rA.skipped)), failed $(length(rA.failed)), $(round(t3; digits = 1)) s",
             "complete with failed chunks; no new failure (the injected chunk $ORD_FAIL was written in session 1, see the manifest check below)", pf(rA.status == "complete with failed chunks" && isempty(rA.failed)))
        # --- B: uninterrupted, same injection, same fixed capacity
        tB = @elapsed rB = campaign(STORE_B, MAN_B; capacity = CAP_FIXED, hooks = EX.RunHooks(mutate = inject()))
        emit("D86y", "uninterrupted campaign B in the main process (warm): status, chunks run, failed, wall s", "$(rB.status): ran $(length(rB.ran)), failed $(rB.failed), $(round(tB; digits = 1)) s",
             "complete with failed chunks; failed = the injected chunk only", pf(rB.status == "complete with failed chunks" && rB.failed == [ORD_FAIL] && length(rB.ran) == NCH))
        before = tree_listing(STORE_A)
        r2 = campaign(STORE_A, MAN_A; capacity = CAP_FIXED, hooks = EX.RunHooks(mutate = inject()))
        emit("D53o/D86t", "second run on the finished campaign: chunks run, files modified", "$(length(r2.ran)), $(before == tree_listing(STORE_A) ? 0 : 1)", "0, 0", pf(isempty(r2.ran) && before == tree_listing(STORE_A)))
        # --- T4: 4 threads and capacity 7 in a child process (execution invariance at campaign scale, I24)
        rt = child_campaign(STORE_T, MAN_T, "T4", 4; cap = 7)
        rt.m === nothing && (say(rt.log); error("child campaign T4 gave no result line"))
        emit("I24/D86q", "campaign on 4 threads with launch capacity 7 (own process): status, failed", "$(rt.m.captures[2]): failed $(rt.m.captures[5])", "same status as B", pf(rt.m.captures[2] == "complete with failed chunks" && parse(Int, rt.m.captures[5]) == 1))
        CAMP["A_t"] = (t1, t2, t3)
    else
        info("D86y", "REUSE: stores A, B, C and T4 taken from the scratch folder")
        m = match(r"CAMPAIGN (\S+) status=(.*?) chunks=(\d+) ran=(\d+) failed=(\d+) threads=(\d+) wall_s=([\d.]+) maxrss_before_bytes=(\d+) maxrss_bytes=(\d+)", read(joinpath(SCR, "child-C.log"), String))
        CAMP["C"] = (wall = parse(Float64, m.captures[7]), rss = parse(Int, m.captures[9]))
    end
    for (tag, mp) in (("A", MAN_A), ("B", MAN_B), ("C", MAN_C), ("T4", MAN_T))
        st = EX.read_json(mp)["status"]
        failed = [c["ordinal"] for c in EX.read_json(mp)["chunks"] if c["status"] == "failed"]
        emit("D86y", "campaign $tag: manifest status and failed chunks (the campaign stops only by the injected cause)", "$st, failed $failed", "complete with failed chunks, [$ORD_FAIL]", pf(st == "complete with failed chunks" && failed == [ORD_FAIL]))
    end
    nfiles = (A = length(chunk_files(STORE_A)), B = length(chunk_files(STORE_B)), C = length(chunk_files(STORE_C)), T4 = length(chunk_files(STORE_T)))
    emit("D86y", "chunk files in stores A, B, C, T4", string(nfiles), "$NCH each", pf(all(==(NCH), values(nfiles))))
    ft = read(joinpath(SCR, "child-C.log"), String)
    info("D86y", "child C result line: " * first(filter(l -> startswith(l, "CAMPAIGN"), split(ft, '\n'))))
end

# ---------------------------------------------------------------------------------------------------
# O20 and D86s (V): the independent reader on every store; the module's read-back against it; I14
function compare_readers(store, R)
    nfile = nattr = ndata = nbad = 0
    bad = String[]
    badd(msg) = (nbad += 1; length(bad) < 8 && push!(bad, msg))
    for fn in chunk_files(store)
        nfile += 1
        ch = EX.read_chunk(joinpath(store, fn))
        for (nm, v) in ch.attrs
            if nm == "provenance"
                for (mem, d) in v, (k, x) in d
                    nattr += 1
                    get(R.pattrs, (fn, mem, k), nothing) == attr_str(x) || badd("$fn provenance/$mem@$k")
                end
                continue
            end
            nattr += 1
            get(R.attrs, (fn, nm), nothing) == attr_str(v) || badd("$fn attr $nm: module $(attr_str(v)) reader $(get(R.attrs, (fn, nm), nothing))")
        end
        id = ch.attrs["method"]
        for o in ch.orbits
            oid = string(o.id)
            for (nm, v) in o.attrs
                nattr += 1
                get(R.oattrs, (fn, oid, nm), nothing) == attr_str(v) || badd("$fn orbit $oid attr $nm")
            end
            sets = Pair{String, Any}["t" => o.t, "state" => o.state]
            o.winding === nothing || push!(sets, "winding" => o.winding)
            for (m, v) in o.lyap; push!(sets, "lyap/$m" => v); end
            for (m, v) in o.diagnostics; push!(sets, "diagnostics/$m" => v); end
            for (m, v) in o.cost; push!(sets, "cost/$m" => v); end
            for (iid, t) in o.indicators
                for f in (:series, :final, :flags, :scalars, :counts)
                    x = getfield(t, f)
                    x === nothing || push!(sets, "indicators/$iid/$f" => x)
                end
                nm = Dict(:series => (:columns, :columns), :final => (:fcolumns, :columns), :flags => (:bits, :bits), :scalars => (:scolumns, :columns))
                for (f, (field, attr)) in nm
                    getfield(t, f) === nothing && continue
                    ndata += 1
                    get(R.dattrs, (fn, oid, "indicators/$iid/$f", String(attr)), nothing) == attr_str(getfield(t, field)) || badd("$fn orbit $oid $f@$attr")
                end
                t.final !== nothing && (ndata += 1; get(R.dattrs, (fn, oid, "indicators/$iid/final", "rows"), nothing) == attr_str(t.rows) || badd("$fn orbit $oid final@rows"))
                t.range !== nothing && (ndata += 1; get(R.dattrs, (fn, oid, "indicators/$iid/scalars", "range"), nothing) == attr_str(t.range) || badd("$fn orbit $oid scalars@range"))
            end
            for (p, A) in sets
                ndata += 1
                d = get(R.data, (fn, oid, p), nothing)
                ok = d !== nothing && d.sha == rawsha(A) && d.eltype == string(eltype(A)) && d.jdims == join(size(A), "x")
                ok || badd("$fn orbit $oid dataset $p differs between readers")
            end
            for (m, thr) in o.thresholds
                ndata += 1
                string(thr) == get(R.dattrs, (fn, oid, "diagnostics/$m", "threshold"), "") || badd("$fn orbit $oid diagnostics/$m threshold")
            end
        end
    end
    return (files = nfile, attrs = nattr, datasets = ndata, bad = nbad, examples = bad)
end
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
    for (k, v) in Ra.pattrs
        get(Rb.pattrs, k, nothing) == v || note("provenance:" * k[3], "$(k[1]) provenance $(k[2]) $(k[3])")
    end
    for (k, v) in Ra.pooled
        get(Rb.pooled, k, nothing) == v || note("pooled", "$(k[1]) pooled/$(k[2])")
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
    for (S, T) in ((Ra.attrs, Rb.attrs), (Rb.attrs, Ra.attrs), (Ra.data, Rb.data), (Rb.data, Ra.data), (Ra.pattrs, Rb.pattrs), (Rb.pattrs, Ra.pattrs))
        for k in keys(S)
            haskey(T, k) || note("key-set", "missing key $(k)")
        end
    end
    return (diffs = diffs, examples = ex, total = sum(values(diffs); init = 0))
end
section("READER") do
    for (tag, store) in (("A", STORE_A), ("B", STORE_B), ("C", STORE_C), ("T4", STORE_T))
        rr = run_reader(store, tag)
        R = parse_reader(rr.path)
        CAMP["R$tag"] = R
        nfiles = length(chunk_files(store))
        nfail = count(l -> occursin("|FAIL|", l), R.checks)
        emit("O20/D86s", "independent reader (HDF5.jl in a separate process, no Chaoslab; checks the layout of R15, E24, E25, E31), store $tag: files read, checks failed, seconds", "$(length(R.files)) of $nfiles files, $nfail failed, $(round(rr.seconds; digits = 1)) s", "all files, 0 failed", pf(rr.ok && length(R.files) == nfiles && nfail == 0))
        foreach(l -> say("   reader: " * l), first([l for l in R.checks if occursin("|FAIL|", l)], 5))
        notes = [l for l in R.checks if occursin("|NOTE|", l)]
        isempty(notes) || info("O20", "store $tag reader notes: $(length(notes)); the first: $(first(notes))")
        cmp = compare_readers(store, R)
        emit("O20", "store $tag: module read_chunk equals the independent reader: files, attributes (root, orbit, provenance, dataset), datasets (sha-256 of the raw bytes, dtype, shape)", "$(cmp.files) files, $(cmp.attrs) attributes, $(cmp.datasets) datasets, $(cmp.bad) differ", "0 differ", pf(cmp.bad == 0 && cmp.files == nfiles))
        foreach(x -> say("   differ: " * x), cmp.examples)
    end
    R = CAMP["RB"]
    # dtype / shape / attribute-list table per kind (D86s)
    for kind in ["tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff"]
        fns = [fn for ((fn, nm), v) in R.attrs if nm == "kind" && v == kind]
        sets = Dict{String, Set{String}}()
        for ((fn, o, p), d) in R.data
            fn in fns || continue
            m = match(r"^(indicators/\w+/\w+|\w+)(/.*)?$", p)
            push!(get!(sets, replace(p, r"indicators/\w+/" => "indicators/<id>/", r"cost/\w+" => "cost/<id>"), Set{String}()), "$(d.eltype) $(d.jdims)")
        end
        info("D86s", "kind $kind: $(length(fns)) files; datasets (eltype and Julia shape; file extents are the reverse): " * join(("$p: " * join(sort!(collect(s))[1:min(end, 2)], " | ") * (length(s) > 2 ? " | ... ($(length(s)) shapes)" : "") for (p, s) in sort(collect(sets); by = first)), "; "))
    end
    emit("D86s", "`lyap` absent from the files of every new kind, present in the L1 files (the reader checks both)", "checked by the reader", "0 failed", pf(all(count(l -> occursin("|FAIL|", l), CAMP["R$t"].checks) == 0 for t in ("A", "B", "C", "T4"))))
    # I14: resumed against uninterrupted; and the execution-invariance comparisons
    df = diff_readers(CAMP["RA"], CAMP["RB"])
    emit("I14/D86t/D86y", "resumed (killed twice) campaign A against uninterrupted campaign B: differing datasets and attributes outside times and host fields (created_utc, cpu, threads, session), via the independent reader", df.total, "0", pf(df.total == 0))
    foreach(x -> say("   diff: " * x), df.examples)
    isempty(df.diffs) || info("I14", "difference classes: $(df.diffs)")
    dt = diff_readers(CAMP["RB"], CAMP["RT4"]; ignore = ("created_utc", "cpu", "threads", "session", "chunk_size", "launches_per_chunk"))
    emit("I24/D86q", "campaign B (23 threads, capacity 1000) against campaign T4 (4 threads, capacity 7): differing datasets and attributes outside times, host fields, chunk_size and launches_per_chunk", dt.total, "0", pf(dt.total == 0))
    foreach(x -> say("   diff: " * x), dt.examples)
    isempty(dt.diffs) || info("I24", "difference classes: $(dt.diffs)")
    dc = diff_readers(CAMP["RB"], CAMP["RC"]; ignore = ("created_utc", "cpu", "threads", "session", "chunk_size", "launches_per_chunk"))
    emit("I24/D86q", "campaign B (capacity 1000) against campaign C (own process, capacity from the free memory): differing datasets and attributes outside times, host fields, chunk_size and launches_per_chunk", dc.total, "0", pf(dc.total == 0))
    # h5dump (the HDF5 library's command-line tool, outside HDF5.jl and outside the module) on indicator datasets
    rng = MersenneTwister(8601)
    files = shuffle(rng, [f for f in chunk_files(STORE_B) if EX.read_chunk_attrs(joinpath(STORE_B, f))["status"] == "ok" && haskey(EX.read_chunk_attrs(joinpath(STORE_B, f)), "kind")])[1:24]
    nbad = ndone = 0
    HDF5_jll.h5dump() do exe
        for fn in files
            ch = EX.read_chunk(joinpath(STORE_B, fn))
            o = ch.orbits[min(2, end)]
            id = ch.attrs["method"]
            t = o.indicators[id]
            for (path, A) in [("/orbits/$(o.id)/indicators/$id/$f", getfield(t, f)) for f in (:series, :final, :flags, :scalars, :counts) if getfield(t, f) !== nothing]
                txt = read(`$exe -d $path -m %.17g -y -w 0 $(joinpath(STORE_B, fn))`, String)
                body = match(r"DATA \{([^}]*)\}", txt)
                body === nothing && (nbad += 1; continue)
                toks = filter(!isempty, split(String(body.captures[1]), r"[,\s]+"))
                vals = eltype(A) <: Integer ? parse.(Int64, toks) : parse.(Float64, toks)
                ndone += 1
                (length(vals) == length(A) && all(i -> isequal(vals[i], Float64(vec(A)[i])) || (isnan(vals[i]) && isnan(vec(A)[i])), eachindex(vals))) || (nbad += 1)
            end
        end
    end
    emit("O20", "HDF5 command-line tool h5dump (-m %.17g; no HDF5.jl, no Chaoslab) on $(length(files)) ok files of the six kinds, every indicator dataset of one orbit each (series, final, flags, scalars, counts): equal to the module's values exactly", "$ndone datasets, $nbad differ", "0 differ", pf(nbad == 0 && ndone > 0))
end

# ---------------------------------------------------------------------------------------------------
# R7 and R19: the manifest
const R7_KEYS = ["study", "schema_version", "plan_sha256", "plan_counts", "status", "reasons", "sessions", "git", "versions", "host", "threads", "backend_table",
                 "device", "device_record", "prior_record", "prereg", "panel", "block", "defaults", "launch", "chunks", "arms", "wall", "overrides", "removed_temporaries"]
const R7_CHUNK = ["ordinal", "key", "hash", "arms", "status", "attempt", "backend", "fallback_reason", "orbits", "steps", "wall_seconds", "first_attempt_reasons", "flags", "setup_check"]
const R7_ARM = ["steps", "wall_seconds", "thread_seconds", "cost_ns_per_orbit_step", "budget_ns", "status"]
const R19_SPEC = ["kind", "spec", "spec_sha256", "chain", "budget", "bytes_per_orbit"]
section("MANIFEST") do
    for (tag, mp, store) in (("B", MAN_B, STORE_B), ("A", MAN_A, STORE_A))
        m = EX.read_json(mp)
        miss = [k for k in R7_KEYS if !haskey(m, k)]
        chunks = m["chunks"]
        missc = unique([k for c in chunks for k in R7_CHUNK if !haskey(c, k)])
        missa = unique([k for (_, a) in m["arms"] for k in R7_ARM if !haskey(a, k)])
        sts = Dict{String, Int}(); for c in chunks; sts[c["status"]] = get(sts, c["status"], 0) + 1; end
        ok7 = isempty(miss) && isempty(missc) && isempty(missa) && length(chunks) == NCH && [c["ordinal"] for c in chunks] == collect(1:NCH) && m["plan_sha256"] == PLAN.hash && m["schema_version"] == 2
        emit("R7", "manifest $tag: R7 keys, chunk entries ($(length(chunks)) of $NCH, statuses $sts), arm entries, plan hash, schema version", "missing $miss / $missc / $missa", "none missing", pf(ok7))
        ind = get(m, "indicators", nothing)
        emit("R19", "manifest $tag: `indicators` table by spec id", ind === nothing ? "absent" : join(sort(collect(keys(ind))), ","), "all six specs", pf(ind !== nothing && sort(collect(keys(ind))) == sort(collect(keys(SPECS)))))
        ind === nothing && continue
        bad = String[]
        for (id, sp) in SPECS
            e = ind[id]
            miss = [k for k in R19_SPEC if !haskey(e, k)]
            isempty(miss) || push!(bad, "$id missing $miss")
            e["kind"] == sp["kind"] || push!(bad, "$id kind")
            e["spec_sha256"] == EX.spec_hash(sp) || push!(bad, "$id spec_sha256")
            e["spec"] == EX.spec_text(sp) || push!(bad, "$id spec text")
            sp["kind"] in ("frequency", "phase_advance", "birkhoff") && (haskey(e, "windows") && haskey(e, "run_length") || push!(bad, "$id windows/run_length"))
            sp["kind"] == "pair" && (haskey(e, "frame_source") && haskey(e, "direction_rule") || push!(bad, "$id frame_source/direction_rule"))
            sp["kind"] == "birkhoff" && (haskey(e, "observable_names") || push!(bad, "$id observable_names"))
            b = e["bytes_per_orbit"]
            b["total"] == b["M"] + b["R"] + b["S"] + b["B"] || push!(bad, "$id bytes terms")
        end
        emit("R19", "manifest $tag: per spec kind, canonical text, spec hash (recomputed from the block), chain record, budget, bytes(orbit) terms (sum equals total), windows and run length, frame source and direction rule, observable names", "$(length(bad)) defects $(first(bad, 3))", "0", pf(isempty(bad)))
        tg = ind["tg"]
        info("R19", "manifest $tag: chain record of tg $(tg["chain"]); bytes(orbit) of tg $(tg["bytes_per_orbit"])")
        specrows = [r for r in m["backend_table"] if get(r, "method", get(r, "spec", "")) in keys(SPECS)]
        emit("R19", "manifest $tag: backend table has a row per (precision, spec id) for the six specs", "$(length(specrows)) of $(length(m["backend_table"])) rows are spec rows", "12 rows (6 specs x 2 precisions, Float64 and Double64)", pf(length(specrows) >= 12))
        # per chunk entries: flags recount, not_applicable, transform seconds
        R = load_store(store)
        nflag = nbadflag = nna = ntr = 0
        for c in chunks
            plan_c = PLAN.chunks[c["ordinal"]]
            id = plan_c.run_type.method
            haskey(SPECS, id) || continue
            c["status"] == "failed" && continue
            nflag += 1
            haskey(c, "indicator_flags") || (nbadflag += 1; continue)
            rec = Dict{String, Int}()
            for o in R[c["ordinal"]].orbits
                t = o.indicators[id]
                t.flags === nothing && continue
                for (i, nm) in enumerate(t.bits)
                    any(w -> (w >> (i - 1)) & 1 == 1, t.flags) && (rec[nm] = get(rec, nm, 0) + 1)
                end
            end
            Dict{String, Int}(k => Int(v) for (k, v) in c["indicator_flags"]) == rec || (nbadflag += 1)
            SPECS[id]["kind"] == "frequency" && (haskey(c, "transform_seconds") ? ntr += 1 : (nbadflag += 1))
            if plan_c.cell == "controls" && id in ("fq", "bk")
                get(get(c, "not_applicable", Dict()), id, nothing) == [883] ? (nna += 1) : (nbadflag += 1)
            end
        end
        emit("R19", "manifest $tag: `indicator_flags` of $nflag indicator chunks equal the recount from the files; `not_applicable` of the controls chunks of fq and bk lists orbit 883 ($nna); `transform_seconds` on $ntr frequency chunks", "$nbadflag defects", "0", pf(nbadflag == 0 && nna == 2 && ntr > 0))
        flagtot = Dict{String, Int}()
        for c in chunks, (k, v) in get(c, "indicator_flags", Dict()); flagtot[k] = get(flagtot, k, 0) + v; end
        info("R19", "manifest $tag: flag counts summed over the chunks: $(sort(collect(flagtot)))")
        arms_ok = all(haskey(a, "specs") for (n, a) in m["arms"] if haskey(first(a for a in PLAN.block.raw["arm"] if a["name"] == n), "indicators"))
        emit("R19", "manifest $tag: per-arm totals have the per-spec measured cost per orbit-step against the spec's budget line", arms_ok, "true", pf(arms_ok))
        for (n, a) in sort(collect(m["arms"]); by = first), (id, c) in sort(collect(get(a, "specs", Dict())); by = first)
            say("COST | arm $n | spec $id | measured $(fmt(c["measured_ns_per_orbit_step"])) ns per orbit-step (thread time) | budget line $(fmt(c["budget_line_ns"])) ns")
        end
        tag == "B" && info("R19", "manifest B: thread count $(m["threads"]), device $(m["device"]), launch $(m["launch"])")
    end
end

# ---------------------------------------------------------------------------------------------------
# the files, orbit by orbit: D86e, D86g, D86h, D86i, D86j, D86k, D86l, D86m (V on real panel orbits), I20, I21, I22, O16, O17
section("FILES") do
    R = load_store(STORE_B)
    chunks_by_id = Dict{String, Vector{Int}}()
    for (ord, ch) in R
        push!(get!(chunks_by_id, ch.attrs["method"], Int[]), ord)
    end
    # ---- D86e: chain records, forms, member lists, spec text, provenance per member
    nbad = nfiles = 0
    expected_forms(members) = ["$m(replace=$(i == length(members)))" for (i, m) in enumerate(members)]
    for id in ("tg", "pr", "st"), ord in chunks_by_id[id]
        ch = R[ord]; a = ch.attrs
        nfiles += 1
        sp = SPECS[id]
        want = sp["kind"] == "stretching" ? ["stretching(replace=true)"] : expected_forms(sp["members"])
        ok = a["chain"] == want && a["spec"] == EX.spec_text(sp) && a["spec_sha256"] == EX.spec_hash(sp) && haskey(a, "provenance") &&
             (sp["kind"] == "stretching" || (sort(collect(keys(a["provenance"]))) == sort(sp["members"]) && length(unique(a["provenance"][m]["indicator"] for m in sp["members"])) == length(sp["members"])))
        ok || (nbad += 1)
    end
    emit("D86e", "chain record equals the forms the module must set (all members but the last non-replacing, the last replacing), the spec text and hash equal the block's, one provenance group per member with its own `indicator` value, on $nfiles files of tg, pr, st", nbad, "0", pf(nbad == 0))
    # ---- D86g/I20: t equals the grid; finiteness; FLI nondecreasing; SALI range
    nt = nfin = nfli = nsali = ntot = 0
    eps64 = eps(Float64)
    for ord in chunks_by_id["tg"], o in R[ord].orbits
        R[ord].attrs["status"] == "ok" || continue
        it = ITEM_OF[o.attrs["item_hash"]]
        t = o.indicators["tg"]
        ntot += 1
        nt += o.t == grid_of(it) ? 0 : 1
        nfin += all(isfinite, t.series) ? 0 : 1
        f = colvec(t, "FLI.fli"); f1 = colvec(t, "FLI.fli1")
        nfli += (issorted(f) && issorted(f1)) ? 0 : 1
        s = colvec(t, "SALI.sali")
        nsali += all(x -> 0 <= x <= sqrt(2) * (1 + 2 * eps64), s) ? 0 : 1
    end
    emit("D86g", "tangent family on $ntot ok orbits (Float64, Double64; wrapped, lifted): t differs from the grid / non-finite column / FLI decreasing / SALI outside [0, sqrt 2 (1 + 2 eps)]", "$nt / $nfin / $nfli / $nsali", "0 / 0 / 0 / 0", pf(nt + nfin + nfli + nsali == 0))
    # ---- I21: the orbit point and the monitor columns equal those of the L1 item of the same orbit (grid of horizon 4000); the window kinds at the common steps
    l1 = Dict{Int, Any}()
    for ord in chunks_by_id["L1"], o in R[ord].orbits
        l1[o.id] = o
    end
    nchk = ndiff = 0
    for id in ("tg", "st"), ord in chunks_by_id[id]
        R[ord].attrs["status"] == "ok" && R[ord].attrs["coordinate_mode"] == "wrapped" && R[ord].attrs["precision"] == "Float64" || continue
        for o in R[ord].orbits
            ref = l1[o.id]
            nchk += 1
            same = isbitseq(o.t, ref.t) && isbitseq(o.state, ref.state) && o.winding !== nothing && isbitseq(o.winding, ref.winding) &&
                   all(m -> isbitseq(o.diagnostics[m], ref.diagnostics[m]), ["det", "symp"])
            ndiff += same ? 0 : 1
        end
    end
    emit("I21", "orbit point (state, winding) and the two monitor columns at the 135 grid steps of tg and st (Float64, wrapped) against the L1 item of the same orbit: orbits compared / differing", "$nchk / $ndiff", "0 differ", pf(ndiff == 0 && nchk > 600))
    nchk = ndiff = 0
    for id in ("fq", "bk", "pa", "pr"), ord in chunks_by_id[id]
        R[ord].attrs["status"] == "ok" && R[ord].attrs["coordinate_mode"] == "wrapped" && R[ord].attrs["precision"] == "Float64" || continue
        for o in R[ord].orbits
            ref = l1[o.id]
            common = intersect(o.t, ref.t)
            ia, ib = [findfirst(==(n), o.t) for n in common], [findfirst(==(n), ref.t) for n in common]
            nchk += 1
            same = isbitseq(o.state[ia, :], ref.state[ib, :]) && isbitseq(o.winding[ia, :], ref.winding[ib, :])
            ndiff += same ? 0 : 1
        end
    end
    emit("I21", "orbit point of fq, bk, pa, pr (Float64, wrapped; windows 2000 and 2001, pair 300 steps) against the L1 item at the grid steps both have: orbits compared / differing", "$nchk / $ndiff", "0 differ", pf(ndiff == 0 && nchk > 1000))
    # ---- sample of orbits for the stand-alone comparisons: 4 per K index (Float64 wrapped), 3 Double64 ones
    rng = MersenneTwister(8602)
    byk = Dict{Int, Vector{Int}}()
    for o in PLAN.items
        o.method == "tg" && o.precision == "Float64" && o.coordinate_mode == "wrapped" && o.design_class != "control" && push!(get!(byk, parse(Int, match(r"^k(\d+):", o.cell).captures[1]), Int[]), o.orbit)
    end
    picks = vcat([shuffle(rng, unique(v))[1:4] for (k, v) in sort(collect(byk); by = first)]...)
    ctrl = [871, 872, 873, 883]
    dbl = [it.orbit for it in PLAN.items if it.method == "tg" && it.precision == "Double64"][1:3:9]
    info("O16", "sample orbits (Float64 wrapped): $picks; controls $ctrl; Double64: $dbl")
    item_for(id, o, prec, mode) = first(it for it in PLAN.items if it.method == id && it.orbit == o && it.precision == prec && it.coordinate_mode == mode)
    function stored(id, o, prec, mode)
        it = item_for(id, o, prec, mode)
        ord = first(c.ordinal for c in PLAN.chunks if it.id in c.items)
        return it, first(x for x in R[ord].orbits if x.id == o), R[ord].attrs["status"]
    end
    # ---- O16, I20, I22: tangent family and pair against the stand-alone replacing reducers at 5 seeded grid steps and the last
    ncomp = nbad = 0
    bad = String[]
    for (o, prec) in vcat([(o, "Float64") for o in vcat(picks, ctrl)], [(o, "Double64") for o in dbl])
        it, rec, st = stored("tg", o, prec, "wrapped")
        st == "ok" || continue
        rows = vcat(sort(shuffle(rng, 1:length(rec.t)-1)[1:5]), length(rec.t))
        for r in rows, m in SPECS["tg"]["members"]
            ref = stand_alone(it, m, rec.t[r])
            names = Dict("SALI" => ["SALI.sali", "SALI.L"], "FLI" => ["FLI.fli", "FLI.fli1"], "MEGNO" => ["MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"])[m]
            got = Float64[colvec(rec.indicators["tg"], nm)[r] for nm in names]
            ncomp += length(names)
            ok = isbitseq(got, ref)
            ok || (nbad += 1; length(bad) < 4 && push!(bad, "orbit $o $prec $m step $(rec.t[r]): stored $got stand-alone $ref"))
        end
    end
    emit("O16/I20/I22", "tg (chain FLI, MEGNO, SALI) against the stand-alone replacing FLI(), MEGNO(), SALI() run directly to the step, on the identical orbit, at 5 seeded grid steps and the last: orbits $(length(picks) + length(ctrl) + length(dbl)), values compared / differing", "$ncomp / $nbad", "0 differ (bitwise)", pf(nbad == 0 && ncomp > 100))
    foreach(x -> say("   differ: " * x), bad)
    # pair
    ncomp = nbad = 0; bad = String[]; ndx = 0; nxi = 0
    for (o, prec) in vcat([(o, "Float64") for o in vcat(picks[1:2:end], [871, 883])], [(o, "Double64") for o in dbl[1:2]])
        it, rec, st = stored("pr", o, prec, "wrapped")
        st == "ok" || continue
        sp = SPECS["pr"]
        xi0 = Float64.(rec.attrs["xi0"])
        rows = vcat(sort(shuffle(rng, 1:length(rec.t)-1)[1:5]), length(rec.t))
        for r in rows, m in sp["members"]
            ref = pair_ref(it, sp, xi0, rec.t[r], m)
            names = m == "RLI" ? ["RLI.DL", "RLI.RLI"] : ["MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"]
            got = Float64[colvec(rec.indicators["pr"], nm)[r] for nm in names]
            ncomp += length(names)
            isbitseq(got, ref.cols) || (nbad += 1; length(bad) < 4 && push!(bad, "pair orbit $o $prec $m step $(rec.t[r])"))
        end
        ps = pair_ref(it, sp, xi0, 1, "RLI").start
        ndx += (Float64(ps.Dx_r) == rec.attrs["Dx_r"] || (hasproperty(ps, :Dx_r) && isequal(Float64(ps.Dx_r), rec.attrs["Dx_r"]))) ? 0 : 1
        nxi += abs(sum(abs2, xi0) - 1) < 1e-14 ? 0 : 1
    end
    emit("O16/D86h", "pr (RLI, MEGNO on the pair map) against the stand-alone RLI(), MEGNO() on the PairMap with the stored xi_0 at 6 grid steps per orbit (horizon 300): values compared / differing; stored Dx_r differing from the realized separation of pair_start; xi_0 not a unit vector", "$ncomp / $nbad; $ndx; $nxi", "0; 0; 0", pf(nbad == 0 && ndx == 0 && nxi == 0 && ncomp > 100))
    foreach(x -> say("   differ: " * x), bad)
    let
        dxr = Float64[]; nout = 0
        for ord in chunks_by_id["pr"], o in R[ord].orbits
            R[ord].attrs["status"] == "ok" || continue
            push!(dxr, o.attrs["Dx_r"] / o.attrs["Dx"])
            nout += (0.5 <= o.attrs["Dx_r"] / o.attrs["Dx"] <= 2) ? 0 : 1
        end
        emit("D86h/R13", "realized separation Dx_r / nominal Dx of all $(length(dxr)) pair orbits: min, max, outside [1/2, 2]", "$(fmt(minimum(dxr))), $(fmt(maximum(dxr))), $nout", "0 outside", pf(nout == 0))
    end
    # RLI of the Y5 control is zero at every grid step (rli RL9a)
    let
        ok = true; nrow = 0
        for ord in chunks_by_id["pr"], o in R[ord].orbits
            o.id == 883 || continue
            t = o.indicators["pr"]
            ok &= all(iszero, colvec(t, "RLI.RLI")) && all(iszero, colvec(t, "RLI.DL")); nrow += length(o.t)
        end
        emit("O18/D86h", "RLI and DL exactly 0 at every grid step for the Y5 control (orbit 883), Float64 wrapped and the Double64 chunk if any ($nrow rows)", ok, "exactly 0", pf(ok && nrow > 0))
    end
    # ---- D86i stretching
    ncomp = nbad = 0; bad = String[]; ncounts = ncb = 0
    ST = SPECS["st"]
    for (o, prec) in vcat([(o, "Float64") for o in vcat(picks, ctrl)], [(o, "Double64") for o in dbl])
        it, rec, st = stored("st", o, prec, "wrapped")
        st == "ok" || continue
        t = rec.indicators["st"]
        rows = vcat(sort(shuffle(rng, 1:length(rec.t)-1)[1:3]), length(rec.t))
        for r in rows
            sf = stretch_final(it, ST, rec.t[r])
            ncomp += 1
            Float64(sf.fin.D) == t.series[r, 1] || (nbad += 1; length(bad) < 4 && push!(bad, "orbit $o $prec step $(rec.t[r]) D stored $(t.series[r, 1]) stand-alone $(Float64(sf.fin.D))"))
        end
        sf = stretch_final(it, ST, rec.t[end])
        f(v) = Float64.(collect(v))
        want = Float64[sf.fin.handoffs; f(sf.fin.sum); f(sf.fin.comp); f(sf.fin.sum2); f(sf.fin.comp2); f(sf.fin.mean); f(sf.fin.std); f(sf.fin.below); f(sf.fin.above)]
        ncomp += 1
        isbitseq(t.scalars[1:17], want) || (nbad += 1; length(bad) < 4 && push!(bad, "orbit $o $prec scalars"))
        isbitseq(t.range, Float64[sf.lo, sf.hi]) || (nbad += 1)
        if t.counts !== nothing
            ncounts += 1
            ncb += (t.counts == sf.fin.counts && all(c -> sum(Int.(t.counts[:, c])) == Int(t.scalars[1]), 1:2)) ? 0 : 1
        end
    end
    emit("D86i/O16", "stretching: stored D at 3 seeded grid steps and the last, the 17 orbit scalars and the range against the stand-alone final that includes the workspace (bitwise); values compared / differing", "$ncomp / $nbad", "0 differ", pf(nbad == 0 && ncomp > 50))
    foreach(x -> say("   differ: " * x), bad)
    nin = nout = nbadc = 0
    for ord in chunks_by_id["st"], o in R[ord].orbits
        R[ord].attrs["status"] == "ok" || continue
        in_sl = ITEM_OF[o.attrs["item_hash"]].in_SL
        (o.indicators["st"].counts !== nothing) == in_sl ? (nin += 1) : (nbadc += 1)
    end
    emit("D86i", "stretching counts are stored for the orbits of the subset S_L only (a count matrix present exactly when the orbit is in S_L): orbits agreeing / disagreeing", "$nin / $nbadc", "0 disagree", pf(nbadc == 0))
    emit("D86i", "stored count matrices of the S_L orbits checked against the stand-alone counts ($ncounts), disagreeing", ncb, "0", pf(ncb == 0 && ncounts > 0))
    # pooled counts: sum over all orbits of the chunk of the stand-alone counts, for 3 chunks
    let
        npool = nbadp = 0
        stch = [ord for ord in chunks_by_id["st"] if R[ord].attrs["status"] == "ok" && R[ord].attrs["precision"] == "Float64"]
        for ord in stch[[1, div(length(stch), 2), length(stch)]]
            tot = zeros(Int64, ST["bins"] + 2, 2)
            for o in R[ord].orbits
                it = ITEM_OF[o.attrs["item_hash"]]
                tot .+= Int64.(stretch_final(it, ST, it.horizon).fin.counts)
            end
            fn = chunkname(ord)
            pooled = h5open(joinpath(STORE_B, fn), "r") do f
                (read(f["pooled/st/counts"]), read_attribute(f["pooled/st"], "orbits"))
            end
            npool += 1
            (pooled[1] == tot && pooled[2] == length(R[ord].orbits)) || (nbadp += 1)
        end
        emit("D86i/E29", "pooled counts of 3 chunks (first, middle, last of the stretching chunks) against the exact integer sum of the stand-alone counts of every orbit of the chunk; chunks checked / differing", "$npool / $nbadp", "0 differ", pf(nbadp == 0 && npool == 3))
    end
    # ---- O17: the window kinds from series built by plain iteration of the map (the Chirikov map written here from the equations, no Chaoslab code in the iteration)
    nfq = nbfq = npa = nbpa = nbk = nbbk = 0; bad = String[]
    FQ, PA, BK = SPECS["fq"], SPECS["pa"], SPECS["bk"]
    wsel = vcat(picks[1:2:end], [871, 872, 874])
    for o in wsel
        it, rec, st = stored("fq", o, "Float64", "wrapped")
        st == "ok" || continue
        ser = plain_series(it, 2 * maximum(FQ["windows"]))
        z = ComplexF64[cis(p.q) for p in ser]
        t = rec.indicators["fq"]
        for (h, W) in enumerate(FQ["windows"])
            r = Chaoslab.frequency_map(z[1:2W], W; order = FQ["order"], positions = Symbol(FQ["positions"]), pad = FQ["pad"], mean_removal = FQ["mean_removal"], mirror_ratio = FQ["mirror_ratio"], floor_multiplier = FQ["floor_multiplier"])
            w1, w2 = r.w1[1], r.w2[1]
            vals = Float64[w1.nu, w2.nu, w1.amplitude, w2.amplitude, w1.mirror_amplitude, w2.mirror_amplitude, w1.second_ratio, w2.second_ratio, r.dist[1], r.delta, r.D]
            nfq += 1
            isbitseq(t.final[h, :], vals) || (nbfq += 1; length(bad) < 4 && push!(bad, "fq orbit $o W=$W"))
        end
        it, rec, st = stored("pa", o, "Float64", "wrapped")
        st == "ok" || continue
        Wm = maximum(PA["windows"])
        ser = plain_series(it, 2Wm + 1)
        zz = ComplexF64[Complex(p.q, p.p) for p in ser]
        inc = Chaoslab.phase_increments(zz)
        t = rec.indicators["pa"]
        for (h, W) in enumerate(PA["windows"])
            nu = Chaoslab.weighted_phase_frequency(inc, W; p = PA["birkhoff_order"])
            delta = Chaoslab.circle_distance(nu[1], nu[2])
            fi = Chaoslab.frequency_index(delta)
            npa += 1
            isbitseq(t.final[h, :], Float64[nu[1], nu[2], delta, fi.D]) || (nbpa += 1; length(bad) < 4 && push!(bad, "pa orbit $o W=$W"))
        end
        it, rec, st = stored("bk", o, "Float64", "wrapped")
        st == "ok" || continue
        Wm = maximum(BK["windows"])
        ser = plain_series(it, 2Wm)
        obs = map(BK["observables"]) do ob
            c = ob["coordinate"]
            ob["kind"] == "cos" ? Chaoslab.CosObservable(c, ob["frequency"]) : ob["kind"] == "sin" ? Chaoslab.SinObservable(c, ob["frequency"]) : Chaoslab.LiftedObservable(c, c + 2)
        end
        S = Matrix{Float64}(undef, 2Wm, length(obs))
        for (j, p) in enumerate(ser), (k, f) in enumerate(obs)
            S[j, k] = Float64(f(SVector(p.q, p.p, p.m, p.n)))
        end
        r = Chaoslab.weighted_birkhoff(S, Tuple(BK["windows"]); order = BK["birkhoff_order"])
        t = rec.indicators["bk"]
        for (h, W) in enumerate(BK["windows"]), k in eachindex(obs)
            ix = Chaoslab.birkhoff_index(r.avg1[k, h], r.avg2[k, h], r.S[k, h])
            vals = Float64[r.avg1[k, h], r.avg2[k, h], r.S[k, h], r.B[h], ix.delta, ix.absdig, ix.reldig, ix.dig]
            nbk += 1
            isbitseq(t.final[(h - 1) * length(obs) + k, :], vals) || (nbbk += 1; length(bad) < 4 && push!(bad, "bk orbit $o W=$W obs $k"))
        end
    end
    emit("O17/D86j(d)", "frequency tables (2 windows) from the plain-iteration series of the Float64 wrapped map and frequency_map: rows compared / differing", "$nfq / $nbfq", "0 differ (bitwise)", pf(nbfq == 0 && nfq > 10))
    emit("O17/D86l", "phase-advance tables from the plain-iteration series, phase_increments and weighted_phase_frequency: rows compared / differing", "$npa / $nbpa", "0 differ (bitwise)", pf(nbpa == 0 && npa > 10))
    emit("O17/D86m", "weighted Birkhoff tables (4 observables, 2 windows) from the plain-iteration series and weighted_birkhoff: rows compared / differing", "$nbk / $nbbk", "0 differ (bitwise)", pf(nbbk == 0 && nbk > 40))
    foreach(x -> say("   differ: " * x), bad)
    # Double64 window tables: the package's forward as the series (the plain Float64 iteration does not apply)
    nd = nbd = 0
    for o in [it.orbit for it in PLAN.items if it.method == "fq" && it.precision == "Double64"][1:4]
        it, rec, st = stored("fq", o, "Double64", "wrapped")
        st == "ok" || continue
        z = with_prec("Double64") do T
            sys, x = sys_x(it, T)
            out = Complex{T}[]
            for _ in 1:(2 * maximum(FQ["windows"]))
                x = Chaoslab.forward(sys, x); push!(out, cis(x[1]))
            end
            out
        end
        for (h, W) in enumerate(FQ["windows"])
            r = Chaoslab.frequency_map(z[1:2W], W; order = FQ["order"], positions = Symbol(FQ["positions"]), pad = FQ["pad"], mean_removal = FQ["mean_removal"], mirror_ratio = FQ["mirror_ratio"], floor_multiplier = FQ["floor_multiplier"])
            nd += 1
            isbitseq(rec.indicators["fq"].final[h, [1, 2, 9, 10, 11]], Float64[r.w1[1].nu, r.w2[1].nu, r.dist[1], r.delta, r.D]) || (nbd += 1)
        end
    end
    emit("O17/D86j(d)", "Double64 frequency rows (nu1, nu2, dist, delta, D) from the series of the package's forward in Double64: rows compared / differing", "$nd / $nbd", "0 differ", pf(nbd == 0 && nd > 0))
    # ---- D86j(b),(f): bytes per orbit and window; no dataset holds the series (a frequency file's dataset lengths are the grid length, the window count or the layout's)
    let
        lay = EX.chunk_layout(joinpath(STORE_B, chunkname(first(chunks_by_id["fq"]))))
        paths = sort(collect(keys(lay)))
        oid = minimum(parse(Int, match(r"^/orbits/(\d+)/", p).captures[1]) for p in paths if occursin(r"^/orbits/\d+/", p))
        orbit1 = [p for p in paths if startswith(p, "/orbits/$oid/")]
        sizes = [prod(lay[p].size) for p in orbit1]
        grid_n = length(grid_of(first(it for it in PLAN.items if it.method == "fq" && it.precision == "Float64")))
        biggest = maximum(sizes)
        emit("D86j(f)", "frequency file: largest dataset of an orbit has $biggest elements (grid length $grid_n, 2 windows x 11 columns = 22): no dataset holds the 2000-sample series", biggest, "<= 2 x grid length (270)", pf(biggest <= 2 * grid_n))
        fq = lay[first(p for p in paths if occursin("indicators/fq/final", p))]
        fl = lay[first(p for p in paths if occursin("indicators/fq/flags", p))]
        bytes_per_window = (prod(fq.size) + prod(fl.size)) * 8 / length(FQ["windows"])
        emit("D86j(b)", "stored bytes per orbit and window for frequency: final ($(fq.eltype) $(fq.size)) + flag word ($(fl.eltype) $(fl.size))", bytes_per_window, "96 (eleven binary64 and one int64 flag word)", pf(bytes_per_window == 96))
    end
    # ---- D86k/flags: counts of flag words over the campaign (degenerate on the fixed points etc.); the module draws no conclusion
    fl = Dict{String, Int}()
    for id in ("fq", "pa", "bk"), ord in chunks_by_id[id], o in R[ord].orbits
        t = o.indicators[id]
        for (i, nm) in enumerate(t.bits), w in t.flags
            (w >> (i - 1)) & 1 == 1 && (fl["$id.$nm"] = get(fl, "$id.$nm", 0) + 1)
        end
    end
    info("D86k", "flag rows set over the campaign (spec.flag => rows): $(sort(collect(fl)))")
    nfp = [(o.id, o.indicators["fq"].flags) for ord in chunks_by_id["fq"] for o in R[ord].orbits if o.id in (871, 872, 873, 875) && R[ord].attrs["precision"] == "Float64"]
    emit("D86k", "fixed-point controls (hyperbolic and elliptic) in the frequency chunk: chunk status ok and the `degenerate` bit set where the series is constant", "$(length(nfp)) orbits; flags $(first(nfp, 3))", "status ok, flags data", pf(all(R[ord].attrs["status"] == "ok" for ord in chunks_by_id["fq"] if PLAN.chunks[ord].cell == "controls")))
end

# ---------------------------------------------------------------------------------------------------
# O18: the closed forms through the recorder (K_std = 0 shear with direction 0, the Y5 control, the K_std = 0 line for the frequency kinds)
section("O18") do
    include(joinpath(@__DIR__, "experiments-indicators-closed.jl"))
end

# ---------------------------------------------------------------------------------------------------
# O19: independent hashing with the platform's command-line SHA-256
section("O19") do
    include(joinpath(@__DIR__, "experiments-indicators-hash.jl"))
end

# ---------------------------------------------------------------------------------------------------
# D86v: CPU rerun of 20 items per arm; the pilot at the own run length of the window kinds
section("RERUN") do
    t = @elapsed r = EX.rerun_check(PP.prereg, PP.panel, STORE_B; threads = NT, thread_override = NT != K19, allow_dirty = true)
    inj = Set(PLAN.chunks[ORD_FAIL].items)
    nb = count(x -> x.bitwise, r)
    arms = unique(x.arm for x in r)
    notb = [x for x in r if !x.bitwise]
    notb_inj = [x for x in notb if x.item in inj]
    emit("D86v/D86y", "CPU reference rerun of rerun_sample = $(PLAN.block.raw["verify"]["rerun_sample"]) items per arm on store B ($(length(arms)) arms, $(round(t; digits = 1)) s): items, bitwise equal, backend of the stored items", "$(length(r)), $nb, $(join(unique(x.backend for x in r), ","))", "all bitwise (CPU)", pf(length(notb) == length(notb_inj)))
    isempty(notb) || info("D86v", "items not bitwise: $(length(notb)), of which in the chunk of the injected failure: $(length(notb_inj)) (a failed chunk stores blank rows from the failure on, so a rerun cannot equal it)")
    kinds = Dict{String, Int}(); for x in r; id = PLAN.items[x.item].method; kinds[id] = get(kinds, id, 0) + 1; end
    info("D86v", "rerun items by spec: $(sort(collect(kinds)))")
end
section("PILOT") do
    out = joinpath(SCR, "pilot")
    rm(out; recursive = true, force = true)
    t = @elapsed res = EX.run_pilot(PP.prereg, PP.panel; out = out, threads = NT, thread_override = NT != K19, allow_dirty = true)
    rep = res.report
    rows = [r for r in rep["cost"] if get(r, "method", "") in keys(SPECS)]
    info("D86v", "pilot report keys: $(sort(collect(keys(rep)))); written $(res.path) in $(round(t; digits = 1)) s")
    txt = read(res.path, String)
    emit("D86v", "pilot report states that the numbers are not study data", occursin("not study data", lowercase(txt)) || occursin("not study", lowercase(txt)), "true", pf(occursin("not study data", lowercase(txt)) || occursin("not study", lowercase(txt))))
    for r in rows
        say("PILOT | spec $(r["method"]) | " * join(("$k=$(fmt(v))" for (k, v) in sort(collect(r); by = first) if k in ("precision", "horizon", "bytes_per_orbit", "capacity", "orbit_pass_seconds", "transform_seconds", "thread_seconds_per_orbit_step", "wall_seconds_per_orbit_step", "flag_shares")), " | "))
    end
    emit("D86v", "pilot rows for the six specs: window kinds at their own run length (fq, bk 2000; pa 2001), series kinds at the pilot horizon (300)", join(sort(unique("$(r["method"]):$(get(r, "horizon", "-"))" for r in rows)), ", "), "fq/bk 2000, pa 2001, others 300; fq reports the orbit pass and the transform separately",
         pf(!isempty(rows) && all(r -> get(r, "horizon", 0) == (r["method"] in ("fq", "bk") ? 2000 : r["method"] == "pa" ? 2001 : 300), rows) && any(r -> haskey(r, "transform_seconds") && haskey(r, "orbit_pass_seconds"), rows)))
    emit("D86y", "study data folder after the pilot and the campaigns: files under data/p1-le-benchmark and entries under data/ unchanged", "$(filecount(STUDY_DATA)) files, $(length(tree_listing(joinpath(ROOT, "data")))) entries", "$(DATA_BEFORE[1]) files, $(length(DATA_BEFORE[2])) entries", pf(filecount(STUDY_DATA) == DATA_BEFORE[1] && tree_listing(joinpath(ROOT, "data")) == DATA_BEFORE[2]))
end

# ---------------------------------------------------------------------------------------------------
# D86n: bytes(orbit) and the launch capacity, with the worked numbers of R17 item 4
section("MEM") do
    include(joinpath(@__DIR__, "experiments-indicators-memory.jl"))
end

# ---------------------------------------------------------------------------------------------------
# D86w (V): the module holds no study value, no threshold, no detection
section("D86W") do
    include(joinpath(@__DIR__, "experiments-indicators-source.jl"))
end

# ---------------------------------------------------------------------------------------------------
section("ISOLATION") do
    now = (filecount(STUDY_DATA), tree_listing(joinpath(ROOT, "data")))
    emit("D86y", "the study's data folder: files under data/p1-le-benchmark and every entry (path, size, modification time) under data/ before and after this script", "$(now[1]) files / $(length(now[2])) entries; listings equal: $(now[2] == DATA_BEFORE[2])", "unchanged (0 files written by this check; the P1 campaign's own $(DATA_BEFORE[1]) files are untouched)", pf(now[1] == DATA_BEFORE[1] && now[2] == DATA_BEFORE[2]))
    gp = joinpath(ROOT, "data", "gpu-benchmark.toml")
    info("D86y", "data/gpu-benchmark.toml untouched (the record is not retaken for D86): size $(filesize(gp)) bytes, listed in the unchanged tree above")
    stray = [p for p in (STORE_A, STORE_B, STORE_C, STORE_T) if startswith(abspath(p), abspath(joinpath(ROOT, "data")))]
    emit("D86y", "every store, manifest and temporary file of this script lies under the scratch folder", isempty(stray) ? "yes ($SCR)" : stray, "yes", pf(isempty(stray)))
end

say("# summary: PASS $(COUNTS["PASS"]) FAIL $(COUNTS["FAIL"]) INFO $(COUNTS["INFO"]), total seconds $(round(time() - T0; digits = 1)), peak RSS of this process $(round(Sys.maxrss() / 2^30; digits = 2)) GiB")
close(LOG)
