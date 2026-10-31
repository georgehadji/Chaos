# Included by verify/experiments-d94.jl: the panels and stores exported for the independent oracle (REDUCED, FIXA, O26, THREADS, COST, README).
# Needs the definitions of the main script (XP, V, HDIR, say, emit, info, section, sha_file) and, for FIXA, O26 and THREADCHILD, the test helpers.

function write_prereg(path, block, what)
    io = IOBuffer()
    println(io, "# Verification input: ", what, "\n\nWritten by verify/experiments-d94.jl. Not a study; the panel block is the only part the module reads.\n\n## Header\n")
    println(io, "| Field | Value |\n|-------|-------|\n| Study slug | d94-verify |")
    println(io, "| SHA-256 of this file (current version) | ", "0"^64, " (rule: SELF) |")
    println(io, "| Amendment log | A1 (2026-10-31, verification input) |\n")
    println(io, "## Panel\n\n```toml")
    TOML.print(io, Dict("panel" => block))
    println(io, "```\n\nEnd of file.")
    text = String(take!(io))
    write(path, replace(text, "0"^64 => XP.prereg_hash(text); count = 1))
    return path
end

# the P1 panel block with the overrides of verify/experiments-d94-reduced.toml (drop the scalar n_c: `quota` replaces it)
function reduced_block(; file = joinpath(@__DIR__, "experiments-d94-reduced.toml"), over = nothing)
    text = replace(read(P1_PREREG, String), "\r\n" => "\n")
    base = TOML.parse(match(r"```toml\n(.*?)\n```"s, text).captures[1])["panel"]
    b = merge(base, over === nothing ? TOML.parsefile(file)["panel"] : over)
    delete!(b, "n_c")
    return b
end

# D and D_B by the library functions, written as TOML for the oracle author
function export_db(panel_path, dir)
    P = TOML.parsefile(panel_path)
    man, orbs = P["manifest"], P["orbit"]
    seed, b, rmax = Int(man["inputs"]["seed"]["value"]), Int(man["db_size"]), Int(man["r_max"])
    out = Dict{String, Any}("seed" => seed, "db_size" => b, "r_max" => rmax)
    for R in 1:rmax
        DR = [(id = Int(o["id"]), k_index = Int(o["k_index"])) for o in orbs if get(o, "in_D", false) && o["round"] <= R]
        db = XP.db_subset(seed, R, DR, b)
        out["R$R"] = Dict{String, Any}("D_size" => length(DR), "D_ids" => sort!([x.id for x in DR]), "DB_size" => length(db), "DB_ids" => db,
                                       "DB_per_k_index" => [count(x -> x.id in db && x.k_index == i, DR) for i in 1:length(man["K_binary64"])])
        ok = length(db) == min(b, length(DR)) && issubset(db, [x.id for x in DR])
        emit("D94q", "D_B for R = $R (library function), |D_R| / |D_B| / per K index $(out["R$R"]["DB_per_k_index"])", "$(length(DR)) / $(length(db))", "|D_B| = min($b, |D_R|), subset of D_R", pf(ok))
    end
    path = joinpath(dir, "p3-reduced-db.toml")
    open(io -> TOML.print(io, out; sorted = true), path, "w")
    return path
end

section("REDUCED") do
    dir = joinpath(V, "reduced"); mkpath(dir)
    pre = write_prereg(joinpath(dir, "p3-reduced-prereg.md"), reduced_block(), "reduced real panel of D94q and D94r")
    out = joinpath(dir, "p3-reduced-panel.toml")
    if get(ENV, "EXV_REUSE", "") == "1" && isfile(out)
        info("D94q", "EXV_REUSE=1: the reduced panel already in the scratch folder is used (its build time is in the first REDUCED lines of this log)")
    else
        foreach(f -> rm(f; force = true), (out, out * ".pools.toml", out * ".provenance.toml"))
        t = @elapsed XP.build_panel(pre; out = out)
        emit("D94q", "reduced real panel built (real K list, seed 20261030, quota and fraction tables, db_size 100, screen_points 4096, r_max 3, controls none), threads $(Threads.nthreads()), wall seconds", t, "completes", "PASS")
    end
    P = TOML.parsefile(out); man, orbs = P["manifest"], P["orbit"]
    tot = man["totals"]
    info("D94q", "orbits $(tot["orbits"]) (drawn $(tot["drawn"]), controls $(tot["controls"])), in D $(count(o -> get(o, "in_D", false), orbs)), panel SHA-256 $(sha_file(out))")
    emit("D94q", "no control orbit with controls = none", tot["controls"], "0", pf(tot["controls"] == 0 && isempty(man["controls"])))
    cells = man["round_cells"]
    for st in ("populated", "short", "empty", "not drawn", "undefined")
        info("D94q", "round cells with status $st: $(count(c -> c["status"] == st, cells)) of $(length(cells))")
    end
    emit("D94q", "K = 0 orbits only in the gate round 1 (orbits at K index 1 outside round 1)", count(o -> o["k_index"] == 1 && o["round"] != 1, orbs), "0", pf(all(o -> o["k_index"] != 1 || o["round"] == 1, orbs)))
    export_db(out, dir)
end

# ---------------------------------------------------------------------------------------------------
# Fixture A of the test helpers (ex_d94_helpers.jl) built into V, and its plan in the four states of the log (open at R_plan 1, 2, 3 and closed)
function fixa_plans(b)
    run = d94_block()
    run["panel_prereg_sha256_accepted"] = [XP.check_prereg(b.prereg).computed_sha256]
    pre = ex_write_prereg_run(ex_path("fixa-run-prereg.md"), run)
    fx = (; panel = b.path, prereg = pre)
    logp = ex_path("fixa.log"); write(logp, "")
    rnd(c) = parse(Int, match(r"cell=\S+:r(\d+)$", c.key).captures[1])
    texts, specs = String[], NamedTuple[]
    for (r, dec) in ((1, "continue"), (2, "continue"), (3, "close"))
        plan = XP.plan_run(pre, b.path; code_version = EX_CODE, round_log = logp)
        push!(texts, "# plan with $(r - 1) log entries\n" * XP.plan_text(plan))
        digest = q -> d94_digest([c.hash for c in plan.chunks if "ref" in c.arms && rnd(c) <= q])
        push!(specs, (round = r, decision = dec))
        d94_write_log(logp, fx, specs, digest)
    end
    push!(texts, "# closed plan\n" * XP.plan_text(XP.plan_run(pre, b.path; code_version = EX_CODE, round_log = logp)))
    return join(texts, "\n=====\n")
end

section("FIXA") do
    HDIR[] = joinpath(V, "fixtureA"); mkpath(HDIR[])
    t = @elapsed b = d94_build("A", D94_FA)
    emit("D94q", "fixture panel A built into $(HDIR[]), wall seconds", t, "completes", "PASS")
    info("D94q", "fixture A: seed $(D94_FA.seed), K $(D94_FA.K), M $(D94_FA.M), r_max $(D94_FA.r_max), orbits $(b.panel["manifest"]["totals"]["orbits"]), panel SHA-256 $(sha_file(b.path))")
    write(ex_path("plan.txt"), fixa_plans(b))
    info("D94q", "plan texts of fixture A (R_plan 1, 2, 3 open and closed) written to plan.txt")
end

# ---------------------------------------------------------------------------------------------------
# O26: a small campaign by rounds (the hand panel and block of the helpers) with the round log written between the runs
section("O26") do
    HDIR[] = joinpath(V, "o26"); mkpath(HDIR[])
    fx = d94_std()
    store = ex_store("camp"); logp = ex_path("camp.log"); write(logp, "")
    digest = r -> d94_digest(d94_ref_hashes(store, r))
    specs, lines = NamedTuple[], String[]
    for (stage, dec) in ((1, "continue"), (2, "continue"), (3, "close"), (4, ""))
        res = ex_campaign(fx; store = store, round_log = logp)
        files = ex_chunk_files(store)
        push!(lines, "after run $stage: " * join(["$(basename(f)) $(d94_attr(f, "chunk_hash"))" for f in files], "; "))
        stage == 4 && break
        push!(specs, (round = stage, decision = dec))
        d94_write_log(logp, fx, specs, digest)
        cp(logp, ex_path("camp-log-after-round$stage.txt"); force = true)
        emit("O26", "run $stage of the campaign finished, chunk files in the store, log entries written", "$(length(files)), $(stage)", "-", "INFO")
    end
    write(ex_path("camp-store-listing.txt"), join(lines, "\n") * "\n")
    n = length(ex_chunk_files(store))
    emit("O26", "campaign through 3 rounds and closure: chunk files in the store", n, "> 0", pf(n > 0))
    for r in 1:3
        emit("O26", "ref arm digest of rounds 1 to $r from the chunk_hash attributes (the helper's HDF5 read)", d94_digest(d94_ref_hashes(store, r)), "-", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# Process-level thread counts: fixture A and its plans built under julia -t 1, 4, 23
section("THREADCHILD") do
    HDIR[] = ENV["EXV_HDIR"]; mkpath(HDIR[])
    b = d94_build("A", D94_FA)
    write(ex_path("plan.txt"), fixa_plans(b))
    say("THREADCHILD | threads $(Threads.nthreads()) panel SHA-256 $(sha_file(b.path)) | -")
end
section("THREADS") do
    files = ("d94-A-panel.toml", "d94-A-panel.toml.pools.toml", "plan.txt")
    dirs = Dict{Int, String}()
    for n in (1, 4, 23)
        dir = joinpath(V, "threads", "t$n"); rm(dir; recursive = true, force = true); mkpath(dir); dirs[n] = dir
        env = merge(Dict(ENV), Dict("EXV_ONLY" => "THREADCHILD", "EXV_HDIR" => dir, "EXV_LOG" => joinpath(dir, "child-log.txt"), "EXV_V" => V))
        cmd = setenv(`$(Base.julia_cmd()) --startup-file=no -t $n --project=$(dirname(Base.active_project())) $(joinpath(@__DIR__, "experiments-d94.jl"))`, env)
        t = @elapsed run(pipeline(cmd; stdout = devnull))
        prov = TOML.parsefile(joinpath(dir, "d94-A-panel.toml.provenance.toml"))
        info("THREADS", "julia -t $n: child $(round(t; digits = 1)) s, provenance threads $(prov["threads"]), ntasks $(prov["ntasks"])")
    end
    for f in files
        shas = [sha_file(joinpath(dirs[n], f)) for n in (1, 4, 23)]
        emit("D94p", "$f from julia -t 1, 4, 23 byte-identical (SHA-256 $(first(shas[1], 12)))", length(unique(shas)), "1 distinct", pf(length(unique(shas)) == 1))
    end
end

# ---------------------------------------------------------------------------------------------------
# D94s cost: one K at the real screen_points (262144) with the real quota and fraction tables; the study screens six K in turn
section("COST") do
    dir = joinpath(V, "cost"); mkpath(dir)
    seed = 20261030
    t_pts = minimum(@elapsed(XP.random_screen_points(seed, i, 262144)) for i in 1:3)
    emit("D94s", "random_screen_points at M = 262144 (one K, minimum of 3), seconds", t_pts, "reported", "INFO")
    rss0 = Sys.maxrss() / 2^20
    over = Dict{String, Any}("K" => ["0.971635406"], "seed" => seed, "screen_source" => "random", "screen_points" => 262144, "r_max" => 1, "gate_rounds" => [1],
                             "controls" => "none", "db_size" => 100,
                             "quota" => Dict("1" => Dict("T" => 40, "R" => 40, "X" => 40, "C" => 100, "S" => 250)),
                             "d_fraction" => Dict("1" => Dict(c => "0.30" for c in ("T", "R", "X", "C", "S"))))
    pre = write_prereg(joinpath(dir, "cost-prereg.md"), reduced_block(over = over), "one K at the real screening size")
    out = joinpath(dir, "cost-panel.toml")
    t = @elapsed XP.build_panel(pre; out = out)
    rss = Sys.maxrss() / 2^20
    emit("D94s", "one K (K_c) at screen_points 262144, N_scr 1e5, threads $(Threads.nthreads()): build_panel wall seconds (Newton search included)", t, "reported", "INFO")
    est_wall = 6 * t / 3600
    emit("D94s", "ESTIMATE (6 K in turn, K = 0 screened too, P71; one K measured x 6): wall hours on $(Threads.nthreads()) threads", est_wall, "<= 2 h", pf(est_wall <= 2))
    extra = 5 * (2 * 8 + 1) * 262144 / 2^20
    emit("D94s", "ESTIMATE of peak memory, MiB (maxrss after one K $(round(rss; digits = 0)) + 5 more K of points and bits $(round(extra; digits = 1)); before $(round(rss0; digits = 0)))", rss + extra, "<= 8192", pf(rss + extra <= 8192))
    P = TOML.parsefile(out)
    info("D94s", "one-K panel: orbits $(P["manifest"]["totals"]["orbits"]), pools $(P["manifest"]["pools"]), status of the round cells $(join(["$(c["class"]) $(c["status"]) pool $(c["pool"]) drawn $(c["drawn"])" for c in P["manifest"]["round_cells"]], "; "))")
    info("D94s", "the real test panel is NOT built: it needs the amended preregistration of the study (P69), which does not exist; D94s is not run")
end

# ---------------------------------------------------------------------------------------------------
# README.txt of the scratch folder: every exported path and what it holds (for the author of the independent oracle)
const README_ITEMS = [
    ("generator.txt", "O23 generator check: 4 seeds (0, 1, 20261012, 2^63+5) x 250 words rand(StableRNG(seed), UInt64) and 250 floats rand(StableRNG(seed)); line: <seed> <u64|f64> <index> <value>; each kind from a fresh stream"),
    ("reduced/p3-reduced-panel.toml", "reduced real panel (D94q, D94r): real K list (0, 0.5, 0.971635406, pi/2, 2pi, 10), seed 20261030, p3 quota and D-fraction tables, db_size 100, screen_points 4096, r_max 3, gate round 1, controls none; [manifest] + [[orbit]]"),
    ("reduced/p3-reduced-panel.toml.pools.toml", "its pool sidecar: [k.<i>] points, bits_sha256, T, C, S as [first, last] runs of point ids"),
    ("reduced/p3-reduced-panel.toml.provenance.toml", "its provenance sidecar (threads, ntasks, backend)"),
    ("reduced/p3-reduced-prereg.md", "the input the panel was built from (the P1 panel block with the overrides of verify/experiments-d94-reduced.toml; the inputs are also in the panel manifest)"),
    ("reduced/p3-reduced-db.toml", "D and D_B per R = 1, 2, 3 from the library functions: [R<r>] D_size, D_ids, DB_size, DB_ids, DB_per_k_index (for comparison with the oracle's D_B; the panel file holds only the in_D flags)"),
    ("fixtureA/d94-A-panel.toml", "fixture panel A of test/experiments/ex_d94_helpers.jl (seed 20261012, K 0, 1.5, 5, 64 points per K, r_max 3, gate [1], controls none, db_size 6)"),
    ("fixtureA/d94-A-panel.toml.pools.toml", "its pool sidecar"),
    ("fixtureA/d94-A-panel.toml.provenance.toml", "its provenance sidecar"),
    ("fixtureA/plan.txt", "plan texts of fixture A with a block of arms ref (P, each_round), dd (D, closure), db (D_B, closure): 0, 1, 2 log entries (R_plan 1, 2, 3, open) and the closed plan; the logs are fixtureA/fixa.log and fixtureA/fixa-run-prereg.md"),
    ("threads/t1, t4, t23", "fixture A panel, pool sidecar and plan.txt built by julia -t 1, 4, 23 (identical bytes; provenance differs in threads)"),
    ("o26/camp", "O26 store: chunk files chunk-NNNNN.h5 of the hand-panel campaign (panel o26/d94-std-panel.toml, run prereg o26/d94-std-prereg.md, arms ref P L1 each_round, dd D L5 closure, db D_B L4 closure), 3 rounds then closure; read root attributes chunk_hash, chunk_key, arms with HDF5 tools"),
    ("o26/camp.log", "the round log in its final state (3 entries: continue, continue, close), format of the test plan (round, decision, counts, panel_sha256, prereg_sha256, ref_digest, prev_sha256); written by the test helper, not by an independent script"),
    ("o26/camp-log-after-round{1,2,3}.txt", "the log after each decision"),
    ("o26/camp-store-listing.txt", "chunk file names and chunk_hash after each of the 4 runs"),
    ("o26/camp-manifest.json*", "run manifests of the campaign"),
    ("p1/", "P1 rebuild (p1-rebuilt-panel*.toml), plan dumps (plan-7f70aa7.txt before D94, plan-new.txt HEAD), recorded keys of the accepted P1 run"),
    ("cost/", "one-K full-size screening build (D94s cost estimate)")]
section("README") do
    io = IOBuffer()
    println(io, "Scratch folder of the stage 04 verification of Amendment D94 (verify/experiments-d94.jl). All paths relative to this folder.\nNothing under data/ was read or written.\n")
    for (p, what) in README_ITEMS
        println(io, p, "\n    ", what, "\n    ", ispath(joinpath(V, split(p, ',')[1])) || contains(p, '{') || contains(p, '*') ? "exists" : "MISSING")
    end
    write(joinpath(V, "README.txt"), String(take!(io)))
    missing_ = [p for (p, _) in README_ITEMS if !(contains(p, '{') || contains(p, '*') || contains(p, ',') || ispath(joinpath(V, p)))]
    emit("README", "paths listed in README.txt / missing", "$(length(README_ITEMS)) / $(length(missing_))", "0 missing", pf(isempty(missing_)))
end
