# D53p manifest, D53u layering and entry point, D53x isolation (proxy). Spec: Amendment D53, R7; Compute; D16.

# a scratch git repository with one commit, one tracked file and an ignored folder
function ex_scratch_repo()
    dir = mktempdir()
    git(args...) = run(pipeline(`git -C $dir -c user.name=ex -c user.email=ex@example.com $args`; stdout = devnull, stderr = devnull))
    git("init")
    write(joinpath(dir, "tracked.txt"), "one\n")
    write(joinpath(dir, ".gitignore"), "ignored/\n")
    git("add", ".")
    git("commit", "-m", "init")
    return dir
end
ex_head(dir) = String(strip(read(`git -C $dir rev-parse HEAD`, String)))

const EX_REC = (threads = 4, cpu_ns_per_step = 8.25, gpu_ns_per_step = 2.47, gpu_faster = true, date = "2026-10-17")
const EX_PRIOR = (threads = 4, cpu_ns_per_step = 8.0, gpu_ns_per_step = 2.5, gpu_faster = true, date = "2026-10-16")

@testset "D53p manifest (R7)" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    repo = ex_scratch_repo()
    store = ex_store("p-clean")
    res = ex_campaign(fx; store = store, repo = repo, allow_dirty = false, device_record = EX_REC, prior_record = EX_PRIOR)
    @test res.status == "complete"
    man = EX.read_json(res.manifest)
    top = ("study", "schema_version", "plan_sha256", "plan_counts", "status", "reasons", "sessions", "git", "versions", "host", "threads",
           "backend_table", "device", "device_record", "prior_record", "prereg", "panel", "block", "defaults", "launch", "chunks",
           "arms", "wall", "overrides", "removed_temporaries")
    @test all(k -> haskey(man, k), top)
    @test man["study"] == "ex-fixture" && man["schema_version"] == 2 && man["plan_sha256"] == plan.hash
    @test man["plan_counts"]["items"] == plan.counts.items && man["plan_counts"]["chunks"] == plan.counts.chunks
    @test man["plan_counts"]["steps_before"] == plan.steps.total_before && man["plan_counts"]["steps_after"] == plan.steps.total_after
    # git: commit of the session, clean-tree flag
    @test man["git"]["sha"] == ex_head(repo) && man["git"]["clean"] === true
    s = only(man["sessions"])
    @test s["git_sha"] == ex_head(repo) && s["clean"] === true
    @test all(k -> occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$", s[k]), ("start_utc", "end_utc"))
    # versions, host, threads
    @test all(k -> haskey(man["versions"], k), ("language", "packages_sha256", "hdf5", "module"))
    @test man["versions"]["language"] == string(VERSION) && all(k -> haskey(man["host"], k), ("cpu", "cores"))
    @test man["threads"]["used"] == 4 && man["threads"]["available"] == Threads.nthreads()
    # backend table of R4: one row per (precision, method) of the plan
    rows = man["backend_table"]
    @test Set((r["precision"], r["method"]) for r in rows) == Set((ex_label(it.precision), it.method) for it in plan.items)
    @test all(r -> all(k -> haskey(r, k), ("request", "backend", "reason", "note")), rows) && all(r -> r["backend"] == "cpu", rows)
    # both K19 records
    @test man["device_record"]["threads"] == 4 && man["device_record"]["gpu_faster"] === true && man["device_record"]["date"] == "2026-10-17"
    @test man["device_record"]["cpu_ns_per_step"] == 8.25 && man["device_record"]["gpu_ns_per_step"] == 2.47
    @test man["prior_record"]["date"] == "2026-10-16" && man["prior_record"]["threads"] == 4
    # inputs: hashes at start and end, header value, build hash of the panel, the block as parsed, defaults, launch inputs
    hdr = EX.check_prereg(fx.prereg).header_sha256
    @test man["prereg"]["sha256_start"] == man["prereg"]["sha256_end"] == man["prereg"]["header"] == hdr
    @test man["panel"]["sha256_start"] == man["panel"]["sha256_end"] == ex_sha256(read(fx.panel)) && man["panel"]["build_hash"] == EX_PANEL_HASH
    @test man["block"]["sha256"] == ex_sha256(man["block"]["text"]) == EX.block_hash(EX.read_run_block(fx.prereg))
    @test all(k -> haskey(man["defaults"], k), ("seed", "bigfloat_bits", "renorm_interval_tau", "precision", "backend",
                                                  "orthogonality_tol", "pairing_residual_tol", "det_residual_tol"))
    @test all(k -> haskey(man["launch"], k), ("free_bytes", "f", "capacity")) && man["launch"]["capacity"] == 64
    # chunk list in plan order; statuses sum to the plan's chunk count (I17)
    ch = man["chunks"]
    @test [c["ordinal"] for c in ch] == collect(1:plan.counts.chunks) && [c["key"] for c in ch] == [c.key for c in plan.chunks]
    need = ("ordinal", "key", "hash", "arms", "status", "attempt", "backend", "fallback_reason", "orbits", "steps", "wall_seconds",
            "first_attempt_reasons", "flags", "setup_check")
    @test all(c -> all(k -> haskey(c, k), need), ch)
    @test all(c -> all(k -> haskey(c["flags"], k), ("diagnostics exceeded", "drift", "excluded")), ch)
    @test count(c -> c["status"] == "ok", ch) == plan.counts.chunks
    # per-arm totals and the wall time against the cap
    for arm in fx.run["arm"]
        a = man["arms"][arm["name"]]
        @test all(k -> haskey(a, k), ("steps", "wall_seconds", "thread_seconds", "cost_ns_per_orbit_step", "budget_ns", "status"))
    end
    @test haskey(man["wall"], "total_seconds") && man["wall"]["cap_hours"] == 0.125
    @test isempty(man["overrides"]) && isempty(man["removed_temporaries"])
    # absent records are written as "absent"
    man0 = EX.read_json(ex_campaign(fx; store = ex_store("p-norec"), repo = repo, allow_dirty = false).manifest)
    @test man0["device_record"] == "absent" && man0["prior_record"] == "absent"
end

@testset "D53p clean tree, override and the thread override" begin
    fx = ex_fix(:small)
    repo = ex_scratch_repo()
    run_it(; kw...) = ex_campaign(fx; store = ex_path("p-dirty-store"), manifest = ex_path("p-dirty-man.json"), repo = repo, allow_dirty = false, kw...)
    clean() = (rm(ex_path("p-dirty-store"); recursive = true, force = true); rm(ex_path("p-dirty-man.json"); force = true))
    # a tracked change refuses
    write(joinpath(repo, "tracked.txt"), "two\n")
    clean()
    @test ex_rejects(() -> run_it(), ["tracked"]) && !ispath(ex_path("p-dirty-store")) && !ispath(ex_path("p-dirty-man.json"))
    res = run_it(allow_dirty = true)
    man = EX.read_json(res.manifest)
    @test man["git"]["clean"] === false && "tracked.txt" in man["overrides"]["dirty"]["paths"]
    @test occursin(r"^[0-9a-f]{64}$", man["overrides"]["dirty"]["diff_sha256"])
    run(`git -C $repo checkout -- tracked.txt`)
    # an untracked, unignored file refuses; an ignored file does not
    mkpath(joinpath(repo, "ignored")); write(joinpath(repo, "ignored", "x.txt"), "x\n")
    clean()
    @test run_it().status == "complete"
    write(joinpath(repo, "new.txt"), "n\n")
    clean()
    @test ex_rejects(() -> run_it(), ["untracked"]) && !ispath(ex_path("p-dirty-store"))
    clean()
    man = EX.read_json(run_it(allow_dirty = true).manifest)
    @test "new.txt" in man["overrides"]["dirty"]["paths"] && only(man["sessions"])["clean"] === false
    # a thread count other than [run.k19].threads runs only with the explicit override, which is recorded
    clean()
    man = EX.read_json(ex_campaign(fx; store = ex_path("p-dirty-store"), threads = 2).manifest)
    @test man["overrides"]["threads"]["used"] == 2 && man["overrides"]["threads"]["k19"] == 4
    clean()
    man = EX.read_json(ex_campaign(fx; store = ex_path("p-dirty-store"), threads = 4).manifest)
    @test !haskey(man["overrides"], "threads")
end

@testset "D53p an input changed between launch groups stops the campaign" begin
    for (what, status) in (("prereg", "stopped: preregistration changed"), ("panel", "stopped: panel changed"))
        fx = ex_fix_with(ex_rb(:small), "p-change-$what")
        path = what == "prereg" ? fx.prereg : fx.panel
        orig = read(path)
        touch_it(pt, info) = (pt === :group_start && info.group == 2 &&
                              write(path, what == "prereg" ? replace(String(orig), "End of file." => "End of file!") : String(orig) * "\n# touched\n"); nothing)
        store = ex_store("p-change-$what")
        res = ex_campaign(fx; store = store, hooks = EX.RunHooks(kill = touch_it))
        @test res.status == status
        @test length(ex_chunk_files(store)) == 4                      # group 1 only: a group in flight finishes, none starts
        man = EX.read_json(res.manifest)
        @test man["status"] == status && man[what]["sha256_start"] != man[what]["sha256_end"]
        @test count(c -> c["status"] == "ok", man["chunks"]) == 4 && count(c -> c["status"] == "pending", man["chunks"]) == 14
        write(path, orig)                                             # restored: the store is resumable
        res2 = ex_campaign(fx; store = store)
        @test res2.status == "complete" && sort(res2.ran) == collect(5:18)
    end
end

@testset "D53u layering and declared dependencies" begin
    root = joinpath(@__DIR__, "..", "..")
    deps = TOML.parsefile(joinpath(root, "Project.toml"))["deps"]
    imported = Set{String}()
    for f in readdir(joinpath(root, "src", "Experiments"); join = true)
        endswith(f, ".jl") || continue
        for line in eachline(f)
            m = match(r"^\s*(?:using|import)\s+(.*)$", line)
            m === nothing && continue
            body = String(split(m.captures[1], '#')[1])
            for part in (occursin(':', body) ? [first(split(body, ':'))] : split(body, ','))       # "using X: a, b" imports X only
                name = String(strip(split(strip(part), r"[:.\s]")[1]))
                isempty(name) || push!(imported, name)
            end
        end
    end
    @test haskey(deps, "HDF5") && imported <= Set(keys(deps))                                 # nothing outside the lower layers (relative imports) and the declared packages
end

@testset "D53u entry point" begin
    root = abspath(joinpath(@__DIR__, "..", ".."))
    script = joinpath(root, "bin", "chaoslab.jl")
    @test isfile(script)
    code = [strip(l) for l in eachline(script) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    @test length(code) <= 12 && any(l -> occursin("Experiments.cli", l), code)
    @test !any(l -> occursin(r"\b(for|while|sqrt|log|exp|sin|cos)\b", l), code)       # it parses arguments and calls the module: no numerics
    # status to exit code: zero only for complete; every stop and pause is nonzero
    @test EX.exit_status("complete") == 0
    @test EX.exit_status("complete with failed chunks") != 0
    for s in ("stopped: control bound", "stopped: failed chunks", "stopped: wall cap", "stopped: preregistration changed",
              "stopped: panel changed", "paused arms")
        @test EX.exit_status(s) != 0
    end
    fx = ex_fix(:small)
    io = IOBuffer()
    common = ["--code-version", EX_CODE, "--threads", "4", "--thread-override", "--allow-dirty", "--capacity", "64"]
    args(f, store) = vcat(["run", f.prereg, f.panel, "--store", store, "--manifest", store * "-manifest.json"], common)
    @test EX.cli(args(fx, ex_store("u-run")); io = io) == 0
    @test length(ex_chunk_files(ex_path("u-run"))) == ex_plan(fx).counts.chunks
    # a stop status gives a nonzero exit: wall cap 1e-12 h stops after the first launch group; cost factor 1e-9 pauses an arm
    cap = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "wall_cap_hours") => 1e-12), "u-cap")
    @test EX.cli(args(cap, ex_store("u-cap")); io = io) != 0
    pause = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "first_chunk_cost_factor") => 1e-9), "u-pause")
    @test EX.cli(args(pause, ex_store("u-pause")); io = io) != 0
    # a refusal and an unknown command: exit status 1 with a message
    bad = ex_path("u-bad-prereg.md"); write(bad, replace(read(fx.prereg, String), "End of file." => "End of file!"))
    io = IOBuffer()
    @test EX.cli(["run", bad, fx.panel, "--store", ex_path("u-bad"), "--manifest", ex_path("u-bad.json")]; io = io) == 1
    @test occursin("SHA-256", String(take!(io))) && !ispath(ex_path("u-bad"))
    io = IOBuffer()
    @test EX.cli(["frobnicate"]; io = io) == 1 && occursin("run", String(take!(io)))
    # pilot and panel through the same entry point
    @test EX.cli(["pilot", fx.prereg, fx.panel, "--out", ex_path("u-pilot"), "--code-version", EX_CODE, "--threads", "4",
                  "--thread-override", "--allow-dirty", "--capacity", "64"]; io = IOBuffer()) == 0
    @test isfile(joinpath(ex_path("u-pilot"), "pilot-report.json"))
    pre = ex_prereg_for()
    @test EX.cli(["panel", pre, "--out", ex_path("u-panel.toml")]; io = IOBuffer()) == 0
    @test read(ex_path("u-panel.toml")) == read(ex_build("base").path)                # same bytes as build_panel
    # as a script: no arguments is a usage error with a nonzero exit status
    p = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --project=$root $script`); stdout = devnull, stderr = devnull))
    @test p.exitcode != 0
end

@testset "D53x isolation: no test writes to the study's data folder" begin
    fx = ex_fix(:small)
    d = EX.default_store(fx.prereg)                                   # <data root>/<study slug>
    @test endswith(replace(d, '\\' => '/'), "/data/ex-fixture")
    @test !isdir(d) || isempty(readdir(d))                            # 0 files after the whole suite; the reduced campaign itself is stage 04
    @test EX.default_store(fx.prereg; root = "R") == joinpath("R", "ex-fixture")
end
