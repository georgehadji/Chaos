# Amendment D53 review fixes (stage 04 independent review: H1, H2, M1 to M9, LOW). Spec: experiments-spec.md Amendment D53, R4 to R8, R10;
# schema/run-store.md. One testset per finding; the earlier assertions are untouched. Helpers: ex_run_helpers.jl, ex_device.jl.

rev_hooks(; kw...) = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0); kw...)

# A session killed after `k` chunks (no manifest reset allowed): returns the fixture, the store and the manifest path.
function rev_killed(tag::AbstractString; k::Integer = 2, fx = ex_fix(:small), kw...)
    store = ex_store("rev-$tag")
    man = store * "-manifest.json"
    rm(man; force = true)
    for suffix in (".1", ".2", ".3")
        rm(man * suffix; force = true)
    end
    @test ex_killed(() -> ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = rev_hooks(kill = ex_kill_after(k)), kw...))
    return (; fx, store, man)
end

# copy of a killed state under another tag (store and manifest files)
function rev_copy(run, tag::AbstractString)
    store = ex_store("rev-$tag")
    cp(run.store, store)
    man = store * "-manifest.json"
    for suffix in ("", ".1", ".2", ".3")
        isfile(run.man * suffix) ? cp(run.man * suffix, man * suffix; force = true) : rm(man * suffix; force = true)
    end
    return (; fx = run.fx, store, man)
end

rev_resume(run; kw...) = ex_campaign(run.fx; store = run.store, manifest = run.man, manifest_reset = false, hooks = rev_hooks(), kw...)
rev_state(store) = Dict(f => (mtime(f), filesize(f)) for f in readdir(store; join = true) if isfile(f))

@testset "D53 review H1: a control-bound stop is not a result" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    c5 = ex_find_chunk(plan; method = "L5", cell = "controls", arm = "main")
    bump(id) = (orbits, info) -> (o = only(o for o in orbits if o.id == id); first(values(o.lyap))[end, 1] += 1.0)
    store = ex_store("rev-h1")
    man = store * "-manifest.json"
    res = ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = rev_hooks(mutate = ex_fault(bump(14); ord = c5.ordinal)))
    @test res.status == "stopped: control bound" && isempty(res.failed)
    f = joinpath(store, "chunk-$(lpad(c5.ordinal, 5, '0')).h5")
    ch = EX.read_chunk(f)
    @test ch.attrs["status"] == "failed" && ch.attrs["control_stop"] === true && ch.attrs["complete"] === true
    @test all(o -> !("excluded" in o.attrs["flags"]), ch.orbits)                       # no orbit is excluded from the analysis
    @test all(f -> EX.read_chunk(f).attrs["control_stop"] === false, filter(!=(f), ex_chunk_files(store)))
    m = EX.read_json(man)
    stop = only(m["control_stops"])
    @test stop["ordinal"] == c5.ordinal && stop["code_version"] == EX_CODE && stop["retried"] === false && any(r -> occursin("orbit 14", r), stop["reasons"])
    @test m["chunks"][c5.ordinal]["control_stop"] === true && all(c -> isempty(c["flags"]["excluded"]), m["chunks"])
    # a plain resume refuses and touches nothing
    before = rev_state(store)
    @test ex_rejects(() -> ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = rev_hooks()), ["control bound", "retry_control_stop"])
    @test rev_state(store) == before && !isdir(joinpath(store, "aside"))
    # the stop survives a changed code version as a stale file (the old rule), and the explicit retry recomputes the chunk
    @test ex_rejects(() -> ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = rev_hooks(), code_version = "ex-code-2"), ["stale"])
    res2 = ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = rev_hooks(), retry_control_stop = true)
    @test res2.status == "complete" && c5.ordinal in res2.ran
    m2 = EX.read_json(man)
    @test m2["overrides"]["retry_control_stop"] == [c5.ordinal] && only(m2["control_stops"])["retried"] === true
    @test EX.read_chunk(f).attrs["status"] == "ok" && EX.read_chunk(f).attrs["control_stop"] === false
    @test any(n -> startswith(n, "chunk-$(lpad(c5.ordinal, 5, '0'))"), readdir(joinpath(store, "aside")))
end

@testset "D53 review H2: the backend table is compared at resume" begin
    fx = ex_fix(:small)
    rec = (device_record = EX_DREC,)
    run = ex_with_device(() -> rev_killed("h2"; fx = fx, rec...), ex_hdev())
    @test EX.read_chunk(ex_chunk_files(run.store)[1]).attrs["backend"] == "gpu"
    before = rev_state(run.store)
    bytes = read(run.man)
    # no device now: the table of the resumed run differs from the manifest's, the run is refused and nothing is written
    @test ex_rejects(() -> rev_resume(run), ["backend table", "backend_override"])
    @test rev_state(run.store) == before && read(run.man) == bytes
    res = rev_resume(run; backend_override = true)
    @test res.status == "complete"
    m = EX.read_json(run.man)
    @test !isempty(m["overrides"]["backend_table"]) && any(p -> occursin("Float64", p), m["overrides"]["backend_table"])
    @test EX.read_chunk(ex_chunk_files(run.store)[4]).attrs["backend"] == "cpu"
    # the same table is accepted
    same = ex_with_device(() -> rev_killed("h2-same"; fx = fx, rec...), ex_hdev())
    res = ex_with_device(() -> rev_resume(same; rec...), ex_hdev())
    @test res.status == "complete" && !haskey(EX.read_json(same.man)["overrides"], "backend_table")
    # without a manifest the kept chunk files alone are checked against the table
    lost = ex_with_device(() -> rev_killed("h2-lost"; fx = fx, rec...), ex_hdev())
    for sfx in ("", ".1", ".2", ".3")
        rm(lost.man * sfx; force = true)
    end
    @test ex_rejects(() -> ex_campaign(lost.fx; store = lost.store, manifest = lost.man, manifest_reset = true, hooks = rev_hooks()), ["was written on gpu"])
end

@testset "D53 review M1: renames retry, nothing is overwritten, the manifest survives a failed replace" begin
    busy(code) = Base.IOError("rename", code)
    # transient EBUSY and EACCES are retried; another error is not
    dir = mktempdir()
    src, dst = joinpath(dir, "a"), joinpath(dir, "b")
    write(src, "x")
    for code in (Base.UV_EBUSY, Base.UV_EACCES)
        calls = Ref(0)
        EX._move_retry(src, dst; pause = 0.0, mover = (a, b) -> (calls[] += 1; calls[] < 3 ? throw(busy(code)) : Base.Filesystem.rename(a, b)))
        @test calls[] == 3 && read(dst, String) == "x"
        mv(dst, src)
    end
    calls = Ref(0)
    @test_throws Base.IOError EX._move_retry(src, dst; pause = 0.0, mover = (a, b) -> (calls[] += 1; throw(busy(Base.UV_ENOENT))))
    @test calls[] == 1
    calls[] = 0
    @test_throws Base.IOError EX._move_retry(src, dst; pause = 0.0, tries = 4, mover = (a, b) -> (calls[] += 1; throw(busy(Base.UV_EBUSY))))
    @test calls[] == 4
    # a failed replace keeps the previous manifest whole, a good one keeps numbered backups (newest first)
    path = joinpath(dir, "man.json")
    EX._write_json(path, Dict("v" => 1); keep = 3)
    @test_throws Base.IOError EX._write_json(path, Dict("v" => 2); keep = 3, mover = (a, b) -> throw(busy(Base.UV_EBUSY)))
    @test EX.read_json(path) == Dict("v" => 1)
    for v in 2:5
        EX._write_json(path, Dict("v" => v); keep = 3)
    end
    @test EX.read_json(path)["v"] == 5 && [EX.read_json("$path.$k")["v"] for k in 1:3] == [4, 3, 2] && !isfile("$path.4")
    # a chunk file whose final name is taken is moved aside first, then renamed without force
    ref = ex_small_run("ref")
    ch = EX.read_chunk(ex_chunk_files(ref.store)[1])
    out = joinpath(dir, "chunk-00001.h5")
    EX.write_chunk(out, ch.attrs, ch.orbits)
    EX.write_chunk(out, merge(ch.attrs, Dict{String, Any}("created_utc" => "2026-01-01T00:00:00Z")), ch.orbits)
    @test EX.read_chunk(out).attrs["created_utc"] == "2026-01-01T00:00:00Z" && EX.read_chunk(out * ".old1").attrs["created_utc"] == ch.attrs["created_utc"]
    @test isempty(filter(endswith(".tmp"), readdir(dir)))
end

@testset "D53 review M2: root attributes in the read-back, flush before the rename" begin
    ref = ex_small_run("ref")
    ch = EX.read_chunk(ex_chunk_files(ref.store)[1])
    dir = mktempdir()
    out = joinpath(dir, "chunk-00001.h5")
    bad = merge(ch.attrs, Dict{String, Any}("note" => "a\0b"))               # HDF5 strings end at the NUL: the read-back differs
    err = try
        EX.write_chunk(out, bad, ch.orbits); nothing
    catch e
        e
    end
    @test err isa ErrorException && occursin("root attributes", err.msg) && occursin("note", err.msg)
    @test !isfile(out)                                                       # no final name for a file that failed its read-back
    @test EX._attr_diff(Dict("a" => 1, "b" => [1.0, NaN]), Dict("a" => 1, "b" => [1.0, NaN])) == String[]
    @test EX._attr_diff(Dict("a" => 1), Dict("a" => 2, "c" => 1)) == ["a", "c"]
    EX.write_chunk(out, ch.attrs, ch.orbits)
    @test EX.read_chunk(out).attrs["complete"] === true
    @test_throws Exception EX._fsync(joinpath(dir, "missing.h5"))
    @test (EX._fsync(out); true)
end

@testset "D53 review M3: the session clock moves at every manifest write" begin
    fx = ex_fix(:small)
    store = ex_store("rev-m3")
    man = store * "-manifest.json"
    t = Ref(0.0)
    log = Tuple{Symbol, Float64}[]
    clock = (ev, info) -> (ev === :manifest_write && (t[] += 1.0); push!(log, (ev, t[])); t[])
    @test ex_killed(() -> ex_campaign(fx; store = store, manifest = man, manifest_reset = false, hooks = EX.RunHooks(clock = clock, kill = ex_kill_after(2))))
    nw = count(x -> x[1] === :manifest_write, log)
    m = EX.read_json(man)
    @test nw >= 3 && only(m["sessions"])["wall_seconds"] == nw && m["wall"]["total_seconds"] == nw      # the last write counts its own time
end

@testset "D53 review M4: chunk files without a manifest" begin
    run = rev_killed("m4")
    # an unreadable manifest and no backup: refused, nothing written
    bad = rev_copy(run, "m4-bad")
    for sfx in (".1", ".2", ".3")
        rm(bad.man * sfx; force = true)
    end
    write(bad.man, "{ not json")
    before = rev_state(bad.store)
    @test ex_rejects(() -> rev_resume(bad), ["manifest_reset", "manifest"])
    @test rev_state(bad.store) == before && read(bad.man, String) == "{ not json"
    # explicit option: goes on, recorded
    res = rev_resume(bad; manifest_reset = true)
    @test res.status == "complete" && EX.read_json(bad.man)["overrides"]["manifest_reset"] === true
    # a numbered backup replaces the missing or unreadable manifest: the sessions of the earlier run are kept
    for (tag, damage) in (("missing", f -> rm(f)), ("garbage", f -> write(f, "][")))
        c = rev_copy(run, "m4-$tag")
        @test isfile(c.man * ".1")
        damage(c.man)
        res = rev_resume(c)
        m = EX.read_json(c.man)
        @test res.status == "complete" && length(m["sessions"]) == 2 && !haskey(m["overrides"], "manifest_reset")
        @test any(r -> occursin("numbered", r) || occursin("backup", r), res.reasons)
    end
    # a first run on an empty or absent store needs no option
    @test ex_campaign(ex_fix(:noref); store = ex_store("rev-m4-first"), manifest_reset = false, hooks = rev_hooks()).status == "complete"
end

@testset "D53 review M5: code version covers the code and nothing else" begin
    root = mktempdir()
    for d in ("src/Sub", "ext", "docs")
        mkpath(joinpath(root, d))
    end
    write(joinpath(root, "src", "A.jl"), "a() = 1\n"); write(joinpath(root, "src", "Sub", "B.jl"), "b() = 2\n")
    write(joinpath(root, "ext", "E.jl"), "e() = 3\n"); write(joinpath(root, "Project.toml"), "name = \"X\"\n"); write(joinpath(root, "Manifest.toml"), "# lock\n")
    v0 = EX.code_version(root)
    for (rel, text) in (("src/CONTEXT.md", "notes\n"), ("src/A.jl.cov", "- a() = 1\n"), ("src/Sub/data.txt", "x\n"), ("ext/readme.md", "r\n"), ("docs/x.jl", "x() = 0\n"))
        write(joinpath(root, rel), text)
        @test EX.code_version(root) == v0                                      # not a *.jl file of src/ or ext/
    end
    write(joinpath(root, "src", "A.jl"), "a() = 1\r\n")                        # CRLF checkout: the same version
    @test EX.code_version(root) == v0
    write(joinpath(root, "ext", "E.jl"), "e() = 4\n")
    v1 = EX.code_version(root)
    @test v1 != v0                                                             # ext/ is code
    write(joinpath(root, "Project.toml"), "name = \"Y\"\n")
    @test EX.code_version(root) != v1                                          # Project.toml is code
    root2 = mktempdir()
    mkpath(joinpath(root2, "src")); write(joinpath(root2, "src", "A.jl"), "a() = 1\n")
    @test EX.code_version(root2) isa String                                    # no ext/, no Project.toml: still defined
end

@testset "D53 review M6: panel-derived values are part of the chunk hash" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    orbits = ex_panel_orbits(:std)
    function plan_with(tag, f)
        os = deepcopy(orbits)
        f(os)
        pp = ex_write_panel(ex_path("rev-m6-$tag-panel.toml"); orbits = os)
        return EX.plan_run(fx.prereg, pp; code_version = EX_CODE)
    end
    hyper = only(o for o in orbits if get(o, "stability", "") == "hyperbolic")
    # the analytic exponent of a control changes its chunk hash, not its item hash
    p = plan_with("exp", os -> (only(o for o in os if o["id"] == hyper["id"])["exponent"] += 1e-3))
    @test [it.hash for it in p.items] == [it.hash for it in plan.items]
    diff = [c.ordinal for c in plan.chunks if p.chunks[c.ordinal].hash != c.hash]
    @test !isempty(diff) && all(o -> hyper["id"] in plan.chunks[o].orbits, diff)
    # a hyperbolic control without a stability key and a positive exponent is the same plan; a zero exponent is refused
    p2 = plan_with("nostab", os -> delete!(only(o for o in os if o["id"] == hyper["id"]), "stability"))
    @test [c.hash for c in p2.chunks] == [c.hash for c in plan.chunks]
    ell = only(o for o in orbits if get(o, "stability", "") == "elliptic")
    @test ex_rejects(() -> plan_with("noell", os -> delete!(only(o for o in os if o["id"] == ell["id"]), "stability")), ["stability", string(ell["id"])])
    # a changed panel file (same orbits) makes the kept files stale
    ref = ex_small_run("ref")
    store = ex_store("rev-m6")
    cp(ref.store, store)
    pp = ex_path("rev-m6-bytes-panel.toml")
    write(pp, read(fx.panel, String) * "\n# bytes differ, orbits do not\n")
    fx2 = (; fx..., panel = pp)
    before = rev_state(store)
    @test ex_rejects(() -> ex_campaign(fx2; store = store), ["stale"]) && rev_state(store) == before
end

@testset "D53 review M7: the thread count must be the thread count of Julia" begin
    n = Threads.nthreads()
    plan(k19) = ex_plan(ex_fix_with(ex_rb_with(ex_rb(:small), ("k19", "threads") => k19), "rev-m7-$k19"))
    p = plan(n)
    @test EX._check_threads(p, n, false) === nothing
    @test ex_rejects(() -> EX._check_threads(plan(n + 1), n + 1, false), ["Threads.nthreads()", "thread_override"])
    @test EX._check_threads(plan(n + 1), n + 1, true) === nothing
    fx = ex_fix_with(ex_rb_with(ex_rb(:small), ("k19", "threads") => n + 1), "rev-m7-run")
    store = ex_store("rev-m7")
    @test ex_rejects(() -> EX.run_campaign(fx.prereg, fx.panel; store = store, manifest = store * ".json", threads = n + 1, allow_dirty = true,
                                           code_version = EX_CODE, capacity = 64), ["nthreads"]) && !ispath(store) && !ispath(store * ".json")
end

@testset "D53 review M8: a cost-cap pause is on disk before the first chunk of its group" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    line = 100.0 * 2 * 1e-9
    cost(factor) = info -> info.steps * line * (("main" in info.arms) ? factor : 1.0)
    main_ords = [c.ordinal for c in plan.chunks if "main" in c.arms]
    store = ex_store("rev-m8")
    man = store * "-manifest.json"
    hooks = EX.RunHooks(clock = ex_clock(cost(2.01)), kill = ex_kill_at(:before_data, main_ords[1]))
    @test ex_killed(() -> ex_campaign(fx; store = store, manifest = man, manifest_reset = false, threads = 1, hooks = hooks))
    @test !isfile(joinpath(store, "chunk-$(lpad(main_ords[1], 5, '0')).h5"))
    @test EX.read_json(man)["arms"]["main"]["status"] == "paused: cost cap"
    res = ex_campaign(fx; store = store, manifest = man, manifest_reset = false, threads = 1, hooks = EX.RunHooks(clock = ex_clock(cost(1.0))))
    @test res.paused == ["main"] && isempty(intersect(res.ran, main_ords))   # the resumed run does not start the paused arm again
end

@testset "D53 review M9: chunk files record the dirty tree and the session" begin
    fx = ex_fix(:small)
    repo = ex_scratch_repo()
    store = ex_store("rev-m9")
    man = store * "-manifest.json"
    run1() = ex_campaign(fx; store = store, manifest = man, manifest_reset = false, repo = repo, allow_dirty = false,
                         hooks = rev_hooks(kill = ex_kill_after(2)))
    @test ex_killed(run1)
    write(joinpath(repo, "tracked.txt"), "changed\n")
    res = ex_campaign(fx; store = store, manifest = man, manifest_reset = false, repo = repo, allow_dirty = true, hooks = rev_hooks())
    @test res.status == "complete"
    for f in ex_chunk_files(store)
        a = EX.read_chunk_attrs(f)
        first_session = ex_ordinal(f) <= 2
        @test a["session"] == (first_session ? 1 : 2) && a["dirty"] === !first_session
    end
    m = EX.read_json(man)
    @test [s["index"] for s in m["sessions"]] == [1, 2] && [s["clean"] for s in m["sessions"]] == [true, false]
    @test "dirty" in ex_schema_attrs() && "session" in ex_schema_attrs() && "control_stop" in ex_schema_attrs()
end

@testset "D53 review LOW: device batch of zero orbits, budget line, sha, stop fraction, cap names, L6" begin
    # csize = 0 is the out-of-memory fallback, not an error
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    rt = first(ex_dev_types(plan))
    g = ex_with_device(() -> EX.device_records(ctx, ex_dev_items(plan, rt)), ex_hdev(; free = 0))
    @test g.recs === nothing && g.prov.backend === :cpu && g.prov.reason == Chaoslab.FALLBACK_REASONS[4]
    # the budget line is the one of the group's precision, not of the arm's first
    run = ex_rb_with(ex_rb(:small), ("arm", 1, "precisions") => ["Double64", "Float64"], ("arm", 1, "orbits") => "S_D", ("arm", 1, "methods") => ["L1"],
                     ("arm", 1, "horizons") => [100])
    run["arm"] = run["arm"][1:1]
    fxm = ex_fix_with(run, "rev-low-line")
    res = ex_campaign(fxm; store = ex_store("rev-low-line"), threads = 1,
                      hooks = EX.RunHooks(clock = ex_clock(info -> info.steps * 1000e-9 * (info.precision == "Float64" ? 1 : 1e-3); start = 0.0)))
    @test res.status == "paused arms" && res.paused == ["main"] && any(r -> occursin("budget line 200.0", r), res.reasons)
    # session hashes come from the plan
    m = EX.read_json(ex_small_run("ref").res.manifest)
    p = ex_plan(fx)
    @test m["prereg"]["sha256_start"] == p.prereg_header && m["panel"]["sha256_start"] == p.panel.sha256
    # the failed fraction is checked at the start of a resumed session
    f1 = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "stop_failed_fraction") => 0.5), "rev-low-stop-a")
    bad(orbits, info) = (first(values(orbits[1].lyap))[1, 1] = NaN)
    store = ex_store("rev-low-stop")
    r1 = ex_campaign(f1; store = store, hooks = rev_hooks(mutate = ex_fault(bad; ord = 3, attempts = (1, 2))))
    @test r1.status == "complete with failed chunks" && r1.failed == [3]
    f2 = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "stop_failed_fraction") => 0.01), "rev-low-stop-b")
    r2 = ex_campaign(f2; store = store, hooks = rev_hooks())
    @test r2.status == "stopped: failed chunks" && isempty(r2.ran) && any(r -> occursin("failed twice", r), r2.reasons)
    # unknown cap override names are refused before anything is written
    st = ex_store("rev-low-cap")
    @test ex_rejects(() -> ex_campaign(fx; store = st, cap_override = ["main", "nosuch"]), ["cap_override", "nosuch"]) && !ispath(st)
    io = IOBuffer()
    code = EX.cli(["run", fx.prereg, fx.panel, "--store", st, "--manifest", st * ".json", "--cap-override", "nosuch", "--thread-override", "--allow-dirty"]; io = io)
    @test code == 1 && occursin("nosuch", String(take!(io)))
    pos, opt = EX._cli_parse(["a", "--retry-control-stop", "--backend-override", "--manifest-reset"])
    @test pos == ["a"] && all(k -> opt[k] === true, ("--retry-control-stop", "--backend-override", "--manifest-reset"))
    # L6 has no extended-precision path: the plan refuses it
    l6 = ex_rb_with(ex_rb(:small), ("arm", 1, "methods") => ["L6"], ("arm", 1, "precisions") => ["Float64", "Double64"])
    @test ex_rejects(() -> ex_plan(ex_fix_with(l6, "rev-low-l6")), ["L6", "Float64 only"])
end
