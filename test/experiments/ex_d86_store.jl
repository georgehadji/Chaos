# Amendment D86, criteria D86s (store layout), D86t (atomic write and resume), D86u (manifest), D86v (pilot and rerun).
# Spec: experiments-spec.md, "Amendment D86": R15, R19, R20, I14, I18, I24. API assumed: experiments-test-plan.md, "API assumed by the D86 tests".

d86_s_run() = d86_run("s-all", d86_all_block())                       # all six kinds, and the arm "mix" with the methods L1 and L5 (old files)

@testset "D86s store layout: attributes, groups, dtypes, shapes, attribute lists" begin
    run = d86_s_run()
    plan = run.plan
    @test run.res.status == "complete"
    seen = Set{String}()
    bad = Dict{String, Int}()
    miss!(tag, ok) = ok || (bad[tag] = get(bad, tag, 0) + 1)
    pooled_seen = Ref(0)
    for file in ex_chunk_files(run.store)
        c = plan.chunks[ex_ordinal(file)]
        it0 = plan.items[first(c.items)]
        raw = d86_raw(file)
        lay = EX.chunk_layout(file)
        miss!("schema version", raw["@schema_version"] == 2 && raw["@complete"] === true)
        if !(it0.method in Set(s["id"] for s in plan.block.raw["indicator"]))
            # an old method: the version 2 file still reads, its lyap group is there
            ch = EX.read_chunk(file)
            miss!("old file", !isempty(first(ch.orbits).lyap) && !haskey(raw, "@kind"))
            continue
        end
        spec = d86_spec_of(plan, it0)
        kind = spec["kind"]
        push!(seen, kind)
        miss!("root", raw["@kind"] == kind && raw["@method"] == it0.method && raw["@vectors"] == it0.vectors && raw["@spec_sha256"] == d86_spec_hash(spec) && raw["@spec"] == d86_spec_text(spec))
        miss!("members/chain", (kind in ("tangent", "pair") ? raw["@members"] == spec["members"] && length(raw["@chain"]) == length(spec["members"]) : haskey(raw, "@chain")))
        miss!("no perturbation", !haskey(raw, "@perturbation"))
        miss!("direction attribute", haskey(raw, "@direction") == (it0.direction !== nothing))
        d86_window_kind(spec) && miss!("windows", raw["@windows"] == spec["windows"] && raw["@run_length"] == d86_runlen(spec))
        kind == "pair" && miss!("pair root", raw["@frame_source"] == spec["frame_source"] && raw["@direction_rule"] == "u = xi_0" && raw["@separation"] == spec["separation"])
        kind == "birkhoff" && miss!("observable_names", length(raw["@observable_names"]) == length(spec["observables"]))
        for i in c.items
            it = plan.items[i]; o = it.orbit; id = it.method
            n = length(raw["orbits/$o/t"])
            L(name) = lay["/orbits/$o/$name"]
            miss!("t", L("t").eltype == "Int64" && L("t").size == (n,))
            miss!("state", L("state").eltype == "Float64" && L("state").size == (n, 2))
            it.coordinate_mode == "wrapped" ? miss!("winding", L("winding").eltype == "Int64" && L("winding").size == (n, 2)) : miss!("no winding", !haskey(lay, "/orbits/$o/winding"))
            miss!("monitor", L("diagnostics/det").size == (n,) && L("diagnostics/symp").size == (n,) && !any(k -> startswith(k, "/orbits/$o/diagnostics/") && !(k in ("/orbits/$o/diagnostics/det", "/orbits/$o/diagnostics/symp")), keys(lay)))
            miss!("no lyap", !any(k -> startswith(k, "/orbits/$o/lyap"), keys(lay)))
            miss!("cost", L("cost/$id").eltype == "Float64" && L("cost/$id").size == (kind == "frequency" ? 6 : 5,))
            if kind in ("tangent", "pair", "stretching")
                ser = L("indicators/$id/series")
                cols = d86_at(raw, o, id, "series", "columns")
                miss!("series", ser.eltype == "Float64" && ser.size == (n, length(cols)))
                miss!("no table", !haskey(lay, "/orbits/$o/indicators/$id/final") && !haskey(lay, "/orbits/$o/indicators/$id/flags"))
            else
                fin, flg = L("indicators/$id/final"), L("indicators/$id/flags")
                rows = length(d86_at(raw, o, id, "final", "rows"))
                miss!("final", fin.eltype == "Float64" && fin.size[1] == rows && fin.size[2] == length(d86_at(raw, o, id, "final", "columns")))
                miss!("flags", flg.eltype == "Int64" && flg.size == (rows,) && !isempty(d86_at(raw, o, id, "flags", "bits")))
                miss!("no series", !haskey(lay, "/orbits/$o/indicators/$id/series"))
                # E25: orbit scalars are `handoffs`, `near_pi`, `complete` (phase advance) and `handoffs`, `complete` (birkhoff); for `frequency` the table says "none (planes = 1)"
                if kind != "frequency"
                    sc = L("indicators/$id/scalars")
                    miss!("scalars", sc.eltype == "Float64" && sc.size == (length(d86_at(raw, o, id, "scalars", "columns")),))
                end
                kind == "phase_advance" && miss!("pa scalars", d86_at(raw, o, id, "scalars", "columns") == ["handoffs", "near_pi", "complete"])
                kind == "birkhoff" && miss!("bk scalars", d86_at(raw, o, id, "scalars", "columns") == ["handoffs", "complete"] && fin.size[1] == length(spec["windows"]) * length(spec["observables"]))
                kind == "frequency" && miss!("fq rows", fin.size == (length(spec["windows"]), 11) && !haskey(lay, "/orbits/$o/indicators/$id/scalars"))
            end
            if kind == "stretching"
                miss!("st scalars", L("indicators/$id/scalars").size == (18,))
                has = haskey(lay, "/orbits/$o/indicators/$id/counts")
                miss!("counts only for the subset", has == plan.panel.orbits[o]["in_SL"])
                has && miss!("counts layout", L("indicators/$id/counts").eltype == "UInt32" && L("indicators/$id/counts").size == (spec["bins"] + 2, 2))
            end
            if kind == "pair"
                miss!("pair orbit attrs", all(k -> haskey(raw, "orbits/$o@$k"), ("Dx", "Dx_r", "mode", "xi0")) && length(raw["orbits/$o@xi0"]) == 2)
            end
        end
        if kind == "stretching"
            miss!("pooled orbits", haskey(raw, "pooled/$(it0.method)@orbits"))
            if haskey(lay, "/pooled/$(it0.method)/counts")                    # (D88) only for a chunk whose orbits share one range
                pl = lay["/pooled/$(it0.method)/counts"]
                miss!("pooled", pl.eltype == "Int64" && pl.size == (spec["bins"] + 2, 2) && get(raw, "pooled/$(it0.method)@range", nothing) isa Vector{Float64})
                pooled_seen[] += 1
            end
        end
    end
    @test isempty(bad)
    @test seen == Set(D86_KINDS) && pooled_seen[] > 0
    # the listing of the store layout and the package read-back equal what was written
    f = first(f for f in ex_chunk_files(run.store) if plan.items[plan.chunks[ex_ordinal(f)].items[1]].method == "tg")
    ch = EX.read_chunk(f)
    out = ex_path("d86-s-roundtrip.h5")
    EX.write_chunk(out, ch.attrs, ch.orbits)
    @test d86_canon_file(out) == d86_canon_file(f)
    # sizes (derived): the series of a three-member tangent item is 8 x grid steps x columns bytes per orbit, a frequency item adds 96 bytes per window
    raw = d86_raw(f)
    o = ch.orbits[1].id
    @test sizeof(d86_ds(raw, o, "tg", "series")) == 8 * length(raw["orbits/$o/t"]) * length(d86_at(raw, o, "tg", "series", "columns"))
end

# ================================================================ D86t
@testset "D86t atomic write killed at 20 seeded points; resume after k chunks equals the uninterrupted campaign" begin
    # The fixture of the kill points: all six kinds, Float64 and Double64, the direction sweep, pooled counts, on the orbit 1 (S_L) and the two extra starts, about 23 chunks
    # instead of 71 (the kill points and resumes of this testset are 31 campaigns; the size of the panel is not what the criterion is about)
    run = d86_all_ind_block()
    foreach(a -> a["orbits"] == "S_D" && (a["orbits"] = "S_L"), run["arm"])
    ref = d86_run("t-ref", run)
    fx = ref.fx
    M = ref.plan.counts.chunks
    refsha = d86_store_sha(ref.store)
    points = (:before_data, :during_data, :after_data, :before_rename, :after_rename)
    rng = StableRNG(d86_seed('t', 1))
    final(store, o) = joinpath(store, "chunk-$(lpad(o, 5, '0')).h5")
    for t in 1:20
        ord, pt = rand(rng, 1:M), points[1 + (t - 1) % 5]
        store = ex_store("d86-t-$t")
        @test ex_killed(() -> ex_campaign(fx; store = store, hooks = EX.RunHooks(kill = ex_kill_at(pt, ord))))
        @test all(f -> EX.read_chunk(f).attrs["complete"] === true, ex_chunk_files(store))          # no incomplete file under a final name
        @test isfile(final(store, ord)) == (pt === :after_rename)
        res = ex_campaign(fx; store = store)
        @test res.status == "complete" && isempty(ex_tmp_files(store))
        @test d86_store_sha(store) == refsha
    end
    # killed after k chunks (10 seeded k), resumed with another thread count: every dataset and attribute named in I14, pooled counts included
    ks = sortperm(rand(rng, M - 1))[1:10]
    for (j, k) in enumerate(ks)
        store = ex_store("d86-t-k$j")
        @test ex_killed(() -> ex_campaign(fx; store = store, threads = 4, hooks = EX.RunHooks(kill = ex_kill_after(k))))
        res = ex_campaign(fx; store = store, threads = 1 + j % 2)
        @test res.status == "complete" && sort(res.ran) == collect((k + 1):M)
        @test d86_store_sha(store) == refsha
        @test all(f -> d86_same(get(d86_raw(f), "pooled/st/counts", nothing), get(d86_raw(joinpath(ref.store, basename(f))), "pooled/st/counts", nothing)), ex_chunk_files(store))
    end
    # a changed spec parameter changes the item and chunk hashes: the old files are reported stale and the run stops
    changed = d86_ed(run, "st", "bins" => 21)
    store = ex_store("d86-t-stale"); cp(ref.store, store)
    err = try
        ex_campaign(d86_fx("t-changed", changed); store = store); nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("stale", err.msg)
    # a final failed pair chunk is not retried on resume
    pr = first(c for c in ref.plan.chunks if ref.plan.items[first(c.items)].method == "pr")
    bad(orbits, info) = (orbits[1].indicators["pr"].series[2, 1] = NaN)
    store = ex_store("d86-t-pair")
    res = ex_campaign(fx; store = store, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), mutate = ex_fault(bad; ord = pr.ordinal, attempts = (1, 2))))
    @test pr.ordinal in res.failed && res.status in ("complete with failed chunks", "stopped: failed chunks")
    fpr = joinpath(store, "chunk-$(lpad(pr.ordinal, 5, '0')).h5")
    t0 = mtime(fpr)
    res2 = ex_campaign(fx; store = store)
    @test pr.ordinal in res2.skipped && !(pr.ordinal in res2.ran)
    @test EX.read_chunk(fpr).attrs["status"] == "failed" && mtime(fpr) == t0
end

# ================================================================ D86u
@testset "D86u manifest: every field of R19, the chain records, the omitted lists, the flag counts, the cost lines" begin
    run = d86_s_run()
    plan = run.plan
    man = EX.read_json(run.res.manifest)
    ind = man["indicators"]
    @test Set(keys(ind)) == Set(s["id"] for s in plan.block.raw["indicator"])
    for s in plan.block.raw["indicator"]
        m = ind[s["id"]]
        it = first(i for i in plan.items if i.method == s["id"])
        @test m["kind"] == s["kind"] && m["spec"] == d86_spec_text(s) && m["spec_sha256"] == d86_spec_hash(s) && m["budget"] == s["budget"]
        if s["kind"] in ("tangent", "pair")
            f = first(f for f in ex_chunk_files(run.store) if plan.items[plan.chunks[ex_ordinal(f)].items[1]].method == s["id"])
            @test m["chain"] == d86_raw(f)["@chain"]                                 # the chain records equal those of the files
        end
        d86_window_kind(s) && @test m["windows"] == s["windows"] && m["run_length"] == d86_runlen(s)
        s["kind"] == "pair" && @test m["frame_source"] == s["frame_source"] && m["direction_rule"] == "u = xi_0"
        s["kind"] == "birkhoff" && @test length(m["observable_names"]) == length(s["observables"])
        b = d86_ref_bytes(plan, it)
        @test m["bytes_per_orbit"] == Dict("M" => b.M, "R" => b.R, "S" => b.S, "B" => b.B, "total" => b.total)
    end
    # backend table rows carry the spec ids; the arm totals carry the measured cost per orbit-step against the budget line
    @test all(id -> any(r -> r["method"] == id, man["backend_table"]), keys(ind))
    @test all(a -> haskey(a, "specs") && all(v -> haskey(v, "measured_ns_per_orbit_step") && haskey(v, "budget_line_ns"), values(a["specs"])), values(man["arms"]))
    # chunk entries: omitted orbits, flag counts recounted from the files, transform seconds for frequency
    omitted = Dict{String, Set{Int}}()
    for ce in man["chunks"]
        ce["status"] == "ok" || continue
        it0 = plan.items[plan.chunks[ce["ordinal"]].items[1]]
        it0.method in keys(ind) || continue
        for (id, os) in get(ce, "not_applicable", Dict())
            union!(get!(omitted, id, Set{Int}()), os)
        end
        raw = d86_raw(joinpath(run.store, "chunk-$(lpad(ce["ordinal"], 5, '0')).h5"))
        count = Dict{String, Int}()
        for i in plan.chunks[ce["ordinal"]].items
            o = plan.items[i].orbit
            haskey(raw, "orbits/$o/indicators/$(it0.method)/flags") || continue
            names = d86_flagset.(raw["orbits/$o/indicators/$(it0.method)/flags"], Ref(raw["orbits/$o/indicators/$(it0.method)/flags@bits"]))
            for nm in reduce(union, names; init = Set{String}())
                count[nm] = get(count, nm, 0) + 1
            end
        end
        got = get(ce, "indicator_flags", Dict())
        @test Dict(k => v for (k, v) in got if v > 0) == count
        ind[it0.method]["kind"] == "frequency" && @test haskey(ce, "transform_seconds")
    end
    @test omitted == Dict(id => Set(os) for (id, os) in plan.not_applicable)
    # a campaign without specs has the keys of R7: no `indicators` key, and no other key is added or lost
    old = EX.read_json(d86_run("b-small", ex_fix(:small).run).res.manifest)
    @test !haskey(old, "indicators") && Set(keys(man)) == union(Set(keys(old)), Set(["indicators"]))
end

# ================================================================ D86v
@testset "D86v pilot at the own run length of the window kinds, the cost cap on the spec line, the CPU reference rerun" begin
    fx = d86_fx("v-all", d86_all_ind_block())
    plan = ex_plan(fx)
    out = ex_path("d86-v-pilot"); rm(out; recursive = true, force = true)
    default_before = isdir(EX.default_store(fx.prereg))
    res = EX.run_pilot(fx.prereg, fx.panel; out = out, threads = 1, thread_override = true, allow_dirty = true, code_version = EX_CODE, capacity = 64,
                       hooks = EX.RunHooks(clock = ex_clock(info -> info.steps * 150e-9; start = 0.0)))
    rep = res.report
    @test occursin("not study data", rep["not_study_data"])
    ids = Set(s["id"] for s in plan.block.raw["indicator"])
    rows = [c for c in rep["cost"] if c["method"] in ids]
    @test Set(c["method"] for c in rows) == ids
    pilot_h = fx.run["pilot"]["horizon"]
    for c in rows
        s = only(s for s in plan.block.raw["indicator"] if s["id"] == c["method"])
        @test all(k -> haskey(c, k), ("thread_seconds_per_orbit_step", "wall_seconds_per_orbit_step", "bytes_per_orbit", "capacity", "flag_shares", "horizon"))
        @test c["horizon"] == (d86_window_kind(s) ? d86_runlen(s) : pilot_h)              # window kinds at their own run length, tau = 1 in the fixture
        s["kind"] == "frequency" && @test haskey(c, "orbit_pass_seconds") && haskey(c, "transform_seconds")
    end
    pfiles = ex_chunk_files(joinpath(out, "store"))
    @test !isempty(pfiles) && all(f -> EX.read_chunk(f).attrs["pilot"] === true, pfiles) && isdir(EX.default_store(fx.prereg)) == default_before
    # the first-chunk cost cap reads the budget line of the spec (safety 3: 300 ns, not the 200 ns of the block): 2.01 times the line pauses the arm, 1.99 continues
    # The arm a has two run types (horizons 100 and 80), so two launch groups: the first group sets the cost, the second must stay pending when the arm is paused (D53l form).
    blk = d86_block([d86_arm("a", "S_D", ["tg"], ["Float64"], [100, 80]), d86_arm("b", "S_D", ["st"], ["Float64"], [100])],
                    [merge(d86_spec("tangent"), Dict("budget" => Dict{String, Any}("Float64" => Dict{String, Any}("ns_cpu_thread" => 100.0, "safety" => 3)))), d86_spec("stretching")])
    fc = d86_fx("v-cap", blk)
    line(arm) = (arm == "a" ? 300.0 : 200.0) * 1e-9
    cost(f) = info -> info.steps * (("a" in info.arms) ? line("a") * f : line("b"))
    hi = ex_campaign(fc; store = ex_store("d86-v-hi"), threads = 1, hooks = EX.RunHooks(clock = ex_clock(cost(2.01))))
    lo = ex_campaign(fc; store = ex_store("d86-v-lo"), threads = 1, hooks = EX.RunHooks(clock = ex_clock(cost(1.99))))
    @test hi.status == "paused arms" && hi.paused == ["a"]
    manh = EX.read_json(hi.manifest)
    # D53l form: the first group of the arm ran, its later chunks are recorded as paused and were not run, the other arm completed
    paused = [c for c in manh["chunks"] if c["status"] == "paused"]
    @test !isempty(paused) && all(c -> c["arms"] == ["a"] && !(c["ordinal"] in hi.ran), paused)
    @test manh["arms"]["a"]["status"] == "paused: cost cap" && manh["arms"]["b"]["status"] == "complete"
    @test any(c -> c["status"] == "ok" && c["arms"] == ["a"], manh["chunks"])
    manl = EX.read_json(lo.manifest)
    @test lo.status == "complete" && isempty(lo.paused) && all(c -> c["status"] == "ok", manl["chunks"])
    # the CPU reference rerun reproduces CPU items bit for bit
    run = d86_all_ind_run()
    chk = EX.rerun_check(run.fx.prereg, run.fx.panel, run.store; code_version = EX_CODE, threads = 2, thread_override = true, allow_dirty = true)
    ci = [r for r in chk if run.plan.items[r.item].method in ids]
    @test !isempty(ci) && all(r -> r.bitwise && r.backend == "cpu", ci)
    @test Set(run.plan.items[r.item].method for r in chk) ⊇ Set(["tg", "fq", "bk"])           # items of several kinds are in the sample
    # a simulated-device item is reported with its difference and does not fail the check
    shas = Dict(s["id"] => d86_spec_hash(s) for s in plan.block.raw["indicator"])
    rec = (threads = 4, cpu_ns_per_step = 8.25, gpu_ns_per_step = 2.47, gpu_faster = true, date = "2026-10-28",
           recorder = Dict{String, Any}(id => (threads = 4, gpu_faster = true, spec_sha256 = sha) for (id, sha) in shas))
    tg1 = first(c for c in plan.chunks if plan.items[first(c.items)].method == "tg")
    shift!(orbits, info) = (orbits[1].indicators["tg"].series[end, 1] += 1e-6)
    store = ex_store("d86-v-dev")
    ex_with_device(() -> ex_campaign(fx; store = store, device_record = rec,
                   hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), mutate = ex_fault(shift!; ord = tg1.ordinal))), ex_hdev())
    dchk = EX.rerun_check(fx.prereg, fx.panel, store; code_version = EX_CODE, threads = 4, thread_override = true, allow_dirty = true, sample = 100)
    drow = [r for r in dchk if r.backend == "gpu" && plan.items[r.item].method == "tg"]
    @test !isempty(drow) && any(r -> !r.bitwise && r.max_abs_diff > 0 && r.gated === false, drow)
end
