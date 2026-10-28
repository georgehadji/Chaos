# Amendment D86, criteria D86n (chunks and memory), D86o (backend, simulated device), D86p (device path on the KernelAbstractions CPU backend),
# D86q (determinism), D86r (health). Spec: experiments-spec.md, "Amendment D86": E28, R16, R17, R18, I24, I25. API assumed: experiments-test-plan.md.

# the capacity a campaign reads from the free memory it is given (manifest, launch.capacity); stopped after the first chunk, which is all the manifest needs
function d86_capacity(fx, free)
    store = ex_store("d86-cap-$(fx.tag)-$free")
    ex_killed(() -> ex_campaign(fx; store = store, capacity = nothing, free_bytes = free, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), kill = ex_kill_after(1))))
    return EX.read_json(store * "-manifest.json")["launch"]["capacity"]
end

function d86_cap_grid(; ppd = 20, first = 100, fr = [0.1, 0.5])
    return Dict{String, Any}("points_per_decade" => ppd, "first_steps" => first, "extra_fractions" => fr)
end
function d86_worked(spec; N, mode = "wrapped", precisions = ["Float64"], orbits = "S_L", grid = d86_cap_grid(), tag)
    blk = d86_one(spec; orbits = orbits, precisions = precisions, horizon = N, mode = mode)
    blk["grid"] = grid
    blk["pilot"]["horizon"] = 40
    return d86_fx("n-" * tag, blk)
end

@testset "D86n sizes of the reducers and the worked numbers of R17" begin
    T = Float64
    rb(r, d = 2, k = 2) = ex_Tg.reducer_bytes(r, T, Val(d), Val(k))
    @test rb(Chaoslab.SALI()) == 72 && rb(Chaoslab.FLI()) == 96 && rb(Chaoslab.MEGNO()) == 200 && rb(Chaoslab.RLI(), 4, 2) == 88
    st = Chaoslab.StretchingNumbers(range = (-1.0, 1.0), bins = 1000, interval = 100, min_per_cell = 1)
    @test sizeof(ex_Tg.reducer_init(st, T, Val(2), Val(2))) == 104 && 4 * ex_Tg.reducer_workspace(st, T, Val(2), Val(2)) == 8016
    wp(H) = Chaoslab.WeightedPhaseAdvance(window = ntuple(i -> 16 + i, H), observable = Chaoslab.PlanarObservable(1, 2))
    @test [rb(wp(H), 2, 1) for H in 1:4] == [96, 160, 224, 288]
    for (K, H) in ((4, 3), (31, 1), (4, 11), (1, 1))
        wb = Chaoslab.WeightedBirkhoff(window = ntuple(i -> 16 + i, H), observables = ntuple(i -> Chaoslab.CosObservable(1, i), K))
        @test rb(wb, 2, 1) == 16 + 8 * H * (2 + 5 * K)
    end
    @test 3 * sizeof(T) == 24
    @test ex_Tg.bytes_per_orbit(Val(2), Val(2), T, 4) == 112 && ex_Tg.bytes_per_orbit(Val(2), Val(1), T, 4) == 80 && ex_Tg.bytes_per_orbit(Val(4), Val(2), T, 8) == 176
    # R17 item 4, 100 000 orbits are not needed: bytes(orbit) is per orbit and the capacity is read from a free memory of 16 GB, f = 1/2
    # (a) tangent FLI, MEGNO, SALI at N = 4000, wrapped: 135 grid steps, c = 17, bytes 18 864, capacity 424 088
    fa = d86_worked(d86_spec("tangent"); N = 4000, tag = "a")
    pa = ex_plan(fa)
    ba = d86_ref_bytes(pa, first(pa.items))
    @test ba.ng == 135 && ba.c == 17 && ba.M == 112 && ba.R == 392 && ba.S == 18360 && ba.total == 18864
    for (free, cap) in ((16 * 10^9, 424088), (4 * 10^9, 106022))
        @test d86_ref_capacity(pa, free) == cap && d86_capacity(fa, free) == cap
    end
    # (b) frequency, windows 250, 1000, 4000 (run length 8000, 141 grid steps): bytes 136 000, capacity 58 823; Double64 B = 256 000
    fb = d86_worked(d86_spec("frequency"; windows = [250, 1000, 4000]); N = 8000, precisions = ["Float64", "Double64"], tag = "b")
    pb = ex_plan(fb)
    bb = d86_ref_bytes(pb, first(it for it in pb.items if it.precision == "Float64"))
    @test bb.ng == 141 && bb.M == 80 && bb.R == 24 && bb.S == 7896 && bb.B == 128000 && bb.total == 136000
    @test d86_ref_bytes(pb, first(it for it in pb.items if it.precision == "Double64")).B == 256000
    @test d86_capacity(d86_worked(d86_spec("frequency"; windows = [250, 1000, 4000]); N = 8000, tag = "b64"), 16 * 10^9) == 58823
    # (c) stretching, 1000 bins, N = 4000: R = 8144, S = 8640, bytes 16 896, capacity 473 484
    fc = d86_worked(d86_spec("stretching"; bins = 1000, interval = 100); N = 4000, tag = "c")
    pc = ex_plan(fc)
    bc = d86_ref_bytes(pc, first(pc.items))
    @test bc.R == 8144 && bc.S == 8640 && bc.total == 16896 && d86_ref_capacity(pc, 16 * 10^9) == 473484 && d86_capacity(fc, 16 * 10^9) == 473484
    # (d) weighted Birkhoff K = 4, H = 3: R = 568, bytes 8 544, capacity 936 329; phase advance with three windows: R = 248
    fd = d86_worked(d86_spec("birkhoff"; windows = [250, 1000, 4000]); N = 8000, tag = "d")
    pd = ex_plan(fd)
    bd = d86_ref_bytes(pd, first(pd.items))
    @test bd.R == 568 && bd.total == 8544 && d86_capacity(fd, 16 * 10^9) == 936329
    # (d) phase advance with three windows: run length 8001 with the fraction 0.1 (an indicator-only arm is exempt from the 0.1 rule, correction 1): 141 grid steps,
    # R = 96 + 128 + 24 = 248, bytes 80 + 248 + 7 896 = 8 224, capacity 972 762
    fe = d86_worked(d86_spec("phase_advance"; windows = [250, 1000, 4000]); N = 8001, tag = "e")
    pe = ex_plan(fe)
    be = d86_ref_bytes(pe, first(pe.items))
    @test be.ng == 141 && be.R == 248 && be.S == 7896 && be.total == 8224
    @test d86_ref_capacity(pe, 16 * 10^9) == 972762 && d86_capacity(fe, 16 * 10^9) == 972762
end

@testset "D86n capacity equals the formula for every kind; chunk keys and stored data do not depend on it; one orbit that does not fit refuses the start" begin
    run = d86_all_ind_block()                                                        # indicator arms only: the earlier device path assumes 256 bytes for the methods
    fx = d86_fx("n-all", run)
    plan = ex_plan(fx)
    keys0 = [(c.key, c.hash) for c in plan.chunks]
    for free in (2 * 10^8, 7 * 10^8, 3 * 10^9)
        @test d86_capacity(fx, free) == d86_ref_capacity(plan, free)
        @test [(c.key, c.hash) for c in ex_plan(fx).chunks] == keys0                  # ids and keys do not depend on the free memory (D53g extended)
    end
    # per kind: the item with the largest bytes(orbit) in each run type matches the independent terms (the capacity of one kind alone)
    for id in ("tg", "st", "pr", "prd", "fq", "pa", "bk")
        sub = deepcopy(run)
        sub["indicator"] = [s for s in sub["indicator"] if s["id"] == id]
        sub["arm"] = [merge(a, Dict("indicators" => [id])) for a in sub["arm"] if id in a["indicators"]]
        isempty(sub["arm"]) && continue
        subfx = d86_fx("n-kind-$id", sub)
        p = ex_plan(subfx)
        @test d86_capacity(subfx, 5 * 10^8) == d86_ref_capacity(p, 5 * 10^8)
    end
    # a spec whose one orbit does not fit: the start is refused naming the spec, bytes(orbit) and the free memory read; nothing is written
    store = ex_store("d86-n-nofit")
    big = d86_fx("n-nofit", d86_one(d86_spec("frequency"; windows = [250]); mode = "wrapped"))
    pbig = ex_plan(big)
    need = d86_ref_bytes(pbig, first(pbig.items)).total
    r = try
        ex_campaign(big; store = store, capacity = nothing, free_bytes = need)       # f x free = need / 2 < bytes(orbit)
        nothing
    catch e
        e
    end
    @test r isa ArgumentError && occursin("fq", r.msg) && occursin(string(need), r.msg) && occursin(string(need), r.msg) && occursin("free", r.msg)
    @test isempty(ex_chunk_files(store))
    # stored datasets are bit-identical for the launch capacities 1, 7 and all orbits, and for chunks run in reverse order
    small = d86_one(d86_spec("tangent"), d86_spec("stretching"); orbits = "S_D")
    a = d86_run("n-c1", small; capacity = 1, order = :reverse)
    b = d86_run("n-c7", small; capacity = 7)
    c = d86_run("n-call", small; capacity = 10^6)
    @test a.res.status == b.res.status == c.res.status == "complete"
    @test d86_store_sha(a.store) == d86_store_sha(b.store) == d86_store_sha(c.store)
end

# ================================================================ D86o
@testset "D86o backend matrix for the six kinds (simulated device)" begin
    R = Chaoslab.FALLBACK_REASONS
    plan = d86_plan("c-all", d86_all_block())
    shas = Dict(s["id"] => d86_spec_hash(s) for s in plan.block.raw["indicator"])
    kinds = Dict(s["id"] => s["kind"] for s in plan.block.raw["indicator"])
    threads = 4
    function expect(req, T, kind, dev, rec, sha)
        req === :cpu && return (:cpu, "")
        T === BigFloat && return (:cpu, R[3])
        kind in ("stretching", "frequency") && return (:cpu, R[6])
        dev.loadable || return (:cpu, R[1])
        dev.functional || return (:cpu, R[2])
        dev.compiles || return (:cpu, R[3])
        if req === :auto
            ok = rec !== nothing && rec.threads == threads && rec.gpu_faster && rec.spec_sha256 == sha
            ok || return (:cpu, R[7])
        end
        return (:gpu, "")
    end
    devs = [ExDev(false, false, false), ExDev(true, false, false), ExDev(true, true, false), ExDev(true, true, true)]
    mismatches = String[]
    used = Set{Symbol}()
    for (id, kind) in kinds, T in (Float64, Chaoslab.Double64, BigFloat), req in (:cpu, :gpu, :auto), dev in devs, state in 1:6
        sha = shas[id]
        entry = state == 1 ? nothing : state == 2 ? (threads = threads, gpu_faster = true, spec_sha256 = sha) : state == 3 ? (threads = threads, gpu_faster = false, spec_sha256 = sha) :
                state == 4 ? (threads = 8, gpu_faster = true, spec_sha256 = sha) : state == 5 ? (threads = threads, gpu_faster = true, spec_sha256 = "f"^64) : nothing
        record = state == 6 ? nothing : (recorder = entry === nothing ? Dict{String, Any}() : Dict{String, Any}(id => entry),)
        got = ex_with_device(() -> EX.resolve_run_backend(req, T, id; threads = threads, record = record, kind = kind, spec_sha256 = sha), dev)
        want = expect(req, T, kind, dev, entry, sha)
        (got.backend, got.reason) == want || push!(mismatches, "$id $T $req $(dev) rec$state: $((got.backend, got.reason)) vs $want")
        push!(used, got.backend)
        reasons_ok = got.reason == "" || got.reason in R[1:9]
        reasons_ok || push!(mismatches, "tenth reason $(got.reason)")
        # a record taken on another chain: the existing seventh reason and the note
        if req === :auto && state == 5 && T !== BigFloat && !(kind in ("stretching", "frequency")) && dev.loadable && dev.functional && dev.compiles
            occursin("another chain", got.note) || push!(mismatches, "note $id")
        end
    end
    @test isempty(mismatches)
    @test used == Set([:cpu, :gpu])
    # the backend table of the plan has one row per (precision label, spec id)
    rec = (recorder = Dict{String, Any}(id => (threads = threads, gpu_faster = true, spec_sha256 = shas[id]) for id in keys(shas)),)
    tbl = ex_with_device(() -> EX.backend_table(plan; threads = threads, record = rec), ExDev(true, true, true))
    @test Set(m for ((p, m), r) in tbl if p == "Float64" && m in keys(shas)) == Set(keys(shas))
    @test all(r -> r.backend == :cpu && r.reason == R[6], [tbl[("Float64", id)] for id in ("st", "fq")])
    @test all(r -> r.backend == :gpu, [tbl[("Float64", id)] for id in ("tg", "pr", "prd", "pa", "bk")])
end

# ================================================================ D86p
d86_rec_same(a, b) = ex_rec_same(a, b) && isequal(a.indicators, b.indicators)

@testset "D86p device path: records equal the CPU loop on the KernelAbstractions CPU backend" begin
    seen = Set{String}()
    bad = String[]
    for plan in (d86_chain_run().plan, d86_win_run().plan)                          # the plans of the two shared campaigns
        ctx = EX.make_ctx(plan)
        for rt in unique(c.run_type for c in plan.chunks)
            (rt.precision == "Float64" && rt.method in keys(Dict(s["id"] => s for s in plan.block.raw["indicator"]))) || continue
            spec = only(s for s in plan.block.raw["indicator"] if s["id"] == rt.method)
            spec["kind"] in ("tangent", "pair", "phase_advance", "birkhoff") || continue
            items = ex_dev_items(plan, rt)
            cpu = EX.compute_records(ctx, items, 1)
            got = ex_with_device(() -> EX.device_records(ctx, items), ex_hdev())
            push!(seen, spec["kind"])
            ok = got.recs !== nothing && got.prov.backend === :gpu && got.prov.reason == "" && length(got.recs) == length(cpu) && all(d86_rec_same.(got.recs, cpu))
            ok || push!(bad, string(rt))
        end
    end
    @test isempty(bad)
    @test seen == Set(["tangent", "pair", "phase_advance", "birkhoff"])
    # stretching and frequency never reach the device recorder (CPU only in this amendment)
    plan = d86_win_run().plan
    ctx = EX.make_ctx(plan)
    fq = ex_dev_items(plan, first(c.run_type for c in plan.chunks if c.run_type.method == "fq" && c.run_type.precision == "Float64"))
    @test ex_with_device(() -> EX.device_records(ctx, fq), ex_hdev()).recs === nothing
end

@testset "D86p injected faults send the group to the CPU with the recorded reason and the data of the CPU; the capacity uses the real reducer bytes" begin
    R = Chaoslab.FALLBACK_REASONS
    blk = d86_one(d86_spec("tangent"; id = "t4"), d86_spec("phase_advance"; id = "pa"); orbits = "S_D")
    blk["arm"] = [d86_arm("a", "S_D", ["t4"], ["Float64"], [100]), d86_arm("b", "S_D", ["pa"], ["Float64"], [65])]
    fx = d86_fx("p-fault", blk)
    plan = ex_plan(fx)
    shas = Dict(s["id"] => d86_spec_hash(s) for s in plan.block.raw["indicator"])
    rec = (threads = 4, cpu_ns_per_step = 8.25, gpu_ns_per_step = 2.47, gpu_faster = true, date = "2026-10-28",
           recorder = Dict{String, Any}(id => (threads = 4, gpu_faster = true, spec_sha256 = sha) for (id, sha) in shas))
    cpu = d86_run("p-cpu", blk)
    ref = d86_store_sha(cpu.store)
    function with_dev(dev, tag)
        store = ex_store("d86-p-$tag")
        res = ex_with_device(() -> ex_campaign(fx; store = store, device_record = rec, hooks = ex_zero_clock()), dev)
        return res, store
    end
    res, store = with_dev(ex_hdev(), "ok")
    @test res.status == "complete" && any(f -> EX.read_chunk(f).attrs["backend"] == "gpu", ex_chunk_files(store))
    @test d86_store_sha(store) == ref                                               # the KA CPU backend gives the bits of the CPU loop
    for (dev, r, tag) in ((ex_hdev((ci, li, c) -> throw(ExOOM())), R[4], "oom"), (ex_hdev((ci, li, c) -> throw(ExCompile())), R[3], "compile"))
        res, store = with_dev(dev, tag)
        att = [EX.read_chunk(f).attrs for f in ex_chunk_files(store)]
        @test res.status == "complete" && all(a -> a["backend"] == "cpu" && a["backend_fallback_reason"] == r, att)
        @test d86_store_sha(store) == ref
    end
    # the capacity of a device batch reads the bytes of the chain: free memory chosen so that 9 orbits fit with a 256 byte reducer allowance and 8 with the 392 bytes of the chain
    tplan = d86_plan("p-t4", d86_one(d86_spec("tangent"); orbits = "S_D"))
    it = first(tplan.items)
    b = d86_ref_bytes(tplan, it)
    f = Rational{BigInt}(ex_Tg.CHUNK_HEADROOM)
    F = 9 * (b.M + b.S + 256) + 600                                                  # f x free
    free = Int(ceil(F / f))
    @test fld(F, b.M + b.S + 256) == 9 && fld(F, b.total) == 8
    ctx = EX.make_ctx(tplan)
    items = ex_dev_items(tplan, first(c.run_type for c in tplan.chunks))
    got = ex_with_device(() -> EX.device_records(ctx, items), ex_hdev(; free = free))
    @test got.recs !== nothing && got.prov.chunk_size == 8
end

# ================================================================ D86q
@testset "D86q determinism: threads, capacity, order, alone and in a group, for every kind in Float64 and Double64" begin
    run = d86_all_ind_block()
    runs = [d86_all_ind_run(), d86_run("q-b", run; threads = 1, capacity = 1, order = :reverse),
            d86_run("q-c", run; threads = Threads.nthreads(), capacity = 7), d86_run("q-d", run; threads = 2, capacity = 64)]
    @test all(r -> r.res.status == "complete", runs)
    shas = [d86_store_sha(r.store) for r in runs]
    @test length(unique(shas)) == 1
    # per file as well, so that a difference is located
    files = [ex_chunk_files(r.store) for r in runs]
    @test all(f -> length(f) == length(files[1]), files)
    for k in 2:4, (f1, fk) in zip(files[1], files[k])
        @test d86_canon_file(f1) == d86_canon_file(fk)
    end
    @test Set(d86_raw(f)["@kind"] for f in files[1]) == Set(D86_KINDS)
    # the manifest records the thread count
    @test [EX.read_json(r.res.manifest)["threads"]["used"] for r in runs] == [4, 1, Threads.nthreads(), 2]
end

# ================================================================ D86r
orbits_id(ch, i) = ch.orbits[i].id
d86_bit(bits, name) = Int64(1) << (findfirst(==(name), bits) - 1)
const D86_FREQ_BITS = ["below_floor", "degenerate", "failed", "tie1", "tie2", "near_zero1", "near_zero2", "mirror1", "mirror2", "refined1", "refined2"]
const D86_PA_BITS = ["below_floor", "undefined window"]
const D86_BK_BITS = ["below_floor", "zero", "reldig_undefined"]

@testset "D86r health: injected cases" begin
    blk = d86_block([d86_arm("a", "S_D", ["t4"], ["Float64"], [100]), d86_arm("f", "S_D", ["fq"], ["Float64"], [64]; coordinate_mode = ["wrapped"]),
                     d86_arm("p", "S_D", ["pa"], ["Float64"], [65]), d86_arm("b", "S_D", ["bk"], ["Float64"], [64]; coordinate_mode = ["wrapped"])],
                    [d86_spec("tangent"; id = "t4"), d86_spec("frequency"; id = "fq"), d86_spec("phase_advance"; id = "pa"), d86_spec("birkhoff"; id = "bk")])
    fx = d86_fx("r-all", blk)
    plan = ex_plan(fx)
    chunk(id) = first(c for c in plan.chunks if plan.items[first(c.items)].method == id)
    function inject(id, f; attempts = (1, 2))
        c = chunk(id)
        store = ex_store("d86-r-$id-$(hash(f))")
        ex_campaign(fx; store = store, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), mutate = ex_fault(f; ord = c.ordinal, attempts = attempts)))
        return EX.read_chunk(joinpath(store, "chunk-$(lpad(c.ordinal, 5, '0')).h5")), store, EX.read_json(store * "-manifest.json")
    end
    ind(o, id) = o.indicators[id]
    # a NaN in a column of the tangent family fails the chunk
    ch, _, _ = inject("t4", (orbits, info) -> ind(orbits[2], "t4").series[3, 2] = NaN)
    @test ch.attrs["status"] == "failed" && any(r -> occursin("non-finite", r) && occursin("orbit $(orbits_id(ch, 2))", r), ch.attrs["failure_reasons"])
    # an Inf in a window table fails it, naming the orbit id and the column
    ch, _, _ = inject("fq", (orbits, info) -> ind(orbits[1], "fq").final[1, 3] = Inf)
    @test ch.attrs["status"] == "failed" && any(r -> occursin("orbit $(ch.orbits[1].id)", r) && occursin("amp1", r), ch.attrs["failure_reasons"])
    # a NaN with its flag passes and is counted per flag in the manifest; the same NaN without the flag fails
    withflag(id, col, name, bits) = (orbits, info) -> (t = ind(orbits[1], id); t.final[1, col] = NaN; t.flags[1] |= d86_bit(bits, name))
    nofl(id, col) = (orbits, info) -> (ind(orbits[1], id).final[1, col] = NaN)
    ch, _, man = inject("fq", withflag("fq", 1, "degenerate", D86_FREQ_BITS))
    @test ch.attrs["status"] == "ok"
    @test only(c for c in man["chunks"] if c["ordinal"] == chunk("fq").ordinal)["indicator_flags"]["degenerate"] >= 1
    @test inject("fq", nofl("fq", 1))[1].attrs["status"] == "failed"
    @test inject("pa", withflag("pa", 1, "undefined window", D86_PA_BITS))[1].attrs["status"] == "ok"
    @test inject("pa", nofl("pa", 1))[1].attrs["status"] == "failed"
    @test inject("bk", withflag("bk", 7, "reldig_undefined", D86_BK_BITS))[1].attrs["status"] == "ok"
    @test inject("bk", nofl("bk", 7))[1].attrs["status"] == "failed"
    zero_(orbits, info) = (t = ind(orbits[1], "bk"); t.final[1, 6:8] .= NaN; t.flags[1] |= d86_bit(D86_BK_BITS, "zero"))
    @test inject("bk", zero_)[1].attrs["status"] == "ok"
    # complete false fails; the frequency status failed fails; a chunk in which every orbit carries below_floor passes
    @test inject("pa", (orbits, info) -> ind(orbits[1], "pa").scalars[3] = 0.0)[1].attrs["status"] == "failed"
    @test inject("fq", (orbits, info) -> ind(orbits[1], "fq").flags[1] |= d86_bit(D86_FREQ_BITS, "failed"))[1].attrs["status"] == "failed"
    allfloor(orbits, info) = for o in orbits
        ind(o, "fq").flags .|= d86_bit(D86_FREQ_BITS, "below_floor")
    end
    @test inject("fq", allfloor)[1].attrs["status"] == "ok"
    # the monitor exceedance rule on an indicator arm: 40 orbits, 2 exceeding (5 percent) pass, 3 fail, the orbits flagged in both cases
    big = d86_one(d86_spec("tangent"); orbits = "P")
    for (n_exceed, want) in ((2, "ok"), (3, "failed"))
        fb = d86_fx("r-big$n_exceed", big; panel = :big)
        exceed!(orbits, info) = for o in orbits[1:n_exceed]
            o.diagnostics["det"][end] = 2 * o.thresholds["det"]
        end
        store = ex_store("d86-r-frac$n_exceed")
        ex_campaign(fb; store = store, hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), mutate = ex_fault(exceed!; ord = 1, attempts = (1, 2))))
        chb = EX.read_chunk(first(ex_chunk_files(store)))
        @test chb.attrs["status"] == want && length(chb.orbits) == 40 && flag_ids(chb, "diagnostics exceeded") == collect(1:n_exceed)
    end
    # the CPU rerun after an injected device fault gives the data of a CPU run
    shas = Dict(s["id"] => d86_spec_hash(s) for s in plan.block.raw["indicator"])
    rec = (threads = 4, cpu_ns_per_step = 8.25, gpu_ns_per_step = 2.47, gpu_faster = true, date = "2026-10-28",
           recorder = Dict{String, Any}(id => (threads = 4, gpu_faster = true, spec_sha256 = sha) for (id, sha) in shas))
    cpu = d86_run("r-cpu", blk)
    c = chunk("t4")
    log = NamedTuple[]
    nan1!(orbits, info) = (ind(orbits[1], "t4").series[2, 1] = NaN)
    store = ex_store("d86-r-rerun")
    hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(nan1!; ord = c.ordinal, log = log))
    res = ex_with_device(() -> ex_campaign(fx; store = store, device_record = rec, hooks = hooks), ex_hdev())
    @test res.status == "complete" && [(l.attempt, l.backend) for l in log if l.ordinal == c.ordinal] == [(1, :gpu), (2, :cpu)]
    @test d86_store_sha(store) == d86_store_sha(cpu.store)
    # the exact-exponent control rule does not read indicator items: a large change in an indicator column of a control orbit stops nothing
    bump(orbits, info) = (o = only(o for o in orbits if o.id == 14); ind(o, "t4").series[end, 1] += 1.0)
    fxc = d86_fx("r-ctrl", d86_block([d86_arm("a", "S_D", ["t4"], ["Float64"], [100]; h1_check = true)], [d86_spec("tangent"; id = "t4")]))
    cc = first(c for c in ex_plan(fxc).chunks if c.cell == "controls")
    resc = ex_campaign(fxc; store = ex_store("d86-r-ctrl"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), mutate = ex_fault(bump; ord = cc.ordinal)))
    @test resc.status == "complete"
end
