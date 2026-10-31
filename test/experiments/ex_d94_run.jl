# Amendment D94, runner: D94n (no peeking before closure) and D94o (resume across rounds). Spec: experiments-spec.md, "Amendment D94":
# E41 to E43, R26 item 4, I32, I33. API assumed: the test plan (RunHooks field `read`, keywords round_log, arms, up_to).
# The campaign is the hand-panel campaign of d94_campaign (4 runs: rounds 1, 2, 3, closure; the log is written by the test between the runs).

# canonical read-back of a store for I14 across plans: every dataset and every attribute, except times, host fields, the session lists, and the plan hash
# (a chunk file carries the hash of the plan that wrote it, and the plan of round r differs from the closed plan by R_plan and the log hash)
const D94_SKIP_ATTRS = (EX_ROOT_SKIP..., "plan_sha256")
function d94_canon_chunk(ch)
    io = IOBuffer()
    for k in sort(collect(keys(ch.attrs)))
        k in D94_SKIP_ATTRS || println(io, k, "=", repr(ch.attrs[k]))
    end
    for o in sort(ch.orbits; by = o -> o.id)
        for k in sort(collect(keys(o.attrs)))
            println(io, "orbit ", o.id, " ", k, "=", repr(o.attrs[k]))
        end
    end
    write(io, ex_canon(ch; attrs = false))
    return take!(io)
end
function d94_canon_store(store)
    io = IOBuffer()
    for f in ex_chunk_files(store)
        println(io, basename(f))
        write(io, d94_canon_chunk(EX.read_chunk(f)))
    end
    return take!(io)
end

@testset "D94n no peeking: no closure arm file before closure, only ref arm datasets and root attributes read" begin
    c = d94_campaign()
    plan = d94_plan(c.fx; log = c.logp)                                  # the closed plan holds every ordinal
    refs = Set(ch.ordinal for ch in d94_ref_chunks(plan))
    clos = Set(ch.ordinal for ch in d94_clo_chunks(plan))
    @test length(refs) == 9 && !isempty(clos) && isdisjoint(refs, clos)
    for (st, nfiles, nran) in zip(c.stages[1:3], (3, 6, 9), (3, 3, 3))
        ords = Set(d94_ordinal(f) for f in st.files)
        @test st.res.status == "complete"
        @test length(ords) == nfiles && isdisjoint(ords, clos)                                  # no closure arm chunk file exists before closure (0 counts)
        @test length(st.res.ran) == 3 && length(st.res.skipped) == nfiles - 3                   # the rounds before are skipped, not recomputed
        @test isempty(st.res.stale) && isempty(st.res.failed)
        ds = [r for r in st.reads if r.kind == :datasets]
        @test !isempty(st.reads) && all(r -> r.kind in (:attrs, :datasets), st.reads)
        @test all(r -> d94_ordinal(r.path) in refs, ds)                                           # the only datasets read are those of ref arm chunks
        @test all(r -> !(d94_ordinal(r.path) in clos), st.reads)                                   # not even the root attributes of a closure chunk file
    end
    @test sort(c.stages[2].res.skipped) == [1, 2, 3] && sort(c.stages[3].res.skipped) == collect(1:6)       # round 1 files are not stale in the later plans
    last = c.stages[4]
    @test last.res.status == "complete" && Set(d94_ordinal(f) for f in last.files) == Set(1:plan.counts.chunks)
    @test sort(last.res.ran) == sort(collect(clos)) && sort(last.res.skipped) == sort(collect(refs))      # at closure the closure arms run, the store holds them only then
    @test read(c.logp, String) == d94_log_text(c.fx, c.specs, c.digest)[1]                         # the module never writes the log (I37)
    # the manifest of a plan with rounds (R27)
    m = EX.read_json(c.store * "-manifest.json")
    r = m["rounds"]
    @test r["closed"] === true && r["R_plan"] == 3 && r["log_entries"] == 3 && r["log_sha256"] == ex_sha256(read(c.logp))
    @test r["D_size"] == 9 && r["DB_size"] == 5 && r["DB_allocation"] == [1, 3, 1]
    @test haskey(r, "release") && r["release"]["ref"] == "each_round" && r["release"]["dd"] == "closure" && r["release"]["db"] == "closure"
    # the first runs of the campaign wrote no manifest chunk entry of a closure arm: the entries of the stage 1 manifest are those of its plan
    @test m["plan_counts"]["chunks"] == plan.counts.chunks
end

@testset "D94o resume across rounds: round by round, killed and resumed, equals one run of the closed plan (I14, I33)" begin
    c = d94_campaign()
    fx = c.fx
    plan = d94_plan(fx; log = c.logp)
    nch = (3, 3, 3, plan.counts.chunks - 9)
    kills = NamedTuple[]
    for (n, seed) in enumerate((d94_seed('o', 1),))                                   # one killed campaign: 8 runs, the seed is recorded in the plan
        rng = StableRNG(seed)
        store = ex_store("d94o-kill-$n")
        logp = ex_path("d94o-kill-$n.log")
        write(logp, "")
        specs = NamedTuple[]
        digest = r -> d94_digest(d94_ref_hashes(store, r))
        for (stage, dec) in ((1, "continue"), (2, "continue"), (3, "close"), (4, nothing))
            k = rand(rng, 1:(nch[stage] - 1))                                                       # killed after k chunks of the stage, at a seeded point
            hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), kill = ex_kill_after(k))
            @test ex_killed(() -> ex_campaign(fx; store = store, round_log = logp, hooks = hooks, threads = 1 + (stage + n) % 4))
            res = ex_campaign(fx; store = store, round_log = logp, threads = 1 + (stage + n + 1) % 4)     # resumed with another thread count
            @test res.status == "complete" && isempty(res.stale)
            @test length(res.ran) == nch[stage] - k && length(res.skipped) >= k
            if dec !== nothing
                push!(specs, d94_spec(stage, dec))
                write(logp, d94_log_text(fx, specs, digest)[1])
            end
        end
        push!(kills, (; store, logp, n))
    end
    # one run of the closed plan with the log already closed (the digests are those of the round by round campaign)
    once = ex_store("d94o-once"); mkpath(once); foreach(f -> cp(f, joinpath(once, basename(f))), filter(f -> "ref" in d94_attr(f, "arms"), ex_chunk_files(kills[1].store)))   # H2: the log needs the ref chunks of its rounds in the store
    ologp = ex_path("d94o-once.log")
    kd = kills[1]
    write(ologp, d94_log_text(fx, [d94_spec(1), d94_spec(2), d94_spec(3, "close")], r -> d94_digest(d94_ref_hashes(kd.store, r)))[1])
    res = ex_campaign(fx; store = once, round_log = ologp, threads = 4)
    @test res.status == "complete" && length(res.ran) == length(d94_clo_chunks(plan)) && length(res.ran) + length(res.skipped) == plan.counts.chunks
    canon_once = d94_canon_store(once)
    @test !isempty(canon_once)
    for kd in kills
        @test [basename(f) for f in ex_chunk_files(kd.store)] == [basename(f) for f in ex_chunk_files(once)]
        canon = d94_canon_store(kd.store)
        @test canon == canon_once                                                                   # 0 differences, closure arms included
        @test ex_sha256(canon) == ex_sha256(canon_once)
        @test isempty(ex_tmp_files(kd.store))
    end
    # the campaign of the shared fixture (round by round without kills) is the same again
    @test d94_canon_store(c.store) == canon_once
    # a changed panel file makes the files of round 1 stale and the run stops (D53o(b))
    orbs = d94_panel_orbits()
    orbs[end]["x0"] += 1e-3
    orbs[end]["x0_dec"] = string(orbs[end]["x0"])
    fx2 = d94_fix("std2", d94_block(); orbits = orbs)
    cp(kills[1].store, ex_store("d94o-stale"))
    e0 = ex_path("d94o-empty.log")
    write(e0, "")
    err = try
        ex_campaign(fx2; store = ex_path("d94o-stale"), round_log = e0)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("stale", err.msg)
end

@testset "D94 references: the readers of the helpers on a campaign of the code that exists (guards)" begin
    fx = D94_FIXES["plain"]                                              # the hand panel without rounds and a block without release keys (grammar file)
    store = ex_store("d94-guard")
    res = ex_campaign(fx; store = store)
    files = ex_chunk_files(store)
    @test res.status == "complete" && length(files) == 4
    @test all(f -> occursin(r"^[0-9a-f]{64}$", d94_attr(f, "chunk_hash")) && "ref" in d94_attr(f, "arms") && occursin("cell=", d94_attr(f, "chunk_key")), files)
    @test all(f -> d94_ordinal(f) in 1:4, files) && isempty(d94_ref_hashes(store, 3))          # no round in these keys: nothing is a ref chunk of a round
    c1 = d94_canon_store(store)
    @test !isempty(c1) && c1 == d94_canon_store(store)
    store2 = ex_store("d94-guard2")
    ex_campaign(fx; store = store2, threads = 1)
    @test d94_canon_store(store2) == c1                                  # thread count does not change the canonical read-back
end
