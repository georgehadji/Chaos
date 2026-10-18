# D53o resume. Spec: Amendment D53, R8; invariants I13, I14. O10: an uninterrupted fixture campaign against campaigns killed at
# seeded points and resumed (with another thread count). Small fixture: 18 chunks, 6 launch groups.

@testset "D53o (a) a second run on a finished campaign does nothing" begin
    ref = ex_small_run("ref")
    store = ex_store("o-a")
    cp(ref.store, store)
    files = ex_chunk_files(store)
    before = Dict(f => (mtime(f), filesize(f)) for f in files)
    names0 = sort(readdir(store))
    res = ex_campaign(ref.fx; store = store)
    @test res.status == "complete" && isempty(res.ran) && length(res.skipped) == ref.plan.counts.chunks
    @test all(f -> (mtime(f), filesize(f)) == before[f], files) && sort(readdir(store)) == names0
    @test isempty(res.failed) && isempty(res.stale)
end

@testset "D53o (b) a stale file stops the run; with the option it is moved and recomputed" begin
    ref = ex_small_run("ref")
    store = ex_store("o-b")
    cp(ref.store, store)
    files = ex_chunk_files(store)
    before = Dict(f => mtime(f) for f in files)
    # a changed code version changes every item hash: every file is stale
    err = try
        ex_campaign(ref.fx; store = store, code_version = "ex-code-2"); nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("stale", err.msg) && occursin("chunk-00001.h5", err.msg)
    @test all(f -> mtime(f) == before[f], files) && !isdir(joinpath(store, "stale"))           # nothing touched
    res = ex_campaign(ref.fx; store = store, code_version = "ex-code-2", stale = :move)
    @test res.status == "complete" && sort(res.ran) == collect(1:ref.plan.counts.chunks)
    @test length(readdir(joinpath(store, "stale"))) == ref.plan.counts.chunks                  # moved, never deleted
    @test all(f -> EX.read_chunk(f).attrs["code_version"] == "ex-code-2", ex_chunk_files(store))
    # a file whose chunk hash differs by one value that decides the status is stale as well
    other = ex_fix_with(ex_rb_with(ref.fx.run, ("health", "fail_fraction") => 0.2), "o-b-ff")
    cp(ref.store, ex_store("o-b2"))
    @test ex_rejects(() -> ex_campaign(other; store = ex_path("o-b2")), ["stale"])
end

@testset "D53o (c) killed after k chunks for 10 seeded k, resumed with another thread count (O10)" begin
    ref = ex_small_run("ref")
    M = ref.plan.counts.chunks
    refcanon = ex_store_canon(ref.store; attrs = true)
    href = ex_sha256(refcanon)
    rng = StableRNG(ex_s53('o', 1))
    ks = sortperm(rand(rng, M - 1))[1:10]
    @test length(unique(ks)) == 10
    for (j, k) in enumerate(ks)
        store = ex_store("o-c-$j")
        @test ex_killed(() -> ex_campaign(ref.fx; store = store, threads = 4, hooks = EX.RunHooks(kill = ex_kill_after(k))))
        @test length(ex_chunk_files(store)) == k
        res = ex_campaign(ref.fx; store = store, threads = 1 + j % 2)                           # 2 or 1 threads instead of 4
        @test res.status == "complete" && sort(res.ran) == collect((k + 1):M)
        canon = ex_store_canon(store; attrs = true)
        @test canon == refcanon                                                                 # 0 differences
        @test ex_sha256(canon) == href                                                          # SHA-256 of the canonical read-back
        man = EX.read_json(res.manifest)
        @test length(man["sessions"]) == 2                                                      # the session list of the killed run is kept
    end
end

@testset "D53o (d) a final failed chunk is not retried on resume" begin
    fx = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "stop_failed_fraction") => 0.5), "o-d")
    ref = ex_small_run("ref")
    plan = ex_plan(fx)
    store = ex_store("o-d")
    bad(orbits, info) = (first(values(orbits[1].lyap))[1, 1] = NaN)
    res = ex_campaign(fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(bad; ord = 3, attempts = (1, 2))))
    @test res.failed == [3] && res.status == "complete with failed chunks"
    f3 = joinpath(store, "chunk-00003.h5")
    a3 = EX.read_chunk(f3).attrs
    @test a3["status"] == "failed" && a3["attempt"] == 2 && a3["complete"] === true && !isempty(a3["failure_reasons"])
    t3 = mtime(f3)
    res2 = ex_campaign(fx; store = store)                                                       # no fault now: the failure is final
    @test isempty(res2.ran) && res2.status == "complete with failed chunks" && 3 in res2.skipped
    @test mtime(f3) == t3 && EX.read_chunk(f3).attrs["status"] == "failed"
    res3 = ex_campaign(fx; store = store, retry_failed = true)                                  # on request only
    @test res3.ran == [3] && res3.status == "complete"
    @test EX.read_chunk(f3).attrs["status"] == "ok"
    @test any(n -> startswith(n, "chunk-00003"), readdir(joinpath(store, "aside")))              # the replaced file is moved aside (D53n)
    @test ex_store_canon(store) == ex_store_canon(ref.store)                                    # datasets equal the clean campaign
end
