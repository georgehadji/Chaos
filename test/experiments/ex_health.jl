# D53j health checks, D53k failure and rerun, D53l caps. Spec: Amendment D53, R5. Faults are injected through the mutate hook of
# RunHooks (called with the chunk's recorded orbits after computing, before the checks); time through the clock hook.

# stop after the first `k` chunks are written (the hook throws; ex_killed catches it): single-chunk checks need no full campaign
first_k(mutate, k) = EX.RunHooks(mutate = mutate, kill = ex_kill_after(k))
chunk_file(store, ord) = joinpath(store, "chunk-$(lpad(ord, 5, '0')).h5")
flag_ids(ch, flag) = sort([o.id for o in ch.orbits if flag in o.attrs["flags"]])

@testset "D53j non-finite values fail the chunk" begin
    fx = ex_fix(:small)
    # NaN in one series, first attempt only: the rerun passes and the file records both
    hit = Ref(0)
    nan!(orbits, info) = (o = orbits[2]; hit[] = o.id; first(values(o.lyap))[3, 1] = NaN)
    store = ex_store("j-nan")
    @test ex_killed(() -> ex_campaign(fx; store = store, hooks = first_k(ex_fault(nan!; ord = 1), 1)))
    a = EX.read_chunk(chunk_file(store, 1)).attrs
    @test a["status"] == "ok" && a["attempt"] == 2 && isempty(a["failure_reasons"])
    @test any(r -> occursin("non-finite", r) && occursin("orbit $(hit[])", r), a["first_attempt_reasons"])
    # Inf in the state, both attempts: the chunk fails, the reason names the orbit, the orbit is flagged excluded
    inf!(orbits, info) = (hit[] = orbits[1].id; orbits[1].state[2, 1] = Inf)
    store = ex_store("j-inf")
    @test ex_killed(() -> ex_campaign(fx; store = store, hooks = first_k(ex_fault(inf!; ord = 1, attempts = (1, 2)), 1)))
    ch = EX.read_chunk(chunk_file(store, 1))
    @test ch.attrs["status"] == "failed" && ch.attrs["attempt"] == 2 && ch.attrs["complete"] === true
    @test any(r -> occursin("non-finite", r) && occursin("orbit $(hit[])", r), ch.attrs["failure_reasons"])
    @test hit[] in flag_ids(ch, "excluded")
end

@testset "D53j diagnostics fraction: strictly more than fail_fraction fails; exceeding orbits flagged either way" begin
    base = ex_rb(:big)
    @test ex_plan(ex_fix(:big)).counts.items == 40
    for (ff, pass_n, fail_n) in ((0.05, 2, 3), (0.1, 4, 5), (0.25, 10, 11))
        fx = ex_fix_with(ex_rb_with(base, ("health", "fail_fraction") => ff), "j-frac-$ff"; panel = :big)
        for (n_exceed, want_status) in ((pass_n, "ok"), (fail_n, "failed"))
            exceed!(orbits, info) = for o in orbits[1:n_exceed]
                o.diagnostics["det"][end] = 2 * o.thresholds["det"]            # the running maximum passes the threshold
            end
            store = ex_store("j-frac-$ff-$n_exceed")
            ex_campaign(fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(exceed!; ord = 1, attempts = (1, 2))))
            ch = EX.read_chunk(chunk_file(store, 1))
            @test ch.attrs["status"] == want_status && length(ch.orbits) == 40
            @test ch.attrs["attempt"] == (want_status == "ok" ? 1 : 2)
            @test flag_ids(ch, "diagnostics exceeded") == collect(1:n_exceed)       # flagged in a passing chunk as well
            want_status == "failed" && @test any(r -> occursin("diagnostics exceeded", r), ch.attrs["failure_reasons"])
        end
    end
end

@testset "D53j exact-exponent controls" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    base = ex_small_run("ref")
    @test base.res.status == "complete"                                              # the unperturbed controls are inside control_bound
    c5 = ex_find_chunk(plan; method = "L5", cell = "controls", arm = "main")
    bump(id) = (orbits, info) -> (o = only(o for o in orbits if o.id == id); first(values(o.lyap))[end, 1] += 1.0)
    # an h1_check item over the bound stops the campaign at once with the reported item
    store = ex_store("j-ctrl")
    res = ex_campaign(fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(bump(14); ord = c5.ordinal)))
    @test res.status == "stopped: control bound" && maximum(res.ran) == c5.ordinal
    @test any(r -> occursin("control bound", r) && occursin("orbit 14", r) && occursin("L5", r), res.reasons)
    man = EX.read_json(res.manifest)
    @test man["status"] == "stopped: control bound" && length(man["chunks"]) == plan.counts.chunks
    @test count(c -> c["status"] == "pending", man["chunks"]) == plan.counts.chunks - length(res.ran)
    # the same control in another arm (not h1_check) does not stop the campaign
    c4 = ex_find_chunk(plan; method = "L4", cell = "controls", arm = "tau")
    res = ex_campaign(fx; store = ex_store("j-ctrl-other"), hooks = EX.RunHooks(mutate = ex_fault(bump(14); ord = c4.ordinal)))
    @test res.status == "complete" && isempty(res.failed)
    # a control with exponent zero (elliptic fixed point, orbit 13) is not checked
    res = ex_campaign(fx; store = ex_store("j-ctrl-zero"), hooks = EX.RunHooks(mutate = ex_fault(bump(13); ord = c5.ordinal)))
    @test res.status == "complete" && isempty(res.failed)
end

@testset "D53j setup check and drift" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    d64 = ex_find_chunk(plan; precision = "Double64")
    tol = 1e-10
    off(x) = (orbits, info) -> (orbits[1].state[1, 1] += x)
    # off by 2 setup_check_tol: first attempt fails with a setup-check reason, the rerun passes
    store = ex_store("j-setup-fail")
    ex_campaign(fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(off(2tol); ord = d64.ordinal)))
    a = EX.read_chunk(chunk_file(store, d64.ordinal)).attrs
    @test a["attempt"] == 2 && any(r -> occursin("setup check", r), a["first_attempt_reasons"])
    # at 0.5 setup_check_tol it passes at the first attempt, and the manifest records the check
    store = ex_store("j-setup-pass")
    res = ex_campaign(fx; store = store, hooks = EX.RunHooks(mutate = ex_fault(off(0.5tol); ord = d64.ordinal)))
    a = EX.read_chunk(chunk_file(store, d64.ordinal)).attrs
    @test a["attempt"] == 1 && a["status"] == "ok"
    chunks = EX.read_json(res.manifest)["chunks"]
    @test chunks[d64.ordinal]["setup_check"] == "ok" && chunks[1]["setup_check"] == "not applicable"
    # no Float64 item of the same orbit, method, tau, direction and mode: "no reference", no failure
    nr = ex_fix(:noref)
    res = ex_campaign(nr; store = ex_store("j-noref"))
    pl = ex_plan(nr)
    cd = ex_find_chunk(pl; precision = "Double64")
    @test res.status == "complete" && EX.read_json(res.manifest)["chunks"][cd.ordinal]["setup_check"] == "no reference"
    # drift: |p| above drift_bound at one grid step flags the orbit (strictly above); a flag, not a failure
    bound = 1e6
    for (p, flagged) in ((bound * (1 + 1e-9), true), (bound, false))
        hit = Ref(0)
        drift!(orbits, info) = (hit[] = orbits[1].id; orbits[1].state[5, 2] = p)
        store = ex_store("j-drift-$flagged")
        @test ex_killed(() -> ex_campaign(fx; store = store, hooks = first_k(ex_fault(drift!; ord = 1), 1)))
        ch = EX.read_chunk(chunk_file(store, 1))
        @test ch.attrs["status"] == "ok" && ch.attrs["attempt"] == 1
        @test flag_ids(ch, "drift") == (flagged ? [hit[]] : Int[])
    end
end

@testset "D53k failure and rerun" begin
    fx = ex_fix(:small)
    nan!(orbits, info) = (first(values(orbits[1].lyap))[2, 1] = NaN)
    # the fault clears: attempt 2 on the CPU with the same inputs, file ok, first attempt's reasons recorded
    log = NamedTuple[]
    store = ex_store("k-transient")
    @test ex_killed(() -> ex_campaign(fx; store = store, hooks = first_k(ex_fault(nan!; ord = 2, log = log), 2)))
    a = EX.read_chunk(chunk_file(store, 2)).attrs
    @test a["status"] == "ok" && a["attempt"] == 2 && !isempty(a["first_attempt_reasons"]) && isempty(a["failure_reasons"])
    @test [(l.attempt, l.backend) for l in log if l.ordinal == 2] == [(1, :cpu), (2, :cpu)]       # the rerun is on CPU
    # the fault stays: file failed, both reasons recorded, the failing orbit counts as excluded
    store = ex_store("k-persistent")
    @test ex_killed(() -> ex_campaign(fx; store = store, hooks = first_k(ex_fault(nan!; ord = 2, attempts = (1, 2)), 2)))
    ch = EX.read_chunk(chunk_file(store, 2))
    @test ch.attrs["status"] == "failed" && ch.attrs["attempt"] == 2
    @test !isempty(ch.attrs["failure_reasons"]) && !isempty(ch.attrs["first_attempt_reasons"])
    @test !isempty(flag_ids(ch, "excluded"))
    # stop rule: strictly more than stop_failed_fraction of the chunks of the plan failed twice
    nanat(ords) = (orbits, info) -> info.ordinal in ords && (first(values(orbits[1].lyap))[1, 1] = NaN)
    for (kind, M, ok_ords, stop_ords, last_ran) in ((:c40, 40, [5, 17], [5, 17, 30], 32), (:c20, 20, [5], [5, 12], 12))
        fx = ex_fix(kind)
        @test ex_plan(fx).counts.chunks == M
        # 2 of 40 and 1 of 20 are exactly 0.05: no stop
        res = ex_campaign(fx; store = ex_store("k-$kind-ok"), hooks = EX.RunHooks(mutate = (info, orbits) -> nanat(ok_ords)(orbits, info)))
        @test res.status == "complete with failed chunks" && sort(res.failed) == ok_ords && sort(res.ran) == collect(1:M)
        # 3 of 40 and 2 of 20 are above it: stop, recorded, resumable
        res = ex_campaign(fx; store = ex_store("k-$kind-stop"), hooks = EX.RunHooks(mutate = (info, orbits) -> nanat(stop_ords)(orbits, info)))
        @test res.status == "stopped: failed chunks" && sort(res.failed) == stop_ords && maximum(res.ran) <= last_ran
        man = EX.read_json(res.manifest)
        @test man["status"] == res.status && any(r -> occursin("failed", r), man["reasons"])
        sts = [c["status"] for c in man["chunks"]]
        @test length(sts) == M && all(s -> s in ("ok", "failed", "pending", "paused"), sts)
        @test count(==("failed"), sts) == length(stop_ords) && count(==("pending"), sts) == M - length(res.ran)    # I17: none silently missing
    end
end

@testset "D53l wall cap with an injected clock" begin
    # a large first-chunk factor so that only the wall cap acts (the cost cap has its own testset)
    fx = ex_fix_with(ex_rb_with(ex_rb(:small), ("health", "first_chunk_cost_factor") => 1e9), "l-wall")   # wall_cap_hours = 0.125 = 450 s; 6 launch groups
    store = ex_store("l-wall")
    res = ex_campaign(fx; store = store, threads = 1, hooks = EX.RunHooks(clock = ex_clock(info -> 150.0)))
    @test res.status == "stopped: wall cap"
    @test sort(res.ran) == collect(1:12)                            # three groups of four chunks: 450 s reach the cap before the fourth
    man = EX.read_json(res.manifest)
    @test man["status"] == "stopped: wall cap" && man["wall"]["total_seconds"] == 450.0 && man["wall"]["cap_hours"] == 0.125
    @test count(c -> c["status"] == "ok", man["chunks"]) == 12 && all(c -> c["status"] == "pending", man["chunks"][13:end])
    res2 = ex_campaign(fx; store = store, threads = 1, hooks = EX.RunHooks(clock = ex_clock(info -> 150.0)))
    man2 = EX.read_json(res2.manifest)
    @test res2.status == "stopped: wall cap" && isempty(res2.ran)    # the sessions' times add; the cap is reached at once
    @test length(man2["sessions"]) == 2 && man2["wall"]["total_seconds"] == 450.0
end

@testset "D53l first-chunk cost cap with an injected clock" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    line = 100.0 * 2 * 1e-9                                         # Float64 budget line: ns_cpu_thread x safety, seconds per orbit-step
    cost(factor) = info -> info.steps * line * (("main" in info.arms) ? factor : 1.0)
    main_ords = [c.ordinal for c in plan.chunks if "main" in c.arms]
    # 2.01 times the line pauses the arm; the other arms continue; the remaining chunks of the arm are paused
    store = ex_store("l-cost-hi")
    res = ex_campaign(fx; store = store, threads = 1, hooks = EX.RunHooks(clock = ex_clock(cost(2.01))))
    @test res.status == "paused arms" && res.paused == ["main"]
    @test sort(res.ran) == setdiff(1:plan.counts.chunks, main_ords[5:end])      # main's first group (4 chunks) ran, the rest did not
    man = EX.read_json(res.manifest)
    @test man["arms"]["main"]["status"] == "paused: cost cap" && man["arms"]["tau"]["status"] == "complete"
    @test all(o -> man["chunks"][o]["status"] == "paused", main_ords[5:end])
    # 1.99 times the line continues
    res = ex_campaign(fx; store = ex_store("l-cost-lo"), threads = 1, hooks = EX.RunHooks(clock = ex_clock(cost(1.99))))
    @test res.status == "complete" && isempty(res.paused)
    # a recorded cap override resumes the paused arm
    res = ex_campaign(fx; store = store, threads = 1, cap_override = ["main"], hooks = EX.RunHooks(clock = ex_clock(cost(1.0))))
    @test res.status == "complete" && sort(res.ran) == main_ords[5:end]
    man = EX.read_json(res.manifest)
    @test man["overrides"]["cap"] == ["main"] && man["arms"]["main"]["status"] == "complete"
end
