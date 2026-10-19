# Amendment D56 and preregistration A5: run blocks of successive amendments merged, per-arm stop exemption with per-orbit
# orthonormalisation-breakdown flags, control scope, L2 arms in wrapped coordinates, budget. Spec: experiments-spec.md Amendment D56.

rev_zero() = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0))

# a preregistration with the run table `run4` (Amendment A4) and, if given, the partial run table `run5` (Amendment A5)
function a5_prereg(path, run4; run5 = nothing)
    text = ex_prereg_text_run(run4)
    run5 === nothing || (text = replace(text, "End of file." => "## Amendment A5\n\n```toml\n" * ex_emit_run(run5) * "```\n\nEnd of file."))
    write(path, replace(text, "0"^64 => ex_sha256(ex_prepare(text)); count = 1))
    return path
end

# fixture with a chosen panel (orbit dictionaries) and run table
function a5_fx(run, tag; orbits = ex_panel_orbits(:std), run5 = nothing)
    pp = ex_write_panel(ex_path("a5-$tag-panel.toml"); orbits = orbits)
    return (; run, panel = pp, prereg = a5_prereg(ex_path("a5-$tag-prereg.md"), run; run5 = run5), tag, orbits)
end

# orbit 13 of the small panel (elliptic fixed point (pi, 0)) becomes a hyperbolic control with the exponent 1: its binary64 start is not exact
function a5_orbits()
    os = ex_panel_orbits(:std)
    o = only(o for o in os if o["id"] == 13)
    o["stability"], o["exponent"] = "hyperbolic", 1.0
    return os
end

@testset "D56 merge of the run tables of successive amendments" begin
    run4 = ex_rb(:small)
    run5 = Dict{String, Any}("health" => Dict{String, Any}("control_scope" => "exact_start"),
                             "budget" => Dict{String, Any}("Float64" => Dict{String, Any}("ns_cpu_thread" => 123.5)),
                             "arm" => [ex_arm("tau", "S_D", ["L4"], ["Float64"], [100]; stop_exempt = true),
                                       ex_arm("extra", "S_D", ["L1"], ["Float64"], [100]; coordinate_mode = ["wrapped"])])
    pr = a5_prereg(ex_path("a5-merge-prereg.md"), run4; run5 = run5)
    rb = EX.read_run_block(pr)
    want = ex_rb_with(run4, ("health", "control_scope") => "exact_start", ("budget", "Float64", "ns_cpu_thread") => 123.5)
    want["arm"] = Any[want["arm"][1], run5["arm"][1], want["arm"][3], want["arm"][4], run5["arm"][2]]   # replaced in place, new arm last
    @test EX.block_hash(rb) == EX.block_hash(EX.run_block(want))
    @test [a["name"] for a in rb.raw["arm"]] == ["main", "tau", "h4f", "h4d", "extra"]
    tau = rb.raw["arm"][2]
    @test !haskey(tau, "tau") && tau["stop_exempt"] === true                  # an arm is replaced as a whole: the old `tau` list is gone
    @test rb.raw["budget"]["Float64"]["ns_gpu"] == 10.0 && rb.raw["budget"]["Float64"]["ns_cpu_thread"] == 123.5 && rb.raw["health"]["control_scope"] == "exact_start"
    @test length(rb.raw["h4"]) == 2 && rb.raw["panel_prereg_sha256_accepted"] == run4["panel_prereg_sha256_accepted"]
    # a file with one block reads as before; a block without the new keys has no `control_scope` or `stop_exempt`
    one = EX.read_run_block(a5_prereg(ex_path("a5-one-prereg.md"), run4))
    @test !haskey(one.raw["health"], "control_scope") && !any(a -> haskey(a, "stop_exempt"), one.raw["arm"])
    # validation after the merge names the key
    bad(r5) = () -> EX.read_run_block(a5_prereg(ex_path("a5-bad-prereg.md"), run4; run5 = r5))
    @test ex_rejects(bad(Dict{String, Any}("health" => Dict{String, Any}("control_scope" => "some"))), ["control_scope"])
    @test ex_rejects(bad(Dict{String, Any}("health" => Dict{String, Any}("bogus" => 1))), ["bogus"])
    @test ex_rejects(bad(Dict{String, Any}("arm" => [ex_arm("tau", "S_D", ["L4"], ["Float64"], [100]; stop_exempt = 1)])), ["stop_exempt"])
end

@testset "D56 the real preregistration: A5 over A4" begin
    chk = EX.check_prereg(EX_REAL_PREREG)
    @test chk.amendments == 5 && chk.header_sha256 == chk.computed_sha256
    rb = EX.read_run_block(EX_REAL_PREREG)
    arms = Dict(a["name"] => a for a in rb.raw["arm"])
    # L2 runs only in wrapped coordinates; the other methods of main and long are unchanged
    for a in values(arms)
        "L2" in a["methods"] && @test get(a, "coordinate_mode", [rb.raw["primary"]["coordinate_mode"]]) == ["wrapped"]
    end
    @test arms["main"]["methods"] == ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"] && arms["main"]["h1_check"] === true
    @test arms["main_l2"]["methods"] == ["L2"] && arms["main_l2"]["orbits"] == "P" && arms["long_l2"]["horizons"] == [100000000]
    @test arms["eps_sweep"]["eps"] == [1e-6, 1e-8, 1e-10, 1e-12]
    @test [n for (n, a) in arms if get(a, "stop_exempt", false)] == ["tau_sweep"] && arms["tau_sweep"]["tau"] == [1, 2, 4, 8, 16]
    @test rb.raw["health"]["control_scope"] == "exact_start" && rb.raw["health"]["control_bound"] == 1e-8
    b = rb.raw["budget"]
    @test (b["Float64"]["ns_cpu_thread"], b["Double64"]["ns_cpu_thread"], b["BigFloat"]["ns_cpu_thread"]) == (207.6, 6634.6, 296311.3)
    @test (b["Float64"]["safety"], b["Double64"]["safety"], b["BigFloat"]["safety"]) == (2.0, 3.0, 3.0) && b["Float64"]["ns_gpu"] == 2.468
    @test rb.raw["panel_prereg_sha256_accepted"] == ["650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46"]
    @test EX.block_hash(rb) == EX.block_hash(EX.run_block(ex_real_run()))          # the module's merge equals the independent one
end

@testset "D56 exactness of a control start is decided at plan time" begin
    fx = a5_fx(ex_rb(:small), "exact"; orbits = a5_orbits())
    plan = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE)
    exact = Dict(it.orbit => it.control_exact for it in plan.items)
    @test exact[14] === true && exact[15] === true                             # (0, 0) is a fixed point in binary64; Y5 is linear
    @test exact[13] === false                                                  # (pi, 0): binary64 pi leaves the fixed point
    # a periodic point with a winding is never exact in the lifted binary64 map (the lifted state moves by the binary64 2 pi)
    wound = a5_orbits()
    o = only(o for o in wound if o["id"] == 13)
    o["K_text"], o["K"], o["x0"], o["p0"], o["wind_m"], o["wind_j"] = "10", 10.0, 0.6793899267509509, 1.4708470335662112e-77, 1, 1
    o["x0_dec"], o["p0_dec"] = "0.6793899267509509", "1.4708470335662112e-77"
    fxw = a5_fx(ex_rb(:small), "exact-wound"; orbits = wound)
    @test !any(it.control_exact for it in EX.plan_run(fxw.prereg, fxw.panel; code_version = EX_CODE).items if it.orbit == 13)
    @test all(it -> !it.control_exact, [it for it in plan.items if it.design_class != "control"])
    # the chunk hash follows the exactness (it is derived from the panel)
    other = a5_orbits()
    only(o for o in other if o["id"] == 13)["x0"] = 0.0
    only(o for o in other if o["id"] == 13)["x0_dec"] = "0.0"
    fx2 = a5_fx(ex_rb(:small), "exact2"; orbits = other)
    p2 = EX.plan_run(fx2.prereg, fx2.panel; code_version = EX_CODE)
    @test [it.control_exact for it in p2.items if it.orbit == 13] == fill(true, count(it -> it.orbit == 13, p2.items))
    @test any(c -> 13 in c.orbits && c.hash != plan.chunks[c.ordinal].hash, p2.chunks)
end

@testset "D56 control scope: the bound for exact starts, INFO for the others" begin
    run = ex_rb(:small)
    scoped = ex_rb_with(run, ("health", "control_scope") => "exact_start")
    plan = ex_plan(a5_fx(scoped, "scope-plan"; orbits = a5_orbits()))
    c5 = ex_find_chunk(plan; method = "L5", cell = "controls", arm = "main")
    # without the key the off-orbit control (exponent 1, measured about 0) stops the campaign as before
    fx0 = a5_fx(run, "scope-all"; orbits = a5_orbits())
    res0 = ex_campaign(fx0; store = ex_store("a5-scope-all"), hooks = rev_zero())
    @test res0.status == "stopped: control bound" && any(r -> occursin("orbit 13", r), res0.reasons)
    # with exact_start the campaign completes and keeps the deviation
    fx = a5_fx(scoped, "scope"; orbits = a5_orbits())
    store = ex_store("a5-scope")
    res = ex_campaign(fx; store = store, hooks = rev_zero())
    @test res.status == "complete" && isempty(res.failed)
    ch = EX.read_chunk(joinpath(store, "chunk-$(lpad(c5.ordinal, 5, '0')).h5"))
    o13, o14 = (only(o for o in ch.orbits if o.id == id) for id in (13, 14))
    @test o13.attrs["control_scope"] == "info" && o13.attrs["control_deviation"] > 0.5
    @test o14.attrs["control_scope"] == "bound" && o14.attrs["control_deviation"] < 0.05
    man = EX.read_json(res.manifest)
    info = [(i["orbit"], i["method"]) for c in man["chunks"] for i in c["control_info"]]
    @test sort(info) == [(13, "L1"), (13, "L5")] && all(c -> haskey(c, "control_info"), man["chunks"])
    # an exact control over the bound still stops
    bump(id) = (orbits, info) -> (o = only(o for o in orbits if o.id == id); first(values(o.lyap))[end, 1] += 1.0)
    res2 = ex_campaign(fx; store = ex_store("a5-scope-stop"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(bump(14); ord = c5.ordinal)))
    @test res2.status == "stopped: control bound" && any(r -> occursin("orbit 14", r), res2.reasons)
    # a resumed store keeps the INFO records (rebuilt from the files)
    man2 = EX.read_json(ex_campaign(fx; store = store, hooks = rev_zero()).manifest)
    @test sort([(i["orbit"], i["method"]) for c in man2["chunks"] for i in c["control_info"]]) == [(13, "L1"), (13, "L5")]
end

@testset "D56 orthonormalisation breakdown in an exempt arm is a per-orbit result" begin
    orb = Dict{String, Any}("id" => 1, "class" => "T", "control" => false, "control_kind" => "none", "k_index" => 1, "K_text" => "10", "K" => 10.0,
        "system" => "standard", "x0" => 4.086524818927347979524711263366043567657470703125, "p0" => 4.9210103675371374976066363160498440265655517578125,
        "x0_dec" => "4.086524818927347979524711263366043567657470703125", "p0_dec" => "4.9210103675371374976066363160498440265655517578125",
        "draw_index" => 1, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, 1)), "in_SD" => true, "in_SL" => true)
    function tau_run(exempt)
        run = ex_rb_with(ex_rb(:small), ("pilot", "horizon") => 32)
        delete!(run, "h4")
        run["arm"] = [ex_arm("tau", "S_D", ["L5"], ["Float64"], [1600]; tau = [16], (exempt ? (; stop_exempt = true) : (;))...)]
        return a5_fx(run, "bd-$exempt"; orbits = [orb])
    end
    # tau = 16 at K = 10 is beyond tau_max(K) = 5: the Householder reorthonormalisation breaks down at a step of the run
    store = ex_store("a5-bd")
    res = ex_campaign(tau_run(true); store = store, hooks = rev_zero())
    ch = EX.read_chunk(joinpath(store, "chunk-00001.h5"))
    o = only(ch.orbits)
    @test res.status == "complete" && isempty(res.failed) && ch.attrs["status"] == "ok" && ch.attrs["attempt"] == 1
    @test o.attrs["flags"] == ["breakdown_tau"] && o.attrs["breakdown_kind"] == "orthonormalisation breakdown"
    @test 0 < o.attrs["breakdown_step"] < 1600 && all(isnan, o.lyap["L5"][end, :]) && !isnan(o.lyap["L5"][1, 1])   # the rows from the failure are blank
    @test EX.read_json(res.manifest)["chunks"][1]["flags"]["breakdown_tau"] == [1]
    # the same run in an arm that is not exempt fails the chunk
    res2 = ex_campaign(tau_run(false); store = ex_store("a5-bd-plain"), hooks = rev_zero())
    @test res2.failed == [1] && res2.status != "complete"
    # the failure text is read, not guessed
    msg = "NumericalHealthError: orthonormalisation breakdown = 3.5e15 (threshold 0.5) at step 48, orbit 1, backend cpu"
    @test EX.breakdown_of(msg) == (kind = "orthonormalisation breakdown", step = 48)
    @test EX.breakdown_of("NumericalHealthError: separation = 1.0 (threshold 0.5) at step 4, orbit 1, backend cpu") === nothing
    @test EX.breakdown_of("") === nothing
end

@testset "D56 an exempt arm is outside the failed-chunk stop rule; other failures still fail" begin
    base = ex_rb_with(ex_rb(:small), ("health", "stop_failed_fraction") => 0.01)
    plan = ex_plan(a5_fx(base, "stop-plan"))
    ctau = ex_find_chunk(plan; method = "L4", arm = "tau")
    @test ctau.arms == ["tau"]
    nan!(orbits, info) = (first(values(orbits[1].lyap))[1, 1] = NaN)
    hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(nan!; ord = ctau.ordinal, attempts = (1, 2)))
    # not exempt: one failed chunk of 18 is above 0.01 and stops the campaign
    res = ex_campaign(a5_fx(base, "stop-plain"); store = ex_store("a5-stop-plain"), hooks = hooks)
    @test res.failed == [ctau.ordinal] && res.status == "stopped: failed chunks"
    # exempt: a failure that is not a breakdown still fails the chunk but is not counted by the rule
    ex = ex_rb_with(base, ("arm", 2, "stop_exempt") => true)
    res = ex_campaign(a5_fx(ex, "stop-exempt"); store = ex_store("a5-stop-exempt"), hooks = hooks)
    @test res.failed == [ctau.ordinal] && res.status == "complete with failed chunks"
    # an exempt arm does not exempt a chunk that another arm also holds
    both = ex_rb_with(base, ("arm", 2, "stop_exempt") => true, ("arm", 2, "methods") => ["L1", "L4"])
    pboth = ex_plan(a5_fx(both, "stop-both"))
    shared = ex_find_chunk(pboth; method = "L1", arm = "tau")
    @test "main" in shared.arms && !EX.chunk_exempt(pboth, shared) && EX.chunk_exempt(pboth, ex_find_chunk(pboth; method = "L4", arm = "tau"))
end
