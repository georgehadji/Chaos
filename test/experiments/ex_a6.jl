# Amendment D57 and preregistration A6: optional per-method control bound `[run.health.control_bound_method]`.
# Spec: experiments-spec.md Amendment D57.

@testset "D57 block: per-method control bound is validated and merged" begin
    run4 = ex_rb(:small)
    withcbm(t) = (r = deepcopy(ex_rb_with(run4, ("health", "control_scope") => "exact_start")); r["health"]["control_bound_method"] = t; r)
    rb = EX.run_block(withcbm(Dict{String, Any}("L5" => 2.0)))
    @test rb.raw["health"]["control_bound_method"] == Dict("L5" => 2.0) && rb.raw["health"]["control_bound"] == run4["health"]["control_bound"]
    # a block without the key reads as before
    @test !haskey(EX.run_block(run4).raw["health"], "control_bound_method")
    # a later run block adds the table to the earlier health table
    pr = a5_prereg(ex_path("a6-merge-prereg.md"), run4; run5 = Dict{String, Any}("health" => Dict{String, Any}("control_bound_method" => Dict{String, Any}("L5" => 3e-8))))
    @test EX.read_run_block(pr).raw["health"]["control_bound_method"] == Dict("L5" => 3e-8)
    # refusals name the key
    @test ex_rejects(() -> EX.run_block(withcbm(Dict{String, Any}("L99" => 1e-8))), ["control_bound_method", "L99"])
    @test ex_rejects(() -> EX.run_block(withcbm(Dict{String, Any}("L5" => -1.0))), ["control_bound_method.L5"])
    @test ex_rejects(() -> EX.run_block(withcbm(Dict{String, Any}("L5" => "big"))), ["control_bound_method.L5"])
    @test ex_rejects(() -> EX.run_block(withcbm(Dict{String, Any}("L2" => 1e-8))), ["control_bound_method.L2", "control_methods"])
    @test ex_rejects(() -> EX.run_block(withcbm(3)), ["control_bound_method", "table"])
end

@testset "D57 health: the bound of a method replaces control_bound for that method only, and enters the chunk hash" begin
    run = ex_rb_with(ex_rb(:small), ("health", "control_scope") => "exact_start")
    with(t) = (r = deepcopy(run); t === nothing || (r["health"]["control_bound_method"] = t); r)
    bump(id, d) = (orbits, info) -> (o = only(o for o in orbits if o.id == id); first(values(o.lyap))[end, 1] += d)
    plan0 = ex_plan(a5_fx(with(nothing), "cbm-plan0"; orbits = a5_orbits()))
    c5 = ex_find_chunk(plan0; method = "L5", cell = "controls", arm = "main")
    c1 = ex_find_chunk(plan0; method = "L1", cell = "controls", arm = "main")
    # a deviation of 0.01 passes the run bound (0.05); a method bound of 1e-3 for L5 stops it, naming the bound
    fxs = a5_fx(with(Dict{String, Any}("L5" => 1e-3)), "cbm-stop"; orbits = a5_orbits())
    res = ex_campaign(fxs; store = ex_store("a6-stop"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(bump(14, 0.01); ord = c5.ordinal)))
    @test res.status == "stopped: control bound" && any(r -> occursin("orbit 14", r) && occursin("method L5", r) && occursin("above 0.001", r), res.reasons)
    fx0 = a5_fx(with(nothing), "cbm-nokey"; orbits = a5_orbits())
    res0 = ex_campaign(fx0; store = ex_store("a6-nokey"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(bump(14, 0.01); ord = c5.ordinal)))
    @test res0.status == "complete"
    # the same deviation in L1, which has no method bound, still passes under the run bound 0.05
    res1 = ex_campaign(fxs; store = ex_store("a6-other"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(bump(14, 0.01); ord = c1.ordinal)))
    @test res1.status == "complete"
    # a method bound above the run bound lets a larger deviation of that method pass
    fxl = a5_fx(with(Dict{String, Any}("L5" => 2.0)), "cbm-loose"; orbits = a5_orbits())
    resl = ex_campaign(fxl; store = ex_store("a6-loose"), hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(bump(14, 1.0); ord = c5.ordinal)))
    @test resl.status == "complete"
    # the table changes the hash of the chunks whose status it decides, and no hash when it is absent
    planl = ex_plan(fxl)
    @test any(c -> c.hash != plan0.chunks[c.ordinal].hash, planl.chunks)
    @test [c.hash for c in ex_plan(a5_fx(with(nothing), "cbm-plan1"; orbits = a5_orbits())).chunks] == [c.hash for c in plan0.chunks]
end

@testset "D57 the real preregistration: A6 over A5" begin
    chk = EX.check_prereg(EX_REAL_PREREG)
    @test chk.amendments == 6 && chk.header_sha256 == chk.computed_sha256
    h = EX.read_run_block(EX_REAL_PREREG).raw["health"]
    @test h["control_bound"] == 1e-8 && h["control_bound_method"] == Dict("L9" => 5e-8) && h["control_scope"] == "exact_start"
    @test all(m -> m in h["control_methods"], keys(h["control_bound_method"]))
end
