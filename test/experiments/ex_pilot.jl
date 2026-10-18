# D53q pilot. Spec: Amendment D53, R9; invariant I15. The pilot uses starts that are not panel orbits, writes outside the study store
# and reports costs. The injected clock gives a deterministic cost, so the projection can be checked exactly (the factor-2 check
# against a measured run of the real arms is stage 04: timing noise).

@testset "D53q pilot" begin
    fx = ex_fix(:small)
    plan = ex_plan(fx)
    run_pilot(out; fx = fx, cost_s = 150e-9, kw...) = EX.run_pilot(fx.prereg, fx.panel; out = out, threads = 1, thread_override = true,
        allow_dirty = true, code_version = EX_CODE, capacity = 64, hooks = EX.RunHooks(clock = ex_clock(info -> info.steps * cost_s)), kw...)
    before = Set(readdir(EX_DIR))
    out = ex_path("q-pilot"); rm(out; recursive = true, force = true)
    res = run_pilot(out)
    rep = res.report
    # isolation (I15): everything is under `out`; the pilot files are marked pilot; no study store or manifest is created
    @test Set(readdir(EX_DIR)) == union(before, Set(["q-pilot"]))
    pfiles = ex_chunk_files(joinpath(out, "store"))
    @test !isempty(pfiles) && all(f -> EX.read_chunk(f).attrs["pilot"] === true, pfiles)
    @test isfile(joinpath(out, "pilot-report.json")) && EX.read_json(joinpath(out, "pilot-report.json"))["not_study_data"] == rep["not_study_data"]
    @test occursin("not study data", rep["not_study_data"])
    # starts: orbits_per_k per K, none within 1e-6 of a panel start (including the extra starts), deterministic, keyed by the seed label
    orbs = fx.orbits
    nK = length(unique(o["k_index"] for o in orbs if o["k_index"] > 0))
    starts = rep["starts"]
    @test length(starts) == nK * fx.run["pilot"]["orbits_per_k"]
    pan = vcat([(o["x0"], o["p0"]) for o in orbs],
               [(parse(Float64, h["x0"]), parse(Float64, h["p0"])) for h in fx.run["h4"]])
    @test minimum(EX.torus_distance((s["x0"], s["p0"]), q) for s in starts for q in pan) > 1e-6
    @test [(s["k_index"], s["x0"], s["p0"]) for s in run_pilot(ex_path("q-pilot-2")).report["starts"]] == [(s["k_index"], s["x0"], s["p0"]) for s in starts]
    other = ex_fix_with(ex_rb_with(fx.run, ("pilot", "seed_label") => "ex-pilot-other"), "q-label")
    srt = [(s["x0"], s["p0"]) for s in run_pilot(ex_path("q-pilot-3"); fx = other).report["starts"]]
    @test srt != [(s["x0"], s["p0"]) for s in starts]
    # the report holds every field of R9 item 3
    for k in ("not_study_data", "starts", "cost", "arms", "k19", "first_departure", "diagnostics_share", "parallel_efficiency", "dirty_tree")
        @test haskey(rep, k)
    end
    @test all(c -> all(k -> haskey(c, k), ("method", "precision", "backend", "thread_seconds_per_orbit_step", "wall_seconds_per_orbit_step")), rep["cost"])
    @test Set((c["method"], c["precision"]) for c in rep["cost"]) == Set((it.method, ex_label(it.precision)) for it in plan.items)
    @test all(c -> isapprox(c["thread_seconds_per_orbit_step"], 150e-9; rtol = 1e-9), rep["cost"])        # the injected cost, one worker
    @test haskey(rep["k19"], "device_record") && haskey(rep["k19"], "prior_record")
    @test sort(collect(keys(rep["first_departure"]))) == string.(1:nK)                                  # per K: Double64 and BigFloat departure steps
    @test all(d -> haskey(d, "Double64") && haskey(d, "BigFloat"), values(rep["first_departure"]))
    @test sort(collect(keys(rep["diagnostics_share"]))) == sort(unique(it.method for it in plan.items))
    # projection: planned unique steps times the measured cost, against the budget line, the first-chunk cap and the wall cap
    for (name, a) in rep["arms"]
        steps = plan.steps.arms[name].steps_after
        @test a["planned_unique_steps"] == steps
        @test isapprox(a["projected_thread_seconds"], steps * 150e-9; rtol = 1e-9)
        @test all(k -> haskey(a, k), ("projected_wall_seconds", "budget_line_ns", "first_chunk_cap_ns", "wall_cap_hours"))
        @test a["first_chunk_cap_ns"] == fx.run["health"]["first_chunk_cost_factor"] * a["budget_line_ns"]
    end
    # a pilot does not change the result of a later study run (bit for bit)
    ref = ex_small_run("ref")
    ex_campaign(fx; store = ex_store("q-after"))
    @test ex_store_canon(ex_path("q-after"); attrs = true) == ex_store_canon(ref.store; attrs = true)
    @test isempty(ex_chunk_files(ex_path("q-after"))) == false
    # the dirty-tree override is recorded in the report; a clean tree gives none
    repo = ex_scratch_repo()
    @test run_pilot(ex_path("q-clean"); repo = repo, allow_dirty = false).report["dirty_tree"]["override"] === false
    write(joinpath(repo, "new.txt"), "n\n")
    @test ex_rejects(() -> run_pilot(ex_path("q-dirty-refused"); repo = repo, allow_dirty = false), ["untracked"])
    d = run_pilot(ex_path("q-dirty"); repo = repo, allow_dirty = true).report["dirty_tree"]
    @test d["override"] === true && "new.txt" in d["paths"]
    # the pilot refuses horizons that are not run types of the plan's grid: pilot horizon above the shortest arm horizon is a block error
    @test ex_rejects(() -> EX.run_block(ex_rb_with(fx.run, ("pilot", "horizon") => 1000)), ["pilot", "horizon"])
end
