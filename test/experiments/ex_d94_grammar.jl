# Amendment D94, criteria D94a (panel input grammar), D94b (absence of change), D94k (run block), D94u (layering and entry point).
# Spec: experiments-spec.md, "Amendment D94": R22, R25, I28. API assumed: the test plan, "Amendment D94 (stage 02)".

@testset "D94 references: hand tables of the helpers (guards, pass before the code exists)" begin
    # D_B allocation of the hand panel (plan, "hand tables"): |D_R| = 9, 6, 3 and b = 5
    @test d94_alloc([1, 6, 2], 5) == [1, 3, 1]                 # floors 0, 3, 1, remainders 5/9, 3/9, 1/9, one extra to the largest
    @test d94_alloc([1, 4, 1], 5) == [1, 3, 1]                 # floors 0, 3, 0, remainders 5/6, 2/6, 5/6, two extras to K 1 and K 3 (tie: smaller K first)
    @test d94_alloc([1, 2, 0], 5) == [1, 2, 0]                 # |D_R| <= b: D_B = D_R
    @test d94_alloc([3, 3], 4) == [2, 2] && d94_alloc([5, 5, 5], 4) == [2, 1, 1]      # ties to the smaller K index
    # D sizes: floor(f n + 1/2)
    @test [d94_dsize(f, n) for (f, n) in ((1 // 2, 3), (1 // 4, 2), (1 // 5, 7), (0 // 1, 5), (1 // 1, 5), (2 // 5, 5), (3 // 10, 5))] == [2, 1, 1, 0, 5, 2, 2]
    @test d94_rat("0.35") == 7 // 20 && d94_rat("1/3") == 1 // 3 && d94_rat("0.5") == 1 // 2 && d94_rat("1") == 1
    @test d94_runs_ids([[1, 3], [7, 7]]) == [1, 2, 3, 7]
    # the hand panel and the block are valid for the code that exists: a plan without release keys on a panel without rounds
    orbs = [Dict{String, Any}(k => v for (k, v) in o if !(k in ("round", "in_D"))) for o in d94_panel_orbits()]
    fx = d94_fix("plain", d94_block([ex_arm("ref", "P", ["L1"], ["Float64"], [20])]); orbits = orbs, r_max = nothing)
    plan = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE)
    @test plan.counts.items == 16 && plan.counts.chunks == 4
    @test length(D94_CELLS) == 9 && sum(c[4] for c in D94_CELLS) == 16 && length(d94_hand_D(3)) == 9
    @test d94_hand_cell(16) == (3, 3, "S") && d94_hand_cell(1) == (1, 1, "T")
end

# ------------------------------------------------------------------------------------------------------------------ D94a
const D94_B0 = Dict{String, Any}("screen_source" => "random", "screen_points" => 64)
d94_inp(over::AbstractDict = Dict{String, Any}(); drop = String[]) = ex_block(over; drop = drop)
d94_b1() = merge(D94_B0, Dict{String, Any}("r_max" => 3))
d94_b2() = merge(D94_B0, Dict{String, Any}("r_max" => 3, "quota" => d94_qtab(D94_FA.quota), "d_fraction" => d94_ftab(D94_FA.dfrac)))
d94_inp2(kw...) = (d = d94_b2(); foreach(((k, v),) -> v === :delete ? delete!(d, k) : (d[k] = v), kw); ex_block(d; drop = ["n_c"]))

@testset "D94a valid inputs are accepted and read back" begin
    df = EX.panel_inputs(ex_block())
    @test df.screen_source == "lattice" && df.r_max == 1 && df.gate_rounds == [1] && df.controls == "standard"          # documented defaults
    @test df.newton_grid_side == df.grid_side == 8                                                                          # P72: the value of grid_side
    r = EX.panel_inputs(d94_inp(D94_B0))
    @test r.screen_source == "random" && r.screen_points == 64
    @test EX.panel_inputs(d94_inp(Dict("controls" => "none"))).controls == "none"
    @test EX.panel_inputs(d94_inp(Dict("controls" => "standard"))).controls == "standard"
    r3 = EX.panel_inputs(d94_inp(d94_b1()))
    @test r3.r_max == 3 && r3.gate_rounds == [1]
    @test EX.panel_inputs(d94_inp(merge(d94_b1(), Dict("gate_rounds" => [2, 3])))).gate_rounds == [2, 3]
    n4 = EX.panel_inputs(d94_inp(Dict("newton_grid_side" => 4)))
    @test n4.newton_grid_side == 4 && n4.grid_side == 8
    q = EX.panel_inputs(d94_inp(merge(D94_B0, Dict("quota" => d94_qtab(D94_FA.quota))); drop = ["n_c"]))
    @test q.quota[1]["T"] == 3 && q.quota[2]["C"] == 3 && q.quota[3]["X"] == 3 && q.quota[1]["R"] == 2
    d = EX.panel_inputs(d94_inp(merge(d94_b1(), Dict("d_fraction" => d94_ftab(D94_FA.dfrac), "db_size" => 5))))
    @test d.d_fraction[2]["C"] == 7 // 20 && d.d_fraction[3]["X"] == 1 // 3 && d.d_fraction[1]["T"] == 2 // 5 && d.db_size == 5
    @test EX.panel_inputs(d94_inp(merge(d94_b1(), Dict("db_size" => 0)))).db_size == 0
    # every key at once, undefined classes of K = 0 given as 0, a closed-form fraction, a TOML float
    all_ = d94_b2()
    all_["quota"]["1"] = Dict{String, Any}("T" => 3, "R" => 2, "C" => 0, "S" => 0, "X" => 0)
    all_["d_fraction"]["1"]["T"] = 0.25
    all_["db_size"] = 4
    all_["gate_rounds"] = [1, 2]
    all_["controls"] = "none"
    all_["newton_grid_side"] = 6
    a = EX.panel_inputs(ex_block(all_; drop = ["n_c"]))
    @test a.r_max == 3 && a.gate_rounds == [1, 2] && a.db_size == 4 && a.newton_grid_side == 6 && a.controls == "none" && a.screen_points == 64
    @test a.d_fraction[1]["T"] == 1 // 4 && a.quota[1]["C"] == 0
    # the raw values stay available for the manifest
    @test a.raw["screen_source"] == "random" && a.raw["quota"]["2"]["T"] == 3
end

@testset "D94a malformed inputs are refused with the key, before any computation" begin
    base(sym) = sym === :b0 ? d94_inp(D94_B0) : sym === :b1 ? d94_inp(d94_b1()) : sym === :lat ? ex_block() : sym === :b2 ? d94_inp2() : error("base")
    function with(sym, f; drop = String[])
        d = deepcopy(base(sym))
        f(d)
        foreach(k -> delete!(d, k), drop)
        return d
    end
    cases = Tuple{String, Symbol, Function, Vector{String}, Vector{String}}[
        ("unknown source", :b0, d -> d["screen_source"] = "grid", ["screen_source"], String[]),
        ("screen_points with the lattice", :lat, d -> d["screen_points"] = 64, ["screen_points"], String[]),
        ("random source without screen_points", :b0, d -> delete!(d, "screen_points"), ["screen_points"], String[]),
        ("screen_points zero", :b0, d -> d["screen_points"] = 0, ["screen_points"], String[]),
        ("screen_points not an integer", :b0, d -> d["screen_points"] = 1.5, ["screen_points"], String[]),
        ("screen_points text", :b0, d -> d["screen_points"] = "64", ["screen_points"], String[]),
        ("controls unknown", :b0, d -> d["controls"] = "some", ["controls"], String[]),
        ("newton_grid_side zero", :b0, d -> d["newton_grid_side"] = 0, ["newton_grid_side"], String[]),
        ("newton_grid_side not an integer", :b0, d -> d["newton_grid_side"] = 2.5, ["newton_grid_side"], String[]),
        ("r_max zero", :b1, d -> d["r_max"] = 0, ["r_max"], String[]),
        ("r_max text", :b1, d -> d["r_max"] = "3", ["r_max"], String[]),
        ("gate_rounds empty", :b1, d -> d["gate_rounds"] = Int[], ["gate_rounds"], String[]),
        ("gate_rounds above r_max", :b1, d -> d["gate_rounds"] = [4], ["gate_rounds"], String[]),
        ("gate_rounds zero entry", :b1, d -> d["gate_rounds"] = [0], ["gate_rounds"], String[]),
        ("gate_rounds repeated", :b1, d -> d["gate_rounds"] = [1, 1], ["gate_rounds"], String[]),
        ("gate_rounds not a list", :b1, d -> d["gate_rounds"] = 1, ["gate_rounds"], String[]),
        ("gate_rounds not integers", :b1, d -> d["gate_rounds"] = [1.5], ["gate_rounds"], String[]),
        ("db_size negative", :b1, d -> d["db_size"] = -1, ["db_size"], String[]),
        ("db_size not an integer", :b1, d -> d["db_size"] = 2.5, ["db_size"], String[]),
        ("d_fraction with r_max absent", :b0, d -> d["d_fraction"] = d94_ftab(D94_FA.dfrac), ["d_fraction", "r_max"], String[]),
        ("db_size with r_max absent", :b0, d -> d["db_size"] = 3, ["db_size", "r_max"], String[]),
        ("quota together with n_c", :b2, d -> d["n_c"] = 3, ["quota", "n_c"], String[]),
        ("quota misses a K index", :b2, d -> delete!(d["quota"], "3"), ["quota"], String[]),
        ("quota misses a class that exists", :b2, d -> delete!(d["quota"]["2"], "X"), ["quota"], String[]),
        ("quota names C at K = 0 with a size above 0", :b2, d -> d["quota"]["1"]["C"] = 1, ["quota"], String[]),
        ("quota size not an integer", :b2, d -> d["quota"]["2"]["T"] = 2.5, ["quota"], String[]),
        ("quota size negative", :b2, d -> d["quota"]["2"]["T"] = -1, ["quota"], String[]),
        ("quota unknown K index", :b2, d -> d["quota"]["7"] = Dict{String, Any}("T" => 1), ["quota"], String[]),
        ("quota K index zero", :b2, d -> d["quota"]["0"] = Dict{String, Any}("T" => 1), ["quota"], String[]),
        ("quota unknown class", :b2, d -> d["quota"]["2"]["Q"] = 1, ["quota"], String[]),
        ("d_fraction above 1", :b2, d -> d["d_fraction"]["2"]["C"] = 1.5, ["d_fraction"], String[]),
        ("d_fraction negative", :b2, d -> d["d_fraction"]["2"]["C"] = -0.1, ["d_fraction"], String[]),
        ("d_fraction not a rational", :b2, d -> d["d_fraction"]["2"]["C"] = "abc", ["d_fraction"], String[]),
        ("d_fraction a bool", :b2, d -> d["d_fraction"]["2"]["C"] = true, ["d_fraction"], String[]),
        ("d_fraction unknown class", :b2, d -> d["d_fraction"]["2"]["Z"] = "0.5", ["d_fraction"], String[]),
        ("seed above the generator range", :b0, d -> d["seed"] = big(2)^64, ["seed"], String[]),
        ("a target is not a key", :b1, d -> d["target_T_s"] = 100, ["target_T_s"], String[]),
        ("a threshold is not a key", :b1, d -> d["threshold_W2"] = 10, ["threshold_W2"], String[]),
        ("a stop value is not a key", :b1, d -> d["stop_targets"] = [1, 2], ["stop_targets"], String[]),
    ]
    @test length(cases) >= 24 && count(c -> startswith(c[1], "quota"), cases) >= 2 && count(c -> startswith(c[1], "d_fraction"), cases) >= 2
    for (label, sym, f, keys, drop) in cases
        d = with(sym, f; drop = sym === :b2 ? String[] : drop)
        @test (label, ex_rejects(() -> EX.panel_inputs(d), keys)) == (label, true)
    end
    # through the builder: the same refusals write no file (panel, provenance sidecar, pool sidecar)
    for (n, (label, sym, f, keys, drop)) in enumerate(cases[[1, 3, 10, 12, 22, 24, 33, 39]])
        d = with(sym, f)
        pre = ex_path("d94a-bad$n.md")
        ex_write_prereg(pre, d)
        out = ex_path("d94a-bad$n-panel.toml")
        @test (label, ex_rejects(() -> EX.build_panel(pre; out = out), keys)) == (label, true)
        @test !any(isfile, (out, out * ".provenance.toml", out * ".pools.toml"))
    end
end

# ------------------------------------------------------------------------------------------------------------------ D94b
@testset "D94b absence of change: lattice panels equal the golden values of the code before the amendment" begin
    G = include(joinpath(@__DIR__, "ex_d94_golden.jl"))                  # commit named in the file; computed by the script of the plan
    bs = d94_lattice_builds()
    for (tag, g) in G.panel_sha
        b = bs[tag]
        @test d94_canon_sha(b.path) == g.sha
        @test length(b.panel["orbit"]) == g.orbits
        @test sort(collect(keys(b.panel["manifest"]))) == g.keys                  # no new manifest key for a panel without keys of the amendment
        @test sort(collect(keys(TOML.parsefile(b.path * ".provenance.toml")))) == g.side
        @test !isfile(b.path * ".pools.toml")                                      # the pool sidecar exists for the random source only
        @test !any(o -> haskey(o, "round") || haskey(o, "in_D") || haskey(o, "point_id"), b.panel["orbit"])
    end
    # the run blocks without a key of the amendment: block text and hash, plan hash and every chunk hash equal those of D86b (code before D86)
    g6 = include(joinpath(@__DIR__, "ex_d86_golden.jl"))
    fx = ex_fix(:small)
    @test EX.block_hash(EX.read_run_block(fx.prereg)) == g6.block_hash
    plan = ex_plan(fx)
    @test plan.hash == g6.plan_hash
    @test [c.hash for c in plan.chunks] == g6.chunk_hashes && [c.key for c in plan.chunks] == g6.chunk_keys
    @test all(c -> !occursin(":r", split(c.key, "cell=")[2]), plan.chunks)       # no round in a cell name of a panel without rounds
    @test !occursin("R_plan", plan.text)
end

# ------------------------------------------------------------------------------------------------------------------ D94k
@testset "D94k valid blocks read, hash and give the canonical text" begin
    run = d94_block()
    rb = EX.run_block(run)
    txt = EX.block_text(rb)
    @test EX.block_hash(rb) == ex_sha256(txt)
    @test occursin("run.arm[1].release = \"each_round\"\n", txt) && occursin("run.arm[2].release = \"closure\"\n", txt)
    @test occursin("run.arm[2].orbits = \"D\"\n", txt) && occursin("run.arm[3].orbits = \"D_B\"\n", txt)
    @test rb.raw["arm"][1]["release"] == "each_round" && rb.raw["arm"][3]["orbits"] == "D_B"
    # an arm without release is not subject to the log; a closure arm may use P, D, D_B; an each_round arm may use P or D
    free = EX.run_block(d94_block(d94_arms(rel = false)))
    @test all(a -> !haskey(a, "release"), free.raw["arm"]) && !occursin(".release", EX.block_text(free))
    for (o, rel) in (("P", "closure"), ("D", "closure"), ("D_B", "closure"), ("P", "each_round"), ("D", "each_round"))
        @test EX.run_block(d94_block([ex_arm("a", o, ["L1"], ["Float64"], [20]; release = rel)])).raw["arm"][1]["orbits"] == o
    end
    # the order of the keys in the file does not matter, and the text of a block with release keys is stable
    rev = d94_fix("rev", run; rev = true)
    @test EX.block_text(EX.read_run_block(rev.prereg)) == txt
    # the old blocks keep their hash (D86b golden): checked in D94b
    @test EX.block_hash(EX.run_block(ex_rb(:small))) == (include(joinpath(@__DIR__, "ex_d86_golden.jl"))).block_hash
end

@testset "D94k malformed blocks are refused with the key path, before any computation" begin
    run = d94_block()
    ed(r, kvs...) = ex_rb_with(r, kvs...)
    blk = [
        (ed(run, ("arm", 1, "release") => "weekly"), ["arm[1].release"]),
        (ed(run, ("arm", 1, "release") => 3), ["arm[1].release"]),
        (ed(run, ("arm", 1, "release") => ""), ["arm[1].release"]),
        (ed(run, ("arm", 1, "release") => ["closure"]), ["arm[1].release"]),
        (ed(run, ("arm", 1, "release") => true), ["arm[1].release"]),
        (ed(run, ("arm", 1, "release") => "Closure"), ["arm[1].release"]),
        (ed(run, ("arm", 3, "release") => "each_round"), ["arm[3]", "D_B"]),                  # D_B does not exist before closure
        (ed(run, ("arm", 1, "orbits") => "D_C"), ["arm[1].orbits"]),
        (ed(run, ("primary", "release") => "closure"), ["primary.release"]),
        (ed(run, ("health", "release") => "closure"), ["health.release"]),
        (ed(run, ("round_log",) => "x.log"), ["round_log"]),
        (ed(run, ("r_max",) => 3), ["r_max"]),
        (ed(d86_one(d86_spec("tangent")), ("indicator", 1, "release") => "closure"), ["indicator[1].release"]),
    ]
    @test length(blk) >= 12
    for (b, texts) in blk
        @test d94_rej(() -> EX.run_block(b), texts...)
    end
    # file-level refusals write nothing: the campaign stops at the plan, before the store is made
    store = ex_store("d94k-refused")
    manifest = store * "-manifest.json"
    # (a) release with a panel of one round
    one = d94_fix("one", run; orbits = d94_panel_orbits()[1:6], r_max = 1)
    emptylog = ex_path("d94k-empty.log")
    write(emptylog, "")
    @test d94_rej(() -> EX.plan_run(one.prereg, one.panel; code_version = EX_CODE, round_log = emptylog), "release")
    @test ex_rejects(() -> ex_campaign(one; store = store, round_log = emptylog), ["release"]) && !ispath(store) && !isfile(manifest)
    # (b) a block with release and no log at run time; the dry plan of --check is allowed and says that no round is released
    fx = d94_std()
    @test d94_rej(() -> EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE), "log")
    @test ex_rejects(() -> ex_campaign(fx; store = store), ["log"]) && !ispath(store) && !isfile(manifest)
    dry = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE, check = true)
    @test dry.rounds.R_plan == 0 && !dry.rounds.closed && occursin("no round released", dry.text)
    # (c) D with a panel that has no in_D
    noD = d94_fix("noD", run; orbits = d94_panel_orbits(with_D = false))
    @test d94_rej(() -> EX.plan_run(noD.prereg, noD.panel; code_version = EX_CODE, round_log = emptylog), "in_D")
    @test ex_rejects(() -> ex_campaign(noD; store = store, round_log = emptylog), ["in_D"]) && !ispath(store)
end

# ------------------------------------------------------------------------------------------------------------------ D94u
@testset "D94u layering and entry point" begin
    src = joinpath(dirname(dirname(@__DIR__)), "src")
    ex_files = [joinpath(src, "Experiments", f) for f in readdir(joinpath(src, "Experiments")) if endswith(f, ".jl")]
    imports(text, name) = occursin(Regex("(?m)^\\s*(?:using|import)\\s[^\\n]*\\b" * name * "\\b"), text)       # code lines only, not the prose of docstrings
    @test all(f -> !imports(read(f, String), "Analysis"), ex_files)                                       # Experiments imports no Analysis code
    for (dir, _, files) in walkdir(src), f in files
        below = first(splitpath(relpath(dir, src))) in ("Foundation", "Systems", "Integrators", "Tangent", "Lyapunov", "Indicators", "Diagnostics", "Shadow", "TimeSeries")
        endswith(f, ".jl") && below && @test !imports(read(joinpath(dir, f), String), "Experiments")
    end
    # the usage text gains exactly one option, the path of the round log, and the parser takes it
    @test count("--round-log", EX.CLI_USAGE) == 1
    @test EX.cli(["run"]; io = IOBuffer()) == 1
    buf = IOBuffer()
    @test EX.cli(["run", "p.md", "q.toml", "--round-log"]; io = buf) == 1 && !occursin("unknown option", String(take!(buf)))     # needs a value, is known
    buf = IOBuffer()
    @test EX.cli(["run", "p.md", "q.toml", "--round-limit", "3"]; io = buf) == 1 && occursin("unknown option", String(take!(buf)))
    # the entry point holds no numerics of the amendment
    cli_src = read(joinpath(src, "Experiments", "run_cli.jl"), String)
    @test !occursin(r"r_max|quota|d_fraction|db_size|screen_points|gate_rounds|newton_grid_side|R_plan", cli_src)
    # a refused log exits with status 1 and names the entry; the statuses of the guard are nonzero
    fx = d94_std()
    bad = ex_path("d94u-bad.log")
    write(bad, d94_log_text(fx, [(round = 1, decision = "continue")], r -> "0"^64; panel = "f"^64)[1])
    store = ex_store("d94u-cli")
    buf = IOBuffer()
    st = EX.cli(["run", fx.prereg, fx.panel, "--store", store, "--manifest", store * "-manifest.json", "--round-log", bad, "--thread-override",
                 "--allow-dirty", "--code-version", EX_CODE, "--threads", "4", "--capacity", "64"]; io = buf)
    msg = String(take!(buf))
    @test st == 1 && occursin("entry 1", msg) && !ispath(store)
    @test all(s -> EX.exit_status(s) != 0, ("complete with failed chunks", "stopped: preregistration changed", "paused: wall cap"))
end
