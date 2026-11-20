# Amendment D125, stage 02: D125b (panel input of the scheme "systems") and D125i (run block grammar of the new kinds).
# Refusal tables are checked as one assertion per table (the list of cases that were not refused as stated is empty), so a case that today's code refuses
# for another reason cannot make a RED run look green case by case. The acceptance testsets come first in each group and fail today.

@testset "D125b valid panel inputs are accepted and read back" begin
    for fx in (D125_FA, D125_FB)
        pin = EX.panel_inputs(d125_block(fx))                                              # RED: ArgumentError (unknown key scheme) today
        @test pin.scheme == "systems"
        @test [c.id for c in pin.cells] == [c["id"] for c in fx.cells]
        @test [c.kind for c in pin.cells] == [c["kind"] for c in fx.cells]
        @test pin.classes.regular == fx.classes[1] && pin.classes.weak == fx.classes[2] && pin.classes.strong == fx.classes[3]
        @test pin.rules.ratio_regular_max == fx.rules["ratio_regular_max"] && pin.rules.weak_fraction == fx.rules["weak_fraction"]
        @test pin.rules.weak_quantile == fx.rules["weak_quantile"]                          # exact text kept (E49: an exact rational)
        @test [c.n_start for c in pin.cells] == [c["n_start"] for c in fx.cells]
        @test [c.id for c in pin.controls] == [c["id"] for c in fx.controls]
        @test pin.seed == fx.seed                                                           # seed given in the block; absent = default of numerics-defaults.toml
    end
    # a panel input without `scheme` is read as before (I40)
    @test EX.panel_inputs(ex_block()).n_c == 3
    # the seed defaults from numerics-defaults.toml when the key is absent (E0)
    @test EX.panel_inputs(d125_edit(d125_block(D125_FA), "seed" => :delete)).seed == EX_SEED
end

# (name, keys the error names, edits); 55 cases, one per refusal class of R29 item 3 (at least three each for quota and fraction tables)
function d125_panel_cases()
    b = d125_block(D125_FA)
    c = (name, keys, edits...) -> (name, keys, d125_edit(b, edits...))
    cases = Any[
        c("no classes table", ["classes"], "classes" => :delete), c("no rules table", ["rules"], "rules" => :delete),
        c("no cell", ["cell"], "cell" => :delete), c("unknown top key", ["threshold_x"], "threshold_x" => 1.0),
        c("base key K", ["K"], "K" => ["0", "1"]), c("base key grid_side", ["grid_side"], "grid_side" => 8), c("base key n_c", ["n_c"], "n_c" => 3),
        c("base key q_max", ["q_max"], "q_max" => 2), c("base key rho", ["rho"], "rho" => 1e-10), c("base key r_max", ["r_max"], "r_max" => 2),
        c("base key gate_rounds", ["gate_rounds"], "gate_rounds" => [1]), c("base key screen_source", ["screen_source"], "screen_source" => "random"),
        c("base key controls", ["controls"], "controls" => "none"), c("base key N_scr", ["N_scr"], "N_scr" => 2000),
        c("class key missing", ["weak"], ("classes", "weak") => :delete), c("class key unknown", ["extra"], ("classes", "extra") => "x1"),
        c("class ids equal", ["classes"], ("classes", "weak") => "reg"), c("class id ill-formed", ["classes"], ("classes", "strong") => "1st"),
        c("rule key missing", ["floor_abs"], ("rules", "floor_abs") => :delete), c("rule key unknown", ["ratio_x"], ("rules", "ratio_x") => 1.0),
        c("regular_max >= chaotic_min", ["ratio_regular_max"], ("rules", "ratio_regular_max") => 6.0),
        c("chaotic_min > weak_min", ["ratio_chaotic_min"], ("rules", "ratio_chaotic_min") => 9.0),
        c("weak_fraction 0", ["weak_fraction"], ("rules", "weak_fraction") => 0.0), c("weak_fraction above 1", ["weak_fraction"], ("rules", "weak_fraction") => 1.5),
        c("weak_quantile 1", ["weak_quantile"], ("rules", "weak_quantile") => "1"), c("weak_quantile text", ["weak_quantile"], ("rules", "weak_quantile") => "most"),
        c("floor_abs 0", ["floor_abs"], ("rules", "floor_abs") => 0.0), c("floor_log negative", ["floor_log"], ("rules", "floor_log") => -1.0),
        c("cell kind missing", ["kind"], ("cell", 1, "kind") => :delete), c("cell kind unknown", ["kind"], ("cell", 1, "kind") => "lorenz"),
        c("chain key in froeschle cell", ["boundary"], ("cell", 1, "boundary") => "open"),
        c("duplicate cell id", ["id"], ("cell", 2, "id") => "y3a"), c("cell id ill-formed", ["id"], ("cell", 1, "id") => "Y 3"),
        c("n_start missing", ["n_start"], ("cell", 1, "n_start") => :delete), c("n_start 0", ["n_start"], ("cell", 1, "n_start") => 0),
        c("scr_back 1", ["scr_back"], ("cell", 1, "scr_back") => 1), c("scr_tau 0", ["scr_tau"], ("cell", 1, "scr_tau") => 0),
        c("scr_steps not a multiple", ["scr_steps"], ("cell", 1, "scr_steps") => 401), c("scr_steps 0", ["scr_steps"], ("cell", 1, "scr_steps") => 0),
        c("parameter not exact", ["nu"], ("cell", 1, "nu") => "abc"), c("parameter not finite", ["kappa"], ("cell", 1, "kappa") => "1/0"),
        c("N = 1", ["N"], ("cell", 2, "N") => 1), c("N above 25", ["N"], ("cell", 2, "N") => 26),
        c("ring with N = 2", ["boundary"], ("cell", 2, "boundary") => "ring"), c("boundary unknown", ["boundary"], ("cell", 2, "boundary") => "closed"),
        c("energy at escape", ["energy"], ("cell", 3, "energy") => "1/6"), c("energy 0", ["energy"], ("cell", 3, "energy") => "0"),
        c("q2 window empty", ["q2_lo"], ("cell", 3, "q2_lo") => "0.4"), c("p2 window reversed", ["p2_lo"], ("cell", 3, "p2_lo") => "0.3"),
        c("integrator on a map cell", ["integrator"], ("cell", 1, "integrator") => Dict("method" => "Leapfrog", "scheme" => "DKD", "step" => "0.1")),
        c("flow cell without integrator", ["integrator"], ("cell", 3, "integrator") => :delete),
        c("quota class missing", ["quota"], ("cell", 1, "quota", "wk") => :delete), c("quota class unknown", ["quota"], ("cell", 1, "quota", "zz") => 1),
        c("quota non-integer", ["quota"], ("cell", 1, "quota", "reg") => 1.5), c("quota negative", ["quota"], ("cell", 1, "quota", "reg") => -1),
        c("fraction class missing", ["d_fraction"], ("cell", 1, "d_fraction", "st") => :delete),
        c("fraction above 1", ["d_fraction"], ("cell", 1, "d_fraction", "reg") => "3/2"), c("fraction negative", ["d_fraction"], ("cell", 1, "d_fraction", "reg") => "-0.1"),
        c("fraction class unknown", ["d_fraction"], ("cell", 1, "d_fraction", "zz") => "1/2"), c("d_min negative", ["d_min"], ("cell", 1, "d_min") => -1),
        c("control odd dimension", ["dimension"], ("control", 1, "dimension") => 3), c("control mu wrong length", ["mu"], ("control", 1, "mu") => ["3"]),
        c("control mu entry 1", ["mu"], ("control", 1, "mu") => ["3", "1"]), c("chain_origin eps nonzero", ["eps"], ("control", 2, "eps") => "0.1"),
        c("seed above generator range", ["seed"], "seed" => 2^64)]
    return cases
end

@testset "D125b malformed panel inputs are refused naming the key, before any computation" begin
    @test EX.panel_inputs(d125_block(D125_FA)).scheme == "systems"                          # the unedited block is valid (RED today)
    cases = d125_panel_cases()
    @test length(cases) >= 40
    bad = [name for (name, keys, blk) in cases if !ex_rejects(() -> EX.panel_inputs(blk), keys)]
    isempty(bad) || @info "D125b: not refused as stated" bad
    @test isempty(bad)
    # 8 through the builder: refused, 0 files (panel, provenance sidecar, pool sidecar)
    pick = ["n_start 0", "scr_steps not a multiple", "N above 25", "energy at escape", "quota negative", "fraction above 1", "base key K", "class ids equal"]
    badb = [name for (name, keys, blk) in cases if name in pick && !d125_refused(blk, keys; tag = "b-" * replace(name, " " => "-"))]
    @test isempty(badb)
end

# ---------------------------------------------------------------- D125i run block grammar
# a valid block with the integrators and the estimators, every key written (I51)
d125_valid_block() = d125_run(Any[d125_arm("maps", "P", ["L1", "L5"], [100]; cells = ["y3a", "ch2"], vectors = "d", frames = ["a"]),
                                  d125_arm("flow", "P", ["L5"], [160]; cells = ["hh1"], vectors = 4, frames = ["a", "b2"], integrators = ["lf", "kl"],
                                           inds = ["tg"], ests = ["cq", "cl"])];
                              estimators = Any[d125_cqr(), d125_clv()], specs = Any[d125_tangent()], pilot_clv = 30)

# (name, key texts, edit of the valid block as a function); 52 cases, at least two per estimator kind
function d125_block_cases()
    v = d125_valid_block
    c = (name, keys, f) -> (name, keys, (d = v(); f(d); d))
    arm = (d, i) -> d["arm"][i]
    est = (d, i) -> d["estimator"][i]
    return Any[
        c("unknown arm key", ["thresh"], d -> arm(d, 1)["thresh"] = 1), c("cells empty", ["cells"], d -> arm(d, 1)["cells"] = String[]),
        c("cells unknown id", ["cells"], d -> arm(d, 1)["cells"] = ["nope"]), c("cells not a list", ["cells"], d -> arm(d, 1)["cells"] = "y3a"),
        c("vectors 1", ["vectors"], d -> arm(d, 1)["vectors"] = 1), c("vectors above d", ["vectors"], d -> arm(d, 1)["vectors"] = 5),
        c("vectors text other than d", ["vectors"], d -> arm(d, 1)["vectors"] = "all"), c("vectors missing for a QR arm", ["vectors"], d -> delete!(arm(d, 1), "vectors")),
        c("frames missing", ["frames"], d -> delete!(arm(d, 1), "frames")), c("frames repeated label", ["frames"], d -> arm(d, 2)["frames"] = ["a", "a"]),
        c("frames ill-formed label", ["frames"], d -> arm(d, 2)["frames"] = ["A-1"]), c("frames empty", ["frames"], d -> arm(d, 1)["frames"] = String[]),
        c("integrators missing on a flow arm", ["integrators"], d -> delete!(arm(d, 2), "integrators")),
        c("integrators unknown id", ["integrators"], d -> arm(d, 2)["integrators"] = ["zz"]),
        c("integrators on a map arm", ["integrators"], d -> arm(d, 1)["integrators"] = ["lf"]),
        c("estimators unknown id", ["estimators"], d -> arm(d, 2)["estimators"] = ["zz"]),
        c("controls arm with an estimator", ["estimators"], d -> (arm(d, 1)["orbits"] = "controls"; arm(d, 1)["estimators"] = ["cq"])),
        c("arm without methods, specs or estimators", ["methods"], d -> (delete!(arm(d, 1), "methods"); delete!(arm(d, 1), "cells"))),
        c("direction in a new-kind arm", ["direction"], d -> arm(d, 1)["direction"] = ["pi/4"]), c("eps in a new-kind arm", ["eps"], d -> arm(d, 1)["eps"] = [1e-8]),
        c("wrapped on a flow cell", ["coordinate_mode"], d -> arm(d, 2)["coordinate_mode"] = ["wrapped"]),
        c("L2 on a new cell", ["L2"], d -> arm(d, 1)["methods"] = ["L2"]), c("L9 on a new cell", ["L9"], d -> arm(d, 1)["methods"] = ["L9"]),
        c("horizon not a multiple of tau", ["horizons"], d -> (arm(d, 1)["tau"] = [3]; arm(d, 1)["horizons"] = [100])),
        c("integrator id duplicate", ["integrator"], d -> (d["integrator"][2]["id"] = "lf"; arm(d, 2)["integrators"] = ["lf"])), c("integrator unknown method", ["method"], d -> d["integrator"][1]["method"] = "Euler9"),
        c("integrator step not exact", ["step"], d -> d["integrator"][1]["step"] = "fast"), c("integrator step zero", ["step"], d -> d["integrator"][1]["step"] = "0"),
        c("integrator scheme missing where required", ["scheme"], d -> delete!(d["integrator"][1], "scheme")),
        c("integrator scheme where refused", ["scheme"], d -> d["integrator"][2] = Dict{String, Any}("id" => "kl", "method" => "Ruth3", "scheme" => "x", "step" => "0.1")),
        c("integrator id ill-formed", ["id"], d -> d["integrator"][1]["id"] = "1lf"),
        c("estimator id equals a method id", ["id"], d -> (est(d, 1)["id"] = "L5"; arm(d, 2)["estimators"] = ["L5", "cl"])),
        c("estimator id equals an integrator id", ["id"], d -> (est(d, 1)["id"] = "KL"; arm(d, 2)["estimators"] = ["KL", "cl"])),
        c("estimator kind unknown", ["kind"], d -> est(d, 1)["kind"] = "qr"), c("estimator budget missing", ["budget"], d -> delete!(est(d, 1), "budget")),
        c("cqr order 3", ["order"], d -> est(d, 1)["order"] = 3), c("cqr form unknown", ["form"], d -> est(d, 1)["form"] = "f99"),
        c("cqr f84 with k != d", ["form"], d -> (est(d, 1)["form"] = "f84"; est(d, 1)["vectors"] = 2)),
        c("cqr guard not positive", ["orthogonality_tol"], d -> est(d, 1)["orthogonality_tol"] = 0.0),
        c("cqr penalty factor text", ["penalty_factor"], d -> est(d, 1)["penalty_factor"] = "big"), c("cqr key of clv", ["window"], d -> est(d, 1)["window"] = 10),
        c("cqr on a map cell", ["estimators"], d -> arm(d, 1)["estimators"] = ["cq"]),
        c("cqr horizon not a multiple of m", ["horizons"], d -> (est(d, 1)["order"] = 4; arm(d, 2)["horizons"] = [161])),
        c("clv vectors != d", ["vectors"], d -> est(d, 2)["vectors"] = 3), c("clv orthonormalizer not a QR method", ["orthonormalizer"], d -> est(d, 2)["orthonormalizer"] = "L9"),
        c("clv checkpoint with block 0", ["block"], d -> est(d, 2)["storage"] = "checkpoint"), c("clv store with block 4", ["block"], d -> est(d, 2)["block"] = 4),
        c("clv window 0", ["window"], d -> est(d, 2)["window"] = 0), c("clv gap_exclude out of range", ["gap_exclude"], d -> est(d, 2)["gap_exclude"] = [4]),
        c("clv unstable_fraction 0", ["unstable_fraction"], d -> est(d, 2)["unstable_fraction"] = 0.0),
        c("clv probability 1", ["probabilities"], d -> est(d, 2)["probabilities"] = ["1/2", "1"]), c("clv probability not exact", ["probabilities"], d -> est(d, 2)["probabilities"] = ["half"]),
        c("clv hist_lo >= hist_hi", ["hist_lo"], d -> est(d, 2)["hist_lo"] = 1.0), c("clv hist_bins 0", ["hist_bins"], d -> est(d, 2)["hist_bins"] = 0),
        c("clv store_angles unknown", ["store_angles"], d -> est(d, 2)["store_angles"] = "all"), c("clv key missing", ["transient_cap"], d -> delete!(est(d, 2), "transient_cap")),
        c("clv with BigFloat", ["precisions"], d -> arm(d, 2)["precisions"] = ["BigFloat(256)"]),
        c("clv pilot window missing", ["clv_window"], d -> delete!(d["pilot"], "clv_window")),
        c("tangent gali_k above d", ["gali_k"], d -> (d["indicator"][1]["members"] = ["GALI"]; d["indicator"][1]["gali_k"] = 5)),
        c("tangent GALI with another member", ["GALI"], d -> (d["indicator"][1]["members"] = ["GALI", "SALI"]; d["indicator"][1]["gali_k"] = 3)),
        c("tangent GALI without gali_k", ["gali_k"], d -> d["indicator"][1]["members"] = ["GALI"]),
        c("tangent gali_k without GALI", ["gali_k"], d -> d["indicator"][1]["gali_k"] = 3),
        c("tangent OFLI on a map cell", ["OFLI"], d -> (d["indicator"][1]["members"] = ["OFLI"]; arm(d, 2)["cells"] = ["y3a"]; delete!(arm(d, 2), "integrators");
                                                      delete!(arm(d, 2), "estimators"))),
        c("tangent SALI not last", ["SALI"], d -> d["indicator"][1]["members"] = ["SALI", "FLI"]),
        c("spec kind pair on a new cell", ["pair"], d -> (d["indicator"] = Any[d86_spec("pair")]; arm(d, 2)["indicators"] = ["pr"])),
        c("unknown key in an estimator", ["thresh"], d -> est(d, 2)["thresh"] = 1.0)]
end

@testset "D125i valid blocks are read, hashed and canonical" begin
    pr = ex_write_prereg_run(ex_path("d125-i-valid-prereg.md"), d125_valid_block())
    rb = EX.read_run_block(pr)                                                              # RED: unknown key integrator today
    txt = EX.block_text(rb)
    @test EX.block_hash(rb) == ex_sha256(txt)
    @test occursin("Leapfrog", txt) && occursin("cqr", txt) && occursin("clv", txt)
    # the merge by id of the integrator and estimator tables (R34 item 5): a later entry replaces, new entries are appended
    a = Dict{String, Any}("integrator" => d125_integrators(), "estimator" => Any[d125_cqr()])
    b = Dict{String, Any}("integrator" => Any[Dict{String, Any}("id" => "lf", "method" => "Leapfrog", "scheme" => "KDK", "step" => "0.03125"),
                                                Dict{String, Any}("id" => "pf", "method" => "PEFRL", "step" => "0.1")],
                          "estimator" => Any[d125_cqr(order = 4), d125_clv()])
    m = EX.merge_run(a, b)
    @test [i["id"] for i in m["integrator"]] == ["lf", "kl", "pf"] && m["integrator"][1]["scheme"] == "KDK"
    @test [e["id"] for e in m["estimator"]] == ["cq", "cl"] && m["estimator"][1]["order"] == 4
    # a block without a key of the amendment is unchanged: the D86b and D94k golden block hashes are asserted in ex_d86_grammar.jl and ex_d94_grammar.jl
    @test EX.block_hash(EX.read_run_block(ex_fix(:small).prereg)) == include(joinpath(@__DIR__, "ex_d86_golden.jl")).block_hash
end

@testset "D125i malformed blocks are refused naming the key path, 0 files" begin
    @test d125_plan(d125_fix("valid", d125_valid_block())) isa EX.RunPlan                    # the unedited block plans (RED today)
    cases = d125_block_cases()
    @test length(cases) >= 50
    bad = [name for (name, keys, run) in cases if !d125_block_refused(run, keys; tag = replace(name, r"[^a-z0-9]+" => "-"))]
    isempty(bad) || @info "D125i: not refused as stated" bad
    @test isempty(bad)
    ek = [length(filter(c -> startswith(c[1], p), cases)) for p in ("cqr", "clv")]
    @test all(>=(2), ek)
end
