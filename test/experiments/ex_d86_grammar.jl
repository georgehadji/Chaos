# Amendment D86, criteria D86a (grammar and validation), D86b (absence of change), D86c (plan), D86d (hashes and keys).
# Spec: experiments-spec.md, "Amendment D86": R11, R12, E22, E26, E27, I18, O19. API assumed: experiments-test-plan.md, "API assumed by the D86 tests".

# an ArgumentError that contains every one of the texts (a message that names only the old unknown key `indicator` is not enough)
function d86_rej(f, texts::AbstractString...)
    try
        f()
    catch e
        return e isa ArgumentError && all(t -> occursin(t, e.msg), texts)
    end
    return false
end

# edit a spec of a block by id: d86_ed(run, "st", "bins" => 0)
function d86_ed(run, id, kvs::Pair...)
    i = findfirst(s -> s["id"] == id, run["indicator"])
    return ex_rb_with(run, (("indicator", i, k) => v for (k, v) in kvs)...)
end
d86_ix(run, id) = findfirst(s -> s["id"] == id, run["indicator"])

@testset "D86a valid blocks read and give the canonical text" begin
    base = d86_all_block()
    rb = EX.run_block(base)
    @test length(rb.raw["indicator"]) == 7 && [s["id"] for s in rb.raw["indicator"]] == [s["id"] for s in base["indicator"]]
    txt = EX.block_text(rb)
    @test EX.block_hash(rb) == ex_sha256(txt)
    for (i, s) in enumerate(base["indicator"])                                   # one line per path, the canonical form of D53
        @test occursin("run.indicator[$i].id = \"$(s["id"])\"\n", txt) && occursin("run.indicator[$i].kind = \"$(s["kind"])\"\n", txt)
        @test occursin("run.indicator[$i].budget.Float64.safety = ", txt)
    end
    @test occursin(r"run\.arm\[1\]\.indicators = \[\"tg\", \"st\", \"pr\", \"prd\"\]\n", txt)
    # one spec of each kind alone
    for k in D86_KINDS
        s = d86_spec(k)
        wk = d86_window_kind(s)
        blk = d86_one(s; mode = k in ("frequency", "birkhoff") ? "wrapped" : nothing)
        r1 = EX.run_block(blk)
        @test r1.raw["indicator"][1]["kind"] == k && r1.raw["arm"][1]["indicators"] == [s["id"]]
        wk && @test r1.raw["arm"][1]["horizons"] == [d86_runlen(s)]
    end
    # an arm with `indicators` non-empty may omit `methods` or give an empty list (correction 7)
    withempty = d86_one(d86_spec("tangent"))
    withempty["arm"][1]["methods"] = String[]
    @test EX.run_block(withempty).raw["arm"][1]["indicators"] == ["tg"]
    # an arm may hold specs only; with methods and specs together both are kept
    @test !haskey(EX.run_block(d86_one(d86_spec("tangent"))).raw["arm"][1], "methods") || isempty(EX.run_block(d86_one(d86_spec("tangent"))).raw["arm"][1]["methods"])
    @test EX.run_block(base).raw["arm"][5]["methods"] == ["L1", "L5"]
    # key order of the file does not matter
    @test EX.block_text(EX.read_run_block(d86_fx("a-rev", base; rev = true).prereg)) == txt
end

@testset "D86a the fraction 0.1: exempt for an indicator-only arm (correction 1), still refused for an arm with a catalog method" begin
    # run length 8001 of phase advance with the windows 250, 1000, 4000, the fraction 0.1 in the grid: accepted, N / 10 is not an integer
    blk = d86_one(d86_spec("phase_advance"; windows = [250, 1000, 4000]))
    @test blk["grid"]["extra_fractions"] == [0.1, 0.5]
    @test EX.run_block(blk).raw["arm"][1]["horizons"] == [8001]
    # an arm that names a catalog method keeps the refusal: horizon 8001 with L1, with tau 1 (a guard, the refusal exists before the amendment)
    old = ex_rb_with(ex_rb(:small), ("arm", 1, "horizons") => [8001])
    @test d86_rej(() -> EX.run_block(old), "extra_fractions")
    # an arm with a method and a spec together is not exempt
    mixed = d86_all_block()
    mixed["arm"][5]["horizons"] = [101]
    @test d86_rej(() -> EX.run_block(mixed), "extra_fractions")
    # an indicator-only arm whose horizon is a multiple of 10 but N / 10 is not a multiple of tau is exempt as well
    tau4 = d86_one(d86_spec("tangent"); horizon = 100, tau = 4)
    @test EX.run_block(tau4).raw["arm"][1]["tau"] == [4]
end

@testset "D86a malformed blocks are refused with the key path" begin
    base = d86_all_block()
    E(id, kvs...) = d86_ed(base, id, kvs...)
    one(spec; kw...) = d86_one(spec; kw...)
    bigobs = [Dict{String, Any}("kind" => "cos", "coordinate" => 1, "frequency" => j) for j in 1:32]
    # (name, kind tag, block, texts that the message must hold)
    cases = Tuple{String, String, Any, Vector{String}}[
        ("missing id", "common", E("tg", "id" => :delete), ["indicator[1].id"]),
        ("id with capitals", "common", E("tg", "id" => "TG"), ["indicator[1].id"]),
        ("id starting with a digit", "common", E("tg", "id" => "1tg"), ["indicator[1].id"]),
        ("duplicate id", "common", E("pr", "id" => "tg"), ["indicator[2].id"]),
        ("id equal to a method id", "common", E("tg", "id" => "l1"), ["indicator[1].id"]),
        ("missing kind", "common", E("tg", "kind" => :delete), ["indicator[1].kind"]),
        ("unknown kind", "common", E("tg", "kind" => "gali"), ["indicator[1].kind"]),
        ("missing budget", "common", E("tg", "budget" => :delete), ["indicator[1].budget"]),
        ("budget without a precision of the arms", "common", ex_rb_with(base, ("indicator", 1, "budget", "Double64") => :delete), ["indicator[1].budget", "Double64"]),
        ("safety below 1", "common", ex_rb_with(base, ("indicator", 1, "budget", "Float64", "safety") => 0.5), ["indicator[1].budget", "safety"]),
        ("non-positive cost", "common", ex_rb_with(base, ("indicator", 1, "budget", "Float64", "ns_cpu_thread") => -1.0), ["indicator[1].budget", "ns_cpu_thread"]),
        ("unknown key threshold", "common", E("tg", "threshold" => 0.5), ["indicator[1].threshold"]),
        ("unknown key class", "common", E("st", "class" => "regular"), ["indicator[4].class"]),
        ("unknown key detect", "common", E("fq", "detect" => true), ["indicator[5].detect"]),
        ("tangent: members missing", "tangent", E("tg", "members" => :delete), ["indicator[1].members"]),
        ("tangent: members empty", "tangent", E("tg", "members" => String[]), ["indicator[1].members"]),
        ("tangent: unknown member", "tangent", E("tg", "members" => ["GALI"]), ["indicator[1].members"]),
        ("tangent: repeated member", "tangent", E("tg", "members" => ["FLI", "FLI", "SALI"]), ["indicator[1].members"]),
        ("tangent: SALI not last", "tangent", E("tg", "members" => ["SALI", "FLI"]), ["indicator[1].members"]),
        ("tangent: key of another kind", "tangent", E("tg", "separation" => 1e-8), ["indicator[1].separation"]),
        ("pair: separation zero", "pair", E("pr", "separation" => 0.0), ["indicator[2].separation"]),
        ("pair: separation negative", "pair", E("pr", "separation" => -1e-8), ["indicator[2].separation"]),
        ("pair: separation not finite", "pair", E("pr", "separation" => NaN), ["indicator[2].separation"]),
        ("pair: separation missing", "pair", E("pr", "separation" => :delete), ["indicator[2].separation"]),
        ("pair: unknown frame source", "pair", E("pr", "frame_source" => "random"), ["indicator[2].frame_source"]),
        ("pair: SALI member", "pair", E("pr", "members" => ["SALI"]), ["indicator[2].members"]),
        ("pair: repeated member", "pair", E("pr", "members" => ["RLI", "RLI"]), ["indicator[2].members"]),
        ("stretching: range lo >= hi", "stretching", E("st", "range" => [1.0, 1.0]), ["indicator[4].range"]),
        ("stretching: range of one entry", "stretching", E("st", "range" => [1.0]), ["indicator[4].range"]),
        ("stretching: bins zero", "stretching", E("st", "bins" => 0), ["indicator[4].bins"]),
        ("stretching: interval zero", "stretching", E("st", "interval" => 0), ["indicator[4].interval"]),
        ("stretching: interval not a multiple of tau", "stretching", one(d86_spec("stretching"; interval = 3); tau = 2), ["indicator[1].interval", "tau"]),
        ("stretching: min_per_cell negative", "stretching", E("st", "min_per_cell" => -1), ["indicator[4].min_per_cell"]),
        ("stretching: counts not a subset", "stretching", E("st", "counts" => "S_X"), ["indicator[4].counts"]),
        ("stretching: key of another kind", "stretching", E("st", "members" => ["SALI"]), ["indicator[4].members"]),
        ("frequency: windows empty", "frequency", E("fq", "windows" => Int[]), ["indicator[5].windows"]),
        ("frequency: windows not increasing", "frequency", E("fq", "windows" => [32, 16]), ["indicator[5].windows"]),
        ("frequency: window below the minimum", "frequency", E("fq", "windows" => [4, 16]), ["indicator[5].windows"]),
        ("frequency: window not an integer", "frequency", E("fq", "windows" => [16.5, 32]), ["indicator[5].windows"]),
        ("frequency: window below 8 (order + 1)", "frequency", E("fq", "order" => 2), ["indicator[5].windows"]),     # correction 7: 16 < 8 x 3
        ("frequency: order 6", "frequency", E("fq", "order" => 6), ["indicator[5].order"]),
        ("frequency: pad 1", "frequency", E("fq", "pad" => 1), ["indicator[5].pad"]),
        ("frequency: mirror_ratio above 1", "frequency", E("fq", "mirror_ratio" => 1.5), ["indicator[5].mirror_ratio"]),
        ("frequency: mirror_ratio zero", "frequency", E("fq", "mirror_ratio" => 0.0), ["indicator[5].mirror_ratio"]),
        ("frequency: floor_multiplier zero", "frequency", E("fq", "floor_multiplier" => 0.0), ["indicator[5].floor_multiplier"]),
        ("frequency: unknown observable", "frequency", E("fq", "observable" => "polar"), ["indicator[5].observable"]),
        ("frequency: unknown positions", "frequency", E("fq", "positions" => "middle"), ["indicator[5].positions"]),
        ("frequency: exp_i_q in lifted mode", "frequency", one(d86_spec("frequency"); mode = "lifted"), ["indicator[1].observable", "lifted"]),
        ("phase advance: window below 16", "phase_advance", E("pa", "windows" => [8, 16]), ["indicator[6].windows"]),
        ("phase advance: birkhoff_order 5", "phase_advance", E("pa", "birkhoff_order" => 5), ["indicator[6].birkhoff_order"]),
        ("phase advance: birkhoff_order 0", "phase_advance", E("pa", "birkhoff_order" => 0), ["indicator[6].birkhoff_order"]),
        ("phase advance: unknown observable", "phase_advance", E("pa", "observable" => "exp_i_q"), ["indicator[6].observable"]),
        ("phase advance: horizon is not tau L", "phase_advance", ex_rb_with(base, ("arm", 4, "horizons") => [64]), ["arm[4].horizons"]),
        ("birkhoff: no observables", "birkhoff", E("bk", "observables" => Dict{String, Any}[]), ["indicator[7].observables"]),
        ("birkhoff: 32 observables", "birkhoff", E("bk", "observables" => bigobs), ["indicator[7].observables"]),
        ("birkhoff: unknown observable kind", "birkhoff", E("bk", "observables" => [Dict{String, Any}("kind" => "tan", "coordinate" => 1, "frequency" => 1)]), ["indicator[7].observables[1]"]),
        ("birkhoff: frequency 0", "birkhoff", E("bk", "observables" => [Dict{String, Any}("kind" => "cos", "coordinate" => 1, "frequency" => 0)]), ["indicator[7].observables[1]"]),
        ("birkhoff: frequency not an integer", "birkhoff", E("bk", "observables" => [Dict{String, Any}("kind" => "cos", "coordinate" => 1, "frequency" => 1.5)]), ["indicator[7].observables[1]"]),
        ("birkhoff: coordinate 3", "birkhoff", E("bk", "observables" => [Dict{String, Any}("kind" => "sin", "coordinate" => 3, "frequency" => 1)]), ["indicator[7].observables[1]"]),
        ("birkhoff: cos in lifted mode", "birkhoff", one(d86_spec("birkhoff"); mode = "lifted"), ["indicator[1].observables", "lifted"]),
        ("birkhoff: window below 16", "birkhoff", E("bk", "windows" => [8, 16]), ["indicator[7].windows"]),
        ("arm with neither methods nor indicators", "arm", ex_rb_with(base, ("arm", 4, "indicators") => :delete), ["arm[4]", "neither"]),
        ("arm with an empty list of ids and no methods", "arm", ex_rb_with(base, ("arm", 4, "indicators") => String[]), ["arm[4]"]),
        ("arm with an unknown id", "arm", ex_rb_with(base, ("arm", 3, "indicators") => ["nosuch"]), ["arm[3].indicators", "nosuch"]),
        ("spec used by no arm", "arm", ex_rb_with(base, ("arm", 2, "indicators") => ["tg", "st", "prd", "pr", "fq"], ("arm", 6, "indicators") => ["tg", "fq"]), ["indicator[7]", "bk"]),
        ("window arm with two tau", "arm", ex_rb_with(base, ("arm", 4, "tau") => [1, 2]), ["arm[4]", "tau"]),
        ("window arm: horizon is not tau L", "arm", ex_rb_with(base, ("arm", 2, "horizons") => [100]), ["arm[2].horizons"]),
        ("BigFloat for an indicator arm", "arm", ex_rb_with(base, ("arm", 1, "precisions") => ["BigFloat"]), ["arm[1].precisions", "BigFloat"]),
        ("indicator id in control_methods", "arm", ex_rb_with(base, ("health", "control_methods") => ["L1", "tg"]), ["control_methods", "tg", "indicator"]),
    ]
    bad = String[]
    for (name, kind, blk, texts) in cases
        d86_rej(() -> EX.run_block(blk), texts...) || push!(bad, name)
    end
    @test isempty(bad)
    @test length(cases) >= 45                                   # the cases cover every refusal class of R11 item 6, several per class (correction 6b)
    for k in ("tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff")
        @test count(c -> c[2] == k, cases) >= 2
    end
    # file level: refused by the campaign before anything is written (0 files in the store, no manifest)
    for (i, (name, blk)) in enumerate([("repeated member", E("tg", "members" => ["FLI", "FLI", "SALI"])), ("threshold key", E("tg", "threshold" => 1.0)),
                                       ("windows not increasing", E("fq", "windows" => [32, 16])), ("unknown id", ex_rb_with(base, ("arm", 3, "indicators") => ["nosuch"])),
                                       ("pair separation", E("pr", "separation" => -1.0)), ("birkhoff order", E("bk", "birkhoff_order" => 7))])
        fx = d86_fx("a-bad$i", blk)
        store = ex_store("d86-a-bad$i")
        r = try
            ex_campaign(fx; store = store, manifest = store * "-manifest.json")
            nothing
        catch e
            e
        end
        @test r isa ArgumentError && occursin("indicator", r.msg) && occursin(r"\[\d\]|nosuch|\.arm|arm\[", r.msg)
        @test (!isdir(store) || isempty(readdir(store))) && !isfile(store * "-manifest.json")
    end
end

@testset "D86a merge of two run blocks by id, unknown keys, control_methods" begin
    base = d86_one(d86_spec("tangent"), d86_spec("stretching"))
    s2 = d86_spec("tangent"; members = ["MEGNO", "SALI"])
    new = d86_spec("pair"; id = "pr2")
    run5 = Dict{String, Any}("indicator" => [s2, new],
                             "arm" => [d86_arm("b", "S_L", ["pr2"], ["Float64"], [100])])
    pr = a5_prereg(ex_path("d86-a-merge-prereg.md"), base; run5 = run5)
    rb = EX.read_run_block(pr)
    ids = [s["id"] for s in rb.raw["indicator"]]
    @test ids == ["tg", "st", "pr2"]                                              # a replaced entry keeps its place, a new one is appended
    @test rb.raw["indicator"][1]["members"] == ["MEGNO", "SALI"]                  # replaced as a whole
    @test rb.raw["indicator"][3]["separation"] == 1e-8
    @test [a["name"] for a in rb.raw["arm"]] == ["a", "b"]
    # a spec replaced as a whole does not keep a key of the old entry that the new one lacks
    run6 = Dict{String, Any}("indicator" => [Dict{String, Any}("id" => "st", "kind" => "tangent", "members" => ["SALI"], "budget" => d86_budget())])
    rb6 = EX.read_run_block(a5_prereg(ex_path("d86-a-merge2-prereg.md"), base; run5 = run6))
    @test !haskey(only(s for s in rb6.raw["indicator"] if s["id"] == "st"), "bins")
end

# ---------------------------------------------------------------- D86b
# Golden values: computed once with the code version before this amendment (commit 38eee93) by the script in the plan, section D86b.
@testset "D86b absence of change: a block without specs" begin
    fx = ex_fix(:small)
    rb = EX.read_run_block(fx.prereg)
    @test !haskey(rb.raw, "indicator") && all(a -> !haskey(a, "indicators"), rb.raw["arm"])
    gold = include(joinpath(@__DIR__, "ex_d86_golden.jl"))
    @test EX.block_hash(rb) == gold.block_hash && ex_sha256(EX.block_text(rb)) == gold.block_hash
    plan = ex_plan(fx)
    @test plan.hash == gold.plan_hash && ex_sha256(EX.plan_text(plan)) == gold.plan_hash
    @test [it.hash for it in plan.items] == gold.item_hashes
    @test [c.key for c in plan.chunks] == gold.chunk_keys && [c.hash for c in plan.chunks] == gold.chunk_hashes
    @test all(it -> !hasproperty(it, :kind) || it.kind === nothing, plan.items) && !occursin("kind=", EX.item_text(first(plan.items)))
    # the chunk files of a campaign: canonical read-back SHA-256 per file equal to the golden values
    run = d86_run("b-small", ex_fix(:small).run)
    @test run.res.status == "complete"
    @test [ex_sha256(d86_canon_file(f)) for f in ex_chunk_files(run.store)] == gold.file_sha
end

# ---------------------------------------------------------------- D86c
@testset "D86c plan of the fixture block with all six kinds" begin
    run = d86_all_block()
    fx = d86_fx("c-all", run)
    plan = d86_plan("c-all", run)
    ref = d86_ref_keys(run, fx.orbits)
    isid(it) = it.method in keys(ref.perspec)
    sp = [it for it in plan.items if isid(it)]
    # 3 controls: the elliptic fixed point (13), the hyperbolic fixed point (14) and Y5 (15) (correction 6a)
    # hand count (S_D = ids 1 2 5 6 9 10 13 14 15, h4 = 16 17): ind1 4 specs x 9 = 36; ind2 tg, st, prd x 3 directions x 9 = 81, pr 9, fq and bk 8 each (Y5 omitted) = 106;
    # ind1dup 9 (all merged into ind1); ind3 9; mix L1 L5 on orbit 1 = 2 (its tg item is merged); dbl 3 specs x 2 starts = 6. Unique: 36 + 106 + 9 + 2 + 6 = 159.
    @test plan.counts.items == length(plan.items) == 159 && length(sp) == 157
    @test length(sp) == sum(v[1] for v in values(ref.perspec))
    @test Dict((it.orbit, it.method, ex_label(it.precision), it.horizon, it.tau, it.direction, it.coordinate_mode) => it.arms for it in sp) == ref.items
    @test [it.id for it in plan.items] == collect(1:159)
    # per arm items before the merge, from the plan's own report
    want_before = Dict("ind1" => 36, "ind2" => 106, "ind1dup" => 9, "ind3" => 9, "mix" => 3, "dbl" => 6)
    @test Dict(n => plan.steps.arms[n].items_before for n in keys(want_before)) == want_before
    @test ref.before == merge(want_before, Dict("mix" => 1))                         # the reference counts the spec items only (mix: 1 tg item)
    # directions: 3 values multiply tangent, stretching and pair-with-direction by 3, the window kinds and pair-with-seeded by 1
    ofspec(id) = [it for it in sp if it.method == id && "ind2" in it.arms]
    @test length(ofspec("tg")) == 27 && length(ofspec("st")) == 27 && length(ofspec("prd")) == 27 && length(ofspec("pr")) == 9 && length(ofspec("fq")) == 8 && length(ofspec("bk")) == 8
    @test all(it -> it.direction === nothing, [it for it in sp if it.method in ("fq", "pa", "bk", "pr")])
    @test all(it -> it.direction !== nothing, [it for it in sp if it.method in ("tg", "st", "prd")])
    @test all(it -> it.vectors == (it.method in ("fq", "pa", "bk") ? 1 : 2), sp) && all(it -> it.eps === nothing, sp)
    @test all(it -> it.kind == only(s["kind"] for s in run["indicator"] if s["id"] == it.method), sp)
    # applicability: the Y5 control (id 15) has no item for exp_i_q and cos / sin specs; the list is in the plan and in its text
    @test plan.not_applicable == Dict("fq" => [15], "bk" => [15]) == ref.omitted
    @test !any(it -> it.orbit == 15 && it.method in ("fq", "bk"), sp) && any(it -> it.orbit == 15 && it.method == "pa", sp)
    txt = EX.plan_text(plan)
    @test occursin(r"fq[^\n]*15", txt) && occursin(r"bk[^\n]*15", txt)
    # identical items of two arms are one item listing both, in arm order
    tg1 = only(it for it in plan.items if it.method == "tg" && it.orbit == 1 && it.horizon == 100)
    @test tg1.arms == ["ind1", "ind1dup", "mix"]
    # order: precision first, then the first arm; methods before specs inside an arm (mix: L1, L5)
    armidx = Dict(a["name"] => i for (i, a) in enumerate(run["arm"]))
    @test issorted(plan.items; by = it -> (EX_PREC_RANK[ex_label(it.precision)], armidx[first(it.arms)]))
    l1, l5 = (only(it for it in plan.items if it.method == m) for m in ("L1", "L5"))
    @test l1.id < l5.id
    groups = Dict{Any, Vector{Int}}()
    for it in plan.items
        push!(get!(groups, (EX_PREC_RANK[ex_label(it.precision)], armidx[first(it.arms)], ex_runtype(it)), Int[]), it.orbit)
    end
    @test all(issorted, values(groups))
    # per-spec step counts (horizon steps per item): tg 9 x 100 + 27 x 64 + 2 x 64; st and prd 9 x 100 + 27 x 64; pr 9 x 100 + 9 x 64; fq and bk 8 x 64 + 2 x 64; pa 9 x 65
    want = Dict("tg" => (38, 2756), "st" => (36, 2628), "prd" => (36, 2628), "pr" => (18, 1476), "fq" => (10, 640), "bk" => (10, 640), "pa" => (9, 585))
    @test Dict(id => (plan.steps.specs[id].items, plan.steps.specs[id].steps) for id in keys(want)) == want == Dict(id => Tuple(v) for (id, v) in ref.perspec)
    @test plan.steps.total_after == 11553 && plan.steps.total_before == 12553 && plan.steps.saving == 1000
    # bytes identical on repeated planning, for 1, 4 and all threads, and for the textual order of the tables
    @test plan.hash == ex_sha256(txt) && !occursin(EX_DIR, txt)
    for nt in (1, 4, Threads.nthreads())
        p2 = ex_plan(fx; ntasks = nt)
        @test EX.plan_text(p2) == txt && p2.hash == plan.hash
    end
    @test EX.plan_text(ex_plan(d86_fx("c-all", run; rev = true))) == txt
    # chunks: one per (run type, cell); the key of a new kind has the second form of R12 item 3
    cell = ex_cell_map(fx.orbits, run)
    @test plan.counts.chunks == length(unique((ex_runtype(it), cell[it.orbit]) for it in plan.items))
    @test sort(vcat((c.items for c in plan.chunks)...)) == collect(1:159)
    for c in plan.chunks
        it = plan.items[first(c.items)]
        it.method in keys(ref.perspec) || continue
        s = only(s for s in run["indicator"] if s["id"] == it.method)
        @test c.key == d86_chunk_key(it, d86_spec_hash(s)) && it.spec == d86_spec_hash(s)
        @test c.hash == ex_sha256(EX.chunk_text(plan, c.ordinal))
    end
end

@testset "D86c refusals: a window arm with another horizon or two tau, an arm that selects no orbit for a spec" begin
    run = d86_all_block()
    @test d86_rej(() -> EX.run_block(ex_rb_with(run, ("arm", 2, "horizons") => [64, 128])), "arm[2].horizons")
    @test d86_rej(() -> EX.run_block(ex_rb_with(run, ("arm", 4, "tau") => [1, 2])), "arm[4]", "tau")
    # tau = 2 with the horizon tau L is accepted
    ok = EX.run_block(d86_one(d86_spec("frequency"); tau = 2, mode = "wrapped"))
    @test ok.raw["arm"][1]["horizons"] == [128]
    # the only orbit of the arm (S_L = orbit 1, a drawn orbit) is fine; a subset in which the spec is undefined for every orbit is an error naming the spec
    blk = d86_one(d86_spec("frequency"); orbits = "S_L", mode = "wrapped")
    @test length(ex_plan(d86_fx("c-sl", blk)).items) == 1
    orbs = ex_panel_orbits(:std)
    for o in orbs
        o["in_SL"] = o["id"] == 15
    end
    pp = ex_write_panel(ex_path("d86-c-y5-panel.toml"); orbits = orbs)
    pr = ex_write_prereg_run(ex_path("d86-c-y5-prereg.md"), blk)
    @test d86_rej(() -> EX.plan_run(pr, pp; code_version = EX_CODE), "fq")
end

# ---------------------------------------------------------------- D86d
@testset "D86d hashes: one case per hashed field, per parameter, and what must not change them" begin
    run = d86_all_block()
    plan = d86_plan("c-all", run)
    base = Dict(k => first(it for it in plan.items if it.method == k) for k in ("tg", "pr", "prd", "st", "fq", "pa", "bk"))
    cases = 0
    for (k, it) in base
        h0 = EX.item_hash(it)
        fields = Pair{Symbol, Any}[:orbit => it.orbit + 100, :system => "linear", :param => "9.5", :x0 => "0.123", :p0 => "0.456", :run_seed => "1", :orbit_seed => "2",
            :method => it.method * "x", :precision => "Double64", :horizon => it.horizon + 1, :tau => it.tau + 1, :coordinate_mode => it.coordinate_mode == "lifted" ? "wrapped" : "lifted",
            :vectors => it.vectors + 1, :grid_sha256 => "f"^64, :code_version => "other", :hash_version => 2, :kind => "other", :spec => "e"^64]
        it.direction === nothing || push!(fields, :direction => "0.7")
        for (f, v) in fields
            @test EX.item_hash(merge(it, NamedTuple{(f,)}((v,)))) != h0
            cases += 1
        end
        @test !occursin("eps=", EX.item_text(it)) && occursin("\nkind=$(it.kind)\nspec=$(it.spec)\n", EX.item_text(it))
        @test endswith(EX.item_text(it), "hash_version=1\nkind=$(it.kind)\nspec=$(it.spec)\n")
        @test (it.direction === nothing) == !occursin("direction=", EX.item_text(it))
    end
    @test cases >= 24
    # every parameter of every kind changes the item hash and the spec hash (24 or more cases); budget, backend, arm name, key order do not
    id(s) = s["id"]
    edits = [("tg", "members" => ["MEGNO", "SALI"]), ("tg", "members" => ["SALI"]), ("pr", "separation" => 2e-8), ("pr", "frame_source" => "direction"), ("pr", "members" => ["MEGNO", "RLI"]),
             ("st", "range" => [-5.0, 5.0]), ("st", "bins" => 21), ("st", "interval" => 8), ("st", "min_per_cell" => 3), ("st", "counts" => "none"),
             ("fq", "observable" => "q_plus_i_p"), ("fq", "windows" => [16, 24, 32]), ("fq", "order" => 0), ("fq", "positions" => "edge"), ("fq", "pad" => 8),
             ("fq", "mean_removal" => false), ("fq", "mirror_ratio" => 0.6), ("fq", "floor_multiplier" => 16.0),
             ("pa", "observable" => "lifted"), ("pa", "windows" => [16, 24]), ("pa", "birkhoff_order" => 3),
             ("bk", "observables" => d86_obs()[1:3]), ("bk", "windows" => [24, 32]), ("bk", "birkhoff_order" => 3)]
    # a changed window moves the run length, so the arm horizon moves with it (R11 item 3)
    function edited(sid, kv)
        r = d86_ed(run, sid, kv)
        if kv.first == "windows"
            s = r["indicator"][d86_ix(r, sid)]
            for (ai, a) in enumerate(r["arm"])
                if sid in get(a, "indicators", String[]) && sid != "tg"
                    r = ex_rb_with(r, ("arm", ai, "horizons") => [d86_runlen(s)])
                end
            end
        end
        return r
    end
    for (j, (sid, kv)) in enumerate(edits)
        r = edited(sid, kv)
        p = ex_plan(d86_fx("d-edit$j", r))
        it0, it1 = (first(it for it in pl.items if it.method == sid) for pl in (plan, p))
        s0 = only(s for s in run["indicator"] if s["id"] == sid); s1 = only(s for s in r["indicator"] if s["id"] == sid)
        @test it1.spec != it0.spec && it1.spec == d86_spec_hash(s1) && it0.spec == d86_spec_hash(s0)
        @test it1.hash != it0.hash
        cases += 1
    end
    @test length(edits) >= 24
    # unchanged: budget, key order, arm name, backend request, thread count, launch capacity
    r2 = ex_rb_with(run, ("indicator", 1, "budget", "Float64", "ns_cpu_thread") => 999.0, ("indicator", 4, "budget", "Float64", "safety") => 7, ("arm", 1, "name") => "renamed",
                    ("backend", "Float64") => "cpu")
    p2 = ex_plan(d86_fx("d-same", r2))
    @test [it.hash for it in p2.items] == [it.hash for it in plan.items] && [it.spec for it in p2.items if it.kind !== nothing] == [it.spec for it in plan.items if it.kind !== nothing]
    @test [c.key for c in p2.chunks] == [c.key for c in plan.chunks]
    prev = ex_plan(d86_fx("c-all", run; rev = true))
    @test [it.hash for it in prev.items] == [it.hash for it in plan.items]
    for nt in (1, 4)
        @test [it.hash for it in ex_plan(d86_fx("c-all", run); ntasks = nt).items] == [it.hash for it in plan.items]
    end
    # the chunk key of an old method is unchanged (the first form of D53)
    c1 = ex_find_chunk(plan; method = "L1")
    it = plan.items[first(c1.items)]
    @test c1.key == "method=L1 precision=Float64 horizon=100 tau=1 eps=- direction=$(it.direction) mode=$(it.coordinate_mode) cell=$(c1.cell)"
end

@testset "D86d independent hashing (O19): 1000 items of every kind built from seeded field values" begin
    run = d86_all_block()
    plan = d86_plan("c-all", run)
    rng = StableRNG(d86_seed('d', 1))
    hex(n) = join(rand(rng, collect("0123456789abcdef"), n))
    # Every field of the item is drawn (correction 6c of the D86 stage 02 checkpoint): the items are not plan items of a fixture, they are built as NamedTuples with the
    # fields of E27 and passed to the module's own item_text and item_hash.
    function seeded_item(k::AbstractString, spec)
        id = spec["id"]
        has_dir = d86_has_dir(spec)
        return (orbit = rand(rng, 1:5000), system = rand(rng, ("standard", "linear")), param = string(round(5 * rand(rng); digits = 3)), x0 = string(6 * rand(rng)),
                p0 = string(6 * rand(rng)), run_seed = string(rand(rng, UInt64)), orbit_seed = string(rand(rng, UInt64)), method = id,
                precision = rand(rng, ("Float64", "Double64")), horizon = rand(rng, 1:2000), tau = rand(rng, (1, 2, 4)), eps = nothing,
                direction = has_dir ? rand(rng, ("pi/4", "0.3", string(round(3 * rand(rng); digits = 4)))) : nothing, coordinate_mode = rand(rng, ("lifted", "wrapped")),
                vectors = k in ("frequency", "phase_advance", "birkhoff") ? 1 : 2, grid_sha256 = hex(64), code_version = "code-" * hex(8), hash_version = 1,
                kind = k, spec = d86_spec_hash(spec))
    end
    for k in ("tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff")
        spec = first(s for s in run["indicator"] if s["kind"] == k)
        items = [seeded_item(k, spec) for _ in 1:1000]
        texts = EX.item_text.(items)
        @test texts == [d86_item_text(it, it.spec) for it in items]
        @test EX.item_hash.(items) == ex_sha256_many(texts)
        # the spec hash of the plan items is the independent one
        @test first(it for it in plan.items if it.method == spec["id"]).spec == d86_spec_hash(spec)
    end
    # the spec hash ignores the budget and the order of the keys (the text of a reversed file is the same)
    s = first(run["indicator"])
    s2 = deepcopy(s); s2["budget"]["Float64"]["ns_cpu_thread"] = 5.0
    @test d86_spec_hash(s) == d86_spec_hash(s2) && !occursin("budget", d86_spec_text(s))
    @test EX.spec_hash(EX.run_block(run), "tg") == d86_spec_hash(s) && EX.spec_text(EX.run_block(run), "tg") == d86_spec_text(s)
end
