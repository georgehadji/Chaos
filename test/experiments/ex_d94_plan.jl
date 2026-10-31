# Amendment D94, plan and guard: D94i (D_B in the plan), D94l (plan with rounds), D94m (round log and guard).
# Spec: experiments-spec.md, "Amendment D94": E38 to E41, R25, R26, I32, I33, I37. API assumed: the test plan.
# Fixture: the hand panel of ex_d94_helpers.jl (16 orbits in 9 round cells, rounds 1 to 3, db_size 5), the block with the arms ref (P, each_round),
# dd (D, closure) and db (D_B, closure). The logs are written by the helpers with the digests of the store of the shared campaign (d94_campaign).

const D94_CONT = (decision = "continue",)
d94_spec(r, dec = "continue"; kw...) = (; round = r, decision = dec, kw...)
const D94_CUTS = [
    ("empty log", NamedTuple[], 1, false),
    ("continue at 1", [d94_spec(1)], 2, false),
    ("continue at 1 and 2", [d94_spec(1), d94_spec(2)], 3, false),
    ("close at 3", [d94_spec(1), d94_spec(2), d94_spec(3, "close")], 3, true),
    ("close at 2", [d94_spec(1), d94_spec(2, "close")], 2, true),
    ("close at 1", [d94_spec(1, "close")], 1, true),
]

function d94_planlog(c, specs; fx = c.fx, tag = "x")
    p = ex_path("d94-$(fx.tag)-$tag-$(abs(hash(specs))).log")
    write(p, d94_log_text(fx, specs, c.digest)[1])
    return p
end

# reference counts of a plan from the panel (hand tables of the fixture): items and chunks of each arm for rounds up to R, closure included or not
function d94_expect(R::Integer, closed::Bool)
    cells = [cl for cl in D94_CELLS if cl[1] <= R]
    ref = (items = sum(cl[4] for cl in cells), chunks = length(cells), per_round = [sum(cl[4] for cl in cells if cl[1] == r) for r in 1:R])
    DR = d94_hand_D(R)
    dd = (items = length(DR), chunks = length(unique(d94_hand_cell(x.id) for x in DR)))
    dbids = d94_db(EX_RUN_SEED, R, DR, 5)
    db = (items = length(dbids), chunks = length(unique(d94_hand_cell(id) for id in dbids)), ids = dbids)
    return (; ref, dd, db, items = ref.items + (closed ? dd.items + db.items : 0), chunks = ref.chunks + (closed ? dd.chunks + db.chunks : 0))
end

d94_sig(plan, ch) = (ch.ordinal, ch.key, ch.hash, [plan.items[i].hash for i in ch.items])

@testset "D94i D_B in the plan: size, subset, allocation per K index (hand tables)" begin
    c = d94_campaign()
    for (name, specs, R, closed) in D94_CUTS
        closed || continue
        plan = d94_plan(c.fx; log = d94_planlog(c, specs))
        DR = d94_hand_D(R)
        dbids = sort(unique(it.orbit for it in plan.items if "db" in it.arms))
        ddids = sort(unique(it.orbit for it in plan.items if "dd" in it.arms))
        @test (name, dbids) == (name, d94_db(EX_RUN_SEED, R, DR, 5))
        @test (name, length(dbids)) == (name, min(5, length(DR)))                       # |D_B| = min(db_size, |D_R|)
        @test issubset(dbids, [x.id for x in DR]) && (name, ddids) == (name, [x.id for x in DR])     # D_B within D_R: no orbit outside D_R
        n = [count(x -> x.k_index == i, DR) for i in 1:3]
        @test (name, [count(id -> d94_hand_cell(id)[2] == i, dbids) for i in 1:3]) == (name, d94_alloc(n, 5))
        @test (name, d94_alloc(n, 5)) == (name, R == 1 ? [1, 2, 0] : [1, 3, 1])           # the hand tables of the plan
    end
    # D_B for R = 2 differs from R = 3 only through R: the same D with another R in the stream text gives another selection
    @test EX.db_subset(EX_RUN_SEED, 2, d94_hand_D(3), 5) == d94_db(EX_RUN_SEED, 2, d94_hand_D(3), 5)
    # an open log releases no D_B: the arm db is planned only at closure
    open_plan = d94_plan(c.fx; log = d94_planlog(c, [d94_spec(1)]))
    @test !any(it -> "db" in it.arms || "dd" in it.arms, open_plan.items)
end

@testset "D94l plan with rounds: stable ordinals, keys and hashes; names; counts; order; determinism" begin
    c = d94_campaign()
    fx = c.fx
    plans = Dict(name => d94_plan(fx; log = d94_planlog(c, specs)) for (name, specs, R, closed) in D94_CUTS)
    for (name, specs, R, closed) in D94_CUTS
        p = plans[name]
        ex = d94_expect(R, closed)
        @test (name, p.rounds.R_plan, p.rounds.closed) == (name, R, closed)
        @test (name, p.counts.items, p.counts.chunks) == (name, ex.items, ex.chunks)            # the hand count per phase and round
        refs, clos = d94_ref_chunks(p), d94_clo_chunks(p)
        @test (name, [sum(length(ch.items) for ch in refs if d94_cell_round(ch.cell) == r) for r in 1:R]) == (name, ex.ref.per_round)
        @test (name, length(refs), sum(length(ch.items) for ch in refs)) == (name, ex.ref.chunks, ex.ref.items)
        @test all(ch -> occursin(r"^k\d+:[TRCSX]:r\d+$", ch.cell) && occursin("cell=" * ch.cell, ch.key), p.chunks)            # the cell names carry the round
        @test all(ch -> d94_cell_round(ch.cell) <= R, p.chunks)                                    # no round beyond R_plan
        if closed
            @test (name, length(clos), sum(length(ch.items) for ch in clos)) == (name, ex.dd.chunks + ex.db.chunks, ex.dd.items + ex.db.items)
            @test minimum(ch.ordinal for ch in clos) > maximum(ch.ordinal for ch in refs)           # the closure chunks come after all each_round chunks
            @test issorted([ch.ordinal for ch in p.chunks]) && [ch.ordinal for ch in p.chunks] == collect(1:length(p.chunks))
        else
            @test isempty(clos) && !any(it -> "dd" in it.arms || "db" in it.arms, p.items)         # no closure arm before closure
        end
        @test issorted([d94_cell_round(ch.cell) for ch in refs])                                    # round ascending inside the each_round phase
        @test occursin("R_plan=$R", p.text) && occursin("closed=$closed", p.text) && occursin(ex_sha256(read(d94_planlog(c, specs))), p.text)
        @test p.rounds.log_sha256 == ex_sha256(read(d94_planlog(c, specs)))
    end
    # stability (I33): ordinal, key and hash of every each_round chunk of the rounds in common are equal in every plan
    base = plans["close at 3"]
    for (name, specs, R, closed) in D94_CUTS
        p = plans[name]
        want = [d94_sig(base, ch) for ch in d94_ref_chunks(base) if d94_cell_round(ch.cell) <= R]
        @test (name, [d94_sig(p, ch) for ch in d94_ref_chunks(p)]) == (name, want)
        @test (name, [it.hash for it in p.items if "ref" in it.arms]) == (name, [it.hash for it in base.items if "ref" in it.arms && d94_hand_cell(it.orbit)[1] <= R])
    end
    # the plan bytes are identical for 1, 4 and all available threads
    for nt in (1, 4, 23, Threads.nthreads())
        @test d94_plan(fx; log = d94_planlog(c, D94_CUTS[4][2]), ntasks = nt).text == base.text
    end
    # a block that lists its keys in another order gives the same plan, but for the line of the log hash (the log carries the header hash of its preregistration)
    rev = d94_fix("rev", d94_block(); rev = true)
    for (name, specs, R, closed) in D94_CUTS[[1, 4]]
        pr = d94_plan(rev; log = d94_planlog(c, specs; fx = rev))
        strip_log(t) = replace(t, r"(?m)^log_sha256=.*\n" => "")
        @test strip_log(pr.text) == strip_log(plans[name].text) && [d94_sig(pr, ch) for ch in pr.chunks] == [d94_sig(plans[name], ch) for ch in plans[name].chunks]
    end
    # the log is read, never written (I37)
    lp = d94_planlog(c, D94_CUTS[3][2])
    before = read(lp)
    d94_plan(fx; log = lp)
    @test read(lp) == before
    # the plan reads a built panel: cell names and counts follow the panel file
    b = d94_build("A", D94_FA)
    run = d94_block([ex_arm("ref", "P", ["L1"], ["Float64"], [20]; release = "each_round"), ex_arm("dd", "D", ["L5"], ["Float64"], [20]; release = "closure")])
    run["panel_prereg_sha256_accepted"] = [b.panel["manifest"]["prereg"]["computed_sha256"]]
    pr = ex_path("d94-built-prereg.md")
    ex_write_prereg_run(pr, run)
    emptylog = ex_path("d94-built.log")
    write(emptylog, "")
    bp = EX.plan_run(pr, b.path; code_version = EX_CODE, round_log = emptylog)
    r1 = [o for o in b.panel["orbit"] if o["round"] == 1]
    @test bp.counts.items == length(r1) && Set(ch.cell for ch in bp.chunks) == Set("k$(o["k_index"]):$(o["class"]):r1" for o in r1)
    @test bp.rounds.R_plan == 1 && !bp.rounds.closed
end

@testset "D94m round log and guard" begin
    c = d94_campaign()
    fx = c.fx
    c1, c2, c3 = d94_spec(1), d94_spec(2), d94_spec(3)
    x1, x2, x3 = d94_spec(1, "close"), d94_spec(2, "close"), d94_spec(3, "close")
    cases = [
        ("empty", NamedTuple[], (1, false)),
        ("continue at 1", [c1], (2, false)),
        ("continue at 1 and 2", [c1, c2], (3, false)),
        ("close at 2", [c1, x2], (2, true)),
        ("close at r_max", [c1, c2, x3], (3, true)),
        ("close at 1", [x1], (1, true)),
        ("a gap", [c1, d94_spec(3)], 2),
        ("a repeated round", [c1, d94_spec(1)], 2),
        ("a changed panel hash", [d94_spec(1; panel = "a"^64)], 1),
        ("a changed preregistration hash", [c1, d94_spec(2; prereg = "b"^64)], 2),
        ("a broken chain", [c1, d94_spec(2; prev = "c"^64)], 2),
        ("a changed ref digest", [d94_spec(1; digest = "d"^64)], 1),
        ("an entry after close", [x1, d94_spec(2)], 2),
        ("continue at r_max", [c1, c2, d94_spec(3)], 3),
        ("an unknown decision", [d94_spec(1, "pause")], 1),
        ("a changed count is data: the digest of the entry text is the chain, not the counts", [d94_spec(1; counts = Dict("n_ts" => 1))], (2, false)),
    ]
    @test length(cases) >= 14
    for (name, specs, want) in cases
        lp = d94_planlog(c, specs; tag = "m")
        if want isa Tuple
            p = d94_plan(fx; log = lp)
            @test (name, p.rounds.R_plan, p.rounds.closed) == (name, want[1], want[2])
        else
            @test (name, d94_rej(() -> d94_plan(fx; log = lp), "entry $want")) == (name, true)
            store = ex_store("d94m-refused")
            @test (name, ex_rejects(() -> ex_campaign(fx; store = store, round_log = lp), ["entry $want"])) == (name, true)
            @test !ispath(store) && !isfile(store * "-manifest.json")                                      # 0 files written
        end
    end
    # a log file that does not exist and a line that is not an entry are refused naming the file or the line
    store = ex_store("d94m-nofile")
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = ex_path("d94-missing.log")), ["log"]) && !ispath(store)
    junk = ex_path("d94m-junk.log")
    write(junk, "this is not a log entry\n")
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = junk), ["entry 1", "log"]) && !ispath(store)
    # a closure arm asked for without a closed log, and a round beyond R_plan, are refused with 0 files
    open_log = d94_planlog(c, [c1]; tag = "open")
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = open_log, arms = ["dd"]), ["dd"]) && !ispath(store)
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = open_log, arms = ["db", "ref"]), ["db"]) && !ispath(store)
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = open_log, up_to = 3), ["round"]) && !ispath(store)
    @test ex_rejects(() -> ex_campaign(fx; store = store, round_log = d94_planlog(c, NamedTuple[]; tag = "e0"), up_to = 2), ["round"]) && !ispath(store)
    # the guard status through the entry point is nonzero (D94u)
    buf = IOBuffer()
    st = EX.cli(["run", fx.prereg, fx.panel, "--store", store, "--manifest", store * "-manifest.json", "--round-log", d94_planlog(c, [d94_spec(1; digest = "d"^64)]; tag = "m"),
                 "--thread-override", "--allow-dirty", "--code-version", EX_CODE, "--threads", "4", "--capacity", "64"]; io = buf)
    @test st == 1 && occursin("entry 1", String(take!(buf))) && !ispath(store)
end

@testset "D94 references: the log writer and the expected counts (guards, pass before the code exists)" begin
    fx = d94_std()
    txt, ents = d94_log_text(fx, [d94_spec(1), d94_spec(2, "close")], r -> "a"^64)
    @test D94_EMPTY_SHA == ex_sha256("") && length(ents) == 2 && txt == join(ents)
    @test occursin("prev_sha256=" * D94_EMPTY_SHA, ents[1]) && occursin("prev_sha256=" * ex_sha256(ents[1]), ents[2])
    @test occursin("decision=close\n", ents[2]) && occursin("count.n_ts=11\ncount.n_w2=4\n", ents[1]) && occursin("panel_sha256=" * ex_sha256(read(fx.panel)), ents[1])
    @test occursin("prereg_sha256=" * ex_header_hash(read(fx.prereg, String)), ents[1]) && startswith(ents[1], "round=1\ndecision=continue\n")
    e3, e2, e1 = d94_expect(3, true), d94_expect(2, true), d94_expect(1, true)
    @test (e3.ref.items, e3.ref.chunks, e3.dd.items, e3.dd.chunks, e3.db.items) == (16, 9, 9, 9, 5)
    @test (e2.ref.items, e2.ref.chunks, e2.dd.items, e2.db.items) == (11, 6, 6, 5) && (e1.ref.items, e1.dd.items, e1.db.items) == (6, 3, 3)
    @test e1.db.ids == [1, 3, 5] && e3.ref.per_round == [6, 5, 5] && d94_expect(2, false).items == 11 && d94_expect(3, false).chunks == 9
    @test [x.id for x in d94_hand_D(3)] == [1, 3, 5, 7, 9, 10, 12, 14, 16] && [x.k_index for x in d94_hand_D(3)] == [1, 2, 2, 2, 2, 3, 2, 2, 3]
    @test d94_digest(["b", "a"]) == ex_sha256("a\nb\n") && d94_cell_round("k2:T:r3") == 3
end
