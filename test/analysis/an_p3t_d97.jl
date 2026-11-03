# Amendment D97, stage 02: the test side of study p3-abstention reads the store of the runner. Spec: analysis-spec.md, Amendment D97, E74 to E82, I47 to I51, O28 to O30, criteria D97a to
# D97j (D97k is stage 04). The toy campaign of an_p3t_helpers.jl is written in the shape of the runner (E82, guarded by an_p3t_shape.jl); the code of Analysis reads the shape of the D95
# toy, so the testsets below are RED because the reader refuses or misreads the store (ArgumentError of the cell form, a missing arm, a missing tg item, wrong values), never because of a
# fixture. The expected values come from the design of the campaign through the reference of the D95 suite (p3t_expected, p3t_counts, p3t_log_entries), with the tags of the orbits
# changed where a testset changes the store (p3t_retag), or from the hand tables stated in each testset (and in the test plan, "Amendment D97").
# Every testset states which E-number it pins. No random draw is made here (the toy campaign draws with p3t_seed('a'); the statistics use the streams of E66).

d97_call(f) = try f() catch e; e isa InterruptException && rethrow(); e end         # an exception as a value, so that one RED call does not hide the assertions after it
d97_ok(r) = !(r isa Exception)
d97_err(r) = r isa Exception ? first(split(sprint(showerror, r), '\n')) : ""
const D97_CACHE = Dict{String, Any}()

# a section of the tables file (Table T<n>), its tables (header and rows, cells stripped), the panel cell of an orbit
d97_sec(tab, n) = (a = findfirst("Table T$n:", tab); b = n == 11 ? findfirst("## F2 data", tab) : findfirst("Table T$(n + 1):", tab); a === nothing ? "" : tab[a.start:(b === nothing ? lastindex(tab) : b.start - 1)])
function d97_tables(text)
    out = Vector{Vector{Vector{String}}}(); cur = Vector{Vector{String}}()
    for l in vcat(split(text, '\n'), [""])
        if startswith(l, "|")
            cells = String[String(strip(c)) for c in split(strip(l, '|'), '|')]
            all(c -> occursin(r"^:?-{2,}:?$", c), cells) || push!(cur, cells)
        elseif !isempty(cur)
            push!(out, cur); cur = Vector{Vector{String}}()
        end
    end
    return out
end
d97_panel_ids(panel, cell) = (m = match(r"^k(\d+):(\w+):r(\d+)$", cell); sort(Int[Int(o["id"]) for o in TOML.parsefile(panel)["orbit"] if o["k_index"] == parse(Int, m.captures[1]) && o["class"] == m.captures[2] && o["round"] == parse(Int, m.captures[3])]))
d97_file(ord) = "chunk-" * lpad(string(ord), 5, '0') * ".h5"
d97_ids(m, k) = Int.(m["sets"]["ids"][k])
d97_drops(m) = [(Int(d["id"]), d["reason"], Int(d["round"]), Int(d["k_index"])) for d in m["sets"]["drops"]]
d97_run(S, I, st, inst; kw...) = d97_call(() -> AN.analyze(; p3t_an_args(I, st, inst)..., threads = 2, kw...))
d97_read(res) = d97_ok(res) ? (m = AX.read_json(res.results), tab = read(res.tables, String)) : (m = nothing, tab = "")
# the rows of the table of the chunks failed twice (T11) that a store must have: the chunks of status failed of the campaign table and `extra` ordinals; orbits from the panel by cell (E80)
d97_failed_rows(S, extra = Int[]) = Set(Vector{String}[[d97_file(c.ordinal), c.cell, join(d97_panel_ids(S.panel, c.cell), ", ")] for c in S.chunks if c.status == "failed" || c.ordinal in extra])
d97_t11_rows(tab) = Set(r for t in d97_tables(d97_sec(tab, 11)) if first(t)[1] == "chunk file" for r in t[2:end])
d97_t3(tab) = (t = d97_tables(d97_sec(tab, 3)); isempty(t) ? Vector{Vector{String}}() : first(t))
d97_t10(tab) = Dict(r[1] => r for t in d97_tables(d97_sec(tab, 10)) for r in t[2:end])

# ---------------------------------------------------------------- D97a: the cell form of E74
@testset "D97a the cells of the store and the expected cells of the panel are those of E74 (the fixture side; passes now)" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    pan = TOML.parsefile(S.panel)
    want(r) = sort(unique(["k$(o["k_index"]):$(o["class"]):r$(o["round"])" for o in pan["orbit"] if o["round"] == r]))
    # by hand from the design (an_p3t_helpers.jl, P3T_LIST and P3T_RCLASS): the orbits 9, 13 and 28 are of class R, so a round and a K index can hold two cells; K index 1 is the K = 0 gate cell
    @test want(1) == ["k1:C:r1", "k2:C:r1", "k2:R:r1", "k3:C:r1", "k4:C:r1"] && want(2) == ["k2:C:r2", "k3:C:r2", "k4:C:r2", "k5:C:r2"] && want(3) == ["k2:C:r3", "k2:R:r3", "k3:C:r3", "k4:C:r3", "k5:C:r3"]
    for r in 1:3
        have = sort(unique([c.cell for c in S.chunks if c.arm == "ref" && c.round == r]))
        @test have == want(r) && length(symdiff(have, want(r))) == 0                           # 0 differences between the cells of the store and the cells of the panel
    end
    @test all(c -> occursin(r"^k\d+:[A-Z]:r\d+$", c.cell) && c.round == parse(Int, match(r":r(\d+)$", c.cell).captures[1]), S.chunks)
    # I51: every orbit of the panel is in exactly one cell of the ref arm, and the cells of an arm hold exactly the orbits of that arm
    for k in ("L1", "L5")
        ids = sort(vcat([c.ids for c in S.chunks if c.arm == "ref" && c.spec == k]...))
        @test ids == 1:31
    end
end

@testset "D97a count-round and analyze read the cells of the runner: the log is written round by round (E74)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    inst = p3t_instance(S, I, "d97a-r1"; rounds = 1, closure = false)
    r = d97_call(() -> AN.count_round(; p3t_cr_args(I, st, inst)...))
    @test d97_ok(r) && r.action === :appended && r.round == 1 && r.counts.orbits == 12 && r.counts.drift == 1 && r.counts.t_s == 2 && r.counts.w2 == 4     # hand counts of round 1 (plan, D95)
    # the cell of the class R holds orbits that the cell of the class C does not: dropping its L1 chunk leaves the round unfinished and the message names the cell of the panel
    inst = p3t_instance(S, I, "d97a-miss"; rounds = 1, closure = false)
    rm(p3t_chunk(S, inst; arm = "ref", role = :L1, cell = "k2:R:r1"))
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["k2:R:r1"]) && !isfile(inst.log)
    # a ref chunk with a cell of the old form (k<i>:r<round>), and one with a cell of no form: refused naming the cell (the discovery run showed that the old form never occurs in a store of the runner)
    for (tag, bad) in (("old", "k3:r1"), ("none", "kk3-x"))
        inst = p3t_instance(S, I, "d97a-$tag"; rounds = 1, closure = false)
        f = p3t_chunk(S, inst; arm = "ref", role = :L1, cell = "k3:C:r1")
        p3t_set_attr(f, "chunk_key", replace(AX.read_chunk_attrs(f)["chunk_key"], "cell=k3:C:r1" => "cell=$bad"))
        @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), [bad]) && !isfile(inst.log)
    end
    c = p3t_clone(p3t_full(), "d97a-an"); f = p3t_chunk(S, c; arm = "ref", role = :L5, cell = "k4:C:r2")
    p3t_set_attr(f, "chunk_key", replace(AX.read_chunk_attrs(f)["chunk_key"], "cell=k4:C:r2" => "cell=k4:r2"))
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["k4:r2"]) && !ispath(c.out)
end

# ---------------------------------------------------------------- D97b: tg at N_c (E75)
@testset "D97b the tg items sit at N_c in ind and d64 (the fixture side; passes now)" begin
    S = p3t_cached()
    tg = [c for c in S.chunks if c.spec == "tg"]
    @test !isempty(tg) && all(c -> c.horizon == P2_NC && c.arm in ("ind", "d64"), tg) && all(c -> c.horizon == P2_NC + 1, [c for c in S.chunks if c.spec == "pa1"]) &&
          all(c -> c.horizon == P2_NC, [c for c in S.chunks if c.spec in ("fq1", "bk1")]) && all(c -> c.horizon == P2_N1, [c for c in S.chunks if c.arm == "ref"])
    # one tg item of the toy decoded by hand: orbit 3 (SALI regular at the budgets 1 to 5 of its template) at the recorded steps of N_c: the last step is 20480 and no step is above it
    f = joinpath(S.store, d97_file(first(c.ordinal for c in S.chunks if c.arm == "ind" && c.spec == "tg" && 3 in c.ids)))
    t = AX.HDF5.h5open(h -> read(h["orbits"]["3"]["t"]), f, "r")
    @test t == P2_GB && t[end] == P2_NC && maximum(t) == 20480
end

@testset "D97b the verdicts of SALI, FLI, MEGNO read the tg item at N_c; |D| is its hand count, not 0" begin
    r = d97_call(p3t_run_main)
    @test d97_ok(r)
    d97_ok(r) || return
    (; X, m) = r
    cmp = Dict((c["set"] == "T_s" ? :Ts : :W2, Symbol(c["indicator"])) => c for c in m["hb"]["comparisons"])
    for (k, x) in X.hb.rows
        k[2] in (:SALI, :FLI, :MEGNO) || continue
        c = cmp[k]
        @test (Int(c["n11"]), Int(c["n10"]), Int(c["n01"]), Int(c["n00"])) == (x.n11, x.n10, x.n01, x.n00)
    end
    cnt(rs) = Dict(string(o) => count(r -> r.out === o, rs) for o in (:chaotic_supported, :regular_supported, :unresolved))
    @test Dict(k => Int(v) for (k, v) in m["outputs"]["E_test"]) == cnt(X.E) && Dict(k => Int(v) for (k, v) in m["outputs"]["T_s"]) == cnt([r for r in X.E if r.Ts])
    # |D| by hand: the 16 orbits of the panel flag in_D, minus the drift orbit 14 = 15; a reader that looks for tg at N_1 finds no Double64 verdict and gives 0
    @test m["sets"]["counts"]["D"] == 15 == length(X.D) && sum(Int(v) for v in values(m["r_ab"]["pairs"])) == 15 && Int.(m["r_ab"]["orbits"]) == X.D
    @test Int(m["hc"]["failures"]) + Int(m["hc"]["nonfailures"]) == 15
end

@testset "D97b a store with tg at N_1 only is refused, naming the item (E75, E61 k)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97b-n1")
    sp0 = first(s for s in p3t_allspecs(S.p) if s.arms == ["ind"] && s.id == "tg")
    for o in p3t_ords(S; arm = "ind", spec = "tg")
        S.chunks[o].status == "ok" || continue
        cell = S.chunks[o].cell
        p3t_alter_chunk(S, c, I, o; sp_edit = sp -> merge(sp, (; horizon = P2_N1, grid = P2_GA)),
                        edit = Dict{String, Any}("horizon" => P2_N1, "chunk_key" => p3t_key(merge(sp0, (; horizon = P2_N1)), cell)))
    end
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["tg"]) && !ispath(c.out)
end

# ---------------------------------------------------------------- D97c: the manifest status (E76)
function d97_healed()
    get!(D97_CACHE, "heal") do
        S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
        failed = [c.ordinal for c in S.chunks if c.status == "failed"]
        S2 = p3t_retag(S, [24, 26], :n)                                                  # the orbits of the failed chunks have their items in the healed store
        inst = p3t_instance(S, I, "d97-heal")
        p3t_heal_chunks(S, inst, I, failed)
        write(inst.log, join(first(p3t_log_entries(S2, I))))
        res = d97_run(S, I, st, inst)
        (; S2, inst, failed, res, X = p3t_expected(S2, I.qt), d97_read(res)...)
    end
end

@testset "D97c both forms of the status are accepted when they say what the entries say" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    R = d97_call(p3t_run_main)
    @test d97_ok(R) && AX.read_json(R.inst.manifest)["status"] == "complete with failed chunks" && count(e -> e["status"] == "failed", AX.read_json(R.inst.manifest)["chunks"]) == 6
    H = d97_healed()
    man = AX.read_json(H.inst.manifest)
    @test man["status"] == "complete" && all(e -> e["status"] == "ok", man["chunks"]) && length(H.failed) == 6
    @test d97_ok(H.res) && H.res.status === :ok
    d97_ok(H.res) || return
    m, X = H.m, H.X
    @test d97_ids(m, "E_test") == [r.id for r in X.E] && 24 in d97_ids(m, "E_test") && 26 in d97_ids(m, "E_test") && d97_drops(m) == [(d.id, d.reason, d.round, d.kidx) for d in X.drops]
    @test [m["decisions"][h]["verdict"] for h in ("HA", "HB", "HC")] == [X.ha.decision, X.hb.decision, X.hc.decision] && isempty(d97_t11_rows(H.tab))
    # the L1 cost of the healed store has no failed entry to exclude; the same hand value as the store with failures (the failed entry is excluded there)
    @test isapprox(parse(Float64, d97_t10(H.tab)["L1"][2]), 250.0; rtol = 1e-12)
end

@testset "D97c the three mismatched forms and any other word are refused, naming the status (E76)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    H = d97_healed()
    set_status(c, s) = (m = AX.read_json(c.manifest); m["status"] = s; AX._write_json(c.manifest, m))
    # `complete` with a failed entry; `complete with failed chunks` with none
    c = p3t_clone(p3t_full(), "d97c-a"); set_status(c, "complete")
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["complete", "failed"]) && !ispath(c.out)
    c = p3t_clone(H.inst, "d97c-b"); set_status(c, "complete with failed chunks")
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["complete with failed chunks"]) && !ispath(c.out)
    # other words, including the other statuses of the runner (D53 item 6) and an absent status
    for (n, w) in enumerate(("running", "stopped: failed chunks", "paused arms", "complete with failed", ""))
        c = p3t_clone(p3t_full(), "d97c-w$n"); set_status(c, w)
        @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["status", w]) && !ispath(c.out)
    end
    # the set of the failed chunks of the manifest is the set of the failed chunk files (the per-chunk check of E61 i, unchanged): a chunk failed in the manifest only
    c = p3t_clone(p3t_full(), "d97c-set")
    m = AX.read_json(c.manifest); e = first(e for e in m["chunks"] if e["status"] == "ok"); e["status"] = "failed"; AX._write_json(c.manifest, m)
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["status", d97_file(e["ordinal"])]) && !ispath(c.out)
end

# ---------------------------------------------------------------- D97d: failed twice in the arms (E77, E80 T11)
function d97_multi()
    get!(D97_CACHE, "multi") do
        S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
        inst = p3t_clone(p3t_full(), "d97-multi")
        # one more failed chunk in ind_pa1 (cell k3:C:r2, orbits 17 18 22), d64 (tg, cell k2:C:r1, orbits 3 4 5), d64_pa1 (cell k3:C:r1, orbits 6 7 8 14), chk (Double64, cell k3:C:r1, orbits 6 8 of D_B);
        # the main store has the ref cell k5:C:r2 (orbit 24) and the ind cells k5:C:r3 (orbit 26) failed
        fails = vcat(p3t_ords(S; arm = "ind_pa1", cell = "k3:C:r2"), p3t_ords(S; arm = "d64", spec = "tg", cell = "k2:C:r1"), p3t_ords(S; arm = "d64_pa1", cell = "k3:C:r1"),
                     p3t_ords(S; arm = "chk", cell = "k3:C:r1", precision = "Double64"))
        p3t_fail_chunks(S, inst, fails)
        S3 = p3t_retag(S, [17, 18, 22], :nov)
        X3 = p3t_expected(S3, I.qt; noD = [3, 4, 5, 6, 7, 8, 14])
        res = d97_run(S, I, st, inst)
        (; S, I, inst, fails, S3, X3, res, d97_read(res)...)
    end
end

@testset "D97d the main store: the ref arm and the ind arms failed twice (orbit 24 excluded, orbit 26 no verdict data), T11 lists the chunks with the orbits of their cells" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; S, X, m, tab) = R
    @test d97_drops(m) == [(14, "drift", 1, 3), (23, "unresolved", 2, 4), (24, "excluded", 2, 5), (26, "no verdict data", 3, 5)] == [(d.id, d.reason, d.round, d.kidx) for d in X.drops]
    rows = d97_t11_rows(tab)
    @test rows == d97_failed_rows(S) && length(rows) == 6
    @test count(r -> r[2] == "k5:C:r3" && r[3] == "26", rows) == 4 && count(r -> r[2] == "k5:C:r2" && r[3] == "24", rows) == 2        # the cells of ind (3 specs) and ind_pa1 (1), of ref (L1 and L5)
end

@testset "D97d one failed chunk in each of ind_pa1, d64, d64_pa1 and chk: the exclusions, the orbits that leave D, T11, T3 (E77, P118)" begin
    R = d97_multi()
    @test d97_ok(R.res) && R.res.status === :ok
    d97_ok(R.res) || return
    (; X3, m, tab, S, fails) = R
    @test length(fails) == 4 && AX.read_json(R.inst.manifest)["status"] == "complete with failed chunks"
    # ind_pa1: the orbits of its failed cell have no verdict data (17 18 22 in the order of the panel); the old exclusions stay
    @test d97_drops(m) == [(d.id, d.reason, d.round, d.kidx) for d in X3.drops] && [d[1] for d in d97_drops(m) if d[2] == "no verdict data"] == [17, 18, 22, 26]
    # d64 and d64_pa1: their orbits leave D (and only D): 3 4 5 and 6 7 8 (14 was dropped as drift); D has 15 - 3 - 3 - 1 = 8 orbits left of the 15, orbit 18 left E_test with ind_pa1
    @test d97_ids(m, "D") == X3.D == [9, 12, 13, 15, 16, 20, 21, 25] && m["sets"]["counts"]["D"] == 8 && all(i -> i in d97_ids(m, "E_test"), [3, 4, 5, 6, 7, 8])
    @test Set([3, 12]) <= Set(d97_ids(m, "D_B")) <= Set(P3T_DB) && Int(m["hc"]["failures"]) == X3.hc.failures && m["decisions"]["HC"]["verdict"] == X3.hc.decision &&
          m["decisions"]["HA"]["verdict"] == X3.ha.decision && m["decisions"]["HB"]["verdict"] == X3.hb.decision
    # T11: one row for each failed chunk, with the orbits of the panel cell
    @test d97_t11_rows(tab) == d97_failed_rows(S, fails) && length(d97_t11_rows(tab)) == 10
    # chk (P118): nothing is refused; the orbits 3 and 12 (their cells did not fail) keep their three values; the orbits 6 and 8 have no Double64 value: their row is absent or shows the missing
    # cell, never a number (reading: the spec says "its T3 row shows the missing cell" and fixes neither how D_B is formed for a failed chunk, since the orbits of the cell are those of the panel
    # cell 6 7 8 14 and not those of D_B, nor whether D_B keeps the orbits; both are accepted here, see the plan)
    t3 = d97_t3(tab)
    row(i) = (k = findfirst(r -> r[1] == string(i), t3[2:end]); k === nothing ? nothing : t3[k + 1])
    @test all(i -> row(i) !== nothing && all(j -> tryparse(Float64, row(i)[j]) !== nothing, 3:5), (3, 12)) && all(i -> row(i) === nothing || (tryparse(Float64, row(i)[4]) === nothing && tryparse(Float64, row(i)[3]) !== nothing), (6, 8))
    @test all(r -> r[1] in string.(P3T_DB), t3[2:end])
end

@testset "D97d a missing item of an orbit whose chunk is not listed as failed stays a refusal (E61 k)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97d-pa1")                         # the pa1 item of orbit 12 dropped from the ind_pa1 chunk of its cell (k4:C:r1), no failure in the manifest
    p3t_alter_chunk(S, c, I, first(p3t_ords(S; arm = "ind_pa1", cell = "k4:C:r1")); drop = [12])
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["pa1", "12"]) && !ispath(c.out)
end

# ---------------------------------------------------------------- D97e and O30: the measured costs (E78)
@testset "D97e the costs of the specs: tg, fq1, bk1 from ind, pa1 from ind_pa1, L1 from the chunk entries of the ref arm (hand values)" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; tab, m, X, I) = R
    t10 = d97_t10(tab)
    num(s) = tryparse(Float64, s)
    # P2_COST of the fixture: tg 237.5, pa1 1810.5, bk1 3000.25, fq1 16850.5 in the Float64 arms (the specs of d64 and d64_pa1 hold ten times these values and are not read, P119); L1: 1e9 s^-1 x
    # (sum of thread_seconds) / (sum of steps) of the ok L1 entries = 250.0 (every ok chunk of the toy runs at 250 ns per step); the failed L1 entry has 7.0 s over 240000 steps and is excluded
    @test num(t10["L1"][2]) !== nothing && isapprox(num(t10["L1"][2]), 250.0; rtol = 1e-12) && isapprox(num(t10["tg"][2]), 237.5; rtol = 1e-12) && isapprox(num(t10["pa1"][2]), 1810.5; rtol = 1e-12) &&
          isapprox(num(t10["bk1"][2]), 3000.25; rtol = 1e-12) && isapprox(num(t10["fq1"][2]), 16850.5; rtol = 1e-12)
    # the failed entry is the one with the large values: with it the cost would be (sum + 7.0) / (steps + 240000)
    es = [e for e in AX.read_json(R.inst.manifest)["chunks"] if startswith(e["key"], "method=L1 ") && e["arms"] == ["ref"]]
    @test count(e -> e["status"] == "failed", es) == 1 && length(es) == 14 && isapprox(1e9 * sum(e["thread_seconds"] for e in es if e["status"] == "ok") / sum(e["steps"] for e in es if e["status"] == "ok"), 250.0; rtol = 1e-12)
    sv = m["sensitivities"]["measured"]
    mv = p3t_hb_variant(X, I.qt; costs = (; cost_l1 = 250.0, cost_tg = 237.5, cost_pa1 = 1810.5, cost_bk1 = 3000.25, cost_fq1 = 16850.5))
    @test sv["decidable"] == mv.decidable && sv["verdict"] == mv.verdict
end

function d97_costs(tag, edit)
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), tag)
    m = AX.read_json(c.manifest); edit(m); AX._write_json(c.manifest, m)
    res = d97_run(S, I, st, c)
    return (; c, res, d97_read(res)...)
end

@testset "D97e a manifest whose measured costs equal the cost keys: the sensitivity (iii) is decidable and equals the base; O30 a: the L1 cost of a three-entry set with a failed entry and zero-step entries" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    X = p3t_expected(S, I.qt)
    # the cost keys of the toy table: tg 2658.0, pa1 1810.5, bk1 3000.6, fq1 16850.4; L1 by hand from the entries: three ok entries (1000 steps 0.5 s, 3000 steps 1.5 s, 2000 steps 1.0 s), every other ok entry has
    # zero steps and zero seconds, the failed entry has 240000 steps and 7.0 s and is excluded: 1e9 x (0.5 + 1.5 + 1.0) / (1000 + 3000 + 2000) = 5.0e5
    R = d97_costs("d97e-ok", function (m)
        for (sp, v) in (("tg", 2658.0), ("fq1", 16850.4), ("bk1", 3000.6)); m["arms"]["ind"]["specs"][sp]["measured_ns_per_orbit_step"] = v; end
        m["arms"]["ind_pa1"]["specs"]["pa1"]["measured_ns_per_orbit_step"] = 1810.5
        hand = [(1000, 0.5), (3000, 1.5), (2000, 1.0)]
        n = 0
        for e in m["chunks"]
            (startswith(e["key"], "method=L1 ") && e["arms"] == ["ref"] && e["status"] == "ok") || continue
            n += 1
            e["steps"], e["thread_seconds"] = n <= 3 ? hand[n] : (0, 0.0)
        end
    end)
    @test d97_ok(R.res)
    d97_ok(R.res) || return
    t10 = d97_t10(R.tab)
    @test isapprox(parse(Float64, t10["L1"][2]), 5.0e5; rtol = 1e-12) && isapprox(parse(Float64, t10["tg"][2]), 2658.0; rtol = 1e-12) && isapprox(parse(Float64, t10["fq1"][2]), 16850.4; rtol = 1e-12)
    mv = p3t_hb_variant(X, I.qt; costs = (; cost_l1 = 5.0e5, cost_tg = 2658.0, cost_pa1 = 1810.5, cost_bk1 = 3000.6, cost_fq1 = 16850.4))
    sv = R.m["sensitivities"]["measured"]
    @test sv["decidable"] == mv.decidable && sv["verdict"] == mv.verdict && Set(Symbol.(sv["failing"] .|> s -> replace(s, "T_s" => "Ts"))) == Set(Symbol("$(k[1]):$(k[2])") for k in mv.failing)
end

@testset "D97e a spec without a value, or an L1 with no steps, makes sensitivity (iii) not decidable naming it; only the Float64 arms are read (P119)" begin
    # tg missing from the specs of ind (the specs of d64 still hold it: they must not be used), and every ok L1 entry with zero steps
    R = d97_costs("d97e-miss", function (m)
        delete!(m["arms"]["ind"]["specs"], "tg")
        for e in m["chunks"]
            (startswith(e["key"], "method=L1 ") && e["arms"] == ["ref"] && e["status"] == "ok") && (e["steps"] = 0)
        end
    end)
    @test d97_ok(R.res)
    d97_ok(R.res) || return
    sv = R.m["sensitivities"]["measured"]
    @test sv["decidable"] == false && sv["verdict"] == "not decidable" && occursin("tg", sv["reason"]) && occursin("L1", sv["reason"]) && !occursin("fq1", sv["reason"]) && !occursin("pa1", sv["reason"])
    t10 = d97_t10(R.tab)
    @test tryparse(Float64, t10["tg"][2]) === nothing && tryparse(Float64, t10["L1"][2]) === nothing && isapprox(parse(Float64, t10["pa1"][2]), 1810.5; rtol = 1e-12)
    @test R.m["decisions"]["HB"]["verdict"] == "not supported"                                  # no sensitivity changes a decision
end

# ---------------------------------------------------------------- D97f: T3 (E79)
@testset "D97f T3: the columns of E79, the three L5 values found by content at the chk horizon, the ladder flip, no column `labels agree`" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; tab, S, X) = R
    t3 = d97_t3(tab); sec3 = d97_sec(tab, 3)
    @test !isempty(t3) && t3[1] == ["orbit", "N_B", "L5 Float64", "L5 Double64", "L5 BigFloat", "ladder would force", "rho"]
    @test !any(c -> occursin("labels agree", lowercase(c)), t3[1]) && occursin("labels agree", lowercase(sec3)) && occursin("omitted", lowercase(sec3)) && occursin("chk", sec3)      # the text states the omission and the cause
    rows = Dict(parse(Int, r[1]) => r for r in t3[2:end])
    @test sort(collect(keys(rows))) == P3T_DB && all(r -> r[2] == string(P3T_NB), values(rows))                # N_B = 1000 of the P1 table, equal to the chk horizon
    D = S.design.f
    for id in P3T_DB
        off = P3T_CHK_OFF[id]
        v = [parse(Float64, rows[id][j]) for j in 3:5]
        @test isapprox(v, [D.lref[id], D.lref[id] + off[1], D.lref[id] + off[2]]; rtol = 1e-12)                  # Float64 from ref at the recorded step N_B, the others from chk by content
        @test rows[id][6] == string(id == 6) && rows[id][7] == string(X.byid[id].rho)                               # the ladder flips on orbit 6 only (0.02 > 5e-3); orbit 8 differs by 0.004 <= 5e-3
    end
    # the count of the orbits on which the ladder would have made rho unresolved, split by rho
    m = match(r"unresolved: C (\d+), R (\d+)", sec3)
    @test m !== nothing && (parse(Int, m.captures[1]), parse(Int, m.captures[2])) == (X.byid[6].rho === :C ? (1, 0) : (0, 1))
end

@testset "D97f a chk horizon different from the key N_B of the P1 table is refused (E79)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97f-hor")
    for o in p3t_ords(S; arm = "chk")
        f = joinpath(c.store, d97_file(o))
        p3t_set_attr(f, "horizon", 2000); p3t_set_attr(f, "chunk_key", replace(AX.read_chunk_attrs(f)["chunk_key"], "horizon=$(P3T_NB)" => "horizon=2000"))
    end
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["horizon", "N_B"]) && !ispath(c.out)
end

# ---------------------------------------------------------------- D97g: the rounds above the closing round (E80, O30 c)
@testset "D97g a panel of 3 rounds closed at round 2: T2 and the counts are over the rounds 1 and 2 (hand counts)" begin
    S, st = p3t_cached(), p3_cached("d95")
    I = p3t_inputs(S, st, "d97g"; target_t_s = 3, target_w2 = 5)                       # decisions continue, close (D95d)
    ents, decs = p3t_log_entries(S, I)
    @test decs == ["continue", "close"]
    inst = p3t_instance(S, I, "d97g"; rounds = 2, clo_rounds = 2); write(inst.log, join(ents))                        # the runner releases the rounds up to the closing round only
    @test maximum(parse(Int, match(r":r(\d+)$", e["key"]).captures[1]) for e in AX.read_json(inst.manifest)["chunks"]) == 2 && length(AX.read_json(inst.manifest)["chunks"]) < S.files
    S2 = merge(S, (; orbits = [o for o in S.orbits if o.round <= 2]))
    X2 = d97_call(() -> p3t_expected(S2, I.qt))
    @test d97_ok(X2)
    res = d97_run(S, I, st, inst)
    @test d97_ok(res) && res.status === :ok
    (d97_ok(res) && d97_ok(X2)) || return
    m, tab = d97_read(res)
    # hand counts: round 1 has 14 orbits (ids 1 to 14: 2 of K = 0, drift 14), round 2 has 10 (ids 15 to 24: U 23, excluded 24); E_test 11 + 8 = 19, the gate 2: 19 + 2 + 3 drops = 24 orbits of the rounds 1 and 2
    @test m["sets"]["counts"]["E_test"] == 19 == length(X2.E) && m["sets"]["counts"]["gate"] == 2 && d97_drops(m) == [(14, "drift", 1, 3), (23, "unresolved", 2, 4), (24, "excluded", 2, 5)] == [(d.id, d.reason, d.round, d.kidx) for d in X2.drops]
    @test length(m["sets"]["drops"]) + m["sets"]["counts"]["gate"] + m["sets"]["counts"]["E_test"] == 24 && all(d -> Int(d["round"]) <= 2, m["sets"]["drops"]) && all(r -> Int(r["round"]) <= 2, m["sets"]["by_k_round"])
    @test [m["sets"]["counts"][k] for k in ("W2", "T_s", "S2", "D", "D_B")] == [count(r -> r.W2, X2.E), count(r -> r.Ts, X2.E), count(r -> r.S2, X2.E), length(X2.D), 4] && d97_ids(m, "E_test") == [r.id for r in X2.E] && !any(i -> i >= 25, d97_ids(m, "E_test"))
    @test !(25 in d97_ids(m, "D")) && length(m["log"]["entries"]) == 2 && m["log"]["entries"][2]["decision"] == "close"
    # T2: one row for each round up to the closing round, by hand: round, orbits, K = 0, E_test (the other columns from the reference)
    t2 = first(d97_tables(d97_sec(tab, 2)))
    @test t2[1][1:4] == ["round", "orbits", "K = 0", "E_test"] && [r[1:4] for r in t2[2:end]] == [["1", "14", "2", "11"], ["2", "10", "0", "8"]]
    for (i, k) in enumerate(1:2)
        rr = [o.id for o in S2.orbits if o.round == k]
        es = [r for r in X2.E if r.id in rr]
        @test (parse(Int, t2[i + 1][5]), parse(Int, t2[i + 1][6]), parse(Int, t2[i + 1][7])) == (count(r -> r.W2, es), count(r -> r.Ts, es), count(r -> r.S2, es))
    end
    @test occursin("drift on L1 or L5 (as the round log counts) 1", d97_sec(tab, 2)) && occursin("excluded 1,", d97_sec(tab, 2))
end

# ---------------------------------------------------------------- D97h: the six arms and the texts (E81, P114, P115, P120)
@testset "D97h the manifest must name the six arms; a manifest of the old four names is refused" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    for a in P3T_ARMS
        c = p3t_clone(p3t_full(), "d97h-$a")
        m = AX.read_json(c.manifest); delete!(m["arms"], a); AX._write_json(c.manifest, m)
        r = d97_call(() -> AN.analyze(; p3t_an_args(I, st, c)...))
        @test r isa ArgumentError && occursin(Regex("\\b$(a)\\b"), sprint(showerror, r)) && occursin("arm", lowercase(sprint(showerror, r))) && !ispath(c.out)      # the whole word (ind is not ind_pa1, nor `index`)
    end
    c = p3t_clone(p3t_full(), "d97h-old")
    m = AX.read_json(c.manifest); for a in ("ind_pa1", "d64_pa1"); delete!(m["arms"], a); end; AX._write_json(c.manifest, m)
    r = d97_call(() -> AN.analyze(; p3t_an_args(I, st, c)...))
    @test r isa ArgumentError && (occursin("ind_pa1", sprint(showerror, r)) || occursin("d64_pa1", sprint(showerror, r))) && !ispath(c.out)
end

@testset "D97h the texts of the readings P114 and P115 state E78, E76, the six arms and E75; the deviation of the dropped column replaces the deviation of tg at N_1 (P120)" begin
    R = d97_call(p3t_run_main)
    @test d97_ok(R)
    d97_ok(R) || return
    (; m, text, tab) = R
    p114 = first(r for r in m["readings"] if startswith(r, "P114"))
    p115 = first(r for r in m["readings"] if startswith(r, "P115"))
    @test all(w -> occursin(w, p114), ["ind_pa1", "L1", "specs", "chunk entries", "complete with failed chunks"]) && all(w -> occursin(w, p115), vcat(P3T_ARMS, ["N_c"]))
    @test occursin(p114, tab) && occursin(p115, tab)                                                                   # the same texts in the tables file
    @test !occursin("tg item of the ind arm is read at the horizon N_1", text) && !occursin("tg item of the ind arm is read at the horizon N_1", tab)
    dev = tab[findfirst("## Deviations", tab).start:findfirst("## Exploratory items", tab).start]
    @test occursin("labels agree", lowercase(dev)) && occursin("chk", dev) && occursin("labels agree", lowercase(text))      # the omitted column and its cause, in both files
end

# ---------------------------------------------------------------- D97i: the fixture is written in the shape of the runner (E82)
@testset "D97i the toy campaign is in the shape of the runner: no cell of the old form, six arms, the items at the horizons of the runner, the keys of the runner (fixture side; passes now)" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    inst = p3t_instance(S, I, "d97i")
    files = p3t_files(inst.store)
    keys_ = [AX.read_chunk_attrs(f) for f in files]
    @test !any(a -> occursin(r"cell=k\d+:r\d+$", a["chunk_key"]), keys_) && all(a -> occursin(r"cell=k\d+:[A-Z]:r\d+$", a["chunk_key"]), keys_)
    @test Set(only(a["arms"]) for a in keys_) == Set(P3T_ARMS) && all(a -> !occursin("  ", a["chunk_key"]) && occursin(" eps=- ", a["chunk_key"]) && occursin(" mode=wrapped ", a["chunk_key"]), keys_)
    @test all(a -> startswith(a["chunk_key"], "method=$(a["method"]) ") && (!(a["method"] in ("L1", "L5")) == occursin(" kind=", a["chunk_key"])) && (occursin(" kind=", a["chunk_key"]) == occursin(" spec=", a["chunk_key"])), keys_)
    chk = [a for a in keys_ if a["arms"] == ["chk"]]
    @test sort(unique(a["precision"] for a in chk)) == ["BigFloat(256)", "Double64"] && all(a -> a["horizon"] == P3T_NB && a["method"] == "L5", chk)
    @test sort(unique(a["precision"] for a in keys_ if a["arms"] == ["d64"] || a["arms"] == ["d64_pa1"])) == ["Double64"] && sort(unique(a["precision"] for a in keys_ if a["arms"] in (["ind"], ["ind_pa1"], ["ref"]))) == ["Float64"]
    man = AX.read_json(inst.manifest)
    @test all(a -> haskey(man["arms"][a], "specs"), ["ind", "ind_pa1", "d64", "d64_pa1"]) && !any(a -> haskey(man["arms"][a], "specs"), ["ref", "chk"]) &&
          Set(keys(man["arms"]["ind"]["specs"])) == Set(["tg", "fq1", "bk1"]) && Set(keys(man["arms"]["ind_pa1"]["specs"])) == Set(["pa1"]) && all(e -> e["steps"] isa Integer && e["thread_seconds"] isa Real && occursin("cell=", e["key"]), man["chunks"])
    # the cell and the arm of every chunk entry are those of its file (the key of the entry is the key of the chunk)
    byord = Dict(e["ordinal"] => e for e in man["chunks"])
    @test all(f -> (a = AX.read_chunk_attrs(f); byord[parse(Int, match(r"chunk-(\d{5})", f).captures[1])]["key"] == a["chunk_key"] && byord[parse(Int, match(r"chunk-(\d{5})", f).captures[1])]["arms"] == a["arms"]), files)
    # a fixture parser finds only cells of E74 in the panel too: every panel orbit is in one cell of each ref item
    pan = TOML.parsefile(inst.panel)
    @test all(o -> haskey(o, "class") && o["class"] in ("C", "R"), pan["orbit"]) && count(o -> o["class"] == "R", pan["orbit"]) == 3
end
