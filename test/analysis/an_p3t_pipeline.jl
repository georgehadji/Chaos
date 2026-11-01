# D95f to D95l, D95n, D95o: `analyze` for study p3-abstention end to end on the toy campaign (31 orbits, 3 rounds, closed log), against the independent reference. Spec: analysis-spec.md,
# Amendment D95, E63 to E73, I38 to I46, criteria D95f to D95l, D95n, D95o. Expected values: p3t_expected (an_p3t_helpers.jl): sets, strata, rule outputs, HA, HB, HC by the reference
# an_p3t_ref.jl (exact rationals, the generator of an_p3_ref.jl) from the designed verdict histories; the sensitivities and the power by the variants below.
# API proposed by the tests: analyze(; prereg, params, test_params, store, panel, round_log, manifest, p2_manifest, p2_prereg, p2_params, calibration_record, p1_prereg, p1_params, out, threads)
# -> (status = :ok, results, tables, f2, f3, decisions) with f2 and f3 (columns, rows, digest); the machine file p3-abstention-test-results.json has the keys read below.

ra(x) = x isa Integer ? P3T_R(x) : (m = match(r"^(-?\d+)(?:/(\d+))?$", String(x)); m === nothing ? error("not a rational: $x") : P3T_R(parse(BigInt, m.captures[1])) // (m.captures[2] === nothing ? 1 : parse(BigInt, m.captures[2])))
ra(x::Nothing) = nothing
const P3T_RUN = Dict{String, Any}()
function p3t_run_main()
    get!(P3T_RUN, "main") do
        S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
        inst = p3t_clone(p3t_full(), "main-run")
        res = AN.analyze(; p3t_an_args(I, st, inst)..., threads = 4)
        (; S, I, st, inst, res, m = AX.read_json(res.results), X = p3t_expected(S, I.qt), text = read(res.results, String), tab = read(res.tables, String))
    end
end

@testset "D95n D95o output: two files, the key order, no clock, no absolute path, the hash of the machine file, the lists" begin
    R = p3t_run_main(); (; I, inst, res, m, text, tab) = R
    @test res.status === :ok && basename(res.results) == "p3-abstention-test-results.json" && basename(res.tables) == "p3-abstention-tables.md" && an_files(inst.out) == ["p3-abstention-tables.md", "p3-abstention-test-results.json"]
    @test text == replace(text, "\r\n" => "\n") && endswith(text, "\n") && !occursin(inst.dir, text) && !occursin(replace(inst.dir, '\\' => '/'), text) && !occursin(r"\d{4}-\d{2}-\d{2}T\d{2}:", text)
    @test !occursin(inst.dir, tab) && !occursin(replace(inst.dir, '\\' => '/'), tab) && !occursin(r"\d{4}-\d{2}-\d{2}T\d{2}:", tab) && !occursin('\r', tab)
    top = ["program", "fingerprints", "parameters", "test_parameters", "log", "sets", "outputs", "decisions", "ha", "hb", "hc", "sensitivities", "power", "flags", "r_ab", "readings", "exploratory"]
    @test all(k -> haskey(m, k), top) && [match(r"^  \"([a-z_0-9]+)\":", l).captures[1] for l in split(text, '\n') if occursin(r"^  \"[a-z_0-9]+\":", l)] == top
    @test Set(keys(m["parameters"])) == Set(string.(P3_PARAM_KEYS)) && Set(keys(m["test_parameters"])) == Set(string.(P3T_KEYS)) && m["exploratory"] == [] && !isempty(m["readings"])
    # the decisions in the order HA, HB, HC; the fingerprints (SHA-256 of the inputs read and of every chunk file read, the commit and the clean flag)
    @test collect(keys(m["decisions"])) == ["HA", "HB", "HC"] || Set(keys(m["decisions"])) == Set(["HA", "HB", "HC"])
    fp = m["fingerprints"]
    @test fp["preregistration"] == I.prereg_hash && fp["p3_parameters"] == p3ref_sha256_file(I.params) && fp["test_parameters"] == p3ref_sha256_file(I.tparams) && fp["calibration_record"] == p3ref_sha256_file(I.record) &&
          fp["panel"] == p3ref_sha256_file(inst.panel) && fp["manifest"] == p3ref_sha256_file(inst.manifest) && fp["round_log"] == p3ref_sha256_file(inst.log) && fp["p2_manifest"] == I.q.p2_manifest_sha256 &&
          haskey(fp, "p1_parameters") && haskey(fp, "git_commit") && fp["git_clean"] isa Bool && length(fp["chunks"]) == 50 && all(v -> occursin(r"^[0-9a-f]{64}$", v), values(fp["chunks"]))
    # the tables file: header, the decisions, T1 to T11, the data blocks of F2 and F3 with digests, the readings, the deviations, the exploratory list; the machine file's hash is stated
    mh = match(r"results SHA-256:? ([0-9a-f]{64})", tab)
    @test mh !== nothing && mh.captures[1] == p3ref_sha256_file(res.results)
    @test all(t -> occursin("Table T$t", tab), 1:11) && occursin("## F2 data", tab) && occursin("## F3 data", tab) && occursin("digest: " * res.f2.digest, tab) && occursin("digest: " * res.f3.digest, tab)
    @test p2_digest(res.f2.columns, res.f2.rows) == res.f2.digest && p2_digest(res.f3.columns, res.f3.rows) == res.f3.digest
    @test occursin(r"(?i)readings", tab) && occursin(r"(?i)deviations", tab) && occursin(r"(?i)exploratory", tab) && occursin(r"(?i)exploratory[^\n]*\n+[^\n]*none", tab)
    for h in ("HA", "HB", "HC"); @test occursin(h, tab); end
end

@testset "D95f labels, drops and sets on the test orbits; the strata; the gate set; I43 and I46" begin
    (; X, m) = p3t_run_main()
    s = m["sets"]
    ids(k) = Int.(s["ids"][k])
    @test ids("E_test") == [r.id for r in X.E] && ids("W2") == [r.id for r in X.E if r.W2] && ids("T_s") == [r.id for r in X.E if r.Ts] && ids("S2") == [r.id for r in X.E if r.S2]
    @test ids("D") == [r.id for r in X.E if r.id in P3T_D] && ids("D_B") == P3T_DB && ids("gate") == [1, 2] && s["counts"]["overlap"] == 0
    @test [s["counts"][k] for k in ("E_test", "W2", "T_s", "S2", "D", "D_B", "gate")] == [25, 10, 5, 7, 15, 4, 2]
    # D is the E_test orbits of the panel flag in_D that carry Double64 verdict data: the drift orbit 14 (in_D, dropped) is not in D
    @test !(14 in ids("D")) && 14 in P3T_D && all(i -> i in ids("E_test"), ids("D"))
    @test Int.(s["strata"]) == X.ha.strata == [5, 5, 7, 1, 7] && sum(Int.(s["strata"])) == s["counts"]["E_test"]
    # drops by reason, round and K index (counted, listed)
    @test [(Int(d["id"]), d["reason"], Int(d["round"]), Int(d["k_index"])) for d in s["drops"]] == [(d.id, d.reason, d.round, d.kidx) for d in X.drops]
    # I43: the drops sum to the panel orbits not in E_test and the gate set; nesting; by-K and by-round rows sum to the pooled rows
    @test length(s["drops"]) + s["counts"]["gate"] + s["counts"]["E_test"] == 31 && all(i -> i in ids("W2"), ids("T_s"))
    rows = s["by_k_round"]
    @test all(k -> sum(Int(r[k]) for r in rows) == s["counts"][k], ("E_test", "W2", "T_s", "S2")) && all(r -> haskey(r, "k_index") && haskey(r, "round"), rows)
    # I46: the sizes of the final sets differ from the counts of the last log entry by exactly the orbits that lost verdict data (orbit 26, a W2 orbit): t_s 0, w2 1
    @test m["sets"]["log_difference"] == Dict("t_s" => 0, "w2" => 1) && Int(m["log"]["entries"][3]["counts"]["w2"]) - s["counts"]["W2"] == 1 && Int(m["log"]["entries"][3]["counts"]["t_s"]) - s["counts"]["T_s"] == 0
    @test [e["decision"] for e in m["log"]["entries"]] == ["continue", "continue", "close"]
end

@testset "D95g verdicts and rule outputs: counts of the three outputs per set; R_a and R_ab" begin
    (; X, m, S) = p3t_run_main()
    cnt(rs) = Dict(string(o) => count(r -> r.out === o, rs) for o in (:chaotic_supported, :regular_supported, :unresolved))
    o = m["outputs"]
    @test Dict(k => Int(v) for (k, v) in o["E_test"]) == cnt(X.E) && Dict(k => Int(v) for (k, v) in o["W2"]) == cnt([r for r in X.E if r.W2]) &&
          Dict(k => Int(v) for (k, v) in o["T_s"]) == cnt([r for r in X.E if r.Ts]) && Dict(k => Int(v) for (k, v) in o["S2"]) == cnt([r for r in X.E if r.S2])
    # R_ab on D: the nine pairs of outputs, the overrides of the ladder (R_a decided, R_ab unresolved), by the reference with the Double64 sequences; never silently R_a
    Dd = S.design.d; D = S.design.f
    ab = Dict{String, Int}(); over = 0
    for r in X.E
        r.id in P3T_D || continue
        a = p3ref_out(D.V[r.id], X.kP, X.th); b = p3ref_out_ab(D.V[r.id], Dd.V[r.id], X.kP, X.th)
        ab["$a|$b"] = get(ab, "$a|$b", 0) + 1
        a !== :unresolved && b === :unresolved && (over += 1)
    end
    @test Dict(k => Int(v) for (k, v) in m["r_ab"]["pairs"]) == ab && Int(m["r_ab"]["overrides"]) == over && sum(values(ab)) == 15 && Int.(m["r_ab"]["orbits"]) == [r.id for r in X.E if r.id in P3T_D]
    # the orbit 14 (drift) and the orbits without Double64 data are not in the R_ab table
    @test !(14 in Int.(m["r_ab"]["orbits"]))
end

@testset "D95h HA end to end: the point values, the bounds and the decision equal the reference" begin
    (; X, m) = p3t_run_main()
    h = m["ha"]; x = X.ha
    @test m["decisions"]["HA"]["verdict"] == x.decision == "not supported" && Int.(m["decisions"]["HA"]["failed"]) == x.failed == [2, 3]
    @test ra(h["e"]) == x.point.e && ra(h["dw"]) == x.point.dw && ra(h["dt"]) == x.point.dt && Int(h["g"]) == x.point.g == 0 && ra(h["a"]["E_test"]) == x.point.aE
    @test ra(h["a"]["W2"]) == x.point.a[:W2] && ra(h["a"]["T_s"]) == x.point.a[:Ts] && ra(h["a"]["S2"]) == x.point.a[:S2]
    for (k, b) in (("e_up", x.bounds.e_up), ("aE_up", x.bounds.aE_up), ("dw_lo", x.bounds.dw_lo), ("dt_lo", x.bounds.dt_lo))
        @test ra(h["bounds"][k]["value"]) == b.value && Int(h["bounds"][k]["dropped"]) == b.dropped
    end
end

# ---------------------------------------------------------------- the variants of HB: matched budgets of a cost table, round-down budgets, without LE
function p3t_hb_variant(X, qt; costs = nothing, down = false)
    q2 = costs === nothing ? qt : merge(qt, costs)
    M = p3ref_matched(q2, X.th.V)
    names = Symbol.(qt.voter_bits)
    bX(Y) = down ? M.X[Y].round_down : (Y === :LE ? M.le_read : M.X[Y].b)
    undecidable = any(Y -> bX(Y) === nothing, vcat(names, :LE))
    undecidable && return (decidable = false, verdict = "not decidable", failing = Tuple{Symbol, Symbol}[])
    D = X.cmp
    cmp = Dict{Tuple{Symbol, Symbol}, Any}()
    for ((A, Y), v) in D
        pos = findfirst(==(bX(Y)), qt.budgets)
        set = sort([r for r in X.E if A === :Ts ? r.Ts : r.W2]; by = r -> r.id)
        cmp[(A, Y)] = (rule = v.rule, x = [X.vfun(r.id, Y)[pos] === :regular for r in set])
    end
    hb = p3tref_hb(cmp, X.cov, X.sizes, qt)
    return (decidable = true, verdict = hb.decision, failing = hb.failing)
end

@testset "D95i HB end to end: the 13 comparisons, the coverage, the flags; the matched budgets" begin
    (; X, m, S, I) = p3t_run_main()
    h = m["hb"]; x = X.hb
    @test m["decisions"]["HB"]["verdict"] == x.decision == "not supported" && h["verdict"] == x.decision && length(h["comparisons"]) == 13
    rows = Dict((c["set"] == "T_s" ? :Ts : :W2, Symbol(c["indicator"])) => c for c in h["comparisons"])
    @test Set(keys(rows)) == Set(keys(x.rows))
    for (k, r) in x.rows
        c = rows[k]
        @test (Int(c["n11"]), Int(c["n10"]), Int(c["n01"]), Int(c["n00"])) == (r.n11, r.n10, r.n01, r.n00) && ra(c["p"]) == r.p && c["margin"] == r.margin && c["pass"] == r.pass && ra(c["fr_rule"]) == r.fr_rule && ra(c["fr_x"]) == r.fr_x
        @test Int(c["n11"]) + Int(c["n10"]) + Int(c["n01"]) + Int(c["n00"]) == (k[1] === :Ts ? 5 : 10) && 0 <= ra(c["p"]) <= 1
    end
    @test Set(Symbol.(h["failing"] .|> s -> replace(s, "T_s" => "Ts"))) == Set(Symbol("$(k[1]):$(k[2])") for k in x.failing) && ra(h["coverage"]["T_s"]) == X.cov[:Ts] && ra(h["coverage"]["W2"]) == X.cov[:W2] && h["coverage_ok"] == x.coverage_ok
    # matched budgets: b_X by the ceiling rule (SALI, FLI, MEGNO at 5120), LE read at N_c and flagged cost-unmatched; the sizes against the targets give the underpowered flag
    @test Int(rows[(:Ts, :SALI)]["b"]) == 5120 && Int(rows[(:W2, :LE)]["b"]) == 20480 && h["cost_unmatched_le"] == true == X.M.le_unmatched && h["underpowered"] == true == x.underpowered
    @test m["flags"]["underpowered_hb"] == true && m["flags"]["cost_unmatched_le"] == true && m["flags"]["underpowered_hc"] == false && isempty(h["defects"])
    # HB without LE: 12 comparisons, reported, no effect on the decision
    cmp12 = Dict(k => v for (k, v) in X.cmp if k != (:W2, :LE))
    w = p3tref_hb(cmp12, X.cov, X.sizes, I.qt)
    @test h["without_le"]["verdict"] == w.decision && length(h["without_le"]["comparisons"]) == 12 && m["decisions"]["HB"]["verdict"] == x.decision
end

@testset "D95j sensitivities (reported, no effect on any decision) and the power at the achieved sizes" begin
    (; X, m, S, I, st, inst) = p3t_run_main()
    qt = I.qt
    sv = m["sensitivities"]
    # (i) the round-down budgets, (ii) the pilot costs, (iii) the costs measured in the test run (the manifest has no value for L1: not decidable, P114), (iv) without LE
    rd = p3t_hb_variant(X, qt; down = true)
    @test sv["round_down"]["decidable"] == rd.decidable && sv["round_down"]["verdict"] == rd.verdict && rd.decidable
    pc = p3t_hb_variant(X, qt; costs = (; cost_l1 = qt.pilot_cost_l1, cost_tg = qt.pilot_cost_tg, cost_pa1 = qt.pilot_cost_pa1, cost_bk1 = qt.pilot_cost_bk1, cost_fq1 = qt.pilot_cost_fq1))
    @test sv["pilot"]["decidable"] == pc.decidable && sv["pilot"]["verdict"] == pc.verdict
    @test sv["measured"]["decidable"] == false && sv["measured"]["verdict"] == "not decidable" && occursin("L1", sv["measured"]["reason"])
    @test sv["without_le"]["verdict"] == p3tref_hb(Dict(k => v for (k, v) in X.cmp if k != (:W2, :LE)), X.cov, X.sizes, qt).decision
    @test Set(Symbol.(sv["round_down"]["failing"] .|> s -> replace(s, "T_s" => "Ts"))) == Set(Symbol("$(k[1]):$(k[2])") for k in rd.failing)
    @test m["decisions"]["HB"]["verdict"] == X.hb.decision                                                                                       # no sensitivity changes a decision
    # the power at the achieved sizes: |W2| = 10 is below target_w2 = 11, so the power of the full HB decision is recomputed from the rates of the record with the streams (pp-achieved-hb, run)
    pw = m["power"]["hb"]
    mod = p3t_model()
    nT, nRW = X.sizes.Ts, X.sizes.W2 - X.sizes.Ts
    @test pw !== nothing && Int(pw["runs"]) == qt.pp_runs && Int(pw["successes"]) == p3tref_power_hb(qt, mod, nT, nRW, qt.pp_runs) && pw["power"] isa Real && abs(pw["power"] - Int(pw["successes"]) / qt.pp_runs) < 1e-12
    @test pw["planned"] == qt.power_hb_planned && 0 <= pw["se"] <= 1
    @test m["power"]["hc"] === nothing                                                                  # failures 7 >= target_f 6: no recomputation (and no type-I check is repeated)
end

@testset "D95j no recomputation when both targets are met; the sensitivity (iii) with the costs of the manifest" begin
    S, I0, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    I = p3t_inputs(S, st, "met"; target_t_s = 5, target_w2 = 10)                 # the final sizes (5 and 10) meet the targets
    inst = p3t_instance(S, I, "met"); write(inst.log, join(first(p3t_log_entries(S, I))))
    m = AX.read_json(AN.analyze(; p3t_an_args(I, st, inst)..., threads = 2).results)
    @test m["power"]["hb"] === nothing && m["flags"]["underpowered_hb"] == false && m["hb"]["underpowered"] == false
    # the measured costs: the manifest gets a value for every spec; tg is measured at 237.5 against 2658, so the matching of SALI, FLI, MEGNO is capped: not decidable (a non-LE comparator capped)
    man = AX.read_json(inst.manifest)
    man["arms"]["ref"]["specs"]["L1"] = Dict{String, Any}("budget_line_ns" => 300.0, "measured_ns_per_orbit_step" => 207.0)
    AX._write_json(inst.manifest, man)
    m2 = AX.read_json(AN.analyze(; p3t_an_args(I, st, inst)..., out = joinpath(inst.dir, "out2")).results)
    X = p3t_expected(S, I.qt)
    mv = p3t_hb_variant(X, I.qt; costs = (; cost_l1 = 207.0, cost_tg = 237.5, cost_pa1 = 1810.5, cost_bk1 = 3000.25, cost_fq1 = 16850.5))
    @test mv.decidable == false && m2["sensitivities"]["measured"]["decidable"] == false && m2["sensitivities"]["measured"]["verdict"] == "not decidable" && !occursin("L1", m2["sensitivities"]["measured"]["reason"])
end

@testset "D95k HC end to end and the power at the achieved failures" begin
    (; X, m) = p3t_run_main()
    h = m["hc"]; x = X.hc
    @test m["decisions"]["HC"]["verdict"] == x.decision == "supported" && h["verdict"] == x.decision
    @test ra(h["auc_i"]) == x.auc_i == 55 // 56 && ra(h["auc_p"]) == x.auc_p == 41 // 56 && ra(h["delta"]) == x.delta == 1 // 4 && ra(h["lower"]) == x.lower && ra(h["upper"]) == x.upper
    @test Int(h["failures"]) == x.failures == 7 && Int(h["nonfailures"]) == x.nonfailures == 8 && Int(h["dropped"]) == x.dropped && h["underpowered"] == false && 0 <= ra(h["auc_i"]) <= 1 && 0 <= ra(h["auc_p"]) <= 1
    # the failures by hand: tied votes (3 4 9 15 16 25) and the wrong plain label (12)
    @test Int.(h["failure_ids"]) == [3, 4, 9, 12, 15, 16, 25]
    # underpowered HC: target_f = 8 > 7 failures: the flag, the planned power next to the recomputed power at the achieved failures, the decision unchanged
    S, st = p3t_cached(), p3_cached("d95")
    I = p3t_inputs(S, st, "hc-under"; target_f = 8)
    inst = p3t_instance(S, I, "hc-under"); write(inst.log, join(first(p3t_log_entries(S, I))))
    m2 = AX.read_json(AN.analyze(; p3t_an_args(I, st, inst)...).results)
    @test m2["hc"]["underpowered"] == true && m2["flags"]["underpowered_hc"] == true && m2["decisions"]["HC"]["verdict"] == "supported" && m2["power"]["hc"]["planned"] == I.t.power_hc_planned &&
          Int(m2["power"]["hc"]["runs"]) == I.qt.pp_runs && 0 <= m2["power"]["hc"]["power"] <= 1 && m2["power"]["hc"]["se"] >= 0
end

@testset "D95l the reported items T1 to T11: sums, the record, the figure data F2 and F3" begin
    (; X, m, res, tab, I, S) = p3t_run_main()
    rec = AX.read_json(I.record)
    sec(n) = (a = findfirst("Table T$n", tab); b = findfirst("Table T$(n + 1)", tab); a === nothing ? "" : tab[a.start:(b === nothing ? lastindex(tab) : b.start - 1)])
    t1 = sec(1)
    @test occursin("36", t1) && occursin("SALI", t1) && occursin("MEGNO", t1) && occursin(string(I.t.target_t_s), t1) && occursin(I.t.calibration_record_sha256, tab) && occursin(I.t.p3_table_sha256, tab)
    @test occursin("round", lowercase(sec(2))) && occursin("25", sec(2)) && occursin("31", sec(2))
    @test occursin("N_B", sec(3)) || occursin("D_B", sec(3))
    t3rows = [l for l in split(sec(3), '\n') if occursin(r"^\|\s*(3|6|8|12)\s*\|", l)]
    @test length(t3rows) == 4                                                                                  # one row per orbit of D_B
    @test occursin("drift", sec(11)) && occursin("failed", sec(11)) && occursin("diagnostics exceeded", sec(11))
    @test occursin("210.3", sec(10)) && occursin("237.5", sec(10)) && occursin("2658", sec(10))              # T10: the pilot keys, the measured values of the manifest and the cost keys
    # T9: the nine pairs of outputs on D sum to |D|
    @test sum(Int(v) for v in values(m["r_ab"]["pairs"])) == 15
    # F2: the rule family around theta* (q_C and q_R moved by one step inside their ranges: 4 points) and the seven single indicators at their matched budgets
    f2 = res.f2
    @test f2.columns == ["curve", "point", "decided_share", "error"] && length(f2.rows) == 4 + 7 && count(r -> r[1] == "rule", f2.rows) == 4 && count(r -> r[1] == "single", f2.rows) == 7
    nE = length(X.E)
    function famrow(qC, qR)
        th = merge(X.th, (; qC, qR))
        outs = [p3ref_out(S.design.f.V[r.id], X.kP, th) for r in X.E]
        dec = count(!=(:unresolved), outs); contra = count(i -> (outs[i] === :chaotic_supported && X.E[i].rho === :R) || (outs[i] === :regular_supported && X.E[i].rho === :C), 1:nE)
        return (dec // nE, dec == 0 ? nothing : contra // dec)
    end
    for (qC, qR) in ((2, 3), (1, 3), (3, 3), (2, 4))
        r = first(r for r in f2.rows if r[1] == "rule" && occursin("q_C=$qC", r[2]) && occursin("q_R=$qR", r[2]))
        want = famrow(qC, qR)
        @test isapprox(r[3], Float64(want[1]); rtol = 1e-12) && (want[2] === nothing ? !isfinite(r[4]) || r[4] === nothing : isapprox(r[4], Float64(want[2]); rtol = 1e-12, atol = 1e-15))
    end
    for Y in ("LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA")
        r = first(r for r in f2.rows if r[1] == "single" && r[2] == Y)
        b = Y == "LE" ? X.M.le_read : X.M.X[Symbol(Y)].b
        pos = findfirst(==(b), I.qt.budgets)
        vs = [S.design.f.V[e.id][Symbol(Y)][pos] for e in X.E]
        dec = count(!=(:undecided), vs); contra = count(i -> (vs[i] === :chaotic && X.E[i].rho === :R) || (vs[i] === :regular && X.E[i].rho === :C), 1:nE)
        @test isapprox(r[3], Float64(dec // nE); rtol = 1e-12) && (dec == 0 || isapprox(r[4], Float64(contra // dec); rtol = 1e-12, atol = 1e-15))
    end
    # F3: the ROC points of s_I and s_P for f on the HC set, for each distinct score in decreasing order (fpr, tpr), with the AUC rows
    f3 = res.f3
    R = p3tref_hc_rows(X.hcrows)
    f = [r.f for r in R]
    for (curve, sc) in (("s_I", [r.sI for r in R]), ("s_P", [r.sP for r in R]))
        want = p3tref_roc(sc, f)
        got = [r for r in f3.rows if r[1] == curve && r[2] != "AUC"]
        @test length(got) == length(want) && all(i -> isapprox(got[i][3], Float64(want[i].fpr); rtol = 1e-12, atol = 1e-15) && isapprox(got[i][4], Float64(want[i].tpr); rtol = 1e-12, atol = 1e-15), eachindex(want))
        a = first(r for r in f3.rows if r[1] == curve && r[2] == "AUC")
        @test isapprox(a[3], Float64(p3tref_auc(sc, f)); rtol = 1e-12)
    end
end

@testset "D95o determinism and read-only: runs at 1, 2 and 4 threads and in other folders give byte-identical files; the inputs, the store and data/ are unchanged" begin
    (; I, st, inst, res) = p3t_run_main()
    c = p3t_clone(p3t_full(), "determinism")
    h0 = p2_datatree_hash(); fp = p3t_fp_all(I, st, c)
    r1 = AN.analyze(; p3t_an_args(I, st, c)..., threads = 1, out = joinpath(c.dir, "o1"))
    @test p3t_fp_all(I, st, c) == fp && p2_datatree_hash() == h0 && an_files(joinpath(c.dir, "o1")) == ["p3-abstention-tables.md", "p3-abstention-test-results.json"]
    io = IOBuffer()
    @test AN.cli(p3t_cli_an(I, st, merge(c, (; out = joinpath(c.dir, "o2"))); extra = ["--threads", "2"]); io = io) == 0
    for f in ("p3-abstention-tables.md", "p3-abstention-test-results.json")
        a = read(joinpath(inst.out, f)); @test a == read(joinpath(c.dir, "o1", f)) == read(joinpath(c.dir, "o2", f))
    end
    @test p3t_fp_all(I, st, c) == fp && p2_datatree_hash() == h0 && r1.status === :ok
    # the command line prints the decisions and exits with 0
    out = String(take!(io))
    @test occursin("HA", out) && occursin("HB", out) && occursin("HC", out)
end
