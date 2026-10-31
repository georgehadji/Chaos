# Amendment D93, stage 02: fixture self-test. Not a test of the module: it does not call Chaoslab.Analysis (except the existing P1 `analyze`, a base function) and it passes
# before the module has any D93 functionality. It shows that a RED result elsewhere in the D93 suite is not caused by a broken fixture, helper or reference implementation
# (stage 02 audit, RED).

@testset "D93 fixture: the parameter table, the preregistration and the mirrors" begin
    @test length(P3_PARAM_KEYS) == 71 && allunique(P3_PARAM_KEYS) && Set(keys(p3_params())) == Set(P3_PARAM_KEYS) && length(P3_CARRIED) == 8
    st = p3_cached("default")
    q = st.q
    tab = TOML.parsefile(st.params)
    @test Set(Symbol.(keys(tab))) == Set(P3_PARAM_KEYS)
    text = read(st.prereg, String)
    @test AX.check_prereg(st.prereg).header_sha256 == st.prereg_hash && occursin("| Study slug | p3-abstention |", text)
    function section(name)
        m = match(Regex("(?m)^#+ " * replace(name, r"([()|.\[\]{}^\$*+?\\])" => s"\\\1") * "\\n(.*?)(?=^## |\\z)", "s"), text)
        m === nothing ? nothing : m.captures[1]
    end
    bad = String[]
    for k in P3_PARAM_KEYS
        e = tab[string(k)]
        sec = section(e["section"])
        sec === nothing && (push!(bad, "$k: section"); continue)
        occursin(e["quote"], sec) || push!(bad, "$k: quote not in section")
        e["text"] isa String && occursin(e["text"], e["quote"]) || push!(bad, "$k: text not in quote")
    end
    @test isempty(bad)
    @test tab["a_cal"]["text"] == "0.30" && tab["p_hb"]["text"] == "0.80" && tab["bootstrap_resamples"]["text"] == "10^4" && tab["delta_tol"]["text"] == "0.05"
    # the mirrors: the p2 table of the study holds the same numbers
    p2 = st.s2.p
    @test q.budgets == p2.budgets && q.b_primary == p2.b_primary && q.N_c == p2.N_c && q.N_1 == p2.N_1 && q.n_s == p2.n_s && q.f_s == p2.f_s && q.min_set == p2.min_set &&
          q.gate_size == p2.gate_size && q.lam_w2 == p2.w2_lambda_max && q.lam_s2 == p2.s2_lambda_min
    @test q.p2_manifest_sha256 == an_sha256(read(st.s2.manifest)) && q.voter_bits == q.indicators_i7[2:end] && q.budgets[end] == q.N_c && q.b_primary in q.budgets
    @test AX.read_json(st.s2.manifest)["study"] == P2_SLUG && AX.read_json(st.p1.manifest)["study"] == P1_SLUG && AX.read_json(st.p1.manifest)["panel"]["sha256_end"] == st.s2.panel_hash
    @test p3_texts(q)[:spec_of] == "LE=l1, SALI=tg, FLI=tg, MEGNO=tg, E1=fq1, E2=pa1, WBA=bk1"
end

@testset "D93 fixture: the store realises the designed verdicts; the P1 reference holds the L5 items of C0" begin
    st = p3_cached("default")
    s2 = st.s2
    ref = p2_ref_all(s2)
    d = st.design
    bad = [(id, X) for ((id, X), v) in ref.V1 if v != d.V[id][X]]
    @test isempty(bad) && length(ref.V1) == sum(length(d.V[id]) for id in 1:44) && s2.chunks == 33 && s2.items == 466
    # the P1 store has, for every panel orbit, the L5 and L1 items at N_1 = 80000 (lifted, pi/4) with the steps N_1 / 10, N_1 / 2 and N_1; the grid of the p2 items has them too
    p1 = st.p1
    items = [(ch, o) for ch in (AX.read_chunk(joinpath(p1.store, f)) for f in sort(readdir(p1.store))) for o in ch.orbits]
    l5 = [(ch, o) for (ch, o) in items if ch.attrs["method"] == "L5" && ch.attrs["horizon"] == P2_N1 && ch.attrs["coordinate_mode"] == "lifted" && ch.attrs["direction"] == "pi/4" && ch.attrs["tau"] == 1 &&
          ch.attrs["precision"] == "Float64"]
    @test sort([o.id for (_, o) in l5]) == collect(1:20) && all(((ch, o),) -> all(n -> n in o.t, (P2_N1 ÷ 10, P2_N1 ÷ 2, P2_N1)), l5)
    @test [o.attrs["flags"] for (_, o) in l5 if o.id == 13] == [["diagnostics exceeded"]] && all(((_, o),) -> isempty(o.attrs["flags"]), [x for x in l5 if x[2].id != 13])
    l1 = [(ch, o) for (ch, o) in items if ch.attrs["method"] == "L1" && ch.attrs["horizon"] == P2_N1 && ch.attrs["coordinate_mode"] == "lifted" && ch.attrs["direction"] == "pi/4" && ch.attrs["tau"] == 1 &&
          ch.attrs["precision"] == "Float64" && "main" in ch.attrs["arms"]]
    lam(ts, id) = (x = first(o for (_, o) in ts if o.id == id); x.lyap[first(keys(x.lyap))][end, 1])
    @test lam(l1, 14) - lam(l5, 14) > 9.9e-7 && abs(lam(l1, 3) - lam(l5, 3)) == 0.0 && lam(l5, 14) == 0.30 + 0.01 * 14
    # the reduced rule on these items (the independent reference of E49 c) gives the designed rho_new of the panel; the p2 L5 items give the designed rho_new of the grid
    qd = (label_lambda_min = 1e-2, label_growth_ratio = 5.0, label_rotation_tol = 1e-7, label_estimator_bound = 1e-8)
    function panel_rho(id)
        c5, o5 = first(x for x in l5 if x[2].id == id)
        o1 = first(o for (_, o) in l1 if o.id == id)
        i(n) = searchsortedfirst(o5.t, n)
        x(n) = o5.state[i(n), 1]
        x0 = parse(Float64, split(o5.attrs["ic"], ',')[1])
        lam5 = o5.lyap["L5"][i(P2_N1), 1]
        rho1, rho2 = (x(P2_N1 ÷ 2) - x0) / (2pi * (P2_N1 ÷ 2)), (x(P2_N1) - x(P2_N1 ÷ 2)) / (2pi * (P2_N1 - P2_N1 ÷ 2))
        return p3ref_label((lam5 = lam5, lam5_tenth = o5.lyap["L5"][i(P2_N1 ÷ 10), 1], lam1 = o1.lyap["L1"][i(P2_N1), 1], rho1 = rho1, rho2 = rho2,
                            diag5 = "diagnostics exceeded" in o5.attrs["flags"], diag1 = "diagnostics exceeded" in o1.attrs["flags"]), qd).rho
    end
    @test all(id -> panel_rho(id) === P3_RHO_NEW[id], 1:16)
    # the existing base pipeline on the P1 reference: the old labels of the panel orbits (what the change table of E49 f compares with)
    res = AN.analyze(; prereg = p1.prereg, manifest = p1.manifest, panel = p1.panel, params = p1.params, store = p1.store, out = joinpath(st.dir, "p1-out"))
    old(id) = (l = res.labels[id]; l.instability == "positive growth supported" && l.structure == "irregular" ? :C : l.instability == "none detected on the tested horizon" && l.structure == "quasiperiodic" ? :R : :U)
    @test old(2) === :U && P3_RHO_NEW[2] === :C && old(13) === :U && P3_RHO_NEW[13] === :U && old(1) === :C && old(9) === :R
end

@testset "D93 fixture: the reference implementation against hand cases and the designed selection" begin
    q = p3_params()
    # E52: 2088 = 720 + 720 + 540 + 108, derived from the keys; the enumeration has exactly that size, every element once, in the order of the amendment
    th = p3ref_theta(q)
    @test length(th) == 2088 == p3ref_count(q) && [count(t -> length(t.V) == k, th) for k in 3:6] == [720, 720, 540, 108] && [t.index for t in th] == 1:2088
    @test th[1] == (index = 1, V = [:SALI, :FLI, :MEGNO], mask = 7, qC = 1, qR = 2, w = 0, v = 0) && th[2].v == 1 && th[3].w == 1 && th[7].qR == 3 && th[end].V == Symbol.(q.voter_bits) && th[end].qC == 6
    @test allunique((t.V, t.qC, t.qR, t.w, t.v) for t in th)
    # E53 on the real cost table: the cost cap removes every V that contains E1 (690 elements pass), LE is always capped
    cap = p3ref_cap(q)
    @test cap == 18105 // 10 * 20480 && count(t -> p3ref_CR(q, t.V, q.b_primary) <= cap, th) == 690
    @test all(t -> (:E1 in t.V) == (p3ref_CR(q, t.V, q.b_primary) > cap), th)
    m = p3ref_matched(q, [:SALI, :FLI, :MEGNO])
    @test m.CR == 2560 * 26580 // 10 && m.X[:SALI].b == 2560 && m.X[:E1].b == 640 && m.X[:E2].b == 5120 && m.X[:WBA].b == 2560 && m.X[:LE].capped && m.le_read == q.N_c && m.le_unmatched
    @test m.X[:E1].round_down == 320 && m.X[:E2].round_down == 2560 && m.X[:SALI].round_down == 2560
    # exact decimals: 2076 // 10, and a tie at the ceiling rule is exact (b x cost == C_R)
    @test p3ref_rat(207.6) == 2076 // 10 && p3ref_rat(1e-4) == 1 // 10000 && p3ref_rat(1.0e10) == 10^10 && p3ref_rat(0.1) == 1 // 10
    # the rule: chaotic precedence, the stability of w = 1, 2 and the voters that are unstable
    V = Dict(:A => p2_sym("rrrrrrcccc"), :B => p2_sym("rrrrrrrrcc"), :C => p2_sym("uuuuuurrrr"))
    th1 = (V = [:A, :B, :C], qC = 1, qR = 2, w = 0, v = 1)
    @test p3ref_out(V, 6, th1) === :chaotic_supported && p3ref_out(V, 6, merge(th1, (qC = 2,))) === :regular_supported && p3ref_out(V, 6, merge(th1, (qC = 2, v = 0))) === :unresolved
    @test p3ref_evidence(V, 6, merge(th1, (w = 1,))) == (nC = 0, nR = 1, nCany = 1) && p3ref_evidence(V, 7, merge(th1, (w = 1,))) == (nC = 1, nR = 2, nCany = 1)
    @test p3ref_evidence(V, 8, merge(th1, (w = 2,))) == (nC = 1, nR = 1, nCany = 2) && p3ref_stable(V[:A], 8, 2) && !p3ref_stable(V[:A], 7, 2)
    @test_throws ErrorException p3ref_stable(V[:A], 1, 2)
    # exact tail of the McNemar test, rounding, ceil5, Hanley and McNeil (O17)
    @test p3ref_mcnemar(0, 1) == 1 // 2 && p3ref_mcnemar(1, 1) == 3 // 4 && p3ref_mcnemar(0, 5) == 1 // 32 && p3ref_mcnemar(2, 3) == 1 // 2 && p3ref_mcnemar(0, 0) == 1
    @test [p3ref_round_away(x) for x in (7125 // 100, 30875 // 100, 2375 // 10, -5 // 2)] == [71, 309, 238, -3] && [p3ref_ceil_step(x // 1, 5) for x in (71, 75, 76)] == [75, 75, 80]
    @test p3ref_hm(1 // 2, 2, 2) == 5 // 48 && p3ref_hm(1, 2, 2) == 0 && p3ref_hm(4 // 5, 10, 20) == 197 // 22500 && p3ref_hm(3 // 5, 5, 100) == 6597 // 350000
    @test p3ref_auc([5, 6], [1, 2]) == 1 && p3ref_auc([1, 2, 3], [1, 2, 3]) == 1 // 2 && p3ref_auc([2, 3], [2, 1]) == 7 // 8 && p3ref_auc(Int[], [1]) === nothing
    # the draw path of the real keys: 17 points, the allowed ones, the sizes quoted in the amendment (1.75 gives 100, 333, 433)
    path = p3ref_path(q)
    @test length(path) == 17 && path[1].rho == 1 && path[end].rho == 5 && path[4].rho == 7 // 4 && (path[4].nT, path[4].nRW, path[4].nW) == (100, 333, 433)
    @test [x.j for x in path if x.allowed] == 3:16 && path[end].nT == 285 && path[end].nW == 1235 && !path[3].allowed
    @test p3ref_sha256("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad" && p3ref_sha256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" &&
          all(t -> p3ref_sha256(t) == an_sha256(t), ("chaoslab-analysis-1|20261030|p3|optimism|1", "x"^55, "y"^56, "z"^64, "w"^119))
    # generators: SplitMix64 and xoshiro256++ against the published first outputs; the stream text rule against an independent SHA-256
    @test p3ref_splitmix(UInt64(0))[1] == 0xe220a8397b1dcdaf && p3ref_splitmix(UInt64(0))[2] == 0x6e789e6aa1b965f4
    g = P3RefGen((UInt64(1), UInt64(2), UInt64(3), UInt64(4)))
    @test [p3ref_next!(g) for _ in 1:2] == UInt64[41943041, 58720359]
    @test p3ref_seed(20261030, "optimism", 1) == UInt64(5692226309264838358) && p3ref_seed(20261030, "pp-hb", (0, 1)) == UInt64(3437253884877108471) && p3ref_seed(20261030, "pp-hc", (20, 1)) == UInt64(15836479185517542896)
    # the designed selection on the 36 orbits of Cal: the sets by design, then a selection exists (the PP chain is run on it)
    st = p3_cached("default")
    rows = p3_ref_rows(st)
    cal = [o for o in rows if o.in_cal]
    @test length(cal) == 36 && count(o -> o.rho === :R, cal) == 12 && count(o -> o.S2, cal) == 8 && count(o -> o.Ts, cal) == 6 && count(o -> o.W2, cal) == 13 && count(o -> o.gate, rows) == 2
    @test sort([o.id for o in cal if o.S2]) == [11, 12, 33, 34, 35, 36, 37, 38] && sort([o.id for o in cal if o.Ts]) == [5, 6, 23, 24, 25, 39] && sort([o.id for o in cal if o.W2]) == [1, 3, 5, 6, 7, 8, 21, 22, 23, 24, 25, 26, 39]
    @test [(r.id, r.V) for r in rows] == [(r.id, r.V) for r in p3_design_rows(st.design)]
    sel = p3ref_select(rows, st.q)
    @test sel.outcome === :selected && sel.alone == (690, 1745, 1044, 1840, 681, 554) && length(sel.feasible) == 62 && length(sel.tol) == 42 && sel.J_min == 19 // 156
    @test all(i -> sel.first_fail[i] == 0, sel.feasible) && all(i -> sel.J[i] <= sel.J_min + p3ref_rat(st.q.delta_tol), sel.tol) && sel.star in sel.tol
    @test sel.star == 147 && sel.thetas[147] == (index = 147, V = [:SALI, :FLI, :E2], mask = 19, qC = 1, qR = 2, w = 1, v = 0) && sel.stats[147].e == 1 // 32 && sel.stats[147].decided == 32
    # hand anchors of the design: three strong orbits (33, 34, 35) have undecided tangent-map indicators, so every V of tg alone abstains on 3 of the 8 strong orbits (F3);
    # the gate orbit 10 is called chaotic at b_P by FLI and E2 only, unstable: a rule that trusts both at w = 0 fires on the gate (F4), at w = 1 it does not
    tg = [i for (i, t) in enumerate(sel.thetas) if all(X -> X in (:SALI, :FLI, :MEGNO), t.V)]
    @test length(tg) == 3 * 2 * 3 * 2 && all(i -> sel.stats[i].u_s2 == 3 // 8 && sel.first_fail[i] in (2, 3), tg)
    fe(qC, w) = findfirst(t -> t.V == [:SALI, :FLI, :E2] && t.qC == qC && t.qR == 2 && t.w == w && t.v == 0, sel.thetas)
    @test sel.stats[fe(2, 0)].g == 1 && sel.stats[fe(2, 1)].g == 0 && !sel.flags[fe(2, 0)][4] && sel.flags[fe(2, 1)][4]
    println("D93 fixture selection: |F| = $(length(sel.feasible)), |Theta_tol| = $(length(sel.tol)), J_min = $(sel.J_min), theta* = $(sel.thetas[sel.star]), alone = $(sel.alone)")
end

@testset "D93 fixture: the generator vectors written from an independent implementation (Python standard library) and the exact power oracle" begin
    # computed with Python 3 (hashlib, statistics.NormalDist, an own SplitMix64 and xoshiro256++), see the test plan, "Stream vectors"
    @test p3ref_splitmix(p3ref_seed(20261030, "optimism", 1)) == (UInt64(15300154784399154531), UInt64(1595416676562361740), UInt64(9699484540660912394), UInt64(17671760824607871275))
    g = p3ref_stream(20261030, "optimism", 1)
    @test [p3ref_next!(g) for _ in 1:4] == UInt64[1960354936839448390, 11115532021063746583, 17181476028355880092, 5824767452774032242]
    g = p3ref_stream(20261030, "optimism", 1)
    @test [p3ref_uniform!(g) for _ in 1:4] == [0.10627105406820064, 0.6025741982784751, 0.931409681822553, 0.3157612763260232]
    # exact power of a single comparison, by the sum over all cells, against the reference simulation (4 standard errors; the same draws as the module must give)
    pp = merge(p3_params(), (; alpha_hb = 0.05, rho_hb = 0.5, c_min = 0.5))
    ex = Float64(p3ref_exact_power(16, 1 // 10, 1 // 5, 1 // 2, 5 // 100, 1 // 2))
    model = p3_model(rule_T = (p = 0.1, d = 1.0), comp = Dict(:SALI => (T = (p = 0.2, s = 0.5), RW = (p = 0.0, s = 0.0))), ts = [:SALI])
    runs = 4000
    succ = p3ref_hb_power(20261030, 0, model, 16, 0, pp, runs)
    se = sqrt(ex * (1 - ex) / runs)
    @test 0 < ex < 1 && abs(succ / runs - ex) <= 4 * se
end
