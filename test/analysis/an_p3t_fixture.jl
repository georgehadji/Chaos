# Amendment D95, stage 02: fixture self-test. Not a test of the module: it calls no D95 function (the existing D93 JSON writer `AN.p3_json` is used to write the record, and the
# existing Experiments readers read the chunk files) and it passes before the module has any D95 functionality. It shows that a RED result elsewhere in the D95 suite is not
# caused by a broken fixture, helper or reference implementation (stage 02 audit, RED). The hand tables below are written from the design (an_p3t_helpers.jl), not from code.

@testset "D95 fixture: the reference helpers on hand cases" begin
    # entry text: 13 lines, the count lines in ascending byte order, the empty-text hash of the first entry
    e = p3tref_entry(; round = 1, decision = "continue", counts = Dict("t_s" => 2, "drift" => 1, "excluded" => 0, "orbits" => 12, "target_t_s" => 5, "target_w2" => 11, "w2" => 4),
                     panel = "a"^64, prereg = "b"^64, digest = "c"^64, prev = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")
    L = split(e, '\n')
    @test length(L) == 14 && L[end] == "" && count(==('\n'), e) == 13
    @test L[1:9] == ["round=1", "decision=continue", "count.drift=1", "count.excluded=0", "count.orbits=12", "count.t_s=2", "count.target_t_s=5", "count.target_w2=11", "count.w2=4"]
    @test L[10:13] == ["panel_sha256=" * "a"^64, "prereg_sha256=" * "b"^64, "ref_digest=" * "c"^64, "prev_sha256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"]
    @test p3ref_sha256("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    # the count rule at its boundaries (O21 g)
    @test p3tref_decide(5, 11, 5, 11, 1, 3) == "close" && p3tref_decide(4, 11, 5, 11, 1, 3) == "continue" && p3tref_decide(5, 10, 5, 11, 2, 3) == "continue" &&
          p3tref_decide(0, 0, 5, 11, 3, 3) == "close" && p3tref_decide(9, 9, 5, 11, 2, 3) == "continue"
    # type 7 percentile by hand: B' = 20, x = 1..20: p = 19/20 gives h = 19.05, 19 + 0.05; p = 1/20 gives h = 1.95, 1 + 0.95; an integer h
    @test p3tref_percentile(1:20, 19 // 20) == 19 + 1 // 20 && p3tref_percentile(1:20, 1 // 20) == 1 + 19 // 20 && p3tref_percentile([1, 2, 3, 4, 5], 1 // 2) == 3 && p3tref_percentile([7], 1 // 2) == 7
    @test p3tref_percentile(Rational{Int}[], 1 // 2) === nothing
    # McNemar: m = 0 gives 1, (n10 = 0, n01 = 5) gives 1/32, (2, 3) gives (C(5,3)+C(5,4)+C(5,5))/32 = 16/32
    @test p3ref_mcnemar(0, 0) == 1 && p3ref_mcnemar(0, 5) == 1 // 32 && p3ref_mcnemar(2, 3) == 1 // 2 && p3ref_mcnemar(5, 0) == 1
    # AUC: separated classes 1, equal scores 1/2, one tied pair one half; ROC of a four-orbit case: scores 3, 2, 2, 1 with failures at 3 and 2
    @test p3tref_auc([3, 3, 1, 1], [true, true, false, false]) == 1 && p3tref_auc([1, 1, 1, 1], [true, true, false, false]) == 1 // 2 && p3tref_auc([2, 1], [true, false]) == 1 &&
          p3tref_auc([2, 2, 1], [true, false, false]) == 3 // 4 && p3tref_auc([1, 2], [false, false]) === nothing
    roc = p3tref_roc([3, 2, 2, 1], [true, true, false, false])
    @test [(r.threshold, r.fpr, r.tpr) for r in roc] == [(3, 0, 1 // 2), (2, 1 // 2, 1), (1, 1, 1)]
    # the stratified resample of two strata of 2 and 3 orbits: the integers come from the independent Python generator (streams.py of the scratch area, SplitMix64 + xoshiro256++)
    @test p3tref_resample(20261030, "boot-ha", 1, [2, 3]) == [[2, 1], [3, 1, 3]] && p3tref_resample(20261030, "boot-ha", 2, [2, 3]) == [[1, 2], [1, 3, 2]]
    # the first integers of the four streams of E66 a (1 to 25) written by the independent generator
    ints(p, m) = (g = p3ref_stream(20261030, p, 1); [p3ref_randint!(g, m) for _ in 1:8])
    @test ints("boot-ha", 25) == [4, 8, 25, 23, 10, 17, 23, 21] && ints("boot-hc", 25) == [7, 6, 12, 4, 21, 23, 16, 15] &&
          ints("pp-achieved-hb", 25) == [18, 11, 19, 7, 5, 25, 15, 24] && ints("pp-achieved-hc", 25) == [21, 24, 23, 24, 4, 4, 15, 5]
    @test p3ref_seed(20261030, "boot-ha", 1) == 6130294642297488880 && p3ref_seed(20261030, "boot-hc", 1) == 15727976254305454371 &&
          p3ref_seed(20261030, "pp-achieved-hb", 1) == 4500577164153853627 && p3ref_seed(20261030, "pp-achieved-hc", 1) == 4279710978371254201
end

@testset "D95 fixture: the preregistration, the tables, the record and the store" begin
    st = p3_cached("d95")
    S = p3t_cached()
    I = p3t_cached_inputs()
    @test length(P3T_KEYS) == 18 && allunique(P3T_KEYS) && length(P3T_PILOT) == 5 && length(S.orbits) == 31 && S.files == 50 && length(readdir(S.store)) == 50
    text = read(I.prereg, String)
    @test AX.check_prereg(I.prereg).header_sha256 == I.prereg_hash && occursin("| Study slug | p3-abstention |", text) && occursin("## $(P3T_A1)", text)
    sec(name) = (m = match(Regex("(?m)^#+ " * replace(name, r"([()|.\[\]{}^\$*+?\\])" => s"\\\1") * "\\n(.*?)(?=^## |\\z)", "s"), text); m === nothing ? nothing : m.captures[1])
    tab = TOML.parsefile(I.tparams)
    @test Set(Symbol.(keys(tab))) == Set(P3T_KEYS)
    bad = String[]
    for k in P3T_KEYS
        e = tab[string(k)]; s = sec(e["section"])
        s === nothing && (push!(bad, "$k: section"); continue)
        occursin(e["quote"], s) || push!(bad, "$k: quote"); occursin(e["text"], e["quote"]) || push!(bad, "$k: text")
    end
    @test isempty(bad)
    # the 71-key table is read by the existing D93 reader against this preregistration (it carries the pilot lines and the section A1 without harm)
    q = AN.p3_read_params(I.params; prereg = I.prereg, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest)
    @test q.r_max == 3.0 && q.n_t_min == 3 && q.bootstrap_resamples == 100 && q.gate_size == 2 && q.min_set == 5
    # the record: parse, the checks a to f of E60 hold by construction
    rec = AX.read_json(I.record)
    @test p3ref_sha256_file(I.record) == I.t.calibration_record_sha256 && p3ref_sha256_file(I.params) == I.t.p3_table_sha256 == rec["fingerprints"]["p3_parameters"]
    @test rec["outcome"]["selection"] == "selected" && rec["parameters"]["test_seed"] == 20261030 && rec["parameters"]["r_max"] == 3.0
    @test rec["pp_hb"]["targets"]["n_T"] == I.t.target_t_s && rec["pp_hb"]["targets"]["n_W"] == I.t.target_w2 && rec["pp_hc"]["target"] == I.t.target_f && rec["pp_hb"]["targets"]["cap"] === false
    @test rec["selection"]["theta"]["V"] == ["SALI", "FLI", "MEGNO", "E2"] && rec["selection"]["theta"]["qC"] == 2 && rec["pp_hc"]["power_expected"]["power"] == I.t.power_hc_planned
    @test rec["pp_hb"]["path"][1]["power"] == I.t.power_hb_planned && rec["pp_hb"]["model"]["rule"]["T"]["p"] == 0.1 && rec["pp_hb"]["model"]["comp"]["SALI"]["RW"]["s"] == 0.3
end

@testset "D95 fixture: the chunk files, the panel, the manifest and the instance" begin
    S = p3t_cached(); I = p3t_cached_inputs()
    inst = p3t_instance(S, I, "fx")
    files = p3t_files(inst.store)
    @test length(files) == 50
    attrs(f) = AX.read_chunk_attrs(f)
    refs = [f for f in files if "ref" in attrs(f)["arms"]]
    @test length(refs) == 24 && all(f -> occursin(r"cell=k\d:r\d$", attrs(f)["chunk_key"]) && attrs(f)["complete"] === true && attrs(f)["study"] == P3_SLUG && attrs(f)["prereg_sha256"] == I.prereg_hash &&
                                         attrs(f)["panel_sha256"] == S.panel_hash && length(attrs(f)["chunk_hash"]) == 64, refs)
    @test count(f -> attrs(f)["status"] == "failed", files) == 6          # the ref cell (2, 5) twice (L1, L5), the ind cell 5 four times
    man = AX.read_json(inst.manifest)
    @test man["status"] == "complete" && man["study"] == P3_SLUG && sort(collect(keys(man["arms"]))) == ["chk", "d64", "ind", "ref"] && man["prereg"]["header"] == I.prereg_hash && length(man["chunks"]) == 50
    pan = TOML.parsefile(inst.panel)
    @test pan["manifest"]["r_max"] == 3 && length(pan["orbit"]) == 31 && [o["round"] for o in pan["orbit"]] == [o.round for o in S.orbits] && count(o -> o["in_D"], pan["orbit"]) == 16 &&
          count(o -> o["K"] == 0.0 && o["round"] == 1, pan["orbit"]) == 2
    # the rounds by hand (see the test plan): counts after round 1, 2, 3
    hand = Dict(1 => (orbits = 12, excluded = 0, drift = 1, t_s = 2, w2 = 4), 2 => (orbits = 22, excluded = 1, drift = 1, t_s = 4, w2 = 7), 3 => (orbits = 29, excluded = 1, drift = 1, t_s = 5, w2 = 11))
    for r in 1:3
        c = p3t_counts(S, I.qt, r); h = hand[r]
        @test (c["orbits"], c["excluded"], c["drift"], c["t_s"], c["w2"], c["target_t_s"], c["target_w2"]) == (h.orbits, h.excluded, h.drift, h.t_s, h.w2, 5, 11)
    end
    ents, decs = p3t_log_entries(S, I)
    @test decs == ["continue", "continue", "close"] && all(e -> count(==('\n'), e) == 13, ents) && occursin("prev_sha256=e3b0c442", ents[1]) && occursin("prev_sha256=$(p3ref_sha256(ents[1]))", ents[2])
    # the instance of 2 rounds holds the ref chunks of rounds 1 and 2 only
    i2 = p3t_instance(S, I, "fx2"; rounds = 2, closure = false)
    @test length(p3t_files(i2.store)) == 16 && all(f -> "ref" in attrs(f)["arms"], p3t_files(i2.store))
end

@testset "D95 fixture: the expected sets, outputs and statistics from the design" begin
    S = p3t_cached(); I = p3t_cached_inputs()
    X = p3t_expected(S, I.qt)
    @test length(X.E) == 25 && X.sizes == (Ts = 5, W2 = 10) && X.ha.strata == [5, 5, 7, 1, 7] && sum(X.ha.strata) == 25
    @test [(d.id, d.reason) for d in X.drops] == [(14, "drift"), (23, "unresolved"), (24, "excluded"), (26, "no verdict data")]
    @test count(r -> r.S2, X.E) == 7 && [r.id for r in X.E if r.Ts] == [3, 4, 15, 16, 25] && [r.id for r in X.E if r.W2 && !r.Ts] == [5, 13, 17, 28, 30] && X.ha.point.g == 0
    @test X.ha.decision == "not supported" && X.ha.failed == [2, 3] && X.ha.point.dw == 1 // 35 && X.ha.point.dt == 3 // 7 && X.ha.point.e == 0 && X.ha.point.aE == 14 // 25
    @test X.hc.failures == 7 && X.hc.nonfailures == 8 && X.hc.decision == "supported" && X.hc.auc_i == 55 // 56 && X.hc.auc_p == 41 // 56 && X.hc.delta == 1 // 4
    @test X.hb.decision == "not supported" && length(X.hb.rows) == 13 && X.M.le_read == 20480 && X.M.X[:SALI].b == 5120 && !X.hb.coverage_ok && X.hb.underpowered
    # the sets are nested: T_s in W2, W2 + S2 + other C = C orbits
    @test all(r -> !r.Ts || r.W2, X.E) && count(r -> r.rho === :C, X.E) == 18 && count(r -> r.rho === :R, X.E) == 7
end
