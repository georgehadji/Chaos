# Amendment D89, stage 02: fixture self-test. Not a test of the module: it does not call Chaoslab.Analysis and passes before the module exists. It shows
# that a RED result elsewhere in the D89 suite is not caused by a broken fixture or a broken reference implementation (stage 02 audit, RED).

@testset "D89 fixture: the parameter table, the preregistration and the run blocks" begin
    @test length(P2_PARAM_KEYS) == 58 && allunique(P2_PARAM_KEYS) && length(P2_MACHINE) == 11 && all(k in P2_PARAM_KEYS for k in P2_MACHINE)
    st = p2_cached("default")
    p = st.p
    @test Set(keys(p)) == Set(P2_PARAM_KEYS)
    tab = TOML.parsefile(st.params)
    @test Set(Symbol.(keys(tab))) == Set(P2_PARAM_KEYS)
    text = read(st.prereg, String)
    @test AX.check_prereg(st.prereg).header_sha256 == st.prereg_hash && occursin("| Study slug | p2-indicators |", text)
    # every non-machine key carries a quote that occurs in its section and contains the value text; machine keys are marked
    function section(name)
        m = match(Regex("(?m)^#+ " * replace(name, r"([()|.\[\]{}^$*+?\\])" => s"\\\1") * "\\n(.*?)(?=^## |\\z)", "s"), text)
        m === nothing ? nothing : m.captures[1]
    end
    bad = String[]
    for k in P2_PARAM_KEYS
        e = tab[string(k)]
        sec = section(e["section"])
        sec === nothing && (push!(bad, "$k: section"); continue)
        if k in P2_MACHINE
            e["machine"] === true || push!(bad, "$k: machine")
        else
            occursin(e["quote"], sec) || push!(bad, "$k: quote not in section")
            v = getfield(p, k)
            vt = v isa AbstractVector ? join(v, ", ") : string(v)
            k in (:gate_classes, :gate_size, :gate_k_index) ? true : (occursin(vt, e["quote"]) || push!(bad, "$k: value text $vt not in quote"))
        end
    end
    @test isempty(bad)
    # the Experiments run-block reader accepts the blocks and holds the values of the machine keys
    rb = AX.read_run_block(st.prereg)
    arm(n) = rb.raw["arm"][findfirst(a -> a["name"] == n, rb.raw["arm"])]
    @test arm("lyap_p")["horizons"] == [p.N_1] && arm("stretch_p")["horizons"] == [p.N_c] && arm("dir_tg")["horizons"] == [p.N_dir]
    @test arm("dir_tg")["direction"] == p.directions && rb.raw["health"]["drift_bound"] == p.drift_bound && rb.raw["panel_prereg_sha256_accepted"] == [p.panel_build_sha256]
    spec(id) = rb.raw["indicator"][findfirst(s -> s["id"] == id, rb.raw["indicator"])]
    @test spec("fq1")["windows"] == p.windows == spec("bk1")["windows"] && 2 .* p.windows == p.budgets && p.budgets[end] == p.N_c == 2 * maximum(p.windows)
    @test p.b_primary in p.budgets && p.budgets[P2_KP] == p.b_primary && length(rb.raw["h4"]) == p.grid_orbit_count == length(p.grid_k_values) * p.grid_i_count * p.grid_j_count
    @test [h["label"] for h in rb.raw["h4"]] == [g.label for g in st.h4] && [h["K"] for h in rb.raw["h4"]] == [g.K for g in st.h4]
    @test p.bootstrap_seed == AN_SEED == TOML.parsefile(joinpath(@__DIR__, "..", "..", "_config", "numerics-defaults.toml"))["seed"]
    # the P1 reference: slug, panel hash, label numbers equal the keys of the P1 table
    p1 = st.p1
    @test AX.read_json(p1.manifest)["study"] == P1_SLUG && AX.read_json(st.manifest)["study"] == P2_SLUG && AX.read_json(p1.manifest)["panel"]["sha256_end"] == st.panel_hash
    q1 = TOML.parsefile(p1.params)
    @test q1["lambda_threshold"] == p.label_lambda_min && q1["growth_ratio"] == p.label_growth_ratio && q1["rotation_tol"] == p.label_rotation_tol && q1["h1a_bound"] == p.label_estimator_bound
    @test length(q1) == 37 && an_sha256(read(st.panel)) == an_sha256(read(p1.panel))
end

@testset "D89 fixture: grids, steps and the store" begin
    st = p2_cached("default")
    p = st.p
    # E26 hand table: the first recorded step not before b_k and the step nearest to n / 10
    @test [p2ref_n(P2_GA, b) for b in p.budgets] == [44, 88, 176, 352, 704, 1408, 2816, 5632, 11264, 22528]
    @test [p2ref_nearest(P2_GA, n / 10) for n in P2_NK] == [4, 9, 18, 35, 70, 141, 282, 563, 1126, 2253]
    @test [p2ref_n(P2_GB, b) for b in p.budgets] == [44, 88, 176, 352, 704, 1408, 2816, 5632, 10240, 20480] == P2_NKB
    @test p2ref_nearest([1, 3], 2) == 3 && p2ref_nearest([1, 2, 4], 3) == 4 && p2ref_nearest([10, 20], 14.9) == 10          # a tie goes to the larger step
    for g in (P2_GA, P2_GAD, P2_GB, P2_GB1)
        @test issorted(g) && allunique(g)
    end
    @test P2_GA[end] == p.N_1 && (p.N_1 ÷ 10) in P2_GA && (p.N_1 ÷ 2) in P2_GA && p.n_s in P2_GA && P2_GAD[end] == p.N_dir && P2_GB[end] == p.N_c && P2_GB1[end] == p.N_c + 1
    @test st.chunks == 33 && st.items == 466 && length(filter(f -> endswith(f, ".h5"), readdir(st.store))) == 33
    man = AX.read_json(st.manifest)
    @test man["status"] == "complete" && man["plan_counts"]["chunks"] == 33 && man["plan_counts"]["items"] == 466 && man["git"]["clean"] === true && isempty(man["control_stops"])
    # one chunk of every kind read back by the real reader
    ch = AX.read_chunk(p2_chunk(st, 2))                               # tg, panel
    @test ch.attrs["kind"] == "tangent" && ch.attrs["method"] == "tg" && ch.attrs["coordinate_mode"] == "wrapped" && isempty(ch.orbits[1].lyap)
    tb = ch.orbits[1].indicators["tg"]
    @test size(tb.series) == (length(P2_GA), 10) && tb.columns[end] == "SALI.L" && ch.orbits[1].t == P2_GA && size(ch.orbits[1].winding) == (length(P2_GA), 2)
    ch = AX.read_chunk(p2_chunk(st, 10))                              # fq1, panel without Y5
    @test length(ch.orbits) == 19 && !any(o -> o.id == 20, ch.orbits) && size(ch.orbits[1].indicators["fq1"].final) == (10, 11) && ch.orbits[1].indicators["fq1"].rows[1] == "W=20"
    ch = AX.read_chunk(p2_chunk(st, 11))                              # bk1
    @test size(ch.orbits[1].indicators["bk1"].final) == (40, 8) && ch.orbits[1].indicators["bk1"].rows[1:2] == ["W=20 cos(1 x[1])", "W=20 sin(1 x[1])"]
    ch = AX.read_chunk(p2_chunk(st, 14))                              # pa1 horizon N_c + 1
    @test ch.attrs["horizon"] == p.N_c + 1 && ch.orbits[1].t[end] == p.N_c + 1 && size(ch.orbits[1].indicators["pa1"].final) == (10, 4)
    ch = AX.read_chunk(p2_chunk(st, 4))                               # L5 on the grid: two exponents
    @test size(ch.orbits[1].lyap["L5"]) == (length(P2_GA), 2) && ch.orbits[1].id == 21 && ch.orbits[1].attrs["ic"] == "$(st.h4[1].x0),$(st.h4[1].p0)" && ch.orbits[1].attrs["params"] == "0.5"
    @test AX.read_chunk(p2_chunk(st, 1)).orbits[17].attrs["control_scope"] == "bound" && AX.read_chunk(p2_chunk(st, 1)).orbits[18].attrs["control_scope"] == "info"
    # flag words: the bit of the j-th name is 2^j; the manifest counts orbits with the flag in at least one row
    fl = [e for e in man["chunks"] if haskey(e, "indicator_flags")]
    @test !isempty(fl) && all(e -> e["indicator_flags"] isa AbstractDict, fl) && man["chunks"][10]["not_applicable"] == Dict("fq1" => [20]) && man["chunks"][11]["not_applicable"] == Dict("bk1" => [20])
    @test st.p1.chunks == 55
end

@testset "D89 fixture: the stored numbers realise the designed verdicts (R1) and the designed sets" begin
    st = p2_cached("default")
    ref = p2_ref_all(st)
    d = st.design
    bad = [(id, X) for ((id, X), v) in ref.V1 if v != d.V[id][X]]
    @test isempty(bad) && length(ref.V1) == sum(length(d.V[id]) for id in 1:44)
    @test count(k -> k[2] === :RLI, keys(ref.V1)) == 4 + 24 && all(!haskey(ref.V1, (20, X)) for X in (:E1, :WBA)) && all(haskey(ref.V1, (i, :LE)) for i in 1:44)
    # the sets, derived by hand from the design (lambda_ref and lambda_ns of the plan, "Fixture design")
    S = p2ref_sets(p2_ref_orbits(st), st.p)
    @test S.E_all == [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]
    @test S.W2 == [1, 3, 5, 7, 8, 21, 22, 23, 24, 26, 39] && S.S2 == [11, 12, 33, 34, 35, 36, 37, 38] && S.T_s == [5, 23, 24, 39] && S.T_s_design == [5, 6]
    @test S.regular == [9, 10, 15, 16, 27, 28, 29, 30, 31, 41, 42, 43] && S.gate == [9, 10] && S.S_D == [1, 5, 7, 9, 11, 12, 15] && S.S_L == [1, 5, 9, 11] && isempty(S.overlap)
    @test S.E8 == vcat([1, 5, 9, 11], 21:31, 33:43) && length(S.E8) == 26
    @test [(x.id, x.reason) for x in S.dropped] == [(2, "unresolved"), (13, "unresolved"), (17, "control"), (18, "control"), (19, "control"), (20, "control"), (32, "unresolved"), (44, "unresolved")]
    # headline values of the default scenario, by hand (arithmetic in the test plan): medians of the censored times on S2, R_X on W2 over S2
    rW, rS = p2_ref_recs(ref, S.W2, P2_X7), p2_ref_recs(ref, S.S2, P2_X7)
    M = p2ref_stat(:H2, (rS,), st.p).M
    @test M == Dict(:LE => 40.0, :SALI => 40.0, :FLI => 40.0, :MEGNO => 80.0, :E1 => 640.0, :E2 => 160.0, :WBA => 320.0) && p2ref_stat(:H2, (rS,), st.p).Q == 16.0
    @test p2ref_stat(:H3, (rW, rS), st.p).R == Dict(:LE => 64.0, :SALI => 16.0, :FLI => 16.0, :MEGNO => 32.0, :E1 => 8.0, :E2 => 16.0, :WBA => 128.0)
    h1 = p2ref_stat(:H1, (rW,), st.p)
    @test h1.n == 11 && h1.used == 11 && h1.U == 1.0 && h1.A_min == 0.0
    F = p2ref_stat(:H5, (rW,), st.p).F
    @test F[:LE] == 4 / 11 && F[:E1] == 7 / 11 && F[:WBA] == 7 / 11 && F[:SALI] == 0.0 && p2ref_stat(:H5, (rW,), st.p).F_med == 0.0
    # detection of two orbits by hand: orbit 5 (T profile), LE "rrrrrrrccc": stable from k = 7 (b = 2560 x 2 = 5120 = N_c / 4), correct; SALI all regular: wrong
    @test ref.Det[(5, :LE)] == (t_det = 5120, outcome = :correct, censored = 5120) && ref.Det[(5, :SALI)] == (t_det = 40, outcome = :wrong, censored = 40960)
    @test ref.Det[(15, :E1)].outcome === :undetected && ref.Det[(15, :E1)].censored == 40960 && ref.Det[(9, :E1)] == (t_det = 320, outcome = :correct, censored = 320)
end

@testset "D89 fixture: reference implementation against hand cases (O7)" begin
    # kappa (O7 a): perfect agreement 1; independent table 0; p_e = 1 gives 0; a = d = 0 with b = c gives -1
    tab(a, b, c, d) = p2ref_table(vcat(fill(:chaotic, a + b), fill(:regular, c + d)), vcat(fill(:chaotic, a), fill(:regular, b), fill(:chaotic, c), fill(:regular, d)))
    @test tab(10, 0, 0, 10).kappa == 1.0 && tab(25, 25, 25, 25).kappa == 0.0 && tab(10, 0, 0, 0).kappa == 0.0 && tab(0, 7, 7, 0).kappa == -1.0
    # a hand table: [[20, 5], [10, 15]]: p_o = 35/50 = 0.7, p_e = (25 x 30 + 25 x 20) / 2500 = 0.5, kappa = 0.2 / 0.5 = 0.4
    t = tab(20, 5, 10, 15)
    @test t.n == 50 && t.A == 0.7 && isapprox(t.kappa, 0.4; atol = 1e-15)
    @test p2ref_table([:chaotic, :na], [:regular, :regular]).n == 1 && p2ref_table(Symbol[], Symbol[]).A === nothing
    # AUC (O7 b): separated 1; identical score sets 0.5; C = {1, 2, 3} against R = {2} gives (1 + 0.5 + 0) / 3 = 0.5; one tied pair counts one half
    @test p2ref_auc([5, 6], [1, 2]) == 1.0 && p2ref_auc([1, 2, 3], [1, 2, 3]) == 0.5 && p2ref_auc([1, 2, 3], [2]) == 0.5 && p2ref_auc([2, 3], [2, 1]) == (0.5 + 1 + 1 + 1) / 4
    @test p2ref_auc([1, 2], [3, 4]) == 0.0 && p2ref_auc([1.0, 2.0], [0.5]; min_set = 2) === nothing && p2ref_auc(fill(1.0, 20), fill(0.0, 19); min_set = 20) === nothing
    @test p2ref_auc(fill(1.0, 20), fill(0.0, 20); min_set = 20) == 1.0 && p2ref_auc([-Inf, 1.0], [-Inf]) == 0.75
    # midranks and Spearman: ties get the mean of their ranks
    @test p2ref_midranks([10.0, 20.0, 20.0, 30.0]) == [1.0, 2.5, 2.5, 4.0] && p2ref_spearman([1.0, 2, 3, 4], [1.0, 2, 3, 4]) == 1.0 && p2ref_spearman([1.0, 2, 3, 4], [4.0, 3, 2, 1]) == -1.0
    @test p2ref_spearman([1.0, 1, 1], [1.0, 2, 3]) === nothing && p2ref_spearman([1.0, 2], [1.0, 2]) === nothing
    # medians and type 7 percentiles
    @test p2ref_median([3, 1, 2]) == 2.0 && p2ref_median([4, 1, 3, 2]) == 2.5 && p2ref_percentile(collect(1.0:11.0), 0.05) == 1.5 && p2ref_percentile([1.0, 2.0], 0.5) == 1.5
    # detection by hand (N_c / 4 = 5120 = b_7): stable from k = 7 allowed, from k = 8 not; a flip at k = 8 gives no earlier stability; undecided at k = 9 is undetected
    p = p2_params()
    sq(s) = p2_sym(s)
    @test p2ref_detection(sq("rrrrrrrccc"), :C, p) == (t_det = 5120, outcome = :correct, censored = 5120)
    @test p2ref_detection(sq("rrrrrrrrcc"), :C, p) == (t_det = nothing, outcome = :undetected, censored = 40960)
    @test p2ref_detection(sq("cccccccrcc"), :C, p) == (t_det = nothing, outcome = :undetected, censored = 40960)
    @test p2ref_detection(sq("ccccccccru"), :C, p).outcome === :undetected && p2ref_detection(sq("cccccccccu"), :C, p).t_det === nothing
    @test p2ref_detection(sq("uurrrrrrrr"), :R, p) == (t_det = 160, outcome = :correct, censored = 160) && p2ref_detection(sq("rrrrrrrrrr"), :C, p) == (t_det = 40, outcome = :wrong, censored = 40960)
    # the hypothesis decisions: the strict and the non-strict comparisons of E32
    @test p2ref_decide(:H2, (n = 20, Q = 10.0, Q_lo = 3.0), p) === :not_supported && p2ref_decide(:H2, (n = 20, Q = 10.0, Q_lo = nextfloat(3.0)), p) === :supported
    @test p2ref_decide(:H2, (n = 4, Q = 10.0, Q_lo = 5.0), p) === :not_decidable && p2ref_decide(:H1, (n = 5, U = 0.2, U_lo = 0.1, A_min = 0.8), p) === :supported
    # the stream rule: SHA-256 of the label text, first 8 bytes big-endian
    @test p2_stream_seed("p2|H1|W2|E_all|R1|all", AN_SEED) == reinterpret(UInt64, reverse(hex2bytes(an_sha256("chaoslab-analysis-1|20261012|p2|H1|W2|E_all|R1|all"))[1:8]))[1]
    @test p2_draws("p2|H2|S2|E_all|R1|all", 8, 3, AN_SEED) == p2_draws("p2|H2|S2|E_all|R1|all", 8, 3, AN_SEED) && size(p2_draws("x", 8, 5, 1)) == (8, 5)
    @test p2_draws("a", 8, 3, 1) != p2_draws("b", 8, 3, 1)
end
