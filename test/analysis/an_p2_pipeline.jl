# D89 pipeline on the synthetic study "p2-indicators" (default scenario): reference class, identities, sets, verdicts, scores, detection, statistics and
# decisions H1 to H6, the extension and ROBUST flags, the R2 version, the gate, AUC and regular range. Spec: analysis-spec.md, Amendment D89.
# Expected values come from the design of an_p2_study.jl through the independent reference of an_p2_ref.jl / an_p2_refstore.jl (hand values in the plan).

const P2_RES = Dict{String, Any}()
p2_res(tag, st; kw...) = get!(() -> p2_run(st; kw...), P2_RES, tag)
ref_u_s2(st, ref, S) = p2ref_stat(:H1, (p2_ref_recs(ref, S.S2, P2_X7),), st.p).U
const P2_REF = Dict{String, Any}()
function p2_refs(tag, st)
    get!(P2_REF, tag) do
        ref = p2_ref_all(st)
        (; ref, S = p2ref_sets(p2_ref_orbits(st), st.p))
    end
end
# fields of the module's decision tuples that differ from the reference tuples
function p2_decision_diff(dec, rf)
    bad = String[]
    for h in (:H1, :H2, :H3, :H5)
        append!(bad, ["$h.$k" for k in p2_diff(getproperty(dec, h), getproperty(rf, h))])
    end
    append!(bad, ["H4.$k" for k in p2_diff(dec.H4, rf.H4)])
    append!(bad, ["H4.design.$k" for k in p2_diff(dec.H4.design, rf.H4design)])
    h6 = dec.H6
    p2_same(h6.dir.n, rf.H6.dir.n) || push!(bad, "H6.dir.n")
    p2_same(h6.dir.delta, rf.H6.dir.delta) || push!(bad, "H6.dir.delta")
    p2_same(h6.dir.dropped, rf.H6.dir.dropped) || push!(bad, "H6.dir.dropped")
    p2_same(h6.prec.n, rf.H6.prec.n) || push!(bad, "H6.prec.n")
    p2_same(h6.prec.delta, rf.H6.prec.delta) || push!(bad, "H6.prec.delta")
    p2_same(h6.prec.dropped, rf.H6.prec.dropped) || push!(bad, "H6.prec.dropped")
    h6.verdict === rf.H6.verdict || push!(bad, "H6.verdict")
    return bad
end

@testset "D89 default scenario: status, outputs, reference class, grid identity, lambda_ref" begin
    st = p2_cached("default")
    res = p2_res("default", st)
    d = st.design
    @test res.status === :ok && res.study == "p2-indicators" && isfile(res.tables) && basename(res.tables) == "p2-indicators-tables.md"
    @test an_files(st.out) == ["p2-indicators-tables.md"]                                   # figures not requested: the tables file and nothing else
    @test Set(keys(res.rho)) == Set(1:44) && all(id -> res.rho[id] === d.rho[id], 1:44)
    # the panel labels are the P1 labels of the base pipeline on the P1 reference inputs (I21)
    p1 = AN.analyze(; p1_args(st.p1)..., out = joinpath(st.dir, "p1-out"))
    @test all(id -> (res.p1_labels[id].instability, res.p1_labels[id].structure) == (p1.labels[id].instability, p1.labels[id].structure), 1:20)
    pos, none, unr = "positive growth supported", "none detected on the tested horizon", "unresolved"
    @test all(id -> res.rho[id] === (p1.labels[id].instability == pos && p1.labels[id].structure == "irregular" ? :C :
                                      p1.labels[id].instability == none && p1.labels[id].structure == "quasiperiodic" ? :R : :U), 1:20)
    # the grid: identity by the table order, labels, K index, i, j; the reduced rule gives the designed class; every grid orbit carries ref_weak
    g = Dict(r.id => r for r in res.grid)
    @test length(res.grid) == 24 && Set(keys(g)) == Set(21:44)
    @test all(r -> g[r.id].label == r.label && g[r.id].kindex == r.kindex && g[r.id].i == r.i && g[r.id].j == r.j && g[r.id].K == parse(Float64, r.K), st.h4)
    @test all(i -> g[i].rho === d.rho[i] && g[i].ref_weak === true, 21:44) && g[21].label == "grid-k1-i00-j00" && g[44].label == "grid-k2-i02-j03" && g[33].kindex == 2
    @test (g[32].instability, g[32].structure) == (none, unr) && (g[44].instability, g[44].structure) == (unr, unr) && (g[27].instability, g[27].structure) == (none, "quasiperiodic")
    @test (g[21].instability, g[21].structure) == (pos, "irregular")
    # lambda_ref and lambda_ns are the values of the P2 L1 item at N_1 and n_s (P13)
    @test all(i -> res.lambda_ref[i] == d.lref[i] && res.lambda_ns[i] == d.lns[i], 1:44)
end

@testset "D89 default scenario: the sets (E29, I16)" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    @test p2_diff(res.sets, S) == Symbol[]
    @test res.sets.W2 == [1, 3, 5, 7, 8, 21, 22, 23, 24, 26, 39] && res.sets.S2 == [11, 12, 33, 34, 35, 36, 37, 38] && res.sets.T_s == [5, 23, 24, 39]
    @test res.sets.T_s_design == [5, 6] && length(res.sets.E_all) == 36 && length(res.sets.E8) == 26 && res.sets.gate == [9, 10] && isempty(res.sets.overlap)
    @test issubset(res.sets.T_s, res.sets.W2) && issubset(res.sets.E8, res.sets.E_all) && isempty(intersect(res.sets.W2, res.sets.S2))
    @test sort(vcat(res.sets.E_all, [x.id for x in res.sets.dropped])) == collect(1:44)
    # T1: every orbit of P and G once, by source and K index; the counts of rho add up; hand values for the rho counts
    T1 = res.items.T1
    col(nm) = findfirst(==(nm), T1.columns)
    row(src, k) = T1.rows[findfirst(r -> r[col("source")] == src && r[col("K index")] == k, T1.rows)]
    @test sum(r[col("orbits")] for r in T1.rows) == 44 && length(T1.rows) == 5 && all(r -> r[col("C")] + r[col("R")] + r[col("U")] == r[col("orbits")], T1.rows)
    cnt(src, k) = Tuple(row(src, k)[col(nm)] for nm in ("orbits", "C", "R", "U"))
    @test cnt("panel", 1) == (13, 8, 2, 3) && cnt("panel", 2) == (6, 3, 2, 1) && cnt("panel", 0) == (1, 1, 0, 0) && cnt("grid", 1) == (12, 6, 5, 1) && cnt("grid", 2) == (12, 8, 3, 1)
    @test sum(r[col("C")] for r in T1.rows) == 26 && sum(r[col("R")] for r in T1.rows) == 12 && sum(r[col("U")] for r in T1.rows) == 6
end

@testset "D89 default scenario: verdicts R1 and R2, scores, steps, detection against the reference" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st); d = st.design
    # the verdict series of every orbit and indicator that has an item; no entry for the others (RLI outside S_L and G, E1 and WBA for the Y5 control)
    keys1 = Set((id, X, :R1) for (id, X) in keys(ref.V1)); keys2 = Set((id, X, :R2) for (id, X) in keys(ref.V2))
    @test Set(keys(res.verdicts)) == union(keys1, keys2)
    @test all(((id, X),) -> res.verdicts[(id, X, :R1)] == ref.V1[(id, X)] == d.V[id][X], keys(ref.V1))
    @test all(((id, X),) -> res.verdicts[(id, X, :R2)] == ref.V2[(id, X)], keys(ref.V2))
    @test count(((id, X),) -> ref.V1[(id, X)] != ref.V2[(id, X)], keys(ref.V1)) > 20                    # R2 is a different rule: the series differ where the thresholds do
    @test !haskey(res.verdicts, (7, :RLI, :R1)) && haskey(res.verdicts, (1, :RLI, :R1)) && haskey(res.verdicts, (30, :RLI, :R1)) && !haskey(res.verdicts, (20, :E1, :R1)) && haskey(res.verdicts, (20, :LE, :R1))
    @test p2_str(res.verdicts[(7, :LE, :R1)]) == "rrrrrccccc" && p2_str(res.verdicts[(10, :SALI, :R1)]) == "uurrrrrrrr" && p2_str(res.verdicts[(15, :WBA, :R1)]) == "uuuuuuuuuu"
    # scores of E31 (nothing where unavailable, minus infinity kept for RLI = 0), steps used (n(b) on the item grid, W for the window items)
    @test all(((id, X),) -> p2_same(res.scores[(id, X)], ref.S[(id, X)]), keys(ref.S)) && Set(keys(res.scores)) == Set(keys(ref.S))
    @test res.steps[(7, :LE)] == P2_NK && res.steps[(7, :SALI)] == P2_NK && res.steps[(1, :RLI)] == P2_NKB && res.steps[(7, :E1)] == P2_WINDOWS && res.steps[(7, :STR)] == P2_NKB
    @test all(((id, X),) -> res.steps[(id, X)] == ref.N[(id, X)], keys(ref.N))
    # the detection of the orbits of E_all (rho C or R): t_det, outcome, censored time; hand values of the plan
    @test Set(keys(res.detection)) == Set((id, X) for id in S.E_all for X in P2_X8 if haskey(ref.V1, (id, X)))
    @test all(((id, X),) -> res.detection[(id, X)] == ref.Det[(id, X)], keys(res.detection))
    @test res.detection[(5, :LE)] == (t_det = 5120, outcome = :correct, censored = 5120) && res.detection[(5, :SALI)] == (t_det = 40, outcome = :wrong, censored = 40960)
    @test res.detection[(15, :E1)] == (t_det = nothing, outcome = :undetected, censored = 40960) && res.detection[(9, :E1)] == (t_det = 320, outcome = :correct, censored = 320)
    @test res.detection[(16, :LE)] == (t_det = nothing, outcome = :undetected, censored = 40960) && res.detection[(28, :MEGNO)].outcome === :undetected
    @test res.detection[(29, :FLI)] == (t_det = nothing, outcome = :undetected, censored = 40960) && res.detection[(41, :SALI)] == (t_det = 160, outcome = :correct, censored = 160)
    # t~ takes only the values b_0 .. b_7 and 2 N_c (derived), the orbits that are not in E_all have no detection
    @test all(v -> v.censored in vcat(P2_BUDGETS[1:8], [2 * st.p.N_c]), values(res.detection)) && all(((id, X),) -> id in S.E_all, keys(res.detection))
end

@testset "D89 default scenario: the six decisions in order, statistics and intervals" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    @test collect(keys(res.decisions)) == [:H1, :H2, :H3, :H4, :H5, :H6]                  # decided in this order (I15)
    rf = p2_ref_hyp(st, ref, S, P2_X7)
    @test p2_decision_diff(res.decisions, rf) == String[]
    h = res.decisions
    # hand values (plan, "Hand-computed expected values"): W2 11 orbits, S2 8, T_s 4, T_s_design 2
    @test (h.H1.verdict, h.H1.n, h.H1.used, h.H1.dropped, h.H1.U, h.H1.U_lo, h.H1.A_min) == (:supported, 11, 11, 0, 1.0, 1.0, 0.0) && length(h.H1.A) == 21
    @test p2_same(h.H1.U_S2, ref_u_s2(st, ref, S))                                      # U on S2 is reported with no rule
    @test h.H2.verdict === :supported && h.H2.n == 8 && h.H2.M == Dict(:LE => 40.0, :SALI => 40.0, :FLI => 40.0, :MEGNO => 80.0, :E1 => 640.0, :E2 => 160.0, :WBA => 320.0) && h.H2.Q == 16.0 && h.H2.Q_lo == 16.0
    @test h.H3.verdict === :supported && (h.H3.n_w, h.H3.n_s) == (11, 8) && h.H3.R == Dict(:LE => 64.0, :SALI => 16.0, :FLI => 16.0, :MEGNO => 32.0, :E1 => 8.0, :E2 => 16.0, :WBA => 128.0) && h.H3.passing == 6
    @test h.H4.verdict === :not_decidable && h.H4.n == 4 && h.H4.dphi == 1.0 && h.H4.design.verdict === :not_decidable && h.H4.design.n == 2    # T_s has 4 < min_set = 5 orbits
    @test h.H5.verdict === :not_supported && h.H5.n == 11 && h.H5.F_med == 0.0 && isapprox(h.H5.F[:LE], 4 / 11) && isapprox(h.H5.F[:E1], 7 / 11) && isapprox(h.H5.F[:WBA], 7 / 11)
    @test h.H6.verdict === :supported && h.H6.dir.n == 7 && h.H6.prec.n == 4 && isapprox(h.H6.dir.delta[:LE], 1 / 7) && isapprox(h.H6.dir.delta[:SALI], 1 / 7) && h.H6.dir.delta[:FLI] == 0.0
    @test h.H6.dir.dropped[:FLI] == 1 && h.H6.prec.delta[:LE] == 0.25 && h.H6.prec.delta[:E2] == 0.25 && h.H6.prec.delta[:E1] == 0.0       # Double64: orbit 11 LE and orbit 5 E2 differ
end
@testset "D89 default scenario: extension on E8 with eight indicators, ROBUST flags, the R2 version" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    ext = p2_ref_hyp(st, ref, S, P2_X8; suffix = "_8", srule = "E8")
    @test p2_decision_diff(res.extension, ext) == String[] && collect(keys(res.extension)) == [:H1, :H2, :H3, :H4, :H5, :H6]
    @test length(res.extension.H1.A) == 28 && !haskey(res.extension.H4.phi, :LE) && length(res.extension.H4.phi) == 7 && res.extension.H1.n == 8 && res.extension.H2.n == 7
    # ROBUST iff the primary and the extension verdicts are both decidable and equal
    @test res.robust == (H1 = :ROBUST, H2 = :ROBUST, H3 = :ROBUST, H4 = :NOT_COMPARABLE, H5 = :ROBUST, H6 = :NOT_COMPARABLE)
    @test collect(keys(res.robust)) == [:H1, :H2, :H3, :H4, :H5, :H6]
    # R2 (E33, owner decision of the stage 03 checkpoint): every statistic and decision of E32 recomputed with the verdicts of R2, the detection times (hence M_X, Q, R_X) and
    # the H6 verdicts of the reference store included; the sets are those of E29; the primary decisions stay those of R1
    r2 = p2_ref_hyp(st, ref, S, P2_X7; rule = :R2)
    @test p2_decision_diff(res.r2, r2) == String[] && res.decisions.H2.verdict === :supported && collect(keys(res.r2)) == [:H1, :H2, :H3, :H4, :H5, :H6]
    @test res.r2.H2.M != res.decisions.H2.M && res.r2.H3.R != res.decisions.H3.R && res.r2.H6.dir.delta != res.decisions.H6.dir.delta          # R2 changes the times and the H6 deltas
    # independent hand computation of a few R2 numbers (design of an_p2_study.jl; R2: LE chaotic iff n lambda >= 3 ln n, SALI iff sali <= 1e-8).
    # LE: the stored n lambda at n_k is r x 0.2 x n'_k, so R2 needs r >= 15 ln(n_k) / n'_k = 14.19, 7.46, 4.31, 2.51, 1.41, 0.77 at k = 0 .. 5: a chaotic row (r >= 8)
    # is chaotic from k = 1 (8.4 < 14.19 at k = 0), a regular row (r about 1.2 to 1.45) from k = 5 (from k = 4 when id mod 5 = 4: 1.44 >= 1.405).
    #   S2 (8 orbits): every LE onset is 1 in R2, t~ = b_1 = 80, M_LE(S2) = 80 (R1: 40). SALI: the regular rows of ids 12, 33, 36 (id mod 3 = 0) hold sali = 1e-9, chaotic under R2:
    #   t~ = 40 for 11, 12, 33, 34, 35, 36 and 80 for 37, 38: the median is 40.
    #   W2 (11 orbits), LE onsets in R2: 5 for 1, 3, 7, 8, 21, 22, 26 (W and W2 profiles), 5 for 5 and 23, 4 for 24 and 39 (id mod 5 = 4): t~ = 1280 nine times and 640 twice: M_LE(W2) = 1280.
    #   R_LE = 1280 / 80 = 16 (R1: 2560 / 40 = 64).
    @test res.r2.H2.M[:LE] == 80.0 && res.r2.H2.M[:SALI] == 40.0 && res.r2.H3.M_w[:LE] == 1280.0 && res.r2.H3.R[:LE] == 16.0 && res.decisions.H3.R[:LE] == 64.0
    @test p2_str(res.verdicts[(11, :LE, :R2)]) == "rccccccccc" && p2_str(res.verdicts[(24, :LE, :R2)]) == "rrrrcccccc" && p2_str(res.verdicts[(36, :SALI, :R2)]) == "cccccccccc"
    # H6 under R2 (verdicts at k = 9 of the direction and Double64 items): the regular LE row of orbit 12 in direction 5 is chaotic under R2 (r >= 1.2 > 0.067), so there is no flip,
    # delta_dir(LE) = 0 (R1: 1/7); the SALI flip of orbit 7 in direction 3 stays (sali = 0.1, id mod 3 = 1): 1/7; the undecided FLI row of orbit 15 in direction 2 holds fli = 3,
    # regular under R2: nothing dropped. Double64: the LE flip of orbit 11 disappears the same way (delta_prec(LE) = 0), the E2 flip of orbit 5 stays (D -1 against -9, same rule): 1/4.
    @test res.r2.H6.dir.delta[:LE] == 0.0 && isapprox(res.r2.H6.dir.delta[:SALI], 1 / 7) && res.r2.H6.dir.dropped[:FLI] == 0 && res.r2.H6.prec.delta[:LE] == 0.0 && res.r2.H6.prec.delta[:E2] == 0.25
    # the R2 version of the extension (I_8 on E8)
    r2e = p2_ref_hyp(st, ref, S, P2_X8; rule = :R2, suffix = "_8", srule = "E8")
    @test p2_decision_diff(res.r2_extension, r2e) == String[] && collect(keys(res.r2_extension)) == [:H1, :H2, :H3, :H4, :H5, :H6]
    # the gate: no gate orbit is chaotic at k = 9; RLI exists for 1 of the 2; no second version
    gate = p2_ref_gate(st, ref, S)
    @test p2_same(collect(res.gate), gate) && all(r -> !r.violated, res.gate) && isempty(res.second) && Set(r.X for r in res.gate) == Set(P2_X8)
    @test (r = first(r for r in res.gate if r.X === :RLI); (r.decided, r.na, r.share) == (1, 1, 0.0)) && (r = first(r for r in res.gate if r.X === :LE); (r.decided, r.chaotic, r.na) == (2, 0, 0))
end

@testset "D89 default scenario: AUC, regular range, leaving budget (E31, M10)" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st); p = st.p
    auc = p2_ref_auc(st, ref, S)
    @test Set(keys(res.auc)) == Set(keys(auc)) && all(k -> p2_same(res.auc[k], auc[k]), keys(auc))
    # hand counts: K index 1 has 13 chaotic and 7 regular orbits of E_all, K index 2 has 11 and 5; by source only the grid at K 1 has both classes at 5 or more
    @test (a = res.auc[(:LE, 6, 1, "pooled")]; (a.n_C, a.n_R) == (13, 7)) && (a = res.auc[(:LE, 6, 2, "pooled")]; (a.n_C, a.n_R) == (11, 5)) && res.auc[(:LE, 6, 1, "pooled")].auc !== nothing
    @test res.auc[(:LE, 6, 1, "panel")].auc === nothing && res.auc[(:LE, 6, 1, "grid")].auc !== nothing && res.auc[(:LE, 6, 2, "grid")].auc === nothing
    # RLI exists for S_L and G: 8 chaotic and 6 regular orbits at K 1 (computed), 9 and 3 at K 2 (fewer than min_set regular orbits: not computed); STR enters the AUC
    @test (a = res.auc[(:RLI, 6, 1, "pooled")]; (a.n_C, a.n_R) == (8, 6) && a.auc !== nothing) && (a = res.auc[(:RLI, 6, 2, "pooled")]; (a.n_C, a.n_R) == (9, 3) && a.auc === nothing)
    @test haskey(res.auc, (:STR, 9, 1, "pooled")) && res.auc[(:STR, 9, 1, "pooled")].auc == 1.0                  # every chaotic D (0.02 ..) is below every regular D (0.5 ..): the score -log10 D is higher
    @test all(a -> a.auc === nothing || 0.0 <= a.auc <= 1.0, values(res.auc))
    # regular range: the 5th and 95th percentile of the scores of the regular set at that K (type 7), for the indicators of I_7 and STR, when at least min_set orbits exist
    orbs = Dict(o.id => o for o in p2_ref_orbits(st))
    bad = String[]
    for X in vcat(P2_X7, :STR), kidx in 1:2, k in 0:9
        sc = Float64[ref.S[(id, X)][k + 1] for id in S.regular if orbs[id].kindex == kidx && haskey(ref.S, (id, X)) && ref.S[(id, X)][k + 1] !== nothing]
        want = length(sc) >= p.min_set ? (lo = p2ref_percentile(sc, 0.05), hi = p2ref_percentile(sc, 0.95), median = p2ref_median(sc), n = length(sc)) : nothing
        got = res.range[(X, kidx, k)]
        (want === nothing ? got === nothing : (got !== nothing && p2_same(got, want))) || push!(bad, "$X $kidx $k")
    end
    @test isempty(bad)
    # leaving budget of the sticky orbits (T_s) and the share below the regular median at k = 9
    badl = String[]
    for id in S.T_s, X in vcat(P2_X7, :STR)
        kidx = orbs[id].kindex
        haskey(ref.S, (id, X)) || continue
        up = [res.range[(X, kidx, k)] === nothing ? nothing : res.range[(X, kidx, k)].hi for k in 0:9]
        sc = ref.S[(id, X)]
        want = nothing
        if all(!isnothing, up)
            for k in 0:9
                all(j -> sc[j + 1] !== nothing && sc[j + 1] > up[j + 1], k:9) && (want = p.budgets[k + 1]; break)
            end
        end
        res.leaving[(id, X)] === want || push!(badl, "$id $X")
    end
    @test isempty(badl)
    for X in vcat(P2_X7, :STR), kidx in 1:2
        ids = [id for id in S.T_s if orbs[id].kindex == kidx && haskey(ref.S, (id, X)) && ref.S[(id, X)][10] !== nothing]
        rg = res.range[(X, kidx, 9)]
        want = (rg === nothing || isempty(ids)) ? nothing : count(id -> ref.S[(id, X)][10] < rg.median, ids) / length(ids)
        @test p2_same(res.below_median[(X, kidx)], want)
    end
end
