# D93d to D93f, D93j to D93l, D93n to D93p end to end: `calibrate` on the synthetic study "p3-abstention" (an_p3_helpers.jl), against the independent reference of an_p3_ref.jl and
# the designed study (hand table in the test plan, "Synthetic study"). Spec: analysis-spec.md, Amendment D93, E49 to E58, I30 to I37.
# The result of `calibrate` is a NamedTuple (the contract of the plan): status, outcome, record, tables, rows, dropped, sets, relabel, verdicts, selection, star, matched, c7, optimism, hb, hc, rab.

const P3_RES = Dict{String, Any}()
p3_res(tag, st; kw...) = get!(() -> p3_run(st; kw...), P3_RES, tag)
const P3_REFROWS = Dict{String, Any}()
p3_refrows(tag, st) = get!(() -> p3_ref_rows(st), P3_REFROWS, tag)

# HB analogue on Cal for theta* (E54 C7): for A in T_s (I_6) and W2 (I_7): b_X, cost, flags, FR_X, FR_rule, n10, n01, the exact McNemar tail, the margin test, the coverage of the rule
function p3ref_c7(rows, q, th)
    k = p3ref_kP(q)
    cal = [o for o in rows if o.in_cal]
    m = p3ref_matched(q, th.V)
    out = Dict{Symbol, Vector{NamedTuple}}()
    for (A, flag, comps) in ((:T_s, o -> o.Ts, Symbol.(q.voter_bits)), (:W2, o -> o.W2, Symbol.(q.indicators_i7)))
        set = [o for o in cal if flag(o)]
        n = length(set)
        rule(o) = p3ref_out(o.V, k, th)
        rfr(o) = rule(o) === :regular_supported
        lst = NamedTuple[]
        for X in comps
            info = m.X[X]
            bX = X === :LE ? m.le_read : info.b
            kx = findfirst(==(bX), q.budgets) - 1
            xfr(o) = o.V[X][kx + 1] === :regular
            frx, frr = count(xfr, set) // n, count(rfr, set) // n
            n10, n01 = count(o -> rfr(o) && !xfr(o), set), count(o -> xfr(o) && !rfr(o), set)
            push!(lst, (X = X, b = bX, cost = bX * info.cost, capped = info.capped, cost_unmatched = X === :LE ? m.le_unmatched : false, fr_x = frx, fr_rule = frr, n10 = n10, n01 = n01,
                        p = p3ref_mcnemar(n10, n01), margin = frx != 0 && frr <= p3ref_rat(q.rho_hb) * frx, coverage = count(o -> rule(o) !== :unresolved, set) // n))
        end
        out[A] = lst
    end
    return out
end

@testset "D93 default scenario: status, the two output files, dispatch result" begin
    st = p3_cached("default")
    res = p3_res("default", st; threads = 2)
    @test res.status === :ok && res.outcome === :selected && isfile(res.record) && isfile(res.tables) && basename(res.record) == "p3-abstention-calibration-record.json" &&
          basename(res.tables) == "p3-abstention-calibration-tables.md"
    @test an_files(st.out) == ["p3-abstention-calibration-record.json", "p3-abstention-calibration-tables.md"]          # the output folder holds only the two files (D93p)
    io = IOBuffer()
    c = p3_copy(st, "pipe-cli")
    @test AN.cli(p3_cli_args(c); io = io) == 0
    @test occursin("p3-abstention-calibration-record.json", String(take!(io))) && an_files(c.out) == ["p3-abstention-calibration-record.json", "p3-abstention-calibration-tables.md"]
end

@testset "D93d relabelling C0 on the stores: rho_new of the panel from the P1 items and of the grid from the p2 items, the diagnostics condition, P45 (P1 L1 item), rotation numbers, change table" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    rl = Dict(r.id => r for r in res.relabel.rows)
    @test Set(keys(rl)) == Set(vcat(1:16, 21:44)) && all(id -> rl[id].rho === P3_RHO_NEW[id], keys(rl))
    @test all(id -> rl[id].source === (id <= 20 ? :panel : :grid) && rl[id].kindex == (id <= 20 ? Int(AN_ORBITS[id]["k_index"]) : st.s2.h4[id - 20].kindex), keys(rl))
    # orbit 13: the L5 item of the P1 store carries "diagnostics exceeded": both axes unresolved. Orbit 14: the P1 L1 item differs from the L5 item by 1e-6 (P45: the P1 L1 item, not the p2 one)
    @test rl[13].instability == AN.UNRESOLVED && rl[13].structure == AN.UNRESOLVED && rl[13].lam5 == 0.30 + 0.01 * 13
    @test rl[14].rho === :U && rl[14].lam1 - rl[14].lam5 > 9.9e-7 && rl[14].instability == AN.UNRESOLVED && rl[3].lam1 == rl[3].lam5
    # the reduced rule has no ladder test, no direction test, no periodic class: orbit 2 (a P1 ladder defect) is chaotic here, regular orbits of the T class are quasiperiodic
    @test rl[2].rho === :C && all(id -> rl[id].rho === :R && rl[id].structure == "quasiperiodic" && rl[id].instability == AN.NONE_DETECTED, (9, 10, 15, 16, 27, 28, 29, 30, 31, 41, 42, 43))
    # grid orbits: lambda_L5(N_1) and lambda_L1(N_1) of the p2 items; the lifted coordinate from the wrapped value and the winding count gives the rotation numbers of the lifted orbit (1e-12)
    d = st.design
    @test all(id -> rl[id].lam5 == d.lref[id] && rl[id].lam1 == d.lref[id], filter(i -> d.rho[i] === :C, 21:43)) && rl[44].lam5 - rl[44].lam1 > 9.9e-7 && rl[44].rho === :U
    @test all(id -> abs(rl[id].rho1 - parse(Float64, st.s2.h4[id - 20].p0) / (2pi)) < 1e-12 && abs(rl[id].rho2 - parse(Float64, st.s2.h4[id - 20].p0) / (2pi)) < 1e-12, setdiff(21:44, [32])) &&
          rl[32].rho === :U && abs(rl[32].rho1 - rl[32].rho2) > 1e-7
    @test all(id -> abs(rl[id].rho1 - Float64(AN_ORBITS[id]["p0"]) / (2pi)) < 1e-12 && abs(rl[id].rho2 - Float64(AN_ORBITS[id]["p0"]) / (2pi)) < 1e-12, setdiff(1:16, [13]))
    # E49 f: the change table of the panel against the labels of the base pipeline on the P1 inputs (read only), the grid against the D89 labels: no change
    p1 = st.p1
    base = AN.analyze(; prereg = p1.prereg, manifest = p1.manifest, panel = p1.panel, params = p1.params, store = p1.store, out = joinpath(st.dir, "p1-out-d93"))
    old(id) = (l = base.labels[id]; l.instability == AN.POSITIVE && l.structure == "irregular" ? :C : l.instability == AN.NONE_DETECTED && l.structure == "quasiperiodic" ? :R : :U)
    exp_panel = Dict{Tuple{Int, Symbol, Symbol}, Int}()
    for id in 1:16
        k = (Int(AN_ORBITS[id]["k_index"]), old(id), P3_RHO_NEW[id]); exp_panel[k] = get(exp_panel, k, 0) + 1
    end
    exp_grid = Dict{Tuple{Int, Symbol, Symbol}, Int}()
    for id in 21:44
        k = (st.s2.h4[id - 20].kindex, P3_RHO_NEW[id], P3_RHO_NEW[id]); exp_grid[k] = get(exp_grid, k, 0) + 1
    end
    @test res.relabel.changes.panel == exp_panel && res.relabel.changes.grid == exp_grid && res.relabel.changes.grid_changed == 0 &&
          res.relabel.changes.panel_changed == count(id -> old(id) !== P3_RHO_NEW[id], 1:16) && res.relabel.changes.panel_changed >= 1 && sum(values(res.relabel.changes.panel)) == 16
    @test old(2) === :U && P3_RHO_NEW[2] === :C && haskey(exp_panel, (1, :U, :C)) && isempty(res.relabel.missing)
    # no orbit of the p2 grid is relabelled by any other rule: the grid labels equal the labels of the p2 pipeline of D89 on the same items
    @test all(id -> rl[id].rho === P3_RHO_NEW[id], 21:44)
end

@testset "D93e, D93f Cal, the sets and the verdicts on the stores equal the independent reference; I35, I33, I37" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    rows = p3_refrows("default", st)
    got = Dict(r.id => r for r in res.rows)
    @test Set(keys(got)) == Set(vcat(1:16, 21:44)) && all(r -> (got[r.id].source, got[r.id].kindex, got[r.id].in_cal, got[r.id].rho, got[r.id].S2, got[r.id].Ts, got[r.id].W2, got[r.id].gate) ==
                                                              (r.source, r.kindex, r.in_cal, r.rho, r.S2, r.Ts, r.W2, r.gate), rows)
    # the sets by hand (plan, "Synthetic study"): S2, T_s, W2, R_W, regular, gate, Cal; the dropped orbits with their reasons; the conservation of counts
    S = res.sets
    @test S.S2 == [11, 12, 33, 34, 35, 36, 37, 38] && S.T_s == [5, 6, 23, 24, 25, 39] && S.W2 == [1, 3, 5, 6, 7, 8, 21, 22, 23, 24, 25, 26, 39] && S.R_W == [1, 3, 7, 8, 21, 22, 26] &&
          S.regular == [9, 10, 15, 16, 27, 28, 29, 30, 31, 41, 42, 43] && S.gate == [9, 10] && isempty(S.overlap) && length(S.cal) == 36
    @test [(d.id, d.reason) for d in S.dropped] == [(13, "unresolved"), (14, "unresolved"), (17, "control"), (18, "control"), (19, "control"), (20, "control"), (32, "unresolved"), (44, "unresolved")] &&
          [(d.id, d.reason) for d in res.dropped] == [(d.id, d.reason) for d in S.dropped]
    @test sort(vcat(S.cal, [d.id for d in S.dropped])) == vcat(1:20, 21:44) && issubset(S.T_s, S.W2) && S.R_W == setdiff(S.W2, S.T_s) && all(i -> got[i].rho === :C, vcat(S.S2, S.W2))
    @test all(d -> d.source === (d.id <= 20 ? :panel : :grid), res.dropped)
    # the verdicts: those of D89 for the same orbit and item (the designed sequences read back by the independent reference), at the ten budgets, for every orbit of Cal and of the gate set; I_7 only
    ref = p2_ref_all(st.s2)
    need = sort(union(S.cal, S.gate))
    @test Set(k[1] for k in keys(res.verdicts)) == Set(need) && all(k -> k[2] in P2_X7, keys(res.verdicts)) && length(res.verdicts) == 7 * length(need)
    bad = [(id, X) for id in need, X in P2_X7 if res.verdicts[(id, X)] != ref.V1[(id, X)]]
    @test isempty(bad) && all(r -> r.V == (r.id in need ? Dict(X => ref.V1[(r.id, X)] for X in P2_X7) : nothing), res.rows)
    # a change of a p3 key never changes a verdict (the verdict rules and the six numbers of R1 come from the p2 table)
    c = p3_copy(st, "pipe-verdict-key")
    q2 = merge(c.q, (; e_cal = 0.25, cost_tg = 2700.0, delta_tol = 0.06))
    t2, _ = p3_prereg_text(q2); write(c.prereg, t2); p3_write_params(c.params, q2)
    r2 = p3_run(c)
    @test r2.status === :ok && r2.verdicts == res.verdicts && r2.sets.cal == S.cal
end

@testset "D93j C3 to C8 on the stores: theta*, F and the counts of each item, J_min, Theta_tol, the C7 items equal the independent reference exactly" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    rows = p3_refrows("default", st)
    ref = p3ref_select(rows, st.q)
    sel = res.selection
    @test same_sel(sel, ref) && res.outcome === :selected && sel.star == 147 && res.star == ref.thetas[147] && length(sel.feasible) == 62 && length(sel.tol) == 42 && sel.alone == (690, 1745, 1044, 1840, 681, 554)
    # the matched budgets and flags of theta* (E53) and the Double64 costs, in the result
    m = p3ref_matched(st.q, res.star.V)
    @test res.matched.CR == m.CR && all(X -> (res.matched.X[X].b, res.matched.X[X].capped, res.matched.X[X].round_down, res.matched.X[X].round_down_none) == (m.X[X].b, m.X[X].capped, m.X[X].round_down, m.X[X].round_down_none), keys(m.X)) &&
          res.matched.le_read == m.le_read && res.matched.le_unmatched && res.matched.X[:LE].capped
    # C7: the first 20 elements of Theta_tol in the order of C6, the HB analogue on Cal (exact McNemar tail, margin test with FR_X = 0 failing, coverage)
    @test res.c7.top == p3ref_order(ref.thetas, ref.stats, ref.tol)[1:20] && length(res.c7.top) == 20 && res.c7.top[1] == 147
    hb = p3ref_c7(rows, st.q, ref.thetas[147])
    for A in (:T_s, :W2)
        @test length(res.c7.hb[A]) == length(hb[A]) && all(i -> res.c7.hb[A][i] == hb[A][i], eachindex(hb[A]))
    end
    @test length(hb[:T_s]) == 6 && length(hb[:W2]) == 7 && hb[:W2][1].X === :LE && hb[:W2][1].cost_unmatched && hb[:W2][1].b == st.q.N_c && all(r -> !r.cost_unmatched, hb[:T_s])
    @test all(r -> r.margin == (r.fr_x != 0 && r.fr_rule <= 1 // 2 * r.fr_x), vcat(hb[:T_s], hb[:W2])) && any(r -> r.fr_x == 0 && !r.margin, vcat(hb[:T_s], hb[:W2])) && all(r -> r.n10 + r.n01 > 0 || r.p == 1, vcat(hb[:T_s], hb[:W2]))
    # the fitted values of theta*: e, the unresolved shares, the gate count, the rule's FR and decided shares (also by source), and J
    s = ref.stats[147]
    @test res.c7.e == s.e && res.c7.u_cal == s.u_cal && res.c7.u_s2 == s.u_s2 && res.c7.g == s.g && res.c7.J == ref.J[147] && res.c7.cost == s.cost
    # the risk-coverage points: the decided share of Cal against e for the rule family around theta* (q_C and q_R each varied by one step inside their ranges, V, w, v fixed; theta* itself
    # included) and for each single indicator of I_7 at its matched budget (its chaotic verdict as chaotic-supported, its regular verdict as regular-supported)
    th = ref.thetas[147]
    fam = [i for (i, t) in enumerate(ref.thetas) if t.V == th.V && t.w == th.w && t.v == th.v && ((t.qC == th.qC && abs(t.qR - th.qR) <= 1) || (t.qR == th.qR && abs(t.qC - th.qC) <= 1))]
    @test length(fam) == 3 && Set((p.qC, p.qR) for p in res.c7.risk_coverage.family) == Set((ref.thetas[i].qC, ref.thetas[i].qR) for i in fam) && length(res.c7.risk_coverage.family) == 3 &&
          all(p -> (i = first(j for j in fam if (ref.thetas[j].qC, ref.thetas[j].qR) == (p.qC, p.qR)); p.coverage == ref.stats[i].decided // 36 && p.e == ref.stats[i].e), res.c7.risk_coverage.family)
    cal = [o for o in rows if o.in_cal]
    for p in res.c7.risk_coverage.single
        bX = p.X === :LE ? m.le_read : m.X[p.X].b
        kx = findfirst(==(bX), st.q.budgets) - 1
        dec = [o for o in cal if o.V[p.X][kx + 1] !== :undecided]
        contra = count(o -> (o.V[p.X][kx + 1] === :chaotic && o.rho === :R) || (o.V[p.X][kx + 1] === :regular && o.rho === :C), dec)
        @test p.b == bX && p.coverage == length(dec) // 36 && p.e == (isempty(dec) ? nothing : contra // length(dec))
    end
    @test Set(p.X for p in res.c7.risk_coverage.single) == Set(Symbol.(st.q.indicators_i7)) && length(res.c7.risk_coverage.single) == 7
end

@testset "D93j outcomes recorded with exit code 0 and no PP: 'no feasible rule' (F empty) and 'sets too small'; the exact fractions in the record" begin
    st = p3_cached("default")
    c = p3_copy(st, "pipe-empty")
    q2 = p3_q_empty(c.q)
    t2, _ = p3_prereg_text(q2); write(c.prereg, t2); p3_write_params(c.params, q2)
    r = p3_run(c)
    @test r.status === :ok && r.outcome === :no_feasible_rule && r.star === nothing && isempty(r.selection.feasible) && r.hb === nothing && r.hc === nothing && r.optimism === nothing &&
          same_sel(r.selection, p3ref_select(p3_refrows("default", st), q2))
    rec = AX.read_json(r.record)
    @test rec["outcome"]["selection"] == "no feasible rule" && rec["selection"]["star"] === nothing && rec["pp_hb"] === nothing && rec["pp_hc"] === nothing && rec["theta_space"]["size"] == 2088
    io = IOBuffer(); @test AN.cli(p3_cli_args(c; extra = ["--out", c.out * "-cli"]); io = io) == 0
    # sets too small: a second designed study with three sticky orbits in Cal (min_set = 5)
    sm = p3_cached("small"; design = p3_design_small())
    rs = p3_run(sm)
    @test rs.status === :ok && rs.outcome === :sets_too_small && rs.star === nothing && rs.hb === nothing && rs.hc === nothing && rs.sets.T_s == [5, 6, 23] && length(rs.sets.T_s) < sm.q.min_set
    @test AX.read_json(rs.record)["outcome"]["selection"] == "sets too small" && AX.read_json(rs.record)["pp_hb"] === nothing
    io = IOBuffer(); @test AN.cli(p3_cli_args(sm; extra = ["--out", sm.out * "-cli"]); io = io) == 0
    # the rationals of the record are texts a/b in lowest terms; the counts are integers
    rec = AX.read_json(p3_res("default", st).record)
    @test rec["selection"]["J_min"] isa String && occursin(r"^-?\d+/\d+$", rec["selection"]["J_min"]) && (a = split(rec["selection"]["J_min"], '/'); gcd(parse(BigInt, a[1]), parse(BigInt, a[2])) == 1) &&
          rec["theta_space"]["size"] isa Integer && rec["selection"]["feasible"] isa Integer && rec["selection"]["tol"] isa Integer && rec["selection"]["J_min"] == "19/156"
end

@testset "D93k, D93l the PP chain on the stores: optimism, rates, path, targets, type-I, AUC_I_cal and the n_f table" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    q = st.q
    rows = p3_refrows("default", st)
    ref = p3ref_select(rows, q)
    # optimism by the reference (the streams (optimism, b), the resamples of Cal by increasing id, C3 to C6 on each)
    opt = p3ref_optimism(rows, q)
    @test res.optimism.used == opt.used && res.optimism.skipped == opt.skipped && res.optimism.thetas == opt.thetas && opt.used + opt.skipped == q.b_opt
    @test isapprox(res.optimism.optimism.T_s, opt.optimism.T_s; rtol = 1e-12, atol = 1e-15) && isapprox(res.optimism.optimism.W2, opt.optimism.W2; rtol = 1e-12, atol = 1e-15)
    # the rates of theta*: r from Cal exactly, the corrected rate r + max(optimism, 0), the decided shares
    h = res.hb.rates
    s = ref.stats[147]
    @test h.r[:T_s] == Float64(s.fr_ts) && h.r[:W2] == Float64(s.fr_w2) && h.r_corr[:T_s] == h.r[:T_s] + max(res.optimism.optimism.T_s, 0) && h.r_corr[:W2] == h.r[:W2] + max(res.optimism.optimism.W2, 0)
    @test isapprox(h.d[:T_s], Float64(s.dec_ts); atol = 1e-15) && all(x -> 0 <= x <= 1, vcat(collect(values(h.p_rule)), collect(values(h.p_comp)), collect(values(h.s_comp)), collect(values(h.d))))
    @test all(p -> h.p_rule[p] >= q.rule_rate_floor, (:T_s, :R_W)) && Set(k[1] for k in keys(h.p_comp)) == Set(Symbol.(q.indicators_i7)) && Set(k[2] for k in keys(h.p_comp)) == Set([:T_s, :R_W])
    # the path: every allowed point simulated at pp_runs runs with the streams (pp-hb, (j, i)); the counts equal the reference simulation of the model of the rates
    md = res.hb.model
    @test md.rule[:T].p == h.p_rule[:T_s] && md.rule[:RW].p == h.p_rule[:R_W] && md.ts == Symbol.(q.voter_bits) && md.w2 == Symbol.(q.indicators_i7)
    refmodel = (rule = md.rule, comp = Dict(X => get(md.comp, X, (T = (p = 0.0, s = 0.0), RW = (p = 0.0, s = 0.0))) for X in P2_X7), ts = md.ts, w2 = md.w2)
    pa = res.hb.path
    @test length(pa) == 17 && [x.allowed for x in pa] == [x.allowed for x in p3ref_path(q)] &&
          all(x -> x.successes == p3ref_hb_power(q.test_seed, x.j, refmodel, x.nT, x.nRW, q, q.pp_runs), [x for x in pa if x.allowed])
    @test res.hb.targets == AN.p3_hb_targets(pa, q) && res.hb.targets.outcome in (:ok, :cap)
    # the type-I check at rho* (or at r_max when the cap binds): for each X the same simulation with the rule rate set to the rate of X; the oracle: the joint rate at most alpha_hb + 3 se
    t = res.hb.targets
    pj = t.outcome === :ok ? t.j : 16
    pt = pa[pj + 1]
    @test all(X -> res.hb.type1.rates[X].joint == p3ref_hb_power(q.test_seed, pj, merge(refmodel, (rule = Dict(:T => (p = refmodel.comp[X][:T].p, d = md.rule[:T].d), :RW => (p = refmodel.comp[X][:RW].p, d = md.rule[:RW].d)),)),
                                                            pt.nT, pt.nRW, q, q.pp_runs) / q.pp_runs, P2_X7)
    @test res.hb.type1.max_joint <= q.alpha_hb + 3 * sqrt(q.alpha_hb * (1 - q.alpha_hb) / q.pp_runs)
    # HC: Cal orbits with rho in {C, R} and K > 0 (the gate and all K = 0 orbits are out); the plain label of the decided Float64 votes of V* at b_P; a tie is a failure; s_I; the AUC
    kP = p3ref_kP(q)
    pop = [o for o in rows if o.in_cal && !(o.source === :panel && o.kindex == 1)]
    votes(o) = Symbol[o.V[X][kP + 1] for X in ref.thetas[147].V if o.V[X][kP + 1] !== :undecided]
    fails = Float64[]; nonf = Float64[]
    for o in pop
        v = votes(o); nC, nR = count(==(:chaotic), v), count(==(:regular), v)
        f = nC == nR ? true : (nC > nR ? :C : :R) !== o.rho                      # a tie is a failure
        score = nC + nR == 0 ? 0.0 : 1 - abs(nC - nR) / (nC + nR)
        push!(f ? fails : nonf, score)
    end
    @test res.hc.n_fail == length(fails) && res.hc.n_nonfail == length(nonf) && length(fails) + length(nonf) == length(pop) && length(pop) == 26 && length(fails) == 2          # failures by design: orbit 29 (two chaotic votes of three on a regular orbit) and orbit 42 (a tie)
    @test res.hc.auc_i_cal == Float64(p3ref_auc(fails, nonf)) && isapprox(res.hc.auc_i_used, 0.5 + q.auc_shrink * (res.hc.auc_i_cal - 0.5); atol = 1e-15)
    hc = AN.p3_pp_hc(q, res.hc.auc_i_cal; runs = q.pp_runs, threads = 1)
    @test [x.successes for x in res.hc.table] == [x.successes for x in hc.table] && res.hc.target == hc.target && res.hc.underpowered == hc.underpowered && res.hc.deliverable_expected == 49 && res.hc.deliverable_worst == 62
end

@testset "D93j R_ab on the Cal orbits that carry Double64 items (P63): reported, equal to the reference" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    ref = p2_ref_all(st.s2)
    q = st.q
    th = res.star
    kP = p3ref_kP(q)
    # the Double64 verdicts of the orbits of S_L (1, 5, 9, 11) at b_P by the independent reference (table R1 on the Double64 items: horizon N_dir for L1 and tg, N_c for the others)
    items = ref.items
    function d64(id)
        V = Dict{Symbol, Vector{Symbol}}()
        for X in th.V
            vals = [p2ref_value(items, st.s2.p, id, X, k; prec = "Double64", horizon = X in (:SALI, :FLI, :MEGNO, :LE) ? P2_NDIR : nothing) for k in 0:9]
            V[X] = [p2ref_verdict(X, :R1, v, st.s2.p) for v in vals]
        end
        return V
    end
    exp = Dict(id => p3ref_out_ab(Dict(X => ref.V1[(id, X)] for X in P2_X7), d64(id), kP, th) for id in (1, 5, 9, 11))
    @test Dict(r.id => r.output for r in res.rab) == exp && Set(r.id for r in res.rab) == Set([1, 5, 9, 11])
end

@testset "D93n output: the record and the tables file hold the items of E58; the digest of the F1 data block; the record hash; no clock, no absolute path" begin
    st = p3_cached("default"); res = p3_res("default", st; threads = 2)
    text = read(res.record, String)
    rec = AX.read_json(res.record)
    @test text == replace(text, "\r\n" => "\n") && endswith(text, "\n") && !occursin(st.dir, text) && !occursin(replace(st.dir, '\\' => '/'), text) && !occursin(r"\d{4}-\d{2}-\d{2}T\d{2}:", text)
    top = ["program", "fingerprints", "parameters", "outcome", "relabelling", "sets", "theta_space", "selection", "c7", "optimism", "pp_hb", "pp_hc", "r_ab", "readings"]
    @test all(k -> haskey(rec, k), top) && p3_top_keys(text) == top                                                                  # fixed order of the keys
    @test Set(keys(rec["parameters"])) == Set(string.(P3_PARAM_KEYS)) && length(rec["parameters"]) == 71 && rec["parameters"]["test_seed"] == st.q.test_seed && rec["parameters"]["pp_runs"] == st.q.pp_runs
    # the fingerprints (C1): SHA-256 of the preregistration under the rule, of the p2 manifest, the P1 manifest, the panel file, the three parameter tables; the commit and the clean flag
    fp = rec["fingerprints"]
    @test fp["p2_manifest"] == st.q.p2_manifest_sha256 == an_sha256(read(st.s2.manifest)) && fp["p1_manifest"] == an_sha256(read(st.p1.manifest)) && fp["panel"] == an_sha256(read(st.s2.panel)) &&
          fp["p3_parameters"] == an_sha256(read(st.params)) && fp["p2_parameters"] == an_sha256(read(st.s2.params)) && fp["p1_parameters"] == an_sha256(read(st.p1.params)) && fp["preregistration"] == st.prereg_hash &&
          haskey(fp, "git_commit") && fp["git_clean"] isa Bool && haskey(rec["program"], "generator") && occursin("xoshiro256++", rec["program"]["generator"]) && occursin("AS241", rec["program"]["quantile"])
    # the tables file: the record hash (a file cannot state its own hash), T1, the F1 data block with its digest recomputed by an independent function
    tab = read(res.tables, String)
    m = match(r"record SHA-256:? ([0-9a-f]{64})", tab)
    @test m !== nothing && m.captures[1] == an_sha256(read(res.record)) && !occursin(st.dir, tab) && !occursin(r"\d{4}-\d{2}-\d{2}T\d{2}:", tab) && occursin("Table T1", tab) && occursin("## F1 data", tab)
    f1 = res.f1
    ref = p3ref_select(p3_refrows("default", st), st.q)
    @test [r[1] for r in f1.rows] == ref.feasible && length(f1.rows) == 62 && f1.columns[1] == "index" && f1.columns[end - 1] == "in Theta_tol" && f1.columns[end] == "theta*" &&
          [r[end - 1] for r in f1.rows] == [i in ref.tol for i in ref.feasible] && [r[end] for r in f1.rows] == [i == ref.star for i in ref.feasible] && count(r -> r[end], f1.rows) == 1
    @test p2_digest(f1.columns, f1.rows) == f1.digest && occursin("digest: " * f1.digest, tab[findfirst("## F1 data", tab).start:end])
    # the counts: Cal plus the dropped orbits are P union G, each once by reason; the by-source rows sum to the pooled rows (record)
    rl = rec["relabelling"]
    @test rl["dropped"]["total"] == 8 && sum(values(rl["dropped"]["by_reason"])) == 8 && rec["sets"]["cal"]["pooled"] == 36 && rec["sets"]["cal"]["panel"] + rec["sets"]["cal"]["grid"] == 36 &&
          rec["sets"]["T_s"]["pooled"] == 6 && rec["sets"]["T_s"]["panel"] + rec["sets"]["T_s"]["grid"] == 6 && rec["sets"]["W2"]["pooled"] == 13 && rec["sets"]["S2"]["pooled"] == 8
    @test rec["theta_space"]["size"] == 2088 && rec["theta_space"]["alone"] == collect(ref.alone) && rec["selection"]["feasible"] == 62 && rec["selection"]["tol"] == 42 && rec["selection"]["index"] == 147 &&
          rec["selection"]["theta"]["V"] == ["SALI", "FLI", "E2"] && (rec["selection"]["theta"]["qC"], rec["selection"]["theta"]["qR"], rec["selection"]["theta"]["w"], rec["selection"]["theta"]["v"]) == (1, 2, 1, 0)
end

@testset "D93o determinism (I31): two runs, 1 and 4 and 2 threads, different output folders give byte-identical files; D93p read-only (I30): inputs and data tree unchanged" begin
    st = p3_cached("default")
    c = p3_copy(st, "pipe-det")
    pre = p3_fingerprint(c)
    h0 = p2_datatree_hash()
    r1 = p3_run(c; out = joinpath(c.dir, "o1"), threads = 1)
    r4 = p3_run(c; out = joinpath(c.dir, "o4"), threads = 4)
    r2 = p3_run(c; out = joinpath(c.dir, "sub", "o2"), threads = 2)
    r2b = p3_run(c; out = joinpath(c.dir, "o2b"), threads = 2)
    for r in (r4, r2, r2b)
        @test read(r1.record) == read(r.record) && read(r1.tables) == read(r.tables)
    end
    @test p3_fingerprint(c) == pre && p2_datatree_hash() == h0 && an_files(joinpath(c.dir, "o1")) == ["p3-abstention-calibration-record.json", "p3-abstention-calibration-tables.md"]
    # the command line with --threads, the same bytes
    io = IOBuffer()
    @test AN.cli(p3_cli_args(c; extra = ["--threads", "3", "--out", joinpath(c.dir, "ocli")]); io = io) == 0
    @test read(joinpath(c.dir, "ocli", "p3-abstention-calibration-record.json")) == read(r1.record)
    # a refused run leaves every input unchanged and creates no folder; regeneration after deleting the outputs gives the same bytes
    bad = p3_copy(st, "pipe-det-bad"); p3_edit_manifest(bad, m -> m["status"] = "running"); p3_resync!(bad)
    pb = p3_fingerprint(bad)
    @test an_rejects(() -> p3_run(bad), ["status"]) && !ispath(bad.out) && p3_fingerprint(bad) == pb
    rm(joinpath(c.dir, "o1"); recursive = true)
    r1b = p3_run(c; out = joinpath(c.dir, "o1"), threads = 1)
    @test read(r1b.record) == read(r4.record)
end

@testset "D93 P62: a preregistration edited with its hash cell updated (Amendment A1) still passes; the recomputed hash is recorded" begin
    st = p3_cached("default")
    c = p3_copy(st, "pipe-a1")
    h = p3_edit_prereg(c, t -> t * "\n## Amendment A1 (fixture)\n\ntheta* = fixture.\n")
    @test h != c.prereg_hash
    r = p3_run(c)
    @test r.status === :ok && AX.read_json(r.record)["fingerprints"]["preregistration"] == h
    # the result of the run on the edited preregistration is the same selection (the preregistration is a check, not an input of a number)
    @test same_sel(r.selection, p3_res("default", st; threads = 2).selection)
end
