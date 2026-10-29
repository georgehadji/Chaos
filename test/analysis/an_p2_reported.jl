# D89 reported items without a decision rule: T11 sensitivities, T15 packages, T16 flags and controls and the K = 0 line, T17 cost (E35). Spec: analysis-spec.md,
# Amendment D89, E35, criterion D89p (content) and the Analysis plan. Expected values: the reference of an_p2_refstore.jl and hand values of the plan.

@testset "D89 T11 sensitivities on S_L*: fq1_raw, fq1_qp, pa_planar and the cos(x) row of WBA against their primary" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    want = p2_ref_sens(st, ref, S)
    @test Set(keys(res.sensitivity)) == Set(keys(want)) && length(want) == 40
    bad = [k for k in keys(want) if !p2_same(res.sensitivity[k], want[k])]
    @test isempty(bad)
    # hand values: S_L* = {1, 5, 9, 11}. fq1_raw has D + 0.5 on every row: |dD| = 0.5, the verdicts are those of fq1 (-9.x + 0.5 stays below -4, -1 + 0.5 above)
    r = res.sensitivity[(:raw_vs_fq1, 3)]
    @test r.n == 4 && r.A == 1.0 && r.kappa == 0.0 && isapprox(r.median_diff, 0.5; atol = 1e-12) && isapprox(r.max_diff, 0.5; atol = 1e-12)       # all four regular: p_e = 1, kappa = 0
    r = res.sensitivity[(:raw_vs_fq1, 5)]
    @test r.n == 4 && r.A == 1.0 && isapprox(r.kappa, 1.0; atol = 1e-12)                                                                           # orbit 11 chaotic, three regular: kappa = 1
    # fq1_qp has D - 0.25 on odd windows only; pa_planar D + 0.1 on odd windows only
    @test isapprox(res.sensitivity[(:qp_vs_fq1, 1)].max_diff, 0.25; atol = 1e-12) && res.sensitivity[(:qp_vs_fq1, 0)].max_diff == 0.0 && res.sensitivity[(:qp_vs_fq1, 0)].median_diff == 0.0
    @test isapprox(res.sensitivity[(:planar_vs_pa1, 1)].max_diff, 0.1; atol = 1e-12) && res.sensitivity[(:planar_vs_pa1, 2)].max_diff == 0.0
    # WBA from the row cos(1 x[1]) alone: orbits with id = 2 mod 3 (orbit 11 of S_L) have the minimum in the second row, so the two rules disagree on them from k = 3
    @test all(k -> res.sensitivity[(:cos_vs_wba, k)].n >= 3, 0:9) && any(k -> res.sensitivity[(:cos_vs_wba, k)].A < 1.0, 0:9)
    @test res.items.T11 isa NamedTuple && length(res.items.T11.rows) == 40
end

@testset "D89 T15 packages: FLI with RLI (on E8) and MEGNO with SALI, joint verdicts at b_P and N_c" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st)
    want = p2_ref_packages(st, ref, S)
    @test Set(keys(res.packages)) == Set(keys(want)) && length(want) == 4 * 2 * 2
    @test all(k -> p2_same(res.packages[k], want[k]), keys(want))
    # hand values: W2 at b_P (k = 6): MEGNO chaotic on W (1, 3, 7) and W2 (8, 21, 22, 26) orbits and regular on the sticky ones; SALI the same: (c, c) 7, (r, r) 4
    @test res.packages[("W2", (:MEGNO, :SALI), 6)] == Dict((c_, c_) => 7, (r_, r_) => 4)
    # FLI with RLI uses the orbits of E8 only: W2 has 8 of them
    @test sum(values(res.packages[("W2", (:FLI, :RLI), 6)])) == 8 && sum(values(res.packages[("W2", (:MEGNO, :SALI), 6)])) == 11
    @test sum(values(res.packages[("regular", (:FLI, :RLI), 9)])) == 9                                                      # the regular orbits of E8: 9, 27 to 31, 41 to 43
end

@testset "D89 T16 flags, omissions, controls and the K = 0 line; T17 cost" begin
    st = p2_cached("default"); res = p2_res("default", st); (; ref, S) = p2_refs("default", st); p = st.p
    man = AX.read_json(st.manifest)
    # flag counts per spec: the sum over the chunks of the manifest's indicator_flags, which the module recomputed from the stored words
    want = Dict{String, Dict{String, Int}}()
    for (ord, sp) in enumerate(st.specs)
        e = man["chunks"][ord]
        haskey(e, "indicator_flags") || continue
        d = get!(want, sp.id, Dict{String, Int}())
        for (k, v) in e["indicator_flags"]; d[k] = get(d, k, 0) + v; end
    end
    @test Set(keys(res.flag_counts)) >= Set(keys(want)) && all(id -> all(((f, n),) -> get(res.flag_counts[id], f, 0) == n, want[id]), keys(want))
    @test res.omissions == Dict("fq1" => [20], "bk1" => [20])
    # controls: lambda_1(N_1) against the analytic exponent with the stored deviation and scope; the recomputed deviation equals the stored one
    rows = Dict(r.id => r for r in res.controls_table)
    @test sort(collect(keys(rows))) == [17, 18, 19, 20] && rows[17].scope == "bound" && rows[18].scope == "info" && rows[19].scope === nothing && rows[20].scope == "bound"
    @test rows[17].lambda == 0.69 && rows[17].deviation == 0.0 && isapprox(rows[18].deviation, 0.3; atol = 1e-12) && isapprox(rows[18].stored_deviation, rows[18].deviation; atol = 1e-12)
    @test rows[20].lambda == 0.96 && rows[20].deviation == 0.0
    @test all(X -> rows[17].verdicts9[X] === :chaotic, P2_X7) && all(X -> rows[19].verdicts9[X] === :regular, P2_X7) && !haskey(rows[20].verdicts9, :E1) && rows[20].verdicts9[:LE] === :chaotic
    # the K = 0 line: the largest circle distance of the stored frequency nu1 to p0 / (2 pi) modulo 1 over the gate orbits, per window (hand: orbit 10 is the largest)
    kl = p2_ref_k0(st, ref, S)
    @test Set(keys(res.k0line)) == Set(p.windows) && all(W -> p2_same((fq1 = res.k0line[W].fq1, pa1 = res.k0line[W].pa1), (fq1 = kl[W].fq1, pa1 = kl[W].pa1); rtol = 1e-9), p.windows)
    @test res.k0line[20].fq1 == 0.0 && isapprox(res.k0line[10240].fq1, 2 * 9 * 1e-10; atol = 1e-14) && isapprox(res.k0line[10240].pa1, 6 * 9 * 1e-10; atol = 1e-14)
    # T17: the measured cost per orbit-step against the line, the bytes per orbit and the transform seconds, read from the manifest (nothing is measured)
    ids = sort(collect(keys(P2_SPEC_INFO)))
    @test Set(keys(res.cost)) == Set(ids)
    for (i, id) in enumerate(ids)
        c = res.cost[id]
        secs = 1.5 * count(sp -> sp.id == id && sp.kind == "frequency", st.specs)
        @test c.measured_ns == P2_COST[id][2] && c.budget_line_ns == P2_COST[id][1] && c.bytes_per_orbit == 9232 + 100 * i && c.transform_seconds == secs
    end
    @test res.cost["fq1"].transform_seconds == 4.5 && res.cost["tg"].transform_seconds == 0.0
end
