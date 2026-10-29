# D89 variant scenario: a flagged (drift) orbit, a fifth sticky orbit, a violating gate for E1, an extension that differs. Spec: analysis-spec.md, Amendment D89,
# E29, E33, criterion D89n (ROBUST, NOT ROBUST, NOT COMPARABLE, the second versions), D89k (flagged orbits dropped and listed).
# Differences from the default scenario (an_p2_study.jl, variant :var): orbit 8 has |p| = 2e6 at every recorded step and the flag drift on all its items; orbit 25 is sticky
# (T_s has 5 orbits = min_set); gate orbit 10 has E1 chaotic at k = 9 (share 1 of 2 > 0.10); RLI of the W2 orbits 1, 21, 22, 26 turns chaotic only at k = 8.

@testset "D89 variant: drift orbit dropped and listed, T_s with five orbits, decisions against the reference" begin
    st = p2_cached("var"; variant = :var)
    res = p2_res("var", st)
    (; ref, S) = p2_refs("var", st)
    @test res.status === :ok && 8 ∉ res.sets.E_all && (8, "drift") in [(x.id, x.reason) for x in res.sets.dropped] && 8 ∉ res.sets.W2
    @test res.sets.T_s == [5, 23, 24, 25, 39] && res.sets.W2 == [1, 3, 5, 7, 21, 22, 23, 24, 25, 26, 39] && p2_diff(res.sets, S) == Symbol[]
    # I16: E_all plus the dropped orbits are P union G; the orbit is counted in T1 once and appears as flagged
    @test sort(vcat(res.sets.E_all, [x.id for x in res.sets.dropped])) == collect(1:44)
    T1 = res.items.T1
    r = p2_rowsof(T1, "source" => "panel", "K index" => 1)[1]
    @test r[p2_tcol(T1, "flagged")] == 1 && sum(x[p2_tcol(T1, "orbits")] for x in T1.rows) == 44
    # the decisions, all six, equal the reference; H4 is decidable with 5 orbits and supported (Delta_phi = 1 on identical sticky orbits)
    rf = p2_ref_hyp(st, ref, S, P2_X7)
    @test p2_decision_diff(res.decisions, rf) == String[]
    @test res.decisions.H4.verdict === :supported && res.decisions.H4.n == 5 && res.decisions.H4.dphi == 1.0 && res.decisions.H4.dphi_lo == 1.0
    @test res.decisions.H5.verdict === :not_supported && res.decisions.H1.verdict === :supported && res.decisions.H4.design.verdict === :not_decidable
    # the extension: H5 on E8 with eight indicators is supported (the mean of the two middle values of eight F_X), the primary is not: NOT ROBUST; H4 ROBUST; H6 NOT COMPARABLE
    ext = p2_ref_hyp(st, ref, S, P2_X8; suffix = "_8", srule = "E8")
    @test p2_decision_diff(res.extension, ext) == String[] && res.extension.H5.verdict === :supported && res.extension.H5.F_med > 0.1
    @test res.robust == (H1 = :ROBUST, H2 = :ROBUST, H3 = :ROBUST, H4 = :ROBUST, H5 = :NOT_ROBUST, H6 = :NOT_COMPARABLE)
end

@testset "D89 variant: a violating gate adds the second versions and leaves every first version unchanged" begin
    st = p2_cached("var"; variant = :var)
    res = p2_res("var", st)
    (; ref, S) = p2_refs("var", st)
    gate = p2_ref_gate(st, ref, S)
    @test p2_same(collect(res.gate), gate)
    viol = [r.X for r in res.gate if r.violated]
    @test viol == [:E1] && (r = first(r for r in res.gate if r.X === :E1); (r.decided, r.chaotic, r.share) == (2, 1, 0.5))
    # the second version without E1: the indicator set I_7 minus E1 (six indicators) and I_8 minus E1 (seven), the other numbers unchanged, labels with the variant minus-E1
    @test collect(keys(res.second)) == [:E1]
    sv = res.second[:E1]
    I6 = [X for X in P2_X7 if X !== :E1]; I7 = [X for X in P2_X8 if X !== :E1]
    s1 = p2_ref_hyp(st, ref, S, I6; var = "minus-E1")
    s2 = p2_ref_hyp(st, ref, S, I7; var = "minus-E1", suffix = "_8", srule = "E8")
    @test p2_decision_diff(sv.decisions, s1) == String[] && p2_decision_diff(sv.extension, s2) == String[]
    @test length(sv.decisions.H1.A) == 15 && length(sv.extension.H1.A) == 21 && !haskey(sv.decisions.H2.M, :E1) && length(sv.decisions.H3.R) == 6
    # the first versions are the full ones (never replaced): the primary decisions equal the reference on the seven indicators
    @test p2_decision_diff(res.decisions, p2_ref_hyp(st, ref, S, P2_X7)) == String[] && haskey(res.decisions.H2.M, :E1) && length(res.decisions.H1.A) == 21
    # the file: the violation is printed with the indicator, the second versions are separate sections
    text = read(res.tables, String)
    @test occursin("rule invalid on the integrable control", text) && occursin("E1", first(split(text[findfirst("rule invalid on the integrable control", text).start:end], '\n')))
    for id in ("T3", "T4", "T5", "T6", "T7", "T8", "T13")
        @test occursin(Regex("(?m)^## $id[^\\n]*confirmatory"), text) && occursin(Regex("(?m)^## $id without E1"), text)
    end
    # the second version uses the H3 count rule of the table (5 of the 6 left): the count keeps h3_count_min
    @test sv.decisions.H3.passing isa Integer
end
