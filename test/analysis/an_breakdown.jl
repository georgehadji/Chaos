# Regression found in stage 04 on the real store: an item of the tau sweep with the flag breakdown_tau (Amendment A5 item 2) stores t = 0 and NaN
# from the first row that was not computed. It is a result of the sweep, not a failure: the run must not be refused, and F7 counts it.

@testset "breakdown_tau items are results, not non-finite refusals (A5 item 2)" begin
    lamh(id, sp, steps, vals) = id == 5 && sp.method == "L4" && sp.tau == 2 ? vcat(vals[1:6], fill(NaN, length(vals) - 6)) : vals
    flagh(id, sp) = id == 5 && sp.method == "L4" && sp.tau == 2 ? ["breakdown_tau"] : String[]
    attrh(id, sp) = id == 5 && sp.method == "L4" && sp.tau == 2 ? Dict("breakdown_step" => 7, "breakdown_kind" => "orthonormalisation") : Dict()
    st = an_cached("breakdown"; lam = lamh, flags = flagh, attrs = attrh)
    # the fixture has the blank rows the runner writes
    i = findfirst(s -> s.method == "L4" && s.tau == 2, st.specs)
    o = AX.read_chunk(an_chunk(st, i)).orbits
    o5 = first(x for x in o if x.id == 5)
    @test o5.t[7:end] == zeros(Int, length(o5.t) - 6) && !any(isfinite, o5.lyap["L4"][7:end, :]) && "breakdown_tau" in o5.attrs["flags"]
    res = an_res("breakdown", st)
    @test res.status === :ok
    row = first(r for r in res.figures["F7"].rows if r[1] == "tau" && r[2] == "2" && r[3] == "L4")
    @test row[4] == 11 && row[7] == 1                      # eleven of the twelve S_D orbits give values, one breakdown is counted
    text = read(res.tables, String)
    @test occursin("orthonormalisation breakdown", text) && occursin("| 5 | L4 | tau | 2 | 7 | orthonormalisation |", text)      # step and kind are tabulated
    # the other items of the same orbit and every decision are those of the default study
    base = an_res("default", an_cached("default"))
    @test all(an_dec(res, k) == an_dec(base, k) for k in (:H1, :H2, :H3, :H4, :H5))
    # a NaN in an item without the flag is still refused
    s2 = an_copy(st, "breakdown-nan")
    j = an_find(s2; method = "L5")
    an_set_value(an_chunk(s2, j), "orbits/1/lyap/L5", (3, 1), NaN)
    @test an_rejects(() -> an_run(s2), ["non-finite"])
end
