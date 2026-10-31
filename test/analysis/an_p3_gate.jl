# D93c, P68 (stage 03 rework of the stage 02 tests, owner decision 3 of the stage 02 checkpoint): a gate orbit that is not in Cal and lacks a verdict item is refused (E48 h). The gate
# orbits of the designed study (9 and 10) are in Cal, so the case needs a fixture the designed study does not have: the L5 item that C0 reads for orbit 9 carries "diagnostics
# exceeded", which makes its reference class unresolved (E49 c): it leaves Cal ("unresolved") but stays in the gate set (E50, not intersected with Cal).

@testset "D93c P68 a gate orbit outside Cal (rho = U) with a missing verdict item is refused, naming the item and the orbit; with the item the same orbit is accepted" begin
    st = p3_copy(p3_cached("default"), "gate-outside-cal")
    p2_set_attr(p3_p1chunk(st, 5), "orbits/9", "flags", ["diagnostics exceeded"])
    ok = p3_run(st; out = joinpath(st.dir, "out-accepted"), threads = 1)
    @test ok.status === :ok && (9, "unresolved") in [(d.id, d.reason) for d in ok.dropped] && 9 in ok.sets.gate && !(9 in ok.sets.cal) && 9 in Set(k[1] for k in keys(ok.verdicts))
    p3_rewrite(st, 14, o -> o.id == 9 ? nothing : o)                                 # the pa1 item of orbit 9 is gone (chunk 14 holds the pa1 items)
    p3_edit_manifest(st, m -> m["plan_counts"]["items"] -= 1); p3_resync!(st)
    pre = p3_fingerprint(st)
    @test an_rejects(() -> p3_run(st; out = joinpath(st.dir, "out-refused")), ["pa1", "9"]) && !ispath(joinpath(st.dir, "out-refused")) && p3_fingerprint(st) == pre
    # the same item missing for an orbit that is neither in Cal nor in the gate set is not a refusal (orbit 13 is unresolved and not a gate orbit)
    st2 = p3_copy(p3_cached("default"), "gate-outside-cal-other")
    p3_rewrite(st2, 14, o -> o.id == 13 ? nothing : o)
    p3_edit_manifest(st2, m -> m["plan_counts"]["items"] -= 1); p3_resync!(st2)
    r2 = p3_run(st2; out = joinpath(st2.dir, "out-accepted"), threads = 1)
    @test r2.status === :ok && !(13 in r2.sets.cal) && !(13 in r2.sets.gate)
end
