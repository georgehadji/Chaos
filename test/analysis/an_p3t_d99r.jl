# Amendment D99, rework after the stage 04 review (finding 1, owner decision 2026-11-04: fix, tests first). Spec: analysis-spec.md, Amendment D99, E83 ("item 3 holds when both rows hold";
# the status of item 3 as a whole stays in the decision line). When HA is not decidable the failed list is empty, so the status of item 3 in the decision line must come from the two rows,
# not from the failed list. Reuses d99_ha, d99_t5_fixture, d99_i3_fails of an_p3t_d99.jl and the d97_* readers. Hand values, delta_a = 0.1; no random draw.

# name => (state, statuses of the two rows of item 3, item 3 as a whole)
const D99R_FIXTURES = [
    ("not decidable, T_s row undefined (point n.a., bound undefined)", d99_ha(decision = "not decidable", failed = Int[], dw = 0.4, dw_lo = 0.1, dt = nothing, dt_lo = nothing), ["ok", "fails"], false),
    ("not decidable, both rows ok (undecided for another reason)", d99_ha(decision = "not decidable", failed = Int[], dw = 0.12, dw_lo = 0.05, dt = 0.2069, dt_lo = 0.158), ["ok", "ok"], true),
    ("not decidable, both points below delta_a", d99_ha(decision = "not decidable", failed = Int[], dw = 0.05, dw_lo = 0.01, dt = 0.08, dt_lo = 0.02), ["fails", "fails"], false)]

@testset "D99 rework: when HA is not decidable, the status of item 3 in the decision line is ok only when both rows of item 3 are ok (hand fixtures)" begin
    for (name, A, rows3, whole) in D99R_FIXTURES
        txt = d99_t5_fixture(A)
        @testset "fixture $name" begin
            @test d97_ok(txt)
            d97_ok(txt) || continue
            rows = [r for t in d97_tables(txt) for r in t[2:end] if length(r) == 6 && r[1] == "3"]
            @test [r[6] for r in rows] == rows3                                                            # guard: the rows follow their own condition (passes now)
            ha = first(l for l in split(txt, '\n') if startswith(l, "HA:"))
            @test startswith(ha, "HA: not decidable.")                                                     # the decision is unchanged (guard; passes now)
            @test d99_i3_fails(ha) == !whole                                                               # item 3 as a whole from the rows
        end
    end
end
