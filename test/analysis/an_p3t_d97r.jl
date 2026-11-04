# Amendment D97, rework after the stage 04 review (2026-11-03): the four MEDIUM findings. E61 k and E77 (items 1 and 2), E74 / E61 e (item 3), E78 (item 4).
# Helpers and fixture functions are those of an_p3t_d97.jl and an_p3t_helpers.jl; every testset names its finding. RED means the reader accepts what it must refuse, or misreads the store.

# ---------------------------------------------------------------- finding 1: a chk item missing for an orbit of D_B whose chunk did not fail is refused (E61 k, E77)
@testset "D97r1 a chk item (Double64 or BigFloat) missing from a chunk that is not listed as failed is refused, naming the orbit and the item" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    # orbits 6 and 8 share the chk cell k3:C:r1; orbit 8 is dropped from one precision only (no failure in the manifest), so it is known to be in D_B by its other item
    for (prec, word) in (("Double64", "Double64"), ("BigFloat(256)", "BigFloat"))
        c = p3t_clone(p3t_full(), "d97r1-$word")
        p3t_alter_chunk(S, c, I, only(p3t_ords(S; arm = "chk", cell = "k3:C:r1", precision = prec)); drop = [8])
        @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["8", word, "chk"]) && !ispath(c.out)
    end
end

@testset "D97r1 a chk chunk failed twice keeps its readings: Double64 failed leaves D_B (no T3 row), BigFloat failed shows n.a. (E77, P118; guard)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97r1-bf")
    p3t_fail_chunks(S, c, p3t_ords(S; arm = "chk", cell = "k3:C:r1", precision = "BigFloat(256)"))
    R = d97_run(S, I, st, c)
    @test d97_ok(R)
    d97_ok(R) || return
    t3 = d97_t3(d97_read(R).tab)
    row(i) = (k = findfirst(r -> r[1] == string(i), t3[2:end]); k === nothing ? nothing : t3[k + 1])
    @test all(i -> row(i) !== nothing && all(j -> tryparse(Float64, row(i)[j]) !== nothing, 3:5), (3, 12))
    @test all(i -> row(i) !== nothing && tryparse(Float64, row(i)[3]) !== nothing && tryparse(Float64, row(i)[4]) !== nothing && tryparse(Float64, row(i)[5]) === nothing, (6, 8))
end

# ---------------------------------------------------------------- finding 2: a D orbit missing a Double64 item of a chunk that did not fail is refused (E77, E61 k)
@testset "D97r2 a Double64 item missing from a d64 or d64_pa1 chunk that is not listed as failed is refused, naming the orbit and the item" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    cell4 = p3t_cell(first(o for o in S.orbits if o.id == 4))                 # orbit 4 is in D; its cell holds the orbits 3 4 5
    for (arm, spec) in (("d64", "tg"), ("d64", "fq1"), ("d64_pa1", "pa1"))
        c = p3t_clone(p3t_full(), "d97r2-$spec")
        p3t_alter_chunk(S, c, I, only(p3t_ords(S; arm, spec, cell = cell4)); drop = [4])
        @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["4", spec]) && !ispath(c.out)
    end
end

@testset "D97r2 orbits of a failed d64 chunk leave D and are counted by |D|: no refusal, D has the orbits of the other cells (E70, P104 b; guard)" begin
    R = d97_multi()
    @test d97_ok(R.res)
    d97_ok(R.res) || return
    @test R.m["sets"]["counts"]["D"] == length(d97_ids(R.m, "D")) == 8 && !any(i -> i in d97_ids(R.m, "D"), [3, 4, 5, 6, 7, 8])
end

# ---------------------------------------------------------------- finding 3: the panel key class (E74, E61 e)
@testset "D97r3 a panel orbit without the key class, or with a class that is not a non-empty string, is an ArgumentError naming the key and the orbit" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    for (tag, edit) in (("none", o -> delete!(o, "class")), ("empty", o -> o["class"] = ""), ("number", o -> o["class"] = 3))
        inst = p3t_instance(S, I, "d97r3-$tag"; panel_edit = d -> edit(first(o for o in d["orbit"] if o["id"] == 5)))
        r = d97_call(() -> AN.p3t_read_panel(inst.panel, 3))
        @test r isa ArgumentError && occursin("class", sprint(showerror, r)) && occursin("orbit 5", sprint(showerror, r))
    end
end

@testset "D97r3 analyze and count-round refuse a panel without class, naming the key and the orbit" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_instance(S, I, "d97r3-run"; panel_edit = d -> delete!(first(o for o in d["orbit"] if o["id"] == 5), "class"))
    write(c.log, join(first(p3t_log_entries(S, I; panel = c.panel_hash))))
    r = d97_call(() -> AN.analyze(; p3t_an_args(I, st, c)..., out = joinpath(c.dir, "out")))
    @test r isa ArgumentError && occursin("class", sprint(showerror, r)) && occursin("orbit 5", sprint(showerror, r)) && !ispath(joinpath(c.dir, "out"))
    r = d97_call(() -> AN.count_round(; p3t_cr_args(I, st, c)...))
    @test r isa ArgumentError && occursin("class", sprint(showerror, r)) && occursin("orbit 5", sprint(showerror, r))
end

# ---------------------------------------------------------------- finding 4: a spec with 0.0 measured ns per orbit-step is no value (E78)
@testset "D97r4 a spec whose measured_ns_per_orbit_step is 0.0 makes sensitivity (iii) not decidable, naming the spec; the other specs stay" begin
    for (arm, spec) in (("ind", "tg"), ("ind_pa1", "pa1"))
        R = d97_costs("d97r4-$spec", m -> m["arms"][arm]["specs"][spec]["measured_ns_per_orbit_step"] = 0.0)
        @test d97_ok(R.res)
        d97_ok(R.res) || continue
        sv = R.m["sensitivities"]["measured"]
        others = filter(!=(spec), ["tg", "fq1", "bk1", "pa1"])
        @test sv["decidable"] == false && sv["verdict"] == "not decidable" && occursin(spec, sv["reason"]) && !any(o -> occursin(o, sv["reason"]), others) && !occursin("L1", sv["reason"])
        t10 = d97_t10(R.tab)
        @test tryparse(Float64, t10[spec][2]) === nothing && all(o -> tryparse(Float64, t10[o][2]) !== nothing, others)
    end
end

# ---------------------------------------------------------------- the LOW findings 5 to 12 (owner decision 2026-11-03: fix them too); 9 is E64 b as written and 11 is a refactor, both without a testset here
@testset "D97r5 analyze without the P1 parameter table refuses when D_B is not empty, naming the P1 parameter table; nothing is written" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97r5")
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)..., p1_params = nothing), ["P1 parameter", "D_B"]) && !ispath(c.out)
end

@testset "D97r6 two chk items of one orbit and one precision without its bit count (BigFloat(256) and BigFloat(512)) are refused, naming the orbit and the precision" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97r6")
    p3t_alter_chunk(S, c, I, only(p3t_ords(S; arm = "chk", cell = "k3:C:r1", precision = "Double64")); edit = Dict{String, Any}("precision" => "BigFloat(512)"),
                    sp_edit = s -> merge(s, (; precision = "BigFloat(512)")))
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), ["6", "BigFloat", "twice"]) && !ispath(c.out)
end

@testset "D97r7 a chunk of an arm outside the six arms of the test run is refused, naming the chunk file and the arm" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97r7")
    o = first(p3t_ords(S; arm = "chk"))
    p3t_set_attr(joinpath(c.store, d97_file(o)), "arms", ["zz"])
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), [d97_file(o), "\"zz\""]) && !ispath(c.out)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, c)...), [d97_file(o), "\"zz\""])
end

@testset "D97r8 a failed chunk whose cell is not of the form of E74 is refused, naming the chunk file and the cell" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "d97r8")
    o = only(p3t_ords(S; arm = "d64", spec = "tg", cell = "k2:C:r1"))
    p3t_fail_chunks(S, c, [o])
    f = joinpath(c.store, d97_file(o))
    p3t_set_attr(f, "chunk_key", replace(AX.read_chunk_attrs(f)["chunk_key"], r"cell=\S+" => "cell=bogus"))
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)...), [d97_file(o), "\"bogus\""]) && !ispath(c.out)
end

@testset "D97r10 the T2 sentence is the sentence of D95 (I50): the orbits of the panel file and the rounds counted" begin
    R = d97_multi()
    @test d97_ok(R.res)
    d97_ok(R.res) || return
    n = length(TOML.parsefile(R.inst.panel)["orbit"])
    @test occursin("The test draw. Orbits in the panel file: $n; rounds ", d97_sec(R.tab, 2)) && !occursin("the panel file holds", R.tab)
end

@testset "D97r12 a manifest chunk entry without a field that analyze reads is refused, naming the field and the ordinal" begin
    S = p3t_cached()
    o = first(c.ordinal for c in S.chunks if c.arm == "ref" && c.spec == "L1" && c.status == "ok")
    for k in ("steps", "thread_seconds", "key", "arms")
        R = d97_costs("d97r12-$k", m -> delete!(first(e for e in m["chunks"] if e["ordinal"] == o), k))
        @test R.res isa ArgumentError && occursin("\"$k\"", sprint(showerror, R.res)) && occursin("$o", sprint(showerror, R.res)) && !ispath(R.c.out)
    end
end
