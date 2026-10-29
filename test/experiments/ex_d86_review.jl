# Stage 04 review fixes (2026-10-29): F1 (Double64 provenance of the stretching reducer reaches HDF5 as text), R3 (a spec that no orbit selects is refused naming the spec).
# R1 (the two refined bits of the frequency flag word) is checked by D86i(d) through D86_FLAGS_CHECKED.

@testset "D86 review F1 Double64 stretching records and stores" begin
    run = d86_run("review-f1", d86_one(d86_spec("stretching"); orbits = "h4", precisions = ["Double64"], horizon = 100, mode = "wrapped"))
    @test !isempty(ex_chunk_files(run.store))
    seen = 0
    d86_each(run) do it, o, raw, c, file
        ks = [k for k in keys(raw) if occursin("provenance/stretching@spectrum_lo", k)]
        seen += length(ks)
        @test all(k -> raw[k] isa AbstractString, ks)
    end
    @test seen > 0
end

@testset "D86 review R3 a spec selected by no orbit is refused by name" begin
    pr = ex_write_prereg_run(ex_path("fx-d86-review-r3-prereg.md"), d86_one(d86_spec("tangent"); orbits = "S_L"))
    pp = ex_write_panel(ex_path("fx-d86-review-r3-panel.toml"); orbits = [merge(o, Dict{String, Any}("in_SL" => false)) for o in ex_panel_orbits(:std)])
    r = try
        EX.plan_run(pr, pp; code_version = EX_CODE)
        nothing
    catch e
        e
    end
    @test r isa ArgumentError && occursin("`tg`", r.msg)
end

@testset "D88 pooled counts carry their range and are not written across ranges" begin
    run = d86_run("review-d88", d86_one(d86_spec("stretching"; id = "st"); orbits = "S_D", horizon = 100))
    plan = run.plan
    shared = mixed = 0
    for file in ex_chunk_files(run.store)
        c = plan.chunks[ex_ordinal(file)]
        raw = d86_raw(file)
        os = [plan.items[i].orbit for i in c.items]
        rs = [raw["orbits/$o/indicators/st/scalars@range"] for o in os]
        @test raw["pooled/st@orbits"] == length(os)
        if all(r -> isequal(r, first(rs)), rs)
            shared += 1
            @test haskey(raw, "pooled/st/counts") && raw["pooled/st@range"] == first(rs)
        else
            mixed += 1
            @test !haskey(raw, "pooled/st/counts") && !haskey(raw, "pooled/st@range")
        end
    end
    @test shared > 0 && mixed > 0
end
