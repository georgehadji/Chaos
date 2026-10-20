# Regressions from the independent review of stage 04 (new tests, added after the approved suite): failed chunks, the common half point of the ladder,
# a control without its scope, parameter ranges, the provenance in the header.

@testset "review: a chunk that failed twice is refused with a named cause" begin
    s = an_copy(an_cached("default"), "rev-failed")
    an_edit_manifest(s, m -> m["chunks"][3]["status"] = "failed")
    @test an_rejects(() -> an_run(s), ["chunk", "failed"]) && !ispath(s.out)
end

@testset "review: the ladder uses one half point for the three precisions" begin
    # a spike of the lifted x at step 500, which only the BigFloat item records (the Float64 and Double64 grids have 501 and not 500 there)
    spike(id, st, steps) = (id == 9 && 500 in steps) ? (s = copy(st); s[findfirst(==(500), steps), 1] += 0.5; s) : st
    st = an_cached("halfpoint"; state3 = spike)
    res = an_res("halfpoint", st)
    @test res.labels[9].ladder_same === true && !(:ladder in res.labels[9].forced)
    @test (res.labels[9].instability, res.labels[9].structure) == (AN_NONE, "quasiperiodic")
end

@testset "review: a control with a positive exponent and no scope is refused" begin
    s = an_copy(an_cached("default"), "rev-scope")
    i = an_find(s; method = "L5")
    AX.HDF5.h5open(an_chunk(s, i), "r+") do h
        AX.HDF5.delete_attribute(h["orbits/17"], "control_scope")
    end
    @test an_rejects(() -> an_run(s), ["control_scope", "17"]) && !ispath(s.out)
end

@testset "review: parameter ranges and a malformed manifest" begin
    for (k, v, word) in ((:bootstrap_level, 1.5, "bootstrap_level"), (:N_1, 0, "N_1"), (:min_set, 0, "min_set"), (:h2_median, NaN, "h2_median"))
        f = an_write_params(an_path("range-$k.toml"), an_params(; (k => v,)...))
        @test an_rejects(() -> AN.read_params(f), [word])
    end
    s = an_copy(an_cached("default"), "rev-manifest")
    an_edit_manifest(s, m -> delete!(m, "chunks"))
    @test an_rejects(() -> an_run(s), ["manifest", "chunk"]) && !ispath(s.out)
end

@testset "review: the header names the module version and the commit" begin
    res = an_res("default", an_cached("default"))
    ver = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["version"]
    @test res.header.module_version == ver && occursin(r"^([0-9a-f]{40}|unknown)$", res.header.git_sha)
    text = read(res.tables, String)
    @test occursin("module version: $ver", text) && occursin("git commit: $(res.header.git_sha)", text)
    @test any(r -> r[1] == "W" && r[2] == 0, res.figures["F2"].rows)               # the pooled series over W
end
