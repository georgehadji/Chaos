# D93a dispatch and the other pipelines unchanged (I36), D93q layering and dependencies. Spec: analysis-spec.md, Amendment D93, E46, I36, criteria D93a and D93q.
# Golden values of the P2 path were taken from the module at HEAD b2bd5cf (D89 and D91 implemented, D93 not), on the P2 fixture study of an_p2_study.jl; those of the P1 path are
# pinned by the existing guard "D89a P1 path unchanged (guard)" in an_p2_dispatch.jl, which this file does not repeat.

p3_norm(path) = join([l for l in split(read(path, String), '\n') if !occursin(r"^- (module version|git commit):", l)], "\n")

@testset "D93a calibrate is refused for the preregistration of any other study, naming the slug; nothing is written" begin
    st = p3_cached("default")
    # the p2 and p1 preregistrations, and a preregistration with an unknown slug (consistent hash rule)
    c = p3_copy(st, "dispatch-p2")
    @test an_rejects(() -> p3_run(c; prereg = c.s2.prereg), ["p2-indicators"]) && !ispath(c.out)
    g = p1_gold()
    @test an_rejects(() -> p3_run(c; prereg = g.prereg), ["p1-le-benchmark"]) && !ispath(c.out)
    text, _ = p3_prereg_text(c.q; slug = "p3-other")
    pre = joinpath(c.dir, "p3-other-preregistration.md"); write(pre, text)
    @test an_rejects(() -> p3_run(c; prereg = pre), ["p3-other"]) && !ispath(c.out)
    # the slug of the preregistration: absent, several, invalid (the rule of D89 E21)
    for (tag, edit) in (("absent", t -> replace(t, r"(?m)^\| Study slug \|[^\n]*\n" => "")), ("several", t -> replace(t, "| Study slug | p3-abstention |\n" => "| Study slug | p3-abstention |\n| Study slug | p3-abstention |\n")),
                        ("invalid", t -> replace(t, "| Study slug | p3-abstention |" => "| Study slug | P3 Abstention |")))
        d = p3_copy(st, "dispatch-slug-$tag")
        write(d.prereg, edit(read(d.prereg, String)))
        @test an_rejects(() -> p3_run(d), ["slug"]) && !ispath(d.out)
    end
    # unknown options and missing arguments give the usage text and exit code 1; --help exits with 0 and names the command
    io = IOBuffer(); @test AN.cli(["calibrate"]; io = io) == 1 && occursin("calibrate", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["calibrate", "--help"]; io = io) == 0 && occursin("calibrate <preregistration> <parameters>", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["calibrate", st.prereg]; io = io) == 1 && occursin("usage", lowercase(String(take!(io))))
    io = IOBuffer(); @test AN.cli(["calibrate", st.prereg, st.params, "--bogus", "x"]; io = io) == 1 && occursin("--bogus", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["calibrate", st.prereg, st.params, "--out"]; io = io) == 1
    io = IOBuffer(); @test AN.cli(["calibrate", st.prereg, st.params, "--threads", "many"]; io = io) == 1
    io = IOBuffer(); @test AN.cli(["calibrate", st.prereg, st.params, "extra-argument"]; io = io) == 1 && !ispath(st.out)
    # the help text names the options of the command (E46)
    io = IOBuffer(); AN.cli(["calibrate", "--help"]; io = io); h = String(take!(io))
    @test all(o -> occursin(o, h), ("--out", "--threads", "--p2-store", "--p2-manifest", "--p2-preregistration", "--p2-parameters", "--panel", "--p1-store", "--p1-manifest", "--p1-preregistration", "--p1-parameters"))
    # the command line with a preregistration of another study: exit code 1, the slug named, nothing written
    io = IOBuffer(); @test AN.cli(vcat(["calibrate", g.prereg], p3_cli_args(c)[3:end]); io = io) == 1 && occursin("p1-le-benchmark", String(take!(io))) && !ispath(c.out)
    io = IOBuffer(); @test AN.cli(vcat(["calibrate", c.s2.prereg], p3_cli_args(c)[3:end]); io = io) == 1 && occursin("p2-indicators", String(take!(io))) && !ispath(c.out)
    # an --out inside data/ refuses (I30), the folder is not created
    bad = joinpath(@__DIR__, "..", "..", "data", "an-p3-refused-out")
    try
        @test an_rejects(() -> p3_run(c; out = bad), ["data"]) && !ispath(bad)
        io = IOBuffer()
        @test AN.cli(p3_cli_args(c; extra = ["--out", bad]); io = io) == 1 && !ispath(bad)
    finally
        ispath(bad) && rm(bad; recursive = true, force = true)         # only this folder, and only if a defective module created it
    end
end

@testset "D93a analyze with the preregistration of this study is refused (exit code 1) naming the slug and the command calibrate" begin
    st = p3_cached("default")
    c = p3_copy(st, "dispatch-analyze")
    kw = (; prereg = c.prereg, manifest = c.s2.manifest, panel = c.s2.panel, params = c.params, store = c.s2.store, out = c.out)
    @test an_rejects(() -> AN.analyze(; kw...), ["p3-abstention", "calibrate"]) && !ispath(c.out)
    io = IOBuffer()
    @test AN.cli(["analyze", c.prereg, c.s2.manifest, c.s2.panel, c.params, "--store", c.s2.store, "--out", c.out]; io = io) == 1
    msg = String(take!(io))
    @test occursin("p3-abstention", msg) && occursin("calibrate", msg) && !ispath(c.out)
    # with the options of the P2 pipeline too
    io = IOBuffer()
    @test AN.cli(vcat(["analyze", c.prereg, c.s2.manifest, c.s2.panel, c.params, "--store", c.s2.store, "--out", c.out, "--p1-store", c.p1.store]); io = io) == 1 && !ispath(c.out)
    # the unknown-slug message of analyze for every other slug is unchanged (D89 E21)
    o = p2_reslug(p2_cached("default"), "slug-other-p3", "p3-other")
    @test an_rejects(() -> p2_run(o), ["p3-other", "p1-le-benchmark", "p2-indicators"]) && !ispath(o.out)
end

@testset "D93a I36 the P1 and P2 pipelines are unchanged: tables file, figure digests, exit codes and messages equal the golden values taken before the amendment" begin
    st = p2_cached("default")
    res = p2_run(st; out = joinpath(st.dir, "gold-out-d93"))
    @test an_sha256(p3_norm(res.tables)) == "e28a1342eea654413142fe46c8cf74d46d480d1c6f64d19a93633acde18d8316" && count(l -> occursin(r"^- (module version|git commit):", l), split(read(res.tables, String), '\n')) == 2
    @test res.digests == Dict("F1" => "306bfe2d7f8e5d70c5a5667e0bc5a139ceffaf1732626e58afab1e43777777ef", "F2" => "730cf8a01a292dba67c198f217c53ee1702311bb22bfcb020bd1d6e66114e5cd", "F3" => "e1709283783531bf5f629d78a548ebdc97392b101634a4c24d0a666cbf38c293", "F4" => "5369fa28d8e321d8133372a3482b46e55b9fab34ecc7551d79a5de55db9051e8", "F5" => "cc0380d0d5237cab6321ca8f4f7ed80264f53a8cdf209f2399975344e27afccc", "F6" => "fbc4ffa18abbeaacc79d5d4fb43b2cb130c6c71632553a43ff5177427bf9b29b")
    io = IOBuffer()
    @test AN.cli(p2_cli_args(st; extra = ["--out", joinpath(st.dir, "gold-out-d93-cli")]); io = io) == 2 && occursin("figures skipped", String(take!(io)))
    c = p2_copy(st, "gold-ref-d93"); p2_edit_manifest(c, m -> m["status"] = "running")
    msg = try p2_run(c); "no error" catch e; replace(replace(sprint(showerror, e), '\\' => '/'), replace(c.dir, '\\' => '/') => "<dir>") end
    @test msg == "ArgumentError: manifest status is \"running\", not complete"
    # the P1 pipeline: the existing guard compares the P1 tables file with its golden hash; here the exit codes and the usage of the P1 command
    q = p1_gold()
    io = IOBuffer(); @test AN.cli(["analyze", q.prereg, q.manifest, q.panel, q.params, "--store", q.store, "--out", joinpath(q.dir, "out-d93")]; io = io) == 2
    io = IOBuffer(); @test AN.cli(["analyze", "--help"]; io = io) == 0
    io = IOBuffer(); @test AN.cli(["analyze", q.prereg]; io = io) == 1
end

@testset "D93q layering and dependencies: Analysis imports Foundation, Experiments and the standard library only; no plotting package; no lower layer imports Analysis; no new package" begin
    root = joinpath(@__DIR__, "..", "..")
    src = joinpath(root, "src")
    ana = filter(f -> endswith(f, ".jl"), readdir(joinpath(src, "Analysis"); join = true))
    lines = [(basename(f), l) for f in ana for l in eachline(f) if occursin(r"^\s*(?:using|import)\s", l)]
    allowed = r"^\s*(?:using|import)\s+(?:\.\.Experiments|\.\.Foundation|SHA|StableRNGs|TOML|LinearAlgebra|Printf|Dates)\b"
    bad = [x for x in lines if !occursin(allowed, x[2])]
    @test isempty(bad)
    @test isempty([x for f in ana for x in [l for l in eachline(f) if occursin(r"Makie|Plots|PyPlot|Gadfly|PlotlyJS|UnicodePlots|VegaLite|Gnuplot", l) && occursin(r"^\s*(?:using|import)\s", l)]])
    @test isempty([(basename(f), l) for f in ana for l in eachline(f) if occursin(r"Base\.require|@eval\s+using|Requires|Pkg\.add", l) && !startswith(strip(l), "#")])
    lower = ["Foundation", "Systems", "Integrators", "Tangent", "Lyapunov", "Indicators", "Diagnostics", "Shadow", "TimeSeries", "Experiments"]
    up = String[]
    for d in lower, f in (isdir(joinpath(src, d)) ? filter(f -> endswith(f, ".jl"), readdir(joinpath(src, d); join = true)) : String[]), l in eachline(f)
        occursin(r"^\s*(?:using|import)\s.*\bAnalysis\b", l) && push!(up, "$d/$(basename(f))")
    end
    @test isempty(up)
    # no new package: the dependency lists of the project are those of the release (the same list the D91 test pins), and the test target is unchanged
    proj = TOML.parsefile(joinpath(root, "Project.toml"))
    @test Set(keys(proj["deps"])) == Set(["Dates", "DoubleFloats", "HDF5", "JSON", "KernelAbstractions", "LinearAlgebra", "PrecompileTools", "SHA", "StableRNGs", "StaticArrays", "TOML"])
    @test Set(keys(proj["weakdeps"])) == Set(["CairoMakie", "CUDA"]) && proj["targets"]["test"] == ["Aqua", "ForwardDiff", "JET", "Test"] && Set(keys(proj["extras"])) == Set(["Aqua", "ForwardDiff", "JET", "Test"])
    # the module has the D93 entry points (the assertions of this testset that fail before the amendment is implemented)
    @test AN.calibrate isa Function && AN.P3_PARAM_KEYS isa Tuple && length(AN.P3_PARAM_KEYS) == 71 && AN.p3_read_params isa Function
end
