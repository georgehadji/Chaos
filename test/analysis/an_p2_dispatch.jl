# D89a dispatch and P1 unchanged (I11), D89t layering. Spec: analysis-spec.md, Amendment D89, E21, I11, criteria D89a and D89t.
# Golden values of the P1 path were taken from the module before Amendment D89, on the P1 fixture with the slug p1-le-benchmark (`p2_p1_study`).

const P1_GOLD_DIR = Ref{Any}(nothing)
p1_gold() = (P1_GOLD_DIR[] === nothing && (P1_GOLD_DIR[] = p2_p1_study(an_path("p1plain"))); P1_GOLD_DIR[])
p1_args(q) = (; prereg = q.prereg, manifest = q.manifest, panel = q.panel, params = q.params, store = q.store)
function p1_copy(q, tag)
    dir = an_path("p1copy-$tag"); rm(dir; recursive = true, force = true); cp(q.dir, dir)
    rm(joinpath(dir, "out"); recursive = true, force = true)
    sub(s) = replace(s, q.dir => dir)
    return merge(q, (; dir, prereg = sub(q.prereg), manifest = sub(q.manifest), panel = sub(q.panel), params = sub(q.params), store = sub(q.store)))
end
# a copy of a P2 study whose slug is replaced everywhere: preregistration (with its hash), manifest and the root attributes of every chunk file
function p2_reslug(st, tag, slug)
    c = p2_copy(st, tag)
    text, h = p2_prereg_text(c.p; slug = slug)
    write(c.prereg, text)
    p2_edit_manifest(c, m -> (m["study"] = slug; m["prereg"] = Dict("header" => h, "sha256_start" => h, "sha256_end" => h)))
    for f in filter(f -> endswith(f, ".h5"), readdir(c.store; join = true))
        p2_set_attr(f, "", "study", slug); p2_set_attr(f, "", "prereg_sha256", h)
    end
    return merge(c, (; prereg_hash = h))
end

@testset "D89a P1 path unchanged (guard: passes before and after the amendment)" begin
    q = p1_gold()
    res = AN.analyze(; p1_args(q)..., out = joinpath(q.dir, "out-gold"), figures = true)
    text = read(res.tables, String)
    # I11: identical to the file before the amendment except the lines "module version" and "git commit"
    ver = r"^- (module version|git commit):"
    norm = join([l for l in split(text, '\n') if !occursin(ver, l)], "\n")
    @test count(l -> occursin(ver, l), split(text, '\n')) == 2
    @test an_sha256(norm) == "0676bfa678c6680d7f3524a1b80db0de606f913295ac979765d29a34ac8acd1e"
    @test res.status === :figures_skipped && basename(res.tables) == "p1-le-benchmark-tables.md" && collect(keys(res.decisions)) == [:H1, :H2, :H3, :H4, :H5]
    gold = Dict("F1" => "525947c4c70b2ce70f80570675e8f64d02af29b4859483c5942c3d538cc27ebf", "F2" => "e9b119437ebee52f70fa066304ce5f5739f4767c626a4f7ee26cbce96fa2c2c7",
                "F3" => "13e6b660a27d05d4599d50562490919a06fc3939cae74f4caa6d8e176b1b21d5", "F4" => "4c93594834b5d099082ed5948571b566b6f3a564b1022c90961449cf9ffa5a39",
                "F5" => "060ded3848689c73173dae7ba7911161d62a59dddd70698ebb57048922380dca", "F6" => "9a4a54b05b97bf4069a887c7bb911df29b9c3df1223c3422441445334fe07366",
                "F7" => "955df07592068c1b11ca03678b4da2d23a1bed464afeef3a9823d2a242f67b49")
    @test res.digests == gold
    # exit codes
    io = IOBuffer()
    @test AN.cli(["analyze", q.prereg, q.manifest, q.panel, q.params, "--store", q.store, "--out", joinpath(q.dir, "out-cli")]; io = io) == 2
    @test occursin("figures skipped", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["analyze", "--help"]; io = io) == 0
    io = IOBuffer(); @test AN.cli(["analyze", q.prereg]; io = io) == 1
    # refusal messages, with the directory replaced
    function refusal(c, f)
        try
            f()
            return "no error"
        catch e
            return replace(replace(sprint(showerror, e), '\\' => '/'), replace(c.dir, '\\' => '/') => "<dir>")
        end
    end
    a = p1_copy(q, "gold-a"); write(a.prereg, read(a.prereg, String) * "x")
    @test refusal(a, () -> AN.analyze(; p1_args(a)..., out = joinpath(a.dir, "o"))) ==
          "ArgumentError: preregistration check failed: <dir>/p1-le-benchmark-preregistration.md: header SHA-256 1bacc65add369700ac81cab785daa4a12f544d8e4cb1eaf3126311c68096a9a5 " *
          "differs from the computed d0368e6e867e18888570a6029ab5505cd7715ba765a202889b4d9d2061cdce1b"
    b = p1_copy(q, "gold-b"); p2_edit_manifest(b, m -> m["status"] = "running")
    @test refusal(b, () -> AN.analyze(; p1_args(b)..., out = joinpath(b.dir, "o"))) == "ArgumentError: manifest status is \"running\", not complete"
    # the P1 parameter table still has exactly its 37 keys and its reader behaves as before
    @test length(AN.read_params(q.params)) == 37 && length(AN.PARAM_KEYS) == 37
    @test an_rejects(() -> AN.read_params(q.params * ".none"), ["does not exist"])
end

@testset "D89a dispatch on the study slug" begin
    st = p2_cached("default")
    # an unknown slug is refused naming the slug and the two known ones; nothing is written
    o = p2_reslug(st, "slug-other", "p3-other")
    @test an_rejects(() -> p2_run(o), ["p3-other", "p1-le-benchmark", "p2-indicators"]) && !ispath(o.out)
    # the header slug and the manifest `study` must agree, for both studies
    m1 = p2_copy(st, "slug-man"); p2_edit_manifest(m1, m -> m["study"] = "p1-le-benchmark")
    @test an_rejects(() -> p2_run(m1), ["study"]) && !ispath(m1.out)
    q = p1_copy(p1_gold(), "slug-man-p1"); p2_edit_manifest(q, m -> m["study"] = "p2-indicators")
    @test an_rejects(() -> AN.analyze(; p1_args(q)..., out = joinpath(q.dir, "o")), ["study"]) && !ispath(joinpath(q.dir, "o"))
    # the P1 parameter table with the P2 preregistration, and the P2 table with the P1 preregistration: refused, the missing and unknown keys named
    @test an_rejects(() -> p2_run(st; params = st.p1.params, out = joinpath(st.dir, "o-x1")), ["budgets", "N_long"]) && !ispath(joinpath(st.dir, "o-x1"))
    g = p1_gold()
    @test an_rejects(() -> AN.analyze(; p1_args(g)..., params = st.params, out = joinpath(g.dir, "o-x2")), ["budgets", "N_long"]) && !ispath(joinpath(g.dir, "o-x2"))
    # the options --p1-* are a usage error for the P1 preregistration (exit code 1, the study named); the keywords are an ArgumentError
    for opt in ("--p1-store", "--p1-preregistration", "--p1-manifest", "--p1-parameters")
        io = IOBuffer()
        @test AN.cli(["analyze", g.prereg, g.manifest, g.panel, g.params, "--store", g.store, "--out", joinpath(g.dir, "o-x3"), opt, g.store]; io = io) == 1
        msg = String(take!(io))
        @test occursin("p2-indicators", msg) && occursin(opt, msg) && !ispath(joinpath(g.dir, "o-x3"))
    end
    @test an_rejects(() -> AN.analyze(; p1_args(g)..., p1_store = g.store, out = joinpath(g.dir, "o-x4")), ["p1"]) && !ispath(joinpath(g.dir, "o-x4"))
    # --out inside data/ refuses (I12); the folder is not created
    bad = joinpath(@__DIR__, "..", "..", "data", "an-p2-refused-out")
    try
        @test an_rejects(() -> p2_run(st; out = bad), ["data"]) && !ispath(bad)
        io = IOBuffer()
        @test AN.cli(p2_cli_args(st; extra = ["--out", bad]); io = io) == 1 && !ispath(bad)
    finally
        ispath(bad) && rm(bad; recursive = true, force = true)         # only this folder, and only if a defective module created it
    end
end

@testset "D89t layering and dependencies" begin
    root = joinpath(@__DIR__, "..", "..")
    src = joinpath(root, "src")
    ana = filter(f -> endswith(f, ".jl"), readdir(joinpath(src, "Analysis"); join = true))
    lines = [(basename(f), l) for f in ana for l in eachline(f) if occursin(r"^\s*(?:using|import)\s", l)]
    @test !isempty(lines)
    allowed = r"^\s*(?:using|import)\s+(?:\.\.Experiments|\.\.Foundation|SHA|StableRNGs|TOML|Statistics|Random|LinearAlgebra|Printf|Dates)\b"
    bad = [x for x in lines if !occursin(allowed, x[2]) && !occursin(r"^\s*(?:using|import)\s+\.\.(?:Experiments|Foundation)\b", x[2])]
    @test isempty(bad)
    @test isempty([x for f in ana for x in [l for l in eachline(f) if occursin(r"Makie|Plots|PyPlot|Figures", l) && occursin(r"^\s*(?:using|import)\s", l)]])
    # no file of a lower layer imports Analysis
    lower = ["Foundation", "Systems", "Integrators", "Tangent", "Lyapunov", "Indicators", "Diagnostics", "Shadow", "TimeSeries", "Experiments"]
    up = String[]
    for d in lower, f in (isdir(joinpath(src, d)) ? filter(f -> endswith(f, ".jl"), readdir(joinpath(src, d); join = true)) : String[]), l in eachline(f)
        occursin(r"^\s*(?:using|import)\s.*\bAnalysis\b", l) && push!(up, "$d/$(basename(f))")
    end
    @test isempty(up)
    # the P2 pipeline is part of the module (this is the one assertion of the testset that fails before the amendment is implemented)
    @test AN.p2_read_params isa Function && AN.P2_PARAM_KEYS isa Tuple && length(AN.P2_PARAM_KEYS) == 58
end
