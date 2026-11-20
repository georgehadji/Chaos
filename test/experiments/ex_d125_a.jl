# Amendment D125, stage 02: D125a (absence of change, I40, O28) and D125u (layering at the package level).
# The golden values are in ex_d125_golden.jl (code before the amendment); these tests PASS before stage 03 and must keep passing after it.
# The three plans take about 5 minutes together (45 148 + 26 964 + 15 569 items, measured: 143 s, 79 s, 197 s with compilation); the rebuild of the real panels
# (37 minutes for p1-le-benchmark) is run only with CHAOSLAB_D125_FULL=1 (the criterion's [full] tag, stage 04).

const D125_GOLD = include(joinpath(@__DIR__, "ex_d125_golden.jl"))
const D125_FULL = get(ENV, "CHAOSLAB_D125_FULL", "0") == "1"

@testset "D125a golden records: plans, chunk hashes, item texts and block hashes of the three studies (I40)" begin
    for study in ("p1-le-benchmark", "p2-indicators", "p3-abstention")
        gold = D125_GOLD.studies[study]
        pre, pan, _ = d125_study_paths(study)
        @test bytes2hex(sha256(read(pre))) == gold.prereg_sha && bytes2hex(sha256(read(pan))) == gold.panel_sha          # the tracked inputs are the ones captured
        log = nothing
        if study == "p3-abstention"
            log = ex_path("d125-p3-log.txt"); write(log, D125_GOLD.p3_log)
        end
        got = d125_golden_capture(study; log = log)
        for f in keys(gold)
            @test getproperty(got, f) == getproperty(gold, f)                                                           # 0 differences, field by field
        end
    end
end

@testset "D125a the real manifests of the three studies agree with the golden plan counts" begin
    out = joinpath(D125_ROOT, "workspaces", "study", "stages", "03-run", "output")
    for study in ("p1-le-benchmark", "p2-indicators", "p3-abstention")
        m = EX.read_json(joinpath(out, study * "-manifest.json"))
        @test m["plan_counts"]["items"] == D125_GOLD.studies[study].items && m["plan_counts"]["chunks"] == D125_GOLD.studies[study].chunks
        @test length(m["chunks"]) == D125_GOLD.studies[study].chunks
    end
end

@testset "D125a the fixture goldens of D86b, D94b and D96b are unchanged and still asserted" begin
    # the guards live in ex_d86_grammar.jl, ex_d94_grammar.jl and ex_d96.jl, which this amendment does not edit (D125a: no existing test changes)
    @test isfile(joinpath(@__DIR__, "ex_d86_golden.jl")) && isfile(joinpath(@__DIR__, "ex_d94_golden.jl"))
    g = include(joinpath(@__DIR__, "ex_d86_golden.jl"))
    @test length(g.item_hashes) == 52 && length(g.chunk_hashes) == 18
    @test sort(collect(keys(include(joinpath(@__DIR__, "ex_d94_golden.jl")).panel_sha))) == ["base", "k0-empty", "k0-full", "noell"]
end

@testset "D125a rebuilt real panels equal the tracked panels (full)" begin
    if !D125_FULL
        @info "D125a panel rebuild skipped (set CHAOSLAB_D125_FULL=1); stage 04 [full] runs it"
        @test true
    else
        skip = ("git_sha", "git_dirty", "package_version", "julia_version", "stablerngs_version")
        canon = p -> (d = TOML.parsefile(p); foreach(k -> delete!(d["manifest"], k), skip); io = IOBuffer(); TOML.print(io, d; sorted = true); ex_sha256(take!(io)))
        for study in ("p1-le-benchmark", "p3-abstention")
            pre, pan, _ = d125_study_paths(study)
            out = ex_path("d125-rebuilt-$study-panel.toml")
            cd(D125_ROOT) do
                EX.build_panel(relpath(pre, D125_ROOT); out = out)
            end
            @test canon(out) == canon(pan)
        end
    end
end

@testset "D125u layering and entry point" begin
    src = joinpath(@__DIR__, "..", "..", "src")
    # Experiments imports the Integrators layer (the F3 one-step map); the layering test of AC28 allows it already
    mod = read(joinpath(src, "Experiments", "Experiments.jl"), String)
    @test occursin(r"(?m)^\s*(?:using|import)\s+\.\.Integrators", mod)                                       # RED: not imported yet
    # the entry point gains no option and holds no numerics and no amendment name
    cli = read(joinpath(src, "Experiments", "run_cli.jl"), String)
    @test !occursin(r"D125|froeschle|henon_heiles", cli)
    # no study value of the draft preregistration in the module: its seed, sizes and energies are not hard-coded (the verify/ script of stage 04 repeats this search)
    code = join((read(joinpath(src, "Experiments", f), String) for f in readdir(joinpath(src, "Experiments")) if endswith(f, ".jl")), "\n")
    @test !occursin("p3-higher-dimensions", code)
    # nothing below Experiments imports it (AC28), and Experiments imports no Analysis code
    @test !occursin(r"(?m)^\s*(?:using|import)\s+\.\.Analysis", code)
end
