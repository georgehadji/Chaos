# AC22 command, AC20 layering and dependencies, AC18 extension wiring. Spec: analysis-spec.md, Compute ("Command", "Layering and dependencies").

an_cli_args(st; extra = String[]) = vcat(["analyze", st.prereg, st.manifest, st.panel, st.params, "--store", st.store, "--out", st.out], extra)

@testset "AC22 command" begin
    st = an_copy(an_cached("default"), "cli")
    # complete inputs without the plotting extension: tables written, exit code 2, the missing figures named
    io = IOBuffer()
    @test AN.cli(an_cli_args(st); io = io) == 2
    out = String(take!(io))
    @test occursin("figures skipped", out) && isfile(joinpath(st.out, "p1-le-benchmark-tables.md"))
    # help
    io = IOBuffer()
    @test AN.cli(["analyze", "--help"]; io = io) == 0
    @test occursin("analyze", String(take!(io)))
    # missing or unknown options, too few arguments: usage text, exit code 1, nothing written
    for (tag, args) in (("value", ["analyze", st.prereg, st.manifest, st.panel, st.params, "--store"]),
                        ("unknown", vcat(["analyze", st.prereg, st.manifest, st.panel, st.params, "--bogus", "x"])),
                        ("few", ["analyze", st.prereg]),
                        ("none", ["analyze"]),
                        ("command", ["frobnicate"]))
        io = IOBuffer()
        out2 = an_path("cli-out-$tag")
        a = copy(args)
        tag in ("value", "unknown") && append!(a, ["--out", out2])
        @test AN.cli(a; io = io) == 1
        msg = String(take!(io))
        @test (occursin("analyze", msg) || occursin("usage", lowercase(msg))) && !ispath(out2)
    end
    # a refusal (bad preregistration): exit code 1, the cause named, nothing written
    bad = an_copy(an_cached("default"), "cli-bad")
    write(bad.prereg, read(bad.prereg, String) * "x")
    io = IOBuffer()
    @test AN.cli(an_cli_args(bad); io = io) == 1
    @test occursin("preregistration", lowercase(String(take!(io)))) && !ispath(bad.out)
    # the entry point routes the command to the module
    entry = read(joinpath(@__DIR__, "..", "..", "bin", "chaoslab.jl"), String)
    @test occursin("analyze", entry) && occursin("Analysis", entry)
end

@testset "AC20 layering and dependencies" begin
    root = joinpath(@__DIR__, "..", "..")
    src = joinpath(root, "src", "Analysis")
    @test isdir(src)
    files = isdir(src) ? filter(f -> endswith(f, ".jl"), readdir(src; join = true)) : String[]
    @test !isempty(files)
    imports = String[]
    plotting = String[]
    for f in files, line in eachline(f)
        occursin(r"^\s*(?:using|import)\s", line) && push!(imports, "$(basename(f)): $line")
        occursin(r"Makie|CairoMakie|GLMakie", line) && push!(plotting, "$(basename(f)): $line")
    end
    @test !isempty(imports)
    bad = [l for l in imports if occursin(r"(?:\.\.|Chaoslab\.)(?:Systems|Integrators|Tangent|Lyapunov|Indicators|Diagnostics|Shadow|TimeSeries|Figures)\b", l)]
    @test !isempty(files) && isempty(bad)
    @test !isempty(files) && isempty(plotting)                                       # Analysis never loads the plotting package
    # the layer list of the layering test
    lay = read(joinpath(@__DIR__, "..", "layering.jl"), String)
    @test occursin(r"\[\"Experiments\"\],\s*\[\"Analysis\"\],\s*\[\"Figures\"\]", lay)
    # the module is part of the package
    @test AN.read_params isa Function
    # Figures: a package extension of CairoMakie (weak dependency), like CUDA
    proj = TOML.parsefile(joinpath(root, "Project.toml"))
    @test haskey(get(proj, "weakdeps", Dict()), "CairoMakie") && get(get(proj, "extensions", Dict()), "ChaoslabCairoMakieExt", nothing) == "CairoMakie"
    @test haskey(get(proj, "compat", Dict()), "CairoMakie")
    @test isfile(joinpath(root, "ext", "ChaoslabCairoMakieExt.jl"))
    # the verify environment has the plotting package
    vproj = TOML.parsefile(joinpath(root, "verify", "Project.toml"))
    @test haskey(vproj["deps"], "CairoMakie")
end
