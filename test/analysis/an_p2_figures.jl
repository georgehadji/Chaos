# D91 figure drawers of study p2-indicators: the module side. Spec: analysis-spec.md, Amendment D91 (E38 to E45, I22 to I29, criteria D91a to D91o).
# CairoMakie is not a test dependency, so the package suite tests only what the analysis layer does around a drawer: it calls the function in
# `Chaoslab.Analysis.FIGURE_DRAWER[]` once, after the tables file is written, with a result that carries `.study` and `.figures` (F1 to F6, each with `.columns` and `.rows`
# like res.figures["F7"] of P1), and it turns the outcome into a status and a command line. The pictures (D91a, D91f to D91m) belong to verify/analysis-p2-figures.jl.
# The drawer here is a FAKE assigned to FIGURE_DRAWER[] (the previous value is restored in `finally`). No random draw: every value is the fixture's own (no seed needed).
# Reused: the P2 fixture (`p2_cached("default")`, `p2_run`, `p2_cli_args`), the P1 fixture (`an_cached("default")`, `an_run`, `an_cli_args`), `p2_digest(cols, rows)` of the
# independent reference (an_p2_ref.jl), `an_files`, `an_path`.
import Logging

# ---------------------------------------------------------------- fakes and runners
const FIG_MSG = "F6 row 12: step value 0 cannot be drawn on a logarithmic axis"
const FIG_P1_FAILED = "figures failed: the tables are written, see the warning above"
const FIG_P1_SKIPPED = "figures skipped: the plotting package is not loaded (run under the verify environment)"

# FIGURE_DRAWER[] = drawer while f() runs; the previous value comes back whatever happens
function fig_with(f, drawer)
    old = AN.FIGURE_DRAWER[]
    AN.FIGURE_DRAWER[] = drawer
    try
        return Logging.with_logger(f, Logging.NullLogger())
    finally
        AN.FIGURE_DRAWER[] = old
    end
end

# what the drawer received, copied at call time: the tables of `res.figures` (nothing when the result has no usable `.figures`), their digests by the module and by the reference
function fig_snapshot(res)
    try
        f = res.figures
        return Dict{String, Any}(String(k) => (columns = String.(collect(f[k].columns)), rows = [collect(r) for r in f[k].rows],
                                               digest = AN.p2_digest(f[k]), digest_ref = p2_digest(f[k].columns, f[k].rows)) for k in sort(collect(keys(f))))
    catch
        return nothing
    end
end

# a fake drawer: writes nothing, records its calls, raises `fail` (a message) when given
function fig_fake(; fail::Union{Nothing, String} = nothing)
    calls = Any[]
    drawer = function (res, out)
        tf = isdir(out) ? filter(endswith("-tables.md"), readdir(out; join = true)) : String[]
        push!(calls, (res = res, out = out, study = try String(res.study) catch; nothing end, snap = fig_snapshot(res), files = an_files(out),
                      tables = length(tf) == 1 ? read(tf[1]) : nothing))
        fail === nothing || throw(ErrorException(fail))
        return nothing
    end
    return drawer, calls
end

fig_lines(text) = filter(!isempty, strip.(split(text, '\n')))
fig_tables(out, slug) = (f = joinpath(out, slug * "-tables.md"); isfile(f) ? read(f) : UInt8[])

# analyze on the P2 fixture, in a folder of its own (the inputs are read only, I12)
function fig_p2(tag; drawer, figures = true)
    st = p2_cached("default"); out = an_path("p2fig-$tag")
    rm(out; recursive = true, force = true)
    res = fig_with(drawer) do
        p2_run(st; out, figures)
    end
    return (; res, out, files = an_files(out), tables = fig_tables(out, P2_SLUG))
end
# the command on the P2 fixture
function fig_cli_p2(tag; drawer)
    st = merge(p2_cached("default"), (; out = an_path("p2fig-cli-$tag")))
    rm(st.out; recursive = true, force = true)
    io = IOBuffer()
    code = fig_with(drawer) do
        AN.cli(p2_cli_args(st); io = io)
    end
    return (; code, text = String(take!(io)), out = st.out, files = an_files(st.out), tables = fig_tables(st.out, P2_SLUG))
end
function fig_p1(tag; drawer, figures = true)
    st = an_cached("default"); out = an_path("p1fig-$tag")
    rm(out; recursive = true, force = true)
    res = fig_with(drawer) do
        an_run(st; out, figures)
    end
    return (; res, out, files = an_files(out), tables = fig_tables(out, P1_SLUG))
end
function fig_cli_p1(tag; drawer)
    st = merge(an_cached("default"), (; out = an_path("p1fig-cli-$tag")))
    rm(st.out; recursive = true, force = true)
    io = IOBuffer()
    code = fig_with(drawer) do
        AN.cli(an_cli_args(st); io = io)
    end
    return (; code, text = String(take!(io)), out = st.out, files = an_files(st.out))
end

const FIG_CACHE = Dict{String, Any}()
fig_get(f, key) = get!(f, FIG_CACHE, key)
function fig_run_p2(tag; fail = nothing, drawer = :fake, figures = true)
    fig_get("p2-$tag-$figures") do
        d, calls = drawer === :fake ? fig_fake(; fail) : (nothing, Any[])
        r = fig_p2("$tag-$figures"; drawer = d, figures)
        merge(r, (; calls))
    end
end
function fig_run_cli_p2(tag; fail = nothing, drawer = :fake)
    fig_get("cli-p2-$tag") do
        d, calls = drawer === :fake ? fig_fake(; fail) : (nothing, Any[])
        merge(fig_cli_p2(tag; drawer = d), (; calls))
    end
end
function fig_run_p1(tag; fail = nothing, drawer = :fake, figures = true)
    fig_get("p1-$tag-$figures") do
        d, calls = drawer === :fake ? fig_fake(; fail) : (nothing, Any[])
        r = fig_p1("$tag-$figures"; drawer = d, figures)
        merge(r, (; calls))
    end
end
function fig_run_cli_p1(tag; fail = nothing, drawer = :fake)
    fig_get("cli-p1-$tag") do
        d, calls = drawer === :fake ? fig_fake(; fail) : (nothing, Any[])
        merge(fig_cli_p1(tag; drawer = d), (; calls))
    end
end

# ---------------------------------------------------------------- tests
@testset "D91b I22 O11: the drawer is called once, after the tables file, with F1 to F6 equal to the items; the digests agree (drawer, result, tables file, reference)" begin
    r = fig_run_p2("ok")
    @test length(r.calls) == 1
    c = isempty(r.calls) ? nothing : r.calls[1]
    @test c !== nothing && c.study == P2_SLUG && abspath(c.out) == abspath(r.out)
    @test c !== nothing && c.tables !== nothing && c.tables == r.tables && !isempty(r.tables)               # the tables file was complete when the drawer was called (E45)
    snap = c === nothing ? nothing : c.snap
    @test snap !== nothing && Set(keys(snap)) == Set("F$i" for i in 1:6)
    text = String(copy(r.tables))
    for i in 1:6
        k = "F$i"
        t = getproperty(r.res.items, Symbol(k))
        s = snap === nothing ? nothing : get(snap, k, nothing)
        @test s !== nothing && s.columns == String.(t.columns) && isequal(s.rows, [collect(x) for x in t.rows])      # the same columns and rows as the item
        @test s !== nothing && s.digest == r.res.digests[k] && s.digest_ref == r.res.digests[k]                       # the digest of the table the drawer received
        j = findfirst("## $k data", text)
        m = j === nothing ? nothing : match(r"digest: ([0-9a-f]{64})", text[j.start:end])
        @test m !== nothing && s !== nothing && m.captures[1] == s.digest                                          # the digest line of its block (E37)
    end
end

@testset "D91d E45: status ok and exit code 0 with the figures written line when the drawer returns" begin
    r = fig_run_p2("ok")
    @test r.res.status === :ok && r.res.study == P2_SLUG
    c = fig_run_cli_p2("ok")
    @test length(c.calls) == 1
    @test c.code == 0
    ls = fig_lines(c.text)
    @test !isempty(ls) && startswith(ls[end], "figures written in") && !occursin("figures skipped", c.text) && !occursin("figures failed", c.text)
end

@testset "D91o I26: a drawer that writes nothing leaves the output folder with the tables file only" begin
    r = fig_run_p2("ok")
    @test r.files == ["p2-indicators-tables.md"]
    c = fig_run_cli_p2("ok")
    @test c.files == ["p2-indicators-tables.md"]
    @test isempty(filter(endswith(".pdf"), c.files))
end

@testset "D91c I23: the tables file with a drawer equals, byte for byte, the file of a run without figures; the drawer is not called without figures" begin
    nofig = fig_run_p2("nofig"; figures = false)
    @test nofig.res.status === :ok && isempty(nofig.calls) && !isempty(nofig.tables)
    ok = fig_run_p2("ok")
    @test ok.tables == nofig.tables && ok.res.digests == nofig.res.digests
    @test fig_run_cli_p2("ok").tables == nofig.tables
end

@testset "D91d E45: a drawer that raises gives status figures_failed, the tables intact, a line that names the cause, exit code 2" begin
    nofig = fig_run_p2("nofig"; figures = false)
    r = fig_run_p2("fail"; fail = FIG_MSG)
    @test length(r.calls) == 1                                                                              # the drawer is reached
    @test r.res.status === :figures_failed
    @test r.files == ["p2-indicators-tables.md"] && r.tables == nofig.tables                                  # the tables stay, no figure file
    c = fig_run_cli_p2("fail"; fail = FIG_MSG)
    @test c.code == 2
    ls = filter(l -> startswith(l, "figures failed"), fig_lines(c.text))
    @test length(ls) == 1 && occursin(FIG_MSG, ls[1])                                                       # names the figure and the row or the cause
    @test !occursin("see the warning above", c.text) && !occursin("figures skipped", c.text)
    @test c.files == ["p2-indicators-tables.md"] && c.tables == nofig.tables
end

@testset "D91d E45: no drawer loaded: figures skipped, the missing plotting package named, not built yet gone, exit code 2, no PDF" begin
    nofig = fig_run_p2("nofig"; figures = false)
    r = fig_run_p2("none"; drawer = :none)
    @test AN.FIGURE_DRAWER[] === nothing                                                                    # the previous value was restored
    @test r.res.status === :figures_skipped && r.files == ["p2-indicators-tables.md"] && r.tables == nofig.tables
    c = fig_run_cli_p2("none"; drawer = :none)
    @test c.code == 2
    ls = filter(l -> startswith(l, "figures skipped"), fig_lines(c.text))
    @test length(ls) == 1
    l = isempty(ls) ? "" : ls[1]
    @test occursin("plotting package", l) || occursin("CairoMakie", l)
    @test occursin("verify", l)                                                                             # the environment that has the package
    @test !occursin("not built yet", c.text)
    @test isempty(filter(endswith(".pdf"), c.files)) && c.files == ["p2-indicators-tables.md"]
    @test fig_run_p2("nofig-none"; drawer = :none, figures = false).res.status === :ok                      # without figures wanted nothing is skipped
end

@testset "D91e I24: P1 with a drawer: F1 to F7 handed over, tables and digests equal to the run without figures" begin
    nofig = fig_run_p1("nofig"; figures = false)
    @test nofig.res.status === :ok && isempty(nofig.calls)
    r = fig_run_p1("ok")
    @test r.res.status === :ok && length(r.calls) == 1
    c = isempty(r.calls) ? nothing : r.calls[1]
    @test c !== nothing && c.study == P1_SLUG && c.snap !== nothing && Set(keys(c.snap)) == Set("F$i" for i in 1:7)      # exactly F1 to F7, no P2 table
    @test c !== nothing && c.tables == r.tables
    @test r.tables == nofig.tables && !isempty(r.tables) && r.res.digests == nofig.res.digests
    @test r.files == ["p1-le-benchmark-tables.md"]
    ok = fig_run_cli_p1("ok")
    @test ok.code == 0 && startswith(fig_lines(ok.text)[end], "figures written in")
end

@testset "D91e I24: P1 messages and exit codes unchanged: failing drawer, no drawer" begin
    f = fig_run_p1("fail"; fail = FIG_MSG)
    @test f.res.status === :figures_failed && length(f.calls) == 1
    @test f.tables == fig_run_p1("nofig"; figures = false).tables
    c = fig_run_cli_p1("fail"; fail = FIG_MSG)
    @test c.code == 2 && FIG_P1_FAILED in fig_lines(c.text)                                                 # the old line, word for word
    s = fig_run_cli_p1("none"; drawer = :none)
    @test s.code == 2 && FIG_P1_SKIPPED in fig_lines(s.text) && !occursin("p2-indicators", s.text)
    @test fig_run_p1("none"; drawer = :none).res.status === :figures_skipped
end

@testset "D91n I29: the analysis layer loads no plotting package, the figures layer imports the analysis layer and the plotting package only" begin
    root = joinpath(@__DIR__, "..", "..")
    ana = filter(f -> endswith(f, ".jl"), readdir(joinpath(root, "src", "Analysis"); join = true))
    @test !isempty(ana)
    code(f) = [(basename(f), strip(first(split(l, '#')))) for l in eachline(f)]                              # comments stripped
    pkgs = r"\b(?:\w*Makie\w*|Plots|PyPlot|Gadfly|PlotlyJS|UnicodePlots|VegaLite|Gnuplot)\b"
    @test isempty([x for f in ana for x in code(f) if occursin(pkgs, x[2])])
    @test isempty([x for f in ana for x in code(f) if occursin(r"\bBase\.require\b|@eval\s+(?:using|import)\b|\bRequires\b|\bPkg\.(?:add|develop)\b", x[2])])
    @test isempty([x for f in ana for x in code(f) if occursin(r"\bFIGURE_DRAWER\[\]\s*=[^=]", x[2])])         # only the extension sets the drawer
    ext = filter(f -> occursin(r"CairoMakie", basename(f)) || occursin(r"CairoMakie", basename(dirname(f))),
                 [joinpath(d, f) for (d, _, fs) in walkdir(joinpath(root, "ext")) for f in fs if endswith(f, ".jl")])
    @test !isempty(ext)
    imp = [(basename(f), l) for f in ext for (_, l) in code(f) if occursin(r"^(?:using|import)\s", l)]
    tops(l) = [String(match(r"\w+", strip(t)).match) for t in split(first(split(replace(l, r"^(?:using|import)\s+" => ""), ':')), ',') if !isempty(strip(t))]
    @test !isempty(imp) && all(x -> all(n -> n in ("Chaoslab", "CairoMakie"), tops(x[2])), imp)
    @test any(x -> occursin("CairoMakie", x[2]), imp) && any(x -> occursin("Chaoslab", x[2]), imp)
end

@testset "D91n I29: the dependency lists are unchanged: weak dependencies CairoMakie and CUDA, the two extensions, the CairoMakie bound 0.15" begin
    proj = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))
    @test Set(keys(proj["weakdeps"])) == Set(["CairoMakie", "CUDA"])
    @test proj["extensions"] == Dict("ChaoslabCUDAExt" => "CUDA", "ChaoslabCairoMakieExt" => "CairoMakie")
    @test proj["compat"]["CairoMakie"] == "0.15" && proj["compat"]["CUDA"] == "6"
    @test !haskey(proj["deps"], "CairoMakie") && !haskey(proj["deps"], "CUDA")                               # weak, not hard
    @test Set(keys(proj["deps"])) == Set(["Dates", "DoubleFloats", "HDF5", "JSON", "KernelAbstractions", "LinearAlgebra", "PrecompileTools", "SHA", "StableRNGs", "StaticArrays", "TOML"])
    @test isfile(joinpath(@__DIR__, "..", "..", "ext", "ChaoslabCairoMakieExt.jl"))
end
