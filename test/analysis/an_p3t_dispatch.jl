# D95a dispatch and D95p layering. Spec: analysis-spec.md, Amendment D95, E59, I44, criteria D95a, D95p.
# The unchanged behaviour of `calibrate` and of the p1 and p2 pipelines (I44) is guarded by the existing golden-value testsets "D93a I36 the P1 and P2 pipelines are unchanged" (an_p3_dispatch.jl),
# the D93 pipeline tests of the record and the tables, and the P1 guard of an_p2_dispatch.jl; they are not repeated. The existing testset "D93a analyze with the preregistration of this study is
# refused" is replaced by this amendment (E59 lifts that refusal): it must be reworked at stage 03 (test plan, "Rework needed at stage 03").

@testset "D95a count-round is refused for any other study, naming the slug; the unknown-slug message of analyze is unchanged" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    inst = p3t_instance(S, I, "disp"; rounds = 1, closure = false)
    # the p2 and p1 preregistrations and a preregistration with an unknown slug
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst; prereg = st.s2.prereg)...), ["p2-indicators"]) && !isfile(inst.log)
    g = p1_gold()
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst; prereg = g.prereg)...), ["p1-le-benchmark"]) && !isfile(inst.log)
    text, _ = p3t_prereg_text(I.q, I.t; slug = "p3-other")
    pre = joinpath(inst.dir, "p3-other-preregistration.md"); write(pre, text)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst; prereg = pre)...), ["p3-other"]) && !isfile(inst.log)
    # analyze with the p3 preregistration runs the test-side pipeline (no refusal naming calibrate): the closed log is required, so an absent log is the refusal here
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, inst)...), ["round log"]) && !ispath(inst.out)
end

@testset "D95a command lines: usage and exit codes, --help, the p3 form of analyze, --out and --round-log inside data/" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    inst = p3t_instance(S, I, "disp-cli"; rounds = 1, closure = false)
    io = IOBuffer(); @test AN.cli(["count-round"]; io = io) == 1 && occursin("count-round", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["count-round", "--help"]; io = io) == 0 && occursin("count-round <preregistration> <parameters> <test-parameters>", String(take!(io)))
    io = IOBuffer(); AN.cli(["count-round", "--help"]; io = io); h = String(take!(io))
    @test all(o -> occursin(o, h), ("--store", "--panel", "--round-log", "--p2-manifest", "--p2-preregistration", "--p2-parameters", "--calibration-record", "--threads"))
    io = IOBuffer(); @test AN.cli(["count-round", I.prereg]; io = io) == 1 && occursin("usage", lowercase(String(take!(io))))
    io = IOBuffer(); @test AN.cli(["count-round", I.prereg, I.params, I.tparams, "--bogus", "x"]; io = io) == 1 && occursin("--bogus", String(take!(io)))
    io = IOBuffer(); @test AN.cli(["count-round", I.prereg, I.params, I.tparams, "--store"]; io = io) == 1
    io = IOBuffer(); @test AN.cli(["count-round", I.prereg, I.params, I.tparams, "--threads", "many"]; io = io) == 1
    io = IOBuffer(); @test AN.cli(["count-round", I.prereg, I.params, I.tparams, "extra-argument"]; io = io) == 1 && !isfile(inst.log)
    # the command line of count-round: exit code 0, the entry, the counts, the targets and the decision are printed, the file is created
    io = IOBuffer(); @test AN.cli(p3t_cli_cr(I, st, inst); io = io) == 0
    out = String(take!(io))
    @test occursin("round=1", out) && occursin("continue", out) && occursin("count.t_s=2", out) && occursin("target", out) && isfile(inst.log)
    c = p3t_clone(p3t_full(), "disp-an")
    io = IOBuffer(); @test AN.cli(["analyze", I.prereg]; io = io) == 1
    io = IOBuffer(); @test AN.cli(p3t_cli_an(I, st, c; extra = ["--bogus", "x"]); io = io) == 1 && occursin("--bogus", String(take!(io)))
    # --out and --round-log inside data/ refuse, nothing is created there
    bad = joinpath(@__DIR__, "..", "..", "data", "an-p3t-refused-cli")
    try
        io = IOBuffer(); @test AN.cli(p3t_cli_an(I, st, merge(c, (; out = bad))); io = io) == 1 && !ispath(bad) && occursin("data", String(take!(io)))
        io = IOBuffer(); @test AN.cli(vcat(p3t_cli_cr(I, st, inst), ["--round-log", bad * ".log"]); io = io) == 1 && !ispath(bad * ".log")
    finally
        ispath(bad) && rm(bad; recursive = true, force = true); ispath(bad * ".log") && rm(bad * ".log"; force = true)         # only these paths, and only if a defective module created them
    end
end

@testset "D95p layering: Analysis imports Foundation, Experiments and the standard library only; the D95 entry points exist; no new package" begin
    root = joinpath(@__DIR__, "..", "..")
    ana = filter(f -> endswith(f, ".jl"), readdir(joinpath(root, "src", "Analysis"); join = true))
    lines = [(basename(f), l) for f in ana for l in eachline(f) if occursin(r"^\s*(?:using|import)\s", l)]
    allowed = r"^\s*(?:using|import)\s+(?:\.\.Experiments|\.\.Foundation|SHA|StableRNGs|TOML|LinearAlgebra|Printf|Dates)\b"
    @test isempty([x for x in lines if !occursin(allowed, x[2])])
    @test isempty([l for f in ana for l in eachline(f) if occursin(r"Makie|Plots|PyPlot|Gadfly|PlotlyJS|UnicodePlots|VegaLite|Gnuplot", l) && occursin(r"^\s*(?:using|import)\s", l)])
    proj = TOML.parsefile(joinpath(root, "Project.toml"))
    @test Set(keys(proj["deps"])) == Set(["Dates", "DoubleFloats", "HDF5", "JSON", "KernelAbstractions", "LinearAlgebra", "PrecompileTools", "SHA", "StableRNGs", "StaticArrays", "TOML"])
    # the entry points of the amendment (the assertions of this testset that fail before the amendment is implemented)
    @test AN.count_round isa Function && AN.p3t_read_params isa Function && length(AN.P3T_PARAM_KEYS) == 18 && AN.p3t_ha isa Function && AN.p3t_hb isa Function && AN.p3t_hc isa Function
end
