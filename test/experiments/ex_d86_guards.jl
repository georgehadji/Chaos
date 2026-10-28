# Amendment D86, criteria D86w (no study value, no threshold, no detection; the (S) part) and D86x (layering and entry point).
# Spec: experiments-spec.md, "Amendment D86": R11 item 6, R21, I25, I27. The (V) part of D86w (a script in verify/ over the p2 values) is stage 04.

# A second fixture block whose specs differ from the first in every value: other kinds subset, chains, windows, bins, orders, observables,
# separation, budgets, tau (D86w). Expected outcomes are computed from the block.
function d86_alt_block()
    bud = Dict{String, Any}("Float64" => Dict{String, Any}("ns_cpu_thread" => 55.0, "safety" => 4))
    specs = [d86_spec("tangent"; id = "ta", members = ["MEGNO", "SALI"], budget = bud),
             d86_spec("pair"; id = "pb", members = ["MEGNO", "RLI"], separation = 1e-7, frame_source = "direction", budget = bud),
             d86_spec("stretching"; id = "sb", range = [-3.0, 3.0], bins = 33, interval = 2, min_per_cell = 1, counts = "S_D", budget = bud),
             d86_spec("frequency"; id = "fb", observable = "q_plus_i_p", windows = [24, 48], order = 2, positions = "edge", pad = 3, mean_removal = false,
                      mirror_ratio = 0.8, floor_multiplier = 4.0, budget = bud),
             d86_spec("birkhoff"; id = "bb", observables = [Dict{String, Any}("kind" => "sin", "coordinate" => 2, "frequency" => 3), Dict{String, Any}("kind" => "lifted", "coordinate" => 1)],
                      windows = [20, 40], birkhoff_order = 3, budget = bud)]
    arms = [d86_arm("ta", "S_D", ["ta", "sb", "pb"], ["Float64"], [90]; tau = [2], coordinate_mode = ["wrapped"], direction = ["0.3"]),
            d86_arm("fa", "S_D", ["fb"], ["Float64"], [192]; tau = [2], coordinate_mode = ["wrapped"]),
            d86_arm("ba", "S_D", ["bb"], ["Float64"], [160]; tau = [2], coordinate_mode = ["wrapped"])]
    return d86_block(arms, specs)
end

# the checks of the D86 suite that need no hand numbers, for any fixture run: the items equal the independent enumeration, and the stored values equal the stand-alone
# reducers and the pure functions on the independent series (O16, O17)
function d86_generic_check(run)
    plan = run.plan
    bad = String[]
    ref = d86_ref_keys(run.fx.run, run.fx.orbits)
    got = Dict((it.orbit, it.method, ex_label(it.precision), it.horizon, it.tau, it.direction, it.coordinate_mode) => it.arms for it in plan.items if it.method in keys(ref.perspec))
    got == ref.items || push!(bad, "items")
    (plan.not_applicable == ref.omitted) || push!(bad, "omitted")
    d86_each(run) do it, o, raw, c, file
        it.method in keys(ref.perspec) || return
        s = d86_spec_of(plan, it)
        if s["kind"] == "tangent"
            ser = d86_ds(raw, o, it.method, "series"); cols = d86_at(raw, o, it.method, "series", "columns"); n = raw["orbits/$o/t"][end]
            for m in s["members"], (name, v) in d86_member_cols(m, d86_stand_alone(it, m, n))
                reinterpret(UInt64, ser[end, findfirst(==(name), cols)]) == reinterpret(UInt64, v) || push!(bad, "tangent $name orbit $o")
            end
        elseif s["kind"] == "stretching"
            n = raw["orbits/$o/t"][end]
            Float64(d86_stretch_final(it, s, n).fin.D) == d86_ds(raw, o, it.method, "series")[end, 1] || push!(bad, "stretching orbit $o")
        elseif s["kind"] == "pair"
            ser = d86_ds(raw, o, it.method, "series"); cols = d86_at(raw, o, it.method, "series", "columns"); n = raw["orbits/$o/t"][end]
            for m in s["members"], (name, v) in d86_member_cols(m, d86_pair_ref(it, s, raw["orbits/$o@xi0"], n, m).fin)
                reinterpret(UInt64, ser[end, findfirst(==(name), cols)]) == reinterpret(UInt64, v) || push!(bad, "pair $name orbit $o")
            end
        elseif s["kind"] == "frequency"
            fin = d86_ds(raw, o, it.method, "final")
            for (r, row) in enumerate(d86_ref_freq(it, s))
                d86_isequal(fin[r, :], row.vals) || push!(bad, "frequency orbit $o")
            end
        elseif s["kind"] == "phase_advance"
            fin = d86_ds(raw, o, it.method, "final")
            for (r, row) in enumerate(d86_ref_pa(it, s))
                d86_isequal(fin[r, :], row.vals) || push!(bad, "phase advance orbit $o")
            end
        else
            fin = d86_ds(raw, o, it.method, "final")
            for (r, row) in enumerate(d86_ref_bk(it, s).rows)
                d86_isequal(fin[r, :], row) || push!(bad, "birkhoff orbit $o")
            end
        end
    end
    return bad
end

@testset "D86w the suite runs with two fixture blocks whose specs differ in every value; no threshold, no class, no detection time is stored" begin
    A = d86_all_ind_run()
    B = d86_run("w-alt", d86_alt_block())
    @test A.res.status == "complete" && B.res.status == "complete"
    @test isempty(d86_generic_check(A))
    @test isempty(d86_generic_check(B))
    # the two blocks share no spec value: the spec hashes of the same kind differ
    ha = Dict(s["kind"] => d86_spec_hash(s) for s in A.fx.run["indicator"])
    hb = Dict(s["kind"] => d86_spec_hash(s) for s in B.fx.run["indicator"])
    @test all(k -> ha[k] != hb[k], intersect(keys(ha), keys(hb))) && length(intersect(keys(ha), keys(hb))) >= 4
    # no dataset or attribute of any kind holds a threshold or a first-crossing step (the residual thresholds of the monitor, R5, are the module's own)
    words = r"threshold|class|detect|crossing|first_cross|sticky|regular|chaotic"i
    bad = String[]
    for run in (A, B), f in ex_chunk_files(run.store)
        raw = d86_raw(f)
        for (k, v) in raw
            occursin("diagnostics/", k) && continue
            if occursin("/indicators/", k) || startswith(k, "pooled/") || k in ("@kind", "@spec", "@spec_sha256", "@chain", "@members", "@observable_names", "@frame_source", "@direction_rule", "@separation", "@windows", "@run_length")
                occursin(words, k) && push!(bad, k)
                v isa AbstractVector{<:AbstractString} && any(x -> occursin(words, x), v) && push!(bad, k)
            end
        end
    end
    @test isempty(bad)
    # the grammar accepts no key that could hold one
    for key in ("threshold", "class", "detect", "first_crossing", "regular_level", "sticky")
        for (kind, base) in (("tangent", d86_spec("tangent")), ("frequency", d86_spec("frequency")))
            blk = d86_one(base; mode = kind == "frequency" ? "wrapped" : nothing)
            blk["indicator"][1][key] = 0.5
            @test d86_rej(() -> EX.run_block(blk), "indicator[1].$key")
        end
    end
end

@testset "D86x layering and entry point" begin
    src = joinpath(@__DIR__, "..", "..", "src")
    dir = joinpath(src, "Experiments")
    allowed = ("Foundation", "Systems", "Integrators", "Tangent", "Lyapunov", "Diagnostics", "Indicators")
    deps = keys(TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["deps"])
    bad = String[]
    imports_indicators = false
    for f in filter(endswith(".jl"), readdir(dir)), line in eachline(joinpath(dir, f))
        occursin(r"^\s*(?:using|import)\s", line) || continue
        occursin(r"\.\.Indicators\b", line) && (imports_indicators = true)
        for m in eachmatch(r"\.\.(\w+)", line)
            m.captures[1] in allowed || push!(bad, "$f: $line")
        end
        for part in split(replace(line, r"^\s*(?:using|import)\s+" => ""), ',')
            head = strip(first(split(strip(part), ':')))
            (isempty(head) || startswith(head, ".")) && continue
            String(first(split(head, '.'))) in deps || push!(bad, "$f: $line")
        end
    end
    @test isempty(bad)
    @test imports_indicators                                           # Indicators is among the declared imports of Experiments
    # nothing below imports Experiments (Analysis sits above it)
    upward = String[]
    for (d, _, fs) in walkdir(src), f in fs
        endswith(f, ".jl") || continue
        owner = first(splitpath(relpath(d, src)))
        owner in ("Experiments", ".", "Analysis") && continue
        for line in eachline(joinpath(d, f))
            occursin(r"^\s*(?:using|import)\s", line) && occursin("Experiments", line) && push!(upward, "$owner/$f")
        end
    end
    @test isempty(upward)
    # the commands `run` and `pilot` gain no option and the entry point contains no numerics
    uio = IOBuffer()
    EX.cli(String[]; io = uio)                                          # no command: the usage text
    usage = String(take!(uio))
    opts = sort(unique(String[m.match for m in eachmatch(r"--[a-z-]+", usage)]))
    @test opts == sort(["--out", "--store", "--manifest", "--threads", "--thread-override", "--allow-dirty", "--repo", "--code-version", "--capacity", "--cap-override", "--stale",
                        "--retry-failed", "--retry-control-stop", "--backend-override", "--manifest-reset", "--prior-record"])
    cli = read(joinpath(dir, "run_cli.jl"), String)
    @test !occursin(r"SALI|FLI|MEGNO|RLI|Stretching|Birkhoff|PhaseAdvance|frequency_map|Indicators|indicator", replace(cli, r"#[^\n]*" => ""))
    # every stop status is nonzero
    statuses = ["stopped: control bound", "stopped: failed chunks", "stopped: wall cap", "paused arms", "complete with failed chunks"]
    @test all(s -> EX.exit_status(s) != 0, statuses) && EX.exit_status("complete") == 0
    # a start refused by the checks before the stream of D86 gives exit status 1 through the entry point
    fx = d86_fx("x-refuse", ex_rb_with(d86_one(d86_spec("tangent")), ("indicator", 1, "members") => ["SALI", "FLI"]))
    io = IOBuffer()
    @test EX.cli(["run", fx.prereg, fx.panel, "--store", ex_store("d86-x-refuse"), "--manifest", ex_path("d86-x-refuse.json"), "--thread-override", "--allow-dirty", "--code-version", EX_CODE], io = io) == 1
    @test occursin("indicator[1].members", String(take!(io)))
end
