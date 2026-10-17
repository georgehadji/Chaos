# AC1 input validation, AC2 preregistration hash, plus the exported-API and layering checks (AC28).

const EX_REAL_PREREG = joinpath(@__DIR__, "..", "..", "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")

@testset "API exported through Chaoslab" begin
    names_needed = (:panel_inputs, :read_prereg, :check_prereg, :prereg_hash, :build_panel, :portrait_content, :periodic_set,
                    :collect_orbits, :classify_trace, :orbit_exponent, :newton_correction, :accept_residual, :periodic_distance,
                    :grid_centres, :torus_delta, :torus_distance, :screen_buffer, :screen_point, :screen_point!, :screen_bytes,
                    :block_median, :stuck_length, :pred_T, :pred_C, :pred_S, :cell_stream_seed, :cell_rng, :draw_candidates,
                    :place_X, :place_R, :place_R0, :orbit_seed, :y5_control, :select_subsets, :git_state, :package_version)
    @test Chaoslab.Experiments isa Module
    @test all(n -> n in names(Chaoslab.Experiments), names_needed)
    @test all(n -> isdefined(Chaoslab, n), names_needed)       # re-exported by the top module like the other submodules
end

@testset "AC1 input validation" begin
    # 35 malformed inputs, each rejected with an ArgumentError that names the key; the valid base passes.
    base = ex_block()
    pin = EX.panel_inputs(base)
    @test pin.n_c == 3 && pin.grid_side == 8 && pin.N_scr == 2000 && pin.w == 100
    @test pin.K == [0.0, 1.5, 5.0]
    mk(f) = (d = deepcopy(base); f(d); d)
    cases = [
        ("missing n_c", mk(d -> delete!(d, "n_c")), ["n_c"]),
        ("missing K", mk(d -> delete!(d, "K")), ["K"]),
        ("unknown key", mk(d -> d["bogus"] = 1), ["bogus"]),
        ("K infinite", mk(d -> d["K"] = ["0", "Inf"]), ["K"]),
        ("K nan", mk(d -> d["K"] = ["NaN"]), ["K"]),
        ("K empty", mk(d -> d["K"] = String[]), ["K"]),
        ("K repeated", mk(d -> d["K"] = ["0.5", "0.5"]), ["K"]),
        ("K repeated as binary64", mk(d -> d["K"] = ["pi/2", "1.5707963267948966"]), ["K"]),
        ("grid_side zero", mk(d -> d["grid_side"] = 0), ["grid_side"]),
        ("N_scr negative", mk(d -> d["N_scr"] = -5), ["N_scr"]),
        ("n_c zero", mk(d -> d["n_c"] = 0), ["n_c"]),
        ("q_max not an integer", mk(d -> d["q_max"] = 1.5), ["q_max"]),
        ("newton_max zero", mk(d -> d["newton_max"] = 0), ["newton_max"]),
        ("N_scr not a multiple of w", mk(d -> d["N_scr"] = 2050), ["N_scr", "w"]),
        ("w not a multiple of tau", mk(d -> d["tau"] = 30), ["tau", "w"]),
        ("q_max beyond the supported range", mk(d -> d["q_max"] = 1_000_000), ["q_max"]),
        ("theta_T equals theta_S", mk(d -> d["theta_T"] = d["theta_S"]), ["theta_T", "theta_S"]),
        ("theta_T above theta_C", mk(d -> (d["theta_T"] = 0.06; d["theta_S"] = 0.07)), ["theta_T", "theta_C"]),
        ("r_T zero", mk(d -> d["r_T"] = 0.0), ["r_T"]),
        ("r_C negative", mk(d -> d["r_C"] = -1.0), ["r_C"]),
        ("delta_R zero", mk(d -> d["delta_R"] = 0.0), ["delta_R"]),
        ("delta_dup zero", mk(d -> d["delta_dup"] = 0.0), ["delta_dup"]),
        ("rho zero", mk(d -> d["rho"] = 0.0), ["rho"]),
        ("rho one", mk(d -> d["rho"] = 1.0), ["rho"]),
        ("f_S zero", mk(d -> d["f_S"] = 0.0), ["f_S"]),
        ("f_S one", mk(d -> d["f_S"] = 1.0), ["f_S"]),
        ("L_S below w", mk(d -> d["L_S"] = 50), ["L_S"]),
        ("L_S above N_scr", mk(d -> d["L_S"] = 3000), ["L_S"]),
        ("l_lo above l_hi", mk(d -> d["log10_r"] = [-2.0, -4.0]), ["log10_r"]),
        ("10^l_hi at least pi", mk(d -> d["log10_r"] = [-4.0, 0.6]), ["log10_r"]),
        ("mu one", mk(d -> d["mu"] = "1"), ["mu"]),
        ("mu nan", mk(d -> d["mu"] = "NaN"), ["mu"]),
        ("s_L above s_D", mk(d -> d["s_L"] = 3), ["s_L", "s_D"]),
        ("cap_D below the number of controls", mk(d -> d["cap_D"] = 1), ["cap_D"]),
        ("alpha nan", mk(d -> d["alpha"] = "NaN"), ["alpha"]),
        ("sn_trace is no longer an input (D51)", mk(d -> d["sn_trace"] = 1e-6), ["sn_trace"]),
    ]
    @test length(cases) >= 24
    for (name, d, keys) in cases
        @testset "$name" begin
            @test ex_rejects(() -> EX.panel_inputs(d), keys)
        end
    end

    @testset "rejected before any computation, no file written" begin
        for (tag, over, drop) in (("bad-LS", Dict{String, Any}("L_S" => 50), String[]),
                                  ("bad-missing", Dict{String, Any}(), ["n_c"]),
                                  ("bad-unknown", Dict{String, Any}("bogus" => 1), String[]))
            pre = ex_write_prereg(ex_path("$tag-prereg.md"), ex_block(over; drop = drop))
            out = ex_path("$tag-panel.toml")
            @test_throws ArgumentError EX.build_panel(pre; out = out)
            @test !isfile(out) && !isfile(out * ".provenance.toml")
        end
    end
end

@testset "AC2 preregistration hash" begin
    # (a) the frozen study file: header hash equals the recomputed one and the independent one
    @testset "frozen p1-le-benchmark" begin
        text = read(EX_REAL_PREREG, String)
        r = EX.check_prereg(EX_REAL_PREREG)
        @test r.header_sha256 == r.computed_sha256 == ex_header_hash(text) == ex_sha256(ex_prepare(text))
        @test r.amendments == length(collect(eachmatch(r"\bA\d+ \(", match(r"(?m)^\| Amendment log \|.*$", text).match)))
        # the machine-readable block (Amendment A1, item 7) is read as written: same keys, same integer values, K as nearest binary64
        blk = TOML.parse(match(r"```toml\n(.*?)\n```"s, text).captures[1])["panel"]
        pin = EX.read_prereg(EX_REAL_PREREG).inputs
        @test pin.n_c == blk["n_c"] && pin.grid_side == blk["grid_side"] && pin.N_scr == blk["N_scr"] && pin.w == blk["w"]
        @test pin.K == [ex_f64(s) for s in blk["K"]]
        @test pin.mu == ex_f64(blk["mu"]) && pin.alpha == ex_f64(blk["alpha"])
        @test pin.theta_T == blk["theta_T"] && pin.rho == blk["rho"]
    end

    pre = ex_write_prereg(ex_path("hash-prereg.md"))
    text = read(pre, String)
    h = ex_header_hash(text)

    @testset "fixture passes; module value equals the independent one (O3)" begin
        r = EX.check_prereg(pre)
        @test r.header_sha256 == r.computed_sha256 == h
        @test r.amendments == 2
        @test EX.prereg_hash(text) == ex_sha256(ex_prepare(text)) == h
        @test EX.read_prereg(pre).inputs.n_c == 3
    end

    @testset "(b) one changed byte at 20 seeded positions refuses" begin
        rng = StableRNG(ex_seed(2, 1))
        bytes = Vector{UInt8}(codeunits(text))
        cell = findfirst("SHA-256 of this file (current version) | ", text)
        digits0 = last(cell) + 1                                  # first digit of the 64-digit value
        note0 = digits0 + 64                                      # the parenthetical is replaced by SELF before hashing: excluded
        slug = first(findfirst("ex-fixture |", text))
        alnum = [i for i in eachindex(bytes) if (UInt8('0') <= bytes[i] <= UInt8('9') || UInt8('a') <= bytes[i] <= UInt8('z') ||
                                                 UInt8('A') <= bytes[i] <= UInt8('Z')) && !(note0 <= i < note0 + 20)]
        positions = vcat([slug, digits0, digits0 + 63], alnum[1 .+ floor.(Int, rand(rng, 17) .* length(alnum))])
        @assert length(positions) == 20
        for (n, i) in enumerate(positions)
            b = copy(bytes)
            b[i] = b[i] == UInt8('x') ? UInt8('y') : UInt8('x')
            f = ex_path("hash-mut$n.md")
            write(f, b)
            @test_throws ArgumentError EX.check_prereg(f)
        end
    end

    @testset "(c) missing row, two rows, value that is not 64 hex digits" begin
        norow = replace(text, r"(?m)^\| SHA-256 of this file.*\n" => "")
        write(ex_path("hash-norow.md"), norow)
        @test_throws ArgumentError EX.check_prereg(ex_path("hash-norow.md"))
        row = match(r"(?m)^\| SHA-256 of this file.*\n", text).match
        write(ex_path("hash-tworows.md"), replace(text, row => row * row))
        @test_throws ArgumentError EX.check_prereg(ex_path("hash-tworows.md"))
        write(ex_path("hash-short.md"), replace(text, h => h[1:63]))
        @test_throws ArgumentError EX.check_prereg(ex_path("hash-short.md"))
        write(ex_path("hash-nonhex.md"), replace(text, h => "g" * h[2:end]))
        @test_throws ArgumentError EX.check_prereg(ex_path("hash-nonhex.md"))
    end

    @testset "(d) CRLF line endings pass" begin
        write(ex_path("hash-crlf.md"), replace(text, "\n" => "\r\n"))
        r = EX.check_prereg(ex_path("hash-crlf.md"))
        @test r.header_sha256 == r.computed_sha256 == h
    end
end

@testset "AC28 layering of src/Experiments" begin
    src = joinpath(@__DIR__, "..", "..", "src")
    dir = joinpath(src, "Experiments")
    allowed_internal = ("Foundation", "Systems", "Integrators", "Tangent", "Lyapunov", "Diagnostics")
    deps = keys(TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["deps"])
    files = isdir(dir) ? filter(endswith(".jl"), readdir(dir)) : String[]
    bad = String[]
    for f in files, line in eachline(joinpath(dir, f))
        occursin(r"^\s*(?:using|import)\s", line) || continue
        for m in eachmatch(r"\.\.(\w+)", line)
            m.captures[1] in allowed_internal || push!(bad, "$f: $line")
        end
        for part in split(replace(line, r"^\s*(?:using|import)\s+" => ""), ',')
            head = strip(first(split(strip(part), ':')))
            (isempty(head) || startswith(head, ".")) && continue
            String(first(split(head, '.'))) in deps || push!(bad, "$f: $line")      # declared dependencies only (standard library packages must be declared too)
        end
    end
    @test !isempty(files) && isempty(bad)
    # nothing below Experiments imports it
    upward = String[]
    for (d, _, fs) in walkdir(src), f in fs
        endswith(f, ".jl") || continue
        owner = first(splitpath(relpath(d, src)))
        owner in ("Experiments", ".") && continue
        for line in eachline(joinpath(d, f))
            occursin(r"^\s*(?:using|import)\s", line) && occursin("Experiments", line) && push!(upward, "$owner/$f")
        end
    end
    @test isdir(dir) && isempty(upward)
end
