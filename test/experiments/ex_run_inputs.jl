# D53a block validation and refusals, D53r (no study value in the code, two fixture blocks). Spec: Amendment D53, R1 item 1.
# Block-level checks go through EX.run_block(dict); file-level refusals through plan_run and run_campaign. All refuse before any
# computation and write nothing.

# the run table of the real preregistration (Amendment A4): the second ```toml block
function ex_real_run()
    text = replace(read(EX_REAL_PREREG, String), "\r\n" => "\n")
    for m in eachmatch(r"```toml\n(.*?)\n```"s, text)
        d = TOML.parse(String(m.captures[1]))
        haskey(d, "run") && return d["run"]
    end
    error("no run block in $EX_REAL_PREREG")
end

@testset "API exported through Chaoslab (D53)" begin
    names_needed = (:run_block, :read_run_block, :block_text, :block_hash, :plan_run, :plan_text, :chunk_text, :item_text, :item_hash,
                    :record_grid, :grid_hash, :deviation_frame, :launch_capacity, :launch_groups, :resolve_run_backend, :backend_table,
                    :cpu_only_method, :code_version, :run_campaign, :RunHooks, :read_chunk, :write_chunk, :chunk_layout, :read_json,
                    :run_pilot, :rerun_check, :exit_status, :cli, :default_store, :hot_loop_bytes)
    @test all(n -> n in names(Chaoslab.Experiments), names_needed)
    @test all(n -> isdefined(Chaoslab, n), names_needed)
end

@testset "D53a block validation" begin
    base = ex_rb(:small)
    rb = EX.run_block(base)                                           # the valid block passes
    @test EX.block_hash(rb) == ex_sha256(EX.block_text(rb))           # canonical text, hashed by the platform tool
    mk(edits...) = ex_rb_with(base, edits...)
    D = :delete
    cases = [
        ("missing primary", mk(("primary",) => D), ["primary"]),
        ("missing arm", mk(("arm",) => D), ["arm"]),
        ("missing grid", mk(("grid",) => D), ["grid"]),
        ("missing backend", mk(("backend",) => D), ["backend"]),
        ("missing health", mk(("health",) => D), ["health"]),
        ("missing budget", mk(("budget",) => D), ["budget"]),
        ("missing k19", mk(("k19",) => D), ["k19"]),
        ("missing verify", mk(("verify",) => D), ["verify"]),
        ("missing pilot", mk(("pilot",) => D), ["pilot"]),
        ("missing accepted list", mk(("panel_prereg_sha256_accepted",) => D), ["panel_prereg_sha256_accepted"]),
        ("missing control_bound", mk(("health", "control_bound") => D), ["control_bound"]),
        ("missing setup_check_tol", mk(("health", "setup_check_tol") => D), ["setup_check_tol"]),
        ("missing points_per_decade", mk(("grid", "points_per_decade") => D), ["points_per_decade"]),
        ("missing direction", mk(("primary", "direction") => D), ["direction"]),
        ("missing safety", mk(("budget", "Float64", "safety") => D), ["safety"]),
        ("missing rerun_sample", mk(("verify", "rerun_sample") => D), ["rerun_sample"]),
        ("unknown top key", mk(("bogus",) => 1), ["bogus"]),
        ("unknown health key", mk(("health", "bogus_h") => 1), ["bogus_h"]),
        ("unknown arm key", mk(("arm", 1, "bogus_a") => 1), ["bogus_a"]),
        ("unknown primary key", mk(("primary", "bogus_p") => 1), ["bogus_p"]),
        ("unknown method", mk(("arm", 1, "methods") => ["L1", "L7"]), ["methods"]),
        ("unknown precision", mk(("arm", 1, "precisions") => ["Float32"]), ["precisions"]),
        ("unknown subset", mk(("arm", 1, "orbits") => "S_X"), ["orbits"]),
        ("unknown backend request", mk(("backend", "Float64") => "tpu"), ["backend"]),
        ("unknown control method", mk(("health", "control_methods") => ["L1", "L77"]), ["control_methods"]),
        ("unknown coordinate mode", mk(("arm", 2, "coordinate_mode") => ["polar"]), ["coordinate_mode"]),
        ("unparsable direction", mk(("arm", 2, "direction") => ["banana"]), ["direction"]),
        ("empty methods", mk(("arm", 1, "methods") => String[]), ["methods"]),
        ("empty horizons", mk(("arm", 1, "horizons") => Int[]), ["horizons"]),
        ("empty precisions", mk(("arm", 1, "precisions") => String[]), ["precisions"]),
        ("h4 subset without h4 starts", mk(("h4",) => D), ["h4"]),
        ("duplicate arm name", mk(("arm", 2, "name") => "main"), ["name"]),
        ("horizon not a multiple of tau", mk(("arm", 2, "horizons") => [101]), ["horizons", "tau"]),
        ("N / 10 not a multiple of tau", mk(("arm", 2, "tau") => [1, 8], ("arm", 2, "horizons") => [40]), ["extra_fractions"]),
        ("extra fraction one", mk(("grid", "extra_fractions") => [0.1, 1.0]), ["extra_fractions"]),
        ("extra fraction zero", mk(("grid", "extra_fractions") => [0.0]), ["extra_fractions"]),
        ("extra fraction negative", mk(("grid", "extra_fractions") => [-0.1]), ["extra_fractions"]),
        ("extra fraction above one", mk(("grid", "extra_fractions") => [0.5, 1.5]), ["extra_fractions"]),
        ("points_per_decade zero", mk(("grid", "points_per_decade") => 0), ["points_per_decade"]),
        ("points_per_decade not an integer", mk(("grid", "points_per_decade") => 2.5), ["points_per_decade"]),
        ("first_steps zero", mk(("grid", "first_steps") => 0), ["first_steps"]),
        ("first_steps not an integer", mk(("grid", "first_steps") => 1.5), ["first_steps"]),
        ("fail_fraction one", mk(("health", "fail_fraction") => 1.0), ["fail_fraction"]),
        ("fail_fraction negative", mk(("health", "fail_fraction") => -0.01), ["fail_fraction"]),
        ("stop_failed_fraction one", mk(("health", "stop_failed_fraction") => 1.0), ["stop_failed_fraction"]),
        ("stop_failed_fraction negative", mk(("health", "stop_failed_fraction") => -0.1), ["stop_failed_fraction"]),
        ("safety below one", mk(("budget", "Float64", "safety") => 0.9), ["safety"]),
        ("cost factor zero", mk(("health", "first_chunk_cost_factor") => 0.0), ["first_chunk_cost_factor"]),
        ("wall cap negative", mk(("health", "wall_cap_hours") => -1.0), ["wall_cap_hours"]),
        ("control_bound zero", mk(("health", "control_bound") => 0.0), ["control_bound"]),
        ("drift_bound zero", mk(("health", "drift_bound") => 0.0), ["drift_bound"]),
        ("setup_check_tol zero", mk(("health", "setup_check_tol") => 0.0), ["setup_check_tol"]),
        ("k19 threads zero", mk(("k19", "threads") => 0), ["threads"]),
        ("rerun_sample zero", mk(("verify", "rerun_sample") => 0), ["rerun_sample"]),
        ("setup_check_steps beyond the shortest horizon", mk(("health", "setup_check_steps", "T") => 10^6), ["setup_check_steps"]),
        ("setup_check_steps class missing", mk(("health", "setup_check_steps", "X") => D), ["setup_check_steps"]),
        ("duplicate h4 label", mk(("h4", 2, "label") => "ex_h4_a"), ["label"]),
        ("pilot horizon above the shortest arm horizon", mk(("pilot", "horizon") => 1000), ["pilot", "horizon"]),
        ("pilot orbits_per_k zero", mk(("pilot", "orbits_per_k") => 0), ["orbits_per_k"]),
        ("accepted hash not hexadecimal", mk(("panel_prereg_sha256_accepted",) => ["xyz"]), ["panel_prereg_sha256_accepted"]),
        ("accepted list empty", mk(("panel_prereg_sha256_accepted",) => String[]), ["panel_prereg_sha256_accepted"]),
        ("h1_check not a boolean", mk(("arm", 1, "h1_check") => "yes"), ["h1_check"]),
        ("eps not positive", mk(("arm", 1, "eps") => [0.0]), ["eps"]),
        ("tau not positive", mk(("arm", 2, "tau") => [0]), ["tau"]),
        ("budget missing for a used precision", mk(("budget", "Double64") => D), ["budget", "Double64"]),
    ]
    @test length(cases) >= 30
    for (name, d, keys) in cases
        @testset "$name" begin
            @test ex_rejects(() -> EX.run_block(d), keys)
        end
    end
end

@testset "D53a refusals before any computation, nothing written" begin
    fx = ex_fix(:small)
    none(store, man) = !ispath(store) && !ispath(man)
    run_it(f; store, man, kw...) = EX.run_campaign(f.prereg, f.panel; store = store, manifest = man, threads = 4, code_version = EX_CODE,
                                                     allow_dirty = true, kw...)
    @testset "preregistration hash differs from its header" begin
        bad = ex_path("a-prereg-bad.md")
        write(bad, replace(read(fx.prereg, String), "End of file." => "End of file!"))
        store, man = ex_path("a-store-1"), ex_path("a-man-1.json")
        @test ex_rejects(() -> run_it((; prereg = bad, panel = fx.panel); store, man), ["SHA-256", "header"]) && none(store, man)
        @test ex_rejects(() -> EX.plan_run(bad, fx.panel; code_version = EX_CODE), ["SHA-256", "header"])
    end
    @testset "panel hash outside the accepted list" begin
        other = ex_fix_with(fx.run, "a-badpanel"; accepted = "d"^64)
        store, man = ex_path("a-store-2"), ex_path("a-man-2.json")
        @test ex_rejects(() -> run_it(other; store, man), ["panel_prereg_sha256_accepted"]) && none(store, man)
    end
    @testset "thread count differs from [run.k19].threads without override" begin
        store, man = ex_path("a-store-3"), ex_path("a-man-3.json")
        @test ex_rejects(() -> EX.run_campaign(fx.prereg, fx.panel; store = store, manifest = man, threads = 3, allow_dirty = true,
                                               code_version = EX_CODE), ["threads", "k19"]) && none(store, man)
    end
    @testset "malformed block in a file" begin
        bad = ex_fix_with(ex_rb_with(fx.run, ("grid", "first_steps") => 0), "a-badblock")
        store, man = ex_path("a-store-4"), ex_path("a-man-4.json")
        @test ex_rejects(() -> run_it(bad; store, man), ["first_steps"]) && none(store, man)
        @test ex_rejects(() -> EX.read_run_block(bad.prereg), ["first_steps"])
    end
    @testset "a plan with zero items is an error" begin
        orbs = ex_panel_orbits()
        foreach(o -> o["in_SL"] = false, orbs)
        pp = ex_write_panel(ex_path("a-nosl-panel.toml"); orbits = orbs)
        run = deepcopy(fx.run); run["arm"] = [ex_arm("only", "S_L", ["L1"], ["Float64"], [100])]
        z = ex_fix_with(run, "a-zero")
        @test ex_rejects(() -> EX.plan_run(z.prereg, pp; code_version = EX_CODE), ["items", "plan"])
    end
end

@testset "D53r the real A4 block is accepted; no study value in the code" begin
    rb = EX.read_run_block(EX_REAL_PREREG)                            # the frozen A4 block as written (also an input check)
    @test EX.block_hash(rb) == ex_sha256(EX.block_text(rb))
    real = ex_real_run()
    @test EX.block_hash(EX.run_block(real)) == EX.block_hash(rb)
    # the two fixture blocks differ in every value of R1 to R5
    a, b = ex_rb(:small), ex_rb(:alt)
    @test EX.block_hash(EX.run_block(a)) != EX.block_hash(EX.run_block(b))
    for k in ("tau", "eps", "direction"); @test a["primary"][k] != b["primary"][k]; end
    for k in ("points_per_decade", "first_steps", "extra_fractions"); @test a["grid"][k] != b["grid"][k]; end
    for k in ("fail_fraction", "control_bound", "control_methods", "setup_check_steps", "setup_check_tol", "drift_bound",
              "stop_failed_fraction", "first_chunk_cost_factor", "wall_cap_hours")
        @test a["health"][k] != b["health"][k]
    end
    @test a["backend"]["Float64"] != b["backend"]["Float64"] && a["k19"] != b["k19"] && a["budget"] != b["budget"]
    @test [x["methods"] for x in a["arm"]] != [x["methods"] for x in b["arm"]]
    @test [x["horizons"] for x in a["arm"]] != [x["horizons"] for x in b["arm"]]
    # source scan: the distinctive numbers and strings of the study block do not occur in the module
    nums = Set{Float64}()
    strs = Set{String}()
    function walk(x)
        x isa AbstractDict && foreach(walk, values(x))
        x isa AbstractVector && foreach(walk, x)
        x isa AbstractString && length(x) >= 12 && occursin(r"\d", x) && push!(strs, x)
        if x isa Real && !(x isa Bool)
            v = abs(Float64(x))
            (v >= 1e5 || (v > 0 && length(replace(string(v), r"[.e+-]|^0+" => "")) >= 4)) && push!(nums, v)
        end
    end
    walk(real)
    @test 33.09 in nums && 1553.8 in nums && 1e8 in nums && !isempty(strs)
    # the files of the first scope (panel builder, D49) are covered by the base spec; every other file of the module is scanned
    old = ("parallel", "geometry", "inputs", "prereg", "periodic", "newton", "orbits", "search", "screening", "predicates", "draw",
           "subsets", "provenance", "cells", "panel_tables", "panel", "portrait")
    viol = String[]
    for f in readdir(joinpath(@__DIR__, "..", "..", "src", "Experiments"); join = true)
        endswith(f, ".jl") && !(replace(basename(f), ".jl" => "") in old) || continue
        text = read(f, String)
        for s in strs
            occursin(s, text) && push!(viol, "$(basename(f)): \"$s\"")
        end
        for m in eachmatch(r"(?<![\w.])(\d[\d_]*(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
            v = tryparse(Float64, replace(m.captures[1], "_" => ""))
            v !== nothing && v in nums && push!(viol, "$(basename(f)): $(m.captures[1])")
        end
    end
    @test isempty(viol)
end
