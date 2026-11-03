# Amendment D97, stage 02, criterion D97i (E82, O28): the shape guard of the toy campaign. The toy store of an_p3t_helpers.jl is written by the real chunk writer, but its SHAPE
# (cells, arm names, keys of the chunks and of the manifest entries, attribute names, item labels) was chosen by the test side; the D95 toy had the shape of the module and not that of the
# runner, which hid eight differences (M1 to M8 of the discovery run of D96) from every D95 test. Here a tiny campaign is run by the real runner (`Chaoslab.Experiments.run_campaign`) with
# the six arms of the test run, round by round through a hand log, then to closure, with one chunk of the ref arm failed twice; the toy store and the real store must have the same
# structure. Nothing here calls Chaoslab.Analysis, so these tests pass before the module is corrected (they are not RED).
# The helpers of the Experiments suite (ex_*, d86_*, d94_*) are loaded into a private module, so that their constants do not meet those of the Experiments suite in `Main`.

module P3TShape
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, TOML
const XT = joinpath(@__DIR__, "..", "experiments")
include(joinpath(XT, "ex_helpers.jl"))
include(joinpath(XT, "ex_run_helpers.jl"))
include(joinpath(XT, "ex_d86_helpers.jl"))
include(joinpath(XT, "ex_d94_helpers.jl"))

# the block of the test run at the scale of a unit test: the six arms of E81 (ref released each round; the closure arms: ind and d64 hold tg, fq1, bk1, the pa1 twins hold pa1 with the
# run length 2 W + 1, chk holds L5 in Double64 and BigFloat(256) on D_B). The horizons are the run lengths of the window specs (64, 65); the ref arm runs 20 steps.
function block()
    tg = d86_spec("tangent"); fq = d86_spec("frequency"; id = "fq1"); bk = d86_spec("birkhoff"; id = "bk1"); pa = d86_spec("phase_advance"; id = "pa1")
    arms = Any[ex_arm("ref", "P", ["L1", "L5"], ["Float64"], [20]; release = "each_round", coordinate_mode = ["wrapped"]),
               d86_arm("ind", "P", ["tg", "fq1", "bk1"], ["Float64"], [64]; release = "closure", coordinate_mode = ["wrapped"]),
               d86_arm("ind_pa1", "P", ["pa1"], ["Float64"], [65]; release = "closure", coordinate_mode = ["wrapped"]),
               d86_arm("d64", "D", ["tg", "fq1", "bk1"], ["Double64"], [64]; release = "closure", coordinate_mode = ["wrapped"]),
               d86_arm("d64_pa1", "D", ["pa1"], ["Double64"], [65]; release = "closure", coordinate_mode = ["wrapped"]),
               ex_arm("chk", "D_B", ["L5"], ["Double64", "BigFloat"], [20]; release = "closure", coordinate_mode = ["wrapped"])]
    d = d94_block(arms)
    d["indicator"] = Any[tg, fq, pa, bk]
    d["health"]["stop_failed_fraction"] = 0.9                      # one failed chunk must not stop the campaign
    return d
end

# the campaign: round 1, rounds 2 and 3 through the hand log, closure; the first ref chunk of round 1 fails twice (NaN in its first item, both attempts). Returns the paths and the plan.
function campaign()
    fx = d94_fix("d97shape", block())
    store = ex_store("d97shape-camp"); logp = ex_path("d97shape.log"); write(logp, "")
    bad(orbits, info) = (first(values(orbits[1].lyap))[1, 1] = NaN)
    r1 = ex_campaign(fx; store = store, round_log = logp, hooks = EX.RunHooks(mutate = ex_fault(bad; ord = 3, attempts = (1, 2))))
    digest = r -> d94_digest(d94_ref_hashes(store, r))
    lspecs = NamedTuple[(round = 1, decision = "continue")]
    d94_write_log(logp, fx, lspecs, digest)
    res = r1
    for (r, dec) in ((2, "continue"), (3, "close"))
        res = ex_campaign(fx; store = store, round_log = logp)
        push!(lspecs, (round = r, decision = dec))
        d94_write_log(logp, fx, lspecs, digest)
    end
    res = ex_campaign(fx; store = store, round_log = logp)
    return (; store, manifest = store * "-manifest.json", status = res.status, fx)
end
end # module

const P3T_REAL = Ref{Any}(nothing)
p3t_real() = (P3T_REAL[] === nothing && (P3T_REAL[] = P3TShape.campaign()); P3T_REAL[])

# the structure of a store: the chunk files with their root attributes (names, value types) and the names of the datasets of one orbit group, by (arm, item)
function p3t_structure(store)
    out = Dict{Tuple{String, String, String}, NamedTuple}()
    for f in p3t_files(store)
        a = AX.read_chunk_attrs(f)
        k = (only(String.(a["arms"])), String(a["method"]), String(a["precision"]))
        names = Set(keys(a))
        sets = String[]
        AX.HDF5.h5open(f, "r") do h
            if haskey(h, "orbits") && !isempty(keys(h["orbits"]))
                g = h["orbits"][first(sort(collect(keys(h["orbits"]))))]
                visit(gr, pre) = for n in keys(gr)
                    o = gr[n]
                    o isa AX.HDF5.Group ? visit(o, pre * n * "/") : push!(sets, pre * n)
                end
                visit(g, "")
            end
        end
        old = get(out, k, nothing)
        out[k] = (attrs = old === nothing ? names : union(old.attrs, names), datasets = old === nothing || isempty(old.datasets) ? Set(sets) : old.datasets, key = String(a["chunk_key"]), a = a)
    end
    return out
end

@testset "D97i the toy campaign has the structure of a campaign of the runner (E82, O28 guard)" begin
    R = p3t_real()
    @test R.status == "complete with failed chunks"                                   # the failed ref chunk of the real campaign; the toy has failed chunks as well
    S = p3t_cached(); I = p3t_cached_inputs()
    toy_man = AX.read_json(S.manifest); real_man = AX.read_json(R.manifest)
    # the manifest: the keys that the module reads are keys of the runner, the arm names are the six of E81, the chunk entries have the same keys
    read_keys = ["study", "status", "arms", "chunks", "prereg", "panel", "plan_sha256", "schema_version"]
    @test all(k -> haskey(real_man, k) && haskey(toy_man, k), read_keys) && issubset(Set(keys(toy_man)), Set(keys(real_man)))
    @test Set(keys(real_man["arms"])) == Set(keys(toy_man["arms"])) == Set(P3T_ARMS)
    arm_keys(m, a) = Set(keys(m["arms"][a]))
    @test all(a -> issubset(arm_keys(toy_man, a), arm_keys(real_man, a)) && ("specs" in arm_keys(real_man, a)) == ("specs" in arm_keys(toy_man, a)), P3T_ARMS)
    @test sort([a for a in P3T_ARMS if haskey(toy_man["arms"][a], "specs")]) == ["d64", "d64_pa1", "ind", "ind_pa1"] == sort([a for a in P3T_ARMS if haskey(real_man["arms"][a], "specs")])
    @test Set(keys(real_man["arms"]["ind"]["specs"])) == Set(["tg", "fq1", "bk1"]) == Set(keys(toy_man["arms"]["ind"]["specs"])) && Set(keys(real_man["arms"]["ind_pa1"]["specs"])) == Set(["pa1"]) == Set(keys(toy_man["arms"]["ind_pa1"]["specs"]))
    @test all(s -> issubset(Set(keys(s)), Set(["budget_line_ns", "measured_ns_per_orbit_step"])) && haskey(s, "measured_ns_per_orbit_step"), (v for a in ("ind", "ind_pa1", "d64", "d64_pa1") for v in values(toy_man["arms"][a]["specs"])))
    ent_keys = Set(["ordinal", "key", "hash", "arms", "status", "steps", "thread_seconds", "orbits", "attempt", "flags"])
    @test all(e -> issubset(ent_keys, Set(keys(e))) && issubset(Set(keys(e)), Set(keys(real_man["chunks"][1])) ∪ Set(["indicator_flags"])), toy_man["chunks"]) && all(e -> issubset(ent_keys, Set(keys(e))), real_man["chunks"])
    @test Set(e["status"] for e in real_man["chunks"]) == Set(["ok", "failed"]) && Set(e["status"] for e in toy_man["chunks"]) == Set(["ok", "failed"])
    @test all(e -> e["arms"] isa AbstractVector && length(e["arms"]) == 1 && e["steps"] isa Integer && e["thread_seconds"] isa Real, vcat(toy_man["chunks"], real_man["chunks"]))
    # the key of a chunk entry: the same fields in the same order, the cell of the form of E74
    keyfields(k) = [first(split(w, "=")) for w in split(k, " ")]
    cellform = r"cell=k\d+:[A-Z]:r\d+$"
    @test all(e -> occursin(cellform, e["key"]), vcat(toy_man["chunks"], real_man["chunks"]))
    real_forms = Set(keyfields(e["key"]) for e in real_man["chunks"]); toy_forms = Set(keyfields(e["key"]) for e in toy_man["chunks"])
    @test toy_forms == real_forms
    # the item with the L1 method (cost of E78): the key of the entries of the ref arm names the method as `method=L1 `
    @test any(e -> startswith(e["key"], "method=L1 ") && e["arms"] == ["ref"], toy_man["chunks"]) && any(e -> startswith(e["key"], "method=L1 ") && e["arms"] == ["ref"], real_man["chunks"])
    # the chunk files: the attribute names of the toy are attribute names of the runner for the same arm and item, the values of the keys that the module reads have the same forms
    Ts, Rs = p3t_structure(S.store), p3t_structure(R.store)
    @test Set(k[1] for k in keys(Rs)) == Set(P3T_ARMS) == Set(k[1] for k in keys(Ts))
    @test Set((k[1], k[3]) for k in keys(Ts)) == Set((k[1], k[3]) for k in keys(Rs))             # arm and precision labels: Float64, Double64, BigFloat(256)
    @test sort(unique(k[3] for k in keys(Ts))) == ["BigFloat(256)", "Double64", "Float64"]
    read_attrs = ["study", "panel_sha256", "prereg_sha256", "chunk_key", "chunk_hash", "arms", "method", "precision", "horizon", "status", "complete", "tau", "direction", "attempt"]
    for k in keys(Ts)
        @test haskey(Rs, k)
        haskey(Rs, k) || continue
        real, toy = Rs[k], Ts[k]
        @test issubset(toy.attrs, real.attrs) && issubset(Set(a for a in read_attrs if a != "direction" || k[2] in ("L1", "L5", "tg")), toy.attrs ∩ real.attrs)
        @test issubset(toy.datasets, real.datasets)
    end
end
