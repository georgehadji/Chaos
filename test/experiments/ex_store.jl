# D53m store layout, D53n atomic write. Spec: Amendment D53, R6, R10; schema/run-store.md (version 2).
# The files are read back through the module (read_chunk, chunk_layout); the independent reader (O13) is stage 04.

const EX_SMALL = Dict{String, Any}()
# uninterrupted small-fixture campaign, cached by tag
function ex_small_run(tag::AbstractString; kind::Symbol = :small, threads::Integer = 4, hooks = nothing, kw...)
    return get!(EX_SMALL, tag) do
        fx = ex_fix(kind)
        store = ex_store("small-$tag")
        res = hooks === nothing ? ex_campaign(fx; store = store, threads = threads, kw...) :
                                  ex_campaign(fx; store = store, threads = threads, hooks = hooks, kw...)
        (; fx, store, res, plan = ex_plan(fx))
    end
end

# root attribute names of schema/run-store.md, from the table
function ex_schema_attrs()
    text = replace(read(joinpath(@__DIR__, "..", "..", "schema", "run-store.md"), String), "\r\n" => "\n")
    sec = match(r"## Root attributes \(provenance\)\n(.*?)\n## Groups"s, text).captures[1]
    names = String[]
    for line in split(sec, '\n')
        startswith(line, "| `") || continue
        append!(names, [String(m.captures[1]) for m in eachmatch(r"`([a-z_0-9]+)`", split(line, '|')[2])])
    end
    return names
end

function ex_design_class(orbits, id, run)
    nP = maximum(o["id"] for o in orbits)
    id > nP && return "h4"
    o = orbits[id]
    return o["control"] ? "control" : o["class"]
end

@testset "D53m store layout" begin
    run = ex_small_run("ref")
    alt = ex_small_run("alt"; kind = :alt, threads = 2)
    @test run.res.status == "complete" && alt.res.status == "complete"
    names = ex_schema_attrs()
    @test length(names) >= 40 && "chunk_hash" in names && "complete" in names
    bad = String[]
    note!(ok, msg) = ok || push!(bad, msg)
    for r in (run, alt)
        files = ex_chunk_files(r.store)
        note!([ex_ordinal(f) for f in files] == collect(1:r.plan.counts.chunks), "file names are the chunk ordinals")
        prereg_header = EX.check_prereg(r.fx.prereg).header_sha256
        block_hash = EX.block_hash(EX.read_run_block(r.fx.prereg))
        for f in files
            ch = EX.read_chunk(f)
            c = r.plan.chunks[ex_ordinal(f)]
            a = ch.attrs
            it0 = r.plan.items[first(c.items)]
            for nm in names
                nm in ("python_version", "perturbation") && continue
                haskey(a, nm) || push!(bad, "$(basename(f)): attribute $nm missing")
            end
            note!(haskey(a, "perturbation") == (it0.method == "L2"), "perturbation exactly for L2")
            it0.method == "L2" && note!(get(a, "perturbation", nothing) == it0.eps, "perturbation text")
            note!(a["schema_version"] == 2 && a["complete"] === true && a["pilot"] === false, "version, complete, pilot")
            note!(a["status"] == "ok" && a["attempt"] == 1 && a["chunk"] == c.ordinal, "status, attempt, chunk")
            note!(a["chunk_key"] == c.key && a["chunk_hash"] == c.hash && a["plan_sha256"] == r.plan.hash, "key, hashes")
            note!(a["panel_sha256"] == ex_sha256(read(r.fx.panel)) && a["prereg_sha256"] == prereg_header && a["config_sha256"] == block_hash,
                  "input hashes")
            note!(a["code_version"] == EX_CODE && a["study"] == "ex-fixture", "code version, study")
            note!(a["method"] == it0.method && a["precision"] == it0.precision && a["horizon"] == it0.horizon && a["tau"] == it0.tau &&
                  a["vectors"] == it0.vectors && a["direction"] == it0.direction && a["coordinate_mode"] == it0.coordinate_mode, "run type")
            note!(a["grid_sha256"] == it0.grid_sha256, "grid hash")
            note!(a["backend"] == "cpu" && a["backend_request"] in ("auto", "cpu") && a["device"] == "cpu", "backend")
            note!(a["first_attempt_reasons"] == String[] && a["failure_reasons"] == String[] && a["arms"] == c.arms, "reasons, arms")
            note!(occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$", a["created_utc"]), "created_utc format")
            note!(occursin(r"^[0-9a-f]{40}$", a["git_sha"]), "git sha")
            note!(all(k -> a[k] isa Integer, ("renorm_interval_tau", "deviation_vectors_k", "seed", "chunk_size", "gpu_chunks",
                                                "steps_per_launch", "launches_per_chunk")), "tangent-core provenance integers")
            note!(a["orthonormalizer"] isa String && a["implementation"] isa String, "gram-schmidt provenance")
            lay = EX.chunk_layout(f)
            note!(!any(k -> startswith(k, "/grid"), keys(lay)), "/grid/ unused")
            for o in ch.orbits
                it = ex_item_of(r.plan, c, o)
                n = length(o.t)
                need = ("system", "params", "ic", "class_label", "seed", "item_hash", "design_class", "flags", "steps_recorded")
                note!(all(k -> haskey(o.attrs, k), need), "orbit $(o.id) attributes")
                note!(o.attrs["item_hash"] == it.hash && o.attrs["steps_recorded"] == n && string(o.attrs["seed"]) == it.orbit_seed, "orbit values")
                note!(o.attrs["design_class"] == ex_design_class(r.fx.orbits, o.id, r.fx.run) && o.attrs["flags"] isa Vector{String}, "design class, flags")
                b = "/orbits/$(o.id)"
                note!(lay["$b/t"] == (eltype = "Int64", size = (n,)), "t dtype")
                note!(lay["$b/state"] == (eltype = "Float64", size = (n, 2)), "state dtype")
                if it.coordinate_mode == "wrapped"
                    note!(get(lay, "$b/winding", nothing) == (eltype = "Int64", size = (n, 2)), "winding dtype")
                else
                    note!(!haskey(lay, "$b/winding"), "no winding in lifted mode")
                end
                note!(lay["$b/lyap/$(it.method)"] == (eltype = "Float64", size = (n, ex_vectors(it.method))), "lyap dtype")
                for d in ex_diag_names(it.method)
                    note!(lay["$b/diagnostics/$d"] == (eltype = "Float64", size = (n,)), "diagnostics dtype")
                end
                note!(lay["$b/cost/$(it.method)"] == (eltype = "Float64", size = (5,)), "cost dtype")
            end
        end
    end
    @test isempty(bad)
    # schema/run-store.md states version 2 with the changes of R10 and keeps the version 1 parts
    doc = replace(read(joinpath(@__DIR__, "..", "..", "schema", "run-store.md"), String), "\r\n" => "\n")
    @test occursin("Run Store Schema (HDF5), version 2", doc) && occursin("## Changes from v1", doc)
    for nm in ("schema_version", "study", "chunk", "git_sha", "prereg_sha256", "julia_version", "cpu", "threads", "backend", "device",
               "backend_fallback_reason", "created_utc", "complete")
        @test nm in names
    end
    @test occursin("/grid/", doc) && occursin("N-body files use", doc)
end

@testset "D53m write and read back bit for bit" begin
    run = ex_small_run("ref")
    f = ex_chunk_files(run.store)[1]
    ch = EX.read_chunk(f)
    # the hook sees the final buffers before the write; the file equals them
    seen = Dict{Int, Any}()
    obs = (info, attrs, orbits) -> (seen[info.ordinal] = (; attrs = deepcopy(attrs), orbits = deepcopy(orbits)); nothing)
    r = ex_campaign(run.fx; store = ex_store("m-observe"), hooks = EX.RunHooks(observe = obs))
    @test sort(collect(keys(seen))) == collect(1:run.plan.counts.chunks)
    for file in ex_chunk_files(ex_path("m-observe"))
        s, back = seen[ex_ordinal(file)], EX.read_chunk(file)
        @test ex_canon(back) == ex_canon(s)
        @test all(i -> back.orbits[i].attrs == s.orbits[i].attrs, eachindex(s.orbits))
    end
    # special values survive: NaN, Inf, -0.0 and a subnormal are stored as written (a failed chunk holds the series that failed)
    orbs = deepcopy(ch.orbits)
    m = first(keys(orbs[1].lyap))
    orbs[1].lyap[m][1, 1] = NaN; orbs[1].lyap[m][2, 1] = Inf; orbs[1].lyap[m][3, 1] = -0.0; orbs[1].lyap[m][4, 1] = 5.0e-324
    p = ex_path("m-special.h5")
    EX.write_chunk(p, ch.attrs, orbs)
    back = EX.read_chunk(p)
    @test ex_canon(back; attrs = true) == ex_canon((; attrs = ch.attrs, orbits = orbs); attrs = true)
    @test back.attrs["complete"] === true && isempty(ex_tmp_files(dirname(p)))
    @test_throws Exception EX.read_chunk(ex_path("m-missing.h5"))
end

@testset "D53n atomic write killed at 20 seeded points" begin
    ref = ex_small_run("ref")
    refcanon = ex_store_canon(ref.store; attrs = true)
    M = ref.plan.counts.chunks
    points = (:before_data, :during_data, :after_data, :before_rename, :after_rename)
    rng = StableRNG(ex_s53('n', 1))
    final(store, o) = joinpath(store, "chunk-$(lpad(o, 5, '0')).h5")
    for t in 1:20
        ord, pt = rand(rng, 1:M), points[1 + (t - 1) % 5]
        store = ex_store("n-$t")
        @test ex_killed(() -> ex_campaign(ref.fx; store = store, hooks = EX.RunHooks(kill = ex_kill_at(pt, ord))))
        # no file with a final name is incomplete (I13)
        @test all(f -> EX.read_chunk(f).attrs["complete"] === true, ex_chunk_files(store))
        @test isfile(final(store, ord)) == (pt === :after_rename)
        tmps = basename.(ex_tmp_files(store))
        res = ex_campaign(ref.fx; store = store)
        @test sort(String.(EX.read_json(res.manifest)["removed_temporaries"])) == sort(tmps)   # listed in the manifest, then removed
        @test isempty(ex_tmp_files(store))
        @test (ord in res.ran) == (pt !== :after_rename)               # the chunk is recomputed unless its file was already final
        @test res.status == "complete"
        @test ex_store_canon(store; attrs = true) == refcanon
    end
    # a temporary file left by an earlier crash in a finished store is listed and removed, nothing is recomputed
    store = ex_store("n-leftover")
    ex_campaign(ref.fx; store = store)
    write(joinpath(store, "chunk-00003.h5.tmp"), UInt8[1, 2, 3])
    res = ex_campaign(ref.fx; store = store)
    @test isempty(res.ran) && isempty(ex_tmp_files(store))
    @test String.(EX.read_json(res.manifest)["removed_temporaries"]) == ["chunk-00003.h5.tmp"]
end
