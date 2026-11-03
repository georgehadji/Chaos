# Amendment D96 (indicator specs held only by closure arms), stage 02. Spec: experiments-spec.md, "Amendment D96": R28, I38, I39, P116.
# D96a: a campaign whose indicator specs are all held by closure arms runs before closure and lists the specs at closure.
# D96b: golden indicator tables (I38), taken from the code before D96 (commit 8c3bac0); these pass now and guard the unchanged cases.
# D96c is the unchanged package suite; D96d (O27) is stage 04. Reuses the helpers of D86 and D94; the only helper here is d96_canon.

# canonical text of a nested table: keys sorted, one `path = repr(value)` line per leaf
function d96_canon(x, path = "")
    io = IOBuffer()
    if x isa AbstractDict
        for k in sort!(collect(string.(keys(x))))
            write(io, d96_canon(x[k], path * "/" * k))
        end
    elseif x isa AbstractVector
        for (i, v) in enumerate(x)
            write(io, d96_canon(v, path * "[$i]"))
        end
    else
        println(io, path, " = ", repr(x))
    end
    return String(take!(io))
end

# p3-shaped block: ref (each_round, L1 and L5 on P); closure arms hold the indicator specs only (tg on D, fq on D_B)
function d96_block_closure_only()
    fq = d86_spec("frequency")                                                      # windows [16, 32], run length 64
    arms = Any[ex_arm("ref", "P", ["L1", "L5"], ["Float64"], [20]; release = "each_round"),
               d86_arm("tc", "D", ["tg"], ["Float64"], [20]; release = "closure"),
               d86_arm("wc", "D_B", ["fq"], ["Float64"], [64]; release = "closure", coordinate_mode = ["wrapped"])]
    d = d94_block(arms)
    d["indicator"] = Any[d86_spec("tangent"), fq]
    return d
end

# tg in an each_round arm and in a closure arm, fq in an each_round arm: every spec has an item before closure
function d96_block_both()
    fq = d86_spec("frequency")
    arms = Any[d86_arm("tr", "P", ["tg", "fq"], ["Float64"], [64]; release = "each_round", coordinate_mode = ["wrapped"]),
               d86_arm("tc", "D", ["tg"], ["Float64"], [64]; release = "closure", coordinate_mode = ["wrapped"])]
    d = d94_block(arms)
    d["indicator"] = Any[d86_spec("tangent"), fq]
    return d
end

const D96_SPEC_IDS = ["tg", "fq"]
const D96_R19 = ["kind", "spec", "spec_sha256", "chain", "budget", "bytes_per_orbit"]

@testset "D96a closure-only indicator specs: the first run completes, the closure run lists the specs (R28, P116)" begin
    fx = d94_fix("d96a", d96_block_closure_only())
    specs = Dict(s["id"] => s for s in fx.run["indicator"])
    store = ex_store("d96a-camp")
    logp = ex_path("d96a-camp.log")
    write(logp, "")
    # before closure: the first run of the campaign completes with no table entry for the specs
    res = ex_campaign(fx; store = store, round_log = logp)                                   # RED: errors in the table (no item of "tg")
    @test res.status == "complete"
    m = EX.read_json(store * "-manifest.json")
    @test isempty(intersect(collect(keys(get(m, "indicators", Dict{String, Any}()))), D96_SPEC_IDS))
    @test !any(f -> "tc" in d94_attr(f, "arms") || "wc" in d94_attr(f, "arms"), ex_chunk_files(store))       # no closure arm chunk before closure
    # rounds 2 and 3, then the closure run
    digest = r -> d94_digest(d94_ref_hashes(store, r))
    lspecs = NamedTuple[(round = 1, decision = "continue")]
    d94_write_log(logp, fx, lspecs, digest)
    for (r, dec) in ((2, "continue"), (3, "close"))
        res = ex_campaign(fx; store = store, round_log = logp)
        @test res.status == "complete"
        push!(lspecs, (round = r, decision = dec))
        d94_write_log(logp, fx, lspecs, digest)
    end
    res = ex_campaign(fx; store = store, round_log = logp)
    @test res.status == "complete"
    m = EX.read_json(store * "-manifest.json")
    ind = m["indicators"]
    @test sort(collect(keys(ind))) == sort(D96_SPEC_IDS)
    for id in D96_SPEC_IDS
        t = ind[id]
        @test all(k -> haskey(t, k), D96_R19)
        @test t["kind"] == specs[id]["kind"] && t["spec"] == d86_spec_text(specs[id]) && t["spec_sha256"] == d86_spec_hash(specs[id])
        @test t["budget"] == Dict("Float64" => specs[id]["budget"]["Float64"])   # the block keeps only the precisions of the arms (D86 item 1)
        @test Set(keys(t["bytes_per_orbit"])) == Set(["M", "R", "S", "B", "total"])
    end
    @test haskey(ind["fq"], "windows") && ind["fq"]["windows"] == [16, 32] && ind["fq"]["run_length"] == 64
    @test !haskey(ind["tg"], "windows") && !haskey(ind["tg"], "run_length")
end

# Golden SHA-256 of d96_canon of the indicator table (EX.indicator_manifest), taken BEFORE D96 at commit 8c3bac0.
const D96_GOLDEN_BOTH = "2d0bb10a68d859cdaa97c026f4d82ea7b20e0a4c1c7887a74df35e0e2e649d9a"
const D96_GOLDEN_NOROUNDS = "331567fafe4b4d3dfbdc92ce6ce10f203daab261d73bbfffe5c4f5a212834f1e"

function d96_tables()
    fx = d94_fix("d96b", d96_block_both())
    logp = ex_path("d96b-open.log")
    write(logp, "")
    pb = d94_plan(fx; log = logp)
    pn = d86_plan("d96-norounds", d86_all_block())
    return (both = d96_canon(EX.indicator_manifest(pb)), norounds = d96_canon(EX.indicator_manifest(pn)))
end

@testset "D96b golden indicator tables are unchanged (I38)" begin
    t = d96_tables()
    @test all(id -> occursin("/" * id * "/", t.both), D96_SPEC_IDS)
    @test ex_sha256(t.both) == D96_GOLDEN_BOTH                       # spec in an each_round and a closure arm
    @test ex_sha256(t.norounds) == D96_GOLDEN_NOROUNDS               # block without rounds, six kinds
end
