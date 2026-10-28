# Included by verify/experiments-indicators.jl (section O19): item text and hash, spec text and hash, chunk key and chunk hash of every work item and chunk of the reduced plan,
# recomputed from the fields of the item and the block with the platform's command-line SHA-256 (sha256sum) and compared with the plan and with the files of campaign B (extends O8).
# The text of an item is written here from E27 (the D53 text of R1 item 3 without `eps`, `direction` only for the items that have one, then `kind=` and `spec=`).
const SHA_DIR = joinpath(SCR, "sha")
function cli_sha256(texts::Vector{String})
    rm(SHA_DIR; recursive = true, force = true); mkpath(SHA_DIR)
    out = String[]
    for chunk in Iterators.partition(eachindex(texts), 800)
        names = String[]
        for i in chunk
            nm = "t$i.txt"
            write(joinpath(SHA_DIR, nm), texts[i])
            push!(names, nm)
        end
        txt = read(Cmd(`sha256sum $names`; dir = SHA_DIR), String)
        append!(out, [String(split(l)[1]) for l in split(strip(txt), '\n')])
    end
    rm(SHA_DIR; recursive = true, force = true)
    return out
end
Sys.which("sha256sum") === nothing && error("no sha256sum on the path")

# canonical text of a spec without its budget, from the TOML values (keys sorted, lists in order, one `path = value` line each)
function canon_spec(spec)
    io = IOBuffer()
    cv(x::AbstractString) = repr(String(x))
    cv(x::AbstractVector) = "[" * join((cv(e) for e in x), ", ") * "]"
    cv(x) = repr(x)
    function rec(x, path)
        if x isa AbstractDict
            foreach(k -> rec(x[k], path * "." * k), sort!(collect(keys(x))))
        elseif x isa AbstractVector && !isempty(x) && first(x) isa AbstractDict
            foreach(((i, e),) -> rec(e, "$path[$i]"), enumerate(x))
        else
            println(io, path, " = ", cv(x))
        end
    end
    rec(Dict{String, Any}(k => v for (k, v) in spec if k != "budget"), "indicator")
    return String(take!(io))
end
function item_text_ind(it)
    io = IOBuffer()
    for f in (:orbit, :system, :param, :x0, :p0, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau, :eps, :direction, :coordinate_mode, :vectors, :grid_sha256, :code_version, :hash_version, :kind, :spec)
        v = hasproperty(it, f) ? getproperty(it, f) : nothing
        (f in (:eps, :direction, :kind, :spec) && v === nothing) && continue
        print(io, f, "=", v, "\n")
    end
    return String(take!(io))
end
const BLOCK_SPECS = Dict(s["id"] => s for s in block_run()["indicator"])
let
    ids = sort(collect(keys(BLOCK_SPECS)))
    shs = cli_sha256([canon_spec(BLOCK_SPECS[i]) for i in ids])
    SPECSH = Dict(ids[i] => shs[i] for i in eachindex(ids))
    ok = all(i -> shs[i] == EX.spec_hash(SPECS[ids[i]]) && canon_spec(BLOCK_SPECS[ids[i]]) == EX.spec_text(SPECS[ids[i]]), eachindex(ids))
    emit("O19/D86d", "spec text (canonical, no budget) and spec hash of the 6 specs: recomputed from the TOML of verify/experiments-indicators-block.toml with sha256sum, equal to the module's", ok, "6 of 6 equal", pf(ok))
    its = PLAN.items
    texts = [item_text_ind(it) for it in its]
    sh = cli_sha256(texts)
    nt = count(i -> texts[i] != EX.item_text(its[i]), eachindex(its))
    nh = count(i -> sh[i] != its[i].hash, eachindex(its))
    nspec = count(i -> its[i].kind !== nothing && its[i].spec != SPECSH[its[i].method], eachindex(its))
    emit("O19/D86d", "item text (written from E27) and item hash (sha256sum) of all $(length(its)) work items of the reduced plan: texts differing from the module's / hashes differing from the plan's", "$nt / $nh", "0 / 0", pf(nt == 0 && nh == 0))
    nk = count(it -> it.kind !== nothing, its)
    emit("O19/D86d", "of which items of the six new kinds / items of L1; the `spec=` line of every new item equals the independently computed spec hash", "$nk / $(length(its) - nk); mismatching spec lines $nspec", "0 mismatching", pf(nspec == 0))
    # chunk keys: rebuilt from the fields of the first item of each chunk
    nkey = 0
    for c in PLAN.chunks
        it = its[first(c.items)]
        key = it.kind === nothing ? "method=$(it.method) precision=$(it.precision) horizon=$(it.horizon) tau=$(it.tau) eps=$(it.eps === nothing ? "-" : it.eps) direction=$(it.direction === nothing ? "-" : it.direction) mode=$(it.coordinate_mode) cell=$(it.cell)" :
              "method=$(it.method) kind=$(it.kind) spec=$(it.spec[1:16]) precision=$(it.precision) horizon=$(it.horizon) tau=$(it.tau) eps=- direction=$(it.direction === nothing ? "-" : it.direction) mode=$(it.coordinate_mode) cell=$(it.cell)"
        nkey += key == c.key ? 0 : 1
    end
    emit("O19/D86d", "chunk key of every one of the $(length(PLAN.chunks)) chunks rebuilt from the item fields (second form of R12 item 3 for the new kinds, the old form for L1): differing from the plan's", nkey, "0", pf(nkey == 0))
    # chunk hash: the module's chunk text hashed by sha256sum
    ctexts = [EX.chunk_text(PLAN, c.ordinal) for c in PLAN.chunks]
    shc = cli_sha256(ctexts)
    nch = count(i -> shc[i] != PLAN.chunks[i].hash, eachindex(shc))
    emit("O19", "chunk hash of all chunks: sha256sum of the chunk text equals the plan's hash; differing", nch, "0", pf(nch == 0))
    # the plan hash
    shp = cli_sha256([PLAN.text])
    emit("O19", "plan hash: sha256sum of the plan text equals the plan's hash", shp[1] == PLAN.hash, "true", pf(shp[1] == PLAN.hash))
    # the files of campaign B: item_hash of every orbit, chunk_key, chunk_hash, spec_sha256 against the plan
    R = load_store(STORE_B)
    nbad = n = 0
    for (ord, ch) in R
        c = PLAN.chunks[ord]
        a = ch.attrs
        n += 1
        ok = a["chunk_key"] == c.key && a["chunk_hash"] == c.hash && a["plan_sha256"] == PLAN.hash
        its_c = [PLAN.items[i] for i in c.items]
        ok &= sort([o.attrs["item_hash"] for o in ch.orbits]) == sort([it.hash for it in its_c])
        haskey(a, "kind") && (ok &= a["spec_sha256"] == its_c[1].spec)
        nbad += ok ? 0 : 1
    end
    emit("O19/D86d", "files of campaign B: chunk_key, chunk_hash, plan_sha256, the item hashes of the orbits and spec_sha256 equal the independently recomputed values; files compared / differing", "$n / $nbad", "0 differ", pf(nbad == 0 && n == NCH))
end
