# Resume (R8): the store is scanned at start. A complete file whose chunk hash equals the plan's is skipped, a final failed file
# is skipped unless a retry is asked for, a file with another hash is stale; temporary files of a killed run are listed and removed.

chunk_path(store, ordinal) = joinpath(store, "chunk-$(lpad(ordinal, 5, '0')).h5")

# A file moved aside is never deleted: `dir/name`, numbered when the name is taken.
function move_aside!(store, path, dir)
    d = joinpath(store, dir)
    mkpath(d)
    dest = joinpath(d, basename(path))
    k = 0
    while ispath(dest)
        k += 1
        dest = joinpath(d, basename(path) * ".$k")
    end
    mv(path, dest)
    return dest
end

_chunk_files(store) = sort!(filter(f -> occursin(r"^chunk-\d{5}\.h5$", basename(f)), readdir(store; join = true)))

# (kind, ordinal) of each file of the store: :keep (complete, hash equal), :stale.
function classify_store(plan::RunPlan, store)
    keep, stale = Dict{Int, String}(), String[]
    for f in _chunk_files(store)
        ord = parse(Int, string(match(r"chunk-(\d{5})\.h5", basename(f)).captures[1]))
        a = try
            read_chunk_attrs(f)
        catch
            Dict{String, Any}()
        end
        ok = 1 <= ord <= plan.counts.chunks && get(a, "complete", false) === true && get(a, "chunk_hash", "") == plan.chunks[ord].hash &&
             get(a, "pilot", false) === false
        ok ? (keep[ord] = f) : push!(stale, f)
    end
    return keep, stale
end

# Stale files stop the run (listing them) or, on request, are moved to stale/ and recomputed.
function handle_stale!(s::Session, stale::Vector{String}, mode::Symbol)
    isempty(stale) && return
    if mode !== :move
        shown = join(basename.(first(stale, 20)), ", ") * (length(stale) > 20 ? ", ... ($(length(stale)) files)" : "")
        throw(ArgumentError("stale chunk files (chunk hash differs from the plan: code version or an input changed): $shown; " *
                            "pass stale = :move to move them to stale/ and recompute"))
    end
    for f in stale
        move_aside!(s.store, f, "stale")
        push!(s.stale, basename(f))
    end
end

"List and remove the temporary files of a killed run (they hold no accepted data)."
function remove_temporaries!(s::Session)
    tmps = sort!(filter(f -> endswith(f, ".tmp"), readdir(s.store; join = true)))
    for f in tmps
        push!(s.removed, basename(f))
        rm(f; force = true)
    end
end

# Fill the manifest entries from the store: files that stay are read back, the others are pending (or paused).
function load_entries!(s::Session, keep::Dict{Int, String}, retry::Bool)
    plan = s.ctx.plan
    s.entries = Dict{String, Any}[]
    for c in plan.chunks
        if haskey(keep, c.ordinal)
            ch = read_chunk(keep[c.ordinal])
            e = chunk_entry(c, ch.attrs, ch.orbits)
            failed = e["status"] == "failed"
            push!(s.entries, e)
            (failed && retry) || push!(s.skipped, c.ordinal)
            failed && retry && (s.entries[end] = pending_entry(c, "pending"))
        else
            push!(s.entries, pending_entry(c, "pending"))
        end
    end
end
