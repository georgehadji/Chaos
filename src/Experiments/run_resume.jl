# Resume (R8): the store is scanned at start. A complete file whose chunk hash equals the plan's is skipped, a final failed file
# is skipped unless a retry is asked for, a file with another hash (or panel file) is stale, a file written by a control-bound stop
# stops the resume; temporary files of a killed run are listed and removed.

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

# Files of the store by kind: `keep` (complete, chunk hash and panel file equal to the plan's), `halted` (kept but written by a
# control-bound stop, H1) and `stale` (everything else). A different panel file is stale even when the chunk hash is equal (M6).
function classify_store(plan::RunPlan, store)
    keep, halted, stale = Dict{Int, String}(), Dict{Int, String}(), String[]
    for f in _chunk_files(store)
        ord = parse(Int, string(match(r"chunk-(\d{5})\.h5", basename(f)).captures[1]))
        a = try
            read_chunk_attrs(f)
        catch
            Dict{String, Any}()
        end
        ok = 1 <= ord <= plan.counts.chunks && get(a, "complete", false) === true && get(a, "chunk_hash", "") == plan.chunks[ord].hash &&
             get(a, "pilot", false) === false && get(a, "panel_sha256", "") == plan.panel.sha256
        ok ? ((get(a, "control_stop", false) === true ? halted : keep)[ord] = f) : push!(stale, f)
    end
    return keep, stale, halted
end

# A control-bound stop is not a result: a plain resume refuses (H1).
function refuse_halted(halted::Dict{Int, String})
    isempty(halted) && return
    ords = sort!(collect(keys(halted)))
    throw(ArgumentError("chunk files $(join((basename(halted[o]) for o in ords), ", ")) were written by a campaign stopped at a control bound " *
                        "(see control_stops in the manifest); resuming repeats the same code: change the code version, or pass " *
                        "retry_control_stop = true to recompute them (recorded)"))
end

# With the retry option the files are moved to aside/ and recomputed; the retry is recorded in the overrides and on the stop records.
function handle_halted!(s::Session, halted::Dict{Int, String})
    isempty(halted) && return
    ords = sort!(collect(keys(halted)))
    foreach(o -> move_aside!(s.store, halted[o], "aside"), ords)
    s.overrides["retry_control_stop"] = ords
    for r in s.control_stops
        r isa AbstractDict && get(r, "ordinal", 0) in ords && (r["retried"] = true)
    end
end

# Backend table of the previous session (precision, method => (request, backend, reason)); empty when the manifest has none.
function _prev_backends(prev)
    prev isa AbstractDict || return Dict{Tuple{String, String}, Tuple{String, String, String}}()
    return Dict{Tuple{String, String}, Tuple{String, String, String}}((String(r["precision"]), String(r["method"])) =>
        (String(r["request"]), String(r["backend"]), String(r["reason"])) for r in get(prev, "backend_table", Any[]))
end

# H2 (spec R4 item 3): a resumed run whose resolved backend table differs from the previous session's, or kept chunks that disagree with
# the table, is refused unless a backend override is given (then recorded).
function check_backends!(s::Session, keep::Dict{Int, String}, override::Bool)
    isempty(keep) && return
    old = _prev_backends(s.prev)
    now = Dict((p, m) => (String(r.request), String(r.backend), r.reason) for ((p, m), r) in s.table)
    problems = String[]
    if !isempty(old)
        for k in sort!(collect(union(keys(old), keys(now))))
            get(old, k, nothing) == get(now, k, nothing) || push!(problems, "$(k[1]) $(k[2]): was $(get(old, k, "absent")), now $(get(now, k, "absent"))")
        end
    end
    for o in sort!(collect(keys(keep)))
        a = read_chunk_attrs(keep[o])
        row = _row(s, s.ctx.plan.items[first(s.ctx.plan.chunks[o].items)])
        b = get(a, "backend", "")
        off = (b == "gpu" && row.backend !== :gpu) || (b == "cpu" && row.backend === :gpu && isempty(get(a, "backend_fallback_reason", "")))
        off && push!(problems, "chunk $o was written on $b, the table now says $(row.backend)")
    end
    isempty(problems) && return
    override || throw(ArgumentError("the resolved backend table differs from the previous session's (device, record or thread count changed): " *
                                    join(first(problems, 10), "; ") * "; pass backend_override = true to resume anyway (recorded)"))
    s.overrides["backend_table"] = problems
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
            e = chunk_entry(c, ch.attrs, ch.orbits; plan = plan)
            failed = e["status"] == "failed"
            push!(s.entries, e)
            (failed && retry) || push!(s.skipped, c.ordinal)
            failed && retry && (s.entries[end] = pending_entry(c, "pending"))
        else
            push!(s.entries, pending_entry(c, "pending"))
        end
    end
end
