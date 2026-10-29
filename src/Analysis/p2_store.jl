# Study p2-indicators (Amendment D89): the reader of the run store of the indicator items (E25, schema/run-store.md "Added by Amendment D86") and the grid starts (E24).
# One item per (orbit, chunk), found by content; the Julia shapes of the schema ("Axis order") are returned.

# The tables of an item that has an indicator spec (`indicators/<id>`); an item of a catalog method has none.
function _p2_tables(o, id)
    t = get(o.indicators, id, nothing)
    t === nothing && return (series = nothing, columns = String[], final = nothing, rows = String[], fcolumns = String[], flagwords = nothing, bits = String[],
                             scalars = nothing, scolumns = String[], range = nothing)
    return (series = t.series, columns = t.columns, final = t.final, rows = t.rows, fcolumns = t.fcolumns, flagwords = t.flags, bits = t.bits,
            scalars = t.scalars, scolumns = t.scolumns, range = t.range)
end

# Items of one chunk file: the root attributes are those of the spec, the first exponent of a catalog method is `lam`.
function _p2_chunk_items(ch, ordinal::Integer)
    a = ch.attrs
    id = String(a["method"])
    return [(; orbit = o.id, id = id, kind = get(a, "kind", nothing), precision = String(a["precision"]), horizon = Int(a["horizon"]), tau = Int(a["tau"]),
             direction = haskey(a, "direction") ? String(a["direction"]) : nothing, mode = String(a["coordinate_mode"]), arms = String.(a["arms"]), chunk = Int(ordinal),
             t = o.t, state = o.state, winding = o.winding, flags = String.(o.attrs["flags"]), attrs = o.attrs,
             lam = haskey(o.lyap, id) ? o.lyap[id][:, 1] : nothing, _p2_tables(o, id)...) for o in ch.orbits]
end

"""
    p2_read_store(manifest, store) -> Vector{NamedTuple}

Every orbit record of every chunk file the manifest lists (E25): `orbit`, `id` (the spec or method id), `kind`, `precision`, `horizon`, `tau`, `direction`
(`nothing` if none), `mode`, `arms`, `chunk`, `t`, `state`, `winding`, `flags`, `attrs`, `lam` (first exponent or `nothing`), `series`, `columns`, `final`, `rows`,
`fcolumns`, `flagwords`, `bits`, `scalars`, `scolumns`, `range`.
"""
function p2_read_store(manifest::AbstractString, store::AbstractString)
    m = read_json(manifest)
    items = NamedTuple[]
    for c in m["chunks"]
        append!(items, _p2_chunk_items(read_chunk(chunk_path(store, c["ordinal"])), c["ordinal"]))
    end
    return items
end

"The item with these properties (`direction = nothing` for an item without one), or `nothing`."
function p2_find_item(items::AbstractVector; orbit, id, precision, horizon, tau, direction)
    i = findfirst(it -> it.orbit == orbit && it.id == id && it.precision == precision && it.horizon == horizon && it.tau == tau && it.direction == direction, items)
    return i === nothing ? nothing : items[i]
end

"Whether the bit of flag `name` is set in the word of row `row` of a window item; the bit of the j-th name has the value 2^(j-1) (P34)."
function p2_flag(item, row::Integer, name::AbstractString)
    j = findfirst(==(name), item.bits)
    j === nothing && throw(ArgumentError("item of orbit $(item.orbit), spec $(item.id) has no flag named $name"))
    return (item.flagwords[row] & (Int64(1) << (j - 1))) != 0
end

# Dictionary of items by content, for the lookups of the pipeline; a repeated content refuses.
p2_key(orbit, id, precision, horizon, tau, direction) = (Int(orbit), String(id), String(precision), Int(horizon), Int(tau), direction === nothing ? nothing : String(direction))

function p2_index(items)
    idx = Dict{Tuple, Any}()
    for it in items
        k = p2_key(it.orbit, it.id, it.precision, it.horizon, it.tau, it.direction)
        haskey(idx, k) && throw(ArgumentError("item stored twice: orbit $(it.orbit), spec $(it.id), $(it.precision), horizon $(it.horizon), direction $(it.direction)"))
        idx[k] = it
    end
    return idx
end

# The first id of the grid orbits: given, or one above the largest id of the panel file, or (no panel known) that of the panel file next to the preregistration, or 884.
function _p2_first_id(prereg, first_id, panel)
    first_id === nothing || return Int(first_id)
    if panel === nothing
        sib = filter(f -> endswith(f, "-panel.toml"), readdir(dirname(abspath(prereg))))
        panel = length(sib) == 1 ? joinpath(dirname(abspath(prereg)), only(sib)) : nothing
    end
    panel === nothing && return 884
    return maximum(Int(o["id"]) for o in TOML.parsefile(panel)["orbit"]) + 1
end

"""
    p2_grid_starts(prereg; first_id = nothing, panel = nothing) -> Vector{NamedTuple}

The starts of grid G in the order of the `[[run.h4]]` tables of the preregistration (E24): `(id, label, kindex, i, j, K, x0, p0)`, the three numbers as the texts of the
tables. The ids follow the panel (the first id is `first_id`, or one above the largest panel id; 884 without a panel). A label that is not `grid-k<K>-i<i>-j<j>` refuses.
"""
function p2_grid_starts(prereg::AbstractString; first_id = nothing, panel = nothing)
    h4 = read_run_block(prereg).raw["h4"]
    id0 = _p2_first_id(prereg, first_id, panel)
    out = NamedTuple[]
    for (n, h) in enumerate(h4)
        m = match(r"^grid-k(\d+)-i(\d+)-j(\d+)$", h["label"])
        m === nothing && throw(ArgumentError("run block h4 label \"$(h["label"])\" is not of the form grid-k<K index>-i<i>-j<j>"))
        push!(out, (id = id0 + n - 1, label = h["label"], kindex = parse(Int, m.captures[1]), i = parse(Int, m.captures[2]), j = parse(Int, m.captures[3]),
                    K = h["K"], x0 = h["x0"], p0 = h["p0"]))
    end
    return out
end
