# Reader of the run store (E1, E2): one item per (orbit, chunk), found by content.

chunk_path(store::AbstractString, ordinal::Integer) = joinpath(store, "chunk-" * lpad(string(ordinal), 5, '0') * ".h5")

_item_eps(a) = haskey(a, "perturbation") ? parse(Float64, a["perturbation"]) : nothing

"Items of one chunk file: the first exponent, the orbit point, the winding, flags, attributes, running maxima and cost of each orbit."
function _chunk_items(ch, ordinal::Integer)
    a = ch.attrs
    m = a["method"]
    return [(orbit = o.id, method = m, precision = a["precision"], horizon = Int(a["horizon"]), tau = Int(a["tau"]), eps = _item_eps(a),
             direction = a["direction"], mode = a["coordinate_mode"], arms = String.(a["arms"]), chunk = Int(ordinal), t = o.t,
             lam = o.lyap[m][:, 1], state = o.state, winding = o.winding, flags = o.attrs["flags"], attrs = o.attrs,
             diag = o.diagnostics, thresholds = o.thresholds, cost = get(o.cost, m, Float64[])) for o in ch.orbits]
end

"""
    read_store(manifest, store) -> Vector

Every item of the chunk files the manifest lists (axis order of `schema/run-store.md`, "Axis order": the exponent of step j is `lam[j]`).
"""
function read_store(manifest::AbstractString, store::AbstractString)
    m = read_json(manifest)
    items = NamedTuple[]
    for c in m["chunks"]
        append!(items, _chunk_items(read_chunk(chunk_path(store, c["ordinal"])), c["ordinal"]))
    end
    return items
end

item_key(orbit, method, precision, horizon, tau, eps, direction, mode) = (Int(orbit), String(method), String(precision), Int(horizon), Int(tau), eps, String(direction), String(mode))
item_key(it) = item_key(it.orbit, it.method, it.precision, it.horizon, it.tau, it.eps, it.direction, it.mode)

"The item with these properties, or `nothing`."
function find_item(items::AbstractVector; orbit, method, precision, horizon, tau, eps, direction, mode)
    k = item_key(orbit, method, precision, horizon, tau, eps, direction, mode)
    i = findfirst(it -> item_key(it) == k, items)
    return i === nothing ? nothing : items[i]
end

"Dictionary of items by key, for lookups in the pipeline."
index_items(items) = Dict{Tuple, Any}(item_key(it) => it for it in items)
