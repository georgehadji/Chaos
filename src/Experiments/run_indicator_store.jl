# Store layout of the indicator specs (Amendment D86, R15; schema/run-store.md, version 2, additive): the datasets under `indicators/<id>/`, the
# provenance groups of the chain, the pooled counts of a chunk and the root attributes of a chunk of a spec.

function _put_ds!(g, path, data, attrs...)
    d = HDF5.create_dataset(g, path, eltype(data), size(data))
    write(d, data)
    for (k, v) in attrs
        HDF5.write_attribute(d, k, v)
    end
    return nothing
end

# One orbit's tables: `series` (n x c) with `columns`; `final` (rows x columns) with `columns` and `rows`, and `flags` with `bits`; `scalars` with `columns`
# (and `range` for stretching); `counts` for the orbits of the named subset.
function _put_indicators!(g, o::OrbitRecord)
    for (id, t) in o.indicators
        base = "indicators/$id/"
        t.series === nothing || _put_ds!(g, base * "series", t.series, "columns" => t.columns)
        if t.final !== nothing
            _put_ds!(g, base * "final", t.final, "columns" => t.fcolumns, "rows" => t.rows)
            _put_ds!(g, base * "flags", t.flags, "bits" => t.bits)
        end
        t.scalars === nothing || _put_ds!(g, base * "scalars", t.scalars, "columns" => t.scolumns, (t.range === nothing ? () : ("range" => t.range,))...)
        t.counts === nothing || _put_ds!(g, base * "counts", t.counts)
    end
    return nothing
end

function _put_provenance!(f, prov)
    for (m, d) in prov
        _put_attrs!(HDF5.create_group(f, "provenance/$m"), d)
    end
    return nothing
end

_read_provenance(f) = (pg = HDF5.open_group(f, "provenance"); Dict{String, Any}(m => _read_attrs(HDF5.open_group(pg, m)) for m in keys(pg)))

# The pooled counts of a chunk (E29, D88): the exact integer sum over the orbits of the count of each cell at the last grid step, stored once per spec with the
# range of the cells. Orbits that do not share one range (derived ranges of different K or of the linear control) have no common cells: the group then holds
# the number of orbits and nothing else.
function pooled_of(orbits)
    out = Dict{String, NamedTuple}()
    ids = sort!(unique(id for o in orbits for (id, t) in o.indicators if t.pool !== nothing))
    for id in ids
        ts = [o.indicators[id] for o in orbits if haskey(o.indicators, id) && o.indicators[id].pool !== nothing]
        shared = all(t -> isequal(t.range, first(ts).range), ts)
        out[id] = (orbits = length(ts), counts = shared ? reduce(+, (Int64.(t.pool) for t in ts)) : nothing, range = shared ? first(ts).range : nothing)
    end
    return out
end

function _put_pooled!(f, orbits)
    for (id, p) in sort!(collect(pooled_of(orbits)); by = first)
        g = HDF5.create_group(f, "pooled/$id")
        HDF5.write_attribute(g, "orbits", p.orbits)
        p.counts === nothing && continue
        _put_ds!(g, "counts", p.counts)
        HDF5.write_attribute(g, "range", p.range)
    end
    return nothing
end

# The pooled groups of a file, in the form of `pooled_of` (write_chunk compares the two after the write).
function _read_pooled(f)
    out = Dict{String, NamedTuple}()
    haskey(f, "pooled") || return out
    pg = HDF5.open_group(f, "pooled")
    for id in keys(pg)
        g = HDF5.open_group(pg, id)
        out[id] = (orbits = Int(HDF5.read_attribute(g, "orbits")), counts = haskey(g, "counts") ? read(HDF5.open_dataset(g, "counts")) : nothing,
                   range = haskey(HDF5.attributes(g), "range") ? HDF5.read_attribute(g, "range") : nothing)
    end
    return out
end

function _read_ds(g, name)
    haskey(g, name) || return nothing, Dict{String, Any}()
    d = HDF5.open_dataset(g, name)
    return read(d), _read_attrs(d)
end

function _read_indicators(g::HDF5.Group)
    out = Dict{String, NamedTuple}()
    haskey(g, "indicators") || return out
    ig = HDF5.open_group(g, "indicators")
    for id in keys(ig)
        s = HDF5.open_group(ig, id)
        series, sa = _read_ds(s, "series")
        final, fa = _read_ds(s, "final")
        flags, ga = _read_ds(s, "flags")
        scalars, ca = _read_ds(s, "scalars")
        counts, _ = _read_ds(s, "counts")
        out[id] = _ind_nt(series = series, columns = get(sa, "columns", String[]), final = final, rows = get(fa, "rows", String[]), fcolumns = get(fa, "columns", String[]),
                          flags = flags, bits = get(ga, "bits", String[]), scalars = scalars, scolumns = get(ca, "columns", String[]), range = get(ca, "range", nothing),
                          counts = counts)
    end
    return out
end

# Bit-for-bit equality of the stored fields of two indicator results (`pool` is in memory only).
_same_field(a::Nothing, b::Nothing) = true
_same_field(a, b::Nothing) = false
_same_field(a::Nothing, b) = false
_same_field(a::AbstractArray, b::AbstractArray) = _same(a, b)
_same_field(a, b) = a == b
function _same_indicators(a::AbstractDict, b::AbstractDict)
    keys(a) == keys(b) || return false
    return all(id -> all(f -> _same_field(getfield(a[id], f), getfield(b[id], f)), (:series, :columns, :final, :rows, :fcolumns, :flags, :bits, :scalars, :scolumns, :range, :counts)), keys(a))
end

# ---------------------------------------------------------------- root attributes of a chunk of a spec
# Only the types HDF5 writes natively go through; anything else (Double64, BigFloat, BigInt, a tuple) is stored as its exact text.
_h5val(v::AbstractString) = String(v)
_h5val(v::Union{Bool, Base.BitInteger, Float32, Float64}) = v
_h5val(v) = string(v)

# The reducers' own provenance, kept apart per member so that no key is overwritten (E23, correction 4 of the D86 stage 02 checkpoint).
function ind_provenance(::Type{T}, item, spec) where {T}
    kind = spec["kind"]
    sys = first(_system_start(T, item))
    if kind in ("tangent", "pair")
        return Dict{String, Any}(m => Dict{String, Any}(String(k) => _h5val(v) for (k, v) in pairs(reducer_provenance(r)))
                                 for (m, r) in zip(spec["members"], chain_members(spec["members"])))
    end
    kind == "frequency" && return Dict{String, Any}("frequency" => Dict{String, Any}("estimator" => "windowed-fourier", "weight" => "cosine", "order" => spec["order"],
                                    "positions" => spec["positions"], "pad" => spec["pad"], "mean_removal" => spec["mean_removal"]))
    r = first(ind_readers(T, item, spec, sys))
    return Dict{String, Any}(ind_members(spec)[1] => Dict{String, Any}(String(k) => _h5val(v) for (k, v) in pairs(reducer_provenance(r))))
end

"Root attributes of a chunk file of an indicator spec (R15), from the first item of the chunk."
function ind_root_attrs(plan::RunPlan, it)
    spec = block_spec(plan.block, it.method)
    kind = spec["kind"]
    a = Dict{String, Any}("kind" => kind, "spec" => spec_text(spec), "spec_sha256" => it.spec, "chain" => chain_record(spec), "members" => ind_members(spec),
                          "provenance" => ind_provenance(_ind_prec_type(it.precision), it, spec))
    kind == "birkhoff" && (a["observable_names"] = ind_observable_names(spec, it))
    ind_window_kind(kind) && (a["windows"] = Int.(spec["windows"]); a["run_length"] = ind_run_length(spec))
    if kind == "pair"
        a["frame_source"], a["direction_rule"], a["separation"] = spec["frame_source"], IND_DIRECTION_RULE, spec["separation"]
    end
    return a
end
