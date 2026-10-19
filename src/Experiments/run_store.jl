# The run store (R6, R10; schema/run-store.md, version 2): one HDF5 file per chunk, written under a temporary name, read back and
# renamed, so that no file with a final name is ever incomplete (I13). Also the small JSON helpers of the manifest and the report.

import HDF5
import JSON

# ---------------------------------------------------------------- JSON
_jsonable(x::AbstractDict) = Dict{String, Any}(String(k) => _jsonable(v) for (k, v) in x)
_jsonable(x::NamedTuple) = Dict{String, Any}(String(k) => _jsonable(v) for (k, v) in pairs(x))
_jsonable(x::AbstractVector) = Any[_jsonable(v) for v in x]
_jsonable(x::Tuple) = Any[_jsonable(v) for v in x]
_jsonable(x::Symbol) = String(x)
_jsonable(x::AbstractFloat) = isfinite(x) ? Float64(x) : string(x)
_jsonable(x::Union{Bool, Integer, AbstractString, Nothing}) = x isa AbstractString ? String(x) : x
_jsonable(x) = string(x)

"Parsed JSON: objects as `Dict{String,Any}`, integers as `Int`, other numbers as `Float64`."
read_json(path::AbstractString) = JSON.parse(read(path, String); dicttype = Dict{String, Any})

# Atomic write of a JSON file: temporary name next to it, then rename.
function _write_json(path::AbstractString, obj)
    mkpath(dirname(abspath(path)))
    tmp = path * ".partial"
    write(tmp, JSON.json(_jsonable(obj)))
    mv(tmp, path; force = true)
    return path
end

# ---------------------------------------------------------------- writing
_put_attrs!(obj, attrs) = foreach(((k, v),) -> HDF5.write_attribute(obj, k, v), pairs(attrs))

function _put_diag!(g, name, v, thr, first_step)
    d = HDF5.create_dataset(g, "diagnostics/" * name, Float64, size(v))
    write(d, v)
    HDF5.write_attribute(d, "threshold", thr)
    HDF5.write_attribute(d, "first_exceedance_step", first_step)
end

function _put_orbit!(f, o::OrbitRecord)
    g = HDF5.create_group(f, "orbits/$(o.id)")
    _put_attrs!(g, o.attrs)
    g["t"] = o.t
    g["state"] = o.state
    o.winding === nothing || (g["winding"] = o.winding)
    for (m, v) in o.lyap
        g["lyap/" * m] = v
    end
    for (nm, v) in o.diagnostics
        _put_diag!(g, nm, v, o.thresholds[nm], o.first_exceedance[nm])
    end
    for (m, c) in o.cost
        g["cost/" * m] = c
    end
end

_bits(x::AbstractArray{Float64}) = reinterpret(UInt64, vec(x))
_bits(x::AbstractArray) = vec(x)

_same(a::AbstractArray, b::AbstractArray) = size(a) == size(b) && _bits(a) == _bits(b)

function _same(a::OrbitRecord, b::OrbitRecord)
    a.id == b.id && a.attrs == b.attrs && _same(a.t, b.t) && _same(a.state, b.state) && keys(a.lyap) == keys(b.lyap) &&
        (a.winding === nothing ? b.winding === nothing : b.winding !== nothing && _same(a.winding, b.winding)) &&
        all(m -> _same(a.lyap[m], b.lyap[m]), keys(a.lyap)) && keys(a.diagnostics) == keys(b.diagnostics) &&
        all(m -> _same(a.diagnostics[m], b.diagnostics[m]) && a.thresholds[m] === b.thresholds[m] && a.first_exceedance[m] == b.first_exceedance[m],
            keys(a.diagnostics)) && all(m -> _same(a.cost[m], b.cost[m]), keys(a.cost))
end

_nohook(point) = nothing

"""
    write_chunk(path, attrs, orbits)

Write a chunk file atomically: under `path * ".tmp"` with `complete` false, orbit groups, close, read back and compare bit for bit,
set `complete`, rename to `path`. `attrs` are the root attributes, `orbits` the orbit records.
"""
function write_chunk(path::AbstractString, attrs::AbstractDict, orbits::AbstractVector; hook = _nohook)
    tmp = path * ".tmp"
    mkpath(dirname(abspath(path)))
    HDF5.h5open(tmp, "w") do f
        _put_attrs!(f, merge(Dict{String, Any}(attrs), Dict{String, Any}("complete" => false)))
        HDF5.create_group(f, "orbits")
        hook(:before_data)
        half = cld(length(orbits), 2)
        for (j, o) in enumerate(orbits)
            _put_orbit!(f, o)
            j == half && hook(:during_data)
        end
    end
    back = read_chunk(tmp)
    length(back.orbits) == length(orbits) && all(_same(a, b) for (a, b) in zip(orbits, back.orbits)) ||
        error("write_chunk: $tmp differs from the buffers after the write")
    hook(:after_data)
    HDF5.h5open(tmp, "r+") do f
        HDF5.delete_attribute(f, "complete")
        HDF5.write_attribute(f, "complete", true)
    end
    hook(:before_rename)
    mv(tmp, path; force = true)
    hook(:after_rename)
    return path
end

# ---------------------------------------------------------------- reading
_read_attrs(obj) = Dict{String, Any}(k => HDF5.read_attribute(obj, k) for k in keys(HDF5.attributes(obj)))

# Diagnostics of one orbit group (none for L2): running maxima, thresholds, first exceedance steps.
function _read_diag(g::HDF5.Group)
    haskey(g, "diagnostics") || return Dict{String, Vector{Float64}}(), Dict{String, Float64}(), Dict{String, Int}()
    dg = HDF5.open_group(g, "diagnostics")
    dsets = Dict{String, HDF5.Dataset}(nm => HDF5.open_dataset(dg, nm) for nm in keys(dg))
    return (Dict{String, Vector{Float64}}(nm => read(d) for (nm, d) in dsets),
            Dict{String, Float64}(nm => Float64(HDF5.read_attribute(d, "threshold")) for (nm, d) in dsets),
            Dict{String, Int}(nm => Int(HDF5.read_attribute(d, "first_exceedance_step")) for (nm, d) in dsets))
end

function _read_orbit(g::HDF5.Group, id::Int)
    lg, cg = HDF5.open_group(g, "lyap"), HDF5.open_group(g, "cost")
    lyap = Dict{String, Matrix{Float64}}(m => read(HDF5.open_dataset(lg, m)) for m in keys(lg))
    diag, thr, fe = _read_diag(g)
    cost = Dict{String, Vector{Float64}}(m => read(HDF5.open_dataset(cg, m)) for m in keys(cg))
    wind = haskey(g, "winding") ? read(HDF5.open_dataset(g, "winding")) : nothing
    return OrbitRecord(id, _read_attrs(g), read(HDF5.open_dataset(g, "t")), read(HDF5.open_dataset(g, "state")), wind, lyap, diag, thr, fe, cost, "")
end

"Read a chunk file: `(attrs, orbits)`, the orbit records in id order."
function read_chunk(path::AbstractString)
    isfile(path) || throw(ArgumentError("no chunk file $path"))
    HDF5.h5open(String(path), "r") do f
        og = HDF5.open_group(f, "orbits")
        ids = sort!([parse(Int, k) for k in keys(og)])
        return (attrs = _read_attrs(f), orbits = OrbitRecord[_read_orbit(HDF5.open_group(og, string(id)), id) for id in ids])
    end
end

# Root attributes only (cheap scan of a store).
read_chunk_attrs(path::AbstractString) = HDF5.h5open(f -> _read_attrs(f), path, "r")

"""
    chunk_layout(path) -> Dict

Dataset path => `(eltype, size)` of every dataset of a chunk file, as stored (element type name, dimensions).
"""
function chunk_layout(path::AbstractString)
    out = Dict{String, NamedTuple}()
    walk(obj, prefix) = for k in keys(obj)
        x = obj[k]
        x isa HDF5.Dataset ? (out[prefix * "/" * k] = (eltype = string(eltype(x)), size = size(x))) : walk(x, prefix * "/" * k)
    end
    HDF5.h5open(f -> walk(f, ""), path, "r")
    return out
end
