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

# A file held for a moment by another process (a reader, a virus scanner) fails with one of these codes on Windows.
const RETRY_CODES = (Base.UV_EBUSY, Base.UV_EACCES, Base.UV_EPERM)

const RETRY_TRIES = 12            # about 6.6 s in all with the pause below (D59)
const RETRY_PAUSE = Ref(0.1)      # seconds, times the attempt number; a test seam
const BACKUP_COPIER = Ref{Any}(cp)   # copies a file for the numbered backups; a test seam

"Call `f()`, trying again when another process holds a file (the codes of `RETRY_CODES`); any other error, and the last failure, are raised."
function _retry(f; tries::Integer = RETRY_TRIES, pause::Real = RETRY_PAUSE[])
    for k in 1:tries
        try
            return f()
        catch e
            (e isa Base.IOError && e.code in RETRY_CODES && k < tries) || rethrow()
            sleep(pause * k)
        end
    end
end

"Rename `src` to `dst` (replacing a file on every platform), trying again when another process holds a file. `mover` is a test seam."
_move_retry(src, dst; mover = Base.Filesystem.rename, tries::Integer = RETRY_TRIES, pause::Real = RETRY_PAUSE[]) =
    _retry(() -> mover(src, dst); tries = tries, pause = pause)

# Rotate `path.1` .. `path.keep`: the newest copy of the file about to be replaced is `path.1`. A copy that stays held after the retries is
# skipped (a backup is a convenience, D59): a warning is printed and the text is pushed to `failures`. Other errors are raised.
function _rotate_backups(path, keep; failures = nothing)
    function copy1(a, b)
        try
            _retry(() -> BACKUP_COPIER[](a, b; force = true))
        catch e
            (e isa Base.IOError && e.code in RETRY_CODES) || rethrow()
            msg = "backup $(basename(b)) skipped: $(sprint(showerror, e))"
            @warn msg
            failures === nothing || push!(failures, msg)
        end
    end
    for k in (keep - 1):-1:1
        isfile("$path.$k") && copy1("$path.$k", "$path.$(k + 1)")
    end
    isfile(path) && copy1(path, "$path.1")
end

# Write of a JSON file: temporary name next to it, then an atomic replace with retries. The previous file stays whole when the
# replace fails; `keep > 0` also keeps that many numbered backups of the previous versions (best effort, D59). `writer` is a test seam.
function _write_json(path::AbstractString, obj; keep::Integer = 0, mover = Base.Filesystem.rename, failures = nothing, writer = write)
    mkpath(dirname(abspath(path)))
    tmp = path * ".partial"
    text = JSON.json(_jsonable(obj))
    _retry(() -> writer(tmp, text))
    keep > 0 && _rotate_backups(path, keep; failures = failures)
    _move_retry(tmp, path; mover = mover)
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
    _put_indicators!(g, o)
end

_bits(x::AbstractArray{Float64}) = reinterpret(UInt64, vec(x))
_bits(x::AbstractArray) = vec(x)

_same(a::AbstractArray, b::AbstractArray) = size(a) == size(b) && _bits(a) == _bits(b)

function _same(a::OrbitRecord, b::OrbitRecord)
    a.id == b.id && isequal(a.attrs, b.attrs) && _same(a.t, b.t) && _same(a.state, b.state) && keys(a.lyap) == keys(b.lyap) &&
        (a.winding === nothing ? b.winding === nothing : b.winding !== nothing && _same(a.winding, b.winding)) &&
        all(m -> _same(a.lyap[m], b.lyap[m]), keys(a.lyap)) && _same_indicators(a.indicators, b.indicators) && keys(a.diagnostics) == keys(b.diagnostics) &&
        all(m -> _same(a.diagnostics[m], b.diagnostics[m]) && a.thresholds[m] === b.thresholds[m] && a.first_exceedance[m] == b.first_exceedance[m],
            keys(a.diagnostics)) && all(m -> _same(a.cost[m], b.cost[m]), keys(a.cost))
end

_nohook(point) = nothing

# Flush the file to the disk before it gets its final name (a power cut must not leave a final name with lost data).
function _fsync(path)
    open(path, "r+") do io
        fd = Base.fd(io)
        fd isa RawFD || (fd = RawFD(fd))   # Julia 1.11 returns an Int, 1.13 a RawFD
        r = Sys.iswindows() ? ccall(:_commit, Cint, (RawFD,), fd) : ccall(:fsync, Cint, (RawFD,), fd)
        r == 0 || error("write_chunk: fsync of $path failed ($r)")
    end
end

# Names of the root attributes whose read-back differs from what was written.
_attr_diff(want, got) = sort!(String[k for k in union(keys(want), keys(got)) if !(haskey(want, k) && haskey(got, k) && isequal(want[k], got[k]))])

# A final name that is taken is moved to `<name>.old<k>` first: the rename itself never overwrites.
function _take_name!(path)
    ispath(path) || return
    k = 1
    while ispath("$path.old$k")
        k += 1
    end
    mv(path, "$path.old$k")
end

"""
    write_chunk(path, attrs, orbits)

Write a chunk file atomically: under `path * ".tmp"` with `complete` false, orbit groups, close, read back and compare bit for bit (root
attributes included), set `complete`, flush to the disk, rename to `path` (a file already there is moved to `<path>.old<k>` first; the
rename retries while another process holds a file). `attrs` are the root attributes, `orbits` the orbit records.
"""
function write_chunk(path::AbstractString, attrs::AbstractDict, orbits::AbstractVector; hook = _nohook)
    tmp = path * ".tmp"
    mkpath(dirname(abspath(path)))
    want = merge(Dict{String, Any}(attrs), Dict{String, Any}("complete" => false))
    prov = pop!(want, "provenance", Dict{String, Any}())      # the reducers' own keys, one group per member (D86, R15)
    HDF5.h5open(tmp, "w") do f
        _put_attrs!(f, want)
        _put_provenance!(f, prov)
        HDF5.create_group(f, "orbits")
        hook(:before_data)
        half = cld(length(orbits), 2)
        for (j, o) in enumerate(orbits)
            _put_orbit!(f, o)
            j == half && hook(:during_data)
        end
        _put_pooled!(f, orbits)
    end
    back = read_chunk(tmp)
    length(back.orbits) == length(orbits) && all(_same(a, b) for (a, b) in zip(orbits, back.orbits)) ||
        error("write_chunk: $tmp differs from the buffers after the write")
    isequal(pooled_of(orbits), back.pooled) || error("write_chunk: pooled counts of $tmp differ from the buffers after the write")
    isequal(prov, get(back.attrs, "provenance", Dict{String, Any}())) || error("write_chunk: provenance of $tmp differs from the buffers after the write")
    bad = _attr_diff(want, Dict{String, Any}(k => v for (k, v) in back.attrs if k != "provenance"))
    isempty(bad) || error("write_chunk: root attributes of $tmp differ from the buffers after the write: " * join(bad, ", "))
    hook(:after_data)
    HDF5.h5open(tmp, "r+") do f
        HDF5.delete_attribute(f, "complete")
        HDF5.write_attribute(f, "complete", true)
    end
    read_chunk_attrs(tmp)["complete"] === true || error("write_chunk: $tmp does not read back complete")
    _fsync(tmp)
    hook(:before_rename)
    _take_name!(path)
    _move_retry(tmp, path)
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
    cg = HDF5.open_group(g, "cost")
    lyap = haskey(g, "lyap") ? (lg = HDF5.open_group(g, "lyap"); Dict{String, Matrix{Float64}}(m => read(HDF5.open_dataset(lg, m)) for m in keys(lg))) :
           Dict{String, Matrix{Float64}}()        # the items of the indicator specs have no lyap group (D86)
    diag, thr, fe = _read_diag(g)
    cost = Dict{String, Vector{Float64}}(m => read(HDF5.open_dataset(cg, m)) for m in keys(cg))
    wind = haskey(g, "winding") ? read(HDF5.open_dataset(g, "winding")) : nothing
    return OrbitRecord(id, _read_attrs(g), read(HDF5.open_dataset(g, "t")), read(HDF5.open_dataset(g, "state")), wind, lyap, diag, thr, fe, cost, "",
                       _read_indicators(g))
end

"Read a chunk file: `(attrs, orbits)`, the orbit records in id order."
function read_chunk(path::AbstractString)
    isfile(path) || throw(ArgumentError("no chunk file $path"))
    HDF5.h5open(String(path), "r") do f
        og = HDF5.open_group(f, "orbits")
        ids = sort!([parse(Int, k) for k in keys(og)])
        attrs = _read_attrs(f)
        haskey(f, "provenance") && (attrs["provenance"] = _read_provenance(f))
        return (attrs = attrs, orbits = OrbitRecord[_read_orbit(HDF5.open_group(og, string(id)), id) for id in ids], pooled = _read_pooled(f))
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
