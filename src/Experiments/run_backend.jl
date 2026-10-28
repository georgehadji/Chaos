# Backend resolution (R4) and launch groups (R3).

"True for the methods whose reducer chain has a CPU-only stage: the library Householder path (L6)."
cpu_only_method(id::AbstractString) = id == "L6"

"Also true for the indicator kinds without a device path in Amendment D86: stretching (the final builds count matrices) and frequency (the series leaves the device)."
cpu_only_method(id::AbstractString, kind) = cpu_only_method(id) || kind in ("stretching", "frequency")

_element_type(label::AbstractString) = first(parse_precision(label == "BigFloat" ? "BigFloat(256)" : label))

"""
    method_record(record, method) -> record or nothing

D54: a record with a `recorder` table (method id => record on the recording chain of that method) is read per method, and a method
without an entry has no record; a record without the field is the D53h form, one record for every method.
"""
method_record(record, method::AbstractString) = record === nothing || !hasproperty(record, :recorder) ? record : get(record.recorder, method, nothing)

"""
    resolve_run_backend(request, T, method; threads, record, kind, spec_sha256) -> (backend, device, reason, note)

Resolve a request (`:auto`, `:gpu`, `:cpu`) for one (element type, method) by the rules of `_config/architecture.md`, "Compute",
with the K19 gate: `:auto` takes the device only when `record` (`nothing` or `(threads, gpu_faster)`) was taken at the thread
count of the run and says the device is faster. No tenth reason is added; a record at another thread count gives the existing
reason "GPU not faster for the element type" and a `note`. For an indicator spec (`kind` given, D86) the record is keyed by the spec id and carries the
`spec_sha256` of the chain it was taken on: a record of another chain gives the same reason with the note "record taken on another chain".
"""
function resolve_run_backend(request::Symbol, ::Type{T}, method::AbstractString; threads::Integer, record = nothing, kind = nothing,
                             spec_sha256 = nothing) where {T}
    record = method_record(record, method)
    other = kind !== nothing && record !== nothing && (hasproperty(record, :spec_sha256) ? record.spec_sha256 : "") != spec_sha256
    note = request === :auto && record !== nothing && record.threads != threads ?
           "device record taken at $(record.threads) threads, the run uses $threads threads (K19)" :
           request === :auto && other ? "record taken on another chain (spec_sha256 differs from the spec hash)" : ""
    function probe()
        p = Foundation.probe(Foundation.DEVICE[], T)
        return merge(p, (gpu_faster = record !== nothing && record.threads == threads && record.gpu_faster && !other,))
    end
    c = resolve_backend(request, T; cpu_only_reducer = cpu_only_method(method, kind), probe = probe)
    return (backend = c.backend, device = c.device, reason = c.reason, note = note)
end

"""
    backend_table(plan; threads, record) -> Dict

One row per (precision label, method) of the plan: `(request, backend, device, reason, note)`; the request is the block's, per precision.
"""
function backend_table(plan::RunPlan; threads::Integer, record = nothing)
    tbl = Dict{Tuple{String, String}, NamedTuple}()
    for it in plan.items
        key = (_label(it.precision), it.method)
        haskey(tbl, key) && continue
        req = Symbol(plan.block.raw["backend"][key[1]])
        r = resolve_run_backend(req, _element_type(key[1]), it.method; threads = threads, record = record, kind = it.kind, spec_sha256 = it.spec)
        tbl[key] = (request = req, r...)
    end
    return tbl
end

"""
    launch_capacity(free_bytes; f, M, reducer, columns, grid_steps) -> Int

Orbits per launch: `floor(f free_bytes / (M + reducer + 8 columns grid_steps))`, exact arithmetic.
"""
function launch_capacity(free_bytes::Integer; f::Real, M::Integer, reducer::Integer, columns::Integer, grid_steps::Integer)
    per = M + reducer + 8 * columns * grid_steps
    r = Rational{BigInt}(f)
    return Int(fld(numerator(r) * big(free_bytes), denominator(r) * per))
end

"""
    launch_groups(plan; capacity) -> Vector{(chunks, orbits, launches)}

Consecutive chunks of one run type, up to `capacity` orbits each; a chunk larger than the capacity is a group of its own run in
`ceil(orbits / capacity)` launches. Grouping never changes a result (I12).
"""
function launch_groups(plan::RunPlan; capacity::Integer)
    cap = max(1, Int(capacity))
    groups = NamedTuple{(:chunks, :orbits, :launches), Tuple{Vector{Int}, Int, Int}}[]
    cur, n = Int[], 0
    flush() = isempty(cur) || push!(groups, (chunks = cur, orbits = n, launches = cld(n, cap)))
    for c in plan.chunks
        m = length(c.items)
        if !isempty(cur) && (plan.chunks[last(cur)].run_type != c.run_type || n + m > cap)
            flush()
            cur, n = Int[], 0
        end
        push!(cur, c.ordinal)
        n += m
    end
    flush()
    return groups
end
