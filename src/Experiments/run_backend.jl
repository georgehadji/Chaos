# Backend resolution (R4) and launch groups (R3).

"True for the methods whose reducer chain has a CPU-only stage: the library Householder path (L6)."
cpu_only_method(id::AbstractString) = id == "L6"

_element_type(label::AbstractString) = first(parse_precision(label == "BigFloat" ? "BigFloat(256)" : label))

"""
    resolve_run_backend(request, T, method; threads, record) -> (backend, device, reason, note)

Resolve a request (`:auto`, `:gpu`, `:cpu`) for one (element type, method) by the rules of `_config/architecture.md`, "Compute",
with the K19 gate: `:auto` takes the device only when `record` (`nothing` or `(threads, gpu_faster)`) was taken at the thread
count of the run and says the device is faster. No tenth reason is added; a record at another thread count gives the existing
reason "GPU not faster for the element type" and a `note`.
"""
function resolve_run_backend(request::Symbol, ::Type{T}, method::AbstractString; threads::Integer, record = nothing) where {T}
    note = request === :auto && record !== nothing && record.threads != threads ?
           "device record taken at $(record.threads) threads, the run uses $threads threads (K19)" : ""
    function probe()
        p = Foundation.probe(Foundation.DEVICE[], T)
        return merge(p, (gpu_faster = record !== nothing && record.threads == threads && record.gpu_faster,))
    end
    c = resolve_backend(request, T; cpu_only_reducer = cpu_only_method(method), probe = probe)
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
        r = resolve_run_backend(req, _element_type(key[1]), it.method; threads = threads, record = record)
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
