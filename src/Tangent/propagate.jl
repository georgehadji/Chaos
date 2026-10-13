"""
    bytes_per_orbit(Val(d), Val(k), T) -> Int

Resident bytes per orbit: `sizeof(T) (d + d k + 2 k) + 16` (spec Compute).
"""
bytes_per_orbit(::Val{d}, ::Val{k}, ::Type{T}) where {d, k, T} = sizeof(T) * (d + d * k + 2k) + 16

"""
    chunk_size(free_bytes, Val(d), Val(k), T; f=0.5, reducer_bytes=0) -> Int

Orbits per chunk: `floor(f free_bytes / (bytes_per_orbit + reducer_bytes))`, at least 1.
"""
function chunk_size(free_bytes::Integer, v::Val{d}, kv::Val{k}, ::Type{T}; f = 0.5, reducer_bytes::Integer = 0) where {d, k, T}
    return max(1, floor(Int, f * free_bytes / (bytes_per_orbit(v, kv, T) + reducer_bytes)))
end

"Result of [`propagate`](@ref): final points `x`, deviation matrices `W`, per-orbit reducer summaries `finals`, run `provenance`."
struct PropagationResult{X, M, F}
    x::X
    W::M
    finals::F
    provenance::NamedTuple
end

# One work item per orbit; the body is the same `_orbit_body` as `propagate_orbit`.
@kernel function _orbit_kernel!(xo, Wo, so, code, hstep, sys, x0, W0, N, tau, reducers)
    i = @index(Global)
    x, W, st, c, h = _orbit_body(sys, x0[i], W0[i], N, tau, reducers)
    xo[i] = x
    Wo[i] = W
    so[i] = st
    code[i] = c
    hstep[i] = h
end

const WORKGROUP = 16   # orbits per CPU work item; the sequential work per orbit is long, so groups stay small

# Initial deviation matrices of one chunk, on the host for every backend (spec, initial vectors).
function _fill_W!(Win, r, seed, W0, ::Val{d}, ::Val{k}, ::Type{T}) where {d, k, T}
    for (i, j) in enumerate(r)
        Win[i] = W0 === nothing ? initial_deviations(seed, j, Val(d), Val(k), T) : SMatrix{d, k, T}(W0[j])
    end
end

function _run_chunk!(out, bufs, r, args, choice)
    xin, Win, xo, Wo, so, code, hstep = bufs
    len = length(r)
    copyto!(xin, 1, args.x0, first(r), len)
    _fill_W!(Win, r, args.seed, args.W0, args.dims..., args.T)
    backend = CPU()
    _orbit_kernel!(backend, WORKGROUP)(xo, Wo, so, code, hstep, args.sys, xin, Win, args.N, args.tau, args.reducers; ndrange = len)
    KernelAbstractions.synchronize(backend)
    i = findfirst(!=(STATUS_OK), view(code, 1:len))   # lowest orbit id first; a failed chunk stores nothing
    i === nothing || code[i] != STATUS_REDUCER || throw(_reducer_error(args.reducers, so[i], first(r) + i - 1, choice.backend))
    i === nothing || throw(_health_error(code[i], xo[i], Wo[i], hstep[i], first(r) + i - 1, choice.backend))
    copyto!(out.x, first(r), xo, 1, len)
    copyto!(out.W, first(r), Wo, 1, len)
    for i in 1:len
        out.finals[first(r) + i - 1] = map(reducer_final, args.reducers, so[i])
    end
end

function _propagate(sys, x0::AbstractVector{SVector{d, T}}, N, tau, ::Val{k}, seed, reducers, backend, chunk, W0) where {d, T, k}
    n = length(x0)
    choice = resolve_backend(backend, T)
    choice.backend === :cpu || error("GPU propagation needs the CUDA package extension")
    s0 = _init_states(reducers, T, Val(d), Val(k))
    csize = chunk !== nothing ? chunk :
            chunk_size(Sys.free_memory(), Val(d), Val(k), T; reducer_bytes = isbitstype(typeof(s0)) ? sizeof(s0) : 0)
    m = min(csize, n)
    bufs = (Vector{SVector{d, T}}(undef, m), Vector{SMatrix{d, k, T, d * k}}(undef, m),
            Vector{SVector{d, T}}(undef, m), Vector{SMatrix{d, k, T, d * k}}(undef, m),
            Vector{typeof(s0)}(undef, m), Vector{Int8}(undef, m), Vector{Int}(undef, m))
    f0 = map(reducer_final, reducers, s0)
    out = (x = Vector{SVector{d, T}}(undef, n), W = Vector{SMatrix{d, k, T, d * k}}(undef, n),
           finals = Vector{typeof(f0)}(undef, n))
    args = (; sys, x0, seed, W0, N, tau, reducers, T, dims = (Val(d), Val(k)))
    for lo in 1:csize:n
        _run_chunk!(out, bufs, lo:min(lo + csize - 1, n), args, choice)
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason,
            precision = precision_label(T), threads = Threads.nthreads(), cpu = Sys.CPU_NAME,
            seed = seed, tau = tau, k = k, chunk_size = csize)
    return PropagationResult(out.x, out.W, out.finals, prov)
end

"""
    propagate(sys, x0, N; k, tau=1, seed, reducers=(), backend=:cpu, chunk=nothing, W0=nothing) -> PropagationResult

Propagate every orbit of `x0` for `N` steps (a positive multiple of `tau`) with `k` deviation vectors, one kernel
launch per chunk of orbits. Orbit ids are `1:length(x0)`; `W0` overrides the seeded initial deviations. Results equal
a serial `map` of [`propagate_orbit`](@ref) bitwise for any chunk size and thread count.
"""
function propagate(sys, x0::AbstractVector{SVector{d, T}}, N::Integer; k::Integer, tau::Integer = 1, seed::Integer,
                   reducers::Tuple = (), backend::Symbol = :cpu, chunk::Union{Int, Nothing} = nothing,
                   W0 = nothing) where {d, T}
    _check_steps(N, tau)
    1 <= k <= d || throw(ArgumentError("k must be in 1:$d, got $k"))
    chunk === nothing || chunk >= 1 || throw(ArgumentError("chunk must be positive, got $chunk"))
    W0 === nothing || length(W0) == length(x0) || throw(DimensionMismatch("W0 needs one matrix per orbit"))
    return _propagate(sys, x0, Int(N), Int(tau), Val(Int(k)), seed, reducers, backend, chunk, W0)
end
