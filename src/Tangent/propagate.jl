"""
    bytes_per_orbit(Val(d), Val(k), T, s = d) -> Int

Resident bytes per orbit: `sizeof(T) (s + d k + 2 k) + 16` for state length `s` (default `d`; spec Compute, amendment D28).
"""
bytes_per_orbit(::Val{d}, ::Val{k}, ::Type{T}, s::Integer = d) where {d, k, T} = sizeof(T) * (s + d * k + 2k) + 16

const CHUNK_HEADROOM = 0.5   # fraction of free memory a chunk may use

"""
    chunk_size(free_bytes, Val(d), Val(k), T; f=0.5, reducer_bytes=0, s=d) -> Int

Orbits per chunk: `floor(f free_bytes / (bytes_per_orbit + reducer_bytes))`, at least 1; `s` is the state length.
"""
function chunk_size(free_bytes::Integer, v::Val{d}, kv::Val{k}, ::Type{T}; f = CHUNK_HEADROOM, reducer_bytes::Integer = 0,
                    s::Integer = d) where {d, k, T}
    return max(1, floor(Int, f * free_bytes / (bytes_per_orbit(v, kv, T, s) + reducer_bytes)))
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

"""
    device_backend(device) -> KernelAbstractions.Backend

Backend a chunk runs on when the run resolved to `device` (`Foundation.DEVICE[]`); the CUDA extension adds `CUDABackend`.
"""
device_backend(::NoDevice) = CPU()

# Orbits per work group; the sequential work per orbit is long, so the CPU groups stay small.
workgroup(::CPU) = 16

# Initial deviation matrices of one chunk, on the host for every backend (spec, initial vectors).
function _fill_W!(Win, r, seed, W0, ::Val{d}, ::Val{k}, ::Type{T}) where {d, k, T}
    for (i, j) in enumerate(r)
        W0 === nothing || promote_type(T, eltype(W0[j])) === T ||
            throw(ArgumentError("W0[$j] has element type $(eltype(W0[j])), wider than $T"))
        Win[i] = W0 === nothing ? initial_deviations(seed, j, Val(d), Val(k), T) : SMatrix{d, k, T}(W0[j])
    end
end

_host(a::Array) = a
_host(a) = Array(a)

# Seven chunk buffers on `backend`: x in, W in, x out, W out, reducer states, status codes, failing step.
function _alloc_bufs(backend, ::Type{SVector{s, T}}, ::Val{d}, ::Val{k}, ::Type{S}, m) where {s, d, k, T, S}
    M = SMatrix{d, k, T, d * k}
    A = KernelAbstractions.allocate
    return (A(backend, SVector{s, T}, m), A(backend, M, m), A(backend, SVector{s, T}, m), A(backend, M, m),
            A(backend, S, m), A(backend, Int8, m), A(backend, Int, m))
end

# One chunk on `backend`; `who` (:cpu or :gpu) goes into the errors. Stores nothing in `out` when a status is not OK.
function _run_chunk!(out, bufs, Wh, r, args, who::Symbol, backend)
    xin, Win, xo, Wo, so, code, hstep = bufs
    len = length(r)
    copyto!(xin, 1, args.x0, first(r), len)
    _fill_W!(Wh, r, args.seed, args.W0, args.dims..., args.T)
    copyto!(Win, 1, Wh, 1, len)
    _orbit_kernel!(backend, workgroup(backend))(xo, Wo, so, code, hstep, args.sys, xin, Win, args.N, args.tau, args.reducers; ndrange = len)
    KernelAbstractions.synchronize(backend)
    code, so, xo, Wo, hstep = _host(code), _host(so), _host(xo), _host(Wo), _host(hstep)
    i = findfirst(!=(STATUS_OK), view(code, 1:len))   # lowest orbit id first; a failed chunk stores nothing
    i === nothing || code[i] != STATUS_REDUCER || throw(_reducer_error(args.reducers, so[i], first(r) + i - 1, who))
    i === nothing || throw(_health_error(code[i], xo[i], Wo[i], hstep[i], first(r) + i - 1, who))
    copyto!(out.x, first(r), xo, 1, len)
    copyto!(out.W, first(r), Wo, 1, len)
    for i in 1:len
        out.finals[first(r) + i - 1] = map(reducer_final, args.reducers, so[i])
    end
end

# FALLBACK_REASONS index that sends a failed device step to the CPU (D29, D33), or 0 for an exception that is not ours:
# 5 a health error (the CPU rerun decides whether it is real), 8 a lost device or context, 4 out of memory, 3 compile failure.
function _fallback_reason(e)
    e isa NumericalHealthError && return 5
    is_device_lost(e) && return 8
    is_oom(e) && return 4
    is_compile_error(e) && return 3
    return 0
end

# A chunk on the device. Returns 0 when it ran, else the reason index; anything unclassified is rethrown.
function _device_chunk!(out, bufs, Wh, r, args, backend)
    try
        _run_chunk!(out, bufs, Wh, r, args, :gpu, backend)
        return 0
    catch e
        why = _fallback_reason(e)
        why == 0 && rethrow()
        return why
    end
end

# Launch budget (D33): one launch must stay under the display watchdog. A launch takes about
# gpu_ns_per_step x N x orbits, and never less than for one wave of LAUNCH_WAVE resident orbits.
const LAUNCH_BUDGET_S = 1.0   # seconds per launch
const LAUNCH_WAVE = 1024      # orbits that run concurrently on the card measured in D29 (the launch time is flat below this)
const DEFAULT_GPU_NS = 100.0  # ns per orbit-step when the machine has no record: about 5 times the value measured on a GTX 1050

# Largest chunk (orbits per launch) whose estimated launch time is within the budget; may be below LAUNCH_WAVE.
_budget_orbits(ns, N) = floor(Int, min(LAUNCH_BUDGET_S / (ns * 1e-9 * N), typemax(Int) / 4))

function _propagate(sys, x0::AbstractVector{SVector{s, T}}, N, tau, ::Val{d}, ::Val{k}, seed, reducers, backend, chunk, W0) where {s, d, T, k}
    n = length(x0)
    s0 = _init_states(reducers, T, Val(d), Val(k))
    rbytes = isbitstype(typeof(s0)) ? sizeof(s0) : 0
    csize_for(free) = chunk !== nothing ? chunk :
                      chunk_size(free, Val(d), Val(k), T; f = CHUNK_HEADROOM, reducer_bytes = rbytes, s = s)
    # chunk buffers on `b` for chunks of `cs` orbits, and the host staging array of the initial deviations
    alloc(b, cs) = (_alloc_bufs(b, SVector{s, T}, Val(d), Val(k), typeof(s0), min(cs, n)),
                    Vector{SMatrix{d, k, T, d * k}}(undef, min(cs, n)))
    choice = resolve_backend(backend, T; cpu_only_reducer = !isbitstype(typeof(s0)))
    dev, csize, gpu_csize, gpu_chunks, bufs, Wh = CPU(), 0, 0, 0, nothing, nothing
    if choice.backend === :gpu
        cap = _budget_orbits(something(gpu_step_ns(DEVICE[], T), DEFAULT_GPU_NS), N)
        if cap < LAUNCH_WAVE   # N alone exceeds the budget
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[8])
        else
            try
                csize = min(csize_for(free_memory(DEVICE[])), cap)
                dev = device_backend(DEVICE[])
                bufs, Wh = alloc(dev, csize)
                gpu_csize = csize
            catch e
                why = _fallback_reason(e)
                why == 0 && rethrow()
                why == 8 && (DEVICE[] = NoDevice())
                choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[why])
            end
        end
    end
    if choice.backend === :cpu
        dev, csize = CPU(), csize_for(Sys.free_memory())
        bufs, Wh = alloc(dev, csize)
    end
    f0 = map(reducer_final, reducers, s0)
    out = (x = Vector{SVector{s, T}}(undef, n), W = Vector{SMatrix{d, k, T, d * k}}(undef, n),
           finals = Vector{typeof(f0)}(undef, n))
    args = (; sys, x0, seed, W0, N, tau, reducers, T, dims = (Val(d), Val(k)))
    lo = 1
    while lo <= n
        r = lo:min(lo + csize - 1, n)
        if choice.backend === :gpu
            why = _device_chunk!(out, bufs, Wh, r, args, dev)
            if why == 0
                gpu_chunks += 1
                lo = last(r) + 1
                continue
            end
            # the rest of the run stays on the CPU, with a chunk size from host memory; the failed chunk stored nothing
            why == 8 && (DEVICE[] = NoDevice())
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[why])
            dev, csize = CPU(), csize_for(Sys.free_memory())
            bufs, Wh = alloc(dev, csize)
            continue
        end
        _run_chunk!(out, bufs, Wh, r, args, :cpu, dev)
        lo = last(r) + 1
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason,
            precision = precision_label(T), threads = Threads.nthreads(), cpu = Sys.CPU_NAME,
            seed = seed, tau = tau, k = k, chunk_size = csize, headroom = CHUNK_HEADROOM,
            gpu_chunks = gpu_chunks, gpu_chunk_size = gpu_csize)
    return PropagationResult(out.x, out.W, out.finals, prov)
end

"""
    propagate(sys, x0, N; k, tau=1, seed, reducers=(), backend=:cpu, chunk=nothing, W0=nothing) -> PropagationResult

Propagate every orbit of `x0` for `N` steps (a positive multiple of `tau`) with `k` deviation vectors, one kernel
launch per chunk of orbits. The state may be longer than `dimension(sys)` (winding counts, D28): deviation matrices
are `dimension(sys) x k`, and every state component is health-checked. Orbit ids are `1:length(x0)`; `W0` overrides the seeded initial deviations. On the CPU backend results equal
a serial `map` of [`propagate_orbit`](@ref) bitwise for any chunk size and thread count. GPU results are not bitwise equal to the CPU: they agree within the tangent-core AC12
tolerances (estimator level; same-orbit |x difference| <= 4 n eps G^n). `provenance` records the backend that finished the run, the fallback reason, `chunk_size` of that backend,
and how many chunks ran on the device first (`gpu_chunks`, `gpu_chunk_size`). A device launch is kept under a 1 s estimated budget (decision D33): the chunk shrinks, and when
N alone exceeds the budget the run uses the CPU with the reason "launch exceeds the display watchdog".
"""
function propagate(sys, x0::AbstractVector{SVector{s, T}}, N::Integer; k::Integer, tau::Integer = 1, seed::Integer,
                   reducers::Tuple = (), backend::Symbol = :cpu, chunk::Union{Int, Nothing} = nothing,
                   W0 = nothing) where {s, T}
    _check_steps(N, tau)
    d = dimension(sys)
    s >= d || throw(ArgumentError("state length $s is below the system dimension $d"))
    1 <= k <= d || throw(ArgumentError("k must be in 1:$d, got $k"))
    chunk === nothing || chunk >= 1 || throw(ArgumentError("chunk must be positive, got $chunk"))
    W0 === nothing || length(W0) == length(x0) || throw(DimensionMismatch("W0 needs one matrix per orbit"))
    return _propagate(sys, x0, Int(N), Int(tau), Val(Int(d)), Val(Int(k)), seed, reducers, backend, chunk, W0)
end
