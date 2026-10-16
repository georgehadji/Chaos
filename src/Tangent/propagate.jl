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

# One work item per orbit and one launch of `steps` steps from global step n0 (D44), in place on the carried buffers;
# the body is the same `_orbit_body` as `propagate_orbit`. The first launch draws fresh reducer states; an orbit that
# failed in an earlier launch keeps its status and failing step and is not advanced.
@kernel function _orbit_kernel!(x, W, st, code, hstep, sys, n0, steps, N, tau, reducers)
    i = @index(Global)
    if n0 == 0 || code[i] == STATUS_OK
        xi, Wi = x[i], W[i]
        s0 = n0 == 0 ? _init_states(reducers, eltype(xi), Val(size(Wi, 1)), Val(size(Wi, 2))) : st[i]
        xn, Wn, sn, c, h = _orbit_body(sys, xi, Wi, s0, n0, steps, N, tau, reducers)
        x[i] = xn
        W[i] = Wn
        st[i] = sn
        code[i] = c
        hstep[i] = h
    end
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

# Five chunk buffers on `backend`, carried between launches (D44): points, deviation matrices, reducer states, status
# codes, failing step.
function _alloc_bufs(backend, ::Type{SVector{s, T}}, ::Val{d}, ::Val{k}, ::Type{S}, m) where {s, d, k, T, S}
    A = KernelAbstractions.allocate
    return (A(backend, SVector{s, T}, m), A(backend, SMatrix{d, k, T, d * k}, m), A(backend, S, m), A(backend, Int8, m),
            A(backend, Int, m))
end

"""
    _before_launch(device, chunk_index, launch, carried)

Test hook (D44), called on the host before every device launch with the device chunk index, the launch index within the
chunk and the live carried buffers (`(x, W, states)` for `propagate`, `(ends,)` for `propagate_shadow`). No-op by default.
"""
_before_launch(::AbstractDevice, ::Int, ::Int, carried) = nothing

# One chunk on `backend` in launches of `args.L` steps (the last may be shorter); `who` (:cpu or :gpu) goes into the
# errors, `ci` is the device chunk index for the hook. Stores nothing in `out` when a status is not OK.
function _run_chunk!(out, bufs, Wh, r, args, who::Symbol, backend, ci::Int = 0)
    x, W, so, code, hstep = bufs
    len = length(r)
    copyto!(x, 1, args.x0, first(r), len)
    _fill_W!(Wh, r, args.seed, args.W0, args.dims..., args.T)
    copyto!(W, 1, Wh, 1, len)
    kernel = _orbit_kernel!(backend, workgroup(backend))
    n0, launch = 0, 0
    while n0 < args.N
        steps = min(args.L, args.N - n0)
        launch += 1
        who === :gpu && _before_launch(DEVICE[], ci, launch, (x = x, W = W, states = so))
        kernel(x, W, so, code, hstep, args.sys, n0, steps, args.N, args.tau, args.reducers; ndrange = len)
        KernelAbstractions.synchronize(backend)
        n0 += steps
    end
    code, so, xo, Wo, hstep = _host(code), _host(so), _host(x), _host(W), _host(hstep)
    i = findfirst(!=(STATUS_OK), view(code, 1:len))   # lowest orbit id first; a failed chunk stores nothing
    i === nothing || code[i] != STATUS_REDUCER || throw(_reducer_error(args.reducers, so[i], first(r) + i - 1, who))
    i === nothing || throw(_health_error(code[i], xo[i], Wo[i], hstep[i], first(r) + i - 1, who))
    copyto!(out.x, first(r), xo, 1, len)
    copyto!(out.W, first(r), Wo, 1, len)
    for i in 1:len
        out.finals[first(r) + i - 1] = map(reducer_final, args.reducers, so[i])
    end
end

# FALLBACK_REASONS index that sends a failed device step to the CPU (D29, D33, D44), or 0 for an exception that is not ours:
# 5 a health error (the CPU rerun decides whether it is real), 9 a lost device or context, 4 out of memory, 3 compile failure.
function _fallback_reason(e)
    e isa NumericalHealthError && return 5
    is_device_lost(e) && return 9
    is_oom(e) && return 4
    is_compile_error(e) && return 3
    return 0
end

# A chunk on the device, every launch inside one guard: returns 0 when it ran, else the reason index; anything
# unclassified is rethrown. A failed chunk is rerun on the CPU from its initial state, never from the carried state (D44).
function _device_chunk!(out, bufs, Wh, r, args, backend, ci::Int)
    try
        _run_chunk!(out, bufs, Wh, r, args, :gpu, backend, ci)
        return 0
    catch e
        why = _fallback_reason(e)
        why == 0 && rethrow()
        return why
    end
end

# Launch budget (D33, D44): one launch must stay under the display watchdog. A launch of L steps takes about
# ns x L x orbits, and never less than for one wave of LAUNCH_WAVE resident orbits; ns is the record's one-wave value.
const LAUNCH_BUDGET_S = 1.0   # seconds per launch
const LAUNCH_WAVE = 1024      # orbits that run concurrently on the card measured in D29 (the launch time is flat below this)
const DEFAULT_GPU_NS = 100.0  # ns per orbit-step when the machine has no record: about 5 times the value measured on a GTX 1050

# Largest chunk (orbits per launch) whose estimated launch of L steps is within the budget; may be below LAUNCH_WAVE.
_budget_orbits(ns, L) = floor(Int, min(LAUNCH_BUDGET_S / (ns * 1e-9 * L), typemax(Int) / 4))

"""
    _budget_steps(ns, tau) -> Int

Steps per launch (D44 item 6): the largest multiple of `tau` whose estimate for one wave of `LAUNCH_WAVE` orbits,
`ns 1e-9 L LAUNCH_WAVE`, is within the 1 s budget; 0 when even `L = tau` exceeds it.
"""
function _budget_steps(ns, tau::Integer)
    fits(m) = ns * 1e-9 * (m * tau) * LAUNCH_WAVE <= LAUNCH_BUDGET_S
    m = floor(Int, min(LAUNCH_BUDGET_S / (ns * 1e-9 * LAUNCH_WAVE * tau), 2.0^52))
    while m > 0 && !fits(m)
        m -= 1
    end
    while fits(m + 1)
        m += 1
    end
    return m * tau
end

"Forced steps per launch for tests (D44 item 5), on the CPU backend and the device path; `nothing` in normal use."
const _LAUNCH_STEPS = Ref{Union{Nothing, Int}}(nothing)

# The forced launch length, checked against tau: nothing, or a positive multiple of tau.
function _forced_launch(tau::Integer)
    L = _LAUNCH_STEPS[]
    L === nothing || (L >= 1 && L % tau == 0) ||
        throw(ArgumentError("forced steps per launch must be a positive multiple of tau, got L = $L, tau = $tau"))
    return L
end

# Device launch plan (D44): steps per launch L (0 when even L = tau exceeds the budget for one wave) and the orbit cap
# at L, never below one wave (item 6). A forced L (tests, the one-wave record run) keeps the wave floor, so the record
# run of 1024 orbits stays one chunk whatever value is on disk.
function _device_plan(forced, ::Type{T}, N::Int, tau::Int) where {T}
    ns = something(gpu_step_ns(DEVICE[], T), DEFAULT_GPU_NS)
    L = forced === nothing ? min(_budget_steps(ns, tau), N) : min(forced, N)
    return L == 0 ? (0, 0) : (L, max(LAUNCH_WAVE, _budget_orbits(ns, L)))
end

# Launches per chunk of N steps in launches of L (0 when L is 0: no launch).
_launches(N, L) = L == 0 ? 0 : cld(N, L)

function _propagate(sys, x0::AbstractVector{SVector{s, T}}, N, tau, ::Val{d}, ::Val{k}, seed, reducers, backend, chunk, W0, forced) where {s, d, T, k}
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
    Lcpu = forced === nothing ? N : min(forced, N)   # the CPU is split only by the test hook (D44 item 5)
    L, launches = forced === nothing ? (0, 0) : (Lcpu, _launches(N, Lcpu))
    Ldev = 0
    if choice.backend === :gpu
        Ldev, cap = _device_plan(forced, T, N, tau)
        if Ldev == 0   # even L = tau exceeds the budget for one wave
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[8])
        else
            try
                csize = min(csize_for(free_memory(DEVICE[])), cap)
                dev = device_backend(DEVICE[])
                bufs, Wh = alloc(dev, csize)
                gpu_csize = csize
                L, launches = Ldev, _launches(N, Ldev)   # kept after a mid-run fallback (D44 item 13)
            catch e
                why = _fallback_reason(e)
                why == 0 && rethrow()
                why == 9 && (DEVICE[] = NoDevice())
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
    args = (; sys, x0, seed, W0, N, tau, reducers, T, dims = (Val(d), Val(k)), L = choice.backend === :gpu ? Ldev : Lcpu)
    lo = 1
    while lo <= n
        r = lo:min(lo + csize - 1, n)
        if choice.backend === :gpu
            why = _device_chunk!(out, bufs, Wh, r, args, dev, gpu_chunks + 1)
            if why == 0
                gpu_chunks += 1
                lo = last(r) + 1
                continue
            end
            # the rest of the run stays on the CPU, with a chunk size from host memory; the failed chunk stored nothing
            why == 9 && (DEVICE[] = NoDevice())
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[why])
            dev, csize = CPU(), csize_for(Sys.free_memory())
            bufs, Wh = alloc(dev, csize)
            args = merge(args, (L = Lcpu,))
            continue
        end
        _run_chunk!(out, bufs, Wh, r, args, :cpu, dev)
        lo = last(r) + 1
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason,
            precision = precision_label(T), threads = Threads.nthreads(), cpu = Sys.CPU_NAME,
            seed = seed, tau = tau, k = k, chunk_size = csize, headroom = CHUNK_HEADROOM,
            gpu_chunks = gpu_chunks, gpu_chunk_size = gpu_csize, launch_steps = L, launches_per_chunk = launches)
    return PropagationResult(out.x, out.W, out.finals, merge(prov, map(reducer_provenance, reducers)...))
end

"""
    propagate(sys, x0, N; k, tau=1, seed, reducers=(), backend=:cpu, chunk=nothing, W0=nothing) -> PropagationResult

Propagate every orbit of `x0` for `N` steps (a positive multiple of `tau`) with `k` deviation vectors, one kernel
launch per chunk of orbits. The state may be longer than `dimension(sys)` (winding counts, D28): deviation matrices
are `dimension(sys) x k`, and every state component is health-checked. Orbit ids are `1:length(x0)`; `W0` overrides the seeded initial deviations. On the CPU backend results equal
a serial `map` of [`propagate_orbit`](@ref) bitwise for any chunk size and thread count. GPU results are not bitwise equal to the CPU: they agree within the tangent-core AC12
tolerances (estimator level; same-orbit |x difference| <= 4 n eps G^n). `provenance` records the backend that finished the run, the fallback reason, `chunk_size` of that backend,
and how many chunks ran on the device first (`gpu_chunks`, `gpu_chunk_size`). A device launch is kept under a 1 s estimated budget (decisions D33, D44): a long run is
split into launches of `launch_steps` steps (a multiple of tau) that carry the state, `launches_per_chunk` per chunk, and the chunk is capped at that length; only when even
one handoff interval exceeds the budget for one wave does the run use the CPU, with the reason "launch exceeds the display watchdog". A lost device gives "device lost".
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
    forced = _forced_launch(tau)
    return _propagate(sys, x0, Int(N), Int(tau), Val(Int(d)), Val(Int(k)), seed, reducers, backend, chunk, W0, forced)
end
