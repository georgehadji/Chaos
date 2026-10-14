# Two-trajectory estimator (spec lyap-two-trajectory, Algorithm 1 to 4): an own per-orbit loop that mirrors Tangent's
# chunk, kernel and fallback pattern (decision D34), without a tangent map.

using KernelAbstractions
using ..Foundation: DEVICE, NoDevice, free_memory, gpu_step_ns
using ..Tangent: device_backend, workgroup, _fallback_reason, CHUNK_HEADROOM, LAUNCH_WAVE, _budget_orbits, _check_steps, _host

# Per-orbit end record. code: 0 ok, 1 non-finite state, 2 separation zero or non-finite, 3 realized separation outside
# [d0/2, 2 d0]; `value` is the offending number and `step` its handoff.
struct ShadowEnd{s, T}
    x::SVector{s, T}
    y::SVector{s, T}
    sum::CompensatedSum{T}
    dprev::T
    code::Int8
    step::Int
    value::T
end

# d0 must widen exactly into T (same rule as `add`), be finite and positive.
function _check_d0(::Type{T}, d0) where {T}
    d0 isa T || d0 isa Integer || (d0 isa AbstractFloat && promote_type(T, typeof(d0)) === T) ||
        throw(ArgumentError("d0 of type $(typeof(d0)) would be narrowed to $T"))
    d = T(d0)
    (isfinite(d) && d > 0) || throw(ArgumentError("d0 must be finite and positive, got $d0"))
    return d
end

# Algorithm 1: shadow at ref + d0 u0, its realized separation, and a code (0, 2 or 3).
@inline function _shadow_start(sys, x::SVector, u0::SVector, d0)
    y = displaced(sys, x, d0 * u0)
    dp = safe_norm(separation(sys, x, y))
    code = !(isfinite(dp) && dp > 0) ? Int8(2) : (dp < d0 / 2 || dp > 2 * d0 ? Int8(3) : Int8(0))
    return y, dp, code
end

@inline function _first_bad_value(x::SVector, y::SVector)
    for v in x
        isfinite(v) || return v
    end
    for v in y
        isfinite(v) || return v
    end
    return zero(eltype(x))
end

# The per-orbit body (Algorithm 1 to 3), non-throwing so the same code runs in the kernel.
@inline function _shadow_body(sys, x0::SVector{s, T}, u0::SVector, N::Int, tau::Int, d0::T) where {s, T}
    x = x0
    y, dprev, code = _shadow_start(sys, x, u0, d0)
    acc = CompensatedSum{T}()
    code == 0 || return ShadowEnd(x, y, acc, dprev, code, 0, dprev)
    n = 0
    for _ in 1:(N ÷ tau)
        for _ in 1:tau
            x = forward(sys, x)
            y = forward(sys, y)
        end
        n += tau
        (all(isfinite, x) && all(isfinite, y)) || return ShadowEnd(x, y, acc, dprev, Int8(1), n, _first_bad_value(x, y))
        delta = separation(sys, x, y)
        d = safe_norm(delta)
        isfinite(d) && d > 0 || return ShadowEnd(x, y, acc, dprev, Int8(2), n, d)
        acc = add(acc, log(d / dprev))
        y = displaced(sys, x, delta * (d0 / d))
        dprev = safe_norm(separation(sys, x, y))
        (dprev >= d0 / 2 && dprev <= 2 * d0) ||
            return ShadowEnd(x, y, acc, dprev, isfinite(dprev) && dprev > 0 ? Int8(3) : Int8(2), n, dprev)
    end
    return ShadowEnd(x, y, acc, dprev, Int8(0), N, dprev)
end

_final(e::ShadowEnd, N::Int) = (n = N, sum = e.sum.s, comp = e.sum.c, X2 = total(e.sum) / N, d_prev = e.dprev)

function _shadow_error(e::ShadowEnd, d0, orbit, backend::Symbol)
    v = Float64(e.value)
    e.code == 1 && return NumericalHealthError("state", v, Inf, e.step, orbit, backend)
    thr = e.code == 3 ? Float64(d0) : (isfinite(v) ? 0.0 : Inf)
    return NumericalHealthError("separation", v, thr, e.step, orbit, backend)
end

"""
    propagate_shadow_orbit(sys, x0, u0, N, tau, d0; orbit=1, backend=:cpu) -> (x, xs, final)

Per-orbit body of [`propagate_shadow`](@ref): reference orbit `x0` and a shadow started at `x0 + d0 u0`, `N` steps,
the shadow brought back to distance `d0` every `tau` steps. Returns the final reference state `x`, the shadow state `xs`
and `final = (n, sum, comp, X2, d_prev)`. Throws `NumericalHealthError` (`"state"` or `"separation"`).
"""
function propagate_shadow_orbit(sys, x0::SVector{s, T}, u0::SVector, N::Integer, tau::Integer, d0;
                                orbit::Integer = 1, backend::Symbol = :cpu) where {s, T}
    _check_steps(N, tau)
    dT = _check_d0(T, d0)
    e = _shadow_body(sys, x0, u0, Int(N), Int(tau), dT)
    e.code == 0 || throw(_shadow_error(e, dT, orbit, backend))
    return (e.x, e.y, _final(e, Int(N)))
end

@kernel function _shadow_kernel!(so, x0, u0, sys, N, tau, d0)
    i = @index(Global)
    so[i] = _shadow_body(sys, x0[i], u0[i], N, tau, d0)
end

# One chunk on `backend`; stores nothing in `ends` when an orbit failed (lowest orbit id reported first).
function _shadow_chunk!(ends, bufs, r, args, who::Symbol, backend)
    xin, uin, so = bufs
    len = length(r)
    copyto!(xin, 1, args.x0, first(r), len)
    copyto!(uin, 1, args.u, first(r), len)
    _shadow_kernel!(backend, workgroup(backend))(so, xin, uin, args.sys, args.N, args.tau, args.d0; ndrange = len)
    KernelAbstractions.synchronize(backend)
    h = _host(so)
    i = findfirst(e -> e.code != 0, view(h, 1:len))
    i === nothing || throw(_shadow_error(h[i], args.d0, first(r) + i - 1, who))
    copyto!(ends, first(r), h, 1, len)
end

# A chunk on the device: 0 when it ran, else the fallback reason index (as Tangent._device_chunk!).
function _shadow_device_chunk!(ends, bufs, r, args, backend)
    try
        _shadow_chunk!(ends, bufs, r, args, :gpu, backend)
        return 0
    catch e
        why = _fallback_reason(e)
        why == 0 && rethrow()
        return why
    end
end

# Host check before any chunk (Algorithm 1): every realized d_0^+ must lie in [d0/2, 2 d0].
function _check_starts(sys, x0, u, d0, ::Type{T}) where {T}
    for i in eachindex(x0)
        _, dp, code = _shadow_start(sys, x0[i], u[i], d0)
        code == 0 || throw(_shadow_error(ShadowEnd(x0[i], x0[i], CompensatedSum{T}(), dp, code, 0, dp), d0, i, :cpu))
    end
end

function _shadow_directions(seed, U0, n, ::Val{d}, ::Type{T}) where {d, T}
    U0 === nothing && return [SVector{d, T}(initial_deviations(seed, j, Val(d), Val(1), T)[:, 1]) for j in 1:n]
    length(U0) == n || throw(DimensionMismatch("U0 needs one vector per orbit"))
    return [(u = SVector{d, T}(U0[j]); u / safe_norm(u)) for j in 1:n]
end

"""
    propagate_shadow(sys, x0, N; d0, tau=1, seed, backend=:cpu, chunk=nothing, U0=nothing) -> PropagationResult

Two-trajectory estimate X2 of the largest exponent for every orbit of `x0` (spec lyap-two-trajectory): a shadow orbit starts
at distance `d0` from the reference orbit along the seeded direction (`U0` overrides it) and is brought back to `d0` every
`tau` steps. `x[i]` is the final reference state, `W[i]` the final shadow state, `finals[i] = (final,)` with
`final = (n, sum, comp, X2, d_prev)`. `d0` must widen exactly into the element type. On the CPU backend the result equals a
serial `map` of [`propagate_shadow_orbit`](@ref) bitwise for any chunk size and thread count.
"""
function propagate_shadow(sys, x0::AbstractVector{SVector{s, T}}, N::Integer; d0, tau::Integer = 1, seed::Integer,
                          backend::Symbol = :cpu, chunk::Union{Int, Nothing} = nothing, U0 = nothing) where {s, T}
    _check_steps(N, tau)
    d = Int(dimension(sys))
    s >= d || throw(ArgumentError("state length $s is below the system dimension $d"))
    chunk === nothing || chunk >= 1 || throw(ArgumentError("chunk must be positive, got $chunk"))
    dT = _check_d0(T, d0)
    u = _shadow_directions(seed, U0, length(x0), Val(d), T)
    _check_starts(sys, x0, u, dT, T)
    return _shadow_run(sys, x0, u, Val(d), Int(N), Int(tau), dT, seed, backend, chunk)
end

function _shadow_run(sys, x0::AbstractVector{SVector{s, T}}, u, ::Val{d}, N, tau, d0, seed, backend, chunk) where {s, T, d}
    n = length(x0)
    csize_for(free) = chunk !== nothing ? chunk :
                      chunk_size(free, Val(d), Val(0), T; f = CHUNK_HEADROOM, reducer_bytes = 3 * sizeof(T), s = 2s)
    alloc(b, cs) = (m = min(cs, n); (KernelAbstractions.allocate(b, SVector{s, T}, m),
                                     KernelAbstractions.allocate(b, SVector{d, T}, m),
                                     KernelAbstractions.allocate(b, ShadowEnd{s, T}, m)))
    choice = resolve_backend(backend, T)
    dev, csize, gpu_csize, gpu_chunks, bufs = CPU(), 0, 0, 0, nothing
    if choice.backend === :gpu
        cap = _budget_orbits(something(gpu_step_ns(DEVICE[], T), 100.0), N)
        if cap < LAUNCH_WAVE   # N alone exceeds the launch budget (D33)
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[8])
        else
            try
                csize = min(csize_for(free_memory(DEVICE[])), cap)
                dev = device_backend(DEVICE[])
                bufs = alloc(dev, csize)
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
        bufs = alloc(dev, csize)
    end
    ends = Vector{ShadowEnd{s, T}}(undef, n)
    args = (; sys, x0, u, N, tau, d0)
    lo = 1
    while lo <= n
        r = lo:min(lo + csize - 1, n)
        if choice.backend === :gpu
            why = _shadow_device_chunk!(ends, bufs, r, args, dev)
            if why == 0
                gpu_chunks += 1
                lo = last(r) + 1
                continue
            end
            why == 8 && (DEVICE[] = NoDevice())
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[why])
            dev, csize = CPU(), csize_for(Sys.free_memory())
            bufs = alloc(dev, csize)
            continue
        end
        _shadow_chunk!(ends, bufs, r, args, :cpu, dev)
        lo = last(r) + 1
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason,
            precision = precision_label(T), threads = Threads.nthreads(), cpu = Sys.CPU_NAME,
            seed = seed, tau = tau, k = 0, chunk_size = csize, headroom = CHUNK_HEADROOM,
            gpu_chunks = gpu_chunks, gpu_chunk_size = gpu_csize,
            estimator = "two-trajectory", d0 = Float64(d0), state_mode = state_mode(sys))
    return PropagationResult([e.x for e in ends], [e.y for e in ends], [(_final(e, N),) for e in ends], prov)
end
