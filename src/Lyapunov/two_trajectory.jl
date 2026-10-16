# Two-trajectory estimator (spec lyap-two-trajectory, Algorithm 1 to 4): an own per-orbit loop that mirrors Tangent's
# chunk, kernel and fallback pattern (decision D34), without a tangent map.

using KernelAbstractions
using ..Foundation: DEVICE, NoDevice, free_memory
using ..Tangent: device_backend, workgroup, _fallback_reason, CHUNK_HEADROOM, _check_steps, _host, _before_launch, _forced_launch,
                 _device_plan, _launches

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

# The per-orbit body (Algorithm 1 to 3), non-throwing so the same code runs in the kernel. One launch (D44): `steps` steps
# (a multiple of tau) from global step n0; at n0 = 0 the shadow starts at x0 + d0 u0, later launches continue from the
# carried end record `e` (reference, shadow, sum, previous separation). A run of N steps is the call with n0 = 0, steps = N.
@inline function _shadow_body(sys, x0::SVector{s, T}, u0::SVector, e::ShadowEnd{s, T}, n0::Int, steps::Int, N::Int, tau::Int,
                              d0::T) where {s, T}
    if n0 == 0
        x = x0
        y, dprev, code = _shadow_start(sys, x, u0, d0)
        acc = CompensatedSum{T}()
        code == 0 || return ShadowEnd(x, y, acc, dprev, code, 0, dprev)
    else
        x, y, acc, dprev = e.x, e.y, e.sum, e.dprev
    end
    n = n0
    for _ in 1:(steps ÷ tau)
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
    return ShadowEnd(x, y, acc, dprev, Int8(0), n, dprev)
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
    e = _shadow_body(sys, x0, u0, _blank_end(x0, dT), 0, Int(N), Int(N), Int(tau), dT)
    e.code == 0 || throw(_shadow_error(e, dT, orbit, backend))
    return (e.x, e.y, _final(e, Int(N)))
end

# One work item per orbit and one launch (D44), in place on the carried end records `so`; an orbit that failed in an
# earlier launch keeps its record and is not advanced.
@kernel function _shadow_kernel!(so, x0, u0, sys, n0, steps, N, tau, d0)
    i = @index(Global)
    if n0 == 0
        so[i] = _shadow_body(sys, x0[i], u0[i], _blank_end(x0[i], d0), n0, steps, N, tau, d0)
    elseif so[i].code == 0
        so[i] = _shadow_body(sys, x0[i], u0[i], so[i], n0, steps, N, tau, d0)
    end
end

# Placeholder end record for the first launch (never read: the body starts the shadow from x0 when n0 = 0).
@inline _blank_end(x::SVector{s, T}, d0::T) where {s, T} = ShadowEnd(x, x, CompensatedSum{T}(), d0, Int8(0), 0, d0)

# One chunk on `backend` in launches of `args.L` steps; stores nothing in `ends` when an orbit failed (lowest orbit id
# reported first). `ci` is the device chunk index for the hook.
function _shadow_chunk!(ends, bufs, r, args, who::Symbol, backend, ci::Int = 0)
    xin, uin, so = bufs
    len = length(r)
    copyto!(xin, 1, args.x0, first(r), len)
    copyto!(uin, 1, args.u, first(r), len)
    kernel = _shadow_kernel!(backend, workgroup(backend))
    n0, launch = 0, 0
    while n0 < args.N
        steps = min(args.L, args.N - n0)
        launch += 1
        who === :gpu && _before_launch(DEVICE[], ci, launch, (ends = so,))
        kernel(so, xin, uin, args.sys, n0, steps, args.N, args.tau, args.d0; ndrange = len)
        KernelAbstractions.synchronize(backend)
        n0 += steps
    end
    h = _host(so)
    i = findfirst(e -> e.code != 0, view(h, 1:len))
    i === nothing || throw(_shadow_error(h[i], args.d0, first(r) + i - 1, who))
    copyto!(ends, first(r), h, 1, len)
end

# A chunk on the device: 0 when it ran, else the fallback reason index (as Tangent._device_chunk!).
function _shadow_device_chunk!(ends, bufs, r, args, backend, ci::Int)
    try
        _shadow_chunk!(ends, bufs, r, args, :gpu, backend, ci)
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
serial `map` of [`propagate_shadow_orbit`](@ref) bitwise for any chunk size and thread count. On the device a long run is split
into launches that carry the state, with the launch budget, fallback reasons and provenance fields of [`propagate`](@ref) (D44).
"""
function propagate_shadow(sys, x0::AbstractVector{SVector{s, T}}, N::Integer; d0, tau::Integer = 1, seed::Integer,
                          backend::Symbol = :cpu, chunk::Union{Int, Nothing} = nothing, U0 = nothing) where {s, T}
    _check_steps(N, tau)
    d = Int(dimension(sys))
    s >= d || throw(ArgumentError("state length $s is below the system dimension $d"))
    chunk === nothing || chunk >= 1 || throw(ArgumentError("chunk must be positive, got $chunk"))
    dT = _check_d0(T, d0)
    u = _shadow_directions(seed, U0, length(x0), Val(d), T)
    forced = _forced_launch(tau)
    _check_starts(sys, x0, u, dT, T)
    return _shadow_run(sys, x0, u, Val(d), Int(N), Int(tau), dT, seed, backend, chunk, forced)
end

function _shadow_run(sys, x0::AbstractVector{SVector{s, T}}, u, ::Val{d}, N, tau, d0, seed, backend, chunk, forced) where {s, T, d}
    n = length(x0)
    csize_for(free) = chunk !== nothing ? chunk :
                      chunk_size(free, Val(d), Val(0), T; f = CHUNK_HEADROOM, reducer_bytes = 3 * sizeof(T), s = 2s)
    alloc(b, cs) = (m = min(cs, n); (KernelAbstractions.allocate(b, SVector{s, T}, m),
                                     KernelAbstractions.allocate(b, SVector{d, T}, m),
                                     KernelAbstractions.allocate(b, ShadowEnd{s, T}, m)))
    choice = resolve_backend(backend, T)
    dev, csize, gpu_csize, gpu_chunks, bufs = CPU(), 0, 0, 0, nothing
    Lcpu = forced === nothing ? N : min(forced, N)   # the CPU is split only by the test hook (D44 item 5)
    L, launches = forced === nothing ? (0, 0) : (Lcpu, _launches(N, Lcpu))
    Ldev = 0
    if choice.backend === :gpu
        Ldev, cap = _device_plan(forced, T, N, tau)   # the tangent stream's budget and default ns value (D44 item 15)
        if Ldev == 0   # even L = tau exceeds the budget for one wave
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[8])
        else
            try
                csize = min(csize_for(free_memory(DEVICE[])), cap)
                dev = device_backend(DEVICE[])
                bufs = alloc(dev, csize)
                gpu_csize = csize
                L, launches = Ldev, _launches(N, Ldev)
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
        bufs = alloc(dev, csize)
    end
    ends = Vector{ShadowEnd{s, T}}(undef, n)
    args = (; sys, x0, u, N, tau, d0, L = choice.backend === :gpu ? Ldev : Lcpu)
    lo = 1
    while lo <= n
        r = lo:min(lo + csize - 1, n)
        if choice.backend === :gpu
            why = _shadow_device_chunk!(ends, bufs, r, args, dev, gpu_chunks + 1)
            if why == 0
                gpu_chunks += 1
                lo = last(r) + 1
                continue
            end
            why == 9 && (DEVICE[] = NoDevice())
            choice = BackendChoice(:cpu, "cpu", FALLBACK_REASONS[why])
            dev, csize = CPU(), csize_for(Sys.free_memory())
            bufs = alloc(dev, csize)
            args = merge(args, (L = Lcpu,))
            continue
        end
        _shadow_chunk!(ends, bufs, r, args, :cpu, dev)
        lo = last(r) + 1
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason,
            precision = precision_label(T), threads = Threads.nthreads(), cpu = Sys.CPU_NAME,
            seed = seed, tau = tau, k = 0, chunk_size = csize, headroom = CHUNK_HEADROOM,
            gpu_chunks = gpu_chunks, gpu_chunk_size = gpu_csize, launch_steps = L, launches_per_chunk = launches,
            estimator = "two-trajectory", d0 = Float64(d0), state_mode = state_mode(sys))
    return PropagationResult([e.x for e in ends], [e.y for e in ends], [(_final(e, N),) for e in ends], prov)
end
