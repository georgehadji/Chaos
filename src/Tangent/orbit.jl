const STATUS_OK, STATUS_STATE, STATUS_W, STATUS_REDUCER = Int8(0), Int8(1), Int8(2), Int8(3)

_status(x, W) = all(isfinite, x) ? (all(isfinite, W) ? STATUS_OK : STATUS_W) : STATUS_STATE

_init_states(rs::Tuple, ::Type{T}, v::Val, k::Val) where {T} = map(r -> reducer_init(r, T, v, k), rs)

# Workspace entries of a chain, and how many of its reducers ask for one (at most one may, so no offsets are needed).
_workspace_len(::Tuple{}, ::Type, ::Val, ::Val) = (0, 0)
function _workspace_len(rs::Tuple, ::Type{T}, v::Val, k::Val) where {T}
    m = reducer_workspace(first(rs), T, v, k)
    m2, c2 = _workspace_len(Base.tail(rs), T, v, k)
    return (m + m2, (m > 0) + c2)
end
function _chain_workspace(rs::Tuple, ::Type{T}, v::Val, k::Val) where {T}
    m, c = _workspace_len(rs, T, v, k)
    c <= 1 || throw(ArgumentError("at most one reducer of a chain may request a workspace, got $c"))
    return m
end

const _NO_WS = UInt32[]   # handed to chains that request no workspace; never written

# Apply the chain in order at one handoff; each reducer sees the W replaced by the one before it and the orbit's workspace.
_apply(::Tuple{}, ::Tuple{}, n, x, W, ws) = ((), W)
function _apply(rs::Tuple, ss::Tuple, n, x, W, ws)
    s, W1 = reducer_step(first(rs), first(ss), n, x, W, ws)
    rest, W2 = _apply(Base.tail(rs), Base.tail(ss), n, x, W1, ws)
    return ((s, rest...), W2)
end

# First failed reducer in the chain, as a host-side error (the chain is a Tuple; recursion keeps it type stable).
_any_failed(::Tuple{}, ::Tuple{}) = false
_any_failed(rs::Tuple, ss::Tuple) = reducer_failed(first(rs), first(ss)) || _any_failed(Base.tail(rs), Base.tail(ss))
function _reducer_error(rs::Tuple, ss::Tuple, orbit, backend)
    reducer_failed(first(rs), first(ss)) && return reducer_error(first(rs), first(ss), orbit, backend)
    return _reducer_error(Base.tail(rs), Base.tail(ss), orbit, backend)
end

# The single per-orbit body (spec I1, I2): tangent at the OLD point, then advance x; handoff every tau steps.
# The state has s >= d components (D28 winding counts ride along), W is d x k with d the system dimension.
# One launch (D44): `steps` steps (a multiple of tau) from the carried point, matrix, reducer states and global step n0;
# a run of N steps is the call with n0 = 0, steps = N and fresh states. The final status check runs only at step N.
# Non-throwing so the same code runs in the kernel: returns (x, W, states, status, global step of the bad handoff).
# `ws` is the orbit's reducer workspace (zeroed by the caller at n0 == 0, carried between launches like the states).
@inline function _orbit_body(sys, x0::SVector, W0::SMatrix, states0::Tuple, n0::Int, steps::Int, N::Int, tau::Int, reducers::Tuple,
                             ws::AbstractVector{UInt32} = _NO_WS)
    x, W, states = x0, W0, states0
    n = n0
    for _ in 1:(steps ÷ tau)
        for _ in 1:tau
            W = fused_product(tangent(sys, x), W)   # explicit fma for Float32 and Float64 SMatrix tangents, else `*` (D115)
            x = forward(sys, x)
        end
        n += tau
        code = _status(x, W)
        code == STATUS_OK || return (x, W, states, code, n)
        states, W = _apply(reducers, states, n, x, W, ws)
        _any_failed(reducers, states) && return (x, W, states, STATUS_REDUCER, n)
    end
    return (x, W, states, n == N ? _status(x, W) : STATUS_OK, n)
end

# First non-finite entry as Float64 (NaN or Inf), for the error record.
_first_bad(v) = Float64(v[findfirst(!isfinite, v)])

function _health_error(code::Integer, x, W, step::Integer, orbit::Integer, backend::Symbol)
    code == STATUS_STATE ? NumericalHealthError("state", _first_bad(x), Inf, step, orbit, backend) :
                           NumericalHealthError("deviation matrix", _first_bad(W), Inf, step, orbit, backend)
end

function _check_steps(N::Integer, tau::Integer)
    (tau >= 1 && N >= 1 && N % tau == 0) ||
        throw(ArgumentError("N must be a positive multiple of tau, got N = $N, tau = $tau"))
    return nothing
end

"""
    propagate_orbit(sys, x0, W0, N, tau, reducers; orbit=1, backend=:cpu, workspace=nothing) -> (x, W, states)

Per-orbit body shared with the kernel: `N` steps of `sys` from `x0` with deviation matrix `W0`, the reducer chain
applied at steps tau, 2 tau, ...; throws `NumericalHealthError` when `x` or `W` is non-finite at a handoff. A chain whose
reducer asks for a workspace ([`reducer_workspace`](@ref)) uses `workspace`, a `UInt32` vector of exactly that length,
zeroed here before the first step and left holding the counts for the three-argument [`reducer_final`](@ref); `nothing`
allocates one. A wrong length, or two workspace reducers, is an `ArgumentError` naming the workspace before any work.
"""
function propagate_orbit(sys, x0::SVector{s, T}, W0::SMatrix{d, k}, N::Integer, tau::Integer, reducers::Tuple;
                         orbit::Integer = 1, backend::Symbol = :cpu, workspace::Union{Nothing, AbstractVector{UInt32}} = nothing) where {s, T, d, k}
    _check_steps(N, tau)
    m = _chain_workspace(reducers, T, Val(d), Val(k))
    ws = workspace === nothing ? (m == 0 ? _NO_WS : zeros(UInt32, m)) : workspace
    length(ws) == m || throw(ArgumentError("workspace has length $(length(ws)), the reducer chain needs $m"))
    for i in eachindex(ws)
        ws[i] = 0x00000000
    end
    s0 = _init_states(reducers, T, Val(d), Val(k))
    x, W, states, code, step = _orbit_body(sys, x0, W0, s0, 0, Int(N), Int(N), Int(tau), reducers, ws)
    code == STATUS_REDUCER && throw(_reducer_error(reducers, states, orbit, backend))
    code == STATUS_OK || throw(_health_error(code, x, W, step, orbit, backend))
    return (x, W, states)
end
