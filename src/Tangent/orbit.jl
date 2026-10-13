const STATUS_OK, STATUS_STATE, STATUS_W, STATUS_REDUCER = Int8(0), Int8(1), Int8(2), Int8(3)

_status(x, W) = all(isfinite, x) ? (all(isfinite, W) ? STATUS_OK : STATUS_W) : STATUS_STATE

_init_states(rs::Tuple, ::Type{T}, v::Val, k::Val) where {T} = map(r -> reducer_init(r, T, v, k), rs)

# Apply the chain in order at one handoff; each reducer sees the W replaced by the one before it.
_apply(::Tuple{}, ::Tuple{}, n, x, W) = ((), W)
function _apply(rs::Tuple, ss::Tuple, n, x, W)
    s, W1 = reducer_step(first(rs), first(ss), n, x, W)
    rest, W2 = _apply(Base.tail(rs), Base.tail(ss), n, x, W1)
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
# Non-throwing so the same code runs in the kernel: returns (x, W, states, status, step of the bad handoff).
@inline function _orbit_body(sys, x0::SVector{s, T}, W0::SMatrix{d, k}, N::Int, tau::Int, reducers::Tuple) where {s, d, k, T}
    x, W = x0, W0
    states = _init_states(reducers, T, Val(d), Val(k))
    n = 0
    for _ in 1:(N ÷ tau)
        for _ in 1:tau
            W = tangent(sys, x) * W
            x = forward(sys, x)
        end
        n += tau
        code = _status(x, W)
        code == STATUS_OK || return (x, W, states, code, n)
        states, W = _apply(reducers, states, n, x, W)
        _any_failed(reducers, states) && return (x, W, states, STATUS_REDUCER, n)
    end
    return (x, W, states, _status(x, W), N)
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
    propagate_orbit(sys, x0, W0, N, tau, reducers; orbit=1, backend=:cpu) -> (x, W, states)

Per-orbit body shared with the kernel: `N` steps of `sys` from `x0` with deviation matrix `W0`, the reducer chain
applied at steps tau, 2 tau, ...; throws `NumericalHealthError` when `x` or `W` is non-finite at a handoff.
"""
function propagate_orbit(sys, x0::SVector, W0::SMatrix, N::Integer, tau::Integer, reducers::Tuple;
                         orbit::Integer = 1, backend::Symbol = :cpu)
    _check_steps(N, tau)
    x, W, states, code, step = _orbit_body(sys, x0, W0, Int(N), Int(tau), reducers)
    code == STATUS_REDUCER && throw(_reducer_error(reducers, states, orbit, backend))
    code == STATUS_OK || throw(_health_error(code, x, W, step, orbit, backend))
    return (x, W, states)
end
