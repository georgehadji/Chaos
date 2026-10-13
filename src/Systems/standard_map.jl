"""
    StandardMap(K)

Chirikov-Taylor standard map on lifted coordinates, state `SVector(x, p)`:
`p' = p + K sin x`, `x' = x + p'` (spec item 1). Tangent `[[1 + c, 1], [c, 1]]` with `c = K cos x`.
`K` must be built in the state's element type (use [`exact`](@ref)).
"""
struct StandardMap{T} <: AbstractMap
    K::T
end

# sin(Inf) and cos(Inf) throw DomainError; non-finite input must give non-finite output (spec Scope).
_sin(x) = isfinite(x) ? sin(x) : x - x
_cos(x) = isfinite(x) ? cos(x) : x - x

dimension(::StandardMap) = 2

function forward(m::StandardMap, z::SVector{2})
    p = z[2] + m.K * _sin(z[1])
    return SVector(z[1] + p, p)
end

# spec item 4: x = x' - p', p = p' - K sin x
function inverse(m::StandardMap, z::SVector{2})
    x = z[1] - z[2]
    return SVector(x, z[2] - m.K * _sin(x))
end

function tangent(m::StandardMap, z::SVector{2})
    c = m.K * _cos(z[1])
    o = one(c)
    return SMatrix{2, 2}(o + c, c, o, o)
end

# c at the preimage x; spec item 4: DT^-1 = [[1, -1], [-c, 1 + c]]
function inverse_tangent(m::StandardMap, z::SVector{2})
    c = m.K * _cos(z[1] - z[2])
    o = one(c)
    return SMatrix{2, 2}(o, -c, -o, o + c)
end

# pi in the element type of x (BigFloat at the precision of x)
_pi(x) = oftype(x, pi)
_pi(x::BigFloat) = BigFloat(pi, precision(x))

# x mod 2pi in [0, 2pi), spec item 2; a result that rounds up to 2pi is the torus point 0.
# BigFloat: arithmetic uses the global precision, so callers wrap in `setprecision`.
function _wrap1(x)
    tp = 2 * _pi(x)
    y = x - tp * floor(x / tp)
    y < 0 && (y += tp)
    return y >= tp ? zero(y) : y
end

"""
    wrap(z)

Each component of `z` modulo 2pi into [0, 2pi), computed in the element type (2pi built in that type). For
display: the absolute error is about `|x| eps_T`, meaningless above `|x| ~ 1/eps_T`.
"""
wrap(z::SVector) = map(_wrap1, z)
