"""
    WrappedStandardMap(K)

Standard map on the torus with winding counts (decision D28), state `SVector(x, p, m, n)`: `x`, `p` in `[0, 2pi_T)`,
integer counts `m`, `n` stored exactly in `T`; the lifted point is `(x + 2pi m, p + 2pi n)`. `dimension` is 2 (the
tangent is 2 x 2). If `|m|` or `|n|` reaches `maxintfloat(T)` a step returns NaN in all four components. `K` must be
built in the state's element type.
"""
struct WrappedStandardMap{T} <: AbstractMap
    K::T
end

dimension(::WrappedStandardMap) = 2

# Reduction into [0, tp) with the count (amendment D28, W-E1): floor, then the guards for a floor that is one off.
@inline function _reduce(u, tp)
    a = floor(u / tp)
    r = u - tp * a
    if r < 0
        r += tp
        a -= 1
    end
    if r >= tp
        r -= tp
        a += 1
    end
    return r, a
end

# Counts that reached the exactness limit (or are NaN) poison the whole state, so the health check of Tangent fires.
@inline function _checked(x, p, m, n)
    lim = maxintfloat(typeof(x))
    (abs(m) < lim && abs(n) < lim) && return SVector(x, p, m, n)
    nan = oftype(x, NaN)
    return SVector(nan, nan, nan, nan)
end

# W-E2: p~ = p + K sin x, x~ = x + p', each reduced with its count; m' = m + n' + b
function forward(m::WrappedStandardMap, z::SVector{4})
    tp = 2 * _pi(z[1])
    p, a = _reduce(z[2] + m.K * _sin(z[1]), tp)
    n = z[4] + a
    x, b = _reduce(z[1] + p, tp)
    return _checked(x, p, z[3] + n + b, n)
end

# W-E4: x~ = x' - p', p~ = p' - K sin x; m = m' - n' + b, n = n' + a
function inverse(m::WrappedStandardMap, z::SVector{4})
    tp = 2 * _pi(z[1])
    x, b = _reduce(z[1] - z[2], tp)
    p, a = _reduce(z[2] - m.K * _sin(x), tp)
    return _checked(x, p, z[3] - z[4] + b, z[4] + a)
end

# W-E3: cos has period 2 pi, so the tangent depends on x only
tangent(m::WrappedStandardMap, z::SVector{4}) = _dt(m.K, z[1])
inverse_tangent(m::WrappedStandardMap, z::SVector{4}) = _dtinv(m.K, z[1] - z[2])

"""
    torus_state(z) -> SVector(x, p, m, n)

State of [`WrappedStandardMap`](@ref) for the lifted point `z = (X, P)`: each coordinate reduced into `[0, 2pi_T)` in `T`
with its integer count (W-E5).
"""
function torus_state(z::SVector{2})
    tp = 2 * _pi(z[1])
    x, m = _reduce(z[1], tp)
    p, n = _reduce(z[2], tp)
    return SVector(x, p, m, n)
end

"""
    lift(s) -> SVector(X, P)

Lifted point `(x + 2pi_T m, p + 2pi_T n)` of the wrapped state `s`, computed in its element type; its absolute error is
about `|m| eps_T 2pi`, so for large counts widen `s` to `BigFloat` first.
"""
function lift(s::SVector{4})
    tp = 2 * _pi(s[1])
    return SVector(s[1] + tp * s[3], s[2] + tp * s[4])
end
