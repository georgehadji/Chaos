"""
    FroeschleMap(nu, kappa, mu)

Froeschle 4D map in the form of Skokos2010 Eq. (16) (coupled-maps spec item 1), lifted coordinates, state
`SVector(x1, x2, p1, p2)` (positions first): `x' = x + p`, `p1' = p1 - nu sin x1' - mu (1 - cos(x1' + x2'))`,
`p2' = p2 - kappa sin x2' - mu (1 - cos(x1' + x2'))`. Tangent `[[I, I], [A, I + A]]` with `A` symmetric at the new
positions. Parameters must be built in the state's element type (use `exact`).
"""
struct FroeschleMap{T} <: AbstractMap
    nu::T
    kappa::T
    mu::T
end

dimension(::FroeschleMap) = 4

# Coupling force common to both momenta; the operation order matches `StandardMap` so that mu = 0 is bitwise Y1 (spec item 6).
@inline _fcoup(mu, x1, x2) = mu * (one(x1) - _cos(x1 + x2))

function forward(m::FroeschleMap, z::SVector{4})
    x1 = z[1] + z[3]
    x2 = z[2] + z[4]
    f = _fcoup(m.mu, x1, x2)
    return SVector(x1, x2, z[3] - m.nu * _sin(x1) - f, z[4] - m.kappa * _sin(x2) - f)
end

# spec item 5: p = p' + nu sin x1' + mu (1 - cos(x1' + x2')), x = x' - p
function inverse(m::FroeschleMap, z::SVector{4})
    f = _fcoup(m.mu, z[1], z[2])
    p1 = z[3] + m.nu * _sin(z[1]) + f
    p2 = z[4] + m.kappa * _sin(z[2]) + f
    return SVector(z[1] - p1, z[2] - p2, p1, p2)
end

# spec item 4: entries a, b, c of A at the new positions
@inline function _abc(m, x1, x2)
    t = m.mu * _sin(x1 + x2)
    return -m.nu * _cos(x1) - t, -t, -m.kappa * _cos(x2) - t
end

# M = [[I, I], [A, I + A]]; column-major SMatrix entries
function _f3_tangent(m, x1, x2)
    a, b, c = _abc(m, x1, x2)
    o = one(a)
    Z = zero(a)
    return SMatrix{4, 4}(o, Z, a, b, Z, o, b, c, o, Z, o + a, b, Z, o, b, o + c)
end

# M^-1 = [[I + A, -I], [-A, I]] with A at the image positions (spec item 5)
function _f3_inverse_tangent(m, x1, x2)
    a, b, c = _abc(m, x1, x2)
    o = one(a)
    Z = zero(a)
    return SMatrix{4, 4}(o + a, b, -a, -b, b, o + c, -b, -c, -o, Z, o, Z, Z, -o, Z, o)
end

tangent(m::FroeschleMap, z::SVector{4}) = _f3_tangent(m, z[1] + z[3], z[2] + z[4])
inverse_tangent(m::FroeschleMap, z::SVector{4}) = _f3_inverse_tangent(m, z[1], z[2])

"""
    WrappedFroeschleMap(nu, kappa, mu)

[`FroeschleMap`](@ref) on the torus with winding counts (decision D28), state `SVector(x1, x2, p1, p2, c1, c2, c3, c4)`:
coordinates in `[0, 2pi_T)` and one exact integer count per coordinate. `dimension` is 4 (the tangent size). If a count
reaches `maxintfloat(T)` a step returns NaN in all eight components. Parameters must be built in the state's element type.
"""
struct WrappedFroeschleMap{T} <: AbstractMap
    nu::T
    kappa::T
    mu::T
end

dimension(::WrappedFroeschleMap) = 4

@inline function _checked8(v::SVector{8})
    lim = maxintfloat(typeof(v[1]))
    (abs(v[5]) < lim && abs(v[6]) < lim && abs(v[7]) < lim && abs(v[8]) < lim) && return v
    return SVector{8}(ntuple(_ -> oftype(v[1], NaN), Val(8)))
end

# spec Algorithm item 3: new positions reduced, then the new momenta; m' = m + n + b, n' = n + a
function forward(m::WrappedFroeschleMap, s::SVector{8})
    tp = 2 * _pi(s[1])
    x1, b1 = _reduce(s[1] + s[3], tp)
    x2, b2 = _reduce(s[2] + s[4], tp)
    f = _fcoup(m.mu, x1, x2)
    p1, a1 = _reduce(s[3] - m.nu * _sin(x1) - f, tp)
    p2, a2 = _reduce(s[4] - m.kappa * _sin(x2) - f, tp)
    return _checked8(SVector(x1, x2, p1, p2, s[5] + s[7] + b1, s[6] + s[8] + b2, s[7] + a1, s[8] + a2))
end

# mirror of the above: momenta first (a), then positions (b); m = m' - n + b with the preimage count n
function inverse(m::WrappedFroeschleMap, s::SVector{8})
    tp = 2 * _pi(s[1])
    f = _fcoup(m.mu, s[1], s[2])
    p1, a1 = _reduce(s[3] + m.nu * _sin(s[1]) + f, tp)
    p2, a2 = _reduce(s[4] + m.kappa * _sin(s[2]) + f, tp)
    n1, n2 = s[7] + a1, s[8] + a2
    x1, b1 = _reduce(s[1] - p1, tp)
    x2, b2 = _reduce(s[2] - p2, tp)
    return _checked8(SVector(x1, x2, p1, p2, s[5] - n1 + b1, s[6] - n2 + b2, n1, n2))
end

tangent(m::WrappedFroeschleMap, s::SVector{8}) = _f3_tangent(m, s[1] + s[3], s[2] + s[4])
inverse_tangent(m::WrappedFroeschleMap, s::SVector{8}) = _f3_inverse_tangent(m, s[1], s[2])

"""
    torus_state(sys::WrappedFroeschleMap, z) -> SVector(x1, x2, p1, p2, c1, c2, c3, c4)

Each coordinate of the lifted point `z` reduced into `[0, 2pi_T)` in `T` with its integer count.
"""
function torus_state(::WrappedFroeschleMap, z::SVector{4})
    tp = 2 * _pi(z[1])
    r = map(v -> _reduce(v, tp), z)
    return SVector{8}(ntuple(i -> i <= 4 ? r[i][1] : r[i - 4][2], Val(8)))
end

"""
    lift(sys::WrappedFroeschleMap, s) -> SVector{4}

Lifted point `z_j + 2pi_T c_j` of the wrapped state `s`, computed in its element type.
"""
function lift(::WrappedFroeschleMap, s::SVector{8})
    tp = 2 * _pi(s[1])
    return SVector{4}(ntuple(j -> s[j] + tp * s[4 + j], Val(4)))
end
