"""
    WrappedStandardMapChain(N, K, eps; boundary = :open)

[`StandardMapChain`](@ref) on the torus with winding counts (decision D28), state a vector of length `4N`:
`(q_1..q_N, p_1..p_N, count q_1..q_N, count p_1..p_N)`, coordinates in `[0, 2pi_T)` and exact integer counts.
`dimension` is `2N` (the tangent size); `tangent` and `inverse_tangent` read the first `2N` entries. If a count reaches
`maxintfloat(T)` a step returns NaN in all `4N` components. The in-place forms allocate nothing.
"""
struct WrappedStandardMapChain{T} <: AbstractMap
    chain::StandardMapChain{T}
end

function WrappedStandardMapChain(N::Integer, K, eps; boundary::Symbol = :open)
    return WrappedStandardMapChain(StandardMapChain(N, K, eps; boundary))
end

dimension(w::WrappedStandardMapChain) = dimension(w.chain)

# Counts out[c0+1:c0+n] at the exactness limit (or NaN) poison the whole state, so the health check fires.
function _guard_counts!(out, c0, n)
    lim = maxintfloat(eltype(out))
    for i in 1:n
        abs(out[c0 + i]) < lim || return fill!(out, oftype(out[1], NaN))
    end
    return out
end

# W-E2 per site: momenta first (force at the pre-step reduced positions), then positions; m' = m + n' + b, n' = n + a
function forward!(out, w::WrappedStandardMapChain, s)
    m = w.chain
    N = m.N
    tp = 2 * _pi(s[1])
    for j in 1:N
        out[N + j], a = _reduce(s[N + j] + _force(m, s, j), tp)
        out[3N + j] = s[3N + j] + a
    end
    for j in 1:N
        out[j], b = _reduce(s[j] + out[N + j], tp)
        out[2N + j] = s[2N + j] + out[3N + j] + b
    end
    return _guard_counts!(out, 2N, 2N)
end

# W-E4: positions x = x' - p' (b), then p = p' - F(x) (a); m = m' - n' + b, n = n' + a with the image count n'
function inverse!(out, w::WrappedStandardMapChain, s)
    m = w.chain
    N = m.N
    tp = 2 * _pi(s[1])
    for j in 1:N
        out[j], b = _reduce(s[j] - s[N + j], tp)
        out[2N + j] = s[2N + j] - s[3N + j] + b
    end
    for j in 1:N
        out[N + j], a = _reduce(s[N + j] - _force(m, out, j), tp)
        out[3N + j] = s[3N + j] + a
    end
    return _guard_counts!(out, 2N, 2N)
end

forward(w::WrappedStandardMapChain, s::AbstractVector) = forward!(similar(s), w, s)
inverse(w::WrappedStandardMapChain, s::AbstractVector) = inverse!(similar(s), w, s)

# W-E3: the tangent depends on the reduced positions only
tangent!(M, w::WrappedStandardMapChain, s) = tangent!(M, w.chain, s)
inverse_tangent!(M, w::WrappedStandardMapChain, s) = inverse_tangent!(M, w.chain, s)
tangent(w::WrappedStandardMapChain, s::AbstractVector) = tangent!(_newmat(w.chain, s), w.chain, s)
inverse_tangent(w::WrappedStandardMapChain, s::AbstractVector) = inverse_tangent!(_newmat(w.chain, s), w.chain, s)

"""
    torus_state(sys::WrappedStandardMapChain, z) -> Vector

State of the wrapped chain for the lifted point `z` (length `2N`): each coordinate reduced into `[0, 2pi_T)` in `T`
with its integer count, laid out `(reduced coordinates, counts)`.
"""
function torus_state(w::WrappedStandardMapChain, z::AbstractVector)
    d = dimension(w)
    tp = 2 * _pi(z[1])
    s = similar(z, 2d)
    for j in 1:d
        s[j], s[d + j] = _reduce(z[j], tp)
    end
    return s
end

"""
    lift(sys::WrappedStandardMapChain, s) -> Vector

Lifted point `z_j + 2pi_T c_j` of the wrapped state `s` (length `4N`), computed in its element type.
"""
function lift(w::WrappedStandardMapChain, s::AbstractVector)
    d = dimension(w)
    tp = 2 * _pi(s[1])
    return [s[j] + tp * s[d + j] for j in 1:d]
end
