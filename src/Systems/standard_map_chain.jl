"""
    StandardMapChain(N, K, eps; boundary = :open)

Chain of `N` coupled standard maps, `2 <= N <= 100` (coupled-maps spec item 8), lifted coordinates, state a vector
`(q_1..q_N, p_1..p_N)` of length `2N`: `p_j' = p_j + F_j(q)`, `q_j' = q_j + p_j'` with
`F_j = K sin q_j + eps sum_k sin(q_j - q_k)` over the neighbours `k` of `j`. `boundary` is `:open` (free ends) or `:ring`
(`N >= 3`). `K` and `eps` are stored in one element type (use `exact`); the state is any `AbstractVector`.
Tangent `[[I + G, I], [G, I]]` with the symmetric banded `G = dF/dq` at the pre-step positions. The in-place forms
`forward!`, `inverse!`, `tangent!`, `inverse_tangent!`, `tangent_action!`, `inverse_tangent_action!` allocate nothing.
"""
struct StandardMapChain{T} <: AbstractMap
    N::Int
    K::T
    eps::T
    ring::Bool
end

function StandardMapChain(N::Integer, K, eps; boundary::Symbol = :open)
    2 <= N <= 100 || throw(ArgumentError("number of sites N must be in 2:100, got $N"))
    boundary in (:open, :ring) || throw(ArgumentError("boundary must be :open or :ring, got :$boundary"))
    boundary === :ring && N < 3 && throw(ArgumentError("a ring needs N >= 3, got $N"))
    K, eps = promote(K, eps)
    return StandardMapChain(Int(N), K, eps, boundary === :ring)
end

dimension(m::StandardMapChain) = 2 * m.N

# Neighbours of site j (0 = none): the ring joins site N to site 1.
@inline _left(m, j) = j > 1 ? j - 1 : (m.ring ? m.N : 0)
@inline _right(m, j) = j < m.N ? j + 1 : (m.ring ? 1 : 0)

# F_j(q) with the operation order of `StandardMap`: K sin q_j + eps (sum of neighbour sines); eps = 0 gives Y1 bitwise.
@inline function _force(m, q, j)
    kl, kr = _left(m, j), _right(m, j)
    s = zero(m.K * q[j])
    kl > 0 && (s += _sin(q[j] - q[kl]))
    kr > 0 && (s += _sin(q[j] - q[kr]))
    return m.K * _sin(q[j]) + m.eps * s
end

# Row j of G(q): neighbours, off-diagonal entries -eps cos(q_j - q_k) (0 if none), diagonal K cos q_j + eps sum cos.
@inline function _grow(m, q, j)
    kl, kr = _left(m, j), _right(m, j)
    g = m.K * _cos(q[j])
    cl = kl > 0 ? m.eps * _cos(q[j] - q[kl]) : zero(g)
    cr = kr > 0 ? m.eps * _cos(q[j] - q[kr]) : zero(g)
    return kl, kr, -cl, -cr, g + cl + cr
end

# Positions of the preimage of the image point z: q = q' - p' (spec item 9); indexable like a vector, no allocation.
struct _Preimage{Z}
    z::Z
    N::Int
end
@inline Base.getindex(p::_Preimage, j::Integer) = p.z[j] - p.z[p.N + j]

_promoted(m::StandardMapChain, z) = promote_type(eltype(z), typeof(m.K), typeof(m.eps))
_new(m::StandardMapChain, z) = similar(z, _promoted(m, z))
_newmat(m::StandardMapChain, z) = zeros(_promoted(m, z), 2 * m.N, 2 * m.N)

function forward!(out, m::StandardMapChain, z)
    N = m.N
    for j in 1:N
        out[N + j] = z[N + j] + _force(m, z, j)
    end
    for j in 1:N
        out[j] = z[j] + out[N + j]
    end
    return out
end

# spec item 9: q = q' - p', p = p' - F(q)
function inverse!(out, m::StandardMapChain, z)
    N = m.N
    for j in 1:N
        out[j] = z[j] - z[N + j]
    end
    for j in 1:N
        out[N + j] = z[N + j] - _force(m, out, j)
    end
    return out
end

forward(m::StandardMapChain, z::AbstractVector) = forward!(_new(m, z), m, z)
inverse(m::StandardMapChain, z::AbstractVector) = inverse!(_new(m, z), m, z)

# G(q) into the block M[r0+1:r0+N, 1:N]
function _fill_G!(M, m, q, r0)
    N = m.N
    for j in 1:N, k in 1:N
        M[r0 + j, k] = zero(eltype(M))
    end
    for j in 1:N
        kl, kr, gl, gr, g = _grow(m, q, j)
        M[r0 + j, j] = g
        kl > 0 && (M[r0 + j, kl] = gl)
        kr > 0 && (M[r0 + j, kr] = gr)
    end
    return M
end

# M = [[I + G, I], [G, I]]
function tangent!(M, m::StandardMapChain, z)
    N = m.N
    o, Z = one(eltype(M)), zero(eltype(M))
    _fill_G!(M, m, z, N)
    for j in 1:N, k in 1:N
        e = j == k ? o : Z
        M[j, k] = M[N + j, k] + e
        M[j, N + k] = e
        M[N + j, N + k] = e
    end
    return M
end

# M^-1 = [[I, -I], [-G, I + G]] with G at the preimage
function inverse_tangent!(M, m::StandardMapChain, z)
    N = m.N
    o, Z = one(eltype(M)), zero(eltype(M))
    _fill_G!(M, m, _Preimage(z, N), N)
    for j in 1:N, k in 1:N
        e = j == k ? o : Z
        g = M[N + j, k]
        M[j, k] = e
        M[j, N + k] = -e
        M[N + j, k] = -g
        M[N + j, N + k] = g + e
    end
    return M
end

tangent(m::StandardMapChain, z::AbstractVector) = tangent!(_newmat(m, z), m, z)
inverse_tangent(m::StandardMapChain, z::AbstractVector) = inverse_tangent!(_newmat(m, z), m, z)

# out = M W, row by row with the banded G: top = W_q + G W_q + W_p, bottom = G W_q + W_p
function tangent_action!(out, m::StandardMapChain, z, W)
    N = m.N
    for j in 1:N
        kl, kr, gl, gr, g = _grow(m, z, j)
        for c in axes(W, 2)
            gw = g * W[j, c]
            kl > 0 && (gw += gl * W[kl, c])
            kr > 0 && (gw += gr * W[kr, c])
            out[N + j, c] = gw + W[N + j, c]
            out[j, c] = W[j, c] + gw + W[N + j, c]
        end
    end
    return out
end

# out = M^-1 W: top = W_q - W_p, bottom = W_p + G (W_p - W_q) with G at the preimage
function inverse_tangent_action!(out, m::StandardMapChain, z, W)
    N = m.N
    q = _Preimage(z, N)
    for j in 1:N
        kl, kr, gl, gr, g = _grow(m, q, j)
        for c in axes(W, 2)
            gd = g * (W[N + j, c] - W[j, c])
            kl > 0 && (gd += gl * (W[N + kl, c] - W[kl, c]))
            kr > 0 && (gd += gr * (W[N + kr, c] - W[kr, c]))
            out[j, c] = W[j, c] - W[N + j, c]
            out[N + j, c] = W[N + j, c] + gd
        end
    end
    return out
end
