"""
    ChainTangent

Banded one-step tangent of the coupled standard-map chain at a fixed-size state (amendment D105): an isbits `2N x 2N`
`AbstractMatrix` holding the entries of `G = dF/dq` (diagonal, left and right neighbour), the flag of the boundary and
the flag of the inverse. Its entries equal the dense `tangent!` / `inverse_tangent!` (bitwise, `+0.0` off the band);
`M * W` for an `SMatrix` frame `W` with `2N` rows applies the band in the operation order of `tangent_action!` /
`inverse_tangent_action!` (cost linear in `N`) and returns an `SMatrix`.
"""
struct ChainTangent{T, N} <: AbstractMatrix{T}
    g::SVector{N, T}
    gl::SVector{N, T}
    gr::SVector{N, T}
    ring::Bool
    inv::Bool
end

Base.size(::ChainTangent{T, N}) where {T, N} = (2N, 2N)

@inline _left(c::ChainTangent{T, N}, j) where {T, N} = j > 1 ? j - 1 : (c.ring ? N : 0)
@inline _right(c::ChainTangent{T, N}, j) where {T, N} = j < N ? j + 1 : (c.ring ? 1 : 0)

# Entry (j, k) of the band G (+0.0 off the band, as the zero-filled dense block).
@inline function _gentry(c::ChainTangent{T}, j, k) where {T}
    k == j && return c.g[j]
    k == _left(c, j) && return c.gl[j]
    k == _right(c, j) && return c.gr[j]
    return zero(T)
end

# Same expressions as `tangent!` ([[I + G, I], [G, I]]) and `inverse_tangent!` ([[I, -I], [-G, I + G]]).
function Base.getindex(c::ChainTangent{T, N}, i::Int, l::Int) where {T, N}
    @boundscheck checkbounds(c, i, l)
    j, k = (i - 1) % N + 1, (l - 1) % N + 1
    G, e = _gentry(c, j, k), j == k ? one(T) : zero(T)
    top, left = i <= N, l <= N
    c.inv && return top ? (left ? e : -e) : (left ? -G : G + e)
    return top ? (left ? G + e : e) : (left ? G : e)
end

# Rows of G X for the SMatrix X (N x k): g X_j + gl X_left + gr X_right, left to right as in the dense-free actions.
function _band(c::ChainTangent{T, N}, X::SMatrix{N, k}) where {T, N, k}
    return SMatrix{N, k}(ntuple(Val(N * k)) do l
        j, col = (l - 1) % N + 1, (l - 1) ÷ N + 1
        kl, kr = _left(c, j), _right(c, j)
        v = c.g[j] * X[j, col]
        kl > 0 && (v += c.gl[j] * X[kl, col])
        kr > 0 && (v += c.gr[j] * X[kr, col])
        v
    end)
end

function Base.:*(c::ChainTangent{T, N}, W::SMatrix{D, k}) where {T, N, D, k}
    D == 2N || throw(DimensionMismatch("frame rows do not match the tangent size 2N"))
    Wq, Wp = _rows(W, 0, Val(N)), _rows(W, N, Val(N))
    if c.inv   # top = W_q - W_p, bottom = W_p + G (W_p - W_q)
        return vcat(Wq - Wp, Wp + _band(c, Wp - Wq))
    end
    GW = _band(c, Wq)   # top = W_q + G W_q + W_p, bottom = G W_q + W_p
    return vcat(Wq + GW + Wp, GW + Wp)
end

# SVector of the entries s[o + 1 .. o + N] of a fixed-size state.
@inline _blk(s, o, ::Val{N}) where {N} = SVector{N}(ntuple(j -> s[o + j], Val(N)))

# The band of G at the positions q (indexable, at least N entries): diagonal K cos q_j + eps cos(left) + eps cos(right),
# off-diagonal minus those cosines, as `_grow` (one entry per site and kind, no heterogeneous tuples).
function _chain_tangent(m, q, ::Val{N}, ::Type{T}, inverse::Bool) where {N, T}
    m.N == N || throw(DimensionMismatch("state sites do not match the chain"))
    g0 = SVector{N, T}(ntuple(j -> m.K * _cos(q[j]), Val(N)))
    cl = SVector{N, T}(ntuple(j -> (kl = _left(m, j); kl > 0 ? m.eps * _cos(q[j] - q[kl]) : zero(T)), Val(N)))
    cr = SVector{N, T}(ntuple(j -> (kr = _right(m, j); kr > 0 ? m.eps * _cos(q[j] - q[kr]) : zero(T)), Val(N)))
    return ChainTangent{T, N}(g0 + cl + cr, -cl, -cr, m.ring, inverse)
end

# Error messages here are constant: string interpolation allocates and does not compile in a GPU kernel (D105d).
# Fixed-size states (SVector of length 2N): the step in blocks, same expressions and order as the in-place forms.
function forward(m::StandardMapChain, z::SVector{S}) where {S}
    N = S ÷ 2
    iseven(S) && m.N == N || throw(DimensionMismatch("state length is not 2N for this chain"))
    T = _promoted(m, z)
    p = SVector{N, T}(ntuple(j -> z[N + j] + _force(m, z, j), Val(N)))
    return vcat(SVector{N, T}(ntuple(j -> z[j] + p[j], Val(N))), p)
end

function inverse(m::StandardMapChain, z::SVector{S}) where {S}
    N = S ÷ 2
    iseven(S) && m.N == N || throw(DimensionMismatch("state length is not 2N for this chain"))
    T = _promoted(m, z)
    q = SVector{N, T}(ntuple(j -> z[j] - z[N + j], Val(N)))
    return vcat(q, SVector{N, T}(ntuple(j -> z[N + j] - _force(m, q, j), Val(N))))
end

tangent(m::StandardMapChain, z::SVector{S}) where {S} = _chain_tangent(m, z, Val(S ÷ 2), _promoted(m, z), false)
function inverse_tangent(m::StandardMapChain, z::SVector{S}) where {S}
    return _chain_tangent(m, _Preimage(z, S ÷ 2), Val(S ÷ 2), _promoted(m, z), true)
end
