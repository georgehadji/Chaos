"""
    LinearSymplecticMap(A)

Linear symplectic map `z' = A z` of dimension `D = 2N`, ordering `(q_1..q_N, p_1..p_N)`; the tangent is `A`
for every state. `Ainv = Omega^-1 A' Omega` is stored (spec item 10).
"""
struct LinearSymplecticMap{T, D, L} <: AbstractMap
    A::SMatrix{D, D, T, L}
    Ainv::SMatrix{D, D, T, L}
end

# Omega^-1 A' Omega with Omega = [0 I; -I 0], written out in blocks: A = [a b; c d] gives [d' -b'; -c' a']
function _omega_inverse(A::SMatrix{D, D, T}) where {D, T}
    iseven(D) || throw(ArgumentError("symplectic map needs an even dimension, got $D"))
    N = div(D, 2)
    return SMatrix{D, D, T}(ntuple(Val(D * D)) do l
        i = (l - 1) % D + 1
        j = div(l - 1, D) + 1
        if i <= N
            j <= N ? A[N + j, N + i] : -A[j - N, N + i]
        else
            j <= N ? -A[N + j, i - N] : A[j - N, i - N]
        end
    end)
end

LinearSymplecticMap(A::SMatrix{D, D, T, L}) where {D, T, L} = LinearSymplecticMap{T, D, L}(A, _omega_inverse(A))

dimension(::LinearSymplecticMap{T, D}) where {T, D} = D
forward(m::LinearSymplecticMap{T, D}, z::SVector{D}) where {T, D} = m.A * z
inverse(m::LinearSymplecticMap{T, D}, z::SVector{D}) where {T, D} = m.Ainv * z
tangent(m::LinearSymplecticMap{T, D}, ::SVector{D}) where {T, D} = m.A
inverse_tangent(m::LinearSymplecticMap{T, D}, ::SVector{D}) where {T, D} = m.Ainv

"""
    rotation_block(T, theta)

2 x 2 rotation `[[cos, sin], [-sin, cos]]` in element type `T`; exponents 0 and 0.
"""
function rotation_block(::Type{T}, theta) where {T}
    s, c = sincos(T(theta))
    return SMatrix{2, 2, T}(c, -s, s, c)
end

"""
    hyperbolic_block(T, mu)

2 x 2 hyperbolic block `[[cosh, sinh], [sinh, cosh]]` in element type `T`; eigenvalues `exp(+-mu)`.
"""
function hyperbolic_block(::Type{T}, mu) where {T}
    m = T(mu)
    return SMatrix{2, 2, T}(cosh(m), sinh(m), sinh(m), cosh(m))
end

# entry l (column-major) of the (q_1..q_N, p_1..p_N) direct sum of 2 x 2 blocks
function _bd_entry(blocks, l, N)
    D = 2N
    i = (l - 1) % D + 1
    j = div(l - 1, D) + 1
    a, r = i <= N ? (i, 1) : (i - N, 2)
    b, c = j <= N ? (j, 1) : (j - N, 2)
    return a == b ? blocks[a][r, c] : zero(eltype(blocks[1]))
end

_block_diagonal(blocks, ::Val{N}) where {N} =
    SMatrix{2N, 2N, eltype(blocks[1])}(ntuple(l -> _bd_entry(blocks, l, N), Val(4N^2)))

"""
    block_diagonal(blocks...)

Direct sum of 2 x 2 blocks acting on the pairs `(q_j, p_j)`, reordered to `(q_1..q_N, p_1..p_N)`.
"""
block_diagonal(b1::SMatrix{2, 2, T, 4}, bs::SMatrix{2, 2, T, 4}...) where {T} =
    _block_diagonal((b1, bs...), Val(1 + length(bs)))

"""
    conjugated(A, B)

Non-normal variant `S^-1 A S` with `S = [I B; 0 I]` and `B` symmetric `N x N`; symplectic when `A` is.
"""
function conjugated(A::SMatrix{D, D}, B::SMatrix{N, N}) where {D, N}
    D == 2N || throw(DimensionMismatch("A is $(D) x $(D) but B is $(N) x $(N)"))
    T = promote_type(eltype(A), eltype(B))
    S(sgn) = SMatrix{D, D, T}(ntuple(Val(D * D)) do l
        i = (l - 1) % D + 1
        j = div(l - 1, D) + 1
        i == j ? one(T) : (i <= N && j > N ? sgn * T(B[i, j - N]) : zero(T))
    end)
    return S(-1) * A * S(1)
end
