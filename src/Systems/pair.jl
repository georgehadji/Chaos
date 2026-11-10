"""
    PairMap(sys)

Product of the map `sys` with itself (rli spec RL5, P1(c)): the state is the concatenation (state A; state B) of two
states of `sys`, `dimension` is `2 dimension(sys)`, and each half evolves by `sys` alone. The one-step tangent is a
[`PairTangent`](@ref): block diagonal, and its action on a frame acts block by block (a dense `2d x 2d` product through
the zero blocks rounds differently in Double64 and would break the bitwise equality with the single-orbit run, D76).
"""
struct PairMap{M <: AbstractMap} <: AbstractMap
    sys::M
end

"""
    PairTangent(A, B)

Block-diagonal `2d x 2d` matrix `diag(A, B)` with `d x d` blocks `A`, `B` (isbits). It is an `AbstractMatrix` (so `size`,
indexing and `Matrix` work); `M * W` for an `SMatrix` `W` with `2d` rows multiplies each block with its own block rows of
`W`, as the `d x d` times `d x k` product of the single-orbit stream; the zero blocks are never multiplied.
"""
struct PairTangent{T, d, L} <: AbstractMatrix{T}
    A::SMatrix{d, d, T, L}
    B::SMatrix{d, d, T, L}
end

Base.size(::PairTangent{T, d}) where {T, d} = (2d, 2d)
function Base.getindex(M::PairTangent{T, d}, i::Int, j::Int) where {T, d}
    @boundscheck checkbounds(M, i, j)
    i <= d && j <= d && return M.A[i, j]
    i > d && j > d && return M.B[i - d, j - d]
    return zero(T)
end

# Rows `r0 + 1 .. r0 + d` of the `2d x k` matrix `W`.
_rows(W::SMatrix{D, k}, r0, ::Val{d}) where {D, k, d} = SMatrix{d, k}(ntuple(l -> W[r0 + (l - 1) % d + 1, (l - 1) ÷ d + 1], Val(d * k)))

function Base.:*(M::PairTangent{T, d}, W::SMatrix{D, k}) where {T, d, D, k}
    D == 2d || throw(DimensionMismatch("tangent of size $(2d) times a matrix with $D rows"))
    RA, RB = fused_product(M.A, _rows(W, 0, Val(d))), fused_product(M.B, _rows(W, d, Val(d)))   # D115
    return SMatrix{D, k}(ntuple(l -> (i = (l - 1) % D + 1; j = (l - 1) ÷ D + 1; i <= d ? RA[i, j] : RB[i - d, j]), Val(D * k)))
end

dimension(P::PairMap) = 2 * dimension(P.sys)

# The two halves of a pair state (any element type, any state length of the base map).
@inline function _halves(z::SVector{S, T}) where {S, T}
    iseven(S) || throw(DimensionMismatch("a pair state needs an even length, got $S"))   # folds on S (D77 item 5)
    return (SVector{S ÷ 2, T}(ntuple(i -> z[i], Val(S ÷ 2))), SVector{S ÷ 2, T}(ntuple(i -> z[S ÷ 2 + i], Val(S ÷ 2))))
end

function forward(P::PairMap, z::SVector)
    a, b = _halves(z)
    return vcat(forward(P.sys, a), forward(P.sys, b))
end

function inverse(P::PairMap, z::SVector)
    a, b = _halves(z)
    return vcat(inverse(P.sys, a), inverse(P.sys, b))
end

function tangent(P::PairMap, z::SVector)
    a, b = _halves(z)
    return PairTangent(tangent(P.sys, a), tangent(P.sys, b))
end

function inverse_tangent(P::PairMap, z::SVector)
    a, b = _halves(z)
    return PairTangent(inverse_tangent(P.sys, a), inverse_tangent(P.sys, b))
end
