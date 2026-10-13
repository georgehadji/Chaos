"""
    symplectic_form(T, d)

`d x d` matrix `Omega = [0 I; -I 0]` in element type `T` (`d` even).
"""
symplectic_form(::Type{T}, d::Integer) where {T} = _omega(T, Val(Int(d)))

function _omega(::Type{T}, ::Val{D}) where {T, D}
    iseven(D) || throw(ArgumentError("symplectic form needs an even dimension, got $D"))
    N = div(D, 2)
    return SMatrix{D, D, T}(ntuple(Val(D * D)) do l
        i = (l - 1) % D + 1
        j = div(l - 1, D) + 1
        i <= N && j == i + N ? one(T) : (i > N && j == i - N ? -one(T) : zero(T))
    end)
end

"""
    det_residual(M)

`abs(det(M) - 1)`; `det` of an `SMatrix` is static and does not allocate.
"""
det_residual(M::AbstractMatrix) = abs(det(M) - one(eltype(M)))

"""
    symplectic_residual(M)

Maximum absolute entry of `M' Omega M - Omega`.
"""
function symplectic_residual(M::StaticMatrix{D, D, T}) where {D, T}
    Om = symplectic_form(T, D)
    return maximum(abs, M' * Om * M - Om)
end
