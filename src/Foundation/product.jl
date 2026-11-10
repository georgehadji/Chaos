"""
    fused_product(A, B) -> SMatrix

Matrix product with a fixed rounding (Amendment D115). For `SMatrix` factors with the same element type `Float32` or `Float64`,
entry `(i, j)` is `A[i, 1] * B[1, j]` followed by `acc = fma(A[i, q], B[q, j], acc)` for `q = 2..p` in ascending order: the order
of the StaticArrays product, but with every fusion explicit, so the result does not depend on the number of columns of `B`, on the
call site or on the CPU (`muladd` lets the compiler decide). Any other element type gets `A * B`. Not exported.
"""
@inline function fused_product(A::SMatrix{M, P, T}, B::SMatrix{P, N, T}) where {M, P, N, T <: Union{Float32, Float64}}
    C = MMatrix{M, N, T}(undef)
    @inbounds for j in 1:N, i in 1:M
        acc = A[i, 1] * B[1, j]
        for q in 2:P
            acc = fma(A[i, q], B[q, j], acc)
        end
        C[i, j] = acc
    end
    return SMatrix(C)
end
@inline fused_product(A, B) = A * B
