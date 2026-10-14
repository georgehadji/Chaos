"""
    safe_norm(v) -> T

Overflow and underflow safe Euclidean norm: scale by the largest absolute entry, then take the square root of the sum
of squares. Zero for a zero vector, NaN or Inf propagate. The one norm used by every estimator (spec L1 Algorithm 1, AC8).
"""
function safe_norm(v::StaticVector{n, T}) where {n, T}
    m = maximum(abs, v)
    (isfinite(m) && m > 0) || return m
    return m * sqrt(sum(abs2, v ./ m))
end
