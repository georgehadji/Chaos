# Pure functions of the Smaller Alignment Index (spec sali, A2 to A4, A14).

"""
    sali_floor(T) -> T

Measurement floor `B_meas(T) = 4 effective_eps(T)` of the direct form of SALI (spec A14): a computed SALI below it is unresolved.
"""
sali_floor(::Type{T}) where {T} = 4 * effective_eps(T)

# Euclidean norm scaled by the largest absolute entry, so it neither overflows nor underflows; zero for a zero vector,
# NaN and Inf propagate (the same rule as `Lyapunov.safe_norm`, which this layer may not import).
function scaled_norm(v::AbstractVector)
    m = maximum(abs, v)
    (isfinite(m) && m > 0) || return m
    return m * sqrt(sum(abs2, v ./ m))
end

# d_- and d_+ of two vectors that are already unit vectors (A2), each evaluated scaled.
@inline unit_indices(u, v) = (scaled_norm(u - v), scaled_norm(u + v))

"""
    alignment_indices(a, b) -> (d_minus, d_plus)

Parallel and antiparallel alignment indices (spec A2) of two nonzero vectors of any dimension: each input is divided by
its scaled norm, then `d_minus = norm(u_a - u_b)` and `d_plus = norm(u_a + u_b)` are taken scaled. Always
`d_minus^2 + d_plus^2 = 4` up to rounding.
"""
function alignment_indices(a::AbstractVector, b::AbstractVector)
    return unit_indices(a / scaled_norm(a), b / scaled_norm(b))
end

"""
    sali_pair(a, b) -> T

`min(alignment_indices(a, b)...)`, the SALI of two nonzero vectors (spec A2): 0 for parallel or antiparallel vectors, `sqrt(2)`
for orthogonal ones. Bitwise unchanged by a sign change of either vector and by exchanging them.
"""
sali_pair(a::AbstractVector, b::AbstractVector) = min(alignment_indices(a, b)...)

"""
    sali_from_log(L) -> T

SALI from the log form `L = ln GALI_2` (spec A4): `sqrt(2) g / sqrt(1 + sqrt(1 - g^2))` with `g = min(exp(L), 1)`. No
cancellation for small `g`, and `exp` underflows to 0 for very negative `L`. Valid for SALI well above `eps`; below that
use `L` itself.
"""
function sali_from_log(L::AbstractFloat)
    g = min(exp(L), one(L))
    return sqrt(2 * one(L)) * g / sqrt(1 + sqrt(1 - g^2))
end

# SALI(0) of a d x 2 frame: 0 when a column has a zero or non-finite norm.
function frame_sali(W::AbstractMatrix{T}) where {T}
    a1, a2 = scaled_norm(W[:, 1]), scaled_norm(W[:, 2])
    (isfinite(a1) && a1 > 0 && isfinite(a2) && a2 > 0) || return zero(T)
    return min(unit_indices(W[:, 1] / a1, W[:, 2] / a2)...)
end

"""
    check_frame(SALI(), W0; orbit = 0) -> SALI(0)

Frame check of a supplied initial deviation frame (spec Algorithm 1, Amendment D67): returns `SALI(0)` of the two columns
if it is at least `10 sali_floor(T)`, else throws `NumericalHealthError("initial frame", ...)` at step 0. Dependent
columns would give an absolute zero (spec A16), which reads as chaos. The caller runs it on every supplied frame before the
stream starts, so a rejected frame never reaches the map; the seeded frame of `propagate` passes by construction.
"""
function check_frame(::SALI, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) == 2 || throw(ArgumentError("SALI needs a frame with 2 columns, got $(size(W0, 2))"))
    s0 = frame_sali(W0)
    threshold = 10 * sali_floor(T)
    s0 >= threshold || throw(NumericalHealthError("initial frame", Float64(s0), Float64(threshold), 0, orbit, :cpu, 0))
    return s0
end
