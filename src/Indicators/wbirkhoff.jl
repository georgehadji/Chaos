# Weighted Birkhoff averages (spec weighted-birkhoff: WB14, WB15, Amendment D84): the observable kinds, the pure function on a stored series, the pair
# index and the shared update. The weight (`_wb`), the compensated sums and the minimum window are those of naff-fma; the reducer in
# `wbirkhoff_reducer.jl` and the function below share one update (`_wb_sample`), so the stream result is bitwise the function on the sampled series.

const WB_MIN_WINDOW = NF_E2_MIN_WINDOW
const WB_FLOOR_MULTIPLIER = 32

"""
    CosObservable(i, k = 1)

Observable `cos(k x[i])` of the orbit point for [`WeightedBirkhoff`](@ref), with an integer frequency `k`: exact on a wrapped coordinate, whose circle is
the rounded `2 pi` (standard-map W-E6). The caller guarantees smoothness; the module cannot check it.
"""
struct CosObservable
    i::Int
    k::Int
end
CosObservable(i::Integer) = CosObservable(i, 1)

"""
    SinObservable(i, k = 1)

Observable `sin(k x[i])` of the orbit point, the companion of [`CosObservable`](@ref).
"""
struct SinObservable
    i::Int
    k::Int
end
SinObservable(i::Integer) = SinObservable(i, 1)

"""
    LiftedObservable(iq, im)

Observable `x[iq] + 2 pi_T x[im]`: the coordinate rebuilt from its winding count with `2 pi` in the element type, bitwise the corresponding entry of
`lift` (D28). For the wrapped standard map `LiftedObservable(2, 4)` is the lifted `p`. It is not periodic: its average is meaningful only for a
bounded lifted coordinate.
"""
struct LiftedObservable
    iq::Int
    im::Int
end

@inline (o::CosObservable)(x) = cos(o.k * x[o.i])
@inline (o::SinObservable)(x) = sin(o.k * x[o.i])
@inline function (o::LiftedObservable)(x)
    tp = 2 * Systems._pi(x[o.iq])
    return x[o.iq] + tp * x[o.im]
end

"""
    observable_name(o) -> String

The identifier of an observable recorded in the provenance: kind, index and frequency for the three kinds, `"function"` for any other callable.
"""
observable_name(o::CosObservable) = "cos($(o.k) x[$(o.i)])"
observable_name(o::SinObservable) = "sin($(o.k) x[$(o.i)])"
observable_name(o::LiftedObservable) = "lifted(x[$(o.iq)], x[$(o.im)])"
observable_name(::Any) = "function"

# Window argument of the constructor and the function: an integer or a nonempty tuple of integers, each at least the minimum window.
function _wb_windows(window)
    ws = window isa Integer ? (Int(window),) : Tuple(Int.(window))
    (!isempty(ws) && all(>=(WB_MIN_WINDOW), ws)) || throw(ArgumentError("window must be at least $WB_MIN_WINDOW, got $window"))
    return ws
end

function _wb_check_order(order)
    order >= 1 || throw(ArgumentError("birkhoff_order must be at least 1, got $order"))
    return Int(order)
end

# The sums of one (observable, horizon): the two window sums and the largest |f| of the 2N samples (the samples of weight zero count for S).
struct WBCell{T}
    A1::CompensatedSum{T}
    A2::CompensatedSum{T}
    S::T
end
WBCell{T}() where {T} = WBCell{T}(CompensatedSum{T}(), CompensatedSum{T}(), zero(T))

# One sample `f` of window `i` (0 or 1) with the weight `w`: A_i += w f (Neumaier), S = max(S, |f|). A weight that is exactly 0 adds nothing.
@inline function _wb_sample(c::WBCell{T}, i::Int, w::T, f::T) where {T}
    S = max(c.S, abs(f))
    iszero(w) && return WBCell{T}(c.A1, c.A2, S)
    return i == 0 ? WBCell{T}(add(c.A1, w * f), c.A2, S) : WBCell{T}(c.A1, add(c.A2, w * f), S)
end

# The unnormalised weights of one horizon (the mirrored-index weight of naff-fma) or the given vector used as it is.
_wb_weights(::Nothing, N::Int, p::Int, ::Type{T}) where {T} = T[_wb(jj, N, p, T) for jj in 0:(N - 1)]
_wb_weights(weights, N::Int, p::Int, ::Type{T}) where {T} = T[T(w) for w in weights]

function _wb_total(w::AbstractVector{T}) where {T}
    B = CompensatedSum{T}()
    for x in w
        iszero(x) || (B = add(B, x))
    end
    return total(B)
end

"""
    weighted_birkhoff(series, window; order = 1, weights = nothing) -> (avg1, avg2, S, B)

Weighted Birkhoff averages (WB14) of a stored series of any float type (Float64, Double64, BigFloat at the ambient precision). `series` is a vector (one
observable) or an `n x K` matrix (a column per observable); `window` is `N` or a tuple of `N` (horizons, all starting at sample 0: the pair `[0, N)`,
`[N, 2N)` of the first `2N` samples, later samples ignored). `avg1`, `avg2`, `S` are `K x H` matrices of the element type, `B` a vector of `H` weight
totals: `avg = A / B` with `A = sum w_jj f`, `B = sum w_jj` (Neumaier sums over the unnormalised weights, window 2 restarting the weight index at 0, one
`B` for both windows) and `S` the largest `abs(f)` of the `2N` samples. The default weight is `birkhoff_weight` of exponent `order` at the mirrored
index (`w_0 = 0`); `weights` is a vector of `N` weights used as given (one window only; equal, `t (1 - t)`, `sin^2` for the recorded lines).
`ArgumentError` naming `window` (an `N` below 16), `birkhoff_order`, `series` (fewer than `2 max(window)` samples) or `weights`.
"""
function weighted_birkhoff(series::AbstractMatrix{T}, window; order::Integer = 1, weights = nothing) where {T <: AbstractFloat}
    ws = _wb_windows(window)
    p = _wb_check_order(order)
    n, K = size(series)
    n >= 2 * maximum(ws) || throw(ArgumentError("series: $n samples given, 2 window = $(2 * maximum(ws)) needed"))
    weights === nothing || (length(ws) == 1 && length(weights) == ws[1]) ||
        throw(ArgumentError("weights: a vector of length window = $(ws[1]) is needed, with a single window"))
    H = length(ws)
    avg1, avg2, S, B = (Matrix{T}(undef, K, H), Matrix{T}(undef, K, H), Matrix{T}(undef, K, H), Vector{T}(undef, H))
    for (h, N) in enumerate(ws)
        w = _wb_weights(weights, N, p, T)
        B[h] = _wb_total(w)
        for k in 1:K
            c = WBCell{T}()
            for j in 0:(2N - 1)
                i = j >= N ? 1 : 0
                c = _wb_sample(c, i, w[j - i * N + 1], series[j + 1, k])
            end
            avg1[k, h], avg2[k, h], S[k, h] = total(c.A1) / B[h], total(c.A2) / B[h], c.S
        end
    end
    return (avg1 = avg1, avg2 = avg2, S = S, B = B)
end
weighted_birkhoff(series::AbstractVector{<:AbstractFloat}, window; kw...) = weighted_birkhoff(reshape(series, :, 1), window; kw...)

"""
    birkhoff_index(avg1, avg2, S; floor_multiplier = 32) -> (delta, absdig, reldig, dig, below_floor, zero)

The index of an average pair (WB15), of the type of its arguments and broadcastable over the `K x H` results: `delta = avg1 - avg2` (never clipped);
`phi = floor_multiplier effective_eps(T) S`; `absdig = -log10 max(abs(delta), phi)`; with `m = (abs(avg1) + abs(avg2)) / 2`,
`reldig = -log10(max(abs(delta), phi) / m)` where `m > phi` and `NaN` elsewhere; `dig = max(absdig, reldig)` where `reldig` is defined, else `absdig`;
`below_floor = abs(delta) <= phi`. For `S = 0` the digits are `NaN` and `zero` (the observable is identically 0) is set, together with `below_floor`.
`absdig` is never above `-log10 phi` and `reldig` never above `-log10(phi / m)`; `dig` can reach `-log10(32 eps)`. `ArgumentError` naming
`floor_multiplier` unless it is positive.
"""
function birkhoff_index(avg1::T, avg2::T, S::T; floor_multiplier::Real = WB_FLOOR_MULTIPLIER) where {T <: AbstractFloat}
    floor_multiplier > 0 || throw(ArgumentError("floor_multiplier must be positive, got $floor_multiplier"))
    delta = avg1 - avg2
    phi = T(floor_multiplier) * effective_eps(T) * S
    a = abs(delta)
    below = a <= phi
    iszero(S) && return (delta = delta, absdig = T(NaN), reldig = T(NaN), dig = T(NaN), below_floor = below, zero = true)
    e = max(a, phi)
    m = (abs(avg1) + abs(avg2)) / 2
    absdig = -log10(e)
    reldig = m > phi ? -log10(e / m) : T(NaN)
    return (delta = delta, absdig = absdig, reldig = reldig, dig = isnan(reldig) ? absdig : max(absdig, reldig), below_floor = below, zero = false)
end
