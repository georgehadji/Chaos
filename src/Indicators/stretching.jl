# Stretching numbers reducer (spec stretching-numbers, SN1 to SN17): cumulative spectra of the per-handoff log growth of two deviation
# vectors and the spectral distance D between them. The counts live in the per-orbit workspace of the Tangent layer (P1 (c)).

const SN_CONVENTION = "cumulative spectra, D at every interval, columns renormalised separately, no orthogonalisation, overflow cells of width h"

"""
    StretchingNumbers(; range, bins = 1000, interval = 100, min_per_cell = 100, replace = true)

Stretching-number spectra reducer (Voglis, Darriba, Maffione; catalog I10 and L13) for `propagate` with `k = 2` (any other
`k` throws `ArgumentError`). At each handoff `n = K tau` it reads the two arriving, unnormalised columns, forms `a_j = log(alpha_j)`,
`alpha_j` the scaled norm of column `j`, adds `a_j` to a compensated sum (the same operations as `Benettin`, [`SALI`](@ref),
[`FLI`](@ref) and [`MEGNO`](@ref), so the sums are bitwise theirs) and `a_j^2` to a second one, and counts the stretching number
`z_j = a_j / tau` in a cell of a cumulative histogram of `bins` cells of equal width `h` over `range = (s_lo, s_hi)`, plus two
overflow cells below and above (never clipped or dropped). The counts are `UInt32` entries of the per-orbit workspace
([`reducer_workspace`](@ref), `2 (bins + 2)` entries, column `j` in entries `(j-1)(bins+2)+1 ..`, cell `c` at offset `c + 1`),
so the isbits state is small and the hot loop allocates nothing.

`range` is required: the reducer cannot see the map ([`stretching_range`](@ref) gives the derived bound for the standard maps
and `LinearSymplecticMap`). `interval` (map steps) and `min_per_cell` are recorded and validated parameters: D is a function of the
counts and is formed by [`reducer_final`](@ref) with the workspace, at whatever handoff the reader asks (the stream has no emission
channel), and the under-sampling flag is `handoffs < min_per_cell * bins`. [`check_run`](@ref) checks `interval` against `tau` and
the run length before the stream. `replace` is the chain position of the one-replacer rule (as `MEGNO`): `true` alone or last in the
chain, `false` before a later replacer. No orthogonalisation, ever.

!!! warning
    A wrong chain position gives wrong numbers without an error (see [`MEGNO`](@ref)): exactly one replacer per chain, placed last.
    At most one reducer of a chain may have a workspace, so `StretchingNumbers` cannot share a chain with another such reducer.

`reducer_final(r, state)` gives the scalar part `(n, handoffs, sum, comp, sum2, comp2, mean, std)` (the `SVector{2}` log sums with
correction terms, the sums of `a^2`, `mean = total / n`, the root of `total2 / (tau n) - mean^2` clamped at 0); with the workspace
it adds `D`, `undersampled`, `counts` (a `(bins + 2) x 2` copy, row `c + 1` is cell `c`) and the overflow counts `below`, `above`.
A zero or non-finite column norm is recorded as a health failure ("deviation matrix", the column) and leaves sums, counts and `W`
unchanged.
"""
struct StretchingNumbers{R <: AbstractFloat}
    lo::R
    hi::R
    bins::Int
    interval::Int
    min_per_cell::Int
    replace::Bool
end

function StretchingNumbers(; range, bins::Integer = 1000, interval::Integer = 100, min_per_cell::Integer = 100, replace::Bool = true)
    lo, hi = range
    R = float(promote_type(typeof(lo), typeof(hi)))
    (isfinite(lo) && isfinite(hi) && lo < hi) || throw(ArgumentError("range must be finite with s_lo < s_hi, got ($lo, $hi)"))
    bins >= 1 || throw(ArgumentError("bins must be at least 1, got $bins"))
    interval >= 1 || throw(ArgumentError("interval must be at least 1 map step, got $interval"))
    min_per_cell >= 0 || throw(ArgumentError("min_per_cell must be nonnegative, got $min_per_cell"))
    return StretchingNumbers{R}(R(lo), R(hi), Int(bins), Int(interval), Int(min_per_cell), replace)
end

"""
    check_run(r::StretchingNumbers, tau, N) -> nothing

The caller's check before the stream (the protocol gives `reducer_init` neither `tau` nor `N`): `ArgumentError` naming `interval`
if it is not a multiple of `tau`, naming `handoffs` if `N / tau` exceeds `10^9` (4-byte counts, exact integer sums of squares).
"""
function check_run(r::StretchingNumbers, tau::Integer, N::Integer)
    r.interval % tau == 0 || throw(ArgumentError("interval $(r.interval) must be a multiple of tau = $tau"))
    N ÷ tau <= 10^9 || throw(ArgumentError("a run of $(N ÷ tau) handoffs exceeds the limit of 10^9 handoffs"))
    return nothing
end

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`; `n` is the last successful map step,
# `K` the number of successful handoffs; s1, s2 the log sums of the columns, q1, q2 those of the squares.
struct StretchState{T}
    s1::CompensatedSum{T}
    s2::CompensatedSum{T}
    q1::CompensatedSum{T}
    q2::CompensatedSum{T}
    value::T
    n::Int
    K::Int
    step::Int
    code::Int8
    column::Int8
end

function reducer_init(::StretchingNumbers{R}, ::Type{T}, ::Val, ::Val{k}) where {R, T, k}
    k == 2 || throw(ArgumentError("StretchingNumbers needs k = 2, got k = $k"))
    promote_type(T, R) === T || throw(ArgumentError("range of type $R would be narrowed to $T"))
    z = CompensatedSum{T}()
    return StretchState{T}(z, z, z, z, zero(T), 0, 0, 0, Int8(0), Int8(0))
end

function reducer_workspace(r::StretchingNumbers, ::Type, ::Val, ::Val{k}) where {k}
    k == 2 || throw(ArgumentError("StretchingNumbers needs k = 2, got k = $k"))
    return 2 * (r.bins + 2)
end

# floor to Int without `floor(Int, ::Double64)`, which allocates 304 bytes: the high part is an integer or at least an ulp from the
# next one, and the low part is at most half an ulp, so it only matters when the high part is an exact integer.
@inline _floor_int(x::AbstractFloat) = floor(Int, x)
@inline function _floor_int(x::Double64)
    f = floor(Int, x.hi)
    return (x.hi == f && x.lo < 0) ? f - 1 : f
end

"""
    spectrum_cell(z, lo, hi, B) -> Int

The cell rule of the spectrum (SN4): `0` below `lo`, `B + 1` above `hi`, else `min(floor((z - lo) / h) + 1, B)` with
`h = (hi - lo) / B`, so `z == lo` is cell 1 and `z == hi` is cell `B` (the last cell is closed). Computed in the type of the
arguments; the reducer calls this function.
"""
@inline function spectrum_cell(z, lo, hi, B::Integer)
    z < lo && return 0
    z > hi && return Int(B) + 1
    return min(_floor_int((z - lo) / ((hi - lo) / B)) + 1, Int(B))
end

@inline function _count!(ws, offset::Int, z, lo, hi, B::Int)
    i = offset + spectrum_cell(z, lo, hi, B) + 1
    ws[i] += 0x00000001
    return nothing
end

function reducer_step(r::StretchingNumbers, st::StretchState{T}, n::Int, x::SVector, W::SMatrix, ws::AbstractVector{UInt32}) where {T}
    st.code == 0 || return st, W
    c1, c2 = W[:, 1], W[:, 2]
    a1, a2 = safe_norm(c1), safe_norm(c2)
    ok1, ok2 = isfinite(a1) && a1 > 0, isfinite(a2) && a2 > 0
    if !(ok1 && ok2)
        column, value = ok1 ? (Int8(2), a2) : (Int8(1), a1)
        return StretchState{T}(st.s1, st.s2, st.q1, st.q2, value, st.n, st.K, n, Int8(1), column), W
    end
    K = st.K + 1
    tau, lo, hi, B = T(n ÷ K), T(r.lo), T(r.hi), r.bins
    l1, l2 = log(a1), log(a2)
    _count!(ws, 0, l1 / tau, lo, hi, B)
    _count!(ws, B + 2, l2 / tau, lo, hi, B)
    st1 = StretchState{T}(add(st.s1, l1), add(st.s2, l2), add(st.q1, l1 * l1), add(st.q2, l2 * l2), st.value, n, K, st.step, st.code, st.column)
    return st1, r.replace ? hcat(c1 / a1, c2 / a2) : W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::StretchingNumbers, st::StretchState, n::Int, x::SVector, W::SMatrix, ws::AbstractVector{UInt32})
    st1, W1 = reducer_step(r, st, n, x, W, ws)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::StretchingNumbers, st::StretchState) = st.code != 0

function reducer_error(::StretchingNumbers, st::StretchState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, Int(st.column))
end

# Mean and standard deviation of z of one column from its two compensated sums (SN2, SN7); tau = n / K is a protocol invariant.
function stretch_moments(s::CompensatedSum{T}, q::CompensatedSum{T}, n::Int, K::Int) where {T}
    K == 0 && return zero(T), zero(T)
    m = total(s) / T(n)
    return m, sqrt(max(total(q) / (T(n ÷ K) * T(n)) - m * m, zero(T)))
end

function reducer_final(::StretchingNumbers, st::StretchState{T}) where {T}
    m1, d1 = stretch_moments(st.s1, st.q1, st.n, st.K)
    m2, d2 = stretch_moments(st.s2, st.q2, st.n, st.K)
    return (n = st.n, handoffs = st.K, sum = SVector(st.s1.s, st.s2.s), comp = SVector(st.s1.c, st.s2.c),
            sum2 = SVector(st.q1.s, st.q2.s), comp2 = SVector(st.q1.c, st.q2.c), mean = SVector(m1, m2), std = SVector(d1, d2))
end

function reducer_final(r::StretchingNumbers, st::StretchState{T}, ws::AbstractVector{UInt32}) where {T}
    B = r.bins
    counts = Matrix{UInt32}(undef, B + 2, 2)
    for j in 1:2, c in 1:(B + 2)
        counts[c, j] = ws[(j - 1) * (B + 2) + c]
    end
    h = (T(r.hi) - T(r.lo)) / B
    return merge(reducer_final(r, st), (D = spectral_distance(view(counts, :, 1), view(counts, :, 2), h), undersampled = st.K < r.min_per_cell * B,
                 counts = counts, below = SVector{2, Int}(counts[1, 1], counts[1, 2]), above = SVector{2, Int}(counts[B + 2, 1], counts[B + 2, 2])))
end

function reducer_provenance(r::StretchingNumbers)
    return (indicator = "stretching", convention = SN_CONVENTION, spectrum_bins = r.bins, spectrum_lo = r.lo, spectrum_hi = r.hi,
            spectrum_cell_width = (r.hi - r.lo) / r.bins, spectrum_interval = r.interval, spectrum_min_per_cell = r.min_per_cell, replace = r.replace)
end

# An integer converted into T without the Float64 detour of `Double64(::Int)` above 2^53 (the sum of squared count differences).
_int_in(::Type{T}, s::Int) where {T} = s < 2^53 ? T(s) : exact(T, s)

"""
    spectral_distance(N1, N2, h) -> typeof(h)

Spectral distance D of two count vectors (SN5): `sqrt(S) / (k sqrt(h))`, `S` the exact `Int` sum of squared differences and `k` the
common total, so `0 <= D <= sqrt(2 / h)`, `D(N, N) = 0` and D is bitwise symmetric. The overflow cells are part of the vectors
and have the width `h` like the others. `ArgumentError` for unequal totals, `DimensionMismatch` for unequal lengths, `0` for
total 0. Counts of any integer type; 64-bit sums are exact up to `10^9` handoffs.
"""
function spectral_distance(N1::AbstractVector{<:Integer}, N2::AbstractVector{<:Integer}, h::AbstractFloat)
    length(N1) == length(N2) || throw(DimensionMismatch("count vectors of lengths $(length(N1)) and $(length(N2))"))
    S, k1, k2 = 0, 0, 0
    for c in eachindex(N1, N2)
        a, b = Int(N1[c]), Int(N2[c])
        k1, k2, S = k1 + a, k2 + b, S + (a - b)^2
    end
    k1 == k2 || throw(ArgumentError("count vectors must have equal totals, got $k1 and $k2"))
    k1 == 0 && return zero(h)
    return sqrt(_int_in(typeof(h), S)) / (k1 * sqrt(h))
end

"""
    spectrum_density(counts, h) -> Vector

Normalised spectrum (SN3) of a count column of length `B + 2` (cells `0..B + 1`): `N_c / (k h)` for `c = 1..B`, `k` the column
total including the two overflow counts.
"""
function spectrum_density(counts::AbstractVector{<:Integer}, h::AbstractFloat)
    k = sum(Int, counts)
    return [k == 0 ? zero(h) : typeof(h)(Int(counts[c + 1])) / (k * h) for c in 1:(length(counts) - 2)]
end

# BigFloat value of a real number, exact for a Double64.
_bigval(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
_bigval(x::AbstractFloat) = BigFloat(x)

# Symmetric derived bound (SN12) padded outward by 8 eps relative, computed at 256 bits and rounded once into T: the stretching numbers
# of Y5 approach ln(mu) to 1e-30, so an unpadded bound would put about half of them into an overflow cell by rounding.
function _padded_range(::Type{T}, L::BigFloat) where {T}
    hi = T(L * (1 + 8 * _bigval(effective_eps(T))))
    return (-hi, hi)
end

"""
    stretching_range(sys) -> (lo, hi)

Derived symmetric cell range `(-hi, hi)` of the spectra (SN12) for `StandardMap`, `WrappedStandardMap` (`hi` just above
`L_sigma = ln sigma(K)`, the one-step singular value bound of L1 E14) and a 2 x 2 `LinearSymplecticMap` (`ln` of the largest
singular value); `hi` is `L (1 + 8 eps)` rounded, in the element type of the map. Every stretching number of the map lies inside
for every `tau`. `ArgumentError` for a map without stretching (largest singular value within `8 eps` of 1 or below) and for a
`LinearSymplecticMap` that is not 2 x 2 (pass the range by hand).
"""
function stretching_range(sys::Union{StandardMap{T}, WrappedStandardMap{T}}) where {T}
    return setprecision(BigFloat, 256) do
        K = _bigval(sys.K)
        S2 = 2 * K^2 + 2 * abs(K) + 5
        _padded_range(T, log((sqrt(S2) + sqrt(S2 - 4)) / 2))
    end
end

function stretching_range(sys::LinearSymplecticMap{T, D}) where {T, D}
    D == 2 || throw(ArgumentError("stretching_range needs a 2 x 2 LinearSymplecticMap, got $D x $D; pass the range explicitly"))
    return setprecision(BigFloat, 256) do
        a, b, c, d = (_bigval(sys.A[i]) for i in 1:4)
        F, det = a^2 + b^2 + c^2 + d^2, a * d - b * c
        L = log(sqrt((F + sqrt(max(F^2 - 4 * det^2, zero(BigFloat)))) / 2))
        L > 8 * _bigval(effective_eps(T)) || throw(ArgumentError("the map has no stretching (largest singular value $(exp(L)) is not above 1)"))
        _padded_range(T, L)
    end
end

"""
    check_frame(r::StretchingNumbers, W0; orbit = 0) -> defect

Frame check of a supplied initial deviation frame (sali Amendment D67 form, fli D72 item 7), in this order: both columns unit,
`abs(norm(w_i) - 1) <= 8 effective_eps(T)` (column 1 first; a non-finite norm fails), then not parallel,
`min(norm(w1 - w2), norm(w1 + w2)) > 8 eps` (a parallel pair gives identical spectra and D = 0 at every handoff). Returns the larger
norm defect, else throws `NumericalHealthError("initial frame", ...)` at step 0 naming the column (column 2 for parallel
columns, with the smaller norm as value). `ArgumentError` for a frame without 2 columns. The caller runs it on every supplied
frame before the stream; the seeded frame passes by construction.
"""
function check_frame(::StretchingNumbers, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) == 2 || throw(ArgumentError("StretchingNumbers needs a frame with 2 columns, got $(size(W0, 2))"))
    threshold = 8 * effective_eps(T)
    fail(value, col) = throw(NumericalHealthError("initial frame", Float64(value), Float64(threshold), 0, orbit, :cpu, col))
    defect = zero(T)
    for i in 1:2
        d = abs(safe_norm(W0[:, i]) - 1)
        (isfinite(d) && d <= threshold) || fail(d, i)
        defect = max(defect, d)
    end
    w1, w2 = W0[:, 1], W0[:, 2]
    para = min(safe_norm(w1 - w2), safe_norm(w1 + w2))
    para > threshold || fail(para, 2)
    return defect
end
