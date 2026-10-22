# MEGNO reducer (spec megno, Algorithm 1 to 4, M1 to M13): weighted sums of the log stretches of two deviation vectors.

"""
    MEGNO(; replace = true)

Mean Exponential Growth factor of Nearby Orbits reducer for `propagate` with `k = 2` (any other `k` throws `ArgumentError`).
At each handoff `n = K tau` it reads the two arriving, unnormalised columns, forms `a_i = log(alpha_i)`, `alpha_i` the scaled
norm of column `i`, and keeps per column five compensated sums (O(1) state per orbit): the log sum `s_i` (the same operations
as [`SALI`](@ref), [`FLI`](@ref) and `Benettin`, so it is bitwise equal to theirs), `R1 = sum n a`, `R2 = sum n^2 a`, and the
running sums `Q = sum Y`, `Q2 = sum 3 R2`. The delivered series are
`Y = 2 R1 / n` (the (1,-1) MEGNO, weights at the end of the window), `Ybar = Q / K` (its mean over the handoffs `1..K`) and
`Ybar20 = tau Q2 / n^3` (the (2,0) form of Maffione2011a, regular level 1/2), per column, column 1 the primary. The weights
are formed in the element type, never as 64-bit integers (`n^3` leaves the integer range at about 2.1e6 handoffs).

`replace` is the chain position of the one-replacer rule: with `replace = true` (alone in a chain, or last) `W` is replaced by
its column-normalised matrix; with `replace = false` the columns pass on unchanged, for a chain in which a later reducer
(`SALI`) replaces them. The reducer cannot see the chain, so the flag is the caller's.

!!! warning
    A wrong position gives wrong numbers without an error: `MEGNO(replace = false)` alone sums unnormalised columns, and a
    replacing `MEGNO()` before `SALI` leaves SALI with unit columns. Use exactly one replacer per chain, placed last (D71).
    The run record does not show `replace` or the chain; the caller (Experiments) records the chain (D72, fli review R1, R2).

`reducer_final` gives `(n, handoffs, Y, Ybar, Ybar20, sum, comp, R1, R2, Q, Q2)`: `n` the last handoff (map step),
`handoffs = K`, the three series and the log sums with compensation terms as `SVector{2}` (column 1 first), and the totals of
the weighted sums. A zero or non-finite column norm is recorded as a health failure ("deviation matrix", the column) and
leaves sums and `W` unchanged. The values are never clipped, capped or replaced; `Y` is unbounded and may be negative.
"""
struct MEGNO
    replace::Bool
end
MEGNO(; replace::Bool = true) = MEGNO(replace)

# The five compensated sums of one column.
struct MegnoSums{T}
    s::CompensatedSum{T}
    R1::CompensatedSum{T}
    R2::CompensatedSum{T}
    Q::CompensatedSum{T}
    Q2::CompensatedSum{T}
end
MegnoSums{T}() where {T} = MegnoSums{T}(CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}())

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`; `n` is the last successful map step,
# `K` the number of successful handoffs.
struct MEGNOState{T}
    c1::MegnoSums{T}
    c2::MegnoSums{T}
    value::T
    n::Int
    K::Int
    step::Int
    code::Int8
    column::Int8
end

function reducer_init(::MEGNO, ::Type{T}, ::Val, ::Val{k}) where {T, k}
    k == 2 || throw(ArgumentError("MEGNO needs k = 2, got k = $k"))
    return MEGNOState{T}(MegnoSums{T}(), MegnoSums{T}(), zero(T), 0, 0, 0, Int8(0), Int8(0))
end

# One column at map step `n` with log stretch `a` (spec M5, weights in the element type): R1, R2, then Y into Q and 3 R2 into Q2.
function megno_add(m::MegnoSums{T}, a::T, n::Int) where {T}
    tn = T(n)
    R1, R2 = add(m.R1, tn * a), add(m.R2, tn * tn * a)
    return MegnoSums{T}(add(m.s, a), R1, R2, add(m.Q, 2 * total(R1) / tn), add(m.Q2, 3 * total(R2)))
end

function reducer_step(r::MEGNO, st::MEGNOState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    c1, c2 = W[:, 1], W[:, 2]
    a1, a2 = safe_norm(c1), safe_norm(c2)
    ok1, ok2 = isfinite(a1) && a1 > 0, isfinite(a2) && a2 > 0
    if !(ok1 && ok2)
        column, value = ok1 ? (Int8(2), a2) : (Int8(1), a1)
        return MEGNOState{T}(st.c1, st.c2, value, st.n, st.K, n, Int8(1), column), W
    end
    st1 = MEGNOState{T}(megno_add(st.c1, log(a1), n), megno_add(st.c2, log(a2), n), st.value, n, st.K + 1, st.step, st.code, st.column)
    return st1, r.replace ? hcat(c1 / a1, c2 / a2) : W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::MEGNO, st::MEGNOState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::MEGNO, st::MEGNOState) = st.code != 0

function reducer_error(::MEGNO, st::MEGNOState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, Int(st.column))
end

# (Y, Ybar, Ybar20) of one column at the last handoff; all zero before the first (spec M5, tau = n / K).
function megno_series(m::MegnoSums{T}, n::Int, K::Int) where {T}
    K == 0 && return zero(T), zero(T), zero(T)
    tn = T(n)
    # n = K tau is a protocol invariant of the stream (handoffs at multiples of tau), so n ÷ K is tau exactly (review R1)
    return 2 * total(m.R1) / tn, total(m.Q) / T(K), T(n ÷ K) * total(m.Q2) / tn^3
end

function reducer_final(::MEGNO, st::MEGNOState{T}) where {T}
    y1, b1, e1 = megno_series(st.c1, st.n, st.K)
    y2, b2, e2 = megno_series(st.c2, st.n, st.K)
    pair(f) = SVector(f(st.c1), f(st.c2))
    return (n = st.n, handoffs = st.K, Y = SVector(y1, y2), Ybar = SVector(b1, b2), Ybar20 = SVector(e1, e2),
            sum = pair(m -> m.s.s), comp = pair(m -> m.s.c), R1 = pair(m -> total(m.R1)), R2 = pair(m -> total(m.R2)),
            Q = pair(m -> total(m.Q)), Q2 = pair(m -> total(m.Q2)))
end

reducer_provenance(::MEGNO) = (indicator = "megno", convention = "weights at the end of the window, mean over the handoffs 1..K")

"""
    check_frame(MEGNO(), W0; orbit = 0) -> defect

Frame check of a supplied initial deviation frame (spec Algorithm 1, Amendment D72 item 7 form): both columns must be unit
vectors, `abs(norm(w_i) - 1) <= 8 effective_eps(T)` (the convention `s(0) = 0`). Returns the larger defect, else throws
`NumericalHealthError("initial frame", ...)` at step 0 naming the column. Independence and orthogonality are not checked. The
caller runs it on every supplied frame before the stream; the seeded frame passes by construction.
"""
function check_frame(::MEGNO, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) == 2 || throw(ArgumentError("MEGNO needs a frame with 2 columns, got $(size(W0, 2))"))
    threshold = 8 * effective_eps(T)
    defect = zero(T)
    for i in 1:2
        d = abs(safe_norm(W0[:, i]) - 1)
        (isfinite(d) && d <= threshold) || throw(NumericalHealthError("initial frame", Float64(d), Float64(threshold), 0, orbit, :cpu, i))
        defect = max(defect, d)
    end
    return defect
end
