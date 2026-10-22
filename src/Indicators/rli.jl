# RLI reducer (spec rli, RL1 to RL13): the Relative Lyapunov Indicator of an orbit and its shadow, run as one orbit of the product map.

"""
    RLI(; replace = true)

Relative Lyapunov Indicator reducer for `propagate` with `k = 2` on a pair system ([`PairMap`](@ref), even dimension `d`).
Column 1 of the frame is `(xi_0, 0)` (orbit A), column 2 is `(0, xi_0)` (orbit B). At each handoff `n = K tau` it reads the
two arriving, unnormalised columns, forms `a_i = log(alpha_i)`, `alpha_i` the scaled norm of column `i`, adds it to the
compensated log sum `s_i` (the same operations as `Benettin`, so `s_i` is bitwise the single-orbit sum), and keeps
`DL = abs(S_B - S_A) / n` with `S_i = s_i + comp_i` (the correction terms included) and the compensated sum `Q` of the `DL`
series; `RLI = Q / K` is the mean over the handoffs `1..K` (O(1) state per pair, never clipped).

`replace` is the chain position of the one-replacer rule: with `replace = true` (alone in a chain, or last) `W` is replaced
by its column-normalised matrix; with `replace = false` the columns pass on unchanged, for a chain in which a later reducer
(`MEGNO`) replaces them. The reducer cannot see the chain, so the flag is the caller's.

!!! warning
    A wrong position gives wrong numbers without an error (see [`MEGNO`](@ref)): use exactly one replacer per chain, placed
    last (D71). The caller builds the pair start with [`pair_start`](@ref) and the block frame, and runs
    [`check_frame`](@ref) on it before the stream.

`reducer_final` gives `(n, handoffs, DL, RLI, sum, comp, Q)`: `n` the last handoff (map step), `handoffs = K`, `DL` and `RLI`
the series values at the last handoff, `sum` and `comp` the log sums with compensation terms as `SVector{2}` (A first), `Q`
the total of the `DL` series. A zero or non-finite column norm is recorded as a health failure ("deviation matrix", the
column) and leaves sums and `W` unchanged. At every handoff the off-block entries of the arriving frame (block B rows of
column 1, block A rows of column 2) must be exactly 0, as the block action of `PairMap` keeps them; otherwise the reducer
records a failure "pairing" (the column, the largest off-block magnitude), so a run on a non-pair stream or without the
block frame fails at the first handoff (spec Amendment D77 item 3). Mismatched `xi_0` between the blocks is not visible
here: [`check_frame`](@ref) rejects it.
"""
struct RLI
    replace::Bool
end
RLI(; replace::Bool = true) = RLI(replace)

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`, 2 = nonzero off-block magnitude `value`; `n` is the last successful map step,
# `K` the number of successful handoffs.
struct RLIState{T}
    sA::CompensatedSum{T}
    sB::CompensatedSum{T}
    Q::CompensatedSum{T}
    value::T
    n::Int
    K::Int
    step::Int
    code::Int8
    column::Int8
end

function reducer_init(::RLI, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k}
    k == 2 || throw(ArgumentError("RLI needs k = 2, got k = $k"))
    iseven(d) || throw(ArgumentError("RLI needs an even dimension, got d = $d"))
    return RLIState{T}(CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}(), zero(T), 0, 0, 0, Int8(0), Int8(0))
end

# DL at map step n from the two compensated log sums, correction terms included (spec RL4, F5); the split difference keeps
# the cancellation of nearby orbits exact where the two sums are close (D77 item 4, review R2).
rli_dl(sA::CompensatedSum{T}, sB::CompensatedSum{T}, n::Int) where {T} = abs((sB.s - sA.s) + (sB.c - sA.c)) / T(n)

# Largest off-block magnitude of each column of a `2h x 2` frame (D77 item 3); the block action keeps both exactly 0.
@inline function _offblock(W::SMatrix{D, 2, T}) where {D, T}
    h = D ÷ 2
    m1, m2 = zero(T), zero(T)
    for i in 1:h
        m1, m2 = max(m1, abs(W[h + i, 1])), max(m2, abs(W[i, 2]))
    end
    return m1, m2
end

function reducer_step(r::RLI, st::RLIState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    c1, c2 = W[:, 1], W[:, 2]
    a1, a2 = safe_norm(c1), safe_norm(c2)
    ok1, ok2 = isfinite(a1) && a1 > 0, isfinite(a2) && a2 > 0
    if !(ok1 && ok2)
        column, value = ok1 ? (Int8(2), a2) : (Int8(1), a1)
        return RLIState{T}(st.sA, st.sB, st.Q, value, st.n, st.K, n, Int8(1), column), W
    end
    m1, m2 = _offblock(W)
    if !(m1 == 0 && m2 == 0)
        column, value = m1 == 0 ? (Int8(2), m2) : (Int8(1), m1)
        return RLIState{T}(st.sA, st.sB, st.Q, value, st.n, st.K, n, Int8(2), column), W
    end
    sA, sB = add(st.sA, log(a1)), add(st.sB, log(a2))
    st1 = RLIState{T}(sA, sB, add(st.Q, rli_dl(sA, sB, n)), st.value, n, st.K + 1, st.step, st.code, st.column)
    return st1, r.replace ? hcat(c1 / a1, c2 / a2) : W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::RLI, st::RLIState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::RLI, st::RLIState) = st.code != 0

function reducer_error(::RLI, st::RLIState, orbit, backend)
    value = Float64(st.value)
    what = st.code == 2 ? "pairing" : "deviation matrix"
    return NumericalHealthError(what, value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, Int(st.column))
end

function reducer_final(::RLI, st::RLIState{T}) where {T}
    dl, rli = st.K == 0 ? (zero(T), zero(T)) : (rli_dl(st.sA, st.sB, st.n), total(st.Q) / T(st.K))
    return (n = st.n, handoffs = st.K, DL = dl, RLI = rli, sum = SVector(st.sA.s, st.sB.s), comp = SVector(st.sA.c, st.sB.c), Q = total(st.Q))
end

reducer_provenance(::RLI) = (indicator = "rli", convention = "same xi_0, shadow along xi_0, mean over the handoffs 1..K")

"""
    check_frame(RLI(), W0; orbit = 0) -> defect

Frame check of a supplied block frame (spec Algorithm 1, sali Amendment D67 form), in this order: both columns unit vectors,
`abs(norm(w_i) - 1) <= 8 effective_eps(T)` (column 1 first); column 1 supported in the first half of the rows and column 2
in the second half (zero elsewhere); the two block vectors bitwise equal. Returns the larger norm defect, else throws
`NumericalHealthError("initial frame", ...)` at step 0 naming the column (column 2 when the block vectors differ; threshold
`8 eps` for a norm defect, 0 for the support and equality checks). `ArgumentError` for a frame without 2 columns or with an
odd number of rows. The caller runs it on every supplied frame before the stream.
"""
function check_frame(::RLI, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) == 2 || throw(ArgumentError("RLI needs a frame with 2 columns, got $(size(W0, 2))"))
    D = size(W0, 1)
    iseven(D) || throw(ArgumentError("RLI needs a frame with an even number of rows, got $D"))
    h, threshold = D ÷ 2, 8 * effective_eps(T)
    fail(value, thr, col) = throw(NumericalHealthError("initial frame", Float64(value), Float64(thr), 0, orbit, :cpu, col))
    defect = zero(T)
    for i in 1:2
        d = abs(safe_norm(W0[:, i]) - 1)
        (isfinite(d) && d <= threshold) || fail(d, threshold, i)
        defect = max(defect, d)
    end
    for (col, rows) in ((1, h + 1:D), (2, 1:h))
        s = maximum(abs, W0[rows, col])
        s == 0 || fail(s, 0.0, col)
    end
    diff = maximum(abs, W0[1:h, 1] .- W0[h + 1:D, 2])
    diff == 0 || fail(diff, 0.0, 2)
    return defect
end

# Minimal image of a coordinate difference for any magnitude; exact for a difference below pi (round gives 0).
@inline _wrap_diff(u, tp) = u - tp * round(u / tp)

# Separation vector of the two halves (A; B) of a pair state: minimal image on the angle coordinates of the standard maps (the
# wrapped one has the counts as components 3, 4, not used), the plain difference otherwise (spec RL13).
@inline function _pair_delta(::Union{StandardMap, WrappedStandardMap}, z::SVector)
    h = length(z) ÷ 2
    tp = 2 * Systems._pi(z[1])
    return SVector(_wrap_diff(z[h + 1] - z[1], tp), _wrap_diff(z[h + 2] - z[2], tp))
end
@inline _pair_delta(::Any, z::SVector{S}) where {S} = SVector(ntuple(i -> z[S ÷ 2 + i] - z[i], Val(S ÷ 2)))

"""
    pair_distance(sys, z) -> T

Euclidean norm (`safe_norm`) of the minimal-image separation of the pair state `z = (state A; state B)` of the base map `sys`:
for `StandardMap` and `WrappedStandardMap` each coordinate difference is reduced modulo 2 pi (any lifted magnitude), for other
maps it is the plain difference.
"""
pair_distance(sys, z::SVector) = safe_norm(_pair_delta(sys, z))

# State B of the pair: A plus the displacement, reduced with the counts of A in wrapped mode (the construction of
# `Lyapunov.displaced`, so the lifted separation is exactly the minimal-image one).
_shadow(::Any, zA::SVector, delta::SVector) = zA + delta
function _shadow(::WrappedStandardMap, zA::SVector{4}, delta::SVector{2})
    tp = 2 * Systems._pi(zA[1])
    x, a = Systems._reduce(zA[1] + delta[1], tp)
    p, b = Systems._reduce(zA[2] + delta[2], tp)
    return SVector(x, p, zA[3] + a, zA[4] + b)
end

# Dx must widen exactly into T (the rule of `add`); finiteness and sign are the separation check of `pair_start`.
function _widen_dx(::Type{T}, Dx) where {T}
    Dx isa T || Dx isa Integer || (Dx isa AbstractFloat && promote_type(T, typeof(Dx)) === T) ||
        throw(ArgumentError("Dx of type $(typeof(Dx)) would be narrowed to $T"))
    return T(Dx)
end

"""
    pair_start(sys, x0, u, Dx; orbit = 0) -> (z, Dx, Dx_r, mode)

Start of a pair on the base map `sys` (spec Algorithm 1, RL13): `z` the pair state (state A; state B) with A = `x0`
(`torus_state(x0)` for `WrappedStandardMap`) and B = A + `Dx u` (reduced with the counts of A in wrapped mode), `Dx` the
nominal separation in the element type of `x0` (an `ArgumentError` if it would be narrowed), `Dx_r` the realized one
(`pair_distance` of `z`) and `mode` `"wrapped"` or `"lifted"`. A `Dx` that is not finite and positive, or a `Dx_r` that is
zero, not finite or outside `[Dx/2, 2 Dx]`, throws `NumericalHealthError("separation", ...)` with the pair id `orbit`
before any map call. `u` is converted to the element type of `x0` (a direction rounded in a narrower type stays rounded).
In lifted mode the separation is the minimal image modulo 2 pi (the L2 rule), so a lifted gap of a multiple of 2 pi reads as
0 and is rejected (spec Amendment D77 item 5).
"""
function pair_start(sys, x0::SVector{s, T}, u::SVector, Dx; orbit::Integer = 0) where {s, T}
    dx = _widen_dx(T, Dx)
    fail(value, thr) = throw(NumericalHealthError("separation", Float64(value), Float64(thr), 0, orbit, :cpu))
    (isfinite(dx) && dx > 0) || fail(dx, isfinite(dx) ? 0.0 : Inf)
    zA = sys isa WrappedStandardMap ? torus_state(x0) : x0
    z = vcat(zA, _shadow(sys, zA, dx * SVector{length(u), T}(u)))
    dr = pair_distance(sys, z)
    (isfinite(dr) && dr > 0) || fail(dr, isfinite(dr) ? 0.0 : Inf)
    (dr >= dx / 2 && dr <= 2 * dx) || fail(dr, dx)
    return (z = z, Dx = dx, Dx_r = dr, mode = sys isa WrappedStandardMap ? "wrapped" : "lifted")
end
