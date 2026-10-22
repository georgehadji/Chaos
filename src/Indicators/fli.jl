# FLI reducer (spec fli-ofli, Algorithm 1 to 5, Q1 to Q6): running maxima of the log norms of two deviation vectors.

"""
    FLI(; replace = true)

Fast Lyapunov Indicator reducer (log form) for `propagate` with `k = 2` (any other `k` throws `ArgumentError`). At each handoff it
reads the two arriving, unnormalised columns, adds `log(alpha_i)`, `alpha_i` the scaled norm of column `i`, to a compensated sum
`S_i` (the same operations as [`SALI`](@ref) and `Benettin`, so the sums are bitwise equal to theirs), and keeps two running
maxima that include the `s = 0` term 0: `FLI_B = max(0, S_1, S_2)` over all handoffs (the basis maximum, primary) and
`FLI_1 = max(0, S_1)` (the first column, secondary). The maximiser of each is the handoff of the last strict increase.

`replace` is the chain position of the one-replacer rule: with `replace = true` (alone in a chain, or last) `W` is replaced by
its column-normalised matrix; with `replace = false` the columns pass on unchanged, for a chain in which a later reducer
(`SALI`) replaces them. The reducer cannot see the chain, so the flag is the caller's.

`reducer_final` gives `(n, fli, fli1, sum, comp, stepB, colB, step1)`: `n` the last handoff, `fli = FLI_B(n)`,
`fli1 = FLI_1(n)`, the sums and compensation terms, and the handoff step (0 is the `s = 0` term) and column at which each
maximum last strictly increased. A zero or non-finite column norm is recorded as a health failure ("deviation matrix", the
column) and leaves sums, maxima and `W` unchanged. The sources cap the unlogged value at `1e20`; no cap is applied here.
"""
struct FLI
    replace::Bool
end
FLI(; replace::Bool = true) = FLI(replace)

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`; `n` is the last successful handoff.
# `f1`, `fB` the running maxima, `step1`, `stepB`, `colB` their maximisers.
struct FLIState{T}
    s1::CompensatedSum{T}
    s2::CompensatedSum{T}
    f1::T
    fB::T
    value::T
    n::Int
    step::Int
    step1::Int
    stepB::Int
    code::Int8
    column::Int8
    colB::Int8
end

function reducer_init(::FLI, ::Type{T}, ::Val, ::Val{k}) where {T, k}
    k == 2 || throw(ArgumentError("FLI needs k = 2, got k = $k"))
    return FLIState{T}(CompensatedSum{T}(), CompensatedSum{T}(), zero(T), zero(T), zero(T), 0, 0, 0, 0, Int8(0), Int8(0), Int8(0))
end

# Running maxima after the sums moved to s1, s2 at handoff n; ties keep the earlier maximiser (strict increase only).
function fli_maxima(st::FLIState{T}, s1::CompensatedSum{T}, s2::CompensatedSum{T}, n::Int) where {T}
    S1, S2 = total(s1), total(s2)
    f1, step1 = S1 > st.f1 ? (S1, n) : (st.f1, st.step1)
    fB, stepB, colB = S1 > st.fB ? (S1, n, Int8(1)) : (st.fB, st.stepB, st.colB)
    S2 > fB && ((fB, stepB, colB) = (S2, n, Int8(2)))
    return f1, step1, fB, stepB, colB
end

function reducer_step(r::FLI, st::FLIState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    c1, c2 = W[:, 1], W[:, 2]
    a1, a2 = safe_norm(c1), safe_norm(c2)
    ok1, ok2 = isfinite(a1) && a1 > 0, isfinite(a2) && a2 > 0
    if !(ok1 && ok2)
        column, value = ok1 ? (Int8(2), a2) : (Int8(1), a1)
        return FLIState{T}(st.s1, st.s2, st.f1, st.fB, value, st.n, n, st.step1, st.stepB, Int8(1), column, st.colB), W
    end
    s1, s2 = add(st.s1, log(a1)), add(st.s2, log(a2))
    f1, step1, fB, stepB, colB = fli_maxima(st, s1, s2, n)
    st1 = FLIState{T}(s1, s2, f1, fB, st.value, n, st.step, step1, stepB, st.code, st.column, colB)
    return st1, r.replace ? hcat(c1 / a1, c2 / a2) : W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::FLI, st::FLIState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::FLI, st::FLIState) = st.code != 0

function reducer_error(::FLI, st::FLIState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, Int(st.column))
end

function reducer_final(::FLI, st::FLIState)
    return (n = st.n, fli = st.fB, fli1 = st.f1, sum = SVector(st.s1.s, st.s2.s), comp = SVector(st.s1.c, st.s2.c),
            stepB = st.stepB, colB = Int(st.colB), step1 = st.step1)
end

reducer_provenance(::FLI) = (indicator = "fli",)

"""
    check_frame(FLI(), W0; orbit = 0) -> defect

Frame check of a supplied initial deviation frame (spec Algorithm 1, Amendment D72 item 7): both columns must be unit vectors,
`abs(norm(w_i) - 1) <= 8 effective_eps(T)` (the convention `S(0) = 0`). Returns the larger defect, else throws
`NumericalHealthError("initial frame", ...)` at step 0 naming the column. Independence and orthogonality are not checked (FLI
does not need them). The caller runs it on every supplied frame before the stream; the seeded frame passes by construction.
"""
function check_frame(::FLI, W0::AbstractMatrix{T}; orbit::Integer = 0) where {T}
    size(W0, 2) == 2 || throw(ArgumentError("FLI needs a frame with 2 columns, got $(size(W0, 2))"))
    threshold = 8 * effective_eps(T)
    defect = zero(T)
    for i in 1:2
        d = abs(safe_norm(W0[:, i]) - 1)
        (isfinite(d) && d <= threshold) || throw(NumericalHealthError("initial frame", Float64(d), Float64(threshold), 0, orbit, :cpu, i))
        defect = max(defect, d)
    end
    return defect
end

"""
    fli_unlogged(F) -> exp(F) or nothing

Display conversion of the log form `F` of the FLI: `exp(F)` in the type of `F`, or `nothing` when the exponential overflows
(above 709.78 in Float64). Never an infinity.
"""
function fli_unlogged(F::AbstractFloat)
    v = exp(F)
    return isfinite(v) ? v : nothing
end
