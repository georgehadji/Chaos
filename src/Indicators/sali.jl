# SALI reducer (spec sali, Algorithm 2 to 4, A7): two deviation vectors, each renormalised to unit length at every handoff.

"""
    SALI()

Reducer for `propagate` with `k = 2` (any other `k` throws `ArgumentError`): at each handoff adds `log(alpha_i)`, `alpha_i` the
scaled norm of column `i`, to a compensated sum `S_i`, replaces `W` by its column-normalised matrix and evaluates
`SALI(n) = min(d_-, d_+)` from the stored unit columns (no second normalisation, so an exactly aligned pair stays at 0,
spec A16). `reducer_final` gives `(n, sali, L, sum, comp, floor)`: `sali` the SALI at the last handoff `n`, `L = -(S_1 + S_2)`
the log form `ln GALI_2` of an orthonormal initial frame (spec A7, Amendment D67; add `ln abs(det W0)` for another frame),
`sum` and `comp` the two sums and their compensation terms, `floor` the measurement floor [`sali_floor`](@ref). A zero or
non-finite column norm is recorded as a health failure ("deviation matrix", the column) and leaves the sums and `W` unchanged.
"""
struct SALI end

# `code` 0 = ok, 1 = zero or non-finite norm `value` of column `column` at handoff `step`; `n` is the last successful handoff.
struct SALIState{T}
    s1::CompensatedSum{T}
    s2::CompensatedSum{T}
    sali::T
    value::T
    n::Int
    step::Int
    code::Int8
    column::Int8
end

function reducer_init(::SALI, ::Type{T}, ::Val, ::Val{k}) where {T, k}
    k == 2 || throw(ArgumentError("SALI needs k = 2, got k = $k"))
    return SALIState{T}(CompensatedSum{T}(), CompensatedSum{T}(), zero(T), zero(T), 0, 0, Int8(0), Int8(0))
end

function reducer_step(::SALI, st::SALIState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    c1, c2 = W[:, 1], W[:, 2]
    a1, a2 = scaled_norm(c1), scaled_norm(c2)
    ok1, ok2 = isfinite(a1) && a1 > 0, isfinite(a2) && a2 > 0
    if !(ok1 && ok2)
        column, value = ok1 ? (Int8(2), a2) : (Int8(1), a1)
        return SALIState{T}(st.s1, st.s2, st.sali, value, st.n, n, Int8(1), column), W
    end
    u1, u2 = c1 / a1, c2 / a2
    d_minus, d_plus = unit_indices(u1, u2)
    st1 = SALIState{T}(add(st.s1, log(a1)), add(st.s2, log(a2)), min(d_minus, d_plus), st.value, n, st.step, st.code, st.column)
    return st1, hcat(u1, u2)
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::SALI, st::SALIState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::SALI, st::SALIState) = st.code != 0

function reducer_error(::SALI, st::SALIState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, Int(st.column))
end

function reducer_final(::SALI, st::SALIState{T}) where {T}
    s, c = SVector(st.s1.s, st.s2.s), SVector(st.s1.c, st.s2.c)
    L = st.n == 0 ? zero(T) : -((s[1] + s[2]) + (c[1] + c[2]))
    return (n = st.n, sali = st.sali, L = L, sum = s, comp = c, floor = sali_floor(T))
end

reducer_provenance(::SALI) = (indicator = "sali",)
