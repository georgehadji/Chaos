# Benettin estimator (spec lyap-benettin, Algorithm 1 to 4, E3, E8): one deviation vector, renormalised at every handoff.

"""
    Benettin()

Reducer for `propagate` with `k = 1`: at each handoff adds `log(alpha)`, `alpha` the norm of the deviation vector, to a
compensated sum and replaces the vector by the unit vector. `reducer_final` gives `(n, sum, comp, lambda)` with
`lambda = (sum + comp) / n`, the finite-time exponent X(n). A zero or non-finite norm is recorded as a health failure
("deviation matrix", column 1) and leaves the sum and the vector unchanged.
"""
struct Benettin end

# `code` 0 = ok, 1 = zero or non-finite norm `value` at handoff `step`; `n` is the last successful handoff.
struct BenettinState{T}
    sum::CompensatedSum{T}
    n::Int
    code::Int8
    value::T
    step::Int
end

function reducer_init(::Benettin, ::Type{T}, ::Val, ::Val{k}) where {T, k}
    k == 1 || throw(ArgumentError("Benettin needs k = 1, got k = $k"))
    return BenettinState{T}(CompensatedSum{T}(), 0, Int8(0), zero(T), 0)
end

function reducer_step(::Benettin, st::BenettinState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    alpha = safe_norm(W[:, 1])
    isfinite(alpha) && alpha > 0 || return BenettinState{T}(st.sum, st.n, Int8(1), alpha, n), W
    return BenettinState{T}(add(st.sum, log(alpha)), n, st.code, st.value, st.step), W / alpha
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::Benettin, st::BenettinState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::Benettin, st::BenettinState) = st.code != 0

function reducer_error(::Benettin, st::BenettinState, orbit, backend)
    value = Float64(st.value)
    return NumericalHealthError("deviation matrix", value, isfinite(value) ? 0.0 : Inf, st.step, orbit, backend, 1)
end

function reducer_final(::Benettin, st::BenettinState{T}) where {T}
    lambda = st.n == 0 ? zero(T) : total(st.sum) / st.n
    return (n = st.n, sum = SVector(st.sum.s), comp = SVector(st.sum.c), lambda = SVector(lambda))
end

reducer_provenance(::Benettin) = (estimator = "benettin",)
