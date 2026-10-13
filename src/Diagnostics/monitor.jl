# Tangent-stream reducer (spec Algorithm, update rules 1 and 2).

"""
    ResidualMonitor(sys, det_tol, orth_tol)

Reducer for `propagate`: at each handoff it checks the det and symplecticity residuals of `tangent(sys, x)`
(`det_tol` is a multiplier of `eps(T)`, see `det_threshold`) and, for two or more columns, the orthogonality
residual of `W` against the absolute tolerance `orth_tol`. Running maxima are kept in an isbits state.
"""
struct ResidualMonitor{S}
    sys::S
    det_tol::Float64
    orth_tol::Float64
end

# Running maxima, the growth maxima G_max and M_max that scale the sum threshold, and the first-exceedance
# record (code 0 = none; 1 det, 2 symplecticity, 3 orthogonality) read by Tangent through `reducer_failed`.
struct MonitorState{T}
    det_max::T
    symp_max::T
    orth_max::T
    sum_max::T
    pair_max::T
    G_max::T
    M_max::T
    code::Int8
    value::T
    threshold::T
    step::Int
end

reducer_init(::ResidualMonitor, ::Type{T}, ::Val, ::Val) where {T} =
    MonitorState{T}(ntuple(_ -> zero(T), Val(7))..., Int8(0), zero(T), zero(T), 0)

const FAILURE_NAMES = ("det", "symplecticity", "orthogonality")

_bad(value, threshold) = !(isfinite(value) && value <= threshold)

# Orthogonality only for a basis of two or more vectors (spec scope).
_orth(W::SMatrix{<:Any, K, T}) where {K, T} = K >= 2 ? eps_orth(W) : zero(T)

function Tangent.reducer_step(r::ResidualMonitor, st::MonitorState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.code == 0 || return st, W
    M = tangent(r.sys, x)
    vals = (eps_det(M), eps_symp(M), _orth(W))
    thrs = (det_threshold(r.det_tol, T, M), symp_threshold(r.det_tol, T, M), T(r.orth_tol))
    code = _bad(vals[1], thrs[1]) ? 1 : _bad(vals[2], thrs[2]) ? 2 : _bad(vals[3], thrs[3]) ? 3 : 0
    v, th = code == 0 ? (zero(T), zero(T)) : (vals[code], thrs[code])
    st = MonitorState{T}(max(st.det_max, vals[1]), max(st.symp_max, vals[2]), max(st.orth_max, vals[3]),
                         st.sum_max, st.pair_max, st.G_max, max(st.M_max, _maxabs(T, M)), Int8(code), v, th, n)
    return st, W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step` and Tangent raises
# the same error with the orbit id.
function reducer_update(r::ResidualMonitor, st::MonitorState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = Tangent.reducer_step(r, st, n, x, W)
    Tangent.reducer_failed(r, st1) && throw(Tangent.reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

Tangent.reducer_failed(::ResidualMonitor, st::MonitorState) = st.code != 0

function Tangent.reducer_error(::ResidualMonitor, st::MonitorState, orbit, backend)
    finite = isfinite(st.value)
    return NumericalHealthError(FAILURE_NAMES[st.code], Float64(st.value), finite ? Float64(st.threshold) : Inf,
                                st.step, orbit, backend)
end

"""
    monitor_exponents(state, lambda, G, n; pairing_tol, tau, M_max, step, orbit) -> state

Feed the exponents `lambda` (all 2N, estimator order) and the window growth `G` to the monitor state: raises
`sum` when `eps_sum` exceeds `sum_threshold` and records `eps_pair`. The pairing check for N >= 2 needs the
finite-time allowance of the spec (O2) and is left to the caller through [`check_residual`](@ref).
"""
function monitor_exponents(st::MonitorState{T}, lambda::StaticVector, G, n; pairing_tol, tau, M_max, step, orbit) where {T}
    G_max = max(st.G_max, T(G))
    M_hi = max(st.M_max, T(M_max))
    s = eps_sum(lambda)
    check_residual("sum", s, sum_threshold(pairing_tol, T, G_max, tau, M_hi); step, orbit)
    return MonitorState{T}(st.det_max, st.symp_max, st.orth_max, max(st.sum_max, s), max(st.pair_max, eps_pair(lambda)),
                           G_max, M_hi, st.code, st.value, st.threshold, st.step)
end

reducer_final(::ResidualMonitor, st::MonitorState) =
    (det_max = st.det_max, symp_max = st.symp_max, orth_max = st.orth_max, sum_max = st.sum_max, pair_max = st.pair_max)
