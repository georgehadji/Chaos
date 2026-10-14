# QR-family spectrum estimator (spec lyap-gram-schmidt Algorithm 1 to 5; API contract D34, "QRSpectrum").

# Largest entry of Q'Q - I: the Diagnostics value eps_orth (D4), duplicated because peer layers cannot import each other.
orth_loss(Q::StaticMatrix{<:Any, <:Any, T}) where {T} = maximum(abs, Q' * Q - one(T) * I)

"""
    QRSpectrum(alg; orthogonality_tol = 8.0)

Reducer for `propagate` with `k >= 1`: at each handoff factors the deviation matrix with `alg` (`CGS()`, `MGS()` or
`RGS()`), adds `log(R_ii)` to compensated sums and replaces `W` by `Q`. It owns the orthogonality check
(`orthogonality_tol` times `effective_eps(T) k g`) and the breakdown rules; failures are recorded in the state and
thrown by `reducer_update` or by `propagate`. `reducer_final` gives `(n, sum, comp, lambda, G_max, orth_max, count_max)`.
"""
struct QRSpectrum{A <: Orthonormalizer}
    alg::A
    orthogonality_tol::Float64
end
QRSpectrum(alg::Orthonormalizer; orthogonality_tol::Real = 8.0) = QRSpectrum(alg, Float64(orthogonality_tol))

# `code` 0 ok, 1 pivot, 2 passes, 3 conditioning (breakdown), 4 orthogonality; `n` is the last successful handoff,
# `step` the failing one.
struct QRState{T, K}
    sums::SVector{K, CompensatedSum{T}}
    n::Int
    G_max::T
    orth_max::T
    count_max::Int
    code::Int8
    column::Int
    value::T
    threshold::T
    step::Int
end

function reducer_init(::QRSpectrum, ::Type{T}, ::Val, ::Val{k}) where {T, k}
    sums = SVector{k, CompensatedSum{T}}(ntuple(_ -> CompensatedSum{T}(), Val(k)))
    return QRState{T, k}(sums, 0, zero(T), zero(T), 0, Int8(0), 0, zero(T), zero(T), 0)
end

_failed(st::QRState{T, K}, code, column, value, threshold, n) where {T, K} =
    QRState{T, K}(st.sums, st.n, st.G_max, st.orth_max, st.count_max, Int8(code), column, T(value), T(threshold), n)

# Failure record of a non-OK orthonormalize: the pivot sits in R[column, column]; passes record 5 against the cap 4.
function _factor_failure(st, R, status, column, n)
    status == QR_PASSES && return _failed(st, 2, column, MAX_PASSES + 1, MAX_PASSES, n)
    pivot = Float64(R[column, column])
    return _failed(st, 1, column, pivot, isfinite(pivot) ? 0.0 : Inf, n)
end

reducer_step(r::QRSpectrum, st::QRState, n::Int, x::SVector, W::SMatrix) =
    st.code == 0 ? _qr_step(r, st, n, W, orthonormalize(r.alg, W)) : (st, W)

# One handoff given the factorisation `fact` of W (shared with the LapackQR state, which owns a workspace).
function _qr_step(r::QRSpectrum, st::QRState{T, K}, n::Int, W::SMatrix, fact) where {T, K}
    st.code == 0 || return st, W
    Q, R, count, status, column = fact
    status == QR_OK || return _factor_failure(st, R, status, column, n), W
    Gw = window_Gw(R)
    eps_ = effective_eps(T)
    if checks_conditioning(r.alg) && !(Float64(Gw * eps_) < 0.5)
        return _failed(st, 3, 0, Gw * eps_, 0.5, n), W
    end
    orth = K >= 2 ? orth_loss(Q) : zero(T)
    thr = orthogonality_threshold(r.alg, r.orthogonality_tol, T, K, Gw)
    if K >= 2 && !(orth <= thr)
        return _failed(st, 4, 0, orth, thr, n), W
    end
    sums = SVector{K, CompensatedSum{T}}(ntuple(i -> add(st.sums[i], log(R[i, i])), Val(K)))
    G = _G_from_Gw(Gw, Val(K))
    new = QRState{T, K}(sums, n, max(st.G_max, G), max(st.orth_max, orth), max(st.count_max, count), st.code, 0,
                        st.value, st.threshold, st.step)
    return new, Q
end

function reducer_update(r::QRSpectrum, st::QRState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::QRSpectrum, st::QRState) = st.code != 0

function reducer_error(::QRSpectrum, st::QRState, orbit, backend)
    kind = st.code == 4 ? "orthogonality" : "orthonormalisation breakdown"
    return NumericalHealthError(kind, Float64(st.value), Float64(st.threshold), st.step, orbit, backend, st.column)
end

function reducer_final(::QRSpectrum, st::QRState{T, K}) where {T, K}
    lambda = SVector{K, T}(ntuple(i -> st.n == 0 ? zero(T) : total(st.sums[i]) / st.n, Val(K)))
    return (n = st.n, sum = SVector{K, T}(ntuple(i -> st.sums[i].s, Val(K))),
            comp = SVector{K, T}(ntuple(i -> st.sums[i].c, Val(K))), lambda = lambda,
            G_max = st.G_max, orth_max = st.orth_max, count_max = st.count_max)
end

reducer_provenance(r::QRSpectrum) =
    (estimator = "qr", orthonormalizer = orthonormalizer_name(r.alg), orthogonality_tol = r.orthogonality_tol)
