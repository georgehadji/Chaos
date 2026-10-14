# Singular-value exponents for d = k = 2 (spec lyap-svd-stabilized V5, V6; API contract D34, "SVDSpectrum").
# The RGS window step of QRSpectrum is reused unchanged, so S_i and lambda_qr are bitwise those of QRSpectrum(RGS()) (SVD AC3a).

const SVD_GUARD_MESSAGE = "singular-value exponents are specified for d = 2 and k = d only; the general-d algorithm is held"

"""
    log_singular_values(S, c) -> SVector{2}

`(ln sigma_1, ln sigma_2)` of `P = [[e^S1, c e^S1], [0, e^S2]]` from the QR totals `S` and the ratio `c` (spec V6): a sum of
squares under the root, so nothing overflows however large `S[1]` is and `sigma_1 = sigma_2` is safe.
"""
function log_singular_values(S::SVector{2, T}, c::T) where {T}
    u = exp(2 * (S[2] - S[1]))
    c2 = c * c
    F = 1 + c2 + u
    D = (1 + c2 - u)^2 + 4 * c2 * u
    l1 = S[1] + log((F + sqrt(D)) / 2) / 2
    return SVector{2, T}(l1, S[1] + S[2] - l1)
end

"""
    SVDSpectrum(; orthogonality_tol = 8.0)

Reducer for `propagate` with `d = k = 2`: the RGS spectrum of `QRSpectrum(RGS())` plus the scalar `c = b / a` of the
triangular product, from which `reducer_final` gives the singular-value exponents `lambda` (decreasing) beside
`lambda_qr`. Other `d` or `k` throw `ArgumentError` in `reducer_init`.
"""
struct SVDSpectrum
    orthogonality_tol::Float64
end
SVDSpectrum(; orthogonality_tol::Real = 8.0) = SVDSpectrum(Float64(orthogonality_tol))

struct SVDState{T}
    qr::QRState{T, 2}
    c::T
end

_qr(r::SVDSpectrum) = QRSpectrum(RGS(), r.orthogonality_tol)

function reducer_init(r::SVDSpectrum, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k}
    (d == 2 && k == 2) || throw(ArgumentError(SVD_GUARD_MESSAGE))
    return SVDState{T}(reducer_init(_qr(r), T, Val(2), Val(2)), zero(T))
end

function reducer_step(r::SVDSpectrum, st::SVDState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    st.qr.code == 0 || return st, W
    q = _qr(r)
    fact = orthonormalize(RGS(), W)     # one factorisation, shared with the QR step (`_qr_step`)
    R = fact[2]
    qr1, W1 = _qr_step(q, st.qr, n, W, fact)
    (qr1.code == 0 && fact[4] == QR_OK) || return SVDState{T}(qr1, st.c), W1
    S1, S2 = total(st.qr.sums[1]), total(st.qr.sums[2])     # totals before window j (V5)
    return SVDState{T}(qr1, st.c + (R[1, 2] / R[1, 1]) * exp(S2 - S1)), W1
end

function reducer_update(r::SVDSpectrum, st::SVDState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(r::SVDSpectrum, st::SVDState) = reducer_failed(_qr(r), st.qr)
reducer_error(r::SVDSpectrum, st::SVDState, orbit, backend) = reducer_error(_qr(r), st.qr, orbit, backend)

function reducer_final(r::SVDSpectrum, st::SVDState{T}) where {T}
    f = reducer_final(_qr(r), st.qr)
    S = SVector{2, T}(total(st.qr.sums[1]), total(st.qr.sums[2]))
    ls = log_singular_values(S, st.c)
    lambda = st.qr.n == 0 ? zero(SVector{2, T}) : ls / st.qr.n
    return (n = f.n, sum = f.sum, comp = f.comp, lambda_qr = f.lambda, lambda = lambda, c = st.c,
            u = exp(2 * (S[2] - S[1])), G_max = f.G_max, orth_max = f.orth_max, count_max = f.count_max)
end

reducer_provenance(r::SVDSpectrum) = (estimator = "svd", orthonormalizer = "RGS", d = 2, orthogonality_tol = r.orthogonality_tol)
