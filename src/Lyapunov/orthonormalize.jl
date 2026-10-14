# Orthonormalisers of the Gram-Schmidt family (spec lyap-gram-schmidt G2, G4, G5) and their window diagnostics (G10).
# Householder, augmented and Givens variants (lyap-householder-qr) add subtypes and `orthonormalize` methods; nothing here
# changes for them. Everything is immutable SMatrix/SVector arithmetic (mutable MMatrix scratch allocated per call).

"Strategy type of the QR family: the factorisation used for `W = Q R` at each handoff."
abstract type Orthonormalizer end
"Classical Gram-Schmidt (G2)."
struct CGS <: Orthonormalizer end
"Modified Gram-Schmidt (G4)."
struct MGS <: Orthonormalizer end
"Repeated Gram-Schmidt, eta = 1/sqrt 2, at most 4 passes per column (G5)."
struct RGS <: Orthonormalizer end

"Status of `orthonormalize`: success, a zero or non-finite pivot, a fifth Gram-Schmidt pass needed."
const QR_OK = Int8(0)
const QR_PIVOT = Int8(1)
const QR_PASSES = Int8(2)

const MAX_PASSES = 4

orthonormalizer_name(::CGS) = "CGS"
orthonormalizer_name(::MGS) = "MGS"
orthonormalizer_name(::RGS) = "RGS"
# Whether the conditioning rule G_w eps >= 1/2 applies (GS Algorithm 5; not for Householder and Givens).
checks_conditioning(::Orthonormalizer) = false
checks_conditioning(::Union{CGS, MGS, RGS}) = true

_col(M, j, ::Val{d}) where {d} = SVector{d}(ntuple(l -> M[l, j], Val(d)))

# Replace column m of the SMatrix M by v.
_setcol(M::SMatrix{d, k, T}, m, v::SVector{d}) where {d, k, T} =
    SMatrix{d, k, T}(ntuple(l -> div(l - 1, d) + 1 == m ? v[mod(l - 1, d) + 1] : M[l], Val(d * k)))

# Reduce column w (number m) against the first m-1 columns of Q; returns (u, r, passes) with r the coefficients
# (R[1:m-1, m]) as a k-vector and passes = 5 when the cap is exceeded.
function _reduce_column(::CGS, Q::SMatrix{d, k, T}, w::SVector{d, T}, m) where {d, k, T}
    u, r = w, zero(SVector{k, T})
    for i in 1:(m - 1)
        qi = _col(Q, i, Val(d))
        c = dot(qi, w)
        r = setindex(r, c, i)
        u = u - c * qi
    end
    return u, r, 1
end

function _reduce_column(::MGS, Q::SMatrix{d, k, T}, w::SVector{d, T}, m) where {d, k, T}
    u, r = w, zero(SVector{k, T})
    for i in 1:(m - 1)
        qi = _col(Q, i, Val(d))
        c = dot(qi, u)
        r = setindex(r, c, i)
        u = u - c * qi
    end
    return u, r, 1
end

function _reduce_column(::RGS, Q::SMatrix{d, k, T}, w::SVector{d, T}, m) where {d, k, T}
    eta = inv(sqrt(T(2)))
    v, r, nold, passes = w, zero(SVector{k, T}), safe_norm(w), 0
    while true
        passes += 1
        vold = v
        for i in 1:(m - 1)
            qi = _col(Q, i, Val(d))
            c = dot(qi, vold)
            r = setindex(r, r[i] + c, i)
            v = v - c * qi
        end
        nnew = safe_norm(v)
        (nnew >= nold * eta || !(isfinite(nnew) && nnew > 0)) && break
        passes == MAX_PASSES && return v, r, MAX_PASSES + 1
        nold = nnew
    end
    return v, r, passes
end

"""
    orthonormalize(alg, W) -> (Q, R, count, status, column)

Factor `W = Q R` with `Q'Q = I`, `R` upper triangular with exact zeros below the diagonal and `R_ii > 0`; column norms
use `safe_norm`. `count` is the number of passes (CGS and MGS 1, RGS the maximum over columns). `status` is `QR_OK`,
`QR_PIVOT` (a pivot is zero or not finite) or `QR_PASSES` (RGS needs a fifth pass); `column` is the failing column.
On failure `R[column, column]` holds the offending pivot and the rest is unspecified.
"""
function orthonormalize(alg::Union{CGS, MGS, RGS}, W::SMatrix{d, k, T}) where {d, k, T}
    Q, R = zero(SMatrix{d, k, T}), zero(SMatrix{k, k, T})
    count, status, column = 1, QR_OK, 0
    for m in 1:k
        u, r, passes = _reduce_column(alg, Q, _col(W, m, Val(d)), m)
        gamma = safe_norm(u)
        R = _setcol(R, m, setindex(r, gamma, m))
        if passes > MAX_PASSES
            status, column = QR_PASSES, m
        elseif !(isfinite(gamma) && gamma > 0)
            status, column = QR_PIVOT, m
        end
        status == QR_OK || break
        Q = _setcol(Q, m, u / gamma)
        count = max(count, passes)
    end
    return Q, R, count, status, column
end

"""
    window_Gw(R) -> G_w

`G_w = ||R||_F ||R^-1||_F` of the upper triangle `R` (spec G10), by back substitution.
"""
function window_Gw(R::SMatrix{k, k, T}) where {k, T}
    Ri = zero(SMatrix{k, k, T})
    for j in 1:k
        Ri = setindex(Ri, inv(R[j, j]), j + (j - 1) * k)
        for i in (j - 1):-1:1
            s = zero(T)
            for l in (i + 1):j
                s += R[i, l] * Ri[l, j]
            end
            Ri = setindex(Ri, -s / R[i, i], i + (j - 1) * k)
        end
    end
    return sqrt(sum(abs2, R)) * sqrt(sum(abs2, Ri))
end

# G from G_w: k = 1 gives 1; k = 2 solves kappa + 1/kappa = G_w; k >= 3 keeps the upper bound G_w (API Q5).
_G_from_Gw(Gw::T, ::Val{1}) where {T} = one(T)
_G_from_Gw(Gw::T, ::Val{2}) where {T} = (Gw + sqrt(max(Gw^2 - 4, zero(T)))) / 2
_G_from_Gw(Gw, ::Val) = Gw

"""
    window_G(R) -> G

Condition number `sigma_max / sigma_min` of the triangle `R` (exact for `k <= 2`), `G_w` as an upper bound for `k >= 3`.
"""
window_G(R::SMatrix{k}) where {k} = _G_from_Gw(window_Gw(R), Val(k))

_gfactor(::Union{CGS, MGS}, Gw) = Gw
_gfactor(::Orthonormalizer, Gw) = one(Gw)

"""
    orthogonality_threshold(alg, tol, T, k, Gw) -> T

`tol * effective_eps(T) * k * g`, `g = Gw` for CGS and MGS and 1 for the other variants (spec Tolerances).
"""
function orthogonality_threshold(alg::Orthonormalizer, tol::Real, ::Type{T}, k::Integer, Gw) where {T}
    return T(tol) * effective_eps(T) * k * _gfactor(alg, T(Gw))
end
