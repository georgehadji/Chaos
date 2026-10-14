# Householder family, hand-written small-d kernel path (spec lyap-householder-qr E4 to E6, E8, E9, Algorithm 1 and 2).
# Everything is immutable SMatrix/SVector arithmetic and loop bounds are type parameters, so the code is isbits and
# kernel-compatible. The diagonal iteration path only (D30). The library path is in lapack_qr.jl.

"Householder reflectors on the deviation matrix (HH E4, E5), hand-written kernel path, provenance \"V-H\"."
struct HouseholderQR <: Orthonormalizer end
"Householder reflectors applied to the augmented matrix [W, I] (HH E6), kernel path, provenance \"V-A\"."
struct AugmentedQR <: Orthonormalizer end
"Givens plane rotations (HH E8), kernel path, provenance \"V-G\"."
struct GivensQR <: Orthonormalizer end

orthonormalizer_name(::HouseholderQR) = "V-H"
orthonormalizer_name(::AugmentedQR) = "V-A"
orthonormalizer_name(::GivensQR) = "V-G"

# Reflections r(d,k) of HH E9: the k-th column of a square matrix has nothing below the diagonal.
_nref(::Val{d}, ::Val{k}) where {d, k} = k < d ? k : k - 1

# Unit u and signed pivot rho = -sign(a_j) c of the reflector for column j (rows j..d); `false` when c is 0 or not finite.
function _reflector(A::SMatrix{d, m, T}, j) where {d, m, T}
    a = SVector{d, T}(ntuple(i -> i >= j ? A[i + (j - 1) * d] : zero(T), Val(d)))
    c = safe_norm(a)
    (isfinite(c) && c > 0) || return zero(a), c, false
    s = a[j] < 0 ? -one(T) : one(T)
    v = setindex(a, a[j] + s * c, j)
    return v / safe_norm(v), -s * c, true
end

# (I - 2 u u') A
function _reflect(A::SMatrix{d, m, T}, u::SVector{d, T}) where {d, m, T}
    w = SVector{m, T}(ntuple(l -> dot(u, _col(A, l, Val(d))), Val(m)))
    return A - (2 * u) * w'
end

# Reflector sweep over the columns of A (width k for V-H, k + d for V-A): the swept A (pivots rho on the diagonal),
# the reflector columns and the number applied.
function _hh_sweep(A::SMatrix{d, m, T}, ::Val{k}) where {d, m, T, k}
    U, count = zero(SMatrix{d, k, T}), 0
    for j in 1:_nref(Val(d), Val(k))
        u, rho, ok = _reflector(A, j)
        if ok
            A = _reflect(A, u)
            U = _setcol(U, j, u)
            count += 1
        end
        A = setindex(A, rho, j + (j - 1) * d)
    end
    return A, U, count
end

# Q = H_1 ... H_r applied to the first k columns of I_d, reflectors in reverse order.
function _form_q(U::SMatrix{d, k, T}) where {d, k, T}
    E = SMatrix{d, k, T}(ntuple(l -> mod(l - 1, d) == div(l - 1, d) ? one(T) : zero(T), Val(d * k)))
    for j in _nref(Val(d), Val(k)):-1:1
        E = _reflect(E, _col(U, j, Val(d)))
    end
    return E
end

# Q from the augmented block S = Q_full' (columns k+1..k+d of A): Q[a, b] = S[b, a].
_q_from_aug(A::SMatrix{d, m, T}, ::Val{k}) where {d, m, T, k} =
    SMatrix{d, k, T}(ntuple(l -> A[div(l - 1, d) + 1 + (k + mod(l - 1, d)) * d], Val(d * k)))

# Sign fix (HH E5) and status. The top-left k x k block of A is R with the pivots on its diagonal. A negative pivot flips
# the row of R and the column of Q; a pivot that is zero or not finite gives QR_PIVOT with the pivot in R[column, column].
function _finish(Q::SMatrix{d, k, T}, A::SMatrix{ld, m, T}, count) where {d, k, ld, m, T}
    sg = SVector{k, T}(ntuple(j -> A[j + (j - 1) * ld] < 0 ? -one(T) : one(T), Val(k)))
    R = SMatrix{k, k, T}(ntuple(l -> (i = mod(l - 1, k) + 1; j = div(l - 1, k) + 1; i > j ? zero(T) : sg[i] * A[i + (j - 1) * ld]), Val(k * k)))
    Qs = SMatrix{d, k, T}(ntuple(l -> Q[l] * sg[div(l - 1, d) + 1], Val(d * k)))
    for j in 1:k
        p = R[j, j]
        (isfinite(p) && p > 0) || return Qs, R, count, QR_PIVOT, j
    end
    return Qs, R, count, QR_OK, 0
end

function orthonormalize(::HouseholderQR, W::SMatrix{d, k, T}) where {d, k, T}
    A, U, count = _hh_sweep(W, Val(k))
    return _finish(_form_q(U), A, count)
end

# [W, I_d]
_augment(W::SMatrix{d, k, T}) where {d, k, T} =
    SMatrix{d, k + d, T}(ntuple(l -> l <= d * k ? W[l] : (mod(l - 1, d) == div(l - 1, d) - k ? one(T) : zero(T)), Val(d * (k + d))))

function orthonormalize(::AugmentedQR, W::SMatrix{d, k, T}) where {d, k, T}
    A, _, count = _hh_sweep(_augment(W), Val(k))
    return _finish(_q_from_aug(A, Val(k)), A, count)
end

# Givens sweep on the rows of [W, I_d], each row an SVector of width k + d. The column j and row i indices are type
# parameters, so every rotation touches two rows at constant positions and the whole sweep unrolls to straight-line code.
# Rotation [c s; -s c] on rows i-1, i zeroes entry (i, j); it is skipped when both entries are zero (HH Algorithm 1).
function _giv_rot(rows::Tuple, count, ::Val{j}, ::Val{i}) where {j, i}
    ra, rb = rows[i - 1], rows[i]
    a, b = ra[j], rb[j]
    rho = safe_norm(SVector(a, b))
    (isfinite(rho) && rho > 0) || return rows, count
    c, s = a / rho, b / rho
    na = setindex(c * ra + s * rb, rho, j)
    nb = setindex(c * rb - s * ra, zero(eltype(rb)), j)
    return ntuple(n -> n == i - 1 ? na : (n == i ? nb : rows[n]), Val(length(rows))), count + 1
end

# The sweep j = 1..J, i = d..j+1 emitted as straight-line code (a generated function: type-level recursion did not infer).
@generated function _giv_sweep(rows::Tuple, ::Val{J}) where {J}
    d = length(rows.parameters)
    body = Any[:(count = 0)]
    for j in 1:J, i in d:-1:(j + 1)
        push!(body, :((rows, count) = _giv_rot(rows, count, Val($j), Val($i))))
    end
    push!(body, :(return rows, count))
    return Expr(:block, body...)
end

function orthonormalize(::GivensQR, W::SMatrix{d, k, T}) where {d, k, T}
    rows = ntuple(i -> SVector{k + d, T}(ntuple(l -> l <= k ? W[i, l] : (l - k == i ? one(T) : zero(T)), Val(k + d))), Val(d))
    rows, count = _giv_sweep(rows, Val(min(k, d - 1)))
    Mt = hcat(rows...)   # Mt[j, i] = rows[i][j]
    R0 = SMatrix{k, k, T}(ntuple(l -> (i = mod(l - 1, k) + 1; j = div(l - 1, k) + 1; i > j ? zero(T) : Mt[j + (i - 1) * (k + d)]), Val(k * k)))
    Q = SMatrix{d, k, T}(ntuple(l -> Mt[k + mod(l - 1, d) + 1 + div(l - 1, d) * (k + d)], Val(d * k)))
    return _finish(Q, R0, count)
end

# Provenance of the whole family (HH AC15e, API contract): the Gram-Schmidt keys plus the iteration path and implementation.
const HHFamily = Union{HouseholderQR, AugmentedQR, GivensQR, LapackQR}
implementation_name(::Union{HouseholderQR, AugmentedQR, GivensQR}) = "kernel"
implementation_name(::LapackQR) = "library"

reducer_provenance(r::QRSpectrum{<:HHFamily}) =
    (estimator = "qr", orthonormalizer = orthonormalizer_name(r.alg), iteration_path = "diagonal",
     implementation = implementation_name(r.alg), orthogonality_tol = r.orthogonality_tol)
