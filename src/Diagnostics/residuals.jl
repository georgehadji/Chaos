# Residual definitions D1 to D5 and thresholds (spec Equations, Algorithm "Tolerances").

"Determinant residual `|det M - 1|` of a one-step tangent (spec D2)."
eps_det(M) = det_residual(M)

"Symplecticity residual: max abs entry of `M' Omega M - Omega` (spec D3)."
eps_symp(M) = symplectic_residual(M)

"Orthogonality loss: max abs entry of `Q'Q - I` (spec D4)."
eps_orth(Q::StaticMatrix{<:Any, <:Any, T}) where {T} = maximum(abs, Q' * Q - one(T) * I)

"Sum residual `|sum_i lambda_i|` (spec D5)."
eps_sum(lambda::StaticVector) = abs(sum(lambda))

"Pairing residual `max_i |lambda_i + lambda_{2N+1-i}|`, estimator order, no sorting (spec D1)."
function eps_pair(lambda::StaticVector{M}) where {M}
    return maximum(i -> abs(lambda[i] + lambda[M + 1 - i]), 1:(M ÷ 2))
end

_maxabs(::Type{T}, M) where {T} = T(maximum(abs, M))

"Det threshold `tol eps_eff(T) (1 + max|M|)`, `eps_eff = effective_eps`."
det_threshold(tol, ::Type{T}, M) where {T} = T(tol) * effective_eps(T) * (1 + _maxabs(T, M))

"Symplecticity threshold `tol eps_eff(T) N^2 (1 + max|M|)^2`, with `2N = size(M, 1)`."
function symp_threshold(tol, ::Type{T}, M) where {T}
    N = size(M, 1) ÷ 2
    return T(tol) * effective_eps(T) * N^2 * (1 + _maxabs(T, M))^2
end

"Sum and pairing threshold `tol eps_eff(T) (G_max / tau + M_max)`, independent of the run length."
sum_threshold(tol, ::Type{T}, G_max, tau, M_max) where {T} = T(tol) * effective_eps(T) * (T(G_max) / tau + T(M_max))

"""
    check_residual(what, value, threshold; step, orbit, backend=:cpu)

Throw a `NumericalHealthError` when `value` is not finite or exceeds `threshold`; the reported threshold is `Inf`
for a non-finite value.
"""
function check_residual(what::String, value, threshold; step::Integer, orbit::Integer, backend::Symbol = :cpu)
    finite = isfinite(value)
    (finite && value <= threshold) && return nothing
    throw(NumericalHealthError(what, Float64(value), finite ? Float64(threshold) : Inf, step, orbit, backend))
end
