# One Newton solve of the periodic-point search with the conditioning test of Amendment D51.

# Largest singular value of the 2 x 2 matrix (p q; r s).
function _opnorm2(p, q, r, s)
    S = p * p + q * q + r * r + s * s
    D = p * s - q * r
    return sqrt((S + sqrt(max(S * S - 4 * D * D, zero(S)))) / 2)
end

function _error_bound(a, b, c, d, f1, f2)
    den = abs(2 - (a + d))
    den == 0 && return typeof(den)(Inf)
    return sqrt(f1 * f1 + f2 * f2) * _opnorm2(a - 1, b, c, d - 1) / den
end

"""
    error_bound(J, F)

Position error bound of a Newton solution (Amendment D51): `||F||_2 / sigma_min(J - I)` for a monodromy matrix `J` with
`det J = 1`, where `sigma_min(J - I) = |2 - tr J| / ||J - I||_2` (since `det(J - I) = 2 - tr J`). Infinite when `tr J = 2`.
"""
error_bound(J, F) = _error_bound(J[1, 1], J[1, 2], J[2, 1], J[2, 2], F[1], F[2])

"True when `error_bound(J, F) <= delta_dup / 10` (Amendment D51)."
well_conditioned(J, F, delta_dup) = error_bound(J, F) <= delta_dup / 10

const _NEWTON_LIMIT = 1e4          # lifted coordinates beyond this are abandoned as diverged

# One Newton solve (E3) at the ambient precision from the seed (x0, p0). The windings are those of the seed's own rounding and
# stay fixed. Returns (code, x, p, den) with den = |2 - tr J| at the end (0 when not reached): 0 accepted (residual <= rho,
# polished towards rho^2, error bound <= delta_dup / 10), 1 singular step before reaching rho, 2 not finite or diverged,
# 3 iteration cap, 4 residual accepted but ill-conditioned (bound too large, Amendment D51).
function _newton_seed(K, x0, p0, q::Integer, rho, nmax::Integer, delta_dup)
    r = _return_map(K, x0, p0, q)
    r === nothing && return 2, x0, p0, 0.0
    tp = 2 * oftype(x0, pi)
    m, j = r[3], r[4]
    x, p = x0, p0
    rho2, prev, res, singular = rho * rho, Inf, Inf, false
    a = b = c = d = fx = fp = zero(x0)
    for it in 0:nmax
        xq, pq, a, b, c, d = _iterate_tangent(K, x, p, q)
        fx, fp = xq - x - tp * m, pq - p - tp * j
        res = Float64(max(abs(fx), abs(fp)))
        isfinite(res) || return 2, x, p, 0.0
        (res <= rho2 || (res <= rho && res >= prev) || it == nmax) && break
        step = _newton_step(a, b, c, d, fx, fp)
        step === nothing && (singular = true; break)
        prev = res
        x, p = x + step[1], p + step[2]
        (abs(x) > _NEWTON_LIMIT || abs(p) > _NEWTON_LIMIT) && return 2, x, p, 0.0
    end
    res <= rho || return (singular ? 1 : 3), x, p, 0.0
    den = Float64(abs(2 - (a + d)))
    return (_error_bound(a, b, c, d, fx, fp) <= delta_dup / 10 ? 0 : 4), x, p, den
end
