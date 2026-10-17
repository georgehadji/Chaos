# Periodic orbits of the lifted standard map: residual, Newton step and stability (E3, E4). Scalar BigFloat arithmetic at the
# ambient precision; callers set the precision with `setprecision` (spec Compute: BigFloat is CPU only).

"""
    classify_trace(t, rho) -> :elliptic | :parabolic | :hyperbolic

Stability from the trace of the monodromy matrix (E4, det = 1): parabolic if `||t| - 2| <= rho`, else elliptic for `|t| < 2`
and hyperbolic for `|t| > 2`.
"""
function classify_trace(t, rho)
    a = abs(t)
    abs(a - 2) <= rho && return :parabolic
    return a < 2 ? :elliptic : :hyperbolic
end

"""
    orbit_exponent(t, q)

Exponent per map step of a generic tangent vector along an orbit of period `q` with trace `t` (E4):
`ln |mu_+| / q` with `|mu_+| = (|t| + sqrt(t^2 - 4)) / 2` for `|t| > 2`, else zero.
"""
function orbit_exponent(t, q::Integer)
    a = abs(t)
    a > 2 || return zero(t)
    return log((a + sqrt(a * a - 4)) / 2) / q
end

"Acceptance rule on the return residual: finite and at most `rho` (E3)."
accept_residual(res, rho) = isfinite(res) && res <= rho

"""
    newton_correction(J, F)

Newton step `Delta = -(J - I)^-1 F` for a monodromy matrix `J` with `det J = 1`, by the closed form of E3 (`det(J - I) = 2 - tr J`).
`nothing` when `tr J = 2` within four units of the element type's epsilon (the step does not exist).
"""
newton_correction(J, F) = _newton_step(J[1, 1], J[1, 2], J[2, 1], J[2, 2], F[1], F[2])

function _newton_step(a, b, c, d, f1, f2)
    den = 2 - (a + d)
    abs(den) <= 4 * eps(typeof(den)) && return nothing
    return (-((d - 1) * f1 - b * f2) / den, -((a - 1) * f2 - c * f1) / den)
end

# Lifted iterate T^q(x, p) and the monodromy matrix (a b; c d) = DT(z_{q-1}) ... DT(z_0) (E1, E3).
function _iterate_tangent(K, x, p, q::Integer)
    a, b, c, d = one(x), zero(x), zero(x), one(x)
    for _ in 1:q
        cc = K * cos(x)
        a, b, c, d = (1 + cc) * a + c, (1 + cc) * b + d, cc * a + c, cc * b + d
        p = p + K * sin(x)
        x = x + p
    end
    return x, p, a, b, c, d
end

# Winding of a displacement: nearest integer, `nothing` when not finite or absurdly large.
function _winding(delta)
    (isfinite(delta) && abs(delta) < 1e9) || return nothing
    return round(Int, delta)
end

# Return residual with the windings of the point itself: (Fx, Fp, m, j, a, b, c, d); `nothing` if not finite.
function _return_map(K, x, p, q::Integer)
    tp = 2 * oftype(x, pi)
    xq, pq, a, b, c, d = _iterate_tangent(K, x, p, q)
    m, j = _winding((xq - x) / tp), _winding((pq - p) / tp)
    (m === nothing || j === nothing) && return nothing
    return xq - x - tp * m, pq - p - tp * j, m, j, a, b, c, d
end

# Residual max-norm of the return map of z with its own windings (E9 divisor test); Inf when not finite.
function _residual(K, x, p, q::Integer)
    r = _return_map(K, x, p, q)
    r === nothing && return Inf
    return Float64(max(abs(r[1]), abs(r[2])))
end

