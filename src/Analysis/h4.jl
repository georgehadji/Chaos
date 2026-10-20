# H4 (E13): the two starts of SalaArtuso2015, and the comparators (Shevchenko2014 fits, Chirikov line) reported without a verdict.

_h4_side(x, ref, tol) = (pass = abs(x.L1 - ref) <= tol && abs(x.L5 - ref) <= tol, L1 = x.L1, L5 = x.L5)

"""
    decide_h4(h4, p) -> (verdict, halfpi, twopi)

`h4` holds lambda_1(T) of L1 and L5 for the start at K = pi/2 (`halfpi`) and K = 2 pi (`twopi`). Supported if both methods are within the
tolerance of the reference at both K, otherwise "did not reproduce" with the numbers.
"""
function decide_h4(h4, p)
    a = _h4_side(h4.halfpi, p.h4_ref_halfpi, p.h4_tol_halfpi)
    b = _h4_side(h4.twopi, p.h4_ref_2pi, p.h4_tol_2pi)
    return (verdict = a.pass && b.pass ? :supported : :did_not_reproduce, halfpi = a, twopi = b)
end

"""
    shevchenko_fit(K)

Fit of lambda_1(K) of Shevchenko2014: 0.1333 K for K < 1 (Eq. 25), 0.469 (K - 1.037)^(1/2) for 1.1 <= K < 4.4 and ln(K/2) + 1/K^2 for
K >= 4.4 (Eq. 26); `nothing` between 1 and 1.1 (no fit).
"""
function shevchenko_fit(K::Real)
    K < 1 && return 0.1333 * K
    K < 1.1 && return nothing
    K < 4.4 && return 0.469 * sqrt(K - 1.037)
    return log(K / 2) + 1 / K^2
end

"The Chirikov line ln(K/2)."
chirikov_line(K::Real) = log(K / 2)
