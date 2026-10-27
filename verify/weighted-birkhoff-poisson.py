"""Criterion 12(a) of weighted-birkhoff, the external form: the Poisson form WB17 evaluated by mpmath quadrature (30 digits), independent of every sum of the module.

    python verify/weighted-birkhoff-poisson.py 50 100 200

For f_j = cos 2 pi (theta0 + j alpha), theta0 = 1/7, alpha = sqrt 2 - 1, the weighted average of the first window of the weight exp(-1/(t (1 - t)))
is  Re[e^(2 pi i theta0) sum_m W(N (alpha - m))] / sum_m W(N m),  W(s) = integral_0^1 exp(-1/(t (1 - t))) e^(2 pi i s t) dt,  the sums over the m with
|N (alpha - m)| <= 500 (resp. |N m| <= 500). Output: one line per N, "N value" with 28 significant digits. Needs mpmath (numpy, scipy, mpmath, sympy are the
installed packages); nothing is downloaded.
"""
import sys

from mpmath import mp, mpf, exp, pi, sqrt, quad, nstr, floor

mp.dps = 30
SMAX = 500


def W(s):
    """integral_0^1 exp(-1/(t(1-t))) exp(2 pi i s t) dt, one subinterval per unit of |s|."""
    K = int(abs(s)) + 4

    def f(t):
        if t <= 0 or t >= 1:
            return mp.mpc(0)
        return exp(-1 / (t * (1 - t))) * exp(2j * pi * s * t)

    return quad(f, [mpf(i) / K for i in range(K + 1)])


def poisson_avg(N, theta0, alpha):
    ms = range(-(SMAX // N + 2), SMAX // N + 4)
    num = sum((W(N * (alpha - m)) for m in ms if abs(N * (alpha - m)) <= SMAX), mp.mpc(0))
    den = sum((W(mpf(N * m)) for m in ms if abs(N * m) <= SMAX), mp.mpc(0))
    return (exp(2j * pi * theta0) * num).real / den.real


if __name__ == "__main__":
    alpha = sqrt(2) - 1
    theta0 = mpf(1) / 7
    for arg in sys.argv[1:]:
        N = int(arg)
        print(N, nstr(poisson_avg(N, theta0, alpha), 28), flush=True)
