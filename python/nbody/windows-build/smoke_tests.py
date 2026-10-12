"""Smoke tests for the Windows builds of REBOUNDx and celmech. Run: python -I -W ignore smoke_tests.py"""
import math
import numpy as np
import rebound, reboundx, celmech

print("rebound", rebound.__version__, "| reboundx", reboundx.__version__, "| celmech", celmech.__version__)
print("reboundx lib:", reboundx.__libpath__)
print("celmech lib :", celmech.clibcelmech._name)
failures = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (" -- " + detail if detail else ""))
    if not ok:
        failures.append(name)


def two_body(integrator, dt=0.01, e=0.1):
    sim = rebound.Simulation()
    sim.G = 1
    sim.add(m=1.)
    sim.add(m=1e-10, a=1., e=e)
    sim.integrator = integrator
    sim.dt = dt
    return sim


# ---- REBOUNDx -----------------------------------------------------------------------------------------
# GR pericentre precession per orbit: 6 pi G M / (c^2 a (1-e^2))
c, e, norb = 100., 0.1, 20
expected = 6 * math.pi / (c ** 2 * (1 - e ** 2)) * norb
for effect, integrator in (("gr_potential", "whfast"), ("gr", "ias15"), ("gr_full", "ias15")):
    sim = two_body(integrator, e=e)
    rebx = reboundx.Extras(sim)
    f = rebx.load_force(effect)
    rebx.add_force(f)
    f.params["c"] = c
    w0 = sim.particles[1].omega
    sim.integrate(2 * math.pi * norb)
    dw = (sim.particles[1].omega - w0 + math.pi) % (2 * math.pi) - math.pi
    check(f"reboundx {effect} ({integrator}) precession", abs(dw / expected - 1) < 0.05,
          f"measured {dw:.5f} rad, expected {expected:.5f} rad, ratio {dw / expected:.4f}")

# operators (src/steppers.c: uses REBOUND's ias15/whfast integrator callbacks) vs plain REBOUND
def final_xy(kind):
    sim = rebound.Simulation()
    sim.G = 1
    sim.add(m=1.)
    sim.add(m=1e-6, a=1., e=0.3)
    if kind == "plain":
        sim.integrator = "ias15"
    else:
        sim.integrator = "none"
        sim.dt = 0.05
        rebx = reboundx.Extras(sim)
        op = rebx.load_operator(kind)
        rebx.add_operator(op, dtfraction=1., timing="post")
    sim.integrate(3.3)
    return np.array([sim.particles[1].x, sim.particles[1].y])

ref = final_xy("plain")
for kind in ("kepler", "ias15"):
    d = np.abs(final_xy(kind) - ref).max()
    check(f"reboundx '{kind}' operator vs plain IAS15", d < 1e-4, f"max |dx| = {d:.2e}")

# natural cubic spline (src/interpolation.c)
rebx = reboundx.Extras(two_body("whfast"))
t = np.linspace(0, 2 * math.pi, 25)
it = reboundx.Interpolator(rebx, t, np.sin(t), "spline")
err = max(abs(it.interpolate(rebx, x) - math.sin(x)) for x in np.linspace(0.1, 6.1, 50))
check("reboundx spline interpolation", err < 1e-4, f"max error {err:.2e}")

# stochastic forces (src/stochastic_forces.c: rand_r shim): deterministic for a fixed seed and non-trivial
def stoch(kappa):
    sim = two_body("ias15")
    sim.rand_seed = 7  # REBOUND seeds from the clock by default
    rebx = reboundx.Extras(sim)
    f = rebx.load_force("stochastic_forces")
    rebx.add_force(f)
    sim.particles[1].params["kappa"] = kappa
    sim.integrate(10.)
    return sim.particles[1].x, sim.particles[1].y
a, b, quiet = stoch(1e-3), stoch(1e-3), stoch(0.)
check("reboundx stochastic_forces reproducible", a == b)
check("reboundx stochastic_forces has an effect", max(abs(np.subtract(a, quiet))) > 1e-8)

# ---- celmech C library ----------------------------------------------------------------------------------
from celmech.miscellaneous import frequency_modified_fourier_transform
dt, N = 0.5, 1024
tt = np.arange(N) * dt
freqs, amps = [0.7, -1.9, 3.1], [1.0, 0.3 * np.exp(0.5j), 0.05]
z = sum(a * np.exp(1j * f * tt) for f, a in zip(freqs, amps))
res = frequency_modified_fourier_transform(tt, z, 3)
got = sorted(res.items(), key=lambda kv: -abs(kv[1]))
ok = all(abs(g[0] - f) < 1e-6 and abs(abs(g[1]) - abs(a)) < 1e-6 for g, (f, a) in zip(got, zip(freqs, amps)))
check("celmech FMFT recovers 3 sinusoids", ok, "; ".join(f"f={k:.8f} amp={abs(v):.6f}" for k, v in got))

# Poisson-series evaluator (src/poisson_series.c, C99 complex arithmetic) against a pure-Python reference
from celmech.c_poisson_series import DFTermSeries
rng = np.random.default_rng(1)
terms = {}
for _ in range(25):
    k = tuple(int(x) for x in rng.integers(-2, 3, size=6))
    zz = tuple(int(x) for x in rng.integers(0, 2, size=4))
    terms[(k, zz)] = float(rng.normal())
ser = DFTermSeries(terms, 1.0, 1.1)
lam = np.array([0.4, 2.2])
xy = rng.normal(size=4) * 0.1 + 1j * rng.normal(size=4) * 0.1
def mypow(x, n):
    return x ** n if n >= 0 else np.conj(x) ** (-n)
def reference():
    # evaluate_term in src/poisson_series.c: exp(i*lam[1])^k0 * exp(i*lam[0])^k1 * prod_i xy_i^k(i+2) * |xy_i|^(2 z_(i+2)%4)
    e = np.exp(1j * lam)
    tot = 0
    for (k, zz), cf in terms.items():
        p = cf * mypow(e[1], k[0]) * mypow(e[0], k[1])
        for i in range(4):
            p *= mypow(xy[i], k[i + 2]) * abs(xy[i]) ** (2 * zz[(i + 2) % 4])
        tot += p
    return tot.real
val, ref_val = ser._evaluate(lam, xy), reference()
check("celmech Poisson series evaluate (C vs Python)", abs(val - ref_val) < 1e-12, f"C={val:.15g} python={ref_val:.15g}")

# Poincare variables / Hamiltonian from a REBOUND simulation
from celmech import Poincare, PoincareHamiltonian
sim = rebound.Simulation()
sim.add(m=1.)
sim.add(m=1e-5, P=1., e=0.02, l=0.3, pomega=0.1)
sim.add(m=1e-5, P=1.52, e=0.03, l=1.1, pomega=2.0)
sim.move_to_com()
pv = Poincare.from_Simulation(sim)
back = pv.to_Simulation()
check("celmech Poincare round trip", all(abs(back.particles[i].a - sim.particles[i].a) < 1e-10 for i in (1, 2)))
Hp = PoincareHamiltonian(pv)
Hp.add_MMR_terms(3, 1, max_order=1, indexIn=1, indexOut=2)
Hp.add_secular_terms(max_order=2, indexIn=1, indexOut=2)
Hp.integrate(5.)
check("celmech PoincareHamiltonian integrates", np.isfinite(Hp.particles[1].e) and Hp.particles[1].e < 0.1,
      f"e1(t=5) = {Hp.particles[1].e:.6f}")

print()
print("ALL PASSED" if not failures else "FAILED: " + ", ".join(failures))
raise SystemExit(1 if failures else 0)
