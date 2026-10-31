# Independent recomputation of the calibration of study p3-abstention (oracle O15 of the analysis spec, Amendment D93).
# Python 3 standard library plus numpy; exact fractions.Fraction for every rational (counts, shares, J, costs, thresholds, p-values).
# Written from analysis-spec.md (Amendment D93: E47 to E57; D89 E23 to E27, E29) and the preregistration, not from the module.
#
# Usage (from the repository root):
#   python verify/analysis-p3-oracle.py <export folder> <result json> [--p3 table] [--p2 table] [--skip-pp] [--selftest]
# Input: <export folder>/p3-export.jsonl written by verify/analysis-p3-export.jl (HDF5 reader, written from schema/run-store.md).
# Output: the result JSON (rationals as "a/b" texts in lowest terms, counts as integers, reals as binary64 numbers).
#
# Readings where the amendment leaves a choice are marked "READING" below and collected in the result under "readings".
import sys, os, json, math, hashlib, bisect, time, itertools, argparse
from fractions import Fraction
from statistics import NormalDist
import numpy as np
import tomllib

T0 = time.time()
def say(*a):
    print("[%6.1f s]" % (time.time() - T0), *a, flush=True)

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
P2PAR_DEFAULT = os.path.join(ROOT, "workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml")
P3PAR_DEFAULT = os.path.join(ROOT, "workspaces/study/stages/01-design/references/p3-abstention-calibration-parameters.toml")

READINGS = []
def reading(txt):
    if txt not in READINGS:
        READINGS.append(txt)

def fr(x):
    """exact rational of the decimal text of a number (the text of the value; floats through their shortest repr)"""
    if isinstance(x, bool):
        raise TypeError("bool")
    if isinstance(x, int):
        return Fraction(x)
    return Fraction(repr(float(x)))

def rtext(q):
    q = Fraction(q)
    return "%d/%d" % (q.numerator, q.denominator)

# ====================================================================================== generators (E55)
M64 = (1 << 64) - 1

def splitmix_next(state):
    state = (state + 0x9E3779B97F4A7C15) & M64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    z ^= z >> 31
    return state, z

def splitmix_state(seed):
    s = seed & M64
    out = []
    for _ in range(4):
        s, z = splitmix_next(s)
        out.append(z)
    return out

def stream_seed(test_seed, purpose, unit):
    txt = "chaoslab-analysis-1|%s|p3|%s|%s" % (test_seed, purpose, unit)
    return int.from_bytes(hashlib.sha256(txt.encode("utf-8")).digest()[:8], "big")

def rotl(x, k):
    return ((x << k) | (x >> (64 - k))) & M64

class Xoshiro:
    """xoshiro256++ seeded from the 8 bytes by SplitMix64 (scalar version)"""
    def __init__(self, state):
        self.s = list(state)
    @classmethod
    def from_seed(cls, seed):
        return cls(splitmix_state(seed))
    def next(self):
        s = self.s
        result = (rotl((s[0] + s[3]) & M64, 23) + s[0]) & M64
        t = (s[1] << 17) & M64
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]
        s[2] ^= t
        s[3] = rotl(s[3], 45)
        return result
    def uniform(self):
        return ((self.next() >> 11) + 0.5) / 9007199254740992.0
    def randint(self, m):
        """integer in 1..m: x mod m + 1 for the first output x <= 2^64 - 1 - (2^64 mod m) (READING 2 of the test plan; the amendment names no rule)"""
        lim = M64 - ((1 << 64) % m)
        while True:
            x = self.next()
            if x <= lim:
                return x % m + 1

def xo_seed_states(seeds):
    """vectorised SplitMix64 expansion of a uint64 array of seeds into four state arrays"""
    s = seeds.astype(np.uint64).copy()
    out = []
    with np.errstate(over="ignore"):
        for _ in range(4):
            s = s + np.uint64(0x9E3779B97F4A7C15)
            z = s
            z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
            z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
            z = z ^ (z >> np.uint64(31))
            out.append(z)
    return out

class XoV:
    """xoshiro256++ for many independent streams at once (one stream per array element)"""
    def __init__(self, seeds):
        self.s0, self.s1, self.s2, self.s3 = xo_seed_states(seeds)
    def next(self):
        s0, s1, s2, s3 = self.s0, self.s1, self.s2, self.s3
        with np.errstate(over="ignore"):
            a = s0 + s3
            result = ((a << np.uint64(23)) | (a >> np.uint64(41))) + s0
            t = s1 << np.uint64(17)
            s2 = s2 ^ s0
            s3 = s3 ^ s1
            s1 = s1 ^ s2
            s0 = s0 ^ s3
            s2 = s2 ^ t
            s3 = (s3 << np.uint64(45)) | (s3 >> np.uint64(19))
        self.s0, self.s1, self.s2, self.s3 = s0, s1, s2, s3
        return result
    def uniform(self):
        k = (self.next() >> np.uint64(11)).astype(np.float64)
        return (k + 0.5) / 9007199254740992.0

# AS241 (PPND16, Wichura 1988), vectorised; the same algorithm as statistics.NormalDist.inv_cdf
_A = [3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3, 1.3731693765509461125e4,
      4.5921953931549871457e4, 6.7265770927008700853e4, 3.3430575583588128105e4, 2.5090809287301226727e3]
_B = [1.0, 4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3,
      2.1213794301586595867e4, 3.9307895800092710610e4, 2.8729085735721942674e4, 5.2264952788528545610e3]
_C = [1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504,
      1.27045825245236838258, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4]
_D = [1.0, 2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1,
      1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9]
_E = [6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 2.96560571828504891230e-1,
      2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7]
_F = [1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
      7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15]

def _poly(c, x):
    r = np.zeros_like(x) + c[7]
    for j in range(6, -1, -1):
        r = r * x + c[j]
    return r

def ppnd16(p):
    p = np.asarray(p, dtype=np.float64)
    q = p - 0.5
    out = np.empty_like(p)
    cen = np.abs(q) <= 0.425
    if cen.any():
        r = 0.180625 - q[cen] * q[cen]
        out[cen] = q[cen] * _poly(_A, r) / _poly(_B, r)
    ten = ~cen
    if ten.any():
        pt = p[ten]; qt = q[ten]
        r = np.where(qt < 0, pt, 1.0 - pt)
        r = np.sqrt(-np.log(r))
        near = r <= 5.0
        val = np.empty_like(r)
        if near.any():
            rr = r[near] - 1.6
            val[near] = _poly(_C, rr) / _poly(_D, rr)
        far = ~near
        if far.any():
            rr = r[far] - 5.0
            val[far] = _poly(_E, rr) / _poly(_F, rr)
        out[ten] = np.where(qt < 0, -val, val)
    return out

def normal_inv(u):
    return ppnd16(u)

# ====================================================================================== exact helpers
def round_half_away(q):
    q = Fraction(q)
    n = math.floor(abs(q) + Fraction(1, 2))
    return n if q >= 0 else -n

def ceil_step(x, step):
    x = Fraction(x)
    return math.ceil(x / step) * step

def mcnemar_tail(m, k):
    """P(Bin(m, 1/2) >= k), exact"""
    if k <= 0:
        return Fraction(1)
    if k > m:
        return Fraction(0)
    return Fraction(sum(math.comb(m, j) for j in range(k, m + 1)), 2 ** m)

def mcnemar_critical_table(mmax, alpha):
    """crit[m] = smallest k with P(Bin(m,1/2) >= k) <= alpha (m + 1 if none, i.e. no k in 0..m; m = 0 gives 1: the tail at k = 0 is 1 and the decision needs m >= 1 anyway)"""
    crit = [1] * (mmax + 1)
    an, ad = alpha.numerator, alpha.denominator
    for m in range(1, mmax + 1):
        # sum from the top: tail(k) * 2^m = sum_{j>=k} C(m,j); want smallest k with tail*ad <= an*2^m
        total = 2 ** m
        c = 1          # C(m, m)
        tail = 0
        k = m
        best = m + 1
        while k >= 0:
            tail += c
            if tail * ad <= an * total:
                best = k
            else:
                break
            if k > 0:
                c = c * k // (m - k + 1)
            k -= 1
        crit[m] = best
    return crit

def hanley_mcneil(A, na, nn):
    Q1 = A / (2 - A)
    Q2 = 2 * A * A / (1 + A)
    return (A * (1 - A) + (na - 1) * (Q1 - A * A) + (nn - 1) * (Q2 - A * A)) / (na * nn)

def mw_auc(fail_scores, non_scores):
    """Mann-Whitney AUC with ties one half, exact; scores are Fractions; None if a class is empty"""
    if not fail_scores or not non_scores:
        return None
    vals = sorted(set(fail_scores) | set(non_scores))
    cf = {v: 0 for v in vals}; cn = {v: 0 for v in vals}
    for s in fail_scores: cf[s] += 1
    for s in non_scores: cn[s] += 1
    below_n = 0
    twice = 0    # 2 * (pairs with fail > non) + (ties)
    for v in vals:
        twice += cf[v] * (2 * below_n + cn[v])
        below_n += cn[v]
    return Fraction(twice, 2 * len(fail_scores) * len(non_scores))

# ====================================================================================== self test of the building blocks (O17 d, e, f)
def selftest():
    ok = True
    def chk(name, cond):
        nonlocal ok
        say(("PASS " if cond else "FAIL ") + name)
        ok = ok and cond
    chk("O17 d rounding half away 71.25->71, 308.75->309, 237.5->238", (round_half_away(Fraction("71.25")), round_half_away(Fraction("308.75")), round_half_away(Fraction("237.5"))) == (71, 309, 238))
    chk("O17 d ceil5 71->75, 75->75, 76->80", (ceil_step(71, 5), ceil_step(75, 5), ceil_step(76, 5)) == (75, 75, 80))
    chk("O17 b Hanley McNeil (1/2,2,2) = 5/48, A=1 gives 0", hanley_mcneil(Fraction(1, 2), 2, 2) == Fraction(5, 48) and hanley_mcneil(Fraction(1), 5, 7) == 0)
    chk("O17 b Hanley McNeil hand cases of the test plan", hanley_mcneil(Fraction(4, 5), 10, 20) == Fraction(197, 22500) and hanley_mcneil(Fraction(3, 5), 5, 100) == Fraction(6597, 350000)
        and hanley_mcneil(Fraction(7, 10), 40, 800) == Fraction(7119, 3200000))
    chk("O17 c AUC separated 1, identical 1/2, one tied pair 1/2",
        mw_auc([Fraction(3), Fraction(4)], [Fraction(1), Fraction(2)]) == 1 and mw_auc([Fraction(1)] * 3, [Fraction(1)] * 2) == Fraction(1, 2) and mw_auc([Fraction(1)], [Fraction(1)]) == Fraction(1, 2))
    sh = lambda A: math.sqrt(2) * NormalDist().inv_cdf(A)
    chk("O17 c binormal shift d(1/2) = 0, d(0.55), d(0.9) against the stored constants", sh(0.5) == 0 and abs(sh(0.55) - 0.1777119809885155) < 1e-12 and abs(sh(0.9) - 1.8123876048736471) < 1e-12)
    crit = mcnemar_critical_table(200, Fraction(1, 20))
    good = True
    for m in range(1, 201):
        best = m + 1
        for k in range(0, m + 1):
            if mcnemar_tail(m, k) * 20 <= 1:
                best = k; break
        good = good and crit[m] == best
    chk("O17 e critical-value table equals the direct sum for m = 1..200 (alpha 1/20)", good)
    chk("O17 e exact McNemar tail m = 5: P(>=4) = 6/32, P(>=0) = 1", mcnemar_tail(5, 4) == Fraction(6, 32) and mcnemar_tail(5, 0) == 1)
    # generator vectors of the test plan, "Stream vectors"
    chk("O17 f SplitMix64 of seed 0 starts 0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4", (lambda s: s[:2] == [0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4])(splitmix_state(0)))
    g = Xoshiro([1, 2, 3, 4])
    chk("O17 f xoshiro256++ state (1,2,3,4) starts 41943041, 58720359, 3588806011781223, 3591011842654386", [g.next() for _ in range(4)] == [41943041, 58720359, 3588806011781223, 3591011842654386])
    s = stream_seed(20261030, "optimism", "1")
    chk("O17 f seed of (optimism, 1) = 5692226309264838358", s == 5692226309264838358)
    chk("O17 f state of (optimism, 1)", splitmix_state(s) == [15300154784399154531, 1595416676562361740, 9699484540660912394, 17671760824607871275])
    g = Xoshiro.from_seed(s)
    chk("O17 f first four outputs of (optimism, 1)", [g.next() for _ in range(4)] == [1960354936839448390, 11115532021063746583, 17181476028355880092, 5824767452774032242])
    g = Xoshiro.from_seed(s)
    u = [g.uniform() for _ in range(4)]
    chk("O17 f first four uniforms of (optimism, 1)", u == [0.10627105406820064, 0.6025741982784751, 0.931409681822553, 0.3157612763260232])
    g = Xoshiro.from_seed(stream_seed(20261030, "pp-hb", "0,1"))
    nrm = [float(normal_inv(np.array([g.uniform()]))[0]) for _ in range(4)]
    chk("O17 f first four normals of (pp-hb, (0,1)) (AS241) to 1e-12", all(abs(a - b) < 1e-12 for a, b in zip(nrm, [0.13039229363167548, -0.36652077601596667, 0.17412024595138367, -0.38487600750808904])))
    chk("O17 f seeds of (pp-hc,(20,1)), (pp-hb,(3,17)), (pp-hc,(49,1)), (optimism,6)",
        (stream_seed(20261030, "pp-hc", "20,1"), stream_seed(20261030, "pp-hb", "3,17"), stream_seed(20261030, "pp-hc", "49,1"), stream_seed(20261030, "optimism", "6"))
        == (15836479185517542896, 4982238177473391303, 16471504700301138262, 982380449943276992))
    g = Xoshiro.from_seed(s)
    chk("O17 f first eight integers in 1..36 of (optimism,1) = 23,20,5,19,31,7,13,10", [g.randint(36) for _ in range(8)] == [23, 20, 5, 19, 31, 7, 13, 10])
    g = Xoshiro.from_seed(s)
    chk("O17 f first eight integers in 1..3 of (optimism,1) = 2,2,2,1,1,1,1,1", [g.randint(3) for _ in range(8)] == [2, 2, 2, 1, 1, 1, 1, 1])
    g = Xoshiro.from_seed(stream_seed(20261030, "optimism", "2"))
    chk("O17 f first eight integers in 1..36 of (optimism,2) = 35,13,30,8,5,24,12,30", [g.randint(36) for _ in range(8)] == [35, 13, 30, 8, 5, 24, 12, 30])
    # vectorised generator equals the scalar one
    seeds = np.array([stream_seed(20261030, "pp-hb", "%d,%d" % (j, i)) for j in range(2) for i in range(1, 4)], dtype=np.uint64)
    xv = XoV(seeds); got = [xv.next() for _ in range(5)]
    sc = [Xoshiro.from_seed(int(s)) for s in seeds]
    chk("vectorised xoshiro256++ equals the scalar generator", all(int(got[t][r]) == sc[r].next() for t in range(5) for r in range(len(seeds))))
    rng = np.random.default_rng(1)
    uu = np.concatenate([rng.random(20000), np.array([1e-12, 1e-300, 1 - 1e-9, 0.01, 0.3, 0.975, 0.5, 0.0749, 0.0751])])
    ref = np.array([NormalDist().inv_cdf(float(x)) for x in uu])
    chk("vectorised AS241 equals statistics.NormalDist.inv_cdf to 1e-12 (and the stored quantiles)", np.max(np.abs(ppnd16(uu) - ref) / np.maximum(1, np.abs(ref))) < 1e-12
        and abs(float(ppnd16(np.array([1e-12]))[0]) - (-7.034483825301132)) < 1e-12 and abs(float(ppnd16(np.array([1 - 1e-9]))[0]) - 5.997807019601638) < 1e-12)
    return ok

# ====================================================================================== the study
I7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA"]
CH, RG, UN = 1, 2, 0     # verdict codes: chaotic, regular, undecided

def load_params(p3path, p2path):
    p3 = tomllib.load(open(p3path, "rb"))
    p2 = tomllib.load(open(p2path, "rb"))
    q = {k: v["value"] for k, v in p3.items()}
    p = {k: v["value"] for k, v in p2.items()}
    return q, p

def recorded_step_ge(ts, b):
    i = bisect.bisect_left(ts, b)
    return i if i < len(ts) else None

def recorded_step_near_tenth(ts, n):
    """index of the recorded step nearest to n / 10, a tie to the larger step (E26)"""
    j = bisect.bisect_left(ts, -(-n // 10))
    best = None
    for i in (j - 1, j, j + 1):
        if 0 <= i < len(ts):
            d = abs(10 * ts[i] - n)
            if best is None or d < best[0] or (d == best[0] and ts[i] > ts[best[1]]):
                best = (d, i)
    return best[1]

def isfin(x):
    return x is not None and math.isfinite(x)

def series_verdict(ts, g, b, thr, lam=None, le_floor=None):
    """ratio rule of E26/E27 on a series: g = growth values per recorded step"""
    i = recorded_step_ge(ts, b)
    if i is None:
        return UN
    ip = recorded_step_near_tenth(ts, ts[i])
    gn, gp = g[i], g[ip]
    if not (isfin(gn) and isfin(gp)) or gp <= 0:
        return UN
    r = gn / gp
    if lam is not None:       # LE: lambda_1(n) >= floor and r >= threshold
        return CH if (lam[i] >= le_floor and r >= thr) else RG
    return CH if r >= thr else RG

def flagbit(flags, bits, idx, name):
    return (flags[idx] >> bits.index(name)) & 1

def window_verdict(item, spec, b, thr):
    rows, c, flags, bits = item["rows"], item["c"], item["flags"], item["bits"]
    W = b // 2
    if spec in ("fq1", "pa1"):
        i = rows.index("W=%d" % W)
        D = c["D"][i]
        bad = flagbit(flags, bits, i, "degenerate" if spec == "fq1" else "undefined window")
        if bad or not isfin(D):
            return UN
        return CH if D >= thr else RG
    idxs = [i for i, r in enumerate(rows) if r.startswith("W=%d " % W)]
    if len(idxs) != 4:
        raise RuntimeError("bk1 window %d has %d rows" % (W, len(idxs)))
    if any(flagbit(flags, bits, i, "zero") for i in idxs) or any(not isfin(c["dig"][i]) for i in idxs):
        return UN
    return CH if min(c["dig"][i] for i in idxs) < thr else RG

def load_export(folder):
    orbits, ser, win = {}, {}, {}
    meta = None
    for ln in open(os.path.join(folder, "p3-export.jsonl"), encoding="utf-8"):
        d = json.loads(ln)
        k = d["k"]
        if k == "orbit":
            orbits[d["id"]] = d
        elif k == "s":
            ser[(d["o"], d["spec"], d["prec"], d["hor"])] = d
        elif k == "w":
            win[(d["o"], d["spec"], d["prec"], d["hor"])] = d
        elif k == "meta":
            meta = d
    return meta, orbits, ser, win

def verdicts_of(oid, prec, ser, win, p, N1, NC, budgets, hor_tg=None):
    """verdicts of I_7 at the ten budgets for one orbit and one precision: dict indicator -> list of ten codes (None when the item is absent)"""
    out = {}
    thr, floor = p["ratio_threshold"], p["le_floor"]
    l1 = ser.get((oid, "L1", prec, N1 if prec == "Float64" else hor_tg))
    tg = ser.get((oid, "tg", prec, N1 if prec == "Float64" else hor_tg))
    if l1 is not None:
        ts, lam = l1["t"], l1["c"]["lam1"]
        g = [n * l for n, l in zip(ts, lam)]
        out["LE"] = [series_verdict(ts, g, b, thr, lam, floor) for b in budgets]
    if tg is not None:
        ts = tg["t"]
        out["SALI"] = [series_verdict(ts, [-x for x in tg["c"]["SALI.L"]], b, thr) for b in budgets]
        out["FLI"] = [series_verdict(ts, tg["c"]["FLI.fli"], b, thr) for b in budgets]
        out["MEGNO"] = [series_verdict(ts, tg["c"]["MEGNO.Ybar20_1"], b, thr) for b in budgets]
    for ind, spec, hor, key in (("E1", "fq1", NC, "e1_d_threshold"), ("E2", "pa1", NC + 1, "e2_d_threshold"), ("WBA", "bk1", NC, "wba_digits_threshold")):
        it = win.get((oid, spec, prec, hor))
        if it is not None:
            out[ind] = [window_verdict(it, spec, b, p[key]) for b in budgets]
    return out

# ====================================================================================== relabelling C0 (E49)
POS, NONE, UNRES = "positive growth supported", "none detected on the tested horizon", "unresolved"

def relabel_orbit(o, p, use_diag=True, alt_ratio=False):
    """the reduced rule of D89 E23 b plus the diagnostics condition of E49 c; returns (instability, structure, rho, rho1, rho2) or None when data are missing"""
    if o["relab_missing"]:
        return None
    rl = o["relab"]
    lam5, lam10, lam1 = rl["lam5_N1"], rl["lam5_N1_10"], rl["lam1_N1_rel"]
    N1 = 10_000_000
    inst = UNRES
    if isfin(lam5) and isfin(lam10) and lam10 > 0:
        r = 10 * (lam5 / lam10) if alt_ratio else 10 * lam5 / lam10
        if lam5 >= p["label_lambda_min"] and r >= p["label_growth_ratio"]:
            inst = POS
        elif lam5 < p["label_lambda_min"] and r < p["label_growth_ratio"]:
            inst = NONE
    h = N1 // 2
    x0, xh, xN = rl["x_0"], rl["x_N1_2"], rl["x_N1"]
    rho1 = (xh - x0) / (2 * math.pi * h)
    rho2 = (xN - xh) / (2 * math.pi * (N1 - h))
    if inst == POS:
        struct = "irregular"
    elif inst == NONE and abs(rho1 - rho2) <= p["label_rotation_tol"]:
        struct = "quasiperiodic"
    else:
        struct = UNRES
    if lam1 is not None and lam5 is not None and abs(lam1 - lam5) > p["label_estimator_bound"]:
        inst, struct = UNRES, UNRES
    if use_diag and (o["diag5"] or o["diag1"]):
        inst, struct = UNRES, UNRES
    rho = "C" if (inst == POS and struct == "irregular") else ("R" if (inst == NONE and struct == "quasiperiodic") else "U")
    return inst, struct, rho, rho1, rho2

def old_rho(o):
    inst, st = o["old"][0], o["old"][1]
    return "C" if (inst == POS and st == "irregular") else ("R" if (inst == NONE and st == "quasiperiodic") else "U")

def ckey(d):
    return {",".join(str(x) for x in (k if isinstance(k, tuple) else (k,))): v for k, v in sorted(d.items(), key=lambda kv: str(kv[0]))}

# ====================================================================================== rule space, rule evaluation, costs (E52, E53)
def popcount(x):
    return bin(x).count("1")

def enumerate_theta(q):
    vb = q["voter_bits"]
    nb = len(vb)
    out = []
    idx = 0
    for size in range(q["min_voters"], nb + 1):
        for mask in range(1, 1 << nb):
            if popcount(mask) != size:
                continue
            for qc in range(1, size + 1):
                for qr in range(size // 2 + 1, size + 1):
                    for w in q["w_range"]:
                        for v in q["v_range"]:
                            idx += 1
                            out.append(dict(index=idx, size=size, mask=mask, V=[vb[j] for j in range(nb) if mask >> j & 1], qC=qc, qR=qr, w=w, v=v))
    return out

def closed_form_size(q):
    nb = len(q["voter_bits"])
    return sum(math.comb(nb, s) * s * (s - s // 2) * len(q["w_range"]) * len(q["v_range"]) for s in range(q["min_voters"], nb + 1))

def theta_text(t):
    return "V=%s qC=%d qR=%d w=%d v=%d" % ("+".join(t["V"]), t["qC"], t["qR"], t["w"], t["v"])

class Costs:
    def __init__(self, q):
        self.q = q
        self.spec = q["spec_of"]
        self.c = {s: fr(q["cost_" + s]) for s in ("l1", "tg", "pa1", "bk1", "fq1")}
        self.d64 = {s: fr(q["ratio_d64_" + s]) for s in ("tg", "fq1", "pa1", "bk1")}
        self.B = list(q["budgets"]); self.bP = q["b_primary"]; self.Nc = q["N_c"]
        self.I6 = list(q["voter_bits"])
    def cx(self, X):
        return self.c[self.spec[X]]
    def cost_per_step(self, V):
        return sum(self.c[s] for s in sorted({self.spec[X] for X in V}))
    def CR(self, V, b=None):
        return (self.bP if b is None else b) * self.cost_per_step(V)
    def cap(self):
        return min(self.cx(X) for X in self.I6) * self.Nc
    def matched(self, V):
        CR = self.CR(V)
        out = {}
        for X in self.I6 + ["LE"]:
            cx = self.cx(X)
            up = next((b for b in self.B if b * cx >= CR), None)
            down = max((b for b in self.B if b * cx <= CR), default=None)
            out[X] = dict(b=up, capped=up is None, down=down)
        bLE = out["LE"]["b"]
        le_read = self.Nc if bLE is None else min(bLE, self.Nc)
        le_unmatched = le_read * self.c["l1"] < CR
        return CR, out, le_read, le_unmatched

# ====================================================================================== the analysis (C0 to C8)
def main_analysis(export_folder, p3path, p2path, skip_pp=False, log=say):
    q, p = load_params(p3path, p2path)
    meta, orbits, ser, win = load_export(export_folder)
    N1, NC, NS = q["N_1"], q["N_c"], q["n_s"]
    B = list(q["budgets"]); bP = q["b_primary"]; kP = B.index(bP)
    res = dict(readings=READINGS)
    log("export loaded: orbits", len(orbits), "series", len(ser), "windows", len(win))
    mir = [(k, q[k], p[k2]) for k, k2 in (("budgets", "budgets"), ("b_primary", "b_primary"), ("N_c", "N_c"), ("N_1", "N_1"), ("n_s", "n_s"), ("f_s", "f_s"), ("min_set", "min_set"),
                                           ("gate_size", "gate_size"), ("lam_w2", "w2_lambda_max"), ("lam_s2", "s2_lambda_min"))]
    bad = [k for k, a, b in mir if a != b]
    res["mirror_check"] = dict(ok=not bad, differing=bad)
    assert not bad, bad
    assert q["p2_manifest_sha256"] == meta["p2_manifest_sha256"]
    # ---------------------------------------------------------------- C0
    P_ids = [i for i in sorted(orbits) if orbits[i]["src"] == "P"]
    G_ids = [i for i in sorted(orbits) if orbits[i]["src"] == "G"]
    assert len(P_ids) == 883 and len(G_ids) == 2304
    relab = {}
    for i, o in orbits.items():
        if o["control"]:
            continue
        relab[i] = relabel_orbit(o, p)
    alt_diff = [i for i in relab if relab[i] is not None and relabel_orbit(orbits[i], p, alt_ratio=True)[2] != relab[i][2]]
    nodiag_diff = [i for i in relab if relab[i] is not None and relabel_orbit(orbits[i], p, use_diag=False)[2] != relab[i][2]]
    res["relabel_sensitivity"] = dict(ratio_as_10_times_quotient_label_differs=alt_diff, diagnostics_condition_off_label_differs=nodiag_diff)
    missing = [(i, orbits[i]["src"], orbits[i]["kidx"], it) for i in sorted(relab) if relab[i] is None for it in orbits[i]["relab_missing"]]
    new_counts = {}
    for i, r in relab.items():
        if r is None:
            continue
        o = orbits[i]
        d = new_counts.setdefault((o["src"], o["kidx"]), {"C": 0, "R": 0, "U": 0})
        d[r[2]] += 1
    axes_counts = {}
    for i, r in relab.items():
        if r is None: continue
        axes_counts[(orbits[i]["src"], r[0], r[1])] = axes_counts.get((orbits[i]["src"], r[0], r[1]), 0) + 1
    chg_panel, chg_grid = {}, {}
    for i, r in relab.items():
        if r is None:
            continue
        o = orbits[i]
        if o["src"] == "P":
            old = old_rho(o)
            key = (o["kidx"], old, r[2])
            chg_panel[key] = chg_panel.get(key, 0) + 1
        else:
            rold = relabel_orbit(o, p, use_diag=False)[2]     # the D89 E23 b label recomputed from the p2 items without the diagnostics condition
            key = (o["kidx"], rold, r[2])
            chg_grid[key] = chg_grid.get(key, 0) + 1
    res["relabelling"] = dict(
        n_panel_noncontrol=sum(1 for i in P_ids if not orbits[i]["control"]), n_grid=len(G_ids),
        missing=[dict(id=a, src=b, kidx=c, item=d) for a, b, c, d in missing],
        rho_new_by_source_kidx={"%s,%d" % k: v for k, v in sorted(new_counts.items())},
        axes_counts={"%s|%s|%s" % k: v for k, v in sorted(axes_counts.items())},
        change_panel=ckey(chg_panel), change_grid=ckey(chg_grid),
        panel_changed=sum(v for (k, a, b), v in chg_panel.items() if a != b), grid_changed=sum(v for (k, a, b), v in chg_grid.items() if a != b),
        panel_old_totals={c: sum(v for (k, a, b), v in chg_panel.items() if a == c) for c in "CRU"},
        panel_new_totals={c: sum(v for (k, a, b), v in chg_panel.items() if b == c) for c in "CRU"},
        note_old_labels="old = the base P1 labelling of D89 E23 a recomputed independently (verify/analysis-p3-export.jl); grid old = the E23 b label recomputed without the diagnostics condition")
    log("C0: panel changed", res["relabelling"]["panel_changed"], "grid changed", res["relabelling"]["grid_changed"], "missing", len(missing))
    # ---------------------------------------------------------------- E50 sets
    reading("drop reasons are applied in the order control, drift, excluded, missing, unresolved (the order of the test plan contract); the flags are those of the p2 Float64 items read (sensitivity: union over all items reported)")
    lamref = lambda i: orbits[i]["lam1_N1"]
    lamns = lambda i: orbits[i]["lam1_ns"]
    dropped = []
    cal_ids = []
    flags_any = {i: set(orbits[i]["fl_p2"]) for i in orbits}
    flags_sens = {i: set(orbits[i]["fl_p2"]) | set(orbits[i]["fl_p2_all"]) | set(orbits[i]["fl_p1"]) for i in orbits}
    for i in sorted(orbits):
        o = orbits[i]
        if o["control"]:
            dropped.append((i, "control"))
        elif "drift" in flags_any[i]:
            dropped.append((i, "drift"))
        elif "excluded" in flags_any[i]:
            dropped.append((i, "excluded"))
        elif relab[i] is None or not isfin(lamref(i)) or not isfin(lamns(i)):
            dropped.append((i, "missing"))
        elif relab[i][2] == "U":
            dropped.append((i, "unresolved"))
        else:
            cal_ids.append(i)
    sens_drop = sorted(i for i in orbits if not orbits[i]["control"] and (flags_sens[i] & {"drift", "excluded"}))
    res["flag_sensitivity_noncontrol_orbits_with_drift_or_excluded_on_any_item"] = sens_drop
    rho = {i: relab[i][2] for i in cal_ids}
    lam_s2, lam_w2, f_s = q["lam_s2"], q["lam_w2"], q["f_s"]
    def sets_for(cal, exact=False):
        S2, Ts, W2 = [], [], []
        for i in cal:
            if rho[i] != "C":
                continue
            lr, ln = lamref(i), lamns(i)
            if exact:
                a, b = Fraction(lr), Fraction(ln)
                s2 = a >= fr(lam_s2); ts = b <= fr(f_s) * a; w2 = (a <= fr(lam_w2)) or ts
            else:
                s2 = lr >= lam_s2; ts = ln <= f_s * lr; w2 = (lr <= lam_w2) or ts
            if s2: S2.append(i)
            if ts: Ts.append(i)
            if w2: W2.append(i)
        return S2, Ts, W2
    S2, Ts, W2 = sets_for(cal_ids)
    S2x, Tsx, W2x = sets_for(cal_ids, exact=True)
    res["set_exact_vs_binary64_differences"] = dict(S2=sorted(set(S2) ^ set(S2x)), T_s=sorted(set(Ts) ^ set(Tsx)), W2=sorted(set(W2) ^ set(W2x)))
    reading("the set conditions on lambda_ref and lambda_1(n_s) are evaluated in binary64 as in D89 E29 (f_s x lambda_ref in binary64); the exact-rational variant is compared and reported")
    regular = [i for i in cal_ids if rho[i] == "R"]
    gate = [i for i in sorted(orbits) if orbits[i]["src"] == "P" and not orbits[i]["control"] and orbits[i]["kidx"] == p["gate_k_index"] and orbits[i]["cls"] in p["gate_classes"]]
    assert len(gate) == q["gate_size"], len(gate)
    tsset = set(Ts)
    RW = [i for i in W2 if i not in tsset]
    setd = dict(Cal=cal_ids, S2=S2, T_s=Ts, W2=W2, R_W=RW, regular=regular, gate=gate)
    def bysrc(lst):
        d = {"pooled": len(lst), "panel": sum(1 for i in lst if orbits[i]["src"] == "P"), "grid": sum(1 for i in lst if orbits[i]["src"] == "G")}
        d["by_source_kidx"] = {"%s,%d" % (s, k): sum(1 for i in lst if orbits[i]["src"] == s and orbits[i]["kidx"] == k) for s in "PG" for k in range(0, 9)
                               if any(orbits[i]["src"] == s and orbits[i]["kidx"] == k for i in lst)}
        return d
    res["sets"] = {k: bysrc(v) for k, v in setd.items()}
    res["sets"]["overlap_S2_W2"] = len(set(S2) & set(W2))
    dr = {}
    for i, why in dropped:
        d = dr.setdefault(why, {"total": 0, "panel": 0, "grid": 0}); d["total"] += 1; d["panel" if orbits[i]["src"] == "P" else "grid"] += 1
    res["dropped"] = dict(total=len(dropped), by_reason=dr, rows=[dict(id=i, reason=w, src=orbits[i]["src"], kidx=orbits[i]["kidx"]) for i, w in dropped],
                          unresolved_by_source_kidx={"%s,%d" % (s, k): sum(1 for i, w in dropped if w == "unresolved" and orbits[i]["src"] == s and orbits[i]["kidx"] == k)
                                                      for s in "PG" for k in range(0, 9) if any(w == "unresolved" and orbits[i]["src"] == s and orbits[i]["kidx"] == k for i, w in dropped)})
    assert len(cal_ids) + len(dropped) == 3187
    log("E50: |Cal|", len(cal_ids), "S2", len(S2), "T_s", len(Ts), "W2", len(W2), "R_W", len(RW), "regular", len(regular), "gate", len(gate), "dropped", len(dropped))
    # ---------------------------------------------------------------- E51 verdicts
    cal_set = set(cal_ids)
    rows_ids = list(cal_ids) + [i for i in gate if i not in cal_set]
    rowpos = {i: n for n, i in enumerate(rows_ids)}
    NR = len(rows_ids)
    V = np.zeros((NR, 7, 10), dtype=np.int8)
    for i in rows_ids:
        vd = verdicts_of(i, "Float64", ser, win, p, N1, NC, B)
        for xi, X in enumerate(I7):
            if X not in vd:
                raise RuntimeError("verdict item missing: orbit %d indicator %s" % (i, X))
            V[rowpos[i], xi, :] = vd[X]
    log("E51: verdict table", V.shape, "undecided share", float((V == 0).mean()))
    res["verdict_summary_at_bP_cal"] = {X: {c: int((V[:len(cal_ids), xi, kP] == code).sum()) for c, code in (("chaotic", CH), ("regular", RG), ("undecided", UN))} for xi, X in enumerate(I7)}
    res["verdict_digest_sha256"] = hashlib.sha256(V.tobytes() + ",".join(map(str, rows_ids)).encode()).hexdigest()
    res["verdicts"] = {str(i): {X: "".join("cru"[0 if c == CH else (1 if c == RG else 2)] for c in V[rowpos[i], xi]) for xi, X in enumerate(I7)} for i in rows_ids}
    # ---------------------------------------------------------------- E52 enumeration; rule output on all rows for all theta
    thetas = enumerate_theta(q)
    T = len(thetas)
    assert T == closed_form_size(q) == 2088, T
    res["theta_space"] = dict(size=T, closed_form=closed_form_size(q), by_size={str(s): sum(1 for t in thetas if t["size"] == s) for s in range(3, 7)})
    Vk = V[:, :, kP]
    SC, SR = {}, {}
    for w in q["w_range"]:
        stab = (Vk != 0)
        for j in range(1, w + 1):
            if kP - j < 0:
                raise ValueError("budget b/2^w outside B")
            stab = stab & (V[:, :, kP - j] == Vk)
        SC[w] = ((Vk == CH) & stab).astype(np.int16)
        SR[w] = ((Vk == RG) & stab).astype(np.int16)
    CHany = (Vk == CH).astype(np.int16)
    OUT = np.zeros((T, NR), dtype=np.int8)
    for t in thetas:
        cols = [1 + j for j in range(6) if t["mask"] >> j & 1]   # voter bit j is I7 index j + 1
        nC = SC[t["w"]][:, cols].sum(1); nR = SR[t["w"]][:, cols].sum(1); nCa = CHany[:, cols].sum(1)
        out = np.where(nC >= t["qC"], 1, np.where((nCa <= t["v"]) & (nR >= t["qR"]), 2, 0))
        OUT[t["index"] - 1] = out
    CSm = (OUT == 1).astype(np.float64); RSm = (OUT == 2).astype(np.float64)
    costs = Costs(q)
    cap = costs.cap()
    CRs = [costs.CR(t["V"]) for t in thetas]
    f1 = np.array([c <= cap for c in CRs])
    lcm = math.lcm(*([c.denominator for c in CRs] + [1]))
    cr_int = np.array([int(c * lcm) for c in CRs], dtype=np.int64)
    sizeA = np.array([t["size"] for t in thetas]); wA = np.array([t["w"] for t in thetas]); vA = np.array([t["v"] for t in thetas]); maskA = np.array([t["mask"] for t in thetas])
    qcA = np.array([t["qC"] for t in thetas]); qrA = np.array([t["qR"] for t in thetas])
    nrows_cal = len(cal_ids)
    def mk(lst):
        m = np.zeros(NR, dtype=np.float64)
        for i in lst: m[rowpos[i]] = 1.0
        return m
    mcal = mk(cal_ids); mS2 = mk(S2); mTs = mk(Ts); mW2 = mk(W2); mgate = mk(gate)
    mrhoC = mk([i for i in cal_ids if rho[i] == "C"]); mrhoR = mk([i for i in cal_ids if rho[i] == "R"])
    thr = {k: fr(q[k]) for k in ("e_cal", "s_cal", "a_cal", "c_min", "delta_tol", "rho_hb", "alpha_hb", "p_hb", "p_hc", "delta_c", "rule_rate_floor", "rule_rate_factor", "auc_shrink", "auc_p")}
    def stats_of(wts):
        cols = np.stack([mcal * wts, mS2 * wts, mTs * wts, mW2 * wts, mgate * wts, mrhoC * wts, mrhoR * wts], axis=1)
        cs = CSm @ cols; rs = RSm @ cols
        n = cols.sum(0)
        r = lambda a: np.rint(a).astype(np.int64)
        cs, rs, n = r(cs), r(rs), r(n)
        d = dict(n_cal=int(n[0]), n_S2=int(n[1]), n_Ts=int(n[2]), n_W2=int(n[3]), n_gate=int(n[4]),
                 cs_cal=cs[:, 0], rs_cal=rs[:, 0], cs_S2=cs[:, 1], rs_S2=rs[:, 1], cs_Ts=cs[:, 2], rs_Ts=rs[:, 2], cs_W2=cs[:, 3], rs_W2=rs[:, 3], g=cs[:, 4],
                 csR=cs[:, 6], rsC=rs[:, 5])
        d["dec"] = d["cs_cal"] + d["rs_cal"]
        d["contra"] = d["csR"] + d["rsC"]
        d["un_cal"] = d["n_cal"] - d["dec"]
        d["un_S2"] = d["n_S2"] - d["cs_S2"] - d["rs_S2"]
        d["dec_Ts"] = d["cs_Ts"] + d["rs_Ts"]; d["dec_W2"] = d["cs_W2"] + d["rs_W2"]
        return d
    def feasibility(st):
        ok = {}
        ok[1] = f1.copy()
        e, s, a, c = thr["e_cal"], thr["s_cal"], thr["a_cal"], thr["c_min"]
        ok[2] = (st["dec"] > 0) & (st["contra"] * e.denominator <= e.numerator * st["dec"])
        ok[3] = np.full(T, st["n_S2"] > 0) & (st["un_S2"] * s.denominator <= s.numerator * st["n_S2"])
        ok[4] = (st["g"] == 0)
        ok[5] = np.full(T, st["n_cal"] > 0) & (st["un_cal"] * a.denominator <= a.numerator * st["n_cal"])
        ok[6] = (np.full(T, st["n_W2"] > 0) & (st["dec_W2"] * c.denominator >= c.numerator * st["n_W2"])
                 & np.full(T, st["n_Ts"] > 0) & (st["dec_Ts"] * c.denominator >= c.numerator * st["n_Ts"]))
        return ok
    def select(st):
        ok = feasibility(st)
        feas = ok[1] & ok[2] & ok[3] & ok[4] & ok[5] & ok[6]
        idx = np.nonzero(feas)[0]
        if len(idx) == 0:
            return dict(ok=ok, feas=idx, J={}, Jmin=None, tol=idx, star=None)
        J = {}
        for t in idx:
            J[int(t)] = (Fraction(int(st["rs_Ts"][t]), st["n_Ts"]) + Fraction(int(st["rs_W2"][t]), st["n_W2"])) / 2
        Jmin = min(J.values())
        tol = [int(t) for t in idx if J[int(t)] <= Jmin + thr["delta_tol"]]
        tol_sorted = sorted(tol, key=lambda t: (int(cr_int[t]), int(sizeA[t]), int(wA[t]), int(vA[t]), int(maskA[t]), int(qcA[t]), int(qrA[t])))
        return dict(ok=ok, feas=idx, J=J, Jmin=Jmin, tol=np.array(tol_sorted, dtype=np.int64), star=int(tol_sorted[0]))
    wfull = np.ones(NR)
    st = stats_of(wfull)
    sel = select(st)
    alone = {k: int(sel["ok"][k].sum()) for k in range(1, 7)}
    cum, run = {}, np.ones(T, dtype=bool)
    for k in range(1, 7):
        run = run & sel["ok"][k]; cum[k] = int(run.sum())
    firstfail = np.zeros(T, dtype=np.int64)
    for t in range(T):
        for k in range(1, 7):
            if not sel["ok"][k][t]:
                firstfail[t] = k; break
    res["theta_space"].update(alone={"F%d" % k: v for k, v in alone.items()}, cumulative_after={"F%d" % k: v for k, v in cum.items()},
                              first_fail_counts={str(k): int((firstfail == k).sum()) for k in range(0, 7)}, F_size=int(len(sel["feas"])))
    res["cost_cap"] = rtext(cap)
    log("C4: F1..F6 alone", alone, "cumulative", cum, "|F|", len(sel["feas"]))
    min_set = q["min_set"]
    sets_small = st["n_Ts"] < min_set or st["n_W2"] < min_set or st["n_S2"] < min_set
    outcome = "sets too small" if sets_small else ("no feasible rule" if len(sel["feas"]) == 0 else "selected")
    res["outcome"] = outcome
    def tstats(t, stt=None):
        stt = st if stt is None else stt
        return dict(index=int(t) + 1, theta=theta_text(thetas[t]), C_R=rtext(CRs[t]), e=(rtext(Fraction(int(stt["contra"][t]), int(stt["dec"][t]))) if stt["dec"][t] > 0 else None),
                    decided=int(stt["dec"][t]), contradicting=int(stt["contra"][t]), u_cal=rtext(Fraction(int(stt["un_cal"][t]), stt["n_cal"])),
                    u_S2=rtext(Fraction(int(stt["un_S2"][t]), stt["n_S2"])), FR_Ts=rtext(Fraction(int(stt["rs_Ts"][t]), stt["n_Ts"])), FR_W2=rtext(Fraction(int(stt["rs_W2"][t]), stt["n_W2"])),
                    dec_Ts=rtext(Fraction(int(stt["dec_Ts"][t]), stt["n_Ts"])), dec_W2=rtext(Fraction(int(stt["dec_W2"][t]), stt["n_W2"])), g=int(stt["g"][t]))
    ctx = dict(q=q, p=p, orbits=orbits, ser=ser, win=win, thetas=thetas, T=T, costs=costs, cap=cap, CRs=CRs, OUT=OUT, V=V, rows_ids=rows_ids, rowpos=rowpos, setd=setd, rho=rho,
               st=st, sel=sel, thr=thr, B=B, kP=kP, tstats=tstats, stats_of=stats_of, select=select, mcal=mcal, mTs=mTs, mW2=mW2, mS2=mS2, mgate=mgate,
               cal_ids=cal_ids, cr_int=cr_int, sizeA=sizeA, wA=wA, vA=vA, maskA=maskA, qcA=qcA, qrA=qrA, NR=NR, nrows_cal=nrows_cal, res=res, log=log, N1=N1, NC=NC)
    res["selection"] = dict(feasible_count=int(len(sel["feas"])))
    if len(sel["feas"]):
        feas_rows = []
        tolset = set(int(x) for x in sel["tol"])
        for t in sel["feas"]:
            t = int(t)
            Jt = sel["J"][t]
            feas_rows.append(dict(index=t + 1, theta=theta_text(thetas[t]), C_R=rtext(CRs[t]), FR_W2=rtext(Fraction(int(st["rs_W2"][t]), st["n_W2"])), FR_Ts=rtext(Fraction(int(st["rs_Ts"][t]), st["n_Ts"])),
                                  J=rtext(Jt), in_tol=bool(t in tolset), star=bool(sel["star"] is not None and t == sel["star"])))
        res["selection"]["F_rows"] = feas_rows
        res["selection"]["J_min"] = rtext(sel["Jmin"])
        res["selection"]["tol_size"] = int(len(sel["tol"]))
        res["selection"]["tol_top20"] = [dict(tstats(t), rank=n + 1) for n, t in enumerate(sel["tol"][:20])]
    log("outcome", outcome, "J_min", res["selection"].get("J_min"), "|Theta_tol|", res["selection"].get("tol_size"))
    if outcome != "selected":
        return res, ctx
    star = sel["star"]
    th = thetas[star]
    ctx["star"] = star
    res["selection"]["star"] = dict(tstats(star), V=th["V"], qC=th["qC"], qR=th["qR"], w=th["w"], v=th["v"], J=rtext(sel["J"][star]))
    log("theta*", theta_text(th), "index", star + 1)
    c7_and_pp(ctx, skip_pp)
    return res, ctx

# ====================================================================================== C7, PP for HB (E56) and HC (E57), R_ab
def fl(q):
    return float(Fraction(q))

OUTNAME = {0: "unresolved", 1: "chaotic-supported", 2: "regular-supported"}

def c7_and_pp(ctx, skip_pp):
    res, q, p, costs, thetas, star, V, OUT, rowpos, setd, orbits, B, kP, thr, log = (ctx[k] for k in ("res", "q", "p", "costs", "thetas", "star", "V", "OUT", "rowpos", "setd", "orbits", "B", "kP", "thr", "log"))
    th = thetas[star]
    Vstar = th["V"]
    CR, matched, le_read, le_unm = costs.matched(Vstar)
    mres = dict(C_R=rtext(CR), le_read=le_read, le_cost_unmatched=bool(le_unm), cost_cap=rtext(ctx["cap"]), indicators={})
    for X in costs.I6 + ["LE"]:
        m = matched[X]
        mres["indicators"][X] = dict(b=m["b"], capped=m["capped"], round_down=m["down"], round_down_none=m["down"] is None, cost=(rtext(m["b"] * costs.cx(X)) if m["b"] is not None else None))
    # over the feasible elements: how many do not carry the LE flag; any non-LE capped comparator
    n_nolf, n_capped_nonLE = 0, 0
    for t in ctx["sel"]["feas"]:
        CRt, mt, ler, unm = costs.matched(thetas[int(t)]["V"])
        n_nolf += (not unm)
        n_capped_nonLE += any(mt[X]["capped"] for X in costs.I6)
    mres["feasible_elements_without_le_cost_unmatched_flag"] = n_nolf
    mres["feasible_elements_with_a_capped_non_le_comparator"] = n_capped_nonLE
    # Double64 costs
    sp = sorted({costs.spec[X] for X in Vstar})
    flt = costs.bP * sum(costs.c[s] for s in sp)
    pair = costs.bP * sum(costs.c[s] + costs.d64[s] * costs.c[s] for s in sp)
    mres["double64"] = dict(pair=rtext(pair), float_votes=rtext(flt), ratio=rtext(flt / pair))
    res["matched"] = mres
    bidx = lambda X: B.index(le_read if X == "LE" else matched[X]["b"])
    pos_of = lambda ids: [rowpos[i] for i in ids]
    out_star = OUT[star]
    Ts_ids, W2_ids, RW_ids = setd["T_s"], setd["W2"], setd["R_W"]
    nT, nW, nRW = len(Ts_ids), len(W2_ids), len(RW_ids)
    # ---------------------------------------------------------------- C7
    c7 = dict(by_source={}, hb={})
    for A, ids in (("T_s", Ts_ids), ("W2", W2_ids)):
        for src in ("P", "G"):
            sub = [i for i in ids if orbits[i]["src"] == src]
            ps = pos_of(sub)
            o = out_star[ps] if sub else np.array([], dtype=np.int8)
            c7["by_source"]["%s,%s" % (A, "panel" if src == "P" else "grid")] = dict(
                n=len(sub), FR=rtext(Fraction(int((o == 2).sum()), len(sub))) if sub else None, decided=rtext(Fraction(int((o != 0).sum()), len(sub))) if sub else None)
    reading("C7 risk-coverage family: theta* with q_C and q_R each varied by one step inside their ranges, as the product grid (up to 9 points, theta* included); the cross (one at a time) is reported as well")
    for A, ids, comps in (("T_s", Ts_ids, costs.I6), ("W2", W2_ids, ["LE"] + costs.I6)):
        ps = pos_of(ids)
        o = out_star[ps]
        rs = (o == 2); dec = (o != 0)
        rows = []
        for X in comps:
            xi = I7.index(X)
            m = matched[X]
            row = dict(X=X, b_ceiling=m["b"], b_read=(le_read if X == "LE" else m["b"]), capped=m["capped"], cost_unmatched=bool(le_unm) if X == "LE" else False)
            if X != "LE" and m["capped"]:
                rows.append(row); continue
            k = bidx(X)
            xr = (V[ps, xi, k] == RG)
            rb = le_read if X == "LE" else m["b"]
            n10 = int((rs & ~xr).sum()); n01 = int((~rs & xr).sum()); n11 = int((rs & xr).sum())
            mm = n10 + n01
            nrule, nx = int(rs.sum()), int(xr.sum())
            row.update(cost=rtext(rb * costs.cx(X)), FR_X=rtext(Fraction(nx, len(ids))), FR_rule=rtext(Fraction(nrule, len(ids))), n_FR_X=nx, n_FR_rule=nrule, n10=n10, n01=n01, n11=n11,
                       p=rtext(mcnemar_tail(mm, n01)) if mm > 0 else "1/1", margin=bool(nx > 0 and nrule * thr["rho_hb"].denominator <= thr["rho_hb"].numerator * nx),
                       coverage=rtext(Fraction(int(dec.sum()), len(ids))))
            row["p_le_alpha"] = bool(mcnemar_tail(mm, n01) <= thr["alpha_hb"]) if mm > 0 else False
            rows.append(row)
        c7["hb"][A] = rows
    # risk-coverage
    stt = ctx["st"]
    calpos = list(range(ctx["nrows_cal"]))
    rhoC = np.array([ctx["rho"][ctx["rows_ids"][i]] == "C" for i in calpos]); rhoR = ~rhoC
    def rc_point(outv):
        dec = outv != 0
        nd = int(dec.sum())
        contra = int(((outv == 1) & rhoR).sum() + ((outv == 2) & rhoC).sum())
        return dict(decided=nd, contradicting=contra, decided_share=rtext(Fraction(nd, len(calpos))), e=(rtext(Fraction(contra, nd)) if nd else None))
    fam, cross = {}, {}
    sz = th["size"]
    for dc in (-1, 0, 1):
        for dr in (-1, 0, 1):
            qc, qr = th["qC"] + dc, th["qR"] + dr
            if 1 <= qc <= sz and sz // 2 + 1 <= qr <= sz:
                idx = next(t["index"] - 1 for t in thetas if (t["mask"], t["qC"], t["qR"], t["w"], t["v"]) == (th["mask"], qc, qr, th["w"], th["v"]))
                pt = rc_point(OUT[idx][:len(calpos)]); pt["theta"] = theta_text(thetas[idx]); pt["index"] = idx + 1
                fam["%d,%d" % (qc, qr)] = pt
                if dc == 0 or dr == 0:
                    cross["%d,%d" % (qc, qr)] = pt
    single = {}
    for X in I7:
        xi = I7.index(X)
        v = V[:len(calpos), xi, bidx(X) if (X == "LE" or not matched[X]["capped"]) else kP]
        outv = np.where(v == CH, 1, np.where(v == RG, 2, 0)).astype(np.int8)
        pt = rc_point(outv); pt["b"] = le_read if X == "LE" else matched[X]["b"]
        single[X] = pt
    c7["risk_coverage"] = dict(family_grid=fam, family_cross=cross, single=single)
    # K > 0 report (P48)
    posK = [i for i in calpos if (orbits[ctx["rows_ids"][i]]["K"] or 0) > 0]
    o = out_star[posK]
    rhoCK = rhoC[posK]
    ndK = int((o != 0).sum()); contraK = int(((o == 1) & ~rhoCK).sum() + ((o == 2) & rhoCK).sum())
    c7["K_positive"] = dict(n=len(posK), decided=ndK, e=rtext(Fraction(contraK, ndK)) if ndK else None, u_cal=rtext(Fraction(len(posK) - ndK, len(posK))))
    c7["theta_star_stats"] = ctx["tstats"](star)
    res["c7"] = c7
    # ---------------------------------------------------------------- R_ab on the Cal orbits that carry Double64 items (P63)
    hd = {k[3] for k in ctx["ser"] if k[2] == "Double64" and k[1] in ("L1", "tg")}
    rab = []
    hor_d = sorted(hd)[0] if hd else None
    specs_need = {costs.spec[X] for X in Vstar}
    for i in ctx["cal_ids"]:
        have = all(((i, "tg", "Double64", hor_d) in ctx["ser"]) if s == "tg" else ((i, s, "Double64", ctx["NC"] + (1 if s == "pa1" else 0)) in ctx["win"]) for s in specs_need)
        if not have:
            continue
        vd = verdicts_of(i, "Double64", ctx["ser"], ctx["win"], p, ctx["N1"], ctx["NC"], B, hor_tg=hor_d)
        outc = int(out_star[rowpos[i]])
        over = False
        for X in Vstar:
            f64 = int(V[rowpos[i], I7.index(X), kP]); d64 = vd[X][kP]
            if f64 != 0 and d64 != 0 and f64 != d64:
                over = True
        rab.append(dict(id=i, R_a=OUTNAME[outc], R_ab=("unresolved" if over else OUTNAME[outc]), override=over))
    res["r_ab"] = dict(rows=rab, n=len(rab), overrides=sum(1 for r in rab if r["override"]))
    log("C7 done; R_ab orbits", len(rab))
    if skip_pp:
        return
    pp_hb(ctx, matched, le_read)
    pp_hc(ctx)

def pp_hb(ctx, matched, le_read):
    res, q, thetas, star, V, OUT, rowpos, setd, orbits, B, kP, thr, log, costs = (ctx[k] for k in ("res", "q", "thetas", "star", "V", "OUT", "rowpos", "setd", "orbits", "B", "kP", "thr", "log", "costs"))
    th = thetas[star]
    st, stats_of, select = ctx["st"], ctx["stats_of"], ctx["select"]
    test_seed = q["test_seed"]
    nrows_cal = ctx["nrows_cal"]
    # ---------------------------------------------------------------- a. optimism
    log("PP-HB: optimism over", q["b_opt"], "resamples")
    sk, used = [], []
    thetas_b = []
    opt_sum = {"T_s": 0.0, "W2": 0.0}
    nvalid = 0
    for b in range(1, q["b_opt"] + 1):
        g = Xoshiro.from_seed(stream_seed(test_seed, "optimism", str(b)))
        draws = [g.randint(nrows_cal) for _ in range(nrows_cal)]
        w = np.zeros(ctx["NR"])
        w += np.concatenate([np.bincount(np.array(draws) - 1, minlength=nrows_cal).astype(np.float64), np.zeros(ctx["NR"] - nrows_cal)])
        stb = stats_of(w)
        sb = select(stb)
        if sb["star"] is None:
            sk.append(b); thetas_b.append(None); continue
        tb = sb["star"]
        thetas_b.append(tb + 1)
        nvalid += 1
        for A, rsk, nk in (("T_s", "rs_Ts", "n_Ts"), ("W2", "rs_W2", "n_W2")):
            cal_fr = float(Fraction(int(st[rsk][tb]), st[nk]))
            res_fr = float(Fraction(int(stb[rsk][tb]), stb[nk]))
            opt_sum[A] += (cal_fr - res_fr)
    optim = {A: (opt_sum[A] / nvalid if nvalid else None) for A in opt_sum}
    res["optimism"] = dict(used=nvalid, skipped=len(sk), skipped_list=sk, theta_b=thetas_b, optimism_T_s=optim["T_s"], optimism_W2=optim["W2"])
    log("optimism: used", nvalid, "skipped", len(sk), optim)
    # ---------------------------------------------------------------- b. strata and rates
    r_ts = float(Fraction(int(st["rs_Ts"][star]), st["n_Ts"]))
    r_w2 = float(Fraction(int(st["rs_W2"][star]), st["n_W2"]))
    r_corr_T = r_ts + max(optim["T_s"], 0.0)
    r_corr_W = r_w2 + max(optim["W2"], 0.0)
    nT, nW = st["n_Ts"], st["n_W2"]
    nRW = nW - nT
    r_corr_RW = (r_corr_W * nW - r_corr_T * nT) / nRW
    r_corr_RW = min(max(r_corr_RW, 0.0), 1.0)
    floor, factor = float(thr["rule_rate_floor"]), float(thr["rule_rate_factor"])
    prule = lambda r: r if r >= floor else max(factor * r, floor)
    p_rule = {"T": prule(r_corr_T), "RW": prule(r_corr_RW)}
    res["pp_hb"] = dict(r_Ts=r_ts, r_W2=r_w2, r_corr_T_s=r_corr_T, r_corr_W2=r_corr_W, r_corr_R_W=r_corr_RW, p_rule_T_s=p_rule["T"], p_rule_R_W=p_rule["RW"])
    out_star = OUT[star]
    ids = {"T": setd["T_s"], "W2": setd["W2"]}
    posT = [rowpos[i] for i in ids["T"]]; posW = [rowpos[i] for i in ids["W2"]]
    def strata_rates(sub=None):
        """rates of C7/E56 b on Cal (or on one source): comparator rates, shared shares, decided shares, strata sizes"""
        selT = [i for i in ids["T"] if sub is None or orbits[i]["src"] == sub]
        selW = [i for i in ids["W2"] if sub is None or orbits[i]["src"] == sub]
        pT = [rowpos[i] for i in selT]; pW = [rowpos[i] for i in selW]
        nt, nw = len(selT), len(selW)
        nrw = nw - nt
        rsT = out_star[pT] == 2; rsW = out_star[pW] == 2
        dT = int((out_star[pT] != 0).sum()); dW = int((out_star[pW] != 0).sum())
        rfrT, rfrW = int(rsT.sum()), int(rsW.sum())
        d = {"T": (Fraction(dT, nt) if nt else None), "RW": (Fraction(dW - dT, nrw) if nrw else None)}
        pX, sX = {"T": {}, "RW": {}}, {"T": {}, "RW": {}}
        for X in I7:
            xi = I7.index(X)
            k = B.index(le_read if X == "LE" else matched[X]["b"])
            xT = V[pT, xi, k] == RG; xW = V[pW, xi, k] == RG
            cT, cW = int(xT.sum()), int(xW.sum())
            shT, shW = int((rsT & xT).sum()), int((rsW & xW).sum())
            pX["T"][X] = Fraction(cT, nt) if nt else None
            pX["RW"][X] = Fraction(cW - cT, nrw) if nrw else None
            sX["T"][X] = Fraction(shT, rfrT) if rfrT else fr(q["shared_default"])
            sX["RW"][X] = Fraction(shW - shT, rfrW - rfrT) if (rfrW - rfrT) else fr(q["shared_default"])
        return dict(n={"T": nt, "RW": nrw}, d=d, pX=pX, sX=sX, rule_fr={"T": rfrT, "RW": rfrW - rfrT})
    rates_pool = strata_rates()
    res["pp_hb"]["rates"] = {"pooled": rates_to_json(rates_pool), "panel": rates_to_json(strata_rates("P")), "grid": rates_to_json(strata_rates("G"))}
    # ---------------------------------------------------------------- the model of one orbit (c)
    def model_from(rates, prule_over=None):
        m = {}
        for S in ("T", "RW"):
            pr = p_rule[S] if prule_over is None else prule_over[S]
            d = fl(rates["d"][S])
            dprob = min(max((d - pr) / (1 - pr), 0.0), 1.0) if pr < 1 else 1.0
            sx = [fl(rates["sX"][S][X]) for X in I7]
            px = [fl(rates["pX"][S][X]) for X in I7]
            qx = [(max(a - b * pr, 0.0) / (1 - pr) if pr < 1 else 0.0) for a, b in zip(px, sx)]
            m[S] = dict(p=pr, dprob=dprob, s=sx, q=qx)
        return m
    model = model_from(rates_pool)
    res["pp_hb"]["model"] = {S: dict(p=model[S]["p"], dprob=model[S]["dprob"], s=model[S]["s"], q=model[S]["q"]) for S in model}
    # ---------------------------------------------------------------- d. draw path
    start, second, rmax = fr(q["rho_path_start"]), fr(q["rho_path_second"]), fr(q["r_max"])
    step = second - start
    npts = int((rmax - start) / step) + 1
    yT, yW = fr(q["y_t"]), fr(q["y_w"])
    path = []
    for j in range(npts):
        rho_r = start + j * step
        nt = round_half_away(rho_r * yT); nrw = round_half_away(rho_r * (yW - yT))
        allowed = nt >= q["n_t_min"] and (nt + nrw) >= q["n_w_min"]
        path.append(dict(j=j, rho=rtext(rho_r), nT=nt, nRW=nrw, nW=nt + nrw, allowed=allowed))
    crit = mcnemar_critical_table(max(pt["nW"] for pt in path) + 1, thr["alpha_hb"])
    crit_arr = np.array(crit, dtype=np.int64)
    runs = q["pp_runs"]
    log("PP-HB: draw path", npts, "points,", sum(1 for pt in path if pt["allowed"]), "allowed; runs", runs)
    for pt in path:
        if not pt["allowed"]:
            pt.update(successes=None, runs=None, power=None, se=None)
            continue
        t1 = time.time()
        r = sim_hb(model, pt["nT"], pt["nRW"], pt["j"], runs, crit_arr, thr)
        pw = r["joint"].sum() / runs
        pt.update(successes=int(r["joint"].sum()), runs=runs, power=float(pw), se=math.sqrt(pw * (1 - pw) / runs))
        log("  point", pt["j"], "rho", pt["rho"], "nT", pt["nT"], "nW", pt["nW"], "successes", pt["successes"], "(%.1f s)" % (time.time() - t1))
    res["pp_hb"]["path"] = path
    pHB = thr["p_hb"]
    star_pt = next((pt for pt in path if pt["allowed"] and Fraction(pt["successes"], runs) >= pHB), None)
    targets = {}
    if not any(pt["allowed"] for pt in path):
        targets = dict(outcome="no allowed point")
    elif star_pt is not None:
        rho_star = fr(star_pt["rho"].split("/")[0]) / fr(star_pt["rho"].split("/")[1])
        targets = dict(outcome="ok", cap=False, j=star_pt["j"], rho_star=star_pt["rho"], n_T=ceil_step(yT * rho_star, q["ceil_step"]), n_W=ceil_step(yW * rho_star, q["ceil_step"]))
    else:
        last = path[-1]
        targets = dict(outcome="cap", cap=True, n_T=int(yT * rmax), n_W=int(yW * rmax), power_at_cap=last["power"], successes_at_cap=last["successes"], j=last["j"], rho_star=None, rho_at_cap=last["rho"])
    res["pp_hb"]["targets"] = targets
    log("PP-HB targets", targets)
    # ---------------------------------------------------------------- f. type-I check
    if targets["outcome"] != "no allowed point":
        pt = path[targets["j"]]
        t1r = {}
        for X in I7:
            xi = I7.index(X)
            over = {"T": fl(rates_pool["pX"]["T"][X]), "RW": fl(rates_pool["pX"]["RW"][X])}
            m2 = model_from(rates_pool, prule_over=over)
            r = sim_hb(m2, pt["nT"], pt["nRW"], pt["j"], runs, crit_arr, thr, single_for=xi)
            jr = r["joint"].sum() / runs
            t1r[X] = dict(joint_successes=int(r["joint"].sum()), joint=float(jr), single_successes=int(r["single"].sum()), single=float(r["single"].sum() / runs), se=math.sqrt(jr * (1 - jr) / runs))
            log("  type-I", X, t1r[X])
        res["pp_hb"]["type1"] = dict(point=pt["j"], per_comparator=t1r, max_joint=max(v["joint"] for v in t1r.values()))
        res["pp_hb"]["type1"]["reading"] = "single = the comparison(s) of X alone (on T_s if X in I_6, and on W2) both pass, coverage and other comparators ignored; same streams as the main point (common random numbers)"

def rates_to_json(r):
    return dict(n=r["n"], d={S: (rtext(v) if v is not None else None) for S, v in r["d"].items()},
                pX={S: {X: (rtext(v) if v is not None else None) for X, v in d.items()} for S, d in r["pX"].items()},
                sX={S: {X: rtext(v) for X, v in d.items()} for S, d in r["sX"].items()}, rule_false_regular=r["rule_fr"])

def sim_hb(model, nT, nRW, j, runs, crit_arr, thr, single_for=None, chunk=64):
    """the full HB decision for `runs` runs; run i uses the stream (pp-hb, "j,i"), i = 1..runs; nine uniforms per orbit in the order rule, seven comparators, decided"""
    seeds = np.array([int.from_bytes(hashlib.sha256(("chaoslab-analysis-1|%s|p3|pp-hb|%d,%d" % (SIM_TEST_SEED[0], j, i)).encode()).digest()[:8], "big") for i in range(1, runs + 1)], dtype=np.uint64)
    xo = XoV(seeds)
    cnt = {S: dict(n11=np.zeros((7, runs), dtype=np.int32), n10=np.zeros((7, runs), dtype=np.int32), n01=np.zeros((7, runs), dtype=np.int32), dec=np.zeros(runs, dtype=np.int32)) for S in ("T", "RW")}
    for S, n in (("T", nT), ("RW", nRW)):
        m = model[S]
        sx = np.array(m["s"])[:, None]; qx = np.array(m["q"])[:, None]
        done = 0
        while done < n:
            c = min(chunk, n - done)
            U = np.empty((c, 9, runs))
            for a in range(c):
                for b in range(9):
                    U[a, b] = xo.uniform()
            rfr = U[:, 0, :] < m["p"]                                   # (c, runs)
            xfr = np.where(rfr[:, None, :], U[:, 1:8, :] < sx[None], U[:, 1:8, :] < qx[None])   # (c, 7, runs)
            decd = rfr | (U[:, 8, :] < m["dprob"])
            R = rfr[:, None, :]
            cnt[S]["n11"] += (R & xfr).sum(0, dtype=np.int32)
            cnt[S]["n10"] += (R & ~xfr).sum(0, dtype=np.int32)
            cnt[S]["n01"] += (~R & xfr).sum(0, dtype=np.int32)
            cnt[S]["dec"] += decd.sum(0, dtype=np.int32)
            done += c
    rn, rd = thr["rho_hb"].numerator, thr["rho_hb"].denominator
    cm, cd = thr["c_min"].numerator, thr["c_min"].denominator
    def comp_ok(n11, n10, n01):
        mm = n10 + n01
        a = n10 + n11; c = n01 + n11
        return (mm >= 1) & (n01 >= crit_arr[np.minimum(mm, len(crit_arr) - 1)]) & (c >= 1) & (a * rd <= rn * c)
    allok = np.ones(runs, dtype=bool)
    singles = np.ones(runs, dtype=bool)
    for xi in range(7):
        okT = comp_ok(cnt["T"]["n11"][xi], cnt["T"]["n10"][xi], cnt["T"]["n01"][xi]) if xi >= 1 else None   # LE is not a comparator on T_s
        s11 = cnt["T"]["n11"][xi] + cnt["RW"]["n11"][xi]; s10 = cnt["T"]["n10"][xi] + cnt["RW"]["n10"][xi]; s01 = cnt["T"]["n01"][xi] + cnt["RW"]["n01"][xi]
        okW = comp_ok(s11, s10, s01)
        if okT is not None:
            allok &= okT
        allok &= okW
        if single_for is not None and xi == single_for:
            singles = (okT & okW) if okT is not None else okW
    covT = cnt["T"]["dec"] * cd >= cm * nT
    covW = (cnt["T"]["dec"] + cnt["RW"]["dec"]) * cd >= cm * (nT + nRW)
    joint = allok & covT & covW
    return dict(joint=joint, single=singles)

SIM_TEST_SEED = [None]

def pp_hc(ctx):
    res, q, thetas, star, V, setd, orbits, B, kP, thr, log, rowpos = (ctx[k] for k in ("res", "q", "thetas", "star", "V", "setd", "orbits", "B", "kP", "thr", "log", "rowpos"))
    th = thetas[star]
    Vs = [I7.index(X) for X in th["V"]]
    # E57: Cal orbits with rho in {C,R} and K > 0
    fails, nonf = [], []
    byfail = {"P": [0, 0], "G": [0, 0]}
    sc_f, sc_n = [], []
    sc_by = {"P": ([], []), "G": ([], [])}
    nK = 0
    for i in ctx["cal_ids"]:
        o = orbits[i]
        if not ((o["K"] or 0) > 0):
            continue
        nK += 1
        v = V[rowpos[i], Vs, kP]
        nC = int((v == CH).sum()); nR = int((v == RG).sum())
        label = "c" if nC > nR else ("r" if nR > nC else "tie")
        rho = ctx["rho"][i]
        f = (label == "tie") or ((label == "c") != (rho == "C"))
        s = Fraction(0) if nC + nR == 0 else 1 - Fraction(abs(nC - nR), nC + nR)
        (sc_f if f else sc_n).append(s)
        (sc_by[o["src"]][0] if f else sc_by[o["src"]][1]).append(s)
    A_cal = mw_auc(sc_f, sc_n)
    hc = dict(n_K_positive=nK, n_fail=len(sc_f), n_nonfail=len(sc_n), by_source={"panel": dict(n_fail=len(sc_by["P"][0]), n_nonfail=len(sc_by["P"][1]), auc=(rtext(mw_auc(*sc_by["P"])) if mw_auc(*sc_by["P"]) is not None else None)),
                                                                                  "grid": dict(n_fail=len(sc_by["G"][0]), n_nonfail=len(sc_by["G"][1]), auc=(rtext(mw_auc(*sc_by["G"])) if mw_auc(*sc_by["G"]) is not None else None))})
    res["pp_hc"] = hc
    if A_cal is None:
        hc["outcome"] = "AUC undefined (no failure or no non-failure on Cal)"
        return
    hc["AUC_I_cal"] = rtext(A_cal); hc["AUC_I_cal_float"] = float(A_cal)
    A_used = Fraction(1, 2) + thr["auc_shrink"] * (A_cal - Fraction(1, 2))
    A_P = thr["auc_p"]
    hc["AUC_I_used"] = rtext(A_used); hc["AUC_P"] = rtext(A_P); hc["true_difference"] = float(A_used - A_P)
    dI = math.sqrt(2) * NormalDist().inv_cdf(float(A_used)); dP = math.sqrt(2) * NormalDist().inv_cdf(float(A_P))
    hc["d_I"] = dI; hc["d_P"] = dP
    log("PP-HC: AUC_I_cal", hc["AUC_I_cal"], "fails", len(sc_f), "nonfails", len(sc_n), "AUC_I_used", float(A_used), "dI", dI, "dP", dP)
    nf0, nf1 = q["nf_start"], q["nf_second"]
    stp = nf1 - nf0
    nn_per = q["nonfail_per_fail"]
    dexp = q["d_expected_total"] // (1 + nn_per); dwor = q["d_worst_total"] // (1 + nn_per)
    cap_nf = nf0 + math.ceil(max(dwor - nf0, 0) / stp) * stp
    lst = list(range(nf0, cap_nf + 1, stp))
    hc["deliverable_expected"] = dexp; hc["deliverable_worst"] = dwor; hc["nf_search_cap"] = cap_nf
    runs = q["pp_runs"]
    zz = float(thr_float(q["z_one_sided"])); dcut = thr["delta_c"]  # exact Fraction (E57 point condition)
    rows = []
    for nf in lst + ([dexp] if dexp not in lst else []):
        t1 = time.time()
        succ = sim_hc(nf, nn_per * nf, dI, dP, runs, dcut, zz)
        pw = succ / runs
        rows.append(dict(nf=nf, successes=succ, runs=runs, power=pw, se=math.sqrt(pw * (1 - pw) / runs), on_list=nf in lst))
        log("  HC n_f", nf, "power", pw, "(%.1f s)" % (time.time() - t1))
    hc["table"] = rows
    first = next((r for r in rows if r["on_list"] and Fraction(r["successes"], runs) >= thr["p_hc"]), None)
    hc["target"] = dict(n_f=first["nf"], power=first["power"], se=first["se"]) if first is not None else "not reached"
    hc["power_at_deliverable_expected"] = next(r["power"] for r in rows if r["nf"] == dexp)
    hc["underpowered"] = bool(first is None or first["nf"] > dexp)
    if first is None:
        hc["power_at_cap"] = rows[len(lst) - 1]["power"]

def thr_float(x):
    return x

def sim_hc(nf, nn, dI, dP, runs, dcut, zz, chunk=2500):
    succ = 0
    for r0 in range(1, runs + 1, chunk):
        r1 = min(runs, r0 + chunk - 1)
        seeds = np.array([int.from_bytes(hashlib.sha256(("chaoslab-analysis-1|%s|p3|pp-hc|%d,%d" % (SIM_TEST_SEED[0], nf, i)).encode()).digest()[:8], "big") for i in range(r0, r1 + 1)], dtype=np.uint64)
        xo = XoV(seeds)
        R = len(seeds)
        nsteps = 2 * (nf + nn)
        U = np.empty((nsteps, R))
        for s in range(nsteps):
            U[s] = xo.uniform()
        Z = ppnd16(U.reshape(-1)).reshape(nsteps, R)
        # order: I failures, I non-failures, P failures, P non-failures
        sI = np.concatenate([Z[:nf] + dI, Z[nf:nf + nn]], axis=0).T          # (R, nf + nn); failures first
        sP = np.concatenate([Z[nf + nn:2 * nf + nn] + dP, Z[2 * nf + nn:]], axis=0).T
        def ucount(S):
            ranks = np.argsort(np.argsort(S, axis=1, kind="stable"), axis=1, kind="stable") + 1
            return ranks[:, :nf].sum(1) - nf * (nf + 1) // 2   # Mann-Whitney U (no ties for continuous scores)
        UI, UP = ucount(sI), ucount(sP)
        AI, AP = UI / (nf * nn), UP / (nf * nn)
        # exact point condition: doubled counts 2U share the denominator 2 nf nn; Delta >= delta_c  <=>  (2UI - 2UP) * den(dc) >= num(dc) * 2 nf nn (int64 numpy, then python-exact bound)
        dc = Fraction(dcut)
        pt = ((UI.astype(object) - UP.astype(object)) * dc.denominator >= dc.numerator * nf * nn)
        def var(A):
            Q1 = A / (2 - A); Q2 = 2 * A * A / (1 + A)
            return (A * (1 - A) + (nf - 1) * (Q1 - A * A) + (nn - 1) * (Q2 - A * A)) / (nf * nn)
        D = AI - AP
        SE = np.sqrt(var(AI) + var(AP))
        succ += int(((pt.astype(bool)) & (D - zz * SE > 0)).sum())
    return succ

def to_jsonable(o):
    if isinstance(o, Fraction):
        return rtext(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, dict):
        return {str(k): to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(x) for x in o]
    return o

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?")
    ap.add_argument("result", nargs="?")
    ap.add_argument("--p3", default=P3PAR_DEFAULT)
    ap.add_argument("--p2", default=P2PAR_DEFAULT)
    ap.add_argument("--skip-pp", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest or a.folder is None:
        sys.exit(0 if selftest() else 1)
    ok = selftest()
    if not ok:
        sys.exit("self test failed")
    q0 = tomllib.load(open(a.p3, "rb"))
    SIM_TEST_SEED[0] = q0["test_seed"]["value"]
    res, ctx = main_analysis(a.folder, a.p3, a.p2, a.skip_pp)
    os.makedirs(os.path.dirname(os.path.abspath(a.result)), exist_ok=True)
    with open(a.result, "w", encoding="utf-8") as f:
        json.dump(to_jsonable(res), f, indent=1, allow_nan=True)
    say("written", a.result)
