"""O23 (Amendment D94): independent reimplementation of the panel draw in rounds, Python standard library only.

    python -I verify/experiments-d94-oracle.py --panel P.toml --pools P.toml.pools.toml [--db EXPORT.toml] [--generator generator.txt]
    python -I verify/experiments-d94-oracle.py --generator-only generator.txt

Written from the text of Amendment D94 (E32 to E38) and the base sections E12 to E16; reads the pools from the sidecar and the inputs
from the panel manifest. Rebuilds the points of every K, the draws of every round, ids, seeds, D and D_B for R = 1 .. r_max.
Output: one line per check `<id> | quantity | value | bound | PASS/FAIL/INFO`, a final count; exit 1 on any FAIL.

Generator (spec Oracle, O23): state = (seed << 1) | 1 mod 2^128, multiplier 0x45a31efc5a35d971261fd0407a968add, each draw multiplies the
state mod 2^128 and returns the high 64 bits; uniform = (low 52 bits of the word) * 2^-52 (confirmed against the installed library's numbers).
Per-orbit seed (E16 points at the tangent-core spec, not read here): the rule is the SplitMix64 mix of (SplitMix64 mix of the run seed) + id,
inferred from the 883 orbits of the P1 panel file (0 differences); listed as INFO.
"""
import argparse, hashlib, sys, tomllib
from decimal import Decimal, getcontext
from fractions import Fraction

M64 = (1 << 64) - 1
MULT = 0x45A31EFC5A35D971261FD0407A968ADD
CLASS_ORDER = ["T", "R", "C", "S", "X"]
RES = []


def chk(id_, qty, val, bound, ok=None, info=False):
    st = "INFO" if info else ("PASS" if ok else "FAIL")
    RES.append(st)
    print("%s | %s | %s | %s | %s" % (id_, qty, val, bound, st))


def diff_check(id_, qty, ndiff, total):
    chk(id_, "%s (of %d)" % (qty, total), ndiff, "0", ndiff == 0)


# ------------------------------------------------------------------ generator and streams
class Gen:
    def __init__(self, seed):
        self.s = ((seed << 1) | 1) % (1 << 128)

    def word(self):
        self.s = (self.s * MULT) % (1 << 128)
        return self.s >> 64

    def uniform(self):
        return (self.word() & ((1 << 52) - 1)) * 2.0 ** -52


def sha_seed(text):
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest()[:8], "big")


def stream(text):
    return Gen(sha_seed(text))


def splitmix(x):
    x = (x + 0x9E3779B97F4A7C15) & M64
    z = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return z ^ (z >> 31)


def orbit_seed(run_seed, oid):
    return splitmix((splitmix(run_seed) + oid) & M64)


def generator_check(path):
    rows = [l.split() for l in open(path, encoding="utf-8") if l.strip() and not l.startswith("#")]
    seeds = sorted({int(r[0]) for r in rows})
    nu = nl = nh = n_u = n_f = 0
    for sd in seeds:
        g = Gen(sd)
        w = [g.word() for _ in range(250)]
        for r in rows:
            if int(r[0]) != sd:
                continue
            i = int(r[2]) - 1
            if r[1] == "u64":
                n_u += 1
                nu += w[i] == int(r[3])
            else:
                n_f += 1
                v = float(r[3])
                nl += (w[i] & ((1 << 52) - 1)) * 2.0 ** -52 == v
                nh += (w[i] >> 12) * 2.0 ** -52 == v
    chk("O23-gen-u64", "64-bit words equal to the library's (state, multiplier, high word)", "%d/%d" % (nu, n_u), "all", nu == n_u and n_u > 0)
    chk("O23-gen-f64", "uniform = low 52 bits * 2^-52 equals the library's Float64", "%d/%d" % (nl, n_f), "all", nl == n_f and n_f > 0)
    chk("O23-gen-alt", "alternative: high 52 bits * 2^-52 equals the library's Float64", "%d/%d" % (nh, n_f), "none expected", info=True)
    return nu == n_u and nl == n_f


# ------------------------------------------------------------------ draws
def fy(cands, n, g, exact_diff=None):
    """partial Fisher-Yates of E14 on the list cands (copied), n draws; returns the drawn items in order.
    j = t + floor(u (M - t + 1)) (1-based). exact_diff: list, appended with 1 when the exact rational floor differs."""
    c = list(cands)
    M = len(c)
    out = []
    for t in range(1, min(n, M) + 1):
        u = g.uniform()
        k = M - t + 1
        j = t + int(u * k)
        if exact_diff is not None:
            ex = t + (Fraction(u) * k).__floor__()
            if ex != j:
                exact_diff.append(1)
        c[t - 1], c[j - 1] = c[j - 1], c[t - 1]
        out.append(c[t - 1])
    return out


getcontext().prec = 110


def dec_pi():
    # Machin: pi = 16 atan(1/5) - 4 atan(1/239)
    def atan_inv(n):
        x = Decimal(1) / n
        s = x
        t = x
        k = 1
        while abs(t) > Decimal(10) ** -105:
            t = -t * x * x
            k += 2
            s += t / k
        return s
    return 16 * atan_inv(5) - 4 * atan_inv(239)


PI = dec_pi()


def dec_sincos(x):
    s = Decimal(0)
    c = Decimal(0)
    t = Decimal(1)
    k = 0
    while abs(t) > Decimal(10) ** -105 or k < 4:
        if k % 4 == 0:
            c += t
        elif k % 4 == 1:
            s += t
        elif k % 4 == 2:
            c -= t
        else:
            s -= t
        k += 1
        t = t * x / k
    return s, c


def twopi_f(u):
    """binary64 product (2 pi as binary64) * u, one binary64 rounding (the reading used by R and X placements if the 256-bit one does not match)"""
    return (2.0 * 3.141592653589793) * u


def twopi_u(u):
    """RN(2 pi u), evaluated far above 256 bits, one rounding to binary64"""
    return float(2 * PI * Decimal(u))


def parse_rational(s):
    s = str(s).strip()
    if "/" in s:
        a, b = s.split("/")
        return Fraction(Decimal(a)) / Fraction(Decimal(b))
    return Fraction(Decimal(s))


def runs_to_ids(runs):
    out = []
    for a, b in runs:
        out.extend(range(a, b + 1))
    return out


def alloc_db(b, counts):
    """counts: dict k_index -> n_i. largest remainder, ties to the smaller K index. returns dict k_index -> n_i'"""
    tot = sum(counts.values())
    a = {i: (b * n) // tot for i, n in counts.items()}
    rem = {i: (b * n) % tot for i, n in counts.items()}
    extra = b - sum(a.values())
    for i in sorted(counts, key=lambda i: (-rem[i], i))[:extra]:
        a[i] += 1
    return a


def build_db(seed, R, b, D_orbits, nK):
    """D_orbits: list of (id, k_index) of D for rounds 1..R; returns (sorted id list, per K index counts list)"""
    D_R = sorted(D_orbits)
    if len(D_R) <= b:
        sel = [i for i, _ in D_R]
    else:
        strata = {}
        for oid, ki in D_R:
            strata.setdefault(ki, []).append(oid)
        al = alloc_db(b, {ki: len(v) for ki, v in strata.items()})
        sel = []
        for ki in sorted(strata):
            if al[ki] > 0:
                sel.extend(fy(strata[ki], al[ki], stream("chaoslab-db-1|%d|%d|%d" % (seed, R, ki))))
    per = [0] * nK
    ks = dict(D_orbits)
    for i in sel:
        per[ks[i] - 1] += 1
    return sorted(sel), per


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel")
    ap.add_argument("--pools")
    ap.add_argument("--db")
    ap.add_argument("--generator")
    ap.add_argument("--generator-only")
    a = ap.parse_args()
    if a.generator_only:
        generator_check(a.generator_only)
        return fin()
    if a.generator:
        generator_check(a.generator)
    panel = tomllib.load(open(a.panel, "rb"))
    pools = tomllib.load(open(a.pools, "rb"))
    man = panel["manifest"]
    inp = man["inputs"]
    seed = int(inp["seed"]["value"])
    Kb = man["K_binary64"]
    nK = len(Kb)
    r_max = int(inp["r_max"]["value"])
    gates = list(inp["gate_rounds"]["value"])
    Mpts = int(inp["screen_points"]["value"])
    bsize = int(inp["db_size"]["value"]) if "db_size" in inp else 0
    quota = {int(i): d for i, d in inp["quota"]["value"].items()}
    dfrac = {int(i): {c: parse_rational(v) for c, v in d.items()} for i, d in inp["d_fraction"]["value"].items()}
    print("# inputs from the manifest: seed %d, %d K, r_max %d, gates %s, screen_points %d, db_size %d" % (seed, nK, r_max, gates, Mpts, bsize))
    orbits = [o for o in panel["orbit"] if "round" in o]
    byid = {o["id"]: o for o in orbits}

    # ---- pools from the sidecar
    pool = {}
    ok_pts = True
    for i in range(1, nK + 1):
        k = pools["k"][str(i)]
        ok_pts &= k["points"] == Mpts
        for c in "TCS":
            pool[(i, c)] = runs_to_ids(k[c])
    chk("O23-pool-points", "sidecar points per K equals screen_points", Mpts, "all K", ok_pts)
    cnt = {(p["k_index"], p["class"]): p["count"] for p in man["pools"]}
    bad = sum(len(pool[(i, c)]) != cnt.get((i, c), -1) for i in range(1, nK + 1) for c in "TCS" if not (i == 1 and Kb[0] == 0.0 and c != "T"))
    diff_check("O23-pool-counts", "manifest pool counts differ from sidecar run lists", bad, 3 * nK)
    # bit-vector hash: try a few encodings (the spec does not fix one) and report
    for i in range(1, nK + 1):
        k = pools["k"][str(i)]
        found = None
        for name, mp in [("T1C2S4", (1, 2, 4)), ("T1S2C4", (1, 4, 2)), ("T4C2S1", (4, 2, 1))]:
            bits = bytearray(Mpts)
            for c, m in zip("TCS", mp):
                for j in pool[(i, c)]:
                    bits[j - 1] |= m
            if hashlib.sha256(bytes(bits)).hexdigest() == k["bits_sha256"]:
                found = name
        chk("O23-bits-sha", "K index %d: bits_sha256 reproduced from the run lists (1 byte per point, mask)" % i, found or "no encoding matched", "info (format not fixed by the spec)", info=True)

    # ---- screening points of every K (E32), checked against drawn pool orbits
    pts = {}
    for i in range(1, nK + 1):
        g = stream("chaoslab-screen-1|%d|%d" % (seed, i))
        lst = []
        for _ in range(Mpts):
            x = twopi_u(g.uniform())
            p = twopi_u(g.uniform())
            lst.append((x, p))
        pts[i] = lst
    inrange = all(0.0 <= x < 2 * 3.141592653589793238 + 1e-15 and 0.0 <= p < 6.2831853071795871 for l in pts.values() for x, p in l)
    chk("O23-pts-range", "all %d points in [0, 2 pi)" % (nK * Mpts), inrange, "true", inrange)
    chk("O23-pts-distinct", "no two screening points equal within a K", all(len(set(l)) == len(l) for l in pts.values()), "true", all(len(set(l)) == len(l) for l in pts.values()))

    # ---- draws (E33 to E36)
    used = {i: set() for i in range(1, nK + 1)}
    cells = {}  # (r,i,c) -> dict(pool,q,drawn,status,draws)
    exact_diff = []
    isK0 = {i: Kb[i - 1] == 0.0 for i in range(1, nK + 1)}
    for r in range(1, r_max + 1):
        for i in range(1, nK + 1):
            for c in ["S", "C", "T"]:
                undefined = isK0[i] and c != "T"
                q = 0 if undefined else int(quota[i].get(c, 0))
                if isK0[i] and r not in gates:
                    q = 0
                cand = [j for j in pool[(i, c)] if j not in used[i]]
                if undefined:
                    cand_n = 0
                else:
                    cand_n = len(cand)
                dr = [] if (q == 0 or undefined) else fy(cand, q, stream("chaoslab-panel-1|%d|%d|%s|%d" % (seed, i, c, r)), exact_diff)
                used[i].update(dr)
                st = "undefined" if undefined else "not drawn" if q == 0 else "populated" if len(dr) == q else "short" if dr else "empty"
                cells[(r, i, c)] = dict(pool=(0 if undefined else cand_n), q=q, drawn=len(dr), status=st, draws=dr)
    chk("O23-fy-exact", "draws where the float floor(u (M-t+1)) differs from the exact rational floor", len(exact_diff), "info", info=True)

    # R and X: counts only here (starts are checked below)
    for r in range(1, r_max + 1):
        for i in range(1, nK + 1):
            for c in "RX":
                undefined = isK0[i] and c == "X"
                q = 0 if undefined else int(quota[i].get(c, 0))
                if isK0[i] and r not in gates:
                    q = 0
                cells[(r, i, c)] = dict(pool=None, q=q, drawn=None, status="undefined" if undefined else None, draws=None)

    # panel orbits per round cell (draw order by draw_index)
    got = {}
    for o in orbits:
        got.setdefault((o["round"], o["k_index"], o["class"]), []).append(o)
    for v in got.values():
        v.sort(key=lambda o: o["draw_index"])

    # drawn counts for R and X from the panel (their sizes are deterministic: q), status check
    for key, cl in cells.items():
        if cl["drawn"] is None:
            n = len(got.get(key, []))
            cl["drawn"] = n
            if cl["status"] is None:
                cl["status"] = "not drawn" if cl["q"] == 0 else "populated" if n == cl["q"] else "short" if n else "empty"

    # ---- ids (E36) and per-orbit comparison
    exp = []
    for r in range(1, r_max + 1):
        for i in range(1, nK + 1):
            for c in CLASS_ORDER:
                cl = cells[(r, i, c)]
                for t in range(1, cl["drawn"] + 1):
                    exp.append((r, i, c, t))
    nnc = [o for o in orbits if not o.get("control")]
    chk("O23-ids-count", "orbits of the panel (non-control) equals the rebuilt count", "%d vs %d" % (len(nnc), len(exp)), "equal", len(nnc) == len(exp))
    ids_ok = [o["id"] for o in nnc] == list(range(1, len(nnc) + 1))
    chk("O23-ids-run", "ids 1..N without a gap, panel order", ids_ok, "true", ids_ok)
    nd = 0
    for oid, (r, i, c, t) in enumerate(exp, 1):
        o = byid.get(oid)
        if o is None or (o["round"], o["k_index"], o["class"], o["draw_index"]) != (r, i, c, t):
            nd += 1
    diff_check("O23-ids-order", "orbits whose (round, K index, class, draw index) differs from the id order of E36", nd, len(exp))

    # pool orbits: point ids, starts
    nd_id = nd_pt = nd_dec = nd_cls = 0
    ndraw = 0
    for (r, i, c), cl in cells.items():
        if c not in "SCT" or cl["draws"] is None:
            continue
        ol = got.get((r, i, c), [])
        if len(ol) != len(cl["draws"]):
            nd_id += abs(len(ol) - len(cl["draws"]))
            continue
        for o, j in zip(ol, cl["draws"]):
            ndraw += 1
            if o.get("point_id") != j:
                nd_id += 1
            x, p = pts[i][j - 1]
            if "%.17g" % o["x0"] != "%.17g" % x or "%.17g" % o["p0"] != "%.17g" % p:
                nd_pt += 1
            if Fraction(Decimal(o["x0_dec"])) != Fraction(x) or Fraction(Decimal(o["p0_dec"])) != Fraction(p):
                nd_dec += 1
            if j not in set(pool[(i, c)]):
                nd_cls += 1
    diff_check("O23-pts-ids", "pool orbits whose point id differs from the rebuilt draw", nd_id, ndraw)
    diff_check("O23-pts-bits", "pool orbit starts differing in the 17-digit decimal text of x0 or p0", nd_pt, ndraw)
    diff_check("O23-pts-dec", "pool orbit starts whose exact decimal text differs from the binary64 pair", nd_dec, ndraw)
    diff_check("O23-pool-class", "drawn pool orbits whose point id is not in the sidecar pool of its class", nd_cls, ndraw)
    allpts = [(o["x0"], o["p0"]) for o in orbits]
    chk("O23-I29-starts", "equal starts (binary64 pair) among all %d orbits" % len(orbits), len(allpts) - len(set(allpts)), "0", len(allpts) == len(set(allpts)))
    twice = 0
    seen = set()
    for o in orbits:
        if "point_id" in o:
            k = (o["k_index"], o["point_id"])
            twice += k in seen
            seen.add(k)
    chk("O23-I29-ids", "point ids drawn twice at a K", twice, "0", twice == 0)
    chk("O23-rounds", "round of every orbit in 1..r_max", all(1 <= o["round"] <= r_max for o in orbits), "true", all(1 <= o["round"] <= r_max for o in orbits))
    gate_bad = sum(1 for o in orbits if isK0[o["k_index"]] and o["round"] not in gates)
    chk("O23-gate", "K = 0 orbits outside the gate rounds", gate_bad, "0", gate_bad == 0)

    # seeds
    nds = ndc = ndr = 0
    for o in nnc:
        r, i, c = o["round"], o["k_index"], o["class"]
        if int(o["cell_seed"]) != sha_seed("chaoslab-panel-1|%d|%d|%s|%d" % (seed, i, c, r)):
            ndc += 1
        if int(o["orbit_seed"]) != orbit_seed(seed, o["id"]):
            nds += 1
        if int(o["run_seed"]) != seed:
            ndr += 1
    diff_check("O23-seed-cell", "cell stream seeds (E33 text) differing", ndc, len(nnc))
    diff_check("O23-seed-run", "run seeds differing from the manifest seed", ndr, len(nnc))
    chk("O23-seed-orbit", "per-orbit seeds differing from the inferred rule (see header; E16's rule lives in another spec)", "%d of %d" % (nds, len(nnc)), "0", nds == 0)

    # round cell records
    mrc = {(c["round"], c["k_index"], c["class"]): c for c in man["round_cells"]}
    nd_rc = {"pool": 0, "q": 0, "drawn": 0, "status": 0}
    pool_info_rx = 0
    for key, cl in cells.items():
        m = mrc.get(key)
        if m is None:
            for k in nd_rc:
                nd_rc[k] += 1
            continue
        for k in ("q", "drawn", "status"):
            if m[k] != cl[k]:
                nd_rc[k] += 1
        if cl["pool"] is not None and m["pool"] != cl["pool"]:
            nd_rc["pool"] += 1
    n_rc = len(cells)
    chk("O23-roundcells-count", "round cell records in the manifest equals rounds x K x 5", "%d vs %d" % (len(mrc), n_rc), "equal", len(mrc) == n_rc)
    for k, v in nd_rc.items():
        diff_check("O23-roundcell-" + k, "round cell records whose %s differs (pool for S, C, T; R and X pools not rebuilt)" % k, v, n_rc)
    # sums per (K, class) against manifest.cells
    mc = {(c["k_index"], c["class"]): c for c in man["cells"]}
    bad = 0
    for (i, c), m in mc.items():
        if m["drawn"] != sum(cells[(r, i, c)]["drawn"] for r in range(1, r_max + 1)):
            bad += 1
    diff_check("O23-cells-sum", "manifest cells whose drawn differs from the sum over rounds", bad, len(mc))
    tot_ok = man["totals"]["orbits"] == len(orbits)
    chk("O23-totals", "manifest totals.orbits equals the orbit count", "%d vs %d" % (man["totals"]["orbits"], len(orbits)), "equal", tot_ok)

    # ---- D (E37)
    expD = set()
    nd_dsize = 0
    nd_dman = 0
    for key in sorted(cells):
        r, i, c = key
        cl = cells[key]
        n = cl["drawn"]
        if n == 0:
            continue
        f = dfrac.get(i, {}).get(c, Fraction(0))
        m = (f * n + Fraction(1, 2)).__floor__()
        sel = fy(list(range(1, n + 1)), m, stream("chaoslab-d-1|%d|%d|%d|%s" % (seed, r, i, c)))
        for t in sel:
            expD.add(got[key][t - 1]["id"])
        have = sum(1 for o in got[key] if o["in_D"])
        if have != m:
            nd_dsize += 1
    nd_flag = sum(1 for o in nnc if o["in_D"] != (o["id"] in expD))
    diff_check("O23-D-size", "round cells whose D count differs from floor(f n + 1/2)", nd_dsize, len(cells))
    diff_check("O23-D-flags", "orbits whose in_D flag differs from the rebuilt D", nd_flag, len(nnc))
    mdf = man.get("d_fraction", {})
    bad = 0
    for i, d in mdf.items():
        for c, v in d.items():
            if parse_rational(v) != dfrac.get(int(i), {}).get(c):
                bad += 1
    chk("O23-D-fracs", "manifest d_fraction tables differ from the input tables", bad, "0", bad == 0)
    chk("O23-D-total", "|D| rebuilt", len(expD), "info", info=True)

    # ---- D_B (E38)
    Dk = [(o["id"], o["k_index"]) for o in nnc if o["id"] in expD]
    exp_db = {}
    if a.db:
        ex = tomllib.load(open(a.db, "rb"))
        chk("O23-DB-seed", "export seed / r_max / db_size equal the manifest", (ex.get("seed"), ex.get("r_max"), ex.get("db_size")), (seed, r_max, bsize), (ex.get("seed"), ex.get("r_max"), ex.get("db_size")) == (seed, r_max, bsize))
    for R in range(1, r_max + 1):
        DR = [(i, k) for i, k in Dk if byid[i]["round"] <= R]
        sel, per = build_db(seed, R, bsize, DR, nK)
        exp_db[R] = (sel, per)
        want = min(bsize, len(DR))
        chk("O23-DB-size-R%d" % R, "|D_B| = min(db_size, |D_R|)", "%d (|D_R| = %d)" % (len(sel), len(DR)), want, len(sel) == want)
        if a.db:
            e = ex["R%d" % R]
            nd = len(set(sel) ^ set(e["DB_ids"]))
            chk("O23-DB-ids-R%d" % R, "D_B ids in the symmetric difference with the exported D_B", nd, "0", nd == 0)
            chk("O23-DB-per-K-R%d" % R, "D_B per K index equals the export", "%s vs %s" % (per, e["DB_per_k_index"]), "equal", per == list(e["DB_per_k_index"]))
            nd = len(set(i for i, _ in DR) ^ set(e["D_ids"]))
            chk("O23-D-ids-R%d" % R, "D_R ids in the symmetric difference with the exported D_R", nd, "0", nd == 0 and len(DR) == e["D_size"])
            chk("O23-DB-n-R%d" % R, "export DB_size equals", "%d vs %d" % (e["DB_size"], len(sel)), "equal", e["DB_size"] == len(sel))
    if not a.db:
        chk("O23-DB-export", "no D_B export given: D_B checked for size only", "-", "-", info=True)

    # ---- R and X starts
    check_rx(panel, man, inp, seed, nK, isK0, got, cells, r_max, quota)
    return fin()


def check_rx(panel, man, inp, seed, nK, isK0, got, cells, r_max, quota):
    """R and X starts (E12, E13). The spec states the 256-bit evaluation for the screening points (E32) and for the R sum; for 2 pi u in
    phi and in the K = 0 x start it states none, so both readings (256-bit RN, and the binary64 product 2 pi * u) are tried and reported."""
    import math
    q_max = int(inp["q_max"]["value"])
    delta_R = inp["delta_R"]["value"]
    lo, hi = inp["log10_r"]["value"]
    ell = []
    for po in man["periodic"]:
        if po["class"] in ("elliptic", "parabolic") and po["q"] <= q_max:
            for pi_, pt in enumerate(po["points"], 1):
                ell.append((po["k_index"], po["id"], pi_, Decimal(pt[0]), Decimal(pt[1])))
    T = lambda v: "%.17g" % v
    var = {}

    def tally(name, ok):
        var[name] = var.get(name, 0) + bool(ok)

    nR0 = nRs = nRx = nXphi = nXx = nR0tot = nRtot = nXtot = nXr0 = nXr1 = 0
    maxulp = 0.0
    nXx2 = nXx_exact = 0
    for (r, i, c), ol in sorted(got.items()):
        if c not in "RX":
            continue
        g = stream("chaoslab-panel-1|%d|%d|%s|%d" % (seed, i, c, r))
        for o in sorted(ol, key=lambda o: o["draw_index"]):
            u1 = g.uniform()
            u2 = g.uniform()
            ph_a, ph_b = twopi_u(u2), twopi_f(u2)
            if c == "R" and isK0[i]:
                nR0tot += 1
                # order of the reduced fractions a/b, 0 <= a/b < 1, 1 <= b <= q_max: by b, then by a (the spec fixes no order; ascending value was tried and fails)
                fr = [Fraction(a, b) for b in range(1, q_max + 1) for a in range(0, b) if math.gcd(a, b) == 1]
                fr_v = sorted(fr)
                f = fr[int(u1 * len(fr))]
                tally("R K=0 fraction index (ascending value order instead)", (o["a"], o["b"]) == (lambda z: (z.numerator, z.denominator))(fr_v[int(u1 * len(fr_v))]))
                p = float(2 * PI * Decimal(f.numerator) / Decimal(f.denominator))
                idx_ok = (o["a"], o["b"]) == (f.numerator, f.denominator)
                tally("R K=0 fraction index (order by b, then a)", idx_ok)
                tally("R K=0 p0 = RN(2 pi a/b) (256-bit)", T(p) == T(o["p0"]))
                tally("R K=0 x0 = RN(2 pi u2) (256-bit)", T(ph_a) == T(o["x0"]))
                tally("R K=0 x0 = binary64 product 2pi*u2", T(ph_b) == T(o["x0"]))
                if not (idx_ok and T(p) == T(o["p0"]) and (T(ph_a) == T(o["x0"]) or T(ph_b) == T(o["x0"]))):
                    nR0 += 1
            elif c == "R":
                nRtot += 1
                pts_i = [e for e in ell if e[0] == i]
                e = pts_i[int(u1 * len(pts_i))]
                nRs += (o["source_orbit"], o["source_point"]) != (e[1], e[2])
                tally("R stored phi = RN(2 pi u2) (256-bit)", T(ph_a) == T(o["phi"]))
                tally("R stored phi = binary64 product 2pi*u2", T(ph_b) == T(o["phi"]))
                okany = False
                for pn, ph in [("phi 256-bit", 2 * PI * Decimal(u2)), ("phi RN", Decimal(ph_a)), ("phi binary64 product", Decimal(ph_b)), ("phi stored", Decimal(o["phi"]))]:
                    s, cc = dec_sincos(ph)
                    ok = T(float(e[3] + Decimal(delta_R) * cc)) == T(o["x0"]) and T(float(e[4] + Decimal(delta_R) * s)) == T(o["p0"])
                    tally("R start (256-bit sum, delta_R binary64) with " + pn, ok)
                    okany |= ok
                nRx += not okany
            else:
                nXtot += 1
                l_d = Decimal(lo) + (Decimal(hi) - Decimal(lo)) * Decimal(u1)
                l_f = lo + (hi - lo) * u1
                tally("X r = 10^l (l 256-bit)", T(float(Decimal(10) ** l_d)) == T(o["r"]))
                tally("X r = 10^l (l binary64, exact pow)", T(float(Decimal(10) ** Decimal(l_f))) == T(o["r"]))
                tally("X r = libm 10.0**l", T(10.0 ** l_f) == T(o["r"]))
                tally("X stored phi = RN(2 pi u2) (256-bit)", T(ph_a) == T(o["phi"]))
                tally("X stored phi = binary64 product 2pi*u2", T(ph_b) == T(o["phi"]))
                if T(ph_a) != T(o["phi"]) and T(ph_b) != T(o["phi"]):
                    nXphi += 1
                okany = False
                for rn, rd in [("r 256-bit", Decimal(10) ** l_d), ("r stored", Decimal(o["r"]))]:
                    for pn, ph in [("phi 256-bit", 2 * PI * Decimal(u2)), ("phi RN", Decimal(ph_a)), ("phi binary64 product", Decimal(ph_b)), ("phi stored", Decimal(o["phi"]))]:
                        s, cc = dec_sincos(ph)
                        ok = T(float(rd * cc)) == T(o["x0"]) and T(float(rd * s)) == T(o["p0"])
                        tally("X start (256-bit product) from %s, %s" % (rn, pn), ok)
                        okany |= ok
                for pn, ph in [("phi RN", ph_a), ("phi binary64 product", ph_b), ("phi stored", o["phi"])]:
                    for rn, rv in [("r stored", o["r"]), ("r libm 10.0**l", 10.0 ** l_f)]:
                        ok = (rv * math.cos(ph) == o["x0"]) and (rv * math.sin(ph) == o["p0"])
                        tally("X start (binary64 r*cos, r*sin) from %s, %s" % (rn, pn), ok)
                        if rn == "r stored" and pn in ("phi binary64 product", "phi stored"):
                            okany_f = ok
                nXx += not okany_f
                ph_exact = Decimal(o["phi"])
                s_, c_ = dec_sincos(ph_exact)
                for val, ex_ in ((o["x0"], Decimal(o["r"]) * c_), (o["p0"], Decimal(o["r"]) * s_)):
                    d_ = abs(Decimal(val) - ex_) / Decimal(math.ulp(val))
                    maxulp = max(maxulp, float(d_))
                    nXx2 += d_ > 2
                nXx_exact += okany_f
                rd = float(Decimal(10) ** l_d)
                nXr1 += abs(rd - o["r"]) > math.ulp(o["r"])
                nXr0 += rd != o["r"]
    if nR0tot:
        chk("O23-R0-start", "K = 0 R orbits unmatched (fraction index, p0, x0 under either reading of 2 pi u2)", nR0, "0 of %d" % nR0tot, nR0 == 0)
    if nRtot:
        chk("O23-R-source", "K > 0 R orbits whose (source orbit, source point) differs from the rebuilt index", nRs, "0 of %d" % nRtot, nRs == 0)
        chk("O23-R-start", "K > 0 R orbits whose start matches under no variant", nRx, "0 of %d" % nRtot, nRx == 0)
    if nXtot:
        chk("O23-X-phi", "X orbits whose stored phi matches under neither reading of 2 pi u2", nXphi, "0 of %d" % nXtot, nXphi == 0)
        chk("O23-X-r", "X orbits whose stored r differs from the correctly rounded 10^l by more than 1 ulp", nXr1, "0 of %d" % nXtot, nXr1 == 0)
        chk("O23-X-r-exact", "X orbits whose stored r differs from the correctly rounded 10^l (any amount; the spec fixes no pow rounding)", nXr0, "info of %d" % nXtot, info=True)
        chk("O23-X-start", "X orbits whose start is farther than 2 ulp from the exact r cos(phi), r sin(phi) of the stored r and phi (E12; start formed in binary64)", nXx2, "0 of %d (max deviation %.2f ulp)" % (nXtot, maxulp), nXx2 == 0)
        chk("O23-X-start-exact", "X orbits whose start equals binary64 r*cos(phi), r*sin(phi) with this platform's cos and sin (the spec fixes no cos rounding)", nXtot - nXx, "info of %d" % nXtot, info=True)
    tot = {"R": nRtot, "X": nXtot, "R K=0": nR0tot}
    for k, v in sorted(var.items()):
        n = tot["R K=0"] if k.startswith("R K=0") else tot[k[0]]
        chk("O23-variant", "%s: matches" % k, "%d of %d" % (v, n), "info (which reading matched)", info=True)


def fin():
    npass, nfail, ninfo = RES.count("PASS"), RES.count("FAIL"), RES.count("INFO")
    print("TOTAL | %d checks | PASS %d | FAIL %d | INFO %d" % (len(RES), npass, nfail, ninfo))
    sys.exit(1 if nfail else 0)


main()
