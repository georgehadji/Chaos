#!/usr/bin/env python
"""Independent decision layer for study p2-indicators (oracle O6 of the analysis spec, Amendment D89).

Written from analysis-spec.md (E23, E26 to E35, O7, O8) and the frozen preregistration; numpy and scipy only; the module
(src/Analysis) was not read.  Input: the table export of verify/analysis-p2-independent.jl (JSON lines, format documented there),
the parameter table (every threshold comes from it), the pilot evidence file (O8).
Usage (repository root):  python verify/analysis-p2-decision.py [export file] [results file]
Output: workspaces/module-build/stages/04-verify/output/analysis-p2-independent-results.txt, one line per quantity: `name | value`.

Names (all names are `|`-free tokens joined by `|` inside the name field; the separator between name and value is ` | `):
  size|<set>[|P|G|K=<K>]            set sizes: Eall E8 W2 S2 Ts Ts_design regular gate SD_star SL_star, W2_8 ... (E8 suffix 8)
  drop|<reason>[|P|G]               controls, flagged, rhoU.     rho|<P|G>|K=<K>|<C|R|U>      rho counts      rhocell|<kidx>|<class>|<C|R|U>
  label|<kidx or 0>|<class or controls>|<instability>|<structure>   panel label pairs per cell (T6)
  step|<LE|SALI|FLI|MEGNO|RLI|STR>|k<k>      the recorded step n(b_k) used           lamref_stat|<name>   few global checks
  vcount|<rule>|<set>|<X>|k<k>|<chaotic|regular|undecided|na>          T2 counts, rule R1 or R2, sets W2 S2 Ts regular
  det|<rule>|<set>|<X>|<correct|wrong|undetected|M|n>                 outcomes, median censored time M_X, orbits with an item
  tdet|<rule>|<set>|<X>|<value>                                       orbits with detection time <value>  (F2)
  <ver>|<item>|<quantity>   with ver = R1.I7 (primary), R1.I8 (extension on E8), R2.I7, R2.I8, and <ver>.minus-<X> for a violating gate indicator
     H1|n_W2 U U_lo U_hi U_drop Amin A|<X>-<Y> U_S2 decision;  H2|n_S2 M|<X> Q Q_lo Q_hi decision;  H3|n_W2 n_S2 MW|<X> MS|<X> R|<X> R_lo|<X> R_hi|<X> npass decision;
     H4|n_Ts phi|<X> dphi dphi_lo dphi_hi decision, H4design|... ;  H5|n_W2 F|<X> Fmed Fmed_lo Fmed_hi decision;  H6|n_SDs n_SLs ddir|<X> ndrop_dir dprec|<X> ndrop_prec max_dir max_prec decision
     robust|H<n> (ROBUST, NOT ROBUST, NOT COMPARABLE) is given once for each rule: <rule>|robust|H<n>
  agree|<set>|k<k>|<X>-<Y>|{A,kappa,N,ncc,ncr,nrc,nrr}               F1 pairs of I_8 for sets W2 S2;    meankappa|<set>|k<k>|I7    spearman|<set>|k<k>|<X>-<Y>
  Fb|<set>|k<k>|<X> (F_X at every budget), Phi|<set>|<X>|<count of flips>|<orbits>   (T7)
  auc|<X>|k<k>|K=<K>|<pooled|P|G>|{auc,nC,nR}      regular range  rng|<X>|k<k>|K=<K>|{lo,hi,median,n}     leave|<X>|K=<K>|{n,defined,below_median,n9}     leavek|<X>|K=<K>|<b>
  gate|<X>|{decided,chaotic,undecided,na,share,violated}     k0line|<fq1|pa1>|W=<W>|maxdist     ctrl|<id>|<X>|v9
  T11|<pair>|k<k>|{n,A,kappa,medabs,maxabs}      T15|<pair>|<set>|k<k>|<vX,vY>      selfcheck|... PASS|FAIL
"""
import sys, json, math, hashlib, itertools, re, time, tomllib
from collections import defaultdict
import numpy as np
from scipy.stats import rankdata, spearmanr

T0 = time.time()
ROOT = __file__.replace("\\", "/").rsplit("/verify/", 1)[0]
SCRATCH = "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/p2ind"
EXPORT = sys.argv[1] if len(sys.argv) > 1 else SCRATCH + "/p2-export.jsonl"
OUT = sys.argv[2] if len(sys.argv) > 2 else ROOT + "/workspaces/module-build/stages/04-verify/output/analysis-p2-independent-results.txt"
PARF = ROOT + "/workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml"
PILOT = ROOT + "/workspaces/study/stages/03-run/output/p2-indicators-pilot-evidence.md"


def say(*a):
    print("[%.1f s]" % (time.time() - T0), *a, flush=True)


with open(PARF, "rb") as fh:
    PT = tomllib.load(fh)
P = {k: v["value"] for k, v in PT.items()}
BUD, NC, N1, NS, NDIR = list(P["budgets"]), int(P["N_c"]), int(P["N_1"]), int(P["n_s"]), int(P["N_dir"])
KP = BUD.index(P["b_primary"])
DIRS = list(P["directions"])
MINSET = int(P["min_set"])
I7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA"]
I8 = I7 + ["RLI"]
FOUR = ["LE", "SALI", "FLI", "MEGNO"]
NB = len(BUD)

RES = []          # (name, value) lines


def out(name, val):
    if isinstance(val, (float, np.floating)):
        val = repr(float(val))
    elif isinstance(val, (bool, np.bool_)):
        val = "true" if val else "false"
    RES.append((name, str(val)))


# ------------------------------------------------------------------------------------------------ load the export
say("loading export", EXPORT)
ORBS, SER, WIN, CHF, NA, META = {}, {}, {}, [], [], []
with open(EXPORT, "r", encoding="utf-8") as fh:
    for line in fh:
        d = json.loads(line)
        k = d["k"]
        if k == "orbit":
            ORBS[d["id"]] = d
        elif k == "s":
            SER[(d["o"], d["spec"], d["prec"], d["hor"], d["dir"])] = (np.array(d["t"], dtype=np.int64), {c: np.array(v, dtype=float) for c, v in d["c"].items()})
        elif k == "w":
            WIN[(d["o"], d["spec"], d["prec"])] = d
        elif k == "chunkflags":
            CHF.append(d)
        elif k == "na":
            NA.append(d)
        else:
            META.append(d)
say("loaded: orbits", len(ORBS), "series", len(SER), "windows", len(WIN))
NORB = len(ORBS)
assert sorted(ORBS) == list(range(1, NORB + 1))
OI = lambda i: i - 1

# ------------------------------------------------------------------------------------------------ orbit table
ids = np.arange(1, NORB + 1)
isG = np.array([ORBS[i]["src"] == "G" for i in ids])
Kval = np.array([ORBS[i]["K"] if ORBS[i]["K"] is not None else np.nan for i in ids])
kidx = np.array([ORBS[i]["kidx"] for i in ids])
isCtrl = np.array([ORBS[i]["control"] for i in ids])
cls = np.array([ORBS[i]["cls"] for i in ids])
inSD = np.array([ORBS[i]["in_SD"] for i in ids])
inSL = np.array([ORBS[i]["in_SL"] for i in ids])
flagged = np.array([("drift" in ORBS[i]["flags"]) or ("excluded" in ORBS[i]["flags"]) for i in ids])

# ------------------------------------------------------------------------------------------------ recorded steps (E26)
def n_of_b(t, b):
    i = int(np.searchsorted(t, b, side="left"))
    if i >= len(t):
        raise RuntimeError("budget above the horizon")
    return i


def nprime(t, n):
    d = np.abs(10 * t - n)               # |t - n/10| scaled by 10, exact in integers
    return int(np.nonzero(d == d.min())[0][-1])   # a tie goes to the larger step


def tidx(t, n):
    i = int(np.searchsorted(t, n, side="left"))
    if i >= len(t) or t[i] != n:
        raise RuntimeError("step %d not recorded" % n)
    return i


# ------------------------------------------------------------------------------------------------ reference class (E23)
POS, NONE_, UNR = "positive growth supported", "none detected on the tested horizon", "unresolved"
rho = np.empty(NORB, dtype="U1")
inst_ax = np.empty(NORB, dtype=object)
str_ax = np.empty(NORB, dtype=object)
lamref = np.full(NORB, np.nan)
lamns = np.full(NORB, np.nan)
for i in ids:
    o = ORBS[i]
    t, c = SER[(i, "L1", "Float64", N1, DIRS[0])]
    lamref[OI(i)] = c["lam1"][tidx(t, N1)]
    lamns[OI(i)] = c["lam1"][tidx(t, NS)]
    if o["src"] == "P":
        ia, sa = o["inst"], o["struct"]
    else:
        t5, c5 = SER[(i, "L5", "Float64", N1, DIRS[0])]
        lam = c5["lam1"][tidx(t5, N1)]
        lam10 = c5["lam1"][tidx(t5, N1 // 10)]
        ia = UNR
        if math.isfinite(lam) and math.isfinite(lam10) and lam10 > 0:
            r = 10 * lam / lam10
            if lam >= P["label_lambda_min"] and r >= P["label_growth_ratio"]:
                ia = POS
            elif lam < P["label_lambda_min"] and r < P["label_growth_ratio"]:
                ia = NONE_
        h = N1 // 2
        xh, xN = c5["xlift"][tidx(t5, h)], c5["xlift"][tidx(t5, N1)]
        x0 = float(o["ic"].split(",")[0])
        rho1 = (xh - x0) / (2 * math.pi * h)
        rho2 = (xN - xh) / (2 * math.pi * (N1 - h))
        if ia == POS:
            sa = "irregular"
        elif ia == NONE_ and abs(rho1 - rho2) <= P["label_rotation_tol"]:
            sa = "quasiperiodic"
        else:
            sa = UNR
        lam1 = c["lam1"][tidx(t, N1)]
        if abs(lam1 - lam) > P["label_estimator_bound"]:
            ia, sa = UNR, UNR
    inst_ax[OI(i)], str_ax[OI(i)] = ia, sa
    rho[OI(i)] = "C" if (ia == POS and sa == "irregular") else ("R" if (ia == NONE_ and sa == "quasiperiodic") else "U")

# ------------------------------------------------------------------------------------------------ verdict and score of one indicator at one budget
# verdict codes: 1 chaotic, 0 regular, -1 undecided; 2 not applicable (no item).  Returns (R1, R2, score) with score nan = none.
def growth(X, t, c, b):
    i = n_of_b(t, b)
    n = int(t[i])
    j = nprime(t, n)
    if X == "LE":
        gn, gp = n * c["lam1"][i], int(t[j]) * c["lam1"][j]
    elif X == "SALI":
        gn, gp = -c["SALI.L"][i], -c["SALI.L"][j]
    elif X == "FLI":
        gn, gp = c["FLI.fli"][i], c["FLI.fli"][j]
    else:
        gn, gp = c["MEGNO.Ybar20_1"][i], c["MEGNO.Ybar20_1"][j]
    ok = math.isfinite(gn) and math.isfinite(gp) and gp > 0
    r = gn / gp if ok else math.nan
    if not ok:
        v1 = -1
    elif X == "LE":
        v1 = 1 if (c["lam1"][i] >= P["le_floor"] and r >= P["ratio_threshold"]) else 0
    else:
        v1 = 1 if r >= P["ratio_threshold"] else 0
    if X == "LE":
        val = gn
        v2 = -1 if not math.isfinite(val) else (1 if val >= P["r2_le_factor"] * math.log(n) else 0)
    elif X == "SALI":
        val = c["SALI.sali"][i]
        v2 = -1 if not math.isfinite(val) else (1 if val <= P["r2_sali_threshold"] else 0)
    elif X == "FLI":
        val = c["FLI.fli"][i]
        v2 = -1 if not math.isfinite(val) else (1 if val >= P["r2_fli_factor"] * math.log(n) else 0)
    else:
        val = c["MEGNO.Ybar20_1"][i]
        v2 = -1 if not math.isfinite(val) else (1 if val >= P["r2_megno_threshold"] else 0)
    sc = r if (ok and math.isfinite(r)) else math.nan
    return v1, v2, sc, n


def pair_rli(t, c, b):
    i = n_of_b(t, b)
    v = c["RLI.RLI"][i]
    if not math.isfinite(v):
        return -1, -1, math.nan, int(t[i])
    vv = 1 if v >= P["rli_threshold"] else 0
    sc = math.log10(v) if v > 0 else (-math.inf if v == 0 else math.nan)
    return vv, vv, sc, int(t[i])


def stretch_score(t, c, b):
    i = n_of_b(t, b)
    v = c["stretching.D"][i]
    return (-math.log10(v) if (math.isfinite(v) and v > 0) else math.nan), int(t[i])


def bit(w, bits, name):
    return (int(w) >> bits.index(name)) & 1 == 1


def rowindex(rows, W):
    return rows.index("W=%d" % W)


def win_e(X, d, b):
    """E1 (fq1) or E2 (pa1) at budget b from a window item d."""
    W = b // 2
    assert b % 2 == 0 and W in P["windows"]
    r = rowindex(d["rows"], W)
    D = d["c"]["D"][r]
    bad = bit(d["flags"][r], d["bits"], "degenerate" if X == "E1" else "undefined window")
    thr = P["e1_d_threshold"] if X == "E1" else P["e2_d_threshold"]
    if bad or not math.isfinite(D):
        return -1, -1, math.nan
    v = 1 if D >= thr else 0
    return v, v, D


def wba_rows(d, W):
    pre = "W=%d " % W
    idx = [k for k, r in enumerate(d["rows"]) if r.startswith(pre)]
    assert len(idx) == 4
    return idx


def win_wba(d, b, single=None):
    W = b // 2
    assert b % 2 == 0 and W in P["windows"]
    idx = wba_rows(d, W)
    if single is not None:
        idx = [k for k in idx if d["rows"][k] == "W=%d %s" % (W, single)]
        assert len(idx) == 1
    digs = [d["c"]["dig"][k] for k in idx]
    if any(bit(d["flags"][k], d["bits"], "zero") for k in idx) or any(not math.isfinite(x) for x in digs):
        return -1, -1, math.nan
    dm = min(digs)
    v = 1 if dm < P["wba_digits_threshold"] else 0
    return v, v, -dm


def indicator(X, o, b, prec="Float64", hor1=N1, direction=DIRS[0], hort=None):
    """(R1, R2, score) for orbit o, indicator X, budget b; items of the primary rows unless hor1/prec/direction given."""
    if X in FOUR:
        key = (o, "L1" if X == "LE" else "tg", prec, hor1, direction)
        if key not in SER:
            return 2, 2, math.nan, 0
        t, c = SER[key]
        return growth(X, t, c, b)
    if X == "RLI":
        key = (o, "pr", "Float64", NC, "")
        if key not in SER:
            return 2, 2, math.nan, 0
        return pair_rli(*SER[key], b)
    spec = {"E1": "fq1", "E2": "pa1", "WBA": "bk1"}[X]
    key = (o, spec, prec)
    if key not in WIN:
        return 2, 2, math.nan, 0
    d = WIN[key]
    r = win_wba(d, b) if X == "WBA" else win_e(X, d, b)
    return r[0], r[1], r[2], 0


# ------------------------------------------------------------------------------------------------ verdict and score tables (primary rows)
say("verdicts and scores")
IND = I8 + ["STR"]
V = {"R1": {X: np.full((NORB, NB), 2, dtype=np.int8) for X in I8}, "R2": {X: np.full((NORB, NB), 2, dtype=np.int8) for X in I8}}
SC = {X: np.full((NORB, NB), np.nan) for X in IND}
STEP = {}
for i in ids:
    for X in I8:
        for k, b in enumerate(BUD):
            v1, v2, sc, n = indicator(X, int(i), b)
            V["R1"][X][OI(i), k] = v1
            V["R2"][X][OI(i), k] = v2
            SC[X][OI(i), k] = sc
            if X in FOUR + ["RLI"] and n:
                STEP.setdefault((X, k), set()).add(n)
    key = (int(i), "st", "Float64", NC, DIRS[0])
    if key in SER:
        t, c = SER[key]
        for k, b in enumerate(BUD):
            sc, n = stretch_score(t, c, b)
            SC["STR"][OI(i), k] = sc
            STEP.setdefault(("STR", k), set()).add(n)
for (X, k), s in sorted(STEP.items()):
    out("step|%s|k%d" % (X, k), "/".join(str(x) for x in sorted(s)))

# ------------------------------------------------------------------------------------------------ sets (E29)
kv = lambda x: ("%r" % float(x)) if not np.isnan(x) else "none"
ctrl = isCtrl
Eall = (~ctrl) & (~flagged) & (rho != "U")
Cm, Rm = Eall & (rho == "C"), Eall & (rho == "R")
S2 = Cm & (lamref >= P["s2_lambda_min"])
Ts = Cm & (lamns <= P["f_s"] * lamref)
W2 = Cm & ((lamref <= P["w2_lambda_max"]) | Ts)
Tsd = Eall & (~isG) & (cls == "S") & (rho == "C")
REG = Rm
inE8 = inSL | isG
E8 = Eall & inE8
SDs, SLs = Eall & inSD, Eall & inSL
GATE = (~isG) & (kidx == P["gate_k_index"]) & (~ctrl) & np.isin(cls, P["gate_classes"])
assert GATE.sum() == P["gate_size"], GATE.sum()
SETS = {"Eall": Eall, "E8": E8, "W2": W2, "S2": S2, "Ts": Ts, "Ts_design": Tsd, "regular": REG, "gate": GATE, "SD_star": SDs, "SL_star": SLs,
        "W2_8": W2 & E8, "S2_8": S2 & E8, "Ts_8": Ts & E8, "Ts_design_8": Tsd & E8, "regular_8": REG & E8, "SD_star_8": SDs & E8, "SL_star_8": SLs & E8}
for nm, m in SETS.items():
    out("size|" + nm, int(m.sum()))
    out("size|%s|P" % nm, int((m & ~isG).sum()))
    out("size|%s|G" % nm, int((m & isG).sum()))
    for K in sorted(set(Kval[m & ~np.isnan(Kval)])):
        out("size|%s|K=%s" % (nm, kv(K)), int((m & (Kval == K)).sum()))
out("size|W2_and_S2", int((W2 & S2).sum()))
out("drop|controls", int(ctrl.sum()))
out("drop|flagged", int((flagged & ~ctrl).sum()))
out("drop|flagged_ids", ",".join(str(i) for i in ids[flagged & ~ctrl]))
for s, sm in (("P", ~isG), ("G", isG)):
    out("drop|rhoU|" + s, int(((rho == "U") & ~ctrl & ~flagged & sm).sum()))
    for r in "CRU":
        out("rho|%s|all|%s" % (s, r), int(((rho == r) & sm).sum()))
        for K in sorted(set(Kval[sm & ~np.isnan(Kval)])):
            out("rho|%s|K=%s|%s" % (s, kv(K), r), int(((rho == r) & sm & (Kval == K)).sum()))
for kx in sorted(set(kidx[~isG])):
    for cl in sorted(set(cls[(~isG) & (kidx == kx)])):
        m = (~isG) & (kidx == kx) & (cls == cl)
        for r in "CRU":
            out("rhocell|%d|%s|%s" % (kx, cl, r), int((m & (rho == r)).sum()))
cnt = defaultdict(int)
for i in ids[~isG]:
    key = (0, "controls") if isCtrl[OI(i)] else (int(kidx[OI(i)]), cls[OI(i)])
    cnt[key + (inst_ax[OI(i)], str_ax[OI(i)])] += 1
for key, v in sorted(cnt.items()):
    out("label|%d|%s|%s|%s" % key, v)
out("sdim|SD_total", int(inSD.sum()))
out("sdim|SL_total", int(inSL.sum()))
say("sets: Eall", Eall.sum(), "W2", W2.sum(), "S2", S2.sum(), "Ts", Ts.sum(), "Ts_design", Tsd.sum(), "regular", REG.sum(), "E8", E8.sum())

# ------------------------------------------------------------------------------------------------ detection (E28)
def detect(vs, rho_o, ctx_b=BUD, nc=NC):
    """vs: verdict series (length 10) of one orbit; returns (tdet or None, outcome, tcens)."""
    c9 = vs[NB - 1]
    tdet = None
    if c9 in (0, 1):
        for k in range(NB):
            if ctx_b[k] <= nc / P["persist_divisor"] and all(vs[j] == c9 for j in range(k, NB)):
                tdet = ctx_b[k]
                break
    if tdet is None:
        return None, "undetected", P["censor_factor"] * nc
    right = (c9 == 1 and rho_o == "C") or (c9 == 0 and rho_o == "R")
    return tdet, ("correct" if right else "wrong"), (tdet if right else P["censor_factor"] * nc)


def outcomes_for(rule):
    O = {}
    for X in I8:
        tc = np.full(NORB, np.nan)
        oc = np.full(NORB, "", dtype=object)
        td = np.full(NORB, np.nan)
        for i in ids[Eall]:
            vs = V[rule][X][OI(i)]
            if vs[0] == 2:
                continue
            a, b_, c_ = detect([int(x) for x in vs], rho[OI(i)])
            tc[OI(i)], oc[OI(i)] = c_, b_
            td[OI(i)] = a if a is not None else np.nan
        O[X] = (tc, oc, td)
    return O


say("outcomes")
OUTC = {r: outcomes_for(r) for r in ("R1", "R2")}


def median(a):
    a = np.sort(np.asarray(a, dtype=float))
    m = len(a)
    if m == 0:
        return math.nan
    return a[m // 2] if m % 2 else (a[m // 2 - 1] + a[m // 2]) / 2


# ------------------------------------------------------------------------------------------------ agreement (E30)
def pair_counts(vx, vy):
    dec = ((vx == 0) | (vx == 1)) & ((vy == 0) | (vy == 1))
    return (int(((vx == 1) & (vy == 1) & dec).sum()), int(((vx == 1) & (vy == 0) & dec).sum()),
            int(((vx == 0) & (vy == 1) & dec).sum()), int(((vx == 0) & (vy == 0) & dec).sum()))


def a_kappa(ncc, ncr, nrc, nrr):
    N = ncc + ncr + nrc + nrr
    if N == 0:
        return math.nan, math.nan
    po = (ncc + nrr) / N
    num = (ncc + ncr) * (ncc + nrc) + (nrc + nrr) * (ncr + nrr)
    N2 = N * N
    if num < N2:
        pe = num / N2
        return po, (po - pe) / (1 - pe)
    return po, 0.0


def spearman_corr(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return math.nan
    rx, ry = rankdata(x[m]), rankdata(y[m])
    if np.all(rx == rx[0]) or np.all(ry == ry[0]):
        return math.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def U_amin(Vr, mask, k, Iset):
    sel = np.nonzero(mask)[0]
    M = np.stack([Vr[X][sel, k] for X in Iset], axis=1)
    dec = np.all((M == 0) | (M == 1), axis=1)
    neq = dec & np.any(M != M[:, :1], axis=1)
    U = neq.sum() / dec.sum() if dec.sum() else math.nan
    As = []
    for X, Y in itertools.combinations(Iset, 2):
        po, _ = a_kappa(*pair_counts(Vr[X][sel, k], Vr[Y][sel, k]))
        As.append(po)
    amin = min(As) if As and not any(math.isnan(a) for a in As) else math.nan
    return U, amin, As, int(mask.sum() - dec.sum())


# ------------------------------------------------------------------------------------------------ bootstrap (E34), own generator
B = int(P["bootstrap_resamples"])
LEVEL = P["bootstrap_level"]
SEED = int(P["bootstrap_seed"])
DROP_MAX = 0.05


def draws(label, m):
    h = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")
    rng = np.random.default_rng([SEED, h])
    return rng.integers(0, m, size=(B, m))


def counts_of(idx):
    m = idx.shape[1]
    flat = (idx + (np.arange(B) * m)[:, None]).ravel()
    return np.bincount(flat, minlength=B * m).reshape(B, m).astype(float)


def interval(vals):
    """vals: replicate values (nan = undefined).  Returns (lo, hi, ndropped) with nan,nan if more than 5 percent dropped."""
    ok = ~np.isnan(vals)
    nd = int((~ok).sum())
    if nd > DROP_MAX * len(vals) or ok.sum() == 0:
        return math.nan, math.nan, nd
    lo, hi = np.percentile(vals[ok], [100 * (1 - LEVEL) / 2, 100 * (1 - (1 - LEVEL) / 2)])   # method linear = type 7
    return float(lo), float(hi), nd


def safe_div(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)
    return r


def med_rows(vals_sel, idx):
    G = vals_sel[idx]
    if np.isnan(G).any():
        return np.nanmedian(G, axis=1)
    return np.median(G, axis=1)


# ------------------------------------------------------------------------------------------------ hypotheses for one version
def dec_word(ok_supported, decidable):
    return "NOT DECIDABLE" if not decidable else ("SUPPORTED" if ok_supported else "NOT SUPPORTED")


def direction_verdict(X, o, d, rule):
    v1, v2, _, _ = indicator(X, o, BUD[-1], hor1=NDIR, direction=d)
    return v1 if rule == "R1" else v2


def double64_verdict(X, o, rule):
    if X in FOUR:
        v1, v2, _, _ = indicator(X, o, BUD[-1], prec="Double64", hor1=NDIR, direction=DIRS[0])
    else:
        v1, v2, _, _ = indicator(X, o, BUD[-1], prec="Double64")
    return v1 if rule == "R1" else v2


H6CACHE = {}


def h6(ver, rule, ind, sdmask, slmask):
    res = {}
    sd = [int(i) for i in ids[sdmask]]
    sl = [int(i) for i in ids[slmask]]
    dd_ok, dp_ok = len(sd) >= MINSET, len(sl) >= MINSET
    res["n_SDs"], res["n_SLs"] = len(sd), len(sl)
    ddir, dprec = {}, {}
    for X in [x for x in FOUR if x in ind]:
        vs = []
        ndrop = 0
        for o in sd:
            key = ("dir", rule, X, o)
            if key not in H6CACHE:
                H6CACHE[key] = [direction_verdict(X, o, d, rule) for d in DIRS]
            row = H6CACHE[key]
            if any(v not in (0, 1) for v in row):
                ndrop += 1
                continue
            vs.append(len(set(row)) > 1)
        ddir[X] = (sum(vs) / len(vs) if vs else math.nan, ndrop)
    for X in [x for x in I7 if x in ind]:
        vs = []
        ndrop = 0
        for o in sl:
            key = ("prec", rule, X, o)
            if key not in H6CACHE:
                H6CACHE[key] = double64_verdict(X, o, rule)
            v64 = H6CACHE[key]
            v32 = int(V[rule][X][OI(o), NB - 1])
            if v64 not in (0, 1) or v32 not in (0, 1):
                ndrop += 1
                continue
            vs.append(v32 != v64)
        dprec[X] = (sum(vs) / len(vs) if vs else math.nan, ndrop)
    mdir = max([v for v, _ in ddir.values() if not math.isnan(v)], default=math.nan)
    mprec = max([v for v, _ in dprec.values() if not math.isnan(v)], default=math.nan)
    sup = (dd_ok and mdir >= P["h6_flip_min"]) or (dp_ok and mprec >= P["h6_flip_min"])
    dec = "SUPPORTED" if sup else ("NOT SUPPORTED" if (dd_ok and dp_ok) else "NOT DECIDABLE")
    for X, (v, nd) in ddir.items():
        out("%s|H6|ddir|%s" % (ver, X), v)
        out("%s|H6|ddir_drop|%s" % (ver, X), nd)
    for X, (v, nd) in dprec.items():
        out("%s|H6|dprec|%s" % (ver, X), v)
        out("%s|H6|dprec_drop|%s" % (ver, X), nd)
    out("%s|H6|n_SDs" % ver, len(sd))
    out("%s|H6|n_SLs" % ver, len(sl))
    out("%s|H6|max_dir" % ver, mdir)
    out("%s|H6|max_prec" % ver, mprec)
    out("%s|H6|decision" % ver, dec)
    return dec


def run_version(ver, rule, ind, setrule, variant="all"):
    """All of E32 for indicator list ind, verdict rule, on the sets of setrule ('Eall' or 'E8')."""
    sfx = "_8" if setrule == "E8" else ""
    sets = {n: SETS[n + sfx] for n in ("W2", "S2", "Ts", "Ts_design", "SD_star", "SL_star")}
    Vr = V[rule]
    tc, oc, _ = (None, None, None)
    O = OUTC[rule]
    decisions = {}
    pairs = list(itertools.combinations(ind, 2))
    # ---- H1
    W = sets["W2"]
    sel = np.nonzero(W)[0]
    n = len(sel)
    U, amin, As, ndr = U_amin(Vr, W, KP, ind)
    out("%s|H1|n_W2" % ver, n)
    out("%s|H1|U" % ver, U)
    out("%s|H1|U_drop" % ver, ndr)
    out("%s|H1|Amin" % ver, amin)
    for (X, Y), a in zip(pairs, As):
        out("%s|H1|A|%s-%s" % (ver, X, Y), a)
    Us2, _, _, ndr2 = U_amin(Vr, sets["S2"], KP, ind)
    out("%s|H1|U_S2" % ver, Us2)
    if n >= 1:
        idx = draws("p2|H1|W2|%s|%s|%s" % (setrule, rule, variant), n)
        Cn = counts_of(idx)
        M = np.stack([Vr[X][sel, KP] for X in ind], axis=1)
        dec = np.all((M == 0) | (M == 1), axis=1)
        neq = (dec & np.any(M != M[:, :1], axis=1)).astype(float)
        Urep = safe_div(Cn @ neq, Cn @ dec.astype(float))
        lo, hi, nd = interval(Urep)
    else:
        lo = hi = math.nan
        nd = 0
    out("%s|H1|U_lo" % ver, lo)
    out("%s|H1|U_hi" % ver, hi)
    out("%s|H1|U_repdrop" % ver, nd)
    decidable = n >= MINSET and not math.isnan(U) and not math.isnan(lo) and not math.isnan(amin)
    d = dec_word(decidable and U >= P["h1_u_min"] and lo >= P["h1_u_lower_min"] and amin <= P["h1_amin_max"], decidable)
    out("%s|H1|decision" % ver, d)
    decisions["H1"] = d
    # ---- H2
    S = sets["S2"]
    selS = np.nonzero(S)[0]
    nS = len(selS)
    MX = {X: median(O[X][0][selS][~np.isnan(O[X][0][selS])]) for X in ind}
    Q = max(MX.values()) / min(MX.values()) if all(not math.isnan(v) for v in MX.values()) else math.nan
    out("%s|H2|n_S2" % ver, nS)
    for X in ind:
        out("%s|H2|M|%s" % (ver, X), MX[X])
        for oc_ in ("correct", "wrong", "undetected"):
            out("%s|H2|share_%s|%s" % (ver, oc_, X), (O[X][1][selS] == oc_).sum() / max((O[X][1][selS] != "").sum(), 1))
    out("%s|H2|Q" % ver, Q)
    if nS >= 1:
        idx = draws("p2|H2|S2|%s|%s|%s" % (setrule, rule, variant), nS)
        Ms = np.stack([med_rows(O[X][0][selS], idx) for X in ind], axis=1)
        Qrep = Ms.max(axis=1) / Ms.min(axis=1)
        lo, hi, nd = interval(Qrep)
    else:
        lo = hi = math.nan
    out("%s|H2|Q_lo" % ver, lo)
    out("%s|H2|Q_hi" % ver, hi)
    decidable = nS >= MINSET and not math.isnan(Q) and not math.isnan(lo)
    d = dec_word(decidable and Q >= P["h2_q_min"] and lo > P["h2_q_lower_gt"], decidable)
    out("%s|H2|decision" % ver, d)
    decisions["H2"] = d
    # ---- H3
    selW = np.nonzero(W)[0]
    nW = len(selW)
    out("%s|H3|n_W2" % ver, nW)
    out("%s|H3|n_S2" % ver, nS)
    MW = {X: median(O[X][0][selW][~np.isnan(O[X][0][selW])]) for X in ind}
    R = {X: MW[X] / MX[X] if (nW and nS) else math.nan for X in ind}
    if nW >= 1 and nS >= 1:
        iw = draws("p2|H3|W2|%s|%s|%s" % (setrule, rule, variant), nW)
        is_ = draws("p2|H3|S2|%s|%s|%s" % (setrule, rule, variant), nS)
        Rl, Rh = {}, {}
        for X in ind:
            rep = med_rows(O[X][0][selW], iw) / med_rows(O[X][0][selS], is_)
            Rl[X], Rh[X], _ = interval(rep)
    else:
        Rl = Rh = {X: math.nan for X in ind}
    npass = 0
    for X in ind:
        out("%s|H3|MW|%s" % (ver, X), MW[X])
        out("%s|H3|MS|%s" % (ver, X), MX[X])
        out("%s|H3|R|%s" % (ver, X), R[X])
        out("%s|H3|R_lo|%s" % (ver, X), Rl[X])
        out("%s|H3|R_hi|%s" % (ver, X), Rh[X])
        if not math.isnan(R[X]) and not math.isnan(Rl[X]) and R[X] >= P["h3_r_min"] and Rl[X] > P["h3_r_lower_gt"]:
            npass += 1
    out("%s|H3|npass" % ver, npass)
    decidable = nW >= MINSET and nS >= MINSET
    d = dec_word(decidable and npass >= P["h3_count_min"], decidable)
    out("%s|H3|decision" % ver, d)
    decisions["H3"] = d
    # ---- H4 (LE left out), and the design ensemble
    ind4 = [x for x in ind if x != "LE"]
    for item, tsname in (("H4", "Ts"), ("H4design", "Ts_design")):
        T = sets[tsname]
        selT = np.nonzero(T)[0]
        nT = len(selT)
        v9 = {X: Vr[X][selT, NB - 1] for X in ind4}
        phi = {X: ((v9[X] == 0).sum() / ((v9[X] == 0) | (v9[X] == 1)).sum()) if (((v9[X] == 0) | (v9[X] == 1)).sum() > 0) else math.nan for X in ind4}
        dphi = max(phi.values()) - min(phi.values()) if all(not math.isnan(v) for v in phi.values()) else math.nan
        out("%s|%s|n_Ts" % (ver, item), nT)
        for X in ind4:
            out("%s|%s|phi|%s" % (ver, item, X), phi[X])
        out("%s|%s|dphi" % (ver, item), dphi)
        if nT >= 1:
            idx = draws("p2|%s|%s|%s|%s|%s" % (item, tsname, setrule, rule, variant), nT)
            Cn = counts_of(idx)
            reps = []
            for X in ind4:
                num = Cn @ (v9[X] == 0).astype(float)
                den = Cn @ ((v9[X] == 0) | (v9[X] == 1)).astype(float)
                reps.append(safe_div(num, den))
            reps = np.stack(reps, axis=1)
            dr = reps.max(axis=1) - reps.min(axis=1)       # nan propagates when any phi_X is undefined
            lo, hi, nd = interval(dr)
        else:
            lo = hi = math.nan
        out("%s|%s|dphi_lo" % (ver, item), lo)
        out("%s|%s|dphi_hi" % (ver, item), hi)
        decidable = nT >= MINSET and not math.isnan(dphi) and not math.isnan(lo)
        d = dec_word(decidable and dphi >= P["h4_dphi_min"] and lo >= P["h4_dphi_lower_min"], decidable)
        out("%s|%s|decision" % (ver, item), d)
        if item == "H4":
            decisions["H4"] = d
    # ---- H5
    out("%s|H5|n_W2" % ver, nW)
    F = {}
    for X in ind:
        a, b_ = Vr[X][selW, KP], Vr[X][selW, NB - 1]
        both = ((a == 0) | (a == 1)) & ((b_ == 0) | (b_ == 1))
        F[X] = ((a != b_) & both).sum() / both.sum() if both.sum() else math.nan
        out("%s|H5|F|%s" % (ver, X), F[X])
    Fmed = median([F[X] for X in ind]) if all(not math.isnan(v) for v in F.values()) else math.nan
    out("%s|H5|Fmed" % ver, Fmed)
    if nW >= 1:
        idx = draws("p2|H5|W2|%s|%s|%s" % (setrule, rule, variant), nW)
        Cn = counts_of(idx)
        reps = []
        for X in ind:
            a, b_ = Vr[X][selW, KP], Vr[X][selW, NB - 1]
            both = (((a == 0) | (a == 1)) & ((b_ == 0) | (b_ == 1))).astype(float)
            flip = (both.astype(bool) & (a != b_)).astype(float)
            reps.append(safe_div(Cn @ flip, Cn @ both))
        reps = np.stack(reps, axis=1)
        fm = np.where(np.isnan(reps).any(axis=1), np.nan, np.median(reps, axis=1))
        lo, hi, nd = interval(fm)
    else:
        lo = hi = math.nan
    out("%s|H5|Fmed_lo" % ver, lo)
    out("%s|H5|Fmed_hi" % ver, hi)
    decidable = nW >= MINSET and not math.isnan(Fmed) and not math.isnan(lo)
    d = dec_word(decidable and Fmed >= P["h5_fmed_min"] and lo >= P["h5_fmed_lower_min"], decidable)
    out("%s|H5|decision" % ver, d)
    decisions["H5"] = d
    # ---- H6
    decisions["H6"] = h6(ver, rule, ind, sets["SD_star"], sets["SL_star"])
    return decisions


# ------------------------------------------------------------------------------------------------ gate (E33)
say("gate")
gate_viol = []
gsel = np.nonzero(GATE)[0]
for X in I8:
    v = V["R1"][X][gsel, NB - 1]
    dec = ((v == 0) | (v == 1)).sum()
    ch = (v == 1).sum()
    share = ch / dec if dec else math.nan
    viol = bool(dec and share > P["gate_max_share"])
    out("gate|%s|decided" % X, int(dec))
    out("gate|%s|chaotic" % X, int(ch))
    out("gate|%s|undecided" % X, int((v == -1).sum()))
    out("gate|%s|na" % X, int((v == 2).sum()))
    out("gate|%s|share" % X, share)
    out("gate|%s|violated" % X, viol)
    if viol:
        gate_viol.append(X)
        out("gate|%s|message" % X, "rule invalid on the integrable control")

# ------------------------------------------------------------------------------------------------ run all versions
DEC = {}
for rule in ("R1", "R2"):
    say("version", rule)
    DEC[(rule, "I7")] = run_version("%s.I7" % rule, rule, I7, "Eall")
    DEC[(rule, "I8")] = run_version("%s.I8" % rule, rule, I8, "E8")
    for h in ("H1", "H2", "H3", "H4", "H5", "H6"):
        a, b_ = DEC[(rule, "I7")][h], DEC[(rule, "I8")][h]
        both = a != "NOT DECIDABLE" and b_ != "NOT DECIDABLE"
        out("%s|robust|%s" % (rule, h), ("ROBUST" if a == b_ else "NOT ROBUST") if both else "NOT COMPARABLE")
for X in gate_viol:
    say("variant minus", X)
    run_version("R1.I7.minus-%s" % X, "R1", [x for x in I7 if x != X], "Eall", "minus-" + X)
    run_version("R1.I8.minus-%s" % X, "R1", [x for x in I8 if x != X], "E8", "minus-" + X)

# ------------------------------------------------------------------------------------------------ reported items
say("reported items")
SETN = {"W2": W2, "S2": S2, "Ts": Ts, "regular": REG}
for rule in ("R1", "R2"):
    for sn, m in SETN.items():
        for X in I8:
            for k in range(NB):
                v = V[rule][X][m, k]
                for nm, code in (("chaotic", 1), ("regular", 0), ("undecided", -1), ("na", 2)):
                    out("vcount|%s|%s|%s|k%d|%s" % (rule, sn, X, k, nm), int((v == code).sum()))
    O = OUTC[rule]
    for sn in ("W2", "S2", "Ts", "regular", "Ts_design"):
        m = SETS[sn]
        sel = np.nonzero(m)[0]
        for X in I8:
            oc = O[X][1][sel]
            for nm in ("correct", "wrong", "undetected"):
                out("det|%s|%s|%s|%s" % (rule, sn, X, nm), int((oc == nm).sum()))
            tcs = O[X][0][sel]
            tcs = tcs[~np.isnan(tcs)]
            out("det|%s|%s|%s|M" % (rule, sn, X), median(tcs))
            out("det|%s|%s|%s|n" % (rule, sn, X), len(tcs))
            tdd = O[X][2][sel]
            for val in sorted(set(tdd[~np.isnan(tdd)])):
                out("tdet|%s|%s|%s|%d" % (rule, sn, X, int(val)), int((tdd == val).sum()))
            # flip count distribution (R1 only to keep the file short)
for sn in ("W2", "S2", "Ts", "regular"):
    sel = np.nonzero(SETN[sn])[0]
    for X in I8:
        Vx = V["R1"][X][sel]
        phi_ = np.zeros(len(sel), dtype=int)
        for k in range(1, NB):
            phi_ += (((Vx[:, k] == 0) | (Vx[:, k] == 1)) & ((Vx[:, k - 1] == 0) | (Vx[:, k - 1] == 1)) & (Vx[:, k] != Vx[:, k - 1]))
        for cnt_v in sorted(set(phi_)):
            out("Phi|%s|%s|%d" % (sn, X, cnt_v), int((phi_ == cnt_v).sum()))
        for k in range(NB):
            a, b_ = Vx[:, k], Vx[:, NB - 1]
            both = ((a == 0) | (a == 1)) & ((b_ == 0) | (b_ == 1))
            out("Fb|%s|k%d|%s" % (sn, k, X), ((a != b_) & both).sum() / both.sum() if both.sum() else math.nan)

# F1: agreement matrices for I_8 at b_P and N_c on W2 and S2
for sn in ("W2", "S2"):
    for k in (KP, NB - 1):
        sel = np.nonzero(SETN[sn])[0]
        kap = []
        for X, Y in itertools.combinations(I8, 2):
            ncc, ncr, nrc, nrr = pair_counts(V["R1"][X][sel, k], V["R1"][Y][sel, k])
            po, kp = a_kappa(ncc, ncr, nrc, nrr)
            pre = "agree|%s|k%d|%s-%s" % (sn, k, X, Y)
            out(pre + "|A", po); out(pre + "|kappa", kp); out(pre + "|N", ncc + ncr + nrc + nrr)
            out(pre + "|ncc", ncc); out(pre + "|ncr", ncr); out(pre + "|nrc", nrc); out(pre + "|nrr", nrr)
            sp = spearman_corr(SC[X][sel, k], SC[Y][sel, k])
            out("spearman|%s|k%d|%s-%s" % (sn, k, X, Y), sp)
            if X in I7 and Y in I7:
                kap.append(kp)
        out("meankappa|%s|k%d|I7" % (sn, k), float(np.nanmean(kap)))

# T9: AUC per K value (pooled, then by source when both counts are large enough)
def auc(sc_c, sc_r):
    nC, nR = len(sc_c), len(sc_r)
    rk = rankdata(np.concatenate([sc_c, sc_r]))
    return (rk[:nC].sum() - nC * (nC + 1) / 2) / (nC * nR)


Kall = sorted(set(Kval[~np.isnan(Kval)]))
for X in IND:
    for k in range(NB):
        sc = SC[X][:, k]
        for K in Kall:
            for src, sm in (("pooled", np.ones(NORB, bool)), ("P", ~isG), ("G", isG)):
                mc = Cm & (Kval == K) & sm & ~np.isnan(sc)
                mr = Rm & (Kval == K) & sm & ~np.isnan(sc)
                if mc.sum() >= MINSET and mr.sum() >= MINSET:
                    pre = "auc|%s|k%d|K=%s|%s" % (X, k, kv(K), src)
                    out(pre + "|auc", auc(sc[mc], sc[mr])); out(pre + "|nC", int(mc.sum())); out(pre + "|nR", int(mr.sum()))

# T14: regular range, leaving budget, share below the regular median
def q7(sorted_vals, p):
    n = len(sorted_vals)
    h = (n - 1) * p / 100.0
    lo = int(math.floor(h))
    fr = h - lo
    if fr == 0 or lo + 1 >= n:
        return float(sorted_vals[lo])
    a, b_ = sorted_vals[lo], sorted_vals[lo + 1]
    if a == b_:
        return float(a)
    if a == -math.inf:
        return -math.inf
    if b_ == math.inf:
        return math.inf
    return float(a + fr * (b_ - a))


for X in I7 + ["STR"]:
    for K in Kall:
        mts = Ts & (Kval == K)
        mreg = Rm & (Kval == K)
        rng_ = {}
        for k in range(NB):
            av = np.sort(SC[X][mreg, k][~np.isnan(SC[X][mreg, k])])
            if len(av) >= MINSET:
                rng_[k] = (q7(av, P["range_lower_percentile"]), q7(av, P["range_upper_percentile"]), float(median(av)), len(av))
                pre = "rng|%s|k%d|K=%s" % (X, k, kv(K))
                out(pre + "|lo", rng_[k][0]); out(pre + "|hi", rng_[k][1]); out(pre + "|median", rng_[k][2]); out(pre + "|n", len(av))
        if mts.sum() == 0 or len(rng_) < NB:
            if mts.sum():
                out("leave|%s|K=%s|n" % (X, kv(K)), int(mts.sum())); out("leave|%s|K=%s|defined" % (X, kv(K)), "no range")
            continue
        n_def, below, n9 = 0, 0, 0
        lv = defaultdict(int)
        for i in np.nonzero(mts)[0]:
            s = SC[X][i]
            above = [(not math.isnan(s[k])) and s[k] > rng_[k][1] for k in range(NB)]
            if above[NB - 1]:
                r_ = NB - 1
                while r_ > 0 and above[r_ - 1]:
                    r_ -= 1
                lv[BUD[r_]] += 1
                n_def += 1
            if not math.isnan(s[NB - 1]):
                n9 += 1
                below += s[NB - 1] < rng_[NB - 1][2]
        out("leave|%s|K=%s|n" % (X, kv(K)), int(mts.sum())); out("leave|%s|K=%s|defined" % (X, kv(K)), n_def)
        out("leave|%s|K=%s|below_median" % (X, kv(K)), below); out("leave|%s|K=%s|n9" % (X, kv(K)), n9)
        for b_, c_ in sorted(lv.items()):
            out("leavek|%s|K=%s|%d" % (X, kv(K), b_), c_)

# T11: sensitivities on S_L*
def verdict_from_sens(spec, o, b, X):
    d = WIN.get((o, spec, "Float64"))
    if d is None:
        return -2, math.nan
    W = b // 2
    r = rowindex(d["rows"], W)
    D = d["c"]["D"][r]
    bad = bit(d["flags"][r], d["bits"], "degenerate" if X == "E1" else "undefined window")
    thr = P["e1_d_threshold"] if X == "E1" else P["e2_d_threshold"]
    if bad or not math.isfinite(D):
        return -1, D
    return (1 if D >= thr else 0), D


sl_ids = [int(i) for i in ids[SLs]]
for pair, (spec, base, X) in {"fq1_raw-fq1": ("fq1_raw", "fq1", "E1"), "fq1_qp-fq1": ("fq1_qp", "fq1", "E1"), "pa_planar-pa1": ("pa_planar", "pa1", "E2")}.items():
    for k, b in enumerate(BUD):
        va, vb, diffs = [], [], []
        for o in sl_ids:
            a, Da = verdict_from_sens(spec, o, b, X)
            c_, Db = verdict_from_sens(base, o, b, X)
            va.append(a); vb.append(c_)
            if math.isfinite(Da) and math.isfinite(Db):
                diffs.append(abs(Da - Db))
        ncc, ncr, nrc, nrr = pair_counts(np.array(va), np.array(vb))
        po, kp = a_kappa(ncc, ncr, nrc, nrr)
        pre = "T11|%s|k%d" % (pair, k)
        out(pre + "|n", ncc + ncr + nrc + nrr); out(pre + "|A", po); out(pre + "|kappa", kp)
        out(pre + "|medabs", median(diffs)); out(pre + "|maxabs", max(diffs) if diffs else math.nan)
for k, b in enumerate(BUD):
    va, vb, diffs = [], [], []
    for o in sl_ids:
        d = WIN[(o, "bk1", "Float64")]
        s1 = win_wba(d, b, single="cos(1 x[1])")
        s4 = win_wba(d, b)
        va.append(s1[0]); vb.append(s4[0])
        if math.isfinite(s1[2]) and math.isfinite(s4[2]):
            diffs.append(abs(s1[2] - s4[2]))
    ncc, ncr, nrc, nrr = pair_counts(np.array(va), np.array(vb))
    po, kp = a_kappa(ncc, ncr, nrc, nrr)
    pre = "T11|wba_cos1-wba|k%d" % k
    out(pre + "|n", ncc + ncr + nrc + nrr); out(pre + "|A", po); out(pre + "|kappa", kp)
    out(pre + "|medabs", median(diffs)); out(pre + "|maxabs", max(diffs) if diffs else math.nan)

# T15: packages
for pname, (X, Y), em in (("FLI+RLI", ("FLI", "RLI"), True), ("MEGNO+SALI", ("MEGNO", "SALI"), False)):
    for sn in ("W2", "S2", "Ts", "regular"):
        m = SETN[sn] & (inE8 if em else True)
        for k in (KP, NB - 1):
            c_ = defaultdict(int)
            for i in np.nonzero(m)[0]:
                c_[(int(V["R1"][X][i, k]), int(V["R1"][Y][i, k]), rho[i])] += 1
            for (a, b_, r_), n_ in sorted(c_.items()):
                out("T15|%s|%s|k%d|%d,%d|rho=%s" % (pname, sn, k, a, b_, r_), n_)

# T16: K = 0 line and control verdicts
for spec in ("fq1", "pa1"):
    for W in P["windows"]:
        mx = 0.0
        for i in ids[GATE]:
            d = WIN[(int(i), spec, "Float64")]
            nu = d["c"]["nu1"][rowindex(d["rows"], W)]
            p0 = float(ORBS[int(i)]["ic"].split(",")[1])
            if math.isfinite(nu):
                dd = abs((nu - p0 / (2 * math.pi) + 0.5) % 1.0 - 0.5)
                mx = max(mx, dd)
        out("k0line|%s|W=%d|maxdist" % (spec, W), mx)
for i in ids[ctrl]:
    for X in I8:
        out("ctrl|%d|%s|v9" % (i, X), int(V["R1"][X][OI(i), NB - 1]))
    out("ctrl|%d|lamref" % i, lamref[OI(i)])
# T16 flag counts per chunk against the manifest (recomputed by the reader) and omissions
nbad = sum(1 for d in CHF if {k: v for k, v in d["recomputed"].items() if v} != {k: v for k, v in d["manifest"].items() if v})
out("T16|flagchunks_differing", nbad)
out("T16|flagchunks", len(CHF))
for d in NA:
    out("T16|omitted|%s" % d["spec"], ",".join(str(x) for x in d["recomputed"]))
    out("T16|omitted_equals_manifest|%s" % d["spec"], d["recomputed"] == d["manifest"])
tot = defaultdict(int)
for d in CHF:
    for k, v in d["recomputed"].items():
        tot["%s|%s" % (d["spec"], k)] += v
for k, v in sorted(tot.items()):
    out("T16|flagsum|" + k, v)

# ------------------------------------------------------------------------------------------------ self-checks O7, O8
say("self checks")
SELF = []


def selfcheck(name, ok, detail=""):
    SELF.append((name, ok))
    out("selfcheck|" + name, ("PASS" if ok else "FAIL") + (" " + detail if detail else ""))
    print(("PASS " if ok else "FAIL "), name, detail)


# O7 (a) kappa
sc_ = lambda t_: a_kappa(*t_)[1]
selfcheck("O7a kappa perfect agreement = 1", abs(sc_((30, 0, 0, 70)) - 1) < 1e-15)
selfcheck("O7a kappa independent [[25,25],[25,25]] = 0", abs(sc_((25, 25, 25, 25))) < 1e-15)
selfcheck("O7a kappa p_o = p_e (equal margin products) = 0", abs(sc_((20, 30, 20, 30))) < 1e-15)   # margins 50/50 and 40/60: p_o = 0.5 = p_e
selfcheck("O7a kappa with p_e = 1 (one verdict only) = 0", sc_((40, 0, 0, 0)) == 0.0 and sc_((0, 0, 0, 40)) == 0.0)
selfcheck("O7a kappa a = d = 0 and b = c = -1", abs(sc_((0, 25, 25, 0)) + 1) < 1e-15)
selfcheck("O7a kappa N = 0 undefined", math.isnan(sc_((0, 0, 0, 0))))
# O7 (b) AUC
selfcheck("O7b AUC separated classes = 1", auc(np.array([3., 4., 5.]), np.array([0., 1.])) == 1.0)
selfcheck("O7b AUC identical score sets = 0.5", auc(np.array([1., 2., 3.]), np.array([1., 2., 3.])) == 0.5)
selfcheck("O7b AUC C={1,2,3} against R={2} = 0.5", auc(np.array([1., 2., 3.]), np.array([2.])) == 0.5)
selfcheck("O7b AUC one tied pair counts one half: C={2,3}, R={2} = 0.75", auc(np.array([2., 3.]), np.array([2.])) == 0.75)
selfcheck("O7b AUC of a score and of its negation sum to 1 (with -inf)", abs(auc(np.array([-math.inf, 1., 2.]), np.array([-math.inf, 0.])) + auc(np.array([math.inf, -1., -2.]), np.array([math.inf, 0.])) - 1) < 1e-15)
rs = np.random.default_rng(1)
okp = True
for _ in range(200):
    c1, r1 = rs.integers(0, 6, 9).astype(float), rs.integers(0, 6, 7).astype(float)
    brute = np.mean([1.0 if x > y else (0.5 if x == y else 0.0) for x in c1 for y in r1])
    okp &= abs(brute - auc(c1, r1)) < 1e-14
selfcheck("O7b AUC rank-sum equals pair counting on 200 random tie-rich samples", bool(okp))

# O7 (c) exact bootstrap: six orbits, all 6^6 resamples of U (three indicators) and Q; compare with the generator by the CDF
def exact_check():
    rs_ = np.random.default_rng(5)
    Vt = rs_.integers(0, 2, size=(6, 3))                 # verdicts of 3 indicators on 6 orbits, all decided
    tt = rs_.choice([200.0, 400.0, 51200.0, 204800.0], size=(6, 3))
    allidx = np.array(list(itertools.product(range(6), repeat=6)))
    neq = (np.any(Vt != Vt[:, :1], axis=1)).astype(float)
    Uex = neq[allidx].mean(axis=1)
    Mex = np.median(tt[allidx], axis=1)                  # (46656, 3) medians per indicator
    Qex = Mex.max(axis=1) / Mex.min(axis=1)
    idx = draws("p2|selfcheck|exact|E_all|R1|all", 6)
    Urep = neq[idx].mean(axis=1)
    Mr = np.stack([np.median(tt[:, j][idx], axis=1) for j in range(3)], axis=1)
    Qrep = Mr.max(axis=1) / Mr.min(axis=1)
    ok_ = True
    for ex, rep in ((Uex, Urep), (Qex, Qrep)):
        for p in (0.025, 0.5, 0.975):
            q = np.quantile(ex, p, method="lower")
            Fex = (ex <= q).mean()
            Fem = (rep <= q).mean()
            ok_ &= abs(Fem - Fex) <= 4 * math.sqrt(Fex * (1 - Fex) / B + 1e-300) + 1e-12
    return bool(ok_)


selfcheck("O7c exact bootstrap distribution of U and Q (6 orbits, 6^6 resamples) matches the generator within 4 standard errors", exact_check())

# Spearman: own rank correlation against scipy.stats.spearmanr (W2, b_P, LE against SALI; STR is not in I_8)
_sel = np.nonzero(W2)[0]
_x, _y = SC["LE"][_sel, KP], SC["FLI"][_sel, KP]
_m = ~(np.isnan(_x) | np.isnan(_y))
selfcheck("Spearman midrank correlation equals scipy.stats.spearmanr (W2, b_P, LE against FLI)", abs(spearman_corr(_x, _y) - float(spearmanr(_x[_m], _y[_m])[0])) < 1e-12)

# O8: pilot cross-check, 24 starts x 8 indicators = 192 R1 verdicts from the tabulated scores
txt = open(PILOT, encoding="utf-8").read()
sec = txt[txt.index("=== budget n = 100000"):]
sec = sec[:sec.index("```")]
rows = []
for ln in sec.splitlines():
    m = re.match(r"^\s*(\d+)\s+(\d+)\s+(.*)\|(.*)$", ln)
    if not m:
        continue
    nums = m.group(3).split()
    vs = m.group(4).split()
    if len(nums) == 10 and len(vs) == 8:
        rows.append((int(m.group(2)), [float(x) for x in nums], vs))
selfcheck("O8 pilot table parsed: 24 starts", len(rows) == 24, str(len(rows)))
nok = ntot = 0
bad = []
for orb, nums, vs in rows:
    lam1, rLE, rS, rF, rM, sali, rli, de1, de2, dig = nums
    tab = {"LE": (lam1 >= P["le_floor"] and rLE >= P["ratio_threshold"]), "SALI": rS >= P["ratio_threshold"], "FLI": rF >= P["ratio_threshold"],
           "MEGNO": rM >= P["ratio_threshold"], "RLI": rli >= P["rli_threshold"], "E1": de1 >= P["e1_d_threshold"], "E2": de2 >= P["e2_d_threshold"],
           "WBA": dig < P["wba_digits_threshold"]}
    for X, tv in zip(["LE", "SALI", "FLI", "MEGNO", "RLI", "E1", "E2", "WBA"], vs):   # column order of the pilot table
        ntot += 1
        if tab[X] == (tv == "C"):
            nok += 1
        else:
            bad.append((orb, X))
selfcheck("O8 pilot: R1 verdicts from the tabulated scores equal the tabulated verdicts", nok == 192 and ntot == 192, "%d of %d%s" % (nok, ntot, " bad " + str(bad) if bad else ""))
out("selfcheck_summary", "%d PASS of %d" % (sum(1 for _, ok in SELF if ok), len(SELF)))

# O9 results from the reader
for d in META:
    if "check" in d:
        out("reader|" + d["check"], ("PASS " if d["ok"] else "FAIL ") + d["detail"])
    else:
        out("reader|counts", json.dumps(d))

# ------------------------------------------------------------------------------------------------ write
with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
    for n_, v_ in RES:
        fh.write("%s | %s\n" % (n_, v_))
say("wrote", OUT, "lines", len(RES))
