# Independent recomputation of the test side of study p3-abstention (oracle O20 of the analysis spec, Amendment D95), adapted to the store of the runner (Amendment D97, O28, O29).
# This file supersedes verify/analysis-p3t-oracle.py for a store written by the runner (the D95 file is kept for the toy store of D95 and is not changed). Changes of D97: the paths
# (fixture layout of the D96 campaign), cell k<i>:<class>:r<round> (E74), six arms (E81), tg at N_c and pa1 at N_c + 1 in the pa1 arms (E75), manifest status (E76), chunks failed twice by
# arm and cell (E77), costs from the manifest specs of ind and ind_pa1 and L1 from the ref chunk entries (E78), T3 by content at the chk horizon without the "labels agree" column (E79),
# counts over the rounds up to the closing round (E80). The sources of the numbers are the chunk files (verify/analysis-p3t-export.jl, unchanged), the panel, the manifest and the parameter tables.
# Python standard library plus numpy; exact fractions.Fraction for every count, share, p-value, bootstrap statistic and percentile.
# Written from analysis-spec.md (Amendment D95 E59 to E73; D93 E47 to E57; D89 E23 to E31, E34), experiments-spec.md (D94 E39 to E43, R23), schema/run-store.md
# and the toy-campaign INPUTS.md. No file of src/, test/, bench/, ext/ and no module output was read before this oracle wrote its result.
# The generators, the verdict functions, the cost model, AUC, McNemar and the power simulators are reused from verify/analysis-p3-oracle.py
# (written independently of the module for D93).
#
# Usage (repository root; <tag> names store-<tag>, manifest-<tag>.json and round-log-<tag>.txt of <out dir>):
#   julia --project=verify verify/analysis-p3t-export.jl <out dir>/store-<tag> <work folder>
#   python -I verify/analysis-d97-oracle.py <fixture dir> <out dir> <tag> <work folder> <result json>
# Places where the spec leaves a choice are marked READING and collected in the result under "readings".
import sys, os, json, math, hashlib, re, importlib.util, time
from fractions import Fraction
import numpy as np
import tomllib

T0 = time.time()
def say(*a):
    print("[%6.1f s]" % (time.time() - T0), *a, flush=True)

HERE = os.path.dirname(os.path.abspath(__file__))
_sp = importlib.util.spec_from_file_location("p3o", os.path.join(HERE, "analysis-p3-oracle.py"))
p3o = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(p3o)
Xoshiro, XoV, stream_seed, ppnd16 = p3o.Xoshiro, p3o.XoV, p3o.stream_seed, p3o.ppnd16
fr, rtext, mcnemar_tail, mcnemar_critical_table, mw_auc, Costs = p3o.fr, p3o.rtext, p3o.mcnemar_tail, p3o.mcnemar_critical_table, p3o.mw_auc, p3o.Costs
series_verdict, window_verdict = p3o.series_verdict, p3o.window_verdict
CH, RG, UN = p3o.CH, p3o.RG, p3o.UN
I7 = p3o.I7
READINGS = []
def reading(t):
    if t not in READINGS:
        READINGS.append(t)

FIXD, OUTD, TAG, WORK, RESULT = sys.argv[1:6]
CP = lambda *a: os.path.join(FIXD, *a)
OP = lambda *a: os.path.join(OUTD, *a)
for _p in (WORK, RESULT):
    assert "/data/" not in os.path.abspath(_p).replace("\\", "/"), "refused: inside data/"
sha_file = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
sha_txt = lambda b: hashlib.sha256(b).hexdigest()

# ------------------------------------------------------------------ inputs
PREREG = CP("inputs", "p3-abstention-preregistration.sealed.md")
q = {k: v["value"] for k, v in tomllib.load(open(CP("inputs", "p3-abstention-calibration-parameters.toml"), "rb")).items()}
qt = {k: v["value"] for k, v in tomllib.load(open(CP("inputs", "p3-abstention-test-parameters.toml"), "rb")).items()}
p2 = {k: v["value"] for k, v in tomllib.load(open(CP("p2", "p2-indicators-analysis-parameters.toml"), "rb")).items()}
p1 = tomllib.load(open(CP("p1", "p1-le-benchmark-parameters.toml"), "rb"))
PANELF = CP("p3-abstention-test-panel.toml")
panel = tomllib.load(open(PANELF, "rb"))
LOGB = open(OP("round-log-%s.txt" % TAG), "rb").read()
RC = LOGB.decode("utf-8").count(chr(10)) // 13          # the closing round: the number of whole entries of the closed log (E62 a, E80)
ORB_ALL = {o["id"]: o for o in panel["orbit"]}
ORB = {i: o for i, o in ORB_ALL.items() if o["round"] <= RC}     # E80: the orbits of rounds above the closing round are not counted anywhere
CELLORB = {}                                                      # E74: the orbits of a cell are those of the panel with its (k_index, class, round)
for _o in ORB.values(): CELLORB.setdefault((_o["k_index"], _o["class"], _o["round"]), []).append(_o["id"])
MANF = OP("manifest-%s.json" % TAG)
MAN = json.load(open(MANF))
REC = json.load(open(CP("inputs", "p3-abstention-calibration-record.json")))
LOGF = OP("round-log-%s.txt" % TAG)
checks = []
def check(name, ok, detail=""):
    checks.append(dict(name=name, ok=bool(ok), detail=str(detail)))
    say(("PASS " if ok else "FAIL ") + name, detail)

N1, NC, NS = q["N_1"], q["N_c"], q["n_s"]
B = list(q["budgets"]); bP = q["b_primary"]; kP = B.index(bP)
RMAX = int(q["r_max"])
TSEED = q["test_seed"]
VS = list(qt["theta_voters"]); QC, QR, WD, VV = qt["theta_q_c"], qt["theta_q_r"], qt["theta_w"], qt["theta_v"]
T_T, T_W, T_F = qt["target_t_s"], qt["target_w2"], qt["target_f"]
thr = {k: fr(q[k]) for k in ("alpha_hb", "rho_hb", "c_min", "delta_c", "alpha_ci", "delta_a", "e_max", "a_max", "p_hb", "p_hc")}

# preregistration hash (SELF rule)
ptxt = open(PREREG, encoding="utf-8", newline="").read().replace("\r\n", "\n")
m = re.search(r"^(\| SHA-256 of this file \(current version\) \|)(.*)(\|)[ \t]*$", ptxt, re.M)
HEADER_HASH = re.search(r"[0-9a-f]{64}", m.group(2)).group()
self_txt = ptxt[:m.start(2)] + " SELF " + ptxt[m.end(2):]
check("E59 preregistration SELF hash equals the header cell", sha_txt(self_txt.encode("utf-8")) == HEADER_HASH, HEADER_HASH)
PANEL_SHA = sha_file(PANELF)

# ------------------------------------------------------------------ export
chunks, items = [], []
for ln in open(os.path.join(WORK, "p3t-export.jsonl"), encoding="utf-8"):
    d = json.loads(ln)
    (chunks if d["k"] == "chunk" else items).append(d)
CHK = {c["fn"]: c for c in chunks}
IT = {}
dup = 0
for it in items:
    c = CHK[it["fn"]]
    if c["status"] == "failed": continue          # E77: a chunk failed twice holds no item read; it is recognised by its arm and cell
    key = (c["arms"][0], it["o"], it["m"], it["prec"])
    if key in IT: dup += 1
    IT[key] = it
check("E61 f no two chunks with the same arm, orbit and item", dup == 0, dup)
check("E61 f every chunk file complete, study is the slug, panel and preregistration hashes equal the files", all(c["complete"] and c["study"] == "p3-abstention" and c["panel_sha256"] == PANEL_SHA and c["prereg_sha256"] == HEADER_HASH for c in chunks),
      sum(1 for c in chunks if not (c["complete"] and c["study"] == "p3-abstention" and c["panel_sha256"] == PANEL_SHA and c["prereg_sha256"] == HEADER_HASH)))
ARMS = ["ref", "ind", "ind_pa1", "d64", "d64_pa1", "chk"]
check("E81/E61 i the six arm names are the keys of the manifest 'arms' object and each has chunks", set(MAN["arms"]) == set(ARMS) and all(any(a in c["arms"] for c in chunks) for a in ARMS), sorted(MAN["arms"]))
reading("E81: the arm names are the keys of the manifest 'arms' object; the arm of an item is the arm of its chunk; chunk cell = text after 'cell=' of the root attribute chunk_key; the round of a ref chunk is the suffix :r<round> of the cell")
cellre = re.compile(r"cell=k(\d+):([A-Z]):r(\d+)$")
REFCELL = {}      # (round, k_index, class, method) -> chunk
for c in chunks:
    mm = cellre.search(c["chunk_key"])
    if mm is None: check("E74 the cell of chunk %s is of the form k<i>:<class>:r<round>" % c["fn"], False, c["chunk_key"][-24:]); continue
    c["cell"] = (int(mm.group(1)), mm.group(2), int(mm.group(3)))
    if c["arms"] == ["ref"]: REFCELL[(c["cell"][2], c["cell"][0], c["cell"][1], c["method"])] = c
def ref_round(c):
    return c["cell"][2]
check("E74 every chunk key ends in a cell of the form k<i>:<class>:r<round>", all("cell" in c for c in chunks), len(chunks))
check("E80 no chunk of a round above the closing round (the runner never released one)", all(c["cell"][2] <= RC for c in chunks), RC)
# manifest chunk hashes equal the chunk file hashes
mh = {ch["key"]: ch["hash"] for ch in MAN["chunks"]}
check("E61 i manifest chunk hashes equal those of the chunk files read", all(mh.get(c["chunk_key"]) == c["chunk_hash"] for c in chunks) and len(mh) == len(chunks), len(chunks))
nfailed_man = sum(1 for ch in MAN["chunks"] if ch["status"] == "failed")
check("E76 manifest study, and status complete (no failed entry) or complete with failed chunks (at least one failed entry)", MAN["study"] == "p3-abstention" and ((MAN["status"] == "complete" and nfailed_man == 0) or (MAN["status"] == "complete with failed chunks" and nfailed_man > 0)), "%s, %d failed entries" % (MAN["status"], nfailed_man))
check("E76 the set of failed chunks of the manifest equals that of the chunk files", sorted(ch["key"] for ch in MAN["chunks"] if ch["status"] == "failed") == sorted(c["chunk_key"] for c in chunks if c["status"] == "failed"), nfailed_man)
check("E61 e panel hash in manifest equals panel file", MAN["panel"]["sha256_start"] == PANEL_SHA == MAN["panel"]["sha256_end"], PANEL_SHA)
check("E59 manifest preregistration hash equals header", MAN["prereg"]["header"] == HEADER_HASH, "")

# panel identity (E61 j)
badid = 0
for (arm, o, mth, prec), it in IT.items():
    po = ORB[o]
    if it["ic"] != po["x0_dec"] + "," + po["p0_dec"] or float(it["params"]) != float(po["K"]):
        badid += 1
check("E61 j identity: stored ic and params equal the panel file for every item read", badid == 0, badid)

# ------------------------------------------------------------------ labels (E64 a) and per-orbit ref data
def item(arm, o, m, prec="Float64"):
    return IT.get((arm, o, m, prec))
def has_flag(it, f):
    return it is not None and f in it["flags"]

POS, NONE, UNRES = p3o.POS, p3o.NONE, p3o.UNRES
def ref_info(o):
    """returns dict: usable (bool), reason, rho, lam_ref, lam_ns, drift, diag"""
    po = ORB[o]
    l1, l5 = item("ref", o, "L1"), item("ref", o, "L5")
    cell1 = REFCELL.get((po["round"], po["k_index"], po["class"], "L1")); cell5 = REFCELL.get((po["round"], po["k_index"], po["class"], "L5"))
    d = dict(o=o, usable=False, reason=None, rho=None, lam_ref=None, lam_ns=None, drift=False)
    if l1 is None or l5 is None:
        failed = any(c is not None and c["status"] == "failed" for c in (cell1, cell5))
        d["reason"] = "excluded, failed twice" if failed else "missing item"
        return d
    if l1["status"] != "ok" or l5["status"] != "ok" or "excluded" in l1["flags"] + l5["flags"]:
        d["reason"] = "excluded, failed twice"; return d
    need5 = [str(N1 // 10), str(N1 // 2), str(N1)]
    miss = [s for s in need5 if s not in l5["lamfull"]] + (["L1@N1"] if str(N1) not in l1["lamfull"] else []) + (["L1@ns"] if str(NS) not in l1["lamfull"] else [])
    if miss:
        d["reason"] = "missing step " + ",".join(miss); return d
    d["usable"] = True
    d["drift"] = ("drift" in l1["flags"]) or ("drift" in l5["flags"])
    diag = ("diagnostics exceeded" in l1["flags"]) or ("diagnostics exceeded" in l5["flags"])
    lam5, lam10 = l5["lamfull"][str(N1)], l5["lamfull"][str(N1 // 10)]
    lam1 = l1["lamfull"][str(N1)]
    inst = UNRES
    if math.isfinite(lam5) and math.isfinite(lam10) and lam10 > 0:
        r = 10 * lam5 / lam10
        if lam5 >= p2["label_lambda_min"] and r >= p2["label_growth_ratio"]:
            inst = POS
        elif lam5 < p2["label_lambda_min"] and r < p2["label_growth_ratio"]:
            inst = NONE
    x0 = float(l5["ic"].split(",")[0]); h = N1 // 2
    xh, xN = l5["xl"][str(h)], l5["xl"][str(N1)]
    rho1 = (xh - x0) / (2 * math.pi * h); rho2 = (xN - xh) / (2 * math.pi * (N1 - h))
    if inst == POS: st = "irregular"
    elif inst == NONE and abs(rho1 - rho2) <= p2["label_rotation_tol"]: st = "quasiperiodic"
    else: st = UNRES
    if abs(lam1 - lam5) > p2["label_estimator_bound"]: inst, st = UNRES, UNRES
    if diag: inst, st = UNRES, UNRES
    d["rho"] = "C" if (inst == POS and st == "irregular") else ("R" if (inst == NONE and st == "quasiperiodic") else "U")
    d["lam_ref"] = lam1; d["lam_ns"] = l1["lamfull"][str(NS)]
    d["axes"] = (inst, st)
    return d

def in_T(lr, ln, mode="exact"):
    if mode == "exact": return Fraction(ln) <= fr(q["f_s"]) * Fraction(lr)
    if mode == "f64": return ln <= q["f_s"] * lr
    return fr(ln) <= fr(q["f_s"]) * fr(lr)          # shortest-decimal reading
def in_S2(lr, mode="exact"):
    return (Fraction(lr) >= fr(q["lam_s2"])) if mode == "exact" else (lr >= q["lam_s2"] if mode == "f64" else fr(lr) >= fr(q["lam_s2"]))
def in_W2(lr, ts, mode="exact"):
    return ts or ((Fraction(lr) <= fr(q["lam_w2"])) if mode == "exact" else (lr <= q["lam_w2"] if mode == "f64" else fr(lr) <= fr(q["lam_w2"])))
reading("E64 c: the comparisons lambda_ns <= f_s x lambda_ref, lambda_ref >= lam_s2, lambda_ref <= lam_w2 are exact rationals of the binary64 data values against the rational of the decimal text of the threshold (I42); the binary64 and shortest-decimal variants are computed and reported as sensitivities")

RI = {o: ref_info(o) for o in ORB}
for o, d in RI.items():
    if d["usable"] and d["rho"] == "C":
        d["T"] = in_T(d["lam_ref"], d["lam_ns"]); d["S2"] = in_S2(d["lam_ref"]); d["W2"] = in_W2(d["lam_ref"], d["T"])
check("I51 every orbit of the panel up to the closing round is in exactly one ref cell: each has its L1 and L5 item or its ref chunk is listed as failed, and no ref item is missing without a failed chunk",
      all(d["reason"] != "missing item" for d in RI.values()) and all(sum(1 for c in chunks if c["arms"] == ["ref"] and c["method"] == "L1" and c["cell"] == (ORB[o]["k_index"], ORB[o]["class"], ORB[o]["round"])) == 1 for o in ORB), len(ORB))
sens_diff = {}
for mode in ("f64", "dec"):
    sens_diff[mode] = sorted(o for o, d in RI.items() if d["usable"] and d["rho"] == "C" and (in_T(d["lam_ref"], d["lam_ns"], mode) != d["T"] or in_S2(d["lam_ref"], mode) != d["S2"]
                                                                                              or in_W2(d["lam_ref"], in_T(d["lam_ref"], d["lam_ns"], mode), mode) != d["W2"]))
say("set-comparison readings: orbits whose set membership differs from the exact reading: binary64", sens_diff["f64"], "decimal", sens_diff["dec"])

# ------------------------------------------------------------------ round log (E62, E63)
def ref_digest(r):
    hs = sorted(c["chunk_hash"] for c in chunks if c["arms"] == ["ref"] and c["complete"] and ref_round(c) <= r)
    return sha_txt("".join(h + "\n" for h in hs).encode())
def counts_round(r):
    ids = [o for o in sorted(ORB) if ORB[o]["K"] > 0 and ORB[o]["round"] <= r]
    n_ex = sum(1 for o in ids if not RI[o]["usable"])
    rest = [o for o in ids if RI[o]["usable"]]
    n_dr = sum(1 for o in rest if RI[o]["drift"])
    ok = [o for o in rest if not RI[o]["drift"]]
    ts = sum(1 for o in ok if RI[o]["rho"] == "C" and RI[o]["T"])
    w2 = sum(1 for o in ok if RI[o]["rho"] == "C" and RI[o]["W2"])
    return dict(drift=n_dr, excluded=n_ex, orbits=len(ids), t_s=ts, target_t_s=T_T, target_w2=T_W, w2=w2)
def render(r, dec, c, prev, pdig):
    return ("round=%d\ndecision=%s\n" % (r, dec) + "".join("count.%s=%d\n" % (k, c[k]) for k in sorted(c)) +
            "panel_sha256=%s\nprereg_sha256=%s\nref_digest=%s\nprev_sha256=%s\n" % (PANEL_SHA, HEADER_HASH, pdig, prev))
logtxt = LOGB.decode("utf-8")
lines = logtxt.split("\n")
assert lines[-1] == ""
lines = lines[:-1]
entries_file = ["\n".join(lines[i:i + 13]) + "\n" for i in range(0, len(lines), 13)]
check("E62 a the log is whole entries of 13 lines", len(lines) % 13 == 0 and len(entries_file) >= 1, len(lines))
LOGR = []
prev = sha_txt(b"")
rounds_with_complete_ref = sorted({ref_round(c) for c in chunks if c["arms"] == ["ref"] and c["complete"]})
for r in range(1, len(entries_file) + 1):
    c = counts_round(r)
    closes = (c["t_s"] >= T_T and c["w2"] >= T_W) or r == RMAX
    e = render(r, "close" if closes else "continue", c, prev, ref_digest(r))
    ent = dict(round=r, counts=c, decision="close" if closes else "continue", ref_digest=ref_digest(r), text=e, prev_sha256=prev,
               equals_file=(e == entries_file[r - 1]))
    LOGR.append(ent)
    check("E62/E63 entry %d text recomputed equals the log file entry (13 lines, counts, ref digest, decision, chain)" % r, ent["equals_file"], "")
    if not ent["equals_file"]:
        fl_ = entries_file[r - 1].split("\n"); el_ = e.split("\n")
        ent["line_differences"] = [(a, b) for a, b in zip(el_, fl_) if a != b]
    prev = sha_txt(e.encode())
check("E63 h the last entry says close", LOGR[-1]["decision"] == "close", LOGR[-1]["decision"])
check("E62 rounds with complete ref chunks all have entries", rounds_with_complete_ref == list(range(1, len(entries_file) + 1)), rounds_with_complete_ref)
check("E64 c the gate set has gate_size orbits", len([o for o in ORB.values() if o["K"] == 0 and o["round"] == 1]) == q["gate_size"], q["gate_size"])

# ------------------------------------------------------------------ verdicts (E65)
def tg_ser(it):
    return dict(t=it["t"], c=it["c"])
def sitem(arm, o, spec, prec):
    """E75, E81: an item of a spec: pa1 is held by the twin arm <arm>_pa1 (horizon N_c + 1), tg, fq1, bk1 by the arm itself (horizon N_c); the horizon of the item found must be the expected one"""
    it = item(arm + "_pa1" if spec == "pa1" else arm, o, spec, prec)
    if it is not None:
        want = NC + 1 if spec == "pa1" else NC
        HORIZON_BAD.append((o, spec, it["hor"])) if it["hor"] != want else None
    return it
HORIZON_BAD = []
def verdicts64(o, prec, arm):
    """verdicts of I_7 at the ten budgets from the items of arm; LE only for Float64 from the ref L1 item. None where an item is absent."""
    out = {}
    thrv, floor = p2["ratio_threshold"], p2["le_floor"]
    if arm == "ind":
        l1 = item("ref", o, "L1")
        if l1 is not None and l1["status"] == "ok":
            ts, lam = l1["t"], l1["lam"]
            g = [n * l for n, l in zip(ts, lam)]
            out["LE"] = [series_verdict(ts, g, b, thrv, lam, floor) for b in B]
    tg = sitem(arm, o, "tg", prec)
    if tg is not None and tg["status"] == "ok" and "c" in tg:
        ts = tg["t"]
        out["SALI"] = [series_verdict(ts, [-x for x in tg["c"]["SALI.L"]], b, thrv) for b in B]
        out["FLI"] = [series_verdict(ts, tg["c"]["FLI.fli"], b, thrv) for b in B]
        out["MEGNO"] = [series_verdict(ts, tg["c"]["MEGNO.Ybar20_1"], b, thrv) for b in B]
    for ind, spec, key in (("E1", "fq1", "e1_d_threshold"), ("E2", "pa1", "e2_d_threshold"), ("WBA", "bk1", "wba_digits_threshold")):
        it = sitem(arm, o, spec, prec)
        if it is not None and it["status"] == "ok" and "rows" in it:
            w = dict(rows=it["rows"], c=it["c"], flags=it["fw"], bits=it["bits"])
            out[ind] = [window_verdict(w, spec, b, p2[key]) for b in B]
    return out
reading("E65: LE is read from the ref L1 item with the budget rule n(b) = smallest recorded step >= b of the whole series (D89 E26, as in D93), not restricted to steps <= N_c; the variant restricted to steps <= N_c is computed as a sensitivity")
def verdicts_LE_capped(o):
    l1 = item("ref", o, "L1")
    ts = [t for t in l1["t"] if t <= NC]; lam = l1["lam"][:len(ts)]
    g = [n * l for n, l in zip(ts, lam)]
    return [series_verdict(ts, g, b, p2["ratio_threshold"], lam, p2["le_floor"]) for b in B]

# ------------------------------------------------------------------ sets and drops (E64)
Kpos = [o for o in sorted(ORB) if ORB[o]["K"] > 0]
GATE = [o for o in sorted(ORB) if ORB[o]["K"] == 0 and ORB[o]["round"] == 1]
drops = {"excluded (ref item incomplete)": [], "drift": [], "rho = U": [], "no verdict data": []}
Vf = {}   # verdicts Float64 for E_test and gate
need_ind = ("tg", "fq1", "pa1", "bk1")
def ind_complete(o):
    return all(sitem("ind", o, s, "Float64") is not None and sitem("ind", o, s, "Float64")["status"] == "ok" and "excluded" not in sitem("ind", o, s, "Float64")["flags"] for s in need_ind)
reading("E64 b: drift is the stored flag drift on any item read for the orbit (ref L1 and L5, the Float64 verdict items, the Double64 items); in the toy campaign the flag occurs on all items of the one orbit that carries it, so the reading is not tested")
def any_drift(o):
    return any("drift" in it["flags"] for (arm, oo, m, pr), it in IT.items() if oo == o and arm in ("ref", "ind", "ind_pa1", "d64", "d64_pa1"))
Etest = []
for o in Kpos:
    d = RI[o]
    if not d["usable"]: drops["excluded (ref item incomplete)"].append(o); continue
    if any_drift(o): drops["drift"].append(o); continue
    if d["rho"] == "U": drops["rho = U"].append(o); continue
    if not ind_complete(o): drops["no verdict data"].append(o); continue
    Etest.append(o)
gate_nodata = [o for o in GATE if not ind_complete(o)]
for o in Etest + [g for g in GATE if g not in gate_nodata]:
    Vf[o] = verdicts64(o, "Float64", "ind")
    Vf[o]["_vcap"] = verdicts_LE_capped(o)
rhoE = {o: RI[o]["rho"] for o in Etest}
S2s = [o for o in Etest if rhoE[o] == "C" and RI[o]["S2"]]
Tss = [o for o in Etest if rhoE[o] == "C" and RI[o]["T"]]
W2s = [o for o in Etest if rhoE[o] == "C" and RI[o]["W2"]]
Cs = [o for o in Etest if rhoE[o] == "C"]
# D: in_D, in E_test, Double64 verdict data of V*
specs_V = sorted({q["spec_of"][X] for X in VS})
def d64_ok(o):
    for s in specs_V:
        if s == "l1": continue
        it = sitem("d64", o, s, "Double64")
        if it is None or it["status"] != "ok" or "excluded" in it["flags"]: return False
    return True
FAILED = sorted((a, c["cell"]) for c in chunks if c["status"] == "failed" for a in c["arms"])
FAILED_ORBITS = {(a, cell): CELLORB[cell] for a, cell in FAILED}          # E77: the orbits of a failed chunk are those of its cell (E74)
check("E77 every failed chunk is of one of the six arms and its cell has orbits in the panel up to the closing round", all(a in ARMS and cell in CELLORB for a, cell in FAILED), FAILED)
Dall = [o for o in Etest if ORB[o]["in_D"]]
Dset = [o for o in Dall if d64_ok(o)]
D_lost = [o for o in Dall if not d64_ok(o)]
DB = sorted({o for (arm, o, m, pr) in IT if arm == "chk" and m == "L5" and pr == "Double64"})
reading("E64 c (D95 clarification M8): D_B is the set of orbits carrying the L5 Double64 item of an arm chk chunk, so an orbit whose Double64 chk chunk failed twice leaves D_B (listed in the result as D_B_failed_chk) and an orbit whose BigFloat chunk failed stays in D_B with a missing BigFloat cell in its T3 row (E77, P118)")
DB_failed = sorted({o for a, cell in FAILED if a == "chk" for o in CELLORB[cell] if ORB[o]["in_D"]} - set(DB))
check("E64 c every orbit of the chk arm has the panel flag in_D", all(ORB[o]["in_D"] for o in DB), DB)
Vd = {o: verdicts64(o, "Double64", "d64") for o in Dset}
check("E75 the horizon of every verdict item read is N_c (N_c + 1 for pa1)", not HORIZON_BAD, HORIZON_BAD[:3])
# E77, E61 k: an item missing for an orbit of a candidate set is allowed only where its chunk is listed as failed (arm and cell)
def failed_cell(arm, o): return (arm, (ORB[o]["k_index"], ORB[o]["class"], ORB[o]["round"])) in FAILED_ORBITS
unlisted = []
for o in Etest:
    for sp in need_ind:
        arm = "ind_pa1" if sp == "pa1" else "ind"
        if sitem("ind", o, sp, "Float64") is None and not failed_cell(arm, o): unlisted.append(("ind", o, sp))
for o in Dall:
    for sp in specs_V:
        if sp == "l1": continue
        arm = "d64_pa1" if sp == "pa1" else "d64"
        if sitem("d64", o, sp, "Double64") is None and not failed_cell(arm, o): unlisted.append(("d64", o, sp))
check("E77/E61 k every missing item of a candidate orbit belongs to a chunk listed as failed (else the module must refuse)", not unlisted, unlisted[:5])

# strata (E64 d)
def stratum(o):
    if o in set(Tss): return 0
    if o in set(W2s): return 1
    if o in set(S2s): return 2
    return 3 if rhoE[o] == "C" else 4
strata = [[], [], [], [], []]
for o in Etest: strata[stratum(o)].append(o)
overlap = [o for o in Tss if o in set(S2s)]
sizes = dict(E_test=len(Etest), W2=len(W2s), T_s=len(Tss), S2=len(S2s), D=len(Dset), D_B=len(DB), gate=len(GATE))

# ------------------------------------------------------------------ R_a output (E52)
def out_Ra(v, b_idx=kP, V=VS, qc=QC, qr=QR, w=WD, vv=VV):
    nC = nR = nCa = 0
    for X in V:
        a = v[X]
        c = a[b_idx]
        stable = c != UN and all(a[b_idx - j] == c for j in range(1, w + 1)) if b_idx - w >= 0 else None
        if stable is None: raise ValueError("budget outside B")
        if c == CH:
            nCa += 1
            if stable: nC += 1
        elif c == RG and stable: nR += 1
    if nC >= qc: return 1
    if nCa <= vv and nR >= qr: return 2
    return 0
OUTN = {0: "unresolved", 1: "chaotic-supported", 2: "regular-supported"}
OUT = {o: out_Ra(Vf[o]) for o in Etest}
GOUT = {o: out_Ra(Vf[o]) for o in GATE if o in Vf}
def contra(o, outv): return (outv == 1 and rhoE[o] == "R") or (outv == 2 and rhoE[o] == "C")

# ------------------------------------------------------------------ bootstrap (E66)
def randint_list(seed, m, n):
    g = Xoshiro.from_seed(seed)
    return [g.randint(m) for _ in range(n)]
def pct7(vals, p):
    n = len(vals); h = (n - 1) * p + 1; j = math.floor(h)
    if h == j: return vals[j - 1]
    return vals[j - 1] + (h - j) * (vals[j] - vals[j - 1])
def bound(values, kind, Bn):
    """values: list of Fraction or None per replicate; kind 'upper' (p = 1 - alpha) or 'lower' (p = alpha)"""
    valid = sorted(v for v in values if v is not None)
    dropped = len(values) - len(valid)
    if dropped * 20 > len(values) or not valid:
        return None, dropped
    p = 1 - thr["alpha_ci"] if kind == "upper" else thr["alpha_ci"]
    return pct7(valid, p), dropped
BR = int(q["bootstrap_resamples"])

# HA
def ha_stats(mult):
    """mult: dict orbit -> multiplicity"""
    dec = sum(m for o, m in mult.items() if OUT[o] != 0)
    con = sum(m for o, m in mult.items() if contra(o, OUT[o]))
    tot = sum(mult.values())
    un = sum(m for o, m in mult.items() if OUT[o] == 0)
    def a(S):
        n = sum(m for o, m in mult.items() if o in S)
        return None if n == 0 else Fraction(sum(m for o, m in mult.items() if o in S and OUT[o] == 0), n)
    W, T, S2_ = set(W2s), set(Tss), set(S2s)
    aW, aT, aS = a(W), a(T), a(S2_)
    return dict(e=(Fraction(con, dec) if dec else None), a_E=Fraction(un, tot) if tot else None, a_W2=aW, a_Ts=aT, a_S2=aS,
                d_W=(aW - aS if aW is not None and aS is not None else None), d_T=(aT - aS if aT is not None and aS is not None else None))
def ha_point():
    return ha_stats({o: 1 for o in Etest})
HA_PT = ha_point()
reps = []
for r in range(1, BR + 1):
    seed = stream_seed(TSEED, "boot-ha", str(r))
    g = Xoshiro.from_seed(seed)
    mult = {o: 0 for o in Etest}
    for st in strata:
        for _ in range(len(st)):
            mult[st[g.randint(len(st)) - 1]] += 1
    reps.append(ha_stats(mult))
HA_B = {}
for key, kind in (("e", "upper"), ("a_E", "upper"), ("d_W", "lower"), ("d_T", "lower")):
    HA_B[key] = bound([x[key] for x in reps], kind, BR)
reading("E66 b: strata drawn in the order T_s, W2 minus T_s, S2, other rho = C, rho = R, each in increasing orbit id; the integer draw is x mod m + 1 for the first generator output x not above 2^64 - 1 - (2^64 mod m) (the rule of the D93 oracle; the amendment names 'unbiased integer draws')")
reading("E66 d: the dropped-replicate limit (more than 5 percent) is applied per statistic: dropped x 20 > number of replicates")
ha_notdec = []
for nm, S in (("E_test", Etest), ("W2", W2s), ("T_s", Tss), ("S2", S2s)):
    if len(S) < q["min_set"]: ha_notdec.append("%s has %d < min_set" % (nm, len(S)))
if HA_PT["e"] is None: ha_notdec.append("e undefined")
for key, (bv, dr) in HA_B.items():
    if bv is None: ha_notdec.append("bound of %s undefined (dropped %d of %d)" % (key, dr, BR))
if gate_nodata: ha_notdec.append("gate orbit without verdict: %s" % gate_nodata)
g_cnt = sum(1 for o in GATE if o in GOUT and GOUT[o] == 1)
HA = dict(point={k: v for k, v in HA_PT.items()}, g=g_cnt, bounds={k: (v[0], v[1]) for k, v in HA_B.items()}, not_decidable=ha_notdec)
if ha_notdec:
    HA["decision"] = "not decidable"; HA["failing"] = []
else:
    fail = []
    if not (HA_B["e"][0] <= thr["e_max"]): fail.append(1)
    if not (HA_B["a_E"][0] <= thr["a_max"]): fail.append(2)
    if not (HA_PT["d_W"] >= thr["delta_a"] and HA_PT["d_T"] >= thr["delta_a"] and HA_B["d_W"][0] > 0 and HA_B["d_T"][0] > 0): fail.append(3)
    if not (g_cnt <= q["g_max"]): fail.append(4)
    HA["failing"] = fail; HA["decision"] = "SUPPORTED" if not fail else "not supported"

# ------------------------------------------------------------------ HB (E68, E69)
I6 = list(q["voter_bits"])
def hb_eval(qq, rounddown=False, noLE=False, caps_le_steps=False):
    cost = Costs(qq)
    CR, matched, le_read, le_unm = cost.matched(VS)
    res = dict(C_R=CR, le_read=le_read, le_cost_unmatched=bool(le_unm), sets={}, defects=[], undecidable_reasons=[])
    allpass, anycap = True, False
    for A, S, comps in (("T_s", Tss, I6), ("W2", W2s, ([] if noLE else ["LE"]) + I6)):
        rows = []
        outv = [OUT[o] for o in S]
        nA = len(S)
        dec = sum(1 for x in outv if x != 0)
        for X in comps:
            m_ = matched[X]
            if X == "LE":
                b = (m_["down"] if rounddown else le_read)
                capped = False
                if b is None: res["undecidable_reasons"].append("LE has no round-down budget"); rows.append(dict(X=X, b=None)); continue
            else:
                b = m_["down"] if rounddown else m_["b"]
                capped = (b is None)
                if capped:
                    if rounddown: res["undecidable_reasons"].append("%s has no round-down budget" % X)
                    else: res["defects"].append(X); res["undecidable_reasons"].append("capped non-LE comparator %s" % X)
                    rows.append(dict(X=X, b=None, capped=True)); continue
            k = B.index(b)
            n11 = n10 = n01 = n00 = 0
            for o, x in zip(S, outv):
                vx = Vf[o]["_vcap"][k] if (X == "LE" and caps_le_steps) else Vf[o][X][k]
                rr = (x == 2); xr = (vx == RG)
                if rr and xr: n11 += 1
                elif rr: n10 += 1
                elif xr: n01 += 1
                else: n00 += 1
            mm = n10 + n01
            p = mcnemar_tail(mm, n01) if mm > 0 else Fraction(1)
            nrule, nx = n11 + n10, n11 + n01
            margin = nx > 0 and nrule * thr["rho_hb"].denominator <= thr["rho_hb"].numerator * nx
            passed = (p <= thr["alpha_hb"]) and margin
            rows.append(dict(X=X, b=b, cost=b * cost.cx(X), capped=False, cost_unmatched=(bool(le_unm) if X == "LE" and not rounddown else False),
                             FR_X=Fraction(nx, nA) if nA else Fraction(0), FR_rule=Fraction(nrule, nA) if nA else Fraction(0), n11=n11, n10=n10, n01=n01, n00=n00, p=p, margin=bool(margin), passed=bool(passed)))
            if not passed: allpass = False
        res["sets"][A] = dict(rows=rows, n=nA, decided=dec, cov=Fraction(dec, nA) if nA else None, coverage_price=(1 - Fraction(dec, nA)) if nA else None)
    covok = all(res["sets"][A]["cov"] is not None and res["sets"][A]["cov"] >= thr["c_min"] for A in ("T_s", "W2"))
    res["all_pass"] = bool(allpass); res["coverage_ok"] = bool(covok)
    ncomp = sum(len(res["sets"][A]["rows"]) for A in res["sets"])
    res["n_comparisons"] = ncomp
    if len(Tss) < q["n_t_min"] or len(W2s) < q["n_w_min"]: res["undecidable_reasons"].append("set below floor")
    res["decision"] = "not decidable" if res["undecidable_reasons"] else ("SUPPORTED" if (allpass and covok) else "not supported")
    return res
HB = hb_eval(q)
HB["underpowered"] = bool(len(Tss) < T_T or len(W2s) < T_W)
HB["cost_unmatched_LE_flag"] = True
HB_noLE = hb_eval(q, noLE=True)
HB_capsteps = hb_eval(q, caps_le_steps=True)
q_pilot = dict(q); q_pilot.update({"cost_" + s: qt["pilot_cost_" + s] for s in ("l1", "tg", "pa1", "bk1", "fq1")})
meas = {}
for _arm in ("ind", "ind_pa1"):       # E78, P119: the Float64 arms only; the arms d64 and d64_pa1 are not read for costs
    for s, v in MAN["arms"][_arm]["specs"].items():
        if "measured_ns_per_orbit_step" in v: meas[s] = v["measured_ns_per_orbit_step"]
_l1 = [ch for ch in MAN["chunks"] if ch["arms"] == ["ref"] and "method=L1 " in ch["key"] and ch["status"] == "ok"]
_sumsteps = sum(ch["steps"] for ch in _l1)
if _sumsteps > 0: meas["l1"] = 1e9 * math.fsum(ch["thread_seconds"] for ch in _l1) / _sumsteps
q_meas = dict(q); q_meas.update({"cost_" + s: meas[s] for s in meas})
missing_meas = [s for s in ("l1", "tg", "pa1", "bk1", "fq1") if s not in meas]
SENS = dict(round_down=hb_eval(q, rounddown=True), pilot=hb_eval(q_pilot), without_LE=HB_noLE)
if missing_meas:
    SENS["measured"] = dict(decision="not decidable", undecidable_reasons=["manifest lacks a measured cost for spec(s): " + ",".join(missing_meas)])
else:
    SENS["measured"] = hb_eval(q_meas)
reading("E69 b (i): LE under the round-down sensitivity is read at the round-down budget of the spec L1 (no round-down budget makes the sensitivity not decidable)")
reading("E78: the measured cost of tg, fq1, bk1 is arms.ind.specs.<spec>.measured_ns_per_orbit_step, of pa1 arms.ind_pa1.specs.pa1.measured_ns_per_orbit_step; that of L1 is 1e9 x (sum of thread_seconds) / (sum of steps) over the manifest chunk entries of the arm ref whose key names method=L1 and whose status is ok")
reading("E68: the coverage condition and the floors do not depend on the budgets; FR_X(A) = 0 or a margin failure is a failed comparison, a failed comparison with m = 0 has p = 1")

# power at achieved sizes (E69 c)
def pp_model(rec):
    m = rec["pp_hb"]["model"]
    out = {}
    for S in ("T", "RW"):
        pr = m["rule"][S]["p"]; d = m["rule"][S]["d"]
        dprob = min(max((d - pr) / (1 - pr), 0.0), 1.0) if pr < 1 else 1.0
        sx = [m["comp"][X][S]["s"] for X in I7]; px = [m["comp"][X][S]["p"] for X in I7]
        qx = [(max(a - b * pr, 0.0) / (1 - pr) if pr < 1 else 0.0) for a, b in zip(px, sx)]
        out[S] = dict(p=pr, dprob=dprob, s=sx, q=qx)
    return out
def sim_hb_unit(model, nT, nRW, runs, crit_arr, purpose):
    seeds = np.array([stream_seed(TSEED, purpose, str(i)) for i in range(1, runs + 1)], dtype=np.uint64)
    xo = XoV(seeds)
    cnt = {S: dict(n11=np.zeros((7, runs), dtype=np.int32), n10=np.zeros((7, runs), dtype=np.int32), n01=np.zeros((7, runs), dtype=np.int32), dec=np.zeros(runs, dtype=np.int32)) for S in ("T", "RW")}
    for S, n in (("T", nT), ("RW", nRW)):
        mo = model[S]; sx = np.array(mo["s"])[:, None]; qx = np.array(mo["q"])[:, None]
        for a in range(n):
            U = np.empty((9, runs))
            for b in range(9): U[b] = xo.uniform()
            rfr = U[0] < mo["p"]
            xfr = np.where(rfr[None, :], U[1:8] < sx, U[1:8] < qx)
            decd = rfr | (U[8] < mo["dprob"])
            R = rfr[None, :]
            cnt[S]["n11"] += (R & xfr); cnt[S]["n10"] += (R & ~xfr); cnt[S]["n01"] += (~R & xfr); cnt[S]["dec"] += decd
    rn, rd = thr["rho_hb"].numerator, thr["rho_hb"].denominator
    cm, cd = thr["c_min"].numerator, thr["c_min"].denominator
    def comp_ok(n11, n10, n01):
        mm = n10 + n01; a = n10 + n11; c = n01 + n11
        return (mm >= 1) & (n01 >= crit_arr[np.minimum(mm, len(crit_arr) - 1)]) & (c >= 1) & (a * rd <= rn * c)
    ok = np.ones(runs, dtype=bool)
    for xi in range(7):
        if xi >= 1: ok &= comp_ok(cnt["T"]["n11"][xi], cnt["T"]["n10"][xi], cnt["T"]["n01"][xi])
        ok &= comp_ok(cnt["T"]["n11"][xi] + cnt["RW"]["n11"][xi], cnt["T"]["n10"][xi] + cnt["RW"]["n10"][xi], cnt["T"]["n01"][xi] + cnt["RW"]["n01"][xi])
    ok &= cnt["T"]["dec"] * cd >= cm * nT
    ok &= (cnt["T"]["dec"] + cnt["RW"]["dec"]) * cd >= cm * (nT + nRW)
    return int(ok.sum())
POW_HB = None
if HB["underpowered"]:
    nT_, nRW_ = len(Tss), len(W2s) - len(Tss)
    crit = np.array(mcnemar_critical_table(len(W2s) + 2, thr["alpha_hb"]), dtype=np.int64)
    runs = q["pp_runs"]
    succ = sim_hb_unit(pp_model(REC), nT_, nRW_, runs, crit, "pp-achieved-hb")
    pw = succ / runs
    POW_HB = dict(nT=nT_, nRW=nRW_, runs=runs, successes=succ, power=pw, se=math.sqrt(pw * (1 - pw) / runs), planned=qt["power_hb_planned"])
    reading("E69 c: rates read from record.pp_hb.model (rule p and d, comparator p and s per stratum, as written, floats of the 12-digit text); the derived probabilities use the D93 oracle's binary64 formulas; unit of the stream = run index 1..pp_runs; orbits of T_s first, then of R_W, nine uniforms per orbit")

# ------------------------------------------------------------------ HC (E70)
def hc_rows(Vs):
    rows = []
    for o in Dset:
        nC = sum(1 for X in VS if Vf[o][X][kP] == CH); nR = sum(1 for X in VS if Vf[o][X][kP] == RG)
        label = "c" if nC > nR else ("r" if nR > nC else "tie")
        f = (label == "tie") or ((label == "c") != (rhoE[o] == "C"))
        sI = Fraction(0) if nC + nR == 0 else 1 - Fraction(abs(nC - nR), nC + nR)
        diff = sum(1 for X in VS if Vf[o][X][kP] != UN and Vd[o][X][kP] != UN and Vf[o][X][kP] != Vd[o][X][kP])
        rows.append(dict(o=o, nC=nC, nR=nR, f=bool(f), sI=sI, sP=Fraction(diff, len(VS))))
    return rows
HCR = hc_rows(None)
reading("E70/P104: s_P denominator |V*|; D orbits lacking Double64 verdicts of V* leave the set (counted); the bootstrap is not stratified by cell")
def auc_of(rows_w, key):
    """rows_w: list of (row, multiplicity); exact Mann-Whitney with ties one half; None if a class is empty"""
    F = [(r[key], m) for r, m in rows_w if r["f"] and m]
    N = [(r[key], m) for r, m in rows_w if (not r["f"]) and m]
    wf, wn = sum(m for _, m in F), sum(m for _, m in N)
    if wf == 0 or wn == 0: return None
    twice = 0
    for sf, mf in F:
        for sn, mn in N:
            twice += mf * mn * (2 if sf > sn else (1 if sf == sn else 0))
    return Fraction(twice, 2 * wf * wn)
pt_I, pt_P = auc_of([(r, 1) for r in HCR], "sI"), auc_of([(r, 1) for r in HCR], "sP")
nfail, nnon = sum(1 for r in HCR if r["f"]), sum(1 for r in HCR if not r["f"])
HC = dict(n=len(HCR), failures=nfail, nonfailures=nnon, AUC_I=pt_I, AUC_P=pt_P, lost_D_orbits=D_lost)
HC["Delta"] = (pt_I - pt_P) if (pt_I is not None and pt_P is not None) else None
dv = []
for r in range(1, BR + 1):
    g = Xoshiro.from_seed(stream_seed(TSEED, "boot-hc", str(r)))
    mult = [0] * len(HCR)
    for _ in range(len(HCR)): mult[g.randint(len(HCR)) - 1] += 1
    wr = list(zip(HCR, mult))
    a, b = auc_of(wr, "sI"), auc_of(wr, "sP")
    dv.append(None if (a is None or b is None) else a - b)
HC["lower"] = bound(dv, "lower", BR); HC["upper"] = bound(dv, "upper", BR)
nd = []
if nfail < q["min_set"]: nd.append("failures %d < min_set" % nfail)
if nnon < q["min_set"]: nd.append("non-failures %d < min_set" % nnon)
if HC["Delta"] is None: nd.append("AUC undefined")
if HC["lower"][0] is None or HC["upper"][0] is None: nd.append("bound undefined")
HC["not_decidable"] = nd
if nd: HC["decision"] = "not decidable"
elif HC["Delta"] >= thr["delta_c"] and HC["lower"][0] > 0: HC["decision"] = "SUPPORTED"
elif HC["Delta"] <= -thr["delta_c"] and HC["upper"][0] < 0: HC["decision"] = "REVERSED"
else: HC["decision"] = "not supported"
HC["underpowered"] = bool(nfail < T_F)
reading("E70: a bound needed by the decision (either one) that is undefined makes HC not decidable, even if only the other would be used")
POW_HC = None
if HC["underpowered"] and nfail >= 1 and nnon >= 1:
    pc = REC["pp_hc"]
    dI = math.sqrt(2) * p3o.NormalDist().inv_cdf(pc["auc_i_used"]); dP = math.sqrt(2) * p3o.NormalDist().inv_cdf(pc["auc_p"])
    zz = q["z_one_sided"]
    def sim_hc_unit(nf, nn, runs):
        seeds = np.array([stream_seed(TSEED, "pp-achieved-hc", str(i)) for i in range(1, runs + 1)], dtype=np.uint64)
        xo = XoV(seeds); R = len(seeds)
        nsteps = 2 * (nf + nn)
        U = np.empty((nsteps, R))
        for s in range(nsteps): U[s] = xo.uniform()
        Z = ppnd16(U.reshape(-1)).reshape(nsteps, R)
        sI = np.concatenate([Z[:nf] + dI, Z[nf:nf + nn]], axis=0).T
        sP = np.concatenate([Z[nf + nn:2 * nf + nn] + dP, Z[2 * nf + nn:]], axis=0).T
        def ucount(S):
            ranks = np.argsort(np.argsort(S, axis=1, kind="stable"), axis=1, kind="stable") + 1
            return ranks[:, :nf].sum(1) - nf * (nf + 1) // 2
        UI, UP = ucount(sI), ucount(sP)
        AI, AP = UI / (nf * nn), UP / (nf * nn)
        dc = thr["delta_c"]
        pt = ((UI.astype(object) - UP.astype(object)) * dc.denominator >= dc.numerator * nf * nn).astype(bool)
        def var(A):
            Q1 = A / (2 - A); Q2 = 2 * A * A / (1 + A)
            return (A * (1 - A) + (nf - 1) * (Q1 - A * A) + (nn - 1) * (Q2 - A * A)) / (nf * nn)
        D = AI - AP; SE = np.sqrt(var(AI) + var(AP))
        return int((pt & (D - zz * SE > 0)).sum())
    runs = q["pp_runs"]
    succ = sim_hc_unit(nfail, nnon, runs); pw = succ / runs
    POW_HC = dict(failures=nfail, nonfailures=nnon, runs=runs, successes=succ, power=pw, se=math.sqrt(pw * (1 - pw) / runs), planned=qt["power_hc_planned"])
    reading("E70 power at the achieved sizes: the D93 E57 simulation with AUC_I_used = record.pp_hc.auc_i_used and AUC_P = record.pp_hc.auc_p, n_f = achieved failures, n_n = achieved non-failures, stream unit = run index")

# ------------------------------------------------------------------ T3, T4, T9, T11, F2, F3
def labelof2(P, itm, N, o, hfix=0):
    if itm is None or str(N // 10) not in itm["lamfull"] or str(N) not in itm["lamfull"]: return (UNRES, UNRES)
    lam, lam10 = itm["lamfull"][str(N)], itm["lamfull"][str(N // 10)]
    inst = UNRES
    if math.isfinite(lam) and math.isfinite(lam10) and lam10 > 0:
        r = 10 * lam / lam10
        if lam >= P["lambda_threshold"] and r >= P["growth_ratio"]: inst = POS
        elif lam < P["lambda_threshold"] and r < P["growth_ratio"]: inst = NONE
    stc = UNRES
    if inst == POS: stc = "irregular"
    elif inst == NONE:
        x0 = float(itm["ic"].split(",")[0])
        if str(hfix) in itm["xl"] and str(N) in itm["xl"]:
            rho1 = (itm["xl"][str(hfix)] - x0) / (2 * math.pi * hfix); rho2 = (itm["xl"][str(N)] - itm["xl"][str(hfix)]) / (2 * math.pi * (N - hfix))
            if abs(rho1 - rho2) <= P["rotation_tol"]: stc = "quasiperiodic"
    return (inst, stc)
def t3(NB):
    """E79: the Double64 and BigFloat L5 values of an orbit of D_B from the items of the arm chk found by content (method L5, precision label without bits, the horizon of the item itself, which must be N_B);
    the Float64 value from the ref L5 item at the recorded step N_B. A value of a failed chunk is absent (None, P118). No label comparison: the chk items have no steps at the scale of N_1 (owner decision 2)."""
    rows = []
    for o in sorted(DB):
        f = item("ref", o, "L5")
        vals = {}
        for prec in ("Double64", "BigFloat"):
            its = [it for (arm, oo, m, pr), it in IT.items() if arm == "chk" and oo == o and m == "L5" and pr == prec]
            assert len(its) <= 1
            vals[prec] = None if not its else (its[0]["lamfull"][str(its[0]["hor"])] if its[0]["hor"] == NB and str(its[0]["hor"]) in its[0]["lamfull"] else "horizon-mismatch")
        lf = f["lamfull"][str(NB)] if f is not None and str(NB) in f["lamfull"] else None
        miss = [k for k, v in (("F64", lf), ("D64", vals["Double64"]), ("BF", vals["BigFloat"])) if v is None]
        if lf is None or None in vals.values() or "horizon-mismatch" in vals.values():
            rows.append(dict(o=o, available=False, missing=miss, lam_F64=lf, lam_D64=vals["Double64"], lam_BF=vals["BigFloat"], ladder_unresolved=None, rho=rhoE.get(o, "dropped"))); continue
        ld, lb = vals["Double64"], vals["BigFloat"]
        lad = abs(lf - ld) > p1["ladder_tol"] or abs(ld - lb) > p1["ladder_tol"]
        rows.append(dict(o=o, available=True, lam_F64=lf, lam_D64=ld, lam_BF=lb, ladder_unresolved=bool(lad), rho=rhoE.get(o, "dropped")))
    return rows
chkh = sorted({c["horizon"] for c in chunks if c["arms"] == ["chk"]})
check("E79 the chk horizon is a single value and equals the key N_B of the P1 parameter table", chkh == [p1["N_B"]], (chkh, p1["N_B"]))
T3 = {"N_B=%d" % p1["N_B"]: t3(p1["N_B"])}
reading("E79: the ladder test is formed from the tolerance rule only (|F64 - D64| or |D64 - BF| above ladder_tol); the label part of the original ladder is dropped with the column (owner decision 2); an orbit with a value absent (failed chunk) has no ladder entry")

def counts_by(S, keyf):
    out = {}
    for o in S:
        d = out.setdefault(keyf(o), {0: 0, 1: 0, 2: 0}); d[OUT[o]] += 1
    return {str(k): {OUTN[a]: n for a, n in v.items()} for k, v in sorted(out.items())}
T4 = {nm: dict(pooled=counts_by(S, lambda o: "all"), by_K=counts_by(S, lambda o: ORB[o]["K"]), by_round=counts_by(S, lambda o: ORB[o]["round"]))
      for nm, S in (("E_test", Etest), ("W2", W2s), ("T_s", Tss), ("S2", S2s))}
T2 = {}
for r in range(1, RC + 1):
    T2[str(r)] = dict(orbits=sum(1 for o in ORB.values() if o["round"] == r), K0=sum(1 for o in ORB.values() if o["round"] == r and o["K"] == 0),
                      E_test=sum(1 for o in Etest if ORB[o]["round"] == r), W2=sum(1 for o in W2s if ORB[o]["round"] == r), T_s=sum(1 for o in Tss if ORB[o]["round"] == r), S2=sum(1 for o in S2s if ORB[o]["round"] == r))
T2k = {}
for o in Etest:
    k = "%d|%d" % (ORB[o]["k_index"], ORB[o]["round"])
    d = T2k.setdefault(k, dict(E_test=0, W2=0, T_s=0, S2=0)); d["E_test"] += 1; d["W2"] += o in W2s; d["T_s"] += o in Tss; d["S2"] += o in S2s
last_counts = LOGR[-1]["counts"]
log_difference = dict(t_s=last_counts["t_s"] - len(Tss), w2=last_counts["w2"] - len(W2s))
T9 = None
pairs = {}
over = 0; errA = [0, 0]; errB = [0, 0]
for o in Dset:
    ra = OUT[o]
    ov = any(Vf[o][X][kP] != UN and Vd[o][X][kP] != UN and Vf[o][X][kP] != Vd[o][X][kP] for X in VS)
    rb = 0 if ov else ra
    over += int(ov and ra != 0)
    k = "%s|%s" % (OUTN[ra], OUTN[rb]); pairs[k] = pairs.get(k, 0) + 1
    for rr, e in ((ra, errA), (rb, errB)):
        if rr != 0: e[1] += 1; e[0] += int(contra(o, rr))
T9 = dict(pairs=pairs, overrides_of_decided=over, decided_Ra=errA[1], contradicting_Ra=errA[0], decided_Rab=errB[1], contradicting_Rab=errB[0])
# T11
fl = {}
reading("T11: the stored flags are counted over the items of the chunks that are not failed (an item of a chunk failed twice is not an item read)")
for (arm, o, m, pr), it in IT.items():
    if CHK[it["fn"]]["status"] == "failed": continue
    for f in it["flags"]:
        k = "%s|%s|%s" % (arm, m, f); fl[k] = fl.get(k, 0) + 1
failed_chunks = [dict(fn=c["fn"], key=c["chunk_key"], arms=c["arms"], orbits=[o for o in sorted(ORB) if False]) for c in chunks if c["status"] == "failed"]
# F2
Cost = Costs(q); CRv, matched, le_read, le_unm = Cost.matched(VS)
def rc_point(outs):
    dec = sum(1 for o in Etest if outs[o] != 0); con = sum(1 for o in Etest if contra(o, outs[o]))
    return dict(decided=dec, contradicting=con, decided_share=Fraction(dec, len(Etest)), e=(Fraction(con, dec) if dec else None))
F2 = []
sz = len(VS)
for dc in (-1, 0, 1):
    for dr in (-1, 0, 1):
        qc, qr = QC + dc, QR + dr
        if 1 <= qc <= sz and sz // 2 + 1 <= qr <= sz:
            outs = {o: out_Ra(Vf[o], qc=qc, qr=qr) for o in Etest}
            F2.append(dict(curve="rule family", point="qC=%d,qR=%d" % (qc, qr), cross=bool(dc == 0 or dr == 0), **rc_point(outs)))
for X in I7:
    if X != "LE" and matched[X]["capped"]: continue
    b = le_read if X == "LE" else matched[X]["b"]
    k = B.index(b)
    outs = {o: (1 if Vf[o][X][k] == CH else (2 if Vf[o][X][k] == RG else 0)) for o in Etest}
    F2.append(dict(curve="single " + X, point="b=%d" % b, **rc_point(outs)))
reading("E72 F2: the rule family is the product grid of q_C and q_R each varied by one step inside their ranges (up to 9 points, theta* included), as in the D93 oracle")
# F3
F3 = []
for key in ("sI", "sP"):
    Fv = [r[key] for r in HCR if r["f"]]; Nv = [r[key] for r in HCR if not r["f"]]
    pts = [(Fraction(0), Fraction(0))] if False else []
    for s in sorted(set(Fv + Nv), reverse=True):
        tpr = Fraction(sum(1 for x in Fv if x >= s), len(Fv)) if Fv else None
        fpr = Fraction(sum(1 for x in Nv if x >= s), len(Nv)) if Nv else None
        F3.append(dict(curve=key, score=s, fpr=fpr, tpr=tpr))
    F3.append(dict(curve=key, auc=auc_of([(r, 1) for r in HCR], key)))
reading("E72 F3: ROC rows for each distinct score value in decreasing order: the pair (FPR, TPR) of the rule 'score >= value'; one AUC row per curve; no (0,0) start row")

# ------------------------------------------------------------------ result
def jf(o):
    if isinstance(o, Fraction): return rtext(o)
    if isinstance(o, dict): return {str(k): jf(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [jf(x) for x in o]
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.floating): return float(o)
    return o
res = dict(
    checks=checks, readings=READINGS,
    fingerprints=dict(prereg_header=HEADER_HASH, panel=PANEL_SHA, log=sha_txt(LOGB), record=sha_file(CP("inputs", "p3-abstention-calibration-record.json")),
                      p3_table=sha_file(CP("inputs", "p3-abstention-calibration-parameters.toml")), test_table=sha_file(CP("inputs", "p3-abstention-test-parameters.toml")),
                      manifest=sha_file(MANF), chunk_files={c["fn"]: c["sha256"] for c in chunks}),
    log=[dict(round=e["round"], counts=e["counts"], decision=e["decision"], ref_digest=e["ref_digest"], prev_sha256=e["prev_sha256"], text=e["text"], equals_file=e["equals_file"],
              line_differences=e.get("line_differences")) for e in LOGR],
    ref_info={str(o): dict(usable=d["usable"], reason=d["reason"], rho=d["rho"], lam_ref=d["lam_ref"], lam_ns=d["lam_ns"], drift=d["drift"], T=d.get("T"), S2=d.get("S2"), W2=d.get("W2"), axes=d.get("axes")) for o, d in RI.items()},
    set_sensitivity=sens_diff,
    sets=dict(D_B_failed_chk=DB_failed, E_test=Etest, W2=W2s, T_s=Tss, S2=S2s, D=Dset, D_B=DB, gate=GATE, rho_C=Cs, strata=strata, overlap_T_S2=overlap, sizes=sizes, D_lost_double64=D_lost, gate_without_verdicts=gate_nodata),
    drops=drops,
    rho_counts={r: sum(1 for o in Kpos if RI[o]["usable"] and RI[o]["rho"] == r) for r in "CRU"},
    verdict_codes=dict(E_test={str(o): {X: "".join("cru"[0 if c == CH else (1 if c == RG else 2)] for c in Vf[o][X]) for X in I7} for o in Etest},
                       gate={str(o): {X: "".join("cru"[0 if c == CH else (1 if c == RG else 2)] for c in Vf[o][X]) for X in I7} for o in GATE if o in Vf}),
    outputs=dict(R_a={str(o): OUTN[OUT[o]] for o in Etest}, gate={str(o): OUTN[GOUT[o]] for o in GOUT},
                 R_ab_D={str(o): ("unresolved" if any(Vf[o][X][kP] != UN and Vd[o][X][kP] != UN and Vf[o][X][kP] != Vd[o][X][kP] for X in VS) else OUTN[OUT[o]]) for o in Dset}),
    HA=HA, HB=HB, HB_without_LE=HB_noLE, HB_le_steps_capped_at_Nc=dict(decision=HB_capsteps["decision"], all_pass=HB_capsteps["all_pass"]), sensitivities=SENS, power_hb_achieved=POW_HB,
    HC=dict(HC, rows=HCR), power_hc_achieved=POW_HC,
    T2=T2, T2_by_k_round=T2k, log_difference=log_difference, T3=T3, T4=T4, T9=T9, T11=dict(flag_counts=fl, failed_chunks=[c["fn"] for c in chunks if c["status"] == "failed"], failed_cells=[dict(fn=c["fn"], arms=c["arms"], cell="k%d:%s:r%d" % c["cell"], orbits=CELLORB[c["cell"]]) for c in chunks if c["status"] == "failed"],
                                               manifest_status=MAN["status"]),
    T10=dict(measured=meas, cost_keys={k: q["cost_" + k] for k in ("l1", "tg", "pa1", "bk1", "fq1")}, pilot={k: qt["pilot_cost_" + k] for k in ("l1", "tg", "pa1", "bk1", "fq1")}),
    F2=F2, F3=F3,
)
os.makedirs(os.path.dirname(os.path.abspath(RESULT)), exist_ok=True)
json.dump(jf(res), open(RESULT, "w", encoding="utf-8"), indent=1)
say("written", RESULT)
say("HA", HA["decision"], HA["failing"], HA["not_decidable"], "HB", HB["decision"], HB["undecidable_reasons"], "HC", HC["decision"], HC["not_decidable"])
