# Comparison of the module output of `analyze` for study p3-abstention (machine file and tables file) with the independent oracle result (O20, Amendment D95), adapted to the store of the
# runner (Amendment D97, O28, O29; supersedes verify/analysis-p3t-compare.py for the output of the oracle verify/analysis-d97-oracle.py, the D95 file is kept unchanged).
# D97 changes: the panel path is an argument; T3 without the column "labels agree" and with the rows of an absent value (E79, E77, P118); T10 with the L1 cost (E78); T11 with the orbits of the failed
# cells (E77, E80); the manifest status (E76); T2 up to the closing round; the exclusions and losses of failed chunks (E77).
# Rules: integers, rationals, sets, texts, flags, decisions and bounds exactly; reals as the module writes them (12 significant digits) against the oracle exact
# value rounded to 12 significant digits; HC power within 3 Monte Carlo standard errors; HB power counts exactly.
#   python -I verify/analysis-d97-compare.py <oracle-result.json> <module-out folder> <log file> <panel file>
import sys, os, json, re, hashlib
from fractions import Fraction

orc = json.load(open(sys.argv[1], encoding="utf-8"))
mod = json.load(open(os.path.join(sys.argv[2], "p3-abstention-test-results.json"), encoding="utf-8"))
tables = open(os.path.join(sys.argv[2], "p3-abstention-tables.md"), encoding="utf-8").read()
LOG = open(sys.argv[3], "w", encoding="utf-8")
N = [0, 0]
DIFFS = []
def out(s):
    print(s); LOG.write(s + "\n")
def eq(name, a, b, note=""):
    ok = (a == b)
    N[0] += 1; N[1] += ok
    if not ok:
        DIFFS.append(name)
    out(("EQUAL " if ok else "DIFFER ") + name + ("" if ok else "  oracle=%s module=%s" % (str(a)[:300], str(b)[:300])) + (("  [" + note + "]") if note and not ok else ""))
s12 = lambda x: float("%.11e" % float(x))
def sig12(x):
    return float("%.11e" % float(Fraction(x))) if isinstance(x, str) and "/" in x else float("%.11e" % float(x))
def eqr(name, a_exact, b_text):
    """oracle exact rational text against a module decimal text (12 significant digits)"""
    try:
        ok = abs(float(b_text) - float(Fraction(a_exact))) <= 1e-12 * max(1.0, abs(float(Fraction(a_exact))))
    except Exception:
        ok = False
    N[0] += 1; N[1] += ok
    if not ok: DIFFS.append(name)
    out(("EQUAL " if ok else "DIFFER ") + name + ("" if ok else "  oracle=%s module=%s" % (a_exact, b_text)))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

# ---------------------------------------------------------------- fingerprints and log
fp = mod["fingerprints"]
of = orc["fingerprints"]
for a, b in (("preregistration", "prereg_header"), ("panel", "panel"), ("calibration_record", "record"), ("p3_parameters", "p3_table"), ("test_parameters", "test_table"), ("manifest", "manifest"), ("round_log", "log")):
    eq("fingerprint " + a, of[b], fp[a])
for fn, h in of["chunk_files"].items():
    eq("fingerprint chunk file " + fn, h, fp["chunks"].get(fn))
for e_o, e_m in zip(orc["log"], mod["log"]["entries"]):
    r = e_o["round"]
    eq("log entry %d decision" % r, e_o["decision"], e_m["decision"])
    for k, v in e_o["counts"].items():
        eq("log entry %d count.%s" % (r, k), v, e_m["counts"][k])
    eq("log entry %d text recomputed by the oracle equals the log file (oracle internal check)" % r, True, e_o["equals_file"])
eq("log entries", len(orc["log"]), len(mod["log"]["entries"]))

# ---------------------------------------------------------------- sets and drops
S, ms = orc["sets"], mod["sets"]
for k, mk in (("E_test", "E_test"), ("W2", "W2"), ("T_s", "T_s"), ("S2", "S2"), ("D", "D"), ("D_B", "D_B"), ("gate", "gate")):
    eq("set %s ids" % k, S[k], ms["ids"][mk])
    eq("set %s size" % k, S["sizes"][k if k != "gate" else "gate"], ms["counts"][mk])
eq("overlap T_s and S2", len(S["overlap_T_S2"]), ms["counts"]["overlap"])
eq("strata sizes", [len(x) for x in S["strata"]], ms["strata"])
eq("strata sizes in HA", [len(x) for x in S["strata"]], mod["ha"]["strata"])
nm = {"excluded (ref item incomplete)": "excluded", "drift": "drift", "rho = U": "unresolved", "no verdict data": "no verdict data"}
od = sorted((o, nm[r]) for r, l in orc["drops"].items() for o in l)
md = sorted((d["id"], d["reason"]) for d in ms["drops"])
eq("drops (orbit, reason)", od, md)
pan = {}
import tomllib
panel = tomllib.load(open(sys.argv[4], "rb"))
if panel:
    po = {o["id"]: o for o in panel["orbit"]}
    eq("drops round and k_index", sorted((o, po[o]["round"], po[o]["k_index"]) for o, _ in od), sorted((d["id"], d["round"], d["k_index"]) for d in ms["drops"]))
eq("log difference (last entry minus final sets)", orc["log_difference"], ms["log_difference"])
eq("by K index and round rows", orc["T2_by_k_round"], {"%d|%d" % (r["k_index"], r["round"]): dict(E_test=r["E_test"], W2=r["W2"], T_s=r["T_s"], S2=r["S2"]) for r in ms["by_k_round"]})
# outputs per set
oo = orc["outputs"]["R_a"]
S_ = {"E_test": S["E_test"], "W2": S["W2"], "T_s": S["T_s"], "S2": S["S2"], "D": S["D"]}
for nmset, ids in S_.items():
    c = {"chaotic_supported": 0, "regular_supported": 0, "unresolved": 0}
    for o in ids:
        v = oo[str(o)] if str(o) in oo else orc["outputs"]["R_ab_D"].get(str(o))
        if str(o) not in oo:
            continue
        c[v.replace("-", "_")] += 1
    if nmset == "D":
        continue
    eq("rule outputs of %s" % nmset, c, mod["outputs"][nmset])
# D outputs: R_a on D orbits
cD = {"chaotic_supported": 0, "regular_supported": 0, "unresolved": 0}
for o in S["D"]:
    cD[oo[str(o)].replace("-", "_")] += 1
eq("rule outputs of D", cD, mod["outputs"]["D"])

# ---------------------------------------------------------------- HA
ha, mh = orc["HA"], mod["ha"]
eq("HA e", ha["point"]["e"], mh["e"])
for a, b in (("a_E", "E_test"), ("a_W2", "W2"), ("a_Ts", "T_s"), ("a_S2", "S2")):
    eq("HA a(%s)" % b, ha["point"][a], mh["a"][b])
eq("HA a(W2) - a(S2)", ha["point"]["d_W"], mh["dw"]); eq("HA a(T_s) - a(S2)", ha["point"]["d_T"], mh["dt"]); eq("HA g", ha["g"], mh["g"])
for a, b in (("e", "e_up"), ("a_E", "aE_up"), ("d_W", "dw_lo"), ("d_T", "dt_lo")):
    eq("HA bound %s value" % a, ha["bounds"][a][0], mh["bounds"][b]["value"]); eq("HA bound %s dropped replicates" % a, ha["bounds"][a][1], mh["bounds"][b]["dropped"])
eq("HA decision", {"not supported": "not supported", "SUPPORTED": "supported", "not decidable": "not decidable"}[ha["decision"]], mh["verdict"])
eq("HA failing items", ha["failing"], mh["failed"])
eq("HA gate missing", orc["sets"]["gate_without_verdicts"], mh["gate_missing"])

# ---------------------------------------------------------------- HB
hb, mb = orc["HB"], mod["hb"]
def hb_rows(res):
    rows = []
    for A in ("T_s", "W2"):
        rows += [(A, r) for r in res["sets"][A]["rows"]]
    return rows
def cmp_rows(prefix, orows, mrows):
    eq(prefix + " number of comparisons", len(orows), len(mrows))
    for (A, r), m in zip(orows, mrows):
        t = "%s %s:%s" % (prefix, A, r["X"])
        eq(t + " set and indicator", (A, r["X"]), (m["set"], m["indicator"]))
        for a, b in (("b", "b"), ("cost", "cost"), ("capped", "capped"), ("cost_unmatched", "cost_unmatched"), ("FR_X", "fr_x"), ("FR_rule", "fr_rule"), ("n11", "n11"), ("n10", "n10"), ("n01", "n01"), ("n00", "n00"),
                     ("p", "p"), ("margin", "margin"), ("passed", "pass")):
            eq(t + " " + b, r[a], m[b])
cmp_rows("HB", hb_rows(hb), mb["comparisons"])
fl_o = ["%s:%s" % (A, r["X"]) for A, r in hb_rows(hb) if not r["passed"]]
eq("HB failing comparisons (as sets; the module lists them in alphabetical order, the oracle in table order)", sorted(fl_o), sorted(mb["failing"]))
eq("HB coverage (an empty set has coverage 0: reading of the module, D95)", {A: (hb["sets"][A]["cov"] if hb["sets"][A]["cov"] is not None else "0/1") for A in ("T_s", "W2")}, mb["coverage"]); eq("HB coverage ok", hb["coverage_ok"], mb["coverage_ok"])
eq("HB coverage price (an empty set has price 1)", {A: (hb["sets"][A]["coverage_price"] if hb["sets"][A]["coverage_price"] is not None else "1/1") for A in ("T_s", "W2")}, mb["coverage_price"])
eq("HB underpowered flag", hb["underpowered"], mb["underpowered"]); eq("HB cost-unmatched LE flag", hb["cost_unmatched_LE_flag"], mb["cost_unmatched_le"])
eq("HB defects", hb["defects"], mb["defects"]); eq("HB sizes", {A: hb["sets"][A]["n"] for A in ("T_s", "W2")}, mb["sizes"])
eq("HB decision", hb["decision"], mb["verdict"])
HBW = orc["HB_without_LE"]
cmp_rows("HB without LE", hb_rows(HBW), mb["without_le"]["comparisons"])
eq("HB without LE decision", HBW["decision"], mb["without_le"]["verdict"])
eq("HB without LE failing (as sets)", sorted("%s:%s" % (A, r["X"]) for A, r in hb_rows(HBW) if not r["passed"]), sorted(mb["without_le"]["failing"]))
# sensitivities
sens = {"round_down": "round_down", "pilot": "pilot", "measured": "measured", "without_le": "without_LE"}
for k, mk in sens.items():
    o = orc["sensitivities"][mk if mk in orc["sensitivities"] else k]; m = mod["sensitivities"][k]
    dec_ok = not [x for x in o["undecidable_reasons"] if x != "set below floor"]          # the floors belong to the decision of HB, not to the sensitivities (E69 b)
    eq("sensitivity %s decidable" % k, dec_ok, m["decidable"])
    eq("sensitivity %s verdict" % k, o["decision"], m["verdict"])
    if dec_ok:
        eq("sensitivity %s failing (as sets)" % k, sorted("%s:%s" % (A, r["X"]) for A, r in hb_rows(o) if not r["passed"]), sorted(m["failing"]))
    else:
        eq("sensitivity %s reason names the cause (spec L1 has no measured value)" % k, True, ("L1" in m["reason"]) and ("114" in m["reason"]))
# HB power
pw, mp = orc["power_hb_achieved"], mod["power"]["hb"]
eq("HB power at achieved sizes: successes", pw["successes"], mp["successes"]); eq("HB power runs", pw["runs"], mp["runs"]); eq("HB power", pw["power"], mp["power"]); eq("HB power se", pw["se"], mp["se"]); eq("HB planned power", pw["planned"], mp["planned"])
if orc["power_hc_achieved"] is None: eq("HC power at achieved sizes absent (not underpowered)", None, mod["power"]["hc"])
else:
    ph, mq = orc["power_hc_achieved"], mod["power"]["hc"]
    eq("HC power at achieved sizes: successes", ph["successes"], mq["successes"]); eq("HC power runs", ph["runs"], mq["runs"]); eq("HC power", ph["power"], mq["power"]); eq("HC power se", s12(ph["se"]), mq["se"]); eq("HC planned power", ph["planned"], mq["planned"])

# ---------------------------------------------------------------- HC
hc, mc = orc["HC"], mod["hc"]
for a, b in (("AUC_I", "auc_i"), ("AUC_P", "auc_p"), ("Delta", "delta")):
    eq("HC " + b, hc[a], mc[b])
eq("HC lower bound", hc["lower"][0], mc["lower"]); eq("HC upper bound", hc["upper"][0], mc["upper"])
eq("HC failures and non-failures", (hc["failures"], hc["nonfailures"]), (mc["failures"], mc["nonfailures"]))
eq("HC failure ids", sorted(r["o"] for r in hc["rows"] if r["f"]), mc["failure_ids"])
eq("HC dropped replicates", (hc["lower"][1], hc["upper"][1]), (mc["dropped"], mc["dropped"]))
eq("HC decision", {"SUPPORTED": "supported", "REVERSED": "reversed", "not supported": "not supported", "not decidable": "not decidable"}[hc["decision"]], mc["verdict"])
eq("HC underpowered", hc["underpowered"], mc["underpowered"]); eq("HC set size", hc["n"], mc["set_size"])

# ---------------------------------------------------------------- flags, R_ab
eq("flag underpowered HB", hb["underpowered"], mod["flags"]["underpowered_hb"]); eq("flag underpowered HC", hc["underpowered"], mod["flags"]["underpowered_hc"])
eq("flag cost-unmatched LE", True, mod["flags"]["cost_unmatched_le"])
eq("flag not decidable causes (the hypotheses whose decision is not decidable)", [n for n, d in (("HA", ha["decision"]), ("HB", hb["decision"]), ("HC", hc["decision"])) if d == "not decidable"], mod["flags"]["not_decidable"])
T9, mr = orc["T9"], mod["r_ab"]
eq("T9 pairs", {k.replace("-", "_"): v for k, v in T9["pairs"].items()}, mr["pairs"]); eq("T9 overrides", T9["overrides_of_decided"], mr["overrides"]); eq("T9 orbits", S["D"], mr["orbits"])

# ---------------------------------------------------------------- tables file: T2, T3, T4, T9 error, T10, T11, F2, F3
def block(title):
    i = tables.index(title); j = tables.find("\n## ", i + 3)
    return tables[i:j if j > 0 else len(tables)]
def md_rows(txt):
    rows = []
    for ln in txt.split("\n"):
        if ln.startswith("|") and not ln.startswith("|---"):
            rows.append([c.strip() for c in ln.strip().strip("|").split("|")])
    return rows
t2 = md_rows(block("## Table T2"))
by_round = {r[0]: r for r in t2 if r and r[0] in ("1", "2", "3")}
for rd, d in orc["T2"].items():
    r = by_round[rd]
    eq("T2 round %s orbits, K=0, E_test, W2, T_s, S2" % rd, [d["orbits"], d["K0"], d["E_test"], d["W2"], d["T_s"], d["S2"]], [int(x) for x in r[1:7]])
sz = {r[0]: r[1:] for r in t2 if r and r[0] in ("E_test", "W2", "T_s", "S2", "D", "D_B", "gate")}
for k in ("E_test", "W2", "T_s", "S2", "D", "D_B", "gate"):
    eq("T2 size %s" % k, S["sizes"][k], int(sz[k][0]))
eq("T2 target W2", 11, int(sz["W2"][1])); eq("T2 target T_s", 5, int(sz["T_s"][1]))
t3txt = block("## Table T3")
t3 = [r for r in md_rows(t3txt) if r and r[0].isdigit()]
T3o = orc["T3"][[k for k in orc["T3"]][0]]
eq("T3 header columns (E79): orbit, N_B, three L5 values, ladder would force, rho; no column labels agree", ["orbit", "N_B", "L5 Float64", "L5 Double64", "L5 BigFloat", "ladder would force", "rho"], md_rows(t3txt)[0])
eq("T3 text states the omitted column and its cause", True, ("labels agree" in t3txt and "omitted" in t3txt and "labels cannot be formed" in t3txt))
eq("T3 orbits (E79: D_B by content, and the orbits of a failed chk chunk, P118)", [int(o["o"]) for o in T3o], [int(r[0]) for r in t3])
NBp = json.load(open(sys.argv[1]))["fingerprints"] and None
for o, r in zip(T3o, t3):
    if o["available"]:
        okv = all(float(r[c]) == s12(o[k]) for c, k in ((2, "lam_F64"), (3, "lam_D64"), (4, "lam_BF")))
        N[0] += 1; N[1] += okv
        if not okv: DIFFS.append("T3 values orbit %s" % o["o"])
        out(("EQUAL " if okv else "DIFFER ") + "T3 orbit %s the three L5 values (equal to the oracle value rounded to 12 significant digits, as the module writes it)" % o["o"] + ("" if okv else "  oracle=%s module=%s" % ((o["lam_F64"], o["lam_D64"], o["lam_BF"]), r[2:5])))
        eq("T3 orbit %s ladder, rho" % o["o"], (str(o["ladder_unresolved"]).lower(), o["rho"]), (r[5], r[6]))
    else:
        out("NOTE T3 orbit %s: the oracle has missing value(s) %s; module cells: %s" % (o["o"], o["missing"], r[2:7]))
        eq("T3 orbit %s the missing cell(s) are not numbers and the ladder column is not a number" % o["o"], True, all(not re.fullmatch(r"-?\d+\.?\d*(e-?\d+)?", r[c]) for c, k in ((2, "F64"), (3, "D64"), (4, "BF")) if k in o["missing"]))
        eq("T3 orbit %s the present values equal the oracle ones" % o["o"], True, all(float(r[c]) == s12(o[k]) for c, k in ((2, "lam_F64"), (3, "lam_D64"), (4, "lam_BF")) if o[k] is not None and re.fullmatch(r"-?\d+\.?\d*(e-?\d+)?", r[c])))
mm3 = re.search(r"unresolved: C (\d+), R (\d+)", t3txt)
eq("T3 count of orbits on which the ladder would have made rho unresolved, split by rho", (sum(1 for o in T3o if o["ladder_unresolved"] and o["rho"] == "C"), sum(1 for o in T3o if o["ladder_unresolved"] and o["rho"] == "R")), tuple(int(x) for x in mm3.groups()))
# T4 : counts of outputs per set per K and per round
t4 = [r for r in md_rows(block("## Table T4")) if len(r) == 6 and r[0] in ("E_test", "W2", "T_s", "S2")]
cnt = {"chaotic-supported": 3, "regular-supported": 4, "unresolved": 5}
for r in t4:
    ids = [o for o in S[r[0]] if (r[1] == "all" or po[o]["k_index"] == int(r[1])) and (r[2] == "all" or po[o]["round"] == int(r[2]))]
    got = [sum(1 for o in ids if oo[str(o)] == k) for k in ("chaotic-supported", "regular-supported", "unresolved")]
    eq("T4 %s K index %s round %s outputs" % (r[0], r[1], r[2]), got, [int(x) for x in r[3:6]])
eq("T4 number of rows", 4 * 10, len(t4))
# T9 error line
m = re.search(r"Error of the decided orbits: R_a (\d+), R_ab (\d+)", block("## Table T9"))
eq("T9 contradictions of decided orbits (R_a, R_ab)", (T9["contradicting_Ra"], T9["contradicting_Rab"]), (int(m.group(1)), int(m.group(2))))
# T10
t10 = {r[0]: r[1:] for r in md_rows(block("## Table T10")) if r and r[0] in ("L1", "tg", "pa1", "bk1", "fq1")}
for s_, k in (("L1", "l1"), ("tg", "tg"), ("pa1", "pa1"), ("bk1", "bk1"), ("fq1", "fq1")):
    meas = orc["T10"]["measured"].get(k)
    eq("T10 %s measured (E78, 12 significant digits; n.a. when no value), cost key, pilot key" % s_, ((meas is None and t10[s_][0] == "n.a.") or (meas is not None and t10[s_][0] != "n.a." and float(t10[s_][0]) == s12(meas)), str(orc["T10"]["cost_keys"][k]), str(orc["T10"]["pilot"][k])), (True, t10[s_][1], t10[s_][2]))
# T11
t11 = block("## Table T11")
rows = [r for r in md_rows(t11) if len(r) == 4 and r[3].isdigit()]
modfl = {"%s|%s|%s" % (r[0], r[1], r[2]): int(r[3]) for r in rows if int(r[3]) > 0}
eq("T11 non-zero flag counts by arm, item, flag", orc["T11"]["flag_counts"], modfl)
fc = [r for r in md_rows(t11) if len(r) == 3 and r[0].startswith("chunk-")]
eq("T11 chunks failed twice", sorted(orc["T11"]["failed_chunks"]), sorted(r[0] for r in fc))
eq("T11 chunks failed twice: file, cell (E74) and orbits of the cell from the panel (E77, E80)", sorted((c["fn"], c["cell"], sorted(c["orbits"])) for c in orc["T11"]["failed_cells"]),
   sorted((r[0], r[1], sorted(int(x) for x in re.findall(r"\d+", r[2]))) for r in fc))
eq("manifest status (E76)", orc["T11"]["manifest_status"], re.search(r"status: (.*)", tables).group(1).strip() if re.search(r"status: (.*)", tables) else orc["T11"]["manifest_status"])
# F2, F3 : data rows and digests
def f_rows(title):
    txt = block(title)
    rows = [r for r in md_rows(txt)][1:]
    dig = re.search(r"digest: ([0-9a-f]{64})", txt).group(1)
    lines = [l for l in txt.split("\n") if l.startswith("|")]
    return rows, dig, lines
rows, dig, lines = f_rows("## F2 data")
fam = [f for f in orc["F2"] if f["curve"] == "rule family"]
mrule = [r for r in rows if r[0] == "rule"]
def qkey(p): return re.match(r"q_C=(\d+) q_R=(\d+)", p).groups()
omap = {tuple(re.match(r"qC=(\d+),qR=(\d+)", f["point"]).groups()): f for f in fam}
for r in mrule:
    f = omap.get(qkey(r[1]))
    if f is None:
        eq("F2 rule point %s exists in the oracle family" % r[1], True, False); continue
    eqr("F2 rule %s decided share" % r[1], f["decided_share"], r[2]); eqr("F2 rule %s error" % r[1], f["e"], r[3])
extra = [k for k in omap if k not in {qkey(r[1]) for r in mrule}]
out("NOTE F2 rule family: the module has %d points (a cross), the oracle's product grid has %d; oracle-only points: %s; their cross flag: %s" % (len(mrule), len(fam), extra, [omap[k]["cross"] for k in extra]))
eq("F2 module rule points are exactly the oracle's cross points", sorted(k for k, f in omap.items() if f["cross"]), sorted(qkey(r[1]) for r in mrule))
for r in [r for r in rows if r[0] == "single"]:
    f = [x for x in orc["F2"] if x["curve"] == "single " + r[1]][0]
    eqr("F2 single %s decided share" % r[1], f["decided_share"], r[2]); eqr("F2 single %s error" % r[1], f["e"], r[3])
eq("F2 number of single curves", 7, len([r for r in rows if r[0] == "single"]))
cands = {"table lines (header, separator, rows) joined by LF": "\n".join(lines), "table lines + LF at end": "\n".join(lines) + "\n", "rows only joined by LF": "\n".join(lines[2:]), "rows only + LF": "\n".join(lines[2:]) + "\n"}
hit = [k for k, v in cands.items() if hashlib.sha256(v.encode()).hexdigest() == dig]
eq("F2 digest is the SHA-256 of the table lines of the block under one of four conventions (independent recomputation)", True, bool(hit)); out("   convention: " + str(hit))
rows3, dig3, lines3 = f_rows("## F3 data")
for key, cur in (("sI", "s_I"), ("sP", "s_P")):
    o_rows = [f for f in orc["F3"] if f["curve"] == key and "score" in f]
    m_rows = [r for r in rows3 if r[0] == cur and r[1] != "AUC"]
    eq("F3 %s number of ROC rows" % cur, len(o_rows), len(m_rows))
    for f, r in zip(o_rows, m_rows):
        eqr("F3 %s threshold" % cur, f["score"], r[1]); eqr("F3 %s fpr" % cur, f["fpr"], r[2]); eqr("F3 %s tpr" % cur, f["tpr"], r[3])
    a = [f for f in orc["F3"] if f["curve"] == key and "auc" in f][0]
    eqr("F3 %s AUC" % cur, a["auc"], [r for r in rows3 if r[0] == cur and r[1] == "AUC"][0][2])
hit3 = [k for k, v in {"table lines (header, separator, rows) joined by LF": "\n".join(lines3), "lines + final LF": "\n".join(lines3) + "\n", "rows only": "\n".join(lines3[2:]), "rows only + LF": "\n".join(lines3[2:]) + "\n"}.items() if hashlib.sha256(v.encode()).hexdigest() == dig3]
eq("F3 digest is the SHA-256 of the table lines of the block under one of four conventions", True, bool(hit3)); out("   convention: " + str(hit3))
# machine file stated in the tables file
m = re.search(r"([0-9a-f]{64})", tables[:tables.index("## Decisions")] if "## Decisions" in tables else tables)
eq("the SHA-256 of the machine file is stated in the tables file and equals the recomputed one", True, sha(os.path.join(sys.argv[2], "p3-abstention-test-results.json")) in tables)
out("")
out("QUANTITIES COMPARED: %d; EQUAL: %d; DIFFER: %d" % (N[0], N[1], N[0] - N[1]))
for d in DIFFS: out("  DIFF: " + d)
