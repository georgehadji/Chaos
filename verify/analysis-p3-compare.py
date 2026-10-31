# Comparison of the independent recomputation (verify/analysis-p3-oracle.py) with the calibration record of the module
# (p3-abstention-calibration-record.json, keys of E58 and the layout of the D93 section of the test plan) and, when given, its tables file.
# Python 3 standard library only.
#
# Usage:  python verify/analysis-p3-compare.py <oracle json> <record json> [<tables md>]
# Rules (D93r): integers, rationals, sets, texts and flags exactly; reals to 1e-12 relative (the record writes 12 significant digits, so a record
# value equal to the oracle value rounded to 12 significant digits is "equal to 12 digits"); HC powers (and the PP quantities that use the
# platform logarithm) within 3 Monte Carlo standard errors.
# Prints one line per quantity: EQUAL, or DIFF with both values. Exit code 0 iff no DIFF.
import sys, json, math, hashlib
from fractions import Fraction

def frac(x):
    if isinstance(x, str):
        return Fraction(x)
    return Fraction(repr(float(x))) if isinstance(x, float) else Fraction(x)

class Cmp:
    def __init__(self):
        self.n = 0; self.bad = 0; self.notes = []
    def _line(self, ok, name, a, b, how="exact"):
        self.n += 1
        if not ok:
            self.bad += 1
        sa, sb = json.dumps(a, default=str)[:160], json.dumps(b, default=str)[:160]
        print(("EQUAL " if ok else "DIFF  ") + name + ("" if ok else "  oracle=%s  record=%s" % (sa, sb)) + ("  [%s]" % how if ok and how != "exact" else ""))
    def exact(self, name, a, b):
        self._line(a == b, name, a, b)
    def rat(self, name, a, b):
        try:
            ok = frac(a) == frac(b)
        except Exception:
            ok = False
        self._line(ok, name, a, b)
    def real(self, name, a, b):
        if a is None or b is None:
            self._line(a is None and b is None, name, a, b); return
        a, b = float(a), float(b)
        rel = abs(a - b) / max(abs(a), 1e-300)
        if a == b or rel <= 1e-12:
            self._line(True, name, a, b, "1e-12 relative")
        elif float("%.12g" % a) == b or abs(a - b) <= 1e-15:
            self._line(True, name, a, b, "equal to 12 significant digits")
        else:
            self._line(False, name, a, b)
    def stat(self, name, succ_a, succ_b, runs):
        pa, pb = succ_a / runs, succ_b / runs
        se = math.sqrt(max(pa * (1 - pa), pb * (1 - pb), 1e-12) / runs)
        if succ_a == succ_b:
            self._line(True, name, succ_a, succ_b, "equal counts")
        else:
            self._line(abs(pa - pb) <= 3 * se, name, succ_a, succ_b, "within 3 standard errors")

def main():
    orc = json.load(open(sys.argv[1], encoding="utf-8"))
    rec = json.load(open(sys.argv[2], encoding="utf-8"))
    tables = open(sys.argv[3], encoding="utf-8").read() if len(sys.argv) > 3 else None
    C = Cmp()
    srcname = {"P": "panel", "G": "grid"}
    # ------------------------------------------------------------------ parameters and fingerprints
    C.exact("outcome.selection", orc["outcome"], {"selected": "selected"}.get(rec["outcome"]["selection"], rec["outcome"]["selection"]))
    # ------------------------------------------------------------------ relabelling (E49)
    rl, rr = orc["relabelling"], rec["relabelling"]
    mine = {}
    for k, v in rl["rho_new_by_source_kidx"].items():
        s, ki = k.split(",")
        for rho, c in v.items():
            if c:
                mine[(srcname[s], int(ki), rho)] = c
    theirs = {(d["source"], d["k_index"], d["rho"]): d["count"] for d in rr["rho_new"] if d["count"]}
    C.exact("relabelling.rho_new by (source, K index, rho)", sorted(map(str, mine.items())), sorted(map(str, theirs.items())))
    ax = {}
    for k, v in rl["axes_counts"].items():
        s, a, b = k.split("|"); ax[(a, b)] = ax.get((a, b), 0) + v
    C.exact("relabelling.axes (instability, structure) pooled", sorted(map(str, ax.items())), sorted(str(((d["instability"], d["structure"]), d["count"])) for d in rr["axes"] if d["count"]))
    for name, key in (("panel", "changes_panel"), ("grid", "changes_grid")):
        mine = {tuple(k.split(",")): v for k, v in rl["change_" + name].items()}
        mine = {(int(a), b, c): v for (a, b, c), v in mine.items()}
        theirs = {(d["k_index"], d["old"], d["new"]): d["count"] for d in rr[key]}
        C.exact("relabelling.change table " + name + " by (K index, old, new)", sorted(map(str, mine.items())), sorted(map(str, theirs.items())))
    C.exact("relabelling.panel_changed", rl["panel_changed"], rr["panel_changed"])
    C.exact("relabelling.grid_changed (D93d: 0 against D89 E23 b)", rl["grid_changed"], rr["grid_changed"])
    C.exact("relabelling.missing data", len(rl["missing"]), len(rr["missing_data"]))
    d1, d2 = orc["dropped"], rr["dropped"]
    C.exact("dropped.total", d1["total"], d2["total"])
    C.exact("dropped.by_reason", {k: v["total"] for k, v in d1["by_reason"].items()}, d2["by_reason"])
    C.exact("dropped rows (id, source, K index, reason)", sorted((r["id"], srcname[r["src"]], r["kidx"], r["reason"]) for r in d1["rows"]), sorted((r["id"], r["source"], r["k_index"], r["reason"]) for r in d2["rows"]))
    # ------------------------------------------------------------------ sets (E50)
    for nm in ("cal", "S2", "T_s", "W2", "R_W", "regular", "gate"):
        mk = {"cal": "Cal"}.get(nm, nm)
        for part in ("pooled", "panel", "grid"):
            C.exact("sets.%s.%s" % (nm, part), orc["sets"][mk][part], rec["sets"][nm][part])
    C.exact("sets.overlap S2 with W2", orc["sets"]["overlap_S2_W2"], rec["sets"]["overlap"])
    mine = {(srcname[k.split(",")[0]], int(k.split(",")[1])): v for k, v in orc["sets"]["Cal"]["by_source_kidx"].items()}
    C.exact("sets.cal by source and K index", sorted(map(str, mine.items())), sorted(str(((d["source"], d["k_index"]), d["count"])) for d in rec["sets"]["cal_by_k"]))
    # ------------------------------------------------------------------ space, F1..F6, selection (E52 to E54)
    ts_o, ts_r = orc["theta_space"], rec["theta_space"]
    C.exact("theta_space.size", ts_o["size"], ts_r["size"])
    C.exact("theta_space.F1..F6 alone", [ts_o["alone"]["F%d" % k] for k in range(1, 7)], ts_r["alone"])
    C.exact("theta_space.first failing item counts (0 = feasible, 1..6)", [ts_o["first_fail_counts"][str(k)] for k in range(0, 7)], ts_r["first_fail"])
    C.rat("theta_space.cost cap of F1", orc["cost_cap"], ts_r["cost_cap"])
    C.exact("theta_space.sizes Cal, S2, T_s, W2", [orc["sets"]["Cal"]["pooled"], orc["sets"]["S2"]["pooled"], orc["sets"]["T_s"]["pooled"], orc["sets"]["W2"]["pooled"]],
            [ts_r["sizes"][k] for k in ("cal", "S2", "Ts", "W2")])
    so, sr = orc["selection"], rec["selection"]
    C.exact("selection.|F|", so["feasible_count"], sr["feasible"])
    C.exact("selection.|Theta_tol|", so["tol_size"], sr["tol"])
    C.rat("selection.J_min", so["J_min"], sr["J_min"])
    C.exact("selection.theta* index", so["star"]["index"], sr["index"])
    C.exact("selection.theta* (V, qC, qR, w, v)", [so["star"]["V"], so["star"]["qC"], so["star"]["qR"], so["star"]["w"], so["star"]["v"]], [sr["theta"]["V"], sr["theta"]["qC"], sr["theta"]["qR"], sr["theta"]["w"], sr["theta"]["v"]])
    C.rat("selection.theta* cost C_R(theta*, b_P)", so["star"]["C_R"], sr["star"]["cost"])
    m_o, m_r = orc["matched"], sr["star"]
    for X in ("SALI", "FLI", "MEGNO", "E1", "E2", "WBA", "LE"):
        a = m_o["indicators"][X]; b = m_r["matched"][X]
        C.exact("matched budget %s (b, capped, round-down, round-down none)" % X, [a["b"], a["capped"], a["round_down"], a["round_down_none"]], [b["b"], b["capped"], b["round_down"], b["round_down_none"]])
    C.exact("LE read budget", m_o["le_read"], m_r["le_read"])
    C.exact("LE cost-unmatched flag", m_o["le_cost_unmatched"], m_r["le_cost_unmatched"])
    for k1, k2 in (("pair", "pair"), ("float_votes", "float"), ("ratio", "ratio")):
        C.rat("Double64 costs " + k1, m_o["double64"][k1], m_r["d64_costs"][k2])
    # F1 data block of the tables file
    if tables is not None:
        i0 = tables.index("## F1 data")
        lines = [l for l in tables[i0:].split("\n") if l.startswith("|")]
        rows = [[c.strip() for c in l.strip("|").split("|")] for l in lines[2:]]
        mine = [(r["index"], r["C_R"], r["FR_W2"], r["FR_Ts"], r["J"], r["in_tol"], r["star"]) for r in so["F_rows"]]
        theirs = [(int(r[0]), r[2], r[3], r[4], r[5], r[6] == "true", r[7] == "true") for r in rows]
        ok = len(mine) == len(theirs) and all(a[0] == b[0] and frac(a[1]) == frac(b[1]) and frac(a[2]) == frac(b[2]) and frac(a[3]) == frac(b[3]) and frac(a[4]) == frac(b[4]) and a[5] == b[5] and a[6] == b[6] for a, b in zip(mine, theirs))
        C._line(ok, "F1 data block rows (index, C_R, FR(W2), FR(T_s), J, in Theta_tol, theta*), %d rows" % len(mine), len(mine), len(theirs))
        dig = [l for l in tables[i0:].split("\n") if l.startswith("digest:")]
        block = [l for l in tables[i0:].split("\n") if l.startswith("|")]
        cand = {"lines joined by LF": "\n".join(block), "lines joined by LF with a final LF": "\n".join(block) + "\n"}
        got = {k: hashlib.sha256(v.encode("utf-8")).hexdigest() for k, v in cand.items()}
        C._line(bool(dig) and dig[0].split()[1] == got["lines joined by LF"], "F1 data block digest = SHA-256 of the table lines (header, separator, rows) joined by LF (independent recomputation)", got["lines joined by LF"], dig[0].split()[1] if dig else None)
    # ------------------------------------------------------------------ C7
    c7o, c7r = orc["c7"], rec["c7"]
    C.exact("c7.top 20 of Theta_tol (indices in the order of C6)", [t["index"] for t in so["tol_top20"]], c7r["top"])
    for A in ("T_s", "W2"):
        ro, rr_ = c7o["hb"][A], c7r["hb"][A]
        C.exact("c7.hb %s comparators" % A, [r["X"] for r in ro], [r["X"] for r in rr_])
        for a, b in zip(ro, rr_):
            tag = "c7.hb %s %s" % (A, a["X"])
            C.exact(tag + " b, capped, cost-unmatched", [a["b_read"], a["capped"], a["cost_unmatched"]], [b["b"], b["capped"], b["cost_unmatched"]])
            if "FR_X" in a:
                C.rat(tag + " cost", a["cost"], b["cost"])
                C.exact(tag + " FR_X, FR_rule, n10, n01, coverage, margin", [frac(a["FR_X"]), frac(a["FR_rule"]), a["n10"], a["n01"], frac(a["coverage"]), a["margin"]],
                        [frac(b["fr_x"]), frac(b["fr_rule"]), b["n10"], b["n01"], frac(b["coverage"]), b["margin"]])
                C.rat(tag + " exact McNemar p", a["p"], b["p"])
    st = so["star"]
    C.rat("c7.e", st["e"], c7r["e"]); C.rat("c7.u_cal", st["u_cal"], c7r["u_cal"]); C.rat("c7.u_s2", st["u_S2"], c7r["u_s2"]); C.exact("c7.g", st["g"], c7r["g"])
    C.rat("c7.J", st["J"], c7r["J"]); C.rat("c7.cost", st["C_R"], c7r["cost"])
    fam_o = c7o["risk_coverage"]["family_cross"]; fam_g = c7o["risk_coverage"]["family_grid"]
    fam_r = {"%d,%d" % (d["qC"], d["qR"]): d for d in c7r["risk_coverage"]["family"]}
    which = "family_cross" if set(fam_r) == set(fam_o) else ("family_grid" if set(fam_r) == set(fam_g) else None)
    C.exact("c7.risk_coverage family points (record = cross: q_C and q_R varied one at a time; oracle also computed the grid)", sorted(fam_r), sorted(fam_o))
    for k in sorted(fam_r):
        src = fam_o.get(k) or fam_g.get(k)
        if src:
            C.exact("c7.risk_coverage family %s (coverage, e)" % k, [frac(src["decided_share"]), frac(src["e"]) if src["e"] else None], [frac(fam_r[k]["coverage"]), frac(fam_r[k]["e"]) if fam_r[k]["e"] else None])
    for d in c7r["risk_coverage"]["single"]:
        a = c7o["risk_coverage"]["single"][d["X"]]
        C.exact("c7.risk_coverage single %s (b, coverage, e)" % d["X"], [a["b"], frac(a["decided_share"]), frac(a["e"]) if a["e"] else None], [d["b"], frac(d["coverage"]), frac(d["e"]) if d["e"] else None])
    nsrc = {"panel": orc["sets"]["Cal"]["panel"], "grid": orc["sets"]["Cal"]["grid"]}
    for s in ("panel", "grid"):
        b = c7r["by_source"][s]
        C.exact("c7.by_source %s (n, FR T_s, FR W2, decided T_s, decided W2)" % s,
                [nsrc[s], frac(c7o["by_source"]["T_s," + s]["FR"]), frac(c7o["by_source"]["W2," + s]["FR"]), frac(c7o["by_source"]["T_s," + s]["decided"]), frac(c7o["by_source"]["W2," + s]["decided"])],
                [b["n"], frac(b["fr_ts"]), frac(b["fr_w2"]), frac(b["dec_ts"]), frac(b["dec_w2"])])
    kp = c7r["by_source"]["k_positive"]
    C.exact("c7.K > 0 (n, e, u)", [c7o["K_positive"]["n"], frac(c7o["K_positive"]["e"]), frac(c7o["K_positive"]["u_cal"])], [kp["n"], frac(kp["e"]), frac(kp["u"])])
    # ------------------------------------------------------------------ R_ab (P63)
    nm = lambda s: s.replace("-", "_")
    C.exact("r_ab rows (id, output) on the Cal orbits with Double64 items", [(r["id"], nm(r["R_ab"])) for r in orc["r_ab"]["rows"]], [(r["id"], r["output"]) for r in rec["r_ab"]])
    # ------------------------------------------------------------------ PP for HB
    oo, orr = orc["optimism"], rec["optimism"]
    C.exact("optimism.used, skipped", [oo["used"], oo["skipped"]], [orr["used"], orr["skipped"]])
    C.exact("outcome.optimism_skipped", oo["skipped"], rec["outcome"]["optimism_skipped"])
    C.exact("optimism.theta_b (index of the element chosen on each resample)", oo["theta_b"], orr["thetas"])
    C.real("optimism(T_s)", oo["optimism_T_s"], orr["optimism"]["T_s"])
    C.real("optimism(W2)", oo["optimism_W2"], orr["optimism"]["W2"])
    ho, hr = orc["pp_hb"], rec["pp_hb"]["rates"]
    C.real("pp_hb.r(theta*, T_s)", ho["r_Ts"], hr["r"]["T_s"]); C.real("pp_hb.r(theta*, W2)", ho["r_W2"], hr["r"]["W2"])
    C.real("pp_hb.r_corr T_s", ho["r_corr_T_s"], hr["r_corr"]["T_s"]); C.real("pp_hb.r_corr W2", ho["r_corr_W2"], hr["r_corr"]["W2"]); C.real("pp_hb.r_corr R_W", ho["r_corr_R_W"], hr["r_corr"]["R_W"])
    C.real("pp_hb.p_rule T_s", ho["p_rule_T_s"], hr["p_rule"]["T_s"]); C.real("pp_hb.p_rule R_W", ho["p_rule_R_W"], hr["p_rule"]["R_W"])
    pool = ho["rates"]["pooled"]
    sname = {"T": "T_s", "RW": "R_W"}
    for S in ("T", "RW"):
        C.real("pp_hb.d(%s)" % sname[S], float(frac(pool["d"][S])), hr["d"][sname[S]])
        for X in ("LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA"):
            C.real("pp_hb.p_X(%s) X=%s" % (sname[S], X), float(frac(pool["pX"][S][X])), hr["p_comp"]["%s|%s" % (X, sname[S])])
            C.real("pp_hb.s_X(%s) X=%s" % (sname[S], X), float(frac(pool["sX"][S][X])), hr["s_comp"]["%s|%s" % (X, sname[S])])
    po, pr = ho["path"], rec["pp_hb"]["path"]
    C.exact("pp_hb.path sizes (j, rho, nT, nRW, nW, allowed), 17 points", [(r["j"], r["rho"], r["nT"], r["nRW"], r["nW"], r["allowed"]) for r in po], [(r["j"], r["rho"], r["nT"], r["nRW"], r["nW"], r["allowed"]) for r in pr])
    for a, b in zip(po, pr):
        if a["allowed"]:
            C.exact("pp_hb.path point %d power counts (successes of %d runs)" % (a["j"], a["runs"]), a["successes"], b["successes"])
    t_o, t_r = ho["targets"], rec["pp_hb"]["targets"]
    C.exact("pp_hb.targets outcome, cap flag", [t_o["outcome"], t_o["cap"]], [t_r["outcome"], t_r["cap"]])
    C.exact("pp_hb.targets n_T, n_W", [t_o["n_T"], t_o["n_W"]], [t_r["n_T"], t_r["n_W"]])
    if t_o.get("cap"):
        C.exact("pp_hb.targets point of the cap (j; rho* is null in the cap case, D93 review fix)", [t_o["j"], None], [t_r["j"], t_r["rho_star"]])
        C.real("pp_hb.targets power at the cap", t_o["power_at_cap"], t_r["power_cap"])
    else:
        C.exact("pp_hb.targets rho*", t_o["rho_star"], t_r["rho_star"])
    ty_o, ty_r = ho["type1"], rec["pp_hb"]["type1"]
    for X, v in ty_o["per_comparator"].items():
        C.real("pp_hb.type-I joint rate X=%s" % X, v["joint"], ty_r["rates"][X]["joint"])
        C.real("pp_hb.type-I single-comparison rate X=%s (reading: comparisons of X alone, coverage ignored)" % X, v["single"], ty_r["rates"][X]["single"])
    C.real("pp_hb.type-I max joint rate", ty_o["max_joint"], ty_r["max_joint"])
    # ------------------------------------------------------------------ PP for HC
    co, cr = orc["pp_hc"], rec["pp_hc"]
    C.exact("pp_hc.failures, non-failures", [co["n_fail"], co["n_nonfail"]], [cr["n_fail"], cr["n_nonfail"]])
    C.real("pp_hc.AUC_I_cal", co["AUC_I_cal_float"], cr["auc_i_cal"]); C.real("pp_hc.AUC_I_used", float(frac(co["AUC_I_used"])), cr["auc_i_used"])
    C.real("pp_hc.AUC_P", float(frac(co["AUC_P"])), cr["auc_p"]); C.real("pp_hc.true difference", co["true_difference"], cr["true_difference"])
    C.real("pp_hc.d(AUC_I_used)", co["d_I"], cr["d_i"]); C.real("pp_hc.d(AUC_P)", co["d_P"], cr["d_p"])
    C.exact("pp_hc.deliverable failures expected, worst", [co["deliverable_expected"], co["deliverable_worst"]], [cr["deliverable_expected"], cr["deliverable_worst"]])
    tab = [r for r in co["table"] if r["on_list"]]
    C.exact("pp_hc.n_f list", [r["nf"] for r in tab], [r["nf"] for r in cr["table"]])
    for a, b in zip(tab, cr["table"]):
        C.stat("pp_hc.power counts n_f=%d (successes of %d runs)" % (a["nf"], a["runs"]), a["successes"], b["successes"], a["runs"])
    C.stat("pp_hc.power at the expected deliverable number n_f=%d" % co["deliverable_expected"], next(r["successes"] for r in co["table"] if r["nf"] == co["deliverable_expected"]), cr["power_expected"]["successes"], co["table"][0]["runs"])
    tgt = co["target"]
    C.exact("pp_hc.target n_f (a check, not exact by the spec)", tgt["n_f"] if isinstance(tgt, dict) else tgt, cr["target"])
    C.exact("pp_hc.underpowered flag", co["underpowered"], cr["underpowered"])
    print("\n%d quantities compared, %d DIFF" % (C.n, C.bad))
    sys.exit(0 if C.bad == 0 else 1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python verify/analysis-p3-compare.py <oracle json> <record json> [<tables md>]")
    main()
