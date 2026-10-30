"""PDF measurements for the Amendment D91 verification (figures F1 to F6 of study p2-indicators).

Written from the spec (workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D91": E38 to E45, O12, O13,
criteria D91a, D91f, D91g, D91h, D91i, D91k, D91e) and from the tables file of the accepted stage 04. It reads nothing from src/, ext/ or test/.
Run (PyMuPDF prints a deprecation warning on import, harmless):
    python -I verify/analysis-p2-figures.py figs  <dir> <tables.md> <params.toml>    D91a, D91f, D91g, D91h, D91i on the six PDFs of <dir>
    python -I verify/analysis-p2-figures.py cmp   <dirA> <dirB> <name> [<name> ...]  pixel comparison of the same-named PDFs at 72 dpi
    python -I verify/analysis-p2-figures.py cmp1  <fileA> <fileB>                    same for two files
    python -I verify/analysis-p2-figures.py unzip <zip> <dir>                         extract a git archive (zip)
    python -I verify/analysis-p2-figures.py sig   <file> [<file> ...]                signature, size, pages, page size
    python -I verify/analysis-p2-figures.py contact <dir> <outdir>                    one PNG per figure at 72 dpi (analysis-d91-F1.png ...)
    python -I verify/analysis-p2-figures.py run <stdout-file> <result-file> <cwd> -- <command ...>   timed run with peak working set (psutil)
Lines for the log are printed as "@@ id | quantity | value | bound | PASS/FAIL/INFO" so the Julia script can count them.
"""
import sys, os, re, json, math, hashlib, time, tomllib, warnings, subprocess
warnings.filterwarnings("ignore")
import numpy as np
import fitz

fitz.TOOLS.mupdf_display_errors(False)
sys.stdout.reconfigure(encoding="utf-8")
IND = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2", "WBA", "RLI"]
TOL = 0.005  # 0.5 percent of the axis span (O12 b)


def emit(i, q, v, b, verdict):
    print(f"@@ {i} | {q} | {v} | {b} | {verdict}", flush=True)


def pf(ok):
    return "PASS" if ok else "FAIL"


# ------------------------------------------------------------------------------------------------------------ blocks (own parse of the tables file)
def typed(s):
    if re.fullmatch(r"-?\d+", s):
        return int(s)
    try:
        return float(s)
    except ValueError:
        return None if s == "nothing" else s


def load_blocks(path):
    text = open(path, "rb").read().decode().replace("\r\n", "\n")
    out = {}
    for m in re.finditer(r"## (F\d) data[^\n]*\n\n((?:\|[^\n]*\n)+)", text):
        lines = m.group(2).rstrip("\n").split("\n")
        cols = [c.strip() for c in lines[0].strip().strip("|").split("|")]
        rows = [[typed(c.strip()) for c in l.strip().strip("|").split("|")] for l in lines[2:]]
        out[m.group(1)] = (cols, rows)
    return out


def as_dicts(blk):
    cols, rows = blk
    return [dict(zip(cols, r)) for r in rows]


def load_params(path):
    t = tomllib.load(open(path, "rb"))
    return dict(budgets=t["budgets"]["value"], b_primary=t["b_primary"]["value"], N_c=t["N_c"]["value"], min_set=t["min_set"]["value"])


# ------------------------------------------------------------------------------------------------------------ page reading
class Page:
    def __init__(s, path):
        s.path = path
        s.doc = fitz.open(path)
        s.pg = s.doc[0]
        s.W, s.H = s.pg.rect.width, s.pg.rect.height
        s.spans = []
        for b in s.pg.get_text("dict")["blocks"]:
            if b["type"] != 0:
                continue
            for l in b["lines"]:
                for sp in l["spans"]:
                    x0, y0, x1, y1 = sp["bbox"]
                    s.spans.append(dict(t=sp["text"].replace("−", "-"), size=sp["size"], bbox=(x0, y0, x1, y1), cx=(x0 + x1) / 2, cy=(y0 + y1) / 2,
                                        rot=abs(l["dir"][1]) > 0.5, bold="Bold" in sp["font"]))
        s.dr = s.pg.get_drawings()

    def text(s):
        return "\n".join(sp["t"] for sp in s.spans)

    def panels(s, minw=100, minh=100):
        out = []
        for d in s.dr:
            r = d["rect"]
            if d["type"] == "f" and d.get("fill") == (1.0, 1.0, 1.0) and r.width > minw and r.height > minh and r.width < s.W - 5:
                out.append((r.x0, r.y0, r.x1, r.y1))
        return sorted(out, key=lambda p: (round(p[1]), p[0]))


def verts(d):
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            a, b = (it[1].x, it[1].y), (it[2].x, it[2].y)
            if not pts or pts[-1] != a:
                pts.append(a)
            pts.append(b)
    return pts


def is_dashed(d):
    return not re.match(r"\s*\[\s*\]", d.get("dashes") or "[]")


def same(c1, c2, tol=0.004):
    return c1 is not None and c2 is not None and all(abs(a - b) <= tol for a, b in zip(c1, c2))


def inside(pt, p, m=1.0):
    return p[0] - m <= pt[0] <= p[2] + m and p[1] - m <= pt[1] <= p[3] + m


def fit(pos, val):
    """pos = a + b val by least squares; returns (a, b, max residual in pt)."""
    b, a = np.polyfit(val, pos, 1)
    return a, b, float(np.max(np.abs(a + b * np.asarray(val) - np.asarray(pos))))


def panel_of(pt, panels):
    for i, p in enumerate(panels):
        if inside(pt, p, 0.5):
            return i
    return None


def xlabels(P, p, vals=None, below=True):
    """numeric tick labels under panel p (rotated or not): list of (value, centre x)."""
    out = []
    for s in P.spans:
        if p[0] - 3 <= s["cx"] <= p[2] + 3 and p[3] < s["cy"] < p[3] + 90 and re.fullmatch(r"-?\d+(\.\d+)?", s["t"]):
            out.append((float(s["t"]), s["cx"]))
    return sorted(out, key=lambda t: t[1])


def ylabels(P, p):
    out = []
    for s in P.spans:
        if not s["rot"] and s["cx"] < p[0] and p[0] - 70 < s["cx"] and p[1] - 12 <= s["cy"] <= p[3] + 12 and re.fullmatch(r"-?\d+(\.\d+)?", s["t"]):
            out.append((float(s["t"]), s["cy"]))
    return sorted(out, key=lambda t: t[1])


def title_above(P, p, k=1):
    """text spans (bold) above the panel p, within its x range, nearest first (top to bottom order returned)."""
    c = [s for s in P.spans if s["bold"] and p[0] - 40 <= s["cx"] <= p[2] + 40 and p[1] - 90 < s["cy"] < p[1]]
    return sorted(c, key=lambda s: s["cy"])


def legend_colour(P, name, startswith=False, maxdx=60):
    """stroke colour of the legend sample drawn just left of a legend label (the page itself gives the palette); the nearest sample to the left wins."""
    cands = [s for s in P.spans if (s["t"].startswith(name) if startswith else s["t"] == name)]
    best = None
    for s in cands:
        for d in P.dr:
            if d["type"] != "s" or d.get("color") is None or len(d["items"]) != 1 or d["items"][0][0] != "l":
                continue
            r = d["rect"]
            gap = s["bbox"][0] - r.x1
            dy = abs((r.y0 + r.y1) / 2 - s["cy"])
            if r.width < 8 or r.width > maxdx or dy > 12 or not (-2 <= gap <= maxdx):
                continue
            score = gap + 3 * dy
            if best is None or score < best[0]:
                best = (score, d["color"], d)
    return (best[1], best[2]) if best else (None, None)


def data_strokes(P, panel, colour, minw=1.5):
    out = []
    for d in P.dr:
        if d["type"] == "s" and same(d.get("color"), colour) and (d.get("width") or 0) >= minw and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), panel):
            out.append(d)
    return out


# ------------------------------------------------------------------------------------------------------------ D91a, D91f
def check_files(dirp, tag=""):
    ok_all = True
    pages = {}
    for k in range(1, 7):
        p = os.path.join(dirp, f"p2-indicators-fig-F{k}.pdf")
        if not os.path.isfile(p):
            emit("D91a", f"F{k} file", "missing", "exists", "FAIL")
            ok_all = False
            continue
        b = open(p, "rb").read()
        d = fitz.open(p)
        sig = b[:5] == b"%PDF-"
        okf = sig and 1000 < len(b) <= 1024 * 1024 and len(d) == 1
        ok_all &= okf
        emit("D91a", f"p2-indicators-fig-F{k}.pdf: signature, size (bytes), pages", f"{sig}, {len(b)}, {len(d)}", "%PDF-, 1000 < size <= 1048576, 1 page", pf(okf))
    extra = sorted(os.listdir(dirp))
    return ok_all


def check_layout(dirp):
    for k in range(1, 7):
        P = Page(os.path.join(dirp, f"p2-indicators-fig-F{k}.pdf"))
        mn = min(s["size"] for s in P.spans) * 504 / P.W
        okf = P.W <= 1000 and P.H <= 1.3 * P.W and mn >= 8 and len(P.doc) == 1
        emit("D91f", f"F{k}: width, height (pt), height/width, smallest span size scaled to 504 pt, spans", f"{P.W:.0f}, {P.H:.0f}, {P.H / P.W:.3f}, {mn:.3f}, {len(P.spans)}", "<= 1000, <= 1.3 x width, >= 8", pf(okf))


# ------------------------------------------------------------------------------------------------------------ D91g text content
def count_t(P, pred):
    return sum(1 for s in P.spans if pred(s["t"]))


def joined(P):
    return " ".join(s["t"] for s in P.spans)


def check_text(dirp, B, prm):
    bud = prm["budgets"]
    bP, Nc = prm["b_primary"], prm["N_c"]
    R = {k: as_dicts(B[k]) for k in B if k.startswith("F")}
    # F1
    P = Page(os.path.join(dirp, "p2-indicators-fig-F1.pdf"))
    need = [(n, 1) for n in IND[:7]] + [("RLI*", 1), ("A_XY", 1), ("kappa_XY", 1)]
    miss = [n for n, _ in need if count_t(P, lambda t, n=n: t == n) < 1]
    rli_n = {}
    for r in R["F1"]:
        if r["Y"] == "RLI":
            rli_n.setdefault(r["set"], set()).add(r["N"])
    ttl = {}
    ok_t = True
    detail = []
    for st in ("W2", "S2"):
        for nm, b in (("b_P", bP), ("N_c", Nc)):
            Ns = [r["N"] for r in R["F1"] if r["set"] == st and r["budget"] == b and r["Y"] != "RLI" and r["X"] != "RLI"]
            Ns += [] if not Ns else []
            allN = [r["N"] for r in R["F1"] if r["set"] == st and r["budget"] == b and "RLI" not in (r["X"], r["Y"])]
            lo, hi = min(allN), max(allN)
            t1 = f"{st}, {nm} = {b}"
            t2 = f"N = {lo}" if lo == hi else f"N = {lo} to {hi}"
            c1, c2 = count_t(P, lambda t, t1=t1: t == t1), 0
            # the second title line is one per panel; it may repeat across panels, so count the pair on neighbouring lines
            tl = [s for s in P.spans if s["t"] == t1]
            c2 = sum(1 for s in tl for u in P.spans if u["t"] == t2 and abs(u["cx"] - s["cx"]) < 40 and 0 < u["cy"] - s["cy"] < 25)
            ok_t &= c1 == 1 and c2 == 1
            detail.append(f"{t1}/{t2}:{c1},{c2}")
    note = [s["t"] for s in P.spans if s["t"].startswith("*RLI")]
    note_ok = len(note) == 1 and all(str(n) in note[0] for v in rli_n.values() for n in v) and all(len(v) == 1 for v in rli_n.values())
    okf = not miss and ok_t and note_ok and "12800" in P.text() and "102400" in P.text()
    emit("D91g", "F1 text: seven names + RLI* + two colour bar names; titles once each (set, budget and name, N range); note under the figure with the N of the RLI pairs",
         f"missing {miss}; titles {detail}; note {note}; RLI N per set {({k: sorted(v) for k, v in rli_n.items()})}", "none missing, titles once, note holds N", pf(okf))
    # F2
    P = Page(os.path.join(dirp, "p2-indicators-fig-F2.pdf"))
    names_ok = all(count_t(P, lambda t, n=n: t == n) >= 1 for n in IND)
    sets = {}
    for r in R["F2"]:
        if r["kind"] == "outcome":
            sets.setdefault(r["set"], {}).setdefault(r["indicator"], 0)
            sets[r["set"]][r["indicator"]] += r["count"]
    detail, ok_t = [], True
    for st in ("S2", "W2", "T_s"):
        Nmax = max(sets[st].values())
        t = f"{st} (N = {Nmax})"
        c = count_t(P, lambda x, t=t: x == t)
        ok_t &= c == 1
        for ind, n in sets[st].items():
            if n != Nmax:
                lab = f"{ind} (N = {n})"
                c2 = count_t(P, lambda x, lab=lab: x == lab)
                ok_t &= c2 == 1
                detail.append(f"{lab}:{c2}")
        detail.append(f"{t}:{c}")
    bud_ok = all(count_t(P, lambda t, b=b: t == str(b)) >= 3 for b in bud)
    out_ok = all(count_t(P, lambda t, w=w: t == w) >= 1 for w in ("correct", "wrong", "undetected"))
    okf = names_ok and ok_t and bud_ok and out_ok
    emit("D91g", "F2 text: eight names, S2 W2 T_s titles with N once each, N labels of the indicators whose N differs, ten budget values (>= 3 each), three outcome names",
         f"names {names_ok}; {detail}; budgets {bud_ok}; outcomes {out_ok}", "all present, titles once", pf(okf))
    # F3
    P = Page(os.path.join(dirp, "p2-indicators-fig-F3.pdf"))
    dec = {}
    for r in R["F3"]:
        dec.setdefault(r["set"], {}).setdefault(r["indicator"], []).append(r["decided"])
    detail, ok_t = [], True
    for st, nm in (("W2", "W2"), ("S2", "S2"), ("T_s", "T_s")):
        mx = max(max(v) for v in dec[st].values())
        t = f"{nm}: up to {mx} orbits decided"
        c = count_t(P, lambda x, t=t: x == t)
        ok_t &= c == 1
        detail.append(f"{t}:{c}")
    mx = max(max(v) for v in dec["regular"].values())
    jt = joined(P)
    reg_ok = "regular set: share of chaotic" in jt and "false-chaotic rate" in jt and f"up to {mx} orbits decided" in jt
    reg_title_count = count_t(P, lambda x: x.startswith("regular set:"))
    rli_leg = [s["t"] for s in P.spans if s["t"].startswith("RLI (")]
    rli_cnt = {st: max(dec[st]["RLI"]) for st in dec}
    rli_ok = len(rli_leg) >= 1 and all(f"{st} {n}" in jt for st, n in rli_cnt.items())
    names_ok = all(count_t(P, lambda t, n=n: t == n or t.startswith(n + " (")) >= 1 for n in IND)
    bud_ok = all(count_t(P, lambda t, b=b: t == str(b)) >= 4 for b in bud)
    bn_ok = f"b_P = {bP}" in jt and f"N_c = {Nc}" in jt
    okf = ok_t and reg_ok and reg_title_count == 1 and rli_ok and names_ok and bud_ok and bn_ok
    emit("D91g", "F3 text: eight names, titles W2 S2 T_s with 'up to N orbits decided' once each, regular set title (false-chaotic rate), RLI own decided counts, budgets (>= 4 each), b_P and N_c",
         f"{detail}; regular {reg_ok}/{reg_title_count}; RLI legend {rli_leg} vs {rli_cnt}; names {names_ok}; budgets {bud_ok}; b_P/N_c {bn_ok}", "all present", pf(okf))
    # F4
    P = Page(os.path.join(dirp, "p2-indicators-fig-F4.pdf"))
    rows = R["F4"]
    Ktexts = []
    for r in rows:
        if r["K"] not in Ktexts:
            Ktexts.append(r["K"])
    # the block text of K is the number as written in the tables file; read the text of the digits again from the file
    names_ok = all(count_t(P, lambda t, n=n: t == n) >= 1 for n in IND + ["STR"])
    detail, ok_t = [], True
    kt_raw = B.get("_Ktext", [])
    for K in kt_raw:
        pr = [r for r in rows if r["source"] == "pooled" and r["K"] == float(K)]
        nC = [r["n_C"] for r in pr]
        nR = [r["n_R"] for r in pr]
        f = lambda a: f"{min(a)}" if min(a) == max(a) else f"{min(a)} to {max(a)}"
        tt = [f"K = {K}", f"n_C {f(nC)}", f"n_R {f(nR)}"]
        cs = [count_t(P, lambda x, t=t: x == t) for t in tt]
        ok_t &= all(c == 1 for c in cs)
        detail.append("/".join(tt) + ":" + "".join(map(str, cs)))
    msg = "AUC not computed: a class has fewer orbits than the minimum set size"
    seq = [sp["t"] for sp in P.spans]
    parts = msg.split(" ")
    want = ["AUC not computed:", "a class has fewer", "orbits than the", "minimum set size"]
    nmsg = sum(1 for i in range(len(seq) - 3) if seq[i:i + 4] == want and " ".join(seq[i:i + 4]) == msg)  # the four lines of the text, one after the other on the page
    chance = count_t(P, lambda t: t.startswith("chance level"))
    okf = names_ok and ok_t and nmsg == 2 and len(kt_raw) == 10 and chance == 1
    emit("D91g", "F4 text: nine names (STR included), ten K titles with n_C and n_R ranges once each, the E42 text for K = 0.0 and 10.0 (blocks of text), chance level in the legend",
         f"names {names_ok}; E42 text blocks {nmsg}; K values {kt_raw}; titles {detail}; chance label {chance}", "all present, E42 text 2 times", pf(okf))
    # F5
    P = Page(os.path.join(dirp, "p2-indicators-fig-F5.pdf"))
    detail = {n: count_t(P, lambda t, n=n: t == n) for n in IND + ["rho"]}
    rows_t = {f"K index {i}": count_t(P, lambda t, i=i: t == f"K index {i}") for i in range(1, 9)}
    g1 = count_t(P, lambda t: t == f"verdict of each indicator at budget b_P = {bP}")
    g2 = count_t(P, lambda t: t == f"verdict of each indicator at budget N_c = {Nc}")
    okf = all(v == 2 for v in detail.values()) and all(v == 2 for v in rows_t.values()) and g1 == 1 and g2 == 1
    emit("D91g", "F5 text: each of the eight names and 'rho' twice (two groups), the eight row titles twice, one group title per budget",
         f"names {detail}; rows {rows_t}; group titles {g1},{g2}", "2 each, 1 title per budget", pf(okf))
    # F6
    P = Page(os.path.join(dirp, "p2-indicators-fig-F6.pdf"))
    cnt = {}
    for r in R["F6"]:
        cnt.setdefault(r["class"], set()).add(r["count"])
    want = {c: sorted(v) for c, v in cnt.items()}
    ok = {c: count_t(P, lambda t, c=c, n=v[0]: t == f"{c} ({n} orbits)") == 1 for c, v in want.items()}
    okf = all(ok.values()) and "T_s is a subset of the class C" in joined(P)
    emit("D91g", "F6 text: class names with their orbit counts once each (counts of the block), the subset note", f"{want}; {ok}", "C 1815, R 1159, T_s 30", pf(okf))


# ------------------------------------------------------------------------------------------------------------ D91h, D91i readers
def segs(paths):
    """all line segments of the paths (cairo drops collinear intermediate vertices, so points are tested against segments)."""
    out = []
    for d in paths:
        for it in d["items"]:
            if it[0] == "l":
                out.append((it[1].x, it[1].y, it[2].x, it[2].y))
    return np.array(out, dtype=float).reshape(-1, 4)


def seg_dist(S, px, py, W, H):
    """smallest distance from a point to the segments, each coordinate scaled by the axis length (so 0.005 = 0.5 percent of the axis span)."""
    if len(S) == 0:
        return None
    x1, y1, x2, y2 = S[:, 0] / W, S[:, 1] / H, S[:, 2] / W, S[:, 3] / H
    x, y = px / W, py / H
    dx, dy = x2 - x1, y2 - y1
    L = dx * dx + dy * dy
    t = np.where(L > 0, np.clip(((x - x1) * dx + (y - y1) * dy) / np.where(L > 0, L, 1), 0, 1), 0)
    return float(np.min(np.hypot(x - (x1 + t * dx), y - (y1 + t * dy))))


def points_on(S, pts, W, H):
    """pts: expected (x pt, y pt). Returns (found, total, max scaled distance of the found ones, misses)."""
    found, maxe, miss = 0, 0.0, []
    for (xp, yp) in pts:
        e = seg_dist(S, xp, yp, W, H)
        if e is not None and e <= TOL:
            found += 1
            maxe = max(maxe, e)
        else:
            miss.append((round(xp, 1), round(yp, 1), e))
    return found, len(pts), maxe, miss


def small_fills(P, panel, colour, maxsz=20):
    return [d for d in P.dr if d["type"] == "f" and same(d.get("fill"), colour) and d["rect"].width < maxsz and d["rect"].height < maxsz
            and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), panel)]


def top_title(P, p, pan):
    """title spans of the panel p (or of the panel above it in the same column)."""
    col = [q for q in pan if abs(q[0] - p[0]) < 2]
    top = min(col, key=lambda q: q[1])
    return title_above(P, top)


def f1_read(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F1.pdf"))
    pan = P.panels()
    info = P.pg.get_image_info(xrefs=True)
    bars = []
    for im in info:
        x0, y0, x1, y1 = im["bbox"]
        pix = fitz.Pixmap(P.doc, im["xref"])
        labs = [(float(s["t"]), s["cy"]) for s in P.spans if not s["rot"] and x1 < s["cx"] < x1 + 60 and y0 - 12 <= s["cy"] <= y1 + 12 and re.fullmatch(r"-?\d+(\.\d+)?", s["t"])]
        v = np.array([a for a, _ in labs])
        yy = np.array([b for _, b in labs])
        b_, a_ = np.polyfit(yy, v, 1)  # value = a_ + b_ y
        bars.append(dict(y0=y0, y1=y1, x=(x0 + x1) / 2, a_=a_, b_=b_, vmin=float(v.min()), vmax=float(v.max()), n=len(labs)))
    S = 4  # render at 4 px per pt and sample the centre pixel of each cell
    pm = P.pg.get_pixmap(matrix=fitz.Matrix(S, S), alpha=False)
    img = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width, 3)
    for bar in bars:  # the colour bar as drawn: the render column at its centre, row by row, with the value of each row from the tick labels
        rows = np.arange(int((bar["y0"] + 0.5) * S), int((bar["y1"] - 0.5) * S))
        bar["arr"] = img[rows, int(bar["x"] * S)].astype(float)
        bar["vals"] = bar["a_"] + bar["b_"] * ((rows + 0.5) / S)
    cells = [d["rect"] for d in P.dr if d["type"] == "f" and d.get("fill") not in (None, (1.0, 1.0, 1.0)) and 10 < d["rect"].width < 40 and 10 < d["rect"].height < 40]
    res = {"A": [], "kappa": []}
    per_panel = {}
    unmatched = 0
    for r in cells:
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        i = panel_of((cx, cy), pan)
        if i is None:
            unmatched += 1
            continue
        p = pan[i]
        pw, ph = (p[2] - p[0]) / 8, (p[3] - p[1]) / 8
        col, row = int((cx - p[0]) // pw), int((cy - p[1]) // ph)
        cl = [s for s in P.spans if s["rot"] and p[0] <= s["cx"] <= p[2] and p[3] < s["cy"] < p[3] + 80]
        cname = min(cl, key=lambda s: abs(s["cx"] - cx))["t"].rstrip("*")
        first = min(q for q in pan if abs(q[1] - p[1]) < 2)
        rl = [s for s in P.spans if not s["rot"] and s["cx"] < first[0] and first[1] - 3 <= s["cy"] <= first[3] + 3]
        rname = min(rl, key=lambda s: abs(s["cy"] - cy))["t"].rstrip("*")
        tt = top_title(P, p, pan)
        m = re.match(r"(\w+), (b_P|N_c) = (\d+)", tt[0]["t"])
        st, bud = m.group(1), int(m.group(3))
        bar = max(bars, key=lambda b: min(b["y1"], p[3]) - max(b["y0"], p[1]))
        c = img[int(cy * S), int(cx * S)].astype(float)
        a, b = bar["arr"][:-1], bar["arr"][1:]
        ab = b - a
        t = np.einsum("ij,ij->i", c - a, ab) / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-9)
        t = np.clip(t, -1.0, 2.0)  # the colours of the end values lie half a sample beyond the first and last bar row
        d2 = np.sum((c - (a + t[:, None] * ab)) ** 2, axis=1)
        s = int(np.argmin(d2))
        val = bar["vals"][s] + t[s] * (bar["vals"][s + 1] - bar["vals"][s])
        val = min(max(val, bar["vmin"]), bar["vmax"])
        key = "A" if bar["vmin"] >= 0 else "kappa"
        res[key].append((st, bud, rname, cname, col, row, val, bar["vmax"] - bar["vmin"], math.sqrt(d2[s])))
        per_panel.setdefault(i, [0, st, bud])[0] += 1
    return P, pan, res, per_panel, unmatched, bars


def check_f1(dirp, B, prm):
    P, pan, res, per_panel, unmatched, bars = f1_read(dirp, B, prm)
    exp = {(r["set"], r["budget"], r["X"], r["Y"]): r for r in as_dicts(B["F1"])}
    err = {"A": [], "kappa": []}
    bad, tri = [], 0
    for key in ("A", "kappa"):
        for (st, bud, X, Y, col, row, val, span, dist) in res[key]:
            e = exp.get((st, bud, X, Y))
            if e is None:
                bad.append((key, st, bud, X, Y))
                continue
            tri += col > row and IND.index(Y) == col and IND.index(X) == row
            err[key].append(abs(val - (e["A"] if key == "A" else e["kappa"])) / span)
    ncell = len(res["A"]) + len(res["kappa"])
    okn = {k: sum(1 for x in v if x <= TOL) for k, v in err.items()}
    emit("D91h", "F1 cell colours decoded with the colour bars of the page (render 4 px/pt, centre of each cell): cells within 0.5 percent of the colour bar span, A_XY and kappa_XY, max error as fraction of the span",
         f"{ncell} cells (rects outside a matrix {unmatched}); A {okn['A']}/{len(err['A'])} max {max(err['A'], default=0):.4f}; kappa {okn['kappa']}/{len(err['kappa'])} max {max(err['kappa'], default=0):.4f}; pairs not in the block {len(bad)}",
         "112 + 112 cells within 0.005", pf(ncell == 224 and okn["A"] == 112 and okn["kappa"] == 112 and not bad))
    emit("D91i", "F1 drawing counts: matrix panels; filled cells per matrix; cells above the diagonal that carry the pair (X, Y) of the block",
         f"{len(pan)} panels; cells per panel {sorted(v[0] for v in per_panel.values())}; consistent {tri}", "8 matrices, 28 cells each, 224 consistent", pf(len(pan) == 8 and len(per_panel) == 8 and all(v[0] == 28 for v in per_panel.values()) and tri == 224))


def check_f2(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F2.pdf"))
    R = as_dicts(B["F2"])
    bud = prm["budgets"]
    pan = P.panels()
    top, bot = pan[:3], pan[3:]
    emit("D91i", "F2 panels (white areas): top row and bottom row", f"{len(top)} + {len(bot)}", "3 + 3 = six panels", pf(len(pan) == 6))
    pal = {n: legend_colour(P, n)[0] for n in IND}
    yl = ylabels(P, top[0])
    ya, yb, yres = fit([c for _, c in yl], [v for v, _ in yl])
    ycal = lambda v: ya + yb * v
    tot = found = 0
    maxe = xres_all = 0.0
    miss_all, nmark = [], []
    for p in top:
        st = title_above(P, p)[0]["t"].split(" ")[0]
        xl = xlabels(P, p)
        xa, xb, xres = fit([c for _, c in xl], [math.log10(v) for v, _ in xl])
        xres_all = max(xres_all, xres)
        xp = [xa + xb * math.log10(b) for b in bud]
        for ind in IND:
            rows = [r for r in R if r["set"] == st and r["indicator"] == ind]
            N = sum(r["count"] for r in rows if r["kind"] == "outcome")
            td = {r["value"]: r["count"] for r in rows if r["kind"] == "t_det"}
            cum, ys = 0, []
            for b in bud:
                cum += td.get(b, 0)
                ys.append(ycal(cum / N))
            # step line: flat at the share of budget k between budget k and k + 1 (midpoints), and the last point
            pts = [((xp[k] + xp[k + 1]) / 2, ys[k]) for k in range(9)] + [(xp[9], ys[9])]
            paths = [d for d in data_strokes(P, p, pal[ind])]
            f, t, e, miss = points_on(segs(paths), pts, p[2] - p[0], p[3] - p[1])
            found += f
            tot += t
            maxe = max(maxe, e)
            miss_all += [(st, ind) + m for m in miss]
            nmark.append(len(small_fills(P, p, pal[ind])))
    emit("D91h", "F2 top row: running share of t_det (step lines) re-read with the tick labels of the page: 9 flat-segment midpoints and the last point per line lie on the drawn path within 0.5 percent of the panel size (3 sets x 8 indicators x 10), max error fraction",
         f"{found}/{tot}, max {maxe:.4f}; label fit residual {xres_all:.2f} pt (x), {yres:.2f} pt (y); misses {miss_all[:3]}", "240 within 0.005 (at least 20)", pf(found == tot and tot == 240))
    emit("D91i", "F2 top row: lines (series paths) and markers per line (one per budget)", f"{len(nmark)} lines; markers per line {sorted(set(nmark))}", "24 lines, 10 markers each", pf(len(nmark) == 24 and set(nmark) == {10}))
    # bottom row: stacked bars; the three colours come from the legend swatches on the page
    sw = {}
    for nm in ("correct", "wrong", "undetected"):
        s = [s for s in P.spans if s["t"] == nm][0]
        sw[nm] = [d["fill"] for d in P.dr if d["type"] == "f" and d.get("fill") and d["rect"].width < 30 and abs((d["rect"].y0 + d["rect"].y1) / 2 - s["cy"]) < 8 and 0 < s["bbox"][0] - d["rect"].x1 < 12][0]
    lum = [sum(sw[n]) for n in ("correct", "wrong", "undetected")]
    totb = okb = nrect = 0
    maxb = 0.0
    for p in bot:
        ti = min(range(3), key=lambda i: abs((top[i][0] + top[i][2]) / 2 - (p[0] + p[2]) / 2))
        st = title_above(P, top[ti])[0]["t"].split(" ")[0]
        xl = [(float(s["t"]), s["cx"]) for s in P.spans if not s["rot"] and p[0] - 3 <= s["cx"] <= p[2] + 3 and p[3] < s["cy"] < p[3] + 25 and re.fullmatch(r"\d(\.\d+)?", s["t"])]
        xa, xb, _ = fit([c for _, c in xl], [v for v, _ in xl])
        labs = [s for s in P.spans if not s["rot"] and p[0] - 120 < s["cx"] < p[0] and p[1] <= s["cy"] <= p[3] and not re.fullmatch(r"-?\d+(\.\d+)?", s["t"])]
        rects = {n: [d["rect"] for d in P.dr if d["type"] == "f" and same(d.get("fill"), sw[n], 0.004) and d["items"][0][0] == "re" and d["rect"].height < 60 and d["rect"].width > 0
                     and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), p)] for n in sw}
        nrect += sum(len(v) for v in rects.values())
        for ind in IND:
            cnt = {r["value"]: r["count"] for r in R if r["set"] == st and r["indicator"] == ind and r["kind"] == "outcome"}
            N = sum(cnt.values())
            lab = [s for s in labs if s["t"].split(" ")[0] == ind][0]
            cum = 0.0
            for oc in ("correct", "wrong", "undetected"):
                share = cnt.get(oc, 0) / N
                rr = [r for r in rects[oc] if abs((r.y0 + r.y1) / 2 - lab["cy"]) < 6]
                if rr:
                    e = max(abs((rr[0].x0 - xa) / xb - cum), abs((rr[0].x1 - xa) / xb - (cum + share)))
                else:
                    e = 0.0 if share == 0 else 1.0
                totb += 1
                okb += e <= TOL
                maxb = max(maxb, e)
                cum += share
    emit("D91h", "F2 bottom row: stacked bars, both edges of every segment (correct, wrong, undetected; colours from the legend swatches) against the cumulative shares of the outcome rows, 3 sets x 8 indicators x 3, max error (axis 0 to 1)",
         f"{okb}/{totb}, max {maxb:.4f}; bar rectangles read {nrect}; swatch luminance order dark to light {lum[0] < lum[1] < lum[2]}", "72 within 0.005, dark to light", pf(okb == totb and totb == 72 and lum[0] < lum[1] < lum[2]))


def check_f3(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F3.pdf"))
    R = as_dicts(B["F3"])
    bud = prm["budgets"]
    pan = P.panels()
    pal = {n: legend_colour(P, n, startswith=(n == "RLI"))[0] for n in IND}
    found = tot = 0
    maxe = 0.0
    nmark, miss_all, rules_ok = [], [], []
    for p in pan:
        st = title_above(P, p)[0]["t"].split(":")[0]
        st = "regular" if st == "regular set" else st  # the title says "regular set", the block says "regular"
        xl = xlabels(P, p)
        xa, xb, xres = fit([c for _, c in xl], [bud.index(int(v)) for v, _ in xl])
        xcal = lambda k, xa=xa, xb=xb: xa + xb * k
        first = min(q for q in pan if abs(q[1] - p[1]) < 2)
        yl = ylabels(P, first)
        ya, yb, yres = fit([c for _, c in yl], [v for v, _ in yl])
        for ind in IND:
            pts = [(xcal(r["k"]), ya + yb * r["fraction"]) for r in R if r["set"] == st and r["indicator"] == ind]
            f, t, e, miss = points_on(segs(data_strokes(P, p, pal[ind])), pts, p[2] - p[0], p[3] - p[1])
            found += f
            tot += t
            maxe = max(maxe, e)
            miss_all += [(st, ind) + m for m in miss]
            nmark.append(len(small_fills(P, p, pal[ind])))
        dv = [d for d in P.dr if d["type"] == "s" and is_dashed(d) and abs(d["rect"].x1 - d["rect"].x0) < 0.1 and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), p)]
        xs = sorted(d["rect"].x0 for d in dv)
        want = [xcal(bud.index(prm["b_primary"])), xcal(bud.index(prm["N_c"]))]
        rules_ok.append(len(xs) == 2 and all(abs(a - b) / (p[2] - p[0]) <= TOL for a, b in zip(xs, want)))
    emit("D91h", "F3 lines: fraction against budget index, 4 panels x 8 indicators x 10 indices re-read with the tick labels of the page: points on the drawn path within 0.5 percent of the panel size, max error fraction",
         f"{found}/{tot}, max {maxe:.4f}; misses {miss_all[:3]}", "320 within 0.005 (at least 20)", pf(found == tot and tot == 320))
    emit("D91h", "F3 dashed rules at b_P (k = 6) and N_c (k = 9): two dashed verticals per panel at the tick positions of 12800 and 102400", f"{rules_ok}", "4 panels True", pf(all(rules_ok) and len(rules_ok) == 4))
    emit("D91i", "F3 drawing counts: panels; lines (series paths) and markers per line (ten points per line)", f"{len(pan)} panels; {len(nmark)} lines; markers per line {sorted(set(nmark))}", "4 panels, 32 lines, 10 markers each", pf(len(pan) == 4 and len(nmark) == 32 and set(nmark) == {10}))


def check_f4(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F4.pdf"))
    R = [r for r in as_dicts(B["F4"]) if r["source"] == "pooled"]
    bud = prm["budgets"]
    pan = P.panels()
    pal = {n: legend_colour(P, n)[0] for n in IND + ["STR"]}
    chance_c, _ = legend_colour(P, "chance level", startswith=True)
    found = tot = gap_ok = gap_tot = 0
    maxe = 0.0
    nlines, nmark, miss_all, chance = [], 0, [], []
    for p in pan:
        K = float(title_above(P, p)[0]["t"].replace("K = ", ""))
        xl = xlabels(P, p)
        xa, xb, xres = fit([c for _, c in xl], [bud.index(int(v)) for v, _ in xl])
        xcal = lambda k, xa=xa, xb=xb: xa + xb * k
        first = min(q for q in pan if abs(q[1] - p[1]) < 2)
        yl = ylabels(P, first)
        ya, yb, yres = fit([c for _, c in yl], [v for v, _ in yl])
        ycal = lambda v, ya=ya, yb=yb: ya + yb * v
        hz = [d for d in P.dr if d["type"] == "s" and same(d.get("color"), chance_c) and is_dashed(d) and abs(d["rect"].y1 - d["rect"].y0) < 0.1 and inside(((d["rect"].x0 + d["rect"].x1) / 2, d["rect"].y0), p)]
        chance.append((K, len(hz), all(abs(d["rect"].y0 - ycal(0.5)) / (p[3] - p[1]) <= TOL for d in hz) if hz else None))
        nl = 0
        for ind in IND + ["STR"]:
            rows = [r for r in R if r["indicator"] == ind and r["K"] == K]
            paths = [d for d in data_strokes(P, p, pal[ind]) if not (ind == "STR" and same(d.get("color"), chance_c))]
            S = segs(paths)
            nl += len(paths) > 0
            pts = [(xcal(r["k"]), ycal(r["AUC"])) for r in rows if r["AUC"] is not None]
            f, t, e, miss = points_on(S, pts, p[2] - p[0], p[3] - p[1])
            found += f
            tot += t
            maxe = max(maxe, e)
            miss_all += [(K, ind) + m for m in miss]
            mk = small_fills(P, p, pal[ind])
            nmark += len(mk)
            # a gap is empty: no marker at the budget and no segment of the series spanning it
            for r in rows:
                if r["AUC"] is None:
                    gap_tot += 1
                    x = xcal(r["k"])
                    hit = any(abs((m["rect"].x0 + m["rect"].x1) / 2 - x) <= 2.0 for m in mk) or any((min(a, c) < x - 2.0 and max(a, c) > x + 2.0) or abs(a - x) <= 2.0 or abs(c - x) <= 2.0 for a, _, c, _ in S)
                    gap_ok += not hit
        nlines.append((K, nl))
    emit("D91h", "F4 lines: pooled AUC against budget index, 10 K x 9 series, every value of the block re-read with the tick labels of the page: points on the drawn path within 0.5 percent of the panel size, max error fraction",
         f"{found}/{tot}, max {maxe:.4f}; misses {miss_all[:3]}", "700 within 0.005 (at least 20)", pf(found == tot and tot == 700))
    emit("D91h", "F4 chance-level rule (K, dashed lines found, at AUC = 0.5)", f"{chance}", "dashed line at 0.5 in every panel that has lines", pf(all(c[2] for c in chance if c[1] > 0) and sum(1 for c in chance if c[1] > 0) >= 8))
    emit("D91i", "F4 chance-level line in the panels where every AUC is nothing (K = 0.0 and 10.0): found / not drawn", f"{[c for c in chance if c[1] == 0]}", "E42: such a panel shows no lines; reported only", "INFO")
    strp = [d for q in pan for d in data_strokes(P, q, pal["STR"]) if not same(d.get("color"), chance_c)]
    emit("D91h", "F4 STR drawn in black and dashed (E38 d): STR paths found; black; dashed", f"{len(strp)} paths; black {all(same(d.get('color'), (0.0, 0.0, 0.0)) for d in strp)}; dashed {all(is_dashed(d) for d in strp)}", "8 panels with a STR path; black; dashed", pf(len(strp) >= 8 and all(same(d.get("color"), (0.0, 0.0, 0.0)) and is_dashed(d) for d in strp)))
    emit("D91h", "F4 gaps: nothing drawn at the budget of a 'nothing' pooled row (no marker within 2 pt, no segment of the series through it), all series, all K", f"{gap_ok}/{gap_tot}", "200 of 200 gaps empty", pf(gap_ok == gap_tot and gap_tot == 200))
    expn = {}
    for r in R:
        expn.setdefault(r["K"], set())
        if r["AUC"] is not None:
            expn[r["K"]].add(r["indicator"])
    expn = {K: len(v) for K, v in expn.items()}
    ok_series = all(expn[K] == n for K, n in nlines)
    nine = sum(1 for _, n in nlines if n == 9)
    emit("D91i", "F4 drawing counts: panels; series drawn per panel (K, series) against the series that have a value in the block; markers of the series colours inside the panels",
         f"{len(pan)} panels; drawn {nlines}; block {sorted(expn.items())}; panels with all nine series {nine} (the spec text says nine lines each; K = 0.0 and 10.0 have none by E42, and RLI has no AUC at two K); {nmark} markers",
         "10 panels; series equal the block; 700 markers", pf(len(pan) == 10 and ok_series and nmark == 700))
    nn = sum(1 for r in R if r["AUC"] is None)
    emit("D91i", "F4 gaps of the block: pooled rows that are 'nothing' / pooled rows", f"{nn} of {len(R)}", "200 of 900", pf(nn == 200 and len(R) == 900))


def check_f5(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F5.pdf"))
    R = as_dicts(B["F5"])
    cells = [d["rect"] for d in P.dr if d["type"] == "f" and d.get("fill") not in (None, (1.0, 1.0, 1.0)) and d["rect"].width < 12 and d["rect"].height < 12 and d["rect"].y1 < 980]
    xs = sorted(set(round((r.x0 + r.x1) / 2, 1) for r in cells))
    ys = sorted(set(round((r.y0 + r.y1) / 2, 1) for r in cells))

    def clusters(v):
        d = np.diff(v)
        pitch = np.median(d)
        return np.split(np.array(v), [i + 1 for i, x in enumerate(d) if x > 1.5 * pitch])

    xc, yc = clusters(xs), clusters(ys)
    emit("D91i", "F5 drawing counts from the page: filled cells; map columns x cell columns; map rows x cell rows (144 maps of 24 x 12 cells: 8 x 9 per group, two groups)",
         f"{len(cells)} cells; {len(xc)} x {sorted(set(len(c) for c in xc))}; {len(yc)} x {sorted(set(len(c) for c in yc))} -> {len(xc) * len(yc)} maps",
         "41472 cells; 9 x 24; 16 x 12 -> 144 maps", pf(len(cells) == 41472 and len(xc) == 9 and len(yc) == 16 and set(len(c) for c in xc) == {24} and set(len(c) for c in yc) == {12}))
    S = 8  # render at 8 px per pt, read the centre of every cell, recover the three colours from the page
    pm = P.pg.get_pixmap(matrix=fitz.Matrix(S, S), alpha=False)
    img = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width, 3)
    cen = np.array([[(r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2] for r in cells])
    pix = img[(cen[:, 1] * S).astype(int), (cen[:, 0] * S).astype(int)].astype(int)
    uniq, cnt = np.unique(pix, axis=0, return_counts=True)
    top3 = uniq[np.argsort(-cnt)[:3]]
    sw = {}
    for nm in ("C chaotic", "R regular", "U undecided"):
        s = [s for s in P.spans if s["t"] == nm][0]
        best = None
        for d in P.dr:
            if d.get("fill") not in (None, (1.0, 1.0, 1.0)) and d["rect"].width < 30 and abs((d["rect"].y0 + d["rect"].y1) / 2 - s["cy"]) < 8 and d["rect"].x1 <= s["bbox"][0] + 2 and s["bbox"][0] - d["rect"].x1 < 40:
                best = d["fill"]
        sw[nm[0]] = np.array([round(c * 255) for c in best])
    near = {k: [tuple(int(x) for x in t) for t in top3 if np.max(np.abs(t - v)) <= 3] for k, v in sw.items()}
    emit("D91h", "F5 colours recovered from the page: distinct cell-centre colours of the render (8 px/pt), the three most frequent, the legend swatches (C chaotic, R regular, U undecided) and the clusters within 3 of 255 of a swatch",
         f"{len(uniq)} distinct; top3 {[tuple(int(x) for x in t) for t in top3]}; swatches { {k: tuple(int(x) for x in v) for k, v in sw.items()} }; matched { {k: len(v) for k, v in near.items()} }",
         "three clusters, each within 3/255 of its swatch", pf(all(len(v) == 1 for v in near.values()) and len(uniq) <= 6))
    keys = list(sw)
    cls = np.full(len(cells), -1)
    for j, k in enumerate(keys):
        cls[np.max(np.abs(pix - sw[k]), axis=1) <= 3] = j
    unk = int(np.sum(cls < 0))
    xl_idx = {round(v, 1): (mi, ii) for mi, c in enumerate(xc) for ii, v in enumerate(c)}
    yl_idx = {round(v, 1): (mi, 11 - jj) for mi, c in enumerate(yc) for jj, v in enumerate(c)}  # the top of a map is j = 11
    colnames = []
    for c in xc:
        cx = (c[0] + c[-1]) / 2
        colnames.append(min([s for s in P.spans if not s["rot"] and s["t"] in IND + ["rho"]], key=lambda s: abs(s["cx"] - cx) + 1000 * (s["cy"] > 80))["t"])
    kidx = []
    for c in yc:
        cy = (c[0] + c[-1]) / 2
        kidx.append(int(min([s for s in P.spans if s["t"].startswith("K index")], key=lambda s: abs(s["cy"] - cy))["t"].split()[-1]))
    gt = sorted([s for s in P.spans if s["t"].startswith("verdict of each indicator at budget")], key=lambda s: s["cy"])
    gbud = ["b_P" if "b_P" in g["t"] else "N_c" for g in gt]
    exp = {(r["K index"], r["i"], r["j"]): r for r in R}
    cmap = {"C": 0, "chaotic": 0, "R": 1, "regular": 1, "U": 2, "undecided": 2}
    bad, seen, bycol = 0, set(), {}
    for n, (cx, cy) in enumerate(cen):
        mi, ii = xl_idx[round(cx, 1)]
        mj, jj = yl_idx[round(cy, 1)]
        grp = mj // 8
        col, K = colnames[mi], kidx[mj]
        seen.add((grp, K, col, ii, jj))
        r = exp[(K, ii, jj)]
        e = cmap[r["rho"]] if col == "rho" else cmap[r[f"{col} {gbud[grp]}"]]
        ok = cls[n] == e
        bad += not ok
        bycol.setdefault(col, [0, 0])
        bycol[col][0] += 1
        bycol[col][1] += not ok
    emit("D91h", "F5 colour of every cell (render 8 px/pt, centre, matched to the legend swatches within 3/255) against rho and the verdict of the block row (K index, i, j); group b_P above, N_c below",
         f"{len(cen)} cells read, {bad} different, {unk} not matched to a swatch; distinct (group, K index, column, i, j) {len(seen)}; groups {gbud}; rows K index {kidx[:8]}", "41472 equal, 0 different", pf(bad == 0 and unk == 0 and len(seen) == 41472 and gbud == ["b_P", "N_c"]))


def verts_poly(d):
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            pts += [(it[1].x, it[1].y), (it[2].x, it[2].y)]
    return pts


def check_f6(dirp, B, prm):
    P = Page(os.path.join(dirp, "p2-indicators-fig-F6.pdf"))
    R = as_dicts(B["F6"])
    pan = P.panels()
    p = pan[0]
    xl = [(float(s["t"]), s["cx"]) for s in P.spans if not s["rot"] and p[0] - 3 <= s["cx"] <= p[2] + 3 and p[3] < s["cy"] < p[3] + 25 and re.fullmatch(r"\d+", s["t"])]
    yl = [(float(s["t"]), s["cy"]) for s in P.spans if not s["rot"] and s["cx"] < p[0] and p[1] - 12 <= s["cy"] <= p[3] + 12 and re.fullmatch(r"\d+(\.\d+)?", s["t"])]
    xa, xb, xres = fit([c for _, c in xl], [math.log10(v) for v, _ in xl])
    ya, yb, yres = fit([c for _, c in yl], [math.log10(v) for v, _ in yl])
    xcal = lambda v: xa + xb * math.log10(v)
    ycal = lambda v: ya + yb * math.log10(v)
    W, H = p[2] - p[0], p[3] - p[1]
    tot = found = nlines = nbands = 0
    maxe = 0.0
    miss_all = []
    for cl in ("C", "R", "T_s"):
        col, _ = legend_colour(P, cl + " (", startswith=True)
        rows = [r for r in R if r["class"] == cl]
        paths = [d for d in data_strokes(P, p, col) if len(verts(d)) > 20]
        nlines += len(paths)
        f, t, e, miss = points_on(segs(paths), [(xcal(r["step"]), ycal(r["p50"])) for r in rows], W, H)
        found += f
        tot += t
        maxe = max(maxe, e)
        miss_all += [(cl, "p50") + m for m in miss]
        polys = sorted([d for d in P.dr if d["type"] == "f" and same(d.get("fill"), col) and len(verts_poly(d)) > 20 and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), p)], key=lambda d: d["rect"].height)
        nbands += len(polys)
        for d, names in zip(polys, (("p25", "p75"), ("p5", "p95"))):
            S = segs([d])
            for nm in names:
                f, t, e, miss = points_on(S, [(xcal(r["step"]), ycal(r[nm])) for r in rows], W, H)
                found += f
                tot += t
                maxe = max(maxe, e)
                miss_all += [(cl, nm) + m for m in miss]
    emit("D91h", "F6 median lines and the two bands (p25 to p75, p5 to p95) of the three classes, every step of the block re-read with the tick labels of the page (both axes logarithmic): points on the drawn path or band outline within 0.5 percent of the panel size, max error fraction",
         f"{found}/{tot}, max {maxe:.4f}; label fit residual {xres:.2f}, {yres:.2f} pt; misses {miss_all[:3]}", "489 x 5 within 0.005 (at least 20)", pf(found == tot and tot == 489 * 5))
    mk = {cl: len(small_fills(P, p, legend_colour(P, cl + " (", startswith=True)[0])) for cl in ("C", "R", "T_s")}
    emit("D91i", "F6 markers on the median lines (E38 d: a marker at every plotted point): markers of the class colour inside the panel per class", f"{mk}", "163 each (one per recorded step)", pf(all(v == 163 for v in mk.values())))
    emit("D91i", "F6 drawing counts: panels; median lines (stroked paths of the class colour); band polygons (two per class)", f"{len(pan)} panel; {nlines} lines; {nbands} bands", "1 panel, 3 lines, 6 bands", pf(len(pan) == 1 and nlines == 3 and nbands == 6))


def guard(f, name, *a):
    try:
        f(*a)
    except Exception as e:
        import traceback
        emit(name, "ERROR while reading the PDF", (type(e).__name__ + ": " + str(e) + " | " + traceback.format_exc().splitlines()[-3].strip())[:400], "no error", "FAIL")


def marker_fp(P, panel, colour):
    for d in P.dr:
        if d["type"] == "f" and same(d.get("fill"), colour) and d["rect"].width < 20 and d["rect"].height < 20 and inside(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2), panel):
            return (tuple(it[0] for it in d["items"]), round(d["rect"].width, 1), round(d["rect"].height, 1))
    return None


def check_encoding(dirp):
    """I28 (E38 c, d): the same indicator has the same colour and marker in F2, F3 and F4; C and R have the colours of F5 in F6; the palette avoids the class colours."""
    pg = {k: Page(os.path.join(dirp, f"p2-indicators-fig-F{k}.pdf")) for k in range(1, 7)}
    pal = {k: {n: legend_colour(pg[k], n, startswith=(n == "RLI" and k == 3))[0] for n in IND} for k in (2, 3, 4)}
    same_col = all(same(pal[2][n], pal[3][n], 1 / 255) and same(pal[2][n], pal[4][n], 1 / 255) for n in IND)
    mk = {}
    for k in (2, 3, 4):
        mk[k] = {}
        for n in IND:
            mk[k][n] = next((f for f in (marker_fp(pg[k], q, pal[k][n]) for q in pg[k].panels()) if f), None)
    same_mk = all(mk[2][n] is not None and mk[2][n][0] == mk[3][n][0] == mk[4][n][0] if mk[4][n] and mk[3][n] and mk[2][n] else False for n in IND)
    P5 = pg[5]
    sw = {}
    for nm in ("C chaotic", "R regular", "U undecided"):
        s = [s for s in P5.spans if s["t"] == nm][0]
        sw[nm[0]] = [d["fill"] for d in P5.dr if d.get("fill") not in (None, (1.0, 1.0, 1.0)) and d["rect"].width < 30 and abs((d["rect"].y0 + d["rect"].y1) / 2 - s["cy"]) < 8 and d["rect"].x1 <= s["bbox"][0] + 2 and s["bbox"][0] - d["rect"].x1 < 40][-1]
    c6 = {cl: legend_colour(pg[6], cl + " (", startswith=True)[0] for cl in ("C", "R", "T_s")}
    cls_ok = same(c6["C"], sw["C"], 2 / 255) and same(c6["R"], sw["R"], 2 / 255)
    dist = min(math.dist(pal[2][n], sw[c]) for n in IND for c in ("C", "R"))
    third = min(math.dist(c6["T_s"], x) for x in list(pal[2].values()) + [sw["C"], sw["R"], sw["U"]])
    emit("I28", "one encoding: the colour of each of the eight indicators equal in the legends of F2, F3 and F4 (within 1/255); the marker shape (item types) equal in the three figures; STR black in F4",
         f"colours equal {same_col}; marker shapes equal {same_mk}; shapes {({n: (mk[2][n][0] if mk[2][n] else None) for n in IND})}", "True, True", pf(same_col and same_mk))
    emit("I28", "classes: C and R colours of F6 (legend) equal the C and R swatches of F5 (within 2/255); indicator palette distance from the two class colours (RGB, 0 to 1.73); T_s colour distance from every other colour",
         f"{cls_ok}; min palette distance from C and R {dist:.3f}; T_s distance {third:.3f}; F5 swatches C {tuple(round(x, 3) for x in sw['C'])} R {tuple(round(x, 3) for x in sw['R'])} U {tuple(round(x, 3) for x in sw['U'])}; F6 C {tuple(round(x, 3) for x in c6['C'])} R {tuple(round(x, 3) for x in c6['R'])} T_s {tuple(round(x, 3) for x in c6['T_s'])}",
         "equal; palette distance > 0.15; T_s distance > 0.15", pf(cls_ok and dist > 0.15 and third > 0.15))


def cmd_figs(dirp, tables, params):
    B = load_blocks(tables)
    # the K values as text of the block (the figure titles use them as written)
    txt = open(tables, "rb").read().decode()
    m = re.search(r"## F4 data[^\n]*\n\n((?:\|[^\n]*\n)+)", txt)
    Kt = []
    for l in m.group(1).split("\n")[2:]:
        if l.startswith("|"):
            k = [c.strip() for c in l.strip("|").split("|")][2]
            if k not in Kt:
                Kt.append(k)
    B["_Ktext"] = Kt
    prm = load_params(params)
    ext = sorted(os.listdir(dirp))
    emit("D91a", "files in the output folder", ext, "tables file and p2-indicators-fig-F1..F6.pdf only", pf(ext == sorted(["p2-indicators-tables.md"] + [f"p2-indicators-fig-F{k}.pdf" for k in range(1, 7)])))
    check_files(dirp)
    guard(check_layout, "D91f", dirp)
    guard(check_text, "D91g", dirp, B, prm)
    for name, f in (("F1", check_f1), ("F2", check_f2), ("F3", check_f3), ("F4", check_f4), ("F5", check_f5), ("F6", check_f6)):
        guard(f, "D91h " + name, dirp, B, prm)


# ------------------------------------------------------------------------------------------------------------ pixel comparison, signature, contact sheets, timed run
def render(path, dpi=72):
    d = fitz.open(path)
    pm = d[0].get_pixmap(dpi=dpi, alpha=False)
    return np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width, pm.n), len(d)


def cmp_files(a, b):
    ia, na = render(a)
    ib, nb = render(b)
    if ia.shape != ib.shape:
        return f"shapes {ia.shape} vs {ib.shape}", -1, False
    nd = int(np.sum(np.any(ia != ib, axis=2)))
    same_bytes = open(a, "rb").read() == open(b, "rb").read()
    return f"{ia.shape[1]} x {ia.shape[0]} px, pages {na}/{nb}, raw bytes equal {same_bytes}", nd, same_bytes


def cmd_cmp(da, db, names):
    for n in names:
        info, nd, sb = cmp_files(os.path.join(da, n), os.path.join(db, n))
        print(f"{n}\t{nd}\t{info}")


def cmd_sig(files):
    for f in files:
        if not os.path.isfile(f):
            print(f"{os.path.basename(f)}\tmissing")
            continue
        b = open(f, "rb").read()
        d = fitz.open(f)
        print(f"{os.path.basename(f)}\t{b[:5] == b'%PDF-'}\t{len(b)}\t{len(d)}\t{d[0].rect.width:.0f}x{d[0].rect.height:.0f}")


def cmd_contact(dirp, outp):
    os.makedirs(outp, exist_ok=True)
    for k in range(1, 7):
        d = fitz.open(os.path.join(dirp, f"p2-indicators-fig-F{k}.pdf"))
        d[0].get_pixmap(dpi=72, alpha=False).save(os.path.join(outp, f"analysis-d91-F{k}.png"))
        print(f"analysis-d91-F{k}.png")


def cmd_run(outf, jsonf, cwd, cmd):
    import psutil
    t0 = time.perf_counter()
    with open(outf, "wb") as fo:
        p = psutil.Popen(cmd, cwd=cwd, stdout=fo, stderr=subprocess.STDOUT)
        peak = 0
        while p.poll() is None:
            try:
                peak = max(peak, p.memory_info().peak_wset)
            except Exception:
                pass
            time.sleep(0.05)
    wall = time.perf_counter() - t0
    open(jsonf, "w").write("exit=%d" % p.returncode + chr(10) + "wall=%r" % wall + chr(10) + "peak_wset=%d" % peak + chr(10))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "figs":
        cmd_figs(*a[1:4])
    elif a[0] == "cmp":
        cmd_cmp(a[1], a[2], a[3:])
    elif a[0] == "cmp1":
        print(cmp_files(a[1], a[2])[1])
    elif a[0] == "sig":
        cmd_sig(a[1:])
    elif a[0] == "unzip":
        import zipfile
        zipfile.ZipFile(a[1]).extractall(a[2])
    elif a[0] == "text":
        P = Page(a[1])
        print(chr(10).join(sp["t"] for sp in P.spans))
    elif a[0] == "enc":
        check_encoding(a[1])
    elif a[0] == "contact":
        cmd_contact(a[1], a[2])
    elif a[0] == "run":
        i = a.index("--")
        cmd_run(a[1], a[2], a[3], a[i + 1:])
    else:
        sys.exit("unknown command")
