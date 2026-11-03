# D95k (V): compares the module's HC power at achieved failures (module.jsonl of verify/analysis-d95k.jl) with the independent simulator
# sim_hc of verify/analysis-p3-oracle.py (own generator and streams): PASS when |p_module - p_oracle| <= 3 sqrt(se_module^2 + se_oracle^2).
# Run: python -I verify/analysis-d95k.py <out folder>
import importlib.util, json, math, os, sys
from fractions import Fraction
from statistics import NormalDist
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("p3oracle", os.path.join(here, "analysis-p3-oracle.py"))
O = importlib.util.module_from_spec(spec); spec.loader.exec_module(O)
O.SIM_TEST_SEED[0] = "d95k-independent"
out = sys.argv[1]
rows = [json.loads(l) for l in open(os.path.join(out, "module.jsonl"))]
nbad = 0
for r in rows:
    dI = math.sqrt(2) * NormalDist().inv_cdf(float(r["auc_i_used"])); dP = math.sqrt(2) * NormalDist().inv_cdf(float(r["auc_p"]))
    s = O.sim_hc(r["nf"], r["nn"], dI, dP, r["runs"], Fraction(r["delta_c"]), float(r["z_one_sided"]))
    n = r["runs"]; pm, po = r["successes"] / n, s / n
    se = math.sqrt(pm * (1 - pm) / n + po * (1 - po) / n)
    ok = abs(pm - po) <= 3 * se if se > 0 else pm == po
    nbad += not ok
    print("D95k | nf=%d nn=%d | module %d, oracle %d of %d | diff %.4f, 3 SE %.4f | %s" % (r["nf"], r["nn"], r["successes"], s, n, pm - po, 3 * se, "PASS" if ok else "FAIL"))
print("D95k | summary: %d sizes, %d outside 3 SE | %s" % (len(rows), nbad, "PASS" if nbad == 0 else "FAIL"))
