"""Compare the golden (module before D95) and head outputs of verify/analysis-p3t-regress.sh, ignoring the lines 'module version' and 'git commit' (I44).
python -I verify/analysis-p3t-regress.py <regress dir>"""
import os, re, sys, hashlib
root = sys.argv[1]
ign = re.compile(r"module version|git commit|git_commit|\"version\"|package_version|git_sha|elapsed|wall", re.I)
tot_bad = 0
for dp, _, fs in sorted(os.walk(os.path.join(root, "golden"))):
    for f in sorted(fs):
        if f.endswith(".out"):
            continue
        g = os.path.join(dp, f); h = g.replace(os.sep + "golden" + os.sep, os.sep + "head" + os.sep)
        if not os.path.exists(h):
            print("REGRESS | %s | missing in head | FAIL" % os.path.relpath(g, root)); tot_bad += 1; continue
        a = open(g, encoding="utf-8", errors="replace").read().replace("\r\n", "\n").split("\n")
        b = open(h, encoding="utf-8", errors="replace").read().replace("\r\n", "\n").split("\n")
        ia = [l for l in a if ign.search(l)]; ib = [l for l in b if ign.search(l)]
        ka = [l for l in a if not ign.search(l)]; kb = [l for l in b if not ign.search(l)]
        same_bytes = hashlib.sha256(open(g, "rb").read()).hexdigest() == hashlib.sha256(open(h, "rb").read()).hexdigest()
        diff = sum(1 for x, y in zip(ka, kb) if x != y) + abs(len(ka) - len(kb))
        tot_bad += diff > 0
        print("REGRESS | %s | lines %d / %d | lines differing outside the ignored version and commit lines: %d | ignored lines golden %d head %d | byte-identical: %s | %s" %
              (os.path.relpath(g, root).replace("golden" + os.sep, ""), len(a), len(b), diff, len(ia), len(ib), same_bytes, "PASS" if diff == 0 else "FAIL"))
        if diff and diff < 6:
            for x, y in zip(ka, kb):
                if x != y: print("   golden:", x[:200]); print("   head  :", y[:200])
print("REGRESS | files differing outside the ignored lines:", tot_bad)
