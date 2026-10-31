"""O26 (Amendment D94, E40): independent writer of the round log.

    python -I verify/experiments-d94-log.py --panel P.toml --prereg-sha HEX --out DIR
        --rounds 3 --decisions continue,continue,close --counts T_s=10,W2=3
        (--store DIR_WITH_chunk-*.h5 [--each-round-arms ref] | --chunks FILE)

Entry text (pinned): lines round=<r>, decision=<continue|close>, count.<name>=<int> in name order, panel_sha256=, prereg_sha256=,
ref_digest=, prev_sha256=; each ends in LF. Log = concatenation. panel_sha256 = SHA-256 of the panel file bytes. prev_sha256 = SHA-256
of the previous entry text (first: SHA-256 of the empty text). ref_digest = SHA-256 of the sorted chunk hashes (hex, each + LF) of the
chunks of the each_round arms of rounds 1..r; the round is the suffix `:r<round>` of the chunk key.
Chunk hashes: from the `chunk_hash` / `chunk_key` / `arms` root attributes of the chunk files read with h5dump when it is installed;
otherwise from --chunks FILE, one line per chunk: `<chunk_key> <chunk_hash> [arm,arm...]` (taken from the plan).
Writes DIR/round-log.txt (valid) and DIR/round-log-bad-<field>-<entry>.txt (one altered byte in one field of one entry).
"""
import argparse, glob, hashlib, os, re, shutil, subprocess, sys

sha = lambda b: hashlib.sha256(b).hexdigest()


def chunks_from_h5(store):
    out = []
    for f in sorted(glob.glob(os.path.join(store, "chunk-*.h5"))):
        t = subprocess.run(["h5dump", "-A", "-g", "/", f], capture_output=True, text=True).stdout
        def attr(n):
            m = re.search(r'ATTRIBUTE "%s" \{.*?DATA \{(.*?)\n\s*\}\s*\n\s*\}' % re.escape(n), t, re.S)
            return re.findall(r'"([^"]*)"', m.group(1)) if m else []
        out.append((attr("chunk_key")[0], attr("chunk_hash")[0], attr("arms")))
    return out


def chunks_from_raw(store):
    """no h5dump: the root attributes chunk_hash and chunk_key are fixed-length ASCII strings stored in the object header; scan the bytes.
    The `arms` attribute is a variable-length string (global heap) and is not read: the arm is taken as the method id of the key
    (pass --each-round-arms with method ids, e.g. L1)."""
    out = []
    for f in sorted(glob.glob(os.path.join(store, "chunk-*.h5"))):
        b = open(f, "rb").read()
        i = b.index(b"chunk_hash\x00")
        h = re.search(rb"[0-9a-f]{64}", b[i:i + 400]).group().decode()
        j = b.index(b"chunk_key\x00")
        k = re.search(rb"method=[ -~]*", b[j:j + 600]).group().decode()
        out.append((k, h, [re.search(r"method=(\S+)", k).group(1)]))
    return out


def chunks_from_file(p):
    out = []
    for ln in open(p, encoding="utf-8"):
        ln = ln.strip()
        if ln and not ln.startswith("#"):
            w = ln.split()
            out.append((w[0], w[1], w[2].split(",") if len(w) > 2 else []))
    return out


def round_of(key):
    m = re.search(r":r(\d+)$", key)
    return int(m.group(1)) if m else None


def entry(r, dec, counts, psha, rsha, ref, prev):
    s = "round=%d\ndecision=%s\n" % (r, dec)
    for k in sorted(counts):
        s += "count.%s=%d\n" % (k, counts[k])
    return s + "panel_sha256=%s\nprereg_sha256=%s\nref_digest=%s\nprev_sha256=%s\n" % (psha, rsha, ref, prev)


def build(rounds, decisions, counts, psha, rsha, chunks, ref_arms):
    es, prev = [], sha(b"")
    for r in range(1, rounds + 1):
        hs = sorted(h for k, h, arms in chunks if round_of(k) is not None and round_of(k) <= r
                    and (not ref_arms or any(a in ref_arms for a in arms)))
        ref = sha("".join(h + "\n" for h in hs).encode())
        e = entry(r, decisions[r - 1], counts[r - 1], psha, rsha, ref, prev)
        es.append(e)
        prev = sha(e.encode())
    return es


def flip(s):
    """alter one byte: hex digit -> another hex digit, digit -> other digit, letter -> other letter"""
    c = s[0]
    n = "0" if c != "0" else "1"
    if c.isalpha() and not c in "abcdef":
        n = "a"
    elif c in "abcdef":
        n = "0"
    return n + s[1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--prereg-sha", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rounds", type=int, required=True)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--counts", default="T_s=10,W2=3", help="name=int,... same at every round, or per round separated by ;")
    ap.add_argument("--store")
    ap.add_argument("--chunks")
    ap.add_argument("--each-round-arms", default="", help="comma list of arm names; empty = every arm in the file counts")
    a = ap.parse_args()
    decisions = a.decisions.split(",")
    assert len(decisions) == a.rounds
    parts = a.counts.split(";")
    cl = [dict((kv.split("=")[0], int(kv.split("=")[1])) for kv in p.split(",") if kv) for p in parts]
    counts = [cl[min(i, len(cl) - 1)] for i in range(a.rounds)]
    psha = sha(open(a.panel, "rb").read())
    if a.store and shutil.which("h5dump"):
        chunks, src = chunks_from_h5(a.store), "h5dump"
    elif a.store and not a.chunks:
        chunks, src = chunks_from_raw(a.store), "raw byte scan of the chunk files (h5dump not installed); arm = method id"
    elif a.chunks:
        chunks, src = chunks_from_file(a.chunks), "plan file " + a.chunks
    else:
        sys.exit("no h5dump and no --chunks file")
    ref_arms = [x for x in a.each_round_arms.split(",") if x]
    os.makedirs(a.out, exist_ok=True)
    es = build(a.rounds, decisions, counts, psha, a.prereg_sha, chunks, ref_arms)
    open(os.path.join(a.out, "round-log.txt"), "w", newline="\n", encoding="utf-8").write("".join(es))
    print("chunk hash source: %s; %d chunks; %d entries; panel_sha256 %s" % (src, len(chunks), len(es), psha))
    fields = ["round", "decision", "panel_sha256", "prereg_sha256", "ref_digest", "prev_sha256"] + ["count." + k for k in sorted(counts[0])]
    for ei in range(len(es)):
        for f in fields:
            lines = es[ei].split("\n")
            for li, ln in enumerate(lines):
                if ln.startswith(f + "="):
                    k, v = ln.split("=", 1)
                    nv = {"decision": "continuf" if v == "continue" else "closf"}.get(f) or (str(int(v) + 1) if v.isdigit() else flip(v))
                    lines[li] = k + "=" + nv
            bad = es[:ei] + ["\n".join(lines)] + es[ei + 1:]
            name = "round-log-bad-%s-%d.txt" % (f.replace(".", "_"), ei + 1)
            open(os.path.join(a.out, name), "w", newline="\n", encoding="utf-8").write("".join(bad))
    print("wrote valid log and %d altered copies to %s" % (len(es) * len(fields), a.out))


main()
