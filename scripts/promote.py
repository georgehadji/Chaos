"""Propose (and with --apply, make) the mechanical merge of a literature topic into /_config/ (literature-review/05).

Usage: python scripts/promote.py TOPIC [--apply] [--root DIR]
Merges: new entries into references.bib (blocked on duplicate keys or DOIs), UNVERIFIED items,
benchmark rows from the notes, and catalog keys marked `?` that now resolve. Annotations for the
bibliography section files are drafted separately by the agent and approved in the same diff.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bibcheck import norm, parse_bib  # noqa: E402

LR = "workspaces/literature-review/stages"
BIB = "_config/bibliography"


def section_lines(text, heading_word):
    """Lines under the first '## ...heading_word...' heading, up to the next '## '."""
    m = re.search(rf"^## [^\n]*{heading_word}[^\n]*\n(.*?)(?=^## |\Z)", text, re.M | re.S | re.I)
    return m.group(1).splitlines() if m else []


def plan(root, topic):
    """Return a dict describing every proposed change, plus blocking problems."""
    out = root / LR
    bib_new_text = (out / f"03-verify/output/{topic}-verified.bib").read_text(encoding="utf-8")
    report = (out / f"03-verify/output/{topic}-verify-report.md").read_text(encoding="utf-8")
    notes_p = out / f"04-synthesize/output/{topic}-notes.md"
    notes = notes_p.read_text(encoding="utf-8") if notes_p.exists() else ""
    existing = parse_bib((root / BIB / "references.bib").read_text(encoding="utf-8"))
    keys = {e["key"] for e in existing}
    dois = {e.get("doi", "").lower() for e in existing if e.get("doi")}
    new = parse_bib(bib_new_text)
    problems = [f"duplicate key {e['key']}" for e in new if e["key"] in keys]
    problems += [f"duplicate DOI {e['doi']} ({e['key']})" for e in new if e.get("doi", "").lower() in dois]
    titles = {norm(e.get("title", "")): e["key"] for e in existing if e.get("title")}
    problems += [f"duplicate title {e['key']} = {titles[norm(e['title'])]}" for e in new
                 if e.get("title") and norm(e["title"]) in titles]
    unverified = [l for l in section_lines(report, "UNVERIFIED") if l.startswith("- ")]
    bench = [l for l in section_lines(notes, "Benchmark") if l.startswith("|") and not re.match(r"^\|[-| ]+\|$", l)][1:]
    resolvable = {e["key"] for e in new}  # only keys this topic adds; `?` on an existing key is a pending human decision
    flips = {}
    for cat in sorted((root / "_config/catalogs").glob("*.md")):
        hits = sorted(set(k for k in re.findall(r"\b([A-Z][A-Za-z]+(?:19|20)\d\d[a-z]?)\?", cat.read_text(encoding="utf-8")) if k in resolvable))
        if hits:
            flips[cat.relative_to(root).as_posix()] = hits
    return {"new": new, "bib_text": bib_new_text, "unverified": unverified, "bench": bench, "flips": flips, "problems": problems}


def show(p):
    print(f"references.bib: +{len(p['new'])} entries: " + ", ".join(e["key"] for e in p["new"]))
    print(f"unverified.md: +{len(p['unverified'])} items")
    print(f"benchmarks.md: +{len(p['bench'])} rows")
    for cat, ks in p["flips"].items():
        print(f"{cat}: resolve {', '.join(k + '?' for k in ks)}")
    for prob in p["problems"]:
        print("BLOCKED:", prob)


def apply(root, topic, p):
    bib = root / BIB
    append(bib / "references.bib", "\n" + p["bib_text"].strip() + "\n")
    if p["unverified"]:
        append(bib / "unverified.md", "\n" + "\n".join(p["unverified"]) + "\n")
    if p["bench"]:
        append(bib / "benchmarks.md", "\n".join(p["bench"]) + "\n")
    for cat, ks in p["flips"].items():
        path = root / cat
        text = path.read_text(encoding="utf-8")
        for k in ks:
            text = re.sub(rf"\b{k}\?", k, text)
        path.write_text(text, encoding="utf-8", newline="\n")
    log = root / LR / f"05-promote/output/{topic}-promote-log.md"
    log.write_text(f"# Promotion log: {topic}\n\nMerged keys: {', '.join(e['key'] for e in p['new'])}\n"
                   f"UNVERIFIED items: {len(p['unverified'])}\nBenchmark rows: {len(p['bench'])}\n"
                   f"Catalog keys resolved: {sum(len(v) for v in p['flips'].values())}\nCommit: record the SHA after committing.\n",
                   encoding="utf-8", newline="\n")
    print(f"applied; log at {log.relative_to(root).as_posix()}")


def append(path, text):
    with open(path, "a", encoding="utf-8", newline="\n") as f:
        f.write(text)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    root = Path(argv[argv.index("--root") + 1]) if "--root" in argv else Path(".")
    p = plan(root, argv[1])
    show(p)
    if p["problems"]:
        return 1
    if "--apply" in argv:
        apply(root, argv[1], p)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
