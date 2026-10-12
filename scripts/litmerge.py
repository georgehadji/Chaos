"""Merge per-source search tables and list identifiers to verify (literature-review/02 and /03).

Usage:
  python scripts/litmerge.py TOPIC              merge data/literature/TOPIC/search-*.md into merged.md
  python scripts/litmerge.py TOPIC --to-verify  write data/literature/TOPIC/to-verify.md from the stage 02 candidate list
Optional: --root DIR (default: current directory).

Merge: one row per paper across sources (rows sharing a DOI or arXiv ID; OpenAlex's 10.48550/arxiv.* DOIs
count as arXiv IDs); a paper is included if any source included it; papers already in references.bib
(by DOI or arXiv ID) are dropped and listed at the end. Narrowing rules are applied afterwards, by hand,
and recorded in the candidate list.
"""
import re
import sys
from pathlib import Path

DATA = "data/literature"
CANDIDATES = "workspaces/literature-review/stages/02-search/output/{}-candidates.md"


def ids(cell):
    """(DOIs, arXiv IDs) named in an identifier cell, lower-cased, version suffixes dropped."""
    c = cell.lower()
    d = {x.rstrip(".") for x in re.findall(r"doi:\s*(10\.[^\s|/]+/[^\s|]+)", c)}
    a = {re.sub(r"v\d+$", "", x) for x in re.findall(r"arxiv:\s*([\w.-]+/?\d*(?:\.\d+)?)", c)}
    a |= {x.split("arxiv.", 1)[1] for x in d if x.startswith("10.48550/arxiv.")}
    return {x for x in d if not x.startswith("10.48550/")}, a


def table_rows(text, width=7):
    """Data rows of the Markdown tables in text with the given column count."""
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == width and cells[0] != "identifier" and not set(cells[0]) <= set("-: "):
            yield cells


def merge(sources, bib_text):
    """sources: {source_name: table_text}. Returns (new_rows, rows_already_in_bib)."""
    bib = bib_text.lower()
    bib_d = set(re.findall(r"doi = \{([^}]+)\}", bib))
    bib_a = set(re.findall(r"eprint = \{([^}]+)\}", bib))
    rows, seen = [], {}
    for src, text in sorted(sources.items()):
        for cells in table_rows(text):
            d, a = ids(cells[0])
            if not d and not a:
                continue
            hit = next((seen[k] for k in [*d, *a] if k in seen), None)
            if hit is not None:
                r = rows[hit]
                r["sources"].add(src)
                r["doi"] |= d
                r["arxiv"] |= a
                if cells[5].lower().startswith("incl"):
                    r["cells"][5] = "include"
                continue
            seen.update({k: len(rows) for k in [*d, *a]})
            rows.append({"cells": cells, "sources": {src}, "doi": set(d), "arxiv": set(a),
                         "known": bool((d & bib_d) or (a & bib_a))})
    return [r for r in rows if not r["known"]], [r for r in rows if r["known"]]


def ident(r):
    return " / ".join([f"doi:{x}" for x in sorted(r["doi"])] + [f"arXiv:{x}" for x in sorted(r["arxiv"])])


def write_merged(root, topic):
    d = root / DATA / topic
    sources = {f.stem.removeprefix("search-"): f.read_text(encoding="utf-8") for f in sorted(d.glob("search-*.md"))}
    rows, known = merge(sources, (root / "_config/bibliography/references.bib").read_text(encoding="utf-8"))
    inc = [r for r in rows if r["cells"][5].lower().startswith("incl")]
    exc = [r for r in rows if r not in inc]
    out = ["| identifier | first author et al. | year | title | SQ | include/exclude | reason | sources |",
           "|---|---|---|---|---|---|---|---|"]
    for r in sorted(inc, key=lambda r: (r["cells"][4], r["cells"][2])) + sorted(exc, key=lambda r: r["cells"][4]):
        out.append("| " + " | ".join([ident(r)] + r["cells"][1:] + [", ".join(sorted(r["sources"]))]) + " |")
    out.append("\nAlready in references.bib (dropped): " + "; ".join(ident(r) for r in known))
    (d / "merged.md").write_text("\n".join(out).replace(chr(0x2014), " - ") + "\n", encoding="utf-8", newline="\n")
    print(f"{topic}: {len(rows) + len(known)} unique, {len(inc)} included, {len(exc)} excluded listed, "
          f"{len(known)} already in references.bib")


def write_to_verify(root, topic):
    text = (root / CANDIDATES.format(topic)).read_text(encoding="utf-8")
    lines = []
    for cells in table_rows(text, width=8):
        if cells[5].lower().startswith("incl"):
            d, a = ids(cells[0])
            lines.append(f"- doi:{sorted(d)[0]}" if d else f"- arXiv:{sorted(a)[0]}")
    (root / DATA / topic / "to-verify.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"{topic}: {len(lines)} identifiers to verify (DOI preferred)")


def main(argv):
    root = Path(argv[argv.index("--root") + 1]) if "--root" in argv else Path(".")
    skip = argv.index("--root") + 1 if "--root" in argv else -1
    args = [a for i, a in enumerate(argv) if i > 0 and i != skip and not a.startswith("--")]
    if not args:
        print(__doc__)
        return 2
    (write_to_verify if "--to-verify" in argv else write_merged)(root, args[0])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
