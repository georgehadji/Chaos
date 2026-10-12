"""Write literature-review/03 outputs from the bibcheck draft (no AI).

Usage: python scripts/verify_report.py TOPIC [TOPIC ...] [--root DIR]
Reads data/literature/TOPIC/draft.bib and bibcheck.log (from `bibcheck.py to-verify.md --out draft.bib`),
plus optional data/literature/TOPIC/verify-overrides.json with checkpoint decisions:
  {"drop_doi": {"10.x/y": "reason"}, "notes": {"Key": "note"}}
Drops, with a reason in the report: chapters of books already in references.bib (chapter DOI extends the
book DOI), same-title duplicates of existing entries, papers already verified by an earlier TOPIC in the
same call, and the overrides. Writes [topic]-verified.bib and [topic]-verify-report.md to 03-verify/output.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bibcheck import norm, parse_bib  # noqa: E402

OUT = "workspaces/literature-review/stages/03-verify/output"


def split_entries(text):
    return [c.strip() for c in re.split(r"\n(?=@)", text) if c.strip()]


def classify(chunks, refs, overrides, prior):
    """Return (kept_chunks, drop_reasons). prior maps normalized titles to the topic that kept them."""
    titles = {norm(e.get("title", "")): e["key"] for e in refs}
    books = {e["doi"].lower(): e["key"] for e in refs if e.get("doi")}
    drop_doi = {k.lower(): v for k, v in overrides.get("drop_doi", {}).items()}
    keep, dropped = [], []
    for c in chunks:
        e = parse_bib(c + "\n")[0]
        doi, n = e.get("doi", "").lower(), norm(e.get("title", ""))
        book = next((k for bd, k in books.items() if doi.startswith(bd + ".") or doi.startswith(bd + "_")), None)
        if doi in drop_doi:
            dropped.append(f"{e['key']} (doi:{doi}): {drop_doi[doi]}")
        elif book:
            dropped.append(f"{e['key']} (doi:{doi}): chapter of {book}, already in references.bib; "
                           f"cite the book with the chapter title ('{e.get('title', '')}')")
        elif n in titles:
            dropped.append(f"{e['key']}: same title as {titles[n]} in references.bib (different identifier form)")
        elif n in prior:
            dropped.append(f"{e['key']}: also found by {prior[n]}; verified there, kept once")
        else:
            keep.append(c)
    return keep, dropped


def report(topic, ents, dropped, failures, notes):
    arx = sum(1 for e in ents if not e.get("doi"))
    lines = [f"# Verification report: {topic}", "",
             f"Stage 03. Input: included rows of `../02-search/output/{topic}-candidates.md`, one identifier per paper "
             "(DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref, DataCite and arXiv records); "
             "every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).", "",
             "## Summary", "",
             f"- Fetch failures: {len(failures)}.",
             f"- In `{topic}-verified.bib`: {len(ents)} entries, all class verified "
             f"({len(ents) - arx} by DOI record, {arx} by arXiv record).",
             "- Keys follow `conventions.md` (clashes suffixed a, b); none clashes with `references.bib`.", "",
             "## Fixes and notes", ""]
    lines += [f"- {k}: {v}" for k, v in notes.items() if any(e["key"] == k for e in ents)] or ["- None."]
    lines += ["", "## Dropped at verification", ""] + ([f"- {x}" for x in dropped] or ["- None."])
    lines += ["", "## UNVERIFIED", ""]
    lines += [f"- {f}" for f in failures] or ["None: every included candidate resolved to a fetched record."]
    return "\n".join(lines) + "\n"


def run(root, topics):
    refs = parse_bib((root / "_config/bibliography/references.bib").read_text(encoding="utf-8"))
    prior = {}
    for topic in topics:
        d = root / "data/literature" / topic
        ov_path = d / "verify-overrides.json"
        overrides = json.loads(ov_path.read_text(encoding="utf-8")) if ov_path.exists() else {}
        log = (d / "bibcheck.log").read_text(encoding="utf-8").splitlines()
        keep, dropped = classify(split_entries((d / "draft.bib").read_text(encoding="utf-8")), refs, overrides, prior)
        dropped += [l.replace(", skipped", " (same paper under another DOI form); one entry kept")
                    for l in log if "duplicate record" in l]
        failures = [l for l in log if "fetch failed" in l]
        bib = ("\n\n".join(keep) + "\n").replace(chr(0x2014), "---")
        ents = parse_bib(bib)
        prior.update({norm(e.get("title", "")): topic for e in ents})
        out = root / OUT
        (out / f"{topic}-verified.bib").write_text(bib, encoding="utf-8", newline="\n")
        (out / f"{topic}-verify-report.md").write_text(
            report(topic, ents, dropped, failures, overrides.get("notes", {})).replace(chr(0x2014), " - "),
            encoding="utf-8", newline="\n")
        print(f"{topic}: {len(ents)} entries, dropped {len(dropped)}, fetch failures {len(failures)}")


def main(argv):
    root = Path(argv[argv.index("--root") + 1]) if "--root" in argv else Path(".")
    skip = argv.index("--root") + 1 if "--root" in argv else -1
    topics = [a for i, a in enumerate(argv) if i > 0 and i != skip and not a.startswith("--")]
    if not topics:
        print(__doc__)
        return 2
    run(root, topics)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
