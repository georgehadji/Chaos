"""Per-slug pipeline status for the ICM workspaces (the `status` and `status <slug>` triggers).

Usage: python scripts/icm_status.py [slug] [--root DIR]
A stage is COMPLETE for a slug when an artifact named in its Outputs table exists for that slug.
"""
import fnmatch
import re
import sys
from pathlib import Path

WORKSPACES = ["literature-review", "module-build", "study", "paper"]


def stage_patterns(context_text):
    """Output artifact patterns of a stage, e.g. '*-search-plan.md' or '*/' for a folder."""
    outs = re.search(r"^## Outputs\n(.*)", context_text, re.M | re.S)
    found = re.findall(r"`output/\[[a-z]+\]([^`]*)`", outs.group(1) if outs else "")
    return [p if p.endswith("/") else "*" + p for p in found]


def slugs_in(output_dir, patterns):
    """Map slug -> True for every file or folder in output_dir matching a pattern."""
    hits = set()
    if not output_dir.is_dir():
        return hits
    for entry in output_dir.iterdir():
        for pat in patterns:
            if pat == "/" and entry.is_dir():
                hits.add(entry.name)
            elif entry.is_file() and fnmatch.fnmatch(entry.name, pat):
                suffix = pat.lstrip("*")
                stem = re.escape(suffix).replace(r"\*", ".*")
                m = re.fullmatch(rf"(.+?){stem}", entry.name)
                if m:
                    hits.add(m.group(1))
    return hits


def workspace_status(root, ws):
    """Return (stage names, {slug: set of complete stages})."""
    stages = sorted(d for d in (root / "workspaces" / ws / "stages").iterdir() if d.is_dir())
    table = {}
    for d in stages:
        pats = stage_patterns((d / "CONTEXT.md").read_text(encoding="utf-8"))
        for slug in slugs_in(d / "output", pats):
            table.setdefault(slug, set()).add(d.name)
    # A run starts in the first stage: a later-stage name without a first-stage artifact (a figure such as
    # lyap-benettin-convergence.png cut at its first hyphen, an amendment plan such as tangent-core-d44) is no slug.
    first = stages[0].name if stages else None
    return [d.name for d in stages], {s: done for s, done in table.items() if first in done}


def catalog_counts(root):
    counts = {}
    for p in sorted((root / "_config" / "catalogs").glob("*.md")):
        for status in re.findall(r"\| (planned|verified)(?: \([^|]*\))? \|$", p.read_text(encoding="utf-8"), re.M):
            counts[status] = counts.get(status, 0) + 1
    return counts


def render(stages, done):
    return "  -->  ".join(f"[{s}] {'COMPLETE' if s in done else 'PENDING'}" for s in stages)


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    root = Path(argv[argv.index("--root") + 1]) if "--root" in argv else Path(".")
    if "--root" in argv:
        args.remove(str(root))
    slug = args[0] if args else None
    for ws in WORKSPACES:
        stages, table = workspace_status(root, ws)
        if slug is None:
            print(f"{ws}: {len(table)} slug(s)")
            for s, done in sorted(table.items()):
                print(f"  {s}: {len(done)}/{len(stages)} stages")
        elif slug in table:
            print(f"{ws} / {slug}\n  {render(stages, table[slug])}")
    if slug is None:
        print("catalogs: " + ", ".join(f"{k} {v}" for k, v in sorted(catalog_counts(root).items())))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
