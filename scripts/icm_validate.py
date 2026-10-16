"""Validate the ICM workspace against _core/CONVENTIONS.md (the Stage 05 checks, rerunnable).

Usage: python scripts/icm_validate.py [repo_root]
Exit code 1 if any check fails. Placeholders left for `setup` and paths built in later
milestones are reported as information, not failures.
"""
import re
import sys
from pathlib import Path

WORKSPACES = ["literature-review", "module-build", "study", "paper"]
CONTRACT_SECTIONS = {"Inputs", "Process", "Checkpoints", "Audit", "Outputs"}
ICM_DIRS = ("workspaces", "_config", "setup", "schema")  # lowercase-hyphen rule applies here
SKIP = {".git", "_core", "data", ".ruff_cache", "__pycache__", "workspace-builder"}
# git-ignored N-body environments and build outputs (.gitignore): third-party code, not ICM content
SKIP_DIRS = tuple(f"python/nbody/windows-build/{d}/" for d in ("venv", "_work", "wheels"))
LATER = ("/scripts", "/bin", "/test", "/src", "/bench", "/CHANGELOG", "/python")
CREATIVE = {("literature-review", "01-scope"), ("literature-review", "02-search"),
            ("literature-review", "04-synthesize"), ("module-build", "01-spec"),
            ("module-build", "02-tests"), ("study", "01-design"), ("study", "04-analyze"),
            ("study", "05-report"), ("paper", "01-outline"), ("paper", "02-draft"),
            ("paper", "03-review")}
BUILD = {("module-build", "03-implement")}
NAME_RE = re.compile(r"_?[a-z0-9][a-z0-9\-]*(\.[a-z0-9]+)*")
FORMATS = {"preregistration": "workspaces/study/stages/01-design/references/preregistration-format.md",
           "spec": "workspaces/module-build/stages/01-spec/references/spec-format.md"}


def section(text, name):
    m = re.search(rf"^## {name}\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else None


def load_stages(root):
    stages = {}
    for ws in WORKSPACES:
        for d in sorted((root / "workspaces" / ws / "stages").iterdir()):
            stages[(ws, d.name)] = d
    return stages


def skipped(rel):
    return (bool(set(rel.parts) & SKIP) or any(part.startswith(".venv") for part in rel.parts)
            or rel.as_posix().startswith(SKIP_DIRS))


def text_files(root):
    for p in root.rglob("*"):
        if p.is_file() and not skipped(p.relative_to(root)) and p.name != "LICENSE":
            if p.suffix in {".md", ".toml", ".jl", ".py", ".yml", ".cff", ".bib"}:
                yield p


def check_structure(root, stages, fails):
    for (ws, st), d in stages.items():
        for need in ("CONTEXT.md", "output/.gitkeep", "references"):
            if not (d / need).exists():
                fails.append(f"structure: {ws}/{st} lacks {need}")


def check_files(root, fails, info):
    placeholders = set()
    for p in text_files(root):
        rel = p.relative_to(root).as_posix()
        t = p.read_text(encoding="utf-8", errors="replace")
        n = t.count("\n")
        if p.name == "CONTEXT.md" and n >= 80:
            fails.append(f"size: {rel} has {n} lines (CONTEXT limit 80)")
        if p.suffix == ".md" and rel.startswith(("_config", "schema", "workspaces")) \
                and p.name not in {"CONTEXT.md", "CLAUDE.md"} and "output" not in p.parts and n >= 200:
            fails.append(f"size: {rel} has {n} lines (reference limit 200)")
        if chr(0x2014) in t:  # em dash
            fails.append(f"em dash: {rel}")
        found = set(re.findall(r"\{\{([A-Z_]+)\}\}", t))
        if found and (p.name == "CLAUDE.md" or rel == "CONTEXT.md"):
            fails.append(f"placeholder in a file that must work before setup: {rel}")
        if found and rel != "setup/questionnaire.md":
            placeholders |= found
        if rel.startswith(ICM_DIRS):
            bad = [part for part in p.relative_to(root).parts
                   if part not in {"CLAUDE.md", "CONTEXT.md", "references.bib"} and not NAME_RE.fullmatch(part.lstrip("."))]
            if bad:
                fails.append(f"naming: {rel} ({bad[0]})")
    return placeholders


def check_placeholders(root, stages, placeholders, fails, info):
    q = (root / "setup/questionnaire.md").read_text(encoding="utf-8")
    asked = set(re.findall(r"\{\{([A-Z_]+)\}\}", q))
    for name in sorted(placeholders - asked):
        fails.append(f"placeholder without a question: {{{{{name}}}}}")
    for (ws, st), d in stages.items():
        t = (d / "CONTEXT.md").read_text(encoding="utf-8")
        if t.count("{{") != (section(t, "Inputs") or "").count("{{"):
            fails.append(f"placeholder outside Inputs: {ws}/{st}")
    if placeholders:
        info.append(f"setup pending: {len(placeholders)} placeholders not yet filled")


def check_inputs(root, stages, fails, pending):
    for (ws, st), d in stages.items():
        inputs = section((d / "CONTEXT.md").read_text(encoding="utf-8"), "Inputs") or ""
        for row in inputs.splitlines()[2:]:
            cells = [c.strip() for c in row.strip("|").split("|")]
            if len(cells) < 3:
                continue
            base = None
            for path in re.findall(r"`([^`]+)`", cells[1]):
                if path.startswith("../"):
                    kind = next((k for k in FORMATS if k in path), None)
                    base = root / FORMATS[kind] if kind else base
                    continue
                target = root / path.lstrip("/") if path.startswith("/") else d / path
                if any(c in path for c in "<[*") or not path.startswith(("/", "references/")):
                    continue  # runtime paths (data/...) and per-run names are not repo files
                if not target.exists():
                    (pending.add(path) if path.startswith(LATER) else fails.append(f"missing path: {ws}/{st} -> {path}"))
                elif target.suffix == ".md":
                    base = target
            quoted = re.findall(r'"([^"]+)"', cells[2])
            if quoted and base is not None and base.exists():
                heads = set(re.findall(r"^#+ (.+?)\s*$", base.read_text(encoding="utf-8"), re.M))
                fails.extend(f'section: {ws}/{st} quotes "{q}", not a heading in {base.relative_to(root).as_posix()}'
                             for q in quoted if q not in heads)


def read_edges(stages):
    edges = []
    pat = re.compile(r"`(?:\.\./(0\d-[a-z-]+)|/workspaces/([a-z-]+)/stages/(0\d-[a-z-]+))/output/([^`]+)`")
    for (ws, st), d in stages.items():
        for m in pat.finditer(section((d / "CONTEXT.md").read_text(encoding="utf-8"), "Inputs") or ""):
            producer = (ws, m.group(1)) if m.group(1) else (m.group(2), m.group(3))
            edges.append(((ws, st), producer, m.group(4).split(",")[0].strip()))
    return edges


def check_graph(root, stages, fails):
    edges = read_edges(stages)
    for consumer, producer, art in edges:
        if consumer[0] == producer[0] and producer[1] >= consumer[1]:
            fails.append(f"cycle risk: {consumer} reads later stage {producer}")
        outs = section((stages[producer] / "CONTEXT.md").read_text(encoding="utf-8"), "Outputs") or ""
        if art.split("-", 1)[-1] not in outs:
            fails.append(f"handoff gap: {consumer} reads {art} not in {producer} Outputs")
    for base in ("_config", "schema"):
        for p in (root / base).rglob("*.md"):
            for m in re.findall(r"workspaces/[a-z-]+[^\s`)]*|[Rr]ead by [a-z-]+/\d\d", p.read_text(encoding="utf-8")):
                fails.append(f"back-reference: {p.relative_to(root).as_posix()} -> {m}")


def check_contracts(stages, fails):
    for (ws, st), d in stages.items():
        t = (d / "CONTEXT.md").read_text(encoding="utf-8")
        fails.extend(f"CONTEXT purity: {ws}/{st} has section {h}"
                     for h in re.findall(r"^## (.+)$", t, re.M) if h not in CONTRACT_SECTIONS)
        steps = len(re.findall(r"^\d+\. ", section(t, "Process") or "", re.M))
        cp, au = section(t, "Checkpoints"), section(t, "Audit")
        if (ws, st) in CREATIVE and not cp:
            fails.append(f"checkpoint: creative stage {ws}/{st} has none")
        fails.extend(f"checkpoint: {ws}/{st} names step {n} of {steps}"
                     for n in re.findall(r"^\| (\d+) \|", cp or "", re.M) if not 1 <= int(n) <= steps)
        if (ws, st) in CREATIVE | BUILD and not au:
            fails.append(f"audit: {ws}/{st} has none")
        if au and not re.search(r"audit", section(t, "Process") or "", re.I):
            fails.append(f"audit: {ws}/{st} process never runs it")


def validate(root):
    fails, info, pending = [], [], set()
    stages = load_stages(root)
    check_structure(root, stages, fails)
    placeholders = check_files(root, fails, info)
    check_placeholders(root, stages, placeholders, fails, info)
    check_inputs(root, stages, fails, pending)
    check_graph(root, stages, fails)
    check_contracts(stages, fails)
    if pending:
        info.append("built in later milestones: " + ", ".join(sorted(pending)))
    return fails, info


def main(argv):
    root = Path(argv[1] if len(argv) > 1 else ".").resolve()
    fails, info = validate(root)
    for line in fails:
        print("FAIL", line)
    for line in info:
        print("INFO", line)
    print("ALL CHECKS PASS" if not fails else f"{len(fails)} FAILURES")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
