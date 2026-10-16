# Stage 05: Promote

Merge the verified results into the shared bibliography after human approval.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 03 | `../03-verify/output/[topic]-verified.bib`, `[topic]-verify-report.md` | Full files | entries and UNVERIFIED list |
| Stage 04 | `../04-synthesize/output/[topic]-notes.md` | Full file, if it exists | annotations and benchmark values |
| Tool | `/scripts/promote.py` | run, do not read | builds the merge diff |

## Process

1. Run `python scripts/promote.py [topic]` (no AI): it proposes appends to `references.bib`, `benchmarks.md` and `unverified.md`, plus catalog `?` keys that now resolve, and blocks on duplicate keys or DOIs. Draft the matching annotation additions to the bibliography section files from the notes (Sonnet subagent).
2. **[Checkpoint]** Present the diff.
3. Apply the approved diff (`promote.py [topic] --apply` plus the approved annotations; a promotion write into `/_config/`), run the audit, commit, push.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Merge diff into `/_config/bibliography/` and catalog key flips | Approve, edit, or reject |

## Audit

| Check | Pass Condition |
|-------|---------------|
| No duplicates | Every key in `references.bib` is unique |
| Size | Every bibliography section file is under 200 lines (split if not) |
| Style | No em dashes |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Updated bibliography | `/_config/bibliography/` | promotion write (final deliverable) |
| Promotion log | `output/[topic]-promote-log.md` | What was merged, the commit SHA |
