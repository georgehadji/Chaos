# Stage 03: Verify

Check every candidate against its fetched record and write verified BibTeX.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../02-search/output/[topic]-candidates.md` | included rows only | what to verify |
| Config | `/_config/bibliography/conventions.md` | Full file | key format, verification classes |
| Tool | `/scripts/bibcheck.py` | run, do not read | resolves identifiers on Crossref, DataCite and arXiv |
| Tool | `/scripts/litmerge.py`, `/scripts/verify_report.py` | run, do not read | identifier list in; verified .bib and report out |

## Process

1. Run `python scripts/litmerge.py [topic] --to-verify`, then `python scripts/bibcheck.py` on that list with `--out` a draft .bib (no AI). It fetches each record and reports failures and duplicate records.
2. Fix mismatches against the fetched record only (Sonnet subagent).
3. Classify each entry: verified, verified-secondary, or UNVERIFIED (per `conventions.md`); `python scripts/verify_report.py [topic]` writes the outputs and drops known duplicates (checkpoint decisions go in `verify-overrides.json` beside the draft).
4. **[Checkpoint]** Present failures and UNVERIFIED items.
5. Run the audit, write the outputs.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Mismatch fixes, failures, UNVERIFIED list | Drop, fix by hand, or keep as UNVERIFIED |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Clean bib | Every entry in the .bib is verified or verified-secondary |
| Keys | Keys follow `conventions.md` and are unique against `/_config/bibliography/references.bib` |
| Report | Every UNVERIFIED item states what failed |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Verified entries | `output/[topic]-verified.bib` | BibTeX |
| Verification report | `output/[topic]-verify-report.md` | Fixes, failures, UNVERIFIED list (read by the human at the checkpoint and by 05-promote) |
