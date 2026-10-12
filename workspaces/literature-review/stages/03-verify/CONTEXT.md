# Stage 03: Verify

Check every candidate against its fetched record and write verified BibTeX.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../02-search/output/[topic]-candidates.md` | included rows only | what to verify |
| Config | `/_config/bibliography/conventions.md` | Full file | key format, verification classes |
| Tool | `/scripts/bibcheck.py` | run, do not read | resolves identifiers on Crossref and arXiv |

## Process

1. Run `python scripts/bibcheck.py` on the candidate list (no AI). It fetches each record and reports mismatches in title, authors, venue, volume, pages.
2. Fix mismatches against the fetched record only (Sonnet subagent).
3. Classify each entry: verified, verified-secondary, or UNVERIFIED (per `conventions.md`).
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
