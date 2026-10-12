# Stage 05: Report

Write the traced, cited results document.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 04 | `../04-analyze/output/[study]-findings.md`, `[study]-tables.md` | Full files | results |
| Stage 03 | `../03-run/output/[study]-manifest.json` | Full file | traceability |
| Config | `/_config/bibliography/` | entries cited | citations |

## Process

1. Write the results document (Sonnet subagent): every number traced to a table or the manifest, every claim cited, limitations stated.
2. **[Checkpoint]** Present the results.
3. Run the audit, write the output, commit.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Results document | Accept, or revise |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Traced | Every number points to a table row or manifest field |
| Cited | Every citation key is verified (no `?`, no UNVERIFIED) |
| Style | No em dashes |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Results | `output/[study]-results.md` | Final deliverable of the study (read by paper) |
