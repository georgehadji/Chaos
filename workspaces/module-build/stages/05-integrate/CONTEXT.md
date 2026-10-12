# Stage 05: Integrate

Flip the catalog status, update docs, and commit.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../04-verify/output/[module]-verification.md` | Full file | accepted evidence |
| Config | `/_config/catalogs/` | this module's rows | status to flip |
| Code | `/src/<Module>/CONTEXT.md`, `/CHANGELOG.md` | Full files | docs to update |

## Process

1. Flip the catalog `Status` of this module's rows to `verified` (promotion write into `/_config/`).
2. Update the folder `CONTEXT.md` and `CHANGELOG.md`.
3. **[Checkpoint]** Present the commit.
4. Commit on a branch. Push or PR only with explicit approval.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Diff and commit message | Approve the commit; push yes or no |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Status flip | `/_config/catalogs/` | promotion write (read by study/01 and later module runs) |
| Integration log | `output/[module]-integration-log.md` | Commit SHA, files changed (final deliverable) |
