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
2. Update the folder `CONTEXT.md` and `CHANGELOG.md`; run the audit.
3. **[Checkpoint]** Present the commit.
4. Commit on a branch and push it. PR only with explicit approval.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Diff and commit message | Approve the commit; PR yes or no |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Status matches evidence | A row is `verified` only if its verification report has no open criterion; a qualifier in parentheses names any deferred part |
| Log | The integration log names the commit SHA and every file changed |
| Changelog | `CHANGELOG.md` has an entry for the module or amendment |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Status flip | `/_config/catalogs/` | promotion write (read by study/01 and later module runs) |
| Integration log | `output/[module]-integration-log.md` | Commit SHA, files changed (final deliverable) |
