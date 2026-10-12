# Stage 04: Verify

Prove the module against oracles and an independent review.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-spec/output/[module]-spec.md` | "Acceptance criteria", "Oracle" | what to prove |
| Stage 03 | `../03-implement/output/[module]-impl-notes.md` | Full file | deviations to check |
| Code | `/src/<Module>/`, `/test/<module>/` | Full | what is verified |

## Process

1. Run the oracle cross-check at each precision the spec names; make the convergence or order plot.
2. Independent review (fresh Sonnet subagent, given the spec and the code, not the notes): math versus code, numerical traps.
3. Write the verification report: evidence for each acceptance criterion, review findings and their fixes.
4. **[Checkpoint]** Present the report.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Evidence per criterion, plot, review findings | Accept, or send back to 01, 02 or 03 |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Evidence | Every acceptance criterion has a number or a plot |
| Review closed | No open CRITICAL or HIGH review finding |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Verification report | `output/[module]-verification.md` | Evidence table, review findings (read by 05) |
| Figures | `output/[module]-*.png` | Convergence and oracle plots |
