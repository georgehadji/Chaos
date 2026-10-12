# Stage 01: Design

Write and freeze the study preregistration.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | study slug, paper it serves | per-run input |
| Config | `/_config/research-question.md` | "Three claims to test", "Test panel", "Outputs retained", the paper's row | what the study must decide |
| Config | `/_config/catalogs/` | rows in scope, `Status` column | only verified methods may run |
| Config | `/_config/value-types.md` | Full file | agree the value the study delivers (Pattern 13) |
| Config | `/_config/numerics-defaults.toml` | Full file | defaults the study may override |
| Config | `/_config/research-ideas.md` | the idea row this study serves | scope |
| Config | `/_config/environment.md` | Full file | compute limits for the budget |
| Reference | `references/preregistration-format.md` | Full file | required sections |

## Process

1. Agree the study scope and its value types with the human.
2. Draft the preregistration (Sonnet subagent): hypotheses, panel, methods, metrics, decision thresholds, compute budget, analysis plan, stopping rules.
3. **[Checkpoint]** Present the preregistration.
4. Run the audit; freeze: record the git SHA and the file's SHA-256 in its header and commit it. Later changes are amendments, logged before any data exists.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Scope and value types | Confirm the study is worth running |
| 2 | Full preregistration | Approve and freeze, or revise |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Verified methods | Every method in scope has catalog `Status` = `verified` |
| Decidable | Every hypothesis has a metric and a numeric decision rule |
| Single home | Thresholds appear only in this file |
| Budget | Estimated CPU hours fit the limits in `/_config/environment.md` |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Preregistration | `output/[study]-preregistration.md` | Frozen; header holds git SHA and SHA-256 |
