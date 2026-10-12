# Stage 01: Outline

Agree the contribution and map every claim to evidence.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | paper slug, target journal | per-run input |
| Config | `/_config/research-question.md` | "Paper sequence", this paper's row | scope and claim |
| Study | `/workspaces/study/stages/05-report/output/[study]-results.md` | each study this paper uses | the evidence |
| Config | `/_config/value-types.md` | Full file | the contribution (Pattern 13) |
| Config | `/_config/publication.md` | Full file | journals, citation style, author data |

## Process

1. Agree the contribution and value types with the human.
2. Draft the outline (Sonnet subagent, ARS `/ars-outline`): one claim per section, each claim mapped to a study result or a citation.
3. **[Checkpoint]** Present the outline.
4. Run the audit, write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Outline with claim-to-evidence map | Approve, reorder, or cut |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Evidence | Every claim maps to a study result row or a verified key |
| Accepted studies | Every study used has an accepted `[study]-results.md` |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Outline | `output/[paper]-outline.md` | Sections, claims, evidence map (read by 02) |
