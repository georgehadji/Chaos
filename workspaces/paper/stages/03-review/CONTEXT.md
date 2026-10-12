# Stage 03: Review

Adversarial review and citation check, then revise.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../02-draft/output/[paper]-draft.md` | Full file | what to review |
| Config | `/_config/bibliography/references.bib` | Full file | citation check |

## Process

1. Adversarial review (fresh Sonnet subagent, ARS `/ars-reviewer`): correctness, overclaiming, missing prior work.
2. Citation check (ARS `/ars-citation-check`) against `references.bib`.
3. **[Checkpoint]** Present the findings and a proposed change list.
4. Apply the approved changes; write the revised draft and the response table.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Review findings, citation report, change list | Which changes to apply |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Closed | Every finding is fixed or answered in the response table |
| Citations | Citation check reports no unresolved key |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Revised draft | `output/[paper]-revised.md` | Manuscript (read by 04) |
| Response table | `output/[paper]-review-response.md` | Finding, action, location (read by the human at the checkpoint) |
