# Stage 02: Draft

Draft the manuscript from results and verified citations only.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-outline/output/[paper]-outline.md` | Full file | structure |
| Study | `/workspaces/study/stages/04-analyze/output/[study]-fig-*.pdf`, `[study]-tables.md` | used figures and tables | figures |
| Study | `/workspaces/study/stages/05-report/output/[study]-results.md` | Full file | numbers |
| Config | `/_config/bibliography/` | entries cited in the outline | citations |

## Process

1. Draft section by section (Sonnet subagent), copying numbers from results files only.
2. Run the audit; list any number or claim that could not be traced.
3. **[Checkpoint]** Present the draft and the untraced items.
4. Write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Draft and untraced items | Proceed to review, or redirect sections |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Traced | Every number matches a results or table row |
| Cited | Every key is verified; no UNVERIFIED item cited |
| Style | No em dashes |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Draft | `output/[paper]-draft.md` | Full manuscript in Markdown with BibTeX keys (read by 03) |
