# Contracts: paper workspace (optional)

DRAFT (Stage 02). One run per paper. `[paper]` is a slug such as `p1-le-algorithms`. Numbers come only from study results; citations only from the verified bibliography. ARS skills are used from the installed plugin (not bundled).

## 01-outline

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | paper slug, target journal | per-run input |
| Config | `/_config/research-question.md` | "Paper sequence", this paper's row | scope and claim |
| Study | `/workspaces/study/stages/05-report/output/[study]-results.md` | each study this paper uses | the evidence |
| Config | `/_config/value-types.md` | Full file | the contribution (Pattern 13) |
| Config | `/_config/publication.md` | Full file | journals, citation style, author data |

### Process

1. Agree the contribution and value types with the human.
2. Draft the outline (Sonnet subagent, ARS `/ars-outline`): one claim per section, each claim mapped to a study result or a citation.
3. **[Checkpoint]** Present the outline.
4. Run the audit, write the output.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Outline with claim-to-evidence map | Approve, reorder, or cut |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Evidence | Every claim maps to a study result row or a verified key |
| Accepted studies | Every study used has an accepted `[study]-results.md` |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Outline | `output/[paper]-outline.md` | Sections, claims, evidence map (read by 02) |

## 02-draft

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-outline/output/[paper]-outline.md` | Full file | structure |
| Study | `/workspaces/study/stages/04-analyze/output/[study]-fig-*.pdf`, `[study]-tables.md` | used figures and tables | figures |
| Study | `/workspaces/study/stages/05-report/output/[study]-results.md` | Full file | numbers |
| Config | `/_config/bibliography/` | entries cited in the outline | citations |

### Process

1. Draft section by section (Sonnet subagent), copying numbers from results files only.
2. Run the audit, write the output. No checkpoint (review is stage 03).

### Audit

| Check | Pass Condition |
|-------|---------------|
| Traced | Every number matches a results or table row |
| Cited | Every key is verified; no UNVERIFIED item cited |
| Style | No em dashes |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Draft | `output/[paper]-draft.md` | Full manuscript in Markdown with BibTeX keys (read by 03) |

## 03-review

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../02-draft/output/[paper]-draft.md` | Full file | what to review |
| Config | `/_config/bibliography/references.bib` | Full file | citation check |

### Process

1. Adversarial review (fresh Sonnet subagent, ARS `/ars-reviewer`): correctness, overclaiming, missing prior work.
2. Citation check (ARS `/ars-citation-check`) against `references.bib`.
3. **[Checkpoint]** Present the findings and a proposed change list.
4. Apply the approved changes; write the revised draft and the response table.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Review findings, citation report, change list | Which changes to apply |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Closed | Every finding is fixed or answered in the response table |
| Citations | Citation check reports no unresolved key |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Revised draft | `output/[paper]-revised.md` | Manuscript (read by 04) |
| Response table | `output/[paper]-review-response.md` | Finding, action, location (read by the human at the checkpoint) |

## 04-finalize

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../03-review/output/[paper]-revised.md` | Full file | final text |
| Config | `/_config/publication.md` | "Journals", "Citation style" | format |
| Guide | `/_config/setup-guides/latex.md` | Full file, first run only | building the PDF |
| Config | `/_config/bibliography/references.bib` | Full file | bibliography |

### Process

1. Convert to the journal LaTeX template (ARS `/ars-format-convert`); build the PDF if LaTeX is installed.
2. Write the AI-use disclosure (ARS `/ars-disclosure`) and the code and data availability statement (repo URL, commit SHA, manifests).
3. **[Checkpoint]** Present the PDF.
4. Commit. Submission and any push only with explicit approval.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | PDF, disclosure, availability statement | Ready to submit, or revise |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Manuscript | `output/[paper]/` | LaTeX sources, figures, .bib, PDF (final deliverable) |
