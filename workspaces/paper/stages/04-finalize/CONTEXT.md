# Stage 04: Finalize

Format for the journal, add disclosures, build the PDF.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../03-review/output/[paper]-revised.md` | Full file | final text |
| Config | `/_config/publication.md` | "Journals", "Citation style" | format |
| Guide | `/_config/setup-guides/latex.md` | Full file, first run only | building the PDF |
| Config | `/_config/bibliography/references.bib` | Full file | bibliography |

## Process

1. Convert to the journal LaTeX template (ARS `/ars-format-convert`); build the PDF if LaTeX is installed.
2. Write the AI-use disclosure (ARS `/ars-disclosure`) and the code and data availability statement (repo URL, commit SHA, manifests).
3. **[Checkpoint]** Present the PDF.
4. Commit and push. Submission only with explicit approval.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | PDF, disclosure, availability statement | Ready to submit, or revise |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Manuscript | `output/[paper]/` | LaTeX sources, figures, .bib, PDF (final deliverable) |
