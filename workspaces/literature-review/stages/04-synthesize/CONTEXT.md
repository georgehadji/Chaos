# Stage 04: Synthesize

Read the verified papers and extract claims and located benchmark values. Optional: skip for small topics.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-scope/output/[topic]-search-plan.md` | sub-questions, catalog rows named | what to extract |
| Previous stage | `../03-verify/output/[topic]-verified.bib` | Full file | which papers to read |
| Sources | `data/literature/[topic]/` | one paper per subagent | full texts |
| Config | `/_config/research-question.md` | "Three claims to test" | what to extract |
| Config | `/_config/catalogs/` | rows named in the search plan | where findings land |

## Process

1. Fetch open-access full texts into `data/literature/[topic]/`; list paywalled papers for the human to add. Then read each paper (one Sonnet subagent per paper): claims, method, benchmark values with page, equation or table number, relevance to catalog rows.
2. Merge into one notes file; list open questions and conflicts between papers.
3. **[Checkpoint]** Present the notes.
4. Run the audit, write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Notes, benchmark values, conflicts | Accept, or send papers back for re-reading |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Located values | Every benchmark value has a key plus page, equation or table number |
| Own words | No quotation longer than one short sentence |
| Verified only | Every key is in `[topic]-verified.bib` or `/_config/bibliography/references.bib` |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Notes | `output/[topic]-notes.md` | Per-paper annotations, benchmark table, open questions |
