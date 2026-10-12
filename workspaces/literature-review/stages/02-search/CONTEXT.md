# Stage 02: Search

Run the search plan and screen candidates by their identifiers.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-scope/output/[topic]-search-plan.md` | Full file | what to search |
| Config | `/_config/bibliography/references.bib` | citation keys and DOIs only | de-duplication |
| Tool | `/scripts/litmerge.py` | run, do not read | merges the per-source tables, drops papers already in the bibliography |

## Process

1. Run each query (Sonnet subagents, one per source, in parallel). Every candidate must carry a DOI, arXiv ID or ADS bibcode found in this run. Nothing from memory.
2. De-duplicate across sources and against existing keys and DOIs: `python scripts/litmerge.py [topic]` (no AI).
3. Screen titles and abstracts with the inclusion rules; record the reason for each exclusion.
4. **[Checkpoint]** Present the screened list.
5. Run the audit, write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Included and excluded candidates with reasons | Add, drop, or re-screen candidates |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Identifiers | Every candidate has a DOI, arXiv ID or ADS bibcode |
| Traceability | Every candidate names the sub-question that found it |
| No duplicates | No candidate matches an existing key or DOI |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Candidate list | `output/[topic]-candidates.md` | Table: identifier, authors, year, title, sub-question, include/exclude, reason |
