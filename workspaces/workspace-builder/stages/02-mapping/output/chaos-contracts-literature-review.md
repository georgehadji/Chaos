# Contracts: literature-review workspace

DRAFT (Stage 02). One run per topic. `[topic]` is a lowercase-hyphen slug (e.g. `lyapunov-algorithms`). Paths starting with `/` are from the repo root. Full texts of papers live in `data/literature/[topic]/` (git-ignored, never committed).

## 01-scope

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | topic slug and the question it serves | per-run input |
| Config | `/_config/research-question.md` | "Primary question", "Three claims to test", "Paper sequence" | keeps the search tied to the claims |
| Config | `/_config/research-ideas.md` | "Prior art to verify before claiming any gap" | known gaps to check |
| Config | `/_config/bibliography/CONTEXT.md` | Full file | what is already covered, per section |

### Process

1. Ask for the topic slug and the question it answers.
2. Read the covered sections of the bibliography for this topic.
3. Draft the search plan: sub-questions, inclusion and exclusion rules, sources (Crossref, arXiv, ADS, INSPIRE), query strings, stop rule.
4. **[Checkpoint]** Present the plan.
5. Run the audit, revise, write the output.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Sub-questions, rules, queries, stop rule | Approve or narrow the scope |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Tied to claims | Every sub-question names a claim (C1-C3), a research line, or a catalog row ID |
| Testable rules | Each inclusion rule can be applied to a title and abstract with a yes/no answer |
| Stop rule | A numeric or source-based stop condition is stated |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Search plan | `output/[topic]-search-plan.md` | Sub-questions, rules, queries, stop rule |

## 02-search

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-scope/output/[topic]-search-plan.md` | Full file | what to search |
| Config | `/_config/bibliography/references.bib` | citation keys and DOIs only | de-duplication |

### Process

1. Run each query (Sonnet subagents, one per source, in parallel). Every candidate must carry a DOI, arXiv ID or ADS bibcode found in this run. Nothing from memory.
2. De-duplicate against existing keys and DOIs.
3. Screen titles and abstracts with the inclusion rules; record the reason for each exclusion.
4. **[Checkpoint]** Present the screened list.
5. Run the audit, write the output.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Included and excluded candidates with reasons | Add, drop, or re-screen candidates |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Identifiers | Every candidate has a DOI, arXiv ID or ADS bibcode |
| Traceability | Every candidate names the sub-question that found it |
| No duplicates | No candidate matches an existing key or DOI |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Candidate list | `output/[topic]-candidates.md` | Table: identifier, authors, year, title, sub-question, include/exclude, reason |

## 03-verify

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../02-search/output/[topic]-candidates.md` | included rows only | what to verify |
| Config | `/_config/bibliography/conventions.md` | Full file | key format, verification classes |
| Tool | `/scripts/bibcheck.py` | run, do not read | resolves identifiers on Crossref and arXiv |

### Process

1. Run `python scripts/bibcheck.py` on the candidate list (no AI). It fetches each record and reports mismatches in title, authors, venue, volume, pages.
2. Fix mismatches against the fetched record only (Sonnet subagent).
3. Classify each entry: verified, verified-secondary, or UNVERIFIED (per `conventions.md`).
4. **[Checkpoint]** Present failures and UNVERIFIED items.
5. Run the audit, write the outputs.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Mismatch fixes, failures, UNVERIFIED list | Drop, fix by hand, or keep as UNVERIFIED |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Clean bib | Every entry in the .bib is verified or verified-secondary |
| Keys | Keys follow `conventions.md` and are unique against `/_config/bibliography/references.bib` |
| Report | Every UNVERIFIED item states what failed |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Verified entries | `output/[topic]-verified.bib` | BibTeX |
| Verification report | `output/[topic]-verify-report.md` | Fixes, failures, UNVERIFIED list (read by the human at the checkpoint and by 05-promote) |

## 04-synthesize (optional, skip for small topics)

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../03-verify/output/[topic]-verified.bib` | Full file | which papers to read |
| Sources | `data/literature/[topic]/` | one paper per subagent | full texts |
| Config | `/_config/research-question.md` | "Three claims to test" | what to extract |
| Config | `/_config/catalogs/` | rows named in the search plan | where findings land |

### Process

1. Read each paper (one Sonnet subagent per paper): claims, method, benchmark values with page, equation or table number, relevance to catalog rows.
2. Merge into one notes file; list open questions and conflicts between papers.
3. **[Checkpoint]** Present the notes.
4. Run the audit, write the output.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Notes, benchmark values, conflicts | Accept, or send papers back for re-reading |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Located values | Every benchmark value has a key plus page, equation or table number |
| Own words | No quotation longer than one short sentence |
| Verified only | Every key is in `[topic]-verified.bib` or `/_config/bibliography/references.bib` |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Notes | `output/[topic]-notes.md` | Per-paper annotations, benchmark table, open questions |

## 05-promote

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 03 | `../03-verify/output/[topic]-verified.bib`, `[topic]-verify-report.md` | Full files | entries and UNVERIFIED list |
| Stage 04 | `../04-synthesize/output/[topic]-notes.md` | Full file, if it exists | annotations and benchmark values |
| Tool | `/scripts/promote.py` | run, do not read | builds the merge diff |

### Process

1. Run `python scripts/promote.py [topic]` (no AI): it proposes appends to `references.bib`, the section files, the benchmark file and the UNVERIFIED file, plus catalog `?` keys that now resolve.
2. **[Checkpoint]** Present the diff.
3. Apply the approved diff (promotion write into `/_config/`), run the audit, commit. Push only with explicit approval.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Merge diff into `/_config/bibliography/` and catalog key flips | Approve, edit, or reject |

### Audit

| Check | Pass Condition |
|-------|---------------|
| No duplicates | Every key in `references.bib` is unique |
| Size | Every bibliography section file is under 200 lines (split if not) |
| Style | No em dashes |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Updated bibliography | `/_config/bibliography/` | promotion write (final deliverable) |
| Promotion log | `output/[topic]-promote-log.md` | What was merged, the commit SHA |
