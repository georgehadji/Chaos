# Verification report: recurrence-indicators

Stage 03. Input: included rows of `../02-search/output/recurrence-indicators-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref, DataCite and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetch failures: 0.
- In `recurrence-indicators-verified.bib`: 26 entries, all class verified (21 by DOI record, 5 by arXiv record).
- Keys follow `conventions.md` (clashes suffixed a, b); one clash with `references.bib` fixed by hand (below).

## Fixes and notes

- Key clash fixed by hand: the fetched record for doi:10.1016/j.physleta.2026.132224 (Recurrence time entropy and weak chaos in Hamiltonian flows) was keyed `Sales2026`, which `references.bib` already holds for a different paper (Linear Dependence, Smaller and Generalized Alignment Indices). It is keyed `Sales2026b`, following the existing pair `Rein2019` / `Rein2019b`, so the existing key and its citations stay unchanged. `verify_report.py` did not flag this clash.

## Dropped at verification

- None.

## UNVERIFIED

None: every included candidate resolved to a fetched record.
