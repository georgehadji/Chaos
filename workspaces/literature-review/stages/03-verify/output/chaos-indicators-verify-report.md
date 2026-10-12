# Verification report: chaos-indicators

Stage 03, 2026-10-12. Input: included rows of `../02-search/output/chaos-indicators-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetched: all; failures: 0.
- In `chaos-indicators-verified.bib`: 35 entries, all class verified (29 by DOI record, 6 by arXiv record).
- Keys follow `conventions.md` (clashes suffixed a, b); none clashes with `references.bib`.
- Every fetched title was compared with its candidate row; no identifier pointed to a different paper.

## Fixes and notes

- None.

## Dropped at verification

- SkokosManos2014: same title as SkokosManos2016 in references.bib (different identifier form)

## UNVERIFIED

None: every included candidate resolved to a fetched record.

## Open for stage 04

- Rows screened on title only (kept because two sources found them) still need an abstract check when the papers are read; see the reason column of the candidate list.
