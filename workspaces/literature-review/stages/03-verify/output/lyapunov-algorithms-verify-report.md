# Verification report: lyapunov-algorithms

Stage 03, 2026-10-12. Input: included rows of `../02-search/output/lyapunov-algorithms-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetched: all; failures: 0.
- In `lyapunov-algorithms-verified.bib`: 35 entries, all class verified (32 by DOI record, 3 by arXiv record).
- Keys follow `conventions.md` (clashes suffixed a, b); none clashes with `references.bib`.
- Every fetched title was compared with its candidate row; no identifier pointed to a different paper.

## Fixes and notes

- None.

## Dropped at verification

- doi:10.1143/ptp/84.4.563: duplicate record of doi:10.1143/ptp.84.563 (same paper under another DOI form); one entry kept

## UNVERIFIED

None: every included candidate resolved to a fetched record.

## Open for stage 04

- Rows screened on title only (kept because two sources found them) still need an abstract check when the papers are read; see the reason column of the candidate list.
