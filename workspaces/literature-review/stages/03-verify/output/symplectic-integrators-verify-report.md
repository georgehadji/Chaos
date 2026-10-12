# Verification report: symplectic-integrators

Stage 03, 2026-10-12. Input: included rows of `../02-search/output/symplectic-integrators-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetched: all; failures: 0.
- In `symplectic-integrators-verified.bib`: 55 entries, all class verified (49 by DOI record, 6 by arXiv record).
- Keys follow `conventions.md` (clashes suffixed a, b); none clashes with `references.bib`.
- Every fetched title was compared with its candidate row; no identifier pointed to a different paper.

## Fixes and notes

- None.

## Dropped at verification

- FarandaMestreTurchetti2012: same title as FarandaMestreTurchetti2012 in references.bib (different identifier form)
- Anon2005 (doi:10.1017/cbo9780511614118.006): chapter of LeimkuhlerReich2004, already in references.bib; cite the book with the chapter title ('The modified equations')
- GerlachSkokos2010: same title as GerlachSkokos2011 in references.bib (different identifier form)
- Danieli2018: same title as Danieli2019 in references.bib (different identifier form)
- arxiv:1702.03354: duplicate record of doi:10.1007/s11075-017-0287-z (same paper under another DOI form); one entry kept

## UNVERIFIED

None: every included candidate resolved to a fetched record.

## Open for stage 04

- Rows screened on title only (kept because two sources found them) still need an abstract check when the papers are read; see the reason column of the candidate list.
