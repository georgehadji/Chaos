# Verification report: standard-map-stickiness

Stage 03, 2026-10-12. Input: included rows of `../02-search/output/standard-map-stickiness-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetched: all; failures: 0.
- In `standard-map-stickiness-verified.bib`: 49 entries, all class verified (48 by DOI record, 1 by arXiv record).
- Keys follow `conventions.md` (clashes suffixed a, b); none clashes with `references.bib`.
- Every fetched title was compared with its candidate row; no identifier pointed to a different paper.

## Fixes and notes

- HomanMeiss2026: published title differs from the arXiv title in the candidate row ('Recurrence and Stickiness in the Noisy Harper Map'); the entry uses the record fetched by DOI.
- EfthymioulosContopoulosVoglis1999: the Crossref record spells the first author 'Efthymioulos'; the same author's 1997 record in this list (Efthymiopoulos1997) spells 'Efthymiopoulos'. Kept as the record gives it; fix by hand at the checkpoint if you agree it is a record typo.

## Dropped at verification

- HarsoulaKaramanosContopoulos2019: also found by lyapunov-algorithms; verified there, kept once
- TomsovicLakshminarayan2007: also found by lyapunov-algorithms; verified there, kept once
- daSilvaMancheinBeims2019: also found by lyapunov-algorithms; verified there, kept once
- Kai2026 (doi:10.21203/rs.3.rs-9004238/v1): Research Square preprint: neither peer-reviewed nor arXiv, fails common inclusion rule (a). Dropped, not UNVERIFIED.
- MancheinBeimsRost2009: also found by lyapunov-algorithms; verified there, kept once

## UNVERIFIED

None: every included candidate resolved to a fetched record.

## Open for stage 04

- Rows screened on title only (kept because two sources found them) still need an abstract check when the papers are read; see the reason column of the candidate list.
