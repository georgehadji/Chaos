# Verification report: lyapunov-addenda

Stage 03. Input: included rows of `../02-search/output/lyapunov-addenda-candidates.md`, one identifier per paper (DOI preferred, else arXiv ID). Tool: `scripts/bibcheck.py` (Crossref, DataCite and arXiv records); every BibTeX field comes from the fetched record (ALL-CAPS names title-cased).

## Summary

- Fetch failures: 0.
- In `lyapunov-addenda-verified.bib`: 8 entries: 7 verified (Crossref DOI records), 1 verified-secondary (Stewart1997, no DOI; see Fixes and notes). Corrected by hand: the script counted the secondary entry as an arXiv record.
- Keys follow `conventions.md`; none clashes with `references.bib` (checked by script). The owner's nine items give 8 new entries; item 9 already exists as DieciVanVleck2002.
- `bibcheck.py` on the final .bib: re-run clean (see command log in the run summary).
- Input exception: Stewart1997 has no identifier, so `litmerge.py --to-verify` cannot list it (it fails with IndexError on a candidate row without DOI or arXiv ID); the other 8 rows were run through the script on a scratch copy and Stewart1997 was added by hand.

## Fixes and notes

- Stewart1997: No DOI or arXiv ID exists (ETNA volume 5 is not in Crossref). Class verified-secondary: (a) Crossref reference list of Oliveira and Stewart 2000 (doi:10.1016/s0024-3795(99)00273-6, entry BIB17) cites Electron. Trans. Numer. Anal. 5 (1997) 29-47 with this title (minus 'and the stability of products'); (b) OpenAlex record W159031619 gives the full title, author David E. Stewart, year 1997; (c) page 1 of the full text carries ETNA Volume 5, pp. 29-47, June 1997. ETNA's own pages were not reachable (host failed, mirror behind a bot challenge).
- DieciVanVleck2008: Record: SIAM J. Numer. Anal. 46(3), 1166-1189 (2008), found by Crossref title search. The supplied qr-error.pdf is the author preprint (20 pages, no SIAM header; title and abstract equal the Crossref record), so page numbers cited from the full text are preprint pages, not journal pages.
- VuNguyenPhan2025: Key follows the Crossref record, which gives family names Vu, Nguyen, Phan (Vietnamese order, matching the printed byline); the owner's brief wrote Linh, Truong, Tuyen. Owner may choose LinhTruongTuyen2025 instead; pending decision at stage 05. Crossref gives no pages or article number. DOI equals the one printed on page 1.
- Shevchenko2024: Record: JETP Letters 120(8), 622-629 (2024). The supplied full text is arXiv:2411.01939v1 (the arXiv record lists the same DOI and title); it is not the journal pagination. eprint added by hand because bibcheck skipped the arXiv ID as a duplicate of the DOI.
- DieciElia2006: First page of the full text reads J. Differential Equations 230 (2006) 502-531, matching the record.
- DieciEliaLopez2004: First page of the full text reads J. Comput. Appl. Math. 164-165 (2004) 255-264, matching the record.
- CincottaGiordanoShevchenko2022: Record only, no full text supplied.
- OliveiraStewart2000: Record only, wanted, no full text supplied. DOI found by Crossref title search in this run.

## Full-text first-page checks

| Item | File | First page against the record |
|---|---|---|
| 1 Stewart1997 | svd-long-product.pdf | ETNA 5, pp. 29-47, June 1997; title, author, abstract (Cholesky LR, hyperbolic Givens, stable products) match |
| 2 DieciVanVleck2008 | qr-error.pdf | title, authors and abstract equal the Crossref record; author preprint (20 pp.) |
| 3 VuNguyenPhan2025 | qr-methods.pdf | printed DOI 10.52825/dae-p.v3i.2632, published 03 Sep 2025, title and three authors match |
| 4 Shevchenko2024 | Shevchenko2024-jetpl.pdf | arXiv:2411.01939v1 (4 Nov 2024), title and author match; preprint, not the JETPL typesetting |
| 7 DieciElia2006 | svd-approximate.pdf | J. Differential Equations 230 (2006) 502-531, title, authors match |
| 8 DieciEliaLopez2004 | smooth-svd-lorentz.pdf | JCAM 164-165 (2004) 255-264, title, authors match |
| 9 DieciVanVleck2002 | spectral-intervals.pdf | title and authors match the existing record; abstract wording differs slightly from Crossref ("during the years" against "over the years"), so it is an author preprint (24 pp. against 27 printed pages 516-542); no SIAM header |
| 5, 6 | none | record only |

## Dropped at verification

- DieciVanVleck2002 (item 9, doi:10.1137/s0036142901392304): the fetched Crossref record equals the existing references.bib entry (same title, volume 40(2), pages 516-542); no new entry. Full text filed under `data/literature/lyapunov-algorithms/DieciVanVleck2002.pdf`.
- arxiv:2411.01939: duplicate record of doi:10.1134/s0021364024602872 (same paper under another DOI form); one entry kept

## UNVERIFIED

None: every included candidate resolved to a fetched record or, for Stewart1997, to two independent secondary records.
