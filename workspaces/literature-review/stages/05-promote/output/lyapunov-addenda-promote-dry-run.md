# Promotion dry run: lyapunov-addenda

NOT APPLIED. Nothing in `/_config/` was written, nothing was committed, there is no commit SHA. This file is a dry-run diff for the owner's stage 05 checkpoint, 2026-10-13. A real promotion log (`lyapunov-addenda-promote-log.md`) exists only after approval and `--apply`.

How it was made: `scripts/promote.py lyapunov-addenda` (no `--apply`) was not run, because the permission classifier refused the call in this session (reason given: modifying shared resources, although without `--apply` the script only reads). The mechanical part below is therefore derived by hand from the same rules and read-only greps: duplicate keys, DOIs and titles against `_config/bibliography/references.bib`, `?` keys in `_config/catalogs/*.md`, the verified .bib, the UNVERIFIED section of the verify report and the benchmark table of the stage 04 notes. Re-run `python scripts/promote.py lyapunov-addenda` at the checkpoint; the expected output is given in section 1. The script only appends, so the removal (section 3) and the catalog edits (section 4) are hand edits listed here.

## 1. Mechanical merge (what promote.py should report)

- `references.bib`: +8 entries, text in `../03-verify/output/lyapunov-addenda-verified.bib`: Stewart1997, DieciVanVleck2008, VuNguyenPhan2025, Shevchenko2024, CincottaGiordanoShevchenko2022, OliveiraStewart2000, DieciElia2006, DieciEliaLopez2004. Checked by grep: none of the 8 keys is in `references.bib`; none of the 7 DOIs is; no title clashes (the ninth DOI, doi:10.1137/s0036142901392304, is DieciVanVleck2002, already present, dropped at stage 03, no new entry). Expected `BLOCKED` lines: none.
- `unverified.md`: +0 items (the verify report's UNVERIFIED section reads "None").
- `benchmarks.md`: +3 rows (section 2).
- Catalog keys with `?` that now resolve: none (no `?` suffix on any of the 8 keys in `_config/catalogs/`).

Entries by class: 7 verified (Crossref records, re-checked by `bibcheck.py` on the final .bib, 0 issues); Stewart1997 verified-secondary (no DOI; two independent registry confirmations, see the verify report), carrying `note = {Verified via secondary sources only; no DOI}` as the existing no-DOI entries do.

Owner decisions pending at this checkpoint:
- D-a. Key for item 3: `VuNguyenPhan2025` (Crossref family names Vu, Nguyen, Phan, Vietnamese order as printed) against `LinhTruongTuyen2025` (the brief's wording). Changing it means editing the key in the verified .bib before applying; nothing else depends on it yet.
- D-b. Stewart1997 accepted as verified-secondary (citable, flagged). The url in the entry is the OpenAlex record, https://openalex.org/W159031619.
- D-c. Shevchenko2024 entry carries `eprint = {2411.01939}`; the full text on disk is the arXiv v1 preprint, not the journal pagination (journal pages 622-629 are from the Crossref record).
- D-d. Whether to cite the two preprint-paginated full texts (DieciVanVleck2008, DieciVanVleck2002) by preprint page in the specs (the notes do), flagged as such.

## 2. benchmarks.md: rows to append (Geist1990 Toda chain)

Source: `data/literature/lyapunov-algorithms/Geist1990.pdf`, checked on the page images, pages as printed. Every number in the old unverified.md bullet and in the catalog D3 paragraph was checked:

| Number in the old entry | Paper | Result |
|---|---|---|
| GS 1.00 | Table III, p.886: GS 1.00 | match |
| RGS 1.02 | Table III, p.886: RGS 1.02 | match |
| H1 1.07 | Table III, p.886: H1 1.07 | match |
| H2 1.17 | Table III, p.886: H2 1.17 | match |
| continuous QR 40 | Table III, p.886: CQR 40.00; text p.886 "up to a factor 40" | match |
| "29-dim" Toda chain (corrections.md, catalog L8) | p.884, Eq. (40): N = 15, state space R^(2N-2) x S^1, [(2N-2)+1] = 29 | match |
| "discrete QR variants agree on lambda_1 and lambda_28 within finite-time fluctuation" (catalog D3) | Table II p.885, and the text beside it: differences between the discrete methods are much smaller than the fluctuations of lambda_i(n) | match; the printed table also gives the lambda_1 and lambda_28 values below |

Two qualifications, not mismatches: lambda_1 on this quasiperiodic attractor is close to zero (about -4e-5 at n = 500), so "lambda_1" is not a chaotic-exponent target; and the paper does not say which attractor the Table III "mean relative computation times" average over.

Rows (columns as in `benchmarks.md`: Quantity | Value | System / parameters | Source key(s) | How confirmed); identical to the three rows of the stage 04 notes, which promote.py harvests:

| Quantity | Value | System / parameters | Source key(s) | How confirmed |
|---|---|---|---|---|
| Mean relative CPU time of Lyapunov-exponent methods (discrete QR variants and continuous QR) | GS 1.00, RGS 1.02, H1 1.07, H2 1.17, CQR 40.00 (GS is the unit) | Damped, driven Toda chain, N = 15, state space R^(2N-2) x S^1, 29-dimensional (Eq. 40, p.884); methods of Table I, p.884 (GS Gram-Schmidt, RGS repeated GS, H1 and H2 Householder by Eqs. 17 and 18, CQR continuous QR). Reproduction target only, never a test oracle | Geist1990 | Read (p.886, Table III; the text on the same page says continuous QR needs "up to a factor 40" more CPU time). Which attractor the "mean" runs over is not stated; the CRAY X-MP/24 with 64-bit arithmetic is named for the CQR runs only (p.886) |
| Finite-time lambda_1 after n = 500 and n = 2000 driving periods, discrete QR variants | n = 500: GS -0.3878e-4, RGS -0.3878e-4, H1 -0.3880e-4, H2 -0.3875e-4 (CQR -0.2218e-4); n = 2000: GS -0.7811e-4, RGS -0.7811e-4, H1 -0.7810e-4, H2 -0.7803e-4 (CQR not run) | Same chain, three-frequency quasiperiodic T^3 solution, a = 3.5, omega_0 = 1.1237, d = 0.1 (Fig. 3 caption, p.885); lambda_1 is close to zero on this attractor. Reproduction target only, never a test oracle | Geist1990 | Read (p.885, Table II). The paper states that differences between the discrete methods are much smaller than the fluctuations of lambda_i(n) (p.885); CQR differs from the discrete values by more than they differ among themselves |
| Finite-time lambda_28 (smallest exponent) after n = 500 and n = 2000 periods | n = 500: -0.1803 for GS, RGS, H1, H2 (CQR -0.1804); n = 2000: -0.1812 for GS, RGS, H1, H2 (CQR not run) | Same chain and attractor. Reproduction target only, never a test oracle | Geist1990 | Read (p.885, Table II, columns lambda_28; the text calls it the smallest exponent). Index 28 against a 29-dimensional state space is consistent with 28 non-trivial exponents plus the time direction S^1, which the paper does not state |

## 3. unverified.md: removal

Delete line 5 of `_config/bibliography/unverified.md`, the whole bullet beginning:

`- **Geist, Parlitz & Lauterborn 1990, Toda-chain relative CPU times (GS 1.00, RGS 1.02, H1 1.07, H2 1.17, continuous QR 40) and lambda_1 values**: Bibliographic record CONFIRMED (Geist1990, ...` and ending `... and are not in the benchmark table.`

Nothing else in that file changes (17 lines remain). The other bullets are untouched.

## 4. Edits that follow from the removal (one home for each fact)

4a. `_config/catalogs/lyapunov.md`, section D3. Heading line 58 `## D3. Candidate benchmark values (user-supplied scoping review, verify before use)` becomes `## D3. Benchmark values (pointer)`. The paragraph on line 61 (the one that says "the numbers themselves are UNVERIFIED") is replaced by (text between the quote marks):

"Geist1990, damped-driven Toda chain: relative CPU times and the finite-time lambda_1 and lambda_28 values were checked against the paper (Tables II and III, pp.885-886, decision D30, O12) and live in _config/bibliography/benchmarks.md, rows with source key Geist1990. A reproduction target only, never a test oracle."

4b. `_config/catalogs/lyapunov.md`, row L8: the phrase `(Geist1990, to verify)` becomes `(Geist1990, Table III p.886)`.

4c. `_config/catalogs/lyapunov.md`, row L9 (sources column `GreeneKim1987, Dieci1997`): proposed addition of `Stewart1997` (general-n discrete product SVD, Alg. (2.1) p.31). Owner's call: it changes what the lyap-svd-stabilized spec may cite, and its P1 is the owner's open point.

4d. `_config/bibliography/lyapunov.md`, line 13 (Geist1990 entry): replace the tail `Bibliographic data confirmed; its CPU-time table could NOT be accessed (see UNVERIFIED).` with `CPU-time table (Table III, p.886) and the lambda_1, lambda_28 values (Table II, p.885) read and entered in benchmarks.md.`

4e. `_config/bibliography/lyapunov.md`, annotations to append (the file has 73 lines; +8 stays under the 200-line limit):

- **Stewart1997** D. E. Stewart (1997). A new algorithm for the SVD of a long product of matrices and the stability of products. *Electronic Transactions on Numerical Analysis* 5, 29-47. https://openalex.org/W159031619 [verified-secondary] -- Discrete product SVD of p square factors without forming the product (LRCH by treppen-iteration, Alg. (2.1) p.31, O(n^3 p) per iteration, backward error Theorem 1 pp.35-36); needs all factors stored; high relative accuracy only for stable products (Def. 1, p.38). Full text on disk.
- **OliveiraStewart2000** S. Oliveira & D. E. Stewart (2000). Exponential splittings of products of matrices and accurately computing singular values of long products. *Linear Algebra and its Applications* 309(1-3), 175-190. doi:10.1016/s0024-3795(99)00273-6 -- Companion to Stewart1997 on accurate singular values of long products. Wanted; no full text, not read.
- **DieciVanVleck2008** L. Dieci & E. S. Van Vleck (2008). On the error in QR integration. *SIAM Journal on Numerical Analysis* 46(3), 1166-1189. doi:10.1137/06067818x -- Global error bounds for the computed Q factor independent of the interval length, under integral separation (Thm 3.2, Thm 4.2); a flow set-up, no map algorithm. Full text is the author preprint.
- **VuNguyenPhan2025** V. H. Linh, N. D. Truong & P. Q. Tuyen as printed (Crossref: Vu, Nguyen, Phan) (2025). On computation of Lyapunov exponents by QR methods with error control for semi-linear DAEs. *DAE Panel* 3. doi:10.52825/dae-p.v3i.2632 -- Discrete (Alg. 3.1) and continuous (Alg. 3.3) QR with embedded Runge-Kutta error control for DAEs; no map algorithm.
- **DieciElia2006** L. Dieci & C. Elia (2006). The singular value decomposition to approximate spectra of dynamical systems. Theoretical aspects. *Journal of Differential Equations* 230(2), 502-531. doi:10.1016/j.jde.2006.08.007 -- Justifies SVD methods for systems with stable, distinct exponents; continuous SVD ODEs only; refers to Stewart1997 for the discrete algorithm. Full text on disk.
- **DieciEliaLopez2004** L. Dieci, C. Elia & L. Lopez (2004). Smooth SVD on the Lorentz group with application to computation of Lyapunov exponents. *Journal of Computational and Applied Mathematics* 164-165, 255-264. doi:10.1016/s0377-0427(03)00644-7 -- Continuous smooth SVD for Lorentz-group fundamental matrices; not for symplectic Jacobians.
- **Shevchenko2024** I. I. Shevchenko (2024). Lyapunov exponents in fundamental models of nonlinear resonance. *JETP Letters* 120(8), 622-629. doi:10.1134/s0021364024602872, arXiv:2411.01939 -- Analytical Lyapunov-time estimates from separatrix-map theory for the pendulum and Andoyer models (kept for the later papers; not annotated further). Full text on disk is arXiv v1.
- **CincottaGiordanoShevchenko2022** P. M. Cincotta, C. M. Giordano & I. I. Shevchenko (2022). Revisiting the relation between the Lyapunov time and the instability time. *Physica D* 430, 133101. doi:10.1016/j.physd.2021.133101 -- Lyapunov time against instability time; kept for the later papers; no full text, not read.

4f. After 4a to 4e the statements "UNVERIFIED" about Geist1990 in `_config/` are gone; the grep that finds them is `grep -n -i geist _config/catalogs/lyapunov.md _config/bibliography/*.md`.

## 5. Files outside `/_config/` changed in this run (already written, nothing to apply)

- Stage outputs under `workspaces/literature-review/stages/01..05/output/lyapunov-addenda-*`.
- `data/literature/lyapunov-addenda/` (full texts under their keys, draft.bib, bibcheck.log, to-verify.md, verify-overrides.json) and `data/literature/lyapunov-algorithms/DieciVanVleck2002.pdf` and `.txt`. `data/` is git-ignored.
- Not committed.
