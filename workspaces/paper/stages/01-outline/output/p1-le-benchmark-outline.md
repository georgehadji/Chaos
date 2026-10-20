# Outline: paper p1-le-benchmark

- Slug: p1-le-benchmark
- Status: approved
- Journal: to be decided
- Studies used: p1-le-benchmark, accepted results document `workspaces/study/stages/05-report/output/p1-le-benchmark-results.md`
- Scope source: `_config/research-question.md`, "Paper sequence", row P1. Claims tested: C1 (asymptotic equivalence) and C2 (finite-time inequivalence).
- Framing (human decision): report as preregistered. No new claims, no softening, no reinterpretation.

Pointer scheme. "RD" means the results document. Pointers inside square brackets such as [T2, header line], [params: h1c_fraction], [findings B, H1], [manifest: status] are the results document's own pointers and are carried over unchanged. Citation keys are written `[@Key]`.

## Contribution

This paper reports one preregistered numerical study of largest-Lyapunov-exponent estimators on the standard map, run on a panel of 883 orbits at five nonzero kick values plus exact-exponent controls. It delivers FALSIFIES: the preregistered rule for C1 was not met (H1 falsified), a test that could have failed and did. It delivers NEGATIVE: H2 (finite-time inequivalence in weak chaos) is not supported, and is reported with the same standing as the supported results. It delivers REPRODUCES as a failed attempt: the reference value of Sala and Artuso at K = pi/2 did not reproduce, while the K = 2 pi value did. The convergence-time regime dependence (H3) and the label informativeness (H5) are supported. TOOL is claimed only in a limited sense: the run used package version 0.0.1, and the tagged release v0.1.0 has not been cut. NOVEL is not claimed.

## Sections

### 1. Introduction (literature only)

- Claim 1.1: The largest-Lyapunov-exponent values of Sala and Artuso at K = pi/2 (0.298) and K = 2 pi (1.172) are the published references against which the study attempts a reproduction.
  - Evidence: [@SalaArtuso2015]; RD, H4 preregistered rule ([params: h4_ref_halfpi, h4_ref_2pi]).
- Claim 1.2: Two published comparators for the chaotic-sea exponent are the fit of Shevchenko and the line ln(K/2) of Chirikov, and the study uses both as comparators that never pass or fail.
  - Evidence: [@Shevchenko2014], [@Chirikov1979]; RD, H4, "Comparators (never pass or fail)".
- Claims C1 and C2 as defined in `_config/research-question.md` are stated here as the questions the study tests. No other literature claim is made in this section.
  - Evidence: `_config/research-question.md` (claims C1 and C2 rows); RD, "Claims C1 and C2".

### 2. Methods and preregistration

- Claim 2.1: The decision rules were frozen on 2026-10-17, before any data existed, with Amendments A1 to A6 all frozen before the first study data file, and the rules are applied without change.
  - Evidence: RD, Header, "What the study tested" [findings A]; preregistration header SHA-256 [manifest: prereg.header].
- Claim 2.2: The panel has 883 orbits, one of which (orbit 881) is dropped as drift, with design classes T, R, C, S and X, and the sets used by the decision rules are E_ch 209, W 355, S 80, S_L 92 (H5b 883, H5c 334).
  - Evidence: RD, Header, "Reading key" and "Panel and sets" [T1], [findings E, reading 8], [T7, rows H5a, H5b, H5c]. Appendix pointer: RD, Deviations and readings, readings 5 and 8.
- Claim 2.3: All reported numbers come from one fully automatic analysis run over the stored series, identified by git commit, module version 0.0.1, panel SHA-256, preregistration header SHA-256, manifest SHA-256, parameter file SHA-256 and bootstrap seed 20261012. The run completed with 1363 chunks, all with status ok.
  - Evidence: RD, Header run-identity table ([findings A, table], [manifest: git.sha], [manifest: versions.module], [manifest: panel.sha256_end], [manifest: prereg.header], [manifest: status], [manifest: plan_counts.chunks]).
- Claim 2.4: The frozen text was applied with the deviations the program states: wrapped-coordinate L2 (A5 item 1), orthonormalisation breakdown recorded as a result in the tau sweep (A5 item 2), control bound scope (A5 item 3), L9 control bound 5.0e-8 and 1.0e-8 for the other methods (A6), and no manifest overrides.
  - Evidence: RD, "Deviations and readings", deviations 1 to 5 ([findings E, Deviations], [T Deviations], [params: control_bound, control_bound_L9], [manifest: overrides]).

### 3. Results: H1 (claim C1), FALSIFIED

- Claim 3.1: H1 is falsified. On E_ch (209 eligible orbits) the pass fractions are 1.0 for H1a (209 of 209) and 1.0 for H1b (209 of 209), and 0.775119617225 for H1c (162 of 209) against the limit of 0.90.
  - Evidence: RD, H1 measured table ([T2, header line; params: h1a_fraction, h1b_fraction, h1c_fraction]).
- Claim 3.2: The shortfall is carried by H1c alone. The 47 failing orbits are in three K cells (K index 2: 27 of 40 pass; K index 3: 13 of 40 pass; K index 4: 33 of 40 pass), K index 5 and 6 pass 40 of 40, and the 9 controls pass.
  - Evidence: RD, H1, "H1c passing orbits by cell" ([T2, rows ...]). Failing-orbit ids and extremes: [T2, Failing orbits].
- Claim 3.3: On the exact-exponent controls under the bound (871, 873, 875, 877, 879, 883), the largest deviations are within the bounds of A6 (L9 control bound 5.0e-8), with 877 and 879 above 1e-8 and below 5.0e-8.
  - Evidence: RD, H1, "Exact-exponent controls" table ([T2, controls, ...], [params: control_bound, control_bound_L9]).
- Reading under FALSIFIES: the preregistered rule for C1 is not met and H1 is reported as falsified.
  - Evidence: RD, H1, final paragraph [findings B, H1].
- EXPLORATORY: D1 (where the H1c failures sit) and D2 (the info controls 878, 880, 881, 882 and their deviations). Evidence: RD, Exploratory results D1, D2. No cause is claimed.

### 4. Results: H2 (claim C2), NOT SUPPORTED (NEGATIVE)

- Claim 4.1: H2 is not supported. Over the weak-chaos set W (355 orbits) the median of R(n_H2) is 2.23588035691e-5 with a 95 percent bootstrap interval of [1.9257235885e-5, 2.31226443693e-5], against the preregistered thresholds of 0.10 for the median and 0.05 for the interval lower end. Both conditions fail.
  - Evidence: RD, H2 ([T3, header line; params: h2_median, h2_interval]).
- Claim 4.2: The result is a null (NEGATIVE) with the same standing as a pass. The frozen rule decides on the pooled W, and the per-cell table is descriptive. Class S at K index 5 has 1 orbit and is not decidable.
  - Evidence: RD, H2, cell table and last paragraph ([T3, rows ...], [findings B, H2], [params: min_set]).

### 5. Results: H3 (claim C2), SUPPORTED

- Claim 5.1: H3 is supported. The median n_conv is 2.511886e6 over W (355 orbits) and 224.0 over the strong-chaos set S (80 orbits), a ratio of 11213.7767857 with a 95 percent bootstrap interval of [7924.465, 14117.3125], against the thresholds of 10 and 3. The share of undefined n_conv is 0.0229885057471.
  - Evidence: RD, H3 ([T4, header line; params: h3_ratio, h3_interval]).
- Claim 5.2: The strong-chaos set S of H3 is the class C orbits at K = 2 pi and K = 10 and is not the panel class S. Undefined n_conv orbits enter the medians at the value N.
  - Evidence: RD, Header "Reading key" and H3 measured paragraph ([findings A], [findings E, reading 8], [findings B, H3]).
- Cell-level medians are descriptive (RD, H3 cell table, [T4, rows ...]).

### 6. Results: H4 (REPRODUCES, [@SalaArtuso2015]), DID NOT REPRODUCE

- Claim 6.1: H4 did not reproduce. At K = pi/2 both L1 and L5 give 0.256066843848 against the reference 0.298 with tolerance 0.01 (fail). At K = 2 pi both give 1.17018723872 against 1.172 with tolerance 0.02 (pass). The reproduction as a whole fails.
  - Evidence: RD, H4 measured table ([T5, row K = pi/2], [T5, row K = 2 pi], [T5, verdict line]); reference [@SalaArtuso2015].
- Claim 6.2: The possibility that the published orbit sits in a different component remains open and the tables do not decide it. Each K rests on one start with no error bar.
  - Evidence: RD, H4 final paragraph ([findings B, H4]).
- Claim 6.3 (comparators, no verdict): the ensemble mean of lambda_L5(N_1) over 40 kept chaotic-sea orbits, the Shevchenko fit and ln(K/2) are reported at the five K values.
  - Evidence: RD, H4 comparators table ([T5, comparators, K index 2 to 6]); [@Shevchenko2014], [@Chirikov1979].
- EXPLORATORY: D4 (H4 at the shorter horizon T = 1e8; cross-check values 0.299256534631 Float64 and 0.299012763029 Double64 at K = pi/2). Evidence: RD, D4 and H4 cross-check table ([T5, cross-check, ...]). No cause is claimed.

### 7. Results: H5, SUPPORTED

- Claim 7.1: H5 is supported. H5a is 1.0 on 92 ladder-subset orbits (limit at least 0.90), H5b is 0.0951302378256 on 883 orbits (limit at most 0.30), and H5c is 0.0568862275449 on 334 orbits (limit at most 0.30).
  - Evidence: RD, H5 table ([T7, rows H5a, H5b, H5c; params: h5a_min, h5b_max, h5c_max]).
- Claim 7.2: H5a is read from the per-precision labels before the ladder test acts, and H5 is supported only if all three sub-rules pass.
  - Evidence: RD, H5 and Deviations and readings ([findings E, reading 5]).
- EXPLORATORY: D3 (H5a is 1.0 while the ladder test forced 35 of 92 labels to unresolved, share 0.380434782609) and D5 (unresolved labels by cell: 65 with both axes unresolved, 19 structure only). Evidence: RD, D3, D5 ([F4, row ladder], [T6, forced table], [T6, label table]).

### 8. Cost and sensitivity (reported without a decision rule)

- Claim 8.1: Relative to L5 in Float64 (168.60570837 ns per orbit-step), L1 costs 0.701549448499, L2 0.503966109187, L4 0.982175001597, L9 1.23876330018 and L6 2.10273197634. L5 costs 4882.6924383 ns in Double64 and 377399.856865 ns in BigFloat(256). Each row is an unweighted mean over a different item set.
  - Evidence: RD, Cost table ([F6, rows ...], [findings C, Cost], [findings E, reading 12]).
- Claim 8.2: Run totals: total wall time 11957.993616580963 seconds against a cap of 12.0 hours, 1363 chunks, 15569 items, 23 threads of 23 available, and all 13 arms complete.
  - Evidence: RD, Cost, run totals ([manifest: wall.total_seconds, wall.cap_hours, plan_counts.chunks, plan_counts.items, threads.used, threads.available, arms.*.status]).
- Claim 8.3: The sensitivity table (31 rows) shows direction, eps, wrapped-versus-lifted and tau-sweep effects. Wrapped versus lifted is 0.0 for all eight methods of V, and the tau sweep shows breakdown from tau = 8 (17 orbits for L4 and L5) and tau = 16 (39 for L4, 37 for L5).
  - Evidence: RD, Sensitivity ([F7, rows ...], [findings C, Sensitivity]).
- Symplectic residuals (figure F5) are shown as a figure only. No numbers are quoted.
  - Evidence: RD, Sensitivity last paragraphs ([findings C]).

### 9. Exploratory results

- All items are EXPLORATORY, not confirmatory, with no frozen rule and no cause claimed. Items D1 to D7 are placed in Sections 3, 6 and 7 where they bear on a hypothesis (D1, D2, D3, D4, D5). D6 (short, empty and not decidable cells) and D7 (cell patterns of R(n_H2) and n_conv) are reported here.
  - Evidence: RD, Exploratory results ([findings D], D6 [T1, manifest cells table], D7 [T3, rows S 3, X 2] and [T4, rows ...]).
- Claim 9.1 (EXPLORATORY): The largest cell median of R(n_H2) is 0.000775501992792 (class S, K index 3) and the smallest is 5.82402750446e-12 (class X, K index 2).
  - Evidence: RD, D7 ([T3, row S 3], [T3, row X 2]).
- Claim 9.2 (EXPLORATORY): Class X shows 0 qualifying grid cells at every K index from 2 to 6 because the preregistered selection rule does not draw class X from the screening.
  - Evidence: RD, D6 ([T1, manifest cells, rows 2 X to 6 X], [findings D6]).

### 10. Limitations

- Claim 10.1: The study rests on one panel draw (a single seed, 883 orbits), so the bootstrap intervals do not reflect a second draw.
  - Evidence: RD, Limitations bullet 1 [T1].
- Claim 10.2: The panel is run in Float64 only, and precision comparisons rest on the ladder subset S_L of 92 orbits.
  - Evidence: RD, Limitations bullet 2 [T7, row H5a].
- Claim 10.3: Cells with fewer than 20 orbits are not decidable (class S at K index 5: 5 drawn, 1 in W) or descriptive (class S at K index 2: 25 drawn, 21 in W). H4 rests on one start per K with no error bar.
  - Evidence: RD, Limitations bullets 5 and 6 ([T1, manifest cells, rows 5 S, 2 S], [T3, row S 5], [T4, row S 5], [findings B, H4]).
- Claim 10.4: Wrapped-coordinate L2 affects the H1c metric, and the tables do not show whether this accounts for any part of the H1c result.
  - Evidence: RD, Deviations and readings, deviation 1 ([findings E, deviation 1], [T Deviations]).
- Claim 10.5: The numbers were checked against the tables file of the run and the analysis program was verified against an independent recomputation (`analysis-verification.md`). The numbers were not recomputed from the stored series for the results document. Figures F1 to F7 were checked only by file size and the `%PDF-` header.
  - Evidence: RD, Limitations last bullet and open-items list ([findings F]).

### 11. Data and code availability

- Claim 11.1: Repository https://github.com/georgehadji/Chaos (GPL-3.0-or-later). The tagged release v0.1.0 is NOT yet cut (cut at milestone M3 exit). The run used package version 0.0.1 at the commit given in the run-identity table, before the tag.
  - Evidence: `_config/publication.md` (repository); RD, Deviations and readings, "Package version" ([findings E, Note], [manifest: versions.module]).
- Claim 11.2: The run manifest and the stored series identify the run, with the SHA-256 values in the run-identity table.
  - Evidence: RD, Header run-identity table.

### 12. Conclusion

- Claim 12.1: Under the preregistered rules, H1 (claim C1) is falsified, H2 (claim C2) is not supported (NEGATIVE), H3 (claim C2) is supported, H4 did not reproduce the published value at K = pi/2, and H5 is supported.
  - Evidence: RD, "Claims C1 and C2" ([findings B, H1, H2, H3]), H4 and H5 sections. H4 and H5 are not tied to C1 or C2 in the findings file.

## Evidence map

| Claim id | Section | Evidence | Status |
|---|---|---|---|
| 1.1 | 1 | [@SalaArtuso2015]; RD H4 [params: h4_ref_halfpi, h4_ref_2pi] | citation + traced |
| 1.2 | 1 | [@Shevchenko2014], [@Chirikov1979]; RD H4 comparators | citation + traced |
| C1, C2 definitions | 1 | `_config/research-question.md` rows C1, C2; RD "Claims C1 and C2" | traced |
| 2.1 | 2 | RD Header [findings A]; [manifest: prereg.header] | traced |
| 2.2 | 2 | RD Header [T1], [findings E, reading 8], [T7] | traced |
| 2.3 | 2 | RD Header run-identity table | traced |
| 2.4 | 2 | RD Deviations 1 to 5 | traced |
| 3.1 | 3 | RD H1 [T2, header line] | traced |
| 3.2 | 3 | RD H1 cell table [T2, rows ...] | traced |
| 3.3 | 3 | RD H1 controls table [T2, controls, ...] | traced |
| 3.D1, 3.D2 | 3 | RD D1, D2 (EXPLORATORY) | traced |
| 4.1 | 4 | RD H2 [T3, header line] | traced |
| 4.2 | 4 | RD H2 [T3, rows ...], [findings B, H2] | traced |
| 5.1 | 5 | RD H3 [T4, header line] | traced |
| 5.2 | 5 | RD Header reading key, [findings E, reading 8] | traced |
| 6.1 | 6 | RD H4 [T5, rows K = pi/2, K = 2 pi]; [@SalaArtuso2015] | citation + traced |
| 6.2 | 6 | RD H4 [findings B, H4] | traced |
| 6.3 | 6 | RD H4 comparators [T5]; [@Shevchenko2014], [@Chirikov1979] | citation + traced |
| 6.D4 | 6 | RD D4 (EXPLORATORY) | traced |
| 7.1 | 7 | RD H5 [T7, rows H5a, H5b, H5c] | traced |
| 7.2 | 7 | RD [findings E, reading 5] | traced |
| 7.D3, 7.D5 | 7 | RD D3, D5 (EXPLORATORY) | traced |
| 8.1 | 8 | RD Cost [F6] | traced |
| 8.2 | 8 | RD Cost run totals [manifest: ...] | traced |
| 8.3 | 8 | RD Sensitivity [F7] | traced |
| 9.1 | 9 | RD D7 [T3, rows S 3, X 2] (EXPLORATORY) | traced |
| 9.2 | 9 | RD D6 [T1, manifest cells] (EXPLORATORY) | traced |
| 10.1 to 10.5 | 10 | RD Limitations, Deviations 1 | traced |
| 11.1 | 11 | `_config/publication.md`; RD "Package version" | traced |
| 11.2 | 11 | RD Header run-identity table | traced |
| 12.1 | 12 | RD "Claims C1 and C2", H1 to H5 | traced |

No claim has status "uncited".

## Unmapped claims

None in the outline. The following are held out because the results document gives no citation key:
- The definition of K_c = 0.971635406 (used as given in [T5, comparators, K index 3; params: K_c in Header]).
- The definitions of the methods L1 to L9 and Q.
These may appear in the Methods text only as "defined in the study design" with a pointer to the repository, not as cited statements.

## Open items for the human

1. Target journal: to be decided; `_config/publication.md` defers the decision to this stage.
2. Tagged release v0.1.0 is not cut (milestone M3 exit); the run used package version 0.0.1, so TOOL cannot be claimed as delivered.
3. The results header lists two git SHAs (d31b655d0c33e2f2f1a67bbbdeba9c8c8a3fd036 from the tables header and findings, 3e405f3324f8d25162580d3380db26c877d7ca63 from [manifest: git.sha]). The study does not relate them.
4. K_c and the L1 to L9 and Q definitions have no citation key in `_config/bibliography/references.bib` as used by the results document.
5. The STRONG_K constant was kept by owner decision.
6. Manifest.toml for Julia 1.10 is unresolved.
7. Figures F1 to F7 were checked only by file size and the `%PDF-` header, not visually.
