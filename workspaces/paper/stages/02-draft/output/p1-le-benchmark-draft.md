# Largest Lyapunov exponent estimators on the standard map: a preregistered benchmark

## Abstract

We report a preregistered numerical study of largest-Lyapunov-exponent estimators on the standard map. <!-- RD: Header, "What the study tested" -->
The panel had 883 orbits at five nonzero kick values, plus exact-exponent controls. <!-- RD: Header, "What the study tested", [T1] -->
The decision rules were frozen before any data existed. <!-- RD: Header, "What the study tested" [findings A] -->
The claim of asymptotic equivalence (H1) was falsified. On 209 eligible orbits the pass fraction for sub-rule H1c was 0.775119617225 against a limit of 0.90. <!-- RD: H1, [T2, header line; params: h1c_fraction] -->
The claim of finite-time inequivalence in weak chaos (H2) was not supported. Over 355 orbits the median of R(n_H2) was 2.23588035691e-5 against a threshold of 0.10. <!-- RD: H2, [T3, header line; params: h2_median] -->
The regime dependence of the convergence time (H3) was supported, with a ratio of medians of 11213.7767857. <!-- RD: H3, [T4, header line] -->
The reproduction of the reference values of Sala and Artuso [@SalaArtuso2015] (H4) failed at K = pi/2, where we obtained 0.256066843848 against 0.298. It passed at K = 2 pi, where we obtained 1.17018723872 against 1.172. <!-- RD: H4, [T5, row K = pi/2], [T5, row K = 2 pi] -->
The labels were informative (H5). <!-- RD: H5, [T7, verdict line] -->
Cost and sensitivity are reported without a decision rule. <!-- RD: Reported without a decision rule -->
The run used package version 0.0.1. The tagged release v0.1.0 has not been cut. <!-- RD: Deviations and readings, "Package version" [findings E, Note; manifest: versions.module] -->

## 1. Introduction

The largest Lyapunov exponent of the standard map is estimated by several methods. Two questions concern how far these estimators agree. These are claim C1 (asymptotic equivalence) and claim C2 (finite-time inequivalence). <!-- research-question.md: rows C1, C2; RD: "Claims C1 and C2" -->

Claim C1 states that the Benettin, QR (Gram-Schmidt or Householder), stabilized SVD and well-conditioned tangent-map methods agree on the largest exponent when they run long enough at adequate precision. <!-- research-question.md: row C1 -->
Claim C2 states that finite-time estimators differ substantially over realistic budgets. The difference is expected to be largest in weak chaos, in sticky layers and near separatrices. <!-- research-question.md: row C2 -->

The study also attempts to reproduce two published values. Sala and Artuso give 0.298 at K = pi/2 and 1.172 at K = 2 pi [@SalaArtuso2015]. <!-- RD: H4, preregistered rule [params: h4_ref_halfpi, h4_ref_2pi] -->

Two published comparators serve for the chaotic-sea exponent. These are the fit of Shevchenko [@Shevchenko2014] and the line ln(K/2) of Chirikov [@Chirikov1979]. The study uses both as comparators that never pass or fail. <!-- RD: H4, "Comparators (never pass or fail)" -->

This paper reports the study as preregistered. It adds no new claim, softens no result and offers no reinterpretation. Section 2 gives the design and the run identity. Sections 3 to 7 give the confirmatory results for hypotheses H1 to H5. Section 8 gives cost and sensitivity. Section 9 gives exploratory results. Sections 10 to 12 give limitations, availability and the conclusion.

## 2. Methods and preregistration

### 2.1 Preregistration

The decision rules were frozen on 2026-10-17, before any data existed. Amendments A1 to A6 were all frozen before the first study data file. The rules are applied here without change. <!-- RD: Header, "What the study tested" [findings A]; [manifest: prereg.header] -->
The preregistered value types were FALSIFIES (C1, C2), REPRODUCES (SalaArtuso2015), TOOL and NEGATIVE allowed. Claim C3 is out of scope. <!-- RD: Header, "What the study tested" [findings A] -->
The preregistration header is identified by the SHA-256 value in Table 1. <!-- RD: Header, run-identity table [manifest: prereg.header] -->

### 2.2 Panel and sets

The study compares estimators on a panel of 883 orbits at five nonzero values of the kick parameter K, plus exact-exponent controls. <!-- RD: Header, "What the study tested" [T1] -->
K index 0 holds the exact-exponent controls. K index 1 is K = 0. K index 2 to 6 are K = 0.5, 0.971635406 (K_c), 1.57079632679 (pi/2), 6.28318530718 (2 pi) and 10. <!-- RD: Header, "Reading key" [T5, rows K index 2 to 6] -->
Design classes are T (KAM torus), R (resonant), C (chaotic sea), S (sticky) and X (near a separatrix or a hyperbolic fixed point). <!-- RD: Header, "Reading key" [findings A] -->
"Class S" is a panel class. It is not the "strong-chaos set S" of H3, which is the class C orbits at K = 2 pi and K = 10. <!-- RD: Header, "Reading key" [findings A] -->
One orbit, orbit 881, is dropped as drift. <!-- RD: Header, "Panel and sets" [T1, last line] -->
The sets used by the decision rules have these sizes: E_ch 209, W 355, S 80 and S_L 92. H5b uses 883 orbits and H5c uses 334. <!-- RD: Header, "Panel and sets" [findings E, reading 8], [T7, rows H5a, H5b, H5c] -->
The program reads the sets as follows. S is the class C orbits at K = 2 pi and K = 10 without a label condition. W takes class C at 0 < K <= K_c. E_ch takes class C at K > 0. <!-- RD: Deviations and readings, reading 8 [findings E, reading 8] -->
The estimator methods L1 to L9 and Q are defined in the study design. They are not cited statements here. The repository (Section 11) holds the definitions. <!-- outline: Unmapped claims -->

### 2.3 Run identity

All reported numbers come from one fully automatic analysis run over the stored series. <!-- RD: Header, run-identity paragraph [findings A] -->
The run completed with 1363 chunks, all with status ok. <!-- RD: Header, run-identity table [manifest: status], [manifest: plan_counts.chunks], [manifest: chunks[].status] -->
Table 1 identifies the run. <!-- RD: Header, run-identity table -->

Table 1. Run identity.

| Item | Value | Pointer |
|---|---|---|
| Git commit (tables header and findings) | d31b655d0c33e2f2f1a67bbbdeba9c8c8a3fd036 | [findings A, table], [T Header: git commit] |
| Git sha recorded by the run manifest | 3e405f3324f8d25162580d3380db26c877d7ca63 | [manifest: git.sha] |
| Module version | 0.0.1 | [T Header: module version], [manifest: versions.module] |
| Panel SHA-256 | 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5 | [findings A, table], [manifest: panel.sha256_end] |
| Preregistration header SHA-256 | 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf | [findings A, table], [manifest: prereg.header] |
| Manifest SHA-256 | ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d | [findings A, table] |
| Parameter file SHA-256 | a544130bc56dc66d8b1f7ec5d1e54197306f14dc19c7cd22fdbd94e4534b7282 | [findings A, table] |
| Bootstrap seed | 20261012 | [findings A, table], [params: bootstrap seed in Header] |
| Run status | complete, 1363 chunks, all with status ok | [manifest: status], [manifest: plan_counts.chunks], [manifest: chunks[].status] |

<!-- RD: Header, run-identity table -->

The tables header and the findings file give the first git commit. The run manifest records the second. The study does not relate the two values, and neither does this paper. <!-- RD: Header, paragraph after run-identity table [findings A, table], [manifest: git.sha] -->

### 2.4 Deviations from the frozen text

The frozen text was applied with the deviations that the program states. <!-- RD: Deviations and readings [findings E, Deviations] -->

1. L2 runs in wrapped coordinates (Amendment A5 item 1). <!-- RD: Deviations and readings, deviation 1 [findings E, deviation 1; T Deviations] -->
2. The tau sweep records an orthonormalisation breakdown as a result, not as a failure (A5 item 2). <!-- RD: Deviations and readings, deviation 2 [findings E, deviation 2] -->
3. The control bound applies to controls 871, 873, 875, 877, 879 and 883 (A5 item 3). Controls 878, 880, 881 and 882 are reported as info. <!-- RD: Deviations and readings, deviation 3 [T Deviations; T2, controls table] -->
4. The control bound of L9 is 5.0e-8 and that of the other methods is 1.0e-8 (A6). <!-- RD: Deviations and readings, deviation 4 [params: control_bound, control_bound_L9; T Deviations] -->
5. The manifest records no overrides. <!-- RD: Deviations and readings, deviation 5 [manifest: overrides] -->

## 3. Results: H1 (claim C1), FALSIFIED

### 3.1 Preregistered rule

The rule applies on E_ch, which is class C orbits at K > 0 plus the exact-exponent controls. The horizon is N_1 = 1e7. <!-- RD: H1, preregistered rule [params: N_1] -->
An orbit passes H1a if the spread over {L1} + Q is at most 1e-8. It passes H1b if |lambda_L9 - lambda_L5| is at most 1e-6. It passes H1c if |lambda_L2 - lambda_L5| is at most 5e-3 at eps_0 = 1e-8. <!-- RD: H1, preregistered rule [params: h1a_bound, h1b_bound, h1c_bound, h1c_eps] -->
H1 is SUPPORTED only if the pass fraction is at least 0.98 for H1a and H1b and at least 0.90 for H1c. It is FALSIFIED if any sub-rule is below its limit. <!-- RD: H1, preregistered rule [params: h1a_fraction, h1b_fraction, h1c_fraction; findings B, H1] -->

### 3.2 Measured

There were 209 eligible orbits. The verdict as printed is falsified. <!-- RD: H1, Measured [T2, header line] -->

Table 2. H1 sub-rules on E_ch.

| Sub-rule | Limit | Passing orbits | Fraction | Result | Pointer |
|---|---|---|---|---|---|
| H1a | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1a_fraction] |
| H1b | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1b_fraction] |
| H1c | 0.90 | 162 of 209 | 0.775119617225 | fail | [T2, header line; params: h1c_fraction] |

<!-- RD: H1, Measured table -->

H1 is falsified. The H1c pass fraction was 162 of 209, which is 0.775119617225, against the limit of 0.90. <!-- RD: H1, [T2, header line] -->
The shortfall is carried by sub-rule H1c alone. H1a and H1b were met in full. <!-- RD: H1, [findings B, H1] -->

Table 3 gives H1c by cell.

Table 3. H1c passing orbits by cell.

| Cell | Orbits | H1c pass | Pointer |
|---|---|---|---|
| controls | 9 | 9 | [T2, row controls] |
| K index 2 | 40 | 27 | [T2, row K index 2] |
| K index 3 | 40 | 13 | [T2, row K index 3] |
| K index 4 | 40 | 33 | [T2, row K index 4] |
| K index 5 | 40 | 40 | [T2, row K index 5] |
| K index 6 | 40 | 40 | [T2, row K index 6] |

<!-- RD: H1, "H1c passing orbits by cell" -->

The 47 failing orbits sit in three K cells. At K index 2, 27 of 40 orbits pass. At K index 3, 13 of 40 pass. At K index 4, 33 of 40 pass. K index 5 and 6 pass 40 of 40, and the 9 controls pass. <!-- RD: H1, "H1c passing orbits by cell" [T2, rows K index 2 to 6, controls] -->
H1a and H1b pass in every cell. <!-- RD: H1, [T2, cell table] -->
The ids of the failing orbits and the extreme deviations are in [T2, Failing orbits]. Every failing orbit has |lambda_L2 - lambda_L5| above 0.005. The smallest is 0.00548462883667 and the largest is 0.0353748064771. <!-- RD: H1, [T2, Failing orbits, orbit 354], [T2, Failing orbits, orbit 368] -->
On the failing orbits the largest Delta_QB is 9.66310365058e-14 and the largest |lambda_L9 - lambda_L5| is 3.58331395889e-7. None of them breaks H1a or H1b. <!-- RD: H1, [T2, Failing orbits, orbit 195], [T2, Failing orbits, orbit 558]; [findings B, H1] -->

### 3.3 Exact-exponent controls

The controls under the bound are 871, 873, 875, 877, 879 and 883. Table 4 gives the largest deviation over the methods. <!-- RD: H1, "Exact-exponent controls" [findings E, deviation 3] -->

Table 4. Exact-exponent controls under the bound.

| Control | Scope | Largest deviation | L9 deviation | Pointer |
|---|---|---|---|---|
| 871 | bound | 5.26802590439e-9 | 5.26802579337e-9 | [T2, controls, 871] |
| 873 | bound | 2.48736364838e-9 | 8.326783706990001e-12 | [T2, controls, 873] |
| 875 | bound | 5.16177345311e-9 | 1.82781789704e-9 | [T2, controls, 875] |
| 877 | bound | 1.79536172595e-8 | 1.79536172595e-8 | [T2, controls, 877] |
| 879 | bound | 2.28260139679e-8 | 2.28260139679e-8 | [T2, controls, 879] |
| 883 | bound | 1.11022302463e-16 | 1.11022302463e-16 | [T2, controls, 883] |

<!-- RD: H1, "Exact-exponent controls" table -->

The largest deviations lie within the bounds of Amendment A6. For controls 877 and 879 the largest deviation is the L9 deviation. It is above 1e-8 and below the L9 control bound of 5.0e-8. <!-- RD: H1, after controls table [params: control_bound, control_bound_L9; findings A] -->
The table does not give the other methods' deviations on these two controls separately. <!-- RD: H1, after controls table [findings B, H1] -->

### 3.4 Reading

The preregistered rule for C1 is not met. H1 is reported as falsified. <!-- RD: H1, final paragraph [findings B, H1] -->
The test could have failed, and it did.

### 3.5 Exploratory items D1 and D2 (EXPLORATORY)

Everything in this subsection is EXPLORATORY. No frozen rule applies and no cause is claimed. <!-- RD: Exploratory results [findings D] -->

D1 (EXPLORATORY): where the H1c failures sit. The 47 H1c failures are in three K cells. K index 2 has 13 of 40 failing (27 pass). K index 3 has 27 of 40 (13 pass). K index 4 has 7 of 40 (33 pass). <!-- RD: D1 [T2, rows K index 2 to 4; findings D1] -->
K index 5 and 6 have no failures (40 of 40 pass each) and the 9 controls all pass. <!-- RD: D1 [T2, rows K index 5, K index 6, controls] -->
The failing ids fall in three ranges: 13 ids between 166 and 196, 27 between 346 and 385, and 7 between 549 and 580. <!-- RD: D1 [T2, Failing orbits; findings D1] -->
These counts match the per-cell failure counts, but the tables do not state the orbit-to-cell mapping. The tables do not show why the failures sit at the three smaller K. <!-- RD: D1 [findings D1] -->

D2 (EXPLORATORY): the info controls. Controls 878, 880, 881 and 882 are reported as info. Their largest deviations are 0.222903869545, 0.443364672685, 0.650205817015 and 0.103506221066. <!-- RD: D2 [T2, controls, 878], [T2, controls, 880], [T2, controls, 881], [T2, controls, 882] -->
These are far above the bound that applies to the other six controls. Under Amendment A5 item 3 these are the controls whose binary64 start is not an exact periodic point. The tables do not say more. <!-- RD: D2 [findings D2] -->

## 4. Results: H2 (claim C2), NOT SUPPORTED (NEGATIVE)

### 4.1 Preregistered rule

R(n) = Delta_V(n) / lambda_ref is evaluated at n_H2 = 1e4 steps over the weak-chaos set W. <!-- RD: H2, preregistered rule [params: n_H2] -->
H2 is SUPPORTED if the median of R(n_H2) is at least 0.10 and the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) is at least 0.05. Otherwise the null is reported. <!-- RD: H2, preregistered rule [params: h2_median, h2_interval, bootstrap_level, bootstrap_resamples; findings B, H2] -->

### 4.2 Measured

W has 355 orbits. The median of R(n_H2) is 2.23588035691e-5. The 95 percent bootstrap interval is [1.9257235885e-5, 2.31226443693e-5]. The verdict as printed is not supported. <!-- RD: H2, Measured [T3, header line] -->
Both conditions fail. The median is below 0.10 and the lower end is below 0.05. <!-- RD: H2, Measured [findings B, H2] -->

Table 5. H2 against the preregistered thresholds.

| Quantity | Measured | Threshold | Result | Pointer |
|---|---|---|---|---|
| Median of R(n_H2) over W (355 orbits) | 2.23588035691e-5 | at least 0.10 | fail | [T3, header line; params: h2_median] |
| Lower end of the 95 percent interval | 1.9257235885e-5 | at least 0.05 | fail | [T3, header line; params: h2_interval] |

<!-- RD: H2, Measured; preregistered rule -->

This is a null (NEGATIVE) result. It has the same standing as a pass. <!-- RD: H2, [findings B, H2] -->
In W, the relative spread among the methods of V at 1e4 steps is far below the preregistered level. <!-- RD: H2, [findings B, H2] -->
The frozen rule decides on the pooled W. The per-cell table in the results document is descriptive. <!-- RD: H2, [findings B, H2] -->
The cell of class S at K index 5 has 1 orbit, below the minimum set size of 20. It is not decidable and is neither pass nor fail. <!-- RD: H2, after cell table [T3, row S 5], [params: min_set], [findings B, H2] -->
Figure F2 shows R(n) by class and K, with the pooled series over W. <!-- RD: Reported without a decision rule, Sensitivity last paragraphs [findings C]; figure F2 -->

## 5. Results: H3 (claim C2), SUPPORTED

### 5.1 Preregistered rule

n_conv is the smallest n after which every method stays within a 0.05 relative band of lambda_ref up to N. Per orbit it is the maximum over the methods of V. <!-- RD: H3, preregistered rule [params: h3_band] -->
H3 is SUPPORTED if the ratio of the median n_conv over W to the median n_conv over the strong-chaos set S is at least 10, and the lower end of the 95 percent bootstrap interval of the ratio is above 3. <!-- RD: H3, preregistered rule [params: h3_ratio, h3_interval] -->
Orbits whose n_conv is undefined are given n_conv = N, and their share is reported. <!-- RD: H3, preregistered rule [findings B, H3] -->

### 5.2 Measured

Table 6. H3 measured values.

| Quantity | Value | Threshold | Pointer |
|---|---|---|---|
| Orbits in W | 355 | | [T4, header line] |
| Orbits in S | 80 | | [T4, header line] |
| Median n_conv over W | 2.511886e6 | | [T4, header line] |
| Median n_conv over S | 224.0 | | [T4, header line] |
| Ratio | 11213.7767857 | at least 10 | [T4, header line; params: h3_ratio] |
| 95 percent bootstrap interval of the ratio | [7924.465, 14117.3125] | lower end above 3 | [T4, header line; params: h3_interval] |
| Share of undefined n_conv | 0.0229885057471 | | [T4, header line] |
| S orbits without positive growth | 0 | | [T4, header line] |

<!-- RD: H3, Measured -->

H3 is supported. The verdict as printed is supported. Both conditions pass. The ratio is at least 10 and the lower end is above 3. <!-- RD: H3, Measured [T4, header line], [findings B, H3] -->
The strong-chaos set S of H3 is the class C orbits at K = 2 pi and K = 10. It is not the panel class S. <!-- RD: Header, "Reading key" [findings A], [findings E, reading 8] -->
The orbits with undefined n_conv enter the medians at the value N. <!-- RD: H3, Measured [findings B, H3] -->
The cell-level medians are descriptive. The class S cell at K index 5 has 1 orbit and is not decidable, as in H2. <!-- RD: H3, cell table [T4, rows ...], [findings B, H3] -->
Figure F3 shows the n_conv distributions. <!-- RD: Reported without a decision rule, Sensitivity last paragraphs [findings C]; figure F3 -->

## 6. Results: H4 (REPRODUCES, [@SalaArtuso2015]), DID NOT REPRODUCE

### 6.1 Preregistered rule

The start is x0 = 1e-3 and p0 = -1e-3, with T = 1e9 iterations. <!-- RD: H4, preregistered rule [params: h4_T; findings B, H4] -->
At K = pi/2 (reference 0.298) both L1 and L5 (Float64, tau = 1) must be within 0.01. At K = 2 pi (reference 1.172) both must be within 0.02. <!-- RD: H4, preregistered rule [params: h4_ref_halfpi, h4_tol_halfpi, h4_ref_2pi, h4_tol_2pi] -->
The reference values are those of Sala and Artuso [@SalaArtuso2015]. H4 is SUPPORTED if both K pass. A failure is reported as "did not reproduce" with the numbers. <!-- RD: H4, preregistered rule [findings A, findings B, H4] -->

### 6.2 Measured

The verdict as printed is did not reproduce. <!-- RD: H4, Measured [T5, verdict line] -->

Table 7. H4 verdict table.

| Start | L1 | L5 | Reference | Tolerance | Pass | Pointer |
|---|---|---|---|---|---|---|
| K = pi/2 | 0.256066843848 | 0.256066843848 | 0.298 | 0.01 | false | [T5, row K = pi/2] |
| K = 2 pi | 1.17018723872 | 1.17018723872 | 1.172 | 0.02 | true | [T5, row K = 2 pi] |

<!-- RD: H4, Measured table -->

H4 did not reproduce. At K = pi/2 both L1 and L5 gave 0.256066843848 against the reference 0.298. The value lies below the reference and outside the tolerance of 0.01. <!-- RD: H4, [T5, row K = pi/2] -->
The K = pi/2 start did not reproduce the published value [@SalaArtuso2015]. <!-- RD: H4, [T5, row K = pi/2] -->
At K = 2 pi the value 1.17018723872 against 1.172 is within the tolerance of 0.02. <!-- RD: H4, [T5, row K = 2 pi] -->
The reproduction as a whole fails. <!-- RD: H4, [findings B, H4] -->
The possibility that the orbit of the published paper sits in a different component remains open. The tables do not decide it. Each K rests on one start, with no error bar. <!-- RD: H4, final paragraph [findings B, H4] -->

### 6.3 Comparators (no verdict)

Table 8 gives three comparators at the five K values: the ensemble mean of lambda_L5(N_1) over the kept chaotic-sea orbits, the fit of Shevchenko [@Shevchenko2014], and the line ln(K/2) [@Chirikov1979]. They never pass or fail. <!-- RD: H4, "Comparators (never pass or fail)" -->
The frozen text designates the line for K = pi/2, 2 pi and 10. The program prints it at all five K. <!-- RD: H4, "Comparators (never pass or fail)" [findings B, H4] -->

Table 8. Comparators.

| K index | K | Orbits | Mean | Shevchenko fit | ln(K/2) | Pointer |
|---|---|---|---|---|---|---|
| 2 | 0.5 | 40 | 0.0575923614667 | 0.06665 | -1.38629436112 | [T5, comparators, K index 2] |
| 3 | 0.971635406 | 40 | 0.091095497804 | 0.12951899962 | -0.721921822156 | [T5, comparators, K index 3] |
| 4 | 1.57079632679 | 40 | 0.298004458001 | 0.342657808664 | -0.24156447527 | [T5, comparators, K index 4] |
| 5 | 6.28318530718 | 40 | 1.17214446011 | 1.17006018176 | 1.14472988585 | [T5, comparators, K index 5] |
| 6 | 10 | 40 | 1.62000488786 | 1.61943791243 | 1.60943791243 | [T5, comparators, K index 6] |

<!-- RD: H4, comparators table -->

### 6.4 Exploratory item D4 (EXPLORATORY)

D4 (EXPLORATORY): H4 at the shorter horizon. The cross-check uses T = 1e8. It has no verdict under the frozen rule. <!-- RD: H4, "Cross-check at the shorter horizon" [params: h4_T_cross]; D4 -->
At K = pi/2 it gives 0.299256534631 (Float64) and 0.299012763029 (Double64). At K = 2 pi it gives 1.17243899416 and 1.17179404609. <!-- RD: D4 [T5, cross-check, row 1.57079632679], [T5, cross-check, row 6.28318530718] -->
The Float64 value at K = pi/2 at T = 1e9 in the verdict table is 0.256066843848. The ensemble mean of lambda_L5(N_1) over the 40 kept chaotic-sea orbits at K = pi/2 is 0.298004458001. <!-- RD: D4 [T5, row K = pi/2], [T5, comparators, K index 4] -->
The tables do not show why the single start at K = pi/2 gives a different value at T = 1e9 than at T = 1e8. No cause is claimed. <!-- RD: D4 [findings D4] -->

## 7. Results: H5, SUPPORTED

### 7.1 Preregistered rule

On the ladder subset S_L, H5a is the fraction of S_L orbits whose two labels are identical at Float64, Double64 and BigFloat. It must be at least 0.90. <!-- RD: H5, preregistered rule [params: h5a_min] -->
H5b is the fraction of all panel orbits with either axis "unresolved". It must be at most 0.30. <!-- RD: H5, preregistered rule [params: h5b_max] -->
H5c is the fraction with structure "unresolved" among orbits with instability "none detected on the tested horizon". It must be at most 0.30. <!-- RD: H5, preregistered rule [params: h5c_max] -->
H5 is supported only if all three sub-rules pass. H5a is read from the per-precision labels before the ladder test acts. <!-- RD: H5, preregistered rule [findings B, H5; findings E, reading 5] -->

### 7.2 Measured

The verdict as printed is supported. <!-- RD: H5, Measured [T7, verdict line] -->

Table 9. H5 sub-rules.

| Sub-rule | Orbits | Fraction | Limit | Result | Pointer |
|---|---|---|---|---|---|
| H5a | 92 | 1.0 | at least 0.90 | pass | [T7, row H5a; params: h5a_min] |
| H5b | 883 | 0.0951302378256 | at most 0.30 | pass | [T7, row H5b; params: h5b_max] |
| H5c | 334 | 0.0568862275449 | at most 0.30 | pass | [T7, row H5c; params: h5c_max] |

<!-- RD: H5, Measured table -->

All three label conditions are met. H5 is supported. <!-- RD: H5, [findings B, H5] -->

### 7.3 Exploratory items D3 and D5 (EXPLORATORY)

D3 (EXPLORATORY): H5a is 1.0 while many ladder labels are forced. H5a is 1.0 on 92 S_L orbits, read before the ladder test acts. <!-- RD: D3 [T7, row H5a] -->
The ladder test then forced the labels of 35 of those 92 orbits to unresolved, a share of 0.380434782609. <!-- RD: D3 [F4, row ladder; T6, forced table] -->
Both hold at once: the per-label agreement across precisions, and the override of 35 labels by the ladder test. <!-- RD: D3 [findings D3] -->
The program prints the unresolved share before the tests as 46 of 883 (0.0520951302378) and after as 84 of 883 (0.0951302378256). The second value is the H5b fraction. <!-- RD: D3 [F4, row unresolved before the tests], [F4, row unresolved after the tests], [T7, row H5b] -->
Other forced tests: direction 7 of 172 (0.0406976744186). Estimator, diagnostics and excluded show 0 of 883 each. <!-- RD: D3 [F4, row direction], [F4, rows estimator, diagnostics, excluded] -->
The tables do not show how the forced counts (35 and 7) relate to the rise from 46 to 84. No cause is claimed. <!-- RD: D3 [findings D3] -->

D5 (EXPLORATORY): unresolved labels by cell. "Both" counts orbits with both axes unresolved (65 in all). "Structure only" counts orbits with instability "none detected on the tested horizon" and structure unresolved (19 in all). The two add to the 84 of H5b. <!-- RD: D5 [findings D5; T6, label table] -->
Cells not listed in Table 10 have none. <!-- RD: D5 -->

Table 10. Unresolved labels by cell (EXPLORATORY).

| K index | Class | Cell orbits | Both | Structure only | Pointer |
|---|---|---|---|---|---|
| 0 | controls | 13 | 3 | 0 | [T1, row 0 controls; T6, rows 0 controls] |
| 2 | C | 40 | 4 | 0 | [T1, row 2 C; T6, rows 2 C] |
| 2 | R | 40 | 5 | 0 | [T1, row 2 R; T6, rows 2 R] |
| 2 | S | 25 | 4 | 0 | [T1, row 2 S; T6, rows 2 S] |
| 2 | T | 40 | 0 | 11 | [T1, row 2 T; T6, rows 2 T] |
| 2 | X | 40 | 4 | 0 | [T1, row 2 X; T6, rows 2 X] |
| 3 | C | 40 | 6 | 0 | [T1, row 3 C; T6, rows 3 C] |
| 3 | R | 40 | 6 | 0 | [T1, row 3 R; T6, rows 3 R] |
| 3 | S | 40 | 3 | 0 | [T1, row 3 S; T6, rows 3 S] |
| 3 | T | 40 | 2 | 5 | [T1, row 3 T; T6, rows 3 T] |
| 3 | X | 40 | 4 | 0 | [T1, row 3 X; T6, rows 3 X] |
| 4 | C | 40 | 2 | 0 | [T1, row 4 C; T6, rows 4 C] |
| 4 | R | 40 | 6 | 1 | [T1, row 4 R; T6, rows 4 R] |
| 4 | S | 40 | 3 | 0 | [T1, row 4 S; T6, rows 4 S] |
| 4 | T | 40 | 2 | 2 | [T1, row 4 T; T6, rows 4 T] |
| 4 | X | 40 | 3 | 0 | [T1, row 4 X; T6, rows 4 X] |
| 5 | S | 5 | 4 | 0 | [T1, row 5 S; T6, rows 5 S] |
| 5 | T | 40 | 4 | 0 | [T1, row 5 T; T6, rows 5 T] |

<!-- RD: D5 table -->

K index 1, K index 6 and K index 5 classes C and X have no unresolved orbit. <!-- RD: D5 [findings D5; T6, rows 1, 6, 5 C, 5 X] -->

## 8. Cost and sensitivity (reported without a decision rule)

### 8.1 Cost

Table 11 gives ns per orbit-step by method and precision, and relative to L5 in the same precision. Each row is an unweighted mean over all items of that method and precision, so the item sets differ between rows. <!-- RD: Cost [findings C, Cost; findings E, reading 12] -->
Relative to L5 in Float64 (168.60570837 ns per orbit-step), L1 costs 0.701549448499, L2 0.503966109187, L4 0.982175001597, L9 1.23876330018 and L6 2.10273197634. <!-- RD: Cost table [F6, rows L5, L1, L2, L4, L9, L6 Float64] -->
L5 costs 4882.6924383 ns in Double64 and 377399.856865 ns in BigFloat(256). <!-- RD: Cost table [F6, row L5 Double64], [F6, row L5 BigFloat(256)] -->

Table 11. Cost per orbit-step.

| Method | Precision | Items | ns per orbit-step | Relative to L5 | Pointer |
|---|---|---|---|---|---|
| L5 | Float64 | 2605 | 168.60570837 | 1.0 | [F6, row L5 Float64] |
| L1 | Float64 | 1917 | 118.285241721 | 0.701549448499 | [F6, row L1 Float64] |
| L2 | Float64 | 1571 | 84.9715628342 | 0.503966109187 | [F6, row L2 Float64] |
| L4 | Float64 | 1915 | 165.600311888 | 0.982175001597 | [F6, row L4 Float64] |
| L9 | Float64 | 1915 | 208.86256373 | 1.23876330018 | [F6, row L9 Float64] |
| L6 | Float64 | 1227 | 354.532614384 | 2.10273197634 | [F6, row L6 Float64] |
| L5 | Double64 | 186 | 4882.6924383 | 1.0 | [F6, row L5 Double64] |
| L5 | BigFloat(256) | 92 | 377399.856865 | 1.0 | [F6, row L5 BigFloat(256)] |

<!-- RD: Cost table -->

The item count of the L1 Float64 row, 1917, is from the table row. The findings file lists that row without a count. <!-- RD: Cost, note after table [findings C, Cost] -->
Run totals: the total wall time was 11957.993616580963 seconds against a cap of 12.0 hours. The plan had 1363 chunks and 15569 items. The run used 23 threads of 23 available. All 13 arms are complete. <!-- RD: Cost, run totals [manifest: wall.total_seconds, wall.cap_hours, plan_counts.chunks, plan_counts.items, threads.used, threads.available, arms.*.status] -->

### 8.2 Sensitivity

The sensitivity table has 31 rows. The column is the median paired difference of lambda at n_H2 against the primary item, with 171 orbits per row except the tau rows. <!-- RD: Sensitivity [findings C, Sensitivity; F7] -->

Direction. L1 and L5 values are of order 1e-9 at offset 0.3 (L5 1.73141166638e-9). They are of order 1e-5 or 1e-6 at the other offsets: L5 is -3.32612943661e-5 at 0.3 + 2pi/5, -3.23438506002e-5 at 0.3 + 3pi/5 and -3.48457834277e-6 at 0.3 + pi/5. L9 is 0.0 in every direction row. <!-- RD: Sensitivity [F7, row direction 0.3 L5], [F7, row direction 0.3 + 2pi/5 L5], [F7, row direction 0.3 + 3pi/5 L5], [F7, row direction 0.3 + pi/5 L5], [F7, rows direction L9] -->

The eps rows for L2. With eps = 1e-10 the value is 2.26005874629e-9. With eps = 1e-12 it is 5.77729472453e-5. With eps = 1e-6 it is 4.04116577712e-10. <!-- RD: Sensitivity [F7, row eps 1.0e-10], [F7, row eps 1.0e-12], [F7, row eps 1.0e-6] -->

Wrapped versus lifted is 0.0 for all eight methods of V. <!-- RD: Sensitivity [F7, rows wrapped] -->
The column median_log10_nconv_ratio is 0.0 in all 31 rows. <!-- RD: Sensitivity [F7, column median_log10_nconv_ratio] -->

Tau sweep. Tau = 2 and 4 show no breakdown, with 171 orbits for L4 and L5. Tau = 8 has 154 orbits for L4 and L5, with breakdown count 17. Tau = 16 has 132 orbits for L4 with breakdown count 39 and 134 orbits for L5 with breakdown count 37. <!-- RD: Sensitivity [F7, rows tau 2, tau 4], [F7, rows tau 8], [F7, rows tau 16] -->
At tau = 16 the breakdown steps are 16 or 32 for L4 and 16 to 288 for L5. <!-- RD: Sensitivity [F7, breakdown listing, tau 16] -->

Symplectic residuals eps_pair(n) and eps_det(n) are shown in Figure F5 only. No numbers for them are quoted. <!-- RD: Sensitivity, last paragraphs [findings C] -->
Figure F1 shows the lambda_m(n) curves. <!-- RD: Sensitivity, last paragraphs [findings C] -->

## 9. Exploratory results

All items in this section are EXPLORATORY. They are not confirmatory, have no frozen rule and carry no claimed cause. <!-- RD: Exploratory results [findings D] -->
Items D1 to D5 appear in Sections 3, 6 and 7, where they bear on a hypothesis. D6 and D7 are reported here. <!-- RD: Exploratory results [findings D] -->

D6 (EXPLORATORY): short, empty and not decidable cells. The panel manifest marks two cells as short ("fewer qualifying grid cells than n_c"): K index 2 class S (25 drawn, 25 qualifying) and K index 5 class S (5 drawn, 5 qualifying). <!-- RD: D6 [T1, manifest cells table, rows 2 S, 5 S] -->
It marks as empty K index 5 class R (no elliptic periodic point) and K index 6 classes R (no elliptic periodic point), T and S (no qualifying grid cell). <!-- RD: D6 [T1, manifest cells, row 5 R], [T1, manifest cells, rows 6 R, 6 T, 6 S] -->
It marks K index 1 classes C, S and X as undefined (class undefined at K = 0). <!-- RD: D6 [T1, manifest cells, rows 1 C, 1 S, 1 X] -->
Class X shows 0 qualifying grid cells at every K index from 2 to 6. The preregistered selection rule does not draw class X from the screening and places the starts near the hyperbolic fixed point. <!-- RD: D6 [T1, manifest cells, rows 2 X to 6 X] [findings D6] -->
In the H2 and H3 tables only class S at K index 5 is not decidable, with 1 orbit. <!-- RD: D6 [T3, row S 5; T4, row S 5] -->

D7 (EXPLORATORY): cell patterns. The largest cell median of R(n_H2) is 0.000775501992792 (class S, K index 3). The smallest is 5.82402750446e-12 (class X, K index 2). <!-- RD: D7 [T3, row S 3], [T3, row X 2] -->
The cell medians of n_conv are between 4.223954e6 and 7.079458e6 at K index 2 and 3. The horizon there is N_1 = 1e7. <!-- RD: D7 [T4, row C 3], [T4, rows C 2, S 2, X 2], [params: N_1] -->
At K index 4 they are 794328.0 (class S) and 63096.0 (class X). <!-- RD: D7 [T4, row S 4], [T4, row X 4] -->
At K index 5 and 6 they lie between 158.0 and 266.5, excluding class S at K index 5. <!-- RD: D7 [T4, row X 6], [T4, row C 5], [findings D7] -->
The tables do not show why. <!-- RD: D7 [findings D7] -->

## 10. Limitations

1. The study rests on one draw of the panel (a single seed) with 883 orbits. The bootstrap intervals resample orbits within that draw and do not reflect a second draw. <!-- RD: Limitations, bullet 1 [T1] -->
2. The panel is run in Float64 only. Precision comparisons come from the ladder subset S_L of 92 orbits, and H5a rests on that subset. <!-- RD: Limitations, bullet 2 [T7, row H5a] -->
3. Cells with fewer than 20 orbits are not decidable or are descriptive. Class S at K index 5 (5 drawn, 1 in W) is not decidable. Class S at K index 2 (25 drawn, 21 in W) is descriptive. H4 rests on one start per K, as the frozen text specifies, with no error bar. <!-- RD: Limitations, bullets 5 and 6 [T1, manifest cells, rows 5 S, 2 S], [T3, row S 5], [T4, row S 5], [findings B, H4] -->
4. L2 runs in wrapped coordinates. This affects the H1c metric, the sub-rule that failed. The tables do not show whether this accounts for any part of the H1c result. <!-- RD: Deviations and readings, deviation 1 [findings E, deviation 1; T Deviations] -->
5. The numbers were checked against the tables file of the run. The analysis program was verified against an independent recomputation (`analysis-verification.md`). The numbers were not recomputed from the stored series for the results document. Figures F1 to F7 were checked only by file size and the `%PDF-` header, not visually. <!-- RD: Limitations, last bullet and open items [findings F] -->

## 11. Data and code availability

The repository is https://github.com/georgehadji/Chaos (GPL-3.0-or-later). <!-- _config/publication.md: repository -->
The tagged release v0.1.0 is not yet cut. It is cut at milestone M3 exit. The run used package version 0.0.1 at the commit given in Table 1, before the tag. <!-- RD: Deviations and readings, "Package version" [findings E, Note], [manifest: versions.module] -->
The run manifest and the stored series identify the run. The SHA-256 values are in Table 1. <!-- RD: Header, run-identity table -->

## 12. Conclusion

Under the preregistered rules, H1 (claim C1) is falsified. H2 (claim C2) is not supported (NEGATIVE). H3 (claim C2) is supported. H4 did not reproduce the published value at K = pi/2. H5 is supported. <!-- RD: "Claims C1 and C2" [findings B, H1, H2, H3]; H4 and H5 sections -->
H4 and H5 are not tied to C1 or C2 in the findings file. <!-- RD: "Claims C1 and C2" [findings B, H4, H5] -->

## Figure list

Figures are in `workspaces/study/stages/04-analyze/output/`. Captions follow the results document. Figures were checked only by file size and the `%PDF-` header.

- Figure F1: `p1-le-benchmark-fig-F1.pdf`. The lambda_m(n) curves.
- Figure F2: `p1-le-benchmark-fig-F2.pdf`. R(n) by class and K, with the pooled series over W.
- Figure F3: `p1-le-benchmark-fig-F3.pdf`. The n_conv distributions.
- Figure F4: `p1-le-benchmark-fig-F4.pdf`. The label tests (forced and unresolved counts, as read in D3).
- Figure F5: `p1-le-benchmark-fig-F5.pdf`. The symplectic residuals eps_pair(n) and eps_det(n).
- Figure F6: `p1-le-benchmark-fig-F6.pdf`. Cost per orbit-step.
- Figure F7: `p1-le-benchmark-fig-F7.pdf`. Sensitivity.

## References

References are resolved from `_config/bibliography/references.bib` at stage 04: [@SalaArtuso2015], [@Shevchenko2014], [@Chirikov1979].

<!--
Untraced items
- Journal: to be decided (outline open item 1); no journal-specific formatting applied.
- K_c = 0.971635406: used as given in RD; no citation key exists.
- Definitions of methods L1 to L9 and Q: no citation key; referred to the study design and repository only.
- Captions of Figure F4, F6 and F7: RD names F4, F6 and F7 only as table sections of the tables file; the captions above are inferred from RD headings (Cost, Sensitivity, D3). Not checked visually.
- Wording of C1 and C2 definitions: taken from _config/research-question.md rows C1 and C2, not from RD.
- Repository URL and GPL-3.0-or-later: taken from the outline (Claim 11.1) and _config/publication.md, not from RD.
- Outline open items 3 to 7 (two git SHAs, STRONG_K, Manifest.toml, figures not checked visually) are not resolved by this draft; only the two-SHA and figure-check facts are stated.
-->
