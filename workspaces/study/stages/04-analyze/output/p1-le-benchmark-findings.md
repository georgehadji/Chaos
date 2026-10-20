# Findings: study p1-le-benchmark (draft for the stage 04 checkpoint)

## A. What this document is

Study p1-le-benchmark belongs to paper P1 and tests claims C1 and C2 (C3 is out of scope, paper P2). It compares estimators of the largest Lyapunov exponent on a panel of orbits at five nonzero values of the kick parameter K, plus exact-exponent controls. The preregistered value types are FALSIFIES (C1, C2), REPRODUCES (SalaArtuso2015), TOOL and NEGATIVE allowed. The decision rules were frozen on 2026-10-17 before any data existed (with Amendments A1 to A6, all before the first study data file), and this document applies them without change.

All numbers below come from one fully automatic analysis run (no human or AI step in it) over the stored series, and from nothing else. The run is identified by these values.

| Item | Value |
|---|---|
| Preregistration header SHA-256 | 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf |
| Panel SHA-256 | 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5 |
| Manifest SHA-256 | ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d |
| Parameter file SHA-256 | a544130bc56dc66d8b1f7ec5d1e54197306f14dc19c7cd22fdbd94e4534b7282 |
| Git commit | d31b655d0c33e2f2f1a67bbbdeba9c8c8a3fd036 |
| Bootstrap seed | 20261012 |
| F4 data digest | 76fe8d67115ab241e37a5aed56f696ba797762139aa38a53da7453b6c4a5ef35 |
| F6 data digest | f13ce56e0deca246c780d3b9f728cb4a586669d7e46934be26a333d785609caf |

Reading key. K index 0 holds the exact-exponent controls. K index 1 is K = 0 (only classes T and R are populated). K index 2 to 6 are, from the comparator table (T5), K = 0.5, 0.971635406 (this is K_c), 1.57079632679 (pi/2), 6.28318530718 (2 pi) and 10. Design classes, as defined in the preregistration: T KAM torus, R resonant (near-)periodic, C chaotic sea, S sticky, X near a separatrix or a hyperbolic fixed point. "Class S" (a panel class) is not the "strong-chaos set S" of H3, which is the class C orbits at K = 2 pi and K = 10. Counts that print as 166.0 or 2.0 are formatting only.

Every threshold of the decision rules in the run's parameter list (H1 bounds and fractions, n_H2, the H2 and H3 limits, the H4 references and tolerances, the H5 limits, the minimum set size 20) equals the number in the frozen text. The L9 control bound of 5.0e-8 comes from Amendment A6 (frozen before any data); the other methods keep 1e-8.

## B. Confirmatory decisions

Order of decision: H1, H2, H3, H4, H5. The panel has 883 orbits; one (881) is dropped as drift (T1). Every eligible set below has at least 20 orbits, except where marked "not decidable".

### H1 (claim C1): FALSIFIED

Frozen rule. On E_ch (class C orbits at K > 0 plus the exact-exponent controls), at horizon N_1 = 1e7: H1a passes an orbit if the spread over {L1} + Q is <= 1e-8; H1b if |lambda_L9 - lambda_L5| <= 1e-6; H1c if |lambda_L2 - lambda_L5| <= 5e-3 at eps_0 = 1e-8. H1 is SUPPORTED only if the pass fraction is >= 0.98 for H1a and H1b and >= 0.90 for H1c. It is FALSIFIED if any sub-rule is below its limit.

Program output: eligible orbits 209. Verdict as printed: **falsified**.

| Sub-rule | Limit | Passing orbits | Fraction | Result |
|---|---|---|---|---|
| H1a | 0.98 | 209 of 209 | 1.0 | pass |
| H1b | 0.98 | 209 of 209 | 1.0 | pass |
| H1c | 0.90 | 162 of 209 | 0.775119617225 | **fail** |

By cell (T2), H1c passing orbits: controls 9 of 9, K index 2 27 of 40, K index 3 13 of 40, K index 4 33 of 40, K index 5 40 of 40, K index 6 40 of 40. H1a and H1b pass in every cell.

Meaning under the value type FALSIFIES. The preregistered rule for C1 is not met, and the shortfall is carried by sub-rule H1c alone. H1a and H1b are met in full. H1 is reported as falsified. The 47 failing orbits (209 minus 162) are listed in T2. Their ids are 166, 169, 171, 172, 173, 175, 183, 189, 190, 192, 194, 195, 196, 346, 348, 349, 350, 351, 352, 354, 355, 356, 357, 358, 359, 360, 361, 362, 365, 366, 367, 368, 370, 372, 373, 375, 377, 378, 381, 385, 549, 553, 558, 566, 568, 575, 580. Every failing orbit has |lambda_L2 - lambda_L5| above 0.005 (smallest listed 0.00548462883667, largest 0.0353748064771). On the same orbits the largest Delta_QB is 9.66310365058e-14 and the largest |lambda_L9 - lambda_L5| is 3.58331395889e-7, so none of them breaks H1a or H1b.

Exact-exponent controls (largest deviation over the methods; scope as in section E):

| Control | Scope | Largest deviation | L9 deviation |
|---|---|---|---|
| 871 | bound | 5.26802590439e-9 | 5.26802579337e-9 |
| 873 | bound | 2.48736364838e-9 | 8.326783706990001e-12 |
| 875 | bound | 5.16177345311e-9 | 1.82781789704e-9 |
| 877 | bound | 1.79536172595e-8 | 1.79536172595e-8 |
| 879 | bound | 2.28260139679e-8 | 2.28260139679e-8 |
| 883 | bound | 1.11022302463e-16 | 1.11022302463e-16 |

For 877 and 879 the largest deviation is the L9 deviation, which is above 1e-8 and below the L9 bound 5.0e-8 (A6). The table does not give the other methods' deviations on these two controls separately. Controls 878, 880, 881 and 882 are reported as info and are not held to the bound (section E).

### H2 (claim C2, finite-time inequivalence): NOT SUPPORTED (NEGATIVE)

Frozen rule. R(n) = Delta_V(n) / lambda_ref, evaluated at n_H2 = 1e4 steps over the weak-chaos set W. H2 is SUPPORTED if the median of R(n_H2) is >= 0.10 and the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) is >= 0.05. Otherwise the null is reported (NEGATIVE).

Program output: W orbits 355; median R(n_H2) 2.23588035691e-5; 95 percent bootstrap interval [1.9257235885e-5, 2.31226443693e-5]. Verdict as printed: **not supported**. Both conditions fail: the median is below 0.10 and the lower end is below 0.05.

Meaning. This is a null (NEGATIVE) result for H2 and has the same standing as a pass. In W, the relative spread among the methods of V at 1e4 steps is far below the preregistered level. By K and class (T3):

| class | K index | W orbits | median | lo | hi | note |
|---|---|---|---|---|---|---|
| C | 2 | 36 | 0.000685409790935 | 1.21772360607e-5 | 0.00158551651977 | |
| S | 2 | 21 | 5.42291035617e-5 | 3.17071908742e-6 | 0.00145378914565 | |
| X | 2 | 36 | 5.82402750446e-12 | 7.73045459959e-13 | 3.49257079589e-11 | |
| C | 3 | 34 | 0.000380436190835 | 7.27740118397e-5 | 0.000658941961158 | |
| S | 3 | 37 | 0.000775501992792 | 0.000343310022299 | 0.000939537300147 | |
| X | 3 | 36 | 2.26586896653e-5 | 2.21983704619e-5 | 2.30966041916e-5 | |
| S | 4 | 37 | 9.63721942695e-5 | 4.60435510795e-5 | 0.000139046731348 | |
| X | 4 | 37 | 2.32406707119e-5 | 2.32155590409e-5 | 2.32870005206e-5 | |
| S | 5 | 1 | 6.01344714582e-5 | 6.01344714582e-5 | 6.01344714582e-5 | not decidable |
| X | 5 | 40 | 1.92365473742e-5 | 1.92334364212e-5 | 1.92413315356e-5 | |
| X | 6 | 40 | 1.6241799957e-5 | 1.62397348555e-5 | 1.62425706568e-5 | |

The cell class S at K index 5 has 1 orbit, below the minimum of 20. It is **not decidable** (count 1) and is neither pass nor fail. The frozen rule decides on the pooled W; the cells are descriptive.

### H3 (claim C2, regime dependence of the convergence time): SUPPORTED

Frozen rule. n_conv is the smallest n after which every method stays within a 0.05 relative band of lambda_ref up to N; per orbit it is the maximum over the methods of V. H3 is SUPPORTED if (median n_conv over W) / (median n_conv over the strong-chaos set S) is >= 10 and the lower end of the 95 percent bootstrap interval of the ratio is > 3. Orbits whose n_conv is undefined are given n_conv = N, and their share is reported.

Program output: orbits W 355, S 80; median n_conv over W 2.511886e6; over S 224.0; ratio 11213.7767857; 95 percent bootstrap interval [7924.465, 14117.3125]; share of undefined n_conv 0.0229885057471; S orbits without positive growth 0. Verdict as printed: **supported**. Both conditions pass (ratio >= 10, lower end > 3).

Meaning. The preregistered regime-dependence criterion is met. A share of 0.0229885057471 of the orbits of W and S has undefined n_conv and enters the medians at the value N. By cell (T4, W and S orbits):

| class | K index | orbits | median | lo | hi | note |
|---|---|---|---|---|---|---|
| C | 2 | 36 | 7.079458e6 | 6.309573e6 | 7.943282e6 | |
| S | 2 | 21 | 7.079458e6 | 5.0e6 | 7.079458e6 | |
| X | 2 | 36 | 7.079458e6 | 5.623413e6 | 7.079458e6 | |
| C | 3 | 34 | 4.223954e6 | 3.03001e6 | 6.309573e6 | |
| S | 3 | 37 | 5.623413e6 | 3.548134e6 | 7.943282e6 | |
| X | 3 | 36 | 5.623413e6 | 3.981072e6 | 6.6945155e6 | |
| S | 4 | 37 | 794328.0 | 707946.0 | 1.0e6 | |
| X | 4 | 37 | 63096.0 | 44668.0 | 89125.0 | |
| C | 5 | 40 | 266.5 | 200.0 | 335.5 | |
| S | 5 | 1 | 501187.0 | 501187.0 | 501187.0 | not decidable |
| X | 5 | 40 | 251.0 | 168.0 | 447.0 | |
| C | 6 | 40 | 212.0 | 178.0 | 251.0 | |
| X | 6 | 40 | 158.0 | 133.5 | 224.0 | |

Class S at K index 5 (1 orbit) is **not decidable**, as in H2.

### H4 (REPRODUCES, SalaArtuso2015): DID NOT REPRODUCE

Frozen rule. Start x0 = 1e-3, p0 = -1e-3, T = 1e9 iterations. At K = pi/2 (reference 0.298) both L1 and L5 (Float64, tau = 1) must be within 0.01. At K = 2 pi (reference 1.172) both must be within 0.02. H4 is SUPPORTED if both K pass. A failure is reported as "did not reproduce" with the numbers, together with the possibility that the paper's orbit sits in a different component (the benchmark row notes a mixed phase space).

Program output. Verdict as printed: **did not reproduce**.

| start | L1 | L5 | reference | tolerance | pass |
|---|---|---|---|---|---|
| K = pi/2 | 0.256066843848 | 0.256066843848 | 0.298 | 0.01 | false |
| K = 2 pi | 1.17018723872 | 1.17018723872 | 1.172 | 0.02 | true |

Meaning. The K = 2 pi value is within its tolerance. The K = pi/2 value is not: L1 and L5 give the same number, which lies below the reference and outside the tolerance 0.01. The reproduction as a whole therefore fails. As the frozen text states, the possibility that the paper's orbit sits in a different component remains open, and the tables do not decide it. Each K rests on one start, with no error bar.

Cross-check at T = 1e8 (no verdict by the frozen rule; L5, Float64 against Double64):

| K | Float64 | Double64 |
|---|---|---|
| 1.57079632679 | 0.299256534631 | 0.299012763029 |
| 6.28318530718 | 1.17243899416 | 1.17179404609 |

Comparators, never pass or fail. The columns are the ensemble mean of lambda_L5(N_1) over the kept chaotic-sea orbits, the Shevchenko2014 fit, and the line ln(K/2) (Chirikov1979). The frozen text designates the line for K = pi/2, 2 pi and 10; the program prints it at all five K.

| K index | K | orbits | mean | Shevchenko2014 | ln(K/2) |
|---|---|---|---|---|---|
| 2 | 0.5 | 40 | 0.0575923614667 | 0.06665 | -1.38629436112 |
| 3 | 0.971635406 | 40 | 0.091095497804 | 0.12951899962 | -0.721921822156 |
| 4 | 1.57079632679 | 40 | 0.298004458001 | 0.342657808664 | -0.24156447527 |
| 5 | 6.28318530718 | 40 | 1.17214446011 | 1.17006018176 | 1.14472988585 |
| 6 | 10 | 40 | 1.62000488786 | 1.61943791243 | 1.60943791243 |

### H5 (are the labels informative?): SUPPORTED

Frozen rule. On the ladder subset S_L: H5a, the fraction of S_L orbits whose two labels are identical at Float64, Double64 and BigFloat, is >= 0.90. H5b, the fraction of all panel orbits with either axis "unresolved", is <= 0.30. H5c, the fraction with structure "unresolved" among orbits with instability "none detected on the tested horizon", is <= 0.30. By the program's reading, H5 is supported only if all three pass. If H5b or H5c failed, decisions D36 or D37 would go back to the owner (a design result, not a study failure).

Program output (T7). Verdict as printed: **supported**.

| Sub-rule | Orbits | Fraction | Limit | Result |
|---|---|---|---|---|
| H5a | 92 | 1.0 | >= 0.90 | pass |
| H5b | 883 | 0.0951302378256 | <= 0.30 | pass |
| H5c | 334 | 0.0568862275449 | <= 0.30 | pass |

Meaning. All three label conditions of the frozen text are met, so no label decision goes back to the owner under these rules. H5a is taken before the ladder test acts (section E); see D3.

## C. Reported without a decision rule

Cost (F6). The table gives ns per orbit-step by method and precision, and relative to L5 in the same precision. Each row is an unweighted mean over all items of that method and precision, so item sets differ between rows (L5 Float64 2605 items, L9 Float64 1915, L2 Float64 1571, L6 Float64 1227, BigFloat rows 92). Float64 rows (ns per orbit-step, relative to L5): L5 168.60570837 (1.0); L1 118.285241721 (0.701549448499); L2 84.9715628342 (0.503966109187); L4 165.600311888 (0.982175001597); L9 208.86256373 (1.23876330018); L6 354.532614384 (2.10273197634). L5 is 4882.6924383 in Double64 and 377399.856865 in BigFloat(256). See figure F6 for the full table.

Sensitivity (F7), 31 rows. The column is the median paired difference of lambda at n_H2 against the primary item (median_dlambda_nH2), 171 orbits per row except the tau rows below.
- Direction: L1 and L5 values are of order 1e-9 at offset 0.3 (L5 1.73141166638e-9) and of order 1e-5 or 1e-6 at the other offsets (L5 -3.32612943661e-5 at 0.3 + 2pi/5, -3.23438506002e-5 at 0.3 + 3pi/5, -3.48457834277e-6 at 0.3 + pi/5); L9 is 0.0 in every direction row.
- L2 with eps = 1e-10 gives 2.26005874629e-9, eps = 1e-12 gives 5.77729472453e-5, eps = 1e-6 gives 4.04116577712e-10.
- Wrapped versus lifted: 0.0 for all eight methods of V.
- median_log10_nconv_ratio is 0.0 in all 31 rows.
- Tau sweep: tau = 2 and 4 show no breakdown (171 orbits for L4 and L5). Tau = 8 has 154 orbits for L4 and L5, with breakdown count 17. Tau = 16 has 132 orbits for L4 (breakdown count 39) and 134 for L5 (breakdown count 37). The breakdown listing of F7 gives step and kind for each item; at tau = 16 the steps are 16 or 32 for L4 and 16 to 288 for L5.

Symplectic residuals eps_pair(n) and eps_det(n). No numbers for these are quoted in this document. See figure F5 (quantiles by method against the tolerance multipliers). Figures F1 (lambda_m(n) curves), F2 (R(n) by class and K, with the pooled series over W), F3 (n_conv distributions) and F5 are not reproduced here.

## D. Exploratory observations

Everything in this section is EXPLORATORY. None is a decision of section B, none has a frozen rule, and none is to be read as confirmatory. The tables do not give causes, and none is claimed.

**D1 (EXPLORATORY, T2): where the H1c failures sit.** The 47 H1c failures are in three K cells. K index 2 has 13 of 40 failing (27 pass), K index 3 has 27 of 40 (13 pass) and K index 4 has 7 of 40 (33 pass). K index 5 and 6 have no failures (40 of 40 pass each), and the 9 controls all pass. The failing ids fall in three ranges: 13 ids between 166 and 196, 27 between 346 and 385, and 7 between 549 and 580. These counts match the per-cell failure counts, but the tables do not state the orbit-to-cell mapping. The tables do not show why the failures sit at the three smaller K.

**D2 (EXPLORATORY, T2): the info controls.** Controls 878, 880, 881 and 882 are reported as info. Their largest deviations are 0.222903869545, 0.443364672685, 0.650205817015 and 0.103506221066, far above the bound that applies to the other six controls. Under Amendment A5 item 3 these are the controls whose binary64 start is not an exact periodic point; the tables do not say more.

**D3 (EXPLORATORY, T6, F4, T7): H5a is 1.0 while many ladder labels are forced.** H5a is 1.0 on 92 S_L orbits, read before the ladder test acts. The ladder test then forced the labels of 35 of those 92 orbits to unresolved (share 0.380434782609). Both hold at once: the per-precision labels agree, and the ladder test overrides 35 of them. The program prints the unresolved share before the tests as 46 of 883 (0.0520951302378) and after as 84 of 883 (0.0951302378256). The second value is the H5b fraction. Other forced tests: direction 7 of 172 (0.0406976744186); estimator, diagnostics and excluded 0 of 883 each. The tables do not show how the forced counts (35 and 7) relate to the rise from 46 to 84.

**D4 (EXPLORATORY, T5): H4 at the shorter horizon.** At T = 1e8 the cross-check gives, at K = pi/2, 0.299256534631 (Float64) and 0.299012763029 (Double64). At K = 2 pi it gives 1.17243899416 and 1.17179404609. The Float64 value at K = pi/2 at T = 1e9 in the verdict table is 0.256066843848. The ensemble mean of lambda_L5(N_1) over the 40 kept chaotic-sea orbits at K = pi/2 is 0.298004458001. The tables do not show why the single start at K = pi/2 gives a different value at T = 1e9 than at T = 1e8.

**D5 (EXPLORATORY, T6): unresolved labels by cell.** "Both" counts orbits with both axes unresolved (65 in all). "Structure only" counts orbits with instability "none detected on the tested horizon" and structure unresolved (19 in all). The two add to the 84 of H5b. Cells not listed have none.

| K index | class | cell orbits | both | structure only |
|---|---|---|---|---|
| 0 | controls | 13 | 3 | 0 |
| 2 | C | 40 | 4 | 0 |
| 2 | R | 40 | 5 | 0 |
| 2 | S | 25 | 4 | 0 |
| 2 | T | 40 | 0 | 11 |
| 2 | X | 40 | 4 | 0 |
| 3 | C | 40 | 6 | 0 |
| 3 | R | 40 | 6 | 0 |
| 3 | S | 40 | 3 | 0 |
| 3 | T | 40 | 2 | 5 |
| 3 | X | 40 | 4 | 0 |
| 4 | C | 40 | 2 | 0 |
| 4 | R | 40 | 6 | 1 |
| 4 | S | 40 | 3 | 0 |
| 4 | T | 40 | 2 | 2 |
| 4 | X | 40 | 3 | 0 |
| 5 | S | 5 | 4 | 0 |
| 5 | T | 40 | 4 | 0 |

K index 1, K index 6 and K index 5 classes C and X have no unresolved orbit.

**D6 (EXPLORATORY, T1): short, empty and not decidable cells.** The panel manifest marks two cells as short ("fewer qualifying grid cells than n_c"): K index 2 class S (25 drawn, 25 qualifying) and K index 5 class S (5 drawn, 5 qualifying). It marks as empty K index 5 class R (no elliptic periodic point) and K index 6 classes R (no elliptic periodic point), T and S (no qualifying grid cell). It marks K index 1 classes C, S and X as undefined (class undefined at K = 0). Class X shows 0 qualifying grid cells at every K index from 2 to 6 because the preregistered selection rule does not draw class X from the screening (it places the starts near the hyperbolic fixed point). In T3 and T4 only class S at K index 5 is not decidable (1 orbit).

**D7 (EXPLORATORY, T3, T4): cell patterns.** The largest cell median of R(n_H2) is 0.000775501992792 (class S, K index 3) and the smallest is 5.82402750446e-12 (class X, K index 2). The cell medians of n_conv are between 4.223954e6 and 7.079458e6 at K index 2 and 3 (N_1 = 1e7), are 794328.0 (class S) and 63096.0 (class X) at K index 4, and lie between 158.0 and 266.5 at K index 5 and 6 (excluding class S at K index 5). The tables do not show why.

## E. Readings and deviations

The program states the readings where the frozen text is silent and the deviations from the frozen text. The texts of the amendments A1 to A6 are in the preregistration; the descriptions below are the program's own.

Readings

1. R(n) and n_conv use the horizon-N_1 series of every orbit; lambda_ref and the labels use the longest horizon of the orbit (P3). Affects H2 and H3 (the series used) and W membership and H5 (the labels).
2. A failing label test forces both axes to unresolved (P4). Affects H5b, H5c and W membership.
3. The second-half rotation number is (x_N - x_h) / (2 pi (N - h)); h is the recorded step nearest N/2, N/2 itself where it is recorded (P5). Enters the label computation, so it bears on T6 and H5.
4. Labels across the three precisions are compared at N_B (P5). Affects H5a and the ladder test.
5. H5a is read from the per-precision labels before the ladder test acts; H5 is supported only if H5a, H5b and H5c all pass (P6). Affects the H5a value (1.0 while 35 labels are forced, D3) and the H5 verdict rule.
6. Percentile bootstrap, orbits resampled with replacement, W and S resampled independently, streams keyed by seed and label (P7). Affects the intervals of H2 and H3.
7. F1 shows the first drawn orbit of each cell, F3 empirical distribution functions, F4 shares per test, F5 quantiles, F7 paired differences (P8). Affects the figures only.
8. S is the class C orbits at K = 2 pi and K = 10 without a label condition; W takes class C at 0 < K <= K_c; E_ch takes class C at K > 0. Fixes the membership of S (80 orbits), W (355) and E_ch (209).
9. Periodic points are the controls of kind fixed or newton; kind y5 is not periodic. Affects how the controls are labelled.
10. The ladder compares the three precisions with one common half point: the recorded step nearest N_B/2 that all three items have. Affects the ladder test and its share in F4.
11. Diagnostics exceeded is read on the L5 item at N_1 and on the L5 item at the label horizon. Affects the diagnostics test (0 forced).
12. F6 is the unweighted mean of thread-nanoseconds per orbit-step over all items of a method and precision (all arms, all sweep values). Affects the cost table only; item sets differ between rows.
13. F2 has a pooled series over W next to the series by class and K; orbits with lambda_ref <= 0 are left out of F2. Affects figure F2 only.
14. Flags drift and excluded are the union over the items of an orbit; diagnostics exceeded is read on the primary L5 item of arm main. Affects which orbits are dropped (881 in T1).

Deviations

1. L2 runs in wrapped coordinates (A5 item 1); H1c compares the wrapped L2 value with the lifted L5 value of the same orbit. Affects the H1c metric, the sub-rule that failed. The tables do not show whether this accounts for any part of the H1c result.
2. The tau sweep records an orthonormalisation breakdown as a result, not as a failure (A5 item 2); orbits with the flag: 40. Affects the tau rows of F7 (fewer than 171 orbits at tau = 8 and 16).
3. Control bound scope (A5 item 3): under the bound 871, 873, 875, 877, 879, 883; reported as info 878, 880, 881, 882. Affects which controls are tested against the bound (T2 controls table).
4. The control bound of L9 is 5.0e-8, of the other methods 1.0e-8 (A6). Affects the control check of 877 and 879, whose largest deviations are L9 values.
5. Overrides recorded in the manifest: none.

Note. The TOOL value type names the tagged release v0.1.0, which is cut at milestone M3 exit; this run was made with package version 0.0.1 at the commit given in section A, before the tag.

## F. Limits

- One draw of the panel (a single seed), 883 orbits. The bootstrap intervals resample orbits within that draw and do not reflect a second draw.
- The panel is run in Float64 only. Precision comparisons come from the ladder subset S_L of 92 orbits (Float64, Double64, BigFloat), and H5a rests on that subset.
- Sets as printed: E_ch 209 (H1), W 355 and S 80 (H2, H3), S_L 92 (H5a), 883 (H5b), 334 (H5c). F7 rows use 171 orbits except the tau rows (132 to 171).
- One orbit (881) is dropped as drift. The controls table lists it as an info control.
- Cells with fewer than 20 orbits: class S at K index 5 (5 drawn, 1 in W) is not decidable in T3 and T4. Class S at K index 2 has 25 drawn and 21 in W. Cell-level medians and intervals are descriptive; the decisions use the pooled sets.
- H4 rests on one start per K, as the frozen text specifies, with no error bar.
- No numbers for the symplectic residual figure F5 are quoted here.
- The numbers were checked against the tables file of the run (and the program itself was verified against an independent recomputation, `analysis-verification.md`); they were not recomputed here from the stored series. Table cosmetics (counts printing as 2.0 or 166.0) are formatting only.
