# Results: study p1-le-benchmark

## Pointer key

Every number carries a pointer in square brackets.

- `[T1]` to `[T7]`, `[F4]`, `[F6]`, `[F7]`: sections of `workspaces/study/stages/04-analyze/output/p1-le-benchmark-tables.md`, with the row named after the comma.
- `[params]`: the parameters line of the Header of the same tables file.
- `[findings X]`: section X of `workspaces/study/stages/04-analyze/output/p1-le-benchmark-findings.md`.
- `[manifest: field]`: a field of `workspaces/study/stages/03-run/output/p1-le-benchmark-manifest.json`.
- In tables, the last column is the pointer for the whole row.

## Header

Study slug: p1-le-benchmark [manifest: study].

Run identity. All numbers in this document come from one fully automatic analysis run over the stored series [findings A]. The run is identified as follows.

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

The findings file gives the first git commit above [findings A, table] and the manifest records the second [manifest: git.sha]. The findings file does not relate the two values, and this document does not either.

What the study tested. Study p1-le-benchmark belongs to paper P1 and tests claims C1 and C2 [findings A]. It compares estimators of the largest Lyapunov exponent on a panel of 883 orbits [T1] at five nonzero values of the kick parameter K, plus exact-exponent controls [findings A]. The preregistered value types were FALSIFIES (C1, C2), REPRODUCES (SalaArtuso2015), TOOL and NEGATIVE allowed [findings A]. The decision rules were frozen on 2026-10-17, before any data existed, with Amendments A1 to A6 all frozen before the first study data file [findings A]. The rules are applied here without change. Claim C3 is out of scope (paper P2) [findings A].

Reading key. K index 0 holds the exact-exponent controls. K index 1 is K = 0. K index 2 to 6 are K = 0.5, 0.971635406 (K_c), 1.57079632679 (pi/2), 6.28318530718 (2 pi) and 10 [T5, rows K index 2 to 6]. Design classes are T (KAM torus), R (resonant), C (chaotic sea), S (sticky) and X (near a separatrix or a hyperbolic fixed point) [findings A]. "Class S" is a panel class. It is not the "strong-chaos set S" of H3, which is the class C orbits at K = 2 pi and K = 10 [findings A]. Counts that print as `166.0` or `2.0` in the tables are formatting only [findings A].

Panel and sets. The panel has 883 orbits and one, orbit 881, is dropped as drift [T1, last line]. The sets used below have these sizes: E_ch 209, W 355, S 80 [findings E, reading 8], S_L 92 [T7, row H5a], 883 for H5b and 334 for H5c [T7, rows H5b, H5c].

## Confirmatory results

### H1 (claim C1): FALSIFIED

Preregistered rule. On E_ch (class C orbits at K > 0 plus the exact-exponent controls), at horizon N_1 = 1e7 [params: N_1], an orbit passes H1a if the spread over {L1} + Q is at most 1e-8 [params: h1a_bound], H1b if |lambda_L9 - lambda_L5| is at most 1e-6 [params: h1b_bound], and H1c if |lambda_L2 - lambda_L5| is at most 5e-3 [params: h1c_bound] at eps_0 = 1e-8 [params: h1c_eps]. H1 is SUPPORTED only if the pass fraction is at least 0.98 for H1a and H1b [params: h1a_fraction, h1b_fraction] and at least 0.90 for H1c [params: h1c_fraction]. It is FALSIFIED if any sub-rule is below its limit [findings B, H1].

Measured. Eligible orbits: 209 [T2, header line]. Verdict as printed: falsified [T2, header line].

| Sub-rule | Limit | Passing orbits | Fraction | Result | Pointer |
|---|---|---|---|---|---|
| H1a | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1a_fraction] |
| H1b | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1b_fraction] |
| H1c | 0.90 | 162 of 209 | 0.775119617225 | fail | [T2, header line; params: h1c_fraction] |

H1c failed. The pass fraction was 162 of 209 = 0.775119617225 against the limit of 0.90 [T2, header line]. This is a negative result for C1 under the preregistered rule. The shortfall is carried by sub-rule H1c alone, and H1a and H1b were met in full [findings B, H1].

H1c passing orbits by cell:

| Cell | Orbits | H1c pass | Pointer |
|---|---|---|---|
| controls | 9 | 9 | [T2, row controls] |
| K index 2 | 40 | 27 | [T2, row K index 2] |
| K index 3 | 40 | 13 | [T2, row K index 3] |
| K index 4 | 40 | 33 | [T2, row K index 4] |
| K index 5 | 40 | 40 | [T2, row K index 5] |
| K index 6 | 40 | 40 | [T2, row K index 6] |

H1a and H1b pass in every cell [T2, cell table]. The 47 failing orbits (209 minus 162) are listed in [T2, Failing orbits]. Their ids are 166, 169, 171, 172, 173, 175, 183, 189, 190, 192, 194, 195, 196, 346, 348, 349, 350, 351, 352, 354, 355, 356, 357, 358, 359, 360, 361, 362, 365, 366, 367, 368, 370, 372, 373, 375, 377, 378, 381, 385, 549, 553, 558, 566, 568, 575, 580 [T2, Failing orbits; findings B, H1]. Every failing orbit has |lambda_L2 - lambda_L5| above 0.005, the smallest being 0.00548462883667 [T2, Failing orbits, orbit 354] and the largest 0.0353748064771 [T2, Failing orbits, orbit 368]. On the failing orbits the largest Delta_QB is 9.66310365058e-14 [T2, Failing orbits, orbit 195] and the largest |lambda_L9 - lambda_L5| is 3.58331395889e-7 [T2, Failing orbits, orbit 558], so none of them breaks H1a or H1b [findings B, H1].

Exact-exponent controls (largest deviation over the methods; the controls under the bound are 871, 873, 875, 877, 879 and 883) [findings E, deviation 3]:

| Control | Scope | Largest deviation | L9 deviation | Pointer |
|---|---|---|---|---|
| 871 | bound | 5.26802590439e-9 | 5.26802579337e-9 | [T2, controls, 871] |
| 873 | bound | 2.48736364838e-9 | 8.326783706990001e-12 | [T2, controls, 873] |
| 875 | bound | 5.16177345311e-9 | 1.82781789704e-9 | [T2, controls, 875] |
| 877 | bound | 1.79536172595e-8 | 1.79536172595e-8 | [T2, controls, 877] |
| 879 | bound | 2.28260139679e-8 | 2.28260139679e-8 | [T2, controls, 879] |
| 883 | bound | 1.11022302463e-16 | 1.11022302463e-16 | [T2, controls, 883] |

For 877 and 879 the largest deviation is the L9 deviation. It is above 1e-8 and below the L9 control bound of 5.0e-8, which comes from Amendment A6 [params: control_bound, control_bound_L9; findings A]. The table does not give the other methods' deviations on these two controls separately [findings B, H1]. Controls 878, 880, 881 and 882 are reported as info and are not held to the bound (see D2).

Reading under value type FALSIFIES: the preregistered rule for C1 is not met, and H1 is reported as falsified [findings B, H1].

### H2 (claim C2, finite-time inequivalence): NOT SUPPORTED (NEGATIVE)

Preregistered rule. R(n) = Delta_V(n) / lambda_ref is evaluated at n_H2 = 1e4 steps [params: n_H2] over the weak-chaos set W. H2 is SUPPORTED if the median of R(n_H2) is at least 0.10 [params: h2_median] and the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) is at least 0.05 [params: h2_interval, bootstrap_level, bootstrap_resamples]. Otherwise the null is reported [findings B, H2].

Measured. W orbits: 355. Median R(n_H2): 2.23588035691e-5. 95 percent bootstrap interval: [1.9257235885e-5, 2.31226443693e-5] [T3, header line]. Verdict as printed: not supported [T3, header line]. Both conditions fail: the median is below 0.10 and the lower end is below 0.05 [findings B, H2].

This is a null (NEGATIVE) result for H2 and has the same standing as a pass [findings B, H2]. In W, the relative spread among the methods of V at 1e4 steps is far below the preregistered level [findings B, H2]. The frozen rule decides on the pooled W. The cells below are descriptive [findings B, H2].

| Class | K index | W orbits | Median | Lo | Hi | Note | Pointer |
|---|---|---|---|---|---|---|---|
| C | 2 | 36 | 0.000685409790935 | 1.21772360607e-5 | 0.00158551651977 | | [T3, row C 2] |
| S | 2 | 21 | 5.42291035617e-5 | 3.17071908742e-6 | 0.00145378914565 | | [T3, row S 2] |
| X | 2 | 36 | 5.82402750446e-12 | 7.73045459959e-13 | 3.49257079589e-11 | | [T3, row X 2] |
| C | 3 | 34 | 0.000380436190835 | 7.27740118397e-5 | 0.000658941961158 | | [T3, row C 3] |
| S | 3 | 37 | 0.000775501992792 | 0.000343310022299 | 0.000939537300147 | | [T3, row S 3] |
| X | 3 | 36 | 2.26586896653e-5 | 2.21983704619e-5 | 2.30966041916e-5 | | [T3, row X 3] |
| S | 4 | 37 | 9.63721942695e-5 | 4.60435510795e-5 | 0.000139046731348 | | [T3, row S 4] |
| X | 4 | 37 | 2.32406707119e-5 | 2.32155590409e-5 | 2.32870005206e-5 | | [T3, row X 4] |
| S | 5 | 1 | 6.01344714582e-5 | 6.01344714582e-5 | 6.01344714582e-5 | not decidable | [T3, row S 5] |
| X | 5 | 40 | 1.92365473742e-5 | 1.92334364212e-5 | 1.92413315356e-5 | | [T3, row X 5] |
| X | 6 | 40 | 1.6241799957e-5 | 1.62397348555e-5 | 1.62425706568e-5 | | [T3, row X 6] |

The cell of class S at K index 5 has 1 orbit [T3, row S 5], below the minimum set size of 20 [params: min_set]. It is not decidable and is neither pass nor fail [findings B, H2].

### H3 (claim C2, regime dependence of the convergence time): SUPPORTED

Preregistered rule. n_conv is the smallest n after which every method stays within a 0.05 relative band [params: h3_band] of lambda_ref up to N, and per orbit it is the maximum over the methods of V. H3 is SUPPORTED if (median n_conv over W) / (median n_conv over the strong-chaos set S) is at least 10 [params: h3_ratio] and the lower end of the 95 percent bootstrap interval of the ratio is above 3 [params: h3_interval]. Orbits whose n_conv is undefined are given n_conv = N, and their share is reported [findings B, H3].

Measured. Orbits in W: 355. Orbits in S: 80. Median n_conv over W: 2.511886e6. Median n_conv over S: 224.0. Ratio: 11213.7767857. 95 percent bootstrap interval: [7924.465, 14117.3125]. Share of undefined n_conv: 0.0229885057471. S orbits without positive growth: 0 [T4, header line]. Verdict as printed: supported [T4, header line]. Both conditions pass: the ratio is at least 10 and the lower end is above 3 [findings B, H3]. The orbits with undefined n_conv enter the medians at the value N [findings B, H3].

| Class | K index | Orbits | Median | Lo | Hi | Note | Pointer |
|---|---|---|---|---|---|---|---|
| C | 2 | 36 | 7.079458e6 | 6.309573e6 | 7.943282e6 | | [T4, row C 2] |
| S | 2 | 21 | 7.079458e6 | 5.0e6 | 7.079458e6 | | [T4, row S 2] |
| X | 2 | 36 | 7.079458e6 | 5.623413e6 | 7.079458e6 | | [T4, row X 2] |
| C | 3 | 34 | 4.223954e6 | 3.03001e6 | 6.309573e6 | | [T4, row C 3] |
| S | 3 | 37 | 5.623413e6 | 3.548134e6 | 7.943282e6 | | [T4, row S 3] |
| X | 3 | 36 | 5.623413e6 | 3.981072e6 | 6.6945155e6 | | [T4, row X 3] |
| S | 4 | 37 | 794328.0 | 707946.0 | 1.0e6 | | [T4, row S 4] |
| X | 4 | 37 | 63096.0 | 44668.0 | 89125.0 | | [T4, row X 4] |
| C | 5 | 40 | 266.5 | 200.0 | 335.5 | | [T4, row C 5] |
| S | 5 | 1 | 501187.0 | 501187.0 | 501187.0 | not decidable | [T4, row S 5] |
| X | 5 | 40 | 251.0 | 168.0 | 447.0 | | [T4, row X 5] |
| C | 6 | 40 | 212.0 | 178.0 | 251.0 | | [T4, row C 6] |
| X | 6 | 40 | 158.0 | 133.5 | 224.0 | | [T4, row X 6] |

Class S at K index 5 (1 orbit) is not decidable, as in H2 [findings B, H3].

### H4 (REPRODUCES, SalaArtuso2015): DID NOT REPRODUCE

Preregistered rule. Start x0 = 1e-3, p0 = -1e-3, T = 1e9 iterations [params: h4_T; findings B, H4]. At K = pi/2 (reference 0.298) both L1 and L5 (Float64, tau = 1) must be within 0.01. At K = 2 pi (reference 1.172) both must be within 0.02 [params: h4_ref_halfpi, h4_tol_halfpi, h4_ref_2pi, h4_tol_2pi]. The reference values are those of Sala and Artuso [@SalaArtuso2015] as named by the value type REPRODUCES [findings A]. H4 is SUPPORTED if both K pass. A failure is reported as "did not reproduce" with the numbers, together with the possibility that the paper's orbit sits in a different component [findings B, H4].

Measured. Verdict as printed: did not reproduce [T5, verdict line].

| Start | L1 | L5 | Reference | Tolerance | Pass | Pointer |
|---|---|---|---|---|---|---|
| K = pi/2 | 0.256066843848 | 0.256066843848 | 0.298 | 0.01 | false | [T5, row K = pi/2] |
| K = 2 pi | 1.17018723872 | 1.17018723872 | 1.172 | 0.02 | true | [T5, row K = 2 pi] |

H4 did not reproduce. At K = pi/2 both L1 and L5 gave 0.256066843848 against the reference 0.298 [T5, row K = pi/2]. The value lies below the reference and outside the tolerance of 0.01. The K = pi/2 start did not reproduce the published value [@SalaArtuso2015]. The K = 2 pi value, 1.17018723872 against 1.172, is within its tolerance of 0.02 [T5, row K = 2 pi]. The reproduction as a whole fails [findings B, H4]. The possibility that the paper's orbit sits in a different component remains open, and the tables do not decide it [findings B, H4]. Each K rests on one start, with no error bar [findings B, H4].

Cross-check at the shorter horizon T = 1e8 [params: h4_T_cross]. This has no verdict under the frozen rule. L5, Float64 against Double64:

| K | Float64 | Double64 | Pointer |
|---|---|---|---|
| 1.57079632679 | 0.299256534631 | 0.299012763029 | [T5, cross-check, row 1.57079632679] |
| 6.28318530718 | 1.17243899416 | 1.17179404609 | [T5, cross-check, row 6.28318530718] |

Comparators (never pass or fail). The columns are the ensemble mean of lambda_L5(N_1) over the kept chaotic-sea orbits, the fit of Shevchenko [@Shevchenko2014], and the line ln(K/2) [@Chirikov1979]. The frozen text designates the line for K = pi/2, 2 pi and 10. The program prints it at all five K [findings B, H4].

| K index | K | Orbits | Mean | Shevchenko fit | ln(K/2) | Pointer |
|---|---|---|---|---|---|---|
| 2 | 0.5 | 40 | 0.0575923614667 | 0.06665 | -1.38629436112 | [T5, comparators, K index 2] |
| 3 | 0.971635406 | 40 | 0.091095497804 | 0.12951899962 | -0.721921822156 | [T5, comparators, K index 3] |
| 4 | 1.57079632679 | 40 | 0.298004458001 | 0.342657808664 | -0.24156447527 | [T5, comparators, K index 4] |
| 5 | 6.28318530718 | 40 | 1.17214446011 | 1.17006018176 | 1.14472988585 | [T5, comparators, K index 5] |
| 6 | 10 | 40 | 1.62000488786 | 1.61943791243 | 1.60943791243 | [T5, comparators, K index 6] |

### H5 (are the labels informative?): SUPPORTED

Preregistered rule. On the ladder subset S_L, H5a is the fraction of S_L orbits whose two labels are identical at Float64, Double64 and BigFloat, and it must be at least 0.90 [params: h5a_min]. H5b is the fraction of all panel orbits with either axis "unresolved", at most 0.30 [params: h5b_max]. H5c is the fraction with structure "unresolved" among orbits with instability "none detected on the tested horizon", at most 0.30 [params: h5c_max]. By the program's reading, H5 is supported only if all three pass [findings B, H5; findings E, reading 5]. If H5b or H5c failed, decisions D36 or D37 would go back to the owner, as a design result and not a study failure [findings B, H5].

Measured. Verdict as printed: supported [T7, verdict line].

| Sub-rule | Orbits | Fraction | Limit | Result | Pointer |
|---|---|---|---|---|---|
| H5a | 92 | 1.0 | at least 0.90 | pass | [T7, row H5a; params: h5a_min] |
| H5b | 883 | 0.0951302378256 | at most 0.30 | pass | [T7, row H5b; params: h5b_max] |
| H5c | 334 | 0.0568862275449 | at most 0.30 | pass | [T7, row H5c; params: h5c_max] |

All three label conditions are met, so no label decision goes back to the owner under these rules [findings B, H5]. H5a is taken before the ladder test acts [findings E, reading 5]. See D3.

## Claims C1 and C2

As the findings file relates them [findings A, findings B]:

- C1 is tested by H1. The preregistered rule for C1 was not met. H1 is falsified [findings B, H1].
- C2 is tested by H2 (finite-time inequivalence) and H3 (regime dependence of the convergence time). H2 is not supported (NEGATIVE) and H3 is supported [findings B, H2, H3].

H4 and H5 are not tied to C1 or C2 in the findings file [findings B, H4, H5].

## Reported without a decision rule

### Cost (F6 and manifest arms)

The cost table gives ns per orbit-step by method and precision, and relative to L5 in the same precision. Each row is an unweighted mean over all items of that method and precision, so the item sets differ between rows [findings C, Cost; findings E, reading 12].

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

The findings file lists the L1 Float64 row without an item count [findings C, Cost]. The count 1917 is from the table row.

Run totals from the manifest: total wall time 11957.993616580963 seconds against a cap of 12.0 hours [manifest: wall.total_seconds, wall.cap_hours]; 1363 chunks and 15569 items in the plan [manifest: plan_counts.chunks, plan_counts.items]; 23 threads used of 23 available on a host with 24 cores [manifest: threads.used, threads.available, host.cores]; device NVIDIA GeForce RTX 5060 Laptop GPU [manifest: device]. For the main arm, the manifest records a cost of 222.15498686174308 ns per orbit-step against a budget of 415.2 ns, with status complete [manifest: arms.main.cost_ns_per_orbit_step, arms.main.budget_ns, arms.main.status]. All 13 arms in the manifest have status complete [manifest: arms.*.status].

### Sensitivity (F7)

The table has 31 rows. The column is the median paired difference of lambda at n_H2 against the primary item (median_dlambda_nH2), with 171 orbits per row except the tau rows [findings C, Sensitivity; F7].

- Direction. L1 and L5 values are of order 1e-9 at offset 0.3 (L5 1.73141166638e-9) [F7, row direction 0.3 L5] and of order 1e-5 or 1e-6 at the other offsets: L5 -3.32612943661e-5 at 0.3 + 2pi/5 [F7, row direction 0.3 + 2pi/5 L5], -3.23438506002e-5 at 0.3 + 3pi/5 [F7, row direction 0.3 + 3pi/5 L5], -3.48457834277e-6 at 0.3 + pi/5 [F7, row direction 0.3 + pi/5 L5]. L9 is 0.0 in every direction row [F7, rows direction L9].
- L2 with eps = 1e-10 gives 2.26005874629e-9 [F7, row eps 1.0e-10], eps = 1e-12 gives 5.77729472453e-5 [F7, row eps 1.0e-12], and eps = 1e-6 gives 4.04116577712e-10 [F7, row eps 1.0e-6].
- Wrapped versus lifted: 0.0 for all eight methods of V [F7, rows wrapped].
- median_log10_nconv_ratio is 0.0 in all 31 rows [F7, column median_log10_nconv_ratio].
- Tau sweep. Tau = 2 and 4 show no breakdown, with 171 orbits for L4 and L5 [F7, rows tau 2, tau 4]. Tau = 8 has 154 orbits for L4 and L5, with breakdown count 17 [F7, rows tau 8]. Tau = 16 has 132 orbits for L4 with breakdown count 39 and 134 orbits for L5 with breakdown count 37 [F7, rows tau 16]. At tau = 16 the breakdown steps are 16 or 32 for L4 and 16 to 288 for L5 [F7, breakdown listing, tau 16].

Symplectic residuals eps_pair(n) and eps_det(n). No numbers for these are quoted in this document. They are shown in figure F5 [findings C].

Figures F1 (lambda_m(n) curves), F2 (R(n) by class and K, with the pooled series over W), F3 (n_conv distributions) and F5 are not reproduced here [findings C].

## Exploratory results

Everything in this section is EXPLORATORY. None is a decision of the confirmatory section, none has a frozen rule, and none is to be read as confirmatory [findings D]. The tables do not give causes, and none is claimed [findings D]. The tables file lists no exploratory items of its own [T Exploratory items: none], so the pointers below name the tables the observations are read from.

**D1 (EXPLORATORY): where the H1c failures sit.** The 47 H1c failures are in three K cells. K index 2 has 13 of 40 failing (27 pass), K index 3 has 27 of 40 (13 pass) and K index 4 has 7 of 40 (33 pass) [T2, rows K index 2 to 4; findings D1]. K index 5 and 6 have no failures (40 of 40 pass each) and the 9 controls all pass [T2, rows K index 5, K index 6, controls]. The failing ids fall in three ranges: 13 ids between 166 and 196, 27 between 346 and 385, and 7 between 549 and 580 [T2, Failing orbits; findings D1]. These counts match the per-cell failure counts, but the tables do not state the orbit-to-cell mapping [findings D1]. The tables do not show why the failures sit at the three smaller K [findings D1].

**D2 (EXPLORATORY): the info controls.** Controls 878, 880, 881 and 882 are reported as info. Their largest deviations are 0.222903869545 [T2, controls, 878], 0.443364672685 [T2, controls, 880], 0.650205817015 [T2, controls, 881] and 0.103506221066 [T2, controls, 882]. These are far above the bound that applies to the other six controls [findings D2]. Under Amendment A5 item 3 these are the controls whose binary64 start is not an exact periodic point. The tables do not say more [findings D2].

**D3 (EXPLORATORY): H5a is 1.0 while many ladder labels are forced.** H5a is 1.0 on 92 S_L orbits, read before the ladder test acts [T7, row H5a]. The ladder test then forced the labels of 35 of those 92 orbits to unresolved, a share of 0.380434782609 [F4, row ladder; T6, forced table]. Both hold at once: the per-precision labels agree, and the ladder test overrides 35 of them [findings D3]. The program prints the unresolved share before the tests as 46 of 883 (0.0520951302378) [F4, row unresolved before the tests] and after as 84 of 883 (0.0951302378256) [F4, row unresolved after the tests]. The second value is the H5b fraction [T7, row H5b]. Other forced tests: direction 7 of 172 (0.0406976744186) [F4, row direction]; estimator, diagnostics and excluded 0 of 883 each [F4, rows estimator, diagnostics, excluded]. The tables do not show how the forced counts (35 and 7) relate to the rise from 46 to 84 [findings D3].

**D4 (EXPLORATORY): H4 at the shorter horizon.** At T = 1e8 the cross-check gives, at K = pi/2, 0.299256534631 (Float64) and 0.299012763029 (Double64) [T5, cross-check, row 1.57079632679]. At K = 2 pi it gives 1.17243899416 and 1.17179404609 [T5, cross-check, row 6.28318530718]. The Float64 value at K = pi/2 at T = 1e9 in the verdict table is 0.256066843848 [T5, row K = pi/2]. The ensemble mean of lambda_L5(N_1) over the 40 kept chaotic-sea orbits at K = pi/2 is 0.298004458001 [T5, comparators, K index 4]. The tables do not show why the single start at K = pi/2 gives a different value at T = 1e9 than at T = 1e8 [findings D4].

**D5 (EXPLORATORY): unresolved labels by cell.** "Both" counts orbits with both axes unresolved (65 in all). "Structure only" counts orbits with instability "none detected on the tested horizon" and structure unresolved (19 in all). The two add to the 84 of H5b [findings D5; T6, label table]. Cells not listed have none.

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

K index 1, K index 6 and K index 5 classes C and X have no unresolved orbit [findings D5; T6, rows 1, 6, 5 C, 5 X].

**D6 (EXPLORATORY): short, empty and not decidable cells.** The panel manifest marks two cells as short ("fewer qualifying grid cells than n_c"): K index 2 class S (25 drawn, 25 qualifying) and K index 5 class S (5 drawn, 5 qualifying) [T1, manifest cells table, rows 2 S, 5 S]. It marks as empty K index 5 class R (no elliptic periodic point) [T1, manifest cells, row 5 R] and K index 6 classes R (no elliptic periodic point), T and S (no qualifying grid cell) [T1, manifest cells, rows 6 R, 6 T, 6 S]. It marks K index 1 classes C, S and X as undefined (class undefined at K = 0) [T1, manifest cells, rows 1 C, 1 S, 1 X]. Class X shows 0 qualifying grid cells at every K index from 2 to 6 [T1, manifest cells, rows 2 X to 6 X] because the preregistered selection rule does not draw class X from the screening and places the starts near the hyperbolic fixed point [findings D6]. In T3 and T4 only class S at K index 5 is not decidable, with 1 orbit [T3, row S 5; T4, row S 5].

**D7 (EXPLORATORY): cell patterns.** The largest cell median of R(n_H2) is 0.000775501992792 (class S, K index 3) [T3, row S 3] and the smallest is 5.82402750446e-12 (class X, K index 2) [T3, row X 2]. The cell medians of n_conv are between 4.223954e6 [T4, row C 3] and 7.079458e6 [T4, rows C 2, S 2, X 2] at K index 2 and 3 (N_1 = 1e7 [params: N_1]), are 794328.0 (class S) [T4, row S 4] and 63096.0 (class X) [T4, row X 4] at K index 4, and lie between 158.0 [T4, row X 6] and 266.5 [T4, row C 5] at K index 5 and 6, excluding class S at K index 5 [findings D7]. The tables do not show why [findings D7].

## Deviations and readings

Deviations from the frozen text, as the program states them [findings E, Deviations]:

1. L2 runs in wrapped coordinates (Amendment A5 item 1). H1c compares the wrapped L2 value with the lifted L5 value of the same orbit. This affects the H1c metric, the sub-rule that failed. The tables do not show whether this accounts for any part of the H1c result [findings E, deviation 1; T Deviations].
2. The tau sweep records an orthonormalisation breakdown as a result, not as a failure (A5 item 2). The number of orbits with the flag is 40 [T Deviations; findings E, deviation 2]. This affects the tau rows of F7, which have fewer than 171 orbits at tau = 8 and 16 [F7, rows tau 8, tau 16].
3. Control bound scope (A5 item 3). The bound applies to controls 871, 873, 875, 877, 879 and 883. Controls 878, 880, 881 and 882 are reported as info [T Deviations; T2, controls table].
4. The control bound of L9 is 5.0e-8 and that of the other methods 1.0e-8 (A6) [params: control_bound, control_bound_L9; T Deviations].
5. Overrides recorded in the manifest: none [manifest: overrides; T Deviations].

Readings where the frozen text is silent are listed in [findings E, Readings 1 to 14] and in the tables file, section Readings. Two readings bear directly on the verdicts above: H5a is read from the per-precision labels before the ladder test acts, and H5 is supported only if H5a, H5b and H5c all pass [findings E, reading 5]; and S is the class C orbits at K = 2 pi and K = 10 without a label condition, W takes class C at 0 < K <= K_c, and E_ch takes class C at K > 0 [findings E, reading 8].

Package version. The TOOL value type names the tagged release v0.1.0, which is cut at milestone M3 exit. This run was made with package version 0.0.1, at the commit given in the header, before the tag [findings E, Note; manifest: versions.module].

## Limitations

From the findings file [findings F]:

- One draw of the panel (a single seed) with 883 orbits [T1]. The bootstrap intervals resample orbits within that draw and do not reflect a second draw.
- The panel is run in Float64 only. Precision comparisons come from the ladder subset S_L of 92 orbits (Float64, Double64, BigFloat), and H5a rests on that subset [T7, row H5a].
- Sets as printed: E_ch 209 [T2, header line], W 355 and S 80 [T3, T4, header lines], S_L 92 [T7, row H5a], 883 for H5b [T7, row H5b], 334 for H5c [T7, row H5c]. F7 rows use 171 orbits except the tau rows, which use 132 to 171 [F7].
- One orbit (881) is dropped as drift [T1, last line]. The controls table lists it as an info control [T2, controls, 881].
- Cells with fewer than 20 orbits: class S at K index 5 (5 drawn, 1 in W) is not decidable in T3 and T4 [T1, manifest cells, row 5 S; T3, row S 5; T4, row S 5]. Class S at K index 2 has 25 drawn and 21 in W [T1, manifest cells, row 2 S; T3, row S 2]. Cell-level medians and intervals are descriptive, and the decisions use the pooled sets.
- H4 rests on one start per K, as the frozen text specifies, with no error bar [findings B, H4].
- No numbers for the symplectic residual figure F5 are quoted in this document.
- The numbers were checked against the tables file of the run. The analysis program was verified against an independent recomputation (`analysis-verification.md`) [findings F]. The numbers were not recomputed from the stored series for this document. Table cosmetics (counts printing as 2.0 or 166.0) are formatting only [findings A, F].

Open items carried over from the stage 04 and 05 hand-off:

- The STRONG_K constant was kept by owner decision.
- Orbit ids and the K index print as `2.0` or `166.0` in some tables (formatting only) [T2, Failing orbits; T5, comparators].
- Manifest.toml for Julia 1.10 is unresolved.
- Figures F1 to F7 were checked only by file size and the `%PDF-` header. They were not checked visually.

## Uncited claims

- The value K_c = 0.971635406 is used as given in the tables [T5, comparators, K index 3; params: K_c in Header]. The findings file gives no citation key for it.
- The definitions of the methods L1 to L9 and Q are not cited here, because the findings file gives no key for them.
