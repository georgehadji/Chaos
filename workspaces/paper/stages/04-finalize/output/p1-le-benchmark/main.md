# Largest Lyapunov exponent estimators on the standard map: a preregistered numerical study



## Abstract

We report a preregistered numerical study of largest-Lyapunov-exponent estimators on the standard map. 
Scope: one panel from a single draw (a single seed), run in Float64, with five nonzero kick values and one start per kick value in the reproduction test (H4). Section 10 gives the limitations. 
The panel had 883 orbits: 13 exact-exponent controls, orbits at K = 0, and orbits at five nonzero kick values. 
The decision rules were frozen before any data existed. 
The claim of asymptotic equivalence (H1) was falsified. The failure is carried by sub-rule H1c alone, with 162 of 209 eligible orbits passing, a pass fraction of 0.775119617225 against a limit of 0.90. Sub-rules H1a and H1b were met in full. H1 is evaluated in Float64 at the horizon N_1 = 1e7, and H1c compares L2 run in wrapped coordinates with L5 in lifted coordinates (limitation 4, Section 10). 
The claim of finite-time inequivalence (H2) was not supported (NEGATIVE). Over the 355 orbits of the preregistered set W the median of R(n_H2) was 2.23588035691e-5 against a threshold of 0.10. Set W is defined in Section 2.2. 
The regime dependence of the convergence time (H3) was supported, with a ratio of medians of 11213.7767857. 
The reproduction of the reference values of Sala and Artuso [@SalaArtuso2015] (H4) did not reproduce at K = pi/2, where we obtained 0.256066843848 against 0.298. It passed at K = 2 pi, where we obtained 1.17018723872 against 1.172. 
The label rule H5 was supported under its three sub-rules (H5a, H5b, H5c). The ladder test forced 35 of the 92 ladder-subset labels to unresolved; this is reported as exploratory item D3. 
Cost and sensitivity are reported without a decision rule. 
The run used package version 0.0.1, before the tagged release v0.1.0. The code is released as Chaoslab v0.1.2 (TOOL), a tagged and documented release with the same computing code (Section 11). 

## 1. Introduction

This paper reports a preregistered numerical study of largest-Lyapunov-exponent estimators on the standard map. The study tests two claims of the project, C1 (asymptotic equivalence) and C2 (finite-time inequivalence). The claim texts below are quoted from the project file `_config/research-question.md`. They are claims of the study. They are not statements supported by the literature cited in this paper. 

Claim C1, as stated in `_config/research-question.md`: "Benettin, QR (Gram-Schmidt / Householder), stabilized SVD and well-conditioned tangent-map methods agree on lambda_1 when run long enough at adequate precision." 
Claim C2, as stated in `_config/research-question.md`: "Finite-time LE estimators differ substantially over realistic budgets, most of all in weak chaos, sticky layers and near separatrices." Here LE stands for Lyapunov exponent and lambda_1 for the largest exponent. 
The tests of C2 in this study run over the preregistered set W of Section 2.2, whose composition is shown in Tables 6 and 8. 

The study also attempts to reproduce two published values. Sala and Artuso give 0.298 at K = pi/2 and 1.172 at K = 2 pi [@SalaArtuso2015]. 

Two published comparators serve for the chaotic-sea exponent. These are the fit of Shevchenko [@Shevchenko2014] and the line ln(K/2) of Chirikov [@Chirikov1979]. The study uses both as comparators that never pass or fail. 

Lyapunov spectra are usually computed by renormalising tangent vectors and orthogonalising them by Gram-Schmidt [@Benettin1980a; @Benettin1980b; @ShimadaNagashima1979]. Gram-Schmidt, Householder and continuous QR variants were compared on a Toda chain; the discrete variants differed by much less than the fluctuations of the finite-time values, and continuous QR took about 40 times the CPU time of Gram-Schmidt [@Geist1990]. Singular-value methods take the singular values of a product of Jacobians without forming the product [@Stewart1997], and their theory holds for systems with stable and distinct exponents [@DieciElia2006].

Finite-time effects are known in particular settings. For standard-like maps the finite-time exponent carries an initial-vector transient that fades like 1/n [@SalaArtuso2015]. For the two-trajectory method on flows, the renormalisation time and the initial distance can produce a spurious largest exponent [@DubeibeBermudezAlmanza2014]. Sticky orbits need long runs: a finite-time exponent in an extreme sticky zone of the standard map had not converged after 5 x 10^7 iterations, against about 10^7 in the chaotic sea [@HarsoulaKaramanosContopoulos2019], and one orbit at k = 0.975 looks confined for 10^6 steps and escapes after 2.3 x 10^9 [@Meiss2015]. Round-off induced divergence of orbits has been related to the reversibility error in discrete-time systems [@FarandaMestreTurchetti2012]. This study compares the estimators above on one map, over five orbit regimes, with a precision ladder on a subset of the orbits.

This paper reports the study as preregistered. It adds no new claim, softens no result and offers no reinterpretation. Section 2 gives the design, the notation and the run identity. Sections 3 to 7 give the confirmatory results for hypotheses H1 to H5. Exploratory subsections inside those sections are marked EXPLORATORY: D1 and D2 in Section 3.5, D4 in Section 6.4, and D3 and D5 in Section 7.3. Section 8 gives cost and sensitivity. Section 9 gives the exploratory items D6 and D7. Sections 10 to 12 give limitations, availability and the conclusion.

## 2. Methods and preregistration

### 2.1 Preregistration

The decision rules were frozen on 2026-10-17, before any data existed. Amendments A1 to A6 were all frozen before the first study data file. 
The rules are applied as frozen, with the deviations listed in Section 2.5. The results document also lists 14 readings where the frozen text is silent (Section 2.5). 
Amendments A1 to A4 are described here by pointer to the frozen text. A1 sets the hash rule, the T predicate, the screening direction, the windowed estimate, the Newton periodic-point search, the subset caps and the machine-readable panel inputs [prereg: Amendment A1, items 1 to 7]. A2 handles a ghost pair of periodic points at K = 2 pi and is superseded by A3 [prereg: Amendment A2, Amendment A3]. A3 sets the acceptance test of the periodic-point search [prereg: Amendment A3]. A4 adds the machine-readable run block, the setup check, the cost cap and pilot, the backend agreement and the scope of the control bound [prereg: Amendment A4, items 1 to 5]. The content of A5 and A6 is given in Section 2.5. 
The preregistered value types were FALSIFIES (C1, C2), REPRODUCES ([@SalaArtuso2015]), TOOL and NEGATIVE allowed. Claim C3 is out of scope. 
TOOL claim: the library Chaoslab, release v0.1.2. Every method the run used has catalog status `verified` (preregistration, Methods). Release v0.1.2 is tagged, carries user documentation and is archived on Zenodo; its package code differs from the analysis commit only in docstrings, in a file-descriptor conversion that is a no-op on Julia 1.13.2 (the Julia version of the run) and in the figure drawing (Section 11). 
The preregistration header is identified by the SHA-256 value in Table 1. The location of the frozen text is given in Section 11. 

### 2.2 Notation and definitions

This subsection restates the notation of the study design. Each item has a pointer to the frozen text. The definition of the standard map and of its tangent map is in `_config/research-question.md` and is not repeated here (pointer only). 

- Variables and horizons. The map variables are (x, p). The H4 start is x0 = 1e-3 and p0 = -1e-3. K is the kick parameter. The horizons are N_1 = 1e7 (all orbits), N_long = 1e8 (designated subset) and N_B = 1e6 (precision ladder). N is the horizon of the series used; R(n) and n_conv use the horizon-N_1 series of every orbit. 
- Finite-time exponent. lambda_m(n) is the finite-time largest Lyapunov exponent of method m after n map steps, in natural log per step, with the Euclidean 2-norm on the lifted (x, p) tangent vector. L9 uses the largest finite-time singular value of the product. lambda_ref of an orbit is lambda_L5 at the longest Float64 horizon of that orbit. 
- Methods. L1 is Benettin. L2 is the two-trajectory method, an approximation of the largest exponent. L3 is CGS, L4 is MGS, L4b is RGS, L5 is Householder (the default reference), L5b is the H2 diagonal iteration path, L6 is Givens, and L9 is the stabilized SVD (d = 2). Their algorithm details are in the study design and the repository, not here. 
- Method sets. Q = {L3, L4, L4b, L5, L5b, L6}. The variational set is V = {L1} + Q + {L9}. The expression "{L1} + Q" means L1 together with the members of Q. 
- Parameters. tau is the fixed reorthonormalization interval of the Q methods (primary tau = 1). eps is the perturbation size of L2 (primary eps_0 = 1e-8). The primary initial deviation direction is the angle pi/4 from the x axis. On the designated subset four more "direction offsets" are used: 0.3, 0.3 + pi/5, 0.3 + 2 pi/5 and 0.3 + 3 pi/5. The precisions are Float64, Double64 and BigFloat. 
- Spreads and ratios. Delta_QB(n) is the maximum minus the minimum of lambda_m(n) over {L1} + Q. Delta_V(n) is the same over V. R(n) = Delta_V(n) / lambda_ref (H2). n_conv is defined in H3 (Section 5.1); the series are stored at recorded values of n, so n_conv takes recorded values of n. 
- Design classes. T (KAM torus), R (resonant), C (chaotic sea), S (sticky) and X (near a separatrix or a hyperbolic fixed point). Design class is the a-priori intended class of an orbit. A cell is a K value times a design class. 
- Controls. The controls are the exact periodic points of the map and one hyperbolic linear map (Y5), whose exponents are known analytically. "Info" controls are defined in Section 2.5. 
- Sets. E_ch is the set of design class C orbits at K > 0 plus the controls with positive analytic exponent (H1). The weak-chaos set W is defined below. The strong-chaos set S is the class C orbits at K = 2 pi and K = 10 (H3). S_D is the designated subset and S_L is the ladder subset, nested in the panel; S_L carries the Double64 and BigFloat runs. The symbol S names three things in this paper: the panel class S (sticky), the strong-chaos set S of H3 and, as S_L, the ladder subset. The context says which. 
- Definition of W (preregistered). W = orbits with instability label "positive growth supported" in the design classes sticky or near-separatrix at any K > 0, plus chaotic-sea orbits at K <= K_c. 
  The name "weak-chaos set" is the preregistration's. By this definition W contains class S and class X orbits at the large K values as well, and Tables 6 and 8 show its cells. The name does not mean that every member of W is in weak chaos. 
- Labels. Each orbit carries two labels. The instability label takes the values "positive growth supported", "none detected on the tested horizon" and "unresolved". The structure label takes the values "periodic", "quasiperiodic", "irregular" and "unresolved". Four disagreement tests force a label to "unresolved": the precision ladder test (on S_L only), the estimators test (any failure of H1a or H1b on that orbit), the direction test (on the designated subset only) and the diagnostics test (a pairing, determinant or orthogonality residual above its tolerance for L5). A label is "forced" when one of these tests sets it to "unresolved". The thresholds of the label rules are in the frozen text. 
- Drift and kept orbits. An orbit is flagged "drift" when its lifted |p| exceeds 1e6 at any recorded time. Drift orbits are removed from E_ch, W and S. "Kept" orbits are the orbits not dropped as drift or excluded. 
- Bootstrap. The bootstrap intervals are percentile intervals at the 95 percent level with 10^4 resamples. Orbits are resampled with replacement, and W and S are resampled independently. 

### 2.3 Panel and sets

The study compares estimators on a panel of 883 orbits. The panel table [T1] gives its cells: 13 controls in the cell of K index 0; 40 orbits each in classes R and T at K index 1; and the cells at K index 2 to 6. At K index 2 the classes are C, R, S, T and X with 40, 40, 25, 40 and 40 orbits. At K index 3 and 4 the classes are C, R, S, T and X with 40 orbits each. At K index 5 the classes are C, S, T and X with 40, 5, 40 and 40 orbits. At K index 6 the classes are C and X with 40 orbits each. Cells that are short, empty or undefined are listed in D6 (Section 9). 
K index 0 holds the exact-exponent controls. K index 1 is K = 0. K index 2 to 6 are K = 0.5, 0.971635406 (K_c, the breakup value of the golden-mean KAM torus, as given in [@FrahmShepelyansky2010; @CarlettiLaskar2000] and fixed in the preregistration), 1.57079632679 (pi/2), 6.28318530718 (2 pi) and 10. 
"Class S" is a panel class. It is not the "strong-chaos set S" of H3, which is the class C orbits at K = 2 pi and K = 10. 
One orbit, orbit 881, is dropped as drift. Under the preregistration, drift orbits are removed from E_ch, W and S. H5b is defined on all panel orbits, and the results document gives 883 for it. 
The sets used by the decision rules have these sizes: E_ch 209, W 355, S 80 and S_L 92. H5b uses 883 orbits and H5c uses 334. 
The program reads S as the class C orbits at K = 2 pi and K = 10 without a label condition. It reads E_ch as class C at K > 0 (plus the controls, Section 3.1). 
For W this paper, the results document and the findings file use the preregistered definition of Section 2.2. The analysis program prints a shorter phrase for W in its list of readings, class C at 0 < K <= K_c [findings E, reading 8]. That phrase is not the frozen definition and does not describe the cells of Tables 6 and 8, which also contain classes S and X (limitation 6, Section 10). 
The estimator methods L1 to L9 and Q are defined in the study design (Section 2.2). The repository (Section 11) holds the details. 

### 2.4 Run identity

All reported numbers come from one fully automatic analysis run over the stored series. 
The run completed with 1363 chunks, all with status ok. 
Table 1 identifies the run. 

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



The tables header and the findings file give the first git commit. The run manifest records the second. The run commit is an ancestor of the analysis commit, and between them only the analysis code, its figure extension, the command-line entry point and the project file change. The stored series therefore come from the run commit and the tables and figures from the analysis commit, with the same series code at both (limitation 8, Section 10). 

### 2.5 Deviations from the frozen text and readings

The frozen text was applied with the deviations that the program states. 

1. L2 runs in wrapped coordinates (Amendment A5 item 1). H1c compares the wrapped L2 value with the lifted L5 value of the same orbit. 
2. The tau sweep records an orthonormalisation breakdown as a result, not as a failure (A5 item 2). 
3. The control bound applies to controls 871, 873, 875, 877, 879 and 883 (A5 item 3). Controls 878, 880, 881 and 882 are reported as info and are not held to the bound. 
4. The control bound of L9 is 5.0e-8 and that of the other methods is 1.0e-8 (A6). 
5. The manifest records no overrides. 

The results document lists 14 readings where the frozen text is silent [findings E, Readings 1 to 14; T Readings]. It names two that bear directly on the verdicts. Reading 5: H5a is read from the per-precision labels before the ladder test acts, and H5 is supported only if H5a, H5b and H5c all pass. Reading 8: the definitions of S, W and E_ch as the program reads them (Section 2.3). The other twelve readings are listed in the tables file, section Readings, and are not described in this paper. 

## 3. Results: H1 (claim C1), FALSIFIED

### 3.1 Preregistered rule

The rule applies on E_ch, which is class C orbits at K > 0 plus the exact-exponent controls. The horizon is N_1 = 1e7. 
An orbit passes H1a if the spread over {L1} + Q is at most 1e-8. It passes H1b if |lambda_L9 - lambda_L5| is at most 1e-6. It passes H1c if |lambda_L2 - lambda_L5| is at most 5e-3 at eps_0 = 1e-8. 
H1 is SUPPORTED only if the pass fraction is at least 0.98 for H1a and H1b and at least 0.90 for H1c. It is FALSIFIED if any sub-rule is below its limit. 

### 3.2 Measured

There were 209 eligible orbits. The verdict as printed is falsified. Table 2 gives the three sub-rules. 

Table 2. H1 sub-rules on E_ch.

| Sub-rule | Limit | Passing orbits | Fraction | Result | Pointer |
|---|---|---|---|---|---|
| H1a | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1a_fraction] |
| H1b | 0.98 | 209 of 209 | 1.0 | pass | [T2, header line; params: h1b_fraction] |
| H1c | 0.90 | 162 of 209 | 0.775119617225 | fail | [T2, header line; params: h1c_fraction] |



H1 is falsified. The H1c pass fraction was 162 of 209, which is 0.775119617225, against the limit of 0.90. 
The shortfall is carried by sub-rule H1c alone. H1a and H1b were met in full. 

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



The 47 failing orbits sit in three K cells. At K index 2, 27 of 40 orbits pass. At K index 3, 13 of 40 pass. At K index 4, 33 of 40 pass. K index 5 and 6 pass 40 of 40, and the 9 controls pass. 
H1a and H1b pass in every cell. 
The ids of the failing orbits and the extreme deviations are in [T2, Failing orbits]. Every failing orbit has |lambda_L2 - lambda_L5| above 0.005. The smallest is 0.00548462883667 and the largest is 0.0353748064771. 
On the failing orbits the largest Delta_QB is 9.66310365058e-14 and the largest |lambda_L9 - lambda_L5| is 3.58331395889e-7. None of them breaks H1a or H1b. 

### 3.3 Exact-exponent controls

The controls under the bound are 871, 873, 875, 877, 879 and 883. Table 4 gives the largest deviation over the methods. 

Table 4. Exact-exponent controls under the bound.

| Control | Scope | Largest deviation | L9 deviation | Pointer |
|---|---|---|---|---|
| 871 | bound | 5.26802590439e-9 | 5.26802579337e-9 | [T2, controls, 871] |
| 873 | bound | 2.48736364838e-9 | 8.326783706990001e-12 | [T2, controls, 873] |
| 875 | bound | 5.16177345311e-9 | 1.82781789704e-9 | [T2, controls, 875] |
| 877 | bound | 1.79536172595e-8 | 1.79536172595e-8 | [T2, controls, 877] |
| 879 | bound | 2.28260139679e-8 | 2.28260139679e-8 | [T2, controls, 879] |
| 883 | bound | 1.11022302463e-16 | 1.11022302463e-16 | [T2, controls, 883] |



The largest deviations lie within the bounds of Amendment A6. For controls 877 and 879 the largest deviation is the L9 deviation. It is above 1e-8 and below the L9 control bound of 5.0e-8. 
The table does not give the other methods' deviations on these two controls separately. 

Accounting for the control ids. The panel table counts 13 controls in the cell of K index 0, and the H1 cell table counts 9 controls in E_ch. The controls table lists 10 ids: the six of Table 4 and the four info controls 878, 880, 881 and 882 (D2). The other 3 controls, orbits 872, 874 and 876, start at the elliptic fixed point (pi, 0) at K = 0.5, K_c and pi/2 and have analytic exponent 0. E_ch and the controls table take only controls with positive analytic exponent, so these three appear in neither; they serve the setup check of the Double64 and BigFloat runs and enter H5b with every panel orbit. 
Orbit 881 is an info control and is the orbit dropped as drift (Section 2.3). 

### 3.4 Reading

The preregistered rule for C1 is not met. H1 is reported as falsified. 
The falsification rests on sub-rule H1c alone. It was evaluated in Float64 at N_1 = 1e7 on one panel, and H1c compares L2 in wrapped coordinates with L5 in lifted coordinates (Section 2.5, item 1; limitation 4 in Section 10). The tables do not show whether the coordinate choice accounts for any part of the H1c result. 

### 3.5 Exploratory items D1 and D2 (EXPLORATORY)

Everything in this subsection is EXPLORATORY. No frozen rule applies and no cause is claimed. 

D1 (EXPLORATORY): where the H1c failures sit. The 47 H1c failures are in three K cells. K index 2 has 13 of 40 failing (27 pass). K index 3 has 27 of 40 (13 pass). K index 4 has 7 of 40 (33 pass). 
K index 5 and 6 have no failures (40 of 40 pass each) and the 9 controls all pass. 
The failing ids fall in three ranges: 13 ids between 166 and 196, 27 between 346 and 385, and 7 between 549 and 580. 
These counts match the per-cell failure counts, but the tables do not state the orbit-to-cell mapping. The tables do not show why the failures sit at the three smaller K. 

D2 (EXPLORATORY): the info controls. Controls 878, 880, 881 and 882 are reported as info. Their largest deviations are 0.222903869545, 0.443364672685, 0.650205817015 and 0.103506221066. 
These are far above the bound that applies to the other six controls. Amendment A5 item 3 designates these four as the controls reported as info; the results document describes them as the controls whose binary64 start is not an exact periodic point. That description is not tested here. The tables do not say more. 

## 4. Results: H2 (claim C2), NOT SUPPORTED (NEGATIVE)

### 4.1 Preregistered rule

R(n) = Delta_V(n) / lambda_ref is evaluated at n_H2 = 1e4 steps over the set W defined in Section 2.2. The preregistration names W the weak-chaos set. 
W = orbits with instability label "positive growth supported" in the design classes sticky or near-separatrix at any K > 0, plus chaotic-sea orbits at K <= K_c. 
H2 is SUPPORTED if the median of R(n_H2) is at least 0.10 and the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) is at least 0.05. Otherwise the null is reported. 

### 4.2 Measured

W has 355 orbits. The median of R(n_H2) is 2.23588035691e-5. The 95 percent bootstrap interval is [1.9257235885e-5, 2.31226443693e-5]. The verdict as printed is not supported. 
Both conditions fail. The median is below 0.10 and the lower end is below 0.05. 

Table 5. H2 against the preregistered thresholds.

| Quantity | Measured | Threshold | Result | Pointer |
|---|---|---|---|---|
| Median of R(n_H2) over W (355 orbits) | 2.23588035691e-5 | at least 0.10 | fail | [T3, header line; params: h2_median] |
| Lower end of the 95 percent interval | 1.9257235885e-5 | at least 0.05 | fail | [T3, header line; params: h2_interval] |



This is a null (NEGATIVE) result. It has the same standing as a pass. 
In W, the relative spread among the methods of V at 1e4 steps is far below the preregistered level. 
The frozen rule decides on the pooled W. Table 6 gives the cells of W; the cell results are descriptive. W is composed of class C at K index 2 and 3, class S at K index 2 to 5 and class X at K index 2 to 6. 

Table 6. H2: R(n_H2) by cell of W (descriptive).

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



No cell median of R(n_H2) reaches the threshold of 0.10. The largest cell median is 0.000775501992792 (class S, K index 3). The direction of the NEGATIVE verdict therefore holds in every cell as well as in the pooled set. This is exploratory item D7 (Section 9) and is not a decision rule. 
The cell of class S at K index 5 has 1 orbit, below the minimum set size of 20. It is not decidable and is neither pass nor fail. 
The bootstrap resamples orbits as if they were independent draws (limitation 1, Section 10). The cell intervals of Table 6 are descriptive. 
Figure F2 shows R(n) by class and K, with the pooled series over W. 

## 5. Results: H3 (claim C2), SUPPORTED

### 5.1 Preregistered rule

n_conv is the smallest n after which every method stays within a 0.05 relative band of lambda_ref up to N. Per orbit it is the maximum over the methods of V. 
H3 is SUPPORTED if the ratio of the median n_conv over W to the median n_conv over the strong-chaos set S is at least 10, and the lower end of the 95 percent bootstrap interval of the ratio is above 3. Here W is the set of Section 2.2 and S is the class C orbits at K = 2 pi and K = 10. 
Orbits whose n_conv is undefined are given n_conv = N, and their share is reported. 

### 5.2 Measured

Table 7 gives the H3 measured values.

Table 7. H3 measured values.

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



H3 is supported. The verdict as printed is supported. Both conditions pass. The ratio is at least 10 and the lower end is above 3. 
The strong-chaos set S of H3 is the class C orbits at K = 2 pi and K = 10. It is not the panel class S. 
The orbits with undefined n_conv enter the medians at the value N. 
Table 8 gives n_conv by cell. The rows of class C at K index 5 and 6 are the orbits of S. The other rows match the cells of Table 6 in class, K index and orbit count. The cell-level medians are descriptive. The class S cell at K index 5 has 1 orbit and is not decidable, as in H2. 

Table 8. H3: n_conv by cell (descriptive).

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



Figure F3 shows the n_conv distributions. 

## 6. Results: H4 (REPRODUCES, [@SalaArtuso2015]), DID NOT REPRODUCE

### 6.1 Preregistered rule

The start is x0 = 1e-3 and p0 = -1e-3, with T = 1e9 iterations. 
At K = pi/2 (reference 0.298) both L1 and L5 (Float64, tau = 1) must be within 0.01. At K = 2 pi (reference 1.172) both must be within 0.02. 
The reference values are those of Sala and Artuso [@SalaArtuso2015]. H4 is SUPPORTED if both K pass. A failure is reported as "did not reproduce" with the numbers. 

### 6.2 Measured

The verdict as printed is did not reproduce. Table 9 gives the verdict table. 

Table 9. H4 verdict table.

| Start | L1 | L5 | Reference | Tolerance | Pass | Pointer |
|---|---|---|---|---|---|---|
| K = pi/2 | 0.256066843848 | 0.256066843848 | 0.298 | 0.01 | false | [T5, row K = pi/2] |
| K = 2 pi | 1.17018723872 | 1.17018723872 | 1.172 | 0.02 | true | [T5, row K = 2 pi] |



H4 did not reproduce. At K = pi/2 both L1 and L5 gave 0.256066843848 against the reference 0.298. The value lies below the reference and outside the tolerance of 0.01. 
The K = pi/2 start did not reproduce the published value [@SalaArtuso2015]. 
At K = 2 pi the value 1.17018723872 against 1.172 is within the tolerance of 0.02. 
The reproduction as a whole fails. 
The possibility that the orbit of the published paper sits in a different component remains open. The tables do not decide it. Each K rests on one start, with no error bar. 

### 6.3 Comparators (no verdict)

Table 10 gives three comparators at the five K values: the ensemble mean of lambda_L5(N_1) over the kept chaotic-sea orbits, the fit of Shevchenko [@Shevchenko2014], and the line ln(K/2) [@Chirikov1979]. They never pass or fail. 
The frozen text designates the line for K = pi/2, 2 pi and 10. The program prints it at all five K. The rows at K = 0.5 and K = K_c lie outside the designated K and are shown as printed; the line takes negative values there (-1.38629436112 and -0.721921822156). The range of applicability of the fit and of the line is not checked in this paper. No verdict is attached to any comparator row. 

Table 10. Comparators.

| K index | K | Orbits | Mean | Shevchenko fit | ln(K/2) | Pointer |
|---|---|---|---|---|---|---|
| 2 | 0.5 | 40 | 0.0575923614667 | 0.06665 | -1.38629436112 | [T5, comparators, K index 2] |
| 3 | 0.971635406 | 40 | 0.091095497804 | 0.12951899962 | -0.721921822156 | [T5, comparators, K index 3] |
| 4 | 1.57079632679 | 40 | 0.298004458001 | 0.342657808664 | -0.24156447527 | [T5, comparators, K index 4] |
| 5 | 6.28318530718 | 40 | 1.17214446011 | 1.17006018176 | 1.14472988585 | [T5, comparators, K index 5] |
| 6 | 10 | 40 | 1.62000488786 | 1.61943791243 | 1.60943791243 | [T5, comparators, K index 6] |



### 6.4 Exploratory item D4 (EXPLORATORY)

D4 (EXPLORATORY): H4 at the shorter horizon. The cross-check uses T = 1e8. It has no verdict under the frozen rule. D4 does not change the H4 verdict, which stands under the frozen rule: H4 did not reproduce. 
At K = pi/2 it gives 0.299256534631 (Float64) and 0.299012763029 (Double64). At K = 2 pi it gives 1.17243899416 and 1.17179404609. 
The Float64 value at K = pi/2 at T = 1e9 in the verdict table is 0.256066843848. The ensemble mean of lambda_L5(N_1) over the 40 kept chaotic-sea orbits at K = pi/2 is 0.298004458001. 
The tables do not show why the single start at K = pi/2 gives a different value at T = 1e9 than at T = 1e8. No cause is claimed. 

## 7. Results: H5, SUPPORTED

### 7.1 Preregistered rule

On the ladder subset S_L, H5a is the fraction of S_L orbits whose two labels are identical at Float64, Double64 and BigFloat. It must be at least 0.90. 
H5b is the fraction of all panel orbits with either axis "unresolved". It must be at most 0.30. 
H5c is the fraction with structure "unresolved" among orbits with instability "none detected on the tested horizon". It must be at most 0.30. 
H5 is supported only if all three sub-rules pass. H5a is read from the per-precision labels before the ladder test acts. 

### 7.2 Measured

The verdict as printed is supported. Table 11 gives the three sub-rules. 

Table 11. H5 sub-rules.

| Sub-rule | Orbits | Fraction | Limit | Result | Pointer |
|---|---|---|---|---|---|
| H5a | 92 | 1.0 | at least 0.90 | pass | [T7, row H5a; params: h5a_min] |
| H5b | 883 | 0.0951302378256 | at most 0.30 | pass | [T7, row H5b; params: h5b_max] |
| H5c | 334 | 0.0568862275449 | at most 0.30 | pass | [T7, row H5c; params: h5c_max] |



All three label conditions are met. H5 is supported under its three sub-rules. 

### 7.3 Exploratory items D3 and D5 (EXPLORATORY)

D3 (EXPLORATORY): H5a is 1.0 while many ladder labels are forced. H5a is 1.0 on 92 S_L orbits, read before the ladder test acts. 
The ladder test then forced the labels of 35 of those 92 orbits to unresolved, a share of 0.380434782609. 
Both hold at once: the per-label agreement across precisions, and the override of 35 labels by the ladder test. 
The program prints the unresolved share before the tests as 46 of 883 (0.0520951302378) and after as 84 of 883 (0.0951302378256). The second value is the H5b fraction. 
Other forced tests: direction 7 of 172 (0.0406976744186). Estimator, diagnostics and excluded show 0 of 883 each. 
The tables do not show how the forced counts (35 and 7) relate to the rise from 46 to 84. No cause is claimed. 

D5 (EXPLORATORY): unresolved labels by cell. "Both" counts orbits with both axes unresolved (65 in all). "Structure only" counts orbits with instability "none detected on the tested horizon" and structure unresolved (19 in all). The two add to the 84 of H5b. 
Cells not listed in Table 12 have none. 

Table 12. Unresolved labels by cell (EXPLORATORY).

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



K index 1, K index 6 and K index 5 classes C and X have no unresolved orbit. 

## 8. Cost and sensitivity (reported without a decision rule)

### 8.1 Cost

Table 13 gives ns per orbit-step by method and precision, and relative to L5 in the same precision. Each row is an unweighted mean over all items of that method and precision, so the item sets differ between rows. 
Relative to L5 in Float64 (168.60570837 ns per orbit-step), L1 costs 0.701549448499, L2 0.503966109187, L4 0.982175001597, L9 1.23876330018 and L6 2.10273197634. These ratios compare rows with different item sets (1917, 1571, 1915, 1915 and 1227 items against 2605 for L5), and they come from a single run that used 23 of 23 available threads, with no repeated timing reported. 
L5 costs 4882.6924383 ns in Double64 and 377399.856865 ns in BigFloat(256). The Double64 and BigFloat rows have 186 and 92 items against 2605 for L5 in Float64, so they are also not like-for-like with the Float64 rows. 

Table 13. Cost per orbit-step.

| Method | Precision | Items | ns per orbit-step | Relative to L5 | Pointer |
|---|---|---|---|---|---|
| L5 | Float64 | 2605 | 168.60570837 | 1.0 | [cost section, row L5 Float64] |
| L1 | Float64 | 1917 | 118.285241721 | 0.701549448499 | [cost section, row L1 Float64] |
| L2 | Float64 | 1571 | 84.9715628342 | 0.503966109187 | [cost section, row L2 Float64] |
| L4 | Float64 | 1915 | 165.600311888 | 0.982175001597 | [cost section, row L4 Float64] |
| L9 | Float64 | 1915 | 208.86256373 | 1.23876330018 | [cost section, row L9 Float64] |
| L6 | Float64 | 1227 | 354.532614384 | 2.10273197634 | [cost section, row L6 Float64] |
| L5 | Double64 | 186 | 4882.6924383 | 1.0 | [cost section, row L5 Double64] |
| L5 | BigFloat(256) | 92 | 377399.856865 | 1.0 | [cost section, row L5 BigFloat(256)] |



The item count of the L1 Float64 row, 1917, is from the table row. The findings file lists that row without a count. 
Run totals: the total wall time was 11957.993616580963 seconds against a cap of 12.0 hours. The plan had 1363 chunks and 15569 items. The run used 23 threads of 23 available. All 13 arms are complete. 

### 8.2 Sensitivity

The sensitivity table has 31 rows. The column is the median paired difference of lambda at n_H2 against the primary item, with 171 orbits per row except the tau rows. 
These rows are evaluated at n_H2 = 1e4. H1c is evaluated at N_1 = 1e7. The rows below are therefore not an evaluation of the H1c result. 

Direction. L1 and L5 values are of order 1e-9 at offset 0.3 (L5 1.73141166638e-9). They are of order 1e-5 or 1e-6 at the other offsets: L5 is -3.32612943661e-5 at 0.3 + 2pi/5, -3.23438506002e-5 at 0.3 + 3pi/5 and -3.48457834277e-6 at 0.3 + pi/5. L9 is 0.0 in every direction row. 

The eps rows for L2. With eps = 1e-10 the value is 2.26005874629e-9. With eps = 1e-12 it is 5.77729472453e-5. With eps = 1e-6 it is 4.04116577712e-10. 

Wrapped versus lifted is 0.0 for all eight methods of V. 
The column median_log10_nconv_ratio is 0.0 in all 31 rows. This is reported as printed. It has not been checked whether a constant column is an artefact. 

Tau sweep. Tau = 2 and 4 show no breakdown, with 171 orbits for L4 and L5. Tau = 8 has 154 orbits for L4 and L5, with breakdown count 17. Tau = 16 has 132 orbits for L4 with breakdown count 39 and 134 orbits for L5 with breakdown count 37. 
At tau = 16 the breakdown steps are 16 or 32 for L4 and 16 to 288 for L5. 

Symplectic residuals eps_pair(n) and eps_det(n) are shown in Figure F5 only. No numbers for them are quoted. 
Figure F1 shows the lambda_m(n) curves. 

## 9. Exploratory results

All items in this section are EXPLORATORY. They are not confirmatory, have no frozen rule and carry no claimed cause. 
Items D1 to D5 appear in Sections 3, 6 and 7, where they bear on a hypothesis. D6 and D7 are reported here. 

D6 (EXPLORATORY): short, empty and not decidable cells. The panel manifest marks two cells as short ("fewer qualifying grid cells than n_c"): K index 2 class S (25 drawn, 25 qualifying) and K index 5 class S (5 drawn, 5 qualifying). 
It marks as empty K index 5 class R (no elliptic periodic point) and K index 6 classes R (no elliptic periodic point), T and S (no qualifying grid cell). 
It marks K index 1 classes C, S and X as undefined (class undefined at K = 0). 
Class X shows 0 qualifying grid cells at every K index from 2 to 6. The preregistered selection rule does not draw class X from the screening and places the starts near the hyperbolic fixed point. 
In the H2 and H3 tables only class S at K index 5 is not decidable, with 1 orbit. 

D7 (EXPLORATORY): cell patterns. The largest cell median of R(n_H2) is 0.000775501992792 (class S, K index 3). The smallest is 5.82402750446e-12 (class X, K index 2). Both are in Table 6. 
The cell medians of n_conv are between 4.223954e6 and 7.079458e6 at K index 2 and 3. The horizon there is N_1 = 1e7. These values are in Table 8. 
At K index 4 they are 794328.0 (class S) and 63096.0 (class X). 
At K index 5 and 6 they lie between 158.0 and 266.5, excluding class S at K index 5. 
Several cells share the same median (for example 7.079458e6 in three cells of Table 8), which is consistent with n_conv taking recorded values of n (Section 2.2). The effect of this resolution on the bootstrap intervals is not examined. 
The tables do not show why. 

## 10. Limitations

1. The study rests on one draw of the panel (a single seed) with 883 orbits. The bootstrap intervals resample orbits within that draw, treat the orbits as independent resampling units, and do not reflect a second draw. The cell tables of H2 and H3 are descriptive. The composition of the 883 orbits is given in Section 2.3. 
2. The panel is run in Float64 only. Precision comparisons come from the ladder subset S_L of 92 orbits, and H5a rests on that subset. 
3. A set or cell with fewer than 20 orbits is not decidable and is neither pass nor fail. In the H2 and H3 tables only class S at K index 5 (5 drawn, 1 in W) is below that minimum. Class S at K index 2 has 25 drawn and 21 in W, which is above the minimum. Cell-level medians and intervals are descriptive, and the decisions use the pooled sets. H4 rests on one start per K, as the frozen text specifies, with no error bar. 
4. L2 runs in wrapped coordinates. This affects the H1c metric, the sub-rule that failed. The tables do not show whether this accounts for any part of the H1c result. 
5. The numbers were checked against the tables file of the run. The analysis program was verified against an independent recomputation (`analysis-verification.md`). The numbers were not recomputed from the stored series for the results document. Figures F1 to F7 were redrawn and viewed at paper stage 04; they plot the figure data of the tables file, unchanged. 
6. W is defined here by the preregistration, as in the results document and the findings file. The analysis program still prints a shorter phrase for W in its list of readings (class C at 0 < K <= K_c), which does not describe the cells of Tables 6 and 8; the membership of W (355 orbits) does not depend on that phrase. 
7. The controls table lists 10 of the 13 controls; the other 3 are the elliptic controls with exponent 0, which no decision rule of H1 uses (Section 3.3). Orbit 881 is an info control and the orbit dropped as drift. 
8. Table 1 gives two git commits, one for the run and one for the analysis. They are related through the repository history only; the study manifests record each commit separately and no single run record names both. 
9. The repository folder `data/` holds run outputs and is git-ignored (project folder map), so it is not part of the repository at the URL of Section 11. The stored series are packed for deposit as a dataset on Zenodo (1363 chunk files with the manifest, the panel, the preregistration, the parameter file, the file-format description and a SHA-256 list); its DOI is 10.5281/zenodo.23534755 (Section 11). 
10. The frozen preregistration is in the repository (Section 11). The results document records no external time-stamped deposit of it. 
11. Cost values come from a single run in which 23 of 23 available threads were used, with no repeated timing reported. The item sets differ between rows, and the Double64 and BigFloat rows have 186 and 92 items. The cost ratios are therefore not controlled comparisons. 
12. n_conv takes recorded values of n, because the series are stored at 20 points per decade, log-spaced, plus the first 100 steps. Its resolution and its effect on the bootstrap intervals are not examined. 
13. The column median_log10_nconv_ratio is 0.0 in all 31 sensitivity rows. It has not been checked whether this is an artefact. 
14. Open items carried over from the study hand-off: the STRONG_K constant was kept by owner decision, and its role is not described in the results document. The package requires Julia 1.11 or later, so no dependency manifest for Julia 1.10 is kept. 

## 11. Data and code availability

The repository is https://github.com/georgehadji/Chaos (GPL-3.0-or-later). 
The frozen preregistration is `workspaces/study/stages/01-design/output/p1-le-benchmark-preregistration.md` in the repository. Its header SHA-256 is in Table 1. The git SHA at freeze recorded in its header is 1e760e229eabf9d8323ad4d8097cc1dbd017f845, the parent commit of the freeze commit. Amendments A1 to A6 are in the same file and in its amendment log. 
The tagged release v0.1.0 was cut at milestone M3 exit, at commit 238e96b, a descendant of both commits of Table 1. The run used package version 0.0.1 at those commits, before the tag. The tagged release v0.1.1, at commit efeac86, adds the user documentation; between the analysis commit and v0.1.1 the package code changes only in docstrings and in a file-descriptor conversion of the chunk writer that is a no-op on Julia 1.13.2, the Julia version of the run. The tagged release v0.1.2, at commit 664866e, changes only the figure drawing (F1 to F7 redrawn, F1 and F2 also drawn in two halves) and the author metadata, and is archived on Zenodo. 
The run manifest and the stored series identify the run. The SHA-256 values are in Table 1. The stored series are deposited as a Zenodo dataset (limitation 9, Section 10); dataset DOI: [10.5281/zenodo.23534755](https://doi.org/10.5281/zenodo.23534755). 

### Availability statement

Repository: <https://github.com/georgehadji/Chaos> (GPL-3.0-or-later). The numbers in this paper come from one analysis run at commit d31b655d0c33e2f2f1a67bbbdeba9c8c8a3fd036. The stored series come from the run at commit 3e405f3324f8d25162580d3380db26c877d7ca63, an ancestor of the analysis commit with the same series code (Section 2.4). The tagged release v0.1.0 is at commit 238e96b; the run used package version 0.0.1, before the tag. The release to cite for the code is v0.1.2 (commit 664866e), archived on Zenodo with DOI [10.5281/zenodo.23535223](https://doi.org/10.5281/zenodo.23535223) and described by `CITATION.cff`. Study manifest: `workspaces/study/stages/03-run/output/p1-le-benchmark-manifest.json` (SHA-256 ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d). Results document: `workspaces/study/stages/05-report/output/p1-le-benchmark-results.md`. The stored series are not part of the repository; they are deposited as a Zenodo dataset (DOI 10.5281/zenodo.23534755).

## AI-use disclosure

Drafting, adversarial review and revision of this manuscript, and the analysis code, were done with Claude models (Anthropic, Claude Sonnet 5.5) operated through Claude Code under the direction of the author. The author decided the scope, the framing and the preregistered rules, and accepted or redirected the output at each stage checkpoint. Numbers in the text were copied from a results document produced by the analysis program, which was checked against an independent recomputation. Citations are limited to entries verified in the project bibliography. The author is responsible for the content.

## 12. Conclusion

Scope: one panel from a single draw (a single seed), run in Float64, with five nonzero kick values and one start per kick value in H4 (Section 10). 
Under the preregistered rules, H1 (claim C1) is falsified, on sub-rule H1c alone, with H1a and H1b met in full. It was evaluated in Float64 at N_1 = 1e7, with L2 in wrapped coordinates (limitation 4). 
Claim C2 has two preregistered tests with different verdicts. H2 (claim C2) is not supported (NEGATIVE). It measures R(n_H2) = Delta_V(n_H2) / lambda_ref over the set W at n_H2 = 1e4 steps. H3 (claim C2) is supported. It compares the median n_conv over W with the median n_conv over the strong-chaos set S, where n_conv is the maximum over V of the smallest n after which a method stays within a 0.05 relative band of lambda_ref. This paper does not combine the two verdicts into one verdict on C2. 
H4 did not reproduce the published value at K = pi/2 and passed at K = 2 pi. 
H5 is supported under its three sub-rules (H5a, H5b, H5c). The ladder test forced 35 of the 92 ladder-subset labels to unresolved (exploratory item D3). 
H4 and H5 are not tied to C1 or C2 in the findings file. 
TOOL: Chaoslab v0.1.2 is a tagged and documented release of the code the run used; every method of the run is `verified` (Section 2.1, Section 11). 

## Figure list

The figure files are in the `figures/` folder next to this manuscript. They were drawn from the figure data of the study's tables file by the figure extension of the repository at commit caaf95a; the figure data are those of the study output, unchanged. In this manuscript an F-label always names a figure. The figures are listed here and not embedded. The captions below describe what each figure draws.

- Figure F1: `p1-le-benchmark-fig-F1.pdf`. lambda_1(n) of the nine methods for the first drawn orbit of each populated cell (design class and K index), with the band between the smallest and the largest of the eight methods of V (L2 excluded).
- Figure F2: `p1-le-benchmark-fig-F2.pdf`. Median of R(n) over the orbits of each cell and over the pooled set W, with the 95 percent bootstrap band.
- Figure F3: `p1-le-benchmark-fig-F3.pdf`. Empirical distribution functions of n_conv per cell; colour by design class, line style by K index.
- Figure F4: `p1-le-benchmark-fig-F4.pdf`. Share of orbits forced to unresolved by each label test, among the orbits the test applies to, and the share of all orbits unresolved before and after the tests.
- Figure F5: `p1-le-benchmark-fig-F5.pdf`. Largest running maximum over the orbits of each stored diagnostic (det, symp, orth, pair, sum) per method, with the median threshold over the orbits; each orbit is tested against its own threshold.
- Figure F6: `p1-le-benchmark-fig-F6.pdf`. Thread nanoseconds per orbit-step by method and precision.
- Figure F7: `p1-le-benchmark-fig-F7.pdf`. Median paired difference of lambda(n_H2) from the primary item for the levels of tau, eps, initial direction and coordinate mode, over the kept S_D orbits.
