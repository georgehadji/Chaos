# Results: study p2-indicators

## Pointer key

Every number carries a pointer in square brackets.

- `[T1]` to `[T17]`, `[F1]` to `[F6]`: sections of `workspaces/study/stages/04-analyze/output/p2-indicators-tables.md`, with the row or statistic after the comma. In T3 to T8 the columns source and K come first (`[T3, all, U]`); in T10 and T13 the hypothesis comes first (`[T13, H3, R E1]`). In tables with one row per indicator X, `[T5, all, row X]` names the statistics of X in that table (T5: M_W2 X, M_S2 X, R X; T6: phi X, design phi X; T7: F X, decided X; T8: delta_dir X, dropped dir X, delta_prec X, dropped prec X).
- `[Decisions]`, `[Header]`: the blocks of that name in the tables file. `[params]`: the parameter list of the Header, key after the colon.
- `[findings X]`: section X of `workspaces/study/stages/04-analyze/output/p2-indicators-findings.md`.
- `[manifest: field]`: a field of `workspaces/study/stages/03-run/output/p2-indicators-manifest.json`.
- `[prereg: part]`: the frozen preregistration `workspaces/study/stages/01-design/output/p2-indicators-preregistration.md`.
- `[git: ...]`: a read-only git command on the repository.
- In tables, the last column is the pointer for the whole row.

## Header

Study slug: p2-indicators [manifest: study].

Run identity. All numbers come from one fully automatic analysis run over the stored series [findings A].

| Item | Value | Pointer |
|---|---|---|
| Preregistration header SHA-256 | a69bddb2a9045f8aac25ba7771fe7e5d90e70f9ff63203646e3684c8e19efff4 | [Header], [manifest: prereg.header] |
| Manifest SHA-256 | cdef3cf423e16a0b54e81e96fd460f4cb32d7b788325ffbf4e0581f1db27acce | [Header], [findings A] |
| Panel file SHA-256 | 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5 | [Header], [manifest: panel.sha256_end] |
| Parameter table SHA-256 (58 keys) | f5278e0f5c499f06674cc23cbe93d7d0a23720f4e4f77ece6ab48c4659063099 | [Header], [findings A] |
| Plan hash | 3823fa4a89fee047ab1f9ca63fee3d4e2dafa372c50b1d63c6de244fa8a18a81 | [Header], [manifest: plan_sha256] |
| Git commit of the analysis (tables header) | 10e438390917ca38b37c2fed329cd36d69070a49 | [Header], [findings A] |
| Git sha recorded by the run manifest | fd6e561fa4d014fbca9fd222fc5436c8a9a8b646 | [manifest: git.sha] |
| Module version | 0.1.2 | [Header], [manifest: versions.module] |
| Bootstrap seed | 20261012 | [Header], [params: bootstrap_seed] |
| Chunks, items | 695, 26964 | [Header], [manifest: plan_counts.chunks, plan_counts.items] |
| Run status | complete, all 695 chunks ok | [manifest: status, chunks[].status] |
| P1 reference labels: preregistration header, manifest and parameter table SHA-256 | 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf; ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d; a544130bc56dc66d8b1f7ec5d1e54197306f14dc19c7cd22fdbd94e4534b7282 | [Header], [findings A] |

The findings file gives the analysis commit and the manifest records the other commit; the findings file does not relate them. In the repository history [git: merge-base --is-ancestor fd6e561 10e4383; diff --name-only fd6e561 10e4383 -- src ext bin Project.toml] fd6e561 is an ancestor of 10e4383, and between them the only changed files under `src`, `ext`, `bin` and `Project.toml` are `src/Analysis/Analysis.jl`, `src/Analysis/CONTEXT.md`, `src/Analysis/cli.jl`, `src/Analysis/pipeline.jl` and new `src/Analysis/p2_*.jl` files. The series come from the run at fd6e561 (one session, 1028.427606344223 seconds [manifest: sessions[0].wall_seconds]); the tables come from the analysis at 10e4383. The commit 10e4383 is in no existing tag [git: tag --contains 10e4383, empty]; the tag for the TOOL claim is the owner's decision and is not set here [findings A].

What the study tested. Study p2-indicators is the P2 part of the joint P1+P2 paper and tests claim C3 (diagnostic inequivalence) with the two P2 outputs, detection latency and classification stability [prereg: Header; findings A]. Value types: FALSIFIES (H1 to H6), NEGATIVE allowed, TOOL; NOVEL and REPRODUCES are not claimed [prereg: Header]. C1 and C2 are not tested; the accepted P1 results give the reference labels. The rules were frozen on 2026-10-29 with Amendments A1 and A2, both before any data file, and are applied here without change [findings A].

Indicators. I_7 = {LE, SALI, FLI, MEGNO, E1, E2, WBA}; I_8 adds RLI [findings A]. LE is the Benettin estimator L1 [@Benettin1980a; @Benettin1980b], the study's reference estimator [prereg: Methods]. SALI [@Skokos2001], FLI [@Froeschle1997], MEGNO [@CincottaSimo2000], RLI [@SandorMaffione2016], E1 is the frequency path of the frequency map analysis [@Laskar2003], WBA is the weighted Birkhoff average [@SanderMeiss2020; @DuignanMeiss2023], E2 the weighted phase advance path, STR the stretching-number distance (threshold-free) [prereg: Scope check]. The indicator pairs FLI+RLI and MEGNO+SALI and the stretching numbers follow [@Maffione2011a] [prereg: Reported without a decision rule, Methods].

Reading key. E_all = panel plus grid orbits, minus controls, minus drift orbit 881, minus orbits with reference class U; E8 = orbits of E_all that carry RLI. rho is the reference class (C chaotic, R regular, U unresolved). W2 is the weak-chaos set (lambda_ref <= 0.1, or sticky), S2 the strong-chaos set (lambda_ref >= 1.0), T_s the sticky ensemble, T_s^design the a-priori design-class-S ensemble [findings A]. b_P = 12800 [params: b_primary]; N_c = 102400 [params: N_c]. "Panel" is the frozen P1 panel of 883 orbits; "grid" is the new grid of 2304 orbits (8 K values, 24 x 12 lattice) [params: grid_orbit_count, grid_i_count, grid_j_count]. Counts printing as `2.0` are formatting only [findings A].

Sets as printed. Primary (I_7 on E_all): W2 214 [T3, all, orbits], S2 448 [T4, all, orbits], T_s 30 [T6, all, orbits], T_s^design 96 [T6, all, design orbits]. Extension (I_8 on E8): W2 90, S2 303 [T13, H3, orbits W2, orbits S2], T_s 15, T_s^design 2 [T13, H4, orbits, design orbits]. Of W2, 125 are panel and 89 grid [T3, panel and grid, orbits]; of S2, 161 panel and 287 grid [T4, panel and grid, orbits]. Every grid orbit carries the flag ref_weak, 288 in each of the 8 grid rows [T1, column ref_weak].

## Confirmatory results

The primary decisions use I_7, the class rules of table R1 and the sets of E_all. The extension repeats each rule with I_8 on E8; a decision is ROBUST if primary and extension give the same verdict. Table R2 (fixed-threshold rules) is a sensitivity, not decisive. A set with fewer than 20 orbits [params: min_set] is not decidable, neither pass nor fail [findings B]. Every threshold of the parameter list equals the frozen text; that comparison was made by eye when the findings were drafted and is not a program output [findings A].

| Hypothesis | Primary (R1, I_7, E_all) | Extension (I_8, E8) | Flag | Table R2 | Pointer |
|---|---|---|---|---|---|
| H1 | SUPPORTED | supported | ROBUST | supported | [Decisions; T13, H1, robust; T10, H1, verdict] |
| H2 | NOT SUPPORTED | not supported | ROBUST | not supported | [Decisions; T13, H2, robust; T10, H2, verdict] |
| H3 | NOT SUPPORTED | not supported | ROBUST | not supported | [Decisions; T13, H3, robust; T10, H3, verdict] |
| H4 | SUPPORTED | not decidable (T_s 15 orbits, design ensemble 2, both below min_set 20) | NOT COMPARABLE | supported | [Decisions; T13, H4, robust; T10, H4, verdict] |
| H5 | SUPPORTED | supported | ROBUST | not supported | [Decisions; T13, H5, robust; T10, H5, verdict] |
| H6 | NOT SUPPORTED | supported | NOT ROBUST | supported | [Decisions; T13, H6, robust; T10, H6, verdict] |

Primary: H1, H4 and H5 supported; H2, H3 and H6 not supported. A not supported result is a result under NEGATIVE and has the same standing as a supported one [findings B].

### H1 (claim C3, classification divergence at a fixed budget): SUPPORTED

Rule. On W2 at b_P: U is the fraction of orbits whose seven I_7 verdicts are not all equal; A_min is the smallest agreement over the 21 pairs. SUPPORTED if U >= 0.20 [params: h1_u_min], the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) >= 0.10 [params: h1_u_lower_min, bootstrap_level, bootstrap_resamples] and A_min <= 0.80 [params: h1_amin_max]; otherwise not supported. U on S2 is a negative control without a rule [findings B, H1].

Measured. W2 orbits 214, decided for every indicator 214, dropped 0 [T3, all, orbits; orbits decided for every indicator; orbits dropped]. Verdict as printed: supported [T3, all, verdict].

| Quantity | Limit | Value | Result | Pointer |
|---|---|---|---|---|
| U | at least 0.20 | 0.303738317757 | pass | [T3, all, U] |
| Lower end of U interval | at least 0.10 | 0.242990654206 | pass | [T3, all, U, lower] |
| Upper end of U interval | none | 0.364485981308 | | [T3, all, U, upper] |
| A_min (pair MEGNO-E1) | at most 0.80 | 0.757009345794 | pass | [T3, all, A_min; A MEGNO-E1] |
| Mean kappa | none | 0.48042891781 | | [T3, all, mean kappa] |
| U on S2 (no rule) | none | 0.0 | | [T3, all, U on S2 (no rule)] |

Of the 21 pairwise agreements A on W2 the smallest are MEGNO-E1 0.757009345794, MEGNO-E2 0.780373831776, LE-E1 0.81308411215 and SALI-E1 0.817757009346, and the largest are LE-SALI 0.995327102804, SALI-FLI 0.971962616822 and LE-FLI 0.967289719626; all 21 are printed [T3, all, A X-Y rows].

By source: panel 125 orbits, U 0.328, A_min 0.752 [T3, panel]; grid 89 orbits, U 0.269662921348, A_min 0.76404494382 [T3, grid]. By K, cells with at least 20 orbits (descriptive; the rule decides on the pooled W2): K = 0.5, 99 orbits, U 0.292929292929, A_min 0.777777777778; K = 0.75, 23 orbits, U 0.304347826087, A_min 0.739130434783; K = 0.971635406, 83 orbits, U 0.253012048193, A_min 0.78313253012 [T3, all, K rows]. The other K cells have 0, 1, 2 or 5 orbits and are not decidable; T3 prints values for some (U 1.0 on 1 orbit at K = 1.25) [T3, all, 1.25].

Extension [T13, H1]: orbits 90, dropped 0; U 0.277777777778, interval [0.188888888889, 0.377777777778]; A_min 0.766666666667; mean kappa 0.474994986519; U on S2 0.0; verdict supported; ROBUST. Table R2 [T10, H1]: U 0.228971962617, interval [0.172897196262, 0.285046728972]; A_min 0.799065420561 against the limit 0.80; mean kappa 0.376311200177; U on S2 0.0; verdict supported.

Reading under FALSIFIES: the preregistered rule is met; H1 is supported and ROBUST [findings B, H1].

### H2 (claim C3, detection-time divergence): NOT SUPPORTED (NEGATIVE)

Rule. M_X is the median over S2 of the censored detection time t_det(X) (the first budget from which the verdict is stable up to N_c and equals rho; a wrong or missing verdict gets 2 N_c = 204800 [params: censor_factor]), for each X in I_7; Q = max_X M_X / min_X M_X. SUPPORTED if Q >= 10 [params: h2_q_min] and the lower end of its 95 percent bootstrap interval > 3 [params: h2_q_lower_gt]; otherwise not supported [findings B, H2].

Measured. S2 orbits 448 [T4, all, orbits]. M_X = 200.0 for each of LE, SALI, FLI, MEGNO, E1, E2 and WBA [T4, all, M LE to M WBA]. Q 1.0, interval [1.0, 1.0] [T4, all, Q]. Verdict as printed: not supported [T4, all, verdict]. Both conditions fail [findings B, H2].

Outcome shares on S2: correct 1.0, wrong 0.0, undetected 0.0 for LE, SALI, FLI, MEGNO, E1 and WBA; E2 correct 0.997767857143, wrong 0.0, undetected 0.00223214285714 [T4, all, share correct, share wrong, share undetected rows], which is 447 correct and 1 undetected orbit [F2, S2 E2 outcome rows]. The tables give no reason for the undetected orbit [findings C, item 6]. By source: panel 161 orbits, Q 1.0; grid 287 orbits, Q 1.0 [T4, panel and grid]. By K: K = 2 pi (368 orbits) and K = 10 (80 orbits), Q 1.0 each; the other K have 0 S2 orbits [T4, all, K rows].

Extension [T13, H2]: S2 orbits 303; M_X = 200.0 for all eight indicators including RLI; Q 1.0, interval [1.0, 1.0]; not supported; ROBUST. Table R2 [T10, H2]: S2 orbits 448; M_X = 200.0 for all seven; Q 1.0, interval [1.0, 1.0]; not supported.

This null has the same standing as a pass [findings B]. The budgets are a doubling lattice from 200 to 102400 [params: budgets]; the H2 medians sit at the smallest budget, 200, and an earlier detection cannot be observed, so the result cannot say whether the indicators differ below 200 [findings E, findings F].

### H3 (detection latency, regime dependence): NOT SUPPORTED (NEGATIVE)

Rule. R_X = M_X(W2) / M_X(S2) for each X in I_7. SUPPORTED if for at least 5 of the 7 indicators [params: h3_count_min] R_X >= 10 [params: h3_r_min] and the lower end of the 95 percent bootstrap interval of R_X > 3 [params: h3_r_lower_gt]; otherwise not supported [findings B, H3].

Measured. W2 orbits 214, S2 orbits 448 [T5, all, orbits W2, orbits S2].

| X | M_W2 | M_S2 | R_X | lower | upper | Meets the test | Pointer |
|---|---|---|---|---|---|---|---|
| LE | 1600.0 | 200.0 | 8.0 | 8.0 | 16.0 | no | [T5, all, row LE] |
| SALI | 1600.0 | 200.0 | 8.0 | 4.0 | 16.0 | no | [T5, all, row SALI] |
| FLI | 3200.0 | 200.0 | 16.0 | 8.0 | 16.0 | yes | [T5, all, row FLI] |
| MEGNO | 3200.0 | 200.0 | 16.0 | 8.0 | 16.0 | yes | [T5, all, row MEGNO] |
| E1 | 2400.0 | 200.0 | 12.0 | 4.0 | 32.0 | yes | [T5, all, row E1] |
| E2 | 400.0 | 200.0 | 2.0 | 1.0 | 4.0 | no | [T5, all, row E2] |
| WBA | 200.0 | 200.0 | 1.0 | 1.0 | 1.0 | no | [T5, all, row WBA] |

Indicators passing: 3 (the program's count); 5 are needed [T5, all, indicators passing]. The "Meets the test" column is the findings file's reading of the frozen test against the printed values, not a program output [findings B, H3]. Verdict as printed: not supported [T5, all, verdict].

By source (point values only): panel W2 125, S2 161, R FLI 16.0, MEGNO 16.0, E1 8.0; grid W2 89, S2 287, R FLI 8.0, MEGNO 16.0, E1 16.0; LE and SALI 8.0, E2 2.0 and WBA 1.0 in both [T5, panel and grid]. By K, T5 prints R = 1024.0 for all seven indicators at K = 2 pi, where W2 has 1 orbit and S2 368 [T5, all, 6.28318530718]; this cell is below min_set and the other K cells print no R [findings B, H3].

Extension [T13, H3]: W2 90, S2 303. R (lower, upper): LE 8.0 (4.0, 16.0); SALI 8.0 (4.0, 16.0); FLI 8.0 (4.0, 32.0); MEGNO 16.0 (8.0, 32.0); E1 16.0 (3.0, 576.0); E2 2.0 (1.0, 4.0); WBA 1.0 (1.0, 1.0); RLI 2.0 (2.0, 4.0). Indicators passing 1; not supported; ROBUST. The count threshold stays at 5 for eight indicators (Reading P21) [findings B, H3]. Table R2 [T10, H3]: W2 214, S2 448; R LE 2.0, SALI 1.0, FLI 2.0, MEGNO 1.0, E1 12.0 (4.0, 32.0), E2 2.0, WBA 1.0; indicators passing 1; not supported. The independent verification counts 0 passing indicators for R2 against the module's 1, a tolerated count with the same verdict [findings C, item 8].

Reading. LE and SALI (R = 8.0) are just under the threshold of 10 on a doubling lattice; E2 and WBA show hardly any regime dependence; under R2 only E1 reaches the level [findings F; T10, H3, R E1].

### H4 (claim C3, sticky orbits): SUPPORTED (primary); extension NOT COMPARABLE

Rule. T_s = orbits with rho = C and lambda_1(n_s) <= f_s lambda_ref, n_s = 1e5 [params: n_s], f_s = 0.2 [params: f_s]. For each X in I_7 minus LE, phi_X is the fraction of T_s orbits with c_X(N_c) = regular (a false regular verdict); Delta_phi = max_X phi_X minus min_X phi_X. LE is left out because T_s is defined from lambda_1. SUPPORTED if Delta_phi >= 0.50 [params: h4_dphi_min] and the lower end of its 95 percent bootstrap interval >= 0.25 [params: h4_dphi_lower_min]; otherwise not supported. T_s^design is reported with the same statistic [findings B, H4].

Measured. T_s orbits 30, all six indicators decided on 30; T_s^design 96 orbits, all decided on 96 [T6, all, orbits; decided rows; design orbits; design decided rows].

| X | phi_X (T_s) | phi_X (T_s^design) | Pointer |
|---|---|---|---|
| SALI | 0.2 | 0.0520833333333 | [T6, all, row SALI] |
| FLI | 0.266666666667 | 0.0520833333333 | [T6, all, row FLI] |
| MEGNO | 0.233333333333 | 0.0625 | [T6, all, row MEGNO] |
| E1 | 0.7 | 0.229166666667 | [T6, all, row E1] |
| E2 | 0.633333333333 | 0.0729166666667 | [T6, all, row E2] |
| WBA | 0.0666666666667 | 0.0104166666667 | [T6, all, row WBA] |

Delta_phi 0.633333333333, interval [0.466666666667, 0.8] [T6, all, Delta_phi]; both conditions pass. Verdict as printed: supported [T6, all, verdict]. In both ensembles the largest phi_X is E1 and the smallest is WBA [findings B, H4]. T_s^design, reported next to the decision and not a decision: Delta_phi 0.21875, interval [0.135416666667, 0.302083333333], verdict as printed not supported [T6, all, design Delta_phi, design verdict].

By source: panel 15 orbits, Delta_phi 0.866666666667; grid 15 orbits, Delta_phi 0.466666666667 [T6, panel and grid]. By K all cells are below min_set (T_s orbits 8, 2, 12, 1, 5, 2 at K = 0.5, 0.75, 0.971635406, 1.25, pi/2, 2) and descriptive [T6, all, K rows].

Extension [T13, H4]: T_s has 15 orbits, below min_set 20. phi: SALI 0.2, FLI 0.266666666667, MEGNO 0.266666666667, E1 0.466666666667, E2 0.533333333333, WBA 0.0666666666667, RLI 0.0; Delta_phi 0.533333333333, interval [0.333333333333, 0.8]; verdict as printed not decidable. The extension design ensemble has 2 orbits (D90: T_s^design is intersected with E8), phi 0.0 for all seven indicators, Delta_phi 0.0 [0.0, 0.0], not decidable. Flag NOT COMPARABLE. These numbers are context, not a verdict [findings B, H4]. Table R2 [T10, H4]: T_s 30 orbits; phi SALI 0.0, FLI 0.0, MEGNO 0.0, E1 0.7, E2 0.633333333333, WBA 0.0666666666667; Delta_phi 0.7, interval [0.566666666667, 0.866666666667]; supported. Design ensemble under R2: 96 orbits, Delta_phi 0.229166666667, interval [0.145833333333, 0.3125], not supported.

Reading. The spread is driven by E1 and E2 calling many sticky orbits regular at N_c, against WBA [T6, all, row E1, row E2, row WBA]. The conclusion depends on the definition of the sticky set: it holds on T_s (30 orbits, from lambda_1 at n_s = 1e5) and under R2, and is not supported on the design ensemble (96 orbits). The paper should say that it rests on one definition [findings F]. The preregistration states a limit of that definition: extreme sticky zones may stay trapped beyond the reference horizon, and such an orbit gets rho = R [prereg: Panel, R4; @HarsoulaKaramanosContopoulos2019].

### H5 (classification stability over the budget): SUPPORTED

Rule. On W2, F_X is the fraction of orbits (decided at b_P and at N_c) with c_X(b_P) different from c_X(N_c); F_med is the median of F_X over I_7. SUPPORTED if F_med >= 0.10 [params: h5_fmed_min] and the lower end of its 95 percent bootstrap interval >= 0.05 [params: h5_fmed_lower_min]; otherwise not supported [findings B, H5].

Measured. W2 orbits 214 [T7, all, orbits].

| X | F_X | Orbits decided | Pointer |
|---|---|---|---|
| LE | 0.177570093458 | 214 | [T7, all, row LE] |
| SALI | 0.182242990654 | 214 | [T7, all, row SALI] |
| FLI | 0.205607476636 | 214 | [T7, all, row FLI] |
| MEGNO | 0.239436619718 | 213 | [T7, all, row MEGNO] |
| E1 | 0.200934579439 | 214 | [T7, all, row E1] |
| E2 | 0.0841121495327 | 214 | [T7, all, row E2] |
| WBA | 0.0280373831776 | 214 | [T7, all, row WBA] |

F_med 0.182242990654, interval [0.130841121495, 0.228971962617] [T7, all, F_med]; both conditions pass. Verdict as printed: supported [T7, all, verdict]. By source: panel 125 orbits, F_med 0.176; grid 89 orbits, F_med 0.191011235955 [T7, panel and grid]. By K, cells with at least 20 orbits: K = 0.5 (99) 0.181818181818; K = 0.75 (23) 0.130434782609; K = 0.971635406 (83) 0.180722891566 [T7, all, K rows]. T7 also prints F_X(b, W2) for every budget b_0 to b_9; at b_9 = N_c it is 0.0 for all seven because the verdict is compared with itself [findings B, H5].

Extension [T13, H5]: 90 orbits; F LE 0.188888888889, SALI 0.188888888889, FLI 0.177777777778, MEGNO 0.233333333333, E1 0.233333333333, E2 0.0777777777778, WBA 0.0333333333333, RLI 0.0222222222222; F_med 0.183333333333, interval [0.111111111111, 0.25]; supported; ROBUST. Table R2 [T10, H5]: 214 orbits; F LE 0.0420560747664, SALI 0.00467289719626, FLI 0.0327102803738, MEGNO 0.00934579439252, E1 0.200934579439, E2 0.0841121495327, WBA 0.0280373831776; F_med 0.0327102803738, interval [0.0140186915888, 0.0607476635514]; not supported.

Reading. That the finite-budget verdict differs from the verdict at N_c on weak chaos depends on the rule family, growth ratio (R1) or fixed threshold (R2). The primary decision stays R1 as preregistered; the R2 result is a sensitivity shown next to it, not instead of it [findings F].

### H6 (classification stability under direction and precision): NOT SUPPORTED (NEGATIVE)

Rule. Both parts at N_c. (a) On the designated subset S_D (orbits of E_all in it), for X in {LE, SALI, FLI, MEGNO}, delta_dir(X) is the fraction of orbits whose verdict is not identical in the five deviation directions. (b) On the ladder subset S_L, for X in I_7, delta_prec(X) is the fraction of orbits whose verdict differs between Float64 and Double64. SUPPORTED if max_X delta_dir(X) or max_X delta_prec(X) >= 0.05 [params: h6_flip_min]; otherwise not supported. No bootstrap interval [findings B, H6].

Measured. Orbits S_D* 109, S_L* 51 [T8, all, orbits S_D*, orbits S_L*].

| X | delta_dir | Dropped (dir) | delta_prec | Dropped (prec) | Pointer |
|---|---|---|---|---|---|
| LE | 0.0 | 7 | 0.0 | 3 | [T8, all, row LE] |
| SALI | 0.045871559633 | 0 | 0.0 | 0 | [T8, all, row SALI] |
| FLI | 0.0 | 0 | 0.0 | 0 | [T8, all, row FLI] |
| MEGNO | 0.0 | 9 | 0.0 | 4 | [T8, all, row MEGNO] |
| E1 | not in the rule | not in the rule | 0.0196078431373 | 0 | [T8, all, row E1] |
| E2 | not in the rule | not in the rule | 0.0 | 0 | [T8, all, row E2] |
| WBA | not in the rule | not in the rule | 0.0 | 0 | [T8, all, row WBA] |

Largest delta_dir 0.045871559633 (SALI); largest delta_prec 0.0196078431373 (E1); both are below 0.05. Verdict as printed: not supported [T8, all, verdict]. An orbit undecided in some direction or precision is dropped and counted (Reading P23) [findings B, H6].

Extension [T13, H6]: S_D* 51, S_L* 51; delta_dir LE 0.0, SALI 0.078431372549, FLI 0.0, MEGNO 0.0 (dropped dir LE 6, MEGNO 7); delta_prec E1 0.0196078431373, all others 0.0 (dropped prec LE 3, MEGNO 4); verdict supported (SALI 0.078431372549 at least 0.05); NOT ROBUST (primary not supported, extension supported). Table R2 [T10, H6]: S_D* 109, S_L* 51; delta_dir SALI 0.0642201834862, LE, FLI and MEGNO 0.0; delta_prec E1 0.0196078431373, others 0.0; no orbit dropped; supported.

Reading. The primary result is a near miss on one indicator (SALI). The paper should not state input-stability of the verdicts as established [findings F].

### Rule-validity gate on the integrable control (T12)

Rule. For each X in I_8 the share of K = 0 gate orbits (80 panel orbits of K index 1, classes T and R, regular by construction [params: gate_k_index, gate_classes, gate_size]) with c_X(N_c) = chaotic must be at most 0.10 [params: gate_max_share]. A violation is reported as "rule invalid on the integrable control" with the indicator named, and a second version of each decision table without that indicator is added, never in place of the first [findings B, T12].

Per indicator (decided, chaotic, undecided, not applicable; share chaotic 0.0 and violated false for every one) [T12, rows by indicator]: LE 80, 0, 0, 0; SALI 80, 0, 0, 0; FLI 80, 0, 0, 0; MEGNO 80, 0, 0, 0; E1 75, 0, 5, 0; E2 80, 0, 0, 0; WBA 75, 0, 5, 0; RLI 8, 0, 0, 72.

No indicator violates the gate [Decisions]; no second version of a decision table is needed or made. The gate set is not intersected with E_all, and orbits undecided at the last budget leave the denominator (Reading P26) [findings B, T12].

## Claim C3

- The preregistration header places both P2 outputs under C3; the headings of H3, H5 and H6 name an output, not the label C3 [findings B].
- H1 (classification divergence at a fixed budget) is supported and ROBUST. H4 (sticky orbits) is supported on its primary set only and NOT COMPARABLE in the extension [Decisions].
- H2 (detection-time divergence) and H3 (its regime dependence) are not supported, both ROBUST [Decisions].
- The stability results depend on the rule family or the set: H5 is supported under R1 and not under R2; H6 is not supported under R1 and supported under R2 and in the extension [Decisions].
- The frozen text gives no rule that combines H1 to H6 into one verdict on C3, and none is applied here. The evidence is mixed, and how to word the claim is the owner's call [findings B, findings F].

## Reported without a decision rule

### Figures F1 to F6

The figure files exist: `p2-indicators-fig-F1.pdf` to `p2-indicators-fig-F6.pdf` in `workspaces/study/stages/04-analyze/output/`, drawn from the figure data blocks F1 to F6 of the tables file by one command (`bin/chaoslab.jl analyze` under `julia --project=verify`, figure drawers of Amendment D91, commit 248f1a1), run on 2026-10-30 with exit code 0. The tables of that run equal the tables file except the git commit line, and the tables file is kept as accepted. The owner visual check of the six figures is accepted [findings E]. No claim of this document rests on how a figure looks.

Data blocks: F1 agreement counts (cc, cr, rc, rr), N, A, kappa and Spearman for the pairs of I_8 on W2 and S2 at b_P and N_c [F1]; F2 censored detection times, detection times and outcomes per indicator on S2 [F2]; F3 fraction chaotic per indicator and budget k [F3]; F4 AUC rows per indicator, budget, K and source [F4]; F5 one row per grid orbit (2304 rows) with the verdicts of the eight indicators at b_P and N_c [F5]; F6 stretching-distance percentiles (p5 to p95) by class C, R and T_s per step [F6].

Two rows as examples. F1, W2 at 12800, pair LE and SALI: cc 185, cr 1, rc 0, rr 28, N 214, A 0.995327102804 [F1, row W2 12800 LE SALI], equal to [T3, all, A LE-SALI]. F2, E2 on S2: 445 orbits have the censored detection time 200 [F2, row S2 E2 t_censored 200].

### Sensitivities (T11)

Verdict agreement A and kappa on S_L (51 orbits per row) between a sensitivity spec and its primary:

| Pair | A at b_P | kappa at b_P | A at N_c | kappa at N_c | Pointer |
|---|---|---|---|---|---|
| cos_vs_wba | 1.0 | 1.0 | 1.0 | 1.0 | [T11, cos_vs_wba, k 6 and 9] |
| planar_vs_pa1 | 0.588235294118 | 0.0 | 0.588235294118 | 0.0 | [T11, planar_vs_pa1, k 6 and 9] |
| qp_vs_fq1 | 1.0 | 1.0 | 1.0 | 1.0 | [T11, qp_vs_fq1, k 6 and 9] |
| raw_vs_fq1 | 0.960784313725 | 0.917874396135 | 0.921568627451 | 0.833333333333 | [T11, raw_vs_fq1, k 6 and 9] |

Median |dD| at b_P: cos_vs_wba 6.75226967672e-8; planar_vs_pa1 5.13320317084; qp_vs_fq1 0.916825918018; raw_vs_fq1 6.2172489378999994e-15 [T11, rows k 6]. The planar phase advance has A 0.607843137255 at k = 3 and 0.588235294118 from k = 4 on [T11, planar_vs_pa1, k 3 to 9]; the tables do not show why [findings D2].

### Cost (T17) and run totals

Measured ns per orbit-step against the budget line [T17, rows by spec]: bk1 3000.59036817 against 15000.0; fq1 16850.3797613 against 84300.0 (transform seconds 3024.198199); fq1_qp 8749.99503139 against 26800.0; fq1_raw 9683.55278019 against 26800.0; pa1 1810.47806122 against 12600.0; pa_planar 248.296322696 against 2000.0; pr 223.405047816 against 26000.0; st 1702.27590881 against 2200.0; tg 2658.01590681 against 18000.0. Every measured value is below its line and st is the closest [findings D6]. Manifest: total wall time 1028.427606344223 seconds against a cap of 12.0 hours [manifest: wall.total_seconds, wall.cap_hours]; 23 threads used of 23 available on a host with 24 cores [manifest: threads.used, threads.available, host.cores]; device NVIDIA GeForce RTX 5060 Laptop GPU [manifest: device]; all 17 arms complete [manifest: arms.*.status]; no control stop and no override [manifest: control_stops, overrides].

### Controls, flags and packages (T15, T16)

Controls held to the bound (871, 873, 875, 877, 879, 883) have stored deviations from 1.11022302463e-16 to 5.26802579337e-9; controls marked info (878, 880, 881, 882) have 0.222018884744, 0.443891109126, 4.08358102888e-9 and 0.103203221967; controls 872, 874 and 876 have no scope (exponent 0.0) [T16, control rows]. The Y5 control (orbit 883) has no torus, so bk1 and fq1 omit it [T16, omission rows]. The K = 0 line (largest circle distance of nu1 to p0/(2 pi)) falls from 8.77500818086e-6 at W = 100 to 5.551115123130001e-17 at W = 12800, 25600 and 51200, with 5 fq1 orbits left out as non-finite [T16, k0line rows]. 

## Exploratory results

Everything in this section is EXPLORATORY. None is a decision of the confirmatory section, none has a frozen rule, and none is to be read as confirmatory. No computation, test or p-value was made for it; the values are read from the tables, and the tables give no causes, and none is claimed [findings D]. The tables file lists no exploratory items of its own [Exploratory items: none], so the pointers name the tables the observations are read from.

**D1 (EXPLORATORY): AUC of the chaos scores [T9, F4].** AUC per K, panel and grid pooled. At b_P (k = 6) the lowest values are STR at K = 0.5 (0.863760302411, n_C/n_R 99/326) and MEGNO at K = 0.971635406 (0.956349206349, 252/172) [T9, F4, rows STR and MEGNO, k 6, pooled]. At N_c (k = 9) they are STR at K = 0.5 (0.864209580467) and MEGNO at K = 0.971635406 (0.912623931813) [T9, F4, rows STR and MEGNO, k 9, pooled]. For LE, SALI, FLI, E1, E2 and WBA at b_P the AUC at the four K values 0.5, 0.971635406, 1.57079632679 and 6.28318530718 lies between 0.952913279133 (SALI, K = 1.57079632679) and 1.0 (FLI, K = 1.57079632679) [T9, F4, rows k 6, pooled]. RLI has 1.0 at K = 0.971635406 (147/129) and K = 1.57079632679 (220/63), and no AUC ("nothing") at K = 0.5 (6/269) and K = 6.28318530718 (295/3), where a class has fewer than 20 orbits [T9, F4, RLI k 6]. T9 prints all K and budgets [findings D1].

**D2 (EXPLORATORY): sensitivities. See the table of T11 above [findings D2].

**D3 (EXPLORATORY): regular range and sticky orbits [T14].** Share of sticky orbits whose score at N_c lies below the median of the regular class, where printed: 0.0 for LE, SALI, FLI, E1, E2 and WBA at K = 0.5, 0.75, 0.971635406, 1.25, 1.57079632679 and 2.0; MEGNO 0.125 at K = 0.5 and 0.5 at K = 2.0; STR 0.625 at K = 0.5 and 0.2 at K = 1.57079632679; 0.0 elsewhere. No share is printed at K = 0, 3, 2 pi and 10 [findings D3; T14].

**D4 (EXPLORATORY): packages cross-tabulated with rho, no combination rule [T15].** On T_s (rho = C), FLI with RLI at b_P: both chaotic 8, first regular and second chaotic 6, both regular 1 (15 orbits); at N_c: both chaotic 11, first regular and second chaotic 4. MEGNO with SALI on T_s at b_P: both chaotic 17, MEGNO chaotic and SALI regular 5, MEGNO regular and SALI chaotic 2, both regular 6 (30 orbits); at N_c: 22, 1, 2, 5 in the same order [T15, T_s rows]. On the regular set (rho = R) at b_P: FLI with RLI both regular 874; MEGNO with SALI both chaotic 4, MEGNO chaotic and SALI regular 27, MEGNO regular and SALI chaotic 3, both regular 1074, MEGNO undecided with SALI chaotic 12, undecided with SALI regular 39; at N_c: 7, 52, 1, 1066, 6, 27 in the same order [T15, regular rows]. The regular-set counts differ between the two packages because FLI+RLI is limited to the orbits that carry RLI [findings D4; @Maffione2011a].

**D5 (EXPLORATORY): flag counts [T16].** bk1 below_floor 1358, reldig_undefined 1346, zero 11; fq1 below_floor 1152, degenerate 14, mirror1 1839, mirror2 1885; pa1 below_floor 1331; pa_planar below_floor 21 [T16, flag rows]. The other flag rows (fq1 near_zero, refined and tie; fq1_qp; fq1_raw) are in T16.

**D6 (EXPLORATORY): cost. See the cost paragraph above [findings D6].

**D7 (EXPLORATORY): cell and budget patterns already in the tables.** (i) H1 per K cell above min_set: U 0.292929292929, 0.304347826087, 0.253012048193 at K = 0.5, 0.75, 0.971635406 [T3, all, K rows]. (ii) H3 by source, point values with R >= 10: FLI and MEGNO on the panel, MEGNO and E1 on the grid, against FLI, MEGNO and E1 pooled; no interval is printed by source [T5, panel and grid]. (iii) F_X(b, W2) falls with the budget for every indicator: LE 0.682692307692 at b_0, 0.177570093458 at b_6, 0.0794392523364 at b_8, 0.0 at b_9; WBA 0.018691588785 at b_0 to b_3 and 0.0 at b_9 [T7, all, F LE W2 b_0, b_6, b_8, b_9; F WBA W2 b_0 to b_3, b_9; findings D7].

## Deviations and readings

Overrides recorded in the manifest: none [manifest: overrides; Deviations]. The items below bear on a confirmatory number, with the decisions each could touch. This document does not say that any of them changed a verdict; the tables do not show it, and no run without them was made [findings C].

1. P13. lambda_ref and lambda_ns of the panel orbits come from the P2 wrapped L1 items; the preregistration says they equal the P1 L1 values and the tables state that they do not [Deviations]. The size of the difference is not in the tables; the findings point to the decision log (D89). It could touch every decision that uses W2, S2 or T_s (H1 to H5, primary, extension and R2); H6 does not use these sets [findings C, item 1].
2. P29. The panel reference labels (rho) are those of the P1 pipeline at the P1 horizon, made with lifted L5 items, while the L1 item of table R1 is wrapped [Deviations]. It could touch the membership of C and R in every set, the correct and wrong outcomes of H2 and H3, the statistic phi of H4, and through E_all H6 [findings C, item 2].
3. E23 b, P14, P19. The reference class of grid orbits uses the reduced L1/L5 rule and every grid orbit carries ref_weak [Deviations; T1, column ref_weak]. It could touch every set with grid orbits (W2 89 of 214, S2 287 of 448, T_s 15 of 30, all of E8 beyond the ladder subset), which is why the tables split H1 to H5 by source [findings C, item 3].
4. P15. The budget step n(b) is the first recorded step not before b, so indicators are read at slightly different steps (up to 1.13 b); the step is printed in T2. It could touch every budget-based verdict, H1 to H5, in R1 and R2 [Deviations; findings C, item 4].
5. P16. g_SALI = -L, because L (ln GALI_2) is negative and, read literally, every SALI verdict would be undecided [Deviations]. It could touch every statistic with SALI (all of H1 to H6; the H6 maximum delta_dir is SALI's) [findings C, item 5].
6. P17. WBA is undecided at a window if any of its four rows has the flag zero or a non-finite dig; E1 and E2 are undecided on their flag or a non-finite D; an orbit undecided for E1 or WBA at every budget is undetected, with the censored time, and stays in the regular set [Deviations]. It could touch H2 and H3 (E1 and WBA) and the E1 and WBA statistics in H1, H4 and H5 [findings C, item 6; T2, T12].
7. D90. In the E8 extension T_s^design is intersected with E8 (owner decision at the stage 04 checkpoint). It touches only the design statistic of H4 in the extension (2 orbits); the H4 flag is NOT COMPARABLE with or without it, because T_s has 15 orbits in the extension [findings C, item 7; T13, H4].
8. D89w end tolerance (verification, not the analysis). The independent recomputation accepts an interval end within max(0.05 of the width, one lattice step 1/n), and one count (R2 H3 indicators passing: module 1, independent 0) is tolerated. No table number changes. The printed lower ends that carry a decision are far from their thresholds: H1 0.242990654206 against 0.10 [T3, all, U, lower], H4 0.466666666667 against 0.25 [T6, all, Delta_phi], H5 0.130841121495 against 0.05 [T7, all, F_med]; this compares printed numbers and is not a recomputation [findings C, item 8].
9. Bootstrap and set readings (P18, P20, P21, P27). One set of draws per label; 5 percent dropped replicates make an interval undefined; intervals for pooled statistics only; an empty denominator makes the decision not decidable (this bears on the H4 extension and on the K and source cells). The H3 count threshold stays at 5 for eight indicators; W2 and S2 may overlap through T_s and the overlap is 0 in every row [findings C, items 9 to 11; T1, column W2 and S2]. The share 0.05 of dropped replicates and the design class S of T_s^design are module constants, not keys of the 58-key table; no threshold is changed by it [findings C, item 12].

The full list of readings is in the tables file, section Readings [Readings].

## Limitations

From the findings file [findings E]:

- One draw of the panel (883 orbits) and one grid (2304 starts on one lattice) [params: grid_orbit_count]. The bootstrap intervals resample orbits within these and do not reflect a second draw; panel and grid orbits are resampled as independent.
- Reduced reference rule and mixed labels. All grid orbits are flagged ref_weak, and the panel labels come from P1 at its horizon with different items [T1, column ref_weak; Deviations].
- RLI on a subset. E8 holds only the ladder-subset and grid orbits, so the extension has W2 90, S2 303, T_s 15 and design ensemble 2 [T13].
- Finite budgets on a doubling lattice from 200 to 102400, so detection-time medians and ratios take lattice values (200, 400, 1600, 3200; ratios 8 and 16). The H2 medians of 200.0 are at the floor of 200 [params: budgets; T4, all, M LE to M WBA].
- Cells with fewer than 20 orbits (K and source cells) are descriptive and not decidable. The K cells other than 0.5, 0.75 and 0.971635406 hold 0 to 5 W2 orbits [T3, T5]; T_s holds 8, 2 and 12 orbits at those three K and 1, 5 and 2 at K = 1.25, pi/2 and 2 [T6, all, K rows].
- T_s size. The primary T_s has 30 orbits and its conclusion rests on one definition (lambda_1 at n_s) [T6, all, orbits; findings F]. The extension of H4 has 15 orbits (design ensemble 2), so H4 is NOT COMPARABLE [T13, H4].
- Orbits with rho = U are outside every set; the dropped orbits are listed in T1 [T1, Orbits dropped from the sets].
- The independent verification of the analysis matched 14146 of 14146 compared quantities after the stage 04 rework (D89w). The pass after the readings was not blind: the first comparison showed 21 differences, and D90 and the interval-end tolerance were settled by the owner after the data were seen. The verification report says it did not re-run either bootstrap generator with other seeds [findings E].
- The numbers were copied from the tables file and were not recomputed from the stored series in this step. One number set of the findings (section C item 1, from the decision log D89) is not in the tables and is not quoted here [findings E].

## Uncited claims

- The definition of the weighted phase advance path E2 is not cited to a key; the findings file and the preregistration give no verified key for it. The definition of STR is cited only to [@Maffione2011a].
- The value K_c = 0.971635406 is used as given in the parameter list [params: grid_k_values]. The preregistration cites [@FrahmShepelyansky2010] for it as the golden-mean breakup [prereg: Grid G]; the findings file gives no key.
- Every threshold of tables R1 and R2 is proposed by this study; the preregistration states that none is taken from the literature [prereg: Analysis plan]. No citation supports the thresholds.
- The standard map as the system of the study is not cited here, because the findings file gives no key.
