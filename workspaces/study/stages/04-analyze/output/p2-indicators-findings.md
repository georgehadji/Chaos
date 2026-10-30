# Findings: study p2-indicators (draft for the stage 04 checkpoint)

## A. What this document is

Study p2-indicators is the P2 part of the joint P1+P2 paper. It tests claim C3 (diagnostic inequivalence) with the two P2 outputs of the paper row, detection latency and classification stability. The preregistered value types are FALSIFIES (H1 to H6 each have a metric and a numeric rule), NEGATIVE allowed (a hypothesis that is not supported is reported with the same figures as a supported one) and TOOL (the run uses only modules whose catalog status is verified; the TOOL claim cites a tagged release whose number is the owner's decision and is not in the tables). NOVEL and REPRODUCES are NOT claimed. C1 and C2 are not tested here; the accepted P1 results are reused as reference labels. The decision rules were frozen on 2026-10-29 (Amendments A1 and A2, both made before any data file), and this document applies them without change.

Sections B to D come from one fully automatic analysis run (`analyze`, no human or AI step in it) over the stored series, and from nothing else. Numbers are copied from `p2-indicators-tables.md` (the tables file) as printed; table cosmetics (counts printing as 2.0) are formatting only. The run is identified by these values.

| Item | Value |
|---|---|
| Preregistration header SHA-256 | a69bddb2a9045f8aac25ba7771fe7e5d90e70f9ff63203646e3684c8e19efff4 |
| Manifest SHA-256 | cdef3cf423e16a0b54e81e96fd460f4cb32d7b788325ffbf4e0581f1db27acce |
| Panel file SHA-256 | 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5 |
| Parameter table SHA-256 (58 keys) | f5278e0f5c499f06674cc23cbe93d7d0a23720f4e4f77ece6ab48c4659063099 |
| Plan hash | 3823fa4a89fee047ab1f9ca63fee3d4e2dafa372c50b1d63c6de244fa8a18a81 |
| Git commit | 10e438390917ca38b37c2fed329cd36d69070a49 |
| Module version | 0.1.2 |
| Bootstrap seed | 20261012 |
| Chunks | 695 |
| Items | 26964 |
| P1 preregistration header SHA-256 (reference labels) | 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf |
| P1 manifest SHA-256 | ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d |
| P1 parameter table SHA-256 | a544130bc56dc66d8b1f7ec5d1e54197306f14dc19c7cd22fdbd94e4534b7282 |

The tables Header also lists `panel_build_sha256` = 650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46. It carries no digest of the figure data blocks.

Reading key. Indicators: I_7 = {LE, SALI, FLI, MEGNO, E1, E2, WBA}; I_8 = I_7 plus RLI. E_all = panel plus grid orbits, minus controls, minus the drift orbit 881, minus orbits with reference class U. E8 = the orbits of E_all that carry RLI. rho is the reference class (C chaotic, R regular, U unresolved). W2 is the weak-chaos set (lambda_ref <= 0.1, or sticky), S2 the strong-chaos set (lambda_ref >= 1.0), T_s the sticky ensemble, T_s^design (printed T_s_design) the a-priori design-class-S ensemble. b_P = 12800 is the primary budget (budget index k = 6), N_c = 102400 the common horizon (k = 9). "Panel" is the frozen P1 panel of 883 orbits; "grid" is the new K x initial-condition grid of 2304 orbits (8 K values, 24 x 12 lattice). The labels S_D* and S_L* are as printed in T8.

Set sizes as printed. Primary (I_7 on E_all): W2 214, S2 448, T_s 30, T_s^design 96 (T3, T5, T6). Extension (I_8 on E8): W2 90, S2 303, T_s 15, T_s^design 2 (T13). Of the W2 orbits 125 are panel and 89 grid; of the S2 orbits 161 are panel and 287 grid (T3, T4). Every grid orbit carries the flag ref_weak (288 in each of the 8 grid rows of T1; see C, item 3).

Every threshold of the decision rules in the run's parameter list equals the number in the frozen text: h1_u_min 0.2, h1_u_lower_min 0.1, h1_amin_max 0.8, h2_q_min 10.0, h2_q_lower_gt 3.0, h3_r_min 10.0, h3_r_lower_gt 3.0, h3_count_min 5, h4_dphi_min 0.5, h4_dphi_lower_min 0.25, h5_fmed_min 0.1, h5_fmed_lower_min 0.05, h6_flip_min 0.05, min_set 20, b_primary 12800, N_c 102400, and the table R1 and R2 numbers (ratio_threshold 5.0, le_floor 0.0001, rli_threshold 1.0e-7, e1_d_threshold and e2_d_threshold -4.0, wba_digits_threshold 6.0, r2_le_factor 3.0, r2_sali_threshold 1.0e-8, r2_fli_factor 3.0, r2_megno_threshold 2.5). This was compared by eye with the frozen text when this document was drafted; it is not a program output.

## B. Confirmatory decisions

Order of decision: H1, H2, H3, H4, H5, H6. The primary decisions use I_7, the class rules of table R1, and the sets built from E_all. The extension repeats each rule with I_8 on E8; a decision is ROBUST if primary and extension give the same verdict. Table R2 (fixed-threshold rules) is a sensitivity, reported, not decisive. A set with fewer than 20 orbits (min_set) is "not decidable": neither pass nor fail.

Summary as printed by the program (Decisions block):

| Hypothesis | Primary (R1, I_7, E_all) | Extension (I_8, E8) | Flag | Table R2 sensitivity |
|---|---|---|---|---|
| H1 | supported | supported | ROBUST | supported |
| H2 | not supported | not supported | ROBUST | not supported |
| H3 | not supported | not supported | ROBUST | not supported |
| H4 | supported | not decidable | NOT COMPARABLE | supported |
| H5 | supported | supported | ROBUST | not supported |
| H6 | not supported | supported | NOT ROBUST | supported |

Primary: H1, H4 and H5 supported; H2, H3 and H6 not supported. A not supported result is a result under the value type NEGATIVE and has the same standing as a supported one.

### H1 (claim C3, classification divergence at a fixed budget): SUPPORTED

Frozen rule. On W2 at b_P = 12800: U(b_P) is the fraction of orbits whose seven I_7 verdicts are not all equal; A_min(b_P) is the smallest agreement over the 21 pairs of I_7. H1 is SUPPORTED if U(b_P) >= 0.20, the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits) is >= 0.10, and A_min(b_P) <= 0.80. Otherwise not supported (NEGATIVE). U on S2 at b_P is reported as a negative control, without a rule.

Program output (T3). Orbits 214; orbits decided for every indicator 214; orbits dropped 0; U 0.303738317757, interval [0.242990654206, 0.364485981308]; A_min 0.757009345794; mean kappa 0.48042891781; U on S2 (no rule) 0.0. Verdict as printed: **supported**. The three conditions read as printed: 0.3037 >= 0.20; 0.2430 >= 0.10; 0.7570 <= 0.80. A_min is the pair MEGNO-E1.

The 21 pairwise agreements A_XY(b_P) on W2 (T3):

| Pair | A |
|---|---|
| LE-SALI | 0.995327102804 |
| LE-FLI | 0.967289719626 |
| LE-MEGNO | 0.92523364486 |
| LE-E1 | 0.81308411215 |
| LE-E2 | 0.836448598131 |
| LE-WBA | 0.878504672897 |
| SALI-FLI | 0.971962616822 |
| SALI-MEGNO | 0.920560747664 |
| SALI-E1 | 0.817757009346 |
| SALI-E2 | 0.841121495327 |
| SALI-WBA | 0.883177570093 |
| FLI-MEGNO | 0.911214953271 |
| FLI-E1 | 0.836448598131 |
| FLI-E2 | 0.859813084112 |
| FLI-WBA | 0.883177570093 |
| MEGNO-E1 | 0.757009345794 |
| MEGNO-E2 | 0.780373831776 |
| MEGNO-WBA | 0.831775700935 |
| E1-E2 | 0.920560747664 |
| E1-WBA | 0.841121495327 |
| E2-WBA | 0.864485981308 |

By source (T3): panel 125 orbits, U 0.328, A_min 0.752; grid 89 orbits, U 0.269662921348, A_min 0.76404494382. By K, cells with at least 20 orbits: K = 0.5, 99 orbits, U 0.292929292929, A_min 0.777777777778; K = 0.75, 23 orbits, U 0.304347826087, A_min 0.739130434783; K = 0.971635406, 83 orbits, U 0.253012048193, A_min 0.78313253012. The other K cells have 0, 1, 2 or 5 orbits (below min_set); T3 prints values for some of them (for example U 1.0 at K = 1.25 on 1 orbit), and they are not decidable. The frozen rule decides on the pooled W2.

Extension (T13, I_8 on E8): orbits 90; decided for every indicator 90; dropped 0; U 0.277777777778, interval [0.188888888889, 0.377777777778]; A_min 0.766666666667; mean kappa 0.474994986519; U on S2 0.0. Verdict **supported**. Flag **ROBUST**.

Table R2 (T10): U 0.228971962617, interval [0.172897196262, 0.285046728972]; A_min 0.799065420561; mean kappa 0.376311200177; U on S2 0.0. Verdict **supported**. The R2 A_min is 0.799065420561 against the limit 0.80.

### H2 (claim C3, detection-time divergence): NOT SUPPORTED (NEGATIVE)

Frozen rule. M_X is the median over S2 of the censored detection time t_det(X) (the first budget from which the verdict is stable up to N_c, provided the stable verdict equals rho; wrong or missing verdicts get 2 N_c = 204800), for each X in I_7. Q = max_X M_X / min_X M_X. H2 is SUPPORTED if Q >= 10 and the lower end of the 95 percent bootstrap interval of Q is > 3. Otherwise not supported (NEGATIVE).

Program output (T4). S2 orbits 448. M_X = 200.0 for each of LE, SALI, FLI, MEGNO, E1, E2 and WBA. Q 1.0, interval [1.0, 1.0]. Verdict as printed: **not supported**. Both conditions fail (Q < 10; the lower end 1.0 is not > 3).

Outcome shares on S2 (T4): correct 1.0, wrong 0.0, undetected 0.0 for LE, SALI, FLI, MEGNO, E1 and WBA; E2 correct 0.997767857143, wrong 0.0, undetected 0.00223214285714. By source: panel 161 orbits, Q 1.0; grid 287 orbits, Q 1.0. By K: K = 2 pi 368 orbits, Q 1.0; K = 10 80 orbits, Q 1.0; the other K have 0 S2 orbits.

Extension (T13): S2 orbits 303; M_X = 200.0 for all eight indicators including RLI; Q 1.0, interval [1.0, 1.0]. Verdict **not supported**. Flag **ROBUST**.

Table R2 (T10): S2 orbits 448; M_X = 200.0 for all seven; Q 1.0, interval [1.0, 1.0]. Verdict **not supported**.

### H3 (detection latency, regime dependence): NOT SUPPORTED (NEGATIVE)

Frozen rule. R_X = M_X(W2) / M_X(S2), the ratio of the censored median detection times in the weak-chaos and strong-chaos sets, for each X in I_7. H3 is SUPPORTED if for at least 5 of the 7 indicators R_X >= 10 and the lower end of the 95 percent bootstrap interval of R_X is > 3. Otherwise not supported (NEGATIVE).

Program output (T5). W2 orbits 214; S2 orbits 448.

| X | M_W2 | M_S2 | R_X | lower | upper | meets R >= 10 and lower > 3 |
|---|---|---|---|---|---|---|
| LE | 1600.0 | 200.0 | 8.0 | 8.0 | 16.0 | no |
| SALI | 1600.0 | 200.0 | 8.0 | 4.0 | 16.0 | no |
| FLI | 3200.0 | 200.0 | 16.0 | 8.0 | 16.0 | yes |
| MEGNO | 3200.0 | 200.0 | 16.0 | 8.0 | 16.0 | yes |
| E1 | 2400.0 | 200.0 | 12.0 | 4.0 | 32.0 | yes |
| E2 | 400.0 | 200.0 | 2.0 | 1.0 | 4.0 | no |
| WBA | 200.0 | 200.0 | 1.0 | 1.0 | 1.0 | no |

Indicators passing: 3 (the program's count; the last column is this document's reading of the frozen test against the printed values). Needed: 5. Verdict as printed: **not supported**.

By source (T5, point values only, no interval): panel W2 125, S2 161: R LE 8.0, SALI 8.0, FLI 16.0, MEGNO 16.0, E1 8.0, E2 2.0, WBA 1.0. Grid W2 89, S2 287: R LE 8.0, SALI 8.0, FLI 8.0, MEGNO 16.0, E1 16.0, E2 2.0, WBA 1.0. By K: T5 prints R = 1024.0 for all seven indicators at K = 2 pi, where W2 has 1 orbit (S2 368), below min_set; the other K cells print no R.

Extension (T13): W2 90, S2 303. R (lower, upper): LE 8.0 (4.0, 16.0); SALI 8.0 (4.0, 16.0); FLI 8.0 (4.0, 32.0); MEGNO 16.0 (8.0, 32.0); E1 16.0 (3.0, 576.0); E2 2.0 (1.0, 4.0); WBA 1.0 (1.0, 1.0); RLI 2.0 (2.0, 4.0). Indicators passing: 1. Verdict **not supported**. Flag **ROBUST**. The count threshold stays at 5 for eight indicators (Reading P21).

Table R2 (T10): W2 214, S2 448. R (lower, upper): LE 2.0 (1.0, 2.0); SALI 1.0 (1.0, 1.0); FLI 2.0 (1.0, 2.0); MEGNO 1.0 (1.0, 2.0); E1 12.0 (4.0, 32.0); E2 2.0 (1.0, 4.0); WBA 1.0 (1.0, 1.0). Indicators passing: 1. Verdict **not supported**.

### H4 (claim C3, sticky orbits): SUPPORTED

Frozen rule. On T_s (orbits with rho = C and lambda_1(n_s) <= f_s lambda_ref, n_s = 1e5, f_s = 0.2), for each X in I_7 minus LE, phi_X is the fraction of T_s orbits with c_X(N_c) = regular (a false regular verdict at N_c); Delta_phi = max_X phi_X minus min_X phi_X. LE is left out because T_s is defined from lambda_1 (circularity). H4 is SUPPORTED if Delta_phi >= 0.50 and the lower end of its 95 percent bootstrap interval is >= 0.25. Otherwise not supported (NEGATIVE). The a-priori ensemble T_s^design is reported with the same statistic.

Program output (T6). T_s orbits 30; all six indicators decided on 30. T_s^design 96 orbits, all decided on 96.

| X | phi_X (T_s) | phi_X (T_s^design) |
|---|---|---|
| SALI | 0.2 | 0.0520833333333 |
| FLI | 0.266666666667 | 0.0520833333333 |
| MEGNO | 0.233333333333 | 0.0625 |
| E1 | 0.7 | 0.229166666667 |
| E2 | 0.633333333333 | 0.0729166666667 |
| WBA | 0.0666666666667 | 0.0104166666667 |

Delta_phi 0.633333333333, interval [0.466666666667, 0.8]. Verdict as printed: **supported** (0.6333 >= 0.50; 0.4667 >= 0.25). Design ensemble: Delta_phi 0.21875, interval [0.135416666667, 0.302083333333], verdict as printed **not supported**. In both ensembles the largest phi_X is E1 and the smallest is WBA.

By source (T6): panel 15 orbits, Delta_phi 0.866666666667; grid 15 orbits, Delta_phi 0.466666666667 (each below min_set, descriptive). By K (orbits, Delta_phi), all below min_set: K = 0.5 (8) 0.75; K = 0.75 (2) 0.5; K = 0.971635406 (12) 0.416666666667; K = 1.25 (1) 1.0; K = pi/2 (5) 1.0; K = 2 (2) 1.0.

Extension (T13, I_8 on E8): T_s has 15 orbits, below min_set 20. phi: SALI 0.2, FLI 0.266666666667, MEGNO 0.266666666667, E1 0.466666666667, E2 0.533333333333, WBA 0.0666666666667, RLI 0.0; Delta_phi 0.533333333333, interval [0.333333333333, 0.8]. Verdict as printed: **not decidable**. The design ensemble in the extension has 2 orbits (D90: T_s^design is intersected with E8 like the other sets), phi 0.0 for all seven indicators, Delta_phi 0.0 [0.0, 0.0], verdict **not decidable**. Flag **NOT COMPARABLE** (primary supported, extension not decidable). The numbers above are printed for context; they are not a verdict.

Table R2 (T10): T_s 30 orbits. phi: SALI 0.0, FLI 0.0, MEGNO 0.0, E1 0.7, E2 0.633333333333, WBA 0.0666666666667. Delta_phi 0.7, interval [0.566666666667, 0.866666666667]. Verdict **supported**. Design ensemble under R2: 96 orbits, Delta_phi 0.229166666667, interval [0.145833333333, 0.3125], **not supported**.

### H5 (classification stability over the budget): SUPPORTED

Frozen rule. On W2, F_X is the fraction of orbits (decided at b_P and at N_c) with c_X(b_P) different from c_X(N_c), for each X in I_7; F_med is the median of F_X over I_7. H5 is SUPPORTED if F_med >= 0.10 and the lower end of the 95 percent bootstrap interval of F_med is >= 0.05. Otherwise not supported (NEGATIVE).

Program output (T7). W2 orbits 214.

| X | F_X | decided |
|---|---|---|
| LE | 0.177570093458 | 214 |
| SALI | 0.182242990654 | 214 |
| FLI | 0.205607476636 | 214 |
| MEGNO | 0.239436619718 | 213 |
| E1 | 0.200934579439 | 214 |
| E2 | 0.0841121495327 | 214 |
| WBA | 0.0280373831776 | 214 |

F_med 0.182242990654, interval [0.130841121495, 0.228971962617]. Verdict as printed: **supported** (0.1822 >= 0.10; 0.1308 >= 0.05).

By source (T7): panel 125 orbits, F_med 0.176; grid 89 orbits, F_med 0.191011235955. By K, cells with at least 20 orbits: K = 0.5 (99) 0.181818181818; K = 0.75 (23) 0.130434782609; K = 0.971635406 (83) 0.180722891566. T7 also prints F_X(b, W2) for every budget b_0 to b_9; at b_9 = N_c it is 0.0 for all seven (the verdict is compared with itself).

Extension (T13): 90 orbits. F: LE 0.188888888889, SALI 0.188888888889, FLI 0.177777777778, MEGNO 0.233333333333, E1 0.233333333333, E2 0.0777777777778, WBA 0.0333333333333, RLI 0.0222222222222. F_med 0.183333333333, interval [0.111111111111, 0.25]. Verdict **supported**. Flag **ROBUST**.

Table R2 (T10): 214 orbits. F: LE 0.0420560747664, SALI 0.00467289719626, FLI 0.0327102803738, MEGNO 0.00934579439252, E1 0.200934579439, E2 0.0841121495327, WBA 0.0280373831776. F_med 0.0327102803738, interval [0.0140186915888, 0.0607476635514]. Verdict **not supported**.

### H6 (classification stability under direction and precision): NOT SUPPORTED (NEGATIVE)

Frozen rule. Both at N_c. (a) On the designated subset S_D (orbits of E_all in it), for X in {LE, SALI, FLI, MEGNO}, delta_dir(X) is the fraction of orbits whose verdict is not identical in the five deviation directions. (b) On the ladder subset S_L (orbits of E_all in it), for X in I_7, delta_prec(X) is the fraction of orbits whose verdict differs between Float64 and Double64. H6 is SUPPORTED if max_X delta_dir(X) >= 0.05 or max_X delta_prec(X) >= 0.05. Otherwise not supported (NEGATIVE). The rule has no bootstrap interval.

Program output (T8). Orbits S_D* 109; orbits S_L* 51.

| X | delta_dir | dropped (dir) | delta_prec | dropped (prec) |
|---|---|---|---|---|
| LE | 0.0 | 7 | 0.0 | 3 |
| SALI | 0.045871559633 | 0 | 0.0 | 0 |
| FLI | 0.0 | 0 | 0.0 | 0 |
| MEGNO | 0.0 | 9 | 0.0 | 4 |
| E1 | not in the rule | not in the rule | 0.0196078431373 | 0 |
| E2 | not in the rule | not in the rule | 0.0 | 0 |
| WBA | not in the rule | not in the rule | 0.0 | 0 |

Largest delta_dir 0.045871559633 (SALI); largest delta_prec 0.0196078431373 (E1). Both are below 0.05. Verdict as printed: **not supported**. An orbit undecided in some direction or precision is dropped from H6 and counted (Reading P23).

Extension (T13): S_D* 51, S_L* 51. delta_dir: LE 0.0, SALI 0.078431372549, FLI 0.0, MEGNO 0.0 (dropped dir LE 6, MEGNO 7); delta_prec: E1 0.0196078431373, all others 0.0 (dropped prec LE 3, MEGNO 4). Verdict **supported** (SALI 0.0784 >= 0.05). Flag **NOT ROBUST** (primary not supported, extension supported).

Table R2 (T10): S_D* 109, S_L* 51. delta_dir SALI 0.0642201834862, LE, FLI and MEGNO 0.0; delta_prec E1 0.0196078431373, others 0.0; no orbit dropped. Verdict **supported** (SALI 0.0642 >= 0.05).

### Rule-validity gate on the integrable control (T12)

Frozen rule. For each X in I_8 the share of K = 0 gate orbits (80 panel orbits of K index 1, classes T and R, regular by construction) with c_X(N_c) = chaotic must be at most 0.10 (parameter gate_max_share 0.1). A violation is reported as "rule invalid on the integrable control" with the indicator named, and a second version of each decision table without that indicator is added (never in place of the first).

Program output (T12): LE decided 80, chaotic 0; SALI 80, 0; FLI 80, 0; MEGNO 80, 0; E1 decided 75, chaotic 0, undecided 5; E2 80, 0; WBA decided 75, chaotic 0, undecided 5; RLI decided 8, chaotic 0, not applicable 72. Share chaotic 0.0 and violated false for every indicator. The program prints "No indicator violates the gate." No second version of a decision table is therefore needed or made. The gate set is not intersected with E_all, and orbits undecided at the last budget leave the denominator (Reading P26).

### Hypothesis-to-claim mapping stated by the preregistration

The header names claim C3 with the two P2 outputs, detection latency and classification stability. The hypothesis headings read: H1 (claim C3, classification divergence at a fixed budget); H2 (claim C3, detection-time divergence); H3 (detection latency, regime dependence); H4 (claim C3, sticky orbits); H5 (classification stability over the budget); H6 (classification stability under direction and precision). The headings of H3, H5 and H6 name an output but not the label C3; the header places both outputs under C3. The frozen text gives no rule that combines H1 to H6 into one verdict on C3, and none is applied here.

## C. Deviations and readings that affect a confirmatory number

The program lists readings (where the frozen text is silent) and deviations (where it departs from the frozen text, or a frozen statement is false). "Overrides recorded in the manifest: none." The items that bear on a confirmatory number are below, with the decisions each could touch. This document does not say that any of them changed a verdict; the tables do not show it, and no run without them was made.

1. P13 (Reading and Deviation): lambda_ref and lambda_ns of the panel orbits come from the P2 wrapped L1 items at N_1 and at n_s. The preregistration says they equal the P1 L1 values; the tables state that they do not. Size of the difference: not in the tables (the decision log, D89, records that 521 of 883 panel orbits differ by more than 1e-6 and 61 by more than 0.01). lambda_ref decides membership of S2 (>= 1.0), W2 (<= 0.1) and T_s (lambda_1(n_s) <= 0.2 lambda_ref). Could touch: every decision that uses W2, S2 or T_s, that is H1 to H5, in the primary, the extension and R2. H6 does not use these sets.
2. P29 (Deviation): the panel reference labels (rho) are those of the P1 pipeline at the P1 horizon N(o), made with lifted L5 items, while the L1 item of table R1 is wrapped. Could touch: rho defines C and R membership in every set, the correct/wrong outcome behind H2 and H3, and the false-regular statistic phi of H4; through E_all it also bears on H6.
3. E23 b (Deviation): the reference class of grid orbits uses the reduced L1/L5 rule of the preregistration, and every grid orbit carries the flag ref_weak (T1: 288 in each of the 8 grid rows). P14 and P19: a grid orbit is identified by the order of the [[run.h4]] tables with the checks of the stored texts, and its structure axis uses the lifted first coordinate with only the L1-against-L5 test. Could touch: every set that holds grid orbits (W2 89 of 214, S2 287 of 448, T_s 15 of 30, all of E8 beyond the ladder subset); the tables split H1 to H5 by source for this reason.
4. P15 (Reading and Deviation): n(b) is the first recorded step not before b, and n'(n) the recorded step nearest to n/10 (ties to the larger step); the indicators are read at slightly different steps (n(b) up to 1.13 b). The step used is printed in T2. Could touch: every budget-based verdict, so H1 (at b_P), H2 to H5 (detection times, F_X), in R1 and R2 alike.
5. P16 (Deviation): g_SALI = -L, because L (ln GALI_2) is negative and, read literally, every SALI verdict would be undecided. Could touch: every statistic that contains SALI (all of H1 to H6: the A_XY pairs with SALI, F_SALI, phi_SALI, M_SALI, delta_dir SALI). The H6 maximum delta_dir is SALI's.
6. P17 (Reading and Deviation): WBA is undecided at a window if any of its four rows has the flag zero or a non-finite dig; E1 and E2 are undecided on their flag or a non-finite D. An orbit undecided for E1 or WBA at every budget is undetected, with the censored time, and stays in the regular set. Could touch: H2 and H3 (detection times of E1 and WBA) and, through the removal of undecided orbit-budgets, the statistics of E1 and WBA in H1, H4 and H5. The tables print the undecided counts (T2, T12). The E2 undetected share 0.00223214285714 in H2 is printed; the tables do not give its reason.
7. D90 (Reading): in the E8 extension, T_s^design is intersected with E8 like the other sets (owner decision at the stage 04 checkpoint). Could touch: the design statistic of H4 in the extension only (2 orbits). The primary decisions are not touched. The H4 flag is NOT COMPARABLE with or without it, because T_s itself has 15 orbits in the extension.
8. D89w end tolerance (verification, not the analysis): the independent recomputation accepts an interval end of U, Delta_phi and F_med within max(0.05 of the interval width, one lattice step 1/n); one count (R2 H3 indicators passing: module 1, independent 0) is a tolerated count. This is a verify-side rule and changes no number in the tables. The tables do not say which interval ends are involved. The printed lower ends that carry a decision are far from their thresholds (H1 0.2430 against 0.10; H4 0.4667 against 0.25; H5 0.1308 against 0.05); this compares printed numbers and is not a recomputation. The R2 H3 verdict (not supported) is the same for 0 or 1 passing indicator.
9. P20 and P27 (Readings): one set of bootstrap draws per label shared by all statistics of the label; 5 percent dropped replicates make an interval undefined; intervals for the pooled statistics only; a statistic with an empty denominator is undefined and makes its decision not decidable. This bears on the extension of H4 (not decidable) and on the K and source cells, which have point values only.
10. P21 (Reading): the H3 count threshold h3_count_min 5 is kept for eight and for six indicators; a second version is made per violating indicator. No indicator violates the gate, so no second version exists.
11. P18 (Reading): W2 and S2 may overlap through T_s; the overlap is counted in T1 ("W2 and S2" is 0 in every row).
12. The share 0.05 of dropped bootstrap replicates and the design class S of T_s^design are constants of the module, not keys of the 58-key parameter table (a gap of the spec). No threshold is changed by it.

## D. Exploratory notes

Everything in this section is exploratory. None of it is a decision of section B and none has a frozen rule. The preregistration lists the contents of T9, T11, T14, T15, T16 and T17 as reported without a decision rule, and the tables file prints "Exploratory items: none". No computation, test or p-value was made for this section; the values are read from the tables. The tables give no causes, and none is claimed.

**D1 (exploratory, T9): AUC of the chaos scores.** AUC(X, b, K) per K, panel and grid pooled, with (n_C/n_R) in brackets. Four K values are shown (T9 prints all K and budgets); "nothing" means no AUC is printed because a class has fewer than 20 orbits. At b_P = 12800 (k = 6):

| X | K = 0.5 | K = 0.971635406 | K = 1.57079632679 | K = 6.28318530718 |
|---|---|---|---|---|
| LE | 0.999595959596 (99/325) | 0.988703859695 (253/183) | 0.974552517703 (328/124) | 0.9972899729 (369/36) |
| SALI | 0.999876061226 (99/326) | 0.978794152707 (253/189) | 0.952913279133 (328/126) | 0.9972899729 (369/36) |
| FLI | 0.999969015306 (99/326) | 0.997594997595 (253/189) | 1.0 (328/126) | 0.99736525143 (369/36) |
| MEGNO | 0.970621036867 (99/317) | 0.956349206349 (252/172) | 0.976480400334 (327/110) | 0.9972899729 (369/36) |
| E1 | 0.999287352048 (99/326) | 0.999330781939 (253/189) | 0.998524003097 (328/126) | 0.997816922614 (369/36) |
| E2 | 0.998946520419 (99/326) | 0.999435347261 (253/189) | 0.999879016647 (328/126) | 0.999548328817 (369/36) |
| WBA | 0.999597198984 (99/326) | 0.999686304034 (253/189) | 0.999927409988 (328/126) | 0.999849442939 (369/36) |
| RLI | nothing (6/269) | 1.0 (147/129) | 1.0 (220/63) | nothing (295/3) |
| STR | 0.863760302411 (99/326) | 0.961122613297 (253/189) | 0.982989740612 (328/126) | 0.998494429389 (369/36) |

At N_c = 102400 (k = 9):

| X | K = 0.5 | K = 0.971635406 | K = 1.57079632679 | K = 6.28318530718 |
|---|---|---|---|---|
| LE | 0.995884773663 (99/324) | 0.991774383079 (253/185) | 0.977845528455 (328/120) | 0.999774164408 (369/36) |
| SALI | 0.991665117432 (99/326) | 0.980927285275 (253/189) | 0.974932249322 (328/126) | 0.999849442939 (369/36) |
| FLI | 1.0 (99/326) | 1.0 (253/189) | 1.0 (328/126) | 0.99736525143 (369/36) |
| MEGNO | 0.931969537648 (99/317) | 0.912623931813 (253/179) | 0.936001892347 (328/116) | 1.0 (368/35) |
| E1 | 0.999814091839 (99/326) | 0.999916347742 (253/189) | 1.0 (328/126) | 1.0 (369/36) |
| E2 | 0.999938030613 (99/326) | 0.999895434678 (253/189) | 1.0 (328/126) | 1.0 (369/36) |
| WBA | 1.0 (99/326) | 0.999958173871 (253/189) | 1.0 (328/126) | 1.0 (369/36) |
| RLI | nothing (6/269) | 1.0 (147/129) | 1.0 (220/63) | nothing (295/3) |
| STR | 0.864209580467 (99/326) | 0.961697722567 (253/189) | 0.983776132404 (328/126) | 0.997816922614 (369/36) |

The lowest values in these rows are STR at K = 0.5 (0.863760302411 at b_P, 0.864209580467 at N_c) and MEGNO at K = 0.971635406 (0.956349206349 at b_P, 0.912623931813 at N_c).

**D2 (exploratory, T11): sensitivities on S_L (51 orbits per row).** Verdict agreement A and kappa between a sensitivity spec and its primary:

| pair | A at b_P | kappa at b_P | A at N_c | kappa at N_c |
|---|---|---|---|---|
| cos_vs_wba | 1.0 | 1.0 | 1.0 | 1.0 |
| planar_vs_pa1 | 0.588235294118 | 0.0 | 0.588235294118 | 0.0 |
| qp_vs_fq1 | 1.0 | 1.0 | 1.0 | 1.0 |
| raw_vs_fq1 | 0.960784313725 | 0.917874396135 | 0.921568627451 | 0.833333333333 |

Median |dD| at b_P: cos_vs_wba 6.75226967672e-8; planar_vs_pa1 5.13320317084; qp_vs_fq1 0.916825918018; raw_vs_fq1 6.2172489378999994e-15. The planar phase advance has A 0.607843137255 at k = 3 and 0.588235294118 from k = 4 on; the tables do not show why.

**D3 (exploratory, T14): regular range and sticky orbits.** Share of sticky orbits whose score at N_c lies below the median of the regular class (below_median), where printed: 0.0 for LE, SALI, FLI, E1, E2 and WBA at K = 0.5, 0.75, 0.971635406, 1.25, 1.57079632679 and 2.0; MEGNO 0.125 at K = 0.5, 0.5 at K = 2.0, 0.0 at the other four; STR 0.625 at K = 0.5, 0.2 at K = 1.57079632679, 0.0 at the other four. No share is printed ("nothing") at K = 0, 3, 2 pi and 10. Leaving-budget rows exist for each of the 30 sticky orbits per indicator; the rows with no leaving budget printed ("nothing") number LE 3, SALI 2, FLI 0, MEGNO 15, E1 0, E2 0, WBA 0, STR 19 (counted from the T14 rows; the program does not print these counts). The tables do not say more about a missing leaving budget than the M10 definition.

**D4 (exploratory, T15): packages FLI+RLI and MEGNO+SALI, cross-tabulated with rho, no combination rule applied.** On T_s (rho = C), FLI with RLI at b_P: both chaotic 8, first regular and second chaotic 6, both regular 1 (15 orbits); at N_c: both chaotic 11, first regular and second chaotic 4. MEGNO with SALI on T_s at b_P: both chaotic 17, MEGNO chaotic and SALI regular 5, MEGNO regular and SALI chaotic 2, both regular 6 (30 orbits); at N_c: 22, 1, 2, 5. On the regular set (rho = R) at b_P: FLI with RLI both regular 874; MEGNO with SALI both chaotic 4, MEGNO chaotic and SALI regular 27, MEGNO regular and SALI chaotic 3, both regular 1074, MEGNO undecided with SALI chaotic 12, MEGNO undecided with SALI regular 39. At N_c on the regular set: 7, 52, 1, 1066, undecided with chaotic 6, undecided with regular 27 (same order). The regular-set counts differ between the two packages because FLI+RLI is limited to the orbits that carry RLI.

**D5 (exploratory, T16): flags, omissions, controls, K = 0 line.** Flag counts: bk1 below_floor 1358, reldig_undefined 1346, zero 11; fq1 below_floor 1152, degenerate 14, mirror1 1839, mirror2 1885, near_zero1 270, near_zero2 303, refined1 11, refined2 14, tie1 35, tie2 38; pa1 below_floor 1331; pa_planar below_floor 21; the fq1_qp and fq1_raw rows are in T16. Omissions: bk1 and fq1, orbit 883 (the Y5 control has no torus). Controls held to the bound (871, 873, 875, 877, 879, 883): stored deviations from 1.11022302463e-16 to 5.26802579337e-9. Controls marked info (878, 880, 881, 882): deviations 0.222018884744, 0.443891109126, 4.08358102888e-9, 0.103203221967. Controls 872, 874 and 876 have no scope (exponent 0.0, lambda_1(N_1) of order 1e-8). The K = 0 line (largest circle distance of nu1 to p0/(2 pi)) falls from 8.77500818086e-6 at W = 100 to 5.551115123130001e-17 at W = 12800, 25600 and 51200; 5 fq1 orbits are left out as non-finite.

**D6 (exploratory, T17): cost.** Measured ns per orbit-step against the budget line ns: bk1 3000.59036817 against 15000.0; fq1 16850.3797613 against 84300.0 (transform seconds 3024.198199); fq1_qp 8749.99503139 against 26800.0; fq1_raw 9683.55278019 against 26800.0; pa1 1810.47806122 against 12600.0; pa_planar 248.296322696 against 2000.0; pr 223.405047816 against 26000.0; st 1702.27590881 against 2200.0; tg 2658.01590681 against 18000.0. Every measured value is below its line; st is the closest to its line.

**D7 (exploratory, T3, T5, T7): cell and budget patterns already in the tables.** (i) H1 per K cell above min_set: U 0.292929292929, 0.304347826087, 0.253012048193 at K = 0.5, 0.75, 0.971635406. (ii) H3 by source, point values with R >= 10: FLI and MEGNO on the panel, MEGNO and E1 on the grid, against FLI, MEGNO and E1 pooled; no interval is printed by source. (iii) F_X(b, W2) falls with the budget for every indicator; for example LE 0.682692307692 at b_0, 0.177570093458 at b_6, 0.0794392523364 at b_8 and 0.0 at b_9; WBA 0.018691588785 at b_0 to b_3 and 0.0 at b_9.

## E. Limits and what is not done

- Figures F1 to F6 are NOT produced. The figure drawers are a follow-up (Reading P11: the command writes the tables and exits with code 2, "figures skipped"). Only the data blocks of F1 to F6 exist in the tables file: F1 agreement of the pairs of I_8, F2 detection times and outcomes, F3 fraction chaotic per indicator and budget, F4 AUC rows, F5 lattice of the grid, F6 stretching distance by class. This document quotes no figure.
- The independent verification of the analysis (module-build analysis stage 04) matched 14146 of 14146 compared quantities after the stage 04 rework (D89w). The pass after the readings was not blind: the first comparison showed 21 differences, and the T_s^design reading (D90) and the interval-end tolerance (max of 0.05 of the width and one lattice step, for U, Delta_phi and F_med) were settled by the owner after the data were seen. One count (R2 H3 indicators passing) is a tolerated count. The verification report says it did not re-run either bootstrap generator with other seeds.
- One draw of the panel (883 orbits) and one grid (2304 starts on one lattice). The bootstrap intervals resample orbits within these and do not reflect a second draw. Panel and grid orbits are resampled as independent.
- All grid orbits are flagged ref_weak (reduced reference rule), and the panel labels come from P1 at its horizon with different items (section C, items 2 and 3).
- The budgets form a doubling lattice from 200 to 102400, so detection-time medians and their ratios take lattice values (for example 200, 400, 1600, 3200; ratios 8 and 16). The smallest budget is 200: the H2 medians of 200.0 are at that floor, and a detection earlier than 200 cannot be observed.
- Cells with fewer than 20 orbits (K and source cells) are descriptive and not decidable. The K cells other than 0.5, 0.75 and 0.971635406 hold 0 to 5 W2 orbits (T3, T5); T_s holds 8, 2 and 12 orbits at those three K and 1, 5 and 2 at 1.25, pi/2 and 2 (T6). The extension of H4 has 15 orbits (design ensemble 2), so H4 is NOT COMPARABLE.
- Orbits with rho = U are outside every set; the dropped orbits are listed in T1.
- The numbers were copied from the tables file of the run and were not recomputed from the stored series in this step. Section C item 1 quotes one number set from the decision log (D89) that is not in the tables.

## F. Interpretation (for the checkpoint)

This section is a draft for the owner to accept or revise. It is not a decision and it applies no rule that is not in section B. A not supported result is a result, not a failure.

- H1 supported and ROBUST: on weak chaos at b_P = 12800 the seven indicators, under the R1 rules, return split verdicts on 30.4 percent of the 214 orbits (lower end of the interval 0.243), and the least agreeing pair, MEGNO-E1, agrees on 75.7 percent. It also holds under R2 (A_min 0.799, close to the limit 0.80) and with RLI. The strong-chaos negative control has U 0.0.
- H2 not supported and ROBUST: the median detection time on strong chaos is the smallest budget, 200, for every indicator, so this budget list shows no detection-time spread. The result cannot say whether the indicators would differ below 200.
- H3 not supported and ROBUST: the weak-to-strong ratio of detection times reaches the preregistered level for 3 of 7 indicators (FLI, MEGNO, E1) where 5 are needed; under R2 only E1 does. LE and SALI (R = 8.0) are just under the threshold of 10 on a doubling budget lattice. Latency depends on the regime for some indicators and hardly for others (E2, WBA).
- H4 supported on T_s (30 orbits) and under R2, but not supported on the design ensemble (96 orbits) and not comparable in the extension (15 orbits). The spread is driven by E1 (0.7) and E2 (0.633333333333) calling many sticky orbits regular at N_c, against WBA at 0.0666666666667. The conclusion therefore depends on the definition of the sticky set; the paper should say that it rests on one definition (T_s from lambda_1 at n_s = 1e5).
- H5 supported under R1 (F_med 0.182242990654) and ROBUST in the extension, but not supported under R2 (F_med 0.0327102803738, lower end 0.0140186915888). The finding that the finite-budget verdict differs from the verdict at N_c on weak chaos depends on the rule family, growth ratio (R1) or fixed threshold (R2). The primary decision stays R1 as preregistered; the R2 result is a sensitivity and should be shown next to it, not instead of it.
- H6 not supported under R1 (max delta_dir 0.045871559633, max delta_prec 0.0196078431373, both under 0.05), supported under R2 and in the extension, flagged NOT ROBUST. The primary result is a near miss on one indicator (SALI); the paper should not state input-stability of the verdicts as established.
- For C3 as a whole the evidence is mixed: classification divergence at a fixed budget (H1) and sticky-orbit divergence (H4, on its primary set only) are supported; detection-time divergence (H2) and its regime dependence (H3) are not; the stability results (H5, H6) depend on the rule family or the set. The preregistration gives no combination rule, so how to word the claim is the owner's call.
- The deviations of section C (lambda_ref from wrapped L1, labels from different items, ref_weak grid) could shift set membership. These tables do not show whether they matter.
