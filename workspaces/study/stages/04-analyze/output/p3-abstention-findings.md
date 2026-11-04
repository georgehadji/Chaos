# Findings: study p3-abstention (draft for the stage 04 checkpoint)

## A. What this document is

This is a draft for the stage 04 checkpoint. The interpretation in section F is for the owner to confirm or revise. Nothing here is a final decision.

Study p3-abstention tests a fixed abstaining rule R_a against single indicators on the standard map. It has three preregistered hypotheses. HA is the error and abstention of the rule. HB is the false-regular rate at matched cost, with a coverage condition. HC compares indicator disagreement with precision disagreement as predictors of label failure. The value types are FALSIFIES and NEGATIVE allowed: a hypothesis that is not supported is reported with the same figures as a supported one. NOVEL, REPRODUCES and TOOL are not claimed.

The decision rules were frozen on 2026-10-30. Amendments A2, A1 and A3 (all dated 2026-11-04) were made before any test data existed. This document applies the frozen rules without change.

Sections B to D come from one automatic run of `bin/chaoslab.jl analyze` (no human or AI step in it). Every number below is copied from `p3-abstention-tables.md` (the tables file) or from the preregistration. Fractions are quoted as the tables print them, with decimals in parentheses added here for reading only. A number that is not in the tables is marked "not in the tables".

| Item | Value |
|---|---|
| Preregistration header SHA-256 | e4809007aea5ad2302b804f2ba34ebcc838c2d09e6af22559144fa81b9072c77 |
| Results SHA-256 (tables file header) | 452eeabb56d090915ae36da1affe2d7c68c3391a887883fb13068b9ba8a6b98a |
| Calibration record SHA-256 | a8cfdbe12a56697fba17db12b26cafced06647fd05dfa0db79af9a1319b3680a |
| Parameter table (71 keys) SHA-256 | 69965ec24cee507d76924423745cf332ecbb6c003880e5fb4b1036a60fc6a5e1 |
| Test seed | 20261030 (register) |
| Frozen rule theta* | V = SALI, FLI, WBA; q_C = 1; q_R = 3; w = 2; v = 0 |
| Rounds | 5 (closed at r_max = 5) |
| Orbits in the panel file | 6670 |

From the stage 03 run (commit message of the closure commit, not from the tables): 1080 of 1080 chunks ok, none failed or retried, no health flag, no control stop, no override. The tables file agrees on the health side: no L5 flag of any kind (T11), and no chunk that failed twice.

Reading key. R_a is the rule form with indicator evidence only, read at b_P = 12800. rho is the reference class (C chaotic, R regular, U unresolved). E_test is the set of test orbits with rho in {C, R}, K > 0, not drift. W2 is the weak-chaos set, S2 the strong-chaos set, T_s the sticky set (T_s is a subset of W2). D is the Double64 subset, D_B the BigFloat check subset. The gate set is the 80 orbits at K = 0.

Test draw (T2, as printed).

| Set | Size | Target |
|---|---|---|
| E_test | 6433 | n.a. |
| W2 | 1379 | 1235 |
| T_s | 172 | 285 |
| S2 | 820 | n.a. |
| D | 1218 | 45 (this is the HC failure target, see section C) |
| D_B | 100 | n.a. |
| gate | 80 | n.a. |

Labels of the orbits with K > 0 read C 5173, R 1260, U 157. Excluded 0. Dropped as drift 0.

## B. Confirmatory decisions

Order of decision: HA, HB, HC. Decisions as printed by the program: HA not supported, HB supported, HC supported. A not supported result is a result under NEGATIVE and has the same standing as a supported one.

### HA (error and abstention of the rule): NOT SUPPORTED (NEGATIVE)

Frozen rule. HA is supported only if all four items hold. (1) The upper bound of e is at most e_max = 0.03. (2) The upper bound of a(E_test) is at most a_max = 0.15. (3) a(W2) - a(S2) >= delta_a and a(T_s) - a(S2) >= delta_a, with delta_a = 0.1, each with the lower bound of the difference above 0. (4) The gate count g is at most g_max = 1. Bounds are one-sided 95 percent bootstrap bounds, resampled within five strata.

Program output (T5). Strata sizes (T_s, W2 minus T_s, S2, other C, R): 172, 1207, 819, 2975, 1260. Overlap of T_s and S2: 1.

| Item | Statistic | Point | Bound | Threshold | Status |
|---|---|---|---|---|---|
| 1 | e (error among decided) | 9/3028 (0.00297) | upper 604305/146095568 (0.00414) | e_max 0.03 | ok |
| 2 | a(E_test) | 377/6433 (0.0586) | upper 8121/128660 (0.0631) | a_max 0.15 | ok |
| 3 | a(W2) - a(S2) | 29781/565390 (0.0527) | lower 200059901/4653173490 (0.0430) | delta_a 0.1 | fails |
| 3 | a(T_s) - a(S2) | 3647/17630 (0.2069) | lower 5575/35303 (0.1579) | delta_a 0.1 | fails (status of item 3 as a whole; this row meets its condition, see section C item 1) |
| 4 | g (gate orbits chaotic-supported) | 0 | none printed | g_max 1 | ok |

Dropped replicates: 0 for every bound.

Unresolved share a(A) per set (T5):

| Set | a |
|---|---|
| E_test | 377/6433 (0.0586) |
| W2 | 76/1379 (0.0551) |
| T_s | 9/43 (0.2093, equal to 36/172) |
| S2 | 1/410 (0.00244, equal to 2/820) |

Failing item: item 3. Items 1, 2 and 4 hold. The W2 row fails on its point value (0.0527 is below 0.1). The T_s row meets its condition (point 0.2069 and lower bound 0.1579, both above 0.1); the tables print the status of item 3 as a whole on both rows (section C).

Rule outputs (T4, all rounds): E_test chaotic supported 5039, regular supported 1017, unresolved 377. W2: 1289, 14, 76. T_s: 129, 7, 36. S2: 818, 0, 2.

### HB (false-regular rate at matched cost, with coverage): SUPPORTED, flagged underpowered and LE cost-unmatched

Frozen rule. HB is supported if, for every one of the 13 comparisons (6 comparators on T_s, 7 on W2), the one-sided exact McNemar p-value is at most alpha_HB = 0.05 and FR_rule(A) <= rho_HB * FR_X(A) with rho_HB = 0.5, and the coverage cov(T_s) and cov(W2) are each at least c_min = 0.5. LE on W2 is always cost-unmatched and is flagged, not excluded. HB is also reported without LE.

Program output (T6). Final sizes: T_s 172, W2 1379. Underpowered: true. Cost-unmatched LE: true. Defects: none. No comparator was capped. Rule false-regular count: T_s 7/172 (0.0407); W2 2/197 (0.0102, equal to 14/1379).

T_s (comparators I_6 = I_7 minus LE):

| Comparator | b_X | FR_X | n10 | n01 | p | Margin | Pass |
|---|---|---|---|---|---|---|---|
| SALI | 51200 | 27/172 (0.157) | 2 | 22 | 1.79e-5 | true | true |
| FLI | 51200 | 8/43 (0.186) | 2 | 27 | 8.12e-7 | true | true |
| MEGNO | 51200 | 33/172 (0.192) | 5 | 31 | 6.46e-6 | true | true |
| E1 | 6400 | 77/86 (0.895) | 0 | 147 | 5.61e-45 | true | true |
| E2 | 51200 | 75/86 (0.872) | 1 | 144 | 3.27e-42 | true | true |
| WBA | 25600 | 29/172 (0.169) | 0 | 22 | 2.38e-7 | true | true |

W2 (all seven comparators):

| Comparator | b_X | FR_X | n10 | n01 | p | Margin | Pass |
|---|---|---|---|---|---|---|---|
| LE (cost-unmatched) | 102400 | 82/1379 (0.0595) | 14 | 82 | 3.62e-13 | true | true |
| SALI | 51200 | 85/1379 (0.0616) | 9 | 80 | 1.15e-15 | true | true |
| FLI | 51200 | 94/1379 (0.0682) | 8 | 88 | 1.84e-18 | true | true |
| MEGNO | 51200 | 153/1379 (0.111) | 12 | 151 | 4.50e-32 | true | true |
| E1 | 6400 | 305/1379 (0.221) | 0 | 291 | 2.51e-88 | true | true |
| E2 | 51200 | 28/197 (0.142) | 4 | 186 | 3.42e-50 | true | true |
| WBA | 25600 | 45/1379 (0.0326) | 2 | 33 | 1.84e-8 | true | true |

The cost column of T6 gives, in thread-nanoseconds: SALI, FLI, MEGNO 136089600; E1 107842560; E2 92697600; WBA 76815360; LE 21258240. The rule's frozen cost is 72430080 (Amendment A1, item 4; not in the tables). No comparator is capped (capped = false in every row).

Coverage (T6). T_s 34/43 (0.791, price 9/43); W2 1303/1379 (0.945, price 76/1379). Condition met: true. Both are above c_min = 0.5.

Without LE (12 comparisons): supported.

Sensitivities (T7, reported, no effect on any decision). Round down: decidable, supported. Pilot costs: decidable, supported. Measured costs: decidable, supported. Without LE: decidable, supported. Failing comparisons: none listed in any.

Power (T7). Power of the full HB decision at the achieved sizes: 3325 of 20000, power 0.16625, standard error 0.00263. Planned power at the cap: 0.30815 (Amendment A1; table T1). The sizes reached the cap-driven situation: T_s 172 of the target 285, W2 1379 of the target 1235. T_s is above its floor of 100 and W2 above its floor of 300, so HB is decidable and flagged underpowered, as the stopping rule says.

### HC (indicator disagreement against precision disagreement): SUPPORTED

Frozen rule. HC is supported if Delta_AUC >= delta_C = 0.10 and the lower end of its one-sided 95 percent paired bootstrap interval is above 0. It is not decidable if D holds fewer than min_set = 20 failures or non-failures. It is flagged underpowered if failures are below the target n_f = 45.

Program output (T8). Set size 1218. Failures 76. Non-failures 1142. Underpowered: false. Dropped replicates: 0.

| Statistic | Value |
|---|---|
| AUC_I | 73889/86792 (0.8513) |
| AUC_P | 119517/173584 (0.6885) |
| Delta_AUC | 28261/173584 (0.1628) |
| lower bound | 54636833/458855600 (0.1191) |
| upper bound | 776737793/3720042160 (0.2088) |

Delta_AUC 0.1628 is at least 0.10 and the lower bound 0.1191 is above 0. Power at the achieved failures: "not recomputed (target reached)", per the tables. The planned HC power was 0.8437 (Amendment A1).

## C. Deviations and readings that affect a confirmatory number

The tables list four deviations (section "Deviations"). The readings chosen (P87 to P115) are in the tables file and are not repeated here. The items below bear on a confirmatory number. This document does not say that any of them changed a verdict.

1. HA item 3, T_s row. The tables print the status "fails" for a(T_s) - a(S2). The printed point is 3647/17630 (0.2069) and the printed lower bound is 5575/35303 (0.1579). Both are above delta_a = 0.1, and the lower bound is above 0. Read as the preregistration is worded, this row would pass. The status column of table T5 is the status of item 3 as a whole, printed on both of its rows (the code marks a row "fails" when item 3 is in the failed list), so the T_s row inherits the failure of the W2 row. The T_s row meets its condition. The HA decision does not change: the W2 row fails on its own (point 0.0527 below 0.1), and item 3 needs both rows.
2. P120 (Deviation). The "labels agree" column of T3 is omitted. The chk items run at their own horizon N_B = 1000000 and have no steps at the scale of N_1, so labels cannot be formed in three precisions. The tables say study amendment A2 states this. T3 is reported only and decides nothing.
3. P95 (Deviation). The HA strata are not disjoint by definition. An orbit can be in both T_s and S2, and the first stratum it qualifies for holds it. The overlap is 1 orbit. This is why the S2 stratum has 819 orbits and the set S2 has 820.
4. P100 (Deviation). A gate orbit without verdict items would make HA not decidable and be listed. HA was decidable and g = 0 is printed, so the tables show no such orbit.
5. P104 (Deviation). D orbits without Double64 verdicts of the voting set leave the HC set. The tables do not print how many left. T2 gives |D| = 1218 and T8 gives a set size of 1218, so by these two numbers none left. This is an inference, not a printed count.
6. The D "target" of 45 in T2 is the HC failure target n_f, not an orbit count. The table puts it in the size column.
7. The preregistration expected the cap to bind for HB, and it did (T_s 172 of 285). The preregistered planning model assumed the yields of the p2 panel. The tables do not say why T_s fell short while W2 exceeded its target.
8. Per-K and per-source tables of HA, HB and HC are listed in the preregistration under "Reported without a decision rule". The tables file has per-K and per-round counts only for rule outputs (T4). Per-K and per-source results for HA, HB and HC: not in the tables.

## D. Exploratory and reported-only results

Everything in this section is reported without a decision rule. The tables file prints "Exploratory items: none"; the items below are the ones the preregistration lists as reported without a rule. No computation or test was added for this document.

**D1 (T3): reference check on D_B.** The three L5 values (Float64, Double64, BigFloat) at N_B = 1000000 and the ladder test, for 100 orbits. Orbits on which the ladder would have made rho unresolved: C 80, R 0 (as printed). Two orbits (404 and 6101) are marked "dropped" in the rho column. No label was changed (D1). The values differ between precisions by several percent on many orbits (for example orbit 312: 0.0361, 0.0737, 0.0693). The tables give no cause.

**D2 (T9): R_a against R_ab on D (1218 orbits).** Chaotic supported under both: 1045. Regular supported under both: 56. Unresolved under both: 49. Chaotic supported under R_a and unresolved under R_ab: 68 (the ladder overrides). All other pairs: 0. Error of the decided orbits: R_a 8/1169 (0.00684), R_ab 8/1101 (0.00727). The ladder removed 68 decided orbits and no error.

**D3 (T10): cost per orbit-step by spec (thread-ns), measured, cost-model key, pilot key.** L1 77.85, 207.6, 168.0. tg 241.46, 2658.0, 300.0. pa1 275.02, 1810.5, 300.0. bk1 450.77, 3000.6, 504.0. fq1 8860.77, 16850.4, 12609.0. Every measured value is below both the frozen cost key and the pilot key.

**D4 (F2 data): risk and coverage.** Digest 8bbef40c258c1f5b016847376c8a04a2e123248def6e4a2307409f7ceb470d48.

| Curve | Point | Decided share | Error |
|---|---|---|---|
| rule | q_C=1 q_R=3 (theta*) | 0.9414 | 0.00297 |
| rule | q_C=2 q_R=3 | 0.8470 | 0.00294 |
| rule | q_C=1 q_R=2 | 0.9770 | 0.00875 |
| single | LE | 0.9868 | 0.01875 |
| single | SALI | 1.0 | 0.03746 |
| single | FLI | 1.0 | 0.02176 |
| single | MEGNO | 0.9559 | 0.07725 |
| single | E1 | 1.0 | 0.11254 |
| single | E2 | 1.0 | 0.03979 |
| single | WBA | 1.0 | 0.00933 |

**D5 (F3 data): ROC points.** Digest d91bb69181466301d50ddf2f270be7de8084e9bbc87592d11ca25df210e80600. Printed AUC: s_I 0.8513, s_P 0.6885. At threshold 2/3: s_I false-positive rate 0.0210, true-positive rate 0.7237; s_P false-positive rate 0.0236, true-positive rate 0.3947.

**D6 (T1): calibration record, transcribed.** Cal sizes: cal 3016 (panel 831, grid 2185), S2 452, T_s 31, W2 228 (R_W 197), regular 1166, gate 80. This is fitted data from p2 and is not evidence.

## E. Limits and what is not done

- No figures were drawn. The analysis writes the data of F2 and F3 with digests only (reading P110). There is no PDF in this folder, so there is no owner visual check to record. F1 (calibration landscape) belongs to the calibration tables, not to this run.
- Power. HB was run at 172 T_s and 1379 W2 orbits. Power at these sizes is 0.16625 (planned 0.30815 at the cap). HB passed all 13 comparisons, but this is not evidence that the test was powerful.
- The orbits are treated as independent in the bootstrap and in the McNemar test. Orbits of one cell and round are not independent. The preregistration states this limit.
- HA, HB and HC use one system, the standard map, and K values that also occur in the calibration data. Nothing is claimed for other maps or for flows.
- Adaptive allocation is not tested. HB compares a fixed-budget rule with single indicators only.
- Reference labels are time-bound (horizon N_1 = 1e7) and use the reduced rule with no ladder and no direction test. The BigFloat check on 100 orbits (D1 above) shows that the ladder would have made rho unresolved on 80 of them. This is a reported number, not a change of any label.
- The test draw is enriched for sticky and weak-chaos orbits by design. The rates a(E_test) and e are not population rates.
- s_I and the label failure f are computed from the same Float64 votes, so AUC_I has a mechanical part (preregistration, Known limits). HC measures the use of that split against precision disagreement. It is not an independent signal.
- The matched cost uses frozen p2 values. The test run measured costs far below them (D3). The sensitivities with pilot and measured costs are supported (T7).
- Calibration is not blind. The fitted rule performance on p2 is optimistic. The confirmatory evidence is the test.
- All numbers were copied from the tables file and the preregistration. They were not recomputed from the stored series in this step. Decimals in parentheses were computed by division for reading only.
- No independent verification report for stage 04 is cited in this draft.

## F. Interpretation (for the checkpoint)

This section is a draft for the owner to accept or revise. It is not a decision and applies no rule that is not in section B. A not supported result is a result, not a failure.

- HA not supported, on item 3. The rule is accurate when it decides: error 9/3028 (0.0030), upper bound 0.0041 against 0.03. It abstains on 5.9 percent of E_test (upper bound 0.063 against 0.15), and it never called a gate orbit chaotic. The failure is about where it abstains. On W2 it abstains on 5.5 percent of orbits against 0.24 percent on S2. That is a gap of 0.053, about half the 0.10 required. So on W2 the rule abstains only a little more than on S2. On T_s the printed abstention is 20.9 percent, a gap of 0.207. That is more than twice delta_a, so the T_s row meets its condition; the "fails" mark on that row is the status of item 3 as a whole (section C item 1). The decision rests on the W2 row alone. One reading for the owner: most W2 orbits are decided, and the rule abstains visibly mainly on the sticky part of W2.
- HB supported, flagged underpowered and LE cost-unmatched. The rule is false-regular on 4.1 percent of T_s and about 1.0 percent of W2. The closest comparator is WBA at 25600 steps, with 16.9 percent on T_s and 3.3 percent on W2. On W2 the cost-unmatched LE (5.9 percent) is also well above the rule. The direction is the same on all 13 comparisons. E1 and E2 are false-regular on 87 to 90 percent of T_s orbits. The coverage price is 20.9 percent of T_s and 5.5 percent of W2 left unresolved, which is the same abstention as in HA. The rule gets its low false-regular rate partly by abstaining on 9 of 43 T_s orbits, and the HB coverage condition (0.5) is far from binding.
- The caveat on power is important. HB passed on all 13 comparisons, but the power at the achieved sizes is 0.166, and T_s is 172 of the 285 planned. The planned power at the cap was already only 0.308. A pass at low power does not show that the test was powerful. It shows that this draw gave large gaps. The p-values are small (the largest is 1.79e-5, on SALI in T_s), but a significant result from an underpowered test can overstate the true gap, and the tables do not give a way to separate the two. HB should be worded as supported with the underpowered flag, and next to the four sensitivities, which agree.
- HC supported. Indicator disagreement (AUC 0.851) predicts label failure better than precision disagreement (AUC 0.689). The gap is 0.163, with a lower bound of 0.119. There were 76 failures against a target of 45, so HC is not underpowered. Two cautions. The planned AUC_P was 0.55 (nearly constant score), but the observed AUC_P is 0.689, so precision disagreement carried more signal than planned. And the Known limits say s_I and f share the same Float64 votes, so part of AUC_I is mechanical. HC should not be worded as an independent signal.
- Reported without a rule. The ladder overrides (T9) turned 68 chaotic-supported outputs into unresolved on D and removed none of the 8 errors. The cost measured in the run is far below the frozen cost key. The BigFloat check shows large differences between precisions in L5 on many D_B orbits.
- How to word the three results together is the owner's call. The preregistration gives no rule that combines HA, HB and HC.
