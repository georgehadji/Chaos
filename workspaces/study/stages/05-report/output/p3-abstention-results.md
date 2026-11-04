# Results: study p3-abstention

## Pointer key

Every number carries a pointer in square brackets.

- `[T1]` to `[T11]`, `[F2]`, `[F3]`: sections of `workspaces/study/stages/04-analyze/output/p3-abstention-tables.md`, with the row or statistic after the comma. In T4 the set, K index and round come first. In T5 the item number and statistic come first. In T6 the set and comparator come first.
- `[Decisions]`, `[Readings]`, `[Deviations]`: the blocks of that name in the tables file.
- `[findings X]`: section X of `workspaces/study/stages/04-analyze/output/p3-abstention-findings.md` (owner-approved).
- `[manifest: field]`: a field of `workspaces/study/stages/03-run/output/p3-abstention-manifest.json`.
- `[prereg: part]`: the frozen preregistration `workspaces/study/stages/01-design/output/p3-abstention-preregistration.md`.
- Decimals in parentheses are division results added for reading only. The fractions are as the tables print them.

## Header

Study slug: p3-abstention [manifest: study].

| Item | Value | Pointer |
|---|---|---|
| Preregistration header SHA-256 | e4809007aea5ad2302b804f2ba34ebcc838c2d09e6af22559144fa81b9072c77 | [Decisions header], [manifest: prereg.header] |
| Preregistration SHA-256 at start and end of the run | equal to the header hash | [manifest: prereg.sha256_start, prereg.sha256_end] |
| Panel file SHA-256 at start and end of the run | d9a365603c9710cc960ebd7a3e17fd7c9068d7255f16c32d0eaad4b67be6f0ce (equal) | [manifest: panel.sha256_start, panel.sha256_end] |
| Results SHA-256 | 452eeabb56d090915ae36da1affe2d7c68c3391a887883fb13068b9ba8a6b98a | [Decisions header] |
| Calibration record SHA-256 | a8cfdbe12a56697fba17db12b26cafced06647fd05dfa0db79af9a1319b3680a | [T1, calibration_record_sha256] |
| Parameter table SHA-256 (71 keys) | 69965ec24cee507d76924423745cf332ecbb6c003880e5fb4b1036a60fc6a5e1 | [T1, p3_table_sha256] |
| Round log SHA-256 | fe7a9697fc9ef163039eca016f17af54dfa0fd3292e5766cd6a923641a4b5a15 | [manifest: rounds.log_sha256] |
| Run git commit (tree clean) | f3adbe133484ee98c38ad697ce9780d864c68a2e | [manifest: git.sha, git.clean] |
| Module version, language | 0.1.2, Julia 1.13.2 | [manifest: versions.module, versions.language] |
| Test seed | 20261030 | [findings A] |
| Run status | complete | [manifest: status] |
| Chunks | 1080, all with status ok | [manifest: plan_counts.chunks, chunks[].status] |
| Items | 45148 | [manifest: plan_counts.items] |
| Sessions | 6, all clean, 23 threads each | [manifest: sessions[].clean, sessions[].threads_used] |
| Total wall time | 2567.266443490982 s (cap 12.0 h) | [manifest: wall.total_seconds, wall.cap_hours] |
| Control stops, overrides | none, none | [manifest: control_stops, overrides] |
| Backend | cpu (CUDA.jl not loadable for the auto requests) | [manifest: device, backend_table] |
| Host | 24 cores, arrowlake-s | [manifest: host] |

What the study tested. Study p3-abstention tests a fixed abstaining rule R_a against single indicators on the standard map. It has three hypotheses: HA (error and abstention of the rule), HB (false-regular rate at matched cost, with a coverage condition) and HC (indicator disagreement against precision disagreement as predictors of label failure) [findings A]. Value types are FALSIFIES and NEGATIVE allowed. NOVEL, REPRODUCES and TOOL are not claimed [findings A]. The decision rules were frozen on 2026-10-30. Amendments A2, A1 and A3 (all dated 2026-11-04) came before any test data existed. The rules are applied here without change [findings A].

Frozen rule theta*: V = SALI, FLI, WBA; q_C = 1; q_R = 3; w = 2; v = 0 [T1, theta_voters to theta_v]. Indicators used: SALI [@Skokos2001], FLI [@Froeschle1997], MEGNO [@CincottaSimo2000], E1 as the frequency path of frequency map analysis [@Laskar2003], WBA as the weighted Birkhoff average [@SanderMeiss2020; @DuignanMeiss2023], and LE as the Benettin estimator [@Benettin1980a]. E2 is the weighted phase advance path. No verified key for E2 is available here, so it is not cited. All thresholds of the rules are set by this study; none comes from the literature [prereg: Bibliography note].

Reading key. rho is the reference class (C chaotic, R regular, U unresolved). E_test is the set of test orbits with rho in {C, R}, K > 0, not drift. W2 is the weak-chaos set, S2 the strong-chaos set, T_s the sticky set (a subset of W2). D is the Double64 subset and D_B the BigFloat check subset. The gate set is the 80 orbits at K = 0 [findings A].

## The test draw

The panel file holds 6670 orbits in 5 rounds [T2 header]. The round log closed at round 5, which is r_max [manifest: rounds.closed, rounds.log_entries; findings A]. Labels of the orbits with K > 0 read C 5173, R 1260, U 157. Excluded 0. Dropped as drift 0 [T2, labels line].

| Set | Size | Target | Pointer |
|---|---|---|---|
| E_test | 6433 | n.a. | [T2, E_test] |
| W2 | 1379 | 1235 | [T2, W2], [T1, target_w2] |
| T_s | 172 | 285 | [T2, T_s], [T1, target_t_s] |
| S2 | 820 | n.a. | [T2, S2] |
| D | 1218 | 45 (the HC failure target, not an orbit count) | [T2, D], [T1, target_f] |
| D_B | 100 | n.a. | [T2, D_B] |
| gate | 80 | n.a. | [T2, gate] |

The W2 target was exceeded and the T_s target was not reached. The tables do not say why [findings C, item 7]. The manifest records |D| = 1232 and |D_B| = 100 for the allocation [manifest: rounds.D_size, rounds.DB_size]. The set D in the tables has 1218 orbits [T2, D]: it is the orbits of E_test with the panel flag in_D (`src/Analysis/p3t_data.jl`, E64 c). Of the 1232 flagged orbits, 4 are at K = 0 (gate orbits, outside E_test) [panel file, flag in_D]; the other 10 are outside E_test (rho U or drift), a number the tables do not print by itself.

## Confirmatory results

Order of decision: HA, HB, HC. Decisions as printed: HA not supported, HB supported, HC supported [Decisions]. A not supported result is a result under NEGATIVE and has the same standing as a supported one [findings B].

| Hypothesis | Decision | Flags | Pointer |
|---|---|---|---|
| HA | NOT SUPPORTED (NEGATIVE) | fails on item 3 | [Decisions; T5] |
| HB | SUPPORTED | underpowered; LE cost-unmatched | [Decisions; T6; T7] |
| HC | SUPPORTED | not underpowered | [Decisions; T8] |

### HA (error and abstention of the rule): NOT SUPPORTED (NEGATIVE)

Rule. HA is supported only if all four items hold. (1) The upper bound of e is at most e_max = 0.03. (2) The upper bound of a(E_test) is at most a_max = 0.15. (3) a(W2) - a(S2) >= delta_a and a(T_s) - a(S2) >= delta_a, with delta_a = 0.1, each with the lower bound of the difference above 0. (4) The gate count g is at most g_max = 1. Bounds are one-sided 95 percent bootstrap bounds, resampled within five strata [findings B].

Measured. Strata sizes (T_s, W2 minus T_s, S2, other C, R) are 172, 1207, 819, 2975, 1260. The overlap of T_s and S2 is 1 orbit [T5 header]. Dropped replicates are 0 for every bound [T5].

| Item | Statistic | Point | Bound | Threshold | Status | Pointer |
|---|---|---|---|---|---|---|
| 1 | e (error among decided) | 9/3028 (0.00297) | upper 604305/146095568 (0.00414) | e_max 0.03 | ok | [T5, 1] |
| 2 | a(E_test) | 377/6433 (0.0586) | upper 8121/128660 (0.0631) | a_max 0.15 | ok | [T5, 2] |
| 3 | a(W2) - a(S2) | 29781/565390 (0.0527) | lower 200059901/4653173490 (0.0430) | delta_a 0.1 | fails | [T5, 3, a(W2) - a(S2)] |
| 3 | a(T_s) - a(S2) | 3647/17630 (0.2069) | lower 5575/35303 (0.1579) | delta_a 0.1 | printed as fails; see below | [T5, 3, a(T_s) - a(S2)] |
| 4 | g (gate orbits chaotic-supported) | 0 | none printed | g_max 1 | ok | [T5, 4] |

Failing item: item 3. The W2 row fails on its point value: 0.0527 is below delta_a = 0.1. The T_s row meets its condition. Its point 0.2069 and its lower bound 0.1579 are both above 0.1. The status column of T5 is the status of item 3 as a whole, printed on both rows, so the T_s row shows "fails" because the W2 row fails [findings C, item 1]. Item 3 needs both rows. The HA decision does not change [findings C, item 1].

Unresolved share a(A) per set: E_test 377/6433 (0.0586), W2 76/1379 (0.0551), T_s 9/43 (0.2093, equal to 36/172), S2 1/410 (0.00244, equal to 2/820) [T5, unresolved share table; findings B].

Rule outputs (chaotic supported, regular supported, unresolved), all rounds: E_test 5039, 1017, 377; W2 1289, 14, 76; T_s 129, 7, 36; S2 818, 0, 2 [T4, set, all, all].

Reading under FALSIFIES. The rule is accurate when it decides: error 0.0030 with an upper bound of 0.0041 against 0.03 [T5, 1]. It abstains on 5.9 percent of E_test [T5, 2] and it never called a gate orbit chaotic [T5, 4]. The failure is where it abstains. On W2 the rule abstains on 5.5 percent of orbits against 0.24 percent on S2, a gap of 0.0527, about half of the 0.1 required [T5, 3]. On T_s the printed abstention is 20.9 percent, a gap of 0.2069 [T5, unresolved share, T_s; T5, 3]. One reading for the owner: most W2 orbits are decided, and the rule abstains visibly mainly on the sticky part of W2 [findings F].

### HB (false-regular rate at matched cost, with coverage): SUPPORTED, flagged underpowered, LE cost-unmatched

Rule. HB is supported if, for each of the 13 comparisons (6 comparators on T_s, 7 on W2), the one-sided exact McNemar p-value is at most alpha_HB = 0.05 and FR_rule(A) <= rho_HB * FR_X(A) with rho_HB = 0.5. In addition cov(T_s) and cov(W2) are each at least c_min = 0.5. LE on W2 is always cost-unmatched. It is flagged and not excluded. HB is also reported without LE [findings B].

Measured. Final sizes: T_s 172, W2 1379. Underpowered: true. Cost-unmatched LE: true. Defects: none [T6 header]. No comparator was capped [T6, capped column]. Rule false-regular count: T_s 7/172 (0.0407), W2 2/197 (0.0102, equal to 14/1379) [T6, FR_rule column].

T_s comparators (n10 and n01 are the discordant counts of the McNemar table; "Margin" is the FR_rule <= 0.5 FR_X test):

| Comparator | b_X | FR_X | n10 | n01 | p | Margin | Pass | Pointer |
|---|---|---|---|---|---|---|---|---|
| SALI | 51200 | 27/172 (0.157) | 2 | 22 | 1.79e-5 | true | true | [T6, T_s, SALI] |
| FLI | 51200 | 8/43 (0.186) | 2 | 27 | 8.12e-7 | true | true | [T6, T_s, FLI] |
| MEGNO | 51200 | 33/172 (0.192) | 5 | 31 | 6.46e-6 | true | true | [T6, T_s, MEGNO] |
| E1 | 6400 | 77/86 (0.895) | 0 | 147 | 5.61e-45 | true | true | [T6, T_s, E1] |
| E2 | 51200 | 75/86 (0.872) | 1 | 144 | 3.27e-42 | true | true | [T6, T_s, E2] |
| WBA | 25600 | 29/172 (0.169) | 0 | 22 | 2.38e-7 | true | true | [T6, T_s, WBA] |

W2 comparators:

| Comparator | b_X | FR_X | n10 | n01 | p | Margin | Pass | Pointer |
|---|---|---|---|---|---|---|---|---|
| LE (cost-unmatched) | 102400 | 82/1379 (0.0595) | 14 | 82 | 3.62e-13 | true | true | [T6, W2, LE] |
| SALI | 51200 | 85/1379 (0.0616) | 9 | 80 | 1.15e-15 | true | true | [T6, W2, SALI] |
| FLI | 51200 | 94/1379 (0.0682) | 8 | 88 | 1.84e-18 | true | true | [T6, W2, FLI] |
| MEGNO | 51200 | 153/1379 (0.111) | 12 | 151 | 4.50e-32 | true | true | [T6, W2, MEGNO] |
| E1 | 6400 | 305/1379 (0.221) | 0 | 291 | 2.51e-88 | true | true | [T6, W2, E1] |
| E2 | 51200 | 28/197 (0.142) | 4 | 186 | 3.42e-50 | true | true | [T6, W2, E2] |
| WBA | 25600 | 45/1379 (0.0326) | 2 | 33 | 1.84e-8 | true | true | [T6, W2, WBA] |

Matched costs in thread-nanoseconds: SALI, FLI and MEGNO 136089600; E1 107842560; E2 92697600; WBA 76815360; LE 21258240 [T6, cost column]. All 13 comparisons pass. The largest p-value is 1.79e-5 (SALI on T_s) [T6, T_s, SALI].

Coverage. T_s 34/43 (0.791, price 9/43); W2 1303/1379 (0.945, price 76/1379). The condition is met and both values are above c_min = 0.5 [T6, coverage line].

Without LE (12 comparisons): supported [T6, without LE line].

Sensitivities (reported, no effect on any decision). Round down, pilot costs, measured costs and without LE are each decidable and supported. No comparison fails in any [T7, sensitivities table].

Power. Power of the full HB decision at the achieved sizes is 3325 of 20000, that is 0.16625, with standard error 0.00263259318449 [T7, power line]. The planned power at the cap was 0.30815 [T1, power_hb_planned; T7, power line]. T_s reached 172 of the target 285 [T2, T_s]. W2 reached 1379 of the target 1235 [T2, W2]. T_s is above its floor of 100 and W2 above its floor of 300, so HB is decidable and flagged underpowered, as the stopping rule says [findings B].

Reading under FALSIFIES. The direction is the same on all 13 comparisons. The nearest comparator is WBA at 25600 steps, false-regular on 29/172 of T_s and 45/1379 of W2 [T6, T_s, WBA; T6, W2, WBA], against 7/172 and 2/197 for the rule [T6, FR_rule]. The cost-unmatched LE is also well above the rule on W2 (82/1379) [T6, W2, LE]. E1 and E2 are false-regular on 77/86 and 75/86 of the T_s orbits they decide on [T6, T_s, E1, E2]. The rule gets part of its low false-regular rate by abstaining on 9 of 43 T_s orbits [T6, coverage line; findings F]. A pass at low power does not show that the test was powerful, and a significant result from an underpowered test can overstate the true gap [findings F]. HB is worded as supported with the underpowered flag, next to the four sensitivities, which agree.

### HC (indicator disagreement against precision disagreement): SUPPORTED

Rule. HC is supported if Delta_AUC >= delta_C = 0.10 and the lower end of its one-sided 95 percent paired bootstrap interval is above 0. It is not decidable if D holds fewer than min_set = 20 failures or non-failures. It is flagged underpowered if failures are below n_f = 45 [findings B].

Measured. Set size 1218, failures 76, non-failures 1142, underpowered false, dropped replicates 0 [T8 header].

| Statistic | Value | Pointer |
|---|---|---|
| AUC_I | 73889/86792 (0.8513) | [T8, AUC_I] |
| AUC_P | 119517/173584 (0.6885) | [T8, AUC_P] |
| Delta_AUC | 28261/173584 (0.1628) | [T8, Delta_AUC] |
| lower bound | 54636833/458855600 (0.1191) | [T8, lower bound] |
| upper bound | 776737793/3720042160 (0.2088) | [T8, upper bound] |

Delta_AUC 0.1628 is at least 0.10 and the lower bound 0.1191 is above 0. There were 76 failures against the target 45 [T8 header; T1, target_f]. Power at the achieved failures was not recomputed because the target was reached [T8]. The planned HC power was 0.8437 [T1, power_hc_planned]. The F3 data give AUC 0.851334224352 for s_I and 0.688525440133 for s_P [F3, AUC rows].

Reading under FALSIFIES. Indicator disagreement predicts label failure better than precision disagreement. Two cautions. The planned AUC_P was 0.55 [findings F]; the observed value is 0.6885, so precision disagreement carried more signal than planned. And s_I and the label failure f are computed from the same Float64 votes, so part of AUC_I is mechanical [prereg: Known limits; findings F]. HC is not worded as an independent signal.

## Reported without a decision rule

These items carry no decision. The tables print "Exploratory items: none" [Exploratory items]. The items below are those the preregistration lists as reported without a rule [prereg: Reported without a decision rule].

**Reference check on D_B [T3].** L5 at N_B = 1000000 in Float64, Double64 and BigFloat for 100 orbits. Orbits on which the ladder would have made rho unresolved: C 80, R 0 [T3, last line]. Two orbits (404 and 6101) are marked "dropped" in the rho column [T3, rows 404, 6101]. No label was changed [T3 header]. The L5 values differ between precisions by several percent on many orbits; for example orbit 312 reads 0.0360795142496, 0.0736856786269 and 0.0693246545405 [T3, row 312]. The tables give no cause. The column "labels agree" is omitted (deviation P120).

**R_a against R_ab on D [T9].** Out of 1218 orbits: chaotic supported under both 1045; regular supported under both 56; unresolved under both 49; chaotic supported under R_a and unresolved under R_ab 68 (ladder overrides); all other pairs 0. Error of the decided orbits: R_a 8/1169, R_ab 8/1101. The ladder removed 68 decided orbits and no error [T9].

**Cost per orbit-step [T10].** Measured, cost-key and pilot-key values in thread-ns: L1 77.85, 207.6, 168.0; tg 241.46, 2658.0, 300.0; pa1 275.02, 1810.5, 300.0; bk1 450.77, 3000.6, 504.0; fq1 8860.77, 16850.4, 12609.0 [T10]. Every measured value is below both keys. Per arm, the manifest records cost_ns_per_orbit_step: ref 113.93, ind 3184.33, ind_pa1 275.02, d64 6990.70, d64_pa1 2004.35, chk 133744.31 [manifest: arms.*.cost_ns_per_orbit_step]. These are arm averages and mix several specs; the per-spec values are in T10.

**Run cost.** Thread-seconds per arm: ref 15198.18, ind 6524.77, ind_pa1 187.84, d64 2645.77, d64_pa1 252.86, chk 26748.86 [manifest: arms.*.thread_seconds]. The sum is 51558.3 thread-seconds (14.3 thread-hours; sum computed here). All six arms have status complete [manifest: arms.*.status]. Steps: ref 133400000000, ind 2049024000, ind_pa1 683014670, d64 378470400, d64_pa1 126158032, chk 200000000 [manifest: arms.*.steps]. Wall seconds: ref 660.79, ind 283.69, ind_pa1 8.17, d64 115.03, d64_pa1 10.99, chk 1325.88 [manifest: arms.*.wall_seconds].

**Risk and coverage [F2].** Digest 8bbef40c258c1f5b016847376c8a04a2e123248def6e4a2307409f7ceb470d48 [F2].

| Curve | Point | Decided share | Error | Pointer |
|---|---|---|---|---|
| rule | q_C=1 q_R=3 (theta*) | 0.9414 | 0.00297 | [F2, row 1] |
| rule | q_C=2 q_R=3 | 0.8470 | 0.00294 | [F2, row 2] |
| rule | q_C=1 q_R=2 | 0.9770 | 0.00875 | [F2, row 3] |
| single | LE | 0.9868 | 0.01875 | [F2, LE] |
| single | SALI | 1.0 | 0.03746 | [F2, SALI] |
| single | FLI | 1.0 | 0.02176 | [F2, FLI] |
| single | MEGNO | 0.9559 | 0.07725 | [F2, MEGNO] |
| single | E1 | 1.0 | 0.11254 | [F2, E1] |
| single | E2 | 1.0 | 0.03979 | [F2, E2] |
| single | WBA | 1.0 | 0.00933 | [F2, WBA] |

**ROC points [F3].** Digest d91bb69181466301d50ddf2f270be7de8084e9bbc87592d11ca25df210e80600. At threshold 2/3, s_I has false-positive rate 0.0210 and true-positive rate 0.7237; s_P has 0.0236 and 0.3947 [F3, threshold 0.666666666667].

**Calibration record [T1].** Pooled sizes: cal 3016, S2 452, T_s 31, W2 228 (R_W 197), regular 1166, gate 80 [T1, pooled column]. This is fitted data from p2 and is not evidence [findings D6].

**Health [T11].** No flag of any kind for the ref arm (diagnostics exceeded 0, drift 0, excluded 0) and no chunk failed twice [T11]. This agrees with the manifest: 1080 chunks, all ok [manifest: chunks[].status].

## Audit trail

- The analysis is the output of one automatic run of `bin/chaoslab.jl analyze` on the closed round log of 5 rounds, with no human or AI step in it [findings A].
- A second analyze run into an empty folder gave byte-identical tables and results [stage 04 audit, commit 984fcf3 message]. The data tree was unchanged.
- The preregistration and panel hashes were equal at the start and end of every session [manifest: prereg.sha256_start, prereg.sha256_end, panel.sha256_start, panel.sha256_end; sessions[].clean].
- The numbers in this document were copied from the tables, the findings and the manifest. They were not recomputed from the stored series in this step [findings E].

## Limitations

- HB is underpowered. Power at the achieved sizes is 0.16625, against 0.30815 planned at the cap [T7, power line; T1, power_hb_planned]. HB passed all 13 comparisons, but this does not show that the test was powerful [findings E].
- The T_s target was not reached: 172 of 285 [T2, T_s]. The tables give no reason [findings C, item 7].
- LE on W2 is cost-unmatched. Its matched cost is 21258240 thread-ns, below the other comparators [T6, W2, LE]. HB is supported without LE [T6, without LE line].
- The matched cost uses frozen p2 values. The measured costs in the run are far below them [T10]. The pilot and measured sensitivities are supported [T7].
- No figure was drawn. The analysis writes the data of F2 and F3 with digests only (reading P110) [Readings, P110]. There is no figure file in the output folder and no owner visual check.
- Per-K and per-source results for HA, HB and HC are not in the tables. The preregistration lists them as reported without a rule [prereg: Reported without a decision rule]. The tables give per-K and per-round counts only for the rule outputs [T4] [findings C, item 8].
- Single machine, CPU backend. The run used one host with 24 cores and 23 threads, on the cpu backend [manifest: host, threads, device].
- Deviations. P120: the "labels agree" column of T3 is omitted [Deviations]. P95: the HA strata are not disjoint by definition, and the first stratum an orbit qualifies for holds it; the overlap is 1 orbit [Deviations; T5 header]. P100: a gate orbit without verdict items would make HA not decidable; HA was decidable and g = 0 [Deviations; T5, 4]. P104: D orbits without Double64 verdicts leave the HC set. The tables do not print how many left. By T2 (|D| = 1218) and T8 (set size 1218), none left; this is an inference and not a printed count [findings C, item 5].
- The T5 status column prints "fails" on the T_s row of item 3 although that row meets its condition (see HA) [findings C, item 1].
- Orbits are treated as independent in the bootstrap and in the McNemar test. Orbits of one cell and round are not independent [prereg: Known limits; findings E].
- One system (the standard map) and K values that also occur in the calibration data. Nothing is claimed for other maps or flows [findings E].
- Adaptive allocation is not tested. HB compares a fixed-budget rule with single indicators only [findings E].
- Reference labels are time-bound (horizon N_1 = 1e7) and use the reduced rule with no ladder and no direction test. The BigFloat check shows the ladder would have made rho unresolved on 80 of 100 orbits [T3, last line]. A similar caveat for sticky orbits is raised for exponents that do not converge over very long runs [@HarsoulaKaramanosContopoulos2019; prereg: Known limits].
- The test draw is enriched for sticky and weak-chaos orbits by design, so a(E_test) and e are not population rates [findings E].
- s_I and f share the same Float64 votes, so AUC_I has a mechanical part [findings E].
- Calibration is not blind. The fitted rule performance on p2 is optimistic. The confirmatory evidence is the test [findings E].
- No independent verification report for stage 04 is cited.

## Uncited claims

- The definition of E2 is not cited, because no verified key is at hand.
- The standard map as the system of the study is not cited here, because the findings file gives no key.
- Every threshold of the rules is set by this study. No citation supports the thresholds [prereg: Bibliography note].
- How to word HA, HB and HC together is the owner's call. The preregistration gives no rule that combines them [findings F].
