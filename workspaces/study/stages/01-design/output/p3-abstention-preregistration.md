# Preregistration: p3-abstention

FROZEN 2026-10-30. Scope and value types were agreed by the owner at checkpoint 1 (2026-10-30); the owner delegated the checkpoint 2 decisions to the coordinator, who set them as listed in "Decisions at the checkpoint (2026-10-30)" after an independent review; the owner then approved the freeze. Later changes are amendments, logged before any test data exist.

## Header

| Field | Value |
|-------|-------|
| Study slug | p3-abstention |
| Paper | P3a (to be numbered by the owner). The paper sequence of `/_config/research-question.md` is not edited here |
| Claims tested | HC tests a consequence of claim C3 (diagnostic inequivalence): that indicator disagreement carries information about label failure. HA and HB test hypotheses of research idea 23 (`/_config/research-ideas.md`, row 23; research lines R1 and R2 of that row) restricted to maps, with preregistered rules. Idea-23 H1 (refinement and indicator disagreement predict label failure better than energy error at matched cost) becomes HC in a changed form: there is no energy error on a map, so HC compares the two disagreement sources with each other. Idea-23 H3 (adaptive allocation with abstention beats uniform refinement) becomes HB in a changed form: a fixed-budget abstaining rule against single indicators; adaptive allocation is NOT tested (Known limits). HA is the accounting of error and abstention that both need |
| Value types | FALSIFIES (HC as a test of a consequence of C3; HA and HB as tests of the idea-23 hypotheses; the decision rules are fixed before any test data exist), NEGATIVE allowed (a hypothesis that is not supported is reported with the same figures and tables). NOT claimed: NOVEL (no prior-art check for this study; `research-ideas.md` says shadowing-based reliability work narrows any novelty claim for idea 23), REPRODUCES, TOOL |
| Git SHA at freeze | 6d8dbe8d600d4338d313d90ed2ea7e56b3e0e6c9 (parent commit of the freeze commit) |
| SHA-256 of this file (current version) | e4809007aea5ad2302b804f2ba34ebcc838c2d09e6af22559144fa81b9072c77 (rule as p2: the whole cell is replaced by the word SELF before hashing, UTF-8, LF) |
| Amendment log | A2 (2026-11-04, section "Amendment A2"; parent commit eae6babd565740ebef0c34e249798c3bc13c8253; replaces hash de9490168d609f6e56ebc12edb2b088a115477b8481c6dd633a45d5723fb303a); A1 (2026-11-04, section "Amendment A1"; parent commit 02994688354738e11405605e8a0074481a7ccb80; replaces hash f614a38449f611a42c56a0aa79bbc2d106dfaab679d1972a2b5349011df1a072); A3 (2026-11-04, the seal of section "Amendment A2", item 4: panel_prereg_sha256_accepted lists the hash of the A1 version; parent commit d56fc34c8aa31e94efc9a70ba6ce1efc9a365b7b; replaces hash b1e1059f9e9553d0d7ba6ff6e0236c4b8c1ac507847b57adc3138823573ad0d5). No amendment follows: the next step is the build of the test panel under the A1 version (build copy p3-abstention-preregistration-build.md in study stage 02 output), the first test data |

Value types, in the words of `value-types.md`:
- FALSIFIES: HA, HB, HC each have a metric and a numeric decision rule below; the numbers are in the register.
- NEGATIVE: the null of each hypothesis is reported with the same standing as a pass.

## Hypotheses

Notation and metric definitions are in "Methods" (the rule R, the cost model, the calibration) and "Metrics". Sets are defined in "Panel". The rule R is used in its frozen form theta* (Amendment A1); until A1 exists nothing here can be evaluated. All numeric thresholds named below (e_max, a_max, delta_a, g_max, c_min, alpha_HB, rho_HB, delta_C, min_set, the floors n_T^min and n_W^min, alpha_CI, the bootstrap settings) are in the register. A set below its minimum size is "not decidable": reported with its count, neither pass nor fail.

### HA (error and abstention of the rule)

Setting: the rule form R_a (indicator evidence only, defined in Methods) at the primary budget b_P, on every test orbit with reference class rho in {C, R}, K > 0, not drift. E_test = that set.

Metrics:
- e = (number of decided orbits whose output contradicts rho) / (number of decided orbits). "Decided" = output chaotic-supported or regular-supported. Contradiction: chaotic-supported with rho = R, or regular-supported with rho = C.
- a(A) = (number of orbits of A with output unresolved) / |A|, for A in {E_test, W2, T_s, S2}.
- g = number of K = 0 gate orbits (80, regular by construction) with output chaotic-supported.

Bounds: every upper bound of HA is a one-sided upper 95 percent bootstrap bound (percentile, resamples within strata); every lower bound of a difference is a one-sided lower 95 percent bootstrap bound. The strata are the disjoint parts of E_test: T_s, W2 minus T_s, S2, the other orbits with rho = C, and the orbits with rho = R.

Decision: HA is SUPPORTED iff all four hold:
1. the upper bound of e is <= e_max (reason: the calibration bound e_cal plus a transfer margin of 0.01; see the register);
2. the upper bound of a(E_test) is <= a_max (reason: the calibration bound a_cal plus a transfer margin of 0.05 to a new draw; see the register);
3. the abstention is larger on the hard sets: a(W2) - a(S2) >= delta_a and a(T_s) - a(S2) >= delta_a, each with the lower bound of the difference > 0;
4. the gate: g <= g_max.
Each of E_test, W2, T_s and S2 must hold at least min_set orbits, else HA is not decidable. Otherwise HA is not supported and the failing item is reported (NEGATIVE).

### HB (false-regular rate at matched cost, with a coverage condition)

Setting: R_a at b_P. Comparators: each single indicator X in I_7 at its own matched budget b_X (Methods, "Cost model and matched budget"), read with the p2 class rules of table R1 as a plain label (chaotic or regular; undecided counts as not false-regular, as an abstention of the rule does). On T_s the comparators are the six indicators of I_6 = I_7 minus LE: LE is circular on T_s, because T_s is built from the L1 series. On W2 the comparators are all seven; LE enters on W2 only, capped at N_c (reason in Methods, "Cost model and matched budget").

Metric: for A in {W2, T_s} (test sets, rho = C by construction), FR_rule(A) = share of orbits of A that R_a labels regular-supported; FR_X(A) = share of orbits of A with c_X(b_X) = regular; cov(A) = decided share of the rule on A (its "coverage price" is 1 - cov(A)). Paired counts per X: n10 = orbits false-regular under the rule only, n01 = under X only.

Decision: HB is SUPPORTED iff all hold:
- for both A in {W2, T_s} and every comparator X of that set (13 comparisons: 6 on T_s, 7 on W2; an intersection-union test, so no multiplicity correction is needed): the one-sided exact McNemar p-value P(Bin(n10 + n01, 1/2) >= n01) <= alpha_HB, and FR_rule(A) <= rho_HB * FR_X(A);
- coverage: cov(T_s) >= c_min and cov(W2) >= c_min (a rule cannot pass HB by abstaining on most orbits; reason in the register).
If FR_X(A) = 0 for some X the conditions fail for that X. A set below its minimum size (the floors n_T^min for T_s and n_W^min for W2) is not decidable. A comparator whose matched budget is capped (the cost of X at its longest horizon is below the rule's cost) makes HB not decidable for that rule, with one exception: LE on W2 is always capped (Methods) and is flagged "cost-unmatched, comparator below the rule's cost" without making HB undecidable; HB is also reported without LE. The calibration excludes rules for which a non-LE comparator is capped (Methods, C4, F1). A cheap comparator (WBA at N_c) may beat the rule; that is a reported result, not a defect (Decision D4).
Otherwise HB is not supported (NEGATIVE).

### HC (indicator disagreement against precision disagreement as predictors of label failure)

Setting: the Double64 subset D of the test draw (Panel), orbits with rho in {C, R}, K > 0, not drift; V* = the voting set of theta*.

Metric: label failure f(o) = 1 iff the plain label of o at b_P is not rho(o), where the plain label is the majority of the decided Float64 votes of V* at b_P (chaotic if more chaotic than regular votes, regular if more regular, a tie counts as a failure). Two scores for f, both at b_P: s_I = vote split 1 - |n_C - n_R| / (n_C + n_R) of V* in Float64 (indicator disagreement; 0 if no decided vote); s_P = share of X in V* whose Float64 and Double64 verdicts are both decided and differ (precision disagreement). Cost: the Float64 votes cost a fraction 1/(1 + r) of the Float64 plus Double64 pair that s_P needs (r = Double64 over Float64 cost, Methods), so s_I is not given more than the cost of s_P. AUC_I and AUC_P = area under the ROC curve of the score for f (ties count one half); Delta_AUC = AUC_I - AUC_P with a paired bootstrap interval.

Decision: HC is SUPPORTED (indicator disagreement predicts label failure better) iff Delta_AUC >= delta_C and the lower end of its one-sided 95 percent bootstrap interval is > 0. It is REVERSED (precision disagreement predicts better) iff Delta_AUC <= -delta_C and the upper end of the interval is < 0. Otherwise it is not supported (NEGATIVE: no demonstrated difference). Not decidable if D holds fewer than min_set failures or fewer than min_set non-failures. If D holds fewer failures than the target n_f of A1, HC is still run and reported as underpowered, with the computed power of A1.

### Reported without a decision rule

Risk-coverage curves (error against decided share) of R_a over the rule family around theta* and of the single indicators at their matched budgets; R_ab on D (rule with the ladder override); tables per K, per source of the calibration and per round; HB recomputed with the costs measured in the test run itself (T17 analogue), HB at the round-down matched budget, with the pilot costs, and without LE; the BigFloat check of the reference labels on D_B; the calibration record (A1) with its p2 comparison.

## Panel

One system: the Chirikov-Taylor standard map (Y1; map and tangent in `/_config/research-question.md`, "Ground truth for variational methods"; lifted coordinates internally, wrapped for the indicator runs as in p2). K regimes of the test: the six K values of P1 (K = 0, 0.5, K_c = 0.971635406 [FrahmShepelyansky2010], pi/2, 2 pi, 10). The calibration data (p2) contain these six K values (panel) and four more in the grid (0.75, 1.25, 2, 3), which the test does not use. The five K > 0 values of the test and K = 0 therefore also occur in the calibration data; the orbits are new, and the rule has no K parameter. Nothing of the p2 grid or of the P1 panel is reused in the test.

### Calibration data (read only)

- The existing p2 store `data/p2-indicators`, its manifest `workspaces/study/stages/03-run/output/p2-indicators-manifest.json` (SHA-256 cdef3cf423e16a0b54e81e96fd460f4cb32d7b788325ffbf4e0581f1db27acce, p2 results, Header) and the analysis tables `workspaces/study/stages/04-analyze/output/p2-indicators-tables.md` (T1 to T17). The P1 store is read for the panel relabelling (C0). Nothing is written under `data/`.
- Relabelling (C0, Decision D1): the panel orbits of p2 carry full P1 labels (ladder and direction tests on subsets), the grid orbits carry the reduced rule. Before calibration every Cal panel orbit is relabelled with the reduced rule of the test (next subsection), read only from the P1 and p2 stores. Calibration uses only the reduced labels, so calibration and test share one label definition. The relabelling changes the panel orbits whose P1 label was unresolved by a ladder or direction test (and any other orbit where the reduced rule differs).
- Calibration orbits Cal = p2 E_all after relabelling (panel and grid pooled; controls, the drift orbit and rho = U removed). The sets W2, T_s, S2, regular and the gate set (80 orbits at K = 0) are rebuilt with the definitions below and lambda_ref from the p2 wrapped L1 items (as p2). For orientation, p2 had W2 = 214, S2 = 448, T_s = 30, regular = 1159 orbits (p2 results, Header, "Sets as printed"); the program prints the exact post-relabelling sizes. An orbit whose store lacks the data for relabelling leaves Cal and is counted.
- The p2 data have been seen. Calibration is therefore exploratory and not blind; this is a stated limit, and the confirmatory evidence is the test.

### Test draw (fresh, held out)

The test draw is built by stage 02 with the P1 selection rule (`p1-le-benchmark-preregistration.md`, "Panel", "Selection rule", classes T, R, C, S, X) with these changes only:

1. New seed test_seed (register). The stream rule of P1 is used with the round number r appended to its text.
2. Fresh screening: 512 x 512 starts per K > 0 at points drawn uniformly from the torus with the test seed (not the P1 256 x 256 lattice), run with L5, Float64, N_scr = 1e5 steps. Reason: the P1 lattice would give the same qualifying cells and the draws would overlap the P1 panel (its sticky pool at K = 0.5 had only 25 qualifying cells, all drawn; panel manifest). The screening points are discarded after assignment; no screening value enters any metric. This needs a builder change ("Before the run").
3. Draw in rounds r = 1, ..., R_max, without replacement from the pools. Per round and K > 0: classes T, R, X: 40 orbits per cell; class C: 100 at K = 0.5 and K_c, 40 at the other K; class S: at most 250 per cell at K_c and pi/2, at most 25 per round at K = 0.5 and 2 pi (their pools are small: about 100 and about 20 over all rounds, scaled from the P1 pool counts by the 4 times larger screening); a short cell is recorded "short" as in P1, never filled from another class. K = 0 (T and R, 40 each) in round 1 only; these 80 orbits are the gate set. A round has at most 600 (T, R, X) + 320 (C) + 550 (S) = 1470 orbits; round 1 adds the 80 gate orbits.
4. No controls (the exact periodic points and Y5 validated the code in P1 and p2).
5. Designated Double64 subset D, fixed when the draw is made from design classes only (never from labels or rule outputs): per cell, a seeded random share of class S and of class C at K = 0.5 and K_c (the larger D fraction of the register), and the smaller D fraction of every other cell (reason: the failures that HC needs live in the sticky and weak-chaos cells). At most about 261 orbits per round (165 + 60 + 36), 1309 over 5 rounds. BigFloat check subset D_B: orbits drawn at random from D, stratified by K (size in the register).
6. No orbit is selected by an indicator under test; the design class is the a-priori class, as in P1.

### Reference class rho (test and calibration)

The reduced P1 label rule, applied identically to every test orbit and, by relabelling, to every Cal orbit (the rule of the p2 grid, `p2-indicators-preregistration.md`, "Reference class rho(o)", grid case): instability axis from L5 at N_1 = 1e7 and structure axis from the rotation numbers at N_1/2 and N_1, both under the P1 "Label rules" (values not restated here); unresolved if |lambda_L1 - lambda_L5| exceeds the P1 H1a bound at N_1 or a diagnostic fails; rho = C iff instability "positive growth supported" and structure "irregular", rho = R iff "none detected on the tested horizon" and "quasiperiodic", rho = U otherwise (U orbits leave every set). lambda_ref = lambda_1(N_1). The P1 ladder test and direction test are NOT part of rho, for any orbit, so that D is not treated differently from the rest. This is Decision D1: one label definition for calibration and test, by relabelling. The BigFloat check on D_B reports how often the ladder (Float64, Double64, BigFloat L5 at 1e6, P1 rule) would have changed rho; it changes no label.

### Sets (this study's own definitions; the parameters are in the register, marked as p2 mirrors)

All sets are subsets of E_test (or Cal) with rho = C.
- S2 (strong chaos): lambda_ref >= lam_S2.
- T_s (sticky): lambda_1(n_s) <= f_s * lambda_ref.
- W2 (weak chaos): lambda_ref <= lam_W2, or in T_s. T_s is a subset of W2.
- Gate set: the 80 K = 0 orbits of round 1.
- Drift orbits (lifted |p| above the P1 drift bound at any recorded time) are removed, as in P1.

### Test draw size and power calculation (procedure PP, frozen in Amendment A1)

Decision D5: no power input is assumed. PP uses calibration data only, never test data. The rule rates, the comparator rates, the coverage and the AUC of the indicator-disagreement score are outputs of the calibration; the targets n_T (T_s), n_W (W2) and n_f (HC failures) are computed by the deterministic procedure PP after theta* is chosen. PP, its inputs, its outputs and the program hash are frozen in Amendment A1 together with theta*, before any test data exist. Sources, per source: PP reports its rates for the panel and the grid parts of Cal separately.

Power targets (register): joint power of the full HB decision (both sets, all comparisons, the coverage condition) >= p_HB; power of the full HC decision >= p_HC.

PP for HB:
1. Optimism correction. B_opt resamples of the Cal orbits (with replacement, seed test_seed). On each, the calibration C3 to C6 is repeated and gives theta_b. For each set A in {T_s, W2}: optimism(A) = mean over resamples of [FR of theta_b on the full Cal set A minus FR of theta_b on the resample's set A]. The corrected rule rate is r_corr(A) = r(theta*, A) + max(optimism(A), 0), r from the Cal set. The rule rate used is p_rule(A) = r_corr(A) if r_corr(A) >= the floor of the register, else max(2 r_corr(A), floor).
2. Strata for the nesting T_s inside W2: T_s and the rest R_W = W2 minus T_s. Rates for R_W come from the Cal counts by difference, with the same floor. For each stratum S, from Cal: comparator rates p_X(S) = FR_X(S) at the matched budget b_X(theta*); shared share s_X(S) = (rule false-regular orbits that are also false-regular under X) / (rule false-regular orbits), s_X = 0.5 if the rule has none; decided share d(S) of theta*.
3. Model of one orbit of stratum S: the rule is false-regular with probability p_rule(S); given that, comparator X is false-regular with probability s_X(S), and otherwise with probability max(p_X(S) - s_X(S) p_rule(S), 0) / (1 - p_rule(S)); comparators are independent given the rule outcome (a stated approximation; positive dependence between comparators would raise the power). The orbit is decided if the rule is false-regular, otherwise with probability (d(S) - p_rule(S)) / (1 - p_rule(S)).
4. Draw path: rounds rho_r in {1, 1.25, ..., R_max} with n_T(rho_r) = round(rho_r y_T) and W2 minus T_s = round(rho_r (y_W - y_T)), where y_T and y_W are the per-round yields below. Only paths with n_T >= n_T^min and n_W >= n_W^min are allowed.
5. For each allowed rho_r, simulate the full HB decision 20000 times (fixed seed test_seed): all 13 comparisons with the exact one-sided McNemar test and the margin, and the coverage condition on both sets. Targets n_T = ceil5(y_T rho*) and n_W = ceil5(y_W rho*), where rho* is the smallest rho_r with joint power >= p_HB. If no rho_r up to R_max reaches p_HB, the targets are the cap, the deliverable yield over R_max rounds (y_T R_max and y_W R_max), A1 records the expected power at the cap, and the study runs anyway (NEGATIVE allowed).
6. Type-I check: the same simulation with the rule rate set equal to each comparator rate; A1 records the rejection rate of the joint decision (expected well below 0.05).

PP for HC:
1. AUC_I_cal = the AUC of s_I for the label failure f, both at b_P, on all Cal orbits with rho in {C, R} and K > 0 (relabelled labels; the design strata of D do not exist in Cal), with the voting set V*. Shrinkage: AUC_I_used = 0.5 + 0.5 (AUC_I_cal - 0.5). AUC_P = the register value (p2 H6: the Float64 and Double64 verdicts differed for 1 of 51 orbits, E1, and for none in the other six indicators, p2 results H6, T8; p2 has too few Double64 orbits to estimate it). True difference in the simulation: AUC_I_used - AUC_P.
2. For n_f in {20, 25, ...}: simulate 20000 binormal equal-variance score pairs for n_f failures and the register's number of non-failures per failure, independent predictors, the full HC decision (the point condition Delta_AUC_hat >= delta_C and the lower bound Delta_AUC_hat - 1.645 SE > 0, SE from the Hanley-McNeil variances of the two AUCs; an analytic stand-in for the bootstrap). Target n_f = the smallest n_f with power >= p_HC.
3. The failures D can deliver: |D| of the expected size times 1 / (1 + non-failures per failure). If the target exceeds it, A1 records the power at the deliverable number, and HC is run and reported as underpowered. No extra round is added for HC (it would need indicator runs before the rounds close).

Yield assumptions (per round; from the p2 panel, p2 tables T1 and the P1 class counts; assumption that larger cells yield the same share): class S gives T_s with shares 7/25 at K = 0.5, 4/40 at K_c and 4/40 at pi/2 (p2 panel T_s per class-S orbit), so about 7 + 25 + 25 + 0 = 57 T_s per round (y_T = 57). Weak-chaos orbits per orbit of classes C, S, X: 93/105 at K = 0.5 (panel K index 2: 93 chaotic, all 93 in W2), 27/120 at K_c (K index 3: 107 chaotic, 27 in W2), 4/120 at pi/2, 1/105 at 2 pi, so about 146 + 88 + 11 + 1 = 247 W2 per round (y_W = 247, T_s included). R_max is in the register; at the worst case R_max rounds deliver at most 285 T_s and 1235 W2.

Other hypotheses at the same size: HA needs about 5000 decided orbits at the expected size; with a true error of 0.02 the upper 95 percent bound is about 0.023 (normal approximation, 1.645 x sqrt(0.02 x 0.98 / 5000) = 0.003), inside e_max.

#### Illustration only (drafting values, NOT inputs; A1 replaces them)

At drafting I ran PP-like simulations with assumed numbers (single comparator, no coverage condition, no optimism correction), to size the study roughly. Assumed (no calibration existed): comparator = WBA at the cap N_c, with rates 0.067 on T_s (2 of 30; p2 results H4, tables T6) and 0.019 on W2 (4 of 214; p2 tables F3, k = 9); rule rates 0.01 on T_s (p10 = p11 = 0.005) and 0.003 on W2 (p10 = 0.001, p11 = 0.002).

| Case | n for power 0.80 | n for power 0.90 | Power at that n |
|------|------------------|------------------|-----------------|
| T_s, as above | 140 | 175 | 0.90 at n = 175 |
| W2, as above | 475 | 590 | 0.90 at n = 590 |
| T_s, rule rate 0.02 | | | 0.72 at n = 175 |
| W2, rule rate 0.006 | | | 0.71 at n = 590 |
| T_s, comparator 0.10 | 85 | | 0.99 at n = 175 |
| W2, comparator 0.047 (WBA at b_P, p2 tables F3, k = 6) | 175 | | 1.00 at n = 590 |

The same code returned a rejection rate of 0.011 for T_s when the rule equals the comparator. For HC, with assumed AUC_I = 0.75 and AUC_P = 0.55, 30 failures gave power 0.80 (20 if AUC_I = 0.80). These per-set powers do not give a joint power; the real procedure is PP above. With the corrected yields (57 T_s per round) 175 T_s would need about 3 to 4 rounds.

The number of rounds follows "Stopping rules, Rounds".

## Methods

Every method of the scope-check table is `verified`. Definitions, tolerances and algorithm identifiers live in the catalogs and the module specs; they are not copied. The runs use the study runner `Chaoslab.Experiments` as in p2 (spec `workspaces/module-build/stages/01-spec/output/experiments-spec.md`, amendments D86 to D88). BigFloat is refused by the runner for indicator kinds; BigFloat is used only for L5 on D_B.

### Scope check against the catalogs

Read on 2026-10-30 from `/_config/catalogs/indicators.md`, `lyapunov.md` and `systems.md` (Status column only).

| ID | Role here | Catalog Status | In scope |
|----|-----------|----------------|----------|
| I1 / L1 Benettin | indicator LE: comparator on W2 only (not a voter, not a comparator on T_s); the L1 half of the label disagreement test | verified | yes |
| L5 Householder | reference label (instability axis) at N_1 = 1e7; BigFloat check subset D_B | verified | yes |
| I2 SALI (GALI_2 = exp(L), p2 D70) | indicator SALI | verified | yes |
| I5 FLI | indicator FLI | verified | yes |
| I8 MEGNO | indicator MEGNO | verified | yes |
| I12 NAFF, path E1 | indicator E1 | verified | yes |
| I12b weighted phase advance, path E2 | indicator E2 | verified | yes |
| I13 weighted Birkhoff averages | indicator WBA | verified | yes |
| Y1 standard map | the only system | verified | yes |
| E precision (Float64, Double64; BigFloat for L5 only), F diagnostics | precision ladder and checks | verified (sections E and F of `lyapunov.md`) | yes |
| I9 RLI | not run: p2 has RLI only on 92 of 883 panel orbits and on the grid, so a rule that needs it on every orbit could not be calibrated on p2 | verified | EXCLUDED (choice) |
| I10 / L13 stretching numbers | not run: no threshold rule exists (p2 table R1) | verified | EXCLUDED (choice) |
| I3 GALI_k (k > 2), I4 | not run | planned | EXCLUDED |
| Y2 to Y4, F1 to F7 flows | not run (owner decision: standard map only) | planned | EXCLUDED |

Indicator sets: I_7 = {LE, SALI, FLI, MEGNO, E1, E2, WBA} (comparators); I_6 = I_7 minus LE (the voters of the rule). LE is removed from the voters because T_s is defined from the L1 series: an LE vote on T_s repeats the definition of the set.

### Indicator verdicts (unchanged from p2)

The verdict c_X(b) in {chaotic, regular, undecided} of an indicator X at a budget b is exactly the p2 rule: budget list B = {200 x 2^k, k = 0 to 9}, read points n(b), growth-ratio rules and the numbers of table R1, as in `p2-indicators-preregistration.md`, Metrics M1 and M2, with its Amendment A2 (text only). No number of table R1 is changed or recalibrated here. The indicator specs (tg for SALI, FLI, MEGNO; fq1 for E1; pa1 for E2; bk1 for WBA; L1 for LE) and their parameters are those of p2, "Indicator specs" and "Choices made explicitly" (items 1 to 8): wrapped coordinates, tau = 1, primary direction pi/4, Float64. Not run: pair (RLI), stretching, the sensitivity specs, the direction arm. The run block is derived from the p2 run block by stage 02/03; the arms and horizons are listed below.

### The rule R

Evidence (a), indicator disagreement and detection time. For an orbit o, a budget b, a voting set V (a subset of I_6) and a stability depth w in the range of the register:
- the verdict of X in V is stable at b if c_X(b) equals c_X(b/2^j) for j = 1 to w and all of these are decided (w = 0: stable if decided). This is the detection-time evidence: a verdict that appeared only recently does not count. It uses only budgets up to b, so it costs nothing extra;
- n_C = number of X in V with c_X(b) = chaotic and stable; n_R = number with c_X(b) = regular and stable; n_Cany = number of X in V with c_X(b) = chaotic, stable or not.

Output (rule form R_a, no precision evidence):
- chaotic-supported iff n_C >= q_C;
- else regular-supported on the tested horizon (b steps) iff n_Cany <= v and n_R >= q_R;
- else unresolved.
Parameters theta = (V, q_C, q_R, w, v): V a subset of I_6 with |V| >= 3, q_C in 1..|V|, q_R in (floor(|V|/2) + 1)..|V|, w in {0, 1, 2}, v in {0, 1}. The search space has 2088 elements (sum over |V| = 3 to 6 of C(6, |V|) x |V| x (|V| - floor(|V|/2)) x 3 x 2 = 720 + 720 + 540 + 108).

Evidence (b), the precision ladder (rule form R_ab, defined only on D, where Double64 runs exist): if for some X in V the Float64 verdict at b and the Double64 verdict at b are both decided and differ, the output is unresolved whatever R_a says. This part is fixed, not tuned. Only R_a is used for HA and HB (Decision D2: p2 holds Double64 on only 51 orbits, S_L*, too few to calibrate R_ab, and the precision evidence exists on D only in the test); R_ab is reported on D.

The rule is a function of stored verdicts only; it reads no orbit data, so the same code evaluates it on the p2 store and on the test store.

### Cost model and matched budget

Cost of a spec in Float64, thread-nanoseconds per orbit-step. The values are frozen here so that the matching does not depend on the test run (Decision D3); the test run's own measurement (p2 T17 analogue) is a reported sensitivity.

| Spec (indicators) | c_s | Source |
|-------------------|-----|--------|
| L1 (LE; comparator only) | 207.6 | p2 preregistration, Budget, Throughput inputs (P1 line) |
| tg (SALI, FLI, MEGNO; one run gives all three) | 2658.0 | p2 results, Cost, T17 row tg |
| pa1 (E2) | 1810.5 | T17 row pa1 |
| bk1 (WBA) | 3000.6 | T17 row bk1 |
| fq1 (E1) | 16850.4 | T17 row fq1 |

The p2 pilot measured much lower values for several specs (tg 300, pa1 300, bk1 504, fq1 12609, L1 168; `workspaces/study/stages/03-run/output/p2-indicators-pilot-evidence.md`, section 1), by a factor up to 9. The T17 values (full run) are used because the brief asks for cost "measured as in p2 T17"; the pilot values are a reported sensitivity.

Cost of the rule: C_R(theta, b) = b x (sum of c_s over the distinct specs needed by V; tg counted once; L1 is never needed). Cost is taken as linear in b (the windows of E1, E2, WBA sum to about 2b samples; a stated approximation).

Matched budget of a single indicator X in I_6: b_X = the smallest element of {200 x 2^k, k = 0 to 9} with b x c_{s(X)} >= C_R(theta, b_P) (the ceiling rule: the comparator gets at least the rule's cost, which is conservative for the claim that the rule wins). If no element reaches C_R, X is "capped" (HB not decidable; the calibration excludes the rule). Sensitivity: the largest element with cost <= C_R (round-down).

LE as a comparator on W2: read at min(matched budget, N_c), that is at the same lattice as the other indicators and never beyond N_c. Reason: the L1 series beyond N_c approaches lambda_ref and rho themselves (both use L1 and L5 at N_1 on the same orbit and code path), so LE at larger budgets would be compared with its own reference. With c = 207.6, LE at N_c costs less than the rule for every theta, so LE is always flagged cost-unmatched (the comparator gets less than the rule's cost, which favours the rule); HB is also reported without LE.

Double64 costs (reported, for the cost of HC): c^D64_s = r_s c_s with r from the pilot (Double64 over Float64 measured thread-ns, pilot evidence section 1): tg 10.5, fq1 1.38, pa1 7.28, bk1 6.56. The precision pair for V* costs b_P x (sum over the distinct specs of V* of c_s + c^D64_s), and the Float64 votes alone cost b_P x (sum of c_s), a fraction 1 / (1 + r) of it.

### Calibration procedure (deterministic; its output is frozen by Amendment A1)

C0. Relabelling (Decision D1): before anything else, every Cal panel orbit is relabelled with the reduced P1 rule (instability axis from L5 at N_1 = 1e7, structure axis from the rotation numbers at N_1/2 and N_1, L1 against L5 test, diagnostics), read only from the P1 and p2 stores. The p2 grid orbits already carry this rule. From here on only the reduced labels are used. The sets W2, T_s, S2, regular, gate are rebuilt.

C1. Inputs, read only: the p2 store and tables above; the verdicts c_X(b_k) for X in I_7, k = 0 to 9, computed by the Analysis module with the table R1 rules; the relabelled rho and sets. Fingerprint recorded in A1: SHA-256 of the p2 manifest, of the P1 store manifest, of the calibration program and of its output.

C2. Search space Theta: the 2088 theta of "The rule R". Budget fixed at b_P (w = 2 needs b_P/2 and b_P/4, in B).

C3. For each theta: the output of R_a on every Cal orbit; the error e_cal(theta) among decided Cal orbits; the unresolved shares on Cal and on S2; FR(T_s; theta) and FR(W2; theta) (share of set orbits labelled regular-supported); the decided shares on Cal T_s and Cal W2; the gate count on the 80 K = 0 orbits; C_R(theta, b_P).

C4. Feasible set F (all must hold; thresholds in the register):
- F1 cost cap: C_R(theta, b_P) <= min over X in I_6 of c_{s(X)} x N_c (the minimum is E2, 1810.5 x 102400 = 1.854e8 ns; no I_6 comparator is capped). By arithmetic this removes E1 from V (fq1 alone costs 16850.4 x b_P = 2.16e8 ns). Reason: HB is a matched-cost claim.
- F2 error: e_cal(theta) <= e_cal.
- F3 no abstention on strong chaos: unresolved share on S2 <= s_cal.
- F4 gate: no K = 0 gate orbit is chaotic-supported.
- F5 coverage overall: unresolved share on Cal <= a_cal.
- F6 coverage on the hard sets: decided share on Cal W2 >= c_min and on Cal T_s >= c_min (the coverage condition of HB).

C5. Objective: J(theta) = (FR(T_s; theta) + FR(W2; theta)) / 2 at b_P. Minimise over F.

C6. Selection (a tie rule against over-fitting a 30-orbit set): J_min = the minimum; Theta_tol = {theta in F : J(theta) <= J_min + delta_tol}. Choose theta* in Theta_tol by the first of: lowest C_R; smallest |V|; smallest w; smallest v; smallest bitmask of V (bit 0 SALI, 1 FLI, 2 MEGNO, 3 E1, 4 E2, 5 WBA); smallest q_C; smallest q_R. The order is total, so theta* is unique.

C7. Reported in A1, not used for selection: the rule's FR on W2 and T_s against every comparator at its matched budget on Cal (the HB analogue), the top 20 elements of Theta_tol, the risk-coverage points, the calibration error and abstention of theta*, and J and all PP rates by source (panel, grid). These are fitted values and optimistic.

C8. If F is empty: A1 records "no feasible rule" and the study stops at calibration, reported as a NEGATIVE result (no test data are produced). Only the owner can change F or Theta; such a change is an amendment before any test data exist and is labelled exploratory (the p2 data have been seen).

C9. Determinism check: the program is run twice and the two outputs must be byte-identical. The program does not exist yet ("Before the run"); it is written, verified and hashed before A1.

C10. After theta* is chosen, procedure PP (Panel, "Test draw size and power calculation") is run and its output goes into A1.

### Precision

Float64 in all main arms. Double64 on D for tg, fq1, pa1, bk1 (horizon 102400, 102401 for pa1). BigFloat(256) for L5 at N_B = 1e6 on D_B only, with Double64 L5 at 1e6 on D_B for the same check. Defaults of `/_config/numerics-defaults.toml` used unchanged: `renorm_interval_tau` (= tau = 1), `bigfloat_bits`, and the three health multipliers. The `seed` default is not used for the test (test_seed is). Overrides: `precision` per arm as stated; `backend` = `auto` (Float64) and `cpu` (Double64), as p2.

### Arms of the test run (draft)

| Arm | Orbits | Items | Horizon |
|-----|--------|-------|---------|
| ref (labels) | all test orbits | L1, L5, Float64 | 1e7 |
| ind | all test orbits | tg, fq1, pa1, bk1, Float64 (LE is read from the ref L1 item, up to N_c) | 102400 (102401 for pa1) |
| d64 | D | tg, fq1, pa1, bk1, Double64 | 102400 (102401 for pa1) |
| chk | D_B | L5, Double64 and BigFloat | 1e6 |

The Float64 L5 value for the check is read from the ref arm at the recorded step 1e6.

## Metrics

Exact definitions; scores are read from the stored series and tables, nothing is recomputed from orbits.

- Verdicts, budgets B, b_P, N_c, N_1 and n_s: as p2 (M1, M2) and the register.
- Rule output o_R(o; theta) in {chaotic-supported, regular-supported, unresolved}, R_a or R_ab (Methods).
- e, a(A), g (HA); FR_rule, FR_X, cov, n10, n01, exact McNemar p (HB); f, s_I, s_P, AUC, Delta_AUC (HC): defined in Hypotheses.
- AUC: the Mann-Whitney statistic, ties one half. Bootstrap: orbits resampled with replacement (HA within the strata of HA; for the paired Delta_AUC the same resample for both scores), the number of resamples and the seed (test_seed) in the register, percentile intervals.
- Cost: thread-ns per orbit-step per spec as the runner reports them (p2 M9, T17 analogue), used only as a sensitivity.
- Risk-coverage: error e against decided share at b_P, for the rule family around theta* (q_C and q_R varied by one step) and for the single indicators at their matched budgets; reported.

## Budget

Environment limits (`/_config/environment.md`): this computer, at most 23 cores, GPU RTX 5060 Laptop (8 GB), Double64 on the CPU, no wall-clock limit stated; the campaign cap of the register applies. Per-spec lines as in "Cost model" (Float64 thread-ns per orbit-step). Safety factors: the register. Thread-hours.

Orbit counts from the per-cell sizes of "Test draw": at most 1470 orbits per round (600 + 320 + 550), plus 80 gate orbits in round 1; D at most about 261 per round. Worst case over R_max = 5 rounds: N = 5 x 1470 + 80 = 7430 orbits, |D| = 1309. Expected: 4 rounds (about 3 to 4 from the yield assumptions), N = 5960, |D| = 1048.

| Item | Worst case orbits | Line | Before safety | With safety, worst case (expected) |
|------|-------------------|------|---------------|------------------------------------|
| Calibration (reading the p2 and P1 stores, 2088 theta, 200 optimism resamples, about 3000 orbits) | 0 | seconds to minutes | about 0 | about 0 |
| Screening, 5 K x 512 x 512 starts, L5, 1e5 steps | 1.31e6 starts | 207.6 ns | 7.6 | 15.1 (15.1) |
| ref: L1 and L5 at 1e7, 4.15 s per orbit | 7430 | 207.6 ns | 8.6 | 17.1 (13.7) |
| ind: tg + fq1 + pa1 + bk1 = 24,319 ns per step, 2.49 s per orbit | 7430 | T17 measured | 5.1 | 10.3 (8.2) |
| d64 on D: tg 6000, fq1 28100, pa1 4200, bk1 5000 (p2 Double64 lines), 4.43 s per orbit | 1309 | p2 lines | 1.6 | 4.8 (3.9) |
| chk: Double64 L5 at 1e6 (6634.6 ns) and BigFloat L5 at 1e6 (18200 ns, P1 line) on 100 orbits | 100 | P1 lines | 0.7 | 2.1 (2.1) |
| Total | | | about 24 | about 49.4 (about 42.9) |

Arithmetic: ref per orbit = 2e7 x 207.6 ns = 4.152 s; 7430 x 4.152 s x 2 = 61,700 s = 17.1 h. ind per orbit = 24,319 ns x 102400 = 2.49 s; 7430 x 2.49 s x 2 = 37,000 s = 10.3 h. d64 per orbit = (6000 + 28100 + 4200 + 5000) ns x 102400 = 4.43 s; 1309 x 4.43 s x 3 = 17,400 s = 4.8 h. Screening = 5 x 262144 x 1e5 x 207.6 ns = 27,200 s = 7.6 h, x 2 = 15.1.

Sensitivity with the indicator arms charged at the p2 T17 line values (tg 18000, fq1 84300, pa1 12600, bk1 15000 ns per step, 129,900 in all, 13.3 s per orbit) and the safety factor 2: ind = 7430 x 13.3 s x 2 = 197,600 s = 54.9 thread-hours (expected 44.0); the total with safety is about 94.0 in the worst case (about 78.7 expected).

Wall time on 23 threads at 50 percent parallel efficiency (11.5 effective threads; the p2 pessimistic bound; P1 took 3.3 h wall and p2 0.29 h): worst case 49.4 / 11.5 = 4.3 h, or 94.0 / 11.5 = 8.2 h with the T17 lines; expected 3.7 h, or 6.8 h with the T17 lines. The worst case stays below the campaign cap of the register, so R_max is not lowered. The p2 pilot suggests the real cost is several times lower for most specs. The GPU is not needed: in the p2 K19 record (`p2-indicators-k19-records.toml`) the device was slower than the CPU for every device-eligible spec, and the test run records its own K19 values before the first run, as p2 did.

Memory and storage: the E1 buffer is 1.65 MB per orbit in Float64 (p2 spec E28), so a round fits one launch at 16 GB free; the store is of the order of the p2 store plus the L1 and L5 series of at most 7430 orbits, far below any limit.

Audit line "Budget": at most about 49.4 thread-hours with safety in the worst case (94.0 with the T17 lines), 4.3 to 8.2 wall hours, core limit 23 respected.

## Analysis plan

Every figure and table and the hypothesis it decides. All use stored series and tables.

| # | Item | Content | Decides |
|---|------|---------|---------|
| T1 | Calibration record | A1 content: fingerprints, relabelling counts, F and its size, J_min, Theta_tol size, theta*, calibration error and abstention, HB analogue on Cal, PP output, by source | none (reported) |
| F1 | Calibration landscape | FR(W2) and FR(T_s) against C_R over F, theta* marked | none |
| T2 | Test draw | rounds, pools, short cells, counts of rho (C, R, U), sizes of E_test, W2, T_s, S2, D, D_B, achieved against the targets n_T and n_W of A1 | none |
| T3 | Reference check | Double64 and BigFloat L5 on D_B: flips of rho under the P1 ladder | none (reported) |
| T4 | Rule outputs | counts of the three outputs per set, per K, per round | context |
| T5 | HA table | e with bound, a per set, differences with bounds, g, set sizes | HA |
| T6 | HB table | per set and comparator: b_X, cost, FR_X, FR_rule, n10, n01, p, margin test, capped flag, and beside each comparison the coverage price (decided share cov(A) of the rule) | HB |
| T7 | HB sensitivities | round-down budget; pilot costs; test-measured costs for matching; without LE | reported |
| F2 | Risk-coverage | e against decided share, rule family and single indicators at their matched budgets | reported |
| T8 | HC table | failures, non-failures, AUC_I, AUC_P, Delta_AUC with interval, power of A1, underpowered flag | HC |
| F3 | ROC curves | s_I and s_P for f | reported |
| T9 | R_ab on D | outputs of R_a against R_ab, number of ladder overrides, error of the decided | reported |
| T10 | Cost | measured thread-ns per orbit-step per spec against the cost model | reported |
| T11 | Flags and health | flag counts, failed chunks | reported |

Rules for the analysis:
- Hypotheses are decided with the numbers of the register and nothing else, in the order HA, HB, HC. No threshold is tuned after a look at test data.
- Any analysis not listed is labelled "exploratory".
- Results are reported by K and by round as well as pooled. Orbits are treated as independent in the bootstrap (within the strata of HA); rounds are not bootstrap units.
- The report cites only keys present in `/_config/bibliography/references.bib` and never an item of `unverified.md`. Keys used in this file: FrahmShepelyansky2010, HarsoulaKaramanosContopoulos2019 (the second in Known limits). Keys the report may add, all checked present: Maffione2011a, SanderMeiss2020, DuignanMeiss2023, Skokos2001, Froeschle1997, CincottaSimo2000, Laskar2003. No threshold of this study comes from the literature: all are set here.

## Stopping rules

Unit of work = chunk (runner): the items of one run type in one cell, resumable.

Chunk failure and handling: as p2, "Stopping rules" (non-finite values, health errors, determinant or symplecticity monitor above the multipliers in more than 5 percent of a chunk, first-chunk cost above 2 times its line, a Double64 item that disagrees with its Float64 item on the orbit coordinates by more than 1e-10 over the first 100 steps for classes T and R and the first 5 steps otherwise). A failed chunk is rerun once on the CPU; a second failure is final and its orbits are "excluded, failed twice". The campaign stops if more than 5 percent of chunks fail twice or the wall cap of the register is reached.

Rounds (the sizing rule): after round r (r = 1 to R_max) only the ref arm (L1, L5 at 1e7) runs for that round's orbits, and only the counts of T_s and W2 are computed from lambda_1 and rho. Stop adding rounds when T_s >= n_T and W2 >= n_W (the targets frozen in Amendment A1), or at r = R_max. Neither the ind arm, the d64 arm nor any rule output is computed or opened until the last round is closed (no peeking: the number of rounds depends only on label-defined set counts). If at R_max the counts are below the targets but at or above the floors n_T^min and n_W^min, HB is decided and flagged "underpowered" with the power recomputed at the achieved sizes; below the floors HB is not decidable for that set.

Calibration: if F is empty the study stops after A1 (Methods, C8).

Early stopping on scientific results: none. No hypothesis is decided from partial data. The only early-stop causes are the failures above, the wall cap and the calibration rule.

## Amendment rule

Frozen after study/01. A change after freeze is an amendment: one numbered entry in the header log with date, what changed, why and the git SHA of the commit that makes it; allowed only before the first test data file exists (the first test screening run counts as test data). Planned and fixed here:
- Amendment A1 (calibration freeze): records the calibration fingerprints, the relabelling counts, F, theta* (V, q_C, q_R, w, v), the cost-model values actually used, the calibration error, abstention and HB analogue, the output of the power procedure PP (optimism corrections, rates by source, shared shares, comparators, coverage, AUC_I_cal and AUC_I_used, targets n_T, n_W, n_f, the expected power at the cap if the cap binds, flags), the calibration program hash and the test seed. Logged before any test screening, draw, label or indicator run exists. After A1 the rule is a fixed function; no parameter changes after test data exist (a later change is not an amendment but an exploratory re-analysis, and the preregistered decision stands).
- Any change of the thresholds in the register after freeze is an amendment of the same kind, allowed only before test data exist and, for the calibration constraints F, labelled exploratory because the p2 data have been seen.

Audit notes (stage 01 contract). "Verified methods": the scope-check table, all `verified`. "Decidable": every hypothesis has a metric and a numeric decision rule; the rule parameters theta* and the targets n_T, n_W, n_f exist only after A1, so "Decidable" holds through PP and C1 to C10, which are deterministic and fully specified here. "Single home": thresholds live in the register; inline text names them. "Budget": the Budget section.

## Known limits

- Reference labels are time-bound: rho is "none detected on the tested horizon" or "positive growth supported" at N_1 = 1e7. A very sticky orbit trapped beyond N_1 gets rho = R and an indicator that calls it regular is counted correct (HarsoulaKaramanosContopoulos2019 report an exponent not converged after 5e7 steps). The reduced label rule has no ladder and no direction test for any orbit, so the labels are weaker than the P1 panel labels; the BigFloat check on D_B measures how much. Calibration labels are made by relabelling with the same rule, so the two share one definition, but the relabelled panel orbits are not the P1 labels that the p2 results used.
- The precision ladder for indicators is Float64 against Double64. BigFloat is refused by the runner for indicator kinds, so the ladder has no top rung for evidence (b); a flip-free ladder does not prove an indicator verdict correct. Double64 runs only on D, so HA and HB use R_a; the full rule R_ab is only an extension on D. With E1 removed from the voters by F1, the precision score of HC is nearly constant: p2 found 1 flip in 51 orbits and it was on E1, so HC may reduce to AUC_I against about 0.5. Moreover s_I and the label failure f are computed from the same Float64 votes, so AUC_I has a mechanical part (a split vote makes a wrong majority more likely); HC measures the usefulness of that split against precision disagreement, not an independent signal.
- One system: the standard map. Nothing is claimed for other maps or for flows. The flow version of idea 23 (timestep refinement) and the adaptive allocation of its H3 are not tested; HB compares a fixed-budget rule with single indicators and cannot say that adaptive allocation helps.
- The calibration is not blind: the rule family, F, J and the rule were designed after the p2 results were read. Its fitted performance on p2 is optimistic and is not evidence; the optimism correction of PP only adjusts the planning, not the test. The confirmatory evidence is the test, which uses new orbits but K values that also occur in the calibration data and the same class rules R1 of p2.
- The numbers of table R1 were set for p2 and not tuned for abstention. The test draw is enriched for sticky and weak-chaos orbits by design; rates such as a(E_test) depend on that mix and are not population rates.
- The matched cost uses p2 measurements (T17) that differ from the p2 pilot by up to a factor 9 for some specs. A cheap comparator such as WBA reaches N_c at the rule's cost, so a rule that looks at several indicators at b_P steps is compared with WBA at 102400 steps; HB can fail for that reason alone, and that is a result, not a defect. LE on W2 is capped at N_c and is cost-unmatched in the rule's favour; HB is also reported without it.
- T_s is a data-defined subset (lambda_1 at n_s); a different definition of stickiness may change HA and HB. The yields used by PP are extrapolated from the p2 panel at larger cell sizes.
- The tests and the bootstrap treat orbits as independent; orbits of one cell and round are not. PP treats comparators as independent given the rule outcome.

## Thresholds register (single home; all set at the checkpoint, 2026-10-30)

"p2 mirror" marks a value copied from p2 so that this study is self-contained; p2's file is its source.

| Name | Value | One-line reason |
|------|-------|-----------------|
| e_cal | 0.02 | below the largest single-indicator false-chaotic rate on the regular set in p2 (MEGNO 31/1159 = 0.027, p2 tables T2 at k = 6) |
| e_max | 0.03 | e_cal plus a transfer margin of 0.01 |
| a_cal | 0.10 | the share of orbits I accept to leave unresolved on the calibration data (judgement; p2's split-vote U is not an abstention rate) |
| a_max | 0.15 | a_cal plus a transfer margin of 0.05 to a new draw |
| s_cal | 0.02 | strong-chaos orbits should almost never be abstained on |
| delta_a | 0.10 | five times s_cal: the rule must abstain visibly more on W2 and T_s than on S2 |
| alpha_CI | 0.05 (one-sided bounds) | the usual level |
| g_max | 1 (of 80 K = 0 orbits) | p2 had 0 of 80 for every indicator; one in 80 keeps the upper bound near 0.06 |
| c_min | 0.5 | a rule that decides fewer than half of the hard-set orbits cannot claim a lower false-regular rate |
| alpha_HB | 0.05 per comparison, intersection-union over 13 | joint level 0.05 without correction |
| rho_HB | 0.5 | a halving: scale-free and not resting on one orbit (comparator rates differ by about 3 to 4 between W2 and T_s, p2 tables T2) |
| n_T^min, n_W^min | 100, 300 | floors; 0.57 and 0.51 of the illustrative drafting targets 175 and 590 |
| min_set | 20 | generic floor for HA sets and HC failures and non-failures; p2 mirror |
| p_HB, p_HC | 0.80, 0.80 | joint power of the full HB decision; power of the full HC decision |
| delta_C | 0.10 | smallest AUC gap worth reporting, above the sampling error at the planned size |
| delta_tol | 0.02 | one T_s orbit moves J by 0.017 |
| B_opt | 200 | optimism bootstrap resamples (cost: seconds) |
| rule-rate floor, factor in PP | 0.01, 2 | margin for the optimism of a fitted rule: the factor 2 only when the corrected rate is below the floor |
| AUC_I shrinkage | 0.5 + 0.5 (AUC_I_cal - 0.5) | halves the fitted AUC advantage |
| AUC_P in PP | 0.55 | p2 H6: 1 flip in 51 orbits, the precision score is nearly constant |
| non-failures per failure in PP | 20 | D is enriched for sticky and weak chaos; D failure rate about 0.05 |
| PP simulation runs | 20000 | enough for power to 0.01 |
| lam_W2, lam_S2, n_s, f_s | 0.1, 1.0, 1e5, 0.2 | p2 mirror (p2 "Panel", W2, S2, T_s), so that calibration and test sets are the same constructs |
| b_P, N_c, N_1, budgets B | 12800, 102400, 1e7, 200 x 2^k (k = 0 to 9) | p2 mirror |
| w range, v range | {0, 1, 2}, {0, 1} | a verdict that held over up to two doublings; at most one stray chaotic vote |
| screening size, R_max, per-cell sizes | 512 x 512 per K, 5, as in "Test draw" | a four times larger pool than P1 at 7.6 thread-hours; the worst case of 5 rounds is 8.2 wall hours at most, below the cap |
| D fractions, size of D_B | 0.30 (S, and C at K = 0.5 and K_c) and 0.05 (other cells), 100 | failures live in the sticky and weak-chaos cells; 100 is the P1 ladder size |
| test_seed | 20261030 | date of drafting; new |
| safety factors, wall cap | 2 (Float64), 3 (Double64, BigFloat), 12 h | p2 mirror |
| bootstrap | 10^4 resamples (HA, HC) | p2 mirror |

## Decisions at the checkpoint (2026-10-30)

Taken by the coordinator on the owner's delegation.

1. D1, reference labels: the reduced P1 label rule (L5 at N_1, rotation-number structure axis, L1 against L5 test, diagnostics; no ladder and no direction test) is applied to every test orbit, plus the BigFloat and Double64 check on D_B. Calibration panel orbits are relabelled with the same rule before calibration (C0). Reason: one label definition for calibration and test, by relabelling.
2. D2, precision stays secondary: HA and HB use R_a; precision enters only HC and R_ab on D. Reason: p2 holds Double64 on only 51 orbits (S_L*), too few to calibrate R_ab; making it primary would add about 27 thread-hours (7430 x 4.43 s x 3 / 3600).
3. D3, cost model: the p2 T17 values are frozen for matching. HB recomputed with the costs measured in the test run itself is a reported item without a decision rule.
4. D4, HB against a cheap comparator (WBA at N_c) is kept as written. A failure of HB is a reported result (NEGATIVE allowed).
5. D5, power: no assumed rates remain. Procedure PP uses calibration data only (rates with a bootstrap optimism correction, coverage, shrunk AUC_I), simulates the full HB and HC decisions and gives n_T, n_W, n_f, frozen in Amendment A1 with theta*. The round cap and the minimums n_T^min, n_W^min are floors or caps. The drafting numbers are an illustration only.
6. D6, thresholds: all proposed values are accepted and are set as listed in the register.

Review fixes applied the same day (independent review, decisions by the coordinator): coverage condition c_min in HB and F6; LE removed from the voters and from the T_s comparisons (search space 2088, LE capped at N_c on W2); Cal panel relabelling and by-source reporting; joint power of the full HB decision in PP; bootstrap optimism correction, shrunk AUC_I and the full HC simulation in PP; s_I and f at the same budget b_P; yields, round arithmetic and budget recomputed (worst case 49.4 thread-hours, 4.3 wall hours; 94.0 and 8.2 with the T17 lines); a_max and e_max justified from calibration bounds; the K claim reworded; value types and audit notes; HA minimum sizes and stratified bootstrap bounds; thresholds named inline and valued only in the register; scope check moved into Methods; C10 after C9.

Notes. The pool sizes of classes C and S at the test K values are unknown until screening; if the targets are not reached in R_max rounds, HB is flagged underpowered or not decidable (Stopping rules). A sixth round or more K values would be an amendment. The paper row name "P3a (to be numbered by the owner)" is used in this text only; `research-question.md` is not edited (in the current sequence P3 is the higher-dimensional paper). TOOL is not claimed.

## Before the run

The study does not start before this is done. Code the study needs and that does not exist yet:
- the relabelling of Cal panel orbits with the reduced rule from the P1 and p2 stores (C0); I did not check that the P1 store holds the rotation-number data and L5 series needed for every panel orbit;
- the rule evaluator (R_a and R_ab, from stored verdicts);
- the calibration program (C0 to C10, including procedure PP and the optimism bootstrap), deterministic, with the byte-identical double run;
- the matched-cost computation (cost model, ceiling and round-down rules, cap flags, LE capped at N_c);
- the HB and HC analysis (exact McNemar, stratified bootstrap bounds, paired bootstrap AUC, the tables and figures of the Analysis plan);
- the panel builder in seeded rounds: fresh screening from the test seed (512 x 512 random starts per K), draw without replacement per round, the stream rule with the round number, the design subsets D and D_B fixed at draw time. I did not check whether the current builder takes a grid size and a point source as parameters.
All of it is built and verified by a module-build amendment before study stage 02 (Pattern 14: the study uses the library as a package). Amendment A1 (calibration freeze) comes after that code is verified and before any test data.

## Checks I ran

- Catalog Status read 2026-10-30: indicators I1, I2, I5, I8, I9, I10, I12, I12b, I13 verified; lyapunov L1, L5, L13 verified; systems Y1, Y5 verified, F1 to F7 planned.
- Bibliography keys grepped in `references.bib` (present) and in `unverified.md` (absent): FrahmShepelyansky2010, HarsoulaKaramanosContopoulos2019, Maffione2011a, SanderMeiss2020, DuignanMeiss2023, Skokos2001, Froeschle1997, CincottaSimo2000, Laskar2003, Darriba2012, Skokos2010, Das2017, SandorMaffione2016, Benettin1980a.
- p2 numbers were read from the p2 preregistration, results, tables and pilot evidence files, not from the HDF5 store (the store was not opened; no Julia was run): set sizes, T_s and W2 false-regular fractions (tables F3, T6), regular-set false-chaotic counts (tables T2, k = 6), per-K and per-class yields (tables T1, P1 class counts), costs (T17), pilot ratios, and the P1 panel pools (manifest cells in `p1-le-benchmark-panel.toml`).
- The search-space size (2088) and the round, orbit-count and budget arithmetic were computed by script. The illustrative power numbers come from Python (numpy, scipy) in the scratch area. None of it is in the repo; the real power calculation is procedure PP, run at calibration, and it has not been run.
- Not checked: that the builder supports fresh screening; that the P1 store allows the relabelling; that the Analysis module can compute the rule; the calibration itself (no calibration was run, so nothing is claimed about the feasibility of F or about what theta* will be); the sizes of the pools at the test K values.
- Nothing was frozen, committed or pushed; this file is the only file written.

## Amendment A2 (2026-11-04)

Approved by the owner at the checkpoint of study amendment P69 of module amendment D94 (2026-11-04), before any test data exist: no test screening, panel, draw, label, indicator run or pilot on a panel orbit exists. Base sections above are left as frozen; where they differ, this section wins. A2 adds the machine-readable panel input block and run block that the library reads (`experiments-spec.md`, "Amendment D94", R22 and R25). It changes no hypothesis, no threshold of the register, no set, no rule, no arm content and no budget line other than the BigFloat line of the chk arm and the screening note, both in item 6. Every value in the two blocks equals the frozen prose (Panel, "Test draw (fresh, held out)"; Methods; Budget; Stopping rules; Thresholds register), the P1 panel block by the rule "P1 selection rule with these changes only", the p2 run block by the rule "the run block is derived from the p2 run block" (Methods, "Indicator verdicts (unchanged from p2)"), or a library default stated in the notes of item 2 and item 3. Section "Amendment A1" (calibration freeze) is still planned and comes after A2.

1. Order of work. A2 (this section); the dry campaign of the round loop (`workspaces/study/stages/03-run/references/round-loop.md`, the gate before the A1 freeze run: a small panel that is not the test panel); the calibration freeze run and Amendment A1; Amendment A3, the seal of item 4; the build of the test panel; the rounds. The build of the test panel is the first test screening and therefore the first test data; no amendment follows it (Amendment rule). The panel is built under the file as A1 leaves it, by the rule of item 4.

2. Panel input block. It is the first TOML block of this file and holds the `[panel]` table, as the library requires. Sources: the P1 panel block (`p1-le-benchmark-preregistration.md`, Amendment A1 item 7, with A2 and A3 as they stand: the key `sn_trace` is absent) for every key of the P1 selection rule; the keys of Amendment D94 from the frozen text of "Test draw (fresh, held out)" and the register. Notes. (a) `n_c` is replaced by the table `quota` (a repeated definition is refused). (b) `quota`: per round and K index; the K = 0 row is drawn in the rounds of `gate_rounds` only. The class S size at K = 10 is 0 (the frozen text lists S sizes for K_c, pi/2, K = 0.5 and 2 pi only; 250 + 250 + 25 + 25 = 550). (c) `d_fraction`: 0.30 for class S and for class C at K = 0.5 and K_c, 0.05 for every other cell including the K = 0 gate cells (2 of 40 per class, 4 in all). With m = floor(f n + 1/2) (exact rationals) the S cells at K = 0.5 and 2 pi (n = 25, f = 0.30) give 8, so a fully populated round has |D| = 262, not the 261 of the frozen text, and the worst case over 5 rounds is 5 x 262 + 4 = 1314; the frozen 261 and 1309 are read as "about" values (the text says "at most about 261"), about 0.4 percent of the d64 cost; the S fraction at K = 10 is inert (size 0). (d) `grid_side` is required by the library and has no role with the random screening source; it keeps the P1 value, and `newton_grid_side` keeps the P1 value so that the periodic set is that of the P1 panel. (e) `seed` is test_seed of the register; `tau` and `bigfloat_bits` are absent and take the defaults of `_config/numerics-defaults.toml` (`renorm_interval_tau` = 1, `bigfloat_bits`), as Methods, "Precision", states. (f) `s_D`, `cap_D`, `s_L`, `cap_L`, `n_ctrl_max`, `mu` are required keys with the P1 values; no arm of this study uses S_D or S_L and there are no controls. (g) K index i counts from 1 in the order of `K`: 1 is K = 0, 2 is K = 0.5, 3 is K_c, 4 is pi/2, 5 is 2 pi, 6 is K = 10.

```toml
[panel]
K = ["0", "0.5", "0.971635406", "pi/2", "2pi", "10"]
grid_side = 256
N_scr = 100000
alpha = "pi/4"
theta_T = 5e-4
g_T = 5.0
r_T = 0.05
theta_C = 0.01
r_C = 0.3
theta_S = 1e-3
w = 1000
f_S = 0.2
L_S = 20000
log10_r = [-4.0, -2.0]
q_max = 6
rho = 1e-12
delta_R = 1e-3
newton_max = 50
delta_dup = 1e-8
n_ctrl_max = 12
mu = "(3 + sqrt(5))/2"
s_D = 7
cap_D = 210
s_L = 4
cap_L = 110
seed = 20261030
screen_source = "random"
screen_points = 262144
r_max = 5
gate_rounds = [1]
db_size = 100
newton_grid_side = 256
controls = "none"

[panel.quota]
"1" = { T = 40, R = 40 }
"2" = { T = 40, R = 40, C = 100, S = 25, X = 40 }
"3" = { T = 40, R = 40, C = 100, S = 250, X = 40 }
"4" = { T = 40, R = 40, C = 40, S = 250, X = 40 }
"5" = { T = 40, R = 40, C = 40, S = 25, X = 40 }
"6" = { T = 40, R = 40, C = 40, S = 0, X = 40 }

[panel.d_fraction]
"1" = { T = 0.05, R = 0.05 }
"2" = { T = 0.05, R = 0.05, C = 0.30, S = 0.30, X = 0.05 }
"3" = { T = 0.05, R = 0.05, C = 0.30, S = 0.30, X = 0.05 }
"4" = { T = 0.05, R = 0.05, C = 0.05, S = 0.30, X = 0.05 }
"5" = { T = 0.05, R = 0.05, C = 0.05, S = 0.30, X = 0.05 }
"6" = { T = 0.05, R = 0.05, C = 0.05, S = 0.30, X = 0.05 }
```

3. Run block. Arms are those of "Arms of the test run (draft)". Notes. (a) The window specs fix the horizon of every arm that uses them (the library refuses an arm that mixes fq1 or bk1, whose run length is 102400, with pa1, whose run length is 102401). The rows "ind" and "d64" of the draft table are therefore two arms each: `ind` and `d64` hold tg, fq1, bk1 at 102400; `ind_pa1` and `d64_pa1` hold pa1 at 102401. No item, orbit set or horizon changes. (b) `ref` is released each round; every other arm is a closure arm (Stopping rules, "Rounds": no peeking). (c) Coordinates are wrapped in every arm, tau = 1, direction pi/4, as the indicator specs of p2 (Methods, "Indicator verdicts"); LE is read from the L1 item of `ref`. (d) There are no controls, so `control_methods`, `control_bound`, `control_scope` and the control classes of `setup_check_steps` are required by the library and inert; they keep the p2 values. (e) The Float64 and Double64 budget lines of the four indicator specs are the lines of "Cost model and matched budget" and Budget; the lines of `ref` and `chk` are the lines of the Budget table, except the BigFloat line (item 6). No `ns_gpu` line is given (Budget: the GPU is not needed); the library then charges the CPU line. (f) The test run records its own K19 values before the first run, at 23 threads. (g) The pilot uses no panel orbit and writes outside the study store (item 5). (h) `panel_prereg_sha256_accepted` is governed by item 4: the entry shown is the placeholder of 64 zeros that Amendment A3 replaces.

```toml
[run]
panel_prereg_sha256_accepted = ["b1e1059f9e9553d0d7ba6ff6e0236c4b8c1ac507847b57adc3138823573ad0d5"]

[run.primary]
tau = 1
eps = 1e-8                    # required by the library; unused (no L2 item in this study)
direction = "pi/4"
coordinate_mode = "wrapped"

[run.grid]
points_per_decade = 20
first_steps = 100
extra_fractions = [0.1, 0.5]

[run.backend]
Float64 = "auto"
Double64 = "cpu"
BigFloat = "cpu"

[run.health]
fail_fraction = 0.05
control_bound = 1e-8          # inert: no controls
control_methods = ["L1"]      # inert: no controls
control_scope = "exact_start" # inert: no controls
setup_check_steps = { T = 100, R = 100, control_elliptic = 100, control_hyperbolic = 5, C = 5, S = 5, X = 5, h4 = 5 }
setup_check_tol = 1e-10
drift_bound = 1e6
rerun_attempts = 1
stop_failed_fraction = 0.05
first_chunk_cost_factor = 2.0
wall_cap_hours = 12

[run.budget.Float64]
ns_cpu_thread = 207.6
safety = 2

[run.budget.Double64]
ns_cpu_thread = 6634.6
safety = 3

[run.budget.BigFloat]
ns_cpu_thread = 341296.4
safety = 3

[run.k19]
threads = 23

[run.verify]
backend_agreement = { per_orbit_sigma = 6.0, mean_standard_errors = 4.0 }
rerun_sample = 20

[run.pilot]
orbits_per_k = 4
horizon = 100000
seed_label = "p3-abstention-pilot"

[[run.indicator]]
id = "tg"
kind = "tangent"
members = ["FLI", "MEGNO", "SALI"]
[run.indicator.budget.Float64]
ns_cpu_thread = 2658.0
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 6000
safety = 3

[[run.indicator]]
id = "fq1"
kind = "frequency"
observable = "exp_i_q"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
order = 1
positions = "half"
pad = 8
mean_removal = true
mirror_ratio = 0.9
floor_multiplier = 32.0
[run.indicator.budget.Float64]
ns_cpu_thread = 16850.4
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 28100
safety = 3

[[run.indicator]]
id = "pa1"
kind = "phase_advance"
observable = "lifted"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
birkhoff_order = 1
[run.indicator.budget.Float64]
ns_cpu_thread = 1810.5
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 4200
safety = 3

[[run.indicator]]
id = "bk1"
kind = "birkhoff"
observables = [
  { kind = "cos", coordinate = 1, frequency = 1 },
  { kind = "sin", coordinate = 1, frequency = 1 },
  { kind = "cos", coordinate = 2, frequency = 1 },
  { kind = "sin", coordinate = 2, frequency = 1 },
]
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
birkhoff_order = 1
[run.indicator.budget.Float64]
ns_cpu_thread = 3000.6
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 5000
safety = 3

[[run.arm]]
name = "ref"
orbits = "P"
methods = ["L1", "L5"]
precisions = ["Float64"]
horizons = [10000000]
release = "each_round"

[[run.arm]]
name = "ind"
orbits = "P"
indicators = ["tg", "fq1", "bk1"]
precisions = ["Float64"]
horizons = [102400]
release = "closure"

[[run.arm]]
name = "ind_pa1"
orbits = "P"
indicators = ["pa1"]
precisions = ["Float64"]
horizons = [102401]
release = "closure"

[[run.arm]]
name = "d64"
orbits = "D"
indicators = ["tg", "fq1", "bk1"]
precisions = ["Double64"]
horizons = [102400]
release = "closure"

[[run.arm]]
name = "d64_pa1"
orbits = "D"
indicators = ["pa1"]
precisions = ["Double64"]
horizons = [102401]
release = "closure"

[[run.arm]]
name = "chk"
orbits = "D_B"
methods = ["L5"]
precisions = ["Double64", "BigFloat"]
horizons = [1000000]
release = "closure"
```

4. Panel validity and the hash of this file. The library accepts a panel only if the hash under which the panel was built, which the panel records, is in `panel_prereg_sha256_accepted`. That list is inside the hashed file, and the hash rule replaces only the header cell by SELF, so a file cannot list the hash of its own final version. Rule. (a) The panel is built under the version of this file that Amendment A1 leaves; call its header hash H_1. In that version the list holds the placeholder of 64 zeros shown above. (b) Amendment A3, the seal, is logged after A1 and before the panel build. It replaces the placeholder by H_1 and changes nothing else except the header hash cell and the amendment log cell; its log entry names H_1 as the hash it replaces. (c) The panel build reads the file as of the commit that makes A1 (a copy kept in the stage 02 output folder next to the panel; its path is the `prereg` path in the panel manifest), not the sealed file. (d) The round log, the chunk files and the run manifest carry the header hash of the sealed file. (e) The sealed file and the build copy differ in exactly three places: the list entry, the header hash cell and the amendment log cell; the check is the diff of the two. A3 is allowed only between A1 and the panel build and its content is a hash; it touches no other value.

5. Additions to "Amendment rule". (a) A pilot run (experiments spec R9) that uses no panel orbit and no test screening, and writes outside the study store, is not test data. (b) The order of amendments is A2, A1, A3 (item 1). (c) The build of the test panel is the first test screening; after it the file is not changed, A3 included, because the round log and the chunk files carry its header hash.

6. Budget note (screening of K = 0). The panel builder screens K = 0 with the same rule and size as the other K (module amendment D94, P71), because the K = 0 class T needs a screening and the P1 lattice would reuse the P1 start cells. Screening is 6 x 262144 x 1e5 x 207.6 ns = 9.07 thread-hours before the safety factor (frozen Budget: 5 K, 7.6), 18.1 with the factor 2. The worst-case total with safety becomes 52.4 thread-hours (97.0 with the T17 lines), 4.6 wall hours (8.4 with the T17 lines) at 11.5 effective threads, below the wall cap of the register. BigFloat line of the chk arm: the frozen Budget charges BigFloat L5 at 18200 ns per orbit-step ("P1 line"), the 4-thread value of P1 Amendment A4, which P1 Amendment A5 replaced after the 23-worker pilot; the pilot table gives BigFloat L5 alone as 341296.4 ns (module-build stage 04, `experiments-runner-verify-log.txt`, PILOT-COST L5 BigFloat cpu). The block uses 341296.4 ns, safety 3 unchanged, so that the first-chunk cost cap (2 x 3 x 18200 = 109,200 ns against the measured cost) does not pause the arm at its first chunk. The chk line becomes 100 x 1e6 x (341296.4 + 6634.6) ns = 9.66 thread-hours before safety, 29.0 with the factor 3 (frozen: 0.7 and 2.1). The worst-case total with safety becomes 79.3 thread-hours (123.9 with the T17 lines), 6.9 wall hours (10.8 with the T17 lines) at 11.5 effective threads, below the wall cap of 12 hours. No other Budget line changes.

## Amendment A1 (2026-11-04)

Approved by the owner at the checkpoint of the calibration freeze (2026-11-04), before any test data exist: no test screening, panel, draw, label, indicator run or pilot on a panel orbit exists. Base sections above and section "Amendment A2" are left as they stand; this section adds the frozen output of the calibration (Methods, "Calibration procedure", C0 to C10) and of procedure PP (Panel, "Test draw size and power calculation"), as section "Amendment rule" plans. It changes no hypothesis, no threshold of the register, no set, no rule, no arm, no block and no budget line. The calibration record is the single home of every number of the calibration; this section states the values the study reads and points to the record for the rest. After this section the rule is a fixed function (Amendment rule). Next: Amendment A3, the seal (section "Amendment A2", item 4).

1. Program and run. The calibration program is the command `calibrate` of the library (`bin/chaoslab.jl calibrate <preregistration> <parameters> --out <dir> --threads N`, module version 0.1.2) at commit 02994688354738e11405605e8a0074481a7ccb80, working tree of the code paths clean (`src`, `ext`, `bin`, `Project.toml`, `Manifest.toml`). It read this file at header hash f614a38449f611a42c56a0aa79bbc2d106dfaab679d1972a2b5349011df1a072 (the version after A2) and the 71-key table `workspaces/study/stages/01-design/references/p3-abstention-calibration-parameters.toml` (SHA-256 69965ec24cee507d76924423745cf332ecbb6c003880e5fb4b1036a60fc6a5e1), and the p2 and P1 stores read only. The input fingerprints (p2 manifest cdef3cf423e16a0b54e81e96fd460f4cb32d7b788325ffbf4e0581f1db27acce, P1 manifest ae6ec80eb08fbd6e98aeef1fbf840a567a03102ab943f2ae0aaf2b3e367f191d, panel file 881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5, and the p2 and P1 parameter tables and preregistrations) are in the record. Test seed: 20261030 (register). C9: the program was run twice, with 23 and with 4 threads, into two empty folders outside the data folder; both runs exited 0 and the two records and the two tables files are byte-identical. The SHA-256 of the data tree was the same before and after the two runs.

2. Outputs. The record `workspaces/study/stages/04-analyze/output/p3-abstention-calibration-record.json` (SHA-256 a8cfdbe12a56697fba17db12b26cafced06647fd05dfa0db79af9a1319b3680a) and the tables `workspaces/study/stages/04-analyze/output/p3-abstention-calibration-tables.md` (SHA-256 14b8f6e096a7fc29b593177531351ff5e075efdabe4b834a059d14f8cce228de; table T1 and the data of figure F1). The folder is the one the commands `count-round` and `analyze` read the record from by default. Outcome: selected (F is not empty, C8 does not apply).

3. Relabelling and sets (C0). 52 panel orbits changed class under the reduced rule, no grid orbit; 171 Cal orbits dropped (13 controls, 158 unresolved); no orbit lacked relabelling data. Cal has 3016 orbits (831 panel, 2185 grid): S2 452, T_s 31, W2 228 (R_W 197), regular 1166, gate 80.

4. Feasible set and theta* (C2 to C6). Theta has 2088 elements; the cost cap of F1 is 185395200 ns; |F| = 382; J_min = 0; |Theta_tol| = 20. theta* (index 371 of the record): V = SALI, FLI, WBA; q_C = 1; q_R = 3; w = 2; v = 0. Its cost is C_R(theta*, b_P) = 72430080 ns with the frozen cost-model values (Methods, "Cost model and matched budget": L1 207.6, tg 2658.0, pa1 1810.5, bk1 3000.6, fq1 16850.4; Double64 ratios tg 10.5, fq1 1.38, pa1 7.28, bk1 6.56). On Cal: calibration error 7/1432 (0.0049) among the decided, unresolved share 19/377 (0.050), unresolved share on S2 0, gate count 0, J = 259/14136 (0.018). Matched budgets (ceiling rule): SALI, FLI, MEGNO and E2 51200, WBA 25600, E1 6400; LE is read at N_c = 102400 and is cost-unmatched.

5. HB analogue on Cal (C7, fitted, optimistic, not evidence). The rule is false-regular on 1 of 31 T_s orbits and 1 of 228 W2 orbits, with decided shares 25/31 and 109/114. The closest comparator is WBA at 25600: 4 of 31 on T_s (exact McNemar p = 3/16) and 6 of 228 on W2 (p = 1/16); every other comparator has p below 0.04 on both sets. The full table, the top 20 of Theta_tol and the risk-coverage points are in the tables file.

6. PP for HB. Optimism over 200 resamples (none skipped): T_s 0.0269, W2 0.0045. Rule rates used: T_s 0.0591 (r_corr), R_W 0.01 (the floor). Comparator rates, shared shares and decided shares by stratum are in the record (`pp_hb.rates`, `pp_hb.model`). No allowed round count up to R_max reaches p_HB: the joint power is at most 0.3602 (rho = 11/4) and 0.30815 at the cap (6163 of 20000 runs), so the cap binds and the targets are the cap, n_T = 285 and n_W = 1235 (procedure PP, step 5). The study runs anyway and HB is expected to be underpowered; HB may be NEGATIVE, and its planned power is the power at the cap (owner decision 2026-10-31 on finding V2 of the module-build verification of D93: run as frozen, no change of the register). Type-I check (step 6): the largest joint rejection rate is 0.0005 (rule rate set to the WBA rate).

7. PP for HC. AUC_I_cal = 0.9146 on 51 failures and 2885 non-failures; AUC_I_used = 0.7073; AUC_P = 0.55. Target n_f = 45 (power 0.8161). Expected deliverable failures 49 (worst case 62), with power 0.8437; HC is not underpowered.

8. Values read by the study. The test-side table `workspaces/study/stages/01-design/references/p3-abstention-test-parameters.toml` quotes the five pilot costs from Methods, "Cost model and matched budget" (the sentence on the p2 pilot), and the following values from this section; the library checks them against the record and the 71-key table.
- calibration_record_sha256 = a8cfdbe12a56697fba17db12b26cafced06647fd05dfa0db79af9a1319b3680a
- p3_table_sha256 = 69965ec24cee507d76924423745cf332ecbb6c003880e5fb4b1036a60fc6a5e1
- target_t_s = 285
- target_w2 = 1235
- target_f = 45
- theta_voters = SALI, FLI, WBA
- theta_q_c = 1
- theta_q_r = 3
- theta_w = 2
- theta_v = 0
- power_hb_planned = 0.30815
- power_hc_planned = 0.8437
- cap_flag = true
