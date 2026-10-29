# Preregistration: p2-indicators

FROZEN 2026-10-29. The scope (claim C3, detection latency, classification stability; value types FALSIFIES, NEGATIVE, TOOL) was confirmed by the owner at checkpoint 1 on 2026-10-29. Every value marked "proposed" was approved by the owner at the study/01 checkpoint (2026-10-29, all as drafted); the one-line rationales stay as the record; the numbered list is in "Open points for the owner at the checkpoint". All thresholds of this study live only in this file. Other facts are pointed to, not copied; the few counts and hashes that this study needs are given with their source.

## Header

| Field | Value |
|-------|-------|
| Study slug | p2-indicators |
| Paper | joint P1+P2 paper (D65, `/_config/research-question.md`, paper sequence); this study is the P2 part |
| Claims tested | C3 (diagnostic inequivalence), with the two P2 outputs of the paper row: detection latency and classification stability. C1 and C2 are not tested here: the accepted P1 results are reused as ground-truth labels (D61, D65) |
| Value types | FALSIFIES (H1 to H6), NEGATIVE allowed, TOOL (`/_config/value-types.md`). NOT claimed: NOVEL (no prior-art check was made for this study, K9) and REPRODUCES (no row of `benchmarks.md` is matched) |
| Git SHA at freeze | 71456146c8aeafa025142b8c2c96f12a6d9cebee (parent commit of the freeze commit) |
| SHA-256 of this file (current version) | f43875763b404fb408c11ea7793ddec812f0ab0a2cd0cdc3cf11896e5afe70df (rule: the whole cell is replaced by the word SELF before hashing, UTF-8, LF; as P1 Amendment A1 item 1) |
| Amendment log | none yet. Planned before any data exists: A1, the generated `[[run.h4]]` tables of the grid (see Panel, "Grid G") |

Value types, in the words of `value-types.md`:
- FALSIFIES: H1 to H6 each have a metric and a numeric decision rule below, fixed before any data exists.
- NEGATIVE: a hypothesis that is not supported is reported with the same figures and tables as a supported one.
- TOOL: the run uses only modules whose catalog `Status` is `verified` and the runner of decisions D53 to D59 and D86 to D88. The TOOL claim cites a tagged release that contains the M4 modules and these amendments; the tag is cut at study exit and its number is the owner's decision (not fixed here).

## Scope check against the catalogs

Read on 2026-10-29 from `/_config/catalogs/indicators.md`, `lyapunov.md` and `systems.md`.

| ID | Role here | Catalog Status | In scope |
|----|-----------|----------------|----------|
| I1 / L1 Benettin | indicator LE: lambda_1(n), reference estimator of the study (see Methods) | verified | yes |
| L5 Householder | reference of the P1 label rule; used only to label the grid | verified | yes |
| I2 SALI (with GALI_2 = exp(L), D70) | indicator SALI | verified | yes |
| I5 FLI | indicator FLI | verified | yes |
| I8 MEGNO | indicator MEGNO | verified | yes |
| I9 RLI | indicator RLI | verified | yes |
| I10 and L13 stretching numbers | indicator STR (spectral distance D_str) | verified (L13: distribution, mean, second moment of one orbit's stretching numbers) | yes |
| I12 NAFF / frequency map, series path E1 | indicator E1 | verified | yes |
| I12b weighted phase advance, path E2 | indicator E2 | verified | yes |
| I13 weighted Birkhoff averages | indicator WBA | verified | yes |
| Y1 standard map, Y5 linear symplectic map | systems (Y5 only as the exact control of P1) | verified | yes |
| E precision (Float64, Double64), F diagnostics | precision and checks | verified | yes (BigFloat is refused by the runner for indicator kinds) |
| I3 GALI_k (k > 2), I4 | not run | planned | EXCLUDED (in 2D, GALI_2 is the same indicator as SALI, D70) |
| I6, I7 OFLI | not run | planned | EXCLUDED |
| I16 recurrence, I12c, I12d, I14, I15 (0-1 test), I17, I18 (ML), I11 | not run | planned | EXCLUDED |
| Y2 Henon map, L14 FTLE fields | not run | planned | EXCLUDED |

Indicator sets used below. I_7 = {LE, SALI, FLI, MEGNO, E1, E2, WBA}. I_8 = I_7 plus RLI. STR has no threshold (catalog I10 row) and enters only the threshold-free statistics. Names: LE is I1/L1, E1 and E2 are the two frequency paths of the naff-fma spec, WBA is I13.

## Hypotheses

Notation and metric definitions are in "Metrics"; the sets W2, S2, T_s are defined in "Panel". Every statistic with a bootstrap interval uses 10^4 resamples of orbits and the seed of `numerics-defaults.toml` (proposed, the P1 convention). A set with fewer than 20 orbits (min_set, proposed, the P1 value) is "not decidable": reported with its count, neither pass nor fail. The primary decisions use I_7 on the orbits of E_all. The same rules with I_8 on the orbits of E8 (the orbits that carry RLI) are an extension: a decision is called ROBUST if both give the same verdict. The class rules are the primary rules of Metrics, table R1; the fixed-threshold rules of table R2 are a sensitivity (the same statistics, reported, not decisive).

### H1 (claim C3, classification divergence at a fixed budget)

Metric: on W2 at the budget b_P = 12800 (proposed: the member of the budget list nearest to the n_H2 = 1e4 of P1), the unanimity-failure fraction U(b_P) of I_7 and the smallest pairwise agreement A_min(b_P) over the 21 pairs of I_7 (definitions in Metrics, M6). Decision: H1 is SUPPORTED (the indicators classify differently) if U(b_P) >= 0.20 (proposed: one weak-chaos orbit in five gets split votes), the lower end of its 95 percent bootstrap interval is >= 0.10 (proposed), and A_min(b_P) <= 0.80 (proposed). Otherwise H1 is not supported and the null (the indicators agree on weak chaos) is reported (NEGATIVE). U on S2 at b_P is reported as a negative control, without a rule.

### H2 (claim C3, detection-time divergence)

Metric: on S2 (strong chaos), M_X = median over S2 of the censored detection time t_det(X) (Metrics, M4), for each X in I_7, and Q = max_X M_X / min_X M_X. Decision: H2 is SUPPORTED if Q >= 10 (proposed: one decade, the P1 H3 convention) and the lower end of the 95 percent bootstrap interval of Q is > 3 (proposed). Otherwise not supported (NEGATIVE).

### H3 (detection latency, regime dependence)

Metric: for each X in I_7, R_X = M_X(W2) / M_X(S2), the ratio of the censored median detection times in the weak-chaos set and the strong-chaos set. Decision: H3 is SUPPORTED if for at least 5 of the 7 indicators R_X >= 10 and the lower end of the 95 percent bootstrap interval of R_X is > 3 (proposed: a majority of indicators, and the P1 H3 numbers 10 and 3). Otherwise not supported (NEGATIVE).

### H4 (claim C3, sticky orbits)

Metric: on T_s (the sticky ensemble of Panel, R4), for each X in I_7 minus LE, phi_X = fraction of T_s orbits with c_X(N_c) = regular (a false regular verdict at the common horizon N_c = 102400, Metrics M1), and Delta_phi = max_X phi_X minus min_X phi_X. LE is left out because T_s is defined from lambda_1(n_s) (circularity). Decision: H4 is SUPPORTED (indicators answer differently on sticky orbits) if Delta_phi >= 0.50 and the lower end of its 95 percent bootstrap interval is >= 0.25 (proposed). Otherwise not supported (NEGATIVE). The a-priori ensemble T_s^design (panel design class S) is reported with the same statistic.

### H5 (classification stability over the budget)

Metric: on W2, F_X = fraction of orbits with c_X(b_P) different from c_X(N_c) (the verdict at the realistic budget is not the verdict at the common horizon), for each X in I_7, and F_med = median over I_7 of F_X. Decision: H5 is SUPPORTED (finite-time verdicts are not stable) if F_med >= 0.10 and the lower end of the 95 percent bootstrap interval of F_med is >= 0.05 (proposed, the P1 H2 numbers). Otherwise not supported (NEGATIVE).

### H6 (classification stability under direction and precision)

Metrics (both at N_c): (a) on the designated subset S_D (orbits of E_all in it), for X in {LE, SALI, FLI, MEGNO}, delta_dir(X) = fraction of orbits whose verdict is not identical in the five deviation directions of the direction arm; (b) on the ladder subset S_L (orbits of E_all in it), for X in I_7, delta_prec(X) = fraction of orbits whose verdict differs between Float64 and Double64. Decision: H6 is SUPPORTED (verdicts are not input-stable) if max_X delta_dir(X) >= 0.05 or max_X delta_prec(X) >= 0.05 (proposed: any change in one orbit in twenty; P1 treated any label change as "unresolved"). Otherwise not supported (NEGATIVE).

### Reported without a decision rule

Agreement matrices and Cohen's kappa for all pairs of I_8 at every budget, the AUC of every indicator including STR (M7), the regular-range tables of the sticky orbits (Metrics, M10), the fixed-threshold rules of table R2, the sensitivities of the observable, of mean removal and of the planar phase advance, the packages of research idea 17 (the Maffione2011a pairs FLI+RLI and MEGNO+SALI against the reference class, no rule), the flag counts, the cost per orbit-step, and the rule-validity gate on K = 0 (Analysis plan). The slope profile and coarse-graining ablation of research idea 25 are not in this study.

## Panel

Three orbit sources: the frozen P1 panel P (reused), a K x initial-condition grid G (new), and the sticky ensemble T_s (a measured subset of P and G, so no extra starts). Stage 02 (panel builder) is not rerun for this study.

### Panel P (reused by reference)

- File: `workspaces/study/stages/02-panel/output/p1-le-benchmark-panel.toml`, SHA-256 `881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5` (the value recorded by the accepted P1 run, `workspaces/study/stages/05-report/output/p1-le-benchmark-results.md`, Header). The panel records the preregistration hash `650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46` as its build hash; the run block below lists it under `panel_prereg_sha256_accepted`. The K values, the cells, the design classes (T, R, C, S, X), the selection rule, the seeds and the subsets are those of the P1 preregistration, section "Panel" and Amendments A1 to A3; they are not restated.
- Counts read from the panel manifest on 2026-10-29 (the panel file is their home): 883 orbits (870 drawn, 13 controls with ids 871 to 883); K index 0 holds the controls, K index 1 to 6 are K = 0, 0.5, K_c, pi/2, 2 pi, 10; the designated subset S_D has 172 orbits (all 13 controls among them), the ladder subset S_L has 92 orbits (no control); design class S has 110 orbits (K = 0.5: 25, K_c: 40, pi/2: 40, 2 pi: 5). The Y5 control is orbit 883 (system `linear`, start (0, 0), eigenvalue pair (mu, 1/mu) with mu = (3 + sqrt(5))/2).
- P1 labels as ground truth: the accepted P1 results (decisions D61 and D65), that is, the two label axes of D36 and D37 under the "Label rules" of the P1 preregistration, as produced by the Analysis module (`bin/chaoslab.jl analyze`) from the P1 run store. The per-orbit label table is not among the tracked outputs (the tables file holds counts per cell, T6); stage 04 of this study reads the labels from the P1 store, never recomputes them with other rules. P1 labels are taken at the longest Float64 horizon of each orbit, as P1 defines them.
- Exclusions: the 13 controls leave every hypothesis set (they enter only the gate, flag and cost tables); orbit 881 stays dropped as drift, as in P1 (lifted |p| above the drift bound at any recorded time, the same flag in the new runs).

### Reference class rho(o)

Each orbit o of P and G has a reference class rho(o) in {C, R, U} (chaotic, regular, unresolved), the ground truth of every correctness statistic:
- Panel orbit: rho = C if instability = "positive growth supported" and structure = "irregular"; rho = R if instability = "none detected on the tested horizon" and structure = "quasiperiodic"; rho = U in every other case (any unresolved axis, structure "periodic").
- Grid orbit: the same two label rules (P1 "Label rules", instability axis from L5 at N_1 = 1e7, structure axis from the rotation numbers of the recorded orbit points at N_1/2 and N_1), with the disagreement tests that the grid arms provide: |lambda_L1 - lambda_L5| above the P1 H1a bound at N_1 gives U. The precision ladder and the direction test of P1 do not exist for the grid; every grid orbit carries the flag ref_weak and all tables split the results by source (panel, grid).
- lambda_ref(o) = lambda_1(N_1) of the L1 item of this study (arms `lyap_p`, `lyap_g`), N_1 = 1e7. For the panel this equals the P1 L1 value (same orbit, same code path); P1 showed L1 and L5 agree to 1e-8 on all of its eligible orbits.

### Grid G (new starts, K x initial-condition)

- K values (8, proposed): 0.5, 0.75, K_c = 0.971635406 (golden-mean breakup, FrahmShepelyansky2010, `benchmarks.md`), 1.25, pi/2, 2, 3, 2 pi. The four values 0.5, K_c, pi/2, 2 pi coincide with the P1 values so that grid and panel cells share K; 0.75, 1.25, 2, 3 fill the transition from the mixed to the mostly chaotic phase space. K = 0 and K = 10 are left to the panel (K = 0 has no chaos, K = 10 has no regular region in P1).
- Initial conditions (proposed): at each K, the cell-centred lattice x0_i = 2 pi (i + 1/2) / 24, i = 0 to 23, and p0_j = pi (j + 1/2) / 12, j = 0 to 11, so p0 lies in (0, pi). The standard map is reversible under (x, p) -> (-x, -p), which maps p0 in (pi, 2 pi) onto (0, pi); the half range avoids pairs of mirror orbits that are the same orbit twice (the bootstrap treats orbits as independent). 288 starts per K, 8 x 288 = 2304 starts. Values are computed in BigFloat(256) and written as 17-digit decimal texts, one rounding.
- Carrier in the runner: the starts are `[[run.h4]]` tables (K, x0, p0, label), labelled `grid-k<K index 1 to 8>-i<00 to 23>-j<00 to 11>`; the runner gives them ids after the last panel id (884 to 3187) and runs them as the subset `h4`. The 2304 tables are generated mechanically from this rule and appended to the run block by Amendment A1 before any data exists; this section is their single home. The rule, the count and the labels are fixed here; a different count is an amendment.
- Reference arms for G: `lyap_g` (L1 and L5 at N_1 = 1e7, the same horizon as the P1 labels).

### Weak-chaos set W2, strong-chaos set S2, sticky ensemble T_s (R4)

All three are subsets of E_all (below) with rho = C. Numbers proposed unless stated.
- S2 (strong chaos): lambda_ref >= 1.0. In P1 this is the chaotic sea at K = 2 pi (1.17) and K = 10 (1.62), the "strong-chaos set" of P1 H3.
- T_s (sticky, R4): lambda_1(n_s) <= f_s lambda_ref, with n_s = 1e5 and f_s = 0.2. Derivation: the only held definition of stickiness is P1's design class S, a stuck stretch of the windowed exponent below 0.2 times its median (P1 Amendment A1 items 2 to 4, f_S = 0.2) found on the screening run of N_scr = 1e5 steps. The standard-map spec holds the map, its tangent, the fixed points and the periodic-orbit data, but no island boundary or cantorus (search for "stick", "island", "cantor" in `standard-map-spec.md`: only the Meiss2015 map form and the Contopoulos island-centre test). A geometric definition would need new literature, so R4 reuses the P1 construction on the quantity this study records for every orbit: an orbit is sticky at the P1 screening horizon n_s if its finite-time exponent has reached less than a fifth of its value at N_1. The P1 results show why N_1 is the minimum reference horizon: the sticky cells converge at 5e5 to 7.1e6 steps (P1 T4, class S rows). Limit (stated, not corrected): extreme sticky zones may stay trapped beyond 1e7 (HarsoulaKaramanosContopoulos2019, benchmarks row: not converged after 5e7); such an orbit gets rho = R and the indicators that call it regular are counted correct. T_s is measured from L1, so LE is excluded from H4. The a-priori ensemble T_s^design = orbits of P with design class S and rho = C is reported next to it.
- W2 (weak chaos): lambda_ref <= 0.1, or the orbit is in T_s. The value 0.1 sits at the K_c chaotic sea (P1 comparators table) and below pi/2, where lambda is about 0.3.
- Regular set: rho = R (used for the AUC and the regular range, M7 and M10).

### Eligible sets

- E_all = (P union G) minus controls (ids 871 to 883, including Y5), minus the drift orbit 881, minus orbits with rho = U. All hypotheses use sets built from E_all.
- E8 = orbits of E_all that carry RLI: members of S_L (the pair arm `pair_l`) and all of G. The extension of each hypothesis uses E8 and I_8.
- Orbits for which a spec is undefined get no item (Y5 has no torus: E1 and WBA are omitted for orbit 883; the plan lists the omission and the plan hash covers it). A verdict that is "undecided" (Metrics, M2) removes that orbit-budget from the statistic of that indicator and is counted in the table.
- K = 0 gate set: the 80 orbits of the panel cells (K index 1, classes T and R); used only for the rule-validity gate (Analysis plan).

## Methods

Every method of the scope-check table is `verified`. Definitions, tolerances and algorithm identifiers live in the catalogs and the module specs; they are not copied. The recording is done by the study runner `Chaoslab.Experiments` under Amendments D86 to D88 (spec `workspaces/module-build/stages/01-spec/output/experiments-spec.md`, sections "Amendment D86", "Deviations accepted at the D86 stage 03 checkpoint (D87)", "Amendment D88"; verification `workspaces/module-build/stages/04-verify/output/experiments-d86-verification.md`; layout `schema/run-store.md`, stated by spec requirement R15). Each indicator is a spec of one of the six kinds `tangent`, `pair`, `stretching`, `frequency`, `phase_advance`, `birkhoff`, listed by the arm key `indicators`; the parameters below are the required keys of the kinds, and the module holds no default and no threshold.

Reference estimator for lambda_1(n): L1 Benettin (single tangent vector). Reasons: it is the cheapest verified estimator; the column-1 log sums of the tangent chain are bitwise those of L1 (spec invariant I22), so LE and the chain see the same quantity; and on the P1 eligible set the spread over L1 and the QR family was at most 1e-8 for every orbit (P1 results, H1a, 209 of 209). The catalog marks L5 Householder as the default reference of the QR family; L5 is used here only to label the grid orbits by the P1 rule, which uses L5.

### Indicator specs

| Spec id | Kind (indicator) | Members or observable | Windows or horizon | Other keys |
|---------|------------------|-----------------------|--------------------|------------|
| tg | tangent (SALI, FLI, MEGNO) | members FLI, MEGNO, SALI; SALI last, replacing form set by the module (E23) | recording grid of the block | columns of E24: `sali`, `L`, `fli`, `fli1`, MEGNO `Y`, `Ybar`, `Ybar20` for both columns |
| pr | pair (RLI, MEGNO) | members RLI, MEGNO | recording grid | separation 1e-12, frame_source `seeded` |
| st | stretching (STR) | | recording grid, D read at each grid step | range `derived`, bins 1000, interval 100, min_per_cell 100, counts `S_D` |
| fq1 | frequency (E1) | observable `exp_i_q` | windows W in {100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200} | order 1, positions `half`, pad 8, mean_removal true, mirror_ratio 0.9, floor_multiplier 32 |
| fq1_raw | frequency (E1, sensitivity) | `exp_i_q` | same windows | as fq1 with mean_removal false |
| fq1_qp | frequency (E1, sensitivity) | `q_plus_i_p` | same windows | as fq1 (mean_removal true) |
| pa1 | phase_advance (E2) | observable `lifted` | same windows | birkhoff_order 1 |
| pa_planar | phase_advance (E2, sensitivity) | `planar` | same windows | birkhoff_order 1 |
| bk1 | birkhoff (WBA) | cos q, sin q, cos p, sin p (coordinate 1 = q, 2 = p; frequency 1) | same windows | birkhoff_order 1 |

### Choices made explicitly (each value "proposed" unless a source is named)

1. mean_removal = true for the primary E1 spec fq1. Cost, known: on the K_std = 0 line the recovered frequency differs from the closed form p_0 / (2 pi) by up to about 3e-9 (22 of 48 values outside 16 eps), while with mean_removal = false all 48 values are within 16 eps (worst ratio 0.012) (verification report, owner notes; the closed form of naff-fma criterion 7(b) is for the second case). Reason for true: the amplitude maximum of exp(i q) is taken over the whole spectrum; for an orbit whose angle dwells near one point (the X class near the hyperbolic fixed point, thin layers, island orbits) the window mean is of order 1, and without removal the zero-frequency line can be the maximum in both windows (naff-fma spec P4 and F5). The distance then sits at the floor and the verdict is regular because of the estimator, not of the dynamics. The bias of 3e-9 gives D of about -8.5, five decades below the decision level of table R1 (D >= -4), so it cannot flip a verdict. The sensitivity spec fq1_raw (mean_removal false) on S_L measures the difference in verdicts and in D (table T11).
2. Windows (E1, E2, WBA): ten windows, doubling from 100 to 51200 (each at least 8 x (order + 1) = 16). The run length is 2 W_max = 102400 handoffs for E1 and WBA and 2 W_max + 1 = 102401 for E2 (spec E22), so the arm horizons are 102400 and 102401 at tau = 1. The pair [0, W), [W, 2W) is evaluated at the budget 2W, which gives the budget list B of Metrics, M1. Smallest budget 200: the strongest chaos of P1 converged at 158 to 266 steps (P1 T4, K = 2 pi and 10). Largest budget 102400: close to the screening horizon n_s = 1e5 of R4, so the sticky horizon is covered. Memory: the E1 transient buffer is 32 W_max = 1.64 MB per orbit in Float64 (spec E28).
3. E1 transform: Hanning window of order 1 (`order` = 1), sample positions (j + 1/2) / W (`half`), coarse-grid factor `pad` = 8, mirror ratio 0.9, index floor 32 eps (`floor_multiplier` = 32). These are the proposals of naff-fma spec P5 and P10 (P10 values are drafting values there, not sourced from a held text). Refinement: the frequency is refined by the stationary-point search in the element type after the coarse search; `pad` is the only refinement key of the grammar, and the flags `refined1`, `refined2` are outputs, not choices.
4. E2: observable `lifted` (increments of the lifted first coordinate, the rotation number of the standard map; the naff-fma default for E2, free of the branch ambiguity of the planar increments, flag `near_pi`); weight w^[1] (`birkhoff_order` = 1; a higher order helps only beyond 15-digit arithmetic, Das2017 as reported in naff-fma P5). The sensitivity spec pa_planar uses the planar observable.
5. WBA: the four observables cos q, sin q, cos p, sin p of the wrapped angles (the label observables of the weighted-birkhoff spec; research idea 25 uses the same set), `birkhoff_order` = 1. K = 4 observables and H = 10 windows give reducer bytes 16 + 8 H (2 + 5 K) = 1776 (spec E28).
6. STR: range `derived` (the module derives the cell range from the system; standard maps and 2 x 2 linear maps only), bins 1000, interval 100, min_per_cell 100 (Maffione2011a pp. 4 and 20 through stretching-numbers spec P4). The horizon 102400 reaches min_per_cell x bins = 1e5 handoffs, so the scalar `first_full_step` is positive for every orbit. counts = `S_D` stores full count matrices for the designated subset (8 kB per orbit). By D88 the pooled counts of a chunk are written only when its orbits share one range: the panel cell chunks (one K) pool; the controls chunk and the single grid chunk (eight K) do not.
7. Pair (RLI with MEGNO): separation 1e-12 and wrapped coordinates (rli spec R8, Float64), frame_source `seeded` (spec P32, the reading of rli P2 that makes the printed regular law possible). Finding F2 of the verification report: on the Y5 control (orbit 883, start at the origin) the shadow orbit of the pair overflows at step 770 and the chunk fails twice. Choice: Y5 is excluded from the pair arms, not the horizon capped. The arms `pair_l` and `pair_g` use the subsets `S_L` and `h4`; orbit 883 belongs to neither (S_L has no control; h4 holds no panel orbit). A horizon cap below 700 steps would remove RLI from every budget of B above 400, so the cap was rejected. Consequence: RLI exists for 92 of the 883 panel orbits (4 per populated cell) and all of G. The runner has no subset "panel without controls"; see the open points.
8. Coordinates and tau: wrapped coordinates in every arm (E1 `exp_i_q`, WBA `cos` and `sin` and the pair need them; wrapped and lifted L1 agree, P1 F7), tau = 1. Primary deviation direction pi/4 as in P1; the direction arm adds the four P1 directions 0.3, 0.3 + pi/5, 0.3 + 2 pi/5, 0.3 + 3 pi/5. `eps` of the block is the P1 key, unused (no L2 item).
9. Precision per arm: Float64 in all main arms. Double64 on S_L (92 orbits) for L1, tg, fq1, pa1 and bk1 (the stability arms d64_*); not for pr and st. BigFloat is refused by the runner for indicator kinds. Double64 is about 47 times the Float64 cost (P1 K20), so it is a stability check, not the main run. The decision levels are far from the Float64 floors (WBA digits floor 14.15 against the level 6; D floor -14.15 against -4).
10. Backend: `auto` for Float64, `cpu` for Double64. By D86 (R16) STR and E1 are CPU-only; tangent, pair, E2 and WBA follow the per-spec device record at the run's thread count. Before the first run the record per spec is taken at 23 threads (P1 K19 rule; the earlier record is kept as `prior_record`, D58). Stored values on the CPU do not depend on threads (spec I24); across backends only the verified tolerance of lambda_1(n) holds and indicators on chaotic orbits are not compared across backends (spec R20 item 2).

Defaults of `/_config/numerics-defaults.toml` used unchanged: `seed`, `renorm_interval_tau` (= tau = 1), `bigfloat_bits` (unused), and the three health multipliers (they apply to L1 and L5; the indicator items carry the determinant and symplecticity monitor). Overrides: `precision` per arm as above; `backend` as above. No key of the file is filled by this study.

### Run block (draft)

The block is the machine-readable form of the arms and specs. Key names are those of the P1 blocks (Amendments A4 to A6) and of experiments spec R11. Syntax was checked with a TOML parser at drafting; the runner has not read it (see "Checks I ran"). The `[[run.h4]]` tables are appended by Amendment A1. Arm names, spec ids and labels are data.

```toml
[run]
panel_prereg_sha256_accepted = ["650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46"]

[run.primary]
tau = 1
eps = 1e-8                    # carried from the P1 block; unused (no L2 item in this study)
direction = "pi/4"
coordinate_mode = "wrapped"

[run.grid]
points_per_decade = 20
first_steps = 100
extra_fractions = [0.1, 0.5]

[run.backend]
Float64 = "auto"
Double64 = "cpu"

[run.health]
fail_fraction = 0.05
control_bound = 1e-8
control_methods = ["L1"]
control_scope = "exact_start"
setup_check_steps = { T = 100, R = 100, control_elliptic = 100, control_hyperbolic = 5, C = 5, S = 5, X = 5, h4 = 5 }
setup_check_tol = 1e-10
drift_bound = 1e6
rerun_attempts = 1
stop_failed_fraction = 0.05
first_chunk_cost_factor = 2.0
wall_cap_hours = 12

[run.budget.Float64]
ns_cpu_thread = 207.6
ns_gpu = 2.468
safety = 2

[run.budget.Double64]
ns_cpu_thread = 6634.6
safety = 3

[run.k19]
threads = 23

[run.verify]
backend_agreement = { per_orbit_sigma = 6.0, mean_standard_errors = 4.0 }
rerun_sample = 20

[run.pilot]
orbits_per_k = 4
horizon = 100000
seed_label = "p2-indicators-pilot"

[[run.indicator]]
id = "tg"
kind = "tangent"
members = ["FLI", "MEGNO", "SALI"]
[run.indicator.budget.Float64]
ns_cpu_thread = 1000
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 6000
safety = 3

[[run.indicator]]
id = "pr"
kind = "pair"
members = ["RLI", "MEGNO"]
separation = 1e-12
frame_source = "seeded"
[run.indicator.budget.Float64]
ns_cpu_thread = 13000
safety = 2

[[run.indicator]]
id = "st"
kind = "stretching"
range = "derived"
bins = 1000
interval = 100
min_per_cell = 100
counts = "S_D"
[run.indicator.budget.Float64]
ns_cpu_thread = 1100
safety = 2

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
ns_cpu_thread = 13400
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 28100
safety = 3

[[run.indicator]]
id = "fq1_raw"
kind = "frequency"
observable = "exp_i_q"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
order = 1
positions = "half"
pad = 8
mean_removal = false
mirror_ratio = 0.9
floor_multiplier = 32.0
[run.indicator.budget.Float64]
ns_cpu_thread = 13400
safety = 2

[[run.indicator]]
id = "fq1_qp"
kind = "frequency"
observable = "q_plus_i_p"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
order = 1
positions = "half"
pad = 8
mean_removal = true
mirror_ratio = 0.9
floor_multiplier = 32.0
[run.indicator.budget.Float64]
ns_cpu_thread = 13400
safety = 2

[[run.indicator]]
id = "pa1"
kind = "phase_advance"
observable = "lifted"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
birkhoff_order = 1
[run.indicator.budget.Float64]
ns_cpu_thread = 1000
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 4200
safety = 3

[[run.indicator]]
id = "pa_planar"
kind = "phase_advance"
observable = "planar"
windows = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200]
birkhoff_order = 1
[run.indicator.budget.Float64]
ns_cpu_thread = 1000
safety = 2

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
ns_cpu_thread = 1000
safety = 2
[run.indicator.budget.Double64]
ns_cpu_thread = 5000
safety = 3

[[run.arm]]
name = "lyap_p"
orbits = "P"
methods = ["L1"]
indicators = ["tg"]
precisions = ["Float64"]
horizons = [10000000]
h1_check = true

[[run.arm]]
name = "lyap_g"
orbits = "h4"
methods = ["L1", "L5"]
indicators = ["tg"]
precisions = ["Float64"]
horizons = [10000000]

[[run.arm]]
name = "stretch_p"
orbits = "P"
indicators = ["st"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "stretch_g"
orbits = "h4"
indicators = ["st"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "pair_l"
orbits = "S_L"
indicators = ["pr"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "pair_g"
orbits = "h4"
indicators = ["pr"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "win_p"
orbits = "P"
indicators = ["fq1", "bk1"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "win_g"
orbits = "h4"
indicators = ["fq1", "bk1"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "phase_p"
orbits = "P"
indicators = ["pa1"]
precisions = ["Float64"]
horizons = [102401]

[[run.arm]]
name = "phase_g"
orbits = "h4"
indicators = ["pa1"]
precisions = ["Float64"]
horizons = [102401]

[[run.arm]]
name = "dir_tg"
orbits = "S_D"
methods = ["L1"]
indicators = ["tg"]
precisions = ["Float64"]
horizons = [1000000]
direction = ["pi/4", "0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 3pi/5"]

[[run.arm]]
name = "dir_st"
orbits = "S_D"
indicators = ["st"]
precisions = ["Float64"]
horizons = [102400]
direction = ["pi/4", "0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 3pi/5"]

[[run.arm]]
name = "sens_win"
orbits = "S_L"
indicators = ["fq1_raw", "fq1_qp"]
precisions = ["Float64"]
horizons = [102400]

[[run.arm]]
name = "sens_phase"
orbits = "S_L"
indicators = ["pa_planar"]
precisions = ["Float64"]
horizons = [102401]

[[run.arm]]
name = "d64_lyap"
orbits = "S_L"
methods = ["L1"]
indicators = ["tg"]
precisions = ["Double64"]
horizons = [1000000]

[[run.arm]]
name = "d64_win"
orbits = "S_L"
indicators = ["fq1", "bk1"]
precisions = ["Double64"]
horizons = [102400]

[[run.arm]]
name = "d64_phase"
orbits = "S_L"
indicators = ["pa1"]
precisions = ["Double64"]
horizons = [102401]
```

Notes on the block. An arm lists `methods` (exponent methods of the catalog) and/or `indicators` (spec ids). Arms keep one subset each (`P`, `S_D`, `S_L`, `h4`), so every spec that covers both the panel and the grid has one arm per source. A spec used by arms of different tau or horizon is allowed only where each arm horizon equals tau x run length (spec R11 item 3): E1 and WBA run 102400 handoffs, E2 runs 102401, hence the separate arms `win_*` and `phase_*`. `h1_check` marks the arm whose L1 items decide the exact-exponent control rule of the runner (control scope `exact_start`, bound 1e-8, as P1 Amendments A5 and A6 item 3; the L9 bound of A6 does not apply, no L9 item). `counts = "S_D"` and `range = "derived"` are explained in Methods, choices 6.

## Metrics

Exact definitions of every reported quantity. Symbols shared with the repo canon are used as in `research-question.md`. Scores are read from the stored series and tables of the run store; nothing is recomputed from orbits.

### M1 Budgets and horizons

- Budget list B = {b_k = 200 x 2^k, k = 0 to 9} = 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400 (proposed). The window kinds (E1, E2, WBA) are read at b = 2 W exactly, the pair [0, W), [W, 2W) (E2 uses one more sample). The series kinds (LE, SALI, FLI, MEGNO, RLI, STR) are read at n(b), the first recorded grid step not before b; with 20 points per decade n(b) <= 1.13 b. All indicators are therefore compared at the same list of budgets.
- Common horizon N_c = b_9 = 102400. Primary budget b_P = b_6 = 12800. Reference horizon N_1 = 1e7 (P1). Screening horizon of R4: n_s = 1e5.

### M2 Verdicts: score and class rule of each indicator

A verdict c_X(b) is "chaotic", "regular" or "undecided". Growth ratio: with n = n(b) and n' the recorded grid step nearest to n/10 (ties to the larger step), r_X(n) = g_X(n) / g_X(n'), where g_LE(n) = n lambda_1(n) (the cumulative log stretching), g_SALI(n) = L(n) (column `L`, ln GALI_2, negative, so the ratio of two negative numbers is positive), g_FLI(n) = column `fli`, g_MEGNO(n) = column `Ybar20_1` (the (2,0) form of column 1, the form of Maffione2011a on maps, megno decision D73). A regular orbit gives r of about 1.2 to 3 (s = ln n + b, L = -2 ln n + c, fli = ln n + b, Ybar = 2), a chaotic orbit about 10 (exponential or linear growth). A non-finite value or a non-positive denominator gives "undecided".

Table R1, the primary rules. "Chaotic iff" is the condition; otherwise "regular". All numbers are proposed and, where the source column names a level, they are placed between the regular and the chaotic level of the catalog row.

| X | Score read at n(b) | Chaotic iff | Where the number comes from |
|---|--------------------|-------------|-----------------------------|
| LE | lambda_1(n), r_LE | lambda_1(n) >= 1e-4 and r_LE >= 5 | the instability rule of P1 ("Label rules") evaluated at n instead of N; same two numbers on purpose, so that LE at N_1 reproduces the P1 instability axis. 5 lies between 1.2 and 10 |
| SALI | r_SALI | r_SALI >= 5 | as LE. In 2D maps regular SALI decays by a power law (catalog I2, sali spec F1), so a growth-law rule is used for the primary decision and the fixed threshold is in table R2 |
| FLI | r_FLI | r_FLI >= 5 | as LE |
| MEGNO | r_MEGNO | r_MEGNO >= 5 | as LE (regular Ybar stays near 2, benchmarks row Hinse2010) |
| RLI | column `RLI` | RLI >= 1e-7 | midpoint in the logarithm between the ordered level (proportional to the separation 1e-12) and the chaotic level, which the catalog I9 row puts 10.5 orders apart at Dx = 1e-12, n = 1e4 |
| E1 | D_E1, column `D` of the frequency table (log10 of the circle distance of the window pair, floored at 32 eps) | D_E1 >= -4 | regular median about -9.8 and chaotic about -0.6 at N = 1000 (catalog I12 row); the level moves with W, so -4 is placed 5 decades above the mean-removal bias (about -8.5, Methods, choice 1) and 3 decades below the chaotic level |
| E2 | D_E2, column `D` of the phase-advance table | D_E2 >= -4 | regular -11.5 and chaotic -0.6 at N = 1000 (catalog I12 row) |
| WBA | dig_min = minimum over the four observables of the column `dig` | dig_min < 6 | regular median 11.2 digits and chaotic 1.2 to 1.6 at N = 1000 (catalog I13 row); a regular torus converges for all four observables, so one unconverged observable is enough to say "chaotic" |
| STR | D_str | no rule | catalog I10: no threshold, the level moves with the number of points Z |

Undecided: a window row with the flag `degenerate`, `undefined window` or `zero` (E1, E2, WBA); a non-finite required value. The flags `below_floor`, `near_zero`, `mirror` and `near_pi` are data and do not change a verdict; their counts are reported.

Table R2, fixed-threshold sensitivity rules ("chaotic iff"; proposed): LE: n lambda_1(n) >= 3 ln n; SALI: SALI(n) <= 1e-8; FLI: fli(n) >= 3 ln n; MEGNO: Ybar20_1(n) >= 2.5 (0.5 above the quasi-periodic value 2). The other indicators have a single rule. R2 rules are horizon dependent by construction for SALI (power-law decay), FLI and LE (logarithmic growth); the table of results shows how much.

### M3 Pilot calibration of table R1

The numbers of table R1 and R2 are marked as calibratable on the pilot (Amendment rule): they may change by an amendment before any panel or grid data exist, on the evidence of the pilot of experiments spec R20, which runs on non-panel starts. Nothing else in this file is calibratable.

### M4 Detection time and outcomes

For an indicator X and an orbit o with rho(o) in {C, R}: the verdict series c_X(b), b in B up to N_c. A verdict is stable from b_k if c_X(b_j) equals c_X(N_c) and is decided for every j >= k, and b_k <= N_c / 4 = 25600 (proposed: the verdict persists over at least two doublings of the budget). The detection time t_det(X, o) is the smallest such b_k: the first budget at which the indicator reaches a stable verdict (research-question.md, "Detection time"). Outcome:
- correct: a stable b_k exists and c_X(N_c) = rho(o);
- wrong: a stable b_k exists and c_X(N_c) differs from rho(o);
- undetected: no stable b_k, or c_X(N_c) undecided.
Censored detection time: t~_det = t_det for correct, and 2 N_c = 204800 for wrong and undetected (a wrong or missing verdict is not a detection of the reference class; the P1 n_conv convention). M_X(A) = median of t~_det over the set A. t_det itself (any stable verdict) and the three outcome shares are reported for every X and set.

### M5 Stability

- Budget stability: F_X(b, A) = fraction of orbits of A, decided at b and at N_c, with c_X(b) different from c_X(N_c). Flip count Phi_X(o) = number of changes of verdict along B (reported).
- Direction stability: delta_dir(X) as in H6, from the arm `dir_tg` (LE, SALI, FLI, MEGNO; five directions, horizon 1e6, read at N_c).
- Precision stability: delta_prec(X) as in H6, from the arms `d64_*` against the Float64 items of the same orbits (read at N_c).

### M6 Agreement and divergence

For a set A, a budget b and two indicators X, Y with decided verdicts on the orbit:
- A_XY(b) = fraction of orbits of A with c_X(b) = c_Y(b); the 2 x 2 contingency table (agreement matrix) is reported for every pair of I_8 at every budget; Cohen's kappa_XY(b) is reported, with the value 0 when the expected agreement is 1.
- U(b) = fraction of orbits of A, decided for all of I_7 at b, whose seven verdicts are not all equal (unanimity failure). Orbits with an undecided verdict are dropped from U and counted.
- A_min(b) = min over the 21 pairs of I_7 of A_XY(b). The mean kappa over the 21 pairs and the Spearman rank correlation of the scores of table M7 between pairs are reported as further cross-method divergence statistics.

### M7 Threshold-free separation (AUC)

Chaos score of each indicator: LE r_LE, SALI r_SALI, FLI r_FLI, MEGNO r_MEGNO, RLI log10 RLI, E1 D_E1, E2 D_E2, WBA minus dig_min, STR minus log10 D_str (regular orbits have the larger D_str, catalog I10 and Darriba2012). AUC(X, b, K) = probability that a random orbit with rho = C has a higher score than a random orbit with rho = R at the same K (ties count one half). Computed per K and per budget, only when both classes have at least 20 orbits; reported, not decisive. STR enters only here.

### M8 Sets, bootstrap, undecidable sets

W2, S2, T_s, E_all, E8 as in Panel. Bootstrap: 10^4 resamples of the orbits of the set (with replacement), the statistic is recomputed on each resample, the 95 percent interval is the percentile interval, seed from `numerics-defaults.toml`. A set with fewer than 20 orbits is "not decidable".

### M9 Cost and flags

Cost per spec: thread seconds and wall seconds per orbit-step, bytes(orbit), launch capacity, transform seconds for E1, as the runner reports them (spec R19, R20). Flag counts per flag and chunk (`below_floor`, `degenerate`, `near_zero`, `mirror`, `tie`, `refined`, `near_pi`, `zero`, `reldig_undefined`, `undefined window`) as the manifest lists them.

### M10 Regular range of sticky orbits (naff-fma R4, weighted-birkhoff R4)

For X in I_7 plus STR and each K value, the regular range at budget b is the interval from the 5th to the 95th percentile (proposed) of the score of M7 over the orbits with rho = R at that K (panel and grid pooled, at least 20 orbits). For a sticky orbit, the leaving budget is the smallest b_k from which its score stays beyond the 95th percentile on the chaotic side for all later budgets up to N_c; the second quantity is the share of sticky orbits whose score at N_c lies below the median of the regular class. Both are reported, with no decision rule.

## Budget

Environment limits (`/_config/environment.md`): this computer, at most 23 cores, GPU RTX 5060 Laptop (8 GB), Double64 on the CPU. No wall-clock limit is stated there; a campaign cap of 12 wall hours is proposed (the P1 cap, about 4.5 times the pessimistic estimate below).

### Throughput inputs

| Quantity | Value | Source |
|----------|-------|--------|
| L1 and L5 Float64, CPU, 23 workers | 207.6 thread-ns per orbit-step (mean of the P1 methods; P1 measured 222 ns on the main arm) | P1 preregistration, Amendment A5 item 4, and P1 results (Cost) |
| L1 Double64, CPU, 23 workers | 6634.6 thread-ns per orbit-step | P1 Amendment A5 item 4 |
| tangent chain tg Float64 | measured 214 to 990 thread-ns per orbit-step at horizon 4000; line 1000 | `experiments-d86-verify-log.txt` and `-fixed.txt`, COST lines (arms lift, ser) |
| tangent chain tg Double64 | measured 2878 to 5928; line 6000 | same, arm ser_d |
| pair (RLI, MEGNO) Float64 | measured 7535 to 12940 at N = 300, where the per-orbit setup dominates; line 13000 (an upper bound for long horizons) | same, arm pair |
| stretching Float64 | measured 298 to 1064 at 100 bins, N = 4000; line 1100 (1000 bins add about 2 operations per step at 200 grid steps) | same, arms lift, ser |
| frequency E1 Float64 | measured 2073 to 6436 thread-ns per orbit-step at windows 250 and 1000 (horizon 2000): first run 4102 and 6436, rerun on the fixed code 2073 and 2436; pilot 1547 to 4764. The transform dominates | same, arm win, and PILOT lines; owner note "4.1 to 6.6 microseconds" |
| frequency E1 Double64 | measured 9776 to 10230 (first run), pilot ratio to Float64 2.1 | same, arm win_d, PILOT |
| phase advance E2 Float64 / Double64 | 524 to 966 / 2259 to 4150 | same, arms pha, pha_d |
| weighted Birkhoff Float64 / Double64 | 397 to 668 / 2176 to 4994 | same, arms win, win_d |
| Float64 GPU, AC13 ensemble | 2.468 ns per orbit-step (the LogNorm record, not the recording chain) | `environment.md` |
| device per spec, 1024 repeated orbits, wall ns per orbit-step against 23 CPU threads | tangent 69 against 37; pair 297 against 696; phase advance 30.6 against 50.7; birkhoff 62 against 44 (takes vary: phase advance 0.6 to 1.1 of CPU, tangent 1.4 to 3.1 times slower) | verification report, owner notes |

Line for E1 at this study's windows, derived: the worst measured cost per transformed sample is 6436 ns x 2000 steps / 2500 samples = 5.15 microseconds (samples per orbit = 2 x (250 + 1000)); the factor 1.3 allows the larger transform sizes (padded size 8 W, log factor 13 at W = 1000 against 17 to 19 at W = 25600 to 51200), giving 6.7 microseconds per sample. One orbit transforms 2 x sum(W) = 2 x 102300 = 204600 samples, so 204600 x 6.7 microseconds = 1.371 s thread time, which is 1.371 s / 102400 steps = 13.4 microseconds per orbit-step (the line of fq1 and its sensitivity specs). Double64 line: 13.4 x 2.1 = 28.1. Safety factors (proposed): 2 for every Float64 line (the lines are single records at short horizons), 3 for Double64 (as P1).

### Orbit-step counts (steps summed over orbits and items; the item count of an arm is orbits x methods and specs x sweep values)

| Arm | Orbits | Items per orbit | Horizon | Steps | Line (thread-s, before safety) |
|-----|--------|-----------------|---------|-------|--------------------------------|
| lyap_p | 883 | 2 (L1, tg) | 1e7 | 1.766e10 | 10663 |
| lyap_g | 2304 | 3 (L1, L5, tg) | 1e7 | 6.912e10 | 32606 |
| stretch_p | 883 | 1 | 102400 | 9.04e7 | 99 |
| stretch_g | 2304 | 1 | 102400 | 2.36e8 | 260 |
| pair_l | 92 | 1 | 102400 | 9.42e6 | 122 |
| pair_g | 2304 | 1 | 102400 | 2.36e8 | 3067 |
| win_p | 882 (Y5 omitted) | 2 (fq1, bk1) | 102400 | 1.81e8 | 1301 |
| win_g | 2304 | 2 | 102400 | 4.72e8 | 3397 |
| phase_p | 883 | 1 | 102401 | 9.04e7 | 90 |
| phase_g | 2304 | 1 | 102401 | 2.36e8 | 236 |
| dir_tg | 172 | 10 (L1, tg, 5 directions) | 1e6 | 1.72e9 | 1039 |
| dir_st | 172 | 5 | 102400 | 8.81e7 | 97 |
| sens_win | 92 | 2 | 102400 | 1.88e7 | 252 |
| sens_phase | 92 | 1 | 102401 | 9.42e6 | 9 |
| Float64 total | | | | 9.017e10 | 53239 |
| d64_lyap | 92 | 2 | 1e6 | 1.84e8 | 1162 |
| d64_win | 92 | 2 | 102400 | 1.88e7 | 312 |
| d64_phase | 92 | 1 | 102401 | 9.42e6 | 40 |
| Double64 total | | | | 2.123e8 | 1514 |

No two arms hold the same item (they differ in subset, spec, horizon or precision), so the runner's merge removes nothing. Each Line entry is steps x the line of the spec in ns (L1 and L5 207.6 or 6634.6; tg 1000 or 6000; pr 13000; st 1100; fq1 and its sensitivities 13400 or 28100; pa1 and pa_planar 1000 or 4200; bk1 1000 or 5000). Example, lyap_g: 2.304e10 x (207.6 + 207.6 + 1000) ns = 2.304e10 x 1415.2 ns = 32606 s.

### CPU-hour arithmetic (thread-hours)

- Float64: 53239 s = 14.79 h; x 2 = 29.6 h.
- Double64: 1514 s = 0.42 h; x 3 = 1.26 h.
- Total with safety: 29.6 + 1.3 = 30.8 thread-hours (all on the CPU).
- Wall time on 23 threads: 30.8 / 23 = 1.3 h at ideal scaling; 2.7 h at 50 percent efficiency (proposed pessimistic bound; the P1 campaign took 3.3 h wall, results Cost). Both are below the proposed 12 h cap and use at most the 23 cores of `environment.md`.

### GPU-hour arithmetic (device-eligible Float64 specs only)

`auto` follows the per-spec record, so the CPU total above is the budget and the GPU can only lower it. Upper bounds if every device-eligible spec ran on the device at the measured per-spec wall rate: tangent chain 3.27e10 x 69 ns = 2258 s = 0.63 GPU-h (more than its 0.34 h on 23 CPU threads at 37 ns, so `auto` is expected to keep it on the CPU); pair 2.45e8 x 297 ns = 73 s = 0.02 GPU-h; phase advance 3.26e8 x 30.6 ns = 10 s = 0.003 GPU-h; birkhoff 3.26e8 x 62 ns = 20 s = 0.006 GPU-h. With the AC13 rate of 2.468 ns per orbit-step on all 9.02e10 Float64 steps the floor is 0.06 GPU-h (0.12 with safety 2); that rate belongs to the LogNorm kernel and is not used for the decision. STR and E1 are CPU-only by D86. Double64 is CPU-only.

### Memory and storage (derived from spec E28; the run reads the real free memory)

- fq1, Float64: bytes(orbit) = M 80 + R 24 + S 8 x 7 columns x about 163 grid steps (9128) + B 32 x 51200 (1638400) = 1.65 MB; at f = 0.5 and 16 GB free the capacity is about 4850 orbits per launch. Double64: B = 64 W_max = 3.28 MB.
- bk1: reducer bytes 1776 plus series 9 kB; tg: bytes(orbit) 112 + 392 + 27880 = about 28 kB at 205 grid steps. All capacities exceed one chunk except the E1 launches of the grid chunk (2304 orbits), which fit one launch at 16 GB free and split otherwise; results do not depend on the capacity (I24).
- Store (orientation, derived): tg series 17 columns x 8 B x 205 steps = 28 kB per item, about 4100 items, so about 0.12 GB; L1, L5 and the tables add about 0.1 GB; total well below 1 GB (the P1 store is 475 MB).

### K19 and K20 for this study

- K19 (backend): before the first run the per-spec device records (tg, pr, pa1, bk1) are taken at 23 threads on the recording chain; the earlier record file is kept as `prior_record` (D58). Verification found the records vary between takes, so the real study takes its own; the CPU total of this section is the budget either way.
- K20 (Double64): CPU only, sized above (1.26 thread-hours with safety); S_L is 92 orbits. Changing S_L usage or the Double64 horizons is an amendment.
- Audit line "Budget" of the stage contract: 30.8 thread-hours with safety, 1.3 to 2.7 wall hours, core limit 23 respected.

## Analysis plan

Every figure and table, and the hypothesis each decides. All plots and tables use the stored series and tables; sets are those of Panel.

| # | Item | Content | Decides |
|---|------|---------|---------|
| T1 | Panel and grid table | orbits by source and K, counts of rho (C, R, U), drops (controls, drift), sizes of E_all, E8, W2, S2, T_s, T_s^design, orbits omitted per spec | none |
| T2 | Verdict table | per indicator, budget and set (W2, S2, T_s, regular): counts chaotic, regular, undecided | context |
| T3 | H1 table | U(b_P) with interval, A_min, the 21 values A_XY(b_P), U on S2 | H1 |
| F1 | Agreement matrices | A_XY and kappa_XY for I_8 at b_P and N_c, sets W2 and S2 | H1 context |
| T4 | H2 table | M_X(S2), outcome shares, Q with interval | H2 |
| F2 | Detection-time distributions | t_det per indicator for S2, W2, T_s | H2, H3 |
| T5 | H3 table | M_X(W2), M_X(S2), R_X with intervals, number of indicators that pass | H3 |
| T6 | H4 table | phi_X, Delta_phi with interval, for T_s and T_s^design | H4 |
| F3 | Verdict against budget | fraction "chaotic" per indicator and budget for W2, S2, T_s, regular set | H4, H5 context |
| T7 | H5 table | F_X, F_med with interval | H5 |
| T8 | H6 table | delta_dir(X), delta_prec(X), set sizes | H6 |
| T9, F4 | AUC | AUC(X, b, K) for I_8 and STR | reported |
| T10 | Fixed-threshold rules | the statistics of H1 to H6 under table R2 | reported (robustness) |
| T11 | Sensitivities | fq1_raw and fq1_qp against fq1, pa_planar against pa1, WBA with cos q alone against dig_min, on S_L: verdict agreement and difference of D | reported |
| T12 | Rule-validity gate | share of K = 0 gate orbits with a chaotic verdict at N_c, per indicator | gate (below) |
| T13 | I_8 extension | H1 to H6 on E8 with I_8, ROBUST flags | extension |
| T14 | Regular range | M10 for the sticky orbits | reported |
| T15 | Packages | FLI+RLI and MEGNO+SALI (research idea 17) against rho at b_P and N_c | reported |
| T16 | Flags and controls | flag counts per spec, control values (K = 0 line, fixed points, Y5), omissions | reported |
| T17 | Cost | thread-ns per orbit-step per spec against its line, bytes(orbit), capacity | reported |
| F5 | Grid maps | for each K, the 24 x 12 lattice coloured by rho and by the verdict of each indicator at b_P and N_c | reported |
| F6 | D_str against n | STR by class | reported |

Rules for the analysis:
- Hypotheses are decided in the order H1 to H6 with the numbers above and nothing else. The primary rules are table R1 on I_7 and E_all. No threshold is tuned after a look at panel or grid data.
- Rule-validity gate (proposed): for each X in I_8 the share of K = 0 gate orbits (80 orbits, regular by construction) with c_X(N_c) = chaotic must be at most 0.10. A violation is reported as "rule invalid on the integrable control" with the indicator named; the decisions stay as preregistered and a second version of each decision table without that indicator is added (never in place of the first).
- Any analysis not listed here (extra figure, different subset, other budget, other rule) is labelled "exploratory".
- A set with fewer than 20 orbits is "not decidable". Results are reported by K and source (panel, grid) as well as pooled.
- The report cites only keys present in `/_config/bibliography/references.bib` and never an item of `unverified.md`. The keys used in this file: Maffione2011a, Darriba2012, HarsoulaKaramanosContopoulos2019, Hinse2010, Das2017, Meiss2015, FrahmShepelyansky2010. The indicator keys the report will add, all checked present: SandorMaffione2016, Skokos2001, Skokos2007, Skokos2010, Froeschle1997, CincottaSimo2000, Mestre2011, Bazzani2023, Laskar2003, SanderMeiss2020, DuignanMeiss2023. No threshold of this study is taken from the literature: `indicators.md` and `benchmarks.md` hold no numeric decision threshold for these indicators, so every number of tables R1 and R2 is proposed by this study, placed between the regular and the chaotic levels that the catalog rows record.

## Stopping rules

Unit of work = chunk: the items of one run type in one panel cell, or all controls, or all `h4` starts (runner R3); resumable (K8, K13).

A chunk fails if any of (runner R5 and R18):
- a non-finite value in a recorded series, or a required value of a window table without its flag (spec E31), `complete` false, or the frequency status `failed`;
- a numerical health error: zero or non-finite column norm, pairing failure, or a pair separation outside [Dx / 2, 2 Dx] (the study chooses wrapped coordinates and Dx = 1e-12 so that this cannot happen);
- the determinant or symplecticity monitor above the tolerance multipliers of `numerics-defaults.toml` in more than 5 percent of the chunk's orbits (the exceeding orbits are flagged in any case);
- an exact-exponent control of arm `lyap_p` (L1) beyond 1e-8 (scope `exact_start`): a code defect, the campaign stops at once;
- the cost of the first completed launch group of an arm above 2 times its line (the spec budget times safety): the arm is paused and resumes only by a recorded cap override that changes no threshold, horizon or subset (P1 Amendment A4 item 2);
- a Double64 item that disagrees with its Float64 item on the orbit coordinates by more than 1e-10 over the first 100 steps (classes T and R, elliptic controls) or the first 5 steps (the other classes, hyperbolic controls, grid starts).

Handling:
- A failed chunk is rerun once with the same inputs on the CPU. A deterministic failure repeats and is final; its orbits are "excluded, failed twice" and leave every set. The campaign stops if more than 5 percent of all chunks have failed twice (proposed, the P1 value) or the wall cap of 12 hours (proposed) is reached.
- Flags are data and never failures: `below_floor`, `degenerate` (a constant series, such as the fixed-point controls), `undefined window`, `near_zero`, `mirror`, `tie`, `near_pi`, `zero`, `reldig_undefined`. The module stops no orbit early and draws no conclusion from a flag (spec I25).
- The omission of orbit 883 from the E1 and WBA specs is planned, listed in the plan and covered by the plan hash; it is not a failure.
- There is no early stopping on scientific results: no hypothesis is decided from partial data, and the rule-validity gate is not a stop. The only early-stop causes are the failures above and the wall cap.
- Interruptions resume from the last stored chunk; the manifest records versions, seed, backend per spec, threads, the device records of K19 and the per-spec chain records.

## Amendment rule

Frozen after study/01. A change after freeze is an amendment: one numbered entry in the header's amendment log, with date, what changed, why and the git SHA of the commit that makes it; allowed only before stage 03-run produces any data file for a panel or grid orbit. Two cases are planned and fixed here:
- Amendment A1 appends the generated `[[run.h4]]` tables of Panel, "Grid G" (a mechanical expansion of the rule in that section; it changes no value of this file).
- Pilot calibration: the pilot of experiments spec R20 uses non-panel starts, writes outside the study store and is not data in this sense (P1 Amendment A4 item 2). On its evidence the owner may amend the numbers of tables R1 and R2 and the budget lines before any panel or grid data exist. After data exists, a changed threshold is not an amendment: it becomes an exploratory analysis and the preregistered decision stands.

## Open points for the owner at the checkpoint

Each value marked "proposed" above needs approval. The list, with the rationale kept in the section named:

1. Decision sets and indicator sets (Hypotheses, Panel): decisions use I_7 on E_all (no RLI); I_8 on E8 is the extension. Reason: RLI covers only S_L and G.
2. Budget list B = 200 x 2^k up to N_c = 102400, primary budget b_P = 12800, persistence rule b_k <= N_c / 4 (Metrics M1, M4).
3. Class-rule numbers of table R1: growth ratio 5, LE floor 1e-4, RLI 1e-7, E1 and E2 D >= -4, WBA 6 digits (minimum over four observables). Table R2: SALI 1e-8, FLI and LE 3 ln n, MEGNO 2.5. Calibratable on the pilot.
4. Censoring value 2 N_c for wrong and undetected verdicts (M4).
5. H1: U >= 0.20, lower bound >= 0.10, A_min <= 0.80.
6. H2: Q >= 10, lower bound > 3.
7. H3: at least 5 of 7 indicators with R_X >= 10 and lower bound > 3.
8. H4: Delta_phi >= 0.50, lower bound >= 0.25.
9. H5: F_med >= 0.10, lower bound >= 0.05.
10. H6: 0.05 for direction and precision flips.
11. Set parameters: lambda_ref = lambda_1(1e7); S2 at lambda_ref >= 1.0; W2 at lambda_ref <= 0.1 or sticky; sticky R4 with n_s = 1e5 and f_s = 0.2; percentile pair 5 and 95 of M10.
12. Grid G: eight K values, 24 x 12 lattice with p0 in (0, pi), 2304 starts, labels, carried as `[[run.h4]]` by Amendment A1.
13. Reference class of grid orbits: P1 label rules on L5 at 1e7 with the L1 against L5 test, flag ref_weak.
14. Horizons: 1e7 for the lyap arms, 102400 and 102401 for the pair, stretching and window arms, 1e6 for the direction and Double64 arms; ten windows 100 to 51200.
15. mean_removal = true for fq1, with the fq1_raw sensitivity; E1 keys order 1, positions `half`, pad 8, mirror_ratio 0.9, floor_multiplier 32 (naff-fma P10 drafting values).
16. E2: observable `lifted`, birkhoff_order 1, planar sensitivity.
17. WBA: four observables (cos q, sin q, cos p, sin p), birkhoff_order 1.
18. STR: range `derived`, bins 1000, interval 100, min_per_cell 100, counts `S_D`.
19. Pair: separation 1e-12, frame_source `seeded`, members RLI and MEGNO, subsets S_L and h4 (Y5 excluded by subset), horizon 102400.
20. Directions (the five of P1) and the Double64 stability arm on S_L.
21. Rule-validity gate: 0.10 on the K = 0 orbits.
22. Budget lines, safety factors (2 and 3), wall cap 12 h, first-chunk cost factor 2, pilot of 4 orbits per K at horizon 1e5.
23. min_set = 20 and the bootstrap settings (1e4 resamples, seed 20261012).
24. Pilot calibration provision (Amendment rule).

Not proposals, unresolved and needing an owner decision or a runner change:
- Owner decision at the checkpoint (2026-10-29): U1 is a separate follow-up (a runner amendment), not part of this freeze; the study stays as drafted, with the grid as `[[run.h4]]` tables (Amendment A1).
- U1. The runner has no grid subset and no subset "panel without controls". The grid goes in as 2304 `[[run.h4]]` tables (one chunk, so one failing orbit sends 2304 orbits to the rerun), and RLI reaches only 92 of the 883 panel orbits. A runner amendment (a `grid` subset built from a rule; a subset without controls) would remove both limits.
- U2. The per-orbit P1 labels are not a tracked file (tables hold counts, T6). Stage 04 must read them from the P1 store through the Analysis module; if the store is absent, the labels must be regenerated first.
- U3. The size of T_s is unknown before data (a data-dependent subset); below 20 orbits H4 is not decidable.
- U4. The run block of this file has been parsed for syntax only; the runner has not read it.

## Checks I ran

- Catalog `Status` read (2026-10-29): indicators.md rows I1, I2, I5, I8, I9, I10, I12, I12b, I13 verified; I3, I4, I6, I7, I11, I12c, I12d, I14 to I18 planned. lyapunov.md rows L1, L5, L13 verified (L13 only for the distribution, mean and second moment of one orbit's stretching numbers; L14 planned). systems.md rows Y1 and Y5 verified, Y2 planned. Sections E and F of lyapunov.md verified.
- Bibliography keys grepped in `references.bib`, all present and none named as an item in `unverified.md`: Maffione2011a, Darriba2012, HarsoulaKaramanosContopoulos2019, Hinse2010, Das2017, Meiss2015, FrahmShepelyansky2010 (cited above), and SandorMaffione2016, Skokos2001, Skokos2007, Skokos2010, Froeschle1997, CincottaSimo2000, Mestre2011, Bazzani2023, Laskar2003, SanderMeiss2020, DuignanMeiss2023 (for the report). The names Skokos2007 and Bazzani2023 occur inside notes of `unverified.md` about other items; they are not listed there as unverified keys.
- Facts re-derived from sources: panel counts (883 orbits, 13 controls, S_D 172, S_L 92, no control in S_L, design class S 110, Y5 = orbit 883 starting at the origin) by parsing the P1 panel file; the arm and item arithmetic by script; the run block by a TOML parser (syntax only).
- Sources not found or not read in full, and the gaps left open: the individual indicator specs (sali, fli-ofli, megno, rli, stretching-numbers, naff-fma, weighted-birkhoff) were read only through targeted excerpts; `schema/run-store.md` was not read (its content is taken through spec requirement R15); `_config/numerics-policy.md` and `_config/architecture.md` were not read; the P1 run store (`data/`) was not opened; no per-orbit P1 label file exists among the tracked outputs (U2); `indicators.md` and `benchmarks.md` of the bibliography hold no numeric threshold for these indicators (all thresholds are proposed here); the cost of the pair spec at long horizons has not been measured (only N = 300), so 13000 ns is a bound; no pilot was run, so all table R1 and R2 numbers are untested on this system.
- Places where the sources contradict the task summary or each other: (a) the summary says pair is CPU only; the spec (R16) makes `pair` device-eligible per record, and the verification report measures the pair on the device at 0.43 of the CPU time; only STR and E1 are CPU-only. (b) The summary quotes 4.1 to 6.6 microseconds per orbit-step for the frequency kind; that is the first verification run (4102 and 6436 ns); the rerun on the fixed code gives 2073 and 2436 ns (arm) and 1547 ns (pilot), so the budget uses the larger values. (c) The summary calls L1 the P1 default; the catalog marks L5 as the default reference of the QR family; this file uses L1 for the indicator LE and L5 only for the grid labels, with the reasons in Methods. (d) The Float64 GPU rate 2.47 ns per orbit-step of `environment.md` is the LogNorm record; D54 and K19 replace it by per-spec records on the recording chain (tangent 69 ns against 37 on 23 CPU threads), so the GPU estimate gives both. (e) `src/Indicators/CONTEXT.md` quotes 1.4 kB for weighted Birkhoff at (K, H) = (4, 11), which matches neither order of the formula 16 + 8 H (2 + 5 K) (spec D86 "Not sourced", item a); the formula is used. (f) Verification criterion D86y mentions an empty study data folder; the P1 store has since filled it (no consequence here, no file under `data/` was touched).
- Not run: the ICM validator (`scripts/icm_validate.py`), any Julia code, the runner's plan or pilot. Nothing was frozen, committed or pushed; this file is the only file written.
