# Preregistration: p1-le-benchmark

FROZEN 2026-10-17. Every value marked "proposed" was approved by the owner at the study/01 checkpoint; the one-line rationales stay as the record. All thresholds live only in this file.

## Header

| Field | Value |
|-------|-------|
| Study slug | p1-le-benchmark |
| Paper | P1 (`/_config/research-question.md`, paper sequence) |
| Claims tested | C1, C2. C3 is out of scope (P2). |
| Value types | FALSIFIES (C1, C2), REPRODUCES (SalaArtuso2015), TOOL (Chaoslab v0.1.0 at M3 exit), NEGATIVE allowed (`/_config/value-types.md`) |
| Git SHA at freeze | 1e760e229eabf9d8323ad4d8097cc1dbd017f845 (parent commit of the freeze commit) |
| SHA-256 of this file (current version) | 650b3558ab9a728747da85b0f46b123be888a6224996baa667543e6e2e1ddf46 (rule: Amendment A1, item 1) |
| Amendment log | A1 (2026-10-17, section "Amendment A1"; parent commit 963345d86570366f5f0ea3ee6950288e20edcfa7; replaces freeze hash fa84e883405a7acbdb4774be2527bbecb8399a4c369fa79fbfc4d563d98bea87); A2 (2026-10-17, section "Amendment A2"; parent commit ceb5227dfcbd258d4f72dfc4ff99602258320b36; replaces hash c9ad3c2cf3ba1220990d118c81e55da2104a7130b0f49e7dcd83cf75461fdcef); A3 (2026-10-17, section "Amendment A3"; parent commit ceb5227dfcbd258d4f72dfc4ff99602258320b36; replaces hash d89d2bf99e4eaa5740be88da436ec28b2d92f7859a278c50fe734ff4a02ae69c) |

Value types, in the words of `value-types.md`:
- FALSIFIES: H1 to H3 each have a decision rule below, fixed before any data exists.
- REPRODUCES: H4 against the benchmark rows with source key SalaArtuso2015.
- TOOL: the run uses only modules whose catalog `Status` is `verified`; the tagged release v0.1.0 is cut at M3 exit.
- NEGATIVE: if H2 or H3 fail, the null is reported with the same figures and tables as a pass.

## Scope check against the catalogs

Read from `/_config/catalogs/lyapunov.md` and `/_config/catalogs/systems.md`.

| ID | Role here | Catalog Status | In scope |
|----|-----------|----------------|----------|
| L1 Benettin | method | verified | yes |
| L2 two-trajectory | method (approximate lambda_1) | verified | yes |
| L3 CGS | method | verified | yes |
| L4 MGS | method | verified | yes |
| L4b RGS | method | verified | yes |
| L5 Householder | method, default reference | verified | yes |
| L5b H2 (diagonal iteration path) | method | verified (H2; H1 deferred) | yes, H2 only |
| L5b H1 (treppen iteration path) | method | deferred | EXCLUDED: not built, so not `verified` |
| L6 Givens | method | verified | yes |
| L9 stabilized SVD | method | verified (d = 2) | yes, d = 2 only, which is all Y1 needs |
| L7 fixed reorthonormalization interval tau | parameter | verified | yes, as a parameter |
| Y1 standard map | system | verified | yes |
| Y5 linear symplectic maps | exact-exponent control | verified | yes |

L7b, L7c, L8, L10 to L18 stay out of scope. Precision backends Float64, Double64, BigFloat: catalog section E, `verified`. Diagnostics (pairing, det, orthogonality, sum residuals): catalog section F, `verified`.

Method sets used below:
- QR family Q = {L3, L4, L4b, L5, L5b(H2), L6}.
- Variational set V = {L1} + Q + {L9}.
- Approximate method: L2.

## Hypotheses

Notation and metric definitions are in "Metrics". Eligible sets are defined in "Panel".

### H1 (claim C1, asymptotic equivalence of lambda_1)

H1 is decided by three sub-rules on the chaotic-eligible set E_ch (orbits with design class chaotic sea, plus the exact-exponent controls). Horizon N_1 = 1e7 steps (proposed: it is the order of 1e7 iterations that HarsoulaKaramanosContopoulos2019 report for a large chaotic-sea exponent to converge, see `benchmarks.md`).

- H1a (same tangent stream, same tau): metric Delta_QB(N_1) = max over {L1} + Q of lambda_m(N_1), minus the min. Orbit passes if Delta_QB(N_1) <= 1e-8 (proposed: about eight orders above Float64 rounding growth over 1e7 steps, so a pass is meaningful and CGS orthogonality loss is allowed for).
- H1b (SVD against QR): metric |lambda_L9(N_1) - lambda_L5(N_1)|. Orbit passes if <= 1e-6 (proposed: 10 times 1/N_1, the order of the finite-time transient between product-SVD values and QR diagonals, see the L9 row in the catalog).
- H1c (two-trajectory): metric |lambda_L2(N_1) - lambda_L5(N_1)| at the primary perturbation size eps_0 = 1e-8 (proposed: between nonlinear saturation and round-off). Orbit passes if <= 5e-3 (proposed: L2 is only an approximation with a bias set by eps).

Decision: H1 is SUPPORTED if, for each sub-rule, the fraction of orbits in E_ch that pass is >= 0.98 for H1a and H1b and >= 0.90 for H1c (proposed: allows rare marginal orbits near islands; the failing orbits are listed). H1 is FALSIFIED if any sub-rule fraction is below its limit. Exact-exponent controls (Y5 hyperbolic linear map and the Y1 hyperbolic fixed points) must satisfy |lambda_m(N_1) - lambda_exact| <= 1e-8 for every method in V (proposed); a violation is a defect report, not a scientific result, and blocks the analysis (see Stopping rules).

### H2 (claim C2, finite-time inequivalence)

Metric: relative finite-time spread R(n) = Delta_V(n) / lambda_ref, where Delta_V(n) = max over V of lambda_m(n) minus the min, and lambda_ref is the orbit's lambda_L5 at its longest horizon (N_1, or N_long for the long subset). Evaluate at n_H2 = 1e4 steps (proposed: a realistic budget, three orders below N_1).

The weak-chaos set W = orbits with instability label "positive growth supported" (D37) in the design classes sticky or near-separatrix at any K > 0, plus chaotic-sea orbits at K <= K_c. Decision: H2 is SUPPORTED if the median of R(n_H2) over W is >= 0.10 (proposed: a ten percent disagreement is substantial for a quantity meant to be a single number), and the lower end of its 95 percent bootstrap interval (10^4 resamples of orbits, seed from `numerics-defaults.toml`) is >= 0.05 (proposed). Otherwise H2 is not supported; the null is reported (NEGATIVE).

### H3 (claim C2, regime dependence of the convergence time)

Metric: convergence time n_conv,m = smallest n such that |lambda_m(k) - lambda_ref| <= 0.05 lambda_ref for all recorded k in [n, N] (proposed 0.05: a five percent band). Per orbit, n_conv = max over V. Strong-chaos set S = chaotic-sea orbits at K = 2 pi and K = 10. Decision: H3 is SUPPORTED if median n_conv over W divided by median n_conv over S is >= 10 (proposed: one decade, matching the "regime-dependent" headline in `research-question.md`) and the lower end of the 95 percent bootstrap interval of the ratio is > 3 (proposed). Orbits whose n_conv is undefined (never inside the band by N) are given n_conv = N and counted; their share is reported. Otherwise H3 is not supported (NEGATIVE).

### H4 (REPRODUCES, SalaArtuso2015)

Orbit and K from `benchmarks.md`, rows with source key SalaArtuso2015: McMillan form, start (x, y) = (1, 2) x 10^-3, T = 1e9 iterations. Converted to the project form p = x - y (as stated in the benchmark row), the start is x0 = 1e-3, p0 = -1e-3.

- K = pi/2: reference 0.298. Pass if L1 and L5 (Float64, tau = 1) both give lambda_1(1e9) within 0.01 of the reference (proposed: the paper gives three digits, a single orbit, no error bar and unstated precision, `benchmarks.md`, K18).
- K = 2 pi: reference 1.172. Pass if L1 and L5 both within 0.02 (proposed: same, plus small accelerator-mode islands at this K).
- Cross-check: Double64, L5, same starts, T = 1e8 (proposed: cost, see Budget). Reported with the Float64 value at the same T; no pass or fail. Chaotic pseudo-orbits differ between precisions, so only the estimator is compared (K16).

H4 is SUPPORTED if both K pass. A failure is reported as "did not reproduce" with the numbers, together with the possibility that the paper's orbit sits in a different component (the benchmark row notes a mixed phase space). Shevchenko2014 fits (rows with that source key) are loose comparators only: reported against the chaotic-sea ensemble mean of lambda_L5(N_1) at each K, never pass or fail. The Chirikov ln(K/2) line (Chirikov1979, via `benchmarks.md`) is reported the same way at K = pi/2, 2 pi and 10.

### H5 (D36 and D37 revisit condition: are the labels informative?)

D36 and D37 name their own failure conditions in `decision-log.md`: the unresolved fraction is too large to be informative, or the structure axis is not decidable on maps without frequency analysis. Decided on the ladder subset S_L (defined in Panel):

- H5a (stability across precision): fraction of S_L orbits whose two labels are identical at Float64, Double64 and BigFloat is >= 0.90 (proposed).
- H5b (informativeness): the fraction of all panel orbits with either axis "unresolved" is <= 0.30 (proposed: a classifier that abstains more than 30 percent of the time is not a useful baseline).
- H5c (structure axis decidable): among orbits whose instability label is "none detected on the tested horizon", the fraction with structure "unresolved" is <= 0.30 (proposed).

If H5b or H5c fails, the label scheme is reported as not informative and decisions D36 or D37 go back to the owner. This is a design result, not a failure of the study.

### Reported without a decision rule

Cost per method (wall time, tangent-vector count, memory, precision), the symplectic pairing residual eps_pair(n), the area-preservation residual eps_det(n), and sensitivities (tau, perturbation size, direction, wrapped vs lifted) are reported as tables and figures. They are used as gates in "Stopping rules" but are not hypotheses.

## Panel

Orbit classes, K regimes and the output list come from `/_config/research-question.md` ("Test panel", "Outputs retained"). The Y1 map definition and tangent map live there too and are not repeated. Raw initial conditions are built by stage 02-panel; this section fixes the selection rule.

### K values

| Regime | K | Note |
|--------|---|------|
| integrable | 0 | exponent exactly 0 |
| small K | 0.5 (proposed: KAM tori plus resonances, well below K_c) | |
| near-critical | K_c = 0.971635406 | source key FrahmShepelyansky2010 (also CarlettiLaskar2000), `benchmarks.md` |
| intermediate | pi/2 | mixed phase space; matches the SalaArtuso2015 row |
| large | 2 pi and 10 (proposed 10: Chirikov regime, ln(K/2) line) | 2 pi matches the SalaArtuso2015 row; small accelerator-mode islands exist there |

### Cells

Cell = K value x design class. Design classes (from `research-question.md`): KAM torus (T), resonant (near-)periodic (R), chaotic sea (C), sticky (S), near separatrix or hyperbolic fixed point (X). At K = 0 only T and R exist; cells that stage 02 cannot populate after its search budget are recorded as empty, never filled from another class. Maximum 27 populated cells (6 K x 5 classes, minus C, S, X at K = 0). Orbits per cell: n_c = 40 (proposed: gives a bootstrap interval for medians per cell and keeps the panel near 1100 orbits).

Design class is the a-priori intended class. All hypotheses and the long, sensitivity and ladder subsets use it for stratification; the measured labels (D36, D37) are separate outputs. Reporting both avoids selecting orbits by the very estimators under test.

### Selection rule (what stage 02 builds)

1. Screening grid: 256 x 256 uniform points of the torus per K, run with L5, Float64, N_scr = 1e5 steps (proposed). The screening points and results are disjoint from the panel and are discarded after assignment; no screening value enters any metric.
2. Predicates on screening data (all proposed):
   - T: lambda_scr <= 1e-4, and at least 0.05 (torus distance) from every periodic point in item 3.
   - C: lambda_scr >= 0.01 and at least 0.3 from every elliptic periodic point in item 3.
   - S: lambda_scr >= 1e-3 and a contiguous stretch of at least 2e4 steps in which the 1e3-step windowed estimate stays below 0.2 times the orbit's median windowed value ("stuck window").
   - X: not from screening. Start at torus distance r from the hyperbolic fixed point (0, 0), with log10 r uniform in [-4, -2] and the angle uniform. This includes starts like the SalaArtuso2015 start (distance about 1e-3).
3. Periodic points: for K > 0, fixed points (0, 0) and (pi, 0) are known exactly (hyperbolic at (0, 0)). Further periodic points of period <= 6 come from Newton's method in BigFloat, accepted when the return residual is <= 1e-12 (proposed). Class R starts at an elliptic periodic point displaced by 1e-3 (proposed), uniform in angle. The exact periodic points themselves are kept as control orbits (structure label "periodic" by construction; their exponent is known analytically from the eigenvalues of the analytic tangent map). At K = 0, class R uses exact rational p0 = 2 pi a / b with b <= 6.
4. Within each predicate set, ICs are drawn uniformly from the qualifying screening cells with a stream derived from the default seed (`numerics-defaults.toml`) and the cell id; the draw order and counts are recorded in the panel manifest.
5. Controls: the exact periodic points of item 3 (at most 12) and one Y5 hyperbolic linear map with eigenvalue pair (mu, 1/mu), mu fixed in stage 02 (proposed: mu = (3 + sqrt(5))/2); exponent ln mu.

### Subsets (nested)

- Full panel P: at most 27 x 40 + 13 controls = 1093 orbits; budget bound N_P = 1100.
- Designated subset S_D: 7 orbits per populated cell plus all controls, capped at 200. Used for the long horizon, tau sweep, perturbation-size sweep, wrapped vs lifted comparison and the multi-direction re-runs.
- Ladder subset S_L: 4 orbits per populated cell from S_D, capped at 100. Used for the precision ladder.

### Eligible sets

- E_ch: design class C (K > 0) plus the controls with positive analytic exponent (hyperbolic fixed points, Y5).
- W and S: defined in H2 and H3, formed from measured labels within design classes.
- Orbits flagged "drift" (lifted |p| exceeds 1e6 at any recorded time, proposed: accelerator-mode growth that degrades the lifted-coordinate arguments of sin and cos) are removed from E_ch, W and S and listed in a separate table.

## Methods

All methods are in the table "Scope check" and are `verified`. Method details, tolerances and algorithm identifiers live in the catalog and the module specs; they are not copied here.

Defaults from `/_config/numerics-defaults.toml` used unchanged: `seed`, `bigfloat_bits`, `orthogonality_tol`, `pairing_residual_tol`, `det_residual_tol`. Overrides:

| Key | Default | This study | Reason |
|-----|---------|------------|--------|
| precision | Float64 | Float64 for the panel; Double64 and BigFloat on S_L and for H4 cross-check | precision ladder (D36, D37); K20 |
| backend | auto | auto for Float64, but resolved by the K19 rerun (see Budget); Double64 and BigFloat always CPU | Double64 does not compile for the device (K20) |
| renorm_interval_tau | 1 | primary tau = 1 (the default); sweep tau in {1, 2, 4, 8, 16} (proposed) on S_D for L4 and L5; sweep stays within tau_max(T) from the tangent-core spec | L7 study |

Other configuration:
- Coordinates: lifted internally (Y1 default); wrapped is a sensitivity factor on S_D (all of V plus L2).
- Perturbation size for L2: primary eps_0 = 1e-8; sweep eps in {1e-6, 1e-8, 1e-10, 1e-12} on S_D (proposed).
- Initial deviation direction: primary angle pi/4 from the x axis in the (x, p) plane (proposed: not aligned with an axis). On S_D, four more angles 0.3, 0.3 + pi/5, 0.3 + 2 pi/5, 0.3 + 3 pi/5 (proposed), giving 5 directions. Methods re-run per direction: L1, L5, L9. The Q-method initial frame is rotated to the same direction.
- Initial deviation magnitude: unit vector for all variational methods.
- Backend choice does not change results beyond the tolerance because estimators are compared, never trajectories (K16); CPU is the reference. Thread count does not change results: stage 04 hashes were equal at 1 and 4 threads for all in-scope methods (`threads-verify-log.txt`).

## Metrics

Exact definitions of every reported quantity. Where a symbol is shared with the repo canon it is used as in `research-question.md`.

- lambda_m(n): finite-time lambda_1 of method m after n map steps, natural log per step, norm = Euclidean 2-norm on the lifted (x, p) tangent vector (this is the "stated canonical norm" of D37; proposed). L9 uses the largest finite-time singular value of the product.
- Recorded n: n = 1..1e3 every step is too costly for storage; the series are stored at 20 points per decade log-spaced, plus the first 100 steps, up to the horizon (proposed; adequate for n_conv at the 5 percent band).
- Horizons: N_1 = 1e7 for all of P; N_long = 1e8 for S_D (proposed: sticky layers need more than 1e7, see `benchmarks.md` rows from HarsoulaKaramanosContopoulos2019 and Meiss2015); N_D = 1e7 and 1e8 and N_B = 1e6 for the ladder (see Budget).
- Delta_QB(n), Delta_V(n): max minus min of lambda_m(n) over the stated set.
- R(n): as in H2. n_conv: as in H3.
- lambda_ref per orbit: lambda_L5 at the orbit's longest Float64 horizon.
- Regular-orbit reference: Skokos2010 (see `benchmarks.md`) gives lambda_1(n) = O(ln n / n) for regular orbits; the instability test below uses this scaling.
- eps_pair(n), eps_det(n), orthogonality loss, sum residual: as in catalog section F and `research-question.md`; flagged against the tolerance multipliers in `numerics-defaults.toml`.
- Cost: wall time, tangent-vector count, memory, precision per orbit-method pair.
- Detection time: not reported in P1 (indicator-based; P2).

### Label rules (axes of D36 and D37)

The thresholds below make the two preregistered axes decidable. They apply to the final Float64 horizon of each orbit (N_1, or N_long).

Instability axis, using L5:
- "positive growth supported": lambda_L5(N) >= 1e-4 (proposed: about two orders above the regular O(ln N / N) floor at 1e7, see Skokos2010 row) and r(N) >= 5, where r(N) = N lambda(N) / [(N/10) lambda(N/10)] (proposed: r tends to about 1.2 for regular O(ln n / n) decay and to 10 for constant lambda).
- "none detected on the tested horizon": lambda_L5(N) < 1e-4 and r(N) < 5.
- "unresolved": all other cases, or any of the disagreement tests below.

Structure axis:
- "periodic": only for the exact periodic points of item 3 in Selection rule, confirmed by a Float64 return residual <= 1e-10 (proposed) after the known period.
- "quasiperiodic": instability "none detected" and the rotation-number estimates from the first and second half of the horizon, rho = (x_n - x_0)/(2 pi n) from the lifted x, agree within 1e-7 (proposed: a coarse test, the 1/n bias of the plain average is about 1e-8 at N_1).
- "irregular": instability "positive growth supported".
- "unresolved": all other cases, or any disagreement test below.

Disagreement tests that force "unresolved" on either axis (D36, D37 definitions):
- Precision ladder (on S_L only): |lambda_Float64 - lambda_Double64| or |lambda_Double64 - lambda_BigFloat| exceeds delta_lad = 5e-3 at N_B = 1e6 (proposed: shadowing and sampling scatter of a chaotic orbit at 1e6 steps is about 1e-3), or the labels differ between precisions.
- Estimators: any failure of H1a or H1b on that orbit.
- Direction (S_D only): the label changes between the 5 initial directions (any change is unresolved).
- Diagnostics: pairing, det or orthogonality residual above the tolerance multipliers of `numerics-defaults.toml` for L5 on that orbit.

Orbits outside S_L and S_D carry labels from estimator disagreement and diagnostics only; the label tables say which tests could be applied to each orbit.

## Budget

Environment limits (`/_config/environment.md`): this computer, max cores 23, GPU RTX 5060 Laptop (8 GB), Double64 backend on CPU. No wall-clock limit is stated there; a campaign cap of 12 wall hours is proposed (about 5 times the estimate below, rationale: leaves room for the repeats of failed chunks).

### Throughput inputs

| Quantity | Value | Source |
|----------|-------|--------|
| Float64 CPU, 4 threads | 8.255 ns per orbit-step (wall), i.e. 33.0 ns of thread time | `environment.md` and `data/gpu-benchmark.toml` record of 2026-10-17 |
| Float64 CPU, 1 thread | 33.09 ns per orbit-step | stage 04 tangent-core report, AC13 record |
| Float64 GPU, AC13 ensemble | 2.468 ns per orbit-step | `environment.md`, same record |
| Float64 GPU, one wave | 9.137 ns per orbit-step | same record |
| Double64 CPU, 4 threads | 388.4 ns per orbit-step (wall), i.e. 1553.8 ns thread time, about 47 times Float64 CPU | `data/gpu-benchmark.toml`, K20 |
| BigFloat(256) / Float64 CPU | about 550 | stage 04 precision report, AC11 (order of magnitude; taken on a loaded machine); thread time about 18200 ns per step |

Safety factors (proposed): 2 for Float64 (the record is one kernel; L2 uses two orbits, the Q methods use k = 2 vectors) and 3 for Double64 and BigFloat (loaded-machine ratios vary by tens of percent).

### Orbit-step counts (steps summed over orbits and methods)

| Arm | Orbits | Methods or runs | Horizon | Steps |
|-----|--------|-----------------|---------|-------|
| Main, Float64 | 1100 | 9 (L1, L2, L3, L4, L4b, L5, L5b H2, L6, L9) | 1e7 | 9.9e10 |
| Long, Float64 | 200 | 9 | 1e8 | 1.8e11 |
| tau sweep, Float64 | 200 | 5 values of tau x 2 methods (L4, L5) | 1e7 | 2.0e10 |
| Wrapped vs lifted, Float64 | 200 | 9 | 1e7 | 1.8e10 |
| Directions, Float64 | 200 | 5 directions x 3 methods (L1, L5, L9) | 1e7 | 3.0e10 |
| L2 eps sweep, Float64 | 200 | 4 values | 1e7 | 8.0e9 |
| H4, Float64 | 2 starts | 2 methods (L1, L5) | 1e9 | 4.0e9 |
| Float64 total | | | | 3.59e11 |
| Ladder Double64, short | 100 | 4 methods (L1, L3, L5, L9) | 1e7 | 4.0e9 |
| Ladder Double64, long | 100 | 1 method (L5) | 1e8 | 1.0e10 |
| H4 cross-check Double64 | 2 starts | L5 | 1e8 | 2.0e8 |
| Double64 total | | | | 1.42e10 |
| Ladder BigFloat | 100 | 3 methods (L1, L5, L9) | 1e6 | 3.0e8 |

### CPU-hour arithmetic (thread-hours)

- Float64: 3.59e11 steps x 33.09 ns = 1.188e13 ns = 11880 s = 3.30 h; x 2 = 6.6 h.
- Double64: 1.42e10 steps x 1553.8 ns = 2.206e13 ns = 22060 s = 6.13 h; x 3 = 18.4 h.
- BigFloat: 3.0e8 steps x 18200 ns = 5.46e12 ns = 5460 s = 1.52 h; x 3 = 4.6 h.
- Total with safety: 6.6 + 18.4 + 4.6 = 29.6 thread-hours (all on CPU).
- Wall time on 23 threads: 29.6 / 23 = 1.3 h at ideal scaling; 2.6 h at 50 percent efficiency (proposed pessimistic bound). Both are below the proposed 12 h cap and use at most the 23 cores of `environment.md`.

### GPU-hour arithmetic (Float64 arm only)

- 3.59e11 steps x 2.468 ns = 886 s = 0.25 h (AC13 ensemble rate); x 2 = 0.5 h.
- 3.59e11 steps x 9.137 ns = 3280 s = 0.91 h (one-wave rate); x 2 = 1.8 h.
- GPU estimate: 0.5 to 1.8 GPU-hours. If the GPU is used for the Float64 arm, the remaining CPU load is Double64 and BigFloat: 18.4 + 4.6 = 23.0 thread-hours (1.0 wall hour at ideal scaling on 23 threads).
- The GPU holds 8 GB; chunk sizes follow the tangent-core spec (one launch per chunk when n = 1e5 fits, D44j) and are not set here.

### K19 and K20

- K19 (backend): before the first run, rerun the AC13 throughput record (`verify/gpu.jl`) at 23 CPU threads (the thread count of the M3 runs). The existing 4-thread record is copied into the run manifest before the rerun overwrites it, and both values (4 threads and 23 threads, CPU and GPU ns per orbit-step) are recorded in the manifest. The Float64 arm then follows `auto` as resolved by the rerun. The CPU-only total of 29.6 thread-hours above is the budget if `auto` picks the CPU; the GPU estimate applies if it picks the GPU. Either choice fits the limits.
- K20 (Double64): Double64 is CPU-only, so the precision arm is sized against its CPU cost above: S_L of 100 orbits, Double64 horizons 1e7 (4 methods) and 1e8 (1 method), all on all 23 CPU threads. Changing the S_L size or the Double64 horizons is an amendment.
- Audit line "Budget" of the stage contract: 29.6 thread-hours, below any limit of `environment.md` (no stated hour limit; core limit 23 respected).

## Analysis plan

Every figure and table, and the hypothesis each decides. All plots use the stored series.

| # | Item | Content | Decides |
|---|------|---------|---------|
| T1 | Panel table | cells, counts per design class, empty cells, orbits dropped as "drift" | none |
| T2 | H1 table | pass fractions of H1a, H1b, H1c on E_ch per K; list of failing orbits; controls table | H1 |
| F1 | lambda_m(n) curves | one orbit per cell, all methods overlaid, log-n axis, Delta_V band | H2 context |
| F2 | R(n) vs n | median and bootstrap band per design class and K | H2 |
| T3 | H2 table | median R(n_H2) and interval over W | H2 |
| F3 | n_conv distributions | by design class and K | H3 |
| T4 | H3 table | median n_conv over W and S, ratio, interval, share of undefined n_conv | H3 |
| T5 | H4 table | lambda per method at T = 1e9 against 0.298 and 1.172; Double64 cross-check; Shevchenko2014 and ln(K/2) comparisons | H4 |
| T6 | Label table | two-axis label counts per cell; labels that changed with direction, precision, estimator | H5 |
| F4 | Label stability | per-test unresolved shares | H5 |
| T7 | H5 table | H5a, H5b, H5c fractions | H5 |
| F5 | Symplectic residuals | eps_pair(n), eps_det(n) against tolerance multipliers by method | gate and reported |
| F6 | Cost | ns per orbit-step by method and precision; relative to L5 | reported |
| F7 | Sensitivity | tau, eps (L2), direction, wrapped vs lifted effect on lambda(n_H2) and on n_conv | reported |

Rules for the analysis:
- Hypotheses are decided in the order H1, H2, H3, H4, H5, with the numbers above and nothing else. No threshold is tuned after a look at the data.
- Any analysis not listed here (extra figure, different subset, different n_H2) is labelled "exploratory" in the report.
- A hypothesis whose eligible set has fewer than 20 orbits (after drops) is reported as "not decidable", with the count, and is neither a pass nor a fail.
- Results for H2 and H3 are always reported by K and design class as well as pooled.
- The report cites only verified keys of `/_config/bibliography/references.bib`; the keys used here are FrahmShepelyansky2010, CarlettiLaskar2000, SalaArtuso2015, Shevchenko2014, Skokos2010, HarsoulaKaramanosContopoulos2019, Meiss2015 and Chirikov1979. None is in `unverified.md`.

## Stopping rules

Unit of work = chunk: one (cell, method or run type, precision, backend) block of orbits at one horizon, resumable (K8, K13).

A chunk fails if any of:
- non-finite value (NaN or Inf) in any recorded series;
- eps_det or orthogonality loss or eps_pair beyond the tolerance multipliers of `numerics-defaults.toml` for that orbit-method pair, in more than 5 percent of the chunk's orbits (proposed); the exceeding orbits are labelled "unresolved" regardless;
- an exact-exponent control exceeds the 1e-8 bound of H1 (this is a code defect, stop and report);
- measured cost of the first completed chunk of an arm exceeds 2 times its Budget line (proposed): the arm is paused and the budget is re-estimated by amendment;
- a Double64 or BigFloat run disagrees in the first 100 steps with the Float64 run beyond 1e-10 on the orbit coordinates (proposed: indicates a setup error, not chaos).

Handling:
- A failed chunk is rerun once with the same seed on CPU. A second failure is recorded; its orbits get both labels "unresolved" and are excluded from the pass fractions.
- The campaign stops if more than 5 percent of all chunks have failed twice (proposed).
- There is no early stopping on scientific results: no hypothesis is decided from partial data. The only early-stop causes are the failures above and the wall cap (12 hours, proposed).
- Interruptions resume from the last stored chunk; the manifest records versions, seed, backend, threads and the throughput values of K19.

## Amendment rule

Frozen after study/01. A change after freeze is an amendment: one numbered entry in the header's amendment log, with date, what changed, why, and the git SHA of the commit that makes it. Amendments are allowed only before stage 03-run produces any data file (including pilot or screening runs that touch panel orbits); the screening grid of stage 02 is not data in this sense because its results are discarded. After data exists, a changed threshold is not an amendment: it becomes an exploratory analysis, labelled as such, and the preregistered decision stands.

## Open points for the owner at the checkpoint

- Confirm the proposed values marked "proposed" above, in particular H1 tolerances, n_H2 = 1e4, the 0.10 and 10x thresholds of H2 and H3, and the label rules (these thresholds make D36 and D37 decidable; D37 itself flags the structure axis as the risk).
- Confirm K = 0.5 and K = 10 as additions to the K values named in the scope, and the 12 h wall cap.
- The mapping from the D37 phrase "stated canonical norm" to the Euclidean 2-norm on lifted (x, p) is proposed here.

## Amendment A1 (2026-10-17)

Approved by the owner at the module-build/01 checkpoint of `experiments` (decision D49), before stage 03-run produced any data. Base sections above are left as frozen; where they differ, this section wins.

1. Hash rule. The header cell "SHA-256 of this file (current version)" holds the SHA-256 of this file with the whole cell (hash and note) replaced by the word SELF. Rule in full: the whole cell is replaced by the word SELF before hashing (UTF-8, LF line ends); each amendment updates this cell and logs the hash it replaces. The freeze hash fa84e883... was computed under the same rule with the freeze-time cell text.
2. T predicate (Selection rule, item 2). Under the frozen threshold class T is empty: at K = 0 every screening point has lambda_scr = ln(sqrt((cos a + N sin a)^2 + sin^2 a)) / N, about 1.12e-4 for a = pi/4 and N = N_scr = 1e5, above 1e-4, and regular orbits at K > 0 sit on the same ln n / n floor (Skokos2010, see `benchmarks.md`). New rule: T requires lambda_scr <= theta_T, the growth ratio g = S(N_scr) / S(N_scr / 10) < g_T with S(N_scr / 10) > 0 (S is the cumulative log stretching, so g = r(N_scr) of "Label rules"), and the frozen distance condition. g tends to about 1.25 for regular decay and to 10 for constant positive lambda. theta_T = 5e-4 (the owner approved 1e-3 at the checkpoint; it is set to 5e-4 so that theta_T stays strictly below the S threshold 1e-3 and T is disjoint from S, as the module requires; the growth ratio carries the separation).
3. Screening deviation direction: the primary direction of Methods (angle pi/4), one vector (k = 1).
4. Windowed estimate (S predicate): the exponent over disjoint blocks of w = 1e3 steps, e_i = (S((i + 1) w) - S(i w)) / w; the median over blocks is the "median windowed value"; a stuck stretch is a run of consecutive blocks with e_i < 0.2 times that median, its length the block count times w.
5. Newton periodic-point search: at most 50 iterations per seed; two accepted points are the same orbit when their torus distance (after reduction to the primitive period) is below delta_dup = 1e-8.
6. Subset caps: S_D cap 210 and S_L cap 110, so that full population is not cut (7 x 27 + 13 = 202, 4 x 27 = 108). The two control slots left after the 10 fixed points go to the Newton orbits of lowest period, ties broken by K index (module rule).
7. Machine-readable panel inputs. The block below is the single home of the panel inputs the module reads; its values equal the frozen prose except where items 2 to 6 change them.

```toml
[panel]
K = ["0", "0.5", "0.971635406", "pi/2", "2pi", "10"]
n_c = 40
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
```

K and mu entries are decimal or closed-form strings; the module evaluates them in BigFloat and rounds once to each working type.

## Amendment A2 (2026-10-17)

Approved by the owner before stage 03-run produced any data (decision D50). The first build of this panel refused at K = 2pi: the binary64 value of 2pi lies 2.4e-16 below the K at which the period-1 accelerator mode is born, and the periodic-point search accepted a ghost pair near x = pi/2 only 1.0e-8 apart, below 10 delta_dup. Two accepted orbits of the same period closer than 10 delta_dup whose traces both lie within sn_trace of +-2 are now one orbit of class parabolic (near saddle-node), listed in the panel manifest; any other close pair still stops the build. sn_trace = 1e-6 (proposed: for the observed pair t - 2 is about 2e-8, and a genuine pair of distinct orbits that close with traces far from 2 is still refused) is added to the panel block of Amendment A1, item 7, marked there. K stays at 2pi, so H4 is unchanged.

## Amendment A3 (2026-10-17)

Approved by the owner before stage 03-run produced any data (decision D51). Supersedes Amendment A2: with sn_trace = 1e-6 the K = 2pi build still refused, and with 1e-5 it completed but kept about 1000 ghost periodic rows near x = pi/2 and two ghost controls. The periodic-point search now accepts a Newton solution only if its position error bound (residual divided by the smallest singular value of J_q - I) is at most delta_dup / 10; residual-accepted points that fail it are rejected as ill-conditioned and counted in the panel manifest. This uses only rho and delta_dup, so the key sn_trace is removed from the panel block of Amendment A1, item 7. K stays at 2pi; H4 is unchanged.
