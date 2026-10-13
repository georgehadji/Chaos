# Spec: diagnostics

Owner checkpoint 2026-10-13 (decision D24): Float32 (binary32) is dropped from M2. Every Float32 item in this spec (element-type lists, the screening tier, the Float32 GPU gate and Float32 test cases) is inactive until revisited at M4; stage 02 writes no Float32 tests.

Library layer Diagnostics (`_config/architecture.md`, "Module map" row Diagnostics). Reducers over the tangent stream: initialise, update at each step or snapshot, finalise. Draft for the M2 checkpoint; status: unreviewed.

## Scope

- Catalog rows: `_config/catalogs/lyapunov.md` section F (built-in checks every estimator reports), items: symplectic pairing residual, area-preservation residual, orthogonality loss, sum of exponents. Retained outputs eps_pair and eps_det: `_config/research-question.md`, "Outputs retained".
- Which residual applies to which estimator row: the det and symplecticity residuals need a one-step tangent (all variational rows L1, L3-L6, L9; not L2, which has no tangent). Sum and pairing need all 2N exponents (full-spectrum rows L3-L6, L9). Orthogonality loss needs a basis of 2 or more vectors (L3-L6).
- Prerequisites (M2 siblings, all `planned` at drafting, `_config/build-order.md` row M2 requires all verified): tangent-core (stream, renormalisation interval tau as recorded parameter), standard-map (Y1) and the linear symplectic maps (Y5). Diagnostics may not import the Lyapunov layer (`_config/architecture.md`, "Layering"): finite-time exponents arrive as data supplied by the estimator that drives the run.
- Deferred, not in this spec: energy and modified-energy errors (M5, flows; key `energy_error_tol` stays with the first integrator spec); the SALI decay versus lambda_1 - lambda_2 cross-check (M4, needs the indicators layer).
- Keys proposed here: `pairing_residual_tol`, `det_residual_tol` (section Algorithm, "Tolerances"). `orthogonality_tol` belongs to lyap-gram-schmidt; this module defines the quantity and compares it with whatever tolerance the run carries.
- Open points for the owner: (O1) both keys are dimensionless multipliers of the machine epsilon of the element type, not absolute numbers (precision is a run parameter, `_config/numerics-policy.md`); the toml comments should say so. (O2) For N >= 2 the pairing health check fires only at finalise and needs a finite-time allowance kappa that has no default until the first N >= 2 system (M3). (O3) The constant c in the pairing bound is an assumption about the orthonormalisation; lyap-* specs should confirm or replace it.

## Equations

Notation. N pairs, phase-space dimension 2N. DT_j is the one-step tangent (Jacobian of the map) at orbit point x_j. Coordinates are ordered (positions, then momenta), as in `_config/research-question.md`, "Ground truth" (x, p). eps is the machine epsilon of the element type T, u = eps/2 the unit roundoff, M_j = max_ik |DT_j[i,k]| (for the standard map always >= 1 since a unit entry is present; for other maps thresholds use max(1, M_j)). lambda_i(n), i = 1..2N, are the finite-time exponents after n map steps in the estimator's own index order (no re-sorting). Snapshots occur every tau steps; a window is the tau steps between snapshots. Q is the 2N x k orthonormal basis after an orthonormalisation, k <= 2N.

Given (source):
- E1. Omega = [[0_N, I_N], [-I_N, 0_N]] (the matrix J_2N), Skokos2010 p. 7-8 (block form printed at the foot of p. 7 and top of p. 8, before Eq. (2)).
- E2. Symplectic map: DT^T Omega DT = Omega, Skokos2010 Eq. (5), p. 8.
- E3. Sum of all exponents: sum_{i=1..2N} chi_i = 0, Skokos2010 Eq. (47), p. 28.
- E4. Pairing: chi_i = -chi_{2N-i+1}, i = 1..N, Skokos2010 Eq. (48), p. 28 (Skokos attributes the proof to his ref. [13]).
- E5. Finite-time form X_1(t) = chi_1 + (b + z(t))/t, b a constant, z zero-mean noise, "valid both for regular and chaotic orbits", Skokos2010 Eq. (63), p. 36; regular orbits give X_1 of order ln t / t, Eq. (62) and the sentence below it, p. 35. Skokos states E5 for the maximal exponent only.
- E6. Standard map p' = p + K sin x, x' = x + p' and DT = [[1 + K cos x, 1], [K cos x, 1]], det = 1: `_config/research-question.md`, "Ground truth" (x, p ordering).

Definitions (the two retained residuals are given by `_config/research-question.md`, "Outputs retained"):
- D1. eps_pair(n) = max_{i=1..N} |lambda_i(n) + lambda_{2N+1-i}(n)|. (Maximum over i = 1..2N gives the same value by symmetry.)
- D2. eps_det(j) = |det DT_j - 1|, per map step. Window maximum W_det(k) = max over the steps of window k; running maximum R_det(n) = max over steps 1..n.
- D3. eps_symp(j) = max_{a,b} |(DT_j^T Omega DT_j - Omega)[a,b]|.
- D4. eps_orth(n) = max_{a,b} |(Q^T Q - I)[a,b]|, Q the basis just after orthonormalisation at the snapshot (before it, the working vectors are not meant to be orthonormal).
- D5. eps_sum(n) = |sum_{i=1..2N} lambda_i(n)|.

Norm choice (D3, D4): maximum absolute entry, for both. Reasons: (i) one running scalar, no square root, no k x k matrix stored; (ii) a maximum is exact and order-independent, so the value is bitwise identical for any reduction order, which makes thread invariance automatic where a Frobenius sum would depend on summation order; (iii) off-diagonal entries are the cosines between basis vectors and diagonal entries are |norm^2 - 1|, directly comparable with eps; (iv) the Frobenius norm is bounded by 2N times the max-abs value (derived: at most (2N)^2 entries), so a Frobenius figure is recoverable.

Derived identities (elementary algebra, steps shown):
- I1 (2D symplecticity equals det). For A = [[a, b], [c, d]] and Omega = [[0, 1], [-1, 0]]: Omega A = [[c, d], [-a, -b]], A^T (Omega A) = [[0, ad - bc], [-(ad - bc), 0]] = det(A) Omega. Hence for 2N = 2, eps_symp = eps_det in exact arithmetic and the symplecticity residual adds nothing in M2; it is informative from 2N = 4 on.
- I2 (sum identity). For a QR-type estimator the window k product W_k = (product of DT over the window) Q_{k-1} = Q_k R_k, with Q_k orthogonal. Then |det R_k| = |det W_k| = prod |det DT| = 1 for symplectic DT (E2 gives det DT = +-1 by taking determinants, and +1 for the maps in scope), so sum_i ln R_k[i,i] = 0 at every window, and lambda_i(n) = (1/n) sum_k ln R_k[i,i] gives sum_i lambda_i(n) = 0 at every finite n, any tau, any initial basis. The sum residual therefore has no finite-time dynamical part: it is round-off plus the det residual of the actual matrices. The same holds for an estimator built on singular values, since the product of the singular values is |det|.
- I3 (pairing for 2N = 2). i = 1 gives lambda_1 + lambda_2, so eps_pair = eps_sum: pure round-off at every n. No finite-time term exists for the standard map.
- I4 (pairing for 2N >= 4). Individual QR-diagonal exponents obey a form like E5 (inference from E5; Skokos2010 states it for chi_1 only). The pair sums need not cancel their (b_i + b_{2N+1-i})/n terms, with orbit-dependent constants that cannot be preregistered. So eps_pair(n) has a part decaying like 1/n besides round-off.
- I5 (pairing is exact for singular values at finite n). From E2, DT^{-1} = Omega^{-1} DT^T Omega, and Omega is orthogonal, so the singular values of DT^{-1} (the reciprocals 1/sigma_i) equal those of DT^T (the sigma_i) as a set; for a product of symplectic matrices the same holds. Hence sigma_i sigma_{2N+1-i} = 1 and a singular-value-based finite-time spectrum (row L9) pairs exactly in exact arithmetic at every n, while QR-diagonal exponents (L3-L6) pair only asymptotically.

Floating-point bounds (derived; model fl(a op b) = (a op b)(1 + delta), |delta| <= u):
- B1 (det of one standard-map step). Exact entries a11 = 1 + K cos x, a21 = K cos x, a12 = a22 = 1. Stored entries m11 = a11 (1 + d1), m21 = a21 (1 + d2) with |d1|, |d2| <= u; this covers shared or independent rounding of the two entries (for example fused versus unfused arithmetic). Then det of the stored matrix minus 1 = a11 d1 - a21 d2, and evaluating m11 m22 - m12 m21 adds one subtraction rounding of at most u |det|. So eps_det <= u (|a11| + |a21| + 1) <= 3 u M = 1.5 eps M to first order. If both entries share one rounded value of K cos x the bound tightens to u (2 + |K|). The error in cos itself cancels, provided both entries use the same evaluation. Scratch check by the spec author (binary64, 10^6 random lifted points per K in {0.5, 0.9716, 2, 5, 10, 100}): max eps_det / (eps M) = 0.5.
- B2 (what the det check does and does not test). Per-step det of an analytic 2 x 2 Jacobian is nearly an algebraic identity: it tests the Jacobian formula and the arithmetic of the element type, not the propagation. The injected-fault criteria (A5) make it non-vacuous.
- B3 (symplecticity residual of one step, general N). Entry (a,b) of DT^T Omega DT is a sum of 2N products of entries (Omega DT only permutes rows with signs, exactly). Each product has relative error <= u, and recursive summation of 2N terms adds <= (2N - 1) u times the sum of absolute terms, which is <= 2N M^2. So eps_symp <= 2 N^2 eps M^2 to first order.
- B4 (sum of exponents per window). For 2N = 2 the second direction is what remains after removing the component along the first: |w_2| is about sigma_1 while R[2,2] is about 1/sigma_1. A backward-stable orthonormalisation perturbs W by about u |W|, so R[2,2] has absolute error about u sigma_1 and relative error about u G_k, with G_k = sigma_1 / sigma_2N of the window product (= sigma_1^2 for symplectic 2N = 2). The window contributes at most c u G_k to ln R[1,1] + ln R[2,2], with c a constant of the orthonormalisation (assumed c <= 4 for 2N = 2, to be confirmed, O3). The det error of the tau actual matrices adds at most 1.5 eps M tau per window. Accumulating the logs with compensated summation (`_config/coding-conventions.md`, "Code") and dividing by n = (number of windows) x tau gives, with G_max and M_max the running maxima:
  eps_sum(n) <= 2 eps G_max / tau + 1.5 eps M_max + (log rounding, at most eps times ln G_max, absorbed) <= 4 eps (G_max / tau + M_max).
  The bound does not grow with n. The compensated accumulation is what keeps the bound free of n; plain accumulation has only a worst-case bound growing linearly in n, so the estimators must follow the convention. Scratch check by the spec author (cat map and standard map K in {0.5, 2, 5, 10}, tau = 1, 2 x 10^4 steps, Householder QR, binary64): measured eps_sum <= 0.2 eps (G_max + 1) in every case, far below the bound.

## Algorithm

Inputs delivered by the tangent stream and the driving estimator (data, not imports):
- Every snapshot (every tau steps, the tangent-core handoff): orbit id, step index n, orbit point x_n; the one-step tangent DT_n at x_n is evaluated through the system (the Systems layer is below Diagnostics, so this is an allowed import). With the default tau = 1 this samples every step; with tau > 1 eps_det and eps_symp are sampled at snapshot points only, and the window maximum of D2 is the maximum over the sampled steps. A per-step hook in the stream is not added (tangent-core spec hands off only every tau steps).
- Every snapshot: the basis Q after orthonormalisation, the exponents lambda_i(n) (all 2N when the estimator is full-spectrum), and G for the window, taken from the same factorisation that produced the exponents.

Per orbit, the reducer keeps a fixed small state (running maxima, current values, first-exceedance record) and appends the retained series to caller-provided storage at each snapshot: window max of eps_det, running max of eps_det, eps_symp, eps_orth, eps_sum, eps_pair, and the thresholds in force. The full series are the retained outputs of `_config/research-question.md`; the sampling schedule is the study's.

Update rules:
1. Each snapshot: compute eps_det at the snapshot point (D2), update window and running maxima, compare with the det threshold.
2. Each snapshot: compute eps_orth, eps_sum, eps_pair, and eps_symp on the last step's tangent (cost O((2N)^3), so snapshot cadence rather than per step); update G_max and M_max; compare.
3. Finalise: report the running maxima, the final residuals and whether any threshold was exceeded.

Tolerances (run parameters, recorded in provenance with the element type; a study overrides them only in its preregistration, `_config/numerics-policy.md`). Each threshold is multiplier x eps(T) x growth model:
- Det: threshold_det(j) = det_residual_tol x eps(T) x (1 + M_j). Proposed det_residual_tol = 4. Justification: the derived worst case is 1.5 eps M (B1), and the threshold-to-bound ratio is 4 (1 + M) / (1.5 M), about 2.7 or more for any M >= 1. For the standard map 1 + M <= 2 + |K|, so at K = 10 and binary64 the threshold is about 1e-14.
- Symplecticity: threshold_symp = det_residual_tol x eps(T) x N^2 x (1 + M)^2, covering B3 (derived 2 N^2 eps M^2) with a margin of at least 2.
- Sum and pairing (N = 1, where they coincide, I3): threshold_sum(n) = pairing_residual_tol x eps(T) x (G_max / tau + M_max), checked at every snapshot, independent of n. Proposed pairing_residual_tol = 16: four times the derived worst-case constant 4 (B4). The scratch measurement sits about 80 times below this, so stage 04 may tighten it; it must not loosen it without redoing B4.
- Pairing for N >= 2: eps_pair(n) is recorded at every snapshot but the health check is applied only at finalise, with threshold = pairing_residual_tol x eps(T) x (G_max / tau + M_max) + kappa / n_end. Decision and reason: I4 shows a part decaying like 1/n with orbit-dependent constants, so a per-snapshot check would need an n-dependent tolerance nobody can preregister, and any fixed value would either fire spuriously early or be blind late; at finalise a single kappa, preregistered per study from the length of the run, is meaningful. kappa has no default and no use in M2, where every system has N = 1 (O2). The sum residual, which has no such part (I2), stays checked at every snapshot for all N.
- Orthogonality: threshold supplied by the run (`orthogonality_tol`, lyap-gram-schmidt spec).

Health error. When a residual exceeds its threshold, or is NaN or infinite (always, regardless of any threshold), a numerical health error is raised with: residual name (det, symplecticity, orthogonality, sum, pairing, orthonormalisation breakdown), value, threshold, step index, orbit id. If several residuals fail at the same snapshot the reported name follows the fixed order det, symplecticity, orthogonality, sum, pairing, orthonormalisation breakdown (the last kind is raised by the orthonormaliser, lyap-gram-schmidt and lyap-householder-qr, not by this module). Over an ensemble the reported orbit is the lowest orbit id that failed, so the report does not depend on thread count or scheduling. A kernel cannot throw on the GPU: it records the first-exceedance record per orbit and the host raises the error at chunk end, which triggers the CPU fallback with the reason recorded (`_config/architecture.md`, "Compute").

## Invariants

- For every step of a symplectic map: eps_det <= threshold_det and eps_symp <= threshold_symp (E2, B1, B3).
- For a full-spectrum QR or singular-value estimator at every n: eps_sum(n) <= threshold_sum(n), with no n-dependent part (I2).
- For 2N = 2: eps_pair(n) = eps_sum(n) and eps_symp = eps_det in exact arithmetic (I1, I3).
- For a singular-value estimator: eps_pair(n) is at round-off level at every n (I5); for a QR estimator at N >= 2 it is not, and the excess is the finite-time part (I4).
- After an orthonormalisation: eps_orth <= the run's orthogonality tolerance.
- Reducers are per orbit and independent: no value of one orbit depends on another, on thread count or on backend choice beyond round-off.
- A residual is never silently clipped, rounded to zero or skipped when not finite.

## Compute

Policy: `_config/architecture.md`, "Compute" (decision D17). The reducers are per-orbit, embarrassingly parallel and fused into the same KernelAbstractions kernel as the orbit ensemble, so the tangent matrices are never copied to the host; one source serves the CPU backend (the verified reference, the only one CI runs) and CUDA. Element types: binary32, binary64 and double-double (DoubleFloats, isbits) run on both; arbitrary precision is not isbits, so it is CPU-only and on a GPU request the recorded fallback reason is "kernel does not compile for the element type" (owner decision, M2 kickoff). Thresholds are expressed in eps of the element type so they transfer between types (the GTX 1050, compute capability 6.1, has weak Float64 throughput, risk K15). Memory per orbit: a fixed state of about 16 scalars of type T (current values, window and running maxima, G_max, M_max, first-exceedance record), independent of n, plus 6 values of type T per recorded snapshot in the series (the sample schedule is the study's; chunk size follows free memory, never a fixed number). The update allocates nothing after warm-up. CPU-vs-GPU agreement: the same pass/fail verdict on the same orbit and threshold, including injected faults; residual values agree to within the threshold itself, never bitwise (FMA contraction and reduction order change round-off residuals at O(1) relative level, `_config/architecture.md`, "Compute", Reproducibility).

## Oracle

- Analytic (package suite): Y5 linear maps with integer entries, for which every residual is exactly zero; the cat map [[2,1],[1,1]] with eigenvalues (3 +- sqrt 5)/2 and exact exponent +-ln((3 + sqrt 5)/2) = +-0.96242365012 (derived: trace 3, det 1, mu = (3 + sqrt 5)/2); a fault map with exact det 1 + delta.
- Independent recomputation in arbitrary precision (256 bit), on the identical orbit: the binary64 orbit points are stored and the entries of DT are recomputed in arbitrary precision at those exact points (never a separate arbitrary-precision trajectory, which would diverge, risk K18); det and sum residuals are recomputed from them. Short runs (n <= 10^3, no external package) belong to the package suite; the long panel runs (n = 10^6) belong to the verification environment driven from module-build stage 04.
- Stage 04 verification environment (not the package suite, owner decision, M2 kickoff): the ChaosTools full-spectrum routine on the same K and start gives an independent sum of exponents, reported next to eps_sum for comparison only (not pass/fail, since ChaosTools has no pairing or det output); arbitrary-precision recomputation over the test panel K in {0.5, 0.9716354 (`_config/bibliography/benchmarks.md`, golden-mean row), 2, 5, 10}; GPU-vs-CPU comparison on the GTX 1050 with the result in the 04-verify report.
- No published lambda(K) value is used as a hard criterion (risk K18).

## Acceptance criteria

Tolerances use theta_det = det_residual_tol x eps(T) x (1 + M) and theta_sum = pairing_residual_tol x eps(T) x (G_max / tau + M_max) as defined above. Element types in the package suite: binary32, binary64, double-double; the arbitrary-precision items run CPU-only.

1. Y5 integer maps: shear [[1,1],[0,1]], cat [[2,1],[1,1]], order-4 rotation [[0,1],[-1,0]], order-6 rotation [[1,-1],[1,0]] (det 1, entries exact), 1000 steps each. Pass: eps_det = 0 and eps_symp = 0 exactly at every step (zero tolerance, all arithmetic exact), in every element type, and no health error.
2. Y5 rotation by angle 0.7311 (entries cos, sin rounded): over 10^6 steps eps_det <= theta_det (that is <= 8 eps(T) since M <= 1) at every step; eps_orth <= 8 eps(T) at every snapshot when the test re-orthonormalises the rotated basis each step with a reference orthonormalisation (scratch measurement 4 eps in binary64).
3. Cat map, tau = 1, n = 10^5, reference QR as estimator: eps_sum(n) <= theta_sum at every snapshot, with G_max = ((3 + sqrt 5)/2)^2 = 6.854 (derived: sigma_1 = mu, G = mu^2); pairing for N = 1 equals eps_sum within 0 (same quantity).
4. Pairing index mapping, synthetic exponents: N = 1, lambda = (L + d, -L), L = ln mu, d = 2^-20: |eps_pair - d| <= eps(T) (|L| + d). N = 2, lambda = (a, b, -b + d2, -a + d1) with d1 = 2^-20, d2 = 2^-18, a = 0.96, b = 0.31: |eps_pair - d2| <= eps(T) x 2 max|lambda_i|, and adjacent-pair indexing (i, i+1) would give a value near 1.3 and is rejected. The same inputs give eps_sum equal to |d1 + d2| to the same tolerance.
5. Injected faults, both sides of every threshold, with exact binary entries so the injected residual is known exactly. Det: A_delta = [[2 + delta, 1], [1, 1]] (det exactly 1 + delta), delta a power of two. With theta the threshold in force: delta in [2 theta, 4 theta) raises a health error at step 1 with name det, value equal to delta (bitwise), the given orbit id and step; delta in (theta/4, theta/2] over 10^4 steps raises none. Orthogonality: a basis with one vector scaled by (1 + delta): eps_orth = 2 delta + delta^2 within eps(T), same two-sided test against the run's orthogonality tolerance. Sum and pairing: add delta to lambda_1 only (N = 1): both residuals equal delta within eps(T), same two-sided test against theta_sum. Symplecticity: a 4 x 4 block map with one block det 1 + delta, two-sided against threshold_symp. Non-finite: a NaN or infinity in DT, in Q or in lambda raises with the right name for any threshold, including a threshold set to infinity.
6. Real-orbit calibration, standard map, K in {0.5, 0.9716354, 2, 5, 10}, 50 seeded random starts each, 10^4 steps each, binary64 and double-double: max eps_det / (eps(T) M) <= 1.5 (derived bound B1), no health error at the proposed keys. Tolerance 1.5 is the derived worst case, not a fit.
7. Stage 04 (long runs): over the panel and 10^6 steps, max eps_sum(n) / (eps(T) (G_max / tau + M_max)) <= 4 (the B4 constant). If violated the key pairing_residual_tol is revised before M2 exit.
8. Arbitrary-precision (256 bit) recomputation, short (package suite): along a 10^3-step binary64 orbit at K = 5, |eps_det(binary64 evaluation) - eps_det(arbitrary-precision evaluation of the same stored entries)| <= eps(T) at every step (they differ only by the final subtraction rounding, at most u |det| < eps).
9. Zero allocation: after warm-up the per-snapshot update allocates 0 bytes, for binary32, binary64, double-double. Arbitrary precision is excluded (it allocates by construction).
10. Thread invariance: the residual series, window and running maxima and the first-exceedance record of 64 orbits are bitwise identical with 1 thread and 4 threads; an injected fault in one chosen orbit is reported with the same orbit id for both, and with two faulty orbits the lower id is reported.
11. Backend agreement: on the CPU backend (CI) the kernel passes items 1, 5 and 6. Locally on the GPU the same items give the same verdicts on identical orbits and thresholds, and residual values differ by at most the threshold in force; the result goes in the 04-verify report. An arbitrary-precision request on a GPU falls back to the CPU with reason "kernel does not compile for the element type" recorded.

## Sources

- Skokos2010 (`data/literature/lyapunov-algorithms/Skokos2010.txt`, arXiv:0811.0882, arXiv PDF pages): block form of J_2N p. 7-8; Eq. (5) p. 8; Eq. (47), Eq. (48) p. 28; Eq. (62) p. 35; Eq. (63) p. 36.
- `_config/research-question.md`: "Ground truth" (map, Jacobian, ordering), "Outputs retained" (definitions of eps_pair and eps_det).
- `_config/catalogs/lyapunov.md`: section F (built-in checks), section D2 (analytic tangent source), rows L1-L9; `_config/catalogs/systems.md` rows Y1, Y5.
- `_config/architecture.md`: "Module map", "Layering", "Compute"; `_config/coding-conventions.md`: "Code", "Tests"; `_config/numerics-policy.md`, `_config/numerics-defaults.toml`.
- `_config/bibliography/benchmarks.md`: row "Golden-mean KAM torus breakup, standard map (Greene K_c)" (K_c = 0.971635406, quoted only as a panel K, not as a criterion).
- Owner decisions, M2 kickoff: double-double as extended precision, arbitrary precision for short CPU oracle runs only, oracle packages in the stage 04 verification environment. Decision D17; risks K15, K18 (`_config/risks.md`).
- Derived in this spec (not cited): identities I1-I5, floating-point bounds B1-B4, and the cat-map exponent; scratch checks quoted in B1, B4 and item 2 were run once by the author and are to be reproduced in stage 04.

## Amendment D30 (notes from the Batch B review)

Owner-approved with decision D30 (questions O1 and O7). No acceptance criterion changes, no renumbering.
- Health kinds: "orthonormalisation breakdown" is appended last to the name list and to the fixed reporting order (Algorithm, "Health error"). This module never raises it; the orthonormalisers do.
- B4 is the home of G: the growth G of a window is the condition number of the window's triangular factor (the sigma ratio of B4), supplied by the estimator that drives the run (lyap-gram-schmidt G10). The bound and the constant c are unchanged.
- The golden-mean panel value is 0.9716354 (benchmarks row, K_c = 0.971635406), in the Oracle and in criterion 6; the stage 04 report already used it.
