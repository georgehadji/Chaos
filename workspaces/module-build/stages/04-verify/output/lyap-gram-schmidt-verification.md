# Verification: lyap-gram-schmidt

Spec: `../../01-spec/output/lyap-gram-schmidt-spec.md` (criteria 1 to 14, Amendment D35). Package tests: `test/lyap-gram-schmidt/` (`gs_*.jl`). Oracle script: `verify/lyap-gram-schmidt.jl`, run as `julia +release -t 2 --project=verify verify/lyap-gram-schmidt.jl` from the repo root; full output in `lyap-gram-schmidt-verify-log.txt` (lines `AC<n> | quantity | value | bound | verdict`; the quantity text names the case that attains the worst value, "[worst at ...]").

Run facts. Julia 1.13.1, 2 threads, base seed 20261012 with offset block +50000 (package tests +20000, standard-map verify 40000+, lyap-two-trajectory verify 45000+), so every sample is fresh. Source: working tree of branch `m2-p1-modules` (`src/Lyapunov/` unchanged by this stage). Wall time 5635 s on a loaded machine (other verify agents running); not a benchmark. Section times (s): AC1 26, AC2 10, AC3 97, AC4 131, AC5 419, AC6 142, AC7 39, AC8 553 (of which the Double64 10^6-step run 297), AC9 1757 (Double64 4096 orbits x 3 variants x 4 runs), AC10 903, AC11 70, AC14 12, oracle 174 + 15, R3 93, figures 644.

Reuse: the helpers of the package tests (`gs_helpers.jl` whole, `gs_s_loop` sliced from `gs_ac06_09.jl`), so a definition error in a helper would not be caught here; the independent pieces are the 256-bit and exact-rational references, the 512-bit Jacobi eigensolver of the sandwich, and ChaosTools. The package suite restricts Double64 to RGS (budget cut); here Double64 runs all three variants everywhere. Values are "worst ratio of error to bound" (pass <= 1) unless stated. eps_eff: Float64 2.22e-16, Double64 7.889e-31 (2^-100, as the precision spec).

## Oracle

ChaosTools 3.5.4 `lyapunovspectrum` (QR of the tangent deviation vectors every step, Householder via `LinearAlgebra.qr`), DynamicalSystemsBase 3.19.3; read from the installed source, hand-coded Jacobian `[[1+c,1],[c,1]]`, c = K cos x, lifted rule (no wrap).

**Finding for the owner (contradicts lyap-benettin Oracle P4).** `lyapunovspectrum(tands, N; u0)` calls `reinit!(tands, u0)`, whose `Q0` defaults to the identity columns, so a seeded `Q0` given to the `TangentDynamicalSystem` constructor is overwritten before the first step. Probe (log, `ORACLE` line): result with seeded Q0 minus result with identity = 0, whereas a manual QR of `J Q0` gives a different lambda_1 (0.7577). ChaosTools therefore cannot start from the identical seeded frame through this routine; both sides start from W_0 = I (tau = 1 on the Chaoslab side). The spec's "identical initial deviation vectors" is not achievable; pathwise agreement from W_0 = I replaces it. The lifted-rule orbits agree bitwise (K = 2, 1000 steps, difference 0).

| Case | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| AC13 regular K = 0.5 (S7 screen, 0 rejected), 1 orbit, n = 10^5, CGS / MGS / RGS | abs diff lambda_1 vs oracle (both 1.04e-4) | 2.342e-15 each | 1e-4 | PASS |
| same | abs sum lambda_1 + lambda_2 over theta_sum (oracle sum -9.5e-19) | 1.27e-5 / 2.31e-5 / 1.59e-5 | 1 | PASS |
| same (beyond spec) | abs diff lambda_2 vs oracle over 2 theta_sum | 0.0687 each | 1 | PASS |
| AC13 chaotic K = 2 pi, 1000 orbits, n = 10^4, identity frame | difference of ensemble means over combined SE, lambda_1 and lambda_2, all variants (means 1.171 / -1.171, SE 8.4e-4) | 1.867e-13 | 4 | PASS |
| same | max pathwise abs diff lambda_1, lambda_2 | 3.31e-14, 3.35e-14 | none (info) | INFO |
| same | abs sum lambda_1 + lambda_2 over theta_sum (M_max = 1+K) | 6.14e-4 | 1 | PASS |
| d = 4, A_S (G12c), fixed point 0, identity frame, n = 60 / 600: ChaosTools vs exact (1/n) ln vol_p, p = 1..4 | abs diff over p theta_sum | 0.0062 / 0.0437 | 1 | PASS |
| same, CGS = MGS = RGS (tau = 1) vs exact | abs diff over p theta_sum | 0.00232 / 0.00144 | 1 | PASS |
| same, variants vs ChaosTools | abs diff of partial sums over 2 p theta_sum | 0.00308 / 0.0216 | 1 | PASS |
| d = 4 n = 600 spectrum (finite time) | lambda (all variants and ChaosTools) | 0.9625, 1.316, -0.9619, -1.317 | asymptotic +-0.9624236501, +-1.3169578969 (G12c) | INFO (offset ln kappa(S)/n expected) |

Recorded line R3 (lifted rule, one seeded orbit, n = 10^6; values ChaosTools sum / Chaoslab abs sum for CGS / theta_sum): K = 0.5: -2.9e-19 / 2.2e-19 / 2.05e-14; K = 0.9716354: 2.8e-17 / 0 / 3.07e-14; K = 2: 5.6e-17 / 0 / 6.37e-14; K = 5: 0 / 0 / 2.45e-13; K = 10: -1.1e-15 / 0 / 8.31e-13 (MGS, RGS in the log). Both sums sit 2 to 5 orders below theta_sum. The "Chaoslab eps_sum of the same run" (diagnostics reducer) was not run; the QRSpectrum sum is reported instead. Non-scored sanity: K = 10 lambda_1 = 1.62 against ln(K/2) = 1.609.

## Evidence table

"Where": T = package test, V = this script. All numbers from `lyap-gram-schmidt-verify-log.txt`; "F" Float64, "D" Double64 (all three variants in both).

| Criterion | Quantity | Value (F / D) | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1 | cat map, 20 frames, tau {1,2,4}, n {40,400}: abs(lambda_1 - G12a) / theta_sum | 3.89e-3 / 1.22e-4 | 1 | PASS | T, V |
| 1 | abs(lambda_1 + lambda_2) / theta_sum | 7.06e-3 / 2.21e-4 | 1 | PASS | T, V |
| 1 | eigenbasis frame abs(lambda_1 - ln mu) / theta_sum | 2.05e-3 / 4.91e-5 | 1 | PASS | T, V |
| 2 | K = 0: abs(S_1 - (1/2) ln(n^2+1)) / (8 n eps), n up to 10^4 | 0.0247 / 9.3e-4 | 1 | PASS | T, V |
| 2 | abs(S_2 + S_1) / (8 n eps) | 0.0313 / 0.00195 | 1 | PASS | T, V |
| 2 | lambda_1(10^4) / (ln(10^4)/10^4 + 1e-12) | 1 - 1e-9 (printed 1) / same | 1 | PASS (marginal by construction: the margin 1e-12 is 1e-9 relative) | T, V |
| 3a | CGS vs MGS bitwise at k = 2, standard map (808 F / 168 D orbit comparisons) | 808 differ / 0 | 0 | **FAIL (F)**, PASS (D) | T (passes), V |
| 3a | random 2x2: CGS vs MGS Q and R (600 matrices) | 134 differ / 0; 4x2 Float64: 0 | 0 | **FAIL (F)**, PASS | T (passes), V |
| 3a | magnitude of the difference: max abs(S_i(CGS) - S_i(MGS)) / eps | 2.97e6 / 0 | information | INFO (see Findings) | V |
| 3b | A_S k = 4 tau = 3: windows of 200 where CGS and MGS Q differ | 200 / 200 | >= 1 | PASS | T, V |
| 4 | 256-bit chain on identical orbit, abs(lambda_i - ref) / theta_sum, K panel, tau {1,2,4} (F 20 orbits x 4 K; D 8 orbits x 3 K) | 0.0369 / 2.5e-4 | 1 | PASS | T, V |
| 4 | scratch constant max abs(lambda_i - ref) / (eps (G_max/tau + M_max)) | 0.590 / 0.0040 | 16 (spec scratch 0.9) | PASS | V |
| 5 | pair bound, K panel x 50 orbits, n = 1000, 36 (variant, tau) pairs | 0.0594 / 0.0048 | 1 | PASS | T, V |
| 5 | A_S k = 4, n = 600, tau {1,2}, 10 frames (F: 5 of 6 runs, CGS tau = 2 trips; D: 6 of 6) | 0.00770 / 9.7e-5 | 1 | PASS | T, V |
| 6a | RGS orth_max / (8 k eps), all runs of AC1, 4, 5, 6, 8, 9 | 0.198 / 0.0199 | 1 | PASS | T, V |
| 6b | CGS, MGS (k = 2) orth_max / (8 k (G_max+1) eps), same runs | 0.0890 / 0.00613 | 1 | PASS | T, V |
| 6b | recorded max loss / (eps G_w) per window, CGS / MGS | 1.26 / 1.29 (F), 0.088 (D) | 16 (spec scratch 1.05) | PASS | V |
| 6c | tau = 3: CGS orthogonality error raised, of 20 frames | 20 | 20 | PASS | T, V |
| 6c | MGS and RGS raised | 0 and 0 | 0 | PASS | T, V |
| 6c | tau = 3 min loss(CGS)/loss(MGS); min loss(MGS)/loss(RGS) | 358; 1032 | >= 100 | PASS | T, V |
| 6c | tau = 1: frames with loss(CGS) > loss(MGS) > loss(RGS); min ratios | 20 of 20; 2.75, 5.96 | 20; > 1 | PASS | T, V |
| 7 | A^20 (F) / A^40 (D, no fallback power needed): breakdown at handoff 1, correct orbit, step, backend, 200 frames x 3 variants | 200 / 200 each | 200 | PASS | T, V |
| 7 | cause split: fifth-pass / conditioning-rule frames; min eps G_w on completed frames | 0 / 200; F 1.31 to 1.96, D 27.8 to 62.2 | >= 1/2 | PASS | V |
| 7 | A^16 (F) / A^32 (D) no breakdown, n = 200: orth_max / bound (CGS = MGS 0.0405, RGS 0.125; D 0.00101, 0.00342) | as listed | 1 | PASS | T, V |
| 7 | I6: CGS and MGS finals and Q bitwise equal at k = 2 on A^16, F | false | true | **FAIL** (same cause as 3a) | T (passes), V |
| 7 | injected zero column 3 and 1, via propagate and reducer_update: kind, column, step, orbit, backend; NaN: propagate fires "deviation matrix" first, reducer_update alone gives breakdown, column 3, threshold Inf; F and D, 3 variants | all true | true | PASS | T (F), V |
| 8 | wrapped state: orbits (of 240) with x, W, S_i, record not bitwise equal to a manual loop, 3 mismatch kinds | 0 / 0 | 0 | PASS | T (F), V |
| 8 | 10^6 steps K = 5, 8 orbits, RGS: abs(lambda_1 + lambda_2) / theta_sum; finite; no health error | 5.0e-8 / 2.5e-5 | 1 | PASS | T (F), V |
| 8 | bytes allocated after warm-up, 3 variants | 0 / 0 | 0 | PASS | T (F), V |
| 9 | 4096 orbits K = 2 pi: propagate (2 threads; chunk 64, 1000, default) vs serial, not bitwise equal | 0 of 4096, all 9 cases / 0 | 0 | PASS | T (F reduced), V |
| 9 | 1 thread vs 4 threads, child processes, QRSpectrum(RGS), 2048 orbits, n = 1e3: hash of x, W, finals, Float64 / Double64 | equal / equal | bitwise | PASS | `verify/threads.jl` |
| 9 | bytes allocated, A_S d = k = 4, 10^4 steps, 3 variants | 0 / 0 | 0 | PASS | T, V |
| 10 | k = 1 vs Benettin, K panel x tau {1,4} x 50 orbits, n = 10^4: abs diff / bound; runs (of 1200) not bitwise equal; realised max diff | 0; 0; 0 (both types) | 1; 0 | PASS | T (F), V |
| 10 | k = 2, first 10 orbits | 0 / 0 | 1 | PASS | T (F), V |
| 11 | abs(sum_{i<=p} lambda_i - (1/n) ln vol_p) / (p theta_sum), A_S n = 60, 10 frames, tau {1,2,3} (CGS tau >= 2 trips, reported) | 1.25e-3 / 2.98e-5 | 1 | PASS | T (F), V |
| 11 | abs(sum_i lambda_i) / theta_sum | 3.52e-3 / 6.7e-5 | 1 | PASS | T (F), V |
| 11 | G8 sandwich, n = 12, 512 bits (Jacobi residual 1e-154): worst of (lower - sum) and (sum - upper) over p theta_sum | 1.76e-3 / 5.4e-5 | 1 | PASS | V |
| 11 | exact ln vol_p inside the bounds: violation (absolute) | 3.6e-17 / 1.1e-33 | 1e-12 | PASS | V |
| 11 | sandwich width over p theta_sum, smallest for p = 1..3 (not vacuous) | 1.6e8 / 4.6e22 | > 1 | INFO | V |
| 12 | RGS and MGS, regular K = 0.5 (S7: 1024 of 3072 regular, 17 replaced), n = 4e4, tau = 10, against a Double64 reference: orbits above 1e-6 CPU / GPU; max CPU-GPU gap, lambda_1 and lambda_2 | 1 / 1; 8.98e-08 (both variants) | count diff <= 3; gap <= 2e-6 | PASS | `verify/gpu.jl` |
| 12 | RGS and MGS, chaotic K = 10, Float64, 1024 orbits, n = 1e4, tau = 2: max per-orbit gap / sigma_FT; difference of means / SE, lambda_1 and lambda_2 | 4.992; 0.6086 (both variants) | 6; 4 | PASS | `verify/gpu.jl` |
| 12 | Double64 kernel | does not compile for the device | compile or reason | FALLBACK(kernel does not compile for the element type) | `verify/gpu.jl` |
| 13 | oracle package | see Oracle | | PASS | V |
| 14 | provenance fields, 126 checks per type (3 variants x 3 tol values) | 0 failed / 0 | 0 | PASS | T (F), V |
| O3 | max abs(sum lambda_i) / (eps (G_max/tau + M_max)), AC5, AC8, AC11, R3 runs | 0.113 / 0.0019 | 4 (spec scratch 0.015 and 0.056) | PASS | V |

Count: 14 of 14 criteria carry a number (criterion 12 from `verify/gpu.jl`). Three lines FAIL: AC3a (two lines) and the I6 line of AC7, one cause.

## Figures

- `lyap-gram-schmidt-orthogonality.png`: per-window orthogonality loss max abs(Q^T Q - I) against G_w for CGS, MGS, RGS; left standard map K = 2 pi, k = 2, tau 1 to 8; right A_S, k = 4, tau 1 to 4; with the RGS bound 8 k eps and the CGS/MGS bound 8 k G_w eps. CGS and MGS follow loss ~ eps G_w (k = 2, they coincide), RGS stays at a few eps independent of G_w; at k = 4 CGS rises above MGS by orders of magnitude and crosses its bound at tau = 3 and 4, as AC6(c) states.
- `lyap-gram-schmidt-oracle.png`: left, abs diff of lambda_1 from ChaosTools against n (K = 2 pi, one orbit, identity frame), three variants, with theta_sum; the curves lie at 1e-15 and below, 2 orders under theta_sum. Right, A_S d = 4: abs error of lambda_1 + lambda_2 against the exact (1/n) ln vol_2 for ChaosTools and the variants, n 10 to 600, with 2 theta_sum.

## Findings

1. **CGS and MGS are not bitwise equal at k = 2 in Float64 (spec I6, AC3a, AC7 I6).** Reproduced in a plain script outside `Pkg.test` (55 of 200 random 2x2 matrices) and in the verify run (134 of 600 matrices; 808 of 808 standard-map orbit comparisons; AC7 I6 false). The difference is in R12 only (1e-17 to 1e-16). Diagnosis from the log: of 2000 random 2x2 matrices, R12 of CGS equals `fma(q2,w2,q1*w1)` in 2000 and R12 of MGS equals the plain two-rounding `q1*w1+q2*w2` in 2000. The cause is StaticArrays' `dot` (`_vecdot` in `StaticArrays/src/linalg.jl` line 273, an `@inbounds @simd` loop): reassociation and fma contraction are the compiler's choice and differ between the CGS and MGS call sites. The package suite passes because `Pkg.test` runs with bounds checks forced on, which changes the generated code. Double64 and d = 4 are unaffected. Consequence: the per-orbit S_i then differ chaotically (max abs diff 2.97e6 eps in S over n = 1200, within the pair bound, which passes everywhere), so the claim "CGS = MGS bitwise at k = 2" holds only per compile mode. Options for the owner: replace `dot` in `_reduce_column` by an explicit ordered sum, or weaken I6 and AC3a to the pair bound. No src, test or tolerance was changed here.
2. **ChaosTools ignores seeded initial deviation vectors** (see Oracle); lyap-benettin P4 and the tangent-core Oracle text should be amended.
3. **AC12 (GPU)** now run in `verify/gpu.jl`: PASS for Float64; Double64 falls back with the recorded reason.
4. **AC7 scratch not reproduced.** The spec expects about 19 of 200 RGS frames needing a fifth pass and a minimum eps G_w of 0.99 (CGS, MGS) and 1.74 (RGS). Here 0 of 200 frames need a fifth pass for RGS (all 200 raise through the conditioning rule) and the minima over completed frames are 1.31 (CGS), 1.31 (MGS), 1.96 (RGS) for Float64 and 27.8, 27.8, 62.2 for Double64 (A^40; no fallback power 44 needed). The criterion passes with margin of 2.6 at least (>= 1/2).
5. AC6(c) min ratios (358 and 1032) are below the spec scratch minima (407 and 2000) but far above the bound 100. AC6(c) tau = 1 min ratios 2.75 and 5.96 (scratch 3.2 and 7.3), strictly greater than 1 as required (margin small, as the spec notes).
6. O3 constant 0.113 (Float64, A_S, CGS) is above the spec's scratch value 0.056 but 35 times below c = 4. AC4 scratch constant 0.59 (spec 0.9).
7. AC9: 2 threads, not 1 versus 4; the check is serial `propagate_orbit` against threaded `propagate` at three chunk sizes (tangent-core AC6 covers 1 versus 4 threads).
8. P1 vacuous-check band (eps^-1/2 < G_w < 1/(2 eps)): A^10 (G_w 2.3e8), A^16 (2.4e13) in Float64 and A^32 (5.6e26) in Double64 are inside it. CGS and MGS loss there is 0.07, 0.26 and 0.004 times eps G_w, while threshold_orth is 8.1e-7, 0.084, 0.0071 (so the check is nearly vacuous but the measured loss is 1e-3 to 1e-2 below it); RGS loss stays at 1 to 2 eps. Recorded, no pass or fail.

## Deviations from the impl notes, checked

- Pivot stored in `R[column, column]` on failure (zero column: 0; NaN column: NaN with threshold Inf): confirmed for all three variants.
- Double64 health error fields are Float64: confirmed.
- G_max, orth_max, count_max, sums, n unchanged after a failing handoff: confirmed (code 1, step of the failing handoff).
- `window_G == _G_from_Gw(window_Gw)` exactly; exact to 2.4e-14 relative against the SVD condition number at k = 2; upper bound at k = 3 (max cond / window_G 0.9999).
- The notes' "Known limits" say `window_Gw` has no overflow guard; the code has a rescale guard (`orthonormalize.jl` lines 111-112) and G_w is scale invariant to 0 relative difference at 1e200 and 1e-200. The note is stale.
- The header comment of `orthonormalize.jl` (line 3) still says "mutable MMatrix scratch allocated per call", while the notes say the MMatrix version was replaced. Stale comment.
- Zero allocation: 0 bytes for `orthonormalize` (4x4, F and D, all variants), and for `propagate_orbit` in AC8 and AC9.
- k = 1: orth_max stays 0.
- Workspace claim "no cut from the plan": the package suite runs Double64 RGS-only; this stage ran all variants in Double64 for criteria 1 to 11, 14 and found no Double64 failure.

## Runtime

Total 5635 s wall (94 min) with 2 threads on a loaded machine; the largest parts are first-call compilation of Double64 specialisations and CairoMakie (figures 644 s) and the Double64 sweeps (AC9 1757 s, AC10 903 s, AC8 553 s). Well above the 10 minute target of the verify house style; reduce AC9 and the Double64 10^6-step run for any rerun.

## Review findings

Independent review: fresh Sonnet subagent given the spec and the code, not the implementation notes. No CRITICAL and no HIGH finding. Disposition: closed (evidence named), fix (proposed change), record (decision row or spec amendment), owner (checkpoint question).

| Severity | Finding | Location | Disposition |
|---|---|---|---|
| MEDIUM | RGS fifth pass (`QR_PASSES`, code 2) untested; `count_max == 2` at tau <= tau_max never asserted | `test/lyap-gram-schmidt` | fix: one test case; see also Findings item 4 (no frame needed a fifth pass here) |
| LOW | `window_Gw` rescale branch untested | `src/Lyapunov/orthonormalize.jl:111-112` | fix: one test case (the guard exists; the impl notes say it does not) |
| LOW | G_max for k >= 3 is the bound G_w, not cond_2(R) (API Q5, deliberate) | `src/Lyapunov/orthonormalize.jl` | record: spec amendment |
| LOW | stale header comment "mutable MMatrix scratch" | `src/Lyapunov/orthonormalize.jl:1-3` | fix: comment |
| unverified | monitor takes G from its own inputs | `src/Diagnostics/monitor.jl:70-86` | owner: confirm intended |
| (verify) | thread count: AC9 ran 2 threads only (Findings item 7) | - | closed: `verify/threads.jl`, child processes `-t 1` and `-t 4`, outputs bitwise equal (`threads-verify-log.txt`) for QRSpectrum(RGS), Float64 and Double64 |
