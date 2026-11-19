# Verification report: clv-ginelli (catalog row L12; spec D122, test amendment D123)

Stage 04 of `workspaces/module-build`, 2026-11-19. Script `verify/clv-ginelli.jl` (environment `verify/`, groups FULL1 to FULL3, FULL1p, FULL2p, AC3bd, Od13a to Od13e, Oe, R1, R2, R3idx, R4, R5, THR, FIG), raw output `clv-ginelli-verify-log.txt` (groups concatenated in this order; most groups ran in parallel processes, so wall-clock times other than R1 are loaded-machine times), figures `clv-ginelli-convergence.png` (left: distance between two seeds of (frame, C) against n_f = n_b, standard map, coupled maps, five F3 chaotic starts and orbit C; right: criterion 3 (a), fixture b, error of columns 2 to 4 against the backward transient with the predicted slopes dashed) and `clv-ginelli-angles.png` (R2, empirical CDF of the minimum principal angle). Code commit 9c0e4d4 (`src/Lyapunov/clv_ginelli.jl`, `clv_backward.jl`). Julia 1.13.2, up to 8 threads. Status: accepted at the stage 04 checkpoint after Amendment D124, 2026-11-19 (owner: "Accept").

## Evidence per criterion

FULL = the package test files at full sizes (`CLV_FULL_RUN = true`, one process per file). FULL1 44514 pass, 2 error (AC3b, AC5); FULL2 9045 pass, 1 error (AC7); FULL3 144 pass, 0 fail. The three errors have one cause, a test-side defect shown at the full size (next paragraph); no assertion about the module failed. FULL1p and FULL2p rerun the two files with that defect removed in the process only (the files are untouched): FULL1p 44858 pass, 14 fail (all AC3b, one start, below); FULL2p 17563 pass, 0 fail, 0 error. No size reduction of any test file.

Test-side defect (not a module defect): with 20 F3 chaotic starts the pilot of `clv_pilot` (300 handoffs) returns lambda_1 < lambda_2 for start 4 (a sticky start: pilot lambda 0.0049, 0.0110, ...), so `clv_nstar` returns n* = -6548 and the module correctly throws `ArgumentError` (n_f must be nonnegative). That ends AC3b, AC5 and AC7 at that start (AC7 part (c) was not reached). With a 1000-handoff pilot n* = 5653 for start 4 (the cap is 6000). This is a calibration weakness of U3 (pilot length), invisible at the package size of 2 starts.

| AC | Result | Evidence |
|----|--------|----------|
| 1 | PASS | worst delta / B_fix 0.195 (vectors), exponents 0.039 of bound, 30 runs, fixtures a to d, Float64 and Double64, tau in {1, 4}; negative controls pass |
| 2 | PASS | (a) 0.026, (b) 0.038, (c) 0.031 of bound; window mean against the delivered spectrum 0.99995 of its allowance (boundary term, as in the plan), 43 cases |
| 3 | PASS (a), (c); (b) FAIL at one of 20 starts | (a) worst slope deviation 0.0064 (bound 0.2). (c) worst ratio 0.247. (b) FULL1p: 20 F3 chaotic starts plus the two maps, 448 pass, 14 fail, all at F3 start 19 (N_b = 776 from the pilot): distance 3.3e-8 at N_b, 2.8e-10 at 2 N_b, 6.8e-14 at 3 N_b (log, AC3bd; the other 19 starts are at or below 6.9e-11 at N_b). The standard map (N_b = 5) and coupled maps (29) pass. Diagnosis: bound calibration, not code (the distance decays to the floor with the transient; n* = 40 / (g_min Dt) from a finite-time pilot has no margin on an orbit that passes a near tangency). R4 below gives the realized transients |
| 4 | PASS | order 2 errors 3.9e-4, 9.8e-5, 2.4e-5; order 4 7.6e-7, 4.7e-8, 3.0e-9 (h = 1/10, 1/20, 1/40; ratios 4.0, 4.0 and 15.9, 16.0) |
| 5 | PASS (FULL1p) | (a) 0.25 of 512 eps, (b) 0.30 of 1e-9, (c) 0.18 of 64 eps; FULL1 errored at start 4 (defect above) |
| 6 | PASS | worst ratio 0.625 of 16 eps (AC6b); AC6a, AC6c pass |
| 7 | PASS (FULL2p) | worst sine and cosine 1.2e-9 (bound 1e-6), 20 chaotic starts plus orbit C; forward exponents pair for 200 of 200 starts (bound 95 percent) |
| 8 | PASS | bitwise store-all against checkpointed, h_b in {1, 7, 64, 500}; call counts exact; bytes 288576 and 24528; sink order and peak bytes |
| 9 | PASS | bitwise first k columns, Householder and RGS (the first reading of U15 holds) |
| 10 | PASS | AC10a to c pass at full size (D123 tolerance 256 bytes); THR group: 1024 orbits, ntasks 1 against 8, chunk 64 against 1000, and a one-thread process against the 8-thread process, all bitwise (C, V, lambda, final state, forward exponents) |
| 11 | PASS | all within B_v (c_v = 1e4). R3idx at other seeds, error / B_v: Float64 F3 orbit C 0.79 (delta, index 2), 0.14 (ln Lambda, index 3); Double64 F3 orbit C 5.1e-3; standard map 3.0e-4 (Float64), 1.0e-3 (Double64) |
| 12 | PASS | AC12a, AC12b pass (FULL3) |
| 13 | PASS at n* except coupled maps Double64 and coupled maps 256 bit (FAIL as written, converges at 2 n*) | see "Criterion 13" |
| 14 | PASS | AC14a to AC14e pass (FULL3). Not encoded and not run: "tau h incompatible with the integrator" (nothing to check in this API, test plan "Not encoded") |

## Criterion 13

The static null-space algorithm (V19, KuptsovParlitz2012 Eq. (55) p.13) in BigFloat(256) on the realized orbit points: forward Gram-Schmidt frame from a random orthonormal start, backward frame by `Q R = M_n^T Q_{n+1}` from another random start, P = (G-)^T G+, null vector of P(1:j-1, 1:j) by the SVD at 256 bits (`svd(...; full = true)`). It uses no back-substitution, no stored R and no C matrix. The same 256-bit tangent products feed the independent reference Ginelli (`clv_ref`, LinearAlgebra only). The 256-bit comparison uses n = 90 / (g_min Dt) (exp(-90) = 1e-39 against the prefactor 1e5) and the module in BigFloat(256) on its own BigFloat orbit, since 1e-30 is not reachable with n* = 40 / (g Dt) (e^-40 = 4e-18). n_s = 50, every sample. The zero block of F3 is compared as the span of (v_2, v_3).

| System | Float64 vs static (bound 1e-8) | Double64 vs static (1e-8) | BigFloat module vs static (1e-30) |
|--------|---------------------|---------------------|-----------------|
| standard map, 11 starts (1 + 10 chaotic), n_f = n_b = 40 | PASS, worst 7.7e-10 (start 2) | PASS, worst 2.96e-19 (start 11) | PASS, first 3 starts at n = 100: 5.4e-73, 2.5e-76, 9.8e-77 |
| coupled maps d = 4, n* = 91 | PASS 1.4e-12 | FAIL 6.0e-5 (2 n* = 182: 2.2e-17) | FAIL 1.3e-20 at n = 204 (2 n: 1.4e-43) |
| F3 orbit C, n* = 1516 (pilot, 1000 handoffs) | PASS 7.1e-12 | PASS 1.7e-15 | PASS 1.5e-33 at n = 3411 |

For the two coupled-maps FAILs the module and the 256-bit reference Ginelli give the same distance to the static algorithm to four digits (6.023e-5 for both; 1.25e-20 for both), and the distance falls to the floor when both transients are doubled. Diagnosis: transient calibration (the same effect as AC3b start 19: the n* rule has no margin on an orbit with a near tangency), not a defect of the module and not a disagreement between the algorithms. Not loosened here; the owner decides at the checkpoint. The static algorithm also cross-checks the 256-bit Ginelli at the 1e-30 level on all three systems at their long transients (the BigFloat rows, as soon as the transient is long enough). Cost ratio (R1): static over 256-bit Ginelli on the same tangent products 0.40 to 0.57 on F3 (the static pass does one backward QR sweep and one SVD per column per sample, but the 256-bit Ginelli stores and solves ns + nb factors; the tangent products, 15 to 25 s at n* on F3, dominate both) and 2.2 on the coupled maps.

## Recorded lines

- R1 (one thread, medians of 3, forward reference = the delivered discrete QR spectrum on the same handoffs; backward time is only inferred by difference, which is within timing noise): store-all total over the forward reference 1.11 (standard map k = 2, Float64), 1.13 (coupled maps k = 4, Float64), 1.09 (F3 orbit C Float64), about 1.0 in Double64 (total smaller than the forward reference in two Double64 cases, so the backward share is below the noise). The backward phase is thus a small share of the total (within noise); the ratio of 2.4 (V18, Ginelli2013 p.15) is not reproduced and a clean forward over backward ratio is not measurable from outside the module. Checkpointed over store-all: 1.40, 1.42, 1.56 (Float64), 1.50, 1.79, 1.78 (Double64), against the reported about 2 (V17). Static algorithm against 256-bit Ginelli: criterion 13 above.
- R2 (1e5 samples each, no size reduction): minimum principal angle between the unstable and stable sets. Standard map: min 2.5e-6, median 0.82, fraction below 1e-2 0.0050, below 1e-3 5.5e-4; CDF slope on 1e-4 to 1e-1 0.86. Coupled maps (sets {1,2} and {3,4}): min 1.0e-5, median 0.47, below 1e-2 0.0100, below 1e-3 8.3e-4, slope 1.07. F3 orbit C (sets {1} and {4}): min 3.1e-4, median 0.68, below 1e-2 0.020, below 1e-3 1.1e-3, slope 1.07. The distributions reach small angles with a power law of exponent near 1 and no gap at zero in these non-hyperbolic maps; the statements of the papers (O-e) are for other systems and are not tested. Pairs inside the zero block of F3 are tangent by construction (pair (2,3): 94 percent below 1e-2). Figure `clv-ginelli-angles.png`.
- R3: ratios in the table above (criteria 1, 2, 3, 5, 6, 11); worst-conditioned vector: in F3 orbit C the vector index 2 for delta (Float64 0.79 of B_v, Double64 5.1e-3) and index 3 for the exponents; in the standard map index 1 (delta) and 1 or 2 (exponent). Vector 4 of F3 is always 2 to 3 orders below the others.
- R4: transient needed for 1e-10 (ladder resolution about 15 percent), against 40 / (g Dt) and the spec's 25 / (g Dt): standard map 11 to 12 (4 starts; predictions 21 and 13, ratio to 25: 0.84 to 0.95); coupled maps 139 (3 starts; 149, 93; ratio to 25: 1.45 to 1.63); F3 chaotic starts 427 to 1417 (ratio to 25: 1.09 to 1.51), start 4 (sticky) 2112 (prediction 5652, 0.60); F3 orbit C 1417 (1516, 947; ratios 0.93 and 1.50). The spec's 25 is exceeded by up to 1.6 times for the maps and F3, the 40 of the plan is met within the noise on 12 of 13 cases but not at AC3b start 19 or in criterion 13 on the coupled maps (it needs 2 n* there). Figure left panel.
- R5: local exponents are intrinsic at common samples: tau = 5 against five tau = 1 sub-intervals, max |ln Lambda difference| 4.9e-9 (standard map), 4.9e-11 (coupled maps); F3 tau = 20 against two tau = 10 sub-intervals 3.1e-14. Householder against RGS: lambda 1.1e-14 (standard map), 1.6e-14 (F3), CLVs 1.6e-15 and 2.9e-14. The standard-map tau value 4.9e-9 is the largest: it is probably rounding amplified along the five sub-steps (growth e^{0.9 x 5}); not diagnosed further.
- R6 (O-e, reported only): test double `Henon3D` (a = 1.76, b = 0.1), 1e6 window steps, k = 3, 3000 transient handoffs each way. Window means 0.2242, 0.1874, -2.714 (start 1), 0.2244, 0.1876, -2.715 (start 2); mean over 8 random starts 0.2243 (standard error 6e-5), 0.1876 (9e-5), -2.714 (1.5e-4), against the printed 0.225, 0.188, -2.716 (KuptsovParlitz2012 p.18, as quoted by the spec): agreement to about 1e-3 (1e-3, 4e-4, 2e-3), not to the printed digits; their sum agrees with ln b (printed digits sum to -2.303). An independent plain-QR run on the same window (LinearAlgebra only) differs from `lambda_fwd` by 4.5e-7. Forward exponents sum minus ln b is 0 to printed digits; the CLV window-mean sum minus ln b is -7e-7 and -1.8e-6 raw (the boundary term of V14 (c), (ln det C_0 - ln det C_ns)/n_s), and 1.4e-16 and -5.7e-16 once that term is removed; the spec's "1e-12" holds only with the boundary term. The FPU and Henon chain statements of Ginelli2013 pp.25 to 26 are not run (systems not in the catalog).

## Size reductions

None for the package files (full sizes). R2 used the full 1e5 samples. The criterion 13 coupled maps and standard map cases use one start for coupled maps and F3, as the spec says. The 256-bit comparison on F3 and the standard map uses its own longer transients (n = 90 / (g Dt)), a change of the transient rule, not of the bound.

## Impl-notes deviations: confirmed or refuted

- Forward exponents `lambda_fwd` per step divided by h (U9): confirmed (AC10b bitwise at full size; R6 independent run 4.5e-7).
- `ClvEnsemble` export and extra provenance keys: not exercised beyond AC10a and AC14e (pass); THR reads `results` and gets bitwise equal orbits.
- Angles all principal angles ascending (U2): used by R2 (the first column is the smallest); consistent.
- Byte counts equal `clv_bytes` (288576, 24528): confirmed (AC8c). The checkpointed `bytes` is an upper bound when no frames are requested: not measured here.
- No orthogonality guard beyond the factorization status: nothing in the runs contradicted it (no health error in any group; the fault injections of AC14 pass).
- BLAS threads not set to 1: no effect seen; 1 against 8 threads bitwise (THR).
- `n_b = 0` convention (NaN last column) and the starting `C0` as sample n_s: not tested here.
- No bound or test changed at stage 04 (D123 as stage 03). The U3 pilot length and the n* rule are the open calibration points (AC3b, AC5, AC7 at the full size, criterion 13 coupled maps).

## Independent review (Sonnet, spec and code, not the notes)

No CRITICAL or HIGH.

| Id | Severity | Finding |
|----|----------|---------|
| M1 | MEDIUM | `h` not validated (`clv_ginelli.jl:25-29`): h = 0 gives lambda = [Inf, -Inf], NaN gives NaN, h < 0 flips the sign, no error (criterion 14 (d) intent) |
| M2 | MEDIUM | `h` stored as Float64 (`clv_ginelli.jl:19, 29, 246`): `lambda_fwd` divides by `T(h)` while the local exponents use `T(tau h)` rounded in Float64; Double64, tau = 20, h = 0.05: relative 5.5e-17 between the two (Double64 eps 2^-100); a BigFloat h is cut to Float64 |
| L1 | LOW | `C0` only size-checked: not upper triangular, unit columns or positive diagonal; with n_b = 0 delivered as sample n_s unnormalised (I1) |
| L2 | LOW | `angle_sets` length not checked: 1 entry gives `BoundsError` after work starts, 3 entries ignore the third |
| L3 | LOW | angles by `acos` of singular values: floor about sqrt(eps); true angles 1e-8 and 1e-10 return 0.0 (spec range 1e-6 met) |
| L4 | LOW | checkpointed `clv_bytes` uses `block`, not `min(block, N)`, and ignores `frames`: matches V17 literally, not the allocated layout |
| L5 | LOW | `x0` length not checked against the map dimension: `MethodError`, not `ArgumentError` |
| L6 | LOW | the sink receives the live recursion buffers: a writing sink corrupts later samples (documented "copy to keep") |
| L7 | LOW | no orthogonality health check of G (QRSpectrum has one); `LapackQR` with Double64 fails at the first handoff, not before work (not run) |
| L8 | LOW | the record keeps `seed`, `cseed` when `W0`, `C0` override them |
| L9 | LOW | a pivot that overflows the division is reported as "backward breakdown" even when the normalised column is representable (AC14c scenario, spec-consistent) |

Checked correct by the review: R_n from M G_n = G_{n+1} R_n and its storage index; C_n = R_n^-1 C_{n+1} with normalisation, V_n = G_n C_n; the triangular solve; Lambda = 1 / ||R^-1 c|| and its interval alignment; the n_b = 0 convention; column 1 exactly e_1; cosines; principal angles in the element type; compensated forward sums; store-all against checkpointed `isequal` over n_b in {0, 1, 5}, n_s in {1, 4, 9}, h_b in {1, 2, 3, 4, 7, 20}, 0 mismatches; block boundaries; health-error step, orbit, column fields; argument errors before work; concurrency (distinct per-orbit buffers, lowest failing orbit rethrown). Not checked by the review: allocation and type inference, `lambda_fwd` bitwise for tau > 1 and non-unit h, thread runs, oracle accuracy, the other orthonormalizers by run.

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number; figures `clv-ginelli-convergence.png`, `clv-ginelli-angles.png` |
| Review closed | no CRITICAL or HIGH finding; M1, M2 fixed under D124 |
| Criteria passing | 11 of 14 clean at the full size; 3 with a failing line, none traced to a defect of the module: 3 (b) (F3 start 19 at n* = 40 / (g_min Dt), 3.3e-8; 6.8e-14 at 3 n*), 5 and 7 (error at the full size from the 300-handoff pilot: lambda_1 < lambda_2 on start 4, negative n*; PASS with a 1000-handoff pilot), 13 coupled maps (Double64 6.0e-5 at n* = 91, 2.2e-17 at 2 n*; 256 bits 1.25e-20 against 1e-30, 1.4e-43 at 2 n; the reference Ginelli gives the same distance) |

## Amendment D124 (owner: "Amend + fix M1-M2")

Spec section "Amendment D124", test plan section of the same name, decision log D124. Test side: pilot spectrum 1000 handoffs (was 300), transient constant c = 120 in n* = c / (g_min Dt) (was 40) for criteria 3 (b) and 13, the 1e-30 comparison of criterion 13 at c = 180 (was 90); AC3c and AC11 use the module's Dt (tau h in the element type); AC5 (a) bound plus 64 eps (1 / |omega_n| + 1 / |omega_{n+1}|), found by the FULL1 rerun (F3 chaotic start 6, |omega(v_2, v_3)| 3.6e-6: ratio 7.99 without the term, 0.049 with it). Code: M1 (`h` finite and positive, `ArgumentError`), M2 (`h` kept as given, Dt = T(tau) T(h)), with the new testset "D124 step length h".

Reruns (`clv-ginelli-d124-verify-log.txt`): package clv suite 6258 pass, 0 fail; Aqua, JET, layering 14 pass; FULL1 44872 pass (AC3b 0 fail, AC5 worst (a) 0.049), FULL2 17563 pass (AC5, AC7 at the full size), FULL3 151 pass, 0 fail. Criterion 13: standard map 11 starts worst 2.3e-11 (Float64), 1.3e-25 (Double64), 256-bit module against static at most 5.4e-73; coupled maps 7.6e-11 (Float64), 8.1e-28 (Double64), 4.0e-46 at 256 bits; F3 orbit C 6.5e-12 (Float64), 5.4e-30 (Double64), 1.1e-72 at 256 bits. FULL1p, FULL2p, AC3bd, Oe, R1 to R5, THR, FIG not rerun (FULL1p, FULL2p are superseded by the helper change; the others do not use the changed lines, except R4's comparison constant, which stays the 40 and 25 of the first run).

Criteria passing after D124: 14 of 14. Open: review L1 to L9 (LOW).
