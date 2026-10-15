# Verification: lyap-householder-qr

Spec: `../../01-spec/output/lyap-householder-qr-spec.md` (criteria 1 to 17, Amendment D35). Test plan: `../../02-tests/output/lyap-householder-qr-test-plan.md` (Amendment D38). Package tests: `test/lyap-householder-qr/`. Oracle script: `verify/lyap-householder-qr.jl`, run as `julia +release -t 2 --project=verify verify/lyap-householder-qr.jl` from the repo root; full output in `lyap-householder-qr-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`).

Run facts. Julia 1.13.1, 2 threads, base seed 20261012 with offset 44000 + 50 ac + sub (the package suite uses 24000 to 24999, so every sample is fresh). Working tree of branch `m2-p1-modules` (clean at start). Wall time 9827 s (2 h 44 min) on a machine shared with other agents; sections: full cross-product 1549 s, AC7 5973 s (almost all of it the 256-bit references at n = 10^6), D35 reproduction 655 s, figures 877 s (first CairoMakie call), the rest under 150 s each. Runtimes are informative only. Values are "worst ratio of error to bound" unless stated; pass at ratio <= 1.

Sample sizes: the measurement functions reuse `hh_helpers.jl` and the definitions of the `hh_*.jl` test files (loaded without their testsets): `hh_randW`, `hh_ref_qr`, `hh_eval`, `hh_ref_std`, `hh_orbit`, `hh_Psi`, `hh_theta`, `hh_theta_psi`, `hh_BS`, `hh_kappa`, `hh_fibbig`, `hh_matpow`, `hh_det`, `hh_regular_starts`, `hh_se` and the rest. A definition error in a shared helper would therefore not be caught by the number sections; the independent pieces are the 256-bit and 2048-bit references, ChaosTools, and the Gram-Schmidt, Benettin, SVD and two-trajectory estimators of sibling modules. The number sections restore the spec sizes that D38 and the package suite reduced: 200 W per shape in AC1 to AC3 (package: 100), tau in {1, 2, 4} in AC7 (cross-product: {1, 2}), 50 Double64 starts in AC7 (package: 10), K = 10 in AC8b, chunk 64 in AC9, all four variants in every section, and Double64 where the spec says T in {Float64, Double64}. AC8 is Float64 only (D35 evidence limit: the unit is not measured for Double64).

## Full cross-product (D38)

`CROSS | full cross-product (D38) assertions pass / fail / error / broken | 2035 / 0 / 0 / 0 | PASS`, 25 min 47 s testset time (1549 s with includes). Largest testsets: AC6 4m46s, AC1 3m47s, AC12 3m07s, AC7 2m46s, AC9 1m38s, AC11 1m30s, AC4 1m29s. The per-testset Test summary table was printed to stdout only (not in the log; the log holds the CROSS line).

## Evidence table

"Where": T = package test, X = cross-product section, V = number section of `verify/lyap-householder-qr.jl`. Every value is from the log.

| AC | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1a | column residual / (8 d eps \|\|w\|\|), worst over 7 shapes x 200 W x variants, Float64 / Double64 | 0.4739 / 0.02231 | 1 | PASS | T, X, V |
| 1b | max\|Q'Q - I\| / (8 k eps), per (d,k) and Float64 (worst over variants): (2,1) 0.9082, (2,2) 0.5304, (3,2) 0.6621, (4,2) 0.4528, (4,4) 0.2407, (8,3) 0.4129, (8,8) 0.1461 | worst 0.9082 (V-H and V-A at (2,1)) | 1 | PASS (marginal) | V |
| 1b | same, Double64 per (d,k): 0.0496, 0.0241, 0.0358, 0.0212, 0.0124, 0.0126, 0.0083 | worst 0.0496 | 1 | PASS | V |
| 1c, 1e | failures of status, shapes, exact zeros below diagonal, R_ii > 0, counts r(d,k) and g(d,k), Float64 / Double64 | 0 / 0 | 0 | PASS | T, X, V |
| 1d | \|R_ii - R_ii^ref\| / (8 d eps kappa \|\|w\|\|), Float64 / Double64 | 0.1574 / 0.004574 | 1 | PASS | T, X, V |
| 2a | max\|Q'Q - I\| / (8 k eps), theta in {2^-10, 2^-20, 2^-40}, all variants: Float64 loss 1.649 eps at every theta (no growth with 1/theta); Double64 0.0345 eps | 0.1031 / 0.002157 | 1 | PASS | T, X, V |
| 2b | \|R_22 - ref\| / (8 d eps \|\|w_2\|\|): Float64 error 0.0082, 0.303, 0.209 eps \|\|w_2\|\| for the three theta; Double64 0.0056, 0.0053, 0.0120 eps | 0.009471 / 0.0003752 | 1 | PASS | T, X, V |
| 3 | variant agreement R_ii, all pairs, Float64 / Double64 | 0.08108 / 0.002257 | 1 | PASS | T, X, V |
| 3 | variant agreement Q, all pairs, Float64 / Double64 | 0.125 / 0.006348 | 1 | PASS | T, X, V |
| 4a | shear, n = 1e4: variants with any nonzero S_i, lambda or \|W - I\|, Float64 (4 variants) / Double64 (3) | 0 / 0 | 0 | PASS (bitwise) | T, X, V |
| 4b | \|S_1 - (1/2) ln(n^2+1)\| / (8 n eps), n = 10, 1e3, 1e4: Float64 0.0814, 0.0076, 0.0491; Double64 0.00083, 0.0034, 0.0032 | 0.0814 / 0.0034 | 1 | PASS | T, X, V |
| 4b | max S_1/n - ln(n)/n at n = 1e4, Float64 / Double64 | 5.001e-13 / 5.0e-13 | 1e-12 | PASS (margin 2) | X, V |
| 5a | \|S_1 - ln\|\|A^n e_1\|\|\| / B_S, tau in {1,2,5}, n 10 to 1e4; tau = 1 max 1.09 n eps (Float64), 0.034 n eps (Double64); spec observed 2.03 | 0.1362 / 0.00366 | 1 | PASS | T, X, V |
| 5b | \|S_1 + S_2\| / (n theta_sum) | 0.01157 / 0.0005421 | 1 | PASS | T, X, V |
| 5c | \|lambda_1 - E10 form\| / (B_S/n) | 0.1464 / 0.003138 | 1 | PASS | T, X, V |
| 5d | \|lambda_1 + lambda_2\| / theta_sum | 0.01059 / 0.0005515 | 1 | PASS | T, X, V |
| 6a | \|C_i - ln vol_i (2048 bit)\| / (i n theta_sum), N = 2, 3, A and A_S, n = 500, Float64 / Double64 | 0.01853 / 0.0004904 | 1 | PASS | T, X, V |
| 6b | \|C_i/n - Lambda_i - B_i/n\| / (i theta_sum), normal A | 0.02064 / 0.0005338 | 1 | PASS | T, X, V |
| 6c | \|sum lambda_i\| / theta_sum | 0.0614 / 0.001713 | 1 | PASS | T, X, V |
| 7a | max_i \|lambda_i - ref256\| / theta_psi, n = 1e3, 5 K, tau in {1,2,4}, 4 variants, 50 starts: Float64 worst 0.02592 (K = 2), Double64 worst 0.002439 (K = 5) | 0.02592 / 0.002439 | 1 | PASS | T, X, V |
| 7a | same at n = 10^6, Float64, first 2 of 50 starts per K (reference 163 to 266 s per K): worst | 0.002688 | 1 | PASS | V |
| 7b | pairs of variants / (theta_psi^A + theta_psi^B): n = 1e3 Float64 0.01344, Double64 0.0004668; n = 10^6 Float64, full 50 starts 0.009818 | 0.01344 | 1 | PASS | T, X, V |
| 7c | \|lambda_1 + lambda_2\| / theta_sum: n = 1e3 Float64 0.003307; n = 10^6 0.0004530 | 0.003307 | 1 | PASS | T, X, V |
| 7c | measured c' of E11 (against the diagnostics reference 4, P8): n = 1e3 worst 0.0529; n = 10^6 worst 0.007248; spec observed <= 0.5 | 0.0529 | < 4 | PASS (margin 75) | V |
| 8a | \|lambda_1^QR - X_Benettin\| / (16 eps Psi_n/n), K in {0.5,2,5,2pi}, 50 orbits, n = 1e4, tau in {1,4}: V-A 0.0547, V-G 0, V-H 0.0547, LAPACK 0.0300 | 0.0547 | 1 | PASS | T, X, V |
| 8b | max_i \|lambda_i^QR - lambda_i^GS\| / (theta_psi^QR + theta_psi^GS), 5 K, tau in {1,2,4}; CGS, MGS, RGS | 0.005522 | 1 | PASS | T, X, V |
| 8b | observed orthogonality loss max\|Q'Q - I\| / (k eps): QR variants 5.5 (bound 8), CGS/MGS/RGS 1.7e8 | - | - | reported | V |
| 8c | SVD V7 sandwich, \|mid - centre\| / half-width (<= 1 inside), 4 K, tau in {1,2,4}, 4 variants | 0.9996 | 1 | PASS (marginal) | T, X, V |
| 8c | \|lambda_1^SVD + lambda_2^SVD\| / theta_sum | 0.001228 | 1 | PASS | T, X, V |
| 8d | ensemble K = 2pi wrapped, 1000 starts, n = 1e4: \|mean diff\| / (4 SE); means 1.17, 1.17 | 3.557e-08 | 1 | PASS | X, V |
| 8d | pathwise K = 0.5, 50 starts, n = 1e3, d0 = 1e-6 | 0.002341 | 1 | PASS | T, X, V |
| 9 | 1 thread vs 4 threads, child processes, 2048 orbits, n = 1e3: hash of x, W, finals; HouseholderQR, GivensQR, AugmentedQR (Float64, Double64), LapackQR (Float64) | all equal | bitwise | PASS | `verify/threads.jl` |
| 9 | non-identical results of `propagate` vs serial `propagate_orbit`, 4096 orbits, n = 1e3, chunks 64, 1000, default: V-H, V-A, V-G; library path | 0, 0, 0; 0 | 0 | PASS (bitwise, 2 threads) | T, X, V |
| 10 | max bytes per `reducer_step` (kernel, Float64 and Double64, 4 shapes) and per `propagate_orbit` n = 100 (30 measurements); library path | 0; 0 | 0 | PASS | T, X, V |
| 11a,b | health cases failing (of 28), kind, orbit 3, step 15, column, value 0, accumulators unchanged | 0 | 0 | PASS | T, X, V |
| 11c | W = s Q_test, s = 1e+-200: max\|R_ii/s - 1\| / (8 d eps), Float64 / Double64 | 0.125 / 0.0009569 | 1 | PASS | T, X, V |
| 11c | max\|Q - Q_test\| / (16 d eps); \|lambda - ln s\| ratio (Float64 / Double64) | 0.09375 / 0.001343; 0 / 0 | 1 | PASS | T, X, V |
| 12 | state length 4 vs 2, lambda, G, final Q differing bitwise (2 T x 2 k x variants) | 0 | 0 | PASS | T, X, V |
| 13 | \|lambda_1 + lambda_2\| / theta_sum(tau), tau in {1,2,4}, Float64 / Double64 | 0.000614 / 0.00008634 | 1 | PASS | T, X, V |
| 13 | \|lambda_i(tau) - lambda_i(1)\| / (theta_sum(tau) + theta_sum(1)) | 0.0001144 / 0.00002052 | 1 | PASS | T, X, V |
| 14a | HouseholderQR, GivensQR, AugmentedQR, regular K = 0.5 (S7: 1024 of 3072 regular, 17 replaced), n = 4e4, tau = 10, against a Double64 reference: orbits above 1e-6 CPU / GPU; max CPU-GPU gap, lambda_1, lambda_2 | 1 / 1; 8.98e-08 (all three) | count diff <= 3; gap <= 2e-6 | PASS | `verify/gpu.jl` |
| 14a | same three, chaotic K = 10, Float64, 1024 orbits, n = 1e4, tau = 2: max per-orbit gap / sigma_FT; difference of means / SE | 4.992; 0.6086 (all three) | 6; 4 | PASS | `verify/gpu.jl` |
| 14a | max over orbits of \|lambda_1 + lambda_2\| / theta_sum, CPU; GPU: K = 0.5 H 4.9e-4; 5.8e-4, G 3.0e-4; 2.3e-4, A as H; K = 10 H 2.4e-4; 2.2e-4, G 1.3e-4; 8.6e-5, A as H | max 5.8e-4 | 1 | PASS | `verify/gpu.jl` |
| 14a | Double64 kernel | does not compile for the device | compile or reason | FALLBACK(kernel does not compile for the element type) | `verify/gpu.jl` |
| 14b | library path QRSpectrum(LapackQR), gpu and auto request on the GTX 1050 machine: backend; device; reason | cpu; cpu; "reducer is CPU-only" (both) | exact | PASS | `verify/gpu.jl` |
| 14c | K15 record, kernel path QRSpectrum(HouseholderQR), 10^6 orbits, d = k = 2, N = 1e4, tau = 2, K = 2 pi: wall s GPU (GTX 1050, CUDA runtime 12.9.0, 202 chunks); CPU 4 threads; CPU 1 thread (skylake); gate GPU <= CPU | 311.3; 1971; 8109; true (GPU 6.3x the 4-thread CPU) | recorded | INFO (gate met) | `verify/gpu.jl` (`gpu-hh-perf-verify-log.txt`) |
| 15 | failed exact-string and provenance checks (of 34) | 0 | 0 | PASS | T, X, V |
| 16 | H1 versus H2 | - | - | EXCLUDED (D30) | |
| 17a | regular K = 0.5, n = 1e5, \|lambda_i - ChaosTools\|, i = 1, 2, four variants (lambda = +-8.826e-05; ln n / n = 1.151e-04) | 7.5e-15 (V-H, V-A), 3.1e-15, 2.7e-15 | 1e-4 | PASS | V |
| 17b | chaotic K = 2pi, 1000 starts, n = 1e4: \|mean lambda_1 - mean ChaosTools\| / (4 SE), four variants; means 1.171, 1.171 | <= 3.114e-14 | 1 | PASS | V |

Recorded lines (no threshold), in the log:
- R4: ensemble mean lambda_1, 1000 starts, n = 1e4: K = 2pi 1.169 +- 0.0021 (benchmark row, single orbit T = 1e9: 1.172); K = 10 1.620 +- 0.0003 (ln(K/2) = 1.609).
- R5: tau = 8: |lambda_1 + lambda_2| = 2.1e-7, theta_sum = 3.4e-3, G_max = 7.7e12; tau = 10: 1.6e-4, theta_sum 3.09, G_max 8.7e15; tau = 16: the reducer raises "orthonormalisation breakdown" (exact zero R_22, step 16) for V-H. The same happens at tau = 10 in a different start (step 2410 of 1e4, separate probe). So above tau_max the Householder family also reaches rank loss in Float64, only later than Gram-Schmidt; lambda_1 at tau = 8 and 10 equals the tau = 1 value to the last bit on this orbit.
- R6: `orthonormalize` time, loaded machine, d = k = 2: G 0.18 us (0.60x of V-H), CGS 0.28, MGS 0.28, V-H 0.30, RGS 0.32, V-A 0.34, LAPACK 2.3 us (7.7x). d = k = 8: CGS 3.6 us, MGS 4.5, V-H 6.4, RGS 6.9, V-A 7.2, V-G 8.4, LAPACK 15.8 us (2.5x). H1 not built, so no H1 : H2 ratio (D30).

## Oracle comparison

- BigFloat reference (spec Oracle): AC1 (d) against the 256-bit factorization of the same stored W (0.157 Float64 worst); AC1 (a), (b) evaluated at 256 bits; AC5 against exact BigInt Fibonacci products; AC6 against 2048-bit Gram determinants and explicit eigenvectors; AC7 against a 256-bit per-step MGS run on the same orbit points (n = 1e3 for all 50 starts and both types; n = 10^6 for 2 starts per K, the subset forced by the reference cost of 163 to 266 s per orbit pair; (b) and (c) at n = 10^6 are on the full panel of 50 starts).
- Across the four variants: AC3 (R_ii 0.081, Q 0.125 of the bound), AC7(b) and AC13. V-H, V-A and LAPACK agree to a ratio that is bitwise in R_ii at (2,1) (AC3 R_ii ratio 0 there).
- Gram-Schmidt QRSpectrum: AC8(b) 0.0055 of the bound (CGS, MGS, RGS each), with an orthogonality loss of 5.5 k eps for the QR variants against 1.7e8 k eps for the Gram-Schmidt variants on the same orbits.
- ChaosTools 3.5.4 `lyapunovspectrum(tands, N; u0)` on a `TangentDynamicalSystem` built with a hand-coded Jacobian and `Q0 = I` (installed API read under `~/.julia/packages/ChaosTools/0SbsT/src/chaosdetection/lyapunovs/lyapunovspectrum.jl`; the returned vector is per column, not sorted by the code). Lifted standard-map rule `p' = p + K sin x, x' = x + p'` (lyap-benettin Oracle). AC17 differences are 3e-15 to 7e-15 pathwise, i.e. rounding level, because ChaosTools evaluates the same floating-point operations on the same lifted orbit; this validates the convention, the sign and column order and the exponent formula, not an independent orbit sample. Lyapunov exponents from ChaosTools both with i = 1 and i = 2 agree to 7.5e-15 against the spec bound 1e-4.
- Also checked: lyap-benettin (AC8a), lyap-svd-stabilized (AC8c) and lyap-two-trajectory (AC8d), all within bounds.

D35 reproduction (6000 chaotic orbits, K in {2, 5}, 3000 starts each, n = 1e3, tau = 1, V-H kernel, 256-bit reference; Figure 2): max_i error / theta_psi = 0.0949 (spec: 0.076); |lambda_1 + lambda_2| / theta_sum = 0.001812 (spec: 0.0035); |lambda_1 - ref| / (eps Psi_n/n) = 1.203 (spec: 0.81; the bound is 8, so 0.150 of the bound); orbits whose error exceeds theta_sum: 5 (spec: 5 of 6000). The spec's figure 0.81 is not reproduced (1.20 on a fresh sample, still 6.6 times inside B_S^Psi); the spec's quoted ratios were taken against a former bound, see the spec note.

## Margins and marginal criteria

- AC1(b) Float64 (2,1): 0.9082 of 8 k eps (realized 7.27 eps, spec observed 2 to 8.5 eps at d = 2) for V-H and V-A; this is the tightest factor criterion and is the bound's constant, not a defect. (2,2) 0.53, (3,2) 0.66.
- AC8(c): SVD V7 sandwich position 0.9996 of the half-width: one orbit sits at the edge of the sandwich of lyap-svd-stabilized V7 (uses theta_psi^QR + theta_psi^SVD); it passes, margin 4e-4. The sandwich is a property of the pair of estimators, not of the QR variants; flag for the review of lyap-svd-stabilized.
- AC4(b): S_1/n - ln(n)/n = 5.0e-13 against 1e-12 (margin 2; the quantity is O(n^-3) away from the exact limit and the bound is the spec's).
- AC7 Double64 K = 5: lambda_1 error 0.7665 eps Psi/n (bound 8 eps Psi/n, spec observed 0.81): the same size as the spec's float observation, in Double64 the trig floor of D35 applies; ratio to theta_psi is 0.0024.
- ALG2 (Algorithm item 2, below): at W = A^20 Q_0 the exponent sum is 0.95 (G_max 2.0e16, theta_sum 71.8), so the criterion holds vacuously there; this is not a pass of accuracy.
- No FAIL on any line of the final log. One earlier FAIL in a smoke run was my script (the D35 lambda_1 error was ratioed to 1 instead of the 8 of B_S^Psi); corrected in the script before the final run, src, tests and tolerances untouched.

## Algorithm item 2 (rank loss; GS AC7 case A^20 in Float64)

At `W = A^20 Q_0` (one handoff, A = [[2,1],[1,1]], exact integer entries, G_max 2.02e16): CGS, MGS and RGS raise the health error "orthonormalisation breakdown" at step 1; V-H, V-A, V-G and LAPACK continue, with max|Q'Q - I| = 1 k eps (V-G 0), |S_1 - ln||A^20 q_1||| = 5.86 eps for all four, and |lambda_1 + lambda_2| = 0.9526 (R_22 carries an O(1) relative error, the conditioning of E10; theta_sum = 71.75 so the sum identity holds only vacuously). The spec statement (Algorithm 2) is confirmed: the Householder family does not raise here, and gives an orthogonal Q but a meaningless R_22. It does raise on exact zero pivots (R5, tau = 16).

## Deviations from the implementation notes, checked

- `G_w` rescale in `orthonormalize.jl` (Gram-Schmidt-owned code): for 40000 random triangles (k = 2, 4) with entries in [1e-100, 1e100] the new `window_Gw` equals a copy of the pre-rescale code bitwise in every draw (0 differences); at `W = s Q B`, s = 1e+-200, all seven reducers give a finite `G_max` (relative difference to the unscaled case 2.9e-16). The change is behaviour-neutral in range and needed outside it; still an owner item (it sits in another slug).
- `LapackQR` BLAS threads: with 2 Julia threads `reducer_init` leaves 1 BLAS thread (the default was already 1 in this process, so the line does not discriminate the setting; the test is weak).
- `LapackQR` pure `orthonormalize` allocates (2743 bytes per 4x4 call) as documented; `reducer_step` allocates 0 bytes (AC10, 4 shapes).
- V-A and V-G form Q from transposed rows of S: covered by AC3 (Q ratios 0.125 at worst, Float64).
- Counted reflections and rotations equal E9 (AC1(e), 0 failures over 7 shapes x 200 W x 2 T).
- Generic BigFloat(128) `HouseholderQR` on a 4x2 matrix: status OK, max|Q'Q - I| = 1.8e-38.
- Health kind: the spec text (AC11) says "deviation matrix"; the code raises "orthonormalisation breakdown" (D34, GS Algorithm item 5); all 28 health cases match the D34 kind, orbit, step 15, column and value 0.
- H1 (treppen path) is not built and AC16 is excluded by D30; no H1 : H2 timing exists.
- Kernel Double64 8x8 timing on the device and GPU compile are not measured here (AC14).

## Figures

- `lyap-householder-qr-factor-error.png`: worst ratio to the bound over 200 W of (a) the column residual, (b) max|Q'Q - I|, (d) the R_ii error against the 256-bit reference, per (d,k), four variants, Float64 (circles) and Double64 (triangles; V-G, V-H, V-A only). All points lie below the bound line (1) with the Float64 (2,1) V-H and V-A at 0.91; Double64 sits one to two decades lower (effective eps 2^-100). The legend hides some low Double64 points in the left panel.
- `lyap-householder-qr-theta-psi.png`: 6000 chaotic orbits (K = 2 red, K = 5 blue), error max_i|lambda_i - ref256| against theta_psi(n) of the orbit with the line error = theta_psi (left), and |lambda_1 + lambda_2| against theta_sum (right). The data cluster two to four decades under the line; the few orbits with large theta_psi are the island-sticking orbits that motivate the Psi term (5 orbits exceed theta_sum, all inside theta_psi).

## Runtime

Whole script 9827 s with 2 threads on a loaded machine (full cross-product 1549 s; AC7 5973 s of which the n = 10^6 reference 1.1 h; D35 655 s; first CairoMakie call 877 s). The `HV_ONLY` switch (`cross, AC1-3, AC4-6, AC7, AC8, AC9-12, AC13, AC15, alg2, dev, D35, AC17, R, figure`) runs the number sections without the cross-product; `HV_N6`, `HV_N6_REF`, `HV_AC7_STARTS`, `HV_D35_PER_K` shrink the long parts. A shorter final run is possible with HV_N6 = 100000 if the 10^6 horizon is not required.

## Review findings

Independent review: fresh Sonnet subagent given the spec and the code, not the implementation notes. No CRITICAL and no HIGH finding. Disposition: closed (evidence named), fix (proposed change), record (decision row or spec amendment), owner (checkpoint question).

| Severity | Finding | Location | Disposition |
|---|---|---|---|
| MEDIUM | LapackState rebuilt per handoff, `info` Refs may allocate; library-path allocation not measured | `src/Lyapunov/lapack_qr.jl` | partly closed: verify measured `orthonormalize` 2743 B, `reducer_step` 0 B; record |
| MEDIUM | non-isbits state gives rbytes = 0 in the chunk budget; LapackWork per orbit not counted | `src/Tangent/propagate.jl:123` | fix: count the work size, or record (small d only) |
| MEDIUM | AC11(b) identical columns raise only for exactly representable data; spec Algorithm 2 (exact zero only) against the AC11(b) wording; kind name "deviation matrix" against D34 "orthonormalisation breakdown" | spec AC11, Algorithm 2 | record: spec fix |
| LOW | `BLAS.set_num_threads(1)` inside the per-orbit `reducer_init`: global, never restored | `src/Lyapunov/lapack_qr.jl:63` | owner: move to call site or record |
| LOW | LAPACK `info` ignored | `src/Lyapunov/lapack_qr.jl:28, 37` | fix: map info != 0 to a health code |
| LOW | AC1(e) reflection count vacuous for LapackQR | `test/lyap-householder-qr` | record |
| LOW | AC9 thread count not varied; library path at one chunk size | `test/lyap-householder-qr` | closed: `verify/threads.jl`, child processes `-t 1` and `-t 4`, outputs bitwise equal (`threads-verify-log.txt`) for HouseholderQR, GivensQR, AugmentedQR (both types) and LapackQR (Float64) |
| LOW | orthogonality health path untested for the HH family; no d > 8 library run | `test/lyap-householder-qr` | fix: test cases |
| LOW | AC7/AC8 Psi_n at tau 2, 4 and Double64 against theta_psi: evidence only (D35 limits) | `verify/lyap-householder-qr.jl` | record |
| LOW | iteration_path always "diagonal"; H1 not implemented (P5/D30 optional) | `src/Lyapunov` | owner: record the decision |
| LOW | "reducer is CPU-only" keyed on `!isbitstype`; a `cpu_only` trait is clearer | `src/Tangent/propagate.jl` | no change before M3 |
