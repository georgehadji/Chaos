# Verification: gali

Spec: `../../01-spec/output/gali-spec.md` (criteria 1 to 21, R1 to R5, Amendments D107 a to l and D108 a to k). Code: `src/Indicators/gali.jl` (stage 03 commit f5e4469, plus the D108j health rule of the rework 2026-11-05: `_first_bad`, `reducer_error`, `gali_values`). Package tests: `test/gali/` (2724 of 2724 at the package sizes after the rework; 2704 of 2704 at stage 03). Oracle script: `verify/gali.jl` (one Julia process per group: `GA_ONLY=AC1-4`, `AC5-8`, `AC9-12`, `AC13-16`, `AC17-20`, `D108`, `DIAG-AC9`, `DIAG-AC11`, `DIAG-AC12`, `DIAG-AC13`, `R3-AC2`, `R1`, `R2`, `R4`, `AC19`, `AC20t`, `figure`; `julia --project=verify -t 8 verify/gali.jl`). GPU: `verify/gpu.jl` with `GPU_ONLY=gali`. Raw logs: `gali-verify-log.txt`, `gpu-gali-verify-log.txt` (each has the run before D108 first and a headed section "D108 rerun (2026-11-05)" appended). Figure: `gali-laws.png`.

Run facts. Julia 1.13.2, 3 to 8 threads per group (the groups ran in parallel on one machine, so timings are loaded), base seed 20261012 + 40000, verify-only seeds 40900 to 40999 (`VSEED`, GPU 40950). The script sets `GA_FULL_RUN = true` and includes the package test files, so the assertions and tolerances are the package ones and the sizes marked D38 there take the spec sizes (1000 sets per pair in 1 and 3, 1e4 pairs in 4, 20 starts in 5, 8, 9, 10, 13, 15, 18, 20, 5 frames in 11 and 12, n = 1e4 for the regular cases of 9, N = 10 in 11, N = 5 in 9 and 12, 200 orbits at n = 1000 in 20). The `@info` lines of the tests (realised ratios, slopes, counts) are copied into the log. Device: NVIDIA GeForce RTX 5060 Laptop GPU, capability 12.0, CUDA.jl 6.4.2.

Result in one line (after the rework D108, 2026-11-05): criteria 1 to 21 PASS at full size, except 21(b) chaotic class, which is reported in units of the plain bound sum (D108i, INFO) and is open for the orbit-divergence term. Before D108 (checkpoint of 2026-11-04): criteria 1 to 8, 10, 16 to 20 PASS; 9, 11, 12, 13, 14, 15(b) and 21(b, chaotic class) FAIL as specified at full size. None of those failures was an error of the reducer against the 256-bit reference: Amendment D108 corrected the bounds, data and references of the criteria (see "Findings" and the rework section), and the one code change is the subnormal health rule D108j (review finding 1), tests first.

Suite totals at full size (assertions), before D108: AC1-4 18 016 pass, 0 fail; AC5-8 917 pass, 0 fail; AC9-12 630 pass, 98 fail (AC9 72, AC11 16, AC12a 10); AC13-16 723 pass, 17 fail (AC13 12, AC14 1, AC15b 4); AC17-20 1527 pass, 0 fail. After D108 (rerun of the changed groups): AC9-12 728 pass, 0 fail; AC13-16 724 pass, 0 fail; D108 20 pass, 0 fail; AC1-4, AC5-8 and AC17-20 not rerun at full size (their package-size testsets pass; the code change touches only the subnormal path). The first run of AC13-16 and AC17-20 had undefined-helper errors (`ga_slopes`, `ga_check_orbit` live in earlier test files); the script now loads the definitions of the earlier files first, and the logs are those of the rerun.

## Evidence table

"T" = package test file at full size run by the script, "V" = extra section of `verify/gali.jl`, "G" = `verify/gpu.jl`. Ratios are error over bound (pass at most 1).

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1 | (a) four routes agree 1e-60; (b) module evaluation, worst error / (2 B_meas), Float64 / Double64, 1000 sets per pair and type | 0.242 / 0.0089 | 1 | PASS | T |
| 1 (c) | negative control, skip-normalisation double fails (b) | fails | fails | PASS | T |
| 2 | worst error / (3 B_meas), constructed Gram matrices, Float64 / Double64 | 0.083 / 0.0030 | 1 | PASS | T, V |
| 3 | range, nesting, sign/permutation/orthogonal invariance, dependent sets, scaling 2^+-498 | all hold | 2 B_meas | PASS | T |
| 4 | k = 2 against SALI, worst ratios (a) to (d), Float64 | 0.25, 0.21, 0.15, 0 | 1 | PASS | T |
| 4 | same, Double64 | 0.0091, 0.0078, 0.0056, 0 | 1 | PASS | T |
| 5 | determinant identity: (a) standard map worst ratio; (c) Y3 and Y4 against the 256-bit product | 0.097; 0.0010 | 1 | PASS | T |
| 5 (d) | negative control (k < d) differs by a factor 10 | holds | 90 percent | PASS | T |
| 6 | floor: m_floor Float64 k = 2, 3, 4; Double64 k = 2, 3, 4; sharper bound violated | 17, 16, 15; 33, 33, 31; 0 of 4000 each | m_floor >= log10(1/(2 B_meas)) - 1; 99 percent | PASS | T |
| 6 (d) | Gram-determinant control fails | 150 of 150 | at least half | PASS | T |
| 7 | argument errors, non-finite, frame check, seeded frames: min G_K(0) / (10 B_meas) Float64 / Double64 | 5.6e13 / 1.6e28 | accepted | PASS | T |
| 8 | nesting bitwise (a), sums bitwise equal to sali (b), K = d = 2 equals sali (c), Y4 (d): bitwise comparisons; G_2 vs B4 conversion worst ratio Float64 / Double64 | 378 of 378; 0.042 / 0.0026; 0 unequal | bitwise; 1 | PASS | T |
| 9 (a) | direct form vs 256-bit reference at resolved handoffs (G_ref >= 1e-70, D108c), worst error / B_G^k: Float64 K = 4 chaotic / regular; Float64 K = 6, 10; Double64 K = 4 chaotic (3546 / 6300; 185, 231; 3663 resolved order-handoffs, the rest skipped and counted) | 0.72 / 0.072; 0.032, 0.076; 0.072 | 1 | PASS (unchanged) | V |
| 9 (a) | the suite's own check at full size (AC9-12, 728 assertions with AC10 to 12) | before D108: 40 assertions NaN (reference exactly 0 beyond 2^-256), 2 more fail; after D108c: 0 fail | 1 | PASS | T |
| 9 (b) | log form, B_L = 32 eps (1 + a_max) sum Psi (D108a): Float64 chaotic Y3, 20 starts, n = 1000 (before D108: B_L = 32 eps sum Psi, worst 60.1) | worst 0.0047 | 1 | PASS | V, T |
| 9 (b) | same, other classes: Float64 regular and orbits C, D; Float64 K = 6, 10; Double64 K = 4 chaotic / regular, K = 6, 10 | 0.0009; 0.00011, 0.000024; 0.000082 / 0.000044, 0.0000027, 0.0000030 | 1 | PASS | V |
| 9 (b) | realised multiplier: worst error / (eps sum Psi), Float64 chaotic K = 4; same over (1 + a_max) (the evidence of D108a) | 1923; 0.15 | info | INFO | V |
| 9 (c) | direct form vs exp(L) at resolved handoffs (G_d,ref >= max(10 B_meas, 1e-70)) | at most 7e-4 | 1 | PASS | V |
| 9 (d) | floor control (D108b): G_ref < B_meas and the value above 2 B_meas + B_G^k (before D108: above 2 B_meas, 3 order-handoffs on Y4 N = 5, k = 2, n = 690 to 730, each inside B_G^k) | 0 violations in every class and type | 0 | PASS | V |
| 9 (e) | negative control (partial-sum log form) fails | fails | 10 B_L | PASS | T |
| 10 | (a) tau in {1, 2, 5, 10}: worst G / L ratio; (b) K = 0 shear to n = 1e6: B_G^k holds, G_3 n and G_4 n^2 | 0.042 / 0.0061; 0.9999999950 and 0.9999999900 at 1e4, 1.0000000000012 and 1.0000000000025 at 1e6 | 1; reported | PASS | T |
| 11 | chaotic law at hyperbolic fixed points, gated window (100, 200) for d <= 10 (D108d), Y3 and N <= 5, BigFloat and Double64, 5 frames: orders within 0.5 / outside / unresolved | 108 / 0 / 172 | 0 outside | PASS | V |
| 11 | same, N = 10 (d = 20), gated window (400, 800) (D108d; before D108 the gated window (100, 200) gave 22 / 28 / 0, FAIL) | 20 / 0 / 30 | 0 outside | PASS | V, T (AC11 0 fail) |
| 11 | same, N = 10, window (100, 200) (not gated after D108d, reported) | 22 / 28 / 0 | reported | INFO | V |
| 11 | same, N <= 5, window (400, 800) (not the gated window) | 36 / 4 / 240 | reported | INFO | V |
| 11 | direct orders failing in the early window (40, 80), all systems | 101 | reported | INFO | T |
| 12 (a) | stable fixed points: Y3 slopes of ln G_k against ln n over [100, 2000], max over 5 frames | 0.184 | 0.2 | PASS | T |
| 12 (a) | Y4 N = 2, 3 (eps = 0.3) and N = 5 (eps = 0.2, D108e), 5 frames: max slope | 0.0087, 0.181, 0.0069 | 0.2 | PASS | T |
| 12 (a) | Y4 N = 5, K = -1: largest mode value K + eps lambda_max at eps = 0.3 (before D108e: +0.0854, hyperbolic) and at eps = 0.2 (D108e) | -0.2764 | < 0 | PASS | V, T |
| 12 (a) | Y4 N = 5, K = -1, eps = 0.2, 5 frames: G_k in range, max slope | 0.0072 | 0.2 | PASS | V |
| 12 (b) | orbit C slopes k = 2, 3, 4; Y4 N = 2 neighbour (reported) | -0.016, -2.0025, -3.998; 0.010, -1.86, -4.14 | [-0.1, 0.1], [-2.3, -1.7], [-4.3, -3.7] | PASS | T |
| 13 | orbit D at n = 3000: l_1, l_1 - l_2 | 0.00765, 0.00676 | [0.005, 0.012], [0.004, 0.02] | PASS | T |
| 13 | offsets, suite at full size (61 gated orbits, none degenerate, 209 pairs checked at valid handoffs, 28 order-values at handoffs with a failed singular-value check counted and not gated, D108f); before D108: 12 assertions failed (values 98 to 215 at n = 800) | 0 fail | 15, 5.0 | PASS | T |
| 13 | offsets at the n where the 256-bit singular values are valid (209 pairs; 6 of 61 orbits have an invalid n): max abs c_k; max abs c_k(n_2) - c_k(100) | 14.83; 4.008 | 15; 5.0 (D108f; was 1.0) | PASS; PASS | V |
| 14 | orbit C slopes (AC12b) | -0.016, -2.0025, -3.998 | bands | PASS | T |
| 14 | other starts in band, 60 starts, 3 excluded as lower-dimensional tori | 34 of 60 (56.7 percent; 59.6 percent of the 57 not excluded) | 50 percent (D108g; was 70 percent, D107i) | PASS | T |
| 15 (a) | median n_floor(Double64) / n_floor(Float64), 99 of 100 starts reach both | k = 3: 1.82; k = 4: 1.87 | [1.5, 2.5] | PASS | T |
| 15 (b) | BigFloat(256) module on its orbit, n = 150, 22 runs (before D108c: 4 of 22 failed by a NaN bound) | 0 fail | 1 | PASS | T |
| 15 (c), (d) | log form means Float64 vs Double64; representation (finite L where exp(L) underflows) | hold | 4 standard errors | PASS | T |
| 16 | health, overflow, equal columns; (e) starts reaching bitwise equal columns | all hold; 479 of 1000 | exact fields; reported | PASS | T |
| 17 | orbit invariance across chains, 1 vs N threads and chunks, 0 bytes after warm-up (Y3, Y4 N = 100), reader first | bitwise; 0 bytes | bitwise; 0 | PASS | T |
| 18 | (a) shear bitwise in both modes; (b) wrapped vs lifted worst ratio | bitwise; 0.0096 | 1 | PASS | T |
| 19 (a) | ChaosTools GALI_2, module / oracle against the exact value at its own orbit, 98 resolved start-handoffs | 0.091 / 0.087 | 1 | PASS | V |
| 19 | oracle orbit equals the stream orbit bitwise (n <= 1000) | 2 of 42 starts (orbit divergence, max coordinate difference 4.4e4 for chaotic starts) | bitwise, else own-orbit references | INFO | V |
| 19 (a) | pairs with abs(G_module - G_oracle) <= B_module + B_oracle | 86 of 98, max absolute difference 1.3e-9 (non-bitwise orbits) | reported | INFO | V |
| 19 (b) | confirmation: ChaosTools GALI_3, GALI_4 on the K = 0 shear (n = 10, 100, 1000), error against closed form; orbit C slopes | at most 4.2e-17; -1.994, -3.948 | 4 B_meas; [-2.3, -1.7], [-4.3, -3.7] | PASS (not withdrawn) | V |
| 19 (b) | module / oracle against the exact value at own orbit, GALI_3 / GALI_4 (89 / 85 resolved) | 0.013, 0.021 / 0.0048, 0.0056 | 1 | PASS | V |
| 19 (b) | pairs within B_module + B_oracle, GALI_3 / GALI_4 | 89 of 89 (max pair ratio 0.0079) / 85 of 85 (0.0015) | 1 | PASS | V |
| 19 (b) | excluded as unresolved (reference below 10 B_meas): regular starts (D108h gate), GALI_2 / GALI_3 / GALI_4; chaotic starts, reported | regular 0 / 0 / 0 percent; chaotic 44.4 / 58.7 / 65.1 percent (pooled before: 29.4 / 32.5 percent for k = 3, 4) | regular at most 10 percent; chaotic reported | PASS (regular); INFO (chaotic, a property of the Float64 floor at n = 1000) | V |
| 20 | Y4 N = 100 open and ring, K = 4, 200 orbits, n = 1000: range, nested, finite; N = 10 against the 256-bit product | true; hold | as stated | PASS | T, V |
| 20 | G_2 of K = 4 vs K = 2 run, identical orbit | bitwise equal, 200 of 200 (both) | within 2 B_meas | PASS | V |
| 20 (V) | time per orbit-step, GALI K = 4 / K = 2 / Benettin k = 1, open; ratio K = 4 to k = 1 | 2697 / 1781 / 1462 ns; 1.85 | reported | INFO | V |
| 20 (V) | same, ring; first-call compile of the K = 4 chain 42.7 s | 2323 / 1764 / 1469 ns; 1.58 | reported | INFO | V |
| 21 (a) | regular class, 1024 orbits, lifted, Float64 on the device: Y3 and Y4 N = 2 at n = 4e4, Y4 N = 5 at n = 1e4; each backend against its own 256-bit reference (first 8, 8, 4 orbits): G, L worst ratio | at most 0.0011 (G), 2.0e-5 (L, B_L of D108a) | 1 | PASS | G |
| 21 (b) | regular class, orbits within 1e-9 (all 1024; max 1.9e-11): worst G, L gap in units of the bound | at most 0.0017, 6.4e-4 | 1 | PASS | G |
| 21 (a) | chaotic class n = 8, 1024 orbits: Y3 G / L worst ratio, CPU and GPU; Y4 N = 5 G | 0.51 / 0.14 (CPU), 0.64 / 0.11 (GPU); 0.049 | 1 | PASS | G |
| 21 (a) | chaotic Y4 N = 2: G and L ratios; floor-control violations (D108b: 2 B_meas + B_G^k below B_meas; G gated for G_ref >= B_meas) | CPU 0.19, 0.0016; GPU 0.13, 0.0018; 0 violations (before D108: 0.19, 0.13; L 0.056, 0.058; 53) | 1; 0 | PASS | G |
| 21 (b) | chaotic class, orbits within 1e-9 (all 1024; max 3.2e-11), gap in units of the plain sum of the single-backend bounds (D108i: reported, the orbit-divergence term of B11 not formed): Y3 G / L; Y4 N = 2 G / L; Y4 N = 5 G (before D108: 0.92 / 1.51; 24.2 / 16.4; 21.1, gated as FAIL) | 0.92 / 0.048; 24.2 / 1.39; 21.1 | reported | INFO (21(a) against each own reference gates the class: PASS) | G |
| 21 | Double64 compile gate, Y3, Y4 N = 2, 5 | FALLBACK, "kernel does not compile for the element type" | compile or reason | FALLBACK (valid result) | G |
| 21 | Float64 reducer kernel (MMatrix stack buffers, Y3, Y4 N = 2 and N = 5, K = 4) | compiles and runs on the device, no fallback | compiles | PASS | G |
| R1 | 1000 chaotic Y3 starts (2, 0.5, 0.1), n <= 300: crossing steps of G_2 (threshold tau_g = B4 map) vs SALI at tau_S = 1e-3, 1e-6, 1e-10: starts with different steps; starts with no crossing | 0, 0, 0; 13, 34, 72 | recorded | INFO | V |
| R2 | crossing of G_3, G_4 at 1e-8 and 1e-16 (G_4 by the log form, G_3 at 1e-16 in Double64): orbit D 470, 120 / 470; orbit C: G_3 none, G_4 8250 at 1e-8 only; 100 chaotic starts: medians G_3 24, 45 / G_4 13, 24; 100 regular starts: G_3 none, G_4 41 crossed at 1e-8 (median 12 420), none at 1e-16 | see log | recorded | INFO | V |
| R3 | realised ratios: criteria 1, 2, 4, 6 (above); 8 to 10 (above); 21 (above) | see rows | recorded | INFO | T, V, G |
| R4 | orbit D, Float64 direct G_4 vs exact at n = 400 / 800 / 1500 / 3000: error | 6e-25 / 6.5e-29 / 5.9e-30 / 9.3e-36 (G_ref 8.6e-16 / 2.1e-23 / 9.6e-29 / 1.1e-38); k = 2: 3e-16 to 8e-16 | recorded | INFO | V |
| R5 | starts classified as lower-dimensional tori | 3 of 60 | reported | INFO | T |
| D107a / D108a | B_L = 32 eps (1 + a_max) sum Psi: Float64 chaotic Y3 (before: 32 eps sum Psi, worst error 1923 eps sum Psi) | worst error 0.15 eps sum Psi (1 + a_max), 0.0047 of the bound (Float64); 0.0026 eps sum Psi (1 + a_max) (Double64) | 1 | PASS | V |
| D107b | L_ref = ln abs(det W_0) + sum ln abs(det M_i) - sum S_i | used in 9, 15, 21; matches the 512-bit product (n = 800 check) | - | PASS | T, V |
| D107c | identity gated by the absolute form (criterion 5) | 0.097; 0.0010 | 1 | PASS | T |
| D107d | tolerance of criterion 4 with the SALI floor | as 4 | 1 | PASS | T |
| D107e | frame check window 2 B_meas; check_frame not in propagate | as 7 | - | PASS | T |
| D107f | Gram control fails in half of the sets at m = 9 to 13 | 150 of 150 | half | PASS | T |
| D107g / D108d | criterion 11 window: (100, 200) for d <= 10, (400, 800) for d > 10 | see 11 | 0.5 | PASS | T, V |
| D107h / D108f | criterion 13 resolved orders only, k = d by the log form, singular-value check, stability bound 5.0 | see 13 | - | PASS | T, V |
| D107i / D108g | criterion 14 full-size pass rate 50 percent | 56.7 percent | 50 percent | PASS | T |
| D107j | criterion 15(a) on the chaotic panel | 1.82, 1.87 | [1.5, 2.5] | PASS | T |
| D107k | criterion 18(b) where wrapped and lifted orbits stay together | 0.0096 | 1 | PASS | T |
| D107l | scale factors, 4D overflow map, trace order, Gaussian frames at d = 200 | used; AC3(e), 16(c), 17, 20 pass | - | PASS | T |
| D108b | floor control 2 B_meas + B_G^k | see 9 (d), 21 (a) | 0 violations | PASS | V, G |
| D108c | reference depth 1e-70, Lambda NaN-safe | see 9 (a), 9 (c), 15 (b) | - | PASS | T, V |
| D108h | 19(b) allowance on regular starts; own-orbit references | regular 0 percent; both sides against their own reference: see 19 | 10 percent | PASS | V |
| D108i | 21(b) chaotic class reported against the plain bound sum | see 21 (b) | - | INFO | G |
| D108j | subnormal column norm is a health failure (review finding 1) | `ga_d108.jl`: RED 6 pass, 14 fail; GREEN 20 of 20 | - | PASS | T |
| D108k | test-side changes | package suite 2724 of 2724; full-size groups above | - | PASS | T, V |

## Figure

`gali-laws.png`. Left: orbit D, ln G_2, ln G_3 (direct form, plotted above 10 B_meas, which ends G_3 at n about 1400) and ln G_4 (log form, exact to any depth) against n, with the predicted slopes -E_k from the exponents of the 256-bit product at n = 3000 (E_2 = 0.0068, E_3 = 0.0153, E_4 = 0.0306), anchored at n = 500. The curves follow the predicted decay up to n about 1000 and then flatten (all three together, with a joint step near n = 2000): orbit D is sticky and its finite-time exponents are not constant, so a single straight line is only a guide; criterion 13 uses the n-dependent exponents l_i(n). Right: orbit C, ln G_k against ln n, slopes 0, -2, -4 (dashed, anchored at n = 1000), 10 B_meas(Float64, 4) dotted.

## Findings and notes for the owner

The findings 1 to 7 below are the state at the checkpoint of 2026-11-04 (before D108); each is resolved by the amendment item named in the rework section. The numbers in them are the old ones.

1. Same-orbit reference saturates at 256 bits. For a strongly chaotic start G_4 falls below 2^-256 (about 1e-77) after 100 to 200 steps, the 256-bit unit columns then give a volume near 2^-256 or exactly 0, `ga_Lambda` forms -Inf - (-Inf) and the bound B_G^k is NaN. The package suite hid this with n = 150 for chaotic starts (D38); at the spec horizon n = 1000 it fails 40 assertions of AC9 and 4 of AC15b by NaN. Gating the direct form where G_ref >= 1e-70 (`DIAG-AC9`, same cases and tolerances) gives worst ratios 0.72 (Float64) and 0.072 (Double64): the reducer agrees with its reference. Remedy for the spec: a reference of 1024 bits for the volume of chaotic starts, or the unit-column volume replaced by the log-singular-value route; or the horizon of criterion 9 stated per class.
2. The log-form bound B_L (D107a, multiplier 32) is exceeded by a factor up to 60 on lifted chaotic Y3 starts in Float64 at n up to 1000, and is not exceeded elsewhere. Cause isolated: the error scales with eps sum Psi (1 + a_max) with a_max the largest lifted coordinate (reached 9e3 to 3e4 at n = 500 to 700): the tangent is evaluated at the rounded sums q + p of a lifted state (the coupled-maps criterion 1 form 16 eps (1 + a_max)). Realised error over eps sum Psi (1 + a_max) is at most 0.15 (Float64), 0.0026 (Double64). The orbit used by the reference is bitwise the module orbit (checked on the chaotic Y3 starts of criterion 9), so this is not orbit divergence. Proposal: B_L = 32 eps (1 + a_max) sum Psi, or wrapped state for long chaotic runs. Criterion 9(d) (3 order-handoffs on one Y4 N = 5 start, errors 4e-15 to 1e-14 at G_ref about 1e-16, each inside B_G^k) fails for the same reason: the chaotic Lambda term in B_G^k exceeds the 2 B_meas floor of 9(d); the control should use the full bound.
3. Criterion 12(a), Y4 N = 5, K = -1, eps = 0.3: the open chain has lambda_max = 3.618, so the mode value is K + eps lambda_max = +0.0854 > 0: a hyperbolic mode, not an elliptic origin, so G_k decays (the suite's 10 assertions fail as 5 NaN slopes and 5 range failures). N = 2, 3 (-0.4, -0.1) are elliptic and pass. With eps = 0.2 (all modes elliptic) N = 5 passes with slopes at most 0.0072. The spec text "elliptic modes" holds for N <= 3 only.
4. Criterion 13. The 256-bit singular values of the long product (suite helper `ga_logsv`) are invalid at n = 800 for 6 of 61 orbits (the sum of the log singular values, which a symplectic product requires to be 0, is up to 215). On one such start the module log form equals the 512-bit reference to the printed digits (L = -2187.9346747370 at n = 800), so the failing offsets (|c| 98 to 215) come from the reference exponents, not from the module. At the valid n the offset bound 15 holds (14.83, tight) and the stability bound 1.0 fails once (4.0, one Y4 N = 3 start with l_1 - l_2 = 0.104 at n = 400, a near-degenerate pair, the slow convergence of the Eq. (44) case). The constant 1.0 is a proposal (P6).
5. Criterion 19(b) is not withdrawn. ChaosTools `gali(tands, T; threshold = 0)` takes the order from the number of deviation vectors and returns GALI_K only (not the lower orders), so GALI_k is checked from a k-vector run; both the closed forms of the K = 0 shear and the regular slopes (-1.994, -3.948) hold, and at all resolved handoffs the pair is within B_module + B_oracle for k = 3, 4. The Y3 hand system does not reproduce the stream orbit bitwise (2 of 42 starts; chaotic starts diverge by e^{l_1 n}); the check therefore uses each side's own orbit against the 256-bit reference at that orbit (equivalent to the orbit-divergence term). The 10 percent exclusion allowance of 19(b) is not met pooled (29 and 33 percent): at n = 1000 the chaotic starts are beyond the Float64 floor (59 and 65 percent of the chaotic start-handoffs), regular starts 0 percent; the exclusion is a property of the horizon.
6. Criterion 21. The reducer kernel (MMatrix buffers, 6.4 kB at d = 200 in the Float64 stream, here d = 4 and 10) compiles for the device and runs for Y3 and Y4 N = 2, 5; Double64 falls back with the recorded reason (accepted). The regular class passes everywhere. In the chaotic class (n = 8) each backend passes against its own reference for Y3 and Y4 N = 5 (floor control for Y4 N = 2 fails by the 9(d) effect, 53 order-handoffs); the CPU-GPU gap on orbits that agree to 3e-11 exceeds the plain sum of the two single-backend bounds by 1.5 (L, Y3) to 24 (G, Y4 N = 2). The orbit-divergence term of B11, G_ref (exp(D_orb) - 1) with D_orb the first-order effect of the 1e-11 orbit gap on the tangents, was not formed in this script: this part of 21(b) is not assessed, and the verdict FAIL is that of the plain sum. Open work: the D_orb machinery of sali criterion 18 transposed to d >= 4.
7. Criterion 14 at full size: 34 of 60 starts in band (56.7 percent). The P6 fallback of 60 percent also fails (59.6 percent of the 57 not excluded); the k = 2 band of 0.6 was not tested.
8. Criterion 20 timings were taken while other groups and a stuck Julia process (first AC20 attempt: it used the default Haar frame at d = 200, whose 200 x 200 QR does not finish, D107l; the script now passes Gaussian frames) were running; the ratios 1.85 and 1.58 are indications only. First call of the K = 4 chain: 42.7 s (stage 03: 40.7 s).
9. Not run: the Double64 rows of 21(a, b) (device refuses Double64), a bitwise 1 against N threads comparison beyond the suite's (AC17b uses the run's threads), criterion 20 N = 10 at n = 1e3 (suite: n <= 50 as specified).

Documentation and code inconsistencies found (for the final reply, CLAUDE.md rule): (a) spec B10 and criteria say B_L = 8 eps sum Psi, D107a says 32, and neither contains (1 + a_max); (b) spec criterion 12 calls Y4 N = 5, K = -1, eps = 0.3 elliptic: it is not; (c) spec criterion 13 says every order is resolved in BigFloat (false at depth) and the test-side singular values fail at n = 800; (d) spec Algorithm item 1 says `propagate` rejects a bad frame, the code leaves `check_frame` to the Experiments layer (D107e); (e) spec criterion 19 and F7 read as if ChaosTools returns all orders up to k, it returns only the order of the vectors; (f) spec criterion 11 gate (n_1 >= 40) fails for N = 10 and holds at (100, 200) only for N <= 5; (g) the test helper `ga_Lambda` is not NaN-safe for an exactly zero reference volume.

## Rework after the review and Amendment D108 (2026-11-05)

Owner decision at the checkpoint: "Amend spec, fix MEDIUM". Spec Amendment D108 (a to k) corrects the bounds, data and references of criteria 9, 10, 11, 12, 13, 14, 15(b), 19(b), 21 and changes the health rule (D108j); the reducer never disagreed with its 256-bit reference. Order of work: tests first, then the code, then the test-side amendments, then the full-size reruns.

1. RED (D108j). New test file `test/gali/ga_d108.jl` (20 assertions, Float64 and Double64; BigFloat has no reachable subnormal range): the subnormal column `[5e-324, 5e-324, 0]` next to `[1, 2, 3]` gives the failure code with column 2 and threshold `floatmin(T)` in `reducer_step` / `reducer_error` / `reducer_update`, no finite value above 1 from `gali_values`, a rejection from `check_frame`; a column just above `floatmin(T)` is accepted. Against the stage 03 code: 6 pass, 14 fail, 0 error (`gali_values` gave 1.1649 for the two columns, true 0.8238).
2. GREEN. `src/Indicators/gali.jl`: `_first_bad` rejects a norm below `floatmin(T)`; `reducer_error` reports the threshold `floatmin(T)` (0 for a zero norm, Inf for a non-finite one); `gali_values` turns a subnormal column norm into NaN (so `check_frame` rejects the frame and the orders containing that column are non-finite; lower orders stay valid). Constant messages in the kernel path, 0 allocations (unchanged), no function over 50 lines; sali not changed. `ga_d108.jl`: 20 of 20.
3. Test-side (D108k, `test/gali/`, test plan section "Rework D108 (2026-11-05)"): `ga_BL` with the a_max factor, `ga_Lambda` NaN-safe and the `G_ref >= 1e-70` gate, floor control with the full bound, criterion 11 window by d, criterion 12 Y4 N = 5 at eps = 0.2, criterion 13 singular-value check and bound 5.0, criterion 14 pass rate 50 percent. Gali package suite: 2724 of 2724 (2704 plus the 20 of `ga_d108.jl`), 0 fail, 2 min 14 s.
4. Whole package suite (`Pkg.test()`, 4 threads, Aqua, JET and layering inside): 1 408 277 pass, 0 fail, 1 existing broken (stage 03: 1 408 257; +20 = `ga_d108.jl`), 32 min, `quality` (Aqua, JET) 12 of 12, `layering` 2 of 2, `gali` 2724 of 2724. A first whole-suite run (31 min) failed one Aqua check, unbound type parameter of `_first_bad(a::NTuple{K, T}) where {K, T}` (K = 0); the signature became `NTuple{K}` with `eltype(a)` (same behaviour, `ga_d108.jl` 20 of 20 again) and the whole suite was run again (log scratch `ga04/full-d108.out`; the first run `full-d108-run1.out`). The full-size groups of step 5 ran on the code before this signature change (identical behaviour).
5. Full-size reruns (`verify/gali.jl`, `verify/gpu.jl` updated to D108a to i; a headed section "D108 rerun (2026-11-05)" is appended to `gali-verify-log.txt` and `gpu-gali-verify-log.txt`, the earlier sections are kept): AC9-12 728 pass, 0 fail (before: 630 / 98); AC13-16 724 pass, 0 fail (before: 723 / 17); D108 20 pass; DIAG-AC9, DIAG-AC11, DIAG-AC12, DIAG-AC13, AC19 as in the table; `GPU_ONLY=gali`: PASS 25, FALLBACK 3 (Double64 gate), INFO 23, FAIL 0 (before: FAIL 2, PASS 23). The check of 21(a) in `verify/gpu.jl` gates the direct bound for all orders with G_ref >= B_meas (was >= 10 B_meas) because D108b states the floor control for G_ref < B_meas only; a first rerun with the old split showed 26 floor-control violations per backend (Y4 N = 2) for orders between B_meas and 10 B_meas, all inside B_G^k.

Verdicts changed by the rework (old to new, key number):

| Criterion | Old | New | Key number |
|---|---|---|---|
| 9 (a) suite | FAIL as run (NaN reference) | PASS | 0 fail; 3546 resolved Float64 order-handoffs, worst 0.72 |
| 9 (b) | FAIL (60.1) | PASS | worst 0.0047 (Float64 chaotic), 0.00008 (Double64) |
| 9 (d) | FAIL (3 handoffs) | PASS | 0 violations of 2 B_meas + B_G^k |
| 10 (a) | PASS | PASS | worst L ratio 0.0002 with the a_max factor |
| 11 (N = 10) | FAIL (28 outside) | PASS | window (400, 800): 20 / 0 / 30 |
| 12 (a) Y4 N = 5 | FAIL (hyperbolic data) | PASS | eps = 0.2, mode value -0.2764, max slope 0.0069 |
| 13 | FAIL (12 assertions; stability 4.008 > 1.0) | PASS | 209 pairs, max abs c_k 14.83, max stability 4.008 <= 5.0 |
| 14 | FAIL (56.7 percent < 70) | PASS | 34 of 60 (56.7 percent) >= 50 |
| 15 (b) | FAIL (4 of 22) | PASS | 0 of 22 fail |
| 19 (b) exclusion | INFO above the allowance | PASS (regular 0 percent), INFO (chaotic 58.7, 65.1 percent) | |
| 21 (a) Y4 N = 2 chaotic | FAIL (53 floor-control violations) | PASS | 0 violations, G 0.19, L 0.0016 |
| 21 (b) chaotic | FAIL (plain sum, up to 24.2) | INFO (reported, D108i) | G 0.92 / 24.2 / 21.1, L 0.048 / 1.39 |

Not run again at full size: AC1-4, AC5-8, AC17-20, R1, R2, R3-AC2, R4, AC20t, figure (unchanged by D108; the AC17-20 health testsets pass in the package suite). Timings of the rerun were taken on a loaded machine (Julia process 62696 of another session and the whole-suite run were active). Open: the orbit-divergence term D_orb of B11 for d >= 4 (21(b) chaotic class, reported), as sali criterion 18 for d = 2; independent review findings 2 to 4 (LOW) unchanged.

## Audit

| Check | Result |
|---|---|
| Evidence | Criteria 1 to 21, R1 to R5, D107 a to l and D108 b, c, h, i, j, k each have a number in the table above; D108 a, d to g are shown in the rows of the criteria they change |
| Review closed | no CRITICAL or HIGH finding; finding 1 (MEDIUM) fixed in the rework (D108j, tests first); findings 2 to 4 (LOW) open, accepted; 5 closed by criterion 21; 6, 7 INFO. Open: 21(b) chaotic class, D_orb term not formed (reported, D108i) |

## Independent review

Fresh Sonnet subagent, given the spec (with D107) and the code (`src/Indicators/gali.jl`, the Tangent protocol, the sali reducer), not the notes or the test plan; read only, short snippets in scratch.

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | `gali.jl:129-141` (`_first_bad`, line 122): only zero and non-finite column norms are rejected; a subnormal column passes and its unit column is rounded to a few denormal ulps. Confirmed: `gali_values([1,2,3],[5e-324,5e-324,0])` gives 1.1649 (true 0.8238, above the bound 1); the sum S_i receives a wrong logarithm. Silent wrong value, pathological input | MEDIUM | fixed, tests first (D108j, rework 2026-11-05): a norm below `floatmin(T)` is a health failure; `ga_d108.jl` RED 6 pass, 14 fail, GREEN 20 of 20 |
| 2 | `gali.jl:105-109`: `gali_values` has no check D >= K; a 2 by 3 input gives G_3 = 6.7e-33 instead of an error. `reducer_init` checks | LOW | open |
| 3 | `gali.jl:109, 183-189`: `gali_values(::AbstractMatrix)` and `check_frame` build a dynamic-size `SMatrix` (type-unstable); first call 0.56 s (Float64), 0.88 s (Double64), 1.85 s (BigFloat) at d = 200, K = 4; later calls 30 us to 0.9 ms | LOW | open; acceptable for one call per frame |
| 4 | `gali.jl:60-69, 183`: integer element types half-supported (exact 0 and 1 entries return an Int; others `InexactError`; `check_frame` fails in `effective_eps`) | LOW | open |
| 5 | `gali.jl:60`: `MMatrix{D,K,T}` buffer of 6.4 kB at D = 200, K = 4 risks the GPU stack | LOW | closed by criterion 21: the Float64 kernel compiles and runs for Y3 and Y4 N = 2, 5 (N = 100 not run on the GPU) |
| 6 | `gali.jl:87, 91`: R_11 is left out of the product; the column is unit within eps, inside 2 k eps | INFO | no change |
| 7 | `gali.jl:168`: L omits ln abs(det W_0); documented, as sali | INFO | no change |

Confirmed correct by the reviewer: the nested modified Gram-Schmidt gives the volume of the leading k unit columns (near-dependent sets in d = 4, K = 4, angles 1e-1 to 1e-20, 1000 sets against 2048 bits: worst 0.14, 0.087, 0.007 of the floor for k = 2, 3, 4); K = d = 2 bitwise equal to sali on W, sums, compensation and L (Float64, Double64, 200 runs of 300 steps); log form sign and compensation; NaN, Inf and zero columns give code 1 with W and sums unchanged (Float64, Double64, BigFloat); columns scaled by 1e+-200 give identical values; `check_frame` rejects a NaN frame and an angle of 1e-15; constant `ArgumentError` text in `reducer_init`; no `@fastmath`; longest function 21 lines, file 190 lines. Not checked by the reviewer: the GPU, the BigFloat reducer over a long run.
