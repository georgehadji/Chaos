# Verification: rli

Spec: `../../01-spec/output/rli-spec.md` (criteria 1 to 15, Amendments D76 and D77). Test plan: `../../02-tests/output/rli-test-plan.md`. Impl notes: `../../03-implement/output/rli-impl-notes.md`. Package tests: `test/rli/` (`rl_*.jl`). Verify script: `verify/rli.jl`, run as `julia -t 2 --project=verify verify/rli.jl` from the repo root; full output in `rli-verify-log.txt` (one line per measurement: `AC<n> | quantity | value | bound | verdict`). GPU: `verify/gpu.jl`, tag `rli` (`GPU_ONLY=rli julia --project=verify -t 4 verify/gpu.jl`), output in `gpu-rli-verify-log.txt`. Steps 1 and 3 (evidence part) of the stage contract; step 2 (independent review) is done by another agent, see "Independent review". Rerun after the fixes of the stage 04 checkpoint (2026-10-22, decision D77): the logs and this report are the fresh ones.

Run facts. Julia 1.13.2, 2 threads (GPU run 4 threads), base seed 20261012, verify offset 59000 (seeds `RL_VSEED0 + 50 ac + sub`; the package tests use `RL_SEED + 36000 .. 36999`, so the verify samples are fresh; the K = 0.5 regular class and the AC14 panels use the package screening and seeds of `rl_case`, as the spec's panel). Working tree of branch `main`, nothing committed. Wall time of the final run 473 s, of which the rli suite at full size 349 s (65 590 assertions passed, 0 failed, 0 errors; the whole package suite at normal size has 1 393 101 pass, 0 fail, 0 error in 12 min 22 s, rli 64 480 including AC11f; after the AC4 change the rli suite alone is 64 460 pass). The measurement helpers are those of `test/rli/rl_helpers.jl` (`rl_ref_rli`, `rl_PairRef`, `rl_BR`, `rl_BDL`, `rl_Psimax`, `rl_pair`); the loops are written in the script, so a definition error in a helper would not be caught by the script; the independent pieces are the ChaosTools oracle line (criterion 14), the closed forms of RL9 and the bitwise identities in the suite. The per-testset counts of the full run are in the log (`SUITE | testset` lines), as are the realized ratios of its `@info` lines.

## Evidence table

"Where": T = package test (full size through the `suite` section), V = `verify/rli.jl`, G = `verify/gpu.jl`. Two entries per row are Float64 / Double64 unless stated. Bound 1 means the ratio error / bound. Rows marked T without a number passed all assertions of the testset at full size; the suite prints realized ratios only for the `@info` lines.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| suite | rli suite at full size (`RL_FULL_RUN`), 34 testsets incl. AC11f | 65 590 pass, 0 fail, 0 error | all pass | PASS | T via V |
| suite | whole package suite, normal size | 1 393 101 pass, 0 fail, 0 error (12 min 22 s; rli 64 480 incl. AC11f; 64 460 after the AC4 change) | all pass | PASS | T |
| 1 | algebra of the definition on the test map T1, tau in {1,3,5}, closed form, mixed signs, doubling, negative controls: realized ratio error / B_R^lin | 0.0166 / 0.0023 | 1 | PASS | T |
| 2 | constant tangent (Y5, K = 0 shear, n to 1e5): every DL_k and RLI exactly 0, sums bitwise equal, negative control fails | all exact (tolerance 0) | 0 | PASS | T |
| 3 | same-orbit reference, K = 0.5 / 2 pi / 10, Dx in {1e-6, 1e-12}, tau in {1,4}, 20 starts to n = 1e4: RLI error over B_R, worst over all cells (Float64); Double64 worst | 0.0185; 0.0008 | 1 | PASS | T |
| 3 | same, DL_k error over its bound, worst (Float64; Double64 0.0026) | 0.0548 | 1 | PASS | T |
| 4 | (a) BigFloat(256) own orbits, 20 starts, n in {100, 1e3, 1e4}, Dx 1e-14 to 1e-6: slope of log RLI against log Dx_r (60 cases) | min 0.9998, max 1.0000; 0 failing | [0.99, 1.01] | PASS | T, V |
| 4 | (a) O(Dx) law of rho (D76 item 3), 60 cases; reported, not gated (D77 item 1), evaluated where Dx_r n <= 1e-3 | below the law in 1 of 60 (start 17, n = 1e4: dev 4.38e-12, 4.39e-10, 8.8e-9, 5.8e-5 at Dx 1e-12, 1e-10, 1e-8, 1e-6; allowance at 1e-10 was 1.77e-10, a first-to-second order cancellation) | reported | INFO (D77) | T, V |
| 4 | (a) realized second-order coefficient max dev_i / Dx_i (the spec's 20 was replaced by D76) | 5884 (start 19, resonance) | reported | INFO | T |
| 4 | (b) Double64 against BigFloat rho, Dx 1e-20 to 1e-8: worst relative gap; slope in range | 9.5e-10; yes | 1e-3; [0.99, 1.01] | PASS | T |
| 4 | (c) direction u = xi_2, ratio rho(xi_0) / rho(xi_2) at Dx = 1e-12 | 0.056 to 13.4 | reported | INFO | T |
| 5 | chaotic independence of Dx, K in {2 pi, 5, 10}, 1000 starts at n = 1e3, 1000 at n = 1e4 (full size), pairs (1e-12, 1e-10), (1e-12, 1e-8): paired difference of means against allowance, largest ratio of difference to allowance | 0.63 (K = 5, n = 1e3, Dx2 = 1e-8: 0.0045 of 0.0071); all others below 0.62 | 1 | PASS | T |
| 5 | R3: ensemble mean RLI at Dx = 1e-12 / 1e-10 / 1e-8, n = 1e4: K = 2 pi; 5; 10; ratio to the regular class at 1e-12 | 0.0213 / 0.0214 / 0.0224; 0.0327 / 0.0330 / 0.0330; 0.0207 / 0.0210 / 0.0208; 4.0e10 to 6.4e10 | reported | INFO | T |
| 6 | tau dependence, tau in {1,2,5,10}, K in {0.5, 2 pi}, n in {100, 1e3}, Dx in {1e-8, 1e-12}: sums (a) and RL10 bound (b) | all pass | 1 | PASS | T |
| 6 | reported: ratio RLI_tau / RLI_1 for the regular class (Float64 and Double64) | 0.728 to 2.366 | reported | INFO | T |
| 7 | (a) Dx_r from stored states, minimal image, 1000 starts, both modes, Dx in {1e-14, 1e-12, 1e-8}; (b) rejection (kind separation, pair id, before any work); D76 items 1, 2 | all pass | 0 | PASS | T |
| 7 | (c) Float64 against BigFloat own orbits, 20 starts, Dx in {1e-8..1e-12}: at least 18 of 20 within 10 percent (all Dx 1e-8 to 1e-12 pass); worst |RLI_F64 / RLI_BF - 1| at Dx = 1e-12, n = 1e3 | 0.083 | 0.1 | PASS | T |
| 7 | (c) realized ratio |RLI_F64 / RLI_BF - 1| / Theta_n, worst over all cells; Dx = 1e-14 (reported, not gated): worst relative error at n = 1e3 | 0.0074; 10.3 | reported | INFO | T |
| 7 | (d) R8 lifted against wrapped, rotational start, n = 1e4: fraction of handoffs with ulp above Dx_r / 8 (lifted); rho lifted / wrapped | 0.924; 1.196 / 0.428 | reported | INFO | T |
| 8 | cross-check with the single-orbit stream and L1: (a) orbit points bitwise, chains {}, {RLI}, {L1}; (b) sums bitwise L1; (d) chain {RLI(f), MEGNO} bitwise | all pass | 0 | PASS | T |
| 8 | (c) BigFloat post-processing of the stored L1 series: DL error over 16 eps tolerance; mean | 0.123 / 0.125 / 0.087 (K = 0.5 / 2 pi / 10); 0.0052 / 0.0050 / 0.0044 | 1 | PASS | T |
| 9 | symmetry (orbit exchange, sign of xi_0) and column independence, both types | all bitwise | 0 | PASS | T |
| 10 | (a) Float64 against Double64, K = 2 pi, 5 (n = 1e3, 1000 starts), K = 5, 10 (n = 1e4, 200 starts): within 4 combined SE; means shown | pass; e.g. 0.05101 / 0.05104 (K = 2 pi), 0.03125 / 0.03152 (K = 5, 1e4) | 4 SE | PASS | T |
| 10 | (b) BigFloat(256) own orbits against 512 bits, K = 2 pi, 50 starts: error over B_R | 0.0055 | 1 | PASS | T |
| 10 | (c) BigFloat against Double64, n = 300, 100 starts; (d) K = 10, n = 1e4 representation (finite, 0 to 2 L_sigma, not clipped) | pass; RLI 0.0073 to 0.032 | 4 SE; 1 | PASS | T |
| 11 | health and frames: (a) zero tangent in block A, B, both; (b) NaN; (c) ten bad frames rejected, seeded ids 1..1e4 accepted; (d) overflow tau 200 / 300; (e) run record | all pass | per testset | PASS | T |
| 11 | (f) pairing check (D77 item 3, review R1), testset "AC11f pairing check (D77)", 11 assertions: pair stream without block frame and plain standard map with identity frame fail at step 3 (kind "pairing"); block frame, Float64 and Double64, wrapped and lifted, both chain positions: no error, equal finals; odd pair length DimensionMismatch | all pass | per testset | PASS | T |
| 12 | (a) orbit points bitwise across chains {}, {RLI}, {RLI(f), MEGNO}, {MEGNO(f), RLI}; (b) `propagate` equals the serial per-pair body, chunk 64 and 1000; (c) 0 bytes after warm-up; (d) ids permuted | all pass | 0 | PASS | T |
| 12 | (b) 1, 2 and 4 threads in child processes, chunk 64 and 1000, chains {RLI}, {RLI(f), MEGNO}, tau in {1,4}: 4096 pairs, K = 2 pi, wrapped, n = 1e3: hashes of x, W and every final | identical (8 cases) | equal | PASS | V |
| 13 | wrapped against lifted: (a) K = 0 bitwise, RLI = 0; (b) K = 0.5, n in {10, 100}, 50 starts | all pass | 0; 1 | PASS | T |
| 14 | ChaosTools 3.5.4 two one-column tangent systems (base and shadow), lifted rule, Dx = 1e-12: oracle orbits against `forward`, bitwise differences: K = 0.5 (20 pairs x 1e4 steps x 2 orbits), K = 2 pi and 5 (1000 x 1e3 x 2) | 0, 0, 0 | 0 | PASS | V |
| 14 | RLI error over B_R + B_oracle (B_oracle = B_R, same orbits bitwise): K = 0.5 (80 comparisons), 2 pi (3000), 5 (3000); medians 9.1e-6, 1.7e-4, 1.7e-4; excluded 0 | 0.0016; 0.0104; 0.0056 | 1; <= 10 percent | PASS | V |
| 14 | DL_n error over B_DL + B_oracle, same panels | 0.0060; 0.0503; 0.0287 | 1 | PASS | V |
| 15 | (a) K = 0.5, Float64, 1024 pairs, n = 4e4, each backend against its own BigFloat(256) reference at both realized orbits (final states reproduced bitwise 1024 of 1024 on both), 65 handoffs: RLI error over B_R, CPU / GPU; DL | 0.104 / 0.098; 0.104 / 0.098 | 1 | PASS | G |
| 15 | (a) chaotic K = 10, 1024 pairs, n = 8, 10, 100: RLI error over B_R, CPU / GPU (all three n); DL, CPU / GPU | 0.154 / 0.154; 0.175 / 0.163 | 1 | PASS | G |
| 15 | (b) regular class, pairs with final states within 1e-9 (558 of 1024): gap over B_R,cpu + B_R,gpu (median 1.2e-5) | 0.0049 | 1 | PASS | G |
| 15 | (b) chaotic n = 8 / 10 (all 1024 enter; largest horizon 6.0e-7 / 2.6e-5): gap over B_R,cpu + B_R,gpu + (1/n) sum_k (D_orb^A + D_orb^B)/k | 0.027 / 0.034 | 1 | PASS | G |
| 15 | (b) looseness (INFO): chaotic gap without the D_orb term, max ratio 4.4e3 (654 of 1024 above 1) at n = 8, 1.4e5 (854 above 1) at n = 10; the term median 2.3e-10 / 4.8e-9, max 2.2e-8 / 4.1e-6 | | reported | INFO | G |
| 15 | (b) item 4 check on the device orbits, n = 8 / 10: max |z_k^cpu - z_k^gpu|_inf / (16 eps M_k Gamma_k), both orbits; against the literal 16 eps M_n Gamma_n | 0.226, 0.276 / 0.226, 0.538 | 1 | PASS | G |
| 15 | (b) statistical, reported: regular class RLI_cpu / RLI_gpu min / 1% / median / 99% / max over 1024 pairs; K = 10, n = 100 (orbits fully separated across backends) | 0.20 / 0.54 / 1 / 1.66 / 3.27; 0.14 / 0.21 / 1.00 / 4.9 / 7.2 | reported | INFO | G |
| 15 | (c) K = 10, 1024 pairs, n = 1e4, tau = 2: max per-pair gap / sigma_FT; difference of means / SE (means 0.02103 / 0.02102) | 5.77; 0.032 | 6; 4 | PASS | G |
| 15 | reported: pairs with final states apart by more than 1e-9, CPU against GPU: regular n = 4e4; chaotic n = 8; 10; 100 | 466; 116; 463; 1024 (of 1024) | reported | INFO | G |
| 15 | Double64 (regular class, chaotic class and the ensemble) | kernel does not compile for the element type (F-C) | recorded fallback reason is the pass | FALLBACK (pass) | G |
| 15 | Float64 pair kernel with the pairing check compiles and runs on the device (probe `compiles = true`); gpu.jl rli section | 22 PASS, 0 FAIL, 3 FALLBACK (Double64) | 0 FAIL | PASS | G |

Recorded lines of the spec. R1 (SandorMaffione2016 Eq. (2), nu = 0.5, ordered (2, 0), chaotic (3, 0), the identity F7 bitwise): orders of magnitude between chaotic and ordered RLI at n = 1e4 (the printed value is 9 to 10, separation and vector not stated there): 12.0 (Dx = 1e-14), 10.5 (1e-12), 8.9 (1e-10); at n = 100 the orbits have not yet separated (2.3 to 2.9). R2 (realized ratios error / bound): at most 0.055 (criterion 3 DL, Float64), 0.0185 (criterion 3 RLI), 0.0166 (criterion 1), 0.0055 (criterion 10(b)), 0.0074 (7(c) over Theta_n), 0.125 (8(c) against 16 eps), 0.0503 (criterion 14 DL against the doubled bound), 0.154 (criterion 15, either backend). R3: the AC5 line above. R4 (threshold demonstration, 500 uniform starts per K, rho = RLI / Dx_r above 1, i.e. RLI above 1e-12 scaled, Dx = 1e-12): fraction labelled chaotic at n = 1e3 / 1e4 / 1e5, K = 0.5: RLI 0.164 / 0.150 / 0.152, L1 (n X above ln n) 0.180 / 0.198 / 0.202, FLI_B 0.61 / 0.62 / 0.62, MEGNO (Ybar_1 above 3) 0.010 / 0.010 / 0.010; K = 2 pi: RLI 1 / 0.996 / 0.996 against 0.996 for L1, FLI, MEGNO (at most 2 disagreements in 500); K = 10: 1.0 for all four (0 disagreements). The thresholds are those of Maffione2011a and the labels are not pre-registered (claim C3 is for the study P2). R5 (20 regular starts x 20 xi_0 x 20 directions, Dx = 1e-9, all 8000 frames accepted by `check_frame`): rho min / median / max at n = 1e3: 0.0012 / 0.234 / 25.2, at n = 1e4: 2.0e-4 / 0.227 / 17.9; the per-start max / min ratio has median 104 and 226 (maximum 1507 and 840): the vector and the direction matter by two orders of magnitude on the same orbit, first order in Dx for each. R6 (resolution horizon, ordered start (2, 0), n to 2e4): RLI exceeds 10 times the worst-case RL11 floor at n = 10 for Dx = 1e-10 and never for Dx = 1e-12 and 1e-14; the measured Float64 error against BigFloat stays below 10 percent for all three (3.7e-5 at Dx = 1e-12 and 1.6e-2 at Dx = 1e-14 at n = 2e4), so the worst-case floor is pessimistic by orders (F5). R7 (cost, K = 2 pi, wrapped, 20 000 pairs x 1000 steps, tau = 1): CPU (2 threads) pair 35.1 ns per orbit-step, L1 single 16.2, MEGNO single 28.0, ratio pair / L1 2.2 and pair / MEGNO 1.26; GPU (RTX 5060, best of 3) pair 7.99 ns (CPU with 4 threads 23.3, speedup 2.9), L1 2.56 (CPU 10.0, 3.9), MEGNO 5.15 (CPU 20.3, 3.9), reducer state 88 bytes Float64 / 144 Double64. R8: the 7(d) row.

## Oracle line (criterion 14), as the spec states it

ChaosTools 3.5.4 has no RLI routine (spec Oracle). The line is the oracle's independent tangent evolution: two `TangentDynamicalSystem` objects with k = 1 and the hand-coded Jacobian of the lifted standard map, one for the base start and one for the shadow start (the stored `displaced` state), both with the same supplied xi_0, one `step!` per map step, the deviation read after each step, its log norm taken in BigFloat(256), then the deviation renormalised in Float64. DL_k and RLI are formed by RL4 in BigFloat (`rl_ref_rli`) from the two per-step log-norm series and compared with the module's pair run (`PairMap` of the lifted `StandardMap`, `RLI()`, supplied block frame, from the suite's `rl_run` path). First both oracle orbits equal the stream's lifted orbits bitwise at every step (0 differences in 400 000 + 4 000 000 steps of both orbits), so no orbit divergence term is needed and B_oracle is the RL11 form with the same tangents. The line tests the stretch series and the pairing against an independent code path, not the definition itself, which no independent implementation has for maps; the definition is tested by the closed forms of RL9 and the printed laws (criteria 1, 4, 5). It was run lifted, because the oracle's rule is the lifted one; the wrapped mode is covered by criterion 13 (wrapped against lifted) and the same-orbit references of criterion 3.

## Decisions taken (D77)

Owner decisions at the stage 04 checkpoint, 2026-10-22, recorded as spec Amendment D77 and decision D77. The findings F-A and F-B of the first run are resolved; F-C is unchanged.

F-A (resolved). The O(Dx) law on rho of criterion 4(a) is reported, not gated. The first proposal, assert the law only where Dx_r n <= 1e-3, still failed 1 of 60 cases at full size: start 17, n = 1e4, a cancellation between first and second order (dev 4.38e-12, 4.39e-10, 8.8e-9, 5.8e-5 at Dx 1e-12, 1e-10, 1e-8, 1e-6; the allowance at 1e-10 was 1.77e-10). The owner rejected a threshold chosen after the data. The slope test (60 of 60, 0.9998 to 1.0) and the BigFloat comparisons 4(b), 4(c) carry the criterion; the cases below the law are reported (1 of 60) with the realized second-order coefficient (5884, start 19). The unrestricted diagnosis in the log (`AC4diag`, 4 starts at n = 1e4) is informational. This was a criterion-form problem, not a code defect: the linear law RL6 holds.

F-B (resolved). Criterion 15(b) is gated in the D68 form (each backend against its own reference; the gap only where the final states agree within 1e-9: regular class 558 of 1024, ratio 0.0049; chaotic class with the D_orb term, 0.027 and 0.034), as for sali, FLI and MEGNO. The statistical form is reported (5.77 sigma per pair and 0.032 SE for the means at K = 10, n = 1e4; the ratio distribution). The regular-class gate is loose (the median B_R is 2.9e-10 against RLI of 1.6e-13), so it is a consistency check, not a sharp one; the relative CPU-GPU gap on agreeing pairs has median 2.5 percent and maximum 160 percent (reported): the regular-class RLI at Dx = 1e-12 is a difference of nearly equal numbers (F5), and the two backends round it differently.

F-C (unchanged). Double64 GPU kernel does not compile (probe `compiles = false`; `kernel does not compile for the element type`): allowed fallback, recorded in three rows. The Float64 pair kernel (PairTangent block action, 88-byte reducer state, 8 isbits state components in wrapped mode, now with the pairing check) compiles and runs on the device with `provenance.backend == :gpu` on every Float64 GPU line.

Review findings R1 to R5 are all fixed in code (D77 items 3 to 5), see "Independent review". Mismatched xi_0 between the blocks remains the caller's `check_frame` (documented).

Other points, no decision needed:
1. R6 and F5: the worst-case floor of RL11 is pessimistic (RLI is within 3.7e-5 of the BigFloat value at Dx = 1e-12, n = 2e4, where the first-order resolution number says unresolved); the Float64 resolution horizon for Dx = 1e-14 is not reached within n = 2e4 either (1.6e-2).
2. Criterion 5 realized differences are at most 0.63 of the allowance (K = 5, n = 1e3); the allowances could be tightened (K18).
3. R4: at K = 0.5 the RLI labels (15 percent) are between L1 (20 percent) and MEGNO (1 percent); FLI_B labels 62 percent at ln n. Not a pass criterion.
4. R1: the chaotic / ordered gap is 12.0, 10.5 and 8.9 orders at n = 1e4 for Dx = 1e-14, 1e-12, 1e-10 against the printed 9 to 10; the gap grows with the smaller Dx because the ordered value scales with Dx and the chaotic one does not.
5. Not measured here: coverage and size claims of the impl notes; the benchmark of `bench/rli.jl`; Experiments chain wiring (the one-replacer rule is the caller's, D72 item 1); `docs/` build.

## Margins and marginal criteria

Largest ratios in the rows with bound 1: 0.154 (criterion 15(a), both backends, chaotic K = 10), 0.125 (criterion 8(c) DL against 16 eps), 0.104 (criterion 15(a), regular class), 0.055 (criterion 3 DL), 0.050 (criterion 14 DL); all others below 0.035. Criterion 5 realized differences are at most 0.63 of their allowance. There are no failures; the O(Dx) law of criterion 4(a) is reported, not gated (D77).

## Deviations from the impl notes, checked

- Zero allocations in the per-pair body: 0 bytes at n = 1e4, tau in {1,4}, both types, both chains (AC12(c), full-size suite).
- Bitwise equality of block A and B with the single-orbit streams and of the sums with the L1 sums (I5): AC8(a), (b), (d) pass at 20 starts to n = 1e4.
- Reducer state 88 bytes Float64 / 144 Double64 (as the notes; the spec proposal was 64 / 112); the pairing check (D77 item 3) adds no state.
- `check_frame` threshold 8 eps: seeded frames for ids 1..1e4 accepted; the ten bad frames rejected; all 8000 frames of R5 accepted.
- Isbits `PairMap`, `PairTangent`, `RLI` and `RLIState`: the Float64 kernel compiles and runs on the device (the notes had not measured it, "Known limits"); Double64 does not compile.
- `Dx_r` is within 4 eps of `safe_norm(separation(...))` of Lyapunov, not bitwise (the notes' deviation): AC7a passes at 1000 starts, both modes.
- The notes' cost remark, "about 1.2 times a single MEGNO orbit": 1.26 here (CPU, 2 threads, K = 2 pi, wrapped).
- Coverage and size claims were not re-measured here.

## Runtime

`verify/rli.jl`: 473 s total (rli suite at full size 349 s, AC4 diagnosis 40 s, AC14 25.4 s, R4 22 s, figures 11.3 s, AC12 with the child processes 9.2 s, R7 6.9 s, R5 4.2 s, R6 1.9 s). GPU run (`GPU_ONLY=rli`, 4 threads): several minutes wall, not timed (65 prefix runs per backend and class at n = 4e4, BigFloat(256) references of 2048 orbits per backend, D_orb for 2 x 1024 chaotic orbits). The machine was not otherwise idle (the CPU and GPU runs overlapped for part of their length), so timings are not benchmarks; the kernel timings above are best of 3 on the RTX 5060 Laptop GPU (compute capability 12.0, CUDA.jl 6.4.2, driver 13.4).

## Figures

- `rli-growth-laws.png`: left, RLI(n) for the ordered start (2, 0) (solid) and the chaotic start (3, 0) (dashed) of SandorMaffione2016 Eq. (2), nu = 0.5, wrapped, Dx = 1e-12, 1e-10, 1e-8: the ordered curves are flat and proportional to Dx, the chaotic curves rise to a plateau near 1e-2 independent of Dx; middle, rho = RLI / Dx_r: the three ordered curves collapse onto one (RL6), the chaotic ones are separated by the factor Dx (RL7); right, RLI(1e3) against Dx_r for eight K = 0.5 regular starts (slope 1, least-squares slopes 0.9967 to 1.0000 over Dx 1e-12 to 1e-6) and eight K = 2 pi starts (flat).
- `rli-chaostools.png`: left, RLI of Chaoslab against RLI of the ChaosTools two-system oracle line (RL4 in BigFloat), K = 0.5, 2 pi, 5, on the diagonal; right, the ratio to B_R + B_oracle against n, all below 0.011 (bound 1).

## Files

- `verify/rli.jl` (new), `verify/gpu.jl` (new tag `rli`, section "AC15 (rli)" before the summary: functions `r15_trace_kernel!`, `r15_trace`, `r15_setup`, `rli_est`, `r15_ref`, `r15_Ms`, `r15_one`, `rli_chaotic_b`, `rli_d15`, `rli_ac12c`, `rli_timing`; it reuses `m15_grid`, `d16_orb` (with its `matsf` keyword), `gpu_why`, `gpu_reason`, `regular_starts`; `gpu_trace` was not changed).
- `rli-verify-log.txt`, `gpu-rli-verify-log.txt`, `rli-growth-laws.png`, `rli-chaostools.png`, this report.
- The first run edited nothing under `src/` or `test/`. The fixes of D77 (pairing check, DL split difference, u conversion, odd pair state, AC4 reported, AC11f) were made afterwards on the owner's decision; the logs and this report are from the rerun. Nothing was committed.

## Independent review

Fresh Sonnet subagent, given the spec (with Amendment D76; the fixes below are D77) and the code (`src/Indicators/rli.jl`, `src/Systems/pair.jl`, with the Tangent protocol and the neighbouring reducers for comparison); not the notes. Probes in the session scratchpad only. Result: no CRITICAL and no HIGH finding.

Confirmed by probes: bitwise orbit equality of the pair blocks with single orbits (Float64, Double64, wrapped, lifted, tau 1 and 3); log sums bitwise those of Benettin for block A; zero blocks stay exactly 0; n = 0 and the first handoff; tau above 1; type stability and 0 allocations after warm-up; failure recorded, not thrown; `pair_start` minimal image, widening and rejection; `check_frame` order. Not verified by the reviewer: the GPU (covered by criterion 15 above), the BigFloat paths, 11(d) overflow at tau = 300 and the RL9 closed forms (covered by the tests).

| # | Severity | Where | Finding | Disposition |
|---|----------|-------|---------|-------------|
| R1 | MEDIUM | `rli.jl:55-67`, `46-50` | The reducer cannot tell a pair stream from any even-d, k = 2 stream: RLI on `PairMap` without W0 (random seeded 4 x 2 frame) gave 0.00999 with no error; on a plain `StandardMap` with an identity frame 0.00805, and `check_frame` accepts it (0.0); mismatched xi_0 with `check_frame` skipped 0.00969. Suggested fix: a failure code in `reducer_step` when W[h+1:D, 1] or W[1:h, 2] is not exactly 0 at a handoff (the block action keeps them exactly 0; isbits, kernel safe) | Fixed (D77 item 3): runtime pairing check in the reducer, off-block entries of the arriving frame must be exactly 0 at every handoff, else health failure kind "pairing" after the norm check; testset "AC11f pairing check (D77)" (pair stream without block frame, plain standard map with identity frame: error at step 3; block frame, Float64 and Double64, wrapped and lifted, both chain positions: no error and equal finals; odd pair length DimensionMismatch); the Float64 GPU kernel compiles with it. Mismatched xi_0 stays the caller's `check_frame` (documented) |
| R2 | LOW | `rli.jl:53` | DL is `total(sB) - total(sA)`; the split form (sB.s - sA.s) + (sB.c - sA.c) is more accurate (7e-8 relative in Float64 at n = 300, Dx = 1e-10; 4e-30 in Double64) | Fixed (D77 item 4): DL from the split compensated difference; criterion 3 DL ratio now 0.055 |
| R3 | LOW | `rli.jl:169-178` | `u` is not converted to T (a Float64 u in a Double64 run is off xi_0 by about 1e-16); a wrong-length x0 for the wrapped map gives a MethodError in `torus_state` | Fixed (D77 item 5): u converted to the element type and documented |
| R4 | LOW | `pair.jl:45`, `:34` | An odd pair state length silently drops the last component | Fixed (D77 item 5): odd pair state rejected |
| R5 | LOW | `rli.jl:122-131` | The lifted separation is reduced mod 2 pi, so a lifted gap of exactly 2 pi reads as 0 and is rejected (the L2 precedent) | Fixed (D77 item 5): lifted 2 pi gap documented |

Audit. Evidence: every acceptance criterion 1 to 15 has a number or a plot above; no row fails (the criterion 4(a) law is reported under D77). Review closed: no open CRITICAL, HIGH or MEDIUM finding; R1 to R5 are fixed in code and covered by the rerun (rli full-size suite 65 590 pass, `verify/rli.jl` 73 PASS, 0 FAIL, gpu.jl rli section 22 PASS, 0 FAIL, 3 FALLBACK).
