# Verification: fli-ofli

Spec: `../../01-spec/output/fli-ofli-spec.md` (criteria 1 to 16, Amendment D72). Test plan: `../../02-tests/output/fli-ofli-test-plan.md`. Impl notes: `../../03-implement/output/fli-ofli-impl-notes.md`. Package tests: `test/fli-ofli/` (`fl_*.jl`). Verify script: `verify/fli.jl`, run as `julia -t 2 --project=verify verify/fli.jl` from the repo root; full output in `fli-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`). GPU: `verify/gpu.jl`, tag `fli` (`GPU_ONLY=fli julia --project=verify -t 4 verify/gpu.jl`), output in `gpu-fli-verify-log.txt`. Steps 1 and 3 (evidence part) of the stage contract; step 2 (independent review) is done by another agent, see "Independent review".

Run facts. Julia 1.13.2, 2 threads (GPU run 4 threads), base seed 20261012, verify offset 58000 (seeds `FL_SEED + 58000 + 50 ac + sub`; the package tests use `FL_SEED + 33000 .. 33999`, so all samples are fresh; the K = 0.5 regular class is screened again with a verify seed, Benettin X(1e5) <= 5e-4 as L1 S7). Working tree of branch `main`, nothing committed. Wall time of the final run 297 s, of which the package suite at full size 143.8 s (1 229 773 assertions, 0 fail, 0 error; the package-size suite has 263 533). The measurement helpers are those of `test/fli-ofli/fl_helpers.jl` (`fl_BFLI` two-point bound, `fl_ref_steps`, `fl_ref_at`, `fl_Psi`) and, for the oracle line only, `test/sali/sl_helpers.jl` (`sl_BSALI`); the loops are written in the script, so a definition error in a helper would not be caught by the script; the independent pieces are the 256-bit closed forms, the L1 `Benettin` estimator, the sali reducer and ChaosTools. Values are "worst ratio of error to bound" (pass at <= 1) unless stated. Final tally of the verify log: 228 PASS, 0 FAIL, 66 INFO. GPU log: 22 PASS, 0 FAIL, 3 FALLBACK (Double64), 30 INFO (55 lines).

File names follow the request of the caller (`fli-verify-log.txt`, `fli-*.png`), not the contract pattern `fli-ofli-*`; stage 05 may rename.

## Evidence table

"Where": T = package test (`test/fli-ofli/`, full size through the `suite` section), V = `verify/fli.jl`, G = `verify/gpu.jl`. Two entries per row are Float64 / Double64.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| suite | package suite at full size (`FL_FULL_RUN`: every multiple of 10, n = 1e6 shear, every n of AC5), assertions | 1 229 773 pass, 0 fail, 0 error | all pass | PASS | T via V |
| 1 | definition algebra, maps (i) to (iii), 40 handoffs: FLI_1, FLI_B against the exact running maxima over B_S(j) | 0.0318 / 0.00355 | 1 | PASS | T, V |
| 1 | FLI_1, FLI_B equal the running maxima of the reducer's own sums with the 0, bitwise; nondecreasing; column exchange (FLI_B bitwise, FLI_1 = other sum's maximum); sign change bitwise; spec values (FLI_B = 3 from j = 10 to 22, FLI_1 = 0 to j = 15, 0.5 at j = 17, (ii) FLI_B = 0.2 j) | 0 failures; worst FLI_1 / B_S on the zero stretches 0.0091 / 0.00036 | 0 | PASS | T, V |
| 1 | negative controls (no memory, no j = 0 term, no log) wrongly satisfying the predicate | 0 of 6 | 0 | PASS | T, V |
| 2 | K = 0 shear, seeded and W0 = I, tau in {1,10}, n in {10,100,1e4,1e6 (Float64)}: both forms against the closed form | 0.0198 / 0.00091 | 1 | PASS | T, V |
| 2 | W0 = I: FLI_1 exactly 0 (nonzero values 0 of 8 / 6); FLI_B against (1/2) ln(n^2 + 1) | 0.0143 / 0.00068 | 1 | PASS | T, V |
| 3 | linear hyperbolic map Y5, tau in {1,5}, n to 200: seeded ids (2 of 100 excluded) | 0.0435 / 0.00167 | 1 | PASS | T, V |
| 3 | supplied frames, bound max(B_S, 8 eps Psi_n) (D72 item 2): abs(c_1) = 1e-3; 1e-6 | 0.0071 / 0.00016; 0.0143 / 0.00020 | 1 | PASS | T, V |
| 3 | reported: abs(c_1) = 1e-6 over the plain B_S (why D72 item 2 was needed) | 21.1 / 6.63 | reported | INFO | V |
| 3 | (b) abs(c_1) = 1e-3, n = 4: FLI_1 == 0 exactly and FLI_B > 0.5 (the initial vector issue, F3); (c) rate over B + mu^(-4n)(1 + 1/c_max^2) | true; 0.0435 / 0.00167 | true; 1 | PASS | T, V |
| 4 | hyperbolic fixed point K = 1 / 2 / 6, 100 ids, tau in {1,5}: both forms against the BigFloat product norm | 0.0460 / 0.0856 / 0.0853 (Double64 0.0013 / 0.0032 / 0.0011) | 1 | PASS | T, V |
| 5 | (a) order-4 rotation, tau in {1,3}, n to 1e4: FLI != 0 at 0 of 13333 handoffs; (b) K = -2 at (0,0), every n in 1..1000: FLI_1, FLI_B against the closed form; (FLI_B - ln((3 + sqrt 5)/2))_+ ; (c) tau = 2: FLI / B_S | 0 failures; 0.0620 / 0.00194; 0; 0.0417 / 0.00195 | 1 | PASS | T, V |
| 6 | same-orbit reference, 20 starts, every multiple of 10 up to 1e4, tau = 1 / 4: FLI_B over B_FLI, Float64 K = 0.5 / 2 pi / 10 | 0.0248 / 0.0548 / 0.0723 (tau = 4: 0.0105 / 0.0515 / 0.0827) | 1 | PASS | T, V |
| 6 | same, Double64 | 0.00066 / **0.383** / 0.297 (tau = 4: 0.00028 / 0.383 / 0.297) | 1 | PASS (largest ratio of the (S) rows) | T, V |
| 6 | FLI_1 over B_FLI: at most 0.0827 (Float64), 0.383 (Double64); pairs with a differing computed and reference maximiser | 0 of 20000 per panel | reported | INFO | V |
| 6 | negative control, Float32 accumulator at n = 1e3, min error / B_FLI over 3 starts, K = 2 pi / 10 | 6.2e6 / 2.5e7 | > 1 | PASS | T, V |
| 7 | basis-free sandwich: (FLI_B - F_*)_+, (F_* - ln sqrt 2 - FLI_B)_+, (FLI_1 - F_*)_+ over B_FLI, all panels, both types | 0 | 1 | PASS | T, V |
| 7 | reported: (F_* - FLI_B) / ln sqrt 2 (min; median; max), F_* - FLI_1 (min; median; max): K = 0.5 / 2 pi / 10 | 7e-11; 0.19; 0.94 / 7.5e-5; 0.42; 0.95 / 0.012; 0.19; 0.997 and 2.4e-11; 0.49; 2.3 / 0.012; 0.45; 4.9 / 0.0041; 0.16; 2.2 | reported | INFO | V |
| 8 | (a) S_1, S_2 bitwise equal to the sali sums (sum and comp) and to the L1 sums of k = 1 runs from columns 1 and 2, 20 starts x 1000 checkpoints per K | 0 differences (K = 0.5, 2 pi, 10) | 0 | PASS | T, V |
| 8 | (b) FLI_1 against the running maximum of n X_L1, over 4 eps max(1, FLI_1), 20000 points per K | 0.228 / 0.225 / 0.177 | 1 | PASS | T, V |
| 8 | (c) FLI_B >= max of (-L/2) minus B_L/2; max(S_1, S_2) = -L/2 + abs(S_1 - S_2)/2 over B_L | 0 and 0.0060 / 0 and 0.0473 / 0 and 0.0689 | 1 | PASS | T, V |
| 8 | R1: drawdown D = FLI_1 - n X_L1 (min; median; max) at n = 1e2, 1e3, 1e4: K = 0.5 | 0; 0.28; 1.05 / 0; 0.54; 2.80 / 0.13; 0.31; 3.13 | reported | INFO | V |
| 8 | R1: K = 2 pi | median 0, max 0.61 / 0.89 / 1.97 | reported | INFO | V |
| 8 | R1: K = 10 | median 0, max 0 / 0 / 0.19 | reported | INFO | V |
| 8 | R1: FLI_B - FLI_1 (median; max), K = 0.5 / 2 pi / 10, n = 1e4 | 0.33; 1.88 / 0.19; 4.93 / 0; 2.19 | reported | INFO | V |
| 8 | R1 label demonstration at n = 1e3 (20 orbits): "FLI_B above ln n" / "n X_L1 above ln n", labels that differ: K = 0.5 / 2 pi / 10 | 13 / 5, 8 differ; 20 / 20, 0; 20 / 20, 0 | reported | INFO | V |
| 9 | (a) regular class, 20 starts, n in [1e4,1e6]: LS slope of FLI_B against ln n in [0.7,1.3] (min 0.9974, median 1.001, max 1.014); transient starts (D72 item 3) | 20 of 20; 0 | >= 16; <= 4 | PASS | T, V |
| 9 | (a) reported: starts with FLI_B below ln n at n = 1e4; FLI_B - ln n averaged per start (min; median; max) | 7 of 20; -2.41; 0.38; 2.68 | reported | INFO | V |
| 9 | (b) 1000 starts, abs(mean(slope over [20,200]) - mean(X(1e4))) / (4 SE + 4/180), K = 2 pi / 5 / 10 (means 1.169 and 1.170; 0.9645 and 0.9536; 1.620 and 1.620) | 0.0289 / 0.212 / 0.0182 | 1 | PASS | T, V |
| 9 | (b) realized abs(difference) of the means (allowance of the spec 0.0222, sali form 0.0111); z = difference / SE | 0.00109 / 0.01096 / 0.00059; z = 0.28 / 1.49 / 0.23 | allowance 0.0222 | INFO | V |
| 9 | (c) K = 5 / 10: FLI_1(1e3)/1e3 and FLI_B(1e3)/1e3 against X(1e4), over 4 SE + 4/n | 0.040, 0.060 / 0.083, 0.152 | 1 | PASS | T, V |
| 9 | R3: Chirikov ln(K/2) (K = 2 pi: 1.172 of the benchmarks row) against the mean slope | 1.172 vs 1.169; 0.916 vs 0.965; 1.609 vs 1.620 | reported (K18) | INFO | V |
| 9 | R4: first handoff at which FLI_B reaches ln 1e20 = 46.05 (min; median; max), starts never reaching it in n = 1e3: K = 2 pi / 5 / 10 | 29; 39; 1000, 3 / 33; 46; 1000, 23 / 22; 29; 43, 0 | reported | INFO | V |
| 10 | tau in {2,5,10} against 1, K = 0.5 / 2 pi, n in {100,1e3}, 20 starts: (a) sums over 16 eps Psi_n | 0.0045 / 0.0233 (Double64 0.00013 / 0.00095) | 1 | PASS | T, V |
| 10 | (b) (FLI_tau - FLI_1)_+ over B^tau + B^1; (FLI_1 - FLI_tau - floor(tau/2) ln sigma)_+ over B^tau + B^1 | 0.0023 / 0.0216 and 0; Double64 at most 0.0010 and 0 | 1 | PASS | T, V |
| 10 | reported: FLI_1 - FLI_tau (median; max) and its max over floor(tau/2) ln sigma(K): K = 0.5, tau = 2 / 5 / 10 | 0.020; 0.25 and 0.34 / 0.057; 0.68 and 0.47 / 0.12; 0.71 and 0.20 | reported | INFO | V |
| 10 | reported: K = 2 pi, tau = 2 / 5 / 10 | 0; 0.89 and 0.39 / 0; 0.89 and 0.19 / 0; 0.89 and 0.078 | reported | INFO | V |
| 11 | (a) Float64 vs Double64, difference of means of FLI_B/n over 4 SE: K = 2 pi, n = 1e3 (z = 1.00); K = 5, n = 1e3; K = 5, n = 1e4; K = 10, n = 1e4 | 0.250; 0.058; 0.064; 0.054 | 1 | PASS | T, V |
| 11 | (b) BigFloat(256) on its own orbit against the 512-bit reference, 50 starts, n in {10,..,60}, 300 rows, both forms | 0.099 | 1 | PASS | T, V |
| 11 | (c) BigFloat vs Double64, n = 300, 100 starts: K = 2 pi / 5 (z = 0.55 / 0.60) | 0.137 / 0.150 | 1 | PASS | T, V |
| 11 | (d) K = 10, n = 1e4, Float64, 5 starts: FLI_B finite, above 709.78 (range 1.60e4 to 1.64e4), `fli_unlogged` gives `nothing`; error over B_FLI; small values give exp(F) | true; 0.0221; true | true; 1 | PASS | T, V |
| 12 | (a) zero tangent (kind, orbit, step, backend, column 1 or 2, value 0, nothing stored, lowest id first; both chains; Double64 per-orbit; k = 1 ArgumentError); (b) NaN not masked; (c) frame check: norms 2, 3, 0, 1 + 1000 eps rejected, unit, dependent and 1 + eps accepted (largest defect 1 eps), seeded ids 1..1e4 / 1..1e3 accepted (largest defect 2 eps / 0.0625 eps against 8 eps), rejection before any map call; (d) diag(e^3, e^-3): tau = 200 no error (FLI_B = 599.8), tau = 300 error at handoff 300; (e) run record on {FLI} | all true | true | PASS | T, V |
| 13 | (a) orbit points bitwise equal across the six chains, 10 cases; (b) {FLI(replace=false), SALI}: FLI and SALI outputs equal the stand-alone values, final W equal; `propagate` equals the serial per-orbit body bitwise, 4096 orbits, tau in {1,4}, chunk 64 and 1000, both chains | 0 of 10; 0 of 10; true (8 of 8) | 0; true | PASS | T, V |
| 13 | (c) outputs hashed in child processes with 1, 2 and 4 threads, and chunk 64 and 1000: x, W and every final of chains {FLI}, {FLI(false), SALI}, tau in {1,4} (8 cases) | identical | identical | PASS | V |
| 13 | (d) bytes allocated by `propagate_orbit` after warm-up, n = 1e4, tau in {1,4}, both chains, both types | 0 | 0 | PASS | T, V |
| 13 | (e) column exchange and sign change, K = 2 pi and 0.5, both types, tau in {1,4}: failures | 0 of 40; 0 of 80 | 0 | PASS | T, V |
| 14 | (a) K = 0, wrapped vs lifted, 100 runs: differences (series, sums, maximisers) | 0 | 0 | PASS | T, V |
| 14 | (b) K = 0.5, n = 10 / 100, 50 starts, abs(FLI_w - FLI_l) over 64 K M_state j^2 eps (rho1 + rho2); excluded comparisons | 0.00017 / 0.00044; 1 / 2 of 100 | 1; <= 10 | PASS | T, V |
| 15 | ChaosTools 3.5.4 `gali`, oracle orbit against `forward`: bitwise differences (K = 0.5 20 x 1e4, K = 2 pi and 5 100 x 12) | 0 of 200000, 1200, 1200 | 0 | PASS | V |
| 15 | one-sided line, K = 0.5 regular, 80 pairs: max (-ln g_oracle/2 - FLI_B) over B_FLI + B_oracle/(2 g_oracle); smallest FLI_B - (-ln g_oracle/2) = 0.0527; pairs below | -207.7; 0.0527; 0 | 1 | PASS | V |
| 15 | K = 2 pi, 3765 resolved pairs (235 of 4000 excluded, 5.9 percent) | 0.0219; smallest margin -4.14e-4 (one pair, g_oracle = 7.35e-14) | 1; <= 10 percent | PASS | V |
| 15 | K = 5, 3946 resolved pairs (54 of 4000 excluded, 1.4 percent) | -0.107; smallest margin 2.8e-4 (at g_oracle = 2.7e-6) | 1; <= 10 percent | PASS | V |
| 16 | (a) regular class K = 0.5, Float64, 1024 orbits, n = 4e4, each backend against its own BigFloat(256) reference (final states reproduced bitwise 1024 of 1024 on both): FLI_B over B_FLI, CPU / GPU; FLI_1 | 5.80e-4 / 8.20e-4; 5.80e-4 / 8.10e-4 | 1 | PASS | G |
| 16 | (a) chaotic class K = 10, Float64, 1024 orbits, n = 8: FLI_B, CPU / GPU; FLI_1 | 0.133 / 0.137; 0.120 / 0.137 | 1 | PASS | G |
| 16 | (a) chaotic class, n = 10: FLI_B, CPU / GPU; FLI_1 | 0.116 / 0.116; 0.116 / 0.132 | 1 | PASS | G |
| 16 | (b) regular class, orbits with final states within 1e-9 (144 of 1024 enter): FLI_B and FLI_1 gap over B_FLI,cpu + B_FLI,gpu; maximiser equal on CPU and GPU for 144 of 144 | 0.00936; 0.00916 | 1 | PASS | G |
| 16 | (b) chaotic class, n = 8, all 1024 enter (largest 16 eps M_n Gamma_n = 1.06e-5): FLI_B and FLI_1 gap over B_FLI,cpu + B_FLI,gpu + max(D_orb at the two maximiser pairs) | 0.0160; 0.0160 | 1 | PASS | G |
| 16 | (b) chaotic class, n = 10, 989 of 1024 enter (35 beyond 1e-4, largest 5.2e-4, reported only, gap ratio at most 0.0017 over the full bound): same | 0.0129; 0.0129 | 1 | PASS | G |
| 16 | (b) chaotic looseness (INFO), n = 8: gap / D_orb term median / max (D_orb term median 3.4e-8, max 4.1e-6; 813 of 1024 nonzero gaps, median 2.3e-12); without the term 659 of 1024 orbits above 1 (max ratio 1.06e5); with the maximisers of the BigFloat series in the D_orb term | 6.2e-5 / 0.0160; 0.0160 | reported | INFO | G |
| 16 | (b) item 4 check on the device orbit, n = 8 / 10: max over orbits and k of abs(z_k^cpu - z_k^gpu)_inf / (16 eps M_k Gamma_k); against the literal 16 eps M_n Gamma_n (orbits above 1: 0 / 0) | 0.118; 0.0488 / 0.097; 0.0609 | 1; 1 | PASS | G |
| 16 | computed maximiser (stepB, colB) differs from the BigFloat series' on either backend; CPU and GPU maximisers differ (chaotic) | 0 of 1024 on each backend; 0 | reported | INFO | G |
| 16 | (c) K = 10, 1024 orbits, n = 1e4, tau = 2, FLI_B/n and FLI_1/n: max per-orbit gap / sigma_FT; difference of means / SE | 4.989; 0.607 and 4.991; 0.607 | 6; 4 | PASS | G |
| 16 | reported: final states apart by more than 1e-9, CPU vs GPU / CPU Float64 vs Double64: regular n = 4e4; chaotic n = 8; chaotic n = 10 | 880 / 876; 74 / 96; 324 / 592 (of 1024) | reported | INFO | G |
| 16 | Double64 (all classes and the ensemble) | kernel does not compile for the element type | recorded fallback reason is the pass | FALLBACK (pass) | G |

Recorded lines of the spec: R1 and R4 in the table above; R2 (realized ratios error / B_FLI over criteria 2 to 8, 11 and 16): at most 0.0856 (Float64, criterion 4), 0.0827 (Float64, criterion 6, K = 10, tau = 4), 0.383 (Double64, criterion 6, K = 2 pi), 0.137 (criterion 16); so the proposed constant 8 eps with the two-point rule holds with a margin of 2.6 or more everywhere, and of 12 or more in Float64; R3 (sanity rates) in the table, loose, K18.

## Oracle line (criterion 15), as the spec states it

The ChaosTools 3.5.4 oracle computes GALI_2, not FLI, so the line is one-sided through sali A7: -ln g_oracle / 2 = (S_1 + S_2) / 2 <= max(S_1, S_2) <= FLI_B. First the oracle orbit equals the stream's lifted orbit bitwise (0 differences in 202 400 steps), so no orbit divergence term is needed. Then the excess (-ln g_oracle/2 - FLI_B) over its bound B_FLI + B_oracle / (2 g_oracle) is at most 0.0219 (K = 2 pi), negative for K = 5 (-0.107) and K = 0.5 (-207.7): no pair has FLI_B below the oracle value by more than 4.1e-4, and that single pair (K = 2 pi, g_oracle = 7.35e-14, where B_oracle / (2 g) is large) is the only one of 7 791 resolved pairs. FLI_B is typically 0.1 to 1 above the lower value (figure `fli-chaostools.png`, left), because S_1 and S_2 differ.

## Findings that need a decision

1. Criterion 6, Double64, K = 2 pi (and K = 10): the realized ratio error / B_FLI is 0.383 (0.297 at K = 10), against 0.055 and 0.072 in Float64. The constant 8 eps with Psi borrowed from L1 holds with a factor 2.6 in Double64; the test plan (item 12) named this as the place a correct implementation could approach 1. Not a failure. Suggestion for the stage 05 spec note: keep 8 eps, record the realized 0.383 (K18).
2. Criterion 3, supplied frames: the plain B_S is exceeded for abs(c_1) = 1e-6 by a factor 21.1 (Float64) and 6.6 (Double64), as the test plan predicted; the D72 item 2 form passes with 0.0143 / 0.0002. The Psi form of the bound is what carries the criterion for ill-conditioned frames.
3. Criterion 9(b): the realized differences of the means (0.0011, 0.0110, 0.0006) are below the allowance 4/180 = 0.0222 by factors 20, 2 and 37; at K = 5 the sali form 2/180 = 0.0111 would be passed by 1.4 percent only (0.01096). D72 item 8 kept 4/180; no change needed, the allowance could be replaced by 0.015 (K18).
4. Criterion 9(a): 7 of 20 regular starts have FLI_B below ln n at n = 1e4 while all 20 have slopes in [0.997, 1.014] over [1e4, 1e6], which confirms the reason for D72 item 3.
5. Criterion 16: Float64 runs on the device (`provenance.backend == :gpu` on every Float64 GPU line; the FLI reducer kernel compiles); the Double64 kernel does not compile (probe `compiles = false`), so the Double64 lines are the recorded fallback reason, the allowed outcome. The chaotic class is gated at n = 8 (all orbits) and n = 10 (989 of 1024; 35 beyond the item 4 horizon 1e-4 are reported). The D_orb term is loose by about four orders (gap / term median 6e-5, max 0.016) and carries the gate: without it 659 of 1024 orbits at n = 8 and 716 of 989 at n = 10 exceed the bound B_FLI,cpu + B_FLI,gpu (max ratio 1.06e5 and 7.9e6). Realized |FLI_ref,cpu - FLI_ref,gpu| (BigFloat series on the two orbits) is itself within the full bound (max ratio 0.016 and 0.013), so the term is not hiding a defect of the reducer. Maximisers: the computed (stepB, colB) equals the BigFloat series' for every orbit on both backends and CPU and GPU maximisers agree on all 144 regular and all chaotic orbits, so the two-point rule did not need its second point in practice.
6. Criterion 13(c): the thread-invariance check (1, 2, 4 threads, child processes, chunk 64 and 1000) passes; it was a stage 04 item by D72 item 6.
7. Criterion 8 R1 (finding F2 of the spec): on the K = 0.5 regular class at n = 1e3, the labels "FLI_B above ln n" and "n X_L1 above ln n" differ for 8 of 20 orbits (13 against 5 above the threshold); on K = 2 pi and K = 10 they never differ and the drawdown median is 0. The memory term matters on the regular class only, and the two thresholds of Maffione2011a Table 1 are not interchangeable there (not a pass criterion; thresholds are pre-registered per study).
8. Not measured here: runtime benchmark on the device (impl notes: GPU not run), coverage and size claims of the impl notes.

## Margins and marginal criteria

Largest ratios in the (S) rows with bound 1: 0.383 (criterion 6, Double64, K = 2 pi), 0.297 (criterion 6, Double64, K = 10), 0.250 (criterion 11(a), K = 2 pi, z = 1.00 of 4), 0.228 (criterion 8(b), K = 0.5), 0.212 (criterion 9(b), K = 5), 0.152 (criterion 9(c)); all others below 0.15. Criterion 3 supplied frames with the plain B_S are not a criterion (INFO). The nonzero-bound structure of criterion 1 (FLI_1 = 0 asserted up to j = 15, D72 item 4) holds with ratio 0.0091. Criterion 9(a) slopes lie in [0.997, 1.014], inside [0.7, 1.3].

## Deviations from the impl notes, checked

- Zero allocations in the per-orbit body (impl notes, API) confirmed: 0 bytes at n = 1e4, tau in {1,4}, both types, both chains.
- Bitwise equality of the sums with sali and with the L1 sums per column (I6, sali I8) confirmed on 3 x 20 x 1000 checkpoints.
- `check_frame` threshold 8 eps (D72 item 7): seeded frames have defects at most 2 eps (Float64); columns of norm 1 + 1000 eps are rejected.
- Isbits reducer and state: the Float64 kernel path compiles and runs on the device (not stated in the notes, which expected a possible fallback); Double64 does not compile.
- The chain rule is not testable from the reducer (D72 item 1); verified here that `FLI(replace = false)` before `SALI` changes nothing in either output (criterion 13(a), (b)).
- Coverage and size claims of the impl notes were not re-measured here.

## Runtime

`verify/fli.jl`: 297 s total (package suite at full size 143.8 s, criteria 6 and 7 51.6 s, criterion 8 35.1 s, criterion 11 12.2 s, criterion 10 10.2 s, criterion 13 8.7 s with the child processes, criteria 3 to 5 7.5 s, criterion 15 3.8 s, the two figures 12.8 s, the rest a few seconds). GPU run (`GPU_ONLY=fli`, 4 threads): about 4 minutes wall (242 s from file creation to last line; BigFloat(256) references for 1024 orbits per backend at n = 4e4, and D_orb for 2 x 1024 chaotic orbits). The machine was not otherwise idle (the two runs overlapped for part of the time), so timings are not benchmarks.

## Figures

- `fli-growth-laws.png`: left, eight regular K = 0.5 orbits, FLI_B against ln n up to n = 1e6, parallel to the slope-1 line ln n (Q8); middle, ensemble means (200 starts) of FLI_B(n) for K = 5, 2 pi, 10 against lambda_1 n with lambda_1 the L1 mean X(1e4) of the identical orbits (0.922, 1.174, 1.62); right, FLI_1 of the reducer against the running maximum of the L1 sum n X(n) (Benettin, k = 1) for K = 0.5 and 2 pi, on the diagonal (worst abs(diff) / (4 eps max(1, FLI_1)) = 0.237, criterion 8(b); the sums themselves are bitwise equal, criterion 8(a)).
- `fli-chaostools.png`: left, FLI_B against -ln(g_oracle)/2 of ChaosTools `gali` on the resolved pairs of criterion 15 (K = 0.5, 2 pi, 5), above the diagonal; right, the margin FLI_B - (-ln g_oracle / 2) and the bound B_FLI + B_oracle / (2 g_oracle) against -ln(g_oracle)/2 (the margin lies far above the bound; the single pair below the diagonal is the point at the bottom).

## Files

- `verify/fli.jl` (new), `verify/gpu.jl` (new tag `fli`, section "AC16 (fli-ofli)" before the summary: functions `d16_one`, `d16_orb`, `d16_chaotic_b`, `fli_d16`, `fli_est`; it reuses `d68_sig2`, `gpu_trace` and `ac12_chaotic` of the sali and shared sections).
- `fli-verify-log.txt`, `gpu-fli-verify-log.txt`, `fli-growth-laws.png`, `fli-chaostools.png`, this report.
- Nothing under `src/` or `test/` was edited; nothing was committed.

## Independent review

Fresh Sonnet subagent, given the spec (with the checkpoint resolutions and Amendment D72) and `src/Indicators/fli.jl`, `Indicators.jl`, plus the Tangent protocol, Foundation, `sali.jl` and `benettin.jl` for comparison; not the notes, the test plan or the tests. Probes in the session scratchpad only. Result: no CRITICAL or HIGH finding.

Confirmed by probes: the sums and compensation terms are bitwise equal to SALI's per column and to the k = 1 Benettin sum (Float64, Double64; K in {0.5, 2 pi, 10}; tau in {1, 4}); the chain {non-replacing FLI, SALI} gives stand-alone results bitwise and leaves the orbit unchanged; the maxima and maximisers match an independent recomputation with zero mismatches (last strict increase, s = 0 term); health errors carry the right kind, column, orbit and step; `@inferred` for Float64, Double64, BigFloat; allocation does not grow with n; a forced relaunch every 50 steps is bitwise equal to one launch; the frame check accepts a defect of exactly 8 eps and rejects NaN and zero columns.

| # | Severity | Where | Finding | Disposition |
|---|----------|-------|---------|-------------|
| R1 | MEDIUM | `fli.jl`, provenance | The run record does not show the replace option or the chain position; in a chain the later reducer's indicator name wins (D34 A2), so nothing records that FLI ran | The chain and its options are recorded by Experiments (D67 item 2, D72 item 1); proposed for stage 05: say so in the docstring and the Indicators CONTEXT limits |
| R2 | MEDIUM | `fli.jl`, replace option | A wrong chain position gives silently wrong numbers: the replacing FLI before SALI leaves SALI near zero sums; the non-replacing FLI alone sums unnormalised columns (5867.4 against the correct 110.7 at K = 2 pi, n = 100) | By design (D72 item 1); a reducer cannot see its position. Proposed for stage 05: a docstring warning and the Indicators CONTEXT limits; the Experiments chain builder must enforce the one-replacer rule (noted for the experiments module, with the D67 frame check) |
| R3 | LOW | `fli.jl`, unlogged helper | NaN input returns `nothing`, read as overflow | The reducer never produces NaN; documented |
| R4 | LOW | `fli.jl`, frame check | Integer matrices raise MethodError | Float frames are the contract; documented |
| R5 | LOW | `Indicators.jl` docstring | Still says the layer carries its own scaled norm; since D69 the norm is in Foundation | Fix the sentence at stage 05 (docstring only) |
| R6 | LOW | spec Compute | State is 96 bytes in Float64 and 152 in Double64, not the proposed 64; chunk sizing reads the real size | Spec note at stage 05, as sali R6 |
| R7 | LOW | `fli.jl`, maximiser | The maximiser is reported as the map step n = j tau, not the handoff index j | Consistent with L1 (steps); state it in the docstring at stage 05 |
| R8 | LOW | `fli.jl`, init | Interpolated message in the k check inside a kernel path | Same pattern as sali and Benettin; the Float64 kernel compiled and ran (criterion 16) |
