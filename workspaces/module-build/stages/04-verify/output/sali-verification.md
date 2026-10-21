# Verification: sali

Spec: `../../01-spec/output/sali-spec.md` (criteria 1 to 18, Amendments D67 and D68). Test plan: `../../02-tests/output/sali-test-plan.md`. Impl notes: `../../03-implement/output/sali-impl-notes.md`. Package tests: `test/sali/` (`sl_*.jl`). Verify script: `verify/sali.jl`, run as `julia -t 2 --project=verify verify/sali.jl` from the repo root; full output in `sali-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`). GPU: `verify/gpu.jl`, tag `sali` (`GPU_ONLY=sali julia --project=verify -t 4 verify/gpu.jl`), output in `gpu-sali-verify-log.txt`. Steps 1 and 3 (evidence part) of the stage contract; step 2 (independent review) is done by another agent, see "Independent review".

Run facts. Julia 1.13.2, 2 threads (GPU run 4 threads), base seed 20261012, verify offset 56000 (seeds `SL_SEED + 56000 + 50 ac + sub`; the package tests use `SL_SEED + 31000 .. 31999`, so all samples are fresh; the K = 0.5 regular class is screened again with a verify seed, Benettin X(1e5) <= 5e-4 as L1 S7). Working tree of branch `main` at HEAD 6253473, nothing committed. Wall time of the final run 168 s, of which the package suite at full size 61.8 s (9210 assertions, 0 fail, 0 error). The measurement helpers are those of `test/sali/sl_helpers.jl` (`sl_BSALI`, `sl_BL`, `sl_Lambda`, `lb_Psi`, the BigFloat reference `sl_ref`); the loops are written in the script, so a definition error in a helper would not be caught by the script; the independent pieces are the 256-bit closed forms, the L1 `Benettin` estimator and ChaosTools. Values are "worst ratio of error to bound" (pass at <= 1) unless stated. Final tally of the verify log: 185 PASS, 0 FAIL, 30 INFO. GPU log (original criterion 18): 5 PASS, 2 FAIL, 2 FALLBACK, 4 INFO. GPU log under Amendment D68, first run (`gpu-sali-d68-verify-log.txt`): 11 PASS, 2 FAIL (both the chaotic-class part (b), item 1 of the findings), 3 FALLBACK (Double64), 10 INFO. GPU log under D68 with the chaotic-class addition of 2026-10-22 (`gpu-sali-d68b-verify-log.txt`, the current record of criterion 18): 14 PASS, 0 FAIL, 3 FALLBACK (Double64), 12 INFO (29 lines).

## Evidence table

"Where": T = package test (`test/sali/`, full size through the `suite` section), V = `verify/sali.jl`, G = `verify/gpu.jl`. Two entries per row are Float64 / Double64.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| suite | package suite at full size (`SL_FULL_RUN`), assertions | 9210 pass, 0 fail, 0 error | all pass | PASS | T via V |
| 1 | (a) d-, d+ against BigFloat(256), over B_meas, 11000 pairs (d = 2 and 4) | 0.75 / 0.04286 | 1 | PASS | T, V |
| 1 | (b) abs(d-^2 + d+^2 - 4) / (16 eps), Amendment D67 | 0.75 / 0.04286 | 1 | PASS | T, V |
| 1 | (b) same against the original 8 eps of the spec text (reported) | 1.5 / 0.08571 | 1 | INFO: the amendment was needed in Float64 | V |
| 1 | (c) GALI_2 relation, over 8 eps | 0 / 0 | 1 | PASS | T, V |
| 1 | (d) abs(SALI - exact) / B_meas, random pairs; 28 constructed pairs | 0.5389 / 0.0211; 0.0228 / 0.000637 | 1 | PASS | T, V |
| 1 | (e) bitwise change under sign change or exchange; (f) exact zeros of 50 | 0 failures; 50 of 50 | 0; 50 | PASS | T, V |
| 1 | (g) orthogonal pairs, abs(SALI - sqrt 2) / (4 eps); negative control d- alone, antiparallel | 0 / 0.003906; 2 | 1; > 1.9 | PASS | T, V |
| 2 | tau in {2,5,10} against tau = 1, n in {10,100}, 20 starts; SALI: K = 0.5 / K = 2 pi (Float64); L: K = 0.5 / 2 pi | 0.0109 / 0.0426; 0.0142 / 0.0481 (Double64 at most 0.0016) | 1 | PASS | T, V |
| 3 | K = 0 shear, SALI against closed form, tau in {1,10}, n in {10,100,1e4,1e6} (Double64 at 1e6 included) | 0.00715 / 0.000113 | 1 | PASS | T, V |
| 3 | L against closed form | 0.00817 / 0.000285 | 1 | PASS | T, V |
| 3 | W0 = I, abs(SALI n - 1) / (3/(8 n^2) + n B_SALI) | 0.9999 / 1 (the bound is a leading-order truncation, ratio is 1 minus O(1/n^2) by construction) | 1 | PASS | T, V |
| 3 | reported: SALI n^2 abs(b1 b2) / abs(det W0) at n = 1e4 (A10: 1 - (a1/b1 + a2/b2)/n) | 1 | reported | INFO | V |
| 4 | linear hyperbolic map: (a) resolved n, 600 points; (b) L, n in {20,100,200}; (c) rate; 100 ids, 0 excluded, both types | (a) 0.0717 / 0.0027; (b) 0.0955 / 0.0024; (c) 0.0955 / 0.0022 | 1 | PASS | T, V |
| 4 | BigFloat product against the A11 eigen-expansion at n = 8 (max L difference) | 2.8e-76 / 1.4e-76 | < 1e-60 | PASS | V |
| 5 | hyperbolic fixed point K = 1 / 2 / 6: (a) SALI, resolved n | 0.0736 / 0.0817 / 0.0879 (Double64 0.0027 / 0.0026 / 0.0037) | 1 | PASS | T, V |
| 5 | (b) L against the A12 closed form with the explicit tail 1.01 (c + 1) mu+^(-2n) (D67 item 4) | 0.0927 / 0.172 / 0.180; Double64 0.0023 / **0.986** / 0.0033 | 1 | PASS (marginal at Double64 K = 2, n = 20: the tail term dominates; frame coefficient c = 241.8) | T, V |
| 6 | (a) order-4 rotation, abs(SALI(n) - SALI(0)) / B_SALI, n up to 9996, tau in {1,3} | 0.00112 / 5.6e-6 | 1 | PASS | T, V |
| 6 | (b) K = -2 at (0,0), all n in 1..1000, 100 ids, tau in {1,2}: against the exact product; even n against SALI(0) | 0.0356 / 0.00154; 0.0190 / 0.00071 | 1 | PASS | T, V |
| 6 | (b) min over n and ids of SALI(n) - (0.146 - B_SALI) (smallest SALI observed 0.2887) | 0.1427 | >= 0 | PASS | T, V |
| 7 | direct form against the BigFloat reference, 20 starts, every multiple of 10: K = 0.5 (n to 1e4) / 2 pi / 10 (n to 1e3); Float64 | 0.0143 / 0.0860 / 0.0576 | 1 | PASS | T, V |
| 7 | same, Double64 | 0.00059 / 0.0036 / 0.0019 | 1 | PASS | T, V |
| 7 | floor control, chaotic: handoffs with SALI_ref < B_meas where SALI_T > 2 B_meas (1978, 1989, 1953, 1970 handoffs below the floor); resolved handoffs exist | 0 | 0 | PASS | T, V |
| 8 | (a) L against the reference, tau in {1,4}: Float64 K = 0.5 / 2 pi / 10 | 0.0184 / 0.0740 / 0.169 (tau = 4: 0.0114 / 0.0766 / 0.114) | 1 | PASS | T, V |
| 8 | (a) Double64 | 0.00047 / 0.0258 / 0.0167 | 1 | PASS | T, V |
| 8 | (b) K = 2 pi, Float64, all handoffs: L finite, min L = -2460, exp(min L) = 0 | finite | finite, < -745 | PASS | T, V |
| 8 | (c) forms agree above the floor, worst ratio over all panels | 0.063 (Float64), 0.0030 (Double64) | 1 | PASS | T, V |
| 8 | (d) negative control, S_1 only, min abs(-S_1 - L_ref) / (10 B_L), n >= 50 | 7.7e11 (Float64) to 2.5e26 (Double64) | > 1 | PASS | T, V |
| 9 | cross-check with L1 on the identical orbit, K = 0.5 / 2 pi / 10: (a), (b) S_i against n X_L1 | 0 (bitwise equal) | 1 | PASS | T, V |
| 9 | (c) L + n (X1 + X2), over 2 B_L | 0.00038 / 0.0240 / 0.0462 | 1 | PASS | T, V |
| 9 | reported: abs(ln(norm w2 / norm w1))/2 at n = 1e4 (median; max) | 0.33; 2.5 / 0.36; 2.8 / 0.41; 1.6 | reported | INFO | V |
| 10 | chaotic rate, 1000 starts: abs(mean(m/2) - mean(X(1e4))) / (4 SE + 2/180), K = 2 pi / 5 / 10 | 0.0465 / 0.111 / 0.207 | 1 | PASS | T, V |
| 10 | realized difference of the means (allowance of the spec 0.0111); z = difference / SE | 0.0014 / 0.0044 / 0.0044; z = 0.29 / 0.62 / 1.74 | allowance 0.0111 | INFO | V |
| 11 | regular decay, 20 starts, slope of L against ln n over [1e4,1e6] in [-2.4,-1.6] (min -2.185, median -1.99, max -1.90); transient starts | 20 of 20; 0 | >= 16; <= 4 | PASS | T, V |
| 11 | R2: n^2 SALI at n = 1e4 / 1e6 (min; median; max) | 0.17; 4.3; 2724 / 0.38; 3.8; 27520 | reported | INFO | V |
| 12 | (a) median n_floor(Double64) / n_floor(Float64), 100 starts, n = 150: K = 5 (98 reach the floor, 2 islands) / 2 pi | 1.882 / 1.879 | [1.5, 2.4] | PASS | T, V |
| 12 | (b) BigFloat(256) against the 512-bit reference, 200 resolved points | 0.0837 | 1 | PASS | T, V |
| 12 | (c) Float64 vs Double64, n = 1e3, 1000 starts, difference over 4 SE: K = 2 pi / 5 | 0.0933 / 0.0046 | 1 | PASS | T, V |
| 12 | (c) BigFloat vs Double64, n = 300, 100 starts: K = 2 pi / 5 | 0.276 / 0.019 | 1 | PASS | T, V |
| 12 | R1: n_floor medians (Float64 / Double64) K = 5: 17 / 34, K = 2 pi: 15 / 29, K = 10: 11 / 21; A14 loose n_res 18/36, 14/28, 10/20 | | reported | INFO | V |
| 13 | (a) zero tangent: kind, orbit, column, step, backend, value, nothing stored, lowest id first; k = 1 ArgumentError. (b) NaN not masked. (c) frame check (dependent, 5 B_meas rejected; 20 B_meas accepted; seeded ids 1..1e4 / 1..1e3 accepted, min SALI(0) / (10 B_meas) = 1.6e14 / 4.5e28), rejection before any map call. (d) overflow tau 200 / 300. (e) run record and floor | all true | true | PASS | T, V |
| 14 | (a) orbit points bitwise equal across the chains, 10 cases; (b) `propagate` equals serial `propagate_orbit` bitwise, 4096 orbits, tau in {1,4}, chunk 64 and 1000; thread counts 1, 2, 4 give identical output hashes (child processes) | 0 differences; true; true | 0; true | PASS | T, V |
| 14 | (c) bytes allocated, n = 1e4, tau in {1,4}, both types | 0 | 0 | PASS | T, V |
| 14 | (d) absorbing zero, 1000 starts: 995 reach exactly 0 (5 regular islands), violations | 0 | 0 | PASS | T, V |
| 15 | (a) K = 0, wrapped vs lifted bitwise, 100 runs | 0 differences | 0 | PASS | T, V |
| 15 | (b) K = 0.5, n = 10 / 100, 50 starts, abs(L_w - L_l) / (64 K M_state n^2 eps (rho1 + rho2)); excluded starts 0 / 1 | 0.00031 / 0.00030 | 1 | PASS | T, V |
| 16 | ChaosTools `gali`, k = 2, K = 0.5 regular, 20 starts, n in {10,100,1e3,1e4}: abs(g(SALI by A4) - g_oracle) / (B_SALI + B_oracle); with g = exp(L) | 0.00832; 0.00809 (largest absolute difference 3.6e-14) | 1 | PASS | V |
| 16 | K = 2 pi, 1000 starts, resolved n in {2,5,8,12}: same ratio (A4 form; exp(L) form 0.102); unresolved pairs excluded 222 of 4000 | 0.110 | 1; <= 10 percent (5.6 percent) | PASS | V |
| 16 | K = 5, 1000 starts, same (exp(L) form 0.0948); excluded 56 of 4000 | 0.0974 | 1; <= 10 percent (1.4 percent) | PASS | V |
| 16 | (c) oracle orbit against `forward(StandardMap)`: steps with a bitwise difference, K = 0.5 (20 x 1e4), 2 pi and 5 (100 x 12) | 0 of 202400 | 0 | PASS | V |
| 17 | Float64 vs Double64 log form, n = 1e4, 200 starts, difference over 4 combined SE: K = 5 / 10 (z = 0.058 / 0.99; means 0.9719, 0.9714 and 1.620, 1.621) | 0.0146 / 0.248 | 1 | PASS | V |
| 18 | Amendment D68 (a), regular class K = 0.5, Float64, 1024 orbits, n = 4e4, each backend against its own BigFloat(256) reference (CPU orbit by `forward`, GPU orbit by a device trace kernel, final states reproduced bitwise for 1024 of 1024 on both): max abs(L_b - L_ref,b) / (8 eps (Psi1 + Psi2)), CPU / GPU | 4.83e-4 / 5.98e-4 | 1 | PASS | G |
| 18 | (a) regular class, direct SALI, max abs(SALI_b - SALI_ref,b) / B_SALI(n), 1023 of 1024 resolved, CPU / GPU | 7.10e-4 / 6.63e-4 | 1 | PASS | G |
| 18 | (a) chaotic class K = 10, Float64, 1024 orbits, n = 8, same references (final states reproduced for 1024 of 1024): L ratio, CPU / GPU | 0.183 / 0.161 | 1 | PASS | G |
| 18 | (a) chaotic class, direct SALI ratio, 953 of 1024 resolved, CPU / GPU | 0.102 / 0.0771 | 1 | PASS | G |
| 18 | (b) regular class, orbits with final states within 1e-9 in every coordinate (144 of 1024 enter): max abs(L_cpu - L_gpu) / (16 eps (Psi1 + Psi2)); max abs(SALI_cpu - SALI_gpu) / (2 B_SALI) | 0.0143; 5.10e-4 | 1; 1 | PASS | G |
| 18 | (b) chaotic class K = 10, n = 8, first D68 run (orbits within 1e-9, 950 of 1024 enter, 903 resolved; no orbit-divergence term): L gap ratio; SALI gap ratio | 6.78e4; 3.60 | 1; 1 | FAIL (superseded by the next three rows; log `gpu-sali-d68-verify-log.txt`) | G |
| 18 | (b) chaotic class, addition of 2026-10-22: entering orbits (16 eps M_n Gamma_n <= 1e-4; all 1024 enter, 0 reported only, largest 16 eps M_n Gamma_n = 1.06e-5): max abs(L_cpu - L_gpu) / (16 eps (Psi1 + Psi2) + D_orb(n)) | 0.0141 | 1 | PASS | G |
| 18 | (b) chaotic class, same orbits, 953 resolved: max abs(SALI_cpu - SALI_gpu) / (2 B_SALI(n) + SALI_ref(n) D_orb(n)) | 0.0598 | 1 | PASS | G |
| 18 | (b) chaotic class, item 4 check on the device orbit (precision spec, Equations 4), 1024 orbits: max over orbits and k = 1..n of abs(z_k^cpu - z_k^gpu)_inf / (16 eps M_k Gamma_k); same against the literal 16 eps M_n Gamma_n | 0.118; 0.0488 | 1; 1 | PASS (orbits above 1: 0 / 0) | G |
| 18 | (b) chaotic class, looseness of the term (INFO): gap / D_orb(n) with gap = abs(L_cpu - L_gpu), median / max over 1024 orbits (853 nonzero gaps, median among them 1.39e-4); D_orb median / max; the gap without D_orb, 663 of 1024 orbits above 16 eps (Psi1 + Psi2), max ratio 1.05e5 | 6.06e-5 / 0.0141; D_orb 8.33e-8 / 1.95e-5 | reported | INFO | G |
| 18 | (c) K = 10, 1024 orbits, n = 1e4, tau = 2, -L/(2n): max per-orbit gap / sigma_FT; difference of means / SE | 4.99; 0.61 | 6; 4 | PASS | G |
| 18 | Reported, not gated: orbits with final states apart by more than 1e-9, CPU vs GPU / CPU Float64 vs Double64 CPU: regular n = 4e4; chaotic n = 8 | 880 / 876; 74 / 96 (of 1024) | reported | INFO | G |
| 18 | Double64 (both classes) | kernel does not compile for the element type | recorded fallback reason is the pass | FALLBACK (pass) | G |

Recorded lines of the spec: R1 and R2 in the table above; R3 (sanity rates): mean rate 2 x mean(m/2) = 2.329 / 1.923 / 3.248 against 2 ln(K/2) = 2.289 / 1.833 / 3.219 at K = 2 pi / 5 / 10 (loose, K18); R4 (realized ratios) are the ratio values of rows 7, 8, 12 and 16: at most 0.169 against B_SALI or B_L in Float64 and 0.026 in Double64 (K = 10 and K = 2 pi), so the proposed constants 8 eps and 4 eps hold with a margin of 6 or more.

## P8 confirmations (before the numbers of criterion 16)

1. Threshold. The default `threshold = 1e-12` of `gali` stops the run: for K = 5, 200 starts, 100 steps requested, 195 starts stop early (after 9 to 14 steps, median 14); the other 5 are regular islands that run the full 100. With `threshold = 0` all 200 runs have the full length (101 values). So threshold = 0 is needed, as the spec says. Confirmed.
2. Prescribed initial frame. `reinit!(tds, u; Q0 = W0)` stores W0 (`current_deviations(tds) == W0` for 200 of 200 starts); `gali` value at step 1 equals our exp(L(1)) to 6.7e-16 (4 eps = 8.9e-16); the default identity frame gives a g(1) that differs by more than 1e-3 for 200 of 200 starts, so the frame really is used. Confirmed; the oracle assumes an orthonormal start (the recorded value 1 at n = 0), which the seeded frame satisfies.
3. Orbit. The ChaosTools lifted rule with a hand-coded Jacobian and the Chaoslab `forward` give bitwise identical orbits at every step: 0 differences in 202400 steps (K = 0.5: 20 starts x 1e4; K = 2 pi and 5: 100 starts x 12). Confirmed; no orbit divergence term is needed in the bound of criterion 16.
Also recorded: the oracle floor. At n = 30, K = 5, 100 starts, `gali` returns values from 0 to 0.206 (median 1.1e-16) where our exp(L) is of order 1e-25: the oracle is blind below about 1e-15, as A14 says; our log form is not.

## Findings that need a decision

1. Criterion 18: resolved by Amendment D68 (owner, 2026-10-21) and its chaotic-class addition (owner, 2026-10-22); no failure remains. The original finding (the CPU and GPU do not produce the same orbit, so the cross-backend same-orbit bound of 16 eps could not hold: regular class, n = 4e4, final x bitwise equal for 0 of 1024 orbits, ratio 1.76e6; chaotic window, SALI gap 3.59 times 2 B_SALI; log `gpu-sali-verify-log.txt`, kept as the record of the original criterion) is answered by the per-backend form. First D68 run (`gpu-sali-d68-verify-log.txt`, all 1024 orbits of both classes, no subsample): (a) each backend against its own BigFloat(256) reference passes with large margins in both classes (L ratio at most 5.98e-4 regular and 0.183 chaotic, SALI ratio at most 7.10e-4 and 0.102); the reference orbit reproduces each backend's final state bitwise for 1024 of 1024 orbits (the GPU orbit is read with a device kernel that runs `forward`). (c) the K = 10 ensemble form passes (4.99 sigma_FT, 0.61 SE). (b) regular class: 144 of 1024 orbits have final states within 1e-9, gap ratios 0.0143 (L) and 5.10e-4 (SALI), PASS. Divergence counts (reported): final states apart by more than 1e-9 for 880 of 1024 regular orbits CPU vs GPU, 876 CPU Float64 vs Double64 CPU (the orbit sensitivity is that of the class, not of the GPU); for the chaotic class 74 and 96 of 1024. The one failure of that run was (b) on the chaotic class, n = 8, K = 10: 950 of 1024 orbits entered (final states within 1e-9), the L gap was 6.78e4 times 16 eps (Psi1 + Psi2) (589 of the 950 above 1) and the SALI gap 3.60 times 2 B_SALI. Cause: the final-state threshold of 1e-9 does not select same-orbit pairs for a chaotic class; the orbit points already differ earlier in the run (largest point difference over all steps, median 3.6e-13, max 3.5e-10 over the agreeing orbits) and the tangent depends on the point through K cos x. Resolved with the D_orb addition (owner, 2026-10-22; spec Amendment D68, "Addition to part (b) for the chaotic class"; log `gpu-sali-d68b-verify-log.txt`, same orbits and seeds, rerun of the whole tag). D_orb(n) is computed in BigFloat(256) on the CPU reference orbit with the seeded frame, per orbit, from M_k, Gamma_k (precision spec, Equations 3), P_{<k} w_i, norm(P_{>k})_2 and norm(P w_i), in `d68_orb` of `verify/gpu.jl`. Result: all 1024 orbits enter (the largest 16 eps M_n Gamma_n is 1.06e-5, below the horizon 1e-4; none is reported only). L gap: max abs(L_cpu - L_gpu) / (16 eps (Psi1 + Psi2) + D_orb) = 0.0141, PASS. SALI gap in the resolved window (953 of 1024): max abs(SALI_cpu - SALI_gpu) / (2 B_SALI + SALI_ref D_orb) = 0.0598, PASS. Item 4 holds for the device: max over orbits and k = 1..n of abs(z_k^cpu - z_k^gpu)_inf / (16 eps M_k Gamma_k) = 0.118, and against the literal 16 eps M_n Gamma_n = 0.0488 (0 orbits above 1 in both), so the term is valid. The term is loose: gap / D_orb has median 6.06e-5 and max 0.0141 (853 of 1024 orbits have a nonzero gap, median 1.39e-4 among them; D_orb median 8.33e-8, max 1.95e-5), because D_orb is built from the worst-case bound 16 eps M_k Gamma_k on the orbit difference while the realized difference is at most 0.118 of it. Without the term the same gap is above 16 eps (Psi1 + Psi2) for 663 of 1024 orbits (max ratio 1.05e5), so the term, not the margin of the old bound, is what carries the gate. (a), (c), the regular class of (b) and the Double64 fallback are unchanged and reproduce the first D68 run line for line.
2. Float64 GPU only. The Double64 kernel does not compile (probe `compiles = false`), so the Double64 lines of criterion 18 are the recorded fallback reason "kernel does not compile for the element type", the allowed outcome. Device: NVIDIA GeForce RTX 5060 Laptop GPU, compute capability 12.0, CUDA.jl 6.4.2, runtime 13.4.0 (`_config/environment.md`, "Local GPU"); every Float64 line has `provenance.backend == :gpu`.
3. Criterion 1(b): the realized Float64 ratio against the original 8 eps is 1.5, so Amendment D67 item 5 (16 eps) was needed; the amended ratio is 0.75.
4. Criterion 12(a): at K = 5 two of 100 starts never reach the floor in 150 steps (regular islands); the median is taken over the starts that reach it in both types (98), and the gate requires at least 80. The spec text does not say how islands are handled.
5. Criterion 5(b): Double64, K = 2, n = 20 is at 0.986 of the bound (the tail term 1.01 (c + 1) mu+^(-2n) dominates; c = 241.8). Marginal, as D67 item 4 anticipated; Double64 at n = 20 only.
6. Criterion 3, W0 = I: the bound is a leading-order truncation, so the ratio is 1 minus O(1/n^2) for any correct implementation (0.9999 Float64); it cannot discriminate small errors. In Double64 at n = 1e6 the check must be formed in Double64 (the Float64 narrowing of SALI n exceeds the 3.75e-13 bound by rounding alone); the script does that.
7. Criterion 16 uses g = SALI sqrt(4 - SALI^2) / 2 (A4) as the gated conversion and also reports g = exp(L): the oracle's absolute floor (about 1e-16) and our direct form's floor are of the same size, so both forms agree with the oracle; the exp(L) form is the more accurate one (example at n = 12, K = 5: oracle 2.92341706e-11, direct 2.92341685e-11, exp(L) 2.92340796e-11, difference 9e-17).

## Margins and marginal criteria

Largest ratios in the (S) rows: 0.986 (criterion 5(b), Double64 K = 2), 0.9999 (criterion 3 identity, by construction), 0.75 (1(a), 1(b) Float64), 0.207 (10), 0.18 (5(b) Float64), 0.169 (8(a) Float64 K = 10); all others below 0.12. Criterion 10 with the finite-window allowance 2/180: realized differences of the means are 0.0014, 0.0044 and 0.0044, below it; the allowance could be replaced by 0.005 (K18, spec P6). Criterion 11 slopes lie in [-2.19, -1.90], inside [-2.4, -1.6]. Criterion 12(a) ratios 1.88 and 1.88 inside [1.5, 2.4].

## Deviations from the impl notes, checked

- `propagate` and `propagate_orbit` do not call `check_frame` (D67 item 1): verified by AC13(c): the caller-side check rejects a dependent frame before any map call and an accepted frame reaches the map. Experiments must call `check_frame` for supplied frames; nothing in this stage changes that.
- Zero allocations in the per-orbit body (impl notes) confirmed: 0 bytes at n = 1e4, tau in {1,4}, both types.
- Isbits state: the Float64 kernel path compiles and runs on the device; the Double64 kernel does not compile (not stated in the impl notes, which expected "Double64 falls back as for the other reducers": consistent).
- Coverage and size claims of the impl notes were not re-measured here.

## Runtime

Final run of `verify/sali.jl`: 168 s total (package suite 61.8 s, AC3 17 s, AC4 to 6 18 s, AC7 and 8 27 s, AC17 8 s, the two figures 16 s, the rest a few seconds; the thread-invariance child processes with 1 and 4 threads are inside AC14, 4.5 s). GPU run (`GPU_ONLY=sali`): 68 s for the original criterion 18; 291 s under D68 (BigFloat(256) references for 1024 orbits per backend, 4 threads); 318 s (5 min 18 s wall) under D68 with the chaotic-class addition (D_orb, same thread count). The machine was not otherwise idle, so timings are not benchmarks.

## Figures

- `sali-decay-laws.png`: left, one chaotic orbit at K = 2 pi: direct SALI in Float64 and Double64 reaching their floors (8.9e-16, 3.2e-30) while the log form continues, against exp(-2 lambda_1 n) with lambda_1 = X(1e4) = 1.182 of the same orbit; middle, six regular K = 0.5 orbits, SALI from the log form against the power law n^-2 (A9) up to n = 1e6; right, ensemble means (200 starts) of -L(n)/(2n) against the L1 mean X(1e4) at K = 5, 2 pi, 10.
- `sali-chaostools.png`: left, g = exp(L) against the ChaosTools `gali` value on the resolved points of criterion 16 (K = 0.5, 2 pi, 5), on the diagonal; right, the difference against its bound B_SALI + B_oracle (differences are 2 to 4 orders below the bound).

## Files

- `verify/sali.jl`, `verify/gpu.jl` (new tag `sali`, section "AC18 (sali)" before the summary; rewritten for D68); the chaotic class of part (b) rewritten on 2026-10-22: functions `d68_sig2`, `d68_orb`, `d68_chaotic_b`, keyword `chaotic` of `sali_d68`).
- `sali-verify-log.txt`, `gpu-sali-verify-log.txt` (original criterion 18, kept as record), `gpu-sali-d68-verify-log.txt` (criterion 18 as amended by D68, first run, kept as record), `gpu-sali-d68b-verify-log.txt` (criterion 18 with the chaotic-class D_orb addition, current), `sali-decay-laws.png`, `sali-chaostools.png`, this report.
- Nothing under `src/` or `test/` was edited; nothing was committed.

## Independent review

Fresh Sonnet subagent, given the spec (with the checkpoint resolutions and Amendment D67) and `src/Indicators/` plus the Tangent stream and `Lyapunov/benettin.jl` for comparison; not the notes, the test plan or the tests. Short probes in the session scratchpad only. Result: no CRITICAL or HIGH finding.

| # | Severity | Where | Finding | Disposition |
|---|----------|-------|---------|-------------|
| R1 | MEDIUM | `alignment.jl`, A4 conversion | The conversion from the log form to SALI is ill-conditioned when the pair is near orthogonal (SALI near sqrt 2): error about eps / cos(phi), up to about 0.35 sqrt(eps) at the limit; the information is not in L. Probe: K = 0.5 regular start, n = 260, direct 1.41199, converted value off by 1640 eps | No code change (not fixable in the conversion). Criterion 8(c) passed (worst 0.169), so the bound with Psi covers the realized cases. Proposed for the stage 05 spec note: state the condition number under A4 and in the docstring |
| R2 | MEDIUM | spec I2 / A3, `alignment.jl` | The computed SALI can exceed sqrt 2 by up to 2 eps (unit columns carry 1 to 2 eps norm error); the converted form is clipped at sqrt 2 | No code change (clipping would hide rounding). Proposed for the stage 05 spec note: range [0, sqrt 2 (1 + 2 eps)] |
| R3 | MEDIUM | `sali.jl`, final | The stream returns only the last handoff; series checks need prefix runs | Contract limit shared with Benettin (D66 keeps the stream unchanged); runs with different n are bitwise prefixes, which the review confirmed and criterion 14 tests. No change |
| R4 | LOW | `sali.jl`, final | L = -((s1 + s2) + (c1 + c2)) rounds s1 + s2 once without compensation (about 2e-12 absolute at K = 10, n = 1e4) | Within B_L in every run (AC8 worst 0.169, AC10 worst 0.207); optional TwoSum, not done |
| R5 | LOW | `alignment.jl`, floor | For BigFloat the floor follows the default precision at call time, not the data | Documented limit: call inside the same `setprecision` block; all verify runs do |
| R6 | LOW | spec Compute | The reducer state (72 bytes Float64, 120 Double64) is added to the memory formula; the spec line "0 bytes" is wrong | Proposed for the stage 05 spec note; chunk sizes are conservative, no risk |
| R7 | LOW | `alignment.jl`, pure functions | A zero or non-finite vector gives NaN without error (docstring says nonzero) | Documented; the stream path raises the health error |
| R8 | LOW | `sali.jl`, norm | Entries above floatmax / sqrt 2 overflow the scaled norm to a health error; subnormal entries lose direction precision | Academic (needs an extreme tau); shows as a health error, not a wrong number |
| R9 | LOW | `sali.jl`, init | Interpolated message in the k check inside a kernel path | The host-side check throws first; same pattern as Benettin; GPU Float64 compiled and ran (criterion 18) |
| R10 | LOW | D67 | Nothing in `src/` calls the frame check, so a dependent supplied frame reads as SALI = 0 | Approved by D67; the Experiments layer must call the frame check on every supplied frame (noted for the experiments module) |

Also noted (LOW, no change): the converted form clips L > 0 to g = 1; the initial state field reads 0 before the first handoff (unreachable through the stream); the frame check on integer matrices raises MethodError; the Foundation health-error docstring does not list the kinds "initial frame" and "deviation matrix".

The review also confirmed: math matches A2, A3, A4, A7 and the D67 items; zero allocation and `@inferred` for Float64 and Double64; W left column-normalised, the orbit untouched, no shared mutable state; A16 held on 100 K = 2 pi orbits; the SALI sums are bitwise equal to a Benettin k = 1 run per column.
