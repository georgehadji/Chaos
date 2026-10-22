# Verification: megno

Spec: `../../01-spec/output/megno-spec.md` (criteria 1 to 15, Amendment D74). Test plan: `../../02-tests/output/megno-test-plan.md`. Impl notes: `../../03-implement/output/megno-impl-notes.md`. Package tests: `test/megno/` (`mg_*.jl`). Verify script: `verify/megno.jl`, run as `julia -t 2 --project=verify verify/megno.jl` from the repo root; full output in `megno-verify-log.txt` (one line per measurement: `AC<n> | quantity | value | bound | verdict`). GPU: `verify/gpu.jl`, tag `megno` (`GPU_ONLY=megno julia --project=verify -t 4 verify/gpu.jl`), output in `gpu-megno-verify-log.txt`. Steps 1 and 3 (evidence part) of the stage contract; step 2 (independent review) is done by another agent, see "Independent review".

Run facts. Julia 1.13.2, 2 threads (GPU run 4 threads), base seed 20261012, verify offset 58000 (seeds `MG_SEED + 58000 + 50 ac + sub`; the package tests use `MG_SEED + 35000 .. 35999`, so these samples are fresh; the K = 0.5 regular class and the AC6 to AC10, AC14 panels use the package screening and seeds of `mg_case`, as the spec's panel). Working tree of branch `main`, nothing committed. Wall time of the final run 487 s, of which the package suite at full size 309.6 s (1 792 980 assertions, 0 fail, 0 error; the package-size suite has 1 000 962). The measurement helpers are those of `test/megno/mg_helpers.jl` (`mg_M5`, `mg_M6`, `mg_ref_s`, `mg_rows`, `mg_Psi`, `mg_BM`); the loops are written in the script, so a definition error in a helper would not be caught by the script; the independent pieces are the closed forms of M8 and M10, the 256-bit norm-form reference, the Abel form M6 against the weighted form M5, and the ChaosTools oracle (own code path for the tangent evolution). The @info lines of the suite (R1 to R6, the realized ratios R2) are forwarded into the log as `SUITE |` lines. The 8 eps in the text of two suite @info lines is stale: the tolerance in the code is 16 eps (Amendment D74 item 3), and the logged ratios are against 16 eps.

## Evidence table

"Where": T = package test (full size through the `suite` section), V = `verify/megno.jl`, G = `verify/gpu.jl`. Two entries per row are Float64 / Double64 unless stated. Bound 1 means the ratio error / bound.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| suite | package suite at full size (`MG_FULL_RUN`: every multiple of 10 as checkpoint, tau = 1 at n = 1e6 in AC2, Double64 at 1e5 in AC1(v), every n in AC5(b)), 26 testsets | 1 792 980 pass, 0 fail, 0 error | all pass | PASS | T via V |
| 1 | (i) a_m = 0.05, tau in {1,5}, both columns: M against the closed forms of M8(a) over 6 B_S | 0.0138 / 0.00045 | 1 | PASS | T, V |
| 1 | (ii) s(k) = ln(1 + k), j <= 200: Y against 2 (ln(j+1) - ln(j!)/j); Y, Ybar, Ybar20 against the BigFloat M5 sums | 0.00322 / 0.00012; 0.00323 / 0.00012 | 1 | PASS | T, V |
| 1 | (iii) alternating +0.2, -0.2: against the M5 sums; abs(Ybar(200)), abs(Ybar20(200)) over 0.1 | 0.033 (0.0033 absolute) | 1 | PASS | T, V |
| 1 | (iv) memory M7: Y_A(40) = 2.05, Y_B(40) = 1.05, equal s(40) (difference 1.3e-15 / 4.9e-32); against M5 sums | 0.0071 / 0.00037 | 1 | PASS | T, V |
| 1 | (v) long run a_m = 0.1, n = 1e5 and 4e6 (Float64), 1e5 (Double64): against M8(a) over 6 B_S; K^3 at 4e6 = 6.4e19 above 2^63, no integer weights | 0.0071 / 0.00042 | 1 | PASS | T, V |
| 1 | negative controls (s/j no memory, no weight, no factor 2, mean over 0..K-1, Int64 n^3 at 4e6 in the suite) wrongly satisfying the predicate | 0 of 16 per type (and the Int64 control in T) | 0 | PASS | T, V |
| 2 | K = 0 shear, seeded and W0 = I, tau in {1,10}, n to 1e6 (tau = 10): three series over 6 B_S; W0 = I column 1 exactly 0; column 2 closed forms | 0.0060 (worst over both types, suite R2) | 1 | PASS | T |
| 3 | linear hyperbolic map Y5, seeded ids (no exclusion, D74 item 2), supplied frames |c_1| = 1e-3, 1e-6 with max(6 B_S, 48 eps Psi_n), tau in {1,5}, n to 200: worst ratio | 0.0189 | 1 | PASS | T |
| 3 | (c) as amended (D74 item 1), |c_1| = 1e-3, n = 4: Y of column 1 and 2 (module and reference agree) | -3.067 and 4.812 (Double64 -3.067 and 4.812) | column 1 < 0, column 2 > 0.5 | PASS | T |
| 4 | hyperbolic fixed point K in {1,2,6}, ids 1..100, tau in {1,5}: three series over 6 B_S | 0.0258 | 1 | PASS | T |
| 5 | rotation (all three series exactly 0, tau in {1,3}); K = -2 elliptic closed forms every n in 1..1000 (tau = 1) and tau = 2: worst ratio over 6 B_S | exact 0; 0.0143 | 0; 1 | PASS | T |
| 6 | same-orbit reference, 20 starts, every multiple of 10 up to 1e4, tau = 1 / tau = 4, K = 0.5 / 2 pi / 10, Float64: error over B_M (three series, both columns) | 0.0064, 0.0165, 0.0260 / 0.0035, 0.0159, 0.0174 | 1 | PASS | T, V |
| 6 | same, Double64 | 0.0003, 0.0835, **0.545** / 0.0001, 0.0836, 0.543 | 1 | PASS (largest ratio of the (S) rows) | T, V |
| 6 | negative control, Float32 accumulator at n = 1e3, min error / B_M over 3 starts: K = 2 pi / 10 | 7.9e5 / 1.7e6 | > 1 | PASS | T, V |
| 6 | reported: the log sums against the borrowed 8 eps Psi_max, worst: Float64 0.026 to 0.100; Double64 0.001, 0.256, **1.637** (K = 0.5, 2 pi, 10) | see finding 1 | reported (D74 item 7) | INFO | T |
| 7 | (a) sums s_1, s_2 of MEGNO bitwise equal to the FLI, sali and L1 sums (20 starts x 1000 checkpoints per K) | 0 differences | 0 | PASS | T |
| 7 | (b) three series against the M6 post-processing of the stored L1 series, and (c) the Mestre relation (A2), with the final S_1 of the FLI reducer, remainder = 2 X_L1: error over 16 eps max(1, max abs s_1), K = 0.5 / 2 pi / 10 | 0.317 / 0.183 / 0.145 (both) | 1 | PASS | T, V (figure panel 3: worst 0.176) |
| 7 | (d) drawdown: D below -tol at 0 of 20000 points (min D / tol = -0.053); Y_1 = 2 (FLI_1 - D - FLIbar) + (2/n)(FLI_1 - D) | 0.317 | 1 | PASS | T, V |
| 7 | negative control (A2 with the running maximum): pairs with D > 10 tol (K = 0.5) all fail the relation; smallest error / tol | 19197 of 19197; 6.4e8 | > 1 | PASS | T, V |
| 7 | R1: Y_1 + 2 D_1 (min; median; max) at n = 1e3: K = 0.5 / 2 pi / 10 (never negative) | 2.23; 2.91; 5.25 / 1118; 1173; 1234 / 1578; 1628; 1700 | reported (>= 0) | INFO | T |
| 7 | R1 labels at n = 1e3 (Ybar20 > 0.5, Ybar > 2, n X_L1 > ln n, FLI_B > ln n): true counts and orbits where the four disagree, K = 0.5 / 2 pi / 10 | (12, 6, 2, 12), 17 disagree / (20 each), 0 / (20 each), 0 | reported | INFO | T |
| 8 | (a) K = 0.5 regular, 20 starts, n = 1e5: starts with Ybar_1 in [1.5,2.5] and Ybar20_1 in [0.4,0.6]; Ybar_1 min; median; max = 1.829; 1.999; 2.031 (all other starts within 1.98 to 2.03); Ybar20_1 = 0.4925; 0.5000; 0.5142 | 20 of 20 | >= 16 | PASS | T, V |
| 8 | (b) |mean(slope) - mean X(1e4)| over 4 SE + allowance, 1000 starts, K = 2 pi / 5 / 10: 2 m_Ybar | 0.074 / 0.075 / 0.023 | 1 | PASS | T, V |
| 8 | same, m_Y; 4 m_20 | 0.060 / 0.0003 / 0.070; 0.061 / 0.023 / 0.075 | 1 | PASS | T, V |
| 8 | realized differences of the means over the allowance alone (4/180, 8/180 for the (2,0) form): at most 0.159 | 0.0007 to 0.159 | reported | INFO | V |
| 8 | R3 regular class n = 1e4: Ybar_1 min; median; max 1.357; 1.994; 2.077, Ybar20_1 0.429; 0.5004; 0.581 (n = 1e5 above); chaotic means: mean X 1.168 / 0.965 / 1.62 against 2 mean m_Ybar 1.171 / 0.968 / 1.625 and the Chirikov ln(K/2) 0.916, 1.609 | in the log | reported (K18) | INFO | T, V |
| 8 | R4 first handoff with Ybar_1 >= 3 (median; never within 200): K = 2 pi / 5 / 10 | 2, 4 / 3, 10 / 1, 0 | reported | INFO | V |
| 8 | R5 random unit vectors (400): Ybar_1(1e4) min; median; max 0.039; 1.999; 2.149; first n within 0.25 of 2: median 40; seeded columns 1.357; 1.994; 2.077 | in the log | reported | INFO | T |
| 8 | R6 Mestre orbits K = 0.9 (n = 1e3, 2e3): chaotic Y_1/2 = 80.6, 197.4 against FLI - FLIbar = 80.44, 197.3 (remainder 0.162, 0.185); regular Y_1/2 = 0.9755, 1.026 against 0.9701, 1.023 (remainder 0.0054, 0.0031) | in the log | reported | INFO | T |
| 9 | tau in {2,5,10} against 1, K = 0.5 / 2 pi, n in {100,1e3}, 20 starts: (a) sums over 16 eps Psi_n, Float64 / Double64 | at most 0.0043 (K = 0.5), 0.0419 (K = 2 pi) / 0.00018, 0.0007 | 1 | PASS | T, V |
| 9 | (b) Y_tau - Y_1 over (tau - 1) L + B_M^tau + B_M^1; (c) Ybar over 3 (tau - 1) L + B_M terms (the B_M terms are negligible, so these are the ratios to the L-term): K = 0.5; K = 2 pi | (b) 0.122; 0.670 (c) 0.089; 0.325 | 1 | PASS | T, V |
| 10 | (a) Float64 against Double64, difference of means of Y_1/n and Ybar_1/n in combined SE: K = 2 pi, 5 (n = 1e3, 1000 starts); K = 5, 10 (n = 1e4, 200 starts) | 1.05; 0.51; 0.29; 0.92 | 4 | PASS | T, V |
| 10 | (b) BigFloat(256) on its own orbit against the 512-bit reference, K = 2 pi, 50 starts, n to 60, 300 rows | 0.0193 | 1 | PASS | T, V |
| 10 | (c) BigFloat against Double64, n = 300, 100 starts: K = 2 pi / 5 | 0.59 / 0.15 SE | 4 | PASS | T, V |
| 10 | (d) K = 10, n = 1e4, Float64, 20 starts: Y_1 = 1.607e4 to 1.637e4 (median 1.620e4), not clipped (differs from the reference by at most 3.6e-12); error over B_M | 0.0077 | 1; band 1.4e4..1.9e4 | PASS | T, V |
| 11 | (a) zero tangent (kind, orbit, step, column 1 or 2, value 0, lowest id first), (b) NaN tangent (column 0, threshold Inf), k = 1 ArgumentError, three chains: failing checks | 0 of 13 | 0 | PASS | T, V |
| 11 | (c) frame check: norms 2, 3, 0, 1 + 1000 eps rejected (4 of 4); unit, dependent, 1 + eps accepted (3 of 3); seeded ids 1..1e4 accepted, largest defect 2.5 eps; rejection before any map call | all | 8 eps | PASS | T, V |
| 11 | (d) diag(e^3, e^-3): tau = 200 no error (Y = Ybar = 2 s, Ybar20 = 3 s to relative 2.5e-16), tau = 300 error at handoff 300, both chains | 0 failures | 0 | PASS | T, V |
| 11 | (e) run record on {MEGNO}: indicator megno, convention tag, k 2, tau 2, precision Float64; final fields | true | true | PASS | T, V |
| 12 | (a) orbit points bitwise equal across the eight chains, 10 cases; (b) stand-alone outputs, final W and sums (SALI, FLI, L1 per column) bitwise equal, 10 cases | 0 of 10; 0 of 10 | 0 | PASS | T, V |
| 12 | (b) `propagate` equals the serial per-orbit body bitwise (x, W, finals), 4096 orbits, chunk 64 and 1000, 3 chains | true (10 of 10) | true | PASS | T, V |
| 12 | (c) outputs hashed in child processes with 1, 2 and 4 threads, chunk 64 and 1000: x, W, every final of {MEGNO}, {MEGNO(f), SALI}, {FLI(f), MEGNO(f), SALI} x tau in {1,4} (12 cases) | identical | identical | PASS | V |
| 12 | (d) bytes allocated by `propagate_orbit` after warm-up, n = 1e4, tau in {1,4}, two chains, both types | 0 (8 of 8) | 0 | PASS | T, V |
| 12 | (e) column exchange and sign change, K = 2 pi and 0.5, both types, tau in {1,4}: failures | 0 of 40; 0 of 80 | 0 | PASS | T, V |
| 13 | (a) K = 0 dyadic, tau in {1,3}: bitwise differences wrapped / lifted | 0 of 100 | 0 | PASS | T, V |
| 13 | (b) K = 0.5, n = 10 / 100, 50 starts, both columns: error over 6 * 64 K M_state n^2 eps rho_max; excluded (rho_max > 50) | 0.00023 / 0.00019; 3 / 6 of 100 | 1; <= 10 | PASS | T, V |
| 14 | ChaosTools 3.5.4 one-column tangent system (hand-coded Jacobian, `step!` loop, own renormalisation, BigFloat log norms, M5 in BigFloat): oracle orbit against `forward`, bitwise differences: K = 0.5 (20 x 1e4 x 2 columns), K = 2 pi and 5 (1000 x 1e3 x 2 columns) | 0, 0, 0 | 0 | PASS | V |
| 14 | |M_Chaoslab - M_oracle| over B_M + B_oracle (B_oracle = B_M: same orbit bitwise, the oracle's own tangents are ours): K = 0.5 (160 pairs); K = 2 pi (6000); K = 5 (6000); excluded 0 of 4000 each | 0.0048; 0.0126; 0.0142 | 1; <= 10 percent | PASS | V |
| 15 | (a) K = 0.5, Float64, 1024 orbits, n = 4e4, each backend against its own BigFloat(256) reference (final states reproduced bitwise 1024 of 1024 on both), Y; Ybar; Ybar20 over B_M, CPU / GPU | 1.36e-4 / 1.58e-4; 4.84e-5 / 4.91e-5; 3.33e-5 / 3.64e-5 | 1 | PASS | G |
| 15 | (a) chaotic K = 10, 1024 orbits, n = 8 / 10: Y, CPU / GPU; Ybar; Ybar20 | 0.031 / 0.032, 0.029 / 0.026; 0.018 / 0.017, 0.018 / 0.019; 0.015 / 0.011, 0.011 / 0.010 | 1 | PASS | G |
| 15 | (b) regular class, orbits with final states within 1e-9 (144 of 1024): gap over B_M,cpu + B_M,gpu, Y; Ybar; Ybar20 | 0.0033; 0.00077; 0.00048 | 1 | PASS | G |
| 15 | (b) chaotic n = 8 (all 1024 enter, largest horizon 1.06e-5) / n = 10 (989 enter, 35 beyond 1e-4, largest 5.2e-4, reported only): gap over B_M,cpu + B_M,gpu + 6 max_g D_orb, Y; Ybar; Ybar20 | 0.0045; 0.00096; 0.0010 / 0.0041; 0.00089; 0.00082 | 1 | PASS | G |
| 15 | (b) chaotic looseness (INFO), n = 8: Y gap without the D_orb term, max ratio 3.5e4 (614 of 1024 above 1); the term 6 max D_orb median 2.5e-7, max 5.9e-5; n = 10: 2.6e6 (685 of 989 above 1), term median 1.2e-5, max 0.0126 | | reported | INFO | G |
| 15 | (b) item 4 check on the device orbit, n = 8 / 10: max of abs(z_k^cpu - z_k^gpu)_inf / (16 eps M_k Gamma_k); against the literal 16 eps M_n Gamma_n (orbits above 1: 0 / 0) | 0.118, 0.049 / 0.097, 0.061 | 1 | PASS | G |
| 15 | (c) K = 10, 1024 orbits, n = 1e4, tau = 2, Y_1/n and Ybar_1/n: max per-orbit gap / sigma_FT; difference of means / SE | 3.94; 0.137 and 5.41; 0.249 | 6; 4 | PASS | G |
| 15 | reported: final states apart by more than 1e-9, CPU against GPU / CPU Float64 against Double64: regular n = 4e4; chaotic n = 8; n = 10 | 880 / 876; 74 / 96; 324 / 592 (of 1024) | reported | INFO | G |
| 15 | Double64 (regular class, chaotic class and the ensemble) | kernel does not compile for the element type | recorded fallback reason is the pass | FALLBACK (pass) | G |

Recorded lines of the spec: R1, R3 to R6 in the table above; R2 (realized ratios error / bound): at most 0.026 (Float64, criteria 2 to 6), 0.545 (Double64, criterion 6, K = 10), 0.317 (criterion 7 against 16 eps), 0.014 (criterion 14, against the doubled bound), 0.032 (criterion 15, Float64 on the device, chaotic class); the factor 48 of M11 holds with a margin of 1.8 or more everywhere.

## Oracle line (criterion 14), as the spec states it

ChaosTools 3.5.4 has no MEGNO routine (spec F7). The line is the oracle's independent tangent evolution: `TangentDynamicalSystem` with k = 1 and the hand-coded Jacobian of the lifted standard map, one `step!` per map step, the deviation read after each step, its log norm taken in BigFloat(256), then the deviation renormalised in Float64 (needed: the unnormalised column overflows Float64 for the chaotic panels, e^900), and the three series formed by M5 in BigFloat from those per-step log norms. Both seeded columns are run (Q0 = column i). First the oracle orbit equals the stream's lifted orbit bitwise at every step (0 differences in 400 000 + 4 000 000 steps), so no orbit divergence term is needed and B_oracle equals B_M (the oracle's tangents are ours). The largest difference over B_M + B_oracle is 0.0048 (K = 0.5), 0.0126 (K = 2 pi) and 0.0142 (K = 5); no pair is excluded. This tests the stretch series and the weighting against an independent code path, not the definition (no independent implementation for maps exists; the definition is tested by the closed forms, the two algebraic forms and the printed laws), as the spec says.

## D74 item 7: the Double64 K = 10 log sums

Measured on the 20 Double64 starts of the AC6 panel, K = 10, n to 1e4, both columns, with Psi_n on the full handoff set g = 1..n (not the grid) against the exact (BigFloat) total sum + comp:

| Quantity | Value |
|---|---|
| max abs(s - s_ref) / (8 eps Psi_n), own Psi_n | 1.637 |
| same with Psi_max over the full handoff set | 1.635 |
| same with Psi_max on the 10-per-decade grid (read at the grid points) | 1.629 |
| same with Psi_max on every multiple of 10 (the suite's R2 line) | 1.637 |
| three series, abs(M - M_ref) / (48 eps Psi_max), Psi_max over the full set | 0.545 |
| columns with the ratio above 1 | 2 of 40 (both columns of one start) |
| max(abs x, abs p) of the lifted orbit at their worst n | 9.85e5 |
| largest ratio over all columns at the n with max(abs x, abs p) below 1e3 | 0.0045 |
| MEGNO sums (sum and comp) against the L1 Benettin sums, same column, n = 1e3 and 1e4, 40 columns | bitwise equal (0 differences) |

Result. The excess is not the grid: the full handoff set gives 1.635, the grid 1.629. It is a real violation of the borrowed log-sum constant (8 eps Psi_n of L1 Amendment D35) by a factor 1.64, on one start of twenty (two columns of forty) at K = 10 in Double64; it is inherited from L1 (the MEGNO sums are bitwise the Benettin sums) and not from the weighted sums. The cause is the Double64 trigonometry floor on the lifted orbit (`PR-trig` finding of the precision spec Amendment D35, lyap-benettin verification section "Double64 trig floor"): the ratio is below 0.005 wherever the orbit coordinates stay below 1e3 and exceeds 1 only after the lifted coordinate reaches about 1e6. The three MEGNO series stay at 0.545 of B_M because the factor 6 of M11 absorbs it. Not a MEGNO matter; see finding 1.

## Findings that need a decision

1. D74 item 7 (above): the log sums of Double64 at K = 10 exceed the borrowed 8 eps Psi_n by 1.64 on a start whose lifted orbit reaches |x| near 1e6 (Double64 sin and cos floor). The criterion as written bounds only the three series (B_M = 48 eps Psi_max, passed at 0.545). Options for stage 05: keep the constant and record 0.545 and 1.64 (K18); or note in the spec that the log-sum bound needs the trig floor term of the precision spec for Double64 at large lifted arguments. This is an L1 matter shared with sali and FLI, not new to MEGNO.
2. Criterion 8(b): the realized differences are below the allowances (4/180, 8/180) by factors 6 to 1500 (largest 0.159 of the allowance, K = 5, 2 m_Ybar); the allowances could be halved (K18).
3. Criterion 7 tolerance: the realized ratio is 0.317 against 16 eps (D74 item 3), so 0.63 against the original 8 eps, as the test plan expected; 16 eps keeps a factor 3.
4. Criterion 3(c) as amended: Y(4) = -3.067 (column 1) and 4.812 (column 2) in the module and the reference; consistent with D74 item 1.
5. Criterion 15: Float64 runs on the device (`provenance.backend == :gpu` on every Float64 GPU line); the Double64 kernel does not compile (probe `compiles = false`), so the Double64 lines are the recorded fallback reason, the allowed outcome. The 200-byte state compiles and runs: the MEGNO kernel takes 4.61 ns per orbit-step on the GPU (20 000 orbits x 1000 steps, K = 2 pi, tau = 1, best of 3) against 22.6 ns on the CPU with 4 threads (speedup 4.9), compared with FLI (96 bytes) 3.33 ns (CPU 16.1 ns, speedup 4.8) and SALI (72 bytes) 4.27 ns (CPU 15.9 ns, speedup 3.7). So the extra state costs 38 percent over FLI on the device and 40 percent on the CPU (the 8 compensated additions and the division per handoff); no register-pressure failure at the extension's stack limit of 4096 B. The D_orb term for the chaotic class is loose by four orders (it carries the gate: without it 614 of 1024 orbits at n = 8 and 685 of 989 at n = 10 exceed B_M,cpu + B_M,gpu); the gate is passed with ratio 0.0045 and 0.0041. The chaotic class is gated at n = 8 (all 1024) and n = 10 (989; 35 beyond the item 4 horizon 1e-4 reported, gap ratio at most 0.0006 over the full bound), as for FLI; n = 4e4 is not possible for the chaotic class (the D_orb sum is quadratic and the orbits separate).
6. Criterion 12(c): the thread-invariance check (1, 2, 4 threads, child processes, chunk 64 and 1000, three chains) passes; it was a stage 04 item by D74 item 8.
7. Criterion 8 R1, finding F2: on the K = 0.5 regular class at n = 1e3 the four labels disagree for 17 of 20 orbits (counts above threshold 12, 6, 2, 12), and agree for all 40 chaotic ones; Y_1 + 2 D_1 is never negative (min 1.37 at n = 100), as M13(b) requires. Not a pass criterion.
8. Not measured here: coverage and size claims of the impl notes; the benchmark of `bench/megno.jl`; Experiments chain wiring (the one-replacer rule is the caller's, D72 item 1).

## Margins and marginal criteria

Largest ratios in the (S) rows with bound 1: 0.545 (criterion 6, Double64, K = 10), 0.317 (criterion 7 against 16 eps), 0.670 (criterion 9(b), K = 2 pi, to the L-term: the realized difference is about 0.6 of the allowed (tau - 1) L), 0.325 (criterion 9(c)), 0.262 (criterion 10(a), z = 1.05 of 4); all others below 0.15 (8(b) 0.075, 14 0.014). Criterion 9(b) is the only one with a ratio above 0.5 besides the Double64 line, and it measures a derived law, not rounding (K18).

## Deviations from the impl notes, checked

- Zero allocations in the per-orbit body: 0 bytes at n = 1e4, tau in {1,4}, both types, both chains.
- Bitwise equality of the sums with sali, FLI and the L1 sums per column (I5): 0 differences on 3 x 20 x 1000 checkpoints (Float64) and on 40 Double64 columns at n = 1e3 and 1e4.
- State size 200 bytes Float64 / 368 Double64 (as the notes), FLI 96, SALI 72.
- `check_frame` threshold 8 eps: seeded frames have defects at most 2.5 eps; columns of norm 1 + 1000 eps rejected.
- Isbits reducer and state: the Float64 kernel compiles and runs on the device (the notes had not measured it); Double64 does not compile.
- The notes' remark that `T(n)^3` is inexact past n about 2e5 in Float64: the three series at n = 1e5 (AC8(a)) and 1e6 (AC2) stay inside B_S and B_M; no growth of the ratio seen.
- Coverage and size claims were not re-measured here.

## Runtime

`verify/megno.jl`: 487 s total (package suite at full size 309.6 s, D74 item 7 86.8 s, criterion 14 34.9 s, the figures 15.6 s, criterion 10 15.2 s, criterion 12 15.2 s with the child processes, the rest under 3 s each). GPU run (`GPU_ONLY=megno`, 4 threads): several minutes wall, not timed (BigFloat(256) references and the 65-point Psi grid for 1024 orbits per backend at n = 4e4, D_orb for 2 x 1024 chaotic orbits). The machine was not otherwise idle, so timings are not benchmarks; the kernel timings above are best of 3 on the RTX 5060 Laptop GPU (compute capability 12.0, CUDA.jl 6.4.2, driver 13.4).

## Figures

- `megno-growth-laws.png`: left, eight regular K = 0.5 orbits, Ybar_1 against n up to 1e5 (log axis), approaching the level 2 (M3, M8(b)), with the slow approach and the transients of F3 visible; middle, ensemble means (200 starts) of Ybar_1(n) for K = 2 pi and 10 against lambda_1 n / 2 with lambda_1 the L1 mean X(1e4) of the identical orbits (1.167, 1.621); right, the Mestre Eq. (A2) residual of Y_1 in units of its tolerance 16 eps max(1, max abs s_1) at every multiple of 10 up to 1e4, three orbits for each of K = 0.5, 2 pi, 10, all below 0.18 (bound 1).
- `megno-chaostools.png`: left, Y of Chaoslab against Y of the ChaosTools one-column oracle line (M5 in BigFloat), K = 0.5, 2 pi, 5, on the diagonal; right, the ratio to B_M + B_oracle against n, all below 0.015 (bound 1).

## Files

- `verify/megno.jl` (new), `verify/gpu.jl` (new tag `megno`, section "AC15 (megno)" before the summary: functions `m15_ref`, `m15_one`, `m15_chaotic_b`, `megno_d15`, `megno_timing`, `megno_est`, `m15_grid`; `d16_orb` of the fli section got a keyword `matsf` so that it takes the megno BigFloat tangent; it reuses `d68_sig2`, `gpu_trace` and `ac12_chaotic` of the sali and shared sections).
- `megno-verify-log.txt`, `gpu-megno-verify-log.txt`, `megno-growth-laws.png`, `megno-chaostools.png`, this report.
- Nothing under `src/` or `test/` was edited; no tolerance was changed; nothing was committed.

## Independent review

Fresh Sonnet subagent, given the spec (with the checkpoint resolutions and Amendment D74), Mestre2011 (Eqs. 5, A1, A2) and `src/Indicators/megno.jl`, `Indicators.jl`, plus the Tangent protocol, Foundation, `fli.jl`, `sali.jl` and `benettin.jl` for comparison; not the notes, the test plan or the tests. Probes in the session scratchpad only. Result: no CRITICAL, HIGH or MEDIUM finding.

Confirmed by probes: the weighted sums follow Mestre Eq. 5 with the weight at the end of the window; Y, Ybar and Ybar20 agree with independent BigFloat(256) sums to 2e-16 relative on chaotic orbits (n up to 1e6, tau 1 and 10) and with the Abel and harmonic forms; Mestre (A2) holds at K = 50; the constant-increment closed forms hold at n = 4e6 with no Int64 overflow; the chain {MEGNO, FLI, SALI} (non-replacing, non-replacing, replacer) gives stand-alone results bitwise, the log sums are bitwise equal to SALI and to Benettin k = 1; the orbit and W are bitwise unchanged; health errors carry the right kind, column, orbit and step; `@inferred` for Float64, Double64, BigFloat; zero allocation; relaunch every 12 steps is bitwise equal to one launch; chunk sizes and threads give bitwise equal finals.

| # | Severity | Where | Finding | Disposition |
|---|----------|-------|---------|-------------|
| R1 | LOW | `megno.jl`, tau | tau is derived as n / K by integer division; only a direct caller with irregular steps could make Ybar20 wrong, the stream cannot | Add a one-line comment at stage 05 that n = K tau is a protocol invariant |
| R2 | LOW | `megno.jl`, provenance | The run record does not show the replace option | Same disposition as fli-ofli R1 and R2: the chain and its options are recorded by Experiments (D72); the docstring already warns |
| R3 | LOW | weighted sums | For bounded or regular orbits the products n a_k carry a rounding the compensated add does not repair: about n eps / 10 absolute on Y (1.0e-13 at n = 1e6 for alternating increments); chaotic orbits stay near 1e-17 relative | Within B_M in every run here (criterion 6 worst 0.545, criterion 9 below 0.67); optional exact-product correction not done; recorded as a limit |
| R4 | LOW | spec Compute | State is 200 bytes in Float64 and 368 in Double64, not 168 and 328; chunk sizing reads the real size | Spec note at stage 05 |
| R5 | LOW | `megno.jl`, add helper | The helper that adds to the five sums is not inlined; on the GPU MEGNO is about 38 percent slower than FLI per orbit-step (4.61 against 3.33 ns) | Performance only; optional `@inline`, not done |

The review also noted, without action: the reducer emits final values only, as FLI and SALI; `propagate` does not run the frame check (caller rule, D67); when both columns fail at one handoff only column 1 is reported, as FLI and SALI.
