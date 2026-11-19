# Verification report: continuous-qr (catalog row L8; spec D120)

Stage 04 of `workspaces/module-build`, 2026-11-18. Script `verify/continuous-qr.jl` (environment `verify/`, groups FULL1 to FULL4, FULL2h, Oa, Oid, Od20, Od20R, Oc9, Oord, AC11d, Od21, AC17T, AC17alloc, AC22, R1, FIG), raw output `continuous-qr-verify-log.txt`, figures `continuous-qr-convergence.png` (exponent error against h, orders 2 and 4: fitted slopes 4.0 and 2.4 on the closed form, 4.2 and 1.94 on orbit C) and `continuous-qr-orthogonality.png` (E(t) with and without the penalty, forms 84 and 92). Code at commit 8850000 (Amendment D121 added 2026-11-19). Status: accepted at the stage 04 checkpoint after Amendment D121, 2026-11-19 (owner: "Accept").

## Evidence per criterion

FULL = the package test files at full sizes (`CQ_FULL_RUN = true`): FULL1 681 pass; FULL2 84 pass, 4 fail (AC6b 2, AC8, AC9); FULL2h (FULL2 at the calibration horizon t = 60 to 100) 86 pass, 2 fail (AC8, AC9); FULL3 67539 pass, 1 fail (AC11); FULL4 96 pass. No size reduction.

| AC | Result | Evidence |
|----|--------|----------|
| 1 | PASS | BigFloat(256): Q = I bitwise at every step, l_i within 16 eps max(1, abs(A_ii) t); FULL1 |
| 2 | PASS | BigFloat: abs(l_i) <= 8 eps t; Q(t) against exp(At) Q0 0.17 of B_lin (worst, order 2) |
| 3 | PASS | degenerate pair individually: error / B_lin 0.044, 0.0016, 0.0052 |
| 4 | FAIL, first ratio only (BigFloat ladder); Float64 FULL1 PASS | fixture b, q = 4: ratios 38.5, 16.4, 15.5, 15.7 against [8, 32], slope 4.25; the first ratio (h = 0.1 to 0.05, Delta norm(A) = 0.70) is pre-asymptotic; without that pair PASS |
| 5 | PASS | error / bound 0.0075 (l_1), 0.075 (Q) |
| 6 | FAIL at full horizon; PASS at t = 60 | AC6b column independence at t = 100: 0.0137, 0.0288 against 5e-3; at t = 60 worst 1.3e-3 |
| 7 | PASS | flow-direction identity on Taylor BigFloat(256) references to about 1e-63; module against the Taylor orbit within the plan bounds (worst about 5e-5, pendulum amplitude 2, q = 4) |
| 8 | FAIL | energy-covector identity, module side: FULL2 (Tmax 300) worst 14.1, 6946 violations; FULL2h (Tmax 100, 83 orbits) worst 1.47 against 3e-2 (stage 02: 1.55e-2 on 11 orbits). The identity itself holds to 1e-63 on the reference; the module against the Taylor orbit for R and C, t <= 60, passes 3e-2. Cause: discretization error amplified along chaotic orbits beyond the pointwise horizon (as henon-heiles U1); 2 starts per class at calibration |
| 9 | FAIL | T = 60, 20 starts per class: lambda error / allowance 2.32 (stage 02: 0.47). T = 1000, lambda form, spec c_d = 5: realized c_d 13.7 (C), 13.5 (R), 299 (chaotic class), 254 (regular class); l form, plan 5e4 h^4: 5.97 (chaotic), 5.08 (regular), R and C pass. Discretization order (orbit C): q = 4 ratios 28.3, 18.9, 16.7, slope 4.36; q = 2 slope 1.99. Spread among discrete variants about 1e-15 (R6) |
| 10 | PASS | AC10a slope -0.929 in [-1.1, -0.8]; AC10b maxima about 1e-3 against 5.1e-3 |
| 11 | FAIL, one assertion | 200 chaotic starts, T = 2000: all lambda_1 > 0; 191 of 200 within 5e4 h^4 / T of the discrete lambda_1 (max 10.7 times, median 0.075); zero pair, pairing, ensemble agreement PASS |
| 12 | PASS | worst abs(sum) / B_sum 4.5e-5 |
| 13 | PASS | worst E / bound 1.0e-3; E_eq ratio h to h/2 4.19 (q = 2), 16.8 (q = 4) |
| 14 | PASS | t = 1e3, 20 starts per class: error / B_C 7.0e-3 (Float64), 1.1e-4 (Double64) |
| 15 | PASS | 51 orbits in agreement, BigFloat error / bound 1.6e-4 |
| 16 | PASS | FULL4 |
| 17 | PASS | 1 against 8 threads bitwise, 4096 starts, t = 1e3, both orders (orbit-level threading: `propagate_cqr` is serial); chunk 64 against 1000 bitwise; 0 bytes (Float64, Double64; k = 1, 2, 4; both orders) |
| 18, 19 | PASS | FULL4 |
| 20 | PASS (plan bounds) | Taylor BigFloat(256) variational reference: F3 C, t = 60: q = 4 5.4e-4 (spec 1e-3), q = 2 0.117 (spec 1e-1 exceeded, plan 1.0 passes); pendulum amplitude 3, q = 4: 0.0129 at t = 100 (plan 0.1); q = 2 skipped (U7, 4.1 reported); orbit R to t = 1000, q = 4: 0.0953 (h = 1/20), 0.0040 (h = 1/40). The discrete RGS spectrum on the same reference: about 1e-12 |
| 21 | PASS | ChaosTools 3.5.4, 100 starts each, t = 500: lambda_1 0.05014 against 0.04787 (difference 0.0023, allowance 0.0195); lambda_2 0.006172 against 0.006355 (1.8e-4, allowance 0.0094) |
| 22 | device NOT RUN; CPU part PASS | FULL4; `:gpu` and `:auto` resolve to the CPU with the reason "reducer is CPU-only" (Float64 and Double64): the module has no kernel. Device: `verify/gpu.jl` writes `data/gpu-benchmark.toml` (owner decision) |

Main finding: the continuous spectrum's error is set by its own Runge-Kutta order (q = 4 with step 2h), not by the integrator. With the order-8 integrator the discrete QR spectrum agrees with the Taylor reference to about 1e-12, while the continuous one is off by 5e-4 at t = 60 and 0.095 at t = 1000 (orbit R, h = 1/20). Criteria 8, 9, 11 measure this gap against bounds calibrated on 2 starts per class at t = 60.

Recorded lines:
- R1: continuous / discrete RGS time, orbit C, N = 20000: Float64 0.36 (q = 2 and 4), 0.34 against Householder; Double64 0.76 (q = 2), 0.69 (q = 4). Geist1990's 40 is not reproduced.
- R2: E_eq, form 84, h = 1/20: 8.5e-5, 1.0e-3, 9.5e-3 for c_nu = 1, 2, 4 (log Od20).
- R3: error / B_C 7.0e-3 (Float64), 1.1e-4 (Double64); error / B_lin at most 0.17.
- R4: orbit C, t = 60: order 2 0.117 against order 4 5.4e-4 at h = 1/20.
- R5: form 92, k = 2, no penalty: E reaches 1e-2 at t = 173.5 (NaN at step 3878); with the penalty never (Skokos2010 p.53 confirmed).
- R6: continuous against discrete about 1e-3 to 1e-1 at T = 1000; among discrete variants about 1e-15.

Implementation-note deviations: confirmed (chunk ignored, CPU only with the "reducer is CPU-only" reason, E_max default 1e-2, AC12 bounded not rounding, equilibrium check test-side, order 2 poor on the pendulum at amplitude 3, R5). Refined: the bench figure of the notes (1.69 ms) was for 2000 steps; the T = 200 failure of the lambda form extends to T = 1000, and the l form is exceeded there too.

## Independent review (Sonnet, spec and code, not the notes)

No CRITICAL or HIGH. Checked correct: Eq. 85 skew part, the form 92 term (spec C10), the rates of Eqs. 86, 87 (penalty only in dQ/dt), Heun and RK4 nodes and weights on the orbit points (observed orders 4.0 and 2.0, capped by the integrator's order), compensated sums, per-unit-time scaling, 0 bytes, guards for k, N, form, equilibrium.

| Id | Severity | Finding |
|----|----------|---------|
| M1 | MEDIUM | `propagate_cqr` does not check the state length against the flow dimension: a 3-vector on the 2D pendulum runs and returns 2 exponents |
| M2 | MEDIUM | the rates use (Q^T A Q)_ii without dividing by q_i . q_i: up to about 1 percent relative bias passes the E_max = 1e-2 guard silently (not measured) |
| M3 | MEDIUM | `:gpu`/`:auto` with Float64 and Double64 record "reducer is CPU-only"; the spec Compute section excludes that reason |
| M4 | MEDIUM | one orbit that throws aborts the whole ensemble, nothing returned |
| L1 | LOW | no validation of `penalty_factor`, `orthogonality_tol`, map type: late, confusing errors |
| L2 | LOW | Q0 not checked against the guard or against `k` (wrong `k` in the record) |
| L3 | LOW | a non-finite state after the start is labelled "deviation matrix", and is missed when A does not depend on it |
| L4 | LOW | non-bits fields (Symbol, Union) block a device kernel |
| L5 | LOW | the flow record is a printed string, lossy for BigFloat parameters |
| L6 | LOW | the spec still says E_max 1e-3 |

Not checked by the review: the size of the M2 bias, the ensemble behaviour of the other Lyapunov modules (M4), implicit integrators.

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number; criterion 22 device part not run (owner decision) |
| Review closed | no CRITICAL or HIGH finding |
| Criteria passing | 16 of 22 clean; 6 with a failing line: 4 (pre-asymptotic first ratio), 6 (horizon), 8, 9, 11 (discretization bounds calibrated on 2 starts per class), none traced to a defect of the equations; M1 to M4 open |

## Amendment D121 (owner: "Amend + fix M1-M4", then "Keep M2, amend")

Spec section "Amendment D121", test plan section of the same name, decision log D121. Criteria restated from the measured numbers: 4 (ratio gate on pairs with Delta ||A||_F <= 0.35, window [2^q / 4, 2 x 2^q]), 6(b) (t <= 60, B_col 5e-3), 8 (0.4 + 8 eps / c at t <= 60), 9 (c_d 20 for q = 4 and 10 for q = 2 at T = 60; T = 1000 reported), 11 (at least 95 percent of starts), 12 (order of the sum reported; the bound stays scored). Code fixes M1 (dimension and frame-shape guard), M2 (rates (Q^T A Q)_ii / (q_i . q_i); bias on 20 random frames with column defects 2.5e-3 to 7.5e-3: 1.66e-2 before, 4.4e-16 after), M3 and M4 (docstrings), with the new testset "D121 dimension guard, normalised rates, orbit-named errors".

M2 trade-off (owner chose to keep M2): the old rates made the trace cancel exactly, so the sum of criterion 12 was a pure discretization term; with the normalised rates the sum carries an O(E) term. Orbit R, q = 4: sum 1.6e-7 before M2, 1.3e-5 after; worst |sum| / B_sum 2.0e-4. Sum ratio h = 1/20 to 1/40 (reported): R 4.23 (q = 2), 11.7 (q = 4); C 0.817 (q = 2), 14.3 (q = 4). Criterion 4, Oa fixture a, q = 2: ratios 5.71, 1.93, 3.22, 3.65, slope 1.67 (other fixtures and orders: slopes 4.07, 2.02, 4.0).

Reruns (`continuous-qr-d121-verify-log.txt`): package continuous-qr suite 1625 pass, 0 fail, before the "Keep M2, amend" lines (not rerun after them; they only widen one window and drop one assertion); Aqua, JET, layering 14 pass; FULL1 680 pass; FULL2 88 pass; FULL3 67536 pass, 0 fail; FULL4 107 pass; Oa no FAIL (criterion 4 all gated pairs in the window); AC11d PASS; Oc9 INFO. Oid, Od20, Od20R, Od21, AC17T, AC17alloc, R1, FIG not rerun after D121.

Criteria passing after D121: 22 of 22, criterion 22 device part not run (owner decision). Open: L1 to L6 (LOW).
