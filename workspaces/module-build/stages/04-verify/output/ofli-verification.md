# Verification report: ofli (catalog row I6; spec D118)

Stage 04 of `workspaces/module-build`, 2026-11-11 (Amendment D119 added 2026-11-18). Script `verify/ofli.jl` (environment `verify/`, groups FULL1 to FULL4, Oa, Ob, Oc, PLATEAUp, PLATEAUf, Od, R, AC16T, AC18r, FIG), raw output `ofli-verify-log.txt`, figures `ofli-plateau.png` (OFLI_B against t on the pendulum and the period-1 orbit, order 4 and 8, h and h/2, with t*) and `ofli-classes.png` (OFLI_B against t for orbits R, C and the period-1 orbit). Code at commit 7eaf9d4. Status: accepted at the stage 04 checkpoint after Amendment D119, 2026-11-18 (owner: "Accept").

## Evidence per criterion

FULL = the package test files at full sizes (`OF_FULL_RUN = true`): FULL1 78 pass, FULL2 839 pass and 2 fail, FULL3 376 pass, FULL4 21 pass and 1 fail.

| AC | Result | Evidence |
|----|--------|----------|
| 1 | FAIL, marginal (d = 2 line); d = 4 PASS | 256 bits, d = 4: worst error / (16 eps max(1, abs g)) 0.36 (tau 1), 0.89 (tau 10). d = 2, tau = 1, N = 1e4: 1.25 (20 eps against 16 eps), growing about as n^0.3 from 4 eps at n = 100 (accumulated rounding of 1e4 renormalisations; 16 eps is a borrowed constant). FULL1 (package sizes, Float64 and Double64) 78 of 78 |
| 2 | PASS | d = 2, tau = 10: 0.5; frame refusal for w_p = 0 in FULL1; the d = 2, tau = 1 line is the one under criterion 1 |
| 3 | PASS | (a) identity on the Taylor reference 1.7e-62 (bound 1e-30); (b) order ratios in FULL1 |
| 4 | (a), (b) PASS; (c) FAIL at amplitude 3 | (a) 1.2e-61, 3.1e-61, including starts at q = 0. (c) start q = 0 (not a turning point), order-4 table, h = 1/20, t <= 10 P: amplitude 2 7.4e-7, amplitude 3 **3.5e-5** against delta 1e-5. The departure peaks early and falls to 8.8e-6 at the end of the window: a handoff-phase offset, not plateau loss. The delta 1e-5 came from the package start, a turning point, which gives 0 error (impl notes) |
| 5 | PASS | period-1 orbit at 256 bits: closure 1.6e-77, T = 6.0756, tr A = -0.9035, det A - 1 = -7.7e-48, T'(E) = -1.347. Exact series, canonical and seeded frames: dG 0.0019, 0.0009 (bound 0.1); dFLI 2.15, 2.22 (at least 1.5). Order 8, h = 1/20: dG 0.0009, dFLI 2.17. Order 4 reported with t* (R1) |
| 6 | PASS | slope for R 1.0022 in [0.7, 1.3]; 19 of 20 regular starts inside (outliers 0.563, 1.09, reported) |
| 7 | PASS | mean slope 0.04925 against mean exponent 0.04976 on the identical orbits; OFLI_B <= FLI_B, 0 violations |
| 8 | PASS | largest spread over frames 0.35 against ln 2 = 0.693 |
| 9 | PASS | BigFloat singular values (Jacobi): 0 violations in 800 rows (R, C) and 1600 rows per class; (F* - OFLI_B) / ln 2 in [0.001, 0.73] |
| 10 | (a), (c), (d) PASS; (b) 2 FAIL (test tolerance) | (a) sums bitwise equal to FLI and SALI, 294 assertions. (b) one start (:reg, 7), both frames: 2.22 and 1.85 against 1. The test bound 4 eps max(1, abs S) has no 1/sin theta term; min sin theta there is 0.024, and the spec's own rounding term 8 eps / sin theta (O12) covers it |
| 11 | PASS | (a) FULL3 0 fail. (c) angle ratios h against h/2: 4.06 (order 2), 16.02 (order 4); order 8 rounding-dominated (0.67 against 256). Candidate B on C: angle 0.1 at t = 831, ratio B/A 4.5e-5 at the end (collapse, F3) |
| 12 | PASS | FULL3 max error / B_G 0.139 over 17136 rows, 0 maximiser differences; Ob 0.042 (R, C), 0.10 (regular class), 0.046 (chaotic class) |
| 13 | PASS | OFLI_1 - OFLI_tau: 4e-4 and 7.7e-3 for R, about 1e-15 for C, against 2 B_G |
| 14 | PASS | (a) to (g) in FULL3 |
| 15 | PASS | BigFloat against its reference 1.05e-3 of B_G; overflow start: G = 3789 finite, error / B_G 7.3e-4 |
| 16 | PASS | FULL4 chunk sizes 64 and 1000 against serial bitwise, 0 bytes per orbit; AC16T 1 against 8 threads, 4096 starts, tau in {1, 20}: final records, states and frames bitwise |
| 17 | PASS | order 8 against Taylor 256 bits: R 8.7e-15, 1.2e-14, 9.3e-15; C 5.3e-15, 3.8e-14, 6.6e-13 (t = 10, 30, 60; bound 1e-6). Order 4 (reported): R 6.8e-6 to 1.9e-5, C 4.7e-6 to 5.5e-4. Orbit R to t = 1000: order 8 9.5e-13, order 4 5.4e-5 |
| 18 | device NOT RUN; CPU part 1 FAIL (gap rule) | each backend against its own reference: 0 rows above B_G (rerun of the same 200 + 200 starts). Gap rule fails on regular start 120 only: abs(OFLI_Float64 - OFLI_Double64) 4.83e-10, 5.19e-10 against B_G,Float64 3.5e-10, 3.3e-10; final states agree within the 1e-9 gate, but the rule has no orbit-divergence term (D68 gives one only to the chaotic class). Device: `verify/gpu.jl` writes `data/gpu-benchmark.toml` (owner decision) |

Recorded lines:
- R1, order-4 plateau-loss time t* (departure above 0.05), ratio t*(h/2) / t*(h) against 2^4 = 16: pendulum amplitude 2: 368, 5737, 9.09e4 at h = 1/2, 1/4, 1/8 (15.6, 15.8); amplitude 3: 4, 2411, 3.77e4, 6.01e5 (603 pre-asymptotic, 15.65, 15.93); F3 period-1 orbit: 2816 at h = 1/4, 4.13e4 at h = 1/8 (14.7), above 2e5 at h = 1/16 and 1/20. Order 8: no t* below the horizon from h = 1/4 (pendulum, 1e6) and h = 1/2 (F3, 2e5); no ratio measurable. The h^-p prediction of O11 holds for order 4.
- R2, FLI_B - OFLI_B (min, q25, median, q75, max): regular class t = 1e3 1.6e-4, 0.088, 0.41, 1.16, 2.77; t = 1e4 1.6e-5, 0.107, 0.416, 1.25, 3.9. Chaotic class t = 1e3 6.5e-7, 0.076, 0.179, 0.357, 1.10; t = 1e4 8.4e-7, 0.050, 0.159, 0.322, 0.88. No negative value. Orbit C at t = 1000: OFLI_B / t 0.0313 against the L1 exponent 0.0317.
- R3: criteria 12 and 15.
- R4, period-1 orbit (order 8 and 4 alike): running supremum 1.265 to 1.267 between 10 T and 100 T; the instantaneous value oscillates in [-0.044, 1.267] and is more than 0.01 below the supremum at 516 of 547 handoffs (spec F2: the sources' "oscillates around a constant" describes the instantaneous series).
- R5: all 200 chaotic starts reach ln 1e20 = 46.05 before t = 1e4; min 400, median 810.5, max 4502.

Size reductions (INFO in the log): Ob at N = 1e4 steps for R and C and N = 2000 steps for 20 starts per class (tau = 10); PLATEAUp horizon 1e6 and PLATEAUf 2e5 time units against Double64 order-10 references; AC18r a 24-start rerun (6 pass). Not repeated: the 10a bitwise sums in Ob (FULL2 covers them), the Taylor 4b slope to 1000 P (FULL1 covers it). AC16T compares final records, not every handoff.

Implementation-note deviations 1 to 4: consistent with the runs. The note on AC4c (0 error at the turning-point start) is confirmed, and the real value is criterion 4c above.

## Independent review (Sonnet, spec and code, not the notes)

No CRITICAL or HIGH.

| Id | Severity | Finding |
|----|----------|---------|
| M1 | MEDIUM | `ofli_initial` (`ofli.jl:134`) has no field health check: a zero or non-finite f(x0) returns -Inf in every column, not the "vector field" error of Algorithm item 2 and criterion 14a |
| M2 | MEDIUM | the reducer's `ofli`, `ofli1` exclude the s = 0 term (U1), and no function returns the combined value; `reducer_final.ofli` read alone can be below the floor of criterion 10d. Spec Algorithm item 1 puts the term inside G_i |
| M3 | MEDIUM | provenance `ofli_flow = string(nameof(typeof(flow)))` drops parameters: two oscillators with different omega record the same name, so P10's mismatch check cannot see them |
| L1 | LOW | cancellation in w - (u . w) u: relative error of the term about eps / sin theta, as O12; below sin theta about 1e-14 the term is noise, and an exactly parallel column gives about S + ln eps instead of -Inf |
| L2 | LOW | `check_frame` for k = 1 uses unscaled `dot(f, f)`: a zero field gives NaN reported as "initial frame"; overflow above about 1e154 and underflow below 1e-162 |
| L3 | LOW | `reducer_step` does not check length(x) against d |
| L4 | LOW | a subnormal column norm is accepted (GALI rejects it, D108j); only with `replace = false` |

Not checked by the review: the device kernel compile, a propagate-level bitwise run against FLI and SALI (FULL2 10a covers it), the 8 eps constant for d above 4.

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number; criterion 18 device part not run (owner decision) |
| Review closed | no CRITICAL or HIGH finding |
| Criteria passing | 14 of 18 clean; 4 with a failing line: 1 (d = 2 marginal, 1.25), 4c (3.5e-5 against 1e-5), 10b (test bound without 1/sin theta), 18 CPU gap rule (one regular orbit); none traced to a defect in `ofli.jl` |

## Amendment D119 (owner: "Amend + fix M1-M3")

Spec section "Amendment D119", test plan section of the same name, decision log D119. Criterion 1 bound 32 eps to n = 1e4 (Oa d = 2 rescored: 20 eps, PASS); criterion 4(c) at the start q = 0 with delta 1e-4 (PLATEAUp rescored: 3.5e-5, 7.4e-7, PASS); criterion 10(b) bound plus 8 eps / sin theta at the maximising column; criterion 18(b) regular-class gap plus the 1e-9 state gate. Code fixes M1 (`ofli_initial` field health check), M2 (`ofli_total`), M3 (flow record `sprint(show, flow)`), with the new testset "D119 field health of ofli_initial, the combined series, the flow record". Reruns (`ofli-d119-verify-log.txt`): package ofli suite 349 pass, 0 fail; Aqua, JET, layering 14 pass; FULL1 78 pass, FULL2 841 pass, FULL3 380 pass, FULL4 22 pass, 0 fail (criterion 18 CPU part PASS under the D119 gap rule). Criteria passing after D119: 18 of 18, criterion 18 device part not run (owner decision). Oa and PLATEAUp not rerun (bounds only). Open: L1 to L4 (LOW).
