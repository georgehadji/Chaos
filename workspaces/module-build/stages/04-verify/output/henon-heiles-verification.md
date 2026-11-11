# Verification report: henon-heiles (catalog row F3; spec D116)

Stage 04 of `workspaces/module-build`, 2026-11-11. Script `verify/henon-heiles.jl` (environment `verify/`, groups FULL1 to FULL4, O1 to O6, AC16T, FIG), raw output `henon-heiles-verify-log.txt`, figures `henon-heiles-energy.png` (energy error against tau, orders 8 and 4) and `henon-heiles-portrait.png` (SALI classification on the section, E = 1/8). Code at commit 4c83095 plus the Amendment D117 fixes. Status: accepted at the stage 04 checkpoint, 2026-11-11 (owner: "Accept, commit").

## Evidence per criterion

FULL = the package test files at full sizes (`HH_FULL_RUN = true`): first run 6621 + 887 + 21 + 65 tests, 2 failures (AC8 det, AC12b); after Amendment D117 FULL2 888 and FULL3 26 pass, 0 fail.

| AC | Result | Evidence |
|----|--------|----------|
| 1 | PASS | 256 bits: max gradient 0, V error 2e-78, Hessian eigenvalues exact (bound 2^-200); the benchmarks.md points have gradient norm 2 |
| 2 | PASS | worst error / (eps max(1, entry)): 1.18 ForwardDiff, 2.6 complex step, 0.03 Double64 (bound 8); tangent step against AD of the one-step map 0.0087 of its bound |
| 3 | PASS | (a) to (d) in FULL1 |
| 4 | PASS | A(E) harmonic limit: defect/E converges to 0.4167; grid area fraction error 6e-6, 4e-6, 6e-6 (bound 1e-3); 10^5 starts and the infeasible status in FULL1 |
| 5 | PASS | 200 starts per energy, T = 1000 |
| 6 | PASS | tau = 1/20: order 8 1.4e-14 (bound 1e-8), order 4 2.4e-7 (bound 1e-4); observed orders 8.03 to 8.08 and 4.00 to 4.4 (plot) |
| 7 | PASS | Taylor reference (order 40, 2^-200, 256 bits), window t <= 60 (U1): crossing counts equal; max time and point error 1e-13, 1.4e-12, 4.6e-11 at E = 1/12, 1/8, 1/6 (bound 1e-8) |
| 8 | PASS after D117 (first run FAIL, det marginal) | D117 rerun (`henon-heiles-d117-verify-log.txt`, M_max with the leg Jacobian): det residual / bound 0.022, 0.042, 0.118; J error 1.6e-13 to 8.3e-13; FULL2 888 pass, 0 fail. First run: J against the 256-bit Taylor finite-difference Jacobian 1.6e-13 to 8.3e-13 (bound 1e-6): PASS. Det residual / bound: 0.097, 0.24, **1.004** at E = 1/12, 1/8, 1/6; FULL2 the same single case (2.2737e-13 against 2.2654e-13). The defect is about eps |J|^2 with max |J| = 23.4; the bound uses M = the one-step tangent norm, not the leg product norm |
| 9 | PASS | signed slope 0.80, absolute slope 1.01 (bound 1.2); mean det defect 4.0e-15 |
| 10 | PASS | 200 x 200, 30396 points: fraction 0.546 (55 +- 10 points); R SALI(2000) 1.107; C SALI(2000) 1.3e-15, 0 by 10^4; S reported (0.69, 8e-5, 0 at 2000, 5000, 10^4); threshold sensitivity 0.529 to 0.555. 500 x 500 not run |
| 11 | PASS | 100 x 100: order 8 0.0057, 0.5455, 0.9654; order 4 0.0054, 0.5459, 0.9661 (gaps > 5 points, table difference < 0.1 point) |
| 12 | PASS after D117 (first run FAIL (b)) | D117 rerun, (b) scored before t = 50: max error 3.0e-9 (E = 0.2), 4.9e-8 (E = 0.25), bound 1e-6, 200/200 references found; [50, 100) reported: 8.5e-5, 6.4e-4; FULL3 26 pass, 0 fail. First run: (a) 0 returns of 1993 escapes continued to 2 R_esc; (c), (d) PASS. (b): against the Taylor reference 7.6e-12 for t < 25, 4.9e-8 for t in [25, 50), above 1e-6 for 7 of 200 escapes, all t > 50 (worst 6.4e-4 at E = 0.25); tau/8 is no closer (chaotic amplification of round-off, as U1) |
| 13 | PASS | (a) to (c), (e) in FULL4 (two methods, T = 1e5); (d) T = 100: Taylor 1.2e-14 (C), 1.3e-10 (E = 1/6); ChaosTools (Tsit5 at 1e-14, `stop_at_Δt = true`) 2.7e-12, 1.9e-9 (bound 1e-6) |
| 14 | PASS | S0 tangent 4.0e-14 (bound 1e-9); O rotation 2.1e-15 (bound 7.8e-12); S0 spectrum (0.9965, 0.0014, -0.0014, -0.9965); O exponents 7e-16 (bound 1e-12) |
| 15 | PASS | all six slopes and rates in FULL4 |
| 16 | CPU PASS; GPU NOT RUN | 8 against 1 thread: states, frames, SALI bitwise, 0 label differences of 30396. GPU: `verify/gpu.jl` writes `data/gpu-benchmark.toml` (owner decision) |
| 17 | PASS | 0 bytes, Float64 and Double64 |

Size reductions (INFO in the log): O3a reference horizon 200 with window 60; O3b legs 2 to 4 of 10 starts per energy; O3c 100 escapes per energy; criterion 11 at 100 x 100 (FULL3). The ChaosTools oracle used Tsit5, not a high-order solver (none in `verify/`); its orbit matches the module to 1.5e-13 over t = 100.

Implementation-note deviations 1 to 7: each confirmed in the code; none contradicted by an oracle (5: dot-product sectors match the atan2 rule on 1993 escapes).

## Independent review (Sonnet, spec and code, not the notes)

No CRITICAL or HIGH.

| Id | Severity | Finding |
|----|----------|---------|
| M1 | MEDIUM | `section.jl` `_leg_jacobian` projects with the exact vector field, not the integrator's d/ds of the sub-step: det J - 1 is O(tau^p) (ForestRuth 8.3e-7, order 8 1.5e-14). AC8's round-off bound holds only for high order |
| M2 | MEDIUM | after a crossing the frame restarts on the sub-step path while the state continues on the grid: J differs from the derivative of the recorded map by O(tau^(p+1)); within the order-8 tolerance |
| L1 | LOW | `energy_start` returns `:ok` with p1 = 0 on the exact boundary (r == 0); AC4 needs p1 > 0 |
| L2 | LOW | `escape_sector` tie at q1 = 0, q2 < 0 goes to label 1 for both signs of zero: reflection equivariance fails on that measure-zero set |
| L3 | LOW | a one-step jump from r <= R to overflow gives `:nonfinite` with no sector (only with huge R_esc or tau) |
| L4 | LOW | grazing test uses `eps(T)`, crossing tolerance `4 effective_eps(T)` (units differ for Double64); code 4 eps against the spec's 16 eps_T |
| L5 | LOW | `tau` not checked to be a float: Int is type-unstable, Rational bisects in exact arithmetic |
| L6 | LOW | `energy_grid` has no `Real` fallback for E (MethodError on 1//8) |
| L7 | LOW | `random_energy_start` fixed 10^4 tries fails for a tiny D(E) in a large window (documented) |
| L8 | LOW | `section_interval` upper root conditioned like sqrt(eps) near E = 1/6 (inherent) |
| L9 | LOW | `max_crossings = 0` returns `:buffer_full` at once |
| L10 | LOW | with `jacobian`, the state is advanced twice per step (about 2x cost) |

Not checked by the review: the tau convergence of crossing and escape times (its command was blocked by a hook); the O3 comparisons are covered by the oracle runs above.

## Audit

| Check | Result |
|-------|--------|
| Evidence | every criterion has a number; criterion 16 GPU part not run (owner decision) |
| Review closed | no CRITICAL or HIGH finding |
| Criteria passing | 17 of 17 after Amendment D117 (first run: AC8 det bound and AC12b failed on bounds) |

## Amendment D117 (owner: "Amend + fix L1,L2,L4")

Spec section "Amendment D117", test plan section of the same name, decision log D117. AC8 M_max includes the leg Jacobian (bound for the order-8 table only, M1 reported); AC12b scored before t = 50. Code fixes L1 (`energy_start` rejects p1 = 0), L2 (`escape_sector` 1-2 tie by the sign bit of q1), L4 (grazing test in `effective_eps`), with the new testset "D117 boundary start and the 1-2 label tie". Package henon-heiles suite after D117: 6904 pass, 0 fail. Reruns of FULL2, FULL3, O3b, O3c: all PASS (`henon-heiles-d117-verify-log.txt`); the other groups are unaffected by the changes (FULL1, FULL4, O1, O2, O3a, O3d, O4, O5, O6, AC16T, FIG not rerun). Open: M2, L3, L5 to L10 (LOW or within the order-8 tolerance).
