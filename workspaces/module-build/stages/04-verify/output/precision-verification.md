# Verification: precision

Spec: `workspaces/module-build/stages/01-spec/output/precision-spec.md` (with Amendments: AC6 5 percent allowance, Double64 effective eps 2^-100; D24: Float32 dropped).
Environment: `verify/` (Julia 1.13.1, 4 threads), script `verify/precision.jl`, log `precision-verify-log.txt`, commit 1361767 (branch m2-p1-modules, working tree with the D28 edits). Reference helpers are the independent ones of `test/precision/helpers.jl`. Seeds: StableRNG from base 20261012 with the same offsets as the package suite. Sample sizes: 200 orbits x 100 steps (AC4, 5, 6), 10^4 points per function (AC3), 10^4 values (AC2), 100 ICs (AC13), 4000 orbits x 500 steps per element type (AC11).

**Timings in AC11 were taken on a loaded machine** (several other Julia processes were running); only ratios are meaningful and they vary by tens of percent between runs (first run: Double64/Float64 57.6 at 1 thread, 56.1 at 4 threads; BigFloat/Float64 364 and 576; the table gives the final run).

## Evidence table

| AC | Quantity | Value | Bound | Verdict | Where |
|----|----------|-------|-------|---------|-------|
| 1 | max err / bound, K in {pi/2, 2pi, 0.1, 1e-3}: Float64 / Double64 / BigFloat(256) | 0.8 / 0.031 / 0.113 | 1 | PASS | suite + verify |
| 2 | round-trip failures, 10^4 values, Double64 / BigFloat(256) | 0 / 0 | 0 | PASS | suite + verify |
| 3 | sin max err / (eps_D max(1,\|f\|)), eps_D = 2^-104 | 0.524 | 16 | PASS | suite + verify |
| 3 | cos, same | 0.662 | 16 | PASS | suite + verify |
| 3 | log, same | 0.755 | 16 | PASS | suite + verify |
| 3 | effective eps (amendment: <= 2^-100): sin / cos / log | 2^-104.93 / 2^-104.60 / 2^-104.41 | 2^-100 | PASS | verify |
| 4 | Float64 n_H min / median / max | 14 / 19 / 100 (cap, 7 orbits) | floor 15 | PASS (7 orbits below floor, at most 10 allowed) | suite + verify |
| 4 | Float64 max \|e_n\| / (eps M_n Gamma_n), n <= n_H | 0.906 | 8 | PASS | suite + verify |
| 4 | Double64 n_H min / median / max | 39 / 55 / 100 (cap, 9 orbits) | floor 40 | PASS (1 orbit below floor) | suite + verify |
| 4 | Double64 max ratio | 0.403 | 8 | PASS | suite + verify |
| 5 | BigFloat 256 vs 512, max ratio, n <= 100 | 0.619 | 8 | PASS | suite + verify |
| 5 | max \|e_10\| | 8.2e-71 | 1e-60 | PASS | suite + verify |
| 6 | horizon min(n_H(Double64), 40) min / median / max | 40 / 40 / 40 | - | INFO | verify |
| 6 | max ratio \|S_n(D64) - S_n(ref)\| / (64 (1+K) eps M_n Gamma_n) | 2.605 | 1 | one orbit over | verify |
| 6 | median / 99th percentile ratio | 3.8e-4 / 0.055 | - | INFO | verify |
| 6 | orbits with ratio > 1 | 1 of 200 | <= 10 (amendment) | PASS | suite + verify |
| 7 | fallback-matrix mismatches (21 cases, Float32 rows dropped) | 0 | 0 | PASS | suite + verify |
| 8 | precision-string failures | 0 | 0 | PASS | suite + verify |
| 9 | GPU vs CPU orbit agreement | not run here | 2 c eps M_n Gamma_n | STAGE04-GPU | verify/gpu.jl |
| 10 | Float32 GPU gate | Float32 dropped from M2 | - | DROPPED D24 | - |
| 11 | Double64/Float64 time ratio, 1 thread / 4 threads | 64.2 / 46.0 | reported | reported | verify |
| 11 | BigFloat(256)/Float64 time ratio, 1 thread / 4 threads | 566 / 547 | reported | reported | verify |
| 11 | 1 to 4 thread speedup Float64 / Double64 / BigFloat | 3.09 / 4.31 / 3.20 | - | INFO | verify |
| 11 | GPU Double64/Float64 | not run here | - | STAGE04-GPU | verify/gpu.jl |
| 12 | closure violations; max allocated bytes in hot loops (Float64, Double64) | 0; 0 | 0 | PASS | suite + verify |
| 13 | median n* (min / max) | 16 (12 / 200) | reported | INFO | suite + verify |
| 13 | fraction of ICs with n* <= 60 | 0.99 | >= 0.8 | PASS | suite + verify |

Counts (log lines): 26 PASS, 0 FAIL, 3 resolved by a companion line ("see count" for the AC6 max ratio, "see exclusions" for the two AC4 horizon floors), 2 STAGE04-GPU (AC9, AC11 GPU part), 1 DROPPED (AC10), the rest INFO.
AC11 ensemble: 4000 orbits x 500 steps (2e6 steps), min of 3 timings; Float64 took 365 ms (1 thread) and 118 ms (4 threads) in the final run. This is a smaller ensemble than the 10^5 x 10^3 of AC10, because BigFloat(256) would take hours on the full one.

## Figures

`precision-orbit-error.png`: left, error against BigFloat(256) versus n for one orbit (index 33, median Float64 horizon), Float64 and Double64, with the bound 8 eps M_n Gamma_n (dashed) and the realised horizons n_H (dotted verticals) and the 1e-4 level; right, ratio error/(eps M_n Gamma_n) for all 200 orbits up to n_H, with the line at c = 8.

## Review findings

Code review (two independent reviewers, stage 04): no CRITICAL or HIGH findings; the MEDIUM findings were fixed in commit 1361767. Stage 02 rework (AC12 allocation helper specialised on the function argument for Julia 1.10) is recorded in the test plan.

Findings on the evidence:
- AC6: one orbit of 200 exceeds the first-order bound (ratio 2.6); the amendment allows 5 percent. The orbit is accelerator-like; the bound has no local-stretch term. No change.
- AC4: n_H is capped at the reference length 100 for 7 (Float64) and 9 (Double64) regular orbits, so their ratio is checked only to n = 100. Accepted: these orbits are regular and the cap is conservative for the horizon count.
- AC13: one initial condition of 100 never separated by 1e-8 within 200 steps; the 0.99 fraction passes.
- AC11: cost ratios (Double64 about 50x, BigFloat(256) about 550x Float64) were taken on a loaded machine and vary by run; treat them as orders of magnitude.
- Effective Double64 eps measured 2^-104.4 to 2^-104.9, inside the amended 2^-100 with a 2^4 margin.

## Notes

- Bound slack: the realised ratios are well below c = 8 (Float64 0.91, Double64 0.40, BigFloat 256 vs 512 0.62), so the constant 8 has a factor of about 9 to 20 of slack on these orbits. The test cannot pass vacuously: the error curve in the figure follows the bound closely from n = 1 to n_H.
- Effective Double64 eps: sin, cos and log all measure between 2^-104.4 and 2^-104.9 against BigFloat(256), at or below the nominal 2^-104, so the amendment value 2^-100 holds with a margin of about 2^4. Nothing goes to the owner.
- AC4 floors: Float64 has 7 orbits with n_H below 15 (minimum 14), within the 10 allowed (5 percent); Double64 has 1 orbit at 39 (floor 40). n_H is capped at the reference length 100 for 7 (Float64) and 9 (Double64) orbits (regular orbits with small Gamma); for these the ratio was checked only up to n = 100, not up to the true horizon.
- AC6 marginal: one orbit exceeds the first-order bound, with ratio 2.605 (the same orbit class the amendment describes; BigFloat 256 vs 512 shows ratio about 10 on one orbit). 199 of 200 are below 1 and the 99th percentile is 0.055, so the bound has no local-stretch term and is loose for typical orbits but violated by an accelerator-like one. Within the amendment (<= 10 allowed).
- AC13: one IC out of 100 never separated by 1e-8 within the 200-step scan limit (max n* = 200); the 0.99 fraction at n* <= 60 passes with margin. The median n* = 16 agrees with the derived 15 to 30.
- AC11: Double64 costs about 46 to 64 times Float64 per step on CPU and BigFloat(256) about 550 times, on a loaded machine. Float64 at 150 to 180 ns per step is itself several times slower than an unloaded core would give (a sin plus a few additions), which compresses or inflates ratios depending on how the load fell on each run; the first run gave Double64/Float64 of 57.6 and 56.1 and BigFloat/Float64 of 364 and 576. Treat the ratios as an order of magnitude. BigFloat scaling to 4 threads (3.2x) is allocation-bound.
- AC1 uses Float32-free type set (D24); AC7 omits the four Float32 rows. AC9 and AC11 GPU parts await `verify/gpu.jl`.
- The package tests for AC1 to AC8, 12, 13 use the same helpers and seeds, so verify and suite measure the same orbits; the verify run adds 10^4 points for AC3 (suite: 3 x 10^3 per function) and reports the realised values.
