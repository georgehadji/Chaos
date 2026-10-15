# Verification: lyap-benettin

Spec: `../../01-spec/output/lyap-benettin-spec.md` (criteria 1 to 19, Amendment D35). Package tests: `test/lyap-benettin/` (`lb_*.jl`). Oracle script: `verify/lyap-benettin.jl`, run as `julia +release -t 2 --project=verify verify/lyap-benettin.jl` from the repo root; full printed output in `lyap-benettin-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`). Figure: `lyap-benettin-convergence.png`.

Run facts. Julia 1.13.1, 2 threads, CUDA loaded (GTX 1050), base seed 20261012 with offsets `41000 + 50 * criterion + sub` (package tests use 21000 to 21999, so every sample is fresh). Wall time 4731 s (79 min) on a machine shared with other verify agents; timings are informative only. Library state: working tree of branch `m2-p1-modules`, no `src/` or `test/` edit. Values are "worst ratio of error to bound" unless stated; a criterion passes at ratio <= 1. Final run: 92 PASS lines, 0 FAIL, 0 ERROR, 29 INFO lines (Y1, PR-trig, some oracle details), 1 FALLBACK (Double64 on the GPU).

Method note (one-home rule). The measurement code is the package's: `lb_helpers.jl` and the driver functions `lb_ac5_data`, `lb_ac6_data`, `lb_regular` are included from `test/` through `include_patched`, which replaces only the seed and the sample size (the script fails if a patch target is missing). Criteria 1 to 4, 7, 8, 9, 10, 11, 17, 19 re-measure the package assertions as worst ratios with new drivers on the package helpers (`lb_run`, `lb_ref_std`, `lb_ref_full`, `lb_err`, `lb_BS`, `lb_Probe`). A definition error in a shared helper would not be caught by this; the independent pieces are the 256-bit references, the closed forms E10 to E13, and ChaosTools.

## Oracle: ChaosTools `lyapunovspectrum` (D27)

Installed: ChaosTools 3.5.4, DynamicalSystemsBase 3.19.3 (`~/.julia/packages/ChaosTools/0SbsT/`). The oracle is the tangent routine `lyapunovspectrum(tands, N; u0)` with a user `TangentDynamicalSystem` (hand-coded Jacobian `[[1 + K cos x, 1], [K cos x, 1]]`, state (x, p), the lifted rule without wrap), k = 1, `Delta t = 1`. It runs the QR of the deviation vector at every step and divides by N.

**Finding for the owner (spec P4 premise).** `lyapunovspectrum(tands, N; u0 = ...)` calls `reinit!(tands, u0)` without `Q0` (`lyapunovspectrum.jl`), and `reinit!` then restores `Q0 = default_deviations(tands)`, the first unit vector e1 (`DynamicalSystemsBase/src/derived_systems/tangent_space.jl`). A `Q0` given to the constructor is discarded: measured, a system built with `Q0 = (0.6, 0.8)` returns exactly the e1 result (5.5833e-5 at K = 0.5, n = 1000, equal to Chaoslab with W0 = e1 to 1e-11, differing from Chaoslab with (0.6, 0.8)). So ChaosTools cannot start from the seeded w0 through this routine. The tangent-core and spec statement that it can (D30, P4) is wrong for the installed version. Consequence used here: the pathwise lines give Chaoslab W0 = e1 (identical start vector on both sides), and the statistical lines (13 loose, 14, 15) compare Chaoslab with its seeded w0 against ChaosTools with e1, which the S5 rule and the 1e-4 bound allow.

**Orbit identity.** The ChaosTools lifted orbit (rule `p += K sin x; x += p`) is bitwise equal to Chaoslab `forward` over 10 starts x 1000 steps at K = 0.5, 2 pi, 5, 10 (10000 of 10000 each). So the oracle sees the identical orbit and the lifted-rule remark of the spec holds.

## Evidence table

"Where": T = package test, V = `verify/lyap-benettin.jl`. Values from `lyap-benettin-verify-log.txt`. Float64 / Double64 where two numbers are given.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1 | (a) and (b), Y5 map, tau {1,5,10}, n {20,100,1e4}, 100 ids (4 excluded, limit 10) | 0.01757 / 0.0008686 (both a and b) | 1 | PASS | T, V |
| 2a | rot4: sum, comp, lambda nonzero cases | 0 / 0 | 0 | PASS | T, V |
| 2b | K = -2: error vs reference; even n; odd n against ln((3+sqrt5)/2) + B_S | 0.08506, 0.08506, 0.9992 / 0.002905, 0.002901, 0.9992 | 1 | PASS (odd marginal, see below) | T, V |
| 3 | K = 0 shear, n up to 1e4 (20 ids) / n = 1e6 (3 ids; Double64 moved from the package suite) | 0.01931 / 0.001012; 0.001296 / 0.0006517 | 1 | PASS | T, V |
| 3 | strict bound X < ln(n+1)/n, max nX / ln(n+1) | 0.9995 / 0.9995 | 1 | PASS (tight by nature) | T, V |
| 3 | w0 = (0,1), n = 1e4 against 8 n eps | 0.000718 / 0.000638 | 1 | PASS | T, V |
| 4 | hyperbolic fixed point, K = 1, 2, 6 | 0.0158, 0.01605, 0.04331 / 0.000883, 0.002473, 0.002693 | 1 | PASS | T, V |
| 5 | same-orbit propagation, K = 0.5 (n 1e4), 2 pi (1e4), 10 (1e3), 20 starts x tau {1,4}; Double64 with 20 starts at every n (moved) | Float64 0.005105, 0.0003142, 7.557e-05; Double64 0.0002186, 1.874e-05, 4.772e-06 | 1 | PASS | T, V |
| 5 | runs with G_first > 1e4 (reported, not dropped): K = 2 pi 15 of 40, K = 10 17 of 40, K = 0.5 none (both types); orbit equals forward, all handoffs read | true | true | PASS | T, V |
| 6 | tau invariance, tau {2,5,10}, K {0.5, 2 pi}, n {100, 1e3}, 240 comparisons | 0.001694 / 6.951e-05 | 1 | PASS | T, V |
| 7 | full BigFloat orbit, K = 0.5 (S7), 50 starts, n = 10, 100 (2 starts with rho_n > 50 excluded at n = 100, limit 5) | n 10: 0.001049 / 2.449e-05; n 100: 0.001938 / 7.167e-06 | 1 | PASS | T, V |
| 7 | same without the M_state factor (superseded printed form, informational) | n 100: 0.2541 / 0.001016 | - | INFO | V |
| 8 | safe_norm relative error / 2 eps (5 scales, 2000 random vectors, |e| <= 900), zero, NaN, Inf | 0.5572 / 0.02056; true | 1 | PASS | T, V |
| 9 | health checks (zero tangent: kind, orbit, step, backend, column, value, nothing stored, lowest id first; NaN tangent; Double64 body; k != 1) | 36 of 36 | all | PASS | T, V |
| 10a | orbit bitwise equal for chains {}, {L1}, {L1, det} | 10 of 10 | all | PASS | T, V |
| 10b | `propagate` (2 threads, chunks 64 and 1000) vs serial `propagate_orbit`, x, W, finals, provenance; Float64 4096 orbits, Double64 1024 orbits (moved) | 4 of 4, 4 of 4 | all | PASS | T (Float64), V |
| 10b | 1 thread vs 4 threads, child processes, 2048 orbits, n = 1e3: hash of x, W, finals, Float64 / Double64 | equal / equal | bitwise | PASS | `verify/threads.jl` (`threads-verify-log.txt`) |
| 10c | bytes allocated per orbit after warm-up, n = 1e4, tau 1 and 4, both types | 0 | 0 | PASS | T, V |
| 11a | K = 0 wrapped vs lifted bitwise | 100 of 100, both types | all | PASS | T (Float64), V |
| 11b | K = 0.5 regular class, n = 10, 100, 50 starts | Float64 0.001209, 0.000934; Double64 7.35e-06, 4.95e-06 | 1 | PASS | T (Float64), V |
| 12a | Float64 vs Double64 lifted, 1000 starts, n = 1e4, K = 2 pi and 5: difference of means / 4 SE | 0.01658, 0.006964 | 1 | PASS | V |
| 12b | Float64 vs BigFloat(256), Double64 vs BigFloat(256), K = 2 pi, 100 starts, n = 1e3 | 0.1423, 0.2070 | 1 | PASS | V |
| 12c | regular class K = 0.5, 100 starts, n = 1e5, max abs(X_64 - X_D64) | 4.601e-08 | 1e-4 | PASS | V |
| 13 | regular class K = 0.5, 20 starts, n = 1e5, max abs(X - lambda_ChaosTools), seeded w0 vs e1 | 4.074e-05 | 1e-4 | PASS | V |
| 13 | same with e1 on both sides | 1.174e-14 | 1e-4 | PASS | V |
| 13 | pathwise, e1 on both sides, K = 0.5, n = 100, 50 starts, tau = 1, against 64 K M_state n^2 eps rho_n (2 excluded) | 3.334e-05 | 1 | PASS | V |
| 13 | pathwise, K = 2 pi, n = 1e3, 20 starts, in units of 8 eps Psi_n | 0.4715 | informational | INFO | V |
| 14 | K = 2 pi, 1000 starts, n = 1e4: difference of means / 4 SE (means 1.1686 vs 1.1687) | 0.00436 | 1 | PASS | V |
| 15 | K = 5 (200 starts): 0.000989; K = 10 (200 starts): 0.01306 (means 0.9537 vs 0.9538; 1.6186 vs 1.6186) | see left | 1 | PASS | V |
| 16 | CPU vs GPU, Float64, chaotic K = 10, 1024 orbits, n = 1e4, tau = 10: max per-orbit gap / sigma; difference of means / SE | 4.127; 1.046 | 6; 4 | PASS | V |
| 16 | regular K = 0.5, 1024 orbits, n = 4e4, tau = 10, against a Double64 reference: orbits above 1e-6 CPU / GPU; max CPU-GPU gap | 17 / 16; 3.66e-07 | count diff <= 3; gap <= 2e-6 | PASS | V |
| 16 | Double64 kernel | does not compile for the device | compile or reason | FALLBACK(kernel does not compile for the element type) | V |
| 17 | regular decay, 20 starts, slope of ln X vs ln n in [1e4, 1e6]: starts inside [-1.1, -0.8] (slopes -0.918 to -0.863, median -0.904) | 20 of 20 | >= 16 | PASS | T, V |
| 18 | alpha range (E14) over the criterion 5 runs, worst per case | Float64 0.9975, 0.99998, 0.9978; Double64 0.9975, 0.9999, 0.9995 | 1 | PASS (marginal by construction, see below) | T, V |
| 18 | alpha range over the criterion 6 runs | 0.9997 / 0.9999 | 1 | PASS | T, V |
| 19a | det residual / 4 eps (2 + K), K = 0.5, 2 pi, 10, tau {1, 4}, nonzero | at most 0.0836 (Float64), 0.0026 (Double64) | 1 | PASS | T, V |
| 19b | injected det fault: raised in [2 theta, 4 theta), quiet in (theta/4, theta/2] | 1 of 1 and 1 of 1, both types | all | PASS | T, V |

All 19 criteria have a number or a result line. Criteria 12 to 16 exist only in `verify/` (and the GPU part on the local device). The Float32 types are out (D24).

## Y1: large-K known answer

Catalog row Y1 (`_config/catalogs/systems.md`): "large-K lambda_1 vs Chirikov estimate ln(K/2)", sources Chirikov1979, Greene1979. It states no bound. The spec (Oracle, "Loose sanity lines") makes it a reported, never pass or fail line (K18) for X(1e6) ensemble means at K = pi/2, 2 pi, 10, against the benchmarks rows (`_config/bibliography/benchmarks.md`): ln(K/2) (Chirikov law, described by Shevchenko2004; the primary source is listed in `unverified.md`, so it counts as secondary-confirmed) and ln(K/2) + 1/K^2 for K >= 4.4 (Shevchenko2014 Eq. 26, "no error bound"). No tolerance exists, so none was invented; all lines are INFO. Measured (Float64, lifted, 100 uniform starts, n = 1e6, plus K = 20 as an extra large-K point):

| K | mean X (SE) | ln(K/2) | ln(K/2) + 1/K^2 | mean - ln(K/2) | mean - (ln(K/2) + 1/K^2) |
|---|---|---|---|---|---|
| pi/2 | 0.2161 (0.0136), median 0.2993 | not asymptotic | - | - | - |
| 2 pi | 1.1723 (0.00013) | 1.1447 | 1.1701 | 0.0276 | 0.0022 |
| 10 | 1.6201 (0.00010) | 1.6094 | 1.6194 | 0.0107 | 0.00067 (6.7 SE) |
| 20 | 2.3052 (0.00009) | 2.3026 | 2.3051 | 0.0026 | 0.00006 (0.7 SE) |

Chirikov ln(K/2) underestimates by 0.011 at K = 10 and 0.0026 at K = 20, with the sign and size of the 1/K^2 correction; the Shevchenko2014 form is within 7e-4 at K = 10 and 6e-5 at K = 20. At K = pi/2 the mean mixes regular and chaotic starts (median 0.2993 is the chaotic-sea value). Single orbit at the E9-converted start (x, p) = (1e-3, -1e-3), n = 1e7 (SalaArtuso2015 used 1e9): K = pi/2 X = 0.2964 against 0.298 (difference -0.0016); K = 2 pi X = 1.1725 against 1.172 (difference 0.0005). The figure (right panel) shows the K = 10 ensemble mean against n with both curves; ChaosTools gives the same curve.

## Double64 trig floor (precision spec Amendment D35, `PR-trig` lines)

Unit: `effective_eps(Double64)` = 2^-100 = 7.889e-31.

- Arguments with random low part, log-uniform in each decade of |x| from 1 to 1e8 (3000 sin and cos evaluations per decade against BigFloat 256 of the exact value): max absolute error 0.030 eps (decade 1 to 10), 0.004 eps for every decade from 10 to 1e8 (2^-105 at worst). No growth with |x|; max relative error 0.037 eps. On such arguments the floor is the rounding level of the type, far below 2^-100.
- Arguments taken from realized lifted Double64 orbits (K = 2 pi, 20 starts, n = 1e4, 28580 points, max |x| = 6.8e6): max absolute error of sin and cos by decade of |x| in units of eps: below 1e3: 0.02 or less; 1e3 to 1e4: 0.45 (0 of 1246 above 1); 1e4 to 1e5: 8.0 (19 of 5274 above 1); 1e5 to 1e6: 64 (62 of 13402); 1e6 to 1e7: 484 (37 of 8335). So on orbit arguments the floor does grow with |x|: abs error at most 1.2e-4 eps |x| over |x| >= 1e3. The decade maxima are only partially linear (heavy tail: under 1 percent of points exceed 1 eps).
- Scratch probe (not in the log; Julia session outside the script): the large errors belong to arguments whose low part is tiny relative to the ulp of the high part. For one orbit point at |x| = 2.1e6 with |lo|/ulp(hi) = 0.004 the error is 229 eps; the same hi with lo multiplied by 10, 100, 200, set to 0, or replaced by a random low part gives 0.001 to 0.004 eps. Orbit points hit this band with probability of order 0.5 percent, random arguments rarely do, which is why a sample of random arguments does not show the floor. The mechanism inside DoubleFloats was not investigated.
- Attribution on 6 starts at K = 2 pi, n = 1e4 lifted Double64: max |stored tangent entry - K cos x at 256 bits| = 1608 eps over the orbit; yet the error of n X is 0.039 of the old linear bound 8 eps (n + G_first) against the 256-bit reference and 0.0069 against the product of the stored tangents. This sample does not reproduce the 1.63 times excess of the Batch B test writers (6 different starts); the sample is small, and the amended B_S (with the 16 K n^2 eps term) covers both.
- Statement for the precision spec: a Double64 trig floor beyond |x| = 1e3 exists on realized orbit arguments (up to about 480 eps at |x| ~ 1e6, bounded by 1.2e-4 eps |x|), not on arbitrary arguments. Below |x| = 1e3 the 2^-100 convention holds (max 0.45 eps). The wrapped state (D28) keeps |x| < 2 pi and is not affected. Per the precision spec, a lifted Double64 bound beyond |x| = 1e3 stays report-only unless the owner adopts a value.

## Margins and marginal criteria

- Criteria 18, 2b (odd n), 3 (strict bound) show ratios of 0.9992 to 0.99998, which is not a loose margin but expected: the bounds are attained by construction. E13 bounds nX at odd n by ln((3 + sqrt 5)/2) = 0.9624, reached when w0 lies along the largest singular direction of M; the strict bound X < ln(n+1)/n of E11 is nearly attained by w0 along the shear; the E14 range [sigma^-tau, sigma^tau] is exact at cos x = 1 with the vector along the extreme singular direction, and over 4e5 handoffs per case the extreme of a chaotic orbit comes within 1e-5 to 3e-3 of it. All stay at or below 1; the log prints five significant figures (K = 2 pi, tau-1 runs 0.99998). A rounding slack of the form 16 eps tau Phi_j is included in the bound. A ratio of 1 or more here would be a defect; none occurred.
- The B_S-based criteria (1, 3, 4, 5, 6) have ratios of 0.0043 to 0.043 for Float64 (criterion 5 K = 0.5: 0.0051) and 5e-6 to 0.0029 for Double64. The Amendment D35 margin estimates (14 to 80 and more) are consistent. As the amendment itself says, for the standard map the n^2 term dominates from n = 20, so criteria 5 and 6 are envelope tests.
- Criterion 7 and 11b: the M_state form has worst 0.0019 (Float64) against 0.2541 for the superseded printed form (n = 100), in line with the amendment (0.018 and 9.95).
- Criterion 13 pathwise: ratio 3.3e-05 for K = 0.5 against the derived bound; for K = 2 pi the two estimators agree to 0.47 of 8 eps Psi_n, and per-start differences of nX at n = 1e4 are 1.5e-14 (AC14 line).
- Criterion 16: Float64 chaotic class max per-orbit gap is 4.13 sigma against 6 (0.69 of the bound); the mean difference is 1.05 SE against 4. It is the closest GPU number; the regular class shows 17 against 16 orbits above 1e-6 (difference 1, limit 3).
- Criterion 12c: the lifted Double64 and Float64 regular series differ by 4.6e-8 only, so the 1e-4 bound is loose by a factor of 2000 (as stated in the spec, finite-time basis ~10/n).
- Criterion 17: slopes -0.92 to -0.86, expected -1 + 1/(b + ln n) of the spec.

## Deviations from the implementation notes, checked

- `safe_norm` scaling by the maximum entry with a division (not a power of two): the 2 eps bound holds with ratio 0.557 (Float64) and 0.021 (Double64) over 2000 random vectors at scales 2^-900 to 2^900 and the five fixed scales; zero, NaN and Inf behave as stated (AC8).
- Zero or non-finite norm recorded with code 1, `value` the norm, `threshold` 0.0 (Inf for non-finite), `column` 1, kind "deviation matrix"; later steps ignored: AC9 36 of 36 checks (error fields, orbit and step, nothing stored, lowest failing id first, NaN with column 0 and threshold Inf).
- `reducer_final` at n = 0 and the reducer protocol: covered by the package testsets "API reducer protocol and provenance" and "API prefix invariance", which the impl notes report green; not re-measured here.
- Runtime of the package suite at the 3 min ceiling (impl notes) and the full-suite state: not re-run in this stage; the owner decision on the ceiling remains open (implementer's note). `lyap-benettin` is not in the D38 move list; the three Double64 items of test plan deviation 3 (AC3 n = 1e6, AC5 with 20 starts at every n, AC10b Double64) were run here and pass.
- GPU: "Double64 `log` compiling on the device is untested": it does not compile for the element type; the run falls back with the recorded reason (AC16).
- Benchmarks: `bench/lyap-benettin.jl` (5 entries) not run (BenchmarkTools is a verify-environment dependency only); timing is not part of the criteria.
- Not in the Benettin log: `verify/gpu.jl` has no Benettin line, so criterion 16 is measured in `verify/lyap-benettin.jl` with the same tolerances and the same pattern as gpu.jl AC12a and AC12c.

## Figure

`lyap-benettin-convergence.png` (three panels, seeds LB_SEED + 43000 + sub): (left) K = 2 pi, 200 starts, ensemble mean of X(n) for Chaoslab (seeded w0) and ChaosTools (e1) at n = 10 to 1e4 with +-4 SE bands; difference of means over 4 SE falls from 0.43 at n = 10 to 0.0015 at n = 1e4 (start-vector term, 1/n). (middle) K = 0.5 regular class, 20 starts, X(n) of both implementations against n (log-log) with ln n / n and the 1e-4 bound; max |X_Chaoslab - X_ChaosTools| is 0.0497 at n = 100, 4.1e-4 at 1e4 and 4.07e-5 at n = 1e5, below 1e-4 from about n = 4e4 on (the criterion is n = 1e5). (right) K = 10 (Y1), 200 starts, mean X(n) for both against ln(K/2) = 1.6094 and ln(K/2) + 1/K^2 = 1.6194; the means settle near 1.620 by n = 1e4.

## Review findings

Independent review: fresh Sonnet subagent given the spec and the code, not the implementation notes. No CRITICAL and no HIGH finding. Disposition: closed (evidence named), fix (proposed change), record (decision row or spec amendment), owner (checkpoint question).

| Severity | Finding | Location | Disposition |
|---|---|---|---|
| MEDIUM | AC10b thread count never varied | `test/lyap-benettin/lb_ac08_11.jl:115-122` | closed: `verify/threads.jl`, child processes `-t 1` and `-t 4`, outputs bitwise equal (`threads-verify-log.txt`), Float64 and Double64 |
| LOW | BenettinState is 48 B against 16 B in the spec Compute section | `src/Lyapunov/benettin.jl:14-20` | record: spec Compute text to the measured layout |
| LOW | Double64 samples of AC1, AC2b, AC3, AC5, AC19 reduced for runtime | `test/lyap-benettin` | record: reported here as reduced samples; full sets in `verify/lyap-benettin.jl` |
| LOW | AC18 probe uses the plain norm, not `safe_norm` | `test/lyap-benettin/lb_helpers.jl:175` | fix: test helper only |
| LOW | no `reducer_step` test at extreme scale (entries near 3e300) | `test/lyap-benettin` | fix: one test case |

## Notes

- Runtime 4731 s with 2 threads; section times in the log: AC16 1637 s (GPU launches on a display GPU shared with other jobs, plus the Double64 reference), AC12 565 s (Double64 lifted at large |x| is slow), AC3 400 s (Double64 n = 1e6), figure 722 s, Y1 290 s, AC10 276 s, AC5 to 18 268 s. Loaded machine; not benchmarks.
- `verify/Manifest.toml` shows as modified in `git status`; it was not edited by this stage's script (other verify agents ran concurrently).
- AC16 used the local GPU; the Double64 line is the recorded fallback, not a kernel pass.
