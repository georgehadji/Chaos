# Verification: lyap-svd-stabilized

Spec: `../../01-spec/output/lyap-svd-stabilized-spec.md` (criteria 1 to 12, Amendment D35). Package tests: `test/lyap-svd-stabilized/`. Oracle script: `verify/lyap-svd-stabilized.jl`, run as `julia +release -t 2 --project=verify verify/lyap-svd-stabilized.jl` from the repo root; full output in `lyap-svd-stabilized-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`). Scope of this report: stage 04 steps 1 and 3 (oracle cross-check, figure, evidence). Step 2 (independent review) is filled in below by the orchestrator.

Run facts. Julia 1.13.1, 2 threads, ChaosTools 3.5.4, DynamicalSystemsBase 3.19.3, base seed 20261012, offsets 50000 + 100 x criterion (package tests use 30100 to 31200, GS 20100 to 21400, standard-map verify 40000+), so every sample is fresh. Effective eps: Float64 2.22e-16, Double64 2^-100 = 7.889e-31 (the value of the precision spec; all Double64 theta_sum use it). Machine loaded by concurrent agents; timings are not benchmarks. The measurement loops are rewritten in the script (the package criteria are `@testset`s); the helpers (`gs_*`, `sv_*`) are included from the test folders, so a helper error would not be caught here. The independent pieces are the 256-bit references (built on each run's own element-type orbit, K18), the closed forms and ChaosTools. Values are worst ratios of error to bound (pass at <= 1) unless stated. The log has 88 PASS lines, 0 FAIL, 0 ERROR (INFO and STAGE04-GPU lines are not verdicts).

Sample sizes beyond the package suite: V6 4 values of c x 20 draws for both types (suite: Double64 one c); AC4 at full spec size for both types (20 orbits, 4 K, tau 1,2,4, n 12,20,100,1000; suite Double64: 4 orbits, 2 K, 2 tau); AC7 standard map in both types (suite: Float64 only); AC3, AC5, AC8, AC9, AC12 in Double64 (stage 04 parts, Amendment D35 Coverage).

## Evidence table

T = package test, V = `verify/lyap-svd-stabilized.jl`. "F" = Float64, "D" = Double64.

| Criterion | Quantity | Value (F; D) | Bound | Verdict | Where |
|---|---|---|---|---|---|
| API V6 | closed form vs 256-bit 2 x 2 singular values, ratio to 16 eps (1 + \|S_1\| + \|S_2\| + ln(1+c^2)) | F 0.0279; D 3.869e-4 | 1 | PASS | T, V |
| API V6 | S_1 = 1e5 finite, ordered, sum within 16 eps 1e5; sigma_1 = sigma_2 finite | 0; true | 1; true | PASS | T, V |
| 1 | cat map, \|lambda_1 - ln mu\|, \|lambda_2 + ln mu\|, QR offset formula G12a, frame/tau pairs (ratios to theta_sum) | F 0.00205, 0.00501, 0.00389, 0.00177; D 1.04e-4, 8.03e-5, 1.22e-4, 8.27e-5 | 1 | PASS | T, V |
| 2 | K = 0: \|lambda_1 - ln((n+sqrt(n^2+4))/2)/n\|, \|lambda_1 + lambda_2\|, SVD-QR offset (frame (e_2, e_1)) | F 0.01115, 0.00432, 0.00315; D 1.23e-4, 5.40e-4, 6.91e-5 | 1 | PASS | T, V |
| 2 | lambda_1(1e4) - ln(1e4)/1e4 | F 1e-12; D 1e-12 | <= 2e-12 | PASS | T, V |
| 3(a) | QR-family output bitwise equal to `QRSpectrum(RGS())`, 600 orbits per type | 0 mismatches (F, D) | 0 | PASS | T (F), V (F, D) |
| 3(b) | V7 sandwich, lower / upper side ratio to theta_sum; ordering violations | F 9.27e-4 / 1.64e-3; D 3.32e-5 / 2.62e-5; 0 | 1; 0 | PASS | T (F), V |
| 3(c) | S_1(500) >= 25 (411 orbits): \|c_1000 - c_500\| / (8 eps max(1, \|c\|)) | 0 (F and D), 0 failing | 1 | PASS | T (F), V |
| 3(c) | 5 <= S_1(500) < 25 (114 orbits): \|dc\| / (0.06 max(1, \|c_500\|)) | 0.1255 (F and D) | 1 | PASS | T (F), V |
| 3(c) | S_1(500) < 5 (75 orbits): worst relative change, reported | 0.1374 (largest in class 0 <= S_1 < 2) | not scored | INFO | V |
| 4 | \|lambda_i - lambda_i^ref\| / theta_sum, 256-bit product at the own orbit, full size; per K (0.5, 2, 5, 2 pi) | F 0.0366 (0.0366, 0.0169, 0.0127, 0.0083); D 8.76e-4 (8.76e-4, 5.33e-4, 2.37e-4, 3.33e-4) | 1 | PASS | T, V |
| 5 | plain Float64 product + `svdvals` at n = 20: \|ln sigma_2^naive - ln sigma_2^ref\| (naive -13.88, ref -25.13) | 11.25 | > 10 | PASS (marginal) | T, V |
| 5 | first non-finite entry of the plain product (F; D) | 606; 589 | < 700 | PASS | T (F), V |
| 5 | estimator finite, ratio to theta_sum at n = 20 and 2000 | F 1.12e-4; D 1.94e-4 | 1 | PASS | T (F), V |
| 6 | rotation 0.7311, n = 1e4, tau 1 and 10: \|lambda_i\| / theta_sum, all finite | F 0.00203, 0.00526; D 2.97e-4, 6.09e-4 | 1 | PASS | T, V |
| 7 | cat map n = 1e5: \|sum\| / theta_sum; \|sum - (S_1+S_2)/n\| / 8 eps | F 0, 0; D 1.10e-4, 3.31e-5 | 1 | PASS | T, V |
| 7 | standard map, K 4 values, n = 1e4, 10 orbits: worst of the same two ratios, eps_pair == eps_sum | F at most 6.1e-4, 0.0226; D at most 1.06e-4, 0.0023; 0 violations | 1; 0 | PASS | T (F), V (F, D) |
| 8 | wrapped map, K {2,5} x tau {1,3}, n = 1200, 80 orbits: record, x, W bitwise vs plain loop | 0 mismatches (F, D) | 0 | PASS | T (F), V |
| 8 | 10^6 steps, K = 5, 8 orbits: finite, \|lambda_1 + lambda_2\| / theta_sum (M_max 6) | F 5.02e-8; D 2.52e-5 | 1 | PASS | T (F), V |
| 8 | bytes allocated after warm-up (also rechecked at end of run) | 0, 0 (F); 0, 0 (D) | 0 | PASS | T (F), V |
| 9 | 4096 orbits, K = 2 pi, chunk 64, 1000, default: x, W, finals (S_i, compensation, c) bitwise vs serial | true x 3 (F), true x 3 (D) | tolerance 0 | PASS | T (F), V |
| 9 | 1 thread vs 4 threads, child processes, 2048 orbits, n = 1e3: hash of x, W, finals, Float64 / Double64 | equal / equal | bitwise | PASS | `verify/threads.jl` |
| 9 | bytes per orbit body after warm-up | 0 (F, D) | 0 | PASS | T (F), V |
| 10 | regular K = 0.5 (S7: 1024 of 3072 regular, 17 replaced), n = 4e4, tau = 10, against a Double64 reference: orbits above 1e-6 CPU / GPU; max CPU-GPU gap, lambda_1, lambda_2 | 1 / 1; 8.98e-08 | count diff <= 3; gap <= 2e-6 | PASS | `verify/gpu.jl` |
| 10 | chaotic K = 10, Float64, 1024 orbits, n = 1e4, tau = 2: max per-orbit gap / sigma_FT; difference of means / SE, lambda_1, lambda_2 | 4.989; 0.6084 | 6; 4 | PASS | `verify/gpu.jl` |
| 10 | c check (spec Amendment D41, ensemble 1-ulp bound): max \|c_cpu - c_gpu\| / max(1, \|c_cpu\|) over the 1024 of 1024 orbits with S_1(500) >= 25 (orbit 441, \|c\| 1.97; median 1.1e-16; absolute max 6.1e-12, max \|c\| 167.9), against the max over the same orbits of the CPU change under a 1-ulp shift of x0[1] | 1.393e-12 | 3.909e-11 | PASS | `verify/gpu.jl` |
| 10 | conditioning of c: CPU change under a 1-ulp shift of x0[1], relative: max; median; at orbit 441 | 3.909e-11; 3.331e-16; 1.315e-12 | reported | INFO | `verify/gpu.jl` |
| 10 | Double64 kernel | does not compile for the device | compile or reason | FALLBACK(kernel does not compile for the element type) | `verify/gpu.jl` |
| 11 | K = 0.5 regular (20 orbits, 0 replaced), n = 1e5: \|lambda_1^SVD - oracle\| / (1e-4 + ln(1+c_n^2)/(2n)); first orbit | 0.1432; 0.127 (SVD 1.143e-4, oracle 9.975e-5, c_n 4.163) | 1 | PASS | V |
| 11 | K = 2 pi, 1000 orbits, n = 1e4: \|difference of ensemble means\| / (4 combined s.e.) (means 1.17, s.e. 0.002407) | 0.001685 | 1 | PASS | V |
| 12 | scope guard: exact message from `propagate` and `reducer_init`, (d,k) in {(4,4),(4,2),(4,1),(2,1)}, before any work | 0 of 8 deviations (F, D) | 0 | PASS | T (F), V |
| 12 | provenance record (estimator, orthonormalizer, d, tol, k, tau, seed, precision label, backend, device, reason; chained monitor) | 0 of 4 mismatches (F, D) | 0 | PASS | T (F), V |
| Figure | exponent error and V7 sandwich | `lyap-svd-stabilized-exponent-error.png` | - | see Figures | V |

## Oracle comparison

- 256-bit recomputation (Oracle section, package-suite oracle repeated at stage 04): AC4 worst ratios above; the reference is built from the 256-bit tangent product at the run's own orbit points with det Y = 1 used analytically. AC5 uses the same reference at n = 20 and 2000. Float64 error in lambda reaches 1.35e5 eps_eff at small n for K = 2 pi (theta_sum is far larger there because G_max enters; the scratch value "8 eps at n = 2000" is a different orbit and not reproduced as such: the ratio to theta_sum is at most 0.037).
- Analytic cases V9: AC1 (cat map), AC2 (K = 0 shear), both types, ratios 1e-4 to 1e-2 of theta_sum.
- ChaosTools `lyapunovspectrum` (AC11), lifted rule `p += K sin x; x += p`, hand-coded Jacobian, `k = 2`, `Delta t = 1`. Facts read from the installed source: `lyapunovspectrum(tands, N; u0)` calls `reinit!(tands, u0)`, whose default `Q0` is the identity columns, so a `Q0` given at construction is overwritten (checked by running with the identity and with a swap frame: results equal). FINDING for the owner: the lyap-benettin Oracle sentence that ChaosTools "accepts the identical initial vectors" holds only through `reinit!(tands, u; Q0 = ...)`, which `lyapunovspectrum` does not expose; the comparison therefore runs both sides with W0 = identity (the identical-frame rule S2 is met with that frame). The pathwise lambda_1^QR against the oracle differs by at most 2.4e-14 (20 regular orbits, n = 1e5) and 1.3e-14 (1000 chaotic orbits, n = 1e4), so the orbit and the QR exponents agree to rounding; the SVD minus oracle difference (1.7e-5 regular, at most 4.1e-4 chaotic, mean 1.6e-5) is the V7 offset ln(1 + c_n^2)/(2n) up to u/(2n). ChaosTools offers no singular-value route; the SVD claim rests on the 256-bit and analytic oracles.
- Published lambda(K) values: not used (K18, spec).

## Margins and marginal criteria

- AC5 negative control: gap 11.25 against the bound 10 (ratio 0.89 of the margin). The control orbit is the first seeded start with ln sigma_1(20) = 25.13, just above the Amendment threshold 25. The D35 60-start table is reproduced: ln sigma_1(20) from 14.22 to 30.96, median 24.5; the gap exceeds 10 (or the plain singular value is exactly 0) for 41 of 60; every orbit with ln sigma_1 >= 24.9 has a gap above 10 (true); smallest ln sigma_1 with gap > 10: 18.87 (exact zero singular value case); largest ln sigma_1 with a finite gap <= 10: 24.44; 29 of 60 starts qualify (D35: 21 of 60, a different seed). The criterion is satisfied but only because the orbit is selected by the Amendment rule; an unselected orbit fails it 19 times out of 60.
- AC5 Double64: the negative control is defined on the Float64 product. A plain Double64 product keeps sigma_2 at n = 20 (gap 0 to display precision, INFO line), so no gap is required; the Double64 part is the estimator side (ratio 1.94e-4) and the overflow (first non-finite entry at n = 589, same exponent range). Interpretation recorded here.
- AC2: lambda_1^SVD(1e4) - ln(1e4)/1e4 = 1.0e-12 against 2e-12 (derived margin 1.0e-12): ratio 0.5.
- AC3(c): class 5 <= S_1 < 25 worst ratio 0.1255 against the amended 0.06 (D35 worst 0.0261 relative, i.e. 0.43 of the constant; here 0.00753 relative = 0.1255 of 0.06): margin 8. Class S_1 >= 25: exactly 0 change in 411 orbits for both types. A concern raised before the run was that the threshold 25 comes from a Float64 estimate (increment below 1e-16) and could fail in Double64 (increment 2^-100 scale needs S_1 above about 40). It did not occur: a scratch re-run of the 200 AC3 starts at tau = 1 (not in the log) found the smallest S_1(500) in the class at 176.6 and none in [25, 100); the increments are about e^(-2 S_1), far below 2^-100. The threshold is conservative, not wrong.
- AC3(c) class counts (same for both types, 600 orbits): S_1 < 0: 0; [0, 2): 3; [2, 5): 72; [5, 25): 114; >= 25: 411. D35 counts were 744 of 1800 below 25; here 189 of 600.
- AC7 K = 2 pi, Float64: sum identity 0.0226 of 8 eps, the largest ratio of the criterion; all others below 0.003.
- AC4 Float64: 0.0366 of theta_sum (K = 0.5), the largest of the measured ratios outside AC1/AC2; no criterion is within a factor of 10 of its bound except AC5 above.
- Two non-claims of the script: the D35 per-orbit bound ratio (largest 0.104 in D35) is derivation only and is not recomputed; H1 to H3 are held and not scored.

## Deviations from impl notes, checked

- `reducer_final` field order equals the contract (n, sum, comp, lambda_qr, lambda, c, u, G_max, orth_max, count_max): confirmed. `lambda` is zero at n = 0 on a fresh state: confirmed.
- "c is not updated when the QR step fails": after two good steps (c = 0.6176) a rank-deficient window gives code 1, c unchanged, W returned unchanged; a further step leaves the failed state untouched; `reducer_failed` routes through the QR state: confirmed.
- Dependency on `_qr_step` in `qr_spectrum.jl`: still defined; the SVD results are bitwise the QR results (AC3a, 1200 orbits over both types).
- Notes statement "0 bytes after warm-up" (AC8, AC9): confirmed in Float64 and Double64, in the section and again at the end of the run (the standard-map W7 measurement was context dependent; here both measurements are 0).
- Notes known limit "GPU compile of exp, log, sqrt for Double64 untested": still untested (AC10 open).
- Notes "coverage not measured": not measured here either.

## Figures

- `lyap-svd-stabilized-exponent-error.png` (seed offset 51300, one seeded start per curve, tau = 1, n up to 1e4 for (b), 40 log-spaced n up to 2000 for (a)). (a) \|lambda_i^SVD(n) - lambda_i^ref(n)\| / theta_sum against n, i = 1 solid and i = 2 dashed, K = 2 pi Float64, K = 2 pi Double64, K = 0.5 Float64; values below 1e-4 are drawn at 1e-4. Worst values over n <= 2000: K = 2 pi Float64 0.00074 and 0.00329; K = 2 pi Double64 1.0e-4 (floor) and 2.6e-4; K = 0.5 0.0263 and 0.0258. All stay far below 1; theta_sum is the unit, so the curves show rounding well inside the tolerance. (b) measured gap lambda_1^SVD - lambda_1^QR with the V7 band ln(1 + c^2)/(2n) to that plus u/(2n), K = 2 pi (red) and K = 0.5 (blue). The gap follows the lower bound and the band is thinner than the line width (u is negligible), so the lines overlap; c settles (K = 2 pi: 2.461 at n = 1000 and 1e4; K = 0.5: 0.2541 then 0.2615, the algebraic convergence of Algorithm item 3). The smallest gap minus lower bound is -2.7e-16 (K = 2 pi) and -7.7e-20 (K = 0.5), rounding below the bound. In (a) the legend covers the right end of the K = 0.5 curves.

## Runtime

Script wall time 1224 s (2 threads, loaded machine, includes about 3 min for the first CairoMakie call and 200 s for the Double64 10^6-step wrapped run). Sections: AC9 457 s (Double64 serial reference of 4096 orbits), AC8 243 s, figure 185 s, AC3 72 s, AC7 41 s, AC4 41 s, AC11 36 s. Above a 10 minute target; most of it is Double64 compilation and long Double64 runs. An earlier full run (before a script fix) took 1252 s.

## Notes

- **Script fix during the run (not a module issue).** The first full run reported AC3(b) Double64 as FAIL (ratios 5.0e8 and 9.6e8): the script evaluated the sandwich term ln(1 + c^2)/(2n) in Float64, whose 1e-19 error sits under a Double64 theta_sum of 1e-28. The term is now evaluated at 256 bits; the final run gives 3.3e-5 and 2.6e-5. The first run's log was overwritten; no src, test or tolerance was changed.
- **AC10 (GPU), run after D40.** `verify/gpu.jl` now has the section; the CPU log line "verify/gpu.jl has no SVDSpectrum section" (line 126 of `lyap-svd-stabilized-verify-log.txt`) predates it and is superseded by `gpu-batchb-verify-log.txt`. The c check failed the absolute 1e-12 (6.1e-12) and a relative 1e-12 (1.39e-12): c is fixed by the first windows and ill-conditioned on some orbits (on orbit 441 a 1-ulp shift of x0[1] moves c by 1.32e-12 on the CPU). Owner's choice: the ensemble 1-ulp bound of Amendment D41, measured in the same run (3.9e-11); it is loose on well-conditioned orbits (median sensitivity 3.3e-16), which the lambda tolerances cover.
- **AC9 thread counts.** Compared: serial against the 2-thread kernel, three chunk sizes. 1 against 4 threads now run in `verify/threads.jl` (bitwise, both types).
- **AC11 spec text.** The oracle uses the identity frame on both sides (see Oracle comparison). The chaotic ensemble difference of means (0.0017 of the allowed 4 s.e.) is small because both sides integrate the same orbits; it tests the estimator output, not statistics of independent samples.
- **AC4 scratch comparison.** The spec's "8 eps at n = 2000" scratch value is not reproduced as an absolute error figure; the criterion is the ratio to theta_sum.
- **Held.** H1 to H3 (general-d algorithm) are not scored in M2.

## Review findings

Independent review: fresh Sonnet subagent given the spec and the code, not the implementation notes. No CRITICAL and no HIGH finding. Disposition: closed (evidence named), fix (proposed change), record (decision row or spec amendment), owner (checkpoint question).

| Severity | Finding | Location | Disposition |
|---|---|---|---|
| MEDIUM | Double64 coverage of AC4 and AC7 thinner than the D35 coverage paragraph states | spec, Amendment D35 Coverage | record: amend the coverage paragraph |
| MEDIUM | no finiteness guard on c and u: exp(S2 - S1) can overflow when det != 1; unreachable for symplectic maps | `src/Lyapunov/svd_spectrum.jl:16-18, 53` | owner: guard (health code) or record as out of scope |
| MEDIUM | AC9 thread count never varied (serial against 2 threads only) | `test/lyap-svd-stabilized` | closed: `verify/threads.jl`, child processes `-t 1` and `-t 4`, outputs bitwise equal (`threads-verify-log.txt`), Float64 and Double64 |
| LOW | AC6 exercises D = 0 only; the V6 test checks finiteness only | `test/lyap-svd-stabilized` | fix: one case with D > 0 |
| LOW | ln sigma_2 could be computed directly (optional, reviewer does not recommend) | `src/Lyapunov/svd_spectrum.jl:12-19` | no change |
| untested | QR health failure inside SVDSpectrum | `src/Lyapunov/svd_spectrum.jl:51` | fix: one test case |
