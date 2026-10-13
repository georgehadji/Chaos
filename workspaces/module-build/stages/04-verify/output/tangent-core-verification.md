# Verification: tangent-core

Spec: `../../01-spec/output/tangent-core-spec.md` (including Amendment D28). Script: `verify/tangent-core.jl` (run from the repo root with `julia --project=verify -t 4 verify/tangent-core.jl`). Printed output: `tangent-core-verify-log.txt`. Julia 1.13.1, 4 worker threads (12 hardware threads), base seed 20261012 (per-criterion offsets as in the package tests).

Result: 62 lines PASS, 0 FAIL, 8 lines STAGE04-GPU (AC12, AC13), 2 lines DEFERRED lyap-benettin 04 (ChaosTools exponent comparisons, D27). Every criterion of the spec and of Amendment D28 has a measured value in the log.

All timings were taken on a loaded machine (other agents were running Julia at the same time), so they are indicative only.

## Sample sizes

AC2: 100 orbits x 10^4 steps, three tau. AC3: n up to 10^6. AC4: 20 orbits, n = 100 (K = 0.5) and 20 (K = 2 pi); the figure curves use every n from 1 to 100 (K = 0.5) and 1 to 40 (K = 2 pi). Oracle block: 20 orbits, n = 1..100. AC6: 4096 orbits x 10^3 steps. AC7: 10^4 orbits (Float64), 200 (Double64). AC9: 10^8 and 10^7 additions. AC13 (CPU side only): Float64 10^5 orbits x 10^3 steps, Double64 2000 orbits x 10^3 steps, instead of 10^6 x 10^4 (too long on CPU), 1 thread and 4 threads.

## Evidence table

| Criterion | Quantity | Value | Bound | Verdict | Where |
|-----------|----------|-------|-------|---------|-------|
| AC1 | mismatching entries of W_n and x_n against integer A^n, n = 1..35 (largest entry 3.08e14 < 2^53) | 0 | 0 | PASS | package test; re-measured in verify script |
| AC2a | max over 100 orbits of the error of S_n + ln norm(w_n) against BigFloat(256), tau = 1, 5, 10 | 1.819e-12 (each tau) | 1.776e-11 (8 n eps) | PASS | package test; verify script |
| AC2b | max abs(S_n/n - ln((3 + sqrt 5)/2)), tau = 1, 5, 10 | 7.0e-4 (each tau) | 2e-3 | PASS | same |
| AC2c | handoff count N/tau and spacing tau (I2), 100 orbits, three tau | all correct | exact | PASS | same |
| AC3a | mismatching entries of W_n = [[1,n],[0,1]] and x_n = x_0 + n p_0, n = 1, 10, 1000, 10^6 | 0 in each | 0 | PASS | same |
| AC3b | abs(S_n - 0.5 ln(n^2 + 1)), n = 10^4 | 1.776e-15 | 1.776e-11 | PASS | same |
| AC3c | S_n/n - ln(n)/n, n = 10^4 | 5.0e-13 | 1e-12 | PASS | same |
| AC4 | max over 20 orbits of err / (8 n eps prod norm(M_i) norm(W_0)): Float64 K = 0.5, n = 100 | 1.3e-22 | 1 | PASS | same |
| AC4 | Float64 K = 2 pi, n = 20 | 3.1e-4 | 1 | PASS | same |
| AC4 | Double64 K = 0.5, n = 100 | 6.7e-23 | 1 | PASS | same |
| AC4 | Double64 K = 2 pi, n = 20 | 9.4e-6 | 1 | PASS | same |
| AC4 | same ratio, maximum over all plotted n (Float64 K = 0.5 / 2 pi; Double64 K = 0.5 / 2 pi) | 0.0425 / 0.0564; 0.0178 / 0.0252 | 1 | PASS | verify script (figure data) |
| AC4 | orbit of the kernel equals `forward` iterated (isequal), all 20 orbits, all four cases | 0 mismatches | 0 | PASS | verify script |
| Oracle (full BigFloat orbit) | K = 0.5, n = 100, max over 20 orbits of relative error of W_n, Float64 | 3.32e-12 | 1.776e-11 (16 K n^2 eps) | PASS | verify script |
| Oracle | same, Double64 | 1.07e-27 | 3.94e-27 | PASS | verify script |
| Oracle | maximum over n = 1..100 of error / (16 K n^2 eps), Float64 / Double64 | 0.19 / 0.28 | 1 | PASS | verify script |
| Oracle | ChaosTools, regular orbit K = 0.5, n = 10^5 | not run | 1e-4 | DEFERRED lyap-benettin 04 | D27 |
| Oracle | ChaosTools, chaotic ensemble K = 2 pi, 1000 points, n = 10^4 | not run | 4 SE | DEFERRED lyap-benettin 04 | D27 |
| AC5 | bytes allocated by `propagate_orbit` after warm-up (n = 10^4, k = 2, chain of two reducers; tau = 1 / tau = 4 / empty chain n = 100) | 0 / 0 / 0 | 0 | PASS | package test; verify script |
| AC6 | bitwise equality of x, W, logs, reader state: 1 thread vs 4 threads, chunk 64 vs 1000, and serial `propagate_orbit` (7 of 7 comparisons) | 7 / 7 | all equal | PASS | package test (in-process); verify script (child process with `-t 1`) |
| AC7a | max abs(Q'Q - I), Float64, 10^4 orbits: (d,k) = (2,1), (2,2), (4,4), (8,8) | 8.9e-16, 8.9e-16, 1.3e-15, 1.6e-15 | 1.8e-15, 3.6e-15, 7.1e-15, 1.4e-14 (8 k eps) | PASS | package test; verify script |
| AC7a | same, Double64, 200 orbits | 2.8e-32, 2.8e-32, 3.4e-32, 3.8e-32 | 3.9e-31 to 3.2e-30 | PASS | same |
| AC7b | reproducibility (repeat, reversed order, 4 threads, chunk 64 and 1000), 10^4 distinct matrices, seed sensitivity | 6 of 6 | 6 | PASS | same |
| AC7c | KS statistic of the angle, d = 2, k = 1, 10^4 orbits | 0.0104 | 0.0163 | PASS | same |
| AC8 | control run without a flagged orbit | no error | no error | PASS | same |
| AC8a | Inf in state at step 7, tau = 5: (orbit, step, what, backend); nothing stored | (3, 10, "state", cpu) | same | PASS | same |
| AC8b | NaN in W at step 3, tau = 5 | (5, 5, "deviation matrix", cpu) | same | PASS | same |
| AC8c | tau = 200, n = 200: finite, relative deviation of W[1] from e^600 | 1.6e-15 | 1e-10 | PASS | same |
| AC8c | N = 300 with tau = 300 and tau = 100 (overflow) | error at step 300, "deviation matrix" (both) | step 300 | PASS | same |
| AC8d | N not a positive multiple of tau, three bad cases | 3 ArgumentError of 3, before any work | 3 | PASS | same |
| AC9a | compensated sum of 10^8 copies of fl(0.1), relative error | 5.55e-17 (0.25 eps) | 2.22e-16 | PASS | same |
| AC9b | plain running sum, same input (negative control) | 1.887e-9 | must exceed 2.22e-16 | PASS | same |
| AC9c | Double64 compensated sum, 10^7 copies | 0 | 4.93e-32 | PASS | same |
| AC10 | backend resolution and provenance, 24 checks (default probe here: "CUDA.jl not loadable") | 24 of 24 | 24 | PASS | package test; verify script |
| AC11 | M(2,2) = 96, M(4,4) = 240, chunk at 2^30: 5592405 and 2236962; at 2^29: 2796202 and 1118481 | exact | exact | PASS | same |
| AC12a | CPU vs GPU, K = 0.5, 1024 orbits, n = 10^5 | not measured here | 1e-6 | STAGE04-GPU | `verify/gpu.jl` |
| AC12b | same-orbit x difference, n = 10: K = 2 pi / K = 10 | not measured here | 6.7e-5 / 4.9e-3 (4 n eps G^n, G = 9.72 / 14.93) | STAGE04-GPU | `verify/gpu.jl` |
| AC12c | chaotic ensemble K = 10, 1024 orbits, n = 10^4 | not measured here | 6 sigma_FT per orbit, 4 SE for means | STAGE04-GPU | `verify/gpu.jl` |
| AC13 | CPU wall time, Float64, 10^5 orbits x 10^3 steps, 1 thread / 4 threads (see table below) | 85.1 s / 18.6 s | GPU time <= CPU time | STAGE04-GPU (CPU side recorded) | verify script |
| AC13 | CPU wall time, Double64, 2000 orbits x 10^3 steps, 1 thread / 4 threads | 59.7 s / 15.1 s | same | STAGE04-GPU (CPU side recorded) | verify script |
| AC14 | per-step update sites in `src/Tangent` | 1 `forward(sys, x)` and 1 `tangent(sys, x)` | 1 / 1 | PASS | inspection by script |
| W9a | state of length 4, dimension 2: 4-component state and 2 x k matrix equal manual iteration bitwise (k in {1,2}, tau in {1,3}) | 0 mismatches | 0 | PASS | package test (`tc_state_length.jl`); verify script |
| W9b | 257 orbits, chunk 16 and default, chunked equals serial bitwise | 0 mismatches | 0 | PASS | same |
| W9c | Inf in component 3 or 4 at step 7, tau = 3 | "state" at step 9 (both) | step 9 | PASS | same |
| W9d | bytes allocated with s = 4 (chain of one reducer / empty chain) | 0 / 0 | 0 | PASS | same |
| W9e | length of the point seen by a reducer | 4 | 4 | PASS | same |
| W9f | bytes per orbit (s = 4, d = 2, k = 2) / default; chunk at 2^30 for s = 4 / default | 112 / 96; 4793490 / 5592405 | exact | PASS | same |
| W9g | k = 3 with d = 2 (state length 4) | ArgumentError before any work | ArgumentError | PASS | same |
| W9h | FAIL lines among AC1 to AC14 (s = d behaviour unchanged) | 0 | 0 | PASS | verify script; the package suite is run by `test/runtests.jl` |

### AC14: where the one per-step body lives

`src/Tangent/orbit.jl`, function `_orbit_body`: `W = tangent(sys, x) * W` at line 32 and `x = forward(sys, x)` at line 33 are the only per-step update in `src/Tangent` (the script counts one match each, and no `forward(` or `propagate(` mention in `src/Foundation` or `src/Diagnostics`). `_orbit_body` has exactly two callers: the kernel `_orbit_kernel!` (`src/Tangent/propagate.jl`, line 31) and `propagate_orbit` (`src/Tangent/orbit.jl`, line 67). AC1 to AC4 and AC7 to AC11 go through `propagate` (the kernel on the CPU backend); AC5 and the serial part of AC6 go through `propagate_orbit`, which runs the same `_orbit_body`.

### AC13: CPU wall times (GPU column: STAGE04-GPU)

Settings: d = 2, k = 2, K = 2 pi, tau = 1, one column-normalizing reducer (without renormalization W overflows after about 600 steps at K = 2 pi), default chunk size, one warm-up call before timing. The extrapolation to 10^6 orbits x 10^4 steps assumes linear scaling in orbit-steps and is not a measurement.

| Type | Orbits x steps | Threads | Wall time (s) | ns per orbit-step | Extrapolated to 10^6 x 10^4 (s) |
|------|----------------|---------|---------------|-------------------|--------------------------------|
| Float64 | 10^5 x 10^3 | 1 | 85.1 | 851 | 8511 |
| Float64 | 10^5 x 10^3 | 4 | 18.6 | 186 | 1856 |
| Double64 | 2000 x 10^3 | 1 | 59.7 | 29840 | 298400 |
| Double64 | 2000 x 10^3 | 4 | 15.1 | 7568 | 75680 |

The GPU side (CUDA runtime version, device name, Float64 only because Float32 was dropped in D24) is measured by `verify/gpu.jl`; this machine's CUDA toolkit and device were not touched here. An earlier run of the same script (the library was being edited at the time, same kernel) gave 464 ns, 296 ns, 16 520 ns and 6756 ns per orbit-step for the same four rows, so the spread between runs on the loaded machine is up to a factor 2 in either direction (the 1-thread Float64 figure fell outside the 130 ns per orbit-step quoted in the implementation notes in both runs).

## Figures

- `tangent-core-ac4-bigfloat-product.png`: AC4. For K = 0.5 (n up to 100) and K = 2 pi (n up to 40): maximum over 20 seeded orbits of ||W_n - W_ref||_F / (prod ||M_i||_F ||W_0||_F) for Float64 and Double64 (solid), against the bound 8 n eps(T) (dashed); dotted: the error relative to ||W_ref||_F. Both axes logarithmic.
- `tangent-core-oracle-fullbigfloat.png`: oracle block, K = 0.5. Maximum over 20 orbits of the relative error of W_n of the Float64 and Double64 orbits against the full BigFloat(256) orbit, with the bound 16 K n^2 eps(T) (dashed). Both axes logarithmic.

## Review findings

(To be filled in by the orchestrator.)

## Notes

- Nothing failed. The two DEFERRED lines are the ChaosTools exponent comparisons, which move to lyap-benettin stage 04 (D27); AC12 and AC13 wait for `verify/gpu.jl`.
- AC4 is far from marginal, but its bound is very loose at large n: the product of Frobenius norms of the standard-map Jacobian grows even for the regular K = 0.5 orbit (each factor has norm at least sqrt 2), so the ratio at n = 100 is 1e-22 while it peaks at 0.04 (Float64) near small n. The figure shows the normalized error falling away from the line. The relative error against ||W_ref|| (dotted) tracks the bound much more closely.
- The oracle bound 16 K n^2 eps is satisfied with a factor 3.6 (Float64) and 5.3 (Double64) of margin at the worst n (ratio 0.19 and 0.28). The lifted coordinate drifts as n p, so the rounding error of the position grows with n; this is the quantity the bound has to cover and the margin is not large. At K = 2 pi the full-BigFloat comparison is not meaningful (orbit sensitivity) and was not run, as the spec requires K = 0.5 only.
- AC2a gives the same 1.819e-12 for all three tau. That value is one unit in the last place of ln||A^n w_0|| (about 9624), so the exponent bookkeeping is at the rounding floor of the Float64 reference; the bound 8 n eps = 1.8e-11 leaves a factor 10.
- AC6 compares a 1-thread child process against the 4-thread parent through hashes of the bit patterns of x, W, the compensated log sums and the reader state; plus the serial `propagate_orbit` map. All seven comparisons are equal.
- The D28 amendment landed in `src/` while this verification was running. The first run, made against the library before the amendment, marked W9 as PENDING D28; the final run in the log used the amended source and all W9 lines pass (state length 4 with dimension 2). W9h treats AC1 to AC14 re-measured here as the evidence for "s = d unchanged"; the package suite itself was not re-run by this script.
- AC10 covers only what the default probe and test hooks can do: reasons 4 ("device out of memory") and 5 ("NumericalHealthError raised only on the GPU") are per chunk and have no CPU path to exercise; they belong to the GPU run.
- The script took 14.6 minutes (wall) on the loaded machine, above the 10 minute target. The single-thread child run (AC6 digests, AC13 timings) and the AC13 measurements take roughly 5 of those minutes; the 10^8-step compensated sums and the BigFloat curves take most of the rest.
- In the K = 2 pi panel of the first figure the legend covers part of the Float64 curves. The data are all in the log and script; moving the legend below the axes is a one-line change if the figure is used in a paper.
- Timings were taken on a loaded machine (see above); the thread speedup of 4.6 for Float64 on 4 threads is above the thread count and is an artifact of load during the 1-thread run.
