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

Code review (two independent reviewers, stage 04): no CRITICAL or HIGH findings; the MEDIUM findings were fixed in commit 1361767 and the suite re-run GREEN on Julia 1.13 and 1.10.

Accepted deviations from the spec, owner to confirm at the checkpoint:
- Reducers do not receive the orbit id. The per-orbit body returns a status code and the host attaches the orbit id when it raises the error (`reducer_error`, `_health_error`), so the error record is complete without passing the id into the kernel.
- The CPU-only reducer trait (fallback reason 6) is deferred to lyap-householder-qr, the first module with a LAPACK reducer. Until then reason 6 fires for a reducer whose state is not isbits (D29 amendment).
- Probe per element type: resolved by the D29 amendment (see the GPU section).

Findings on the evidence:
- AC10 covers reasons 4 (device out of memory) and 5 (health error only on the GPU) only through the GPU path; see the GPU section for how each was exercised.
- Runtime 14.6 min against the 10 min target, on a loaded machine; the single-thread child run and the 10^8-step sums dominate. Accepted for a stage 04 script, not for CI.
- AC4: the Frobenius-product bound is loose at large n (ratio 1e-22 at n = 100); the check that bites is at small n (0.04). Recorded, no change.

## Notes

- Nothing failed. The two DEFERRED lines are the ChaosTools exponent comparisons, which move to lyap-benettin stage 04 (D27); AC12 and AC13 wait for `verify/gpu.jl`.
- AC4 is far from marginal, but its bound is very loose at large n: the product of Frobenius norms of the standard-map Jacobian grows even for the regular K = 0.5 orbit (each factor has norm at least sqrt 2), so the ratio at n = 100 is 1e-22 while it peaks at 0.04 (Float64) near small n. The figure shows the normalized error falling away from the line. The relative error against ||W_ref|| (dotted) tracks the bound much more closely.
- The oracle bound 16 K n^2 eps is satisfied with a factor 5.2 (Float64) and 3.6 (Double64) of margin at the worst n (ratio 0.19 and 0.28). The lifted coordinate drifts as n p, so the rounding error of the position grows with n; this is the quantity the bound has to cover and the margin is not large. At K = 2 pi the full-BigFloat comparison is not meaningful (orbit sensitivity) and was not run, as the spec requires K = 0.5 only.
- AC2a gives the same 1.819e-12 for all three tau. That value is one unit in the last place of ln||A^n w_0|| (about 9624), so the exponent bookkeeping is at the rounding floor of the Float64 reference; the bound 8 n eps = 1.8e-11 leaves a factor 10.
- AC6 compares a 1-thread child process against the 4-thread parent through hashes of the bit patterns of x, W, the compensated log sums and the reader state; plus the serial `propagate_orbit` map. All seven comparisons are equal.
- The D28 amendment landed in `src/` while this verification was running. The first run, made against the library before the amendment, marked W9 as PENDING D28; the final run in the log used the amended source and all W9 lines pass (state length 4 with dimension 2). W9h treats AC1 to AC14 re-measured here as the evidence for "s = d unchanged"; the package suite itself was not re-run by this script.
- AC10 covers only what the default probe and test hooks can do: reasons 4 ("device out of memory") and 5 ("NumericalHealthError raised only on the GPU") are per chunk and have no CPU path to exercise; they belong to the GPU run.
- The script took 14.6 minutes (wall) on the loaded machine, above the 10 minute target. The single-thread child run (AC6 digests, AC13 timings) and the AC13 measurements take roughly 5 of those minutes; the 10^8-step compensated sums and the BigFloat curves take most of the rest.
- In the K = 2 pi panel of the first figure the legend covers part of the Float64 curves. The data are all in the log and script; moving the legend below the axes is a one-line change if the figure is used in a paper.
- Timings were taken on a loaded machine (see above); the thread speedup of 4.6 for Float64 on 4 threads is above the thread count and is an artifact of load during the 1-thread run.

## GPU

Measured by `verify/gpu.jl` (decisions D29 and D33), raw output `gpu-verify-log.txt`. Machine: GeForce GTX 1050 (sm_61, 3 GB, display GPU, WDDM driver 581.83), CUDA.jl 6.4.2 with runtime 12.9.0, Julia 1.13.1, 4 CPU threads. Windows TDR registry: TdrDelay 8 s, TdrLevel 3 (the default delay is 2 s). Every GPU line below was checked for `provenance.backend == gpu`; an element type or run that ended on the CPU is printed as FALLBACK(reason), which is a result of the fallback path and not a kernel pass. Since D33 a launch is kept under a 1 s estimated budget (gpu_ns_per_step x N x orbits, at least one wave of 1024 orbits); a run beyond it goes to the CPU with the reason "launch exceeds the display watchdog".

| Criterion | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| AC12b | K = 2 pi, n = 10, same orbit, max over 1024 orbits, \|x_cpu - x_gpu\| | 1.231e-08 | 6.702e-05 | PASS |
| AC12b | K = 10, n = 10, same orbit, max over 1024 orbits | 7.248e-07 | 0.004898 | PASS |
| AC12c | K = 10, 1024 orbits, n = 1e4, tau = 10: max per-orbit \|S_cpu/n - S_gpu/n\| / sigma_FT (sigma_FT = 0.009283); \|difference of means\| / SE | 4.599; 1.936 | 6; 4 | PASS |
| AC12a | K = 0.5, 1024 orbits, n = 1e5, tau = 10 (symmetric form, D33) | - | count diff <= 3; gap <= 2e-6 | FALLBACK(launch exceeds the display watchdog) |
| AC12a | same at n = 4e4 (largest n whose one-wave launch fits the budget; launch 1.4 s): orbits above 1e-6 against the Double64 reference, CPU / GPU; max CPU-GPU gap over the 1006 orbits where both agree with the reference | 18 / 18 (diff 0); 3.016e-07 | count diff <= 3; gap <= 2e-6 | PASS |
| AC12a | same, maxima over all 1024 orbits (chaotic layer included): \|S_cpu - S_ref\| 0.05609, \|S_gpu - S_ref\| 0.06607, \|S_cpu - S_gpu\| 0.04523 | 0.04523 | - | INFO |
| AC13 | Float64, 1e5 orbits x 1e3 steps, GPU 20.2 ns/orbit-step (chunk 8000) vs CPU 4 threads 102.6 ns/orbit-step; speedup 5.08; extrapolated to 1e6 x 1e4: 202 s GPU, 1026 s CPU | 2.02 s vs 10.26 s | GPU <= CPU | PASS |
| AC13 | Float64, vs CPU 1 thread 368.3 ns/orbit-step (child run, 2e4 orbits x 1e3 steps); speedup 18.2 | 2.02 s vs 36.83 s | GPU <= CPU | PASS |
| AC13 | Double64, 2000 orbits x 1e3 steps; CPU 4 threads 3006 ns/orbit-step, 1 thread 2.26e+04 ns/orbit-step | - | GPU <= CPU | FALLBACK(kernel does not compile for the element type) |

Findings.
- AC12a. The D29 form compared CPU and GPU on a CPU-selected class and failed at 1.162e-06 on one orbit. The D33 form compares each backend with a Double64 reference on all orbits. At n = 1e5 the single launch of 1024 orbits is estimated at 2.1 s, above the budget, so `propagate` runs it on the CPU and the FALLBACK line is the result; at n = 4e4 the criterion holds: 18 orbits above 1e-6 for each backend (chaotic layer, about 2 percent) and a CPU-GPU gap of 3.0e-07 on the orbits where both agree with the reference. The maxima over all orbits are chaotic-layer orbits, where no backend pair agrees.
- AC13. On an idle machine the GPU is faster than 4 CPU threads for Float64 (5.1x in this run, 4.3x in the earlier one), so the record has `gpu_faster = true` for Float64 and `auto` resolves to the GPU for Float64 on this machine once CUDA is loaded. The CPU times in the CPU section of this report (185.6 and 851.1 ns/orbit-step) were taken while other jobs ran. The record (`data/gpu-benchmark.toml`, git-ignored) holds the numbers of the last run.
- Double64 does not compile for the GPU: the kernel reaches a lazily initialised `__subtf3` call from Quadmath (reached through DoubleFloats) and allocations of objects with references (`Vector{Base.Padding}`). The probe returns `compiles = false`; `resolve_backend` gives CPU with reason 3 for `:gpu` and `:auto`, and an explicit `:gpu` run completes on the CPU with that reason in its provenance.
- Launch time. A launch is bounded by the number of steps per orbit, not by the number of orbits: 256 and 1024 orbits x 1e4 steps took 0.38 s and 0.30 s, 1024 x 1e5 steps 1.5 to 2.9 s. The launch budget of D33 therefore refuses such a run instead of relying on the watchdog delay.

## Amendment D44 (2026-10-16)

### Evidence

- **Package criteria:**
  - D44a, b, d, e, f, g, the package part of h, and the CPU-backend part of c: `test/tangent-core/tc_d44.jl` and `tc_device.jl`, 783 pass.
  - Two-trajectory D44k(a, b, d): `test/lyap-two-trajectory/l2_d44.jl`, 477 pass.
  - Full `Pkg.test()`: 50197 pass (stage 03 notes).
- **Hardware criteria:** `verify/gpu.jl`, raw output `gpu-d44-verify-log.txt` (all sections except hh-perf) and `gpu-d44c-verify-log.txt` (the D44c rerun).
- **Machine:** as in the GPU section above (GTX 1050, TdrDelay 8 s).
- **Record:** one-wave value 28.06 ns per orbit-step (item 9: 1024 orbits x 1e4 steps, one launch, slower of two runs). AC13 value 23.08 ns; GPU faster than the CPU for Float64. The budget then gives L = 34800 for tau = 10.

| Criterion | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| D44c [local] | Device bitwise, K = 2 pi, N = 2e3, forced L in {tau, 20 tau, 130 tau} against L = N and the unforced run: Benettin, QRSpectrum(RGS), SVDSpectrum at tau = 2 (1024 orbits each) and tau = 10 (Benettin 1024 orbits; RGS and SVD 58 orbits, see findings) | equal (6 rows) | bitwise | PASS |
| D44h [local] | The extension's lost-device codes 702, 719, 700, 709, 4, 100, 214, 714, 715, 716, 201, 710, 999: reason of each; codes 2 (out of memory) and 1 (invalid value) are not lost-device codes | "device lost" for every code | "device lost", 700 and 702 included | PASS |
| D44i | One stepping body: sites of `forward(sys` and `tangent(sys` in `src/Tangent`; callers of `_orbit_body` (`propagate.jl:36` kernel, `orbit.jl:69` `propagate_orbit`) | 1 / 1; 2 | 1 / 1; 2 | PASS |
| D44j [local] AC12a | tangent-core LogNorm S/n, 1024 orbits, n = 1e5, tau = 10: orbits above 1e-6 vs the Double64 reference, CPU / GPU; gap; per-launch estimate; measured full launch (ratio) | 20 / 20; 9.444e-08; 1.0 s; 1.508 s (1.51) | diff <= 3; <= 2e-6; <= 1 s; <= 4 s | PASS |
| D44j [local] lyap-benettin AC16 | Benettin, same form | 18 / 17; 1.754e-07; 1.0 s; 1.334 s (1.33) | same | PASS |
| D44j [local] lyap-gram-schmidt AC12 | QRSpectrum(RGS), lambda_1 and lambda_2 (S7 regular class) | 1 / 1; 5.434e-07; 1.0 s; 1.341 s (1.34) | same | PASS |
| D44j [local] lyap-gram-schmidt AC12 | QRSpectrum(MGS) | 1 / 1; 5.434e-07; 1.0 s; 1.265 s (1.27) | same | PASS |
| D44j [local] lyap-householder-qr AC14a | HouseholderQR; \|lambda_1 + lambda_2\| / theta_sum, CPU and GPU | 1 / 1; 5.434e-07; 1.0 s; 1.415 s (1.42); 3.1e-4 and 2.6e-4 | same; <= 1 | PASS |
| D44j [local] lyap-householder-qr AC14a | GivensQR; theta_sum ratio | 1 / 1; 5.434e-07; 1.0 s; 2.009 s (2.01); 2.1e-4 and 2.2e-4 | same; <= 1 | PASS |
| D44j [local] lyap-householder-qr AC14a | AugmentedQR; theta_sum ratio | 1 / 1; 5.434e-07; 1.0 s; 2.983 s (2.98); 3.1e-4 and 2.6e-4 | same; <= 1 | PASS |
| D44j [local] lyap-svd-stabilized AC10a | SVDSpectrum, lambda_1 and lambda_2 | 1 / 1; 5.434e-07; 1.0 s; 2.679 s (2.68) | same | PASS |
| D44k(c) [local] | propagate_shadow device bitwise, lifted and wrapped, K = 0.5 and 2 pi, tau = 2 and 10, 256 orbits, N = 2e3, forced L in {tau, 20 tau, 130 tau} against L = N | equal (8 rows) | bitwise | PASS |
| AC12a (rerun) | n = 1e5, the original criterion, now on the device in 3 launches of 34800 steps (GPU run 2.46 s in all) | 20 / 20 (diff 0); 9.444e-08 | count diff <= 3; gap <= 2e-6 | PASS (was FALLBACK) |

Every D44j row ran on the device: backend gpu, empty reason, 1 device chunk of 1024 orbits, L = 34800, 3 launches per chunk. Each full launch was measured in one chunk of one launch, as the provenance confirms.

The Batch B rows at n = 4e4, the AC12b and c rows and the chaotic K = 10 rows all pass again with the relaunch code, at the same values as `gpu-batchb-verify-log.txt`. The 13 FALLBACK lines are all Double64, "kernel does not compile for the element type", unchanged.

### Findings

- **D44c at tau = 10.**
  - The first run failed for RGS and SVD. Every run of those two rows, the single launch included, ended on the CPU with reason 5: on some orbits only the device raised "orthonormalisation breakdown". The comparison was therefore CPU against CPU (`d44c_diag.jl`: all equal, all CPU).
  - The procedure now also drops orbits whose single-launch device run raises; this extends the D44a correction to the device side. After that, 58 of 1024 orbits remain: 780 break down on the CPU and 186 more only on the device. All rows pass.
  - K = 2 pi with tau = 10 is far outside the usable RGS window, so the bitwise check at tau = 10 rests on 58 orbits. The tau = 2 rows keep all 1024.
  - Owner decision: accept 58 orbits, or move the tau = 10 rows to a smaller K.
- **The launch estimate is low for reducers heavier than the record's.**
  - The one-wave record uses the LogNorm reducer. Measured full launches at L = 34800 take 1.27 to 2.98 s against the 1 s estimate (AugmentedQR is the slowest).
  - All are within the 4 s gate and the 8 s TdrDelay, but the margin to the watchdog is 2.7x, not 8x.
  - Candidate for `_config/risks.md`: a per-reducer record, or a cost factor, if heavier reducers or slower cards follow.
- **CPU timing.** The 4-thread CPU value of this run (218 ns per orbit-step) is twice the earlier 102.6. Other jobs ran during part of the run. `gpu_faster` stays true either way.

### Review findings (independent Sonnet review: spec and code only)

None CRITICAL or HIGH.

| # | Severity | Finding | Disposition |
|---|---|---|---|
| 1 | MEDIUM | The provenance fields take the device L as soon as device setup succeeds. A first chunk that fails in launch 1 still reports L, although item 13 says 0 when "the device was never used". | Kept. The owner-approved test D44k expects L after a lost device in launch 3 of chunk 1. Read with item 13: the device was used once launches ran. The only case left open is a failure in launch 1 of chunk 1. Owner decision. |
| 2 | LOW | `_budget_steps` upward loop unbounded for a tiny, corrupt record value | Fixed: m capped at `typemax(Int) / (4 tau)` |
| 3 | LOW | A forced L with a setup fallback records the forced L on the CPU | As specified: item 13 says a forced CPU step chunk records its L (test hook only) |
| 4 | LOW | `Val(size(Wi, 1))` in the kernel may not constant-fold on the GPU | Refuted on hardware: the probe compiles and every Float64 row ran on the device |
| 5 | LOW | The shadow start is computed again on the device at launch 1, besides the host check | Kept: it keeps the single-launch path unchanged (item 4), and the device needs the start code |
| 6 | LOW | An OOM or lost device inside `compiles()` escapes the probe | Pre-existing (D29), unchanged by D44 |
| 7 | LOW | The record run aborts the script on a fallback | Intended |

## Amendment D46: stage 04 on the RTX 5060 (2026-10-17)

### Evidence

- **Machine:** NVIDIA GeForce RTX 5060 Laptop GPU (compute capability 12.0, 7.96 GiB), Intel Core Ultra 9 275HX, TdrDelay 8 s, CUDA.jl 6.4.2 (verify Manifest), Julia 1.13.2, 4 threads. From D46 on, this is the machine of every [local] criterion; the GTX 1050 sections above stay on record.
- **Record (`data/gpu-benchmark.toml`, per machine, git-ignored):** Float64 GPU 2.468 ns per orbit-step (AC13 ensemble), one-wave value 9.137 ns (item 9; GTX 1050: 28.06), CPU 8.255 ns with 4 threads and 33.09 ns with 1 thread; GPU faster (speedup 3.3 against 4 threads, 13.4 against 1). Double64 falls back, "kernel does not compile for the element type", as on the GTX 1050.
- **Full run before D46:** `gpu-rtx5060-verify-log.txt` (893 s): PASS 64, FAIL 8, FALLBACK 13 (all Double64), INFO 5. The 8 FAIL lines are the D44j rows, which failed only their precondition "more than 1 launch": the budget fits n = 1e5 in one launch on this card (L = 100000). Everything else of Batch A, Batch B and D44 passed.
- **D46 rerun (tag d44):** `gpu-rtx5060-d46-verify-log.txt`: PASS 33 of 33, after the D46 corrections (spec Amendment D44, corrections at stage 04 on the RTX 5060) and the item 13 code change.
- **Package suite:** 49621 pass, 0 fail (D46 tests: 25 in `test/tangent-core/tc_d44.jl`, 16 in `test/lyap-two-trajectory/l2_d44.jl`), Aqua and JET clean.

| Criterion | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| D44c [local] | Device bitwise, N = 2e3, forced L in {tau, 20 tau, 130 tau} against L = N and the unforced run: Benettin, QRSpectrum(RGS), SVDSpectrum; tau = 2 at K = 2 pi, tau = 10 at K = 0.5 (D46) | equal (6 rows), 1024 of 1024 orbits kept in every row | bitwise | PASS |
| D44h [local] | Lost-device codes of the extension table | "device lost" for every code | as D44 | PASS |
| D44i | One stepping body | 1 / 1; 2 | 1 / 1; 2 | PASS |
| D44j [local] AC12a | tangent-core LogNorm, n = 1e5, tau = 10, L = 1e5, 1 launch per chunk: orbits above 1e-6 CPU / GPU; gap; estimate; measured full launch | 20 / 20; 9.444e-08; 1.0 s; 0.345 s | diff <= 3; <= 2e-6; <= 1 s; <= 4 s | PASS |
| D44j [local] lyap-benettin AC16 | Benettin, same form | 18 / 17; 1.754e-07; 1.0 s; 0.370 s | same | PASS |
| D44j [local] lyap-gram-schmidt AC12 | QRSpectrum(RGS) and (MGS), lambda_1 and lambda_2 | 1 / 1; 5.434e-07; 1.0 s; 0.582 s and 0.500 s | same | PASS |
| D44j [local] lyap-householder-qr AC14a | Householder, Givens, Augmented; theta_sum ratio CPU and GPU | 1 / 1; 5.434e-07; 0.497 s, 0.473 s, 0.492 s; at most 3.1e-4 | same; <= 1 | PASS |
| D44j [local] lyap-svd-stabilized AC10a | SVDSpectrum | 1 / 1; 5.434e-07; 1.0 s; 0.617 s | same | PASS |
| D44k(c) [local] | propagate_shadow device bitwise, lifted and wrapped, K = 0.5 and 2 pi, tau = 2 and 10 | equal (8 rows) | bitwise | PASS |
| Item 13 (D46) | Provenance after a fault in launch 1 of chunk 1 (lost device, out of memory, compile error): 0 and 0 unforced, the CPU values forced; device values kept after launch 3 or chunk 2 | as specified | exact | PASS (package suite) |

### Findings

- The D44j values equal the GTX 1050 values to every printed digit (same seeds, same Float64 kernel), so the relaunch at three launches and the single launch here agree with the same Double64 reference.
- With one launch per chunk at n = 1e5, the relaunch path on this card is exercised by D44c and D44k at N = 2e3 with forced L (owner, D46).
- The measured full launches (0.34 to 0.62 s) are below the 1 s estimate, so the launch estimate is conservative on this card; the GTX 1050 finding about heavier reducers (1.3 to 3.0 s against 1 s) does not reproduce here.
