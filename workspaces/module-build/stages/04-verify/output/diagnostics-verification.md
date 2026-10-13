# Verification: diagnostics

Spec: `workspaces/module-build/stages/01-spec/output/diagnostics-spec.md` (no Amendments section). Script: `verify/diagnostics.jl`, run as `julia --project=verify -t 4 verify/diagnostics.jl` (Julia 1.13.1, 4 threads, base seed 20261012). Raw output: `diagnostics-verify-log.txt` (83 lines: 74 PASS, 0 FAIL, 1 DEFERRED, 1 STAGE04-GPU, 2 INFO). Timings were taken on a loaded machine (total wall time 14 min; the run competed with other agents). Element types: Float64 and Double64 (Float32 dropped, D24).

Sample sizes: AC1 4 maps x 1000 steps; AC2 10^6 steps; AC3 10^5 steps; AC6 5 K x 50 starts x 10^4 steps (plus the 10^6-step AC7 orbits); AC7 5 K x 10^6 steps, one orbit per K, full size (not reduced); AC8 10^3 steps at K = 5, 10^5 steps per K for the 256-bit det comparison, 10^4 steps per K for the 256-bit estimator; AC10 64 orbits x 100 steps.

## Evidence table

Value is the worst case over the rows listed in the log. "Pkg" is the package suite (`test/diagnostics/runtests.jl`), "verify" is `verify/diagnostics.jl`.

| AC | Quantity | Value | Bound | Verdict | Where |
|----|----------|-------|-------|---------|-------|
| 1 | max eps_det, eps_symp, monitor maxima, 4 integer maps x 1000 steps, health errors | 0 (both types) | 0 | PASS | Pkg, re-measured in verify |
| 2 | max eps_det / eps, rotation 0.7311, 10^6 steps | 0 (Float64), 0.25 (Double64) | 8 | PASS | Pkg, re-measured in verify |
| 2 | max eps_orth / eps (reference MGS) | 2 (Float64), 1 (Double64) | 8 | PASS | same |
| 3 | G_max, cat map, n = 10^5 | 6.8541 = mu^2 (rtol 1e-12) | mu^2 | PASS | same |
| 3 | max_n eps_sum / (eps (G_max/tau + M_max)), 0 violations of theta_sum | 0.282 | 16 (theta_sum) | PASS | same |
| 3 | eps_pair equals eps_sum (N = 1), lambda_1 within 1e-3 of ln mu | gap 0 | 0 | PASS | same |
| 4 | N = 1 and N = 2 synthetic pairing, error / tolerance; adjacent pairing would give 1.27, measured eps_pair 3.8e-6 | 0 (both types) | 1 | PASS | same |
| 5 | det fault, both sides of theta (bitwise value, orbit, step) | 1/1 raise, 1/1 quiet (both types) | all | PASS | same |
| 5 | orthogonality fault, max abs(eps_orth - (2d + d^2)) / eps | 0.004 (Float64), 0 (Double64); 1/1, 1/1 | all | PASS | same |
| 5 | sum and pairing fault, both sides of theta_sum | 1/1, 1/1 (both types) | all | PASS | same |
| 5 | symplecticity fault (4 x 4 block) both sides, det-1 non-symplectic map named "symplecticity" | 1/1, 1/1, yes (both types) | all | PASS | same |
| 5 | non-finite (NaN, Inf in threshold, DT, Q, lambda; threshold Inf), right name; never clipped | 16/16, yes (both types) | all | PASS | same |
| 6 | max eps_det / (eps M), 50 starts x 10^4 steps, K = 0.5, 0.9716354, 2, 5, 10 | 0.500 (Float64, all K), 0.242 to 0.250 (Double64) | 1.5 | PASS | Pkg, re-measured in verify |
| 6 | same quantity on the 10^6-step AC7 orbits | 0 to 0.500 (Float64), 0.125 to 0.250 (Double64) | 1.5 | PASS | verify |
| 6 | threshold violations at the proposed keys | 0 in every case | 0 | PASS | verify |
| 7 | max_n eps_sum / (eps (G_max/tau + M_max)), 10^6 steps, K panel, Float64 | 0.196 (K=0.5), 0.187 (0.9716354), 0.043 (2), 0.107 (5), 0.086 (10) | 4 | PASS | verify |
| 7 | same, Double64 | 0.072, 0.044, 0.051, 0.051, 0.032 | 4 | PASS | verify |
| 7 | snapshots above 4, above 16 | 0, 0 in all 10 cases | 0 | PASS | verify |
| 7 | ChaosTools full-spectrum sum next to eps_sum | n/a | n/a | DEFERRED lyap-benettin 04 | D27 |
| 8 | max abs(eps_det(binary64) - eps_det(256 bit)) / eps, K = 5, 10^3 steps | 0.5 | 1 | PASS | Pkg, re-measured in verify |
| 8 | same, K panel, 10^5 steps per K | 0.5 in every K | 1 | PASS | verify |
| 8 | 256-bit QR estimator on entries recomputed at the stored points: max eps_sum / (eps_256 (G + M)), 10^4 steps | 0.059 to 0.209 | 16 | PASS | verify |
| 9 | bytes allocated per call after warm-up (propagate_orbit, reducer_update, monitor_exponents) | 0/0/0 (both types) | 0 | PASS | Pkg, re-measured in verify |
| 10 | 64 orbits, x, W, residual maxima bitwise equal, 4 threads vs 1 thread (separate process); equal to serial map | true, true | bitwise | PASS | verify |
| 10 | reported orbit id, fault targets (33), (17,40), (40,17), 4 threads and 1 thread | 33, 17, 17 both | 33, 17, 17 | PASS | verify |
| 11 | CPU kernel: item 1 (4 integer maps, both types), item 5 (det fault both sides), item 6 (panel through `propagate`, 0 health errors, det_max bitwise equal to a direct loop on the same starts) | max residual 0; 1/1, 1/1; gap 0 | as items | PASS | verify |
| 11 | `resolve_backend(:gpu, BigFloat)`: cpu, reason "kernel does not compile for the element type" | as stated | as stated | PASS | verify |
| 11 | GPU-vs-CPU on the GTX 1050 (items 1, 5, 6) | see `verify/gpu.jl` | <= threshold in force | STAGE04-GPU | later run |

Informational: `reducer_update` time per call 304 ns (Float64) and 5.19 us (Double64), loaded machine (AC9 INFO rows). Per-case AC7 runtimes 20 to 117 s (Double64 K = 5 the slowest at 117 s).

## Figures

- `diagnostics-det-residual.png`: per-start max eps_det / (eps M) (50 starts x 10^4 steps, dots) and the 10^6-step orbit (stars) per K, Float64 and Double64, with the B1 bound 1.5.
- `diagnostics-sum-residual.png`: eps_sum(n) / (eps (G_max/tau + M_max)) against n, running maximum per log bin (20 bins per decade), all K and both types, with the bound 4 (AC7) and the key value 16. Values of exactly 0 are clipped to 1e-3 for the log axis (Float64 at K = 0.5 and 0.9716354 reach exactly 0 after about 3 x 10^4 steps).

## Review findings

Code review (two independent reviewers, stage 04): no CRITICAL or HIGH findings; the MEDIUM findings (pending failure wins in `monitor_exponents`, one shared `_failure`) were fixed in commit 1361767. The Tangent patch adding `reducer_step`, `reducer_failed` and `reducer_error` lets the monitor report its failure with the orbit id.

Findings on the evidence:
- Criterion 7: max 0.196 against the pass bound 4, so `pairing_residual_tol = 16` sits about 80 times above the measurement. The spec allows tightening; no change until lyap-* runs give more orbits.
- Criterion 11, kernel path: the monitor's `reducer_final` does not export M_max, so the script checks det_max / eps <= threshold instead of eps_det / (eps M). With M >= 1 this is the stricter check. Accepted.
- Criterion 7 uses one orbit per K, not an ensemble; recorded as a sample-size limit of this stage, revisited in lyap-gram-schmidt stage 04.

## Notes

- Nothing failed. The criterion 7 pass bound (4) is met with a factor of 20 or more: the measured maximum is 0.196 (Float64, K = 0.5), so `pairing_residual_tol = 16` sits 82 times above the measurement and 4 times above the derived constant; the spec allows tightening it, not loosening it. No revision of the key is needed on this evidence.
- Criterion 6 is not marginal: the Float64 value is 0.5 (one half, the rounding of the subtraction), matching the spec author's scratch figure, and Double64 stays at 0.25. The bound 1.5 is 3 to 6 times larger. The B1 first-order derivation assumes the double-double operations obey the same u bound; the measured Double64 values give no sign of a violation.
- Criterion 7 uses one orbit per K (seed SEED + 700 + K index, lifted coordinates, tau = 1, identity start basis), not an ensemble; K = 0.5 and 0.9716354 are regular on these starts (lambda_1 below 1e-4), so G_max there (2.0 and 4.6) is small. G is taken from the R factor of the same factorisation.
- The 256-bit estimator (AC8 sum row) shows the exact-matrix sum residual at the 256-bit rounding level (0.06 to 0.21 in units of eps_256 (G + M)), as identity I2 predicts. The binary64 lambda_1 differs from the 256-bit value by up to 3e-15 over 10^4 steps (K = 0.5, 0.9716354); this is information, not a criterion (the binary64 matrices are rounded versions of the 256-bit ones).
- Criterion 10 compares 4 threads with 1 thread through a child Julia process started by the script (`--child`); the digest covers x, W and the five residual maxima of all 64 orbits.
- Criterion 11: the GPU part (`STAGE04-GPU`) is not run here. The monitor's `reducer_final` does not export M_max, so the kernel-path panel check reports det_max / eps (max 4.0 over the panel, K up to 10 so M up to 11) and a bitwise match against a direct loop instead of eps_det / (eps M).
- D28 wrapped-state names are not used by this module; no criterion is PENDING D28.

## GPU

Measured by `verify/gpu.jl` (decisions D29 and D33), raw output `gpu-verify-log.txt`. Machine: GeForce GTX 1050 (sm_61, 3 GB, display GPU, WDDM driver 581.83), CUDA.jl 6.4.2 with runtime 12.9.0, Julia 1.13.1, 4 CPU threads. Windows TDR registry: TdrDelay 8 s, TdrLevel 3 (the default delay is 2 s). Every GPU line below was checked for `provenance.backend == gpu`; an element type or run that ended on the CPU is printed as FALLBACK(reason), which is a result of the fallback path and not a kernel pass. Since D33 a launch is kept under a 1 s estimated budget (gpu_ns_per_step x N x orbits, at least one wave of 1024 orbits); a run beyond it goes to the CPU with the reason "launch exceeds the display watchdog".

| Criterion | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| AC11 item 1 | integer maps x 1000 steps, Float64: max residual, GPU / CPU | 0 / 0 | 0 (exact) | PASS |
| AC11 item 1 | same, Double64 | - | 0 (exact) | FALLBACK(kernel does not compile for the element type) |
| AC11 item 5 | det fault, theta = 2.665e-15: raise on both backends / quiet on the GPU with det_max within theta of the CPU | 1/1, 1/1 | all | PASS |
| AC11 item 6 | standard-map panel K in {0.5, 0.9716354, 2, 5, 10}, 50 starts x 1e4 steps: health errors; max \|det_max GPU - CPU\| / det threshold in force | 0; 0 | 0; 1 | PASS |

Item 5, raise side: a determinant fault on the GPU is recorded by the monitor in its state, the host sees a non-OK status at the end of the chunk, the chunk reruns on the CPU (reason 5) and the CPU raises the error that reaches the caller, so that error carries `backend = cpu` and the run returns no provenance; the verdict (raise on both sides of the threshold) matches the CPU. The quiet side is checked from the GPU provenance (no fallback) and the det maxima. Item 6: the det maxima on the GPU equal those on the CPU for every start (ratio 0 to the threshold). The monitor state compiles for the GPU for Float64; Double64 does not (cause in the tangent-core GPU section).
