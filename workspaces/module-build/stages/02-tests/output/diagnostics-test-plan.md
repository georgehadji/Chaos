# Test plan: diagnostics

Spec: `workspaces/module-build/stages/01-spec/output/diagnostics-spec.md`. Tests: `test/diagnostics/runtests.jl` (helpers prefixed `dg_`) and `test/diagnostics/faultmaps.jl` (fault map and renormalising reducer, loaded lazily inside the testsets that need them so a missing `Chaoslab.AbstractMap` errors only those). Element types: Float64 and Double64 (Float32 dropped, D24). Keys hard-coded from `_config/numerics-defaults.toml`: `det_residual_tol = 4`, `pairing_residual_tol = 16`, `seed = 20261012`. Seeds: base + criterion number (AC6 adds 600/700 + panel index per type).

## Criterion to test

| AC | Testset | Samples | Tolerances |
|----|---------|---------|-----------|
| 1 | AC1 integer maps exact | 4 integer maps x 1000 steps x {Float64, Double64}; direct functions and ResidualMonitor reducer | eps_det = eps_symp = orth = 0 exactly, no health error. Orbit wrapped mod 2 pi each step (tangent of a linear map does not depend on it; the lift overflows for the cat map near n = 740) |
| 2 | AC2 rotation 0.7311 | 10^6 steps x {Float64, Double64}; MGS reference re-orthonormalisation each step | eps_det <= det_threshold and <= 8 eps(T); eps_orth <= 8 eps(T) |
| 3 | AC3 cat map reference QR | n = 10^5, tau = 1, Householder `qr` of A*Q, `CompensatedSum` of log abs R_ii | eps_sum(n) <= sum_threshold(16, Float64, G = mu^2 = 6.854, 1, M = 2) at every n; eps_pair = eps_sum exactly; lambda_1 within 1e-3 of ln mu (estimator sanity) |
| 4 | AC4 pairing index mapping | N = 1 and N = 2 synthetic exponents x {Float64, Double64} | N = 1: eps(T)(abs L + d); N = 2: 2 eps(T) max abs lambda; adjacent pairing rejected (value < 0.5) |
| 5 | AC5 injected faults / det | A_delta = [[2 + delta, 1], [1, 1]], deltas = powers of two in [2 theta, 4 theta) (raise: name det, value = delta bitwise, orbit 7, step 1) and (theta/4, theta/2] over 10^4 steps (quiet), via `propagate_orbit` | exact |
| 5 | AC5 injected faults / orthogonality | one basis vector scaled by 1 + delta, residual 2 delta + delta^2, both sides of an absolute run tolerance (2^-30 Float64, 2^-60 Double64), via `reducer_update` | within eps(T) |
| 5 | AC5 injected faults / sum and pairing | delta added to lambda_1 (N = 1), `monitor_exponents` (name sum, orbit 9, step 3) and direct `check_residual("pairing", ...)`, both sides of sum_threshold | within eps(T)(abs L + 1) |
| 5 | AC5 injected faults / symplecticity | 4 x 4 `block_diagonal` with one block det 1 + delta, both sides of symp_threshold (N = 2); plus diag(2, 1/2, 1, 1) (det 1, not symplectic) through the monitor | within 16 eps(T); exact name |
| 5 | AC5 injected faults / non-finite | NaN and Inf, direct `check_residual` for all five names with threshold Inf; NaN/Inf in DT, in W (as Q) and in lambda through the monitor with Inf tolerances; eps_* functions never return a finite value for non-finite input | name must match, threshold field Inf |
| 6 | AC6 standard-map calibration | K in {0.5, 0.9716, 2, 5, 10}, 50 starts, 10^4 steps, Float64 and Double64 (5 x 50 x 10^4 per type) | max eps_det/(eps M) <= 1.5; no step above det_threshold |
| 7 | stage 04 | long runs 10^6 steps | not a package test |
| 8 | AC8 256-bit recomputation | 10^3 steps, K = 5, one seeded start, `setprecision(BigFloat, 256) do` on the stored entries | <= eps(Float64) |
| 9 | AC9 zero allocation | `@allocated` of `propagate_orbit` (100 steps, k = 1), `reducer_update` (k = 2), `monitor_exponents`, after warm-up, in functions; Float64 and Double64 | == 0 |
| 10 | AC10 thread invariance | 64 seeded orbits, K = 2, 100 steps: `propagate` vs serial `map` over `propagate_orbit`, bitwise (x, W, reducer finals). Fault: det 1 + 2^-20 on orbit 37 only; on orbits 41 and 17 | bitwise `===`; reported orbit 37, then 17 (lowest id), same as the serial first error |
| 11 | AC11 CPU backend items 1, 5, 6 | item 1: four integer maps x 1000 steps through `propagate(backend = :cpu)`, zero start (fixed point, lift cannot overflow), `dg_Renorm` keeps W finite; item 5: det fault through `propagate`; item 6: panel K, 50 starts, 10^4 steps (Float64), `dg_Renorm` keeps W finite | det_max = symp_max = 0; fault verdicts as AC5a; det_max <= 1.5 eps (1 + K) |
| 11 | GPU part | locally on the GTX 1050, 04-verify report | stage 04 |

Testset count: 10 criterion testsets (AC1, AC2, AC3, AC4, AC5, AC6, AC8, AC9, AC10, AC11); AC5 nests five fault testsets (det, orthogonality, sum and pairing, symplecticity, non-finite), 14 error units in the RED run. Deferred to stage 04: AC7, GPU part of AC11.

Estimated runtime with the real implementation: about 55 to 60 s on CPU, dominated by AC6 in Double64 (about 5 s per K, measured with a throwaway stub) and AC2 in Double64 (about 5 s). AC11 item 6 (2.5e6 kernel steps in Float64) is unmeasured, since the stub kernel allocates per step. No sample size was reduced.

## Spec ambiguities resolved

- `ResidualMonitor.orth_tol` is treated as an absolute threshold on eps_orth (spec: "threshold supplied by the run"); `det_tol` is a multiplier as in `det_threshold`. AC5b depends on this.
- Exponent vectors passed to `eps_sum`, `eps_pair`, `monitor_exponents` are `SVector{2N,T}`.
- `propagate` is assumed not to renormalise (tangent-core: only reducers replace W), so tests that run long use the test reducer `dg_Renorm`; ResidualMonitor orthogonality is exercised only through direct `reducer_update` calls (k >= 2), never through unnormalised propagation.
- For N = 1 the monitor reports "sum" before "pairing" (fixed order); the "pairing" name is tested through `check_residual` directly, since the API has no kappa and no N >= 2 finalise check.
- The injected-fault orbit id comes from the `orbit` keyword of `propagate_orbit` and from `propagate`; `reducer_update` itself has no orbit argument.

## RED log

Command (repo root): `julia --project=<scratchpad>\m2\testenv -e "using Test, Chaoslab; include(raw\"test/diagnostics/runtests.jl\")"`. Result: 0 passed, 0 failed, 14 errored (the 9 flat testsets plus the 5 nested AC5 ones). Every testset Errors on the first missing Chaoslab name.

| Testset | Error | Fail | First error |
|---------|-------|------|-------------|
| AC1 integer maps exact | 1 | 0 | UndefVarError `LinearSymplecticMap` |
| AC2 rotation 0.7311 | 1 | 0 | UndefVarError `LinearSymplecticMap` |
| AC3 cat map reference QR | 1 | 0 | UndefVarError `sum_threshold` |
| AC4 pairing index mapping | 1 | 0 | UndefVarError `eps_pair` |
| AC5 / det | 1 | 0 | UndefVarError `AbstractMap` in `Chaoslab` (loading faultmaps.jl) |
| AC5 / orthogonality | 1 | 0 | UndefVarError `AbstractMap` in `Chaoslab` |
| AC5 / sum and pairing | 1 | 0 | UndefVarError `sum_threshold` |
| AC5 / symplecticity | 1 | 0 | UndefVarError `AbstractMap` in `Chaoslab` |
| AC5 / non-finite | 1 | 0 | UndefVarError `AbstractMap` in `Chaoslab` |
| AC6 standard-map calibration | 1 | 0 | UndefVarError `StandardMap` |
| AC8 256-bit recomputation | 1 | 0 | UndefVarError `StandardMap` |
| AC9 zero allocation | 1 | 0 | UndefVarError `StandardMap` |
| AC10 thread invariance | 1 | 0 | UndefVarError `StandardMap` |
| AC11 CPU backend items 1, 5, 6 | 1 | 0 | UndefVarError `AbstractMap` in `Chaoslab` |

Helper check: the suite was also run against a throwaway stub implementation of the whole API kept outside the repo (scratchpad `m2/dg_fake`, never committed). All pass except the four AC9 allocation checks, which fail only because the stub is not isbits; the helpers, thresholds and window selection are therefore consistent.
