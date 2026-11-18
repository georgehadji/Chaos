# Implementation notes: continuous-qr (catalog row L8; spec D120)

Code: `src/Lyapunov/continuous_qr.jl` (182 lines), include and exports in `src/Lyapunov/Lyapunov.jl`, `src/Lyapunov/CONTEXT.md` (75 lines), `docs/src/api/lyapunov.md` (one phrase; the autodocs block picks up the docstrings), `bench/continuous-qr.jl`. Spec: `../01-spec/output/continuous-qr-spec.md`. Test plan: `../02-tests/output/continuous-qr-test-plan.md` (its API proposal, U1 to U13 and the E_max default 1e-2 of U6 are the operative reading). Tests are untouched.

## Status

| Run | Result |
|-----|--------|
| New suite `test/continuous-qr/` (driver `red_cq.jl`, scratch environment `env05`, Julia 1.13.2, 1 thread) | 1615 pass, 0 fail, 0 error, 82 s (run twice, once before and once after the Aqua fix below; same totals) |
| Aqua, JET (`test/quality.jl`), layering (`test/layering.jl`), one process | first run: 1 failure (Aqua unbound type parameter in `_totals`, fixed); rerun 14 pass, 0 fail, 2 min 31 s |
| Regression, one process: `lyap-benettin`, `lyap-gram-schmidt`, `lyap-householder-qr`, `lyap-svd-stabilized`, `lyap-two-trajectory`, `tangent-core`, `integrator-catalog`, `henon-heiles` | 66754 pass, 0 fail, 0 error, 2 broken (the pre-existing broken marks), 9 min 58 s |
| Coverage, `--code-coverage=@src/Lyapunov/continuous_qr.jl`, new suite alone | 85 of 85 executable lines, 100 percent; the `.cov` file was deleted |
| Whole package suite `Pkg.test()` (main session) | 1428377 pass, 0 fail, 0 error, 3 broken (integrator-catalog 2, weighted-birkhoff 1, as before), 37 min 47 s |
| Not run | `verify/` scripts, `CQ_FULL_RUN` sizes, device runs, `--check-bounds=yes` (stage 04) |

## Design

- One file, 12 functions, longest `_cqr_run` (18 lines). `ContinuousQR{M, F}` holds `map`, `flow`, `order`, `penalty_factor`, `orthogonality_tol`, `form`. It is a per-orbit loop, not a reducer; the orbit is advanced by `forward(map, z)` only and the Jacobian is `vector_field_jacobian(flow, z)`, so it works for any map with a `tau` field (the provenance additionally reads `map.method`).
- `cqr_rhs(Q, A, nu; form)` is the pure penalty form: form 84 `Q S - nu Q E`, form 92 `A Q - Q U + Q S (I + E) - nu Q E`, `U = Q'AQ`, `S` from the strict lower triangle of `U` (`_skew_lower`), `E = Q'Q - I`; `nu = c_nu ||A||_F` with `c_nu` 1 (`k = d`) or 2 (`k < d`). `cqr_growth(Q, A)` gives `(Q'AQ)_ii`.
- Order 2: Heun on one orbit step (2 right-hand sides, 1 new Jacobian). Order 4: RK4 on two orbit steps with Jacobian nodes A_n, A_{n+1}, A_{n+1}, A_{n+2} (4 right-hand sides, 2 new Jacobians). Total per run: N + 1 Jacobian calls, 2N right-hand sides.
- Exponents: the increment `(h/6) sum w_i u_i` of each Q-step goes into one `CompensatedSum` per column (tuple of accumulators, isbits); `l` is the total, `lambda = l / t`.
- Guard: after each Q-step `E = max |Q'Q - I|`; non-finite or above `orthogonality_tol` (default `CQ_EMAX = 1e-2`, U6) throws `NumericalHealthError("orthogonality", ...)` with the step. Before the loop: "state" (non-finite x0) and "deviation matrix" (non-finite Q0 or Jacobian, step 0; a Jacobian at orbit step j reports step j, U4). All argument errors (`order`, `form`, `k`, odd N for order 4, N <= 0, `:f84` with `k != d`, `every < 1`) are `ArgumentError` before any work.
- Series: `_record!` dispatches on the sink (`nothing` or a vector), so `propagate_cqr_orbit` has no `Union` in its path and allocates nothing; `propagate_cqr_series` shares `_cqr_run`, so the series is a bitwise prefix of the run by construction.
- `propagate_cqr`: serial loop over orbits, `finals[j]` is the orbit `final`; provenance as listed in the plan. `resolve_backend(backend, T; cpu_only_reducer = true)` gives the backend record, so every request resolves to the CPU (`:gpu` and `:auto` carry the reason "reducer is CPU-only"); `chunk` is accepted and unused.
- Layering: Lyapunov now imports `..Integrators: method_name, nominal_order` (layer 3 below layer 5, allowed) and `..Tangent: initial_deviations, PropagationResult`.

## Deviations from the spec and the test plan

| Item | Deviation |
|------|-----------|
| Exports | `CQ_EMAX` is exported beside the seven names of the plan (the test suite reads it as `CQ_EMAX`; the helper defines its own constant of that name, and both hold 1e-2) |
| `propagate_cqr_series`, `propagate_cqr_orbit` | take an extra keyword `backend` (series) for the error record; the plan lists `every, orbit` only |
| GPU | no device path (spec criterion 22 device part is stage 04); a `:gpu` request falls back to the CPU with the reason of `resolve_backend` for a CPU-only reducer, not a kernel-specific reason |
| `chunk` | accepted, ignored (the CPU loop has no chunks; chunk invariance of AC17a is trivial) |
| Arithmetic | the expressions follow the scratch calibration implementation of stage 02 operation by operation, so the calibrated bounds carry over; no bound was changed |

No test was edited and no bound loosened.

## Calibration margins (realized on this implementation, h = 1/20, t = 60 unless said)

Realized values equal those of the stage 02 scratch run (same arithmetic); they were read again from the log of this run for AC3, AC5, AC6a, AC6b, AC8, AC9, AC12, AC13a.

| Test | Bound | Realized | Margin |
|------|-------|----------|--------|
| AC3 worst error / B_lin | 1 | 0.017 | 59x |
| AC5 worst error / bound | 1 | 0.138 | 7.2x |
| AC6a worst error / bound | 1 | 0.389 | 2.6x |
| AC6b column independence | 5e-3 | 1.3e-3 | 3.8x |
| AC7 error (q = 4) pend 2 / pend 3 / F1 / F3 R / F3 C | 5e-3, 1e-1, 5e-3 | 8.7e-4, 2.3e-2, 7.5e-5, 1.2e-4, 4e-4 | 5.7x, 4.3x, 67x, 42x, 12x |
| AC7 error (q = 2) pend 2 / F1 / F3 R / F3 C | 1.0 | 0.235, 0.012, 0.205, 0.338 | 4.3x, 83x, 4.9x, 3.0x |
| AC8 worst error | 3e-2 | 1.55e-2 | **1.9x (under 2x)** |
| AC9 worst lambda error / allowance | 1 | 0.469 | **2.1x (just above 2x)** |
| AC9 worst l error at the long horizon | 5e4 h^4 | 1.70e4 h^4 | 2.9x |
| AC12 worst \|sum\| / B_sum | 1 | 3.4e-5 | 30000x |
| AC13a worst E / bound | 1 | 1.0e-3 (bound 1e-2 for q = 4) | 10x |
| AC14 worst error / B_C | 1 | 3.8e-4 (Float64), 7.1e-6 (Double64) | 2600x |
| AC20 worst error / B | 1 | 0.14 | 7x |
| Order ratios h to h/2 | [2^q/2, 2 x 2^q] | 3.85 to 4.7 (q = 2), 8.6 to 20 (q = 4) | inside; AC7 pend 2 q = 4 at t = 10 gave 16.8 |

Flags: AC8 (1.9x) is the thinnest bound, and the plan already records it as U8 (discretization allowance 3e-2 against the spec's 1e-3); AC9 at 2.1x is borderline and the T = 1000 allowance is not calibrated (stage 04). R5 (no penalty form 92) 173.5 time units to E = 1e-2 against Inf with the penalty, as in the plan.

## Performance (`bench/continuous-qr.jl`, benchmark environment, 1 thread, medians; indicative, BenchmarkTools resolution 100 ns)

| Entry | Time | Allocations |
|-------|------|-------------|
| `cqr_rhs`, form 84, d = k = 4 | 0.1 us | 0 |
| `cqr_growth` | 0.1 us | 0 |
| one Q-step, order 4 (2 orbit steps) / order 2 (1 orbit step), with the orbit steps and Jacobians | 1.8 / 1.2 us | 0 |
| orbit C, 100 time units (2000 orbit steps, KahanLi8, h = 1/20), order 4 / order 2 | 1.69 / 2.14 ms | 0 |
| discrete RGS spectrum, same orbit, same length | 4.49 ms | 0 |
| 256 orbits, 100 time units, k = 4, order 4, CPU | 441 ms | 1071 (the result vectors and the frames) |

Cost ratio to the discrete RGS spectrum (R1): 0.38 for order 4 and 0.48 for order 2 in this benchmark (the test suite's own timing reported 0.30 and 0.63; both are noisy, the order-4 variant is cheaper because it takes half as many Q-steps).

## Known limits

- CPU only, serial over orbits; no kernel path. `chunk` is ignored.
- The guard default 1e-2 is not tight for order 2 (E_max about 8e-3 at h = 1/20 on F3, plan U6); the tests pass `5e-2` for q = 2 through `orthogonality_tol`.
- Order 2 on the pendulum of amplitude 3 has an error of 3.4 in the direction identity (AC7 skips it); use order 4 for large amplitudes.
- The lambda form of the long-horizon comparison with the discrete variants fails on orbit C at T = 200 (differences grow in events); the l form is the one tested (plan U9).
- The module accepts f(z0) = 0 (U1); rejection of an equilibrium start for the direction identity is the caller's frame builder.
- Float32 and BigFloat were not run beyond the precision ladder of AC15 (Float64, Double64, BigFloat(256)).

## Audit

| Check | Result |
|-------|--------|
| GREEN | new suite 1615 pass, 0 fail, 0 error; regression suites 66754 pass, 0 fail, 2 broken (pre-existing); whole package suite 1428377 pass, 0 fail, 3 broken |
| Quality | Aqua and JET clean, layering passes (14 of 14, after fixing one Aqua unbound type parameter); coverage of `continuous_qr.jl` 85 of 85 lines, 100 percent |
| Size | longest function `_cqr_run` 18 lines; `continuous_qr.jl` 182 lines, `Lyapunov.jl` 31 lines; `src/Lyapunov/CONTEXT.md` 75 lines |
| Tests untouched | `git diff 165c3ed -- test/` is empty |
| Numerics | no `@fastmath`; the exponent sums are `CompensatedSum` accumulators (one per column); the `.cov` file was deleted |

## Documentation and code inconsistencies

- The spec's default `E_max` is 1e-3 (P5); the approved test plan (U6) sets 1e-2, and `CQ_EMAX = 1e-2` is the code default. The spec text is not updated.
- The spec says the equilibrium rejection and "initial frame" belong to the module (C-items); the plan (U1) moved them to the test-side frame builder, and the code follows the plan.
- The spec's "within rounding" for the AC12 sum is not reachable (U5); the code gives the bounded sum of the plan.
