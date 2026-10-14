# Implementation notes: diagnostics

Code: `src/Diagnostics/` (`Diagnostics.jl`, `residuals.jl`, `monitor.jl`); bench `bench/diagnostics.jl` (BenchmarkTools not installed here, so the suite was not run; numbers below come from a plain timing loop).

## Status

`test/diagnostics/` in the shipped tree: 187 pass, 3 fail, 1 error (all four are orbit-id assertions, see "Open point"). With the Tangent patch below applied in a scratch copy of the tree: 195 pass, 0 fail, 0 error; `test/tangent-core/` still 423/423 in that copy.

## Open point: reducer cannot signal failure with an orbit id (needs a Tangent change)

`reducer_update` has no orbit argument and `_orbit_body` neither catches nor inspects reducer states, so a reducer can only (a) throw from inside the update (orbit unknown, and inside the kernel it surfaces as a task failure, not as a `NumericalHealthError`, so AC10 cannot be satisfied) or (b) poison `W` (then Tangent reports "deviation matrix", not "det"). The shipped monitor does (a) with `orbit = 0`. Failing tests: AC5/det (`e.orbit == 7`, 2 assertions), AC10 (error from the kernel), AC11 item 5 (`e.orbit == 1`).

Scratch fix, verified GREEN (files in the session scratchpad `m2/`: `dg_tangent_reducer_step.patch` for `src/Tangent/{Tangent,orbit,propagate}.jl`, `dg_monitor_with_step.jl` to replace `src/Diagnostics/monitor.jl`):
- Tangent gains `reducer_step` (default calls `reducer_update`; `_apply` calls it), `reducer_failed(r, state)` (default `false`) and `reducer_error(r, state, orbit, backend)`. The body stops an orbit with a new status `STATUS_REDUCER` when any state reports failure; `propagate_orbit` and `_run_chunk!` (lowest orbit id first, as for the other statuses) build the error with the true orbit id.
- The monitor overrides `reducer_step` (non-throwing: it records the first exceedance, code/value/threshold/step, in its isbits state) and keeps `reducer_update` as a thin throwing wrapper, which the tests call directly (AC5b, AC5d, AC5e).
This also matches the spec ("a kernel records the first-exceedance record per orbit, the host raises").

## Deviations from the spec

- Pairing health check for N >= 2 at finalise (needs kappa, O2) is not implemented: `monitor_exponents` checks `sum` at every call and records `eps_pair`; the pairing check is left to the caller via `check_residual("pairing", ...)`. For N = 1 pairing equals sum (I3), so only `sum` can fire (as the tests expect).
- `monitor_exponents` ignores its positional `n` (the contract has it plus `step`; errors use `step`).
- Names in the fixed order are checked det, symplecticity, orthogonality (update) and sum (exponents); only the first failing name is reported per call.
- `det_threshold` and friends return `T`; the error stores `Float64(value)` and `Float64(threshold)` (`Inf` threshold when the value is non-finite).
- eps_det and eps_symp are sampled at the handoff points (spec amendment). Orthogonality is checked on the `W` handed to the reducer, only when it has 2 or more columns.

## Performance (Float64 and Double64, standard map K = 2, 1 thread, -O2)

- `reducer_update` (2 x 2): 163 ns Float64, 4.0 us Double64 (dominated by `sin`/`cos` of the tangent), 0 allocations in both (AC9 passes).
- `propagate_orbit` k = 1: 293 ns/step with the monitor versus 155 ns/step bare.
- Suite runtime about 100 s (AC6 and AC2 in Double64 dominate).

## Known limits

- For 2N >= 4 `eps_det` and `eps_symp` come from Systems (`det_residual` builds a dense matrix for D > 3), so the 0-allocation claim is verified for 2 x 2 only.
- Coverage, Aqua and JET not run in this stage run (orchestrator runs them).

## Review fixes (stage 04)

- `monitor_exponents` raises a failure already recorded in the state (det, symplecticity, orthogonality) before the sum and pairing checks; error construction shared with `reducer_error` through `_failure`.
- `det_residual` for D > 3 no longer allocates a dense matrix (see standard-map notes); the D > 3 allocation caveat above no longer applies.

## Amendment D34 (2026-10-14)

Code: `src/Diagnostics/residuals.jl` (A1b: `det_threshold`, `symp_threshold` and `sum_threshold` use `effective_eps(T)` instead of `eps(T)`; docstrings say `eps_eff`), `src/Diagnostics/monitor.jl` (A3: `reducer_final(::ResidualMonitor, st)` adds `M_max` and `G_max`; the state already tracked both).

### Result
`test/diagnostics/`: 245 pass, 0 fail, 0 error (195 before, 50 new in `dg_d34.jl`), Julia 1.13.1, 2 threads, 2 min 10 s. No Batch A testset needed a change (none compares the final NamedTuple as a whole).

### Deviations and notes
- Float64 and BigFloat thresholds are bitwise unchanged; Double64 thresholds are 16 times looser (2^-100 instead of 2^-104), which removes the false det alarm of a Double64 chain at a residual that meets the rounding model.
- `G_max` in the final is 0 in a plain run (only `monitor_exponents` raises it); the caller forms the end-of-run sum check as `sum_threshold(pairing_tol, T, f.G_max, tau, f.M_max)` (Q1, decision D34).

### Known limits
- None new.
