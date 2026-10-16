# Implementation notes: lyap-two-trajectory (M2 Batch B, stage 03)

Branch m2-p1-modules, base HEAD 9382616. Not committed.

## Files

- `src/Lyapunov/separation.jl` (36 lines): `separation`, `displaced` (default methods and `WrappedStandardMap` methods, spec T2, T3), unexported `state_mode`.
- `src/Lyapunov/two_trajectory.jl` (about 215 lines): `propagate_shadow`, `propagate_shadow_orbit`, kernel `_shadow_kernel!`, end record `ShadowEnd`.
- `src/Lyapunov/Lyapunov.jl`: two `include` lines and four export names (`propagate_shadow`, `propagate_shadow_orbit`, `separation`, `displaced`) only.
- `bench/lyap-two-trajectory.jl`: four entries.

## Result

`test/lyap-two-trajectory/` is GREEN: 568 pass, 0 fail, 0 error (Julia 1.13.1, 2 threads, suite time 127 s compile included, budget 3 minutes, D35). No cuts from the plan were needed. Tests untouched.

## Design (as the contract, D34)

- Own per-orbit loop; chunk, kernel launch, D33 launch budget and CPU fallback duplicated from `Tangent._propagate` (Q7). Imported unexported Tangent names: `device_backend`, `workgroup`, `_fallback_reason`, `CHUNK_HEADROOM`, `LAUNCH_WAVE`, `_budget_orbits`, `_check_steps`, `_host`. No Tangent change.
- The body is non-throwing and returns an isbits `ShadowEnd` (final x, shadow y, compensated sum, d_prev, code, step, value); host throws after the chunk and stores nothing for a failed chunk (lowest orbit id first). Codes: 1 "state" (non-finite component in either orbit), 2 "separation" (d zero or non-finite, threshold 0.0 or Inf), 3 "separation" (realized d^+ outside [d0/2, 2 d0], threshold d0).
- The start check (Algorithm 1) runs on the host for all orbits before any chunk (step 0, backend `:cpu`); the body repeats the same computation, so `propagate_shadow_orbit` also rejects at step 0.
- `chunk_size(...; s = 2s, reducer_bytes = 3 sizeof(T))` reproduces the spec formula `sizeof(T)(2s + 3) + 16`; the real `ShadowEnd` is larger (failure record), Q12.
- `U0` is normalised with `safe_norm`; the seeded direction is used as is (already unit), so `propagate_shadow` equals a serial `map` of `propagate_shadow_orbit` bitwise.
- `d0` must widen exactly into T (else `ArgumentError`), finite and positive.
- Wrapped `separation` is the minimal image into (-pi, pi] of the wrapped difference; `displaced` wraps with `Systems._reduce` and adds the counts to those of the reference (Q10).
- Provenance: Tangent's keys with `k = 0`, plus `estimator`, `d0::Float64`, `state_mode`.

## Deviations from the contract

None in behaviour. The only addition: `state_mode` is a small unexported function in Lyapunov (dispatch on the system) to fill the provenance key.

## Performance (loaded machine, other agents running; Julia 1.13.1, one thread, Float64, K = 5)

- `propagate_shadow_orbit`, 1e5 steps: about 730 ns per step wrapped (two `forward` steps plus separation and displacement per handoff at tau = 1), about 665 ns lifted. 0 bytes allocated after warm-up (AC11, tested for Float64 and Double64, wrapped and lifted).
- Suite time with 2 threads, compile included: 127 s on Julia 1.13.1 (machine quiet), 194 s on Julia 1.10.12 while two other full suites ran in parallel (budget 3 minutes, D35). No plan cut applied; if stage 04 times it on a quiet 1.10 and it exceeds 3 minutes, the plan's cuts apply in their listed order.

## Other suites (Julia 1.13.1 and 1.10.12, 2 threads, SKIP standard-map)

All pass on both versions: quality.jl (Aqua + JET, 12/12), layering.jl (2/2), lyap-benettin (15009), lyap-gram-schmidt (19555), tangent-core (565), precision (150), diagnostics (245). The runner script printed "FAILED MethodError" on 1.13 for each; that is its `Test.get_test_counts` call, not a test failure (each testset summary shows all passed). The lyap-householder-qr and lyap-svd-stabilized suites belong to the other agent and were not run here. Aqua and JET ran over the whole package including my files (uncommitted state of the other agent's files at the time).

## Known limits

- GPU path untested here (no CUDA in this stage): isbits `ShadowEnd` and the Double64 `log` and `sqrt` must compile for the device, otherwise the chunk falls back with reason 3 as in Tangent.
- `separation`/`displaced` exist for `WrappedStandardMap` and the plain default only; a new system with extra state components must define both methods.
- The start check and the body both compute `_shadow_start` (one extra separation per orbit, negligible).
- Series values X2(n) at several n come from re-running with different N (prefix invariance, Q9); there is no emission channel.

## Amendment D44 (2026-10-16)

Code (`src/Lyapunov/two_trajectory.jl`):
- `_shadow_body` takes the carried end record and `(n0, steps, N)`. At `n0 = 0` it starts the shadow inside the body, as before; later launches continue from the carried reference, shadow, sum and previous separation.
- `_shadow_kernel!` works in place on the `ShadowEnd` buffer and skips orbits that failed earlier.
- `_shadow_chunk!` loops over launches and calls `Tangent._before_launch(device, chunk, launch, (ends = buf,))`.
- The budget, the default ns value, the forced launch length, lost device (reason 9) and the provenance fields come from Tangent (`_device_plan`, `_forced_launch`, `_launches`). The literal 100.0 is gone.

### Result
- `test/lyap-two-trajectory/`: 477 pass (D44k a, b, d and the device path on the CPU backend).
- Full suite green, as noted in tangent-core-impl-notes.md, Amendment D44.

### Notes
- `propagate_shadow_orbit` passes a placeholder end record (`_blank_end`), which the body never reads at `n0 = 0`.
- D44k(c) on the GTX 1050 is stage 04.
