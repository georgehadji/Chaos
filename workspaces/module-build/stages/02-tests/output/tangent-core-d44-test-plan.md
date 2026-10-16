# tangent-core, Amendment D44 (step-chunked relaunch): stage 02 test plan

Spec: `01-spec/output/tangent-core-spec.md`, "Amendment D44" (items 1 to 15, criteria D44a to D44k). Stage 02 RED run, 2026-10-16.
Interface fixed by the orchestrator: `Chaoslab.Tangent._LAUNCH_STEPS` (Ref{Union{Nothing, Int}}), `Tangent._before_launch(device, chunk_index, launch, carried)`,
`Tangent._budget_steps(ns, tau)`, provenance `launch_steps` and `launches_per_chunk`, `FALLBACK_REASONS[9] == "device lost"`, `propagate_shadow` on `Tangent.DEFAULT_GPU_NS`.

Files: `test/tangent-core/tc_d44.jl` (new), `test/lyap-two-trajectory/l2_d44.jl` (new), registered in the two `runtests.jl`; helper `tc44_budget_L` in `tc_helpers.jl`.
Seeds: tangent-core `TC_SEED + 4400 + sub` (sub 1 to 8, one per testset), lyap-two-trajectory `l2_seed(18, sub)` (sub 1 to 5, inside the slug's range 25000 to 25999). Every random draw is a `StableRNG(seed)`.

## Criterion to test

| Criterion | Test (file:testset) |
|---|---|
| Interface (items 1, 2 names; L not a positive multiple of tau throws ArgumentError, CPU and device path) | `tc_d44.jl:D44 interface`, `l2_d44.jl:D44k(a)` (last loop) |
| D44a CPU bitwise, carried state (k = 1 Benettin, k = 2 QRSpectrum(RGS()), k = 2 SVDSpectrum; tau 1, 2, 10; L = tau, 20 tau, 130 tau; N = 2000; x, W, finals `isequal`) | `tc_d44.jl:D44a CPU bitwise, carried state` |
| D44b health error step (Inf at step 12, tau 5, L 5: step 15, orbit id, backend, kind "state"; also NaN in W; reducer failure in launch 2 reports the single-launch step; device path reruns on the CPU, same step) | `tc_d44.jl:D44b health error step` |
| D44c, CPU-backend device (FakeDevice, forced L < N equals L = N and the unforced run bitwise) | `tc_d44.jl:D44c device path on the CPU backend, bitwise` |
| D44c [local] GTX 1050 part | stage 04, `verify/gpu.jl` |
| D44d budget (exact `_budget_steps` values, estimate from provenance at most 1 s, N = 1e6 gives 26 launches, even L = tau too long gives the CPU with the budget reason) | `tc_d44.jl:D44d budget`; changed: `tc_device.jl:D29 / launch budget (display watchdog)` |
| D44e lost device in launch 3 (first, middle, last chunk; result equals CPU bitwise; earlier chunks counted; reason "device lost"; DEVICE reset; carried-state corruption at launch 2 never reaches the result; lost device in setup) | `tc_d44.jl:D44e lost device in launch j > 1`; changed: `tc_device.jl:D29 / lost device` |
| D44f OOM (reason 4), compile error (reason 3) and a device-only health error (reason 5, flaky map armed by the hook) in launch 3, first, middle, last chunk | `tc_d44.jl:D44f other classified faults in launch j > 1` |
| D44g provenance (single launch N and 1; split L and cld(N, L); 0 and 0 for a CPU request, setup fallbacks, budget fallback; device values kept after a mid-run fallback; forced CPU step chunk keeps backend :cpu) | `tc_d44.jl:D44g provenance`; also `D44d budget`, `D44e`, `tc_device.jl:D29 / launch budget` |
| D44h ninth reason (9 entries, 9th "device lost", lost gives it, budget gives only the budget string; setup loss) | `tc_d44.jl:D44h ninth reason`, `D44e`; changed list tests (below). The extension code table (700, 702, ...) is stage 04. |
| D44i inspection (one stepping body) | stage 04, `verify/gpu.jl` and the 04-verify report |
| D44j [local] payoff at the original n | stage 04, `verify/gpu.jl` |
| D44k(a) two-trajectory CPU bitwise (K = 0.5, 2 pi, 256 orbits, d0 = 1e-8, tau 1, 2, 10, L = tau, 20 tau, 130 tau; x, W, sum, comp, d_prev, X2) | `l2_d44.jl:D44k(a) CPU bitwise, carried state` |
| D44k(b) error step is the global step (separation collapse and Inf at step 12, tau 5, L 5 to 20; device path reruns on the CPU) | `l2_d44.jl:D44k(b) error step is the global step` |
| D44k(c) [local] two-trajectory device bitwise on the GTX 1050 | stage 04, `verify/gpu.jl` |
| D44k(c) on the CPU backend (FakeDevice, forced L < N equals L = N and the CPU run) | `l2_d44.jl:D44k device path on the CPU backend, bitwise` |
| D44k(d) shared default ns (no record: L = 9765 for tau = 1, 9760 for tau = 10, N = 20000 splits into 3 launches; record 25 ns; budget fallback 0 and 0; forced CPU step chunk) | `l2_d44.jl:D44k(d) budget and the shared default ns value` |
| Item 15 lost device, OOM, provenance for the estimator (carried `ends` corrupted at launch 2, fault at launch 3, all three chunks) | `l2_d44.jl:D44k lost device and classified faults in launch j > 1` |

Stage 04, verify/gpu.jl (local, GTX 1050): D44c hardware part, D44i (inspection), D44j, D44k(c) hardware part, the extension code table of D44h.

## Notes on choices

- D44a uses 256 seeded orbits, not 1024 (the spec allows a reduction; at least 64 are kept, asserted). The full 1024 belongs in the stage 04 `verify/` script. The package suite is compile-bound (D38); the 256-orbit form runs the 36 bitwise comparisons in a few seconds each once compiled.
- K = 2 pi with tau = 10 makes the repeated Gram-Schmidt reducer raise "orthonormalisation breakdown" (orthogonality 0.73 above the 0.5 tolerance, found at orbit 1 step 1030 of the 1024-orbit draw) on some orbits in the unforced run. The test therefore drops the orbits whose single-launch `propagate_orbit` raises a health error (`tc44_clean`) and passes the seeded frames as `W0`, so the orbit ids keep their frames. D44c uses the same filter. This is a property of the existing reducer, not of D44.
- The hook devices (`tc44_HookDevice`, `l2_HookDevice`) overload `_before_launch` under `@static if isdefined(...)`, so a missing hook shows up as failing assertions and as the explicit `isdefined` test, not as a load error of the file.
- `carried` is mutated in place by the corruption test (`carried.x[1]`, `carried.W[1]`; `carried.ends[1]` for the estimator). The tests need these to be the live vectors the next launch reads; if the implementation hands out copies, the corruption test still passes (nothing to leak) but proves nothing; the hook plan says plain Vectors on the CPU backend.
- "Chunk stores nothing" (D44b) is observed as: the call throws and returns no result; the function-level `out` is not reachable from a test.
- `chunk_size` expectation of the existing budget case is kept at 1100 (L = N fits, so the cap is the D33 cap at L = N). See "Interface notes" in the report.

## RED run

Command: no suite selection mechanism exists in `test/runtests.jl` (it includes every `test/<module>/runtests.jl`), so the affected sets were run one by one with a wrapper that does `using Test, Chaoslab; @testset ... include("test/<set>/runtests.jl")` (scratchpad scripts, not committed): tangent-core, lyap-two-trajectory, precision. The full `Pkg.test()` was not run (about 26 min, D38).

| Set | Pass | Fail | Error | Total |
|---|---|---|---|---|
| tangent-core | 499 | 26 | 46 | 571 |
| lyap-two-trajectory | 312 | 7 | 13 | 332 |
| precision | 149 | 1 | 0 | 150 |

All Fail and Error lines belong to the new or changed tests; every other testset of the three suites passes unchanged (AC1 to AC11, W9, D34 and the older D29 testsets).

tangent-core, by testset:

| Testset | Result | Cause |
|---|---|---|
| D44 interface | 11 Fail, 3 Error | `_LAUNCH_STEPS not defined in Chaoslab.Tangent`; `isdefined(:_before_launch)` false; `@test_throws ArgumentError` got UndefVarError |
| D44a | 2 Pass (orbit filter, backend), 1 Error | `_LAUNCH_STEPS not defined` at the first forced run |
| D44b | 1 Pass (single-launch step 15), 1 Error | `_LAUNCH_STEPS not defined` |
| D44c | 1 Error | `_LAUNCH_STEPS not defined` |
| D44d | 7 Pass, 6 Fail, 30 Error | `UndefVarError: _budget_steps not defined in Chaoslab.Tangent` (the 19 budget-function assertions); the device runs with N beyond one launch go to the CPU in the current code (`backend === :gpu` Fail) and the provenance has no `launch_steps` (`FieldError: type NamedTuple has no field launch_steps`) |
| D44e, D44f | 1 Error each | `_LAUNCH_STEPS not defined` (the first forced run) |
| D44g | 3 Error | two `FieldError: type NamedTuple has no field launch_steps` (CPU request, device single launch); the forced-L cases stop at `UndefVarError: _LAUNCH_STEPS not defined` |
| D44h | 2 Pass, 3 Fail, 1 Error | `length(FALLBACK_REASONS) == 9` (it is 8), `unique == 9`, the lost-device hook never fires so `DEVICE[]` is not reset; `_LAUNCH_STEPS not defined` |
| D29 / defaults and hooks | 1 Fail | 8 reasons instead of 9 |
| D29 / launch budget | 3 Fail, 4 Error | the budget(500) and the no-record N = 10_000 cases still fall back to :cpu; `launch_steps` key missing |
| D29 / lost device | 1 Fail | reason is still the budget string |
| AC10 (tc_ac09_11.jl) | 1 Fail | reason list has 8 entries |

lyap-two-trajectory:

| Testset | Result | Cause |
|---|---|---|
| AC10 health and recording | 3 Fail | provenance lacks `launch_steps`, `launches_per_chunk` (extended haskey list, 3 element types) |
| D44k(a), (b), device bitwise, lost/faults | 1 Error each (plus 1 Pass in (a), (b)) | `_LAUNCH_STEPS not defined in Chaoslab.Tangent` |
| D44k(d) | 4 Fail, 9 Error | the estimator refuses N = 20000 (own constant 100.0, no split): `backend :cpu`; `type NamedTuple has no field launch_steps`; `_LAUNCH_STEPS not defined` |

precision AC7: 1 Fail (nine reasons).

No parse errors, no errors inside helpers other than the named missing symbols (`tc_with_launch_steps` / `l2_with_launch_steps` touch `_LAUNCH_STEPS` by design). D44d first ran with unreduced W (K = 5 overflows within about 1000 steps at N = 156251); it now uses `tc_LogNorm()`. Because a testset aborts at the first `UndefVarError`, the assertions behind it (hook-driven cases, bitwise comparisons) are not exercised until stage 03; their single-launch halves (the error steps of D44b/D44k(b), the ZeroTan reducer failure at step 10, the K = 2 pi, tau = 10 reference runs of the estimator) were run separately on the current code and behave as written.

## Verified tests changed, owner review

| File:line | Old | New |
|---|---|---|
| `test/tangent-core/tc_device.jl:93` | `length(FALLBACK_REASONS) == 8 && [8] == "launch exceeds the display watchdog"` | 9 entries, `[8]` unchanged, `[9] == "device lost"` |
| `tc_device.jl` "launch budget", budget cap case | chunk_size == 1100 only | unchanged, plus `launch_steps == N`, `launches_per_chunk == 1` (one launch fits) |
| `tc_device.jl` "launch budget", `budget_ns(500)` case | CPU, reason "launch exceeds the display watchdog", `gpu_chunks == 0` | device split: backend :gpu, reason "", `gpu_chunks == 3`, `launch_steps == L`, `launches_per_chunk == cld(N, L)` with L from the independent formula `tc44_budget_L` (48 for tau 4), result equals the CPU reference |
| `tc_device.jl` "launch budget", no-record N = 10_000, tau = 10 | CPU, budget reason | device: result equals the CPU run, backend :gpu, reason "", `launch_steps == 9760`, `launches_per_chunk == 2`; also N = 100 gives N and 1 |
| `tc_device.jl:191` (lost device) | reason "launch exceeds the display watchdog" | "device lost" |
| `test/precision/runtests.jl:163` | eight reasons | nine, "device lost" last |
| `test/tangent-core/tc_ac09_11.jl:34` | eight reasons | nine, "device lost" last |
| `verify/tangent-core.jl:709` | eight reasons | nine, "device lost" last (not run by the package suite) |
| `test/lyap-two-trajectory/l2_runtime.jl:42` | haskey list of 8 keys | adds `:launch_steps`, `:launches_per_chunk` |

Searched and left unchanged: no test asserts an exact provenance key set (all use `haskey` lists or `keys(base)`); `test/lyap-householder-qr/hh_backend.jl:49` and `verify/lyap-householder-qr.jl:654` carry the same eight-key `haskey` list and still hold (superset not required). `verify/gpu.jl:561` (comment and FALLBACK line of AC12a at n = 1e5, which D44j removes) is a stage 04 item. No other hit for "display watchdog", "FALLBACK_REASONS" or `== 8`. Docs that still say eight reasons (`src/Foundation/CONTEXT.md`, `_config/architecture.md`, the `FALLBACK_REASONS` docstring) are stage 03 and stage 05 items.

## Stage 02 audit

| Check | Result |
|---|---|
| Coverage of spec | PASS. Every criterion D44a to D44k maps to a package test or, where marked [local] or inspection (D44c hardware, D44i, D44j, D44k(c) hardware, the extension code table), to a named stage 04 item in `verify/gpu.jl`. |
| RED | PASS. Every new or changed assertion fails or errors for missing functionality (UndefVarError for `_LAUNCH_STEPS`, `_budget_steps`; missing provenance keys; 8 vs 9 reasons; wrong reason string; :cpu where a relaunch is expected). No parse errors; no error inside a helper other than those names. |
| Reproducible | PASS. All draws use `StableRNG` seeds `TC_SEED + 44xx` and `l2_seed(18, sub)`; the other inputs are fixed grids. |
| No memory values | PASS. No regression values from memory: the numbers in D44d (39062, 39060, 9765, 9760, 97656.25) are the spec's own, cross-checked by an independent formula in the tests (`tc44_budget_L`); no bibliography benchmark is used. |
