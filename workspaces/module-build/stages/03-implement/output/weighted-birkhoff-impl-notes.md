# Implementation notes: weighted-birkhoff (stage 03)

Spec: `workspaces/module-build/stages/01-spec/output/weighted-birkhoff-spec.md` (with Amendment D84). Tests: `test/weighted-birkhoff/` (untouched since `29dc943`; `git diff 29dc943 -- test` is empty). Plan: `workspaces/module-build/stages/02-tests/output/weighted-birkhoff-test-plan.md`. The API is the plan's proposal, unchanged.

## What was built

| File | Lines | Holds |
|------|-------|-------|
| `src/Indicators/wbirkhoff.jl` | 157 | `CosObservable`, `SinObservable`, `LiftedObservable`, `observable_name`; the shared update `_wb_sample` on the cell `WBCell`; `weighted_birkhoff`; `birkhoff_index` |
| `src/Indicators/wbirkhoff_reducer.jl` | 143 | `WeightedBirkhoff`, state `WBState` (`WBHorizon`, `WBCell`), the protocol methods, `reducer_final`, the `check_run` method, `reducer_provenance` |
| `src/Indicators/Indicators.jl` | +3 | seven exports (`weighted_birkhoff`, `birkhoff_index`, `WeightedBirkhoff`, `CosObservable`, `SinObservable`, `LiftedObservable`, `observable_name`), two includes after `naff_reducer.jl`, module docstring; `src/Chaoslab.jl` re-exports unchanged |
| `src/precompile.jl` | +4 | one K = 2, H = 2 `propagate` and one `weighted_birkhoff` with `birkhoff_index` per element type |
| `bench/weighted-birkhoff.jl` | 44 | BenchmarkTools suite with the budget in its header |
| `docs/src/api/indicators.md`, `src/Indicators/CONTEXT.md` (79 lines), `CHANGELOG.md` | | docs sentence, context, entry (catalog status not flipped: stage 05) |
| `_config/` | 0 | no keys (the spec declares none: floor multiplier 32 and minimum window 16 are constants, echoed in the provenance) |

No change in Tangent, Systems, Foundation or `ext/ChaoslabCUDAExt.jl` (the extension is generic over reducers; the reader needed nothing). E2 (`naff_reducer.jl`) is not edited. Reused as they are: `_wb` (the mirrored-index weight), `birkhoff_weight`, `CompensatedSum`, `add`, `total`, `effective_eps`, `NF_E2_MIN_WINDOW` (aliased `WB_MIN_WINDOW`).

Design: the pure function and the reducer share one update (`_wb_sample` on a `WBCell` of two Neumaier sums and S), and the weight of sample `j` is `_wb(jj, N, p, T)` in both, so the reducer is bitwise the function on the sampled series by construction (AC6). The state is `NTuple{H, WBHorizon{T, K}}` (B once per horizon, a cell per observable) plus `n::Int`, `handoffs::Int32`, `code::Int16`. Nothing reachable from `reducer_init` or `reducer_step` throws or allocates: the observables are mapped over the tuple, a non-finite value sets `code` to its index and returns the state otherwise unchanged; all argument checks are in the constructor and `check_run`. `birkhoff_index` (with its `ArgumentError`) is called only from `reducer_final`, which is host code.

## Results

- Folder `test/weighted-birkhoff/`, `julia --project=. --check-bounds=yes -t 2`, default sizes: **611 passed, 0 failed, 0 errors, 1 broken** (30 testsets; the broken one is the skipped device testset AC13), 24 to 30 s warm. First run after the code was written, no iteration needed on the tests.
- With `WB_FULL_RUN = true` (spec sizes): the criteria testsets AC1 to AC12 and the API pass (**594 passed, 0 failed, 1 broken**), but **5 errors** in the recorded-line testsets R1, R3, R6, R7, R8: `BoundsError: attempt to access 50-element Vector{Int64} at index [1:1000]` in the test helper `wb_case` (`wb_helpers.jl:168`, called from `wb_records.jl:7`): the helper takes 1000 starts from a list of 50 (`lb_regular` candidates). A test defect at the spec size, not the module (the stack has no `src/` frame). See "Requests back to stage 02".
- Full package suite (all 19 folders, quality and layering first, scratch environment with Aqua, JET, ForwardDiff on `LOAD_PATH`, `--check-bounds=yes`, 2 threads): first pass **1 398 062 passed, 1 failed, 1 broken**, 14 min 16 s. The one failure was Aqua "unbound type parameters": `_first_nonfinite(vals::NTuple{K, T})` is unbound for `K = 0`; fixed by taking `vals::Tuple`. After the fix quality (12), layering (2) and the weighted-birkhoff folder were rerun (625 passed, 1 broken); the other 18 folders had passed and the fix touches only this file, so the suite is 1 398 063 passed, 0 failed, 1 broken by the sum (not rerun whole). Aqua and JET (`quality`, 12 tests) clean; layering test (2) passes; `python scripts/icm_validate.py` passes.
- Coverage (`--code-coverage=user`, folder only, default sizes): **98.4 percent** (125 of 127 executable lines; `wbirkhoff.jl` 66/66, `wbirkhoff_reducer.jl` 59/61). Uncovered: the success return of `reducer_update` (line 89, the tests call it on the failure path only, as for E2) and the `throw` for an observable that throws in `check_run` (line 132; the bad observables of AC9 return a non-finite value or a tuple rather than throw). The `.cov` files were deleted.
- Audit: no `@fastmath`; compensated (Neumaier) sums for A and B; no em dash in any touched file; files 157 and 143 lines; the longest function (`weighted_birkhoff`) is about 25 lines.

## Reducer bytes (measured, `sizeof` of the state, the test asserts `<=` the P7 formula)

`16 + sizeof(T) H (2 + 5 K)` with equality: Float64 72, 192, 544, 1952 and Double64 128, 368, 1072, 3888 for (K, H) = (1, 1), (4, 1), (4, 3), (4, 11). Workspace 0, `isbits`, `propagate_orbit` allocates 0 bytes after warm-up (AC8, and the bench).

## Performance (`bench/weighted-birkhoff.jl`, `-t 1`, median of the default BenchmarkTools run, other jobs may run on the 24-thread machine so read as orders of magnitude)

| Item | Time |
|------|------|
| reducer step Float64, (K, H) = (1, 1) / (4, 1) / (4, 11) | 21 ns / 61 ns / 254 ns, 0 allocations |
| reducer step Double64, (4, 1) | 0.53 microseconds, 0 allocations |
| `reducer_final` Float64, (1, 1) / (4, 1) / (4, 11) | 2.6 ns / 16 ns / 145 ns, 0 allocations |
| 1000 steps of one orbit, `propagate_orbit`: L1 / E2 lifted / wb (1, 1) / wb (4, 1) / wb (4, 11) / wb (4, 1) + MEGNO chain | 46 / 47 / 42 / 75 / 352 / 109 microseconds, all 0 allocations |
| `propagate`, 1024 orbits x 1000 steps, wb (4, 1), window 400 | 81 ms (about 80 ns per orbit-step; 3155 allocations of the driver) |
| `weighted_birkhoff`, K = 4, N = 1000 / 10^4 (Float64), N = 1000 (Double64) | 37 microseconds / 0.44 ms / 0.46 ms (11 allocations) |
| `birkhoff_index` | constant-folded in the bench (0.001 ns), no meaningful number |

Budget of the bench header (orientation, the spec gates no time): step at (4, 1) under 0.5 microsecond: met (61 ns); one orbit at (4, 1) within 3 times L1: met (1.6); at (4, 11) within 8 times: met barely (7.6; the eleven exponentials per handoff). The reader at (1, 1) costs nothing over the stream (within noise of L1).

## GPU smoke (local device, NVIDIA GeForce RTX 5060 Laptop GPU, scratch scripts outside the repo; `verify/` not touched)

Reducer `WeightedBirkhoff(window = (250, 1000), observables = (cos q, sin q, p lifted, cos 2q))`, chain with `Benettin`, 256 orbits, N = 2000, `propagate(...; backend = :gpu)`; run after the last code change.

- Float64: provenance backend `:gpu`, no fallback; the kernel compiles and runs (no `InvalidIRError`); every final is complete and finite on K_std = 0.5 and 10.
- Identical samples (K_std = 0: addition, product and wrap only; 256 of 256 final orbit points bitwise equal on both backends): device reducer against the CPU reducer, the largest difference over orbits and both windows is 1.71 eps S in the averages (gate of AC13(b): 4 eps S) and 0.081 eps S in Delta (gate 8 eps S); the finals are not bitwise equal (device `sin`, `cos` and `exp` differ in the last bits).
- Orbit level, K_std = 0.5 and 10 (reported, not gated, 13(c)): the largest difference in `avg1` is 0.79 and 461 (the lifted `p` observable grows on the chaotic orbits; the orbits themselves diverge).
- Double64: the kernel does not compile; `propagate` falls back to the CPU with the recorded reason "kernel does not compile for the element type" (the E2 precedent, 13(d)). The result is then the CPU result (bitwise equal).

## Deviations from the spec and the plan

None in the API or the numbers. Readings the code fixes (all from the plan's Items 5, 10, 11, 17, all in the spec by D84):
1. `birkhoff_index` for `S = 0` sets `zero` and `below_floor` together (Item 5); the digits are `NaN`.
2. `S` includes the weight-zero samples (Item 10); `reducer_final` before a horizon completes gives `NaN` for the digits and `delta` and clear flags (`below_floor = false`, `zero = false`), and a running `S`.
3. The constructor rejects `observables` that are not a nonempty tuple of at most 32767 isbits callables, naming `observables` (the plan names the empty tuple and the type; the isbits and size rules follow Item 17 and the reducer's `isbits` requirement).
4. `check_run` calls each observable at `x0` inside `try`; a throw becomes an `ArgumentError` naming `observable` with the exception text (host only). `check_run` does not check observable indices against the dimension (no `d` keyword, unlike E2): the `x0` evaluation catches a bad index as a throw.
5. `WB_MIN_WINDOW` aliases `NF_E2_MIN_WINDOW` (one home); `WB_FLOOR_MULTIPLIER = 32` is new and is the default of `birkhoff_index` and the provenance value.
6. The Aqua fix above: `_first_nonfinite` takes `::Tuple` (a signature fact, no numeric change).

## Requests back to stage 02 (no test was edited)

1. `WB_FULL_RUN = true`: `wb_case` (`test/weighted-birkhoff/wb_helpers.jl:168`) indexes 1000 starts from the 50 candidates of `lb_regular()`, so R1, R3, R6, R7, R8 error at the spec size (`BoundsError ... [1:1000]`, 5 errors; criteria testsets all pass). Either cap the number of starts at what the class provides, or draw the missing starts as the plan's record lines intend. This blocks nothing at the default sizes and must be fixed before `verify/weighted-birkhoff.jl` runs the records at spec size (stage 04).
2. Coverage gap in the tests only: no test makes an observable throw inside `check_run` (the `ArgumentError` path of `check_run` for a throwing observable, line 132 of `wbirkhoff_reducer.jl`) and none calls `reducer_update` on the success path. Optional.

## Known limits and untestable here

- AC13 (device testset) is skipped in the package suite (no device in the test process); the GPU smoke above covers compile and the identical-sample comparison; the D68 form with samples read back from the device kernel (trace kernel of `verify/gpu.jl`) is stage 04.
- The observable's smoothness is the caller's; a discontinuous or non-periodic observable is accepted (R3).
- The state is `16 + sizeof(T) H (2 + 5 K)` bytes per orbit (1.4 kB at (4, 11) in Float64); the retention decision is the study's (P10). Double64 does not run on the device.
- `K` is limited to 32767 (an `Int16` failure code); the handoff count saturates at the `Int32` maximum, harmless since horizons freeze after `2N` handoffs.
- A chain with E2 has `2N + 1` handoffs; `check_run` is exact for the stand-alone reducer (D84 item 1).
- Not run here (stage 04): criterion 12 external form in `verify/`, the D68 device form of 13(a) and (b), 1 versus N threads in child processes (AC11), the records at spec size, the cost harness R5 on the device.
