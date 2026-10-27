# Implementation notes: naff-fma (stage 03)

Spec: `workspaces/module-build/stages/01-spec/output/naff-fma-spec.md` (with Amendments D80 and D81). Tests: `test/naff-fma/` (untouched since `3b49469`; `git diff 3b49469 -- test` is empty). Plan: `workspaces/module-build/stages/02-tests/output/naff-fma-test-plan.md`.

## What was built

| File | Lines | Holds |
|------|-------|-------|
| `src/Indicators/naff_fft.jl` | 31 | radix-2 FFT (no new dependency), used only to locate coarse-search candidates |
| `src/Indicators/naff.jl` | 240 | `window_weights`, `WindowFrequency`, `window_frequency` and its internals (E1, NF2, NF17 to NF20) |
| `src/Indicators/naff_map.jl` | 117 | `circle_distance`, `frequency_distance`, `frequency_index`, `FrequencyMap`, `frequency_map`, `check_frequency_map` |
| `src/Indicators/naff_reducer.jl` | 224 | `birkhoff_weight(s)`, `phase_increments`, `weighted_phase_frequency`, `PlanarObservable`, `LiftedAngleObservable`, `WeightedPhaseAdvance`, protocol methods, `check_run` method, provenance |
| `src/Indicators/Indicators.jl` | +8 | exports (the 14 names of the plan plus the types `WindowFrequency`, `FrequencyMap`) and includes; `src/Chaoslab.jl` re-exports them unchanged |
| `src/precompile.jl` | +3 | one E2 `propagate` and one `frequency_map` per element type |
| `bench/naff.jl` | 48 | BenchmarkTools suite |
| `docs/src/api/indicators.md`, `src/Indicators/CONTEXT.md` (78 lines), `CHANGELOG.md` | | docs, context, entry (catalog status not flipped: stage 05) |
| `_config/numerics-defaults.toml` | 0 | the seven P10 keys were added at stage 03 and struck by the owner at the stage 04 start (the module does not read them; the keys stay keyword defaults, as for sali, megno and stretching) |

No change in Tangent, Systems or Foundation. The pure E2 function and the reducer share one update (`_accumulate`), one weight (`_wb`, mirrored index, D81 item 10), one increment (`_phase_inc`) and one frequency formula, so the planar stream result is bitwise the function on the series sampled at the handoffs (AC10c). The pure function works in BigFloat (`pi` at the precision of the data, plan item 14).

## Results

- Folder `test/naff-fma/`, `julia --project=. --check-bounds=yes -t 2`: **848 passed, 0 failed, 0 errors** (35 testsets, 44 s warm; the plan's stub took 87 s). With `NF_FULL_RUN = true` (spec sizes: 20 starts, N up to 4096 and 10001, N = 16000 in AC8(c)): **858 passed**, 9 min.
- Full package suite (all 19 folders, quality and layering first, scratch environment with Aqua, JET, ForwardDiff): **1 397 405 passed, 0 failed**, 14 min 18 s. Aqua and JET (`quality`, 12 tests) clean; layering test (2) passes; `python scripts/icm_validate.py` passes.
- Coverage of the four new files with `--code-coverage=user` (folder only): **99.3 percent** (304 of 306 lines; `naff.jl` 150/151, `naff_fft.jl` 26/26, `naff_map.jl` 43/43, `naff_reducer.jl` 85/86). Uncovered: the 20-iteration non-convergence return of `_stationary` (`naff.jl:177`) and the success return of `reducer_update` (`naff_reducer.jl:191`; the tests call it only on the failure path).
- Audit: no `@fastmath`; no em dash in any touched file; files 31, 240, 117, 224 lines; every function under 50 lines; the largest, `_coarse`, is about 30.
- Realized values reported by the tests (stage 04 will state them again): AC10d error over B_E2 at most 0.056 (Float64) and 0.0018 (Double64); AC7b leakage ratio 0.99 (Float64) and 1.33 (Double64) of one bound; D81a largest error times 2N 0.6259 in both types; no tie and no non-unique reference in the AC2, AC9, AC13 windows at the default sizes.

## Performance (`bench/naff.jl`, default flags, one thread, median, other jobs were running on the 24-thread machine so read as orders of magnitude)

| Item | Time |
|------|------|
| `window_frequency`, Float64, N = 10^3 / 10^4 / 10^5 | 0.33 ms / 5.3 ms / 77 ms (36 allocations: the FFT and spectrum buffers, `L >= 8 N` complex, 16 MB at N = 10^5) |
| `frequency_map`, pair of windows, N = 10^3 / 10^4 | 0.67 ms / 10.8 ms (N = 10^5: twice the window value, about 0.15 s); Double64 at N = 10^3: 5.1 ms |
| `window_weights` / `birkhoff_weights`, N = 10^4 | 0.14 ms / 0.27 ms |
| `phase_increments` (2 x 10^4 samples) / `weighted_phase_frequency` (2 windows of 9999) | 0.18 ms / 0.37 ms |
| reducer step, planar Float64 / lifted Float64 / planar Double64 | 22 ns / 9 ns / 0.72 microseconds, 0 allocations |
| 1000 steps of one orbit, `propagate_orbit`: L1 / E2 lifted (window 400) / MEGNO / stretching / E2 + MEGNO chain | 40 / 36 / 50 to 83 / 57 / 112 microseconds, all 0 allocations |
| `propagate`, 1024 orbits x 1000 steps, E2 lifted | 41 ms (about 40 ns per orbit-step, 3128 allocations of the driver) |

So the E2 reader is cheap beside the stream (one orbit of 1000 steps with E2 alone is within the noise of the L1 run), and the series path costs about 0.3 microseconds per sample of a window at N = 10^3 to 0.8 microseconds at N = 10^5 (the FFT, `log N` growth, and the Newton iterations on the full window). The ratio of the host transform to the orbit pass (spec R5) is a stage 04 measurement.

## Deviations from the spec and the plan

1. **Coarse search (E1 step 3).** The spec allows a transform "only for the starting bin"; implemented as a Float64 radix-2 FFT of the zero-padded window (length the next power of two of at least `max(pad N, 8 N)`), whose local maxima within 1.5 percent of the top (at most 64) name the candidates; the grid points `m0 - 2 .. m0 + 2` of the true grid `-1/2 + m/(pad N)` around each candidate are evaluated by direct sums (Float64, exactly reduced phases, blocks of 64) and decide the bin, the tie flag (64 eps of Float64, D81 item 7) and the lowest-index rule. A first version with a 0.1 percent threshold and `L >= pad N` lost a near-tie of two peaks 0.16 percent apart under a frequency shift (AC4(c), scalloping of the FFT grid); the 8 N padding bounds the scalloping loss of a Hann peak under 0.5 percent (derived: `0.5 |phi''(0)| (pi/16)^2`, `|phi''(0)| = 0.219`), the threshold has a 3x margin. Not covered: more than 64 local maxima within 1.5 percent of the top (a broad-band window), where the maximum of the grid could be missed; none seen in the suite.
2. **`phase_increments` at the origin.** A sample exactly 0 gives `NaN` for the two increments that touch it, in the pure function as in the reducer (one shared `_phase_inc`). The spec states the origin rule for the reducer only ("stage 03 behaviour, a proposal"); the pure function sharing it keeps the bitwise equality of AC10c.
3. **`reducer_final` has no status field** (plan Items 19): an undefined window shows as `NaN` with `complete == true`; an unfinished one as `NaN` with `complete == false`.
4. **Types exported:** `WindowFrequency`, `FrequencyMap` (not in the plan's list; needed to name the result in documentation).
5. **Keys of P10 added** to `_config/numerics-defaults.toml` (spec P10: "filled when the owner approves", approved by D80; plan item 17; amendment D81: "filled by stage 03 once the owner approves"). Seven keys: `frequency_min_window = 16`, `frequency_pad = 8`, `frequency_floor_multiplier = 32`, `frequency_mirror_ratio = 0.9`, `frequency_weight = "hanning"`, `frequency_window_order = 1`, `frequency_mean_removal = true`. The module does not read them (keyword defaults, checked only through the tests' `settings` echo); `frequency_mean_removal = true` is the spec's value for the angle-like observable while the function default is `false` as the plan fixes (the choice is recorded per run). No earlier module put its keys in this file (the stretching keys of its spec are not there), so the owner may prefer to strike them: they are a separate, easily removed hunk. The experiments and analysis suites were rerun after the edit (see the report). Struck by the owner on 2026-10-27: the file is back to its earlier content.
6. **FFT dependency.** None added; a 31-line radix-2 transform serves. A threaded transform library, if later adopted, brings the thread-count obligation of the spec (set to 1 when orbit threads exceed 1); the in-house FFT is single-threaded, so nothing is owed now.
7. **`CONTEXT.md`.** The five per-run spec bullets were merged into a path pattern plus one label line, and "Deferred" became a bullet of "Limits" to stay under 80 lines; no fact was dropped.

## Known limits and untestable here

- The series feed of Experiments (criterion 16) does not exist (moved to the study spec, D81 item 2). E1 holds 32 N bytes per orbit in Float64 plus about 24 N bytes of transform buffers inside `window_frequency` at pad 8; the chunk capacity formula must count both.
- E1 element types: Float64 and Double64 are tested; BigFloat is not (the coarse search converts to Float64; nothing in the code forbids it). E2: Float64, Double64, BigFloat (pure function).
- Tie decisions of a Double64 series rest on a Float64 coarse value (D81 item 7); the refined frequency of either tied bin is in its bracket.
- Reducer state: `n` is an `Int` since the stage 04 follow-up (below); `handoffs` is an `Int32` that saturates; the `2N + 1` handoffs requirement is the caller's `check_run` (the reducer cannot see the run).
- The mirror rule, the second-line ratio and the near-zero flag use the transform's local maxima; the second-line ratio is therefore resolution-limited by the 8 N padding.
- Not run here (stage 04 or later): criteria 14 (independent BigFloat orbit, ChaosTools period routine, Python NAFF packages), 15 (GPU; the reducer is a plain isbits reader that should compile for Float64, Double64 untried), 1 versus N threads in child processes, the cost harness (R5), the panel lines R1, R3 to R6, the Eq. (79) part of R2, R9.
- No test appeared wrong. Observation for the record: AC4(c) (modulation by 0.1) is sensitive to near-ties: two spectral peaks 0.16 percent apart changed the bin of a version of the coarse search whose FFT grid did not cover the pair. The test is right; it is a good guard of Deviation 1.

## Stage 04 follow-up (independent review; items fixed in `src/Indicators/naff*.jl`, tests in `test/naff-fma/nf_review.jl`)

Folder: 895 passed (848 committed + 47 new in 5 testsets, one per fixed item, each reproducing the reviewer's scenario; all fail on the stage 03 code: 11 failures and 6 errors). The committed tests are untouched.

| Item | Change |
|------|--------|
| M1 increments near +-pi | `reducer_final` gains `near_pi`: the number of increments that entered a window with abs(Delta) > 0.9 pi (planar observable only; 0 for the lifted one; saturating `Int16`). The value of E2 is unchanged (reviewer scenario: E1 0.495, E2 0.3205, now with `near_pi > N`). Spec I9 asks for a recorded flag where this is a property of the data. |
| M2 tie at large N | the direct sums of the coarse search (`_abs_phi64`) are Neumaier-compensated and restart the exponential every 32 samples (was 64, plain sum); the real-cosine tie holds at N = 8000, 20000 and 100000. |
| M3 step counts | `n` is stored as `Int`; `handoffs` saturates at the `Int32` maximum instead of wrapping or throwing; no `InexactError` is reachable. |
| M4 provenance keys | NOT changed: `test/naff-fma/nf_api.jl` (testset "API reducer constructor and observables", lines 104 to 108) reads `indicator`, `estimator`, `window` and the other old names, so the orchestrator must amend the spec and that test (suggested names `naff_indicator`, `naff_estimator`, `naff_window`, or a nested `naff` field). The collision is real: in a chain with `Benettin` the later reducer wins. |
| L magnitudes | `window_frequency` scales the window by an exact power of two (largest component in [1, 2)) and scales the amplitudes back, so 1e-300 to 1e307 samples work and the scale-by-2^k covariance stays bitwise; `_phase_inc` falls back to power-of-two normalised samples when `z conj(zprev)` leaves 1e-280 to 1e280 (the normal path is unchanged bit for bit). |
| L non-finite transform | `_coarse` returns `nothing` for a non-finite transform and the window is `:failed`, so `maximum` of an empty set cannot throw (not reachable for finite data after the scaling). |
| L zero weight | `_accumulate` skips a weight that is exactly 0 (shared by the pure function and the reducer, so still bitwise equal): a NaN increment with the weight 0 (the first increment of each window) no longer spoils the window; AC11b planar origin stays green (the origin inside the window still gives NaN). |
| L bracket | the converged stationary point must lie within 1/M of the coarse bin as the spec says (NF2); iterates may wander to 2/M on the way. The docstring and the spec agree now; no committed test changed. |
| L pad | `pad >= 2` (ArgumentError naming `pad`) in `window_frequency` and `frequency_map`; no committed test uses `pad = 1`. |
| L indices | the observable indices are checked against 1:2d (d entries for lifted maps, 2d for the wrapped ones; the committed hand tests pass `Val(2)` with a state of 4 entries, so a bound of d is not possible) by the caller-side `check_run(r, tau, N; d = d)`, ArgumentError naming `observable`. NOT in `reducer_init`: a run-time message there made the Float64 device kernel fail to compile (InvalidIRError); no string-building throw and no `exponent` call is reachable from `reducer_init`, `reducer_step` or `reducer_final` (the binary exponent is read from the bits). The Float64 device kernel compiles again (lifted and planar probes on the local GPU). |
| not changed | a constant series with one ulp of noise is `:ok` under mean removal (spec-conformant, D81 item 8); noted in `CONTEXT.md` limits. |

State size: unchanged. Layout `hor, prev1, prev2, n::Int, handoffs::Int32, nwrap::Int16, code::Int8` is 96 bytes in Float64 and 176 in Double64 for one horizon, plus 64 and 128 per extra horizon (checked by `sizeof` in the new tests). The places that state the size and stay correct: spec Compute, I10, criterion 10(b); `reducer_final` field list (spec Algorithm, plan API table, `CONTEXT.md` row, docstring) gains `near_pi`; the spec text "`n::Int32`" appears only in the plan (Items 12), not in the spec.
