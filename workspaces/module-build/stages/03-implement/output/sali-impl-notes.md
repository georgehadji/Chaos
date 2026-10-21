# Implementation notes: sali

Code: `src/Indicators/` (`Indicators.jl` 21 lines, `sali.jl` 68, `alignment.jl` 73, `CONTEXT.md`), a new layer included and re-exported in `src/Chaoslab.jl` after Lyapunov. Imports Foundation, Systems and Tangent only (layering test passes; the scaled norm is the layer's own, since Indicators may not import Lyapunov). Also changed: `src/precompile.jl` (one `propagate` call with `SALI()`), `docs/make.jl` (layer list), `docs/src/api.md` (new "Indicators" section), `docs/src/index.md` (layer table), `CHANGELOG.md`. Bench: `bench/sali.jl` (8 entries). Tests untouched (`git diff -- test` empty).

## Status

Whole suite, `julia --project=. -t 2 -e "using Pkg; Pkg.test()"`, Julia 1.13.2: 64126 pass, 0 fail, 0 error, 7 min 06 s. Suites: quality 12 (Aqua and JET), layering 2, analysis 434, diagnostics 245, experiments 4861, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, precision 150, sali 9210 (26 s in the full run, 23 testsets), standard-map 474, tangent-core 808. All 9210 sali assertions passed at the first run, including AC1(b) at the amended 16 eps and the absorbing zero (AC14(d)).

## Design

- `SALI` is an empty isbits struct; `SALIState{T}` holds two `CompensatedSum`, the SALI at the last handoff, the failing value, `n`, the failing step, a status code and the failing column (isbits for Float64 and Double64, so the kernel path is the same as for Benettin).
- Step: scaled norms of both columns; a zero or non-finite norm sets the failure record (column 1 first, else 2) and returns state and `W` otherwise unchanged. Otherwise `log(alpha_i)` goes into the Neumaier sums, the columns are divided by `alpha_i`, and `d_-`, `d_+` are taken as scaled norms of the difference and the sum of the two stored unit columns (no second normalisation, spec A16; this is the point of the stage 02 stub run). `reducer_update`, `reducer_failed`, `reducer_error` as `Benettin`.
- `reducer_final`: `L = -((S_1 + S_2) + (c_1 + c_2))`, zero at `n = 0`; `floor = sali_floor(T)`.
- Pure functions use one `scaled_norm` (same rule as `safe_norm`: scale by the largest entry, NaN and Inf propagate). `alignment_indices` accepts any dimension and any `AbstractVector`. `sali_from_log` clips `g` at 1. `check_frame` returns `SALI(0)` of the two columns (0 if a column norm is zero or non-finite) or throws `NumericalHealthError("initial frame", value, 10 B_meas, 0, orbit, :cpu, 0)`; it needs exactly 2 columns (ArgumentError otherwise).

## Deviations from the spec

None beyond Amendment D67 (frame check run by the caller before the stream, no stream hook; `L = -(S_1 + S_2)`; frame source and chain order recorded by Experiments). `propagate` and `propagate_orbit` do not call `check_frame` themselves (D67 item 1: the caller does), so a supplied frame reaches the map unchecked unless the caller checks it; Experiments must call it when it builds supplied frames. The unit-column handling of `reducer_step` is the only place where the arithmetic differs from a naive "normalise inside the index": `W` is replaced by the same unit columns the index is computed from.

## Performance (bench/sali.jl, scratch use of the `verify` environment, 2 threads, BenchmarkTools, machine not otherwise idle, indicative)

| Entry | Time | Allocations |
|-------|------|-------------|
| `propagate_orbit`, K = 6, 1000 steps, SALI (k = 2) | 37.7 us (37.7 ns per step) | 0 |
| same, Benettin (k = 1) | 40.1 us | 0 |
| `reducer_step` Float64, d = k = 2 | about 100 ns (timer resolution) | 0 |
| `reducer_step` Double64 | 300 ns | 0 |
| `check_frame` Float64 | about 100 ns | 0 |
| `sali_pair`, `sali_from_log` Float64 | at or below the timer resolution (100 ns) | 0 |
| `propagate`, 1024 orbits, 1000 steps, CPU | 22 ms (21.7 ns per orbit-step) | 3140 (chunk bookkeeping) |

The step of the pair costs about the same as the single-vector Benettin orbit, because the standard-map tangent product dominates. The pure-function entries are wrapped in `Ref` to defeat constant folding; sub-resolution numbers are not meaningful. There is no per-module budget file in `bench/`; the earlier modules record no budget either.

## Audit

| Check | Result |
|-------|--------|
| GREEN | 64126 pass, 0 fail, 0 error (all suites above) |
| Quality | Aqua and JET clean (quality 12/12, the OWN_MODULES list picks up `Indicators` automatically); layering 2/2; coverage of `src/Indicators/` 52 of 52 executable lines (100 percent; `Indicators.jl` has none), measured by `--code-coverage=user` on the sali suite and parsing the `.cov` files, which were deleted afterwards |
| Size | longest function (`reducer_step`) 17 lines, `check_frame` 8; files 21, 68, 73 lines; all under 50 and 400 |
| Tests untouched | `git diff -- test` empty |
| Numerics | no `@fastmath`; both log sums are `CompensatedSum` (Neumaier); `L` added as `(S_1 + S_2) + (c_1 + c_2)` |
| Docs | `julia --project=docs docs/make.jl` exits 0 with `checkdocs = :exports`; no missing-docstring warning. One size warning only: `api.html` 148.6 KiB over the 100 KiB warn limit (limit 200) |

## Known limits

- k = 2 only (`ArgumentError` otherwise); the direct form is blind below `sali_floor(T)` (use `L`), and `sali_from_log` is valid only for SALI well above eps.
- `check_frame` is not called by the stream; callers with supplied frames must call it (D67).
- `L` omits `ln abs(det W0)`; for a non-orthonormal supplied frame the caller adds it.
- GPU: not run here (AC18 and `verify/gpu.jl` are stage 04). The reducer is isbits and uses only `log`, `sqrt`, `min`, `maximum`; Double64 falls back as for the other reducers.
- `api.html` is near the 200 KiB hard limit of Documenter (148.6 KiB now); the next indicator modules should split the API page.
- Realized error ratios against B_SALI and B_L are not reported here; stage 04 records them (R4).
