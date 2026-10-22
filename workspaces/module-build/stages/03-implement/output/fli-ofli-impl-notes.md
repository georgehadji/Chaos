# Implementation notes: fli-ofli (FLI only, catalog I5)

Code: `src/Indicators/fli.jl` (new, 122 lines), the existing layer `Chaoslab.Indicators` (`Indicators.jl` now 22 lines: one export line and one `include`). Sali files are unchanged. Also changed: `src/precompile.jl` (two `propagate` calls, `{FLI(replace = false), SALI}` and `{FLI}`), `src/Indicators/CONTEXT.md` (file table, exports table, invariants, tests, deferred), `docs/src/api.md` (one sentence of the Indicators section; the `@autodocs` block picks the new names up), `CHANGELOG.md` (Unreleased, Added). Bench: `bench/fli.jl` (9 entries). Tests untouched (`git diff -- test` empty, `test/fli-ofli/` is the stage 02 state). Nothing committed.

## Status

Whole suite, `julia --project=. -t 2 -e "using Pkg; Pkg.test()"`, Julia 1.13.2: 327659 pass, 0 fail, 0 error, 6 min 55 s (7 min 26 s wall including start). Suites: quality 12 (Aqua and JET), layering 2, analysis 434, diagnostics 245, experiments 4861, fli-ofli 263533 (45 s in the full run, 27 testsets), lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, precision 150, sali 9210, standard-map 474, tangent-core 808. The fli-ofli suite passed at the first run (263533 of 263533, the count of the stage 02 stub dry run).

## Design

- `FLI(; replace = true)` is an isbits struct holding the chain-position flag (D72 item 1): `replace = true` replaces `W` by the column-normalised matrix (the same `c / alpha` operations as `SALI`); `replace = false` returns `W` unchanged, for a chain with SALI after it. The reducer cannot see the chain, so the one-replacer rule is the caller's (Experiments).
- `FLIState{T}` (isbits for Float64, Double64, BigFloat excepted from allocation-freedom): two `CompensatedSum`, the two running maxima `f1`, `fB`, the failing value, `n`, failing step, the maximisers `step1`, `stepB`, `colB`, a status code and the failing column. About 3 compensated pairs of `T` plus 4 values of `T` and small integers.
- Step: scaled norms (`safe_norm`, D69) of the two arriving columns; a zero or non-finite norm sets the failure record (column 1 first, else 2) and returns state and `W` unchanged (as SALI). Else `log(alpha_i)` goes into the Neumaier sums (the same operations as SALI and Benettin, so the sums are bitwise equal, criterion 8a), then `fli_maxima` updates `FLI_1` and `FLI_B` by strict comparison on the compensated totals, both starting at 0 (the s = 0 term, P4). Ties keep the earlier maximiser (D72 item 5), column 1 before column 2.
- `reducer_final`: `(n, fli, fli1, sum, comp, stepB, colB, step1)`; `fli` is `FLI_B` (primary), `fli1` is `FLI_1`. The maximisers are part of the reducer state, so the D44 relaunch carries them (P8) with no extra mechanism: the state is what is carried.
- `check_frame(FLI(), W0; orbit)`: accepts a column whose norm defect is at most `8 effective_eps(T)` (D72 item 7), returns the larger defect, else `NumericalHealthError("initial frame", defect, threshold, 0, orbit, :cpu, column)`. `ArgumentError` for a frame without 2 columns. A method of the existing `check_frame`.
- `fli_unlogged(F)`: `exp(F)` or `nothing` on overflow.
- Log form never clipped, no cap, no early stop.

## Deviations from the spec

None beyond Amendment D72 (items 1, 5 and 7 are implemented as stated: replacer option, last-strict-increase maximiser, 8 eps frame defect). API names follow the test plan: `FLI`, `fli_unlogged`, `check_frame(FLI(), ...)`, final fields `n, fli, fli1, sum, comp, stepB, colB, step1`. `reducer_provenance` is `(indicator = "fli",)`; in a chain the later reducer's name wins in the run record (D34 A2), so chain position and frame source are recorded by Experiments.

## Performance (bench/fli.jl, `verify` environment, 2 threads, BenchmarkTools medians, machine not otherwise idle, indicative)

| Entry | Time | Allocations |
|-------|------|-------------|
| `propagate_orbit`, K = 6, 1000 steps, `FLI()` (k = 2) | 32.6 us (32.6 ns per step) | 0 |
| same, `SALI()` | 37.5 us | 0 |
| same, chain `{FLI(replace = false), SALI()}` | 65.4 us | 0 |
| `reducer_step` Float64, replacer and reader | about 100 ns (timer resolution) | 0 |
| `reducer_step` Double64 | 200 ns | 0 |
| `propagate`, 1024 orbits, 1000 steps, CPU, `FLI()` | 20.1 ms (19.6 ns per orbit-step) | 3143 (chunk bookkeeping) |
| `check_frame`, `fli_unlogged` Float64 | below the timer resolution (not meaningful) | 0 |

The chain costs about the sum of the two reducers' extra work (two scaled norms and two logs per handoff, done twice, P1 (b)) on top of the common tangent product. No per-module budget file exists in `bench/`.

## Audit

| Check | Result |
|-------|--------|
| GREEN | 327659 pass, 0 fail, 0 error, every earlier suite included |
| Quality | Aqua and JET clean (quality 12/12; the OWN_MODULES list picks `Indicators` up); layering 2/2 (Indicators imports Foundation, Systems and Tangent only); coverage of the new code: `src/Indicators/fli.jl` 48 of 48 executable lines (100 percent), measured with `--code-coverage=user` on the fli-ofli suite and parsing the `.cov` files, which were deleted afterwards (alignment.jl 2 of 2, sali.jl 22 of 24 in that run: only the sali failure branch, covered by the sali suite) |
| Size | longest function `reducer_step` 14 lines, `check_frame` 9, `fli_maxima` 6; files 122, 22, 68, 65 lines; all under 50 and 400 |
| Tests untouched | `git diff -- test` empty |
| Numerics | no `@fastmath`; both log sums are `CompensatedSum` (Neumaier); maxima compared on `total(S_i)` |
| Docs | `julia --project=docs docs/make.jl` exits 0 with `checkdocs = :exports`; no missing-docstring warning; one size warning only: `api.html` 152.9 KiB over the 100 KiB warn limit (limit 200; was 148.6 KiB) |

## Known limits

- k = 2 only (`ArgumentError` otherwise); map dimension d >= 2, the closed forms of the spec use d = 2.
- The reducer does not know its chain position; a chain with two replacers (or none where a later reducer needs unit columns) is a caller error that no test or runtime check can see (D72 item 1). Experiments must pass `replace = false` before SALI.
- `check_frame` is not called by the stream; callers with supplied frames must call it (as for SALI, D67).
- FLI depends on `tau` (spec Q10): the maximum runs over handoffs only; the default `tau = 1` matches the sources.
- The unlogged value overflows Float64 above 709.78; use `fli_unlogged`, which returns `nothing` there.
- GPU: not run here (criterion 16 is stage 04). The reducer and state are isbits and use only `log`, `max`, comparisons and the scaled norm; BigFloat is CPU-only as for the other reducers.
- Criteria 13(c) over 1 versus N threads in child processes and 15 (ChaosTools `gali` line) are stage 04; realized error ratios against B_FLI are reported by the suite as `@info` lines (realized in the stage 03 run: max error / B_FLI 0.02 to 0.11 in Float64, 0.0005 to 0.66 in Double64, the largest at K = 2 pi; BigFloat 0.074) and recorded at stage 04.
- `api.html` is at 152.9 KiB of the 200 KiB hard limit; the API page should be split before the next indicator module.
