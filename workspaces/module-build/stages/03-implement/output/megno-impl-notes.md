# Implementation notes: megno (catalog I8)

Code: `src/Indicators/megno.jl` (new, 132 lines), the existing layer `Chaoslab.Indicators` (`Indicators.jl` now 23 lines: `MEGNO` added to the export line, one `include`). SALI, FLI files unchanged. Also changed: `src/precompile.jl` (two `propagate` calls, `{MEGNO(replace = false), SALI}` and `{MEGNO}`), `src/Indicators/CONTEXT.md` (file, exports, invariants, tests, deferred, limits), `src/Experiments/CONTEXT.md` (one clause on the chain rule), `docs/src/api/indicators.md` (one sentence; the `@autodocs` block picks the new names up), `CHANGELOG.md` (Unreleased, Added). Bench: `bench/megno.jl` (8 entries). Tests untouched (`git diff -- test` empty, `test/megno/` is the stage 02 state). Nothing committed.

## Status

Whole suite, `julia --project=. -t 2 -e "using Pkg; Pkg.test()"`, Julia 1.13.2: 1328621 pass, 0 fail, 0 error, 10 min 10 s in the testset summary (10 min 47 s wall including start). Suites: quality 12 (Aqua and JET), layering 2, analysis 434, diagnostics 245, experiments 4861, fli-ofli 263533, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, megno 1000962 (1 min 52 s in the full run, 26 testsets), precision 150, sali 9210, standard-map 474, tangent-core 808. The megno suite passed at the first run (1000962 of 1000962).

## Design

- `MEGNO(; replace = true)` is an isbits struct holding the chain-position flag (D72 item 1), exactly as `FLI`: `replace = true` replaces `W` by the column-normalised matrix (the same `c / alpha` operations as `SALI`); `replace = false` returns `W` unchanged for a chain with SALI after it. The reducer cannot see the chain; the one-replacer rule is the caller's (Experiments).
- `MEGNOState{T}`: two `MegnoSums{T}` (one per column, each five `CompensatedSum`: `s`, `R1`, `R2`, `Q`, `Q2`), the failing value, `n` (last successful map step), `K` (successful handoffs), failing step, status code, failing column. This is the whole reducer state: the D44 relaunch (P9) carries exactly this, no extra mechanism. Size: 200 bytes in Float64, 368 in Double64 (spec proposed 168 and 328: the failure record and the `step`, `n`, `K` fields are extra; FLI is 96 bytes in Float64).
- Step: scaled norms (`safe_norm`, D69) of the two arriving columns; a zero or non-finite norm sets the failure record (column 1 first, else 2) and returns state and `W` unchanged. Else per column `a = log(alpha)`, `s += a` (the same Neumaier operations as SALI, FLI, Benettin, so `sum` and `comp` are bitwise equal), `R1 += n a`, `R2 += n^2 a`, `Q += 2 R1 / n`, `Q2 += 3 R2`, all in the element type with `n` converted by `T(n)`. The weights are never integers (spec Algorithm item 2): `n^3` is formed as `T(n)^3` at the final only.
- `K` is counted by the reducer and `tau = n / K` (exact integer division), so the reducer needs no `tau` field and none is stored (Amendment D74 item 4).
- `reducer_final`: `(n, handoffs, Y, Ybar, Ybar20, sum, comp, R1, R2, Q, Q2)`; `Y = 2 R1 / n`, `Ybar = Q / K`, `Ybar20 = tau Q2 / n^3` from the totals (sum plus compensation); all zero at `n = 0`. Series are `SVector{2, T}`, column 1 primary. Never clipped, capped or replaced; no early stop.
- `reducer_provenance`: `(indicator = "megno", convention = "weights at the end of the window, mean over the handoffs 1..K")`. In a chain the later reducer's name wins in the run record (D34 A2).
- `check_frame(MEGNO(), W0; orbit)`: a method of the existing `check_frame`, the same rule as FLI (unit columns to 8 effective eps), with its own `ArgumentError` text. Written as a copy of the FLI method (9 lines) rather than a shared helper, to leave the verified fli file untouched.
- The weighted-sum form of M5 is used (not the Abel form of M6): it needs no history, and the tests check both against the BigFloat references.

## Deviations from the spec

None. Amendment D74 item 4 (final fields: handoff count and weighted sums; `tau` and chain position are not fields) is implemented as stated. The state size differs from the spec proposal (above), which the spec leaves to this stage.

## Performance (bench/megno.jl, `verify` environment, 2 threads, BenchmarkTools medians, machine not otherwise idle, indicative)

| Entry | Time | Allocations |
|-------|------|-------------|
| `propagate_orbit`, K = 6, 1000 steps, `MEGNO()` (k = 2) | 53.2 us (53 ns per step) | 0 |
| same, `SALI()` | 37.5 us | 0 |
| same, chain `{MEGNO(replace = false), SALI()}` | 125.7 us (noisy; about the sum of the two plus one tangent product) | 0 |
| `reducer_step` Float64, replacer and reader | about 100 ns (timer resolution) | 0 |
| `reducer_step` Double64 | 400 ns | 0 |
| `propagate`, 1024 orbits, 1000 steps, CPU, `MEGNO()` | 42.7 ms (41.7 ns per orbit-step) | 3143 (chunk bookkeeping) |
| `check_frame` Float64 | below the timer resolution (not meaningful) | 0 |

Against FLI (32.6 ns per step, 19.6 ns per orbit-step in `propagate`) the extra cost is the 8 compensated additions and the division per handoff. No per-module budget file exists in `bench/`.

## Realized error ratios (suite `@info` lines, stage 03 run, to be re-recorded at stage 04)

- AC2 shear, max error / 6 B_S: 0.006. AC3 linear hyperbolic: 0.019. AC4 fixed point: 0.026. AC5 elliptic and rotation, tau = 1: 0.014. AC10(b) BigFloat(256) on its own orbit: max error / B_M 0.019.
- AC6 same-orbit, max error / B_M (three series), `worstS` is the log sums against 8 eps Psi_max (reported):

| K_std | Float64 tau 1 | Float64 tau 4 | Double64 tau 1 | Double64 tau 4 | worstS Float64 | worstS Double64 |
|-------|---------------|---------------|----------------|----------------|----------------|-----------------|
| 0.5 | 0.0064 | 0.0035 | 0.0003 | 0.0001 | 0.026 | 0.001 |
| 2 pi | 0.0165 | 0.0159 | 0.058 | 0.058 | 0.049 / 0.053 | 0.199 |
| 10 | 0.026 | 0.017 | 0.446 | 0.446 | 0.100 / 0.089 | 1.63 |

  The Double64 K = 10 log sums exceed the borrowed 8 eps Psi_max by 1.63 (Amendment D74 item 7: reported, not asserted); the three series stay at 0.45 of B_M.
- AC7 (b) and (c), error / (8 eps max(1, max |s_1|)): 0.317 (K = 0.5), 0.183 (2 pi), 0.145 (10).
- AC8 (growth laws, ensemble means, Float64): K = 2 pi: mean X(1e4) 1.1716, 2 m_Ybar 1.1734, m_Y 1.1782, 4 m_20 1.1758; K = 5: 0.9488, 0.9448, 0.9505, 0.9475; K = 10: 1.6196, 1.6250, 1.6224, 1.6241. Regular class K = 0.5 at n = 1e5: no start outside the intervals (Ybar_1 about 1.995, Ybar20_1 about 0.4999). R5: random unit vectors, Ybar_1(1e4) min, median, max = 0.039, 1.999, 2.149.

## Audit

| Check | Result |
|-------|--------|
| GREEN | 1328621 pass, 0 fail, 0 error, every earlier suite included |
| Quality | Aqua and JET clean (quality 12/12; the OWN_MODULES list picks `Indicators` up); layering 2/2 (Indicators imports Foundation, Systems and Tangent only); coverage of the new code: `src/Indicators/megno.jl` 52 of 52 executable lines (100 percent), measured with `--code-coverage=user` on the megno suite and parsing the `.cov` files, which were deleted afterwards |
| Size | longest function `reducer_step` 11 lines, `reducer_final` 8, `check_frame` 9; files 132, 23, 68, 65, 127 lines; all under 50 and 400 |
| Tests untouched | `git diff -- test` empty |
| Numerics | no `@fastmath`; all five sums per column are `CompensatedSum` (Neumaier), the log sum with the same operations as L1, SALI, FLI; no integer weights, so no overflow of `n^3` (AC1(v) at 4e6 passes) |
| Docs | `julia --project=docs docs/make.jl` exits 0, no warnings (the per-layer API pages are small, no size warning) |

## Known limits

- k = 2 only (`ArgumentError` otherwise); map dimension d >= 2, the closed forms of the spec use d = 2.
- The reducer does not know its chain position; a chain with two replacers (or none where a later reducer needs unit columns) is a caller error that no test or runtime check can see (D72 item 1). Experiments must pass `replace = false` before SALI.
- `check_frame` is not called by the stream; callers with supplied frames must call it (D67 form).
- MEGNO depends on `tau` (M12); default `tau = 1` matches the sources. The convention (weights at the end of the window, mean over the handoffs 1..K) is P3, recorded in the provenance.
- The regular level 2 (0.5 for `Ybar20`) holds for a deviation with a component across the invariant curve; approached as O(ln K / K). A single-column result near tangency is slow (F3, R5).
- Rounding: the weighted sums grow like `n^3` (Q2); in Float64 `T(n)^3` is inexact past `n` about 2e5 (relative error one rounding), inside the B_M of AC6 and AC10 (checked to n = 1e5 in the suite; stage 04 reports the full-run sizes).
- GPU: not run here (criterion 15 is stage 04). Reducer and state are isbits and use only `log`, `+`, `*`, `/`, comparisons and the scaled norm; BigFloat is CPU-only as for the other reducers. The GPU register pressure of a 200-byte state is not measured.
- Criteria 12(c) over 1 versus N threads in child processes, 14 (ChaosTools one-column line) and 15 (CPU vs GPU) are stage 04.
