# Implementation notes: rli (catalog I9)

Code: `src/Systems/pair.jl` (new, 65 lines: `PairMap`, `PairTangent`), `src/Indicators/rli.jl` (new, 179 lines: `RLI`, `RLIState`, `check_frame` method, `pair_start`, `pair_distance`). Edited: `src/Systems/Systems.jl` (export line, one `include`, module docstring), `src/Indicators/Indicators.jl` (export line, one `include`, docstring), `src/Systems/CONTEXT.md` and `src/Indicators/CONTEXT.md` (Files, Exports, invariant, tests, limit, spec pointer lines). Bench: `bench/rli.jl` (10 entries). The Tangent layer is unchanged. Tests untouched (`git diff --stat -- test` empty). Nothing committed. Not done here (outside the files named for this stage; the megno run touched them): `src/precompile.jl`, `docs/src/api/*.md`, `CHANGELOG.md`, `src/Experiments/CONTEXT.md`.

## Status

Whole suite, `julia --project=. -t 2 -e "using Pkg; Pkg.test()"`, Julia 1.13: 1393090 pass, 0 fail, 0 error, 9 min 31 s in the testset summary. Suites: quality 12 (Aqua and JET), layering 2, analysis 434, diagnostics 245, experiments 4861, fli-ofli 263533, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, megno 1000962, precision 150, rli 64469 (33 testsets, 39 s in the full run), sali 9210, standard-map 474, tangent-core 808. The rli suite passed at the first run (64469 of 64469, about 1 min including first compile).

## Design

- `PairMap(sys)`: an `AbstractMap` holding the base map. `dimension` doubles. `forward`, `inverse` split the state into its halves (`S / 2` components each, so a wrapped base map has state length 4 per half), apply the base map to each and `vcat` them; bitwise the two single orbits. Pure Systems code, no Foundation import (D23).
- Block action (D76 item 5): `tangent(P, z)` returns `PairTangent(A, B)`, an isbits `AbstractMatrix{T}` of size `2d x 2d` holding the two `d x d` tangents. `size` and `getindex` (zeros off the diagonal blocks) make `Matrix(...)` work. `Base.:*(::PairTangent, ::SMatrix)` takes the `d` top rows and the `d` bottom rows of `W` as `SMatrix{d,k}`, forms `A * top` and `B * bottom` (the same `d x d` times `d x k` operation as the single-orbit stream) and assembles the `2d x k` result. The zero blocks are never multiplied. `src/Tangent/orbit.jl` forms `tangent(sys, x) * W`, which dispatches to this method, so the Tangent layer needed no change.
- `RLI(; replace = true)` follows `MEGNO`: isbits struct with the chain-position flag; `RLIState{T}` holds three `CompensatedSum` (`sA`, `sB`, `Q`), the failure value, `n`, `K`, failing step, code, column. Step: `safe_norm` of the two arriving columns (the zero block of each column adds exactly 0 to the sum of squares, so the norm is bitwise the single-orbit one); on a zero or non-finite norm the failure is recorded (column 1 when both fail) and state and `W` stay unchanged; else `s_i += log(alpha_i)` (the Benettin operations), `DL = abs(total(sB) - total(sA)) / T(n)`, `Q += DL`. `reducer_final` recomputes `DL` from the stored sums and takes `RLI = total(Q) / T(K)`; both zero at `K = 0`. The `DL` of the last handoff is not stored (same expression, so bitwise the value used for `Q`).
- `reducer_init` rejects `k != 2` and odd `d` with the pinned `ArgumentError` messages.
- `check_frame(RLI(), W0; orbit)`: norm defect first (column 1 then 2, `8 effective_eps`), then support (column 1 zero in the second half, column 2 zero in the first half, threshold 0), then bitwise equality of the block vectors (column 2, threshold 0). `ArgumentError` for a frame without 2 columns or with an odd row count.
- `pair_distance(sys, z)`: for `StandardMap` and `WrappedStandardMap` the coordinate differences of the two halves reduced by `u - 2pi * round(u / 2pi)` (exact for a difference below pi, any magnitude otherwise), `safe_norm`; other maps the plain difference. `pair_start` builds A (`torus_state` for the wrapped map), B = A + `Dx u` (wrapped: `Systems._reduce` with the counts of A, the construction of `Lyapunov.displaced`, copied because Indicators may not import Lyapunov), computes `Dx_r` by `pair_distance` and rejects before any map call. `Dx` is widened with the `add` rule (`ArgumentError` if it would narrow); non-finite or non-positive `Dx` and a zero, non-finite or out-of-`[Dx/2, 2 Dx]` `Dx_r` are `NumericalHealthError("separation")` with the pair id.
- Reducer state size: 88 bytes in Float64, 144 in Double64 (spec proposal 64 and 112: the `n`, `K`, `step` and failure fields are extra, as for MEGNO).

## Deviations from the spec and test plan

- None in behaviour. Placement as D76 item 4: `pair_start` and `pair_distance` are in Indicators, `PairMap` in Systems. `PairTangent` is a second exported name in Systems (the API table asked only for "a block-aware object").
- `Dx_r` is not bitwise `safe_norm(separation(...))` of Lyapunov (the minimal image uses `round`, not the two comparisons); the tests ask 4 eps and the suite passes (AC7a). Both reduce a small difference with no operation.

## Performance (`bench/rli.jl`, `verify` environment, 2 threads, BenchmarkTools medians, machine not otherwise idle, indicative)

| Entry | Time | Allocations |
|-------|------|-------------|
| `propagate_orbit`, pair of the standard map K = 6, 1000 steps, `RLI()` (k = 2, d = 4) | 66.0 us (66 ns per step) | 0 |
| same, chain `{RLI(replace = false), MEGNO()}` | 104.3 us | 0 |
| single orbit, `MEGNO()`, 1000 steps (reference) | 53.4 us | 0 |
| `reducer_step` Float64, replacer and reader | 25.4 ns | 0 |
| `reducer_step` Double64 | 274 ns | 0 |
| `propagate`, 1024 pairs, 1000 steps, CPU, `RLI()` | 36.4 ms (35.5 ns per pair-step) | 71 (chunk bookkeeping) |
| `PairTangent` times a 4 x 2 frame | 12.5 ns | 0 |
| `pair_start` Float64 | 5.6 ns (likely constant-folded, indicative) | 0 |
| `check_frame` Float64 | 49 ns | 10 (400 bytes; range slices of the frame, not in the stream) |

The spec states no time budget for RLI (memory formula only), and `bench/` holds no budget file. The pair costs about 1.2 times a single MEGNO orbit here because the standard-map step dominates.

## Audit

| Check | Result |
|-------|--------|
| GREEN | 1393090 pass, 0 fail, 0 error, every earlier suite included (above) |
| Quality | Aqua and JET clean (quality 12/12, `Indicators` and `Systems` are in the OWN_MODULES list); layering 2/2. Coverage: `src/Indicators/rli.jl` 74 of 75 executable lines (98.7 percent; the uncovered line is the `return` after the throw in `reducer_update`, the direct-use success path), `src/Systems/pair.jl` 26 of 26 (100 percent); `--code-coverage=user` on the rli suite, `.cov` files parsed and deleted afterwards |
| Size | longest function `check_frame` 20 lines, `pair_start` 12, `reducer_step` 11; files 179 and 65 lines; all under 50 and 400 |
| Tests untouched | `git diff --stat -- test` empty |
| Numerics | no `@fastmath`; the three sums per pair are `CompensatedSum` (Neumaier), the log sums with the operations of L1; `DL` formed from the totals including the correction terms (F5); `T(n)`, `T(K)` weights, no integer arithmetic on large counts |
| Git | `git status` shows only: `src/Systems/pair.jl`, `src/Indicators/rli.jl`, `bench/rli.jl` (new), `src/Systems/Systems.jl`, `src/Indicators/Indicators.jl`, the two `CONTEXT.md` files (edited), this notes file |

## Known limits

- k = 2 only, even dimension, the block frame supplied by the caller; the stream does not run `check_frame`, callers must (D67 form). A frame whose block vectors differ is rejected, but a frame of the right shape built from a different `xi_0` than the one used for the direction `u` is the caller's choice (P3 uses `u = xi_0`).
- The reducer does not know its chain position; one replacer per chain, last (D71, D72).
- `reducer_provenance` carries the indicator name and convention only; `Dx`, `Dx_r`, mode, direction rule, seed, frame source and chain are recorded by Experiments from the `pair_start` record (D76 item 6).
- `pair_distance` ignores the winding counts of the wrapped state (the minimal image of the stored coordinates is the lifted separation modulo the periodicity, exact for pairs closer than pi).
- GPU: not run here (criterion 15 is stage 04). `PairMap`, `PairTangent`, `RLI` and `RLIState` are isbits and use only `log`, `+`, `*`, `/`, `abs`, comparisons and the scaled norm; the block-aware tangent action must still be shown to compile for the kernel at stage 04. Criteria 12(c) (1 versus N threads in child processes), 14 and 15 are stage 04.
- Docs: no `docs/` build was run, and the `docs/src/api` pages were not edited (the `@autodocs` blocks should pick up the new exports; the megno run added a sentence to `indicators.md`).
