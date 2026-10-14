# Implementation notes: lyap-svd-stabilized

Code: `src/Lyapunov/svd_spectrum.jl` (74 lines: `SVDSpectrum`, `SVDState`, `log_singular_values`, the reducer methods, `reducer_provenance`). `Lyapunov.jl` only gained the include and the two export names. Bench: `bench/lyap-svd-stabilized.jl` (3 entries, not run). No function over 50 lines. Tests untouched (`git status test` clean).

## Status

`test/lyap-svd-stabilized/` alone, 2 threads: 12915 pass, 0 fail, 0 error on Julia 1.13.1 (1 min 26 s testset, 87 s wall, 63% compilation; an earlier run on a loaded machine took 2 min 56 s) and on 1.10.12 (1 min 44 s). Both under the 3 minute ceiling (D35), so no cut from the plan was applied.

Other suites, 2 threads, same tree (with the parallel Householder files present):

| Suite | 1.13.1 | 1.10.12 |
|---|---|---|
| lyap-benettin | 15009 pass | 15009 pass |
| lyap-gram-schmidt | 19555 pass | 19555 pass |
| tangent-core | 565 pass | 565 pass |
| precision | 150 pass | 150 pass |
| diagnostics | 245 pass | 245 pass |
| standard-map | 774 pass | skipped (slow on 1.10, 100 min measured earlier) |
| layering | 2 pass | 2 pass |
| quality (Aqua, JET) | 11 pass, 1 fail | 11 pass, 1 fail |

The one quality failure is Aqua "unbound type parameters" on `_giv_sweep` and `_giv_rot` in `src/Lyapunov/householder.jl` (the parallel Householder slug); nothing from `svd_spectrum.jl`. JET and the other Aqua checks pass. The runs used throwaway environments `scratchpad/m2/svdenv113` and `svdenv110` (Chaoslab developed from this tree plus Aqua, JET, ForwardDiff), because `testenv` has no Aqua.

## Design

The state is `SVDState{T}(qr::QRState{T,2}, c::T)`. The step factors `W` once with `orthonormalize(RGS(), W)`, passes that factorisation to `_qr_step` (the shared step of `QRSpectrum(RGS())`, factored out by the Householder slug as `_qr_step(r, st, n, W, fact)`), and on success adds `(R12/R11) exp(S2 - S1)` to `c` with the totals before the window (V5); `c` is a plain sum. `S_i`, `lambda_qr` and all health records are therefore the QR step's own (AC3a bitwise, tested). `log_singular_values` is V6 verbatim. The guard is in `reducer_init` (exact message in `SVD_GUARD_MESSAGE`), so `propagate` rejects other `d` or `k` before any work.

## Deviations

- `reducer_final` is `(n, sum, comp, lambda_qr, lambda, c, u, G_max, orth_max, count_max)` as the contract; `lambda` is zero at `n = 0`.
- Dependency: the step relies on `_qr_step` in `qr_spectrum.jl`. If that helper is renamed, this file must follow. (A first version called `reducer_step` and `orthonormalize` separately, factoring twice; same results, about 1.6 us per step.)
- `c` is not updated when the QR step fails (the state and `W` are left as the QR step leaves them).

## Performance (1.13.1, 1 thread, machine loaded by parallel agents, indicative)

`reducer_step` for d = k = 2: 0 bytes allocated (AC8, AC9 zero-allocation checks pass). Time per step about 1.1 us for the QR part plus one `exp` and one divide; the double-factoring version measured 1.56 us under load.

## Known limits

- d = k = 2 only (D32). General d is held.
- `u = exp(2 (S2 - S1))` underflows to 0 for large `S1 - S2`, which is the intended behaviour.
- GPU compile of `exp`, `log`, `sqrt` for Double64 is untested here (stage 04, per-chunk fallback).
- Coverage not measured in this run.
- AC3(c) class counts: at K = 0.5 the largest relative change of `c` in the unscored class is 9.9 (reported only, as the amended test intends).
