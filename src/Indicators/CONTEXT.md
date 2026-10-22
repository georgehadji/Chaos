# Indicators

Chaos indicators as reducers of the tangent stream and the pure functions behind them. Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer), so the shared norm `safe_norm` comes from Foundation (D69). Built by the sali run (catalog I2) and the fli-ofli run (catalog I5, FLI only; OFLI deferred to M5); gali (I3) was deferred to M5 (D70); the megno run (catalog I8) added the MEGNO reducer.

Spec: `workspaces/module-build/stages/01-spec/output/sali-spec.md` (A1 to A16, criteria 1 to 18; decision D66, Amendments D67 and D68). Verified: `stages/04-verify/output/sali-verification.md`. Test plan: `stages/02-tests/output/sali-test-plan.md`. Notes: `stages/03-implement/output/sali-impl-notes.md`.

FLI spec: `workspaces/module-build/stages/01-spec/output/fli-ofli-spec.md` (Q1 to Q13, criteria 1 to 16; decisions D71, D72). Test plan: `stages/02-tests/output/fli-ofli-test-plan.md`. Notes: `stages/03-implement/output/fli-ofli-impl-notes.md`.

MEGNO spec: `workspaces/module-build/stages/01-spec/output/megno-spec.md` (M1 to M13, criteria 1 to 15; decisions D73, D74). Test plan: `stages/02-tests/output/megno-test-plan.md`. Notes: `stages/03-implement/output/megno-impl-notes.md`.

## Files

| File | Holds |
|------|-------|
| `Indicators.jl` | module, exports, include order |
| `sali.jl` | `SALI` reducer and its isbits state `SALIState` (`k = 2`) |
| `alignment.jl` | `alignment_indices`, `sali_pair`, `sali_from_log`, `sali_floor`, `check_frame` (SALI method) |
| `fli.jl` | `FLI` reducer and its isbits state `FLIState` (`k = 2`), `check_frame` (FLI method), `fli_unlogged` |
| `megno.jl` | `MEGNO` reducer and its isbits state `MEGNOState` (`k = 2`, five compensated sums per column), `check_frame` (MEGNO method) |

## Exports

| Name | What it does |
|------|--------------|
| `SALI()` | reducer for `k = 2`: compensated `ln` of both column norms, columns renormalised, `SALI(n) = min(d_-, d_+)` from the stored unit columns; final `(n, sali, L, sum, comp, floor)` with `L = -(S_1 + S_2)` |
| `alignment_indices(a, b)` | `(d_minus, d_plus)` of two nonzero vectors of any dimension |
| `sali_pair(a, b)` | `min(alignment_indices(a, b)...)` |
| `sali_from_log(L)` | SALI from `L = ln GALI_2` (spec A4), no cancellation |
| `sali_floor(T)` | `4 effective_eps(T)`, the measurement floor of the direct form |
| `check_frame(SALI(), W0; orbit)` | `SALI(0)` of a supplied frame, or `NumericalHealthError("initial frame")` below `10 sali_floor(T)` |
| `FLI(; replace = true)` | reducer for `k = 2`: compensated `ln` of both column lengths read before normalisation, running maxima `FLI_B` (basis, primary) and `FLI_1` (first column) including the `s = 0` term, maximisers; `replace = false` leaves `W` to a later reducer; final `(n, fli, fli1, sum, comp, stepB, colB, step1)` |
| `check_frame(FLI(), W0; orbit)` | larger norm defect of a supplied frame, or `NumericalHealthError("initial frame")` above `8 eps` |
| `MEGNO(; replace = true)` | reducer for `k = 2`: per column the compensated log sum (bitwise that of SALI, FLI, L1) and the weighted sums `R1 = sum n a`, `R2 = sum n^2 a`, `Q = sum Y`, `Q2 = sum 3 R2`; delivers `Y`, `Ybar` (standard (1,-1) form) and `Ybar20` (the (2,0) mean, regular level 1/2) per column, column 1 primary; `replace = false` leaves `W` to a later reducer; final `(n, handoffs, Y, Ybar, Ybar20, sum, comp, R1, R2, Q, Q2)`, the series `SVector{2}` |
| `check_frame(MEGNO(), W0; orbit)` | larger norm defect of a supplied frame, or `NumericalHealthError("initial frame")` above `8 eps` |
| `fli_unlogged(F)` | `exp(F)`, or `nothing` where the exponential overflows |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| SALI is computed from the stored unit columns, so an exactly aligned pair stays at 0 for every later handoff | spec A16, AC14(d) |
| The orbit is unchanged by the reducer; the CPU path allocates nothing after warm-up | AC14(a), AC14(c) |
| A zero or non-finite column norm is recorded in the state with its column and raised as `NumericalHealthError` by the stream, never thrown inside a kernel | AC13(a) |
| `S_1`, `S_2` are Neumaier sums; `L` has no `1/g` amplification and stays finite where `exp(L)` underflows | spec A7, A15, AC8 |
| `SALI(0)` comes from `check_frame`, not from the reducer; the caller runs it before the stream on every supplied frame | Amendment D67 item 1, AC13(c) |
| For a supplied non-orthonormal frame the caller adds `ln abs(det W0)` to `L` | Amendment D67 item 3 |
| In a chain of column reducers exactly one replaces the columns, the last; `FLI(replace = false)` goes before `SALI`, `FLI()` alone replaces. The reducer cannot check this | D71, D72 item 1 |
| `FLI_1`, `FLI_B` are nondecreasing, `FLI_B >= FLI_1 >= max(0, S_1)`; the log form is never clipped (no cap, no early stop) | spec I1, Representation |
| FLI depends on `tau` (the maximum runs over the handoffs only) | spec Q10 |
| MEGNO is a linear functional of the log sums `s(0..K)` (Abel form, spec M6), not a function of `s(K)` alone; a constant added to `s` changes nothing; `Y`, `Ybar`, `Ybar20` are never clipped (Y is unbounded and may be negative) | spec I1, Representation, M7 |
| MEGNO depends on `tau` (weights at the end of the window, mean over the handoffs `1..K`); the weights are formed in the element type, no 64-bit integer `n^3` | spec M12, criterion 1(v) |
| MEGNO state is 200 bytes in Float64 and 368 in Double64 (reducer bytes for the tangent-core memory formula); a step-chunked relaunch carries it whole | spec P9 |

## Tests and verification

- Package: `test/sali/` (23 testsets, AC1 to AC15 and the API); `test/fli-ofli/` (27 testsets); `test/megno/` (26 testsets).
- Stage 04: `verify/sali.jl` (ChaosTools GALI_2 oracle, AC16 and AC17) and `verify/gpu.jl` (AC18).
- Bench: `bench/sali.jl`, `bench/fli.jl`, `bench/megno.jl`.

## Deferred

| Item | Where |
|------|-------|
| GALI_k (D70), OFLI (flows, D71) | M5 |
| MEGNO stage 04 (ChaosTools one-column tangent line, CPU vs GPU) | `megno` run, stage 04 |

## Limits

- Callers that supply an initial frame must run `check_frame` on it before the stream (D67); the stream does not, and a dependent frame would read as SALI = 0.
- `sali_from_log` loses accuracy near SALI = sqrt 2 (condition number about 1 / cos(phi), spec Amendment D68, R1); use the direct form there.
- The computed SALI lies in [0, sqrt 2 (1 + 2 eps)] (R2). For BigFloat the floor follows the default precision at call time.
- One replacer per chain, placed last (D71): `FLI(replace = false)` and `MEGNO(replace = false)` before `SALI`; `FLI()` or `MEGNO()` alone. A wrong position gives wrong numbers without an error (fli review R2), and the run record does not show `replace` or the chain (R1). Experiments must build and record chains under this rule.
- MEGNO needs a supplied frame to be unit-length (`check_frame(MEGNO(), W0)` before the stream), else `s(0) = 0` is wrong. The regular level 2 (0.5 for `Ybar20`) needs a deviation with a component across the invariant curve (spec F3, M9); `Y` of a vector along it is 0. Level and growth laws are approached slowly (O(ln K / K)); thresholds are pre-registered per study, the module never stops an orbit.
