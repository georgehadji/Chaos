# Indicators

Chaos indicators as reducers of the tangent stream and the pure functions behind them. Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer), so the shared norm `safe_norm` comes from Foundation (D69). Built by the sali run (catalog I2); gali (I3), fli-ofli and megno follow in M4.

Spec: `workspaces/module-build/stages/01-spec/output/sali-spec.md` (A1 to A16, criteria 1 to 18; decision D66, Amendments D67 and D68). Verified: `stages/04-verify/output/sali-verification.md`. Test plan: `stages/02-tests/output/sali-test-plan.md`. Notes: `stages/03-implement/output/sali-impl-notes.md`.

## Files

| File | Holds |
|------|-------|
| `Indicators.jl` | module, exports, include order |
| `sali.jl` | `SALI` reducer and its isbits state `SALIState` (`k = 2`) |
| `alignment.jl` | `alignment_indices`, `sali_pair`, `sali_from_log`, `sali_floor`, `check_frame` |

## Exports

| Name | What it does |
|------|--------------|
| `SALI()` | reducer for `k = 2`: compensated `ln` of both column norms, columns renormalised, `SALI(n) = min(d_-, d_+)` from the stored unit columns; final `(n, sali, L, sum, comp, floor)` with `L = -(S_1 + S_2)` |
| `alignment_indices(a, b)` | `(d_minus, d_plus)` of two nonzero vectors of any dimension |
| `sali_pair(a, b)` | `min(alignment_indices(a, b)...)` |
| `sali_from_log(L)` | SALI from `L = ln GALI_2` (spec A4), no cancellation |
| `sali_floor(T)` | `4 effective_eps(T)`, the measurement floor of the direct form |
| `check_frame(SALI(), W0; orbit)` | `SALI(0)` of a supplied frame, or `NumericalHealthError("initial frame")` below `10 sali_floor(T)` |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| SALI is computed from the stored unit columns, so an exactly aligned pair stays at 0 for every later handoff | spec A16, AC14(d) |
| The orbit is unchanged by the reducer; the CPU path allocates nothing after warm-up | AC14(a), AC14(c) |
| A zero or non-finite column norm is recorded in the state with its column and raised as `NumericalHealthError` by the stream, never thrown inside a kernel | AC13(a) |
| `S_1`, `S_2` are Neumaier sums; `L` has no `1/g` amplification and stays finite where `exp(L)` underflows | spec A7, A15, AC8 |
| `SALI(0)` comes from `check_frame`, not from the reducer; the caller runs it before the stream on every supplied frame | Amendment D67 item 1, AC13(c) |
| For a supplied non-orthonormal frame the caller adds `ln abs(det W0)` to `L` | Amendment D67 item 3 |

## Tests and verification

- Package: `test/sali/` (23 testsets, AC1 to AC15 and the API).
- Stage 04: `verify/sali.jl` (ChaosTools GALI_2 oracle, AC16 and AC17) and `verify/gpu.jl` (AC18).
- Bench: `bench/sali.jl`.

## Deferred

| Item | Where |
|------|-------|
| GPU run (kernel path is isbits and compiles on the CPU; CUDA run is AC18, stage 04) | sali spec criterion 18 |
| Catalog row I2 and benchmarks row "SALI, regular orbit" wording for 2D maps | sali spec P3, stage 05 |
| GALI_k with k > 2, FLI/OFLI, MEGNO | later M4 modules |

## Limits

- Callers that supply an initial frame must run `check_frame` on it before the stream (D67); the stream does not, and a dependent frame would read as SALI = 0.
- `sali_from_log` loses accuracy near SALI = sqrt 2 (condition number about 1 / cos(phi), spec Amendment D68, R1); use the direct form there.
- The computed SALI lies in [0, sqrt 2 (1 + 2 eps)] (R2). For BigFloat the floor follows the default precision at call time.
