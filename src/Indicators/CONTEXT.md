# Indicators

Chaos indicators as reducers of the tangent stream and the pure functions behind them. Imports Foundation, Systems and Tangent; never a Lyapunov module (a peer layer), so the shared norm `safe_norm` comes from Foundation (D69). Runs: sali (I2), fli-ofli (I5, FLI only; OFLI to M5), megno (I8), rli (I9), stretching-numbers (I10, L13; added the reducer workspace to Tangent), naff-fma (I12, I12b: frequency map analysis of an orbit series, the one module that reads the orbit, not the deviation matrix; no Tangent change).

- Per run, in `workspaces/module-build/stages/`: spec `01-spec/output/<run>-spec.md`, test plan `02-tests/output/<run>-test-plan.md`, notes `03-implement/output/<run>-impl-notes.md`; sali also `04-verify/output/sali-verification.md`.
- Spec labels and decisions: sali A1 to A16, criteria 1 to 18 (D66, D67, D68); fli-ofli Q1 to Q13, 1 to 16 (D71, D72); megno M1 to M13, 1 to 15 (D73, D74); rli RL1 to RL13, 1 to 15 (D75 to D77); stretching-numbers SN1 to SN17, 1 to 16 (D78, D79); naff-fma NF1 to NF34, criteria 1 to 15, D80a to D80c, D81a, D81b (D80, D81).

## Files

| File | Holds |
|------|-------|
| `Indicators.jl` | module, exports, include order |
| `sali.jl`, `alignment.jl` | `SALI` reducer and its isbits state `SALIState` (`k = 2`); `alignment_indices`, `sali_pair`, `sali_from_log`, `sali_floor`, `check_frame` (SALI method) |
| `fli.jl` | `FLI` reducer and its isbits state `FLIState` (`k = 2`), `check_frame` (FLI method), `fli_unlogged` |
| `megno.jl` | `MEGNO` reducer and its isbits state `MEGNOState` (`k = 2`, five compensated sums per column), `check_frame` (MEGNO method) |
| `rli.jl` | `RLI` reducer and its isbits state `RLIState` (`k = 2`, pair frame), `check_frame` (RLI method), `pair_start`, `pair_distance` |
| `naff_fft.jl`, `naff.jl`, `naff_map.jl` | series path (E1): radix-2 FFT for the coarse search; `window_weights`, `window_frequency` (`WindowFrequency`); `frequency_map` (`FrequencyMap`), `check_frequency_map`, `circle_distance`, `frequency_distance`, `frequency_index` |
| `naff_reducer.jl` | E2: `birkhoff_weight(s)`, `phase_increments`, `weighted_phase_frequency`, `WeightedPhaseAdvance` with `PlanarObservable`, `LiftedAngleObservable`, state `WPAState`, `check_run` method |
| `stretching.jl` | `StretchingNumbers` and its isbits state `StretchState` (`k = 2`, counts in the Tangent workspace), `spectral_distance`, `spectrum_cell`, `spectrum_density`, `stretching_range`, `check_run`, `check_frame` (method) |

## Exports

| Name | What it does |
|------|--------------|
| `SALI()` | reducer for `k = 2`: compensated `ln` of both column norms, columns renormalised, `SALI(n) = min(d_-, d_+)` from the stored unit columns; final `(n, sali, L, sum, comp, floor)` with `L = -(S_1 + S_2)` |
| `alignment_indices(a, b)`, `sali_pair(a, b)`, `sali_from_log(L)`, `sali_floor(T)` | `(d_minus, d_plus)` of two nonzero vectors; their `min`; SALI from `L = ln GALI_2` without cancellation (spec A4); the floor `4 effective_eps(T)` |
| `check_frame(SALI(), W0; orbit)` | `SALI(0)` of a supplied frame, or `NumericalHealthError("initial frame")` below `10 sali_floor(T)` |
| `FLI(; replace = true)` | reducer for `k = 2`: compensated `ln` of both column lengths read before normalisation, running maxima `FLI_B` (basis, primary) and `FLI_1` (first column) including the `s = 0` term, maximisers; `replace = false` leaves `W` to a later reducer; final `(n, fli, fli1, sum, comp, stepB, colB, step1)` |
| `MEGNO(; replace = true)` | reducer for `k = 2`: per column the compensated log sum (bitwise that of SALI, FLI, L1) and the weighted sums `R1 = sum n a`, `R2 = sum n^2 a`, `Q = sum Y`, `Q2 = sum 3 R2`; delivers `Y`, `Ybar` (standard (1,-1) form) and `Ybar20` (the (2,0) mean, regular level 1/2) per column, column 1 primary; `replace = false` leaves `W` to a later reducer; final `(n, handoffs, Y, Ybar, Ybar20, sum, comp, R1, R2, Q, Q2)`, the series `SVector{2}` |
| `check_frame(FLI() or MEGNO(), W0; orbit)` | larger norm defect of a supplied frame, or `NumericalHealthError("initial frame")` above `8 eps` |
| `RLI(; replace = true)` | reducer for `k = 2` on a `PairMap` (even `d`, block frame): per column the compensated log sum (bitwise L1), `DL = abs(S_B - S_A) / n` with the correction terms, `RLI = Q / K` the mean of `DL` over the handoffs; `replace = false` leaves `W` to a later reducer; final `(n, handoffs, DL, RLI, sum, comp, Q)` |
| `check_frame(RLI(), W0; orbit)` | norm defect of a supplied block frame, or `NumericalHealthError("initial frame")` for a non-unit column, support outside its block, or unequal block vectors |
| `pair_start(sys, x0, u, Dx; orbit)` | `(z, Dx, Dx_r, mode)`: pair state (A; B = A + `Dx u`, wrapped with the counts of A), realized separation, `"wrapped"` or `"lifted"`; `NumericalHealthError("separation")` for a `Dx_r` outside `[Dx/2, 2 Dx]`, `ArgumentError` for a narrowing `Dx` |
| `pair_distance(sys, z)` | `safe_norm` of the minimal-image separation of a pair state (plain difference for other maps) |
| `fli_unlogged(F)` | `exp(F)`, or `nothing` where the exponential overflows |
| `StretchingNumbers(; range, bins, interval, min_per_cell, replace)` | reducer for `k = 2`: compensated log sums (bitwise L1) and of squares per column, a cumulative histogram of `z = ln(alpha)/tau` in `2 (bins + 2)` workspace counts (overflow cells never clipped); `reducer_final(r, st, ws)` adds `D`, `undersampled`, `counts`, `below`, `above` |
| `window_frequency(z; order, positions, pad, mean_removal, mirror_ratio)`, `frequency_map(series, N; ...)` | E1: fundamental frequency of a window (cosine window of order 0 to 5, stationary point of `abs2(phi)`), of the pair `[0, N)`, `[N, 2N)` per plane; `delta`, `D = log10 max(delta, 32 eps)`, flags; nothing thrown for data, `check_frequency_map` raises |
| `WeightedPhaseAdvance(; window, observable, birkhoff_order)` | O(1) reader for any `k >= 1` (columns pass unchanged): weighted mean phase advance of one or several horizons, `2 max(window) + 1` handoffs checked by `check_run(r, tau, N)`; final `(n, handoffs, complete, nu1, nu2, delta, D, below_floor)` |
| `spectral_distance`, `spectrum_cell`, `spectrum_density`, `stretching_range`, `check_run`, `check_frame(r, W0)` | D from two count vectors (exact `Int` sums); the cell rule; `N_c/(k h)`; the derived range padded by `8 eps` (standard maps, 2 x 2 linear); interval and 10^9-handoff check; unit and non-parallel columns |

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
| RLI is the mean of `DL_k` over the handoffs `1..K` (not a function of `DL` at the last handoff), `0 <= DL <= 2 L_sigma`; `DL` and `RLI` are never clipped; constant tangents give exactly 0; it depends on `tau` | rli spec RL4, RL8, RL10, I1, I2 |
| Stretching: counts sum to `k` in every column, D bitwise symmetric and in `[0, sqrt(2/h)]`, the log sums bitwise L1; a state of 104 (Float64) or 176 (Double64) bytes plus `4 (bins + 2)` per column; no allocation after warm-up; D depends on `tau`, bins and the column pair and has no threshold; at most one workspace reducer per chain | stretching spec SN4 to SN6, I1, I4, I8 |
| naff-fma: E1 is a pure function of the series (Experiments owns the feed, spec P1 (a2), criterion 16 is the study's); E2 shares one update with `weighted_phase_frequency` (bitwise equal for the planar observable); state 96 (Float64) or 176 (Double64) bytes plus 64 or 128 per extra horizon, no workspace, no allocation after warm-up; tie tolerance 64 eps of Float64 for every type | naff-fma spec Amendments D80, D81 |
| MEGNO state is 200 bytes in Float64 and 368 in Double64 (reducer bytes for the tangent-core memory formula); a step-chunked relaunch carries it whole | spec P9 |

## Tests and verification

- Package: `test/sali/` (23 testsets, AC1 to AC15 and the API); `test/fli-ofli/` (27 testsets); `test/megno/` (26 testsets); `test/rli/` (34 testsets); `test/stretching/` (42 testsets); `test/naff-fma/` (35 testsets, `NF_FULL_RUN` selects the spec sizes).
- Stage 04: `verify/sali.jl` (ChaosTools GALI_2 oracle, AC16 and AC17) and `verify/gpu.jl` (AC18).
- Bench: `bench/sali.jl`, `bench/fli.jl`, `bench/megno.jl`, `bench/rli.jl`, `bench/stretching.jl`, `bench/naff.jl`.

## Limits

- Deferred: GALI_k (D70) and OFLI (flows, D71) go to M5; L14 (FTLE grids) is an Experiments product (D78).
- Callers that supply an initial frame must run `check_frame` on it before the stream (D67); the stream does not, and a dependent frame would read as SALI = 0.
- `sali_from_log` loses accuracy near SALI = sqrt 2 (condition number about 1 / cos(phi), spec Amendment D68, R1); use the direct form there.
- The computed SALI lies in [0, sqrt 2 (1 + 2 eps)] (R2). For BigFloat the floor follows the default precision at call time.
- One replacer per chain, placed last (D71): `FLI(replace = false)` and `MEGNO(replace = false)` before `SALI`; `FLI()` or `MEGNO()` alone. A wrong position gives wrong numbers without an error (fli review R2), and the run record does not show `replace` or the chain (R1). Experiments must build and record chains under this rule.
- MEGNO needs a supplied frame to be unit-length (`check_frame(MEGNO(), W0)` before the stream), else `s(0) = 0` is wrong. The regular level 2 (0.5 for `Ybar20`) needs a deviation with a component across the invariant curve (spec F3, M9); `Y` of a vector along it is 0. Level and growth laws are approached slowly (O(ln K / K)); thresholds are pre-registered per study, the module never stops an orbit.
- MEGNO in Double64 over long lifted orbits: the log sums inherit the Double64 sin and cos floor at large lifted arguments from L1 (a factor 1.64 over 8 eps Psi_n at |x| near 1e6, K = 10); use the wrapped state there. For bounded or regular orbits Y carries about n eps / 10 absolute from the products of the weights and the increments (megno spec, stage 04 notes, R3).
- RLI needs the caller to build the pair (`pair_start`), the block frame (column 1 = `(xi_0, 0)`, column 2 = `(0, xi_0)`) and to run `check_frame(RLI(), W0)` before the stream: only it sees a mismatched `xi_0` between the blocks. The reducer checks the pairing at every handoff (off-block entries exactly 0, else health error "pairing", D77), so a non-pair stream or a missing block frame fails at the first handoff. Use the wrapped state (lifted Float64 loses Dx = 1e-12 near |x| = 1e4, rli R8) and Dx = 1e-12 in Float64 (rli R6). `RLI(replace = false)` before `MEGNO()` or `RLI()` alone, same chain rule as above. The direction rule, seed, frame source and chain are recorded by Experiments (D76 item 6).
- Stretching numbers: `range` is the caller's (`stretching_range(sys)` for the standard maps), the caller runs `check_run(r, tau, N)` and `check_frame(r, W0)` before the stream, and the series of D is read by calling `reducer_final(r, st, ws)` at the evaluation handoffs (no emission channel). The counts are 8 (bins + 2) bytes per orbit in every `propagate` final (8 kB at 1000 cells): the recorder must consume chunk by chunk. The 2D panel does not show the printed regular constant of D; it is not a chaos certificate (spec F5). GPU: the Float64 kernel with the workspace verified on the local GPU, Double64 kernel falls back (D68 form). Edge ties follow the rounding of the cell width, the reducer records a failure past 10^9 handoffs, and no Experiments path allocates the workspace yet.
- naff-fma: E1 needs the series of the window pair (32 N bytes per orbit in Float64) from a feed Experiments does not have yet (obligations: naff-fma test plan, "What Experiments must supply"); the closed-form interpolation, the multi-term decomposition and the 4D map of Bazzani are not built (spec P7, NF28, NF32). A genuine tie closer than 1.4e-14 in a Double64 series is decided by the Float64 transform (D81 item 7). The index separates regular from chaotic orbits only per window pair (sticky orbits, spec F6); thresholds are pre-registered per study. A planar sample exactly at the origin makes its window undefined (`NaN`, no failure).
