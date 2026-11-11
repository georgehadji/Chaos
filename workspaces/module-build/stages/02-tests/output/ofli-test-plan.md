# Test plan: ofli

Status: approved at the stage 02 checkpoint, 2026-11-11 (owner: "Approve tests"); RED rerun 6 testsets pass, 26 error (all `UndefVarError` for `OFLI`). Points U1 to U13 are the operative reading of the spec criteria they narrow.

Spec: `workspaces/module-build/stages/01-spec/output/ofli-spec.md` (approved at the stage 01 checkpoint, D118, open points P1 to P12 as proposed; 18 acceptance criteria). Suite: `test/ofli/` (`runtests.jl` includes `of_helpers.jl`, `of_ac01_05.jl`, `of_ac06_10.jl`, `of_ac11_15.jl`, `of_ac16_18.jl`). Tags: S package suite, V verification environment (stage 04), G local GPU.

Conventions copied from henon-heiles: every testset begins by naming the new symbols it uses (`_ = (OFLI,)`), so a missing module errors that testset only; no top-level call and no helper calls new API at include time; package sizes under the D38 3-minute ceiling, full spec sizes when `verify/ofli.jl` sets `OF_FULL_RUN = true` before including `runtests.jl`; regression numbers come only from `_config/bibliography/benchmarks.md`, everything else is computed inside the test (closed forms, BigFloat references at the realized orbit points, same-run controls).

## API proposal for stage 03

Layer `Chaoslab.Indicators`, file `src/Indicators/ofli.jl`, the reducer protocol of the Tangent layer (D71/D72: the reducer sees W unnormalised at each handoff n = j tau; exactly one reducer of a chain replaces W, and it is the last).

| Name | Meaning |
|------|---------|
| `OFLI(flow; replace = true)` | reducer. `flow` supplies `vector_field` and `dimension`. With `replace = true` it replaces W by unit columns (the same renormalisation as FLI) and must be last in a chain; `replace = false` leaves W to a later replacer |
| `reducer_init(r::OFLI, T, Val(d), Val(k))` | isbits state (compensated sums S_i, running maxima G_i of the projected and G_i of the unprojected form, step and column of the maximisers, sin theta at the maximisers, minimum sin theta, step count, failure code); throws `ArgumentError` if `dimension(flow) != d` |
| `reducer_step(r::OFLI, st, n, x, W)` | non-throwing. g_i(n) = S_i(old) + ln norm(P_{f(x_n)} w_arr,i), with f the exact field at the numerical state x_n (candidate A, P1); strict increase only, ties keep the earlier maximiser. Failure code 2 (field zero or non-finite) or 1 (deviation matrix, with the column) |
| `reducer_update`, `reducer_failed`, `reducer_error` | the throwing form builds `NumericalHealthError(what, value, threshold, step, orbit, backend, column)` with `what` = "vector field" or "deviation matrix" |
| `reducer_final(r::OFLI, st)` | NamedTuple `(n, ofli, ofli1, inst, fli, fli1, sum, comp, stepB, colB, step1, sinB, sin1, sinmin, nminf)`: `ofli` = OFLI_B, `ofli1` = OFLI_1 (handoffs only), `inst` the instantaneous projected value, `fli`/`fli1` the companion unprojected maxima, `sum`/`comp` the compensated sum, the rest diagnostics |
| `reducer_provenance(r::OFLI)` | `(indicator = "ofli", ofli_flow = string(nameof(typeof(flow))), ofli_replace = replace)` merged into the `propagate` run record |
| `ofli_initial(flow, x0, W0)` | pure function, `SVector{k}` of the s = 0 terms ln norm(P_{f(x0)} w_i) (U1) |
| `check_frame(r::OFLI, W0, x0; orbit = 0)` | new method of the exported `check_frame`: for k = 1, `NumericalHealthError` with kind "initial frame" if the column is parallel to f(x0) within 8 eps; a basis (k >= 2) is accepted |

`fli_unlogged` is reused for the display form (criterion 15b). Exports: `OFLI`, `ofli_initial` (plus the new `check_frame` method). Nothing in `src/` is edited at stage 02.

## Criterion to testset

| AC | Tag | Testset (file) | Sample | Tolerance |
|----|-----|----------------|--------|-----------|
| 1 | S | AC1 algebra on the free particle (`of_ac01_05.jl`); AC1 negative controls (same file) | d = 4, h = 2^-4, N = 1e4, tau in {1, 10}, Float64 and Double64 | 16 eps max(1, abs g); controls "unprojected", "parallel/gradH" must fail |
| 2 | S | AC2 free particle d = 2 | d = 2, frame refusal included | closed form, 16 eps |
| 3 | S | AC3a, AC3b, AC3c (same file) | oscillator, 256 bits; order-2 and order-4 tables; controls | 256-bit identity; ratio 4 and 16 within slack; controls ":lagged", ":no_s0" fail |
| 4 | S, V | AC4a, AC4b, AC4c | pendulum, Double64 order-10 reference | identity 1e-18 (Taylor 1e-30 is V); continuous minimum periodic 1e-12; numerical vs exact delta 1e-5 |
| 5 | S | AC5 | Newton fixed point of the return map via `section_map`, periodic orbit of F3 | bounded OFLI, FLI slope margin 1.5 |
| 6 | S | AC6 (`of_ac06_10.jl`) | R and the regular class | slope vs ln t in [0.7, 1.3] |
| 7 | S | AC7 | chaotic class | slope vs mean exponent, FLI minus OFLI median |
| 8 | S | AC8 | two seeded orthonormal frames | spread under ln sqrt(4) |
| 9 | S | AC9 | sandwich F* - ln sqrt(4) <= OFLI_B <= F* | exact |
| 10 | S | AC10a, AC10b | chains (FLI(replace=false), OFLI), (SALI, OFLI) | bitwise sums; identities |
| 11 | S | AC11a, AC11c (`of_ac11_15.jl`) | candidate B vs A on the oscillator, pendulum, R, C | rigorous bound sin(angle) norm(w), relative slack 1e-9; angle ratios reported |
| 12 | S | AC12; AC12 negative control | 1584 rows, Float64 and Double64, both forms | error <= B_G; the single-precision double must exceed B_G |
| 13 | S | AC13 | tau in {1, 20}, R and C | OFLI_tau <= OFLI_1 + 2 B_G |
| 14 | S | AC14a to AC14e | equilibrium, zero tangent, frames, chains, run record, dimension guard | exact error fields; bitwise orbit |
| 15 | S | AC15a, AC15b | BigFloat 256 vs its own 512-bit reference; overflow G = 777 at N = 26000, h = 0.25 | error <= B_G |
| 16 | S | AC16a, AC16b (`of_ac16_18.jl`) | chunk sizes vs serial; allocation | bitwise; 0 bytes |
| 17 | S, V | AC17 | Double64 order-10 at 1/80 as stand-in | delta_8 = 1e-6 for the order-8 table; Taylor form is V |
| 18 | S, G | AC18 CPU part | regular and chaotic class, Float64 and Double64 as two backends, `:gpu` fallback | per-backend B_G; gap rule on the regular class only; device part is G |

Stage 04 only: criterion 17 Taylor (TaylorIntegration reference in `verify/`), criterion 18 device runs and the D_orb chaotic-class gap, the full spec sizes (`OF_FULL_RUN`), the thread-count comparison (U11).

## Tests that pass already

6 of the 32 testsets use no new API and pass on the code as it is: AC1 negative controls, AC3a, AC3c, AC4a, AC4b, AC12 negative control (17 assertions). They test the test-side doubles and the oracles themselves, so a defect in a double or an oracle shows now, before the module exists.

## Spec points found while writing the tests

- U1: the s = 0 term cannot be formed in the reducer (`reducer_init` sees no x0 or W0). Resolved by the pure function `ofli_initial`; `ofli` and `ofli1` cover handoffs only, and the tests combine both.
- U2: the criterion 1 controls "s = 0 dropped" and "field at x_(n-1)" cannot fail on the free particle (constant field, nondecreasing g). They run on the oscillator (AC3c, start (0, 0.3) for `:no_s0`).
- U3: the 16 eps of criterion 1 holds because h = 2^-4 (dyadic), so the tangent is exact.
- U4: criterion 10(a) "{OFLI, FLI}, FLI not replacing" is read as the chain (FLI(replace=false), OFLI()), since exactly one replacer sits last. The chain (SALI, OFLI) is tested for orbit points only.
- U5: criterion 4(b) "exact OFLI constant after the first period within 1e-12" holds for the continuous minimum only; the handoff supremum converges from below. Tested as the continuous minimum per period.
- U6: the 1e-30 identity (3a, 4a) is tested at 256 bits for the oscillator; in the package the pendulum is tested at 1e-18 against the Double64 order-10 reference (the Taylor BigFloat reference is stage 04).
- U7: criterion 14(g): h and the map's flow id are not in the reducer's record; the test checks tau, k, precision, `ofli_flow`, `ofli_replace`. The experiment layer compares flow identifiers (P10).
- U8: criterion 7 comparison OFLI_B <= FLI_B uses 4 eps max(1, abs FLI), not B_G (same S_i, so tighter).
- U9: criterion 9's sigma_1 is Float64 (from the BigFloat product).
- U10: criterion 18's B_G gap rule applies to the regular class only; the chaotic class needs the D68 orbit-divergence term D_orb (sali-spec Amendment D68), which needs the CPU reference tangent products and is stage 04. The package test reports the number of agreeing orbits.
- U11: the thread-count comparison (1 vs N) cannot be set in-process; `verify/ofli.jl` reruns with `JULIA_NUM_THREADS`.
- U12: the AC11 inequality uses a relative slack of 1e-9.
- U13: the package CPU part of 18 uses OF_NCLASS starts per class (full run: 200).

## Calibration (throwaway scratch implementation outside the repo, deleted)

Package run against the scratch implementation: 345 assertions, all pass, 53.5 s (AC1-5 6.7 s, AC6-10 11.9 s, AC11-15 25.9 s, AC16-18 8.9 s; under the 3-minute ceiling). Measured margins:

| Item | Measured | Bound |
|------|----------|-------|
| AC1 constant column (renormalisation rounding) | 8.9e-16 | 16 eps (relaxed from bitwise) |
| AC3b order 2, h = 0.1 / 0.05 | 2.11e-3 / 5.03e-4 (ratio 4.2) | predicted ratio 4 |
| AC3b order 4, h = 0.1 / 0.05 | 1.32e-5 / 8.39e-7 (ratio 15.8) | predicted ratio 16 |
| AC4c t* (reported) | Inf at h = 1/4 and 1/8 for amplitude 2; 2364 at h = 1/4 for amplitude 3 | reported |
| AC5 period, trace, T'(E) | 6.0756, -0.9035, -1.347 | reported |
| AC5c dG, dFLI (order-4 table, two orbits) | 0.0019, 0.00054; 2.15, 1.88 | dFLI margin against 1.5 is the tight one |
| AC6 slope R, regular class | 1.0032, about 1.007 | [0.7, 1.3] |
| AC7 mean slope vs mean exponent | 0.05592 vs 0.05528; median FLI minus OFLI 0.257 | reported with SE |
| AC8 spread R / C (t = 100, 1000) | 0.320, 0.351 / 0.121, 0.292 | 0.693 (ln sqrt 4) |
| AC9 (F* - OFLI_B) / ln sqrt 4 | 0.0346 to 0.732, median 0.333 | [0, 1] |
| AC11b candidate B: R | first angle above 0.1 never, angle 3e-13 | reported |
| AC11b candidate B: C | first angle above 0.1 at t = 831, ratio 0.323 (collapse) | reported |
| AC11c angle ratios | 4.06 (order 2), 16.02 (order 4); order 8 rounding dominated | reported |
| AC12 max error / B_G | 0.040 (median 4.7e-4), 1584 rows, 0 maximiser differences | 1 |
| AC13 differences | R: -1.8e-15, 3.7e-4, 7.7e-3; C: about 1e-15 | 2 B_G |
| AC15a BigFloat error / B_G | 1.05e-3 | 1 |
| AC15b error / B_G at G = 777 | 8.4e-4 | 1 |
| AC16b allocation | 0 bytes, Float64 and Double64, including the chain with SALI | 0 |
| AC17 order-4 errors R / C (t = 10, 30, 60) | 1.3e-5 to 1.7e-5 / 5.2e-6 to 5.8e-4 | reported; order-8 error passes 1e-6 (value not printed, to record at stage 03) |
| AC18 B_G values | about 1e-11 (Float64), 1e-26 (Double64), ratios 1e-3 to 1e-5; 4 agreeing orbits | 1 |

Not printed in the package run (pass within the bound, margin not recorded): AC1 worst ratio, AC3a and AC4a identities, AC4b spread, AC4c error, AC17 order 8. Record them at stage 03 from the first green run.

## Seeds

All random draws use StableRNGs with `of_seed(ac, sub) = 20261012 + 43000 + 50 ac + sub`. Offsets used: 43000 to 43999 (the range for ofli; criterion number ac selects the 50-wide block). No value comes from memory: regression numbers come from `_config/bibliography/benchmarks.md` rows or are computed in the test.

## RED log

Run: `test/ofli/` only (files `of_ac01_05.jl` to `of_ac16_18.jl` after `of_helpers.jl`), driver in a temp dir, environment with Test, Chaoslab (developed), StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, real Chaoslab with no OFLI.

- 32 testsets: 6 pass (17 assertions), 26 Error, 0 Fail.
- Totals: Pass 17, Error 26, Total 43, 5.7 s.
- All 26 errors are `UndefVarError: OFLI not defined in Main`. No parse error, no MethodError, no helper error; no helper touches new API at include time.
- Errored testsets: AC1 algebra, AC2, AC3b, AC4c, AC5, AC6 to AC9, AC10a, AC10b, AC11a, AC11c, AC12, AC13, AC14a to AC14e, AC15a, AC15b, AC16a, AC16b, AC17, AC18.
- Against the scratch implementation: 345 pass, 0 fail, 0 error, 53.5 s.

## Stage 04 note

`verify/ofli.jl` sets `OF_FULL_RUN = true` before including `test/ofli/runtests.jl`, adds the Taylor reference for criteria 3a, 4a, 17, the D_orb chaotic-class gap and device runs for 18, and the thread-count rerun (U11).
