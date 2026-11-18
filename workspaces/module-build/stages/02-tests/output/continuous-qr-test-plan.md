# Test plan: continuous-qr

Status: approved at the stage 02 checkpoint, 2026-11-18 (owner: "Approve tests"), with the E_max default 1e-2 of U6. RED run of the new suite: 42 pass, 35 error, 77 total, all errors `UndefVarError` (`ContinuousQR` 34, `cqr_rhs` 1). Run against a throwaway scratch implementation outside the repo (session scratch folder `cq-cal`, never in the repo): 1615 pass, 0 fail, 0 error, 71.8 s.

Spec: `workspaces/module-build/stages/01-spec/output/continuous-qr-spec.md` (approved at the stage 01 checkpoint, D120; 22 acceptance criteria, reproduction lines R1 to R6). Suite: `test/continuous-qr/` (`runtests.jl`, `cq_helpers.jl`, `cq_ac01_05.jl`, `cq_ac06_10.jl`, `cq_ac11_15.jl`, `cq_ac16_22.jl`, 1104 lines, 42 testsets).

Conventions copied from ofli: every testset begins by naming the new symbols it uses (`_ = (ContinuousQR, ...)`), so a missing module errors that testset only; no top-level call and no helper calls new API at include time; helpers carry the prefix `cq_`; package sizes stay under the D38 ceiling of 3 minutes (full run 72 s with the scratch); the full spec sizes sit behind `CQ_FULL_RUN`; test-side doubles act as negative controls. Nothing in `src/` is edited at stage 02.

## API proposal for stage 03

Layer `Chaoslab.Lyapunov`, file `src/Lyapunov/continuous_qr.jl`. The module owns a per-orbit loop (not a tangent reducer): the orbit is advanced by the delivered `IntegratorMap`, the Jacobian comes from `vector_field_jacobian(flow, z)`.

| Name | Meaning |
|------|---------|
| `ContinuousQR(map, flow; order = 4, penalty_factor = nothing, orthogonality_tol = CQ_EMAX, form = :auto)` | `map` is an `AbstractMap` with `forward` and a `.tau` field. `order` in {2, 4} (Heun on one orbit step, RK4 on two orbit steps with nodes A_n, A_{n+1}, A_{n+1}, A_{n+2}); anything else throws `ArgumentError`. `form` in `:auto` (k = d gives 84, k < d gives 92), `:f84`, `:f92`; `:f84` with k != d throws `ArgumentError`. `penalty_factor` c_nu defaults to 1 for k = d and 2 for k < d, nu = c_nu \|\|A\|\|_F |
| `cqr_rhs(Q, A, nu; form)` | right-hand side G of the penalty form, pure |
| `cqr_growth(Q, A)` | u_i = (Q^T A Q)_ii, pure |
| `propagate_cqr_orbit(cq, x0, Q0, N; orbit = 1, backend = :cpu)` | one orbit; returns `(x, Q, final)`, zero allocation after warm-up. `final` fields `n, h, t, l, lambda, E, Emax, nrhs, njac` |
| `propagate_cqr_series(cq, x0, Q0, N; every, orbit = 1)` | same run, also returns `series`, records `(n, l, E, Q)` every `every` Q-steps; the series is a bitwise prefix of the run |
| `propagate_cqr(cq, xs, N; k, seed, backend, chunk, Q0)` | many orbits; `PropagationResult(x, W, finals, provenance)`, `finals[j]` equals the orbit `final`, `W` holds Q. Provenance keys `cqr_h, cqr_integrator, cqr_integrator_order, cqr_order, cqr_form, cqr_penalty_factor, cqr_orthogonality_tol, cqr_flow` plus `backend, device, backend_fallback_reason, precision, seed, k` |
| `cqr_bytes(cq, Val(d), Val(k), T)` | memory formula of C21: `sizeof(T) (d + 4dk + 2k + 3 (d/2)^2) + 16` (720 bytes for Float64, d = k = 4) |

Errors: `NumericalHealthError(what, value, threshold, step, orbit, backend[, column])` with `what` = "state" (non-finite orbit state), "deviation matrix" (non-finite Q or Jacobian), "orthogonality" (guard E_max). `ArgumentError` for odd N with order 4, N <= 0, k = 0 or k > d, all before any work. Exports: `ContinuousQR`, `cqr_rhs`, `cqr_growth`, `propagate_cqr_orbit`, `propagate_cqr_series`, `propagate_cqr`, `cqr_bytes`.

## Criterion to testset

| AC | Testset (file) | Sample and bound |
|----|----------------|------------------|
| 1 | ALG, AC1, AC1 controls (`cq_ac01_05.jl`) | upper triangular A, Q = I stays bitwise; rounding-only growth; controls: upper-triangle U, flipped S, wrong sign of nu |
| 2 | AC2 (`cq_ac01_05.jl`) | skew A, l_i = 0 to rounding, Q(t) = exp(At) Q0 against 256-bit exponential, capped B_lin |
| 3 | AC3, AC3 control, AC3 oracle (`cq_ac01_05.jl`) | hyperbolic fixtures, individual l_i against Gram-determinant reference; control: a sorted-exponent double misses 5 of 20 frames |
| 4 | AC4 (`cq_ac01_05.jl`) | ratios of consecutive h in [2^q / 2, 2 x 2^q] above the rounding floor |
| 5 | AC5, AC5 control (`cq_ac01_05.jl`) | rotated frame, l_1, sum, Q(t) against closed form; control checks l_1 and Q only (U3) |
| 6 | AC6a, AC6b (`cq_ac06_10.jl`) | one column against ln \|\|Y q\|\| on F1, inverted oscillator, pendulum; column independence k = 1, 2, 4 on F3 |
| 7 | AC7 (equilibrium, identity) (`cq_ac06_10.jl`) | flow-direction identity C17 on pendulum A = 2, 3, F1, F3 R and C; h to h/2 ratios reported |
| 8 | AC8 (`cq_ac06_10.jl`) | energy-covector identity C18 on 11 orbits up to T_valid |
| 9 | AC9 (`cq_ac06_10.jl`) | same-orbit agreement with Householder, RGS and Benettin; lambda form at T = 60, l form at the long horizon (U9) |
| 10 | AC10a, AC10b (`cq_ac06_10.jl`) | long-time slope of lambda_1 on a regular orbit; bounded maxima |
| 11 | AC11 (`cq_ac11_15.jl`) | statistics over F3 classes: pairing, zero pair, comparison with the discrete variant B |
| 12 | AC12, AC12 control (`cq_ac11_15.jl`) | sum of exponents against volume growth; polar-factor double (U5) |
| 13 | AC13a to AC13d (`cq_ac11_15.jl`) | E_max bound, defect scaling Delta^q (R2), injected fault 1.05 Q raises "orthogonality", stability limits 1.1x and 1.3x |
| 14 | AC14, AC14 control (`cq_ac11_15.jl`) | same-orbit rounding against BigFloat(256) reference, Float64 and Double64; control: single-precision Q^T A Q |
| 15 | AC15 (`cq_ac11_15.jl`) | precision ladder Float64, Double64, BigFloat(256) |
| 16 | AC16 (`cq_ac16_22.jl`) | orbit bitwise equal to the integrator alone, to the discrete spectrum and to Benettin |
| 17 | AC17a, AC17b, AC17c (`cq_ac16_22.jl`) | chunk invariance, zero allocation, series is a prefix of the run |
| 18 | AC18 (`cq_ac16_22.jl`) | N + 1 Jacobian calls, N orbit steps, 2N right-hand sides, memory formula |
| 19 | AC19a to AC19e (`cq_ac16_22.jl`) | equilibrium accepted, health-error kinds, argument errors, provenance and default E_max |
| 20 | AC20 (`cq_ac16_22.jl`) | Double64 order-10 variational reference as package stand-in; the Taylor form is stage 04 |
| 21 | none | ChaosTools comparison is stage 04 only |
| 22 | AC22 CPU part (`cq_ac16_22.jl`) | Float64 and Double64 as two backends, `:gpu` request falls back to the CPU with its reason; device part is stage 04 |
| R1, R2, R5 | R1, R2, R5 (`cq_ac16_22.jl`) | reported (cost ratio, equilibrium defect, no-penalty growth of form 92); R3, R4, R6 are reported from AC3, AC5, AC6a, AC9, AC12, AC14, AC20 |

## Tests that pass already

Seven testsets use no new name and pass in the RED run (42 tests): AC1 negative controls (10), AC3 control (21), AC3 oracle (3), AC5 negative control (4), AC7 equilibrium refusal (2), AC12 polar control (1), AC14 control (1). They test the oracles and the test-side doubles, so a wrong helper fails now rather than at stage 04.

## Spec points found while writing the tests

- U1: the equilibrium rejection ("initial frame") is done by the test-side direction-frame builder; the module accepts f(z0) = 0 (AC19a).
- U2: the series is delivered by a separate `propagate_cqr_series`; a `Union{Nothing, Vector}` field in the return tuple of one function allocates 640 bytes or more per call and breaks AC17b.
- U3: the sum l_1 + l_2 cannot detect a wrong S sign (the trace is invariant under an orthogonal Q), so the AC5 control checks l_1 and Q only.
- U4: Jacobian-NaN step convention. Order 2: call 1 is the initial one, call j + 1 is orbit step j, so the 10th call reports step 9.
- U5: the polar-factor double does not give a sum within rounding (the stage values are not orthonormal, C15). It passes the AC12 bound without any penalty, so the sum alone does not test the penalty; the spec wording "within rounding" is not reachable.
- U6: the E_max default of 1e-3 (P5) is too tight. Measured E_max at t = 100 on F3 (10 starts, order-8 table): h = 1/20, q = 2: 8.3e-3 (k = 1, 2), 6.4e-3 (k = 4); h = 1/20, q = 4: 2.3e-3, 2.2e-3, 1.2e-4; h = 1/40, q = 2: 1.76e-3, 1.76e-3, 1.53e-3; h = 1/40, q = 4: 1.06e-4, 1.02e-4, 7.3e-6. The defect scales as Delta^q and E_eq is about Delta^2 \|\|S\|\|^2 / 2 for q = 2 (penalty evaluated on off-manifold stages). The tests use `CQ_EMAX = 1e-2` for q = 4 and 5e-2 for q = 2 (AC19e checks the default equals `CQ_EMAX`). OWNER DECISION needed on the default.
- U7: B_col 5e-3 (realized 1.3e-3 at t = 60). B_dir: q = 4 5e-3 and 1e-1 for pendulum amplitude 3; q = 2 1.0, pendulum amplitude 3 skipped for q = 2 (realized error 3.4).
- U8: B_en discretization allowance 3e-2 (spec 1e-3; realized worst 1.55e-2).
- U9: lambda form at T = 60: c_d = 10 (q = 4), 2 (q = 2) instead of 5 and 20 (realized ratio 0.47 of the allowance). At the long horizon the lambda form of the spec fails by about 17x on orbit C at T = 200 (differences grow in events along the orbit); the l form is used, \|l^C - l^D\| <= 5e4 h^4 (realized 1.7e4 h^4 on C, 2e3 to 3.5e3 elsewhere). T = 1000 is not calibrated.
- U10: AC13d stability limits: 1.1x and 1.3x the limit raise "orthogonality", 0.5x and 0.9x do not; the check is final E <= E0 because E_max rises transiently to about 4e-3.
- U11: B_lin is vacuous at large t, so AC2, AC5 and AC6a use the cap `min(B, 0.05)` (AC2 and AC5 control 0.25).
- U12: the AC3 control uses h = 0.025 (B_lin 0.028 against a frame-dependent spread up to 0.49) and requires at least one miss; 5 of 20 frames miss.
- U13: `CountingFlow` counts force calls only, so Jacobian counting needs the test-side wrapper `cq_CountFlow`; `IntegratorMap`'s tau is a field, not a method.

## Calibration (throwaway scratch implementation outside the repo, session scratch folder `cq-cal`)

The scratch implemented the proposal in about 90 lines of Julia and ran the whole suite once green. Realized values, h = 1/20, t = 60 unless said:

| Test | Bound | Realized |
|------|-------|----------|
| AC3 worst error / B_lin | 1 | 0.017 |
| AC5 worst error / bound | 1 | 0.14 |
| AC6a worst error / bound | 1 | 0.39; absolute 1.6e-4 (q = 4), 1.9e-2 (q = 2) |
| AC6b column independence | 5e-3 | 1.3e-3 |
| AC7 error (q = 4) pend 2 / pend 3 / F1 / F3 R / F3 C | 5e-3, 1e-1, 5e-3 | 8.7e-4, 2.3e-2, 7.5e-5, 1.2e-4, 4e-4 |
| AC7 error (q = 2) pend 2 / F1 / F3 R / F3 C | 1.0 | 0.235, 0.012, 0.205, 0.338 (pend 3: 3.4, skipped) |
| AC8 worst error | 3e-2 | 1.55e-2, T_valid = 100 on 11 orbits |
| AC9 worst lambda error / allowance | 1 | 0.47; l form worst 1.7e4 h^4 against 5e4 |
| AC12 worst \|sum\| / B_sum | 1 | 3.4e-5 |
| AC13a worst E / bound | 1 | 1.0e-3 (the bound is 1e-2 / 5e-2, so about 10x margin to the guard) |
| AC14 worst error / B_C | 1 | 3.8e-4 (Float64), 7.1e-6 (Double64) |
| AC20 worst error / B | 1 | 0.14 (q = 4); q = 2 excluding pend 3 passes at 1.0 |
| Order ratios h to h/2 (AC4, AC7, AC13b) | [2^q/2, 2 x 2^q] | 4.0 to 4.7 (q = 2), 8.6 to 20 (q = 4, AC7); AC7 F3 C q = 4 at t = 10 gave 0.17 (rounding-limited, reported only) |
| R1 cost ratio q = 4 / q = 2 to discrete RGS | reported | 0.30, 0.32 |
| R5 time to E = 1e-2, form 92, no penalty / with penalty | reported | 173.5 / Inf |

Total scratch run: Pass 1615, Fail 0, Error 0, 71.8 s (RED 18.2 s).

## Seeds

StableRNGs. `cq_seed(ac, sub) = 20261012 + 44000 + 50 ac + sub` (`CQ_BASE = CQ_SEED + 44000`); frames and starts through `initial_deviations(cq_seed(ac, sub), orbit, Val(d), Val(k), T)`. Henon-Heiles starts: orbit R at (0.1, 0), orbit C at (-0.015, 0), energy `CQ_E`, h = `CQ_H` = 1/20.

## RED log

`julia --project=<scratch env> red_cq.jl` (suite only, no implementation): Pass 42, Error 35, Total 77, 18.2 s. Errors: `UndefVarError: ContinuousQR` 34, `UndefVarError: cqr_rhs` 1. No parse error, no helper error. The 42 passes are the seven testsets of "Tests that pass already".

## Stage 04 note

Criterion 21 (ChaosTools), the device part of 22 (GPU, `verify/gpu.jl`), the thread-count comparison of 17, the Taylor-form external reference of 20 and the full spec sizes (`CQ_FULL_RUN`) run at stage 04. The T = 1000 allowance of AC9 must be calibrated there. Stage 03 must not loosen a bound in this file without an amendment.
