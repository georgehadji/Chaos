# Test plan: integrator-catalog

Spec: `workspaces/module-build/stages/01-spec/output/integrator-catalog-spec.md` (approved at the stage 01 checkpoint, D111, open points P1 to P17 as proposed). Suite: `test/integrator-catalog/` (`runtests.jl` includes `ic_helpers.jl`, `ic_flows.jl`, `ic_ac01_03.jl`, `ic_ac04_08.jl`, `ic_ac09_12.jl`, `ic_ac13_16.jl`, `ic_ac17_23.jl`; picked up by `test/runtests.jl` automatically). 43 testsets. Nothing under `src/` or `data/` was touched.

## API proposal for stage 03

Layers: new flow types in `Chaoslab.Systems`; new submodule `Chaoslab.Integrators` (`src/Integrators/`, between Systems and Tangent in `src/Chaoslab.jl`, imports Foundation and Systems only). States are `SVector{2n}` in the layout (q_1..q_n, p_1..p_n) (F2, T2: (q1, q2, p1, p2)); Tao: (q; p; x; y). All functions are generic in the element type of the state (Float64, Double64, BigFloat, and `ForwardDiff.Dual`); method objects are built in an element type `T` (as the maps: `FroeschleMap(T(a), ...)`), tables isbits (tuples) so one kernel serves CPU and GPU.

### Systems: Hamiltonian flows (P1)

| Name | Meaning |
|------|---------|
| `abstract type AbstractFlow`, `AbstractSeparableFlow <: AbstractFlow` | supertypes (architecture doc names; the code has none yet, Inconsistencies (a)) |
| `HarmonicOscillator(w::Vector{T})` (F1), `Pendulum(eps::T)` (F4, H = p^2/2 - eps cos q), `Kepler()` (F2, planar, GM = 1), `HenonHeilesFlow()` (T2, spec item 18), `TaoNonseparable()` (T1, H = (Q^2 + 1)(P^2 + 1)/2; not separable) | test flows |
| `dimension(flow)`, `hamiltonian(flow, z)` | 2n; H in the element type |
| `drift(flow, z, t)`, `kick(flow, z, t)` | exact partial flows of T = \|p\|^2/2 (q += t p) and of V (p -= t grad V(q), for F4 including eps), either sign of `t`; separable flows only |
| `drift_jacobian(flow, z, t)`, `kick_jacobian(flow, z, t)` | elementary Jacobians (SMatrix) |
| `hamiltonian_gradient(flow, z)`, `vector_field(flow, z)` | (dH/dq, dH/dp) and J grad H; all flows including T1 |
| `exact_flow(flow, z, t)`, `exact_flow_matrix(flow, t)` | F1 only, rotation in each plane |
| `angular_momentum(flow, z)`, `runge_lenz(flow, z)` | F2; Runge-Lenz convention A = p x L - q/\|q\| (planar 2-vector) |
| `perturbation_scale(flow)`, `corrector_force(flow, z)`, `corrector_hessian(flow, z)` | eps of the flow (1 for flows without one); gradient and Hessian in q of W = \|grad V\|^2 with V the unscaled perturbation (F4: W = sin^2 q, force sin 2q), for the SABAC corrector |
| `CountingFlow(flow)`, `kick_count(cf)`, `reset_counts!(cf)` | wrapper counting `kick` calls; forwards everything else |

### Integrators

| Name | Meaning |
|------|---------|
| `advance(m, flow, z, tau)` | one step of size `tau` (either sign); a flow integrated by a method |
| `advance_n(m, flow, z, tau, n; compensated = false)` | n steps; adjacent kicks of consecutive steps share one evaluation (KDK leapfrog: n + 1 evaluations); `compensated = true`: update carried as value plus correction term (Hairer2008 Sect. 2) |
| `advance_to(m, flow, z, tau, ns; compensated = false)` | stored states after the step counts `ns` (increasing), one run |
| `advance_checked(m, flow, z, tau) -> (state, status)` | `status.ok::Bool`, `.iterations::Int`, `.increment`, `.stages` (implicit); explicit: ok, 0 iterations. A failed implicit step returns an all-NaN state and `ok = false`; nothing throws |
| `tangent_advance(m, flow, z, W, tau) -> (znew, Wnew)` | the deviation frame `W` (2n x k SMatrix or Matrix) through the elementary Jacobians in stage order (implicit: implicit differentiation); no renormalization |
| `IntegratorMap(m, flow, tau) <: AbstractMap` | `forward`, `inverse` (adjoint table at -tau), `tangent` (SMatrix, = `tangent_advance` with the identity), `inverse_tangent`, `dimension`; consumed by `Tangent.propagate` like any map |
| `advance_ensemble(m, flow, x0, tau, n; backend = :cpu, chunk = nothing) -> (states, ok)` | KernelAbstractions kernel over orbits; per-orbit status; result of orbit j independent of the chunking |
| `SplittingMethod(c, d; order, palindromic, source, name, digits = Inf)` | generic table (drift weights `c`, kick weights `d`, k stages, drift then kick, first stage first); `ArgumentError` at construction for k = 0, sum conditions off by more than 4 k eps_tab, a table declared palindromic without `c_i = c_{k+1-i}`, `d_i = d_{k-i}` |
| `coefficients(m) -> (c, d)`, `stage_count`, `nominal_order`, `is_palindromic`, `table_digits` (`Inf` for algebraic tables), `force_evaluations` (per single step), `method_name`, `method_source` (comma separated keys of `references.bib`), `adjoint_method(m)` | metadata. Palindromic tables end with `d_k = 0` and are stored with possibly zero end drifts (KDK: c = (0, 1, 0), d = (1/2, 1/2, 0)); S9 stored kick first (`c_1 = 0`, cumulative drifts are the gamma_i of Eq. (2.8)) |
| `SymplecticEuler(T, :DK or :KD)` S1, `Leapfrog(T, :DKD or :KDK)` S2, `McLachlan2(T)` S3, `Ruth3(T)` S4, `ForestRuth(T)` S5, `triple_jump(m, order)` S6 and `suzuki_jump(m, order)` S7 (composition builder, raise by 2 per level from `m`, default element type of `m`), `PEFRL(T)` S8, `CalvoSanzSerna4(T)` S9, `Yoshida(T, 6 or 8, :A or :D)` S11, `KahanLi(T, 6 or 8 or :s7odr6)` S12, `SofroniouSpaletta10(T)` S13 | explicit rows; every constructor also callable without `T` (Float64) |
| `SABA(T, n)`, `SBAB(T, n)`, `SABAC(T, n)`, `SBABC(T, n)`, `ABA864(T)`, `ABAH864(T)` | P rows; the corrected rows add the corrector kick before and after; `nominal_order` is the order at eps = 1 (2 for SABA/SBAB, 4 for SABAC/SBABC n >= 2 and ABA864/ABAH864) |
| `ImplicitMidpoint(T)` G1, `Gauss(T, s; maxiter = 50, tol = nothing, standard = false)` G2, `gauss_coefficients(T, s) -> (A, b, c)` | fixed-point solve to the Hairer2008 Sect. 3 stopping rule; `tol` replaces the acceptance bound 8 eps max(1, \|Y\|) by an absolute bound (used for 256-bit references); `standard = true` is the Hairer2008 Fig. 2.1 implementation (negative control, stage 04) |
| `TaoExtended(T, omega; order = 2)` G3, `extend_state(z)` = (q, p, q, p), `physical_state(zext)` | Strang formula, higher order by the triple jump with exponent 1/(l - 1); `omega` required |
| `taylor_reference(flow, z, t)` | declared (no methods) in Integrators, methods in the extension `ChaoslabTaylorExt` (weak dependency `TaylorIntegration`, R1) |

## Criterion to testset

Tags: S package suite, V verification environment (stage 04), G local GPU. Test-side constants in `ic_helpers.jl`: `ic_B` = B(n) of the spec, `ic_etab` = eps_tab, catalog `IC_S`, `IC_P`, `IC_PC`, `IC_GAUSS`, test-side polynomial step matrices for F1, test-side Hamiltonians and the BigFloat reference `ic_ref` (oracle O1).

| AC | Tag | Testset (file) | Sample | Tolerance |
|----|-----|----------------|--------|-----------|
| flows (P1) | S | FL1 to FL7 (ic_flows) | 5 separable flows + T1, 20 seeded points, 3 element types | 8 to 64 eps; Jacobians against ForwardDiff |
| 1 | S | AC1 | every explicit table, 256 bits, plus Float64/Double64 | 4 k eps_tab, bitwise symmetry, digits, source keys in `references.bib` and not in `unverified.md` |
| 2 | S | AC2 | every explicit table (S1 to S13, P rows), polynomial step matrix of F1 | 64 eps_tab S^p; higher linear order reported |
| 3(a) | S | AC3a | builders for orders 4, 6, 8 (triple jump 3, 9, 27 pieces; Suzuki 5, 25, 125), closed forms of Yoshida1990 Eq. (2.11) / Omelyan2002 Eq. (8) | 8 eps; polynomial equality to 1e-60 |
| 3(b), (b3), (b4), (b5) | S | AC3b, AC3b3, AC3b4, AC3b5 | Yoshida A against s7odr6, stage counts, z of S3 and the ratio, S13 digits/count/norm/slope 11, S9 RKN form against the splitting form | 1e-14; 5e-5; 0.6; 1e-5, 0.3; 16 eps |
| 3(c), 3(d), 3(b2) | V | none (stage 04) | O6 package solvers, O8 REBOUND arrays, McLachlan prints | max(16 eps, 1e-13); 1e-30; 1e-14, 1e-19 |
| 4 | S | AC4 | S2, S3, S5, S6 (4, 6, 8), S7 (4, 6, 8), S8, S9, S11, S12, S13 on F4 and F2, T = 8, tau0 chosen, J = 4 | slope within 0.2 (p <= 4) / 0.3 of p; 2 resolved halvings at least |
| 5 | S | AC5 | S1 (both orderings), S4 on F4 and F2; local exponent p + 1 for S1, S4, S9 on F4 and T2 in 256 bits | 0.2; 0.3 |
| 6 | S | AC6a, AC6b, AC6c, AC6 (ABA), AC6d | pendulum q = 1, p = 0.3, T = 8, eps 1e-8, 1e-3, 1e-2; SABA/SBAB n <= 4 (+ n = 6, 10 in 256 bits), ABA864, ABAH864, SABAC/SBABC n = 2..4 | 0.5 (1.0 for n = 6, 10); 0.3 and eps-slope 0.3 (n = 3, 4); ratio < 1; 1.0; ratio >= 10 and slope 0.5 / 1.5 |
| 7 | S | AC7 | G2 s = 1..4 on F4 and F2; coefficients; stage residual; iterations at h L rowsum = 0.1; divergence at 4; stopping rule | 0.2; 64 eps; 8 eps max(1, \|Y\|); <= 25; 100 percent; >= 90 percent |
| 8 | S | AC8a, AC8b,c | T1, omega = 10 and 100, T = 20; l = 2 and 4; negative control of the printed exponent | 0.2, 0.3; 1e-3; coupling monotone |
| 9 | S | AC9 | every method (explicit, corrected, Gauss s = 1..4, Tao l = 2, 4 in the extended space), F1 n = 3, F4, F2, 10 points per case, Float64 and Double64 | 8 d eps \|\|M\|\|^2 (+ 64 eps d implicit) |
| 10 | S | AC10 | every palindromic method n in {1, n_H}; adjoint property for the others; S1 adjoint; Gauss; Tao Strang | 16 B(n) (times 4 for M_n), 16 eps |
| 11 | S | AC11 | ForwardDiff (Float64, Double64) all methods on F4 and F2; n = 1, 10 step products; 256-bit central difference (2^-100) all methods and Gauss | 16 eps (1 + \|M\|) k; 2^-130 |
| 12 | S | AC12a, AC12b, AC12c | F1 exact exponential, all methods; F2 e = 0.5 tangent of S2 (tau = 2^-16) and S12-6 against finite differences of reference orbits (2 columns, 192 bits); frame against two-trajectory finite difference (S2, S12-6, S5) | B(n) + 10 truncation; 2^-30 relative |
| 13 | S | AC13 | S2, S5 to S9, S11 to S13 on F1 n = 3, F2 e = 0.1, 0.5, F4; P rows on F4 eps = 1e-2; 1e4 steps, 1000 samples (spec 1e5) | 2x first half + 16 eps M; slope x length <= 0.1 max |
| 14 | S | AC14 | every explicit method incl. corrected, F2 e = 0.1, 0.5, 1e4 steps (spec 1e5) | B(n) |
| 15 | S/V | AC15 | T2, H0 = 1/8, 100 starts (spec 200), S12 order 8, step 2/240, n up to 1e5 (spec 3e6), energy in 256 bits | package: slope without compensation within 0.2 (0.1 full), end ratio >= 10, mean within 3 SE; (a) and (d) full size only |
| 16 | S | AC16 | allocation of `advance` and `tangent_advance`, all explicit methods, F1 n = 1, 3 and F4, Float64 and Double64; counts through `CountingFlow` against metadata and the table; long-run form | 0 bytes; equality |
| 17 | S | AC17 | all explicit methods, F4 and F2, 2 orbits per method (spec 200), n_H steps; Gauss s = 1, 2; effective eps of Double64 | B(n) of the lower precision; 2^-100 |
| 18 | S/V/G | AC18 | S2, S5, S8, S12-6, SABA3, F4 and F2, 500 orbits (spec 1e5), CPU kernel against a loop, chunks 7 and 1000 bitwise; Gauss ensemble with a failing orbit. GPU against CPU and 1 against 8 threads: stage 04 / G | B(n); bitwise; statuses |
| 19 | S/V | AC19 | declaration of `taylor_reference`, weak dependency and extension file in the repo; the analytic checks run only when the extension is loaded (`@test_skip` otherwise) | 2^-180, 2^-100 |
| 20 | - | reserved (R2 deferred) | no test | |
| 21 | S/V | AC21 | method against R1 on F4, only when the extension is loaded | B(n) |
| 22 | S | AC22 | construction rejections, non-finite input (every explicit method, 4 flows; Gauss; Tao), ensemble with a NaN orbit | ArgumentError; non-finite output, no throw |
| 23 | S | AC23 | source parse of `src/Integrators` and `src/Systems`, position in `Chaoslab.jl`, negative control of the checker | empty set of upward imports |
| API | S | API | hand values (dyadic) for S1, S2; `IntegratorMap` as a map; `Tangent.propagate` consumption bitwise | exact |

Stage 04 only: criterion 3(c), 3(d), 3(b2); the GPU parts of 18; criterion 15(a) and (d) at full size (and the full sizes of 4, 9, 11, 13, 14, 17 and 18 through `IC_FULL_RUN = true` set by `verify/integrator-catalog.jl`); the Taylor oracle (R1, O7) under the verify environment.

## Deviations and spec items not testable as written

1. **API shape.** Method objects are built in an element type (`Leapfrog(T, :DKD)`), not converted at the step. Reason: the spec requires ForwardDiff derivatives of the stepper (O3) and decimal table constants "held in the element type"; a table converted inside `advance` could not follow a `Dual`.
2. **Spec 3(a)/(b) wording.** Closed forms are compared as kick weights and cumulative drift positions (invariant under merging of adjacent drifts) and as the F1 step polynomial; an implementation may store merged or unmerged tables.
3. **Criterion 2 for P rows.** The nominal order of a P row is its order at eps = 1 on F1 (SABA/SBAB 2, ABA864/ABAH864 4); SABAC/SBABC are not tables and are not in AC1 to AC3.
4. **Criterion 3(b5)** requires S9 to be stored kick first (`c_1 = 0`); the RKN form is written out on the test side from the stored table. Whether S9 is palindromic is not asserted (spec silent); reversibility is tested for whatever `is_palindromic` reports.
5. **Criterion 5.** Eq. (30) of Ruth1983 is not transcribed in the spec; the five nonlinear order conditions are tested through the local error exponent p + 1 in 256 bits on the pendulum and T2. The equality with ForestRuth1990 Eqs. (4.3), (4.4) is stage 04.
6. **Criterion 6.** (a) n = 6, 10 use tolerance 1.0 (transition to the eps^2 tau^2 branch at tau = 1/2, 1e-8 is close; spec 0.5). (b) n = 1, 2 are excluded from the eps^2 scaling: their leading error is eps tau^(2n) comparable with eps^2 tau^2 at tau <= 1/8, eps = 1e-2 (n = 1 slope in eps is 1). (d) n = 1 excluded (its remainder is eps tau^2); slope tolerance 0.5 for n = 2 and 1.5 for n = 3, 4 (transition towards eps tau^(2n)). ABA864/ABAH864: slope 4 branch at eps = 1e-2, tau in {1/16, 1/32, 1/64}, 256 bits, tolerance 1.0. All provisional (P17); reported with `@info`.
7. **Criterion 7.** The iteration bound "at most 25 at h L max sum|a| <= 0.5" cannot hold at the end of its range (a contraction factor 0.5 gains one bit per iteration; 53 iterations for binary64): it is tested at 0.1 (about 16 iterations), and reported. The reference orbits use `Gauss(...; maxiter = 500, tol = 1e-40)` because the default cap (50) and tolerance (8 eps at 256 bits) are too tight for a 256-bit reference (at step 1/64 on F2 e = 0.5 the contraction is 0.15). The spec's O1 for G2 itself is the Taylor reference; here G2 s <= 4 is compared with Gauss s = 8, an independent method of different order (the extension is not loaded in the package suite).
8. **Criterion 11.** Complex step omitted (the Kepler and pendulum flows use `abs`, `sqrt` and `sin` on real data; ForwardDiff plus 256-bit central difference are kept). The binary64/double-double bound 16 eps (1 + \|M\|) is multiplied by the stage count k (126 stages of Suzuki 8 round differently in ForwardDiff and in the tangent step; k is the metadata). The 256-bit bound is 2^-130 (1 + \|M\|) instead of "64 eps": the rounding of a central difference with step 2^-100 is 2^-155.
9. **Criterion 12.** Package size: two columns of the reference tangent (192 bits, Gauss s = 8 orbits at perturbation 2^-40 instead of 2^-60: the cost of a 256-bit orbit pair is about 5 s per orbit); tau = 2^-16 for S2 (Float64) and 1/128 for the 6th order table so that the truncation is below 2^-30 relative (my estimate, provisional); the truncation term of 12(a) is measured in double-double.
10. **Criterion 13** at 1e4 steps and 1000 samples (least squares slope noise needs about 1000 samples); S13 at tau = 1/50 is round-off dominated and may fail the "2x first half + 16 eps M" bound by chance (provisional, P17).
11. **Criterion 15** in the package suite: 100 starts, n up to 1e5, slope tolerance 0.2 and the compensated end-ratio only; the slope of (a) needs n >= 1e4 where the rms is at the evaluation floor, so (a) and (d) run at full size only.
12. **Criterion 18** needs the API `advance_ensemble`; the 1 against 8 threads comparison is tested through the chunking (the thread count is the run's); GPU parts are stage 04.
13. **Criteria 19 and 21** (R1): the extension arrives with a weak dependency. Needed in `Project.toml` (owner): `[weakdeps] TaylorIntegration` and `[extensions] ChaoslabTaylorExt = "TaylorIntegration"` plus `ext/ChaoslabTaylorExt.jl` (AC19 tests that they exist). The analytic checks run only when `Base.get_extension(Chaoslab, :ChaoslabTaylorExt)` is not `nothing` and are `@test_skip` otherwise; the verify environment loads the package and runs them.
14. **Criterion 16.** The force count is tested on the table methods (`IC_EXPL`); the corrector kick of SABAC/SBABC is a second force (of W) and is not counted by `kick_count`.
15. **Spec 8c/9 data.** Regression numbers used are the ones the spec prints with a source: 0.1932 and `z` formula (McLachlan1995 Table 2), 6.18133 and 35 stages (SofroniouSpaletta2005 Table 1, Eq. (17)), slope 11 (derived), digit counts of the tables (spec Algorithm item 3). `_config/bibliography/benchmarks.md` has no integrator row (Inconsistencies (g)); no number comes from memory.
16. Runge-Lenz sign convention (`A = p x L - q/|q|`) is fixed by AC test FL6 because the spec does not define it. T1/T2 start values are test choices (`IC_START`), recorded in `ic_helpers.jl`.

## Runtime

Package suite (D38, 3 minute ceiling). Reductions are marked "D38" in the test files and listed in the criterion table. Dominant costs: the 256-bit reference orbits (six: F4 and F2 at T = 8, the pendulum at eps 1e-8, 1e-3, 1e-2, T1 at T = 20), measured at 5 to 8 s each with a throwaway allocating stub of Gauss (a static-array implementation should be several times faster); AC15 (about 100 starts x 1e5 steps x 17 stages, twice); AC12b (4 orbits in 192 bits); AC13 (16 tables x 4 flows x 1e4 steps). The module-dependent runtime cannot be measured before stage 03; estimate 2 to 4 minutes. If stage 03 exceeds the ceiling, the first items to move to `verify/integrator-catalog.jl` are AC12b, AC15 and the 256-bit part of AC11 (same assertions, `IC_FULL_RUN`).

## Seeds

Base 20261012 + 41000 (range 41000..41999), `ic_seed(ac, sub) = base + 40 ac + sub` (sub < 40), StableRNGs only. Used: FL 0 (sub 1 to 6), AC3 (sub 1), AC7 (1), AC9 (1), AC10 (1), AC11 (1), AC15 (1), AC17 (1), AC18 (1), API 24 (sub 1, 2). Every draw comes from one of these generators; no unseeded randomness.

## Dry run of the helpers

A throwaway allocating stub in the session scratchpad (`scratchpad/ic02/stub.jl`: flows F2/F4, a splitting stepper, a Gauss solver; nothing of it is in the repository) was used to run the helpers `ic_steppoly`, `ic_rotpoly`, `ic_composepoly`, `ic_tj_weights`, `ic_ref`, `ic_pick_tau0`, `ic_scan`, `ic_slopes`, `ic_check_order`: leapfrog polynomial differs from the rotation first at tau^3, triple-jump polynomial equals the composition of fractions to 4e-77, measured slopes S2 on F4 and F2 2.00 and triple jump 3.83 to 4.00. The test bodies themselves were not run against an implementation.

## RED log

Command (driver `scratchpad/ic02/red.jl`: `using Test, Chaoslab` and the include of `test/integrator-catalog/runtests.jl`; environment `scratchpad/cm02/testenv`): `julia --project=<testenv> -t 2 red.jl`, 7 s wall. Result: 11 Pass, 3 Fail, 40 Error, 43 testsets. The 40 errors are `UndefVarError` on the first statement of each module-dependent testset (the line naming the new symbols: `HarmonicOscillator`, `SymplecticEuler`, `Gauss`, `SABA`, `coefficients`, `advance_n`, `advance_to`, `tangent_advance`, `IntegratorMap`, `TaoExtended`, `taylor_reference`, `SplittingMethod`, ... one per testset); no MethodError, no parse error, no error inside a helper, no include-time error. The 3 failures and 11 passes are all in AC23 (`isdir(src/Integrators)`, non-empty file list, position of `Integrators/Integrators.jl` in `src/Chaoslab.jl` fail; the parse of the existing Systems sources and the negative control of the checker pass). `test/layering.jl` passes (2 of 2). The full `Pkg.test()` was started in the background and killed after 20 minutes without a summary (all suites run before the summary prints), so it is not part of this evidence; `test/quality.jl` (Aqua, JET) was not run in this stage because the scratch test environment lacks both packages.

## Audit

| Check | Result | Evidence |
|---|---|---|
| Coverage of spec | PASS | criteria 1 to 19, 21 to 23 map to testsets (table above); 20 is a reserved placeholder; the V parts (3(c), 3(d), 3(b2), GPU, full sizes, R1 oracle) are named with their stage 04 home |
| RED | PASS | 40 of 40 module-dependent testsets error with `UndefVarError`; AC23 fails on the missing source directory; no parse or helper error |
| Reproducible | PASS | seeds listed above, StableRNGs only |
| No memory values | PASS | no regression number outside the spec's cited rows (deviation 15); coefficient identities are computed test-side in 256 bits |

## Amendment D112 (stage 02, 2026-11-06)

Spec: `integrator-catalog-spec.md`, criteria 3(b3), 4, 5, 6(a), 7, 8(c), 12, 13 (the changed sentences carry the mark D112). Owner decision "Amend tests" after the stage 03 run (impl notes, "Tests that do not hold": 119 test-side failures in the new suite and 3 in the analysis suite). Tests were fixed first, then rerun to GREEN; `src/` is unchanged (no test fix revealed a code defect). Files: `ic_helpers.jl`, `ic_flows.jl`, `ic_ac01_03.jl`, `ic_ac04_08.jl`, `ic_ac09_12.jl`, `ic_ac13_16.jl`, `test/analysis/an_p2_figures.jl`, `an_p3_dispatch.jl`. Before/after are the stage 03 run (7503 pass, 119 fail, 3 min 12 s) and the run after the amendment (7180 pass, 0 fail, 0 error, 2 broken = the extension skips of AC19 and AC21, 2 min 16 s).

| # | Item | Change | Why | Before / after |
|---|------|--------|-----|----------------|
| 1 | FL1, FL2 `:f1_3` (`ic_H`, `ic_gradV`) | `w^2` computed in the element type of the state (`T(w)^2`) | the Float64 square carried 1e-16 relative error into a Double64 or BigFloat reference that no implementation can meet | 60 failures (40 + 20) to 0: FL1 738 of 738, FL2 1000 of 1000 |
| 2 | AC3b3 (S3 against leapfrog) | the pendulum ratio is replaced by the ratio on F1 with n = 1 of the tau^3 coefficients of the step matrices (operator 2-norm, 256 bits), times 4 for equal work; asserted below 0.5; the pendulum value is no longer computed | the pendulum ratio depends on the start, the ordering and T (1.48, 0.93, 0.80, 0.51, 0.62, 0.58 measured); a linear problem has no such dependence | 1 failure to 0; ratio 1.48 (test configuration) to 0.318 (bound 0.5; McLachlan1995 Table 2 effective constants 0.026 and 0.070, 0.37) |
| 3 | AC4, AC5, AC7 start rule | `ic_asym`: the order is asserted on the two finest resolved halvings of a scan j = 0..8 (AC7: 0..5, implicit steps cost an iteration), at least three halvings resolved; earlier slopes only reported; orders 3 and 4 measured in double-double (`ic_slope_type`); same rule for every method | the first error below 1e-2 is outside the asymptotic range for the tables of order 6 and 8, S13 and S4 (first slopes 5.5 to 7.6, 11.08, 2.34); on F2 the slopes of S4 settle only near tau 2^-10 (4.1, 4.4, 5.0, 3.0, 2.1, 2.7, 2.87, 2.94) and a second-halving start would still fail | 20 failures to 0 (AC4 38 of 38, AC5 14 of 14, AC7 709 of 709); last two slopes within 0.01 to 0.2 of nominal, S13 10.01, S4 2.94 (F2), 2.99 (F4) |
| 4 | AC6a | n = 1, 2 unchanged (binary64, eps 1e-8); n = 3 in 256 bits at eps 1e-8; n = 4 in 256 bits at eps 1e-16; n = 6, 10: slope 2 within 0.3 at eps 1e-8 (the eps^2 tau^2 branch) plus the quadrature property of the kicks (moments k <= 2n - 1 to 64 eps, k = 2n above 1e-14) for n = 1 to 4, 6, 10 | binary64 for n = 3, 4 is at the round-off floor (3.2e-16, noise); at eps 1e-8, n = 4 and 6, 10 are in the eps^2 tau^2 branch even in 256 bits (a tau^12 is about 1e-20 eps; n = 6 at eps 1e-20 gives 12.0 then 9.6 then 1.8, below the reference resolution beyond that) | 8 failures to 0; n = 3: 6.0 (SBAB3 6.31), n = 4: 8.00, n = 6, 10: 2.04 to 2.05; moments at most 8.6e-78, k = 2n at least 1.5e-12 |
| 5 | AC8c | the law is asserted: slope 2 within 0.2 between tau 0.01 and 0.005 for omega 10 and 100, and the constants for omega 10 and 100 agree to 10 percent at tau 0.01; "decreasing with omega" dropped | the coupling is the Strang error c(omega) tau^2 with c 1.27 to 1.35 for omega 3 to 1000 (measured), not monotone in omega | 1 failure to 0; 1.272e-4 (omega 10) and 1.295e-4 (omega 100), +1.8 percent, slopes 2.00 |
| 6 | AC12b | S2: slope 2 within 0.3 between tau 2^-15 and 2^-16, and the Richardson extrapolate (4 W(tau/2) - W(tau))/3 within 2^-30; the 6th order table unchanged | the truncation of S2 is C tau^2 (3.9e-9 at 2^-16) against 2^-30 = 9.3e-10; tau 2^-18 would need 2.6 million steps | 1 failure to 0; 1.56e-8, 3.91e-9 (slope 2.001), extrapolate 4.0e-12; S12-6 1.17e-12 |
| 7 | AC13 | tables of order 6 or higher and ABA864, ABAH864 run in double-double (eps_T of that type), the others in binary64; same two assertions | in binary64 the energy error of 1e4 steps is the round-off random walk (1e-14 to 6e-14), not truncation-bounded oscillation; AC15 (Brouwer's law) covers that regime | 28 failures to 0 (164 of 164), 3.8 s |
| 8 | analysis tests | `an_p2_figures.jl` (D91n I29) and `an_p3_dispatch.jl`: weak dependencies CairoMakie, CUDA, TaylorIntegration; extensions the three; compat TaylorIntegration 0.18; hard dependencies, extras and targets exact as before | Project.toml gained the weak dependency and extension of AC19 | 3 failures to 0 (whole suite) |
| 9 | runtime | AC11 (ForwardDiff compile per method and flow) runs 15 representative methods in the package suite, all 41 explicit methods when `IC_FULL_RUN` is set (`IC_AC11_IDS`, same assertions) | one closure and tag per method: 0.3 s each | AC11 85 s (688 assertions) to 33 s (263); suite 3 min 12 s to 2 min 16 s against the 3 minute ceiling (D38) |

Weakening, stated: (2) the nonlinear equal-work check of S3 is gone (only the linear one is scored); (3) the early, pre-asymptotic slopes are reported, not scored; (4) the order 12 and 20 of SABA6, SABA10, SBAB6, SBAB10 is asserted through the quadrature property, not through an error slope (the slope is not observable); (5) the monotone decrease of the coupling in omega is not asserted; (7) the tables of order 6 or higher are not run in binary64 here; (9) AC11 covers a subset in the package run.
