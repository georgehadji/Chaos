# Test plan: henon-heiles

Status: approved at the stage 02 checkpoint, 2026-11-11 (owner: "Approve tests"), with the test-side readings U1 to U10 (section "Spec points found while writing the tests") as written; RED rerun independently: 32 pass, 26 error (all `UndefVarError` for proposed names), 0 fail.

Spec: `workspaces/module-build/stages/01-spec/output/henon-heiles-spec.md` (approved at the stage 01 checkpoint, D116, open points P1 to P10 as proposed). Suite: `test/henon-heiles/` (`runtests.jl` includes `hh_helpers.jl`, `hh_ac01_05.jl`, `hh_ac06_09.jl`, `hh_ac10_12.jl`, `hh_ac13_17.jl`; picked up by `test/runtests.jl` automatically). 30 testsets. Nothing under `src/`, `_config/` or `data/` was touched. Seeds: `20261012 + 42000 + 40 n + sub` (`hh_seed`, n the criterion), generator StableRNGs. Package sizes follow D38; `verify/henon-heiles.jl` sets `HH_FULL_RUN = true` before including `runtests.jl` for the full sizes of the spec.

## API proposal for stage 03

Reused as they are: `HenonHeilesFlow` (promoted in place to F3, P1; docstring and `src/Systems/CONTEXT.md` lose "T2 test-only"), `hamiltonian`, `force`, `force_jacobian`, `vector_field`, `advance`, `advance_n`, `tangent_advance`, `IntegratorMap`, `propagate`, `propagate_orbit`, the estimators `QRSpectrum`, `Benettin`, `SALI`, `GALI`, `eps_pair`, `eps_sum`, `sum_threshold`, `initial_deviations`. Element types: every new function is generic in the element type (Float64, Double64, BigFloat); states are `SVector{4}` (q1, q2, p1, p2); methods and step sizes are built in the element type (D113).

Systems (new, in `src/Systems/`):

| Name | Meaning |
|------|---------|
| `critical_points(flow::HenonHeilesFlow, ::Type{T})` | `NTuple{4, SVector{2,T}}`: O, S0, S1, S2 = (0,0), (0,1), (+sqrt(3)/2, -1/2), (-sqrt(3)/2, -1/2), correctly rounded (Equation 4) |
| `escape_energy(flow, ::Type{T})` | `T(1)/6` |
| `d3_elements(::Type{T})` | `NTuple{6, SMatrix{2,2,T,4}}`: identity, rotation by 2 pi/3, rotation by 4 pi/3, reflection diag(-1, 1), then the two rotations times the reflection (Equation 5); entries correctly rounded |
| `section_interval(flow, E)` | `(a, b)`, the roots of u(q2) = E bounding the origin component of the section line (Equation 7); `(-1/2, 1)` at E = 1/6 |
| `energy_start(flow, E, q2, p2)` | `(z, status)`: z = (0, q2, p1, p2) with p1 from Equation 8, status `:ok`; `(nothing, :outside)` if p2^2/2 + u(q2) > E (or, for E above 1/6, r > 1); `(nothing, :nonfinite)` for a non-finite argument. Never a NaN state |
| `random_energy_start(flow, E, q2window, p2window, rng; maxtries = 10_000)` | rejection draw (Algorithm 2) inside the origin component for E at most 1/6, inside r <= 1 above; `(z, :ok)` or `(nothing, :no_admissible_point)` (`:nonfinite` for a non-finite E) |
| `energy_grid(flow, E, n)` | `Vector{SVector{4,T}}`: cell centres of the n x n grid over [a, b] x [-sqrt(2E), sqrt(2E)], kept where inside D(E), p1 from Equation 8; E at most 1/6; deterministic (the portrait starts, Algorithm 4; P2 puts the classification with Experiments or the study, so no portrait function is proposed) |

Integrators (new, in `src/Integrators/`; section and escape tools, P2):

| Name | Meaning |
|------|---------|
| `locate_crossing(m, flow, z, tau)` | for z with q1 < 0 <= q1 of the next step: `(zc, s)`, the root s in (0, tau] of q1(Phi_s z) = 0 by the method's own step (Equation 11), `zc = advance(m, flow, z, s)`, `abs(zc[1]) <= 16 eps`; allocation free for isbits states |
| `section_map(m, flow, z0, tau, n; jacobian = false, max_crossings = typemax(Int))` | record with fields `t::Vector{T}`, `z::Vector{SVector{4,T}}` (full crossing states, q1 about 0, p1 > 0), `J::Vector{SMatrix{2,2,T,4}}` (empty unless `jacobian`; `J[i]` is the Jacobian of the leg from section point i-1 to i, Equation 12; `J[1]` the identity when the start lies on the section), `grazing::Int` (bracketed crossings with abs(p1) < 8 eps at the root, flagged and not recorded), `status` in `:ok`, `:nonfinite`, `:buffer_full`. A start with q1 = 0 and p1 > 0 is the first point at t = 0 (Algorithm 3). Never throws on non-finite values |
| `escape_sector(flow, q)` | `Int` in 0, 1, 2: the exit with centre angle pi/2, 7 pi/6, 11 pi/6 nearest to the polar angle of q (Equation 14); allocation free |
| `locate_escape(m, flow, z, tau, R)` | for z with r <= R < r of the next step: `(zc, s)`, root of r(Phi_s z) = R |
| `escape_orbit(m, flow, z0, tau, n; R_esc = 10)` | `(t, sector, status, z)`: first passage of r = R_esc located as above (`status = :escaped`, `z` the state at r = R_esc, `sector` of Equation 14); `:bound` (no escape in n steps; `t` NaN, `sector` -1); `:nonfinite` (non-finite state or overflow). Never throws |

Spectrum, GALI and SALI of a flow orbit use the existing reducers on `IntegratorMap`; the existing estimators report exponents per step (the sum of the log growth over the step count), so the tests divide by the step size to compare with the per-time statements of the spec (spec point U4).

## Criterion to testset

Tags: S package suite, V verification environment (stage 04), G local GPU. "Existing" = passes on the code as it is today (listed below).

| AC | Tag | Testset (file) | Sample | Tolerance |
|----|-----|----------------|--------|-----------|
| 1 | S | AC1a (existing), AC1b (`hh_ac01_05`) | O, S0, S1, S2 in 256 bits; module points; benchmarks.md points (+-sqrt(3)/2, +1/2) | 2^-200; gradient norm at least 1; module data within eps |
| 2 | S | AC2 (existing) | 1000 seeded points in [-1.5, 1.5]^2, Float64 and Double64; ForwardDiff both, complex-step Float64 | 8 eps max(1, abs(entry)) |
| 3(a) | S | AC3a | 6 elements of `d3_elements` at 1000 points, 256 bits, and binary64 | 8 eps max(1, abs(V)) (binary64 adds 4 eps abs(q)_1 abs(grad V)_1, spec point U6) |
| 3(b) | S | AC3b (existing) | order-8 table, 1000 starts | bitwise |
| 3(c) | S | AC3c | 20 starts at E = 1/8, n = 1000, rotation in 256 bits | 16 B(n) exp(lambda_1 n tau), lambda_1 from Benettin on the orbit |
| 3(d) | S | AC3d (existing) | order-8, order-4, leapfrog KDK; 5 regular starts (E = 1/20) each, n = 10^4 | B(n) (regular orbits only, spec point U10) |
| 4 | S | AC4a, AC4b, AC4c | 10^5 starts at E = 1/12, 1/8, 1/6; 1000 x 1000 grid; 1000 infeasible trials per E (four energies) | 4 eps max(1, E); p1 > 0; domain; area relative 1e-3; status 100 percent |
| 5 | S | AC5 | E = 1/12, 1/8, 1/6, 40 starts each (full 200), T = 10^3, order-8, tau = 1/20 | r <= r_b (1 + 1e-6), abs(p)^2 <= 2E (1 + 1e-6), `:bound` |
| 6 | S | AC6 | E = 1/8, 10 starts (full 200), T = 10^3, tau, tau/2, tau/4, both tables | 1e-8 (order-8), 1e-4 (order-4); scaling reported |
| 7 | S | AC7a, AC7b | E = 1/12, 1/8, 1/6, 4 starts each (full 20), T = 100 (full 1000), Double64 order-10 reference at tau/4 | I6 invariants for all crossings; count equal and location 1e-8 up to T_AGREE = 60 (spec point U1); grazing flagged |
| 8 | S | AC8 | first three legs of 2 starts per E (full 10), reference central differences, step 1e-6, Double64 | 1e-6 relative to max(1, abs(entry)); det residual from `det_residual_tol` |
| 9 | S | AC9 | orbit C (E = 1/8), 1000 consecutive crossings | log-log slope of the running defect at most 1.2, signed and absolute sums (spec point U5) |
| 10 | S (24 x 24), V (200 x 200) | AC10 | E = 1/8, T_f = 2000, SALI < 1e-8; orbits R, S, C of Equation 9 (p1 rebuilt, printed p1 within 1e-6) | 55 percent +- 10 points; SALI(R) in [0.5, 2]; SALI(C) < 1e-8 and at most 1e-12 by t = 10^4; S reported |
| 11 | S (14 x 14), V (100 x 100) | AC11 | E = 1/12, 1/8, 1/6, order-8 and order-4 tables, T_f = 2000 | gaps at least 5 points; tables differ by at most 3 points |
| 12(a) | S | AC12a | E = 0.2, 0.25, 100 starts each (full 1000), T = 200, R_esc = 10, continuation to 2 R_esc | r at the flag at least R_esc (1 - 16 eps); sector equals the test-side polar-angle rule; 0 returns below R_esc / 2; default R_esc = 10 |
| 12(b) | S | AC12b | 8 escapes per energy before t = 100 (full 100), Double64 reference | 1e-6 |
| 12(c) | S | AC12c | all starts; reflection bitwise (label k to (-k) mod 3, same t); rotation by 2 pi/3 for escapes before t = 50 | exact labels; at least 90 percent of early escapes escape again |
| 12(d) | S | AC12d | NaN, Inf, overflow starts for `escape_orbit`, `section_map`, `random_energy_start` | status, no exception |
| 13(a)-(c) | S | AC13a to c | orbit C and a seeded chaotic start at E = 1/6 (first of 200 draws with lambda_1(500) > 0.02), T = 2 10^4 (full 10^5), RGS and MGS | `eps_pair` and `eps_sum` against `sum_threshold` plus 10 ln(T)/T allowance (spec point U8); abs(lambda_2), abs(lambda_3) <= 10 ln(T)/T |
| 13(d) | S (package reference), V (Taylor, ChaosTools) | AC13d | orbit C, T = 100, seeded Haar frame, Double64 order-10 reference at tau/4 | 1e-6 (spec point U7) |
| 13(e) | S (30 starts, T = 2000), V (100, 10^4) | AC13e | E = 1/8; tau 1/20 against 1/40; order-8 against order-4 | 3 standard errors |
| 14 | S | AC14 | S0 at T = 5 (closed form), T = 100 spectrum; O at T = 10 and 10^3; points from `critical_points` | 1e-9; 0.05 on sorted exponents; B(n); 1e-12 |
| 15 | S | AC15a, AC15b | R: GALI_2..4 at 10 times in [100, 2000], two Haar frames; C: 100 times in [10, 1000], two frames | slope within 0.25 of 0, -2, -4 (values above 10 floor); rate within 25 percent of E_k times the exponents of the same stretch (spec points U2, U3) |
| 16 | S, G | AC16 | E = 1/8 portrait 10 x 10 (full 40 x 40), default chunk against chunk 7 against the serial `propagate_orbit` map | bitwise states and labels; the GPU part is stage 04 (G) |
| 17 | S | AC17 | `advance`, `tangent_advance`, `locate_crossing`, `locate_escape`, `escape_sector`; Float64 and Double64 | 0 bytes after warm-up |

Coverage: all 17 criteria map to at least one testset. Stage-04-only parts: the TaylorIntegration and ChaosTools spectra of AC13(d) (O3, O4; the package environment has neither; a Double64 order-10 reference replaces them in the package suite), the GPU comparison of AC16, and the full-size grids and horizons (flag `HH_FULL_RUN`).

Regression values: only the escape energy 1/6 (benchmarks.md "Henon-Heiles escape energy", Zotos2017), the section fraction 55 percent and the printed p1 of orbits R, S, C (spec Equation 9 and O5, AntonopoulosManosSkokos2010 Sect. 3.1, read per the spec), and the closed forms derived in the spec (Equations 4 to 8, 16). No value from memory. The benchmarks.md row already carries the corrected ordinate -1/2 (D116, 2026-11-10), so Inconsistency (a) of the spec is resolved in the row.

## Tests that pass already (existing code)

AC1a (14 tests: the critical points, levels and Hessian eigenvalues from `force`, `hamiltonian`, `force_jacobian` of the existing flow), AC2 (force and Hessian against ForwardDiff and complex-step in Float64 and Double64), AC3b (reflection commutes bitwise with one order-8 step), AC3d (time reversal returns within B(n) for three palindromic tables). They guard the existing flow and become regression tests of the promotion (P1).

## Spec points found while writing the tests

- U1 (criteria 7, 8). A pointwise agreement of 1e-8 with the reference to T = 10^3 cannot hold for a chaotic start: two integrators differing by 1e-14 diverge like exp(lambda t) (lambda about 0.1 at E = 1/6). The tests compare count and location up to T_AGREE = 60 and check only the invariants (I6) beyond it.
- U2 (criterion 15). The regular-orbit laws (G_3 of order t^-2, G_4 of order t^-4) need generic deviation vectors: with the frame (1, 0, 0, 0), (0, 1, 0, 0), ... the drafting run on the existing code gave slopes -0.85 and -1.75 for orbit R over t in [100, 3200] (two tangent-like vectors, spec B7 with m > 0), with seeded Haar frames -1.9 to -2.0 and -3.8 to -3.9. The tests use seeded frames.
- U3 (criterion 15). E_k times lambda_1 of the long run (0.039 to 0.041 for C) is not the rate of GALI_k on the early resolved stretch of the sticky orbit C (G_3 decays at 0.044 to 0.053 against 2 lambda_1 = 0.08). The tests compare with the exponents of the same stretch (ratios 0.96 to 1.23 in the drafting run, bound 25 percent) and report the literal reference.
- U4 (criteria 13, 14). The existing estimators return exponents per step (`QRSpectrum`, `Benettin`: sum of the log growth over the step count), not per time; for tau = 1/20 the factor is 20. The criteria are per time; the tests convert.
- U5 (criterion 9). "Cumulative |det J - 1|" is read as the running sum, in two forms: signed (slope 0.5 for a random walk to 1 for a bias) and absolute (slope 1 by construction); both against 1.2. A product of 10^3 Jacobians of a chaotic orbit has no usable determinant (entries about e^300).
- U6 (criterion 3a). In binary64 the rounded matrix and product add 4 eps abs(q)_1 abs(grad V)_1 to the bound; in 256 bits the bound is the spec's as written.
- U7 (criterion 13d). Taylor and ChaosTools references are stage 04; the package check uses the same estimator, frame and renormalization time on the Double64 order-10 table at tau/4 (agreement 1e-6).
- U8 (criterion 13a). The finite-time allowance of the pairing has no default (diagnostics-spec O2); the tests use 10 ln(T)/T per time unit, the shape of criterion 13(b).
- U9 (placement). Criteria 10, 11 and 16 need a portrait of grids; P2 leaves it with Experiments or the study. The tests build it from `energy_grid` and the existing reducers, so no Experiments function is proposed; stage 03 may add one without changing the tests.
- U10 (criterion 3d). B(n) is the round-off bound for regular orbits; at E = 1/8 the chaotic starts would amplify it exponentially, so the test uses E = 1/20 (regular).
- Drafting numbers (existing code, order-8 table KahanLi at tau = 1/20): energy error 5e-15 over T = 10^3 (order-4 table 1.9e-7, 7.8e-10 at tau/4); S0 tangent error at T = 5 is 3e-13; exponents at O are 7e-16 per time; SALI of R at 2000 is 1.11; of C 1.9e-12 at 1000 and 1.3e-15 at 2000; of S 0.69 at 2000 and 0 at 10^4; chaotic fractions at E = 1/12, 1/8, 1/6: 0.0, 0.52, 1.0 (14 x 14) and 0.004, 0.539, 0.968 (24 x 24), order-4 table within 1.3 points. These are drafting measurements, not regression values of the tests.

## RED log

Run: `julia --project=<scratch env05> run_hh.jl <repo>` from the repo root (`using Test, Chaoslab, StaticArrays, StableRNGs; include(test/henon-heiles/runtests.jl)`), whole new folder only, Julia 1.13.2, 1 thread, 6.2 s after precompilation.

```
Test Summary:   | Pass  Error  Total  Time
henon-heiles    |   32     26     58  3.7s
```

30 testsets: 4 pass (AC1a, AC2, AC3b, AC3d: 32 tests), 26 error, 0 fail, 0 broken. Every error is an `UndefVarError` for a proposed name, raised by the first line of the testset (the tuple naming the new symbols): `random_energy_start` 12 testsets, `energy_start` 5, `energy_grid` 3, `critical_points` 2 (AC1b, AC14), `d3_elements` 1, `section_interval` 1, `section_map` 1, `escape_orbit` 1. No parse error, no error in a helper, no failure of an existing-code test.

Check of the test logic (not RED): the same suite against a throwaway implementation of the proposed API (scratch only, not delivered, bisection crossings, 1 thread) gives 6898 pass, 2 fail, 0 error in 61 s; the 2 failures are the allocation of the scratch `escape_sector` (it builds a vector), a defect of the scratch code. The Jacobian convention of AC8 and the section-map tests agree with the Double64 reference on that implementation.

## Amendment D117 (2026-11-11)

Spec section "Amendment D117". Test changes: `hh_ac06_09.jl` AC8, M includes `maximum(abs, J)` (item 1); `hh_ac10_12.jl` AC12b, escapes before t = 50 scored, reference horizon 51 (item 2); new testset "D117 boundary start and the 1-2 label tie" in `hh_ac10_12.jl` (item 3: `energy_start` at p1 = 0 is `:outside`; `escape_sector` at q1 = +0 and -0 gives 2 and 1). `verify/henon-heiles.jl` O3b and O3c follow items 1 and 2.
