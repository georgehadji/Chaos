# Test plan: standard-map

Spec: `../../01-spec/output/standard-map-spec.md`. Suite: `test/standard-map/` (`runtests.jl` includes `sm_helpers.jl`, `sm_ac01_06.jl`, `sm_ac07_13.jl`, `sm_ac14.jl`, `sm_y5.jl`, `sm_ac18.jl`). API names from the Batch A contract only; helpers are prefixed `sm_`. Seeds: `StableRNG(20261012 + offset)` with the offset written in each testset (1000 x criterion + index; 14001, 15001, 16001, 17001, 18001 for the later ones). `eps_T`: Float64 2^-52; Double64 2^-100 (see Spec ambiguities). All BigFloat work is in `setprecision(BigFloat, 256) do ... end`. No Float32 tests (D24). The outer `@testset "standard-map"` lets the file run on its own without a failing criterion aborting the rest.

## Criterion to testset

| AC | Testset | Sample sizes | Tolerances |
|----|---------|--------------|-----------|
| 1 | "AC1 tangent vs derivative" | Float64: 1000 points per listed K (17) plus 1000 with K uniform in [-20,20]; Double64: 200 per K plus 200 uniform | each DT entry vs ForwardDiff (Float64, Double64) and complex step h=1e-30 (Float64 only): 8 eps (1+abs K) |
| 2 | "AC2 determinant and symplecticity" | Float64 1000 per K (+ uniform); Double64 200 per K (+ uniform) | `det_residual`, `symplectic_residual`, and a plain recomputation of both: 8 eps (1+abs K) |
| 3 | "AC3 inverse round trip" | 200 orbits per listed K, Float64 and Double64; n in {1, n_H(K)}; both orders | plain difference <= tau(n) = 16 eps M_n n Gamma_n (M_n on the orbit of the first map applied) |
| 4 | "AC4 inverse tangent" | as AC2 | max entry of DT(z) DT_inv(forward z) - I <= 8 eps (1+abs K)^2 |
| 5 | "AC5 lifted vs wrapped" | 200 per listed K, Float64 and Double64; wrapped run starts from `wrap(z0)` and applies `wrap` after every step; checked at every k <= n_H | torus distance <= tau(k) |
| 6 | "AC6 K to -K symmetry" | 200 per listed K, Float64 and Double64, `exact(T, pi)`; every k <= n_H | torus distance <= tau(k) |
| 7 | "AC7 K=0 exactness" | orbits: 1000 (Float64), 200 (Double64), grid data 2^-10, 1000 steps; tangent products: 20 starts each, n in {1,2,10,100,1000,10^4} | exact equality (`==`) |
| 8 | "AC8 K=0 finite-time exponent" | Float64 and Double64, n in {1,10,100,10^3,10^4,10^6}; product of one-step tangents followed along the orbit up to 10^4 steps, then the (state-independent) K=0 tangent is reused | abs error vs 256-bit ln(sqrt(n^2+1))/n <= 8 eps; abs(10^4 X/ln 10^4 - 1) <= 1e-9; w0=(1,0) gives X == 0 |
| 9 | "AC9 fixed points and eigenvalues" | K in {-6,-2,1,3,4.5,6} at (0,0) and (pi,0), Float64 and Double64; Jordan powers n in {10,100,1000} for K=4 at (pi,0) and K=0 | class from the exact trace; trace 8 eps (1+abs K); abs lambda_+ or angle 32 eps (1+abs K); max entry of DT^n / n = 4 or 1 within 1e-12 |
| 10 | "AC10 fixed point is exact" | every listed K, 10^4 steps, Float64 and Double64 | exact equality |
| 11 | "AC11 benchmark trace row" | K in {1,3,3.9,4.1,5} (`exact(T, "3.9")`), x = 2 pi (1/2) built in T, Float64 and Double64 | trace 2-K at (pi,0) and 2+K at (0,0): 8 eps (1+abs K); class elliptic for K<4, hyperbolic for K>4. Benchmark row "Stability of the central fixed point at (1/2, 0)" |
| 12 | "AC12 Datseris eigenvectors at K=1" | one matrix, Float64 | derived ratios (s-K)/2 and -(K+s)/2 to 1e-12; eigen equation residual 1e-12; printed 0.62 within 0.005. Benchmark row "Eigenvectors of the fixed point (0,0)" |
| 13 | "AC13 Contopoulos orbit data" | O1 (2 steps), P1 (18 steps, and not 9), K=5, Float64 | torus distance <= 2 pi 1e-4 (and > 2 pi 1e-4 at 9 steps). Benchmark row "Island centre and unstable 4/9 orbit point" |
| 14 | "AC14 element-type agreement" | 200 orbits per listed K, n = n_H(K); Float64 vs 256-bit, Double64 vs 256-bit; also Chaoslab 256-bit vs an independent plain-formula 256-bit orbit | lifted and torus difference <= tau(n_H); Double64 uses the measured effective eps = max(2^-104, largest relative sin/cos error on the orbit x values), pass requires it <= 2^-100 |
| 15 | "AC15 Y5 symplecticity" | N in {1,2,3,5}, 10 draws each, Float64 and Double64, A and A_S (B symmetric in [-1,1]); also `dimension`, `tangent == A`, `symplectic_form` block form | max abs(A' Omega A - Omega) <= 16 eps normF(A)^2; abs(det A - 1) <= 16 eps normF(A)^N (module residuals and plain recomputation) |
| 16 | "AC16 Y5 exact exponents" | N in {1,2,3,5}, 10 draws each, n in {10,50}, Float64 only (sigma_max from an SVD) | 16 eps for block-diagonal A; 16 eps + ln(kappa(S))/n for A_S; A^-1 from `LinearSymplecticMap(A).Ainv` |
| 17 | "AC17 Y5 inverse" | N in {1,2,3,5}, 10 draws each, Float64 and Double64 | A^-1 A, A A^-1 vs I, `Ainv` vs inverse blocks built with `rotation_block(T,-theta)` and `hyperbolic_block(T,-mu)`, `inverse(forward z)`, `inverse_tangent`: max entry <= 16 eps normF(A)^2 |
| 18 | "AC18 map-level reproducibility" (map-level part only) | 10^5 seeded orbits, K in {1,5}, n = n_H(K), Float64: serial loop vs `Threads.@threads` loop, `isequal`; zero allocation of `forward`, `inverse`, `tangent`, `inverse_tangent` (Float64, Double64) | bitwise equal; `@allocated == 0` after warm-up inside a function |
| 19 | stage 04 | | ChaosTools cross-check in the separate environment |
| scope | "SC1 non-finite input" (extra, from spec Scope) | NaN, Inf, -Inf in x and in p | non-finite output, no exception |

## Deferred to stage 04 or tangent-core (no test here)

- AC18, CPU-backend kernel vs plain loop and 1 vs 8 threads of `propagate`: needs `Tangent.propagate`; belongs to the tangent-core suite (thread invariance, zero allocation of `propagate_orbit`) and stage 04.
- AC18, GPU vs CPU on 10^5 orbits and the Double64 GPU compile gate with reason "kernel does not compile for the element type": stage 04 (local GPU).
- AC14 heavy part (O7: larger BigFloat sample; effective Double64 sin/cos accuracy on a larger set): stage 04. The package test already measures the effective eps on the 200-orbit sample.
- AC1/O6 ForwardDiff of the n-step map vs the product of tangents at 10^4 points: stage 04.
- AC19 (ChaosTools standard map): stage 04.
- Float32 items of AC14 and AC18: dropped by D24.

## Spec ambiguities resolved

1. Double64 eps. The spec's nominal 2^-104 cannot be met by double-double arithmetic: with DoubleFloats, `symplectic` residuals of DT in Double64 reach about 1.4 x the 8 eps (1+abs K) bound at nominal eps (K = 20), and `sin`/`cos` measure 7.2e-31 (about 2^-100.4) relative on a 2000-point sample. All Double64 tolerances therefore use eps = 2^-100, the limit criterion 14 already accepts for the effective value. AC14 still measures the effective value (floor 2^-104) and requires it <= 2^-100. Owner should confirm.
2. "1000 points per check": read as 1000 per listed K for Float64. For Double64 the sample is 200 per K (3600 points in total, still at least 1000) to keep runtime down.
3. "Wrapped-evolution orbit" (AC5): built from the API as `wrap(forward(sys, .))` after every step, starting at `wrap(z0)`.
4. ||A|| in AC15 and AC17 is the Frobenius norm (spec leaves the norm open; it bounds the 2-norm from above).
5. AC16 is Float64 only: sigma_max(A^n) needs an SVD, which is not available for Double64 in the test environment.
6. AC13 coordinate order: unit-torus pairs read as (x_u, y), mapped to (2 pi x_u, 2 pi y).
7. Extra scope item SC1 asserts that Inf input to `forward` gives a non-finite result instead of the DomainError that `sin(Inf)` raises; the implementation has to guard it.

## Runtime

Measured against a scratch stub of the Systems API (outside the repo, not committed): about 105 s in total on this machine with 2 threads, most of it compilation (Y5 with 10x10 `SMatrix` and Double64 about 55 s, AC1 9 s, AC14 8 s). The compile cost of the real implementation adds to this; the sample sizes are already reduced for Double64 as listed above. Steady-state evaluation is a few seconds. The stub run passed all 541 assertions, which validates the helpers (tolerance ratios: AC2 Float64 0.80, AC4 0.27, AC1 complex step 0.14, others below 0.1).

## RED log

Command (repo root): `julia --project=<scratchpad>/m2/testenv -e 'using Test, Chaoslab; include(raw"test/standard-map/runtests.jl")'` (wall time 49 s). Result: 0 Pass, 0 Fail, 107 Error; every error is an `UndefVarError` on a missing Chaoslab name (the module `Chaoslab` is empty). Counts are `@test` records plus one record per testset aborted by an exception.

| Testset | Error/Fail | First error |
|---------|-----------|-------------|
| AC1 tangent vs derivative | 1 Error | UndefVarError `StandardMap` |
| AC2 determinant and symplecticity | 2 Error | UndefVarError `StandardMap` |
| AC3 inverse round trip | 2 Error | UndefVarError `StandardMap` |
| AC4 inverse tangent | 2 Error | UndefVarError `StandardMap` |
| AC5 lifted vs wrapped | 2 Error | UndefVarError `StandardMap` |
| AC6 K to -K symmetry | 2 Error | UndefVarError `exact` |
| AC7 K=0 exactness | 4 Error | UndefVarError `StandardMap` |
| AC8 K=0 finite-time exponent | 1 Error | UndefVarError `StandardMap` |
| AC9 fixed points and eigenvalues | 1 Error | UndefVarError `exact` |
| AC10 fixed point is exact | 1 Error | UndefVarError `StandardMap` |
| AC11 benchmark trace row | 1 Error | UndefVarError `exact` |
| AC12 Datseris eigenvectors at K=1 | 1 Error | UndefVarError `tangent` |
| AC13 Contopoulos orbit data | 1 Error | UndefVarError `StandardMap` |
| AC14 element-type agreement | 2 Error | UndefVarError `StandardMap` |
| AC15 Y5 symplecticity | 1 Error | UndefVarError `block_diagonal` |
| AC16 Y5 exact exponents | 1 Error | UndefVarError `block_diagonal` |
| AC17 Y5 inverse | 80 Error | UndefVarError `block_diagonal` |
| AC18 map-level reproducibility | 1 Error | UndefVarError `StandardMap` |
| SC1 non-finite input | 1 Error | UndefVarError `StandardMap` |

No parse errors, no errors inside helpers.

## Amendment D28 (2026-10-13)

Test file `test/standard-map/sm_wrapped.jl` (included last from `runtests.jl`, inside the outer testset). One testset per criterion W1 to W8 of the spec amendment. Helpers carry the prefix `smw_`; they reuse `sm_tau`, `sm_nH`, `sm_tdist`, `sm_absdiff`, `sm_orbit`. Seeds: `SM_SEED` plus 5001 to 5062 (distinct per criterion and element type). K in {0.5, 2, 5, 10}.

| Criterion | Testset | Sample sizes | Tolerance |
|---|---|---|---|
| W1 torus agreement | W1 torus agreement with the lifted map | 200 orbits per K, n = 1..n_H(K), Float64 and Double64 | torus distance <= tau(n); also x, p in [0, 2 pi), m, n integers |
| W2 lifted reconstruction | W2 lifted reconstruction | same orbits as W1; torus_state and lift checks on 100 seeded points in [-100, 100]^2 and one known point | lift at 256 bits vs lifted orbit <= tau(n); lift(torus_state(z)) - z <= 16 eps(100) |
| W3 counts exact and bounded | W3 counts exact and bounded | 1 orbit, 10^6 steps, K = 5; first 50 steps teacher-forced; free run n_H(5) = 6 steps | range, isinteger, below 2^53; increments equal; torus <= tau(n) and counts equal |
| W4 K = 0 | W4 K = 0 known answer | 20 seeded dyadic orbits (i / 1024, i in 0:6433), k <= 1000 | p, n unchanged; m_k equal to the 256-bit floor |
| W5 tangent agreement | W5 tangent agreement | 100 (Float64) and 50 (Double64) orbits per K, n <= n_H(K) | AC4 bound against the 256-bit product at the same points; vs lifted product: 2 x bound + n abs(K) tau(n) Gamma_n sqrt 2 |
| W6 inverse round trip | W6 inverse round trip | 200 orbits per K, n in {1, n_H(K)}, both orders, Float64 and Double64 | torus <= tau(n); counts exactly 0 |
| W7 zero allocation | W7 zero allocation | forward, inverse, tangent, inverse_tangent, propagate_orbit (100 steps), Float64 and Double64 | 0 bytes |
| W8 overflow guard | W8 overflow guard | m and n at maxintfloat - 1 (forward), m at -(maxintfloat - 1) (inverse), m two below the limit | NaN in all four components; NumericalHealthError "state" from propagate_orbit and propagate; exact count below the limit |

Deviation from the orchestrator's W3 wording (design point that does not hold): the first 50 steps cannot be compared by a free-running 256-bit orbit, because the orbit is chaotic. Measured at K = 5 from (1, 2): torus distance between the binary64 and 256-bit orbits is 2.7e-12 at step 10, 7e-7 at step 20, 6e-3 at step 30, 1.6 at step 40; first distance above 1e-3 at step 27, first count mismatch at step 34. So W3 recomputes the count increments at 256 bits from the exact binary64 state of each of the first 50 steps (teacher-forced), and compares the free-running 256-bit orbit only up to n_H(5) = 6 steps.

Worst observed ratios (error over bound, 1 = bound): W1 0.026 (Float64), 2.4e-4 (Double64); W2 0.024, 3.7e-4; W5 against 256 bits 0.069, 2.3e-3, against lifted 1.1e-3, 9.9e-6; W6 0.063, 2.8e-3; W3 free run 2.3e-3. W4: the closest approach of any x0 + k p0 to a multiple of 2 pi over the 20 seeded orbits is 4.9e-4, far above the 1e-12 scale where a floor could differ.

### RED log (Amendment D28)

Run: helpers plus `sm_wrapped.jl` only, before any source change. 13 Error, 0 Pass, 0 Fail; no parse errors, no errors inside helpers.

| Testset | Result | First error |
|---|---|---|
| W1 torus agreement with the lifted map | 1 Error | UndefVarError `WrappedStandardMap` |
| W2 lifted reconstruction | 3 Error | UndefVarError `WrappedStandardMap` (two), `torus_state` |
| W3 counts exact and bounded | 1 Error | UndefVarError `WrappedStandardMap` |
| W4 K = 0 known answer | 1 Error | UndefVarError `WrappedStandardMap` |
| W5 tangent agreement | 1 Error | UndefVarError `WrappedStandardMap` |
| W6 inverse round trip | 1 Error | UndefVarError `WrappedStandardMap` |
| W7 zero allocation | 1 Error | UndefVarError `WrappedStandardMap` |
| W8 overflow guard | 4 Error | UndefVarError `WrappedStandardMap` |
