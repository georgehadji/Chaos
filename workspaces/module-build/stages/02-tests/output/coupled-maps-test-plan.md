# Test plan: coupled-maps

Spec: `../../01-spec/output/coupled-maps-spec.md` (approved, D103, open points P1 to P12 as proposed; Amendment D104 at the stage 02 checkpoint: R4 and R13). Test plan approved by the owner at the stage 02 checkpoint, 2026-11-05. Suite: `test/coupled-maps/` (`runtests.jl` includes `cm_helpers.jl`, `cm_ac01_04.jl`, `cm_ac05_06.jl`, `cm_ac07_12.jl`, `cm_ac13_16.jl`, `cm_ac17_19.jl`). Helpers are prefixed `cm_` (wrapped-count helpers `cmw_`). Seeds: `StableRNG(20261012 + offset)`, offset = 1000 x criterion + index (written in each testset; 17001, 18001, 18002, 19001 + N for the later ones). `eps_T`: Float64 2^-52; Double64 2^-100 (reading R1). All BigFloat work is in `setprecision(BigFloat, 256) do ... end`. No Float32 tests (D24). The outer `@testset "coupled-maps"` lets the file run on its own without a failing criterion aborting the rest. No regression value from the literature is used in the package suite: every exact answer is a closed form derived in the spec (items 6, 7, 9, 10), a K = 0 shear, or a 256-bit recomputation; the published numbers (criterion 20) go to stage 04.

## Criterion to testset

Sample sizes are the D38 reduction (package suite carries a representative case, the full stated size runs at stage 04 in `verify/coupled-maps.jl` with identical assertions and tolerances). "Y3 sets" are the five listed parameter sets plus one uniform draw in [-5, 5]^3 per point; "Y4 panel" is the fixed list of eight (N, K, eps, boundary) configurations in `CM_Y4_CONFIGS` (N in {2, 3, 5, 10}, K in {0, 0.5, 1, 3, 5}, eps in {0.2, 0.6102, 1, -0.5}, open and ring), points drawn per configuration.

| AC | Testset (file) | Sample sizes | Tolerances |
|----|----------------|--------------|-----------|
| 1 | "AC1 tangent vs derivative" (`cm_ac01_04.jl`) | Float64: 10 points per Y3 set and per Y4 config; Double64: 3 points per case, N <= 5 | every entry of M vs ForwardDiff, complex step h = 1e-30 (Float64), and the plain-formula tangent `cm_ref_tangent` (items 4, 9): 16 eps (1 + a_max) |
| 2 | "AC2 inverse tangent" | as AC1 | max abs(M(z) M_inv(T z) - I) <= 16 eps (1 + a_max)^2; M_inv vs ForwardDiff of `inverse` and vs the plain-formula inverse tangent (items 5, 9): 16 eps (1 + a_max) |
| 3 | "AC3 determinant and symplecticity" | as AC1, Y4 N <= 10 | 256-bit det and M' Omega M - Omega from the rounded matrix: abs(det - 1) <= 8 d eps normF(M) normF(M^-1) (inverse in 256 bits), entries <= 8 d eps normmax(M)^2; blocks A and G equal to the bit (lower-left block symmetric); module `det_residual` and `symplectic_residual` (also on a dense `Matrix`) within 64 d and 32 d times the same scales (reading R3) |
| 4 | "AC4 inverse round trip" | 10 orbits per Y3 set and Y4 config Float64, 4 Double64; n in {1, n_H}; both orders; plus a 256-bit check of the module orbit and inverse against the independent plain formulas (1 orbit per case) | plain lifted difference <= tau(n) = 32 eps M_n n sigma^n |
| 5 | "AC5 lifted vs wrapped"; "AC5 torus_state and lift" | 5 orbits per case Float64 (N <= 5), 2 Double64 (N <= 3), every k <= n_H; counts: 2 orbits Float64 and 1 Double64, 50 teacher-forced steps; torus_state/lift: 20 points in [-100, 100]^d for three systems plus one known Y3 point | torus distance of wrapped state vs lifted orbit <= tau(k); `lift` at 256 bits vs lifted orbit <= tau(k); coordinates in [0, 2 pi_T), counts integers; counts equal the 256-bit increments computed from the exact state of each step (`cmw_forced_ok`, rules of spec Algorithm item 3); `lift(torus_state(z))` vs z <= 16 eps(100) |
| 6 | "AC6 wrapped inverse"; "AC6 overflow guard" | 5 orbits per case Float64, 2 Double64, n in {1, n_H}, both orders; ten overflow states | coordinates <= tau(n), counts exactly equal; overflow: a count at maxintfloat - 1 that must increase gives NaN in all 2d components (position count, momentum count, inverse; Y3 and chain), at maxintfloat - 2 the step stays finite with the exact count |
| 7 | "AC7 K=0 exactness" | orbits: 20 Float64 and 4 Double64 per system (Y3, chain N = 2, 3 and 100 ring), 1000 steps, data multiples of 2^-10 below 8; products: one start for Y3, chain N = 3, N = 5 ring (Float64), Y3 and N = 3 (Double64), n in {1, 2, 10, 100, 1000, 10^4} | exact equality |
| 8 | "AC8 K=0 finite-time exponent" | Y3, chain N = 3, N = 5 ring (Float64), Y3 and N = 3 (Double64), n in {1, 10, 100, 10^3, 10^4, 10^6}; the K = 0 tangent is state independent (checked exactly on 100 orbit points), the n-step product is M0^n | product == [[I, nI], [0, I]]; abs(X(n) - ln(sqrt(n^2 + 1))/n) <= 8 eps (256-bit reference); position-only vector: X == 0 exactly |
| 9 | "AC9 decoupled limit bitwise" | (a) 3 orbits per (nu, kappa) of five sets Float64, 1 Double64, 1000 steps; (b) N in {2, 5, 10}, K in {1, 3}, open and ring, 2 orbits Float64, 1 Double64, 1000 steps; (c) Y3 and chain N in {3, 10} open and ring, 3 Float64, 1 Double64; plus a control that the coupling does propagate a perturbation | zero difference (`==`) at every step: Y3 (mu = 0) mapped by C(x, p) = (x + p, p) equals the Y1 module with K = -nu and K = -kappa; every chain site equals Y1 with the same K; perturbed pair or site leaves the others unchanged |
| 10 | "AC10 parameter symmetries" | 10 orbits per case Float64, 3 Double64 (N <= 5), every k <= n_H | torus distance <= tau(k) after shifting positions by pi built in T (Y3: (nu, kappa) -> -(nu, kappa), mu fixed; Y4: K -> -K) |
| 11 | "AC11 Y3 fixed point" | five sets plus two extra, 10^4 steps, Float64 and Double64; pairs (nu_i, kappa_{i+2}) from {-6, -2, 1, 3, 4.5, 6}; nu = 4 Jordan block powers n in {10, 100, 1000} | origin stays (zero difference); trace of the pair blocks 2 - nu, 2 - kappa to 8 eps (1 + abs(nu) + abs(kappa)); decoupled blocks (off-pair entries exactly 0); class; angle or ln lambda_+ to 32 eps (1 + abs(par)) times the condition factor max(1, 1/sqrt(abs(tr^2 - 4))); max entry of M^n / n within 1e-12 of 4 |
| 12 | "AC12 Y4 fixed point and exact spectrum" | origin stays 10^4 steps (N = 3 ring, N = 100, Double64 N = 5); N in {2, 3, 5, 10, 50, 100}, open and ring (ring N >= 3), K in {-2, 0.5, 1, 3}, eps in {0.2, 0.6102, -0.5} | eigenvalues of G(0) (Float64) vs 256-bit K + eps l_k to 16 eps (abs(K) + 4 abs(eps)); off-block entries of (U+U)' M(0) (U+U) <= 16 N eps (1 + a_max) (reading R4); trace per mode to 32 eps (1 + a_max); class and exponent outside the parabolic band (reading R4) |
| 13 | "AC13 pairing and sum at finite n" | Y3 five sets and Y4 N in {2, 3, 5} configs, 2 orbits each, Float64, n = 5 and n = n_H; 256-bit product of the rounded tangents, singular values by a cyclic Jacobi solver in 256 bits (`cm_svdvals`, no extra dependency) | abs(lambda_i + lambda_{d+1-i}) <= B(n)/n and abs(sum lambda_i) <= B(n)/n, B(n) = 8 n eps prod normF(M_k) sqrt(d) |
| 14 | "AC14 total momentum for K = 0" | N in {2, 3, 5, 10}, eps in {0.2, 0.6102, 1, -0.5}, open and ring, 5 orbits Float64, 2 Double64, every k <= n_H | abs(sum p(k) - sum p(0)) <= 8 N eps M_k k (256-bit sums) |
| 15 | "AC15 conversion of the printed form" | N = 2, beta in {0.1, 0.3051}, K in {1, 3}, 20 orbits each, Float64 | printed unit-torus lines (symplectic reading, test-side Float64) converted by q = 2 pi x_u, eps = 2 beta vs the chain: torus distance <= tau(k), k <= n_H |
| 16 | "AC16 element-type agreement" | 5 orbits per case Float64, 2 Double64 (N <= 5), plus N = 100 for K in {1, 5}, eps in {0.2, 1}: 2 orbits Float64 and 1 Double64 each | Float64 and Double64 vs 256-bit orbit at n = n_H: lifted and torus difference <= tau(n_H) with the eps of the lower precision; Double64 uses the measured effective eps (largest relative sin/cos error on the trig arguments of the sample, floor 2^-104), pass requires <= 2^-100 |
| 17 | "AC17 large N and bandwidth" (package part) | N in {10, 50, 100}, open and ring, (K, eps) in {(3, 0.6102), (5, 1)}, frame of 4 seeded vectors in [-1, 1] | M W and M^-1 W by `tangent_action!` and `inverse_tangent_action!` vs the plain-formula dense product and the module dense product, 16 eps (1 + a_max) per entry |
| 18 | "AC18 map-level reproducibility"; "AC18 CPU kernel vs plain loop" (package part) | Y3 at (0.5, 0.1, 1e-3) and (2, 0.5, 0.1): 10^5 orbits lifted, 10^4 wrapped, n = n_H; chain N in {2, 5}: 10^4 orbits; `propagate` (CPU backend) vs serial loop: 2000 orbits lifted and wrapped Y3, k = 4 | serial vs `Threads.@threads` loop `isequal`; `propagate` final states `isequal` to the plain loop |
| 19 | "AC19 zero allocation" | Y3 lifted and wrapped (pure functions) Float64 and Double64; chain N in {2, 10, 100}, Float64 and Double64: in-place forms lifted and wrapped, frame actions | `@allocated == 0` after warm-up inside a function; in-place results equal the pure forms bitwise |
| 20 | stage 04 | | published behaviour of orbits C and D (SALI, finite-time exponents), reported not gated (P10); needs the SALI and Lyapunov estimators |
| 21 | stage 04 | | ChaosTools hand-written discrete system in `verify/` (D27), 100 orbits per system |
| scope | "API dimension and constructors"; "SC1 non-finite input" | | dimensions (4, 2N, wrapped), supertypes, field names, constructor errors (reading R5); NaN, Inf, -Inf in positions and momenta give non-finite output with no exception (lifted forward, inverse, tangent; wrapped chain forward) |

## Deferred to stage 04 (no test in the package suite)

- AC17 (V), timing: tangent action on a frame of 4 vectors, time at N = 100 at most 2.6 times the time at N = 50, reported, not a CI check.
- AC18, GPU parts: GPU vs CPU on 10^5 Y3 orbits (and 10^4 Y4 orbits at N in {2, 5}) with the tau(n) bound, the Double64 GPU compile gate with reason "kernel does not compile for the element type", and the "1 and 8 CPU threads" comparison (the process thread count cannot be changed inside the suite; the package suite compares a serial loop with a `Threads.@threads` loop). `propagate` on the chain with `Vector` states depends on the P8 representation decision of stage 03 and is not tested here.
- AC20 (V, heavy, loose, K18): SALI of C and D, finite-time exponents of C and D at 10^6, from the Indicators and Lyapunov modules; reported to the owner, not a build failure.
- AC21 (V): ChaosTools oracle in the separate environment; O7, O9 (ForwardDiff of the n-step map vs the product of tangents at 10^4 points), O10 (256-bit reference orbits at the larger sample, measured double-double sin/cos accuracy), and the full D38 sample sizes of criteria 1 to 6, 9, 10, 14, 16.

## Assumed API (stage 03 implements exactly this)

Systems module, exported through `Chaoslab`. Order of positions and momenta: the first d/2 components are positions, the last d/2 momenta (Y3: `(x1, x2, x3, x4)`, spec item 1; chain: `(q_1..q_N, p_1..p_N)`).

| Name and signature | Contract |
|--------------------|----------|
| `FroeschleMap(nu, kappa, mu)` (`FroeschleMap{T} <: AbstractMap`, fields `nu`, `kappa`, `mu`) | Y3, state `SVector{4}`; parameters built in the state's element type like `StandardMap`; `dimension == 4`; `forward`, `inverse`, `tangent`, `inverse_tangent` return `SVector{4}` or `SMatrix{4, 4}` and accept any element type (Dual, Complex, Double64, BigFloat) |
| `StandardMapChain(N, K, eps; boundary = :open)` (`StandardMapChain{T} <: AbstractMap`, field `N`) | Y4; `boundary` in `:open`, `:ring`; `ArgumentError` for `N < 2`, `N > 100`, `:ring` with `N < 3` and unknown boundary symbols; `dimension == 2N` |
| `forward(sys, z)`, `inverse(sys, z)` for the chain | `z` any `AbstractVector{<:Number}` of length 2N (the tests pass `Vector{T}`; no StaticArrays for N > 4, P8); result a vector of the promoted element type; the protocol of `Systems.jl` otherwise unchanged |
| `tangent(sys, z)`, `inverse_tangent(sys, z_image)` for the chain | an `AbstractMatrix` convertible by `Matrix(...)`, 2N x 2N, `M = [[I + G, I], [G, I]]` (G at the pre-step positions) and `M^-1 = [[I, -I], [-G, I + G]]` (G at the preimage positions `q' - p'`) |
| `forward!(out, sys, z)`, `inverse!(out, sys, z)`, `tangent!(M, sys, z)`, `inverse_tangent!(M, sys, z)` for the chain, lifted and wrapped | in-place forms, return the filled argument, 0 bytes allocated after warm-up (AC19), results equal the pure forms bitwise; `out` and `z` do not alias |
| `tangent_action!(out, sys, z, W)`, `inverse_tangent_action!(out, sys, z_image, W)` for the chain | `out = M W` and `out = M^-1 W` for a `2N x k` frame `W` using the banded structure (spec item 9, cost linear in N), `out::Matrix{T}`, 0 bytes allocated (AC19) |
| `WrappedFroeschleMap(nu, kappa, mu)` | wrapped Y3, state `SVector{8}` = `(x1..x4, c1..c4)` (reduced coordinates in `[0, 2 pi_T)`, one exact integer count per coordinate), `dimension == 4` (tangent size, as `WrappedStandardMap`); `tangent` is the lifted one at the reduced state; a count at `maxintfloat(T)` gives NaN in all 8 components; usable by `propagate` (state longer than `dimension`) |
| `WrappedStandardMapChain(N, K, eps; boundary = :open)` | wrapped chain, state `AbstractVector` of length 4N, layout `(q_1..q_N, p_1..p_N, count q_1..q_N, count p_1..p_N)`, `dimension == 2N`; same overflow rule on 4N components; `tangent`, `inverse_tangent`, `tangent!`, `inverse_tangent!` read only the first 2N entries |
| `torus_state(sys, z)`, `lift(sys, s)` with `sys` a wrapped Y3 or chain system | two-argument forms (the length of the state does not identify the system: the wrapped chain at N = 2 and wrapped Y3 both have length 8); `torus_state` reduces each coordinate into `[0, 2 pi_T)` with its count, `lift` returns `z_j + 2 pi_T c_j` in the element type of `s` (length d) |
| `symplectic_residual(M::AbstractMatrix)` | extension of the existing `StaticMatrix` method to a dense `Matrix` (chain tangents are dense); `det_residual` already accepts any matrix |

Counts of the wrapped Y3 follow the derivation of the spec (Algorithm item 3): `m_j' = m_j + n_j + b_j` with the old momentum count `n_j`, `n_j' = n_j + a_j`; the wrapped chain follows W-E2/W-E4 of standard-map (inverse: `m = m' - n' + b`, `n = n' + a`, with the image count `n'`).

## Readings of the spec

- R1. Double64 eps. As in standard-map, DoubleFloats `sin`/`cos` do not reach the nominal 2^-104. All Double64 tolerances use 2^-100, the limit AC16 already accepts. AC16 still measures the effective value (floor 2^-104); on the sample the stub run measured 2^-104 (4.9e-32). Owner should confirm.
- R2. Sample sizes. The stated sizes (1000 points, 200 orbits) are reduced as in the table (D38 pattern), runtime of the whole suite on a stub implementation 20 s of evaluation after compile. The full sizes belong to `verify/coupled-maps.jl`.
- R3. Module residuals in AC3. The spec fixes the bound for a 256-bit evaluation of the rounded matrix, not for the module residual functions (evaluated in T). The package test applies the 256-bit bounds exactly and bounds the module functions at 64 d eps normF normF (determinant) and 32 d eps normmax^2 (symplectic), generous multiples that catch a wrong formula, not rounding. The symplectic residual needs a method for a dense `Matrix` (Assumed API).
- R4. AC12 off-block and parabolic modes. (a) The bound 16 eps (1 + a_max) for the off-block entries of (U+U)' M(0) (U+U) cannot hold at large N for a LAPACK eigensolver: measured on the real formulas the ratio is 0.02 (N = 2), 0.6 (N = 10), 20 (N = 50), 44 (N = 100). The test uses 16 N eps (1 + a_max) (ratio 0.44 at N = 100). (b) The specified grid contains exactly parabolic modes (g_k = 0 for K = 1, eps = -0.5 and N even; g_k = -4 for K = -2, eps = -0.5 on a ring with N even), where the class is decided by rounding. The class is asserted only where abs(abs(tr) - 2) > 64 eps (1 + a_max), the exponent only beyond 1e-6 and with the condition factor max(1, 1/sqrt(abs(tr^2 - 4))); the trace of every mode is always checked.
- R5. Constructor validation. The spec states the ranges (2 <= N <= 100, ring for N >= 3) but not the error type; the tests expect `ArgumentError` (as `symplectic_form` for an odd dimension). Field names `nu`, `kappa`, `mu`, `N` are part of the Assumed API.
- R6. AC13. B(n) is used literally with the 256-bit product of the Float64 tangents of the realised Float64 orbit; the measured ratio is about 1e-4, so the bound is conservative. Singular values come from a Jacobi eigenvalue solver in 256 bits applied to Y' Y (the package has no GenericLinearAlgebra).
- R7. AC5 counts. Teacher forcing recomputes, from the exact state of each step, the floor of the position sum (old momentum value) and of the momentum update with the force evaluated at the reduced new positions (Y3) or the pre-step reduced positions (chain), as the spec Algorithm item 3 states; a floor within 1e-15 of a tie would differ, which the seeds avoid.
- R8. AC9 bitwise. The comparison uses `==` (so -0.0 equals 0.0) and assumes no fused multiply-add contraction in either code path (Julia does not contract without `@fastmath`); the fallback in P9 is not implemented.
- R9. AC11 pairs. (nu, kappa) pairs are (v_i, v_{i+2}) over the list {-6, -2, 1, 3, 4.5, 6}, so every listed value is used as both nu and kappa; "max entry of M^n equal to 4n" is taken on the (1, 3) pair block of M(0), and tested as `abs(max/n - 4) <= 1e-12`.
- R10. AC15. The printed unit-torus map is evaluated test-side in Float64 in the symplectic reading; the literal-reading determinant remark of the spec is a finding, not scored, and has no test.
- R11. AC16 effective eps. The measured quantity is the relative error of Double64 `sin`/`cos` on the arguments the maps use (positions, neighbour differences, x1 + x2 for Y3), not a full orbit comparison against a double-double reference.
- R12. AC17/AC19 names. The spec asks for "the tangent action on a frame" and "a preallocated workspace" without names; the in-place and action functions of the Assumed API are my choice. The spec's frame action for the wrapped chain is not tested.
- R13. AC18 threads. 1 vs 8 CPU threads cannot be set inside one process, so the package test compares a serial loop with the threaded loop of the running process (at least 2 threads under `Pkg.test` on 2 threads); the 1 and 8 thread statement is stage 04.

## Inconsistencies between spec, docs and code

- Operation names: `_config/architecture.md` ("Module map") names `step`, `tangent`, `invariants` and StaticArrays for dimension <= 8; the code (`src/Systems/`) has `forward`, `inverse`, `tangent`, `inverse_tangent`, `dimension` and residual functions on `SVector` states. Followed the code (spec Finding (g)).
- `src/Systems/Systems.jl` docstrings say `forward` takes an `SVector`; a chain with N > 4 needs a general vector (P8). The tests accept any `AbstractVector`; the Assumed API only requires `SVector` for Y3.
- `symplectic_residual` exists only for `StaticMatrix` and `symplectic_form` builds an `SMatrix` (`Val` dimension); chain tangents up to 200 x 200 are dense. Needs a dense method (R3) and the spec says the residuals are "exposed for any state".
- `torus_state(z::SVector{2})` and `lift(s::SVector{4})` dispatch on the length of the state; the new wrapped states have length 8 for both wrapped Y3 and the wrapped chain at N = 2, so the new methods take the system as first argument (Assumed API); the Y1 one-argument methods stay.
- `src/Systems/CONTEXT.md` says Systems imports nothing from Foundation; the tests use `exact(T, pi)` from Foundation only through `using Chaoslab` (as the standard-map suite does).
- Spec criterion 12 bound for the rotated off-block entries is not attainable at N >= 50 (R4); criterion 11 reads "ratio to n within 1e-12" which is `abs(max/n - 4)`; criterion 18 "10^5 seeded orbits" and "1 and 8 CPU threads" cannot be run in the package suite (R13). Constants of criteria 1 to 3, 11, 12 are the P9 first-order proposals; stage 04 reports realized ratios.

## Runtime

Measured against a scratch stub of the Systems API (outside the repo, not committed; `scratchpad/cm02/stub.jl`, the stub allocates in the chain `inverse_tangent!` and frame actions, so 18 of the 40 AC19 assertions fail there by design): 556 assertions, 538 pass, the 18 AC19 ones fail; about 20 s of test evaluation after compilation on 2 threads, whole script 23 s. Compile cost of the real implementation adds to this. Worst observed ratios (error over bound, stub): AC1 0.05 (complex step 0.05), AC2 0.06, AC3 0.009, AC4 0.010, AC5 0.006, AC6 0.015, AC10 0.016, AC12 eigenvalues 0.28, off-block 0.44, mode trace 0.12, exponent 0.13, AC13 1e-4, AC14 0.021, AC15 0.008, AC16 2e-6, AC8 residuals below 5e-17 (bound 1.8e-15). The ratios in the first rows are partly zero because the stub uses the same formulas as the test references.

## RED log

Command (repo root): `julia --project=<scratchpad>/cm02/testenv -t2 <scratchpad>/cm02/red.jl`, where `red.jl` does `using Test, Chaoslab` and includes `test/coupled-maps/runtests.jl` inside a testset (testenv: Chaoslab by path, Test, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, ForwardDiff); wall time 8 s. Result: 0 Pass, 5 Fail, 69 Error, 0 Broken (74 records: `@test` records plus one record per testset aborted by an exception). Every error is `UndefVarError` on a missing name (`FroeschleMap` 34, `StandardMapChain` 24, `WrappedFroeschleMap` 6, `WrappedStandardMapChain` 5). The 5 Fail are the `@test_throws ArgumentError` constructor checks, which saw `UndefVarError` instead. No parse errors, no errors inside helpers, no warnings; no guard test passes (the suite has no assertion that exercises only existing code).

| Testset | Fail | Error | First error |
|---------|------|-------|-------------|
| API dimension and constructors | 5 | 19 | UndefVarError `FroeschleMap` |
| SC1 non-finite input | | 1 | UndefVarError `FroeschleMap` |
| AC1 tangent vs derivative | | 1 | UndefVarError `FroeschleMap` |
| AC2 inverse tangent | | 1 | UndefVarError `FroeschleMap` |
| AC3 determinant and symplecticity | | 1 | UndefVarError `FroeschleMap` |
| AC4 inverse round trip | | 3 | UndefVarError `FroeschleMap` |
| AC5 lifted vs wrapped | | 1 | UndefVarError `FroeschleMap` |
| AC5 torus_state and lift | | 1 | UndefVarError `WrappedFroeschleMap` |
| AC6 wrapped inverse | | 1 | UndefVarError `FroeschleMap` |
| AC6 overflow guard | | 1 | UndefVarError `WrappedFroeschleMap` |
| AC7 K=0 exactness | | 13 | UndefVarError `FroeschleMap` / `StandardMapChain` |
| AC8 K=0 finite-time exponent | | 1 | UndefVarError `FroeschleMap` |
| AC9 decoupled limit bitwise | | 7 | UndefVarError `FroeschleMap` / `StandardMapChain` |
| AC10 parameter symmetries | | 2 | UndefVarError `FroeschleMap` |
| AC11 Y3 fixed point | | 1 | UndefVarError `FroeschleMap` |
| AC12 Y4 fixed point and exact spectrum | | 1 | UndefVarError `StandardMapChain` |
| AC13 pairing and sum at finite n | | 1 | UndefVarError `FroeschleMap` |
| AC14 total momentum for K = 0 | | 2 | UndefVarError `StandardMapChain` |
| AC15 conversion of the printed form | | 1 | UndefVarError `StandardMapChain` |
| AC16 element-type agreement | | 2 | UndefVarError `StandardMapChain` |
| AC17 large N and bandwidth | | 1 | UndefVarError `StandardMapChain` |
| AC18 map-level reproducibility | | 1 | UndefVarError `FroeschleMap` |
| AC18 CPU kernel vs plain loop | | 1 | UndefVarError `FroeschleMap` |
| AC19 zero allocation | | 5 | UndefVarError `FroeschleMap` / `StandardMapChain` |

Because a missing name aborts most testsets at their first helper call, the bodies of the later assertions are validated by the stub run above, not by the RED run.
