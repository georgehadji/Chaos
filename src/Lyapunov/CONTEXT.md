# Lyapunov

Exponent estimators as reducers of the tangent stream, plus the two-trajectory estimator that runs its own loop. Imports Foundation, Systems and Tangent. Built by the lyap-benettin (catalog L1), lyap-gram-schmidt (L3 CGS, L4 MGS, L4b RGS), lyap-householder-qr (L5 Householder, L6 Givens, L5b augmented), lyap-svd-stabilized (L9) and lyap-two-trajectory (L2) runs.

Spec: `workspaces/module-build/stages/01-spec/output/<slug>-spec.md` for the five slugs above (amendments D35 for all five bounds, D41 svd-stabilized, D42 two-trajectory; owner answers D30; API contract D34). Reports: `stages/04-verify/output/<slug>-verification.md`, same five slugs.

## Files

| File | Holds |
|------|-------|
| `Lyapunov.jl` | module, exports, include order |
| `norm.jl` | `safe_norm`, the one norm every estimator uses |
| `benettin.jl` | `Benettin` reducer and its isbits state (`k = 1`) |
| `orthonormalize.jl` | `Orthonormalizer` types `CGS`, `MGS`, `RGS`, status codes, `orthonormalize` for them, `window_Gw`, `window_G`, `orthogonality_threshold` |
| `qr_spectrum.jl` | `QRSpectrum` reducer, its isbits state, the shared handoff step `_qr_step` |
| `householder.jl` | kernel-path `HouseholderQR`, `AugmentedQR`, `GivensQR`; the sign fix `_finish`; provenance of the whole Householder family |
| `lapack_qr.jl` | `LapackQR`, the LAPACK library path (non-isbits state with a preallocated workspace) |
| `svd_spectrum.jl` | `SVDSpectrum` reducer, `log_singular_values` |
| `separation.jl` | `separation` and `displaced` hooks, with `WrappedStandardMap` methods (minimal image) |
| `two_trajectory.jl` | `propagate_shadow`, `propagate_shadow_orbit`, the per-orbit body, kernel, chunking and device fallback |

## Exports

| Name | What it does |
|------|--------------|
| `safe_norm(v)` | scaled Euclidean norm; no overflow or underflow, zero for a zero vector, NaN and Inf propagate |
| `Benettin()` | reducer for `k = 1`: adds `log` of the deviation norm to a compensated sum, renormalises; final `(n, sum, comp, lambda)` |
| `QRSpectrum(alg; orthogonality_tol)` | reducer for `k >= 1`: factors `W = Q R` at each handoff, sums `log R_ii`, replaces `W` by `Q`; owns the orthogonality and conditioning checks |
| `CGS`, `MGS`, `RGS` | Gram-Schmidt `Orthonormalizer` strategies (classical, modified, repeated with a pass cap) |
| `HouseholderQR`, `AugmentedQR`, `GivensQR` | kernel-path strategies (isbits, run on the GPU), provenance V-H, V-A, V-G; diagonal iteration path only |
| `LapackQR` | library strategy (`QRSpectrum(LapackQR())`): LAPACK `geqrf`/`orgqr`, Float64 only, CPU only |
| `Orthonormalizer`, `orthonormalize(alg, W)` | strategy supertype; returns `(Q, R, count, status, column)` with `R_ii > 0` |
| `QR_OK`, `QR_PIVOT`, `QR_PASSES` | status codes of `orthonormalize` |
| `window_Gw`, `window_G`, `orthogonality_threshold` | window conditioning of `R` and the orthogonality bound built from it |
| `SVDSpectrum(; orthogonality_tol)` | reducer for `d = k = 2`: the RGS spectrum plus singular-value exponents `lambda` beside `lambda_qr` |
| `log_singular_values(S, c)` | log singular values of the triangular product from the QR totals and the ratio `c` |
| `propagate_shadow(sys, x0, N; d0, tau, seed, backend, chunk, U0)` | two-trajectory estimate X2 for every orbit; returns a `PropagationResult` |
| `propagate_shadow_orbit` | the same body for one orbit, serial; throws `NumericalHealthError` |
| `separation(sys, ref, shadow)`, `displaced(sys, ref, delta)` | separation vector and shadow placement; plain difference by default, minimal image for `WrappedStandardMap` |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| `orthonormalize` returns `Q'Q = I`, `R` upper triangular with exact zeros below and `R_ii > 0`; a zero or non-finite pivot gives `QR_PIVOT` with the pivot in `R[column, column]` | lyap-gram-schmidt spec G-items, lyap-householder-qr spec E5 |
| Every estimator takes norms through `safe_norm` | lyap-benettin spec Algorithm 1, AC8 |
| Sums are compensated; `reducer_final` gives `lambda = (sum + comp) / n` | lyap-benettin spec E3 |
| A zero or non-finite norm, pivot, state or separation is recorded in the state and raised as `NumericalHealthError` with orbit, step and backend, never thrown inside a kernel | D30, `Tangent` failure path |
| `SVDSpectrum` reuses the RGS step unchanged: its `lambda_qr` equals that of `QRSpectrum(RGS())` bitwise | lyap-svd-stabilized spec AC3a |
| `SVDSpectrum` refuses any `d` or `k` other than 2 at `reducer_init`; general `d` is held | D32, `SVD_GUARD_MESSAGE` |
| CPU outputs are bitwise identical for 1 and 4 threads | `verify/threads.jl` |
| `propagate_shadow` on the CPU equals a serial `map` of `propagate_shadow_orbit` bitwise, for any chunk size and thread count | lyap-two-trajectory spec, stage 04 report |
| `d0` must widen exactly into the element type; a realized start separation outside `[d0/2, 2 d0]` raises before any chunk | lyap-two-trajectory spec Algorithm 1 |
| Kernel-path estimators (Benettin, MGS, RGS, the three Householder variants, SVDSpectrum, `propagate_shadow`) run on the CUDA GPU in Float64; CGS shares the kernel but was not run on the GPU | `verify/gpu.jl` |
| Double64 on the GPU falls back to the CPU with reason "kernel does not compile for the element type" | `verify/gpu.jl` |
| `QRSpectrum(LapackQR())` falls back to the CPU with reason "reducer is CPU-only" and needs element type Float64 | `lapack_qr.jl`, stage 04 report |
| `SVDSpectrum` needs the 4 KiB device stack that the CUDA extension sets at load | D40 |
| Provenance names the estimator, orthonormalizer, iteration path (`diagonal`), implementation (`kernel` or `library`) and tolerance | lyap-householder-qr spec AC15e |

## Tests and verification

- Package: `test/lyap-benettin/`, `test/lyap-gram-schmidt/`, `test/lyap-householder-qr/`, `test/lyap-svd-stabilized/`, `test/lyap-two-trajectory/`.
- Oracle runs: `verify/<slug>.jl` for each of the five slugs, `verify/threads.jl` (thread identity), `verify/gpu.jl` (GPU lines).
- Bounds, tolerances and measured ratios live in the specs and verification reports; this file points to them.

## Deferred

| Item | Where |
|------|-------|
| H1 (Householder, second iteration path) | lyap-householder-qr spec, held; only H2 (diagonal iteration) exists |
| General-`d` singular-value exponents | held for M3, D32 |
| lyap-two-trajectory criterion 13 (external two-trajectory oracle), sent back to stage 01 | D43 |
| Float32 | D24 |
