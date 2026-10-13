# Implementation notes: standard-map

Code: `src/Systems/` (`Systems.jl`, `standard_map.jl`, `linear_symplectic.jl`, `residuals.jl`, 241 lines in total, no function over 50 lines). Bench: `bench/standard-map.jl` (4 entries, not run by the tests).

## Result
`test/standard-map/runtests.jl` with the full Chaoslab loaded (Foundation `exact` present): 541 pass, 0 fail, 0 error, about 3 min with 2 threads.

## Deviations from the spec
- None in the maps. `wrap` maps a result that rounds up to exactly 2pi to 0 (the same torus point), so the range stays [0, 2pi).
- Non-finite guard: `sin` and `cos` are replaced by `x - x` (NaN) when `x` is not finite, so Inf gives NaN instead of a DomainError (SC1). The guard is one `isfinite` branch per call.
- `det_residual` uses the StaticArrays closed form up to 3 x 3 and LU on a dense copy above that (allocates for D > 3).
- `LinearSymplecticMap(A).Ainv` is `Omega^-1 A' Omega` written out in blocks (index moves and sign flips, no multiplications), so it adds no rounding.
- Systems imports nothing from Foundation (D23). `pi` in the element type comes from `oftype(x, pi)` (BigFloat: `BigFloat(pi, precision(x))`); Double64 gives the correctly rounded double-double pi.

## Performance (Float64 and Double64, K = 5, one orbit, single thread, 10^6 steps)
| Quantity | Float64 | Double64 |
|---|---|---|
| forward | 193 ns/step | 4.8 us/step |
| forward plus tangent | 328 ns/step | 9.6 us/step |
| allocations | 0 | 0 |
The orbit is chaotic and unbounded in the lifted coordinates, so `sin`/`cos` run on large arguments; steps near the origin are faster. No BenchmarkTools numbers (not a dependency).

## Known limits
- JET and Aqua not run here (not in the scratch environment); the orchestrator runs them.
- `symplectic_residual` and `det_residual` take static matrices (`symplectic_residual` needs `StaticMatrix`; `det_residual` any `AbstractMatrix`).
- `symplectic_form(T, d)` is type-stable only when `d` is a compile-time constant.
- Double64 `sin`/`cos` accuracy is the library's (about 2^-100), which is why the tests use eps 2^-100 (test plan, ambiguity 1).

## Review fixes (stage 04)

- `wrap`/`_wrap1` never return a negative value (`y < 0 && (y += tp)` before the `>= tp` guard); docstring states wrap is for display (absolute error about |x| eps_T, meaningless above |x| ~ 1/eps_T); BigFloat uses the global precision, so callers wrap in `setprecision`.
- `StandardMap` docstring: K must be built in the state's element type (use `exact`).
- `conjugated(A, B)` throws `ArgumentError` unless `B == B'`.
- `det_residual` (src/Systems/residuals.jl) uses `det(M)` on the SMatrix for every size; the dense `Matrix(M)` path and `_det` are removed.

## Amendment D28 (2026-10-13)

Code: `src/Systems/wrapped_standard_map.jl` (new, 81 lines; `WrappedStandardMap`, `torus_state`, `lift`, internal `_reduce` and `_checked`), `standard_map.jl` (the two tangent formulas moved into `_dt` and `_dtinv`, shared by both systems; bit-identical), `Systems.jl` (include and three exports). Bench: two wrapped entries added to `bench/standard-map.jl`.

### Result
`test/standard-map/runtests.jl` with the full Chaoslab: 774 pass, 0 fail, 0 error (541 before, 233 new in W1 to W8). Machine heavily loaded by other Julia processes, so the run took 9 min instead of 3 min.

### Deviations
- `_pi(x)` (internal, used by `wrap` and the wrapped system) is now `oftype(x, 3.141592653589793) + oftype(x, 1.2246467991473532e-16)` instead of `oftype(x, pi)`. Same value for Float64 (hi + lo rounds to hi) and Double64 (it is the correctly rounded double-double pi); the old form allocated 160 bytes per call for Double64 (BigFloat detour), which broke W7. BigFloat keeps its own method.
- Reduction order: the guard for a negative remainder runs before the guard for a remainder at or above 2 pi (spec W-E1). The reverse order can return exactly 2 pi for a tiny negative remainder.
- Exactness limit is `maxintfloat(typeof(x))` (2^53 Float64, 2^106 Double64, 2^precision BigFloat); reaching it (`>=`) gives NaN in all four components, also for NaN counts.
- W3 as written cannot hold for a free-running 256-bit orbit over 50 steps (chaos): the 256-bit counts are recomputed teacher-forced from the exact binary64 state, and the free run is compared for n_H(5) = 6 steps. Evidence in the test plan (first torus distance above 1e-3 at step 27, first count mismatch at step 34).

### Performance (K = 5, one orbit, one thread, best of 3, machine loaded, indicative only)
| Quantity | Float64 | Double64 |
|---|---|---|
| forward, lifted | 222 ns/step | 6.2 us/step |
| forward, wrapped | 208 ns/step | 1.9 us/step |
| propagate_orbit k = 2, lifted | 225 ns/step | 10.1 us/step |
| propagate_orbit k = 2, wrapped | 260 ns/step | 3.9 us/step |
Allocations: 0 for all of them. The lifted orbit is unbounded, so its sin and cos run on large arguments (slow range reduction in Double64); the wrapped state stays in [0, 2 pi), which is why Double64 wrapped is about 3 times faster. In Float64 the two floors cost about as much as the saved range reduction.

### Known limits
- Counts are exact only below `maxintfloat(T)` (2^53, about 9e15 in Float64). The count m is the running sum of n, so it grows faster than the number of steps: for a diffusive orbit at K = 5 (P wanders about 3.5 sqrt(k), n about 0.56 sqrt(k)) m is about 1e13 at 10^9 steps and reaches 2^53 near 10^10.5 steps; for an accelerator-mode orbit n grows linearly and m quadratically, so the limit is reached after about sqrt(2^53) = 1e8 steps. The guard turns that into a health error, not a wrong number; Double64 and BigFloat have no practical limit here. (Order-of-magnitude estimates, not measured.)
- `lift` in T loses absolute accuracy for large counts (|m| eps 2 pi); widen to BigFloat first.
- The wrapped system lives on a circle of circumference `2 * pi_T`, not 2 pi (spec W-E6): the shift per wrap is below pi eps_T.
