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
