# Systems

Maps with analytic one-step tangents. Imports nothing from Foundation (D23). Built by the standard-map run (catalog Y1 and Y5, `_config/catalogs/systems.md`).

Spec: `workspaces/module-build/stages/01-spec/output/standard-map-spec.md` (SM-; items are its Equations numbers, W- its Amendment D28). Reports: `stages/04-verify/output/standard-map-verification.md`.

## Files

| File | Holds |
|------|-------|
| `Systems.jl` | module, exports, the five generic functions of the map protocol |
| `standard_map.jl` | `StandardMap` (lifted), `wrap`, shared tangent formulas, non-finite guard |
| `wrapped_standard_map.jl` | `WrappedStandardMap` (torus state with winding counts), `torus_state`, `lift` |
| `linear_symplectic.jl` | `LinearSymplecticMap` and block builders for Y5 |
| `residuals.jl` | `symplectic_form`, `det_residual`, `symplectic_residual` |

## Exports

| Name | What it does |
|------|--------------|
| `AbstractMap`; `dimension`, `forward`, `inverse`, `tangent`, `inverse_tangent` | the protocol: dimension, one step, inverse step, Jacobian at the pre-step point, Jacobian of the inverse at the image |
| `StandardMap(K)` | `p' = p + K sin x`, `x' = x + p'`, state `SVector(x, p)`, lifted coordinates |
| `WrappedStandardMap(K)` | same map on the torus, state `(x, p, m, n)` with integer counts; dimension 2 |
| `torus_state`, `lift` | lifted point to wrapped state and back (exact counts, `T` arithmetic) |
| `wrap` | reduce to `[0, 2pi)` for display only |
| `LinearSymplecticMap(A)` | `z' = A z`, stored inverse `Omega^-1 A' Omega` |
| `rotation_block`, `hyperbolic_block`, `block_diagonal`, `conjugated` | build Y5 test maps, including the non-normal variant |
| `symplectic_form(T, d)`, `det_residual(M)`, `symplectic_residual(M)` | `Omega`, `abs(det M - 1)`, max entry of `M' Omega M - Omega` |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| `tangent` is evaluated at the pre-step point; analytic tangent matches ForwardDiff and complex step | SM Equations 3, SM-AC1 |
| Det and symplecticity residuals stay within the stated multiple of `eps(T)` | SM-I1, SM-I2, SM-AC2 |
| `inverse(forward(z))` returns `z` within the orbit budget; the inverse tangent matches | SM-I4, SM-AC3, SM-AC4 |
| Lifted orbit wrapped equals the torus orbit; `K` to `-K` symmetry holds | SM-I3, SM-I5, SM-AC5, SM-AC6 |
| `K = 0` and Y5 results are exact in exact arithmetic | SM-AC7, SM-AC8, SM-AC16 |
| Non-finite input gives non-finite output (no `DomainError` from `sin`) | SM Scope ("Domain"), W-I5 |
| Wrapped state: coordinates in `[0, 2pi_T)`, counts exact integers, lift agrees with the lifted orbit | W-I1 to W-I3, W1, W2 |
| At a count of `maxintfloat(T)` a step returns NaN in all four components, so the health check fires | W-I5, W8 |
| The orbit step allocates nothing after warm-up | W7 |
| `K` must be built in the state's element type (use `exact`) | `StandardMap` docstring |

## Tests and verification

- Package: `test/standard-map/` (`sm_ac01_06.jl`, `sm_ac07_13.jl`, `sm_ac14.jl`, `sm_ac18.jl`, `sm_y5.jl`, `sm_wrapped.jl`).
- Oracle run: `verify/standard-map.jl` (256-bit orbits, ForwardDiff, ChaosTools map convention); log and figure in `stages/04-verify/output/`.
- Notes: `stages/03-implement/output/standard-map-impl-notes.md`; test plan `stages/02-tests/output/standard-map-test-plan.md`.

## Deferred

| Item | Where |
|------|-------|
| Y1 stays `planned` until lyap-benettin verifies the large-K `lambda_1` known-answer test | catalog `systems.md`, D22 |
| Exponent comparisons with ChaosTools | lyap-benettin stage 04, D27 |
| Reset of winding counts at each run checkpoint (the rule is set, no checkpoint writer exists yet) | D31 |
| Counts beyond `maxintfloat(T)` (about 1e8 steps on an accelerator-mode orbit in Float64) | health error today, D28, D31 |
| Float32 | D24 |
| Flows and integrators | M5 |
