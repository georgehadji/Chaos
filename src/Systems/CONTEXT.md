# Systems

Maps with analytic one-step tangents. Imports nothing from Foundation (D23). Built by the standard-map run (catalog Y1 and Y5) and the coupled-maps run (Y3, Y4), `_config/catalogs/systems.md`.

Spec: `workspaces/module-build/stages/01-spec/output/standard-map-spec.md` (SM-; items are its Equations numbers, W- its Amendment D28). Coupled maps: `coupled-maps-spec.md` (items; Amendment D104), tests `test/coupled-maps/`. Reports: `stages/04-verify/output/standard-map-verification.md`. `PairMap`: `workspaces/module-build/stages/01-spec/output/rli-spec.md` (RL5, P1(c), Amendment D76); tests `test/rli/rl_pair.jl`.

## Files

| File | Holds |
|------|-------|
| `Systems.jl` | module, exports, the five generic functions of the map protocol |
| `standard_map.jl` | `StandardMap` (lifted), `wrap`, shared tangent formulas, non-finite guard |
| `wrapped_standard_map.jl` | `WrappedStandardMap` (torus state with winding counts), `torus_state`, `lift` |
| `linear_symplectic.jl` | `LinearSymplecticMap` and block builders for Y5 |
| `pair.jl` | `PairMap` (product of a map with itself) and its block-diagonal `PairTangent` |
| `froeschle_map.jl` | `FroeschleMap` (Y3, lifted `SVector{4}`) and `WrappedFroeschleMap` (`SVector{8}`), `torus_state`, `lift` |
| `standard_map_chain.jl` | `StandardMapChain` (Y4, `Vector` state), in-place forms, banded frame actions |
| `chain_tangent.jl` | `ChainTangent` (banded tangent, isbits) and the `SVector`-state `forward`, `inverse`, `tangent`, `inverse_tangent` of the chain (D105) |
| `wrapped_standard_map_chain.jl` | `WrappedStandardMapChain` (state length `4N`), `torus_state`, `lift`, `SVector`-state forms (D105) |
| `residuals.jl` | `symplectic_form`, `det_residual`, `symplectic_residual` |
| `flow_core.jl` | `AbstractFlow`, `AbstractSeparableFlow`; the exact partial flows `drift`, `kick` (and increments, Jacobians), `hamiltonian_gradient`, `vector_field`, `CountingFlow` (integrator-catalog, D111) |
| `flow_models.jl` | `HarmonicOscillator` (F1), `Pendulum` (F4), `Kepler` (F2), `HenonHeilesFlow` (T2), `TaoNonseparable` (T1) with Hamiltonian, force, Hessians and the SABAC corrector pair |

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
| `PairMap(sys)`, `PairTangent` | state (A; B) of two states of `sys`, dimension `2 d`; `forward`, `inverse` act on each half, `tangent` is `diag(M(A), M(B))` whose product with a frame acts block by block (rli, D76) |
| `FroeschleMap(nu, kappa, mu)`, `WrappedFroeschleMap` | Y3 in the form of Skokos2010 Eq. (16): state `(x1, x2, p1, p2)`, force at the new positions; wrapped state adds one count per coordinate; dimension 4 |
| `StandardMapChain(N, K, eps; boundary)`, `WrappedStandardMapChain` | Y4, `2 <= N <= 100`, `:open` or `:ring` (`N >= 3`), state `(q, p)` of length `2N`; wrapped `(q, p, count q, count p)`; dimension `2N`; `ArgumentError` outside the ranges |
| `ChainTangent` | banded `2N x 2N` tangent at an `SVector` chain state (length `2N`, wrapped `4N`): entries equal the dense tangent, `M * W` for an `SMatrix` frame is the banded action; chain values are isbits, steps and products allocate 0 bytes (D105) |
| `forward!`, `inverse!`, `tangent!`, `inverse_tangent!`, `tangent_action!`, `inverse_tangent_action!` | chain only: in-place forms (0 bytes after warm-up), dense `2N x 2N` tangent, `M W` and `M^-1 W` for a frame in cost linear in `N` |
| `torus_state(sys, z)`, `lift(sys, s)` | two-argument forms for the wrapped Y3 and chain (state length does not identify the system) |
| `HarmonicOscillator(w)`, `Pendulum(eps)`, `Kepler()`, `HenonHeilesFlow()`, `TaoNonseparable()` | Hamiltonian flows, states `(q, p)`: `hamiltonian`, `force`, `force_jacobian`, `drift`, `kick` and their Jacobians, `vector_field`, `exact_flow` (F1), `angular_momentum`, `runge_lenz` (F2), `perturbation_scale`, `corrector_force`, `corrector_hessian`; `CountingFlow` counts force evaluations |
| `symplectic_form(T, d)`, `det_residual(M)`, `symplectic_residual(M)` | `Omega`, `abs(det M - 1)`, max entry of `M' Omega M - Omega` (also for a dense `Matrix`) |

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
| Y3 at `mu = 0` and the chain at `eps = 0` equal Y1 bitwise (same operation order, no `@fastmath`); a coupled step is plain sines, no FMA | coupled-maps criterion 9 |
| Chain count at `maxintfloat(T)` gives NaN in all `4N` components; Y3 in all 8 | coupled-maps criterion 6 |
| `K` must be built in the state's element type (use `exact`) | `StandardMap` docstring |

## Tests and verification

- Package: `test/coupled-maps/` (`cm_ac01_04.jl` to `cm_ac17_19.jl`) and `test/standard-map/` (`sm_ac01_06.jl`, `sm_ac07_13.jl`, `sm_ac14.jl`, `sm_ac18.jl`, `sm_y5.jl`, `sm_wrapped.jl`).
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
| Chain on the GPU: the CPU stream works with `SVector` states (D105); the GPU comparison is criterion D105d | stage 04 of the rework |
| Flows and integrators | M5 |
