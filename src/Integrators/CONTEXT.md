# Integrators

Symplectic and reference integrators for Hamiltonian flows (catalog `_config/catalogs/integrators.md`, rows S1 to S9, S11 to S13, P1 to P3, G1 to G3, R1, tangent dynamics by the same splitting). Imports Foundation and Systems only (layering test). Built by the integrator-catalog run (M5).

Spec: `workspaces/module-build/stages/01-spec/output/integrator-catalog-spec.md` (D111); section and escape tools: `henon-heiles-spec.md` (D116), tests `test/henon-heiles/`. Tests `test/integrator-catalog/`, plan `stages/02-tests/output/integrator-catalog-test-plan.md`, notes `stages/03-implement/output/integrator-catalog-impl-notes.md`.

## Files

| File | Holds |
|------|-------|
| `Integrators.jl` | module, exports, `taylor_reference` (declared; methods from `ext/ChaoslabTaylorExt.jl`) |
| `splitting.jl` | `SplittingMethod` (table of drift weights `c` and kick weights `d`, isbits, name and source as type parameters), validation, metadata, `adjoint_method`, `triple_jump`, `suzuki_jump` |
| `tables.jl` | S1 to S9, S11 to S13: `SymplecticEuler`, `Leapfrog`, `McLachlan2`, `Ruth3`, `ForestRuth`, `PEFRL`, `CalvoSanzSerna4`, `Yoshida`, `KahanLi`, `SofroniouSpaletta10` (constants copied from the cited source as strings, each with source key and equation or table number) |
| `tables_near.jl` | P rows: `SABA`, `SBAB` (Gauss-Legendre and Gauss-Lobatto quadrature), `ABA864`, `ABAH864`, `CorrectedMethod`, `SABAC`, `SBABC` |
| `stepper.jl` | `advance`, `advance_n` (merged kicks, optional compensated update), `advance_to`, `advance_checked`, `tangent_advance` for the table and corrected methods |
| `gauss.jl` | `Gauss`, `ImplicitMidpoint`, `gauss_coefficients`: fixed-point solver of Hairer2008 Sect. 3 with compensated stage and update sums (D113), status, tangent by implicit differentiation |
| `tao.jl` | `TaoExtended`, `extend_state`, `physical_state`: Tao2016 extended phase space, Strang and triple jump |
| `integrator_map.jl` | `IntegratorMap`: the one-step map as a Systems `AbstractMap` |
| `ensemble.jl` | `advance_ensemble`: one KernelAbstractions kernel over orbits |
| `section.jl` | `locate_crossing` (sub-step root solver of the method's own step), `section_map` (upward crossings of `q1 = 0`, optional Jacobian of the leg in `(q2, p2)`) (henon-heiles, D116) |
| `escape.jl` | `escape_sector`, `locate_escape`, `escape_orbit` (first passage of `r = R_esc`, exit label) |

## Rules a caller can rely on

| Rule | Where |
|------|-------|
| States are `SVector{2n}` (q, then p; Tao: q, p, x, y); a method is built in the element type of the state (`Leapfrog(Double64)`) | all |
| `advance(m, flow, z, tau)` accepts either sign of `tau`; the inverse of a step is the adjoint table at `-tau` (the table itself when palindromic) | `adjoint_method`, `IntegratorMap` |
| Table methods allocate nothing after warm-up, also under `--check-bounds=yes`: immutable `SVector` state, tuple tables, no `MArray` (decision D110) | `stepper.jl` |
| Non-finite input gives non-finite output; an implicit step that fails returns an all-NaN state and `ok = false`, nothing throws | `gauss.jl`, `ensemble.jl` |
| A table with a violated sum condition, a zero stage count or a false palindromic claim is rejected at construction (`ArgumentError`) | `SplittingMethod` |
| The method's coefficient type, the state's element type and a floating step size must agree, else `ArgumentError` (D113); `force_evaluations` is the cost per step of a long run, a corrected method adds its 2 corrector gradients (D113) | `stepper.jl`, `splitting.jl` |
| No renormalization of deviation vectors here (Tangent, Lyapunov own that) | `tangent_advance` |

## Limits

Fixed step only. Tao's tangent frame lives in the extended space (open point P11). The Gauss tangent builds a dense `s * 2n` system (heap, not for kernels). `advance_ensemble` takes a KernelAbstractions backend (`ensemble_backend` hook, defined for the CUDA device in `ext/ChaoslabCUDAExt.jl`, D113); the GPU device runs (`verify/gpu.jl`) were not run at stage 04. Static arrays only: the full tangent frame is measured up to n = 50 degrees of freedom, a `Vector` path is open (risk K22).
