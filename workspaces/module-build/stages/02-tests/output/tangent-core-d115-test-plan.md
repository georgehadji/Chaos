# tangent-core, Amendment D115 (frame product with explicit fma): stage 02 test plan

Spec: `01-spec/output/tangent-core-spec.md`, "Amendment D115" (items 1 to 4, criteria D115a to D115h). Interface: `Chaoslab.Foundation.fused_product(A, B)`, not exported.

Files: `test/tangent-core/tc_d115.jl` (new), registered in `test/tangent-core/runtests.jl`; `test/runtests.jl` logs the CPU of the job (D115g). Seeds `TC_SEED + 11500 + sub` (sub 1 to 4), every draw a `StableRNG(seed)`.

## Criterion to test

| Criterion | Test (file:testset) |
|---|---|
| D115a definition (Float64, Float32; m = p in {2, 4, 6, 8}, n in {1, 2, 3, 4, 8}; 20 draws each, B scaled by 10^(-3..3) per entry; equals the scalar fma reference bitwise; type `SMatrix{m, n, T}`) | `tc_d115.jl:D115a definition` |
| D115b column independence (leading c columns and single columns equal the product of the sliced frame, bitwise) | `tc_d115.jl:D115b column independence` |
| D115c other element types (Double64, BigFloat: `isequal` to `A * B`) | `tc_d115.jl:D115c other element types` |
| D115d stream (standard map K = 6, k = 1, 2; Froeschle map (0.5, 0.1, 0.001), k = 2, 3, 4; `PairMap` of the standard map, k = 1, 2; N = 500, tau = 1, no reducer: `propagate_orbit` equals the fused loop bitwise, leading columns of the widest run equal the narrower runs) | `tc_d115.jl:D115d stream` |
| D115e no allocation (4 x 4 Float64; skipped under forced bounds checks, D110) | `tc_d115.jl:D115e no allocation` |
| D115f no golden moves | whole package suite, stage 04 |
| D115g CI matrix | CI after the push: gali AC8a, AC8b, AC17b on all four Julia jobs; the "test machine" line of `test/runtests.jl` |
| D115h GPU | not run (`verify/gpu.jl` writes under `data/`) |

## RED expectation

D115a to D115e error with `UndefVarError: fused_product` (D115d errors in its first reference loop). The rest of the tangent-core suite is unchanged.
