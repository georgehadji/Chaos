# Implementation notes: lyap-gram-schmidt

Code: `src/Lyapunov/orthonormalize.jl` (147 lines: `Orthonormalizer`, `CGS`, `MGS`, `RGS`, `QR_*`, `orthonormalize`, `window_Gw`, `window_G`, `orthogonality_threshold`) and `src/Lyapunov/qr_spectrum.jl` (92 lines: `QRSpectrum`, `QRState`, `orth_loss`, the reducer methods and `reducer_provenance`). `Lyapunov.jl` includes both and exports the new names. Bench: `bench/lyap-gram-schmidt.jl` (6 entries; not run, BenchmarkTools is not a dependency). No function over 50 lines. `effective_eps`, the `column` field of `NumericalHealthError` and `reducer_provenance` already existed (D34 amendments in Foundation and Tangent), so nothing outside `src/Lyapunov/` changed.

## Status

`test/lyap-gram-schmidt/` on Julia 1.13.1, 2 threads: 19553 pass, 0 fail, 2 error. Both errors are the same test line (`gs_api.jl:66`, API2): `HouseholderQR` is not defined. That line belongs to the Householder slug's type; it is left undefined on purpose (see "Test dependency"). Everything else in the suite is green, including the zero-allocation checks and AC10 bitwise equality with Benettin.

Full package suite (`Pkg.test`, 2 threads), Julia 1.13.1: quality (Aqua and JET) 12/12, layering 2/2, diagnostics 245, precision 150, standard-map 774, tangent-core 565, lyap-benettin 15009, lyap-gram-schmidt as above. The remaining non-pass results are the RED suites of lyap-householder-qr, lyap-svd-stabilized and lyap-two-trajectory (101 fail, 42 error in total including the 2 above); the errors are `UndefVarError` (HouseholderQR, LapackQR, SVDSpectrum, log_singular_values, propagate_shadow, propagate_shadow_orbit), the 101 failures are all inside lyap-two-trajectory (its `separation` and wrapped-state checks, which cannot pass without that module).

Julia 1.10.12: `test/lyap-gram-schmidt/` alone in a fresh environment (Chaoslab developed from this tree): 19553 pass, 2 error (the same `HouseholderQR` line), 2 min 13 s testset time. The 1.10 `Pkg.test` of the whole package was started but stopped before its summary (standard-map alone took 100 min on 1.10 in the Benettin run and is untouched by this change); its log showed no quality or layering failure before the lyap-gram-schmidt block, and Aqua and JET are clean on 1.13 (12/12). A full 1.10 run is left to the orchestrator.

## Runtime against the 3 minute ceiling (D35)

lyap-gram-schmidt on `julia -t 2`, test code compiled in the run, package precompiled: 2 min 31 s testset time (2 min 45 s wall including start-up and package load) in a single-suite run; 3 min 15 s inside the full-suite run, which shared the machine with several other Julia processes. So the folder is under the ceiling when run alone and at or over it on a busy machine. No cut from the test plan was applied. If the owner wants margin, the plan's cut order is: API1 to API3 units, Double64 parts of AC1 and AC4, chunk 1000 for RGS in AC9, the 10^6-step run of AC8. Largest testsets (single-suite run): AC1 19 s, AC10 11 s, AC8 15 s, AC9 17 s, AC3 14 s, API1 and API2 about 14 s each.

## Deviations from the contract

- `orthonormalize` returns `R[column, column]` = the offending pivot on a failed status. The contract says Q and R are unspecified then; the step reads the pivot from there to fill `value` of the health error (the contract gives no other channel). A Householder `orthonormalize` should do the same or the generic failure record needs a value argument.
- The state stores `value` and `threshold` as `T` (as the contract), converted back to `Float64` in `reducer_error`.
- On any failure `G_max` and `orth_max` are also left unchanged (the contract only names sums and W); on success they are updated every handoff, `orth_max` stays 0 for `k = 1`.
- `window_Gw` is computed once per handoff and `G` derived from it (`_G_from_Gw`); `window_G(R)` is the public wrapper.
- Hot-path scratch is immutable `SMatrix`/`SVector` built with `setindex` and a tuple-based column setter. A first version with `MMatrix` scratch allocated 288 bytes per call (and 432 bytes per step in `propagate_orbit`) and failed the zero-allocation tests; it was replaced.

## Design for the Householder slug (additive)

`QRSpectrum{A}` and `QRState{T,K}` do not depend on the algorithm: the step calls `orthonormalize(r.alg, W)`, `checks_conditioning(alg)` (false by default, true for CGS, MGS, RGS) and `orthogonality_threshold` (`g = 1` by default through `_gfactor`). A new variant adds a subtype of `Orthonormalizer`, an `orthonormalize` method, `orthonormalizer_name` and (for the extra keys) a more specific `reducer_provenance` method. `LapackQR` needs its own `reducer_init` returning a non-isbits state. `QR_PASSES` is only produced by RGS.

## Performance (Julia 1.13.1, 1 thread, top level, machine heavily loaded by other agents, indicative only)

- `orthonormalize` 4 x 4 Float64: CGS 1.9 us, MGS 1.4 us, RGS 1.6 us; 0 bytes.
- `reducer_step` RGS d = k = 2: 1.1 us; `propagate_orbit` on the standard map K = 2, k = 2: 0.9 us per step (CGS and RGS alike), 0 bytes in `propagate_orbit` (AC8, AC9).
- These are well above a tight kernel (about 3 to 5 times Benettin per step); the cost is the `setindex` and column-setter copies. A hand-unrolled `k = 2` path would be the first optimisation if a sweep needs it.

## Known limits

- `window_G` is exact for `k <= 2`; for `k >= 3` it returns `G_w`, an upper bound (API Q5), so `sum_threshold` is only looser.
- `window_Gw` uses plain squared sums for the Frobenius norms (no overflow guard); entries beyond about 1e150 would overflow.
- Test dependency: `gs_api.jl:66` calls `orthogonality_threshold(HouseholderQR(), ...)`. The Gram-Schmidt suite can only be fully green once the Householder slug defines `HouseholderQR`. Declaring the empty type here would make the GS suite green but would change the failure mode of the Householder suite from `UndefVarError` to `MethodError`; the owner should decide (the type is one line, `struct HouseholderQR <: Orthonormalizer end`).
- GPU: the step uses `dot`, `sqrt`, `log`, division and `setindex` on isbits values; compilation for the device is untested here (stage 04, per-chunk fallback with reason 3).
- Coverage not measured in this run.
