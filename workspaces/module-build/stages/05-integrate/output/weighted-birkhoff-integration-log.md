# Integration log: weighted-birkhoff

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I13 (weighted Birkhoff averages, convergence rate) | planned -> verified; main-issue column states the pair index (log10 of the difference of the weighted averages on [0, N) and [N, 2N), floor 32 eps S), the R1 levels for cos q (regular median 11.2 digits at N = 1000, the Float64 floor 14.15 by N = 4000, chaotic 1.2 to 1.6; a level that moves with N, not a threshold), the Double64 level 24.9, the misreading of the wrapped and lifted angles, the device slower than the CPU at 4 observables and 11 horizons, and that weighted FLI, RRE and the parameterization method are not built; sources Das2016, SanderMeiss2020, DuignanMeiss2023, Das2017 |

I12c (power spectrum) and I12d (time-frequency) stay planned (M6, D80); I14 (REM, M5) stays planned; I17 carries the arclength candidate (46de9e6, unrelated to this module's evidence).

Evidence: `stages/04-verify/output/weighted-birkhoff-verification.md`, accepted at the stage 04 checkpoint (2026-10-28) with Amendment D85: criteria 1 to 13 pass (12(b) in the D85 quantile form, 12(c) nothing to test, 13(d) Double64 FALLBACK); folder suite 716 checks at default and spec sizes (105 of them in `wb_review.jl`); full package suite 1,398,168 pass, 0 fail, 1 broken (AC13 without CUDA in the process); `verify/weighted-birkhoff.jl` 78 PASS, 0 FAIL; GPU Float64 36 PASS, 0 FAIL, Double64 3 FALLBACK; independent review with no CRITICAL or HIGH finding, M1, M2, L1, L2, L3, L5 fixed, L4, L6, L8, L9 recorded.

Open: none for the module. Recorded limits: order 5 and above refused (the same limit in naff-fma E2 is unguarded: follow-up); at most 31 observables; Float64 digits plateau at 14 to 17; provenance keys collide with E2 in a chain and there is no `convention` key (L6, Experiments records the chain); the unnormalised sum overflows above floatmax / B (L9); the device kernel is slower than the CPU at many horizons. Not built: weighted FLI, RRE, the parameterization method, the multiscale profile (idea 25). R4 (sticky orbits) and criterion 16 go into the study spec. Follow-up for Experiments (obligation in `src/Experiments/CONTEXT.md`): the observable tuple, the horizons, the chain record. Follow-up modules: `rem` (I14, M5), `henon-map` (Y2), I12c and I12d (M6); then study p2 and the M4 exit (risks review, minor release). The full package suite was last run in full at stage 04 (1,398,168 checks); after it only the catalog, CONTEXT, changelog and spec notes changed.

## Commits

- 34d31f9 module-build weighted-birkhoff stage 01: spec approved (D83)
- 29dc943 module-build weighted-birkhoff stage 02: tests (30 testsets, 29 RED) and test plan; spec Amendment D84
- 6cc5554 module-build weighted-birkhoff stage 03: WeightedBirkhoff reducer, weighted_birkhoff and birkhoff_index in Indicators (GREEN, coverage 98.4%), bench, docs, impl notes
- 5ebe919 module-build weighted-birkhoff stage 02 correction: wb_case regular class beyond 50 starts (helper only, owner approved)
- 46de9e6 catalog: arclength indicator candidate on row I17 (unrelated to the module)
- 95041c2 module-build weighted-birkhoff stage 04: verification, review fixes, Amendment D85
- 5a886ef module-build weighted-birkhoff stage 05: catalog I13 verified; spec notes, Indicators limits, Experiments obligation, changelog

## Files changed by the module

Short hashes of the weighted-birkhoff commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (6cc5554, 95041c2, 5a886ef)
- `_config/build-order.md` (34d31f9)
- `_config/catalogs/indicators.md` (5a886ef; 46de9e6 for I17)
- `_config/decision-log.md` (34d31f9, 29dc943, 95041c2)
- `bench/weighted-birkhoff.jl` (6cc5554)
- `docs/src/api/indicators.md` (6cc5554)
- `src/Experiments/CONTEXT.md` (5a886ef)
- `src/Indicators/CONTEXT.md` (6cc5554, 95041c2, 5a886ef)
- `src/Indicators/Indicators.jl` (6cc5554)
- `src/Indicators/wbirkhoff.jl` (6cc5554, 95041c2)
- `src/Indicators/wbirkhoff_reducer.jl` (6cc5554, 95041c2)
- `src/precompile.jl` (6cc5554)
- `test/weighted-birkhoff/` `runtests.jl` (29dc943, 95041c2), `wb_helpers.jl` (29dc943, 5ebe919), `wb_review.jl` (95041c2), `wb_selftest.jl`, `wb_ac01_03.jl`, `wb_ac04_06.jl`, `wb_ac07_09.jl`, `wb_ac10_13.jl`, `wb_api.jl`, `wb_records.jl` (29dc943)
- `verify/gpu.jl`, `verify/weighted-birkhoff.jl`, `verify/weighted-birkhoff-poisson.py` (95041c2)
- `01-spec/output/weighted-birkhoff-spec.md` (34d31f9, 29dc943, 95041c2, 5a886ef)
- `02-tests/output/weighted-birkhoff-test-plan.md` (29dc943, 5ebe919)
- `03-implement/output/weighted-birkhoff-impl-notes.md` (6cc5554, 95041c2)
- `04-verify/output/` `weighted-birkhoff-verification.md`, `weighted-birkhoff-verify-log.txt`, `gpu-weighted-birkhoff-verify-log.txt`, `weighted-birkhoff-oracle.png`, `weighted-birkhoff-digits.png` (95041c2)
- `05-integrate/output/weighted-birkhoff-integration-log.md` (this commit)
