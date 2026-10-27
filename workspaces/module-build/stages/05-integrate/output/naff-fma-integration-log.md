# Integration log: naff-fma

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I12 (frequency map analysis, NAFF, index D) | planned -> verified; main-issue column states the index definition and floor (32 eps), the 9-decade separation of the regular and chaotic medians at N = 1000 in 2D (a level that moves with N, not a threshold) and that the E1 series feed is Experiments' (criterion 16, study spec); Laskar1990, LaskarFroeschleCelletti1992, Laskar1993 marked not held, Laskar2003 and Bazzani2023 added |
| `_config/catalogs/indicators.md` | I12b (refined Fourier / NAFF variants) | planned -> verified; main-issue column states the Hanning default with error of order T^-4, window orders 0 to 5, refinement in the element type, E2 as an O(1) reducer, and that the closed-form interpolation, multi-term decomposition and collocation are not built; sources Gomez2010a, Bartolini1996, Laskar2003, LuqueVillanueva2014, Das2017, Gomez2010b (not held) |

I12c (power spectrum) and I12d (time-frequency) stay planned (deferred to M6, D80). I13 (weighted Birkhoff, own module `weighted-birkhoff`) and I14 (REM, M5) stay planned.

Evidence: `stages/04-verify/output/naff-fma-verification.md`, accepted at the stage 04 checkpoint (2026-10-27) with Amendment D82: criteria 1 to 14 and 15 (in the D82 form) pass, criterion 16 is in the study spec; folder suite 905 checks at spec sizes (895 at package sizes), package suite 1,397,452 pass, 0 fail; `verify/naff-fma.jl` 14 PASS; GPU Float64 12 PASS, 0 FAIL in the D82 form, Double64 2 FALLBACK; independent review with no CRITICAL or HIGH finding, M1 to M3 fixed, M4 recorded.

Open: none for the module. Recorded limits: one-ulp-noise constant series are `:ok` under mean removal; a broad-band window with more than 64 near-top local maxima could miss the true maximum; provenance keys collide in chains (M4); E2 on the planar observable is branch dependent near |Delta| = pi (`near_pi` counts them). Not run: Python NAFF packages (14(c), not installed), the amplitude-error slope, the regular-orbit error law against a true frequency, the rest of R2, R8, R9; R6 used the standard map's elliptic point (pi, 0) because Y2 (Henon) is not available. Sources not held: Laskar1990, LaskarFroeschleCelletti1992, Laskar1993, Laskar1999, Gomez2010b. Follow-up for Experiments (obligation in `src/Experiments/CONTEXT.md`): the series feed for E1 and the chain record; criterion 16 goes into the study spec. Follow-up modules: `weighted-birkhoff` (I13), `rem` (I14, M5), `henon-map` (Y2), I12c and I12d (M6). The full package suite was last run in full at stage 04 (1,397,452 checks); after it only the catalog, CONTEXT, changelog and spec notes changed.

## Commits

- db3992c module-build naff-fma stage 01: spec approved (D80): series estimator for NAFF, O(1) stream reducer for the weighted phase advance
- bd66adc module-build naff-fma stage 01: Amendment D80 sources the NAFF iteration from Laskar2003 (NF16 to NF23); Hanning error order T^-4
- 2dbe702 module-build naff-fma stage 01: owner confirms D80c
- feaaa80 bibliography: BazzaniCapoaniGiovannozzi2024 and Calleja2021 (kept by the owner for later topics)
- f870742 bibliography: ParkHowellStewart2025 (kept by the owner for later topics)
- 6e3a3cf bibliography: LuqueVillanueva2014 (supplied by the owner for naff-fma I12b)
- 3b49469 module-build naff-fma stage 02: tests (35 testsets, RED) and test plan; spec Amendments D80 sources, D81
- fb6886f module-build naff-fma stage 03: series path and WeightedPhaseAdvance reducer in Indicators (GREEN, coverage 99.3%), bench, docs, impl notes
- 7eddf23 module-build naff-fma stage 03: strike the P10 keys from numerics-defaults.toml (owner)
- 07d80dc module-build naff-fma stage 04: verification (criteria 1 to 16, CPU and GPU), review fixes M1 to M3, near_pi, check_run observable check; spec Amendment D82
- 93b9ba1 module-build naff-fma stage 05: catalog I12 and I12b verified; spec notes, Indicators limits, Experiments obligation, changelog

## Files changed by the module

Short hashes of the naff-fma commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (fb6886f, 07d80dc, 93b9ba1)
- `_config/bibliography/indicators.md` (f870742, 6e3a3cf)
- `_config/bibliography/maps.md` (feaaa80)
- `_config/bibliography/references.bib` (feaaa80, f870742, 6e3a3cf)
- `_config/catalogs/indicators.md` (93b9ba1)
- `_config/decision-log.md` (db3992c, bd66adc, 3b49469, 07d80dc)
- `_config/numerics-defaults.toml` (fb6886f, 7eddf23; net unchanged, the P10 keys were added and struck)
- `bench/naff.jl` (fb6886f)
- `docs/src/api/indicators.md` (fb6886f)
- `src/Experiments/CONTEXT.md` (93b9ba1)
- `src/Indicators/CONTEXT.md` (fb6886f, 07d80dc, 93b9ba1)
- `src/Indicators/Indicators.jl` (fb6886f)
- `src/Indicators/naff.jl` (fb6886f, 07d80dc)
- `src/Indicators/naff_fft.jl` (fb6886f)
- `src/Indicators/naff_map.jl` (fb6886f, 07d80dc)
- `src/Indicators/naff_reducer.jl` (fb6886f, 07d80dc)
- `src/precompile.jl` (fb6886f)
- `test/naff-fma/nf_ac01_03.jl` (3b49469)
- `test/naff-fma/nf_ac04_06.jl` (3b49469)
- `test/naff-fma/nf_ac07_09.jl` (3b49469)
- `test/naff-fma/nf_ac10_13.jl` (3b49469)
- `test/naff-fma/nf_api.jl` (3b49469)
- `test/naff-fma/nf_d80.jl` (3b49469)
- `test/naff-fma/nf_d81.jl` (3b49469)
- `test/naff-fma/nf_d81_helpers.jl` (3b49469)
- `test/naff-fma/nf_helpers.jl` (3b49469)
- `test/naff-fma/nf_records.jl` (3b49469)
- `test/naff-fma/nf_review.jl` (07d80dc)
- `test/naff-fma/nf_selftest.jl` (3b49469)
- `test/naff-fma/runtests.jl` (3b49469, 07d80dc)
- `verify/gpu.jl` (07d80dc)
- `verify/naff-fma.jl` (07d80dc)
- `stages/01-spec/output/naff-fma-spec.md` (db3992c, bd66adc, 2dbe702, 3b49469, 07d80dc, 93b9ba1)
- `stages/02-tests/output/naff-fma-test-plan.md` (3b49469)
- `stages/03-implement/output/naff-fma-impl-notes.md` (fb6886f, 7eddf23, 07d80dc)
- `stages/04-verify/output/gpu-naff-verify-log.txt` (07d80dc)
- `stages/04-verify/output/naff-fma-oracle.png` (07d80dc)
- `stages/04-verify/output/naff-fma-spectra.png` (07d80dc)
- `stages/04-verify/output/naff-fma-verification.md` (07d80dc)
- `stages/04-verify/output/naff-fma-verify-log.txt` (07d80dc)
- `stages/05-integrate/output/naff-fma-integration-log.md` (the commit that adds this log)
