# Integration log: fli-ofli

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`. The module delivers FLI (catalog I5) only; OFLI (I6) moved to M5 (D71).

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I5 (FLI) | planned -> verified; main-issue column notes the relation to L1 (running maximum of the L1 log sum), the label difference on regular orbits and the caller-set chain position (D72); sources Skokos2010, Maffione2011a added |
| `_config/catalogs/indicators.md` | I6 (OFLI) | M4 -> M5 (D71, stage 01) |

Evidence: `stages/04-verify/output/fli-ofli-verification.md`, accepted at the stage 04 checkpoint (2026-10-22): criteria 1 to 16 pass on CPU and on the local GPU in the D68 form; independent review with no CRITICAL or HIGH finding. Decision rows: D65, D66, D69, D71, D72.

Open: none for the module. Follow-ups for Experiments (recorded in `src/Experiments/CONTEXT.md`, "Obligations from the Indicators layer"): run `check_frame` on supplied frames (D67) and build and record chains under the one-replacer rule (D71, D72; review R1, R2). The docs API page is near its size limit (153 KiB of 200).

## Commits

- f90115f module-build fli-ofli stage 05: catalog I5 verified; chain-position warnings (review R1, R2), Experiments obligations, docstring and spec notes (R5, R6, R7)
- 1bea614 module-build fli-ofli stage 04: verification accepted (criteria 1 to 16, CPU and GPU in the D68 form)
- f7492e6 module-build fli-ofli stage 03: FLI reducer in the Indicators layer (GREEN, coverage 100%)
- d579bca module-build fli-ofli stage 02: RED tests (27 testsets) and test plan; spec Amendment D72
- afd5cca module-build fli-ofli stage 01: spec approved (D71: FLI only, OFLI to M5, one column replacer per chain)

## Files changed by the module

Short hashes of the fli-ofli commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (f90115f, f7492e6)
- `_config/build-order.md` (afd5cca)
- `_config/catalogs/indicators.md` (f90115f, afd5cca)
- `_config/decision-log.md` (d579bca, afd5cca)
- `bench/fli.jl` (f7492e6)
- `docs/src/api.md` (f7492e6)
- `src/Experiments/CONTEXT.md` (f90115f)
- `src/Indicators/CONTEXT.md` (f90115f, f7492e6)
- `src/Indicators/Indicators.jl` (f90115f, f7492e6)
- `src/Indicators/fli.jl` (f90115f, f7492e6)
- `src/precompile.jl` (f7492e6)
- `test/fli-ofli/fl_ac01_05.jl` (d579bca)
- `test/fli-ofli/fl_ac06_08.jl` (d579bca)
- `test/fli-ofli/fl_ac09_11.jl` (d579bca)
- `test/fli-ofli/fl_ac12_14.jl` (d579bca)
- `test/fli-ofli/fl_api.jl` (d579bca)
- `test/fli-ofli/fl_helpers.jl` (d579bca)
- `test/fli-ofli/runtests.jl` (d579bca)
- `verify/fli.jl` (1bea614)
- `verify/gpu.jl` (1bea614)
- `stages/01-spec/output/fli-ofli-spec.md` (f90115f, d579bca, afd5cca)
- `stages/02-tests/output/fli-ofli-test-plan.md` (d579bca)
- `stages/03-implement/output/fli-ofli-impl-notes.md` (f7492e6)
- `stages/04-verify/output/fli-chaostools.png` (1bea614)
- `stages/04-verify/output/fli-growth-laws.png` (1bea614)
- `stages/04-verify/output/fli-ofli-verification.md` (1bea614)
- `stages/04-verify/output/fli-verify-log.txt` (1bea614)
- `stages/04-verify/output/gpu-fli-verify-log.txt` (1bea614)
