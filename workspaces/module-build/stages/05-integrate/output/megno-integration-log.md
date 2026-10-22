# Integration log: megno

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I8 (MEGNO) | planned -> verified; main-issue column states the map relation to FLI and L1 (Mestre2011 Eq. A2), the label difference on regular orbits and the caller-set chain position (D72); Maffione2011a added to the sources |

Evidence: `stages/04-verify/output/megno-verification.md`, accepted at the stage 04 checkpoint (2026-10-22): criteria 1 to 15 pass on CPU and on the local GPU in the D68 form; independent review with no CRITICAL, HIGH or MEDIUM finding. Decision rows: D65, D66, D69, D71, D72, D73, D74.

Open: none for the module. Recorded limits: the Double64 log sums at large lifted arguments inherit the trig floor of L1 (spec stage 04 notes, D74 item 7); the product rounding of R3. Follow-ups for Experiments as for fli-ofli (`src/Experiments/CONTEXT.md`, "Obligations from the Indicators layer").

## Commits

- 455f8cf module-build megno stage 05: catalog I8 verified; spec notes (D74 item 7 Double64 trig floor, R1, R3, R4), Indicators limits
- a33e2c1 module-build megno stage 04: verification accepted (criteria 1 to 15, CPU and GPU in the D68 form)
- 3fef0f3 module-build megno stage 03: MEGNO reducer in the Indicators layer (GREEN, coverage 100%)
- 7d375fc module-build megno stage 02: RED tests (26 testsets) and test plan; spec Amendment D74
- 11aaa1d module-build megno stage 01: spec approved (D73)

## Files changed by the module

Short hashes of the megno commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (455f8cf, 3fef0f3)
- `_config/catalogs/indicators.md` (455f8cf)
- `_config/decision-log.md` (7d375fc, 11aaa1d)
- `bench/megno.jl` (3fef0f3)
- `docs/src/api/indicators.md` (3fef0f3)
- `src/Experiments/CONTEXT.md` (3fef0f3)
- `src/Indicators/CONTEXT.md` (455f8cf, 3fef0f3)
- `src/Indicators/Indicators.jl` (3fef0f3)
- `src/Indicators/megno.jl` (455f8cf, 3fef0f3)
- `src/precompile.jl` (3fef0f3)
- `test/megno/mg_ac01_05.jl` (7d375fc)
- `test/megno/mg_ac06_08.jl` (7d375fc)
- `test/megno/mg_ac09_10.jl` (7d375fc)
- `test/megno/mg_ac11_13.jl` (7d375fc)
- `test/megno/mg_api.jl` (7d375fc)
- `test/megno/mg_helpers.jl` (7d375fc)
- `test/megno/runtests.jl` (7d375fc)
- `verify/gpu.jl` (a33e2c1)
- `verify/megno.jl` (a33e2c1)
- `stages/01-spec/output/megno-spec.md` (455f8cf, 7d375fc, 11aaa1d)
- `stages/02-tests/output/megno-test-plan.md` (7d375fc)
- `stages/03-implement/output/megno-impl-notes.md` (3fef0f3)
- `stages/04-verify/output/gpu-megno-verify-log.txt` (a33e2c1)
- `stages/04-verify/output/megno-chaostools.png` (a33e2c1)
- `stages/04-verify/output/megno-growth-laws.png` (a33e2c1)
- `stages/04-verify/output/megno-verification.md` (a33e2c1)
- `stages/04-verify/output/megno-verify-log.txt` (a33e2c1)
