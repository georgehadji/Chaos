# Integration log: sali

Milestone M4 (D65), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I2 (SALI) | planned -> verified; regular-orbit column corrected for 2D maps (power-law decay, spec F1, A9) and the Float64 window note added (F2) |
| `_config/bibliography/benchmarks.md` | "SALI, regular orbit" | scope narrowed to flows and maps of dimension >= 4; new row "SALI, regular orbit, 2D maps" (Skokos2007, Skokos2003a) |

Evidence: `stages/04-verify/output/sali-verification.md`, accepted at the stage 04 checkpoint (2026-10-21; criterion 18 closed 2026-10-22 with the D68 orbit-divergence term). Decision rows: D65, D66, D67, D68.

Open: none for the module. Follow-ups: the Experiments layer must call `check_frame` on every supplied frame (D67); `scaled_norm` repeats `Lyapunov.safe_norm` until both move to a lower layer (when gali or megno need it); the API page of the docs is near its size limit.

## Commits

- b4c775e module-build sali stage 05: catalog I2 verified; I2 and benchmarks rows corrected for 2D maps (regular SALI decays by a power law)
- 2928402 module-build sali stage 04: verification accepted (criteria 1 to 18; spec Amendment D68: criterion 18 per-backend references and orbit-divergence term)
- 6253473 module-build sali stage 03: new Indicators layer with the SALI reducer (GREEN, coverage 100%)
- 9889b72 module-build sali stage 02: RED tests (23 testsets) and test plan; spec Amendment D67
- 6600705 module-build sali stage 01: spec approved (D66: column-wise indicators share a chain, QR family in a separate pass)

## Files changed by the module

Short hashes of the sali commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (b4c775e, 6253473)
- `_config/bibliography/benchmarks.md` (b4c775e)
- `_config/catalogs/indicators.md` (b4c775e)
- `_config/decision-log.md` (2928402, 9889b72, 6600705)
- `bench/sali.jl` (6253473)
- `docs/make.jl` (6253473)
- `docs/src/api.md` (6253473)
- `docs/src/index.md` (6253473)
- `src/Chaoslab.jl` (6253473)
- `src/Indicators/CONTEXT.md` (b4c775e, 6253473)
- `src/Indicators/Indicators.jl` (6253473)
- `src/Indicators/alignment.jl` (6253473)
- `src/Indicators/sali.jl` (6253473)
- `src/precompile.jl` (6253473)
- `test/sali/runtests.jl` (9889b72)
- `test/sali/sl_ac01_03.jl` (9889b72)
- `test/sali/sl_ac04_06.jl` (9889b72)
- `test/sali/sl_ac07_09.jl` (9889b72)
- `test/sali/sl_ac10_12.jl` (9889b72)
- `test/sali/sl_ac13_15.jl` (9889b72)
- `test/sali/sl_helpers.jl` (9889b72)
- `verify/gpu.jl` (2928402)
- `verify/sali.jl` (2928402)
- `stages/01-spec/output/sali-spec.md` (2928402, 9889b72, 6600705)
- `stages/02-tests/output/sali-test-plan.md` (9889b72)
- `stages/03-implement/output/sali-impl-notes.md` (6253473)
- `stages/04-verify/output/gpu-sali-d68-verify-log.txt` (2928402)
- `stages/04-verify/output/gpu-sali-d68b-verify-log.txt` (2928402)
- `stages/04-verify/output/gpu-sali-verify-log.txt` (2928402)
- `stages/04-verify/output/sali-chaostools.png` (2928402)
- `stages/04-verify/output/sali-decay-laws.png` (2928402)
- `stages/04-verify/output/sali-verification.md` (2928402)
- `stages/04-verify/output/sali-verify-log.txt` (2928402)
