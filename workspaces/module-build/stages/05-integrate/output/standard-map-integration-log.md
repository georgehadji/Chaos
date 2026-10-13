# Integration log: standard-map

Branch `m2-p1-modules`, M2 Batch A (decision D22). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

The flips were written by the orchestrator (promotion write into `_config/catalogs/`); this run did not touch the catalogs.

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/systems.md` | Y5 (linear symplectic maps) | planned -> verified |
| `_config/catalogs/systems.md` | Y1 (standard map) | stays planned, as agreed: lyap-benettin must first verify the large-K `lambda_1` known-answer test |

Evidence: `stages/04-verify/output/standard-map-verification.md` (criteria 1 to 19 and W1 to W8; Y5 criteria 15 to 17 pass). Decision rows: D22, D23, D28, D31.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/standard-map.jl` (e93ce8f, 55e5987)
- `src/Systems/Systems.jl` (e93ce8f, 55e5987)
- `src/Systems/linear_symplectic.jl` (e93ce8f, 1361767)
- `src/Systems/residuals.jl` (e93ce8f, 1361767)
- `src/Systems/standard_map.jl` (e93ce8f, 1361767, 55e5987)
- `src/Systems/wrapped_standard_map.jl` (55e5987)
- `test/standard-map/runtests.jl` (08fb0f9, 55e5987)
- `test/standard-map/sm_ac01_06.jl` (08fb0f9)
- `test/standard-map/sm_ac07_13.jl` (08fb0f9)
- `test/standard-map/sm_ac14.jl` (08fb0f9)
- `test/standard-map/sm_ac18.jl` (08fb0f9)
- `test/standard-map/sm_helpers.jl` (08fb0f9)
- `test/standard-map/sm_wrapped.jl` (55e5987)
- `test/standard-map/sm_y5.jl` (08fb0f9)
- `verify/standard-map.jl` (64b5f92)
- `stages/01-spec/output/standard-map-spec.md` (2a21173, 55e5987)
- `stages/02-tests/output/standard-map-test-plan.md` (08fb0f9, 55e5987)
- `stages/03-implement/output/standard-map-impl-notes.md` (e93ce8f, 1361767, 55e5987)
- `stages/04-verify/output/standard-map-agreement.png` (64b5f92)
- `stages/04-verify/output/standard-map-verification.md` (64b5f92, 72565b7, ace6f60)
- `stages/04-verify/output/standard-map-verify-log.txt` (64b5f92, 72565b7)

- Shared files this module also touched: `src/Chaoslab.jl` (e93ce8f), `test/quality.jl` (e93ce8f); `src/Systems/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- Folder `CONTEXT.md` files for Foundation, Systems, Tangent, Diagnostics and `ext/`; `CHANGELOG.md` section "Added (M2 Batch A)", "Changed", "Deferred"; this log.

Commit SHA: abef540
