# Integration log: diagnostics

Branch `m2-p1-modules`, M2 Batch A (decision D22). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

The flips were written by the orchestrator (promotion write into `_config/catalogs/`); this run did not touch the catalogs.

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | section F (built-in checks every estimator reports), status line | added: verified |

Evidence: `stages/04-verify/output/diagnostics-verification.md` (74 PASS, 0 FAIL; item 7 ChaosTools sum deferred to lyap-benettin 04). Decision rows: D22, D24, D27.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/diagnostics.jl` (e93ce8f)
- `src/Diagnostics/Diagnostics.jl` (e93ce8f)
- `src/Diagnostics/monitor.jl` (e93ce8f, 1361767)
- `src/Diagnostics/residuals.jl` (e93ce8f)
- `test/diagnostics/faultmaps.jl` (08fb0f9)
- `test/diagnostics/runtests.jl` (08fb0f9)
- `verify/diagnostics.jl` (64b5f92)
- `stages/01-spec/output/diagnostics-spec.md` (2a21173)
- `stages/02-tests/output/diagnostics-test-plan.md` (08fb0f9)
- `stages/03-implement/output/diagnostics-impl-notes.md` (e93ce8f, 1361767)
- `stages/04-verify/output/diagnostics-det-residual.png` (64b5f92)
- `stages/04-verify/output/diagnostics-sum-residual.png` (64b5f92)
- `stages/04-verify/output/diagnostics-verification.md` (64b5f92, 72565b7, ace6f60)
- `stages/04-verify/output/diagnostics-verify-log.txt` (64b5f92)

- Patch to Tangent made for this module in stage 04: `reducer_step`, `reducer_failed`, `reducer_error` (`src/Tangent/Tangent.jl`, `orbit.jl`, `propagate.jl`; e93ce8f, 1361767).
- New in stage 05: `src/Diagnostics/CONTEXT.md`; Amendment D30 in `diagnostics-spec.md`.

## Docs written in stage 05

- Folder `CONTEXT.md` files for Foundation, Systems, Tangent, Diagnostics and `ext/`; `CHANGELOG.md` section "Added (M2 Batch A)", "Changed", "Deferred"; this log.

Commit SHA: (filled by orchestrator)
