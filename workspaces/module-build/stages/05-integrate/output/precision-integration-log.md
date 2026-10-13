# Integration log: precision

Branch `m2-p1-modules`, M2 Batch A (decision D22). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

The flips were written by the orchestrator (promotion write into `_config/catalogs/`); this run did not touch the catalogs.

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | section E (precision variants), status line | added: verified (Float64, Double64, BigFloat; Float32 dropped, D24) |

Evidence: `stages/04-verify/output/precision-verification.md` (26 PASS, 0 FAIL; AC10 dropped by D24). Decision rows: D20, D24, D22.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/precision.jl` (e93ce8f)
- `test/precision/helpers.jl` (08fb0f9, e93ce8f)
- `test/precision/runtests.jl` (08fb0f9, e93ce8f, ace6f60)
- `verify/precision.jl` (64b5f92)
- `stages/01-spec/output/precision-spec.md` (2a21173, 08fb0f9)
- `stages/02-tests/output/precision-test-plan.md` (08fb0f9, e93ce8f)
- `stages/03-implement/output/precision-impl-notes.md` (e93ce8f, 1361767)
- `stages/04-verify/output/precision-orbit-error.png` (64b5f92)
- `stages/04-verify/output/precision-verification.md` (64b5f92, 72565b7, ace6f60)
- `stages/04-verify/output/precision-verify-log.txt` (64b5f92)

- Source files of this module live in `src/Foundation/Foundation.jl` (`precision_label`, `parse_precision`, `exact`, `Double64` re-export; e93ce8f, 1361767) and are described in `src/Foundation/CONTEXT.md` (new, stage 05). Orbit comparisons use `StandardMap` from `src/Systems/`.
- Environment fill: `_config/environment.md` "High precision" (ace6f60).

## Docs written in stage 05

- Folder `CONTEXT.md` files for Foundation, Systems, Tangent, Diagnostics and `ext/`; `CHANGELOG.md` section "Added (M2 Batch A)", "Changed", "Deferred"; this log.

Commit SHA: (filled by orchestrator)
