# Integration log: tangent-core

Branch `m2-p1-modules`, M2 Batch A (decision D22). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

The flips were written by the orchestrator (promotion write into `_config/catalogs/`); this run did not touch the catalogs.

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L7 (reorthonormalization interval, parameter not module) | planned -> verified |
| `_config/catalogs/lyapunov.md` | section D2, status line (analytic tangent source) | added: verified |

Evidence: `stages/04-verify/output/tangent-core-verification.md` (62 PASS, 0 FAIL in the CPU run; GPU section from `verify/gpu.jl`). Decision rows: D17, D22, D23, D25, D28, D29, D33.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/tangent-core.jl` (e93ce8f)
- `ext/ChaoslabCUDAExt.jl` (ace6f60)
- `src/Tangent/Tangent.jl` (e93ce8f, ace6f60)
- `src/Tangent/initial.jl` (e93ce8f)
- `src/Tangent/orbit.jl` (e93ce8f, 55e5987)
- `src/Tangent/propagate.jl` (e93ce8f, 1361767, 55e5987, ace6f60)
- `test/tangent-core/runtests.jl` (08fb0f9, 55e5987, ace6f60)
- `test/tangent-core/tc_ac01_04.jl` (08fb0f9)
- `test/tangent-core/tc_ac05_08.jl` (08fb0f9)
- `test/tangent-core/tc_ac09_11.jl` (08fb0f9, ace6f60)
- `test/tangent-core/tc_device.jl` (ace6f60)
- `test/tangent-core/tc_helpers.jl` (08fb0f9, ace6f60)
- `test/tangent-core/tc_state_length.jl` (55e5987)
- `verify/gpu.jl` (ace6f60)
- `verify/tangent-core.jl` (64b5f92, ace6f60)
- `stages/01-spec/output/tangent-core-spec.md` (2a21173, 55e5987, ace6f60)
- `stages/02-tests/output/tangent-core-test-plan.md` (08fb0f9, 55e5987, ace6f60)
- `stages/03-implement/output/tangent-core-impl-notes.md` (e93ce8f, 1361767, 55e5987, ace6f60)
- `stages/04-verify/output/gpu-verify-log.txt` (ace6f60)
- `stages/04-verify/output/tangent-core-ac4-bigfloat-product.png` (64b5f92)
- `stages/04-verify/output/tangent-core-oracle-fullbigfloat.png` (64b5f92)
- `stages/04-verify/output/tangent-core-verification.md` (64b5f92, 72565b7, ace6f60)
- `stages/04-verify/output/tangent-core-verify-log.txt` (64b5f92)

- Foundation pieces this module created (health error, compensated sums, backend resolution, device hooks): `src/Foundation/Foundation.jl` (e93ce8f, 1361767, ace6f60). Shared: `Project.toml`, `Manifest.toml`, `Manifest-v1.13.toml`, `verify/Project.toml`, `verify/Manifest.toml` (33f49bc, ace6f60). Docs this module's decisions changed: `_config/architecture.md`, `_config/decision-log.md`, `_config/environment.md`, `_config/numerics-policy.md` (see `git log main..HEAD`).
- New in stage 05: `src/Tangent/CONTEXT.md`, `src/Foundation/CONTEXT.md`, `ext/CONTEXT.md`, Amendment D30 in `tangent-core-spec.md`.

## Docs written in stage 05

- Folder `CONTEXT.md` files for Foundation, Systems, Tangent, Diagnostics and `ext/`; `CHANGELOG.md` section "Added (M2 Batch A)", "Changed", "Deferred"; this log.

Commit SHA: (filled by orchestrator)
