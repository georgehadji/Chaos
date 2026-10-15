# Integration log: lyap-gram-schmidt

Branch `m2-p1-modules`, M2 Batch B (decision D43). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L3 (CGS) | planned -> verified |
| `_config/catalogs/lyapunov.md` | L4 (MGS) | planned -> verified |
| `_config/catalogs/lyapunov.md` | L4b (RGS) | planned -> verified |

Evidence: `stages/04-verify/output/lyap-gram-schmidt-verification.md`, accepted at the Batch B stage 04 checkpoint (2026-10-16). Decision rows: D30, D34, D35, D43.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/lyap-gram-schmidt.jl` (d2ca183)
- `src/Lyapunov/orthonormalize.jl` (d2ca183, 0d6e276)
- `src/Lyapunov/qr_spectrum.jl` (d2ca183, 0d6e276)
- `test/lyap-gram-schmidt/gs_ac01_05.jl` (629be54)
- `test/lyap-gram-schmidt/gs_ac06_09.jl` (629be54)
- `test/lyap-gram-schmidt/gs_ac10_14.jl` (629be54)
- `test/lyap-gram-schmidt/gs_api.jl` (629be54)
- `test/lyap-gram-schmidt/gs_helpers.jl` (629be54)
- `test/lyap-gram-schmidt/runtests.jl` (629be54)
- `verify/lyap-gram-schmidt.jl` (ece423c)
- `stages/01-spec/output/lyap-gram-schmidt-spec.md` (9c41595, 629be54, 6ad3423)
- `stages/02-tests/output/lyap-gram-schmidt-test-plan.md` (629be54)
- `stages/03-implement/output/lyap-gram-schmidt-impl-notes.md` (d2ca183)
- `stages/04-verify/output/lyap-gram-schmidt-oracle.png` (ece423c)
- `stages/04-verify/output/lyap-gram-schmidt-orthogonality.png` (ece423c)
- `stages/04-verify/output/lyap-gram-schmidt-verification.md` (ece423c)
- `stages/04-verify/output/lyap-gram-schmidt-verify-log.txt` (ece423c)

- Shared files touched on this branch (all slugs, all commits): `src/Chaoslab.jl` (e93ce8f, 76f3929, a2cc052), `src/Lyapunov/Lyapunov.jl` (76f3929, d2ca183, 0d6e276), `src/precompile.jl` (a2cc052), `verify/gpu.jl` (ace6f60, ece423c), `verify/threads.jl` (ece423c), `ext/ChaoslabCUDAExt.jl` (ace6f60, 6e2af63); `src/Lyapunov/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- `src/Lyapunov/CONTEXT.md` (all five Batch B slugs); `CHANGELOG.md` section "Added (M2 Batch B)" and "Deferred"; spec Amendments D43; this log.

Commit SHA: 6ad3423
