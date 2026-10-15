# Integration log: lyap-svd-stabilized

Branch `m2-p1-modules`, M2 Batch B (decision D43). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L9 (d = 2) | planned -> verified |

Evidence: `stages/04-verify/output/lyap-svd-stabilized-verification.md`, accepted at the Batch B stage 04 checkpoint (2026-10-16). Decision rows: D30, D34, D35, D43, D41, D40.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/lyap-svd-stabilized.jl` (0d6e276)
- `src/Lyapunov/svd_spectrum.jl` (0d6e276)
- `test/lyap-svd-stabilized/runtests.jl` (629be54)
- `test/lyap-svd-stabilized/svd_ac01_05.jl` (629be54)
- `test/lyap-svd-stabilized/svd_ac06_12.jl` (629be54)
- `test/lyap-svd-stabilized/svd_helpers.jl` (629be54)
- `verify/lyap-svd-stabilized.jl` (ece423c)
- `stages/01-spec/output/lyap-svd-stabilized-spec.md` (9c41595, ed77bf5, 629be54, ece423c)
- `stages/02-tests/output/lyap-svd-stabilized-test-plan.md` (629be54, ece423c)
- `stages/03-implement/output/lyap-svd-stabilized-impl-notes.md` (0d6e276)
- `stages/04-verify/output/lyap-svd-stabilized-exponent-error.png` (ece423c)
- `stages/04-verify/output/lyap-svd-stabilized-verification.md` (ece423c)
- `stages/04-verify/output/lyap-svd-stabilized-verify-log.txt` (ece423c)

- Shared files touched on this branch (all slugs, all commits): `src/Chaoslab.jl` (e93ce8f, 76f3929, a2cc052), `src/Lyapunov/Lyapunov.jl` (76f3929, d2ca183, 0d6e276), `src/precompile.jl` (a2cc052), `verify/gpu.jl` (ace6f60, ece423c), `verify/threads.jl` (ece423c), `ext/ChaoslabCUDAExt.jl` (ace6f60, 6e2af63); `src/Lyapunov/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- `src/Lyapunov/CONTEXT.md` (all five Batch B slugs); `CHANGELOG.md` section "Added (M2 Batch B)" and "Deferred"; spec Amendments D43; this log.

Commit SHA: 6ad3423
