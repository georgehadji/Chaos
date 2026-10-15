# Integration log: lyap-two-trajectory

Branch `m2-p1-modules`, M2 Batch B (decision D43). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L2 (two-trajectory) | planned -> verified |

Evidence: `stages/04-verify/output/lyap-two-trajectory-verification.md`, accepted at the Batch B stage 04 checkpoint (2026-10-16). Decision rows: D30, D34, D35, D43, D42.

Open: criterion 13 sent back to stage 01 (spec Amendment D43); criterion 17 regular part accepted in the D42 configuration, recorded as not discriminating.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/lyap-two-trajectory.jl` (0d6e276)
- `src/Lyapunov/separation.jl` (0d6e276)
- `src/Lyapunov/two_trajectory.jl` (0d6e276)
- `test/lyap-two-trajectory/l2_analytic.jl` (629be54, f9dc0b5)
- `test/lyap-two-trajectory/l2_helpers.jl` (629be54, f9dc0b5)
- `test/lyap-two-trajectory/l2_runtime.jl` (629be54, f9dc0b5)
- `test/lyap-two-trajectory/l2_sweeps.jl` (629be54, f9dc0b5)
- `test/lyap-two-trajectory/runtests.jl` (629be54, f9dc0b5)
- `verify/lyap-two-trajectory.jl` (f9dc0b5, ece423c)
- `stages/01-spec/output/lyap-two-trajectory-spec.md` (9c41595, 629be54, ece423c, 68d7172, 6ad3423)
- `stages/02-tests/output/lyap-two-trajectory-test-plan.md` (629be54, f9dc0b5)
- `stages/03-implement/output/lyap-two-trajectory-impl-notes.md` (0d6e276)
- `stages/04-verify/output/lyap-two-trajectory-oracle.png` (ece423c)
- `stages/04-verify/output/lyap-two-trajectory-ushape.png` (ece423c)
- `stages/04-verify/output/lyap-two-trajectory-verification.md` (ece423c, 68d7172)
- `stages/04-verify/output/lyap-two-trajectory-verify-log.txt` (ece423c)

- Shared files touched on this branch (all slugs, all commits): `src/Chaoslab.jl` (e93ce8f, 76f3929, a2cc052), `src/Lyapunov/Lyapunov.jl` (76f3929, d2ca183, 0d6e276), `src/precompile.jl` (a2cc052), `verify/gpu.jl` (ace6f60, ece423c), `verify/threads.jl` (ece423c), `ext/ChaoslabCUDAExt.jl` (ace6f60, 6e2af63); `src/Lyapunov/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- `src/Lyapunov/CONTEXT.md` (all five Batch B slugs); `CHANGELOG.md` section "Added (M2 Batch B)" and "Deferred"; spec Amendments D43; this log.

Commit SHA: 6ad3423
