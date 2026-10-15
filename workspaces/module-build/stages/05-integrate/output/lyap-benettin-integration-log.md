# Integration log: lyap-benettin

Branch `m2-p1-modules`, M2 Batch B (decision D43). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L1 (Benettin) | planned -> verified |
| `_config/catalogs/systems.md` | Y1 (standard map) | planned -> verified, large-K known answer reported (report section "Y1: large-K known answer") |

Evidence: `stages/04-verify/output/lyap-benettin-verification.md`, accepted at the Batch B stage 04 checkpoint (2026-10-16). Decision rows: D30, D34, D35, D43.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/lyap-benettin.jl` (76f3929)
- `src/Lyapunov/benettin.jl` (76f3929)
- `src/Lyapunov/norm.jl` (76f3929)
- `test/lyap-benettin/lb_ac01_04.jl` (629be54)
- `test/lyap-benettin/lb_ac05_07.jl` (629be54)
- `test/lyap-benettin/lb_ac08_11.jl` (629be54)
- `test/lyap-benettin/lb_ac17_19.jl` (629be54)
- `test/lyap-benettin/lb_helpers.jl` (629be54)
- `test/lyap-benettin/lb_psi.jl` (629be54)
- `test/lyap-benettin/runtests.jl` (629be54)
- `verify/lyap-benettin.jl` (ece423c)
- `stages/01-spec/output/lyap-benettin-spec.md` (9c41595, 629be54, 6ad3423)
- `stages/02-tests/output/lyap-benettin-test-plan.md` (629be54)
- `stages/03-implement/output/lyap-benettin-impl-notes.md` (76f3929)
- `stages/04-verify/output/lyap-benettin-convergence.png` (ece423c)
- `stages/04-verify/output/lyap-benettin-verification.md` (ece423c)
- `stages/04-verify/output/lyap-benettin-verify-log.txt` (ece423c)

- Shared files touched on this branch (all slugs, all commits): `src/Chaoslab.jl` (e93ce8f, 76f3929, a2cc052), `src/Lyapunov/Lyapunov.jl` (76f3929, d2ca183, 0d6e276), `src/precompile.jl` (a2cc052), `verify/gpu.jl` (ace6f60, ece423c), `verify/threads.jl` (ece423c), `ext/ChaoslabCUDAExt.jl` (ace6f60, 6e2af63); `src/Lyapunov/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- `src/Lyapunov/CONTEXT.md` (all five Batch B slugs); `CHANGELOG.md` section "Added (M2 Batch B)" and "Deferred"; spec Amendments D43; this log.

Commit SHA: 6ad3423
