# Integration log: lyap-householder-qr

Branch `m2-p1-modules`, M2 Batch B (decision D43). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L5 (Householder) | planned -> verified |
| `_config/catalogs/lyapunov.md` | L6 (Givens) | planned -> verified |
| `_config/catalogs/lyapunov.md` | L5b (H2; H1 deferred) | planned -> verified |

Evidence: `stages/04-verify/output/lyap-householder-qr-verification.md`, accepted at the Batch B stage 04 checkpoint (2026-10-16). Decision rows: D30, D34, D35, D43.

## Files changed by the module across M2

Short hashes are the commits on `main..m2-p1-modules` that touched the file (`git log --name-only`). Stage paths are relative to `workspaces/module-build/`.

- `bench/lyap-householder-qr.jl` (0d6e276)
- `src/Lyapunov/householder.jl` (0d6e276)
- `src/Lyapunov/lapack_qr.jl` (0d6e276)
- `test/lyap-householder-qr/hh_backend.jl` (629be54, f9dc0b5)
- `test/lyap-householder-qr/hh_cross.jl` (629be54, 0d6e276, f9dc0b5)
- `test/lyap-householder-qr/hh_exact.jl` (629be54, f9dc0b5)
- `test/lyap-householder-qr/hh_factor.jl` (629be54, 0d6e276, f9dc0b5)
- `test/lyap-householder-qr/hh_helpers.jl` (629be54, f9dc0b5)
- `test/lyap-householder-qr/hh_panel.jl` (629be54, 0d6e276, f9dc0b5)
- `test/lyap-householder-qr/hh_runtime.jl` (629be54, 0d6e276, f9dc0b5)
- `test/lyap-householder-qr/runtests.jl` (629be54)
- `verify/lyap-householder-qr.jl` (f9dc0b5, ece423c)
- `stages/01-spec/output/lyap-householder-qr-spec.md` (9c41595, 629be54)
- `stages/02-tests/output/lyap-householder-qr-test-plan.md` (629be54, f9dc0b5)
- `stages/03-implement/output/lyap-householder-qr-impl-notes.md` (0d6e276)
- `stages/04-verify/output/lyap-householder-qr-factor-error.png` (ece423c)
- `stages/04-verify/output/lyap-householder-qr-theta-psi.png` (ece423c)
- `stages/04-verify/output/lyap-householder-qr-verification.md` (ece423c, 68d7172)
- `stages/04-verify/output/lyap-householder-qr-verify-log.txt` (ece423c)

- Shared files touched on this branch (all slugs, all commits): `src/Chaoslab.jl` (e93ce8f, 76f3929, a2cc052), `src/Lyapunov/Lyapunov.jl` (76f3929, d2ca183, 0d6e276), `src/precompile.jl` (a2cc052), `verify/gpu.jl` (ace6f60, ece423c), `verify/threads.jl` (ece423c), `ext/ChaoslabCUDAExt.jl` (ace6f60, 6e2af63); `src/Lyapunov/CONTEXT.md` (new, stage 05).

## Docs written in stage 05

- `src/Lyapunov/CONTEXT.md` (all five Batch B slugs); `CHANGELOG.md` section "Added (M2 Batch B)" and "Deferred"; spec Amendments D43; this log.

Commit SHA: 6ad3423
