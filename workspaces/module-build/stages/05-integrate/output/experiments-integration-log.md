# Integration log: experiments

Branch `module/experiments` (decision D47: one branch per change), milestone M3 (decision D49). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

None. `experiments` is not a catalog method or system (spec, Scope); its state lives in this log, the verification report and `build-order.md` (M3 row).

Evidence: `stages/04-verify/output/experiments-verification.md`, accepted at the stage 04 checkpoint (2026-10-17) after the D52 fixes; rerun PASS 93, FAIL 0. Decision rows: D49, D50 (withdrawn), D51, D52. Study panel accepted at the study/02 checkpoint (2026-10-17).

Open: the runner and run store (later amendment of this slug, before study/03, D49); `Manifest.toml` for Julia 1.10 needs `Pkg.resolve()` on 1.10 (CI `lts` job).

## Files changed on the branch

Short hashes are the commits on `main..module/experiments` that touched the file. Stage paths are relative to `workspaces/module-build/`.

- `Manifest-v1.13.toml` (3991155)
- `Project.toml` (3991155)
- `_config/build-order.md` (9279977)
- `_config/decision-log.md` (9279977, 3991155, 7106c42)
- `_config/research-ideas.md` (7106c42)
- `bench/experiments.jl` (3991155)
- `src/Chaoslab.jl` (3991155)
- `src/Experiments/Experiments.jl` (3991155)
- `src/Experiments/cells.jl` (3991155, 7106c42)
- `src/Experiments/draw.jl` (3991155, 7106c42)
- `src/Experiments/geometry.jl` (3991155)
- `src/Experiments/inputs.jl` (3991155, 7106c42)
- `src/Experiments/newton.jl` (3991155)
- `src/Experiments/orbits.jl` (3991155, 7106c42)
- `src/Experiments/panel.jl` (3991155, 7106c42)
- `src/Experiments/panel_tables.jl` (3991155, 7106c42)
- `src/Experiments/parallel.jl` (3991155)
- `src/Experiments/periodic.jl` (3991155)
- `src/Experiments/portrait.jl` (3991155)
- `src/Experiments/predicates.jl` (3991155)
- `src/Experiments/prereg.jl` (3991155, 7106c42)
- `src/Experiments/provenance.jl` (3991155)
- `src/Experiments/screening.jl` (3991155)
- `src/Experiments/search.jl` (3991155, 7106c42)
- `src/Experiments/subsets.jl` (3991155)
- `test/experiments/ex_determinism.jl` (ceb5227)
- `test/experiments/ex_draw.jl` (ceb5227)
- `test/experiments/ex_geometry.jl` (ceb5227)
- `test/experiments/ex_helpers.jl` (ceb5227)
- `test/experiments/ex_inputs.jl` (ceb5227, 3991155, 7106c42)
- `test/experiments/ex_panel.jl` (ceb5227, 3991155, 7106c42)
- `test/experiments/ex_periodic.jl` (ceb5227, 3991155, 7106c42)
- `test/experiments/ex_predicates.jl` (ceb5227)
- `test/experiments/ex_screening.jl` (ceb5227)
- `test/experiments/ex_subsets.jl` (ceb5227)
- `test/experiments/runtests.jl` (ceb5227)
- `verify/Manifest.toml` (0e30b63)
- `verify/experiments.jl` (0e30b63, 7106c42, 9ef8bd7)
- `stages/01-spec/output/experiments-spec.md` (9279977, 3991155, 7106c42)
- `stages/02-tests/output/experiments-test-plan.md` (ceb5227, 3991155, 7106c42)
- `stages/03-implement/output/experiments-impl-notes.md` (3991155, 7106c42)
- `stages/04-verify/output/experiments-portraits.png` (0e30b63, f627688)
- `stages/04-verify/output/experiments-verification.md` (0e30b63)
- `stages/04-verify/output/experiments-verify-log.txt` (0e30b63, f627688)
- `workspaces/study/stages/01-design/output/p1-le-benchmark-preregistration.md` (9279977, 3991155)
- `workspaces/study/stages/02-panel/output/p1-le-benchmark-panel.toml` (050988e)
- `workspaces/study/stages/02-panel/output/p1-le-benchmark-portraits.png` (050988e)
- `src/Experiments/CONTEXT.md`, `CHANGELOG.md` (this stage's commit, named in the commit message)
