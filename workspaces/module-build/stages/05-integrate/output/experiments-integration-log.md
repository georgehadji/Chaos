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

## Runner and run store (second scope, D53 to D57)

Date 2026-10-20. Branch `module/experiments-runner` (decision D47: one branch per change), milestone M3 (decision D49). Stage 05 of `workspaces/module-build`, second integration of the slug.

### Catalog rows flipped

None. `experiments` is not a catalog method or system; its state lives in this log, the verification reports and `build-order.md` (M3 row, unchanged: it closes when the study results are accepted).

Evidence: `stages/04-verify/output/experiments-runner-verification.md`, accepted at the stage 04 checkpoint (2026-10-20) after D55 (review fixes), preregistration Amendments A5 and A6 and decisions D56 and D57. Verification run: PASS 80, FAIL 6 in the main run; the report traces each FAIL row (three stale script assumptions, corrected and rerun clean; the L9 control bound, closed by A6; and three rows with no effect on a decision: the L2 recorder record, the L2 per-orbit device agreement and the tau_sweep pilot ratio). Decision rows: D53, D54, D55, D56, D57.

### Checks at integration

| Check | Result |
|-------|--------|
| `Pkg.test()` (Julia 1.13.2, 4 threads) | 54445 of 54445 pass in 5 min 19 s |
| experiments suite inside `Pkg.test` | 4824 of 4824 in 3 min 16 s; the D35 ceiling is 3 min (2 min 50 s before D57 on the same machine; the cause of the 16 s over the ceiling is not diagnosed: D57 adds five small campaigns, and machine load was not controlled); the owner decides whether to raise the ceiling for this suite or trim it |
| ICM validator | ALL CHECKS PASS |
| Preregistration | `check_prereg` accepts 6 amendments, header hash 62011ba5af514da4f1238aef9ae6cb35ca21093c4343b71dd696df360b9994bf; panel hash of A4 still accepted |
| Study data | `data/p1-le-benchmark/` does not exist; no study chunk has been run |

Open: `Manifest.toml` for Julia 1.10 still needs `Pkg.resolve()` on 1.10 (CI `lts` job); the L2 recorder record is refreshed by the K19 rerun at the start of study/03.

### Files changed on the branch

Short hashes are the commits on `main..module/experiments-runner` that touched the file. Stage paths are relative to `workspaces/module-build/`. This stage also changes `CHANGELOG.md`, `src/Experiments/CONTEXT.md` and this log (the stage commit, named in its message).

- `CLAUDE.md` (77452d6)
- `Manifest-v1.13.toml` (8b5dae8)
- `Project.toml` (8b5dae8)
- `_config/architecture.md` (9cd8df6)
- `_config/decision-log.md` (d6b78ca, 4219a42, d787d25, 9cd8df6, 77452d6)
- `_config/risks.md` (9cd8df6)
- `bench/experiments.jl` (8b5dae8)
- `bin/chaoslab.jl` (8b5dae8)
- `schema/run-store.md` (4219a42, d787d25, 8b5dae8, 77452d6)
- `src/Experiments/CONTEXT.md` (9cd8df6, 4a4873b, 8b5dae8)
- `src/Experiments/Experiments.jl` (4a4873b, 8b5dae8)
- `src/Experiments/run_backend.jl` (9cd8df6, 8b5dae8)
- `src/Experiments/run_block.jl` (d6b78ca, 4219a42, 8b5dae8)
- `src/Experiments/run_campaign.jl` (d787d25, 9cd8df6, 8b5dae8)
- `src/Experiments/run_cli.jl` (d787d25, 8b5dae8)
- `src/Experiments/run_device.jl` (d787d25, 4a4873b)
- `src/Experiments/run_grid.jl` (d787d25, 8b5dae8)
- `src/Experiments/run_groups.jl` (4219a42, d787d25, 4a4873b, 8b5dae8)
- `src/Experiments/run_health.jl` (d6b78ca, 4219a42, 8b5dae8)
- `src/Experiments/run_manifest.jl` (4219a42, d787d25, 9cd8df6, 4a4873b, 8b5dae8)
- `src/Experiments/run_pilot.jl` (4219a42, d787d25, 4a4873b, 8b5dae8)
- `src/Experiments/run_plan.jl` (d6b78ca, 4219a42, d787d25, 8b5dae8)
- `src/Experiments/run_record.jl` (4a4873b, 8b5dae8)
- `src/Experiments/run_resume.jl` (d787d25, 8b5dae8)
- `src/Experiments/run_session.jl` (d787d25, 8b5dae8)
- `src/Experiments/run_store.jl` (d787d25, 8b5dae8)
- `test/experiments/ex_a5.jl` (d6b78ca, 4219a42)
- `test/experiments/ex_a6.jl` (d6b78ca)
- `test/experiments/ex_chunks.jl` (8b5dae8, 0082ef8)
- `test/experiments/ex_device.jl` (4a4873b)
- `test/experiments/ex_health.jl` (8b5dae8, 0082ef8)
- `test/experiments/ex_inputs.jl` (0082ef8)
- `test/experiments/ex_manifest.jl` (8b5dae8, 0082ef8)
- `test/experiments/ex_pilot.jl` (8b5dae8, 0082ef8)
- `test/experiments/ex_plan.jl` (4219a42, 0082ef8)
- `test/experiments/ex_record.jl` (9cd8df6)
- `test/experiments/ex_resume.jl` (0082ef8)
- `test/experiments/ex_review.jl` (d787d25)
- `test/experiments/ex_run_helpers.jl` (d787d25, 8b5dae8, 0082ef8)
- `test/experiments/ex_run_inputs.jl` (4219a42, 0082ef8)
- `test/experiments/ex_store.jl` (0082ef8)
- `test/experiments/runtests.jl` (d6b78ca, 4219a42, d787d25, 9cd8df6, 4a4873b, 0082ef8)
- `verify/Manifest.toml` (05792a3, 4a4873b)
- `verify/Project.toml` (05792a3)
- `verify/experiments-device.jl` (9cd8df6, 4a4873b)
- `verify/experiments-runner-reader.jl` (fa84391, 05792a3)
- `verify/experiments-runner.jl` (d6b78ca, 968ab29, fa84391, 05792a3)
- `stages/01-spec/output/experiments-spec.md` (d6b78ca, 4219a42, 9cd8df6, 77452d6)
- `stages/02-tests/output/experiments-test-plan.md` (d6b78ca, 4219a42, d787d25, 9cd8df6, 4a4873b, 8b5dae8, 0082ef8)
- `stages/03-implement/output/experiments-impl-notes.md` (d6b78ca, 4219a42, d787d25, 9cd8df6, 4a4873b, 8b5dae8)
- `stages/04-verify/output/experiments-runner-cost.png` (b0b51ce, 05792a3)
- `stages/04-verify/output/experiments-runner-verification.md` (d6b78ca, 968ab29, 05792a3)
- `stages/04-verify/output/experiments-runner-verify-log.txt` (d6b78ca, 968ab29, 05792a3)
- `stages/04-verify/output/gpu-rtx5060-d53-verify-log.txt` (9cd8df6, 4a4873b)
- `workspaces/study/stages/01-design/output/p1-le-benchmark-preregistration.md` (d6b78ca, 4219a42, 77452d6)
