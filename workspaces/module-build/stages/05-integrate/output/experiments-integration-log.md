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

## Amendment D86 (stage 05)

Date 2026-10-29. Third scope of the slug: recording of the six indicator kinds by the study runner (decisions D86, D87, D88; spec `experiments-spec.md`, "Amendment D86" and the sections after it). Worked and committed on `main` (owner decision at the stage 05 checkpoint, 2026-10-29: no branch, as for stages 01 to 04 of this amendment; commit f731b16). Milestone M3 (decision D49), stage 05 of `workspaces/module-build`, third integration of the slug. Commit SHA: (set when committed).

### Catalog rows flipped

None. The D86 recording adds no catalog row: `experiments` is not a catalog method or system, and the rows it records (I2, I5, I8, I9, I10, I12, I12b, I13) are already `verified` on their own module evidence, with no qualifier that this amendment changes (`_config/catalogs/indicators.md`). `_config/risks.md` is reviewed at the milestone exit, not here (D86z).

Evidence: `stages/04-verify/output/experiments-d86-verification.md`, accepted at the stage 04 checkpoint (2026-10-29) after the review fixes F1, R1, R2 (D88), R3 and the owner decisions of D88. Decision rows: D86, D87, D88.

### Evidence summary

Quoted from the verification report; the numbers are those of the final code unless noted.

| Item | Result |
|------|--------|
| CPU oracles, `verify/experiments-indicators.jl`, 23 threads | first full run PASS 113, FAIL 0, INFO 63 (on 5b1f91d with the F1 workaround); rerun on the fixed code PASS 113, FAIL 1 (the line that expected the F1 abort), INFO 62, 343 s, peak RSS 1.21 GiB; the reader accepts the D88 layout |
| Device, `verify/experiments-indicators-device.jl`, RTX 5060 Laptop | PASS 32, FAIL 0, FALLBACK 4 (the Double64 passes), INFO 14; peak RSS 2.62 GiB |
| Reduced campaign (D86y) | 2793 items, 352 chunks, 9 arms; 108.2 s for `run_campaign` (limit 900 s); peak host memory 1.01 GiB (limit 8 GB); failed chunks exactly the injected one; resume equals an uninterrupted run (0 differences); CPU rerun 180 of 180 bitwise |
| Independent reader (O20) | 352 of 352 files in each of 4 stores, 0 failed checks; the module's `read_chunk` equals it on 47812 attributes and 33977 datasets, 0 differ; `h5dump` agrees on 49 datasets |
| Plan stability (D86b) | plan hash 81b04cdf5fd6cfd2 equal for the code of ba21bcc and the current code; 15569 of 15569 item hashes and 1363 of 1363 chunk keys and hashes equal |
| Package suites at the stage 04 checkpoint | experiments folder 6844 pass, 0 fail, 0 error; full `Pkg.test()` 1,400,156 pass, 0 fail, 0 error, 1 broken (outside experiments), 15 min |
| Coverage | `src/Experiments` 96.1 percent of 3440 lines (quoted from the orchestrator, taken before the review fixes; new `run_indicator*.jl` files 94 to 100 percent) |
| Review | 0 CRITICAL, 0 HIGH, 3 MEDIUM (R1, R2, R3: fixed), 9 LOW (R4 settled as D88; R5, R8, R9 fixed after the verification report; R6, R7, R10, R11, R12 open) |
| Criteria D86a to D86y | all PASS or covered by the suite ((S) rows D86f, D86r, D86x); D86z is this stage |

### Checks at integration

| Check | Result |
|-------|--------|
| `Pkg.test()` after the fixes R5, R8, R9 (f731b16) | 1,400,166 pass, 0 fail, 0 error, 1 broken (outside experiments); quality 12 of 12; 17 min |
| experiments folder after the fixes R5, R8, R9 | 6854 pass, 0 fail, 0 error |
| ICM validator | ALL CHECKS PASS (`python scripts/icm_validate.py`) |
| Audit | status matches evidence: no catalog row changes and no open criterion beyond the LOW findings below; this log names the files (the commit SHA is set at the commit); `CHANGELOG.md` has the entry under Unreleased, Added |

### Documents changed by this stage (D86z)

- `schema/run-store.md`: section "Added by Amendment D86 (version stays 2)": root attributes of a spec file, `provenance/<member>`, the datasets under `indicators/<id>/`, `pooled/<id>` with `range` (D88), `lyap` optional, `kind` marks the new files; version decision P34 (additive, precedent D56).
- `_config/architecture.md`, "Compute", GPU record: the per-spec entry `[<device>."Float64".recorder.<spec id>]` and its `spec_sha256` field.
- `src/Experiments/CONTEXT.md`: section "Indicators in the runner (D86)" says what is implemented and in which files, the CPU-only kinds (sixth fallback reason) and the limits that remain. The heading "Obligations from the Indicators layer" named by D86z no longer exists there: stage 03 had already replaced it.
- `src/Indicators/CONTEXT.md`: the Limits lines that named missing Experiments paths (chain recording, stretching workspace, naff-fma series feed) point to `src/Experiments/CONTEXT.md`.
- `stages/01-spec/output/naff-fma-spec.md`: criterion 16 carries a status note, gated by D86j.
- `CHANGELOG.md`: one entry under Unreleased, Added.

### Open items

- LOW review findings, not fixed: R6 and R7 (the per-spec cost report and the cost cap judge the first group of an arm only and pool precisions), R10 (the pilot report keeps the last group when a spec has several), R11 (`_stretch_D` restates the D formula of the reducer; guarded by the bitwise test D86i), R12 (the stretching workspace is allocated per orbit, not per slot: garbage, not wrong data).
- Fixed by the orchestrator after the verification report and still uncommitted: R5 (plan-time check that the derived stretching range exists), R8 (the write-time read-back compares `/pooled/<id>`), R9 (the indicator manifest table is cached per plan); their suites were running when this log was written.
- F2 (design limit): the pair on the Y5 control overflows near step 770 at N = 4000; the p2 preregistration must keep the pair horizon below about 700 for Y5 or leave Y5 out of the pair arm.
- Double64 on the device is untested: all four device-capable kinds fall back to the CPU with reason 3, so only the fallback is exercised.
- The device path was measured on the local RTX 5060 only; per-spec cost records vary between takes (tangent 1.4 to 3.1 times slower than the CPU, phase advance 0.6 to 1.1), so `auto` for phase advance is not stable and the study takes its own record.
- Coverage was not rerun after the review fixes.
- The independent reader is HDF5.jl in a separate process plus `h5dump`, not h5py; the command line entry (`bin/`) was not run for the new kinds.
- Full-suite numbers of the final working tree: 1,400,166 pass, 1 broken (above).
- `_config/risks.md` is reviewed at the milestone exit.

### Files changed since the previous log

Short hashes are the commits on `main` from 38eee93 to 366677f that touched the file: 38eee93 stage 01, ba21bcc stage 02, 2e1c100 and 5b1f91d stage 03, 62fd6c2 and 366677f stage 04. "working tree" means uncommitted at the time of writing (the review fixes R5, R8, R9 and this stage's documents); its commit is f731b16 (stage 05). Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (working tree)
- `_config/architecture.md` (working tree)
- `_config/decision-log.md` (38eee93, 5b1f91d, 366677f)
- `bench/experiments.jl` (2e1c100)
- `schema/run-store.md` (working tree)
- `src/Experiments/CONTEXT.md` (2e1c100, working tree)
- `src/Experiments/Experiments.jl` (2e1c100)
- `src/Experiments/run_backend.jl` (2e1c100)
- `src/Experiments/run_block.jl` (2e1c100)
- `src/Experiments/run_campaign.jl` (2e1c100)
- `src/Experiments/run_device.jl` (2e1c100)
- `src/Experiments/run_groups.jl` (2e1c100)
- `src/Experiments/run_health.jl` (2e1c100)
- `src/Experiments/run_indicator.jl` (2e1c100)
- `src/Experiments/run_indicator_chain.jl` (2e1c100, working tree)
- `src/Experiments/run_indicator_device.jl` (2e1c100)
- `src/Experiments/run_indicator_health.jl` (2e1c100)
- `src/Experiments/run_indicator_record.jl` (2e1c100, 62fd6c2)
- `src/Experiments/run_indicator_report.jl` (2e1c100, working tree)
- `src/Experiments/run_indicator_store.jl` (2e1c100, 62fd6c2, 366677f, working tree)
- `src/Experiments/run_manifest.jl` (2e1c100)
- `src/Experiments/run_pilot.jl` (2e1c100)
- `src/Experiments/run_plan.jl` (2e1c100, 62fd6c2)
- `src/Experiments/run_record.jl` (2e1c100)
- `src/Experiments/run_resume.jl` (2e1c100)
- `src/Experiments/run_session.jl` (2e1c100)
- `src/Experiments/run_store.jl` (2e1c100, working tree)
- `src/Indicators/CONTEXT.md` (38eee93, working tree)
- `test/experiments/ex_d86_chains.jl` (ba21bcc, 2e1c100)
- `test/experiments/ex_d86_golden.jl` (ba21bcc)
- `test/experiments/ex_d86_grammar.jl` (ba21bcc, 2e1c100)
- `test/experiments/ex_d86_guards.jl` (ba21bcc)
- `test/experiments/ex_d86_helpers.jl` (ba21bcc, 2e1c100, 62fd6c2)
- `test/experiments/ex_d86_review.jl` (62fd6c2, 366677f, working tree)
- `test/experiments/ex_d86_run.jl` (ba21bcc, 2e1c100)
- `test/experiments/ex_d86_store.jl` (ba21bcc, 2e1c100, 366677f)
- `test/experiments/ex_d86_windows.jl` (ba21bcc, 2e1c100, 366677f)
- `test/experiments/ex_inputs.jl` (ba21bcc)
- `test/experiments/runtests.jl` (ba21bcc, 62fd6c2)
- `verify/experiments-indicators-block.toml` (62fd6c2)
- `verify/experiments-indicators-campaign.jl` (62fd6c2)
- `verify/experiments-indicators-closed.jl` (62fd6c2)
- `verify/experiments-indicators-device.jl` (62fd6c2)
- `verify/experiments-indicators-hash.jl` (62fd6c2)
- `verify/experiments-indicators-memory.jl` (62fd6c2)
- `verify/experiments-indicators-plan.jl` (62fd6c2)
- `verify/experiments-indicators-prep.jl` (62fd6c2)
- `verify/experiments-indicators-reader.jl` (62fd6c2, 366677f)
- `verify/experiments-indicators-series.jl` (62fd6c2)
- `verify/experiments-indicators-source.jl` (62fd6c2)
- `verify/experiments-indicators.jl` (62fd6c2)
- `stages/01-spec/output/experiments-spec.md` (38eee93, ba21bcc, 2e1c100, 5b1f91d, 366677f)
- `stages/01-spec/output/naff-fma-spec.md` (working tree)
- `stages/02-tests/output/experiments-test-plan.md` (ba21bcc, 2e1c100)
- `stages/03-implement/output/experiments-d86-impl-notes.md` (2e1c100)
- `stages/03-implement/output/weighted-birkhoff-impl-notes.md` (38eee93)
- `stages/04-verify/output/experiments-d86-device.png` (62fd6c2, 366677f)
- `stages/04-verify/output/experiments-d86-review.md` (62fd6c2)
- `stages/04-verify/output/experiments-d86-verification.md` (62fd6c2, 366677f)
- `stages/04-verify/output/experiments-d86-verify-log-fixed.txt` (62fd6c2, 366677f)
- `stages/04-verify/output/experiments-d86-verify-log.txt` (62fd6c2)
- `stages/04-verify/output/gpu-rtx5060-d86-verify-log-fixed.txt` (62fd6c2, 366677f)
- `stages/04-verify/output/gpu-rtx5060-d86-verify-log.txt` (62fd6c2)
- `stages/05-integrate/output/experiments-integration-log.md` (working tree)

## Amendment D94 (panel builder and run blocks in rounds, study p3-abstention), stage 05

Branch: `main` (commits pushed, no force-push). Decision D94 (stage 01; owner decisions P69 to P86 approved as proposed).

### Catalog rows flipped

None: the amendment adds no catalog method or indicator.

Evidence: `stages/04-verify/output/experiments-d94-verification.md`, accepted at the stage 04 checkpoint (2026-11-01). Criteria D94a to D94r, D94t, D94u and O26 pass; D94s not run (the real test panel needs the study amendment P69; cost estimated 0.18 h, 1.6 GiB); D94v closed by this stage. Full `Pkg.test()` after the review fixes 1403999 pass, 0 fail, 1 broken (existing); Aqua and JET clean; independent review: C1, H1, H2, M1, M2 and five LOW fixed with 16 probes, M3 and four LOW recorded.

### D94v at integration

| Item | Where | Result |
|------|-------|--------|
| Cell naming with rounds, no layout change, P70 (no orbit attribute `round`) | `schema/run-store.md`, "Added by Amendment D94" | done |
| Panel format and "Files" | `src/Experiments/CONTEXT.md` (fourth scope, five file rows, section "Rounds (D94)", 78 lines) | done |
| Study stage `02-panel` notes the new keys | `workspaces/study/stages/02-panel/references/panel-rounds.md`, input row in its `CONTEXT.md` | done |
| Decision log row D94 | `_config/decision-log.md` (7f70aa7) | present |
| `_config/risks.md` review | due at the M4 exit, not at this amendment | not done here |

Spec text (owner decisions M3 and V3): `experiments-spec.md`, "Clarifications at stage 05": the recorded pool size, the binary64 product for `2 pi u_2`, the order of the K = 0 fractions, the class-bit byte encoding, P70, and the recorded deviations L1, L2, L3, L8. No behaviour or hash changes. `python scripts/icm_validate.py`: all checks pass.

### Open (not verified or deferred)

- D94s and O23 to O26 at full size on the real test panel: after the study amendment P69 (machine-readable panel and run blocks).
- V1: a count altered in the last log entry is accepted; the test-side Analysis amendment recomputes the counts of every entry from the ref store.
- V2: the per-orbit seed rule part of O23 is not independent. V4: REV probes run on the fixed code only.
- Byte identity across machines. The command line entry was run for `--round-log` only on fixtures.
- Not done here: the test-side Analysis amendment (HA, HB, HC, P81 counts, the round-log writer), P69, the A1 freeze run and Amendment A1; CI for recent commits not checked.

### Commits

- `7f70aa7` stage 01: D94 spec, decision log row, owner decisions
- `25f3458` stage 02: failing tests (7 files), plan, golden values
- `ab61826` stage 03: five new files, 14 changed, GREEN, recorded test rework (D86x), impl notes
- `5b36b36` stage 04: verification, review fixes, oracle and log scripts, recorded test rework (D94o)
- (this log) `experiments integration log (D94)`: the commit that adds this section and the documentation edits.

### Files changed on `main` for D94 (short hashes of the commits above that touched the file)

- `_config/decision-log.md` (7f70aa7)
- `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (7f70aa7, this log)
- `test/experiments/ex_d94_golden.jl`, `ex_d94_grammar.jl`, `ex_d94_helpers.jl`, `ex_d94_panel.jl`, `ex_d94_plan.jl`, `ex_d94_streams.jl`, `runtests.jl` (25f3458)
- `test/experiments/ex_d94_run.jl` (25f3458, 5b36b36)
- `test/experiments/ex_d86_guards.jl` (ab61826)
- `workspaces/module-build/stages/02-tests/output/experiments-test-plan.md` (25f3458, ab61826)
- `src/Experiments/Experiments.jl`, `cells.jl`, `inputs.jl`, `panel.jl`, `panel_tables.jl`, `rounds_panel.jl`, `rounds_streams.jl`, `run_block.jl`, `run_campaign.jl`, `run_cli.jl`, `run_groups.jl`, `run_manifest.jl`, `run_resume.jl`, `run_store.jl` (ab61826)
- `src/Experiments/inputs_rounds.jl`, `round_log.jl`, `run_plan.jl`, `run_rounds.jl`, `run_session.jl` (ab61826, 5b36b36)
- `src/Experiments/CONTEXT.md` (this log)
- `bench/experiments.jl` (ab61826)
- `workspaces/module-build/stages/03-implement/output/experiments-d94-impl-notes.md` (ab61826)
- `verify/experiments-d94.jl`, `experiments-d94-panels.jl`, `experiments-d94-screen.jl`, `experiments-d94-review.jl`, `experiments-d94-reduced.toml`, `experiments-d94-oracle.py`, `experiments-d94-log.py` (5b36b36)
- `workspaces/module-build/stages/04-verify/output/experiments-d94-verification.md`, `experiments-d94-verify-log.txt`, `experiments-d94-verify-log-fixed.txt`, `experiments-d94-oracle-log.txt` (5b36b36)
- `schema/run-store.md` (this log)
- `workspaces/study/stages/02-panel/CONTEXT.md`, `references/panel-rounds.md` (this log)
- `CHANGELOG.md` (this log)
- `workspaces/module-build/stages/05-integrate/output/experiments-integration-log.md` (this log)

### Audit

| Check | Result |
|-------|--------|
| Status matches evidence (no catalog row; D94s deferred to P69 and named; findings V1 to V4 recorded) | pass |
| Log names the commit SHAs and every file changed | pass |
| `CHANGELOG.md` has an entry for the amendment | pass |

## Amendment D96 (indicator specs held only by closure arms), stage 05

Branch: `main` (commits pushed, no force-push). Decision D96 (stage 01; owner decision P116 approved as proposed).

### Catalog rows flipped

None: the amendment adds no catalog method or indicator.

Evidence: `stages/04-verify/output/experiments-d96-verification.md`, accepted at the stage 04 checkpoint (2026-11-03). Criteria D96a to D96d and O27 pass; full `Pkg.test()` 1404618 pass, 0 fail, 1 broken (existing); Aqua and JET clean; independent review: no CRITICAL, HIGH or MEDIUM, one LOW recorded (empty table before closure, consistent with R28). Test correction at stage 03 (D96a budget line, owner decision), recorded in the test plan. `python scripts/icm_validate.py`: all checks pass.

### Open (not verified or deferred)

- The count side of O27 ran with the scratch Analysis patches of the discovery (cell form M1 and the tg horizon); the Analysis amendment D97 reruns `verify/experiments-d96*.jl` without them, which closes O24 on the Analysis side.
- Byte identity across machines; one machine only. `_config/risks.md` review: due at the M4 exit.

### Commits

- `8c3bac0` stage 01: R28, I38, I39, O27, D96a to D96d, P116; decision log row
- `ee66bc6` stage 02: failing tests (`ex_d96.jl`), test plan section
- `c41d085` stage 03: the guard in `_indicator_manifest`, D96a budget line corrected, impl notes
- `8c34f37` stage 04: O27 scripts, verification report, verify log
- (this log) `experiments integration log (D96)`: the commit that adds this section and the documentation edits.

### Files changed on `main` for D96

- `_config/decision-log.md` (8c3bac0)
- `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (8c3bac0)
- `test/experiments/ex_d96.jl` (ee66bc6, c41d085), `test/experiments/runtests.jl` (ee66bc6)
- `workspaces/module-build/stages/02-tests/output/experiments-test-plan.md` (ee66bc6, c41d085)
- `src/Experiments/run_indicator_report.jl` (c41d085)
- `workspaces/module-build/stages/03-implement/output/experiments-d96-impl-notes.md` (c41d085)
- `verify/experiments-d96.jl`, `verify/experiments-d96-cli.sh` (8c34f37)
- `workspaces/module-build/stages/04-verify/output/experiments-d96-verification.md`, `experiments-d96-verify-log.txt` (8c34f37)
- `src/Experiments/CONTEXT.md`, `CHANGELOG.md`, `workspaces/module-build/stages/05-integrate/output/experiments-integration-log.md` (this log)

### Audit

| Check | Result |
|-------|--------|
| Status matches evidence (no catalog row; the patched count side named and deferred to D97) | pass |
| Log names the commit SHAs and every file changed | pass |
| `CHANGELOG.md` has an entry for the amendment | pass |
