# Integration log: analysis

Branch `module/analysis` (decision D47: one branch per change), milestone M3 (decision D60). Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

None. `analysis` is not a catalog method or system (spec, Scope); its state lives in this log, the verification report and `build-order.md` (M3 row: modules experiments (D49), analysis (D60)).

Evidence: `stages/04-verify/output/analysis-verification.md`, accepted at the stage 04 checkpoint (2026-10-20): verify log PASS 40, FAIL 0, INFO 2; independent recomputation equal on 883 of 883 orbits; full `Pkg.test` 54916 pass, 0 fail, 0 error; independent review closed (2 CRITICAL and 2 HIGH fixed or closed by a refusal). Decision row: D60.

Open: K values of the strong-chaos set S are a constant of the module (owner decision, spec corrections item 7); chunks that failed twice are refused, not analysed (the study has none); `Manifest.toml` for Julia 1.10 needs `Pkg.resolve()` on 1.10 (CI `lts` job); the experiments suite is over its 3 minute ceiling (D35), as at the end of `experiments`; study stage 04 (running `analyze`, findings) has not been started and needs the final manifest from branch `study/p1-le-benchmark-run` and the parameter file of this branch on one branch.

## Commits on `main..module/analysis`

- `3a9836c` analysis spec (D60): study analysis module between Experiments and Figures; spec with the stage 04 corrections; decision, architecture and build-order rows
- `2ef703b` analysis tests (D60): 434 tests on a synthetic study written through the store writer; layering list gains Analysis; experiments AC28 skips the new layer
- `eedf2c6` analysis implementation (D60): src/Analysis, CairoMakie figures extension, analyze command, bench
- `687dfdf` analysis verification (D60): verify scripts, independent recomputation, report and log; parameter file of study 04
- `e05652a` analysis integration (D60): module CONTEXT, changelog, extension notes
- (this log) `analysis integration log (D60)`: the commit that adds this file.

## Files changed on the branch

Short hashes are the commits above that touched the file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (e05652a)
- `Project.toml` (eedf2c6)
- `_config/architecture.md` (3a9836c)
- `_config/build-order.md` (3a9836c)
- `_config/decision-log.md` (3a9836c)
- `bench/analysis.jl` (eedf2c6)
- `bin/chaoslab.jl` (eedf2c6)
- `ext/CONTEXT.md` (e05652a)
- `ext/ChaoslabCairoMakieExt.jl` (eedf2c6)
- `src/Analysis/Analysis.jl` (eedf2c6)
- `src/Analysis/CONTEXT.md` (e05652a)
- `src/Analysis/checks.jl` (eedf2c6)
- `src/Analysis/cli.jl` (eedf2c6)
- `src/Analysis/figdata.jl` (eedf2c6)
- `src/Analysis/h1.jl` (eedf2c6)
- `src/Analysis/h4.jl` (eedf2c6)
- `src/Analysis/labels.jl` (eedf2c6)
- `src/Analysis/orbits.jl` (eedf2c6)
- `src/Analysis/params.jl` (eedf2c6)
- `src/Analysis/pipeline.jl` (eedf2c6)
- `src/Analysis/sets.jl` (eedf2c6)
- `src/Analysis/stats.jl` (eedf2c6)
- `src/Analysis/store.jl` (eedf2c6)
- `src/Analysis/tables.jl` (eedf2c6)
- `src/Chaoslab.jl` (eedf2c6)
- `test/analysis/an_bootstrap.jl` (2ef703b)
- `test/analysis/an_breakdown.jl` (2ef703b)
- `test/analysis/an_cli.jl` (2ef703b)
- `test/analysis/an_fixture.jl` (2ef703b)
- `test/analysis/an_h1.jl` (2ef703b)
- `test/analysis/an_h4.jl` (2ef703b)
- `test/analysis/an_h5.jl` (2ef703b)
- `test/analysis/an_helpers.jl` (2ef703b)
- `test/analysis/an_labels.jl` (2ef703b)
- `test/analysis/an_params.jl` (2ef703b)
- `test/analysis/an_pipeline.jl` (2ef703b)
- `test/analysis/an_review.jl` (2ef703b)
- `test/analysis/an_sets.jl` (2ef703b)
- `test/analysis/runtests.jl` (2ef703b)
- `test/experiments/ex_inputs.jl` (2ef703b)
- `test/layering.jl` (2ef703b)
- `verify/Manifest.toml` (687dfdf)
- `verify/analysis-independent.jl` (687dfdf)
- `verify/analysis.jl` (687dfdf)
- `workspaces/study/stages/04-analyze/references/p1-le-benchmark-analysis-parameters.toml` (687dfdf)
- `stages/01-spec/output/analysis-spec.md` (3a9836c)
- `stages/02-tests/output/analysis-test-plan.md` (2ef703b)
- `stages/03-implement/output/analysis-impl-notes.md` (eedf2c6)
- `stages/04-verify/output/analysis-independent-results.txt` (687dfdf)
- `stages/04-verify/output/analysis-k0-closed-form.png` (687dfdf)
- `stages/04-verify/output/analysis-verification.md` (687dfdf)
- `stages/04-verify/output/analysis-verify-log.txt` (687dfdf)
- `stages/05-integrate/output/analysis-integration-log.md` (this log)

## Amendment D89 (fourth scope: study p2-indicators), stage 05

Branch: `main` (commits pushed, no force-push). Decisions D89 (spec, stage 01) and D90 (stage 04 checkpoint).

### Catalog rows flipped

None, as for the base module: `analysis` is not a catalog method or system.

Evidence: `stages/04-verify/output/analysis-verification.md` (section "Amendment D89 (stage 04)"), accepted at the stage 04 checkpoint (2026-10-30) with the two owner readings of D90. After the rework: verify log PASS 137, FAIL 0, INFO 88; 14,146 of 14,146 quantities equal to the independent HDF5-only implementation (one count, R2 I_7 H3 `npass`, differs by the one lower end that moved and is a tolerated count); analysis folder suite 1214 pass; full `Pkg.test()` passes; real run 48.8 s at 23 threads, 1.88 GiB peak (limits 15 min and 8 GB); the P1 `analyze` tables file differs from the accepted one in the two header lines only; independent review closed (0 CRITICAL, 1 HIGH fixed, 7 MEDIUM, 4 LOW).

Updated in this stage: `CHANGELOG.md` (entry for the amendment and the schema note), `_config/decision-log.md` (D90), `schema/run-store.md` (P34: the value of a flag bit is 2^(j-1) for the name at position j of `bits`; no change of format), `src/Analysis/CONTEXT.md` (stage 03).

### Open (not verified or deferred)

- The figure drawers for the P2 figures F1 to F6 (P11): `analyze` writes the figure data and exits with code 2, status `figures_skipped`, until the extension has them (decision D89, follow-up).
- Not verified: the input changed during the run (E36 o) beyond a before and after hash; cold-cache timing; draw-by-draw equality of the bootstrap generators (they differ by design); the pilot-store half of D89y (P30: that store is not in the repository); the second versions of T3 to T8, per-orbit values, most by-source and by-K point values, F5, F6, T17 and the layout and digests of the tables file on the independent side.
- The comparison tolerance and the T_s^design reading were settled after the first comparison showed the differences, so the pass after them is not a blind test of those two.
- Unexplained: an `InterruptException` appeared in one analysis folder-suite run of the module-build session and once in a subagent's run; reruns were clean. Nothing in the code raises interrupts.
- Follow-ups for the owner: U1 (runner accepts a subset of the grid), U2 (the P1 labels are not a tracked file), U3 (T_s size risk), the P2 figure drawers; M4 exit tasks (results accepted, `_config/risks.md` review, minor release); the LOW findings of D86 (R6, R7, R10, R11, R12, F2); CI for recent commits has not been checked.

### Commits (stage paths relative to `workspaces/module-build/`)

- `f7f0df2` stage 01: D89 in the decision log, owner decisions of the checkpoint appended to the spec
- `99905e9` stage 02: tests (19 files), plan with criterion map, owner decisions A1 and A2
- `a3ca6dd` stage 03: implementation (p2 files in `src/Analysis`, dispatch on the study slug, p2 CLI options), A1 slug rename in 3 P1 test files
- `9ed8855` stage 02/03 rework: A2 `text` field per parameter key, E33 R2 version recomputed with R2 verdicts
- `b792724` stage 04 (part): independent review fixes (K pooling by K value and others), review tests, independent implementation and results, draft parameter table
- `c3288a3` stage 04: verification (real-store run, transcription of the 58-key table, comparison, report)
- `03a2519` stage 04 rework: T_s^design intersected with E8 in the extension, D89w tolerance, re-verification
- (this log) `analysis integration log (D89)`: the commit that adds this section, with the changelog, decision-log and schema edits.

### Files changed on `main` for D89 (short hashes of the commits above that touched the file)

- `_config/decision-log.md` (f7f0df2, this log)
- `bench/analysis.jl` (a3ca6dd)
- `src/Analysis/Analysis.jl` (a3ca6dd)
- `src/Analysis/CONTEXT.md` (a3ca6dd)
- `src/Analysis/cli.jl` (a3ca6dd)
- `src/Analysis/p2_agree.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_boot.jl` (a3ca6dd)
- `src/Analysis/p2_check.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_dispatch.jl` (a3ca6dd)
- `src/Analysis/p2_hitems.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_hyp.jl` (a3ca6dd, 9ed8855, 03a2519)
- `src/Analysis/p2_items.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_orbits.jl` (a3ca6dd)
- `src/Analysis/p2_params.jl` (a3ca6dd, 9ed8855, b792724)
- `src/Analysis/p2_pipeline.jl` (a3ca6dd, 9ed8855, b792724, 03a2519)
- `src/Analysis/p2_render.jl` (a3ca6dd)
- `src/Analysis/p2_report.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_ritems.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_rules.jl` (a3ca6dd, b792724)
- `src/Analysis/p2_sets.jl` (a3ca6dd)
- `src/Analysis/p2_stat.jl` (a3ca6dd)
- `src/Analysis/p2_store.jl` (a3ca6dd)
- `src/Analysis/pipeline.jl` (a3ca6dd)
- `test/analysis/an_cli.jl` (a3ca6dd)
- `test/analysis/an_helpers.jl` (a3ca6dd)
- `test/analysis/an_p2_agreement.jl` (99905e9)
- `test/analysis/an_p2_bootstrap.jl` (99905e9)
- `test/analysis/an_p2_build.jl` (99905e9, b792724)
- `test/analysis/an_p2_detection.jl` (99905e9)
- `test/analysis/an_p2_dispatch.jl` (99905e9)
- `test/analysis/an_p2_fixture.jl` (99905e9, 9ed8855)
- `test/analysis/an_p2_helpers.jl` (99905e9, 9ed8855, b792724)
- `test/analysis/an_p2_output.jl` (99905e9, b792724)
- `test/analysis/an_p2_params.jl` (99905e9, 9ed8855, b792724)
- `test/analysis/an_p2_pipeline.jl` (99905e9, 9ed8855, b792724, 03a2519)
- `test/analysis/an_p2_ref.jl` (99905e9)
- `test/analysis/an_p2_refstore.jl` (99905e9, 9ed8855, b792724, 03a2519)
- `test/analysis/an_p2_refusals.jl` (99905e9)
- `test/analysis/an_p2_reported.jl` (99905e9, b792724)
- `test/analysis/an_p2_review.jl` (b792724)
- `test/analysis/an_p2_sets.jl` (99905e9)
- `test/analysis/an_p2_stats.jl` (99905e9)
- `test/analysis/an_p2_study.jl` (99905e9, b792724)
- `test/analysis/an_p2_variant.jl` (99905e9)
- `test/analysis/an_p2_verdicts.jl` (99905e9)
- `test/analysis/an_pipeline.jl` (a3ca6dd)
- `test/analysis/runtests.jl` (99905e9, b792724)
- `verify/analysis-p2-decision.py` (b792724)
- `verify/analysis-p2-independent.jl` (b792724)
- `verify/analysis-p2.jl` (c3288a3, 03a2519)
- `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (f7f0df2, 9ed8855, 03a2519)
- `workspaces/module-build/stages/02-tests/output/analysis-test-plan.md` (99905e9, 9ed8855, b792724, 03a2519)
- `workspaces/module-build/stages/03-implement/output/analysis-d89-impl-notes.md` (a3ca6dd, 9ed8855, b792724, 03a2519)
- `workspaces/module-build/stages/04-verify/output/analysis-p2-independent-results.txt` (b792724)
- `workspaces/module-build/stages/04-verify/output/analysis-p2-verify-log.txt` (c3288a3, 03a2519)
- `workspaces/module-build/stages/04-verify/output/analysis-verification.md` (c3288a3, 03a2519)
- `workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml` (b792724)
- `CHANGELOG.md` (this log)
- `schema/run-store.md` (this log)
- `workspaces/module-build/stages/05-integrate/output/analysis-integration-log.md` (this log)

## Amendment D91 (figure drawers F1 to F6 of study p2-indicators), stage 05

Branch: `main` (commits pushed, no force-push). Decision D91 (stage 01; owner decisions P35 to P42 approved as proposed).

### Catalog rows flipped

None, as for the base module.

Evidence: `stages/04-verify/output/analysis-verification.md` (section "Amendment D91: figure drawers F1 to F6 of study p2-indicators (stage 04)"), accepted at the stage 04 checkpoint (2026-10-30) including the owner visual check D91l of the six figures. Verify log PASS 143, FAIL 6, INFO 18: the six fails are the comparison with the tracked study-folder files `p1-le-benchmark-fig-F1` to `F6.pdf`, which date from before the P1 figure redraw (91abcfa, 2026-10-20); a copy of the tree before D91 gives the same P1 pictures as the current module, so D91 does not change P1 (criterion D91e met in substance). Analysis folder suite 1287 pass; full `Pkg.test()` passed at b2d7ebe (stage 03); real P2 command exit 0, 94.7 s, 3.21 GB, figure step 24.0 s cold; independent code review closed (0 CRITICAL, 0 HIGH, 2 MEDIUM fixed, LOW findings fixed or kept as listed in the report).

Updated in this stage: `CHANGELOG.md` (entry for the amendment), `src/Analysis/CONTEXT.md` (exit code 2, result fields, the open P11 reading), `ext/CONTEXT.md` (the P2 drawers).

### Open (not verified or deferred)

- The full `Pkg.test()` was not rerun after the stage 04 review fixes (17 min); they are confined to `ext/p2_figures.jl` and `ext/p2_figures_b.jl`, which the package suite does not load except for the layering checks (folder suite 1287 pass).
- Aqua, JET, coverage and bench for the extension files.
- The reading P11 in the accepted P2 tables file still says that the drawers do not exist (owner decision; the file is accepted).
- The tracked P1 figure files `p1-le-benchmark-fig-F1` to `F6.pdf` in `workspaces/study/stages/04-analyze/output` are stale (V1 of the report); refreshing them is a P1 change for the owner to decide.
- Still to do in the study: draw the P2 figures with the real command into `workspaces/study/stages/04-analyze/output` (study stage 04, audit "regenerable"), then study stage 05.
- Follow-ups for the owner: U1, U2, U3, M4 exit tasks (results accepted, `_config/risks.md` review, minor release), the LOW findings of D86, CI for recent commits has not been checked.

### Commits

- `c2f01e6` stage 01: D91 spec, decision log row, owner decisions
- `72e4463` stage 02: tests (module side, fake drawer, 10 testsets), plan with criterion map
- `b2d7ebe` stage 03: module side, drawers F1 to F6, one test import fixed
- `8ece640` stage 04: verification (independent script, log, contact sheets, report), review fixes in the drawers
- (this log) `analysis integration log (D91)`: the commit that adds this section and the documentation edits.

### Files changed on `main` for D91 (short hashes of the commits above that touched the file)

- `_config/decision-log.md` (c2f01e6)
- `src/Analysis/cli.jl` (b2d7ebe)
- `src/Analysis/p2_pipeline.jl` (b2d7ebe)
- `src/Analysis/pipeline.jl` (b2d7ebe)
- `src/Analysis/CONTEXT.md` (this log)
- `ext/ChaoslabCairoMakieExt.jl` (b2d7ebe)
- `ext/p2_figures.jl` (b2d7ebe, 8ece640)
- `ext/p2_figures_b.jl` (b2d7ebe, 8ece640)
- `ext/CONTEXT.md` (this log)
- `test/analysis/an_p2_figures.jl` (72e4463, b2d7ebe)
- `test/analysis/runtests.jl` (72e4463)
- `verify/analysis-p2-figures.jl` (8ece640)
- `verify/analysis-p2-figures.py` (8ece640)
- `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (c2f01e6)
- `workspaces/module-build/stages/02-tests/output/analysis-test-plan.md` (72e4463)
- `workspaces/module-build/stages/03-implement/output/analysis-d91-impl-notes.md` (b2d7ebe)
- `workspaces/module-build/stages/03-implement/output/analysis-d91-drawers-notes.md` (b2d7ebe)
- `workspaces/module-build/stages/04-verify/output/analysis-d91-verify-log.txt` (8ece640)
- `workspaces/module-build/stages/04-verify/output/analysis-d91-F1.png` to `analysis-d91-F6.png` (8ece640)
- `workspaces/module-build/stages/04-verify/output/analysis-verification.md` (8ece640)
- `CHANGELOG.md` (this log)
- `workspaces/module-build/stages/05-integrate/output/analysis-integration-log.md` (this log)

### Audit

| Check | Result |
|-------|--------|
| Status matches evidence (no catalog row; no open D91 criterion except the stated one: D91e fails only against stale tracked P1 files, see above) | pass, with that note |
| Log names the commit SHAs and every file changed | pass |
| `CHANGELOG.md` has an entry for the amendment | pass |

## Amendment D93 (calibration of study p3-abstention), stage 05

Branch: `main` (commits pushed, no force-push). Decision D93 (stage 01; owner decisions P43 to P68 approved as proposed, P61 accepted as an audit exception).

### Catalog rows flipped

None, as for the base module.

Evidence: `stages/04-verify/output/analysis-verification.md` (section "Amendment D93: calibration of study p3-abstention (stage 04)"), accepted at the stage 04 checkpoint (2026-10-31). Criteria D93a to D93r pass. Real runs on the P1 and p2 stores (outputs in the scratchpad, not in the repository): O15 233 of 233 quantities equal to an independent Python oracle, O19 pass, byte-identical record and tables at 1, 4 and 23 threads, 105 to 222 s, at most 2.45 GB, data tree listing hash unchanged; O16 49 of 50 (finding V1). Full `Pkg.test()` after the review fixes 1401609 pass, 0 fail, 1 broken (existing); independent code review: 0 CRITICAL, 0 HIGH, MEDIUM and LOW fixed or recorded.

Updated in this stage: `CHANGELOG.md` (entry for the amendment), `src/Analysis/CONTEXT.md` (D93 section rewritten as prose: file list, invariants of `calibrate`, measured cost, limits V1 to V4, the D35 overrun; the file is back under the 80-line limit of `icm_validate.py`, which it exceeded since 076ae64).

### Open (not verified or deferred)

- V1: the plug-in HC test rejects at 0.0527 at AUC 0.5 and n_f = 30; recorded, no change (owner decision at stage 04).
- V2: PP for HB at most 0.360, targets at the cap (285 T_s, 1235 W2); for the owner at Amendment A1, before any test data.
- V3: no refusal test of the chunk-file guard. V4: verdicts on the real stores checked indirectly.
- No run on a clean committed tree yet: the A1 freeze run happens on committed code.
- Byte identity across machines (P58); the analysis folder suite is over the D35 ceiling.
- Not done here: Amendment A1, the Experiments amendment (panel builder, run blocks), the test-side analysis HA to HC; `_config/risks.md` review is due at the M4 exit, not at this amendment; CI for recent commits not checked.

### Commits

- `b2bd5cf` stage 01: D93 spec, decision log row, owner decisions
- `94bea67` stage 02: failing tests (11 files), plan, owner decisions
- `076ae64` stage 03: `calibrate`, 17 `p3_*` files, GREEN, recorded test rework, impl notes
- `71b11d4` stage 04: verification on the real stores, review fixes, oracle and O16 scripts, parameter table
- (this log) `analysis integration log (D93)`: the commit that adds this section and the documentation edits.

### Files changed on `main` for D93 (short hashes of the commits above that touched the file)

- `_config/decision-log.md` (b2bd5cf)
- `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (b2bd5cf)
- `test/analysis/an_p3_dispatch.jl`, `an_p3_fixture.jl`, `an_p3_params.jl`, `an_p3_pipeline.jl`, `an_p3_pp.jl`, `an_p3_ref.jl`, `an_p3_rule.jl`, `an_p3_select.jl`, `an_p3_sets.jl` (94bea67)
- `test/analysis/an_p3_helpers.jl`, `an_p3_refusals.jl`, `runtests.jl` (94bea67, 076ae64)
- `test/analysis/an_p3_gate.jl` (076ae64)
- `workspaces/module-build/stages/02-tests/output/analysis-test-plan.md` (94bea67, 076ae64)
- `src/Analysis/p3_c7.jl`, `p3_json.jl`, `p3_label.jl`, `p3_params.jl`, `p3_pp.jl`, `p3_select.jl`, `p3_stat.jl`, `p3_tables.jl`, `p3_verdict.jl`, `p3_hb.jl`, `p3_rule.jl` (076ae64)
- `src/Analysis/p3_hc.jl`, `p3_pipeline.jl`, `p3_read.jl`, `p3_record.jl`, `p3_rng.jl`, `p3_sets.jl` (076ae64, 71b11d4)
- `src/Analysis/Analysis.jl`, `cli.jl`, `p2_dispatch.jl` (076ae64)
- `src/Analysis/CONTEXT.md` (076ae64, this log)
- `bin/chaoslab.jl` (076ae64)
- `bench/analysis.jl` (076ae64)
- `verify/analysis-p3.jl` (076ae64)
- `verify/analysis-p3-export.jl`, `analysis-p3-oracle.py`, `analysis-p3-compare.py`, `analysis-p3-o16.jl` (71b11d4)
- `workspaces/module-build/stages/03-implement/output/analysis-d93-impl-notes.md` (076ae64)
- `workspaces/module-build/stages/04-verify/output/analysis-d93-verify-log.txt`, `analysis-d93-o16-log.txt`, `analysis-d93-suite-full-log.txt` (71b11d4)
- `workspaces/module-build/stages/04-verify/output/analysis-verification.md` (71b11d4)
- `workspaces/study/stages/01-design/references/p3-abstention-calibration-parameters.toml` (71b11d4)
- `CHANGELOG.md` (this log)
- `workspaces/module-build/stages/05-integrate/output/analysis-integration-log.md` (this log)

### Audit

| Check | Result |
|-------|--------|
| Status matches evidence (no catalog row; no open D93 criterion; findings V1 to V4 recorded) | pass |
| Log names the commit SHAs and every file changed | pass |
| `CHANGELOG.md` has an entry for the amendment | pass |

## Amendment D95 (test side of study p3-abstention), stage 05

Branch: `main` (commits pushed, no force-push). Decision D95 (stage 01; owner decisions P87 to P115 approved as drafted, P107 the audit exception of D93 P61 extended).

### Catalog rows flipped

None, as for the base module.

Evidence: `stages/04-verify/output/analysis-verification.md` (section "Amendment D95: test side of study p3-abstention (stage 04)" and its addendum), accepted at the stage 04 checkpoint (2026-11-03). Toy campaign of O22 (e): O20 independent Python oracle 645 of 645 quantities equal before and after the review fixes; O23 every altered log refused by `count-round` and `analyze` (1065 cases; the D94 finding V1 on the last entry is closed); D95a golden outputs of 263717c equal except the commit lines; byte-identical outputs at 1, 4 and 23 threads, `data/` listing unchanged; D95k HC power within 3 SE of an independent simulator at 8 sizes; `count-round` about 18 s, `analyze` about 29 s, at most 1.6 GB. Independent code review: 1 CRITICAL, 3 HIGH, 6 MEDIUM, 7 LOW, all fixed with 17 verify probes (owner decision "fix all"), one recorded test rework (three lines of the C1 fixture record). Full `Pkg.test()` after the fixes 1404598 pass, 0 fail, 1 broken (existing).

Updated in this stage: `analysis-spec.md` ("Clarifications at stage 05" of D95: the readings of the review fixes, the spec points found by the oracle author, O24, the recorded tg deviation), `src/Analysis/CONTEXT.md` (D95 cost and limits; 78 lines), `CHANGELOG.md` (entry), `workspaces/study/stages/03-run/references/round-loop.md` (new: the run and `count-round` loop and the O24 gate before the A1 freeze run) with its Inputs row in `workspaces/study/stages/03-run/CONTEXT.md`.

### Open (not verified or deferred)

- O24 runner acceptance on a runner-made store: a study gate after P69, before the A1 freeze run (owner decision 2026-11-03).
- No real test store yet (P115): arm names, item horizons and manifest fields rest on the pinned readings; time and memory at study size not measured.
- Toy coverage (V4): no false-regular rule, no tie at a set boundary, no dropped HA replicate.
- The M6 probe does not discriminate; three "before" probe failures (M9, L12, L17) are missing helpers, not reproduced wrong results.
- Byte identity across machines. Aqua and JET at stage 04 ran in a scratch environment; the full suite run includes them.
- Not done here: the study amendment P69, D94s, the A1 freeze run, A1; `_config/risks.md` review is due at the M4 exit.

### Commits

- `edb523e` stage 01: D95 spec, decision log row, owner decisions
- `88cfad0` stage 02: failing tests (10 files), plan
- `8a6adff` stage 03: `count-round` and the p3 test side of `analyze`, 12 `p3t_*` files, GREEN, recorded test reworks, impl notes
- `77770b8` stage 04: verification, review fixes, verify and oracle scripts, logs
- `2f9a36c`, `eb2af04` stage 04 addendum: full suite, D95k (V), O24 decision
- (this log) the commit that adds this section and the documentation edits

### Files changed on `main` for D95 (short hashes of the commits above that touched the file)

- `_config/decision-log.md` (edb523e)
- `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (edb523e, this log)
- `test/analysis/an_p3t_dispatch.jl`, `an_p3t_pipeline.jl`, `an_p3t_ref.jl`, `an_p3t_refusals.jl`, `an_p3t_stats.jl`, `runtests.jl` (88cfad0)
- `test/analysis/an_p3t_log.jl` (88cfad0, 8a6adff)
- `test/analysis/an_p3t_fixture.jl`, `an_p3t_helpers.jl`, `an_p3t_params.jl` (88cfad0, 77770b8)
- `test/analysis/an_p3_dispatch.jl` (8a6adff)
- `workspaces/module-build/stages/02-tests/output/analysis-test-plan.md` (88cfad0, 8a6adff)
- `src/Analysis/p3t_figs.jl`, `p3t_json.jl`, `p3t_stats.jl`, `Analysis.jl` (8a6adff)
- `src/Analysis/cli.jl`, `p2_dispatch.jl`, `p3t_count.jl`, `p3t_data.jl`, `p3t_inputs.jl`, `p3t_log.jl`, `p3t_params.jl`, `p3t_pipeline.jl`, `p3t_run.jl`, `p3t_store.jl`, `p3t_tables.jl` (8a6adff, 77770b8)
- `src/Analysis/CONTEXT.md` (8a6adff, this log)
- `bin/chaoslab.jl`, `bench/analysis.jl` (8a6adff)
- `workspaces/module-build/stages/03-implement/output/analysis-d95-impl-notes.md` (8a6adff)
- `verify/analysis-p3t.jl`, `analysis-p3t-alter.jl`, `analysis-p3t-cli.sh`, `analysis-p3t-run.py`, `analysis-p3t-runner.jl`, `analysis-p3t-power.jl`, `analysis-p3t-regress.sh`, `analysis-p3t-regress.py`, `analysis-p3t-export.jl`, `analysis-p3t-oracle.py`, `analysis-p3t-compare.py`, `analysis-d95-review.jl` (77770b8)
- `verify/analysis-d95k.jl`, `analysis-d95k.py` (2f9a36c)
- `workspaces/module-build/stages/04-verify/output/analysis-d95-verify-log.txt`, `analysis-d95-verify-log-fixed.txt`, `analysis-d95-oracle-log.txt`, `analysis-d95-oracle-log-fixed.txt`, `analysis-d95-review-before.txt`, `analysis-d95-review-after.txt` (77770b8)
- `workspaces/module-build/stages/04-verify/output/analysis-d95k-log.txt` (2f9a36c)
- `workspaces/module-build/stages/04-verify/output/analysis-verification.md` (77770b8, 2f9a36c, eb2af04)
- `CHANGELOG.md`, `workspaces/study/stages/03-run/CONTEXT.md`, `workspaces/study/stages/03-run/references/round-loop.md` (this log)
- `workspaces/module-build/stages/05-integrate/output/analysis-integration-log.md` (this log)

### Audit

| Check | Result |
|-------|--------|
| Status matches evidence (no catalog row; no open D95 criterion; O24 deferred to a study gate by owner decision, findings V2 to V4 recorded) | pass |
| Log names the commit SHAs and every file changed | pass |
| `CHANGELOG.md` has an entry for the amendment | pass |
