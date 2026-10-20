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
