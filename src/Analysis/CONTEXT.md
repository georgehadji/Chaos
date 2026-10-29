# Analysis

Study analysis (decision D60): reads a finished study run (manifest, run store, panel file, preregistration, parameter table) and computes the preregistered tables T1 to T7, the numbers behind the figures F1 to F7, the labels of D36 and D37 and the decision of H1 to H5 by the frozen rules. It reads only, owns no scientific threshold, and has no plotting dependency (the figures are a package extension). Layer between `Experiments` and `Figures` (`_config/architecture.md`).

Spec: `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (including "Corrections found in stage 04"). Tests: `test/analysis/`, plan `stages/02-tests/output/analysis-test-plan.md`. Report: `stages/04-verify/output/analysis-verification.md`. Command: `julia --project=verify bin/chaoslab.jl analyze <preregistration> <manifest> <panel> <parameters> [--store dir] [--out dir] [--threads N]`; exit code 0 done, 1 refusal or usage, 2 tables written but figures skipped or failed.

## Files

| File | Holds |
|------|-------|
| `Analysis.jl` | module, exports `analyze`, `read_params`, `read_store`, include order |
| `params.jl` | `read_params`: the parameter table E0 (37 keys, types and ranges checked) |
| `stats.jl` | `stream_seed`, `bootstrap_rng`, `bootstrap_interval`, `bootstrap_ratio`, `bootstrap_band`: percentile bootstrap, streams keyed by seed and label |
| `h1.jl` | `h1_orbit`, `decide_h1`, `control_deviations`, `check_controls` |
| `sets.jl` | `orbit_sets` (E_ch, W, S), `relative_spread`, `convergence_time`, `orbit_convergence`, `decide_h2`, `decide_h3` |
| `h4.jl` | `decide_h4`, `shevchenko_fit`, `chirikov_line` |
| `labels.jl` | `growth_ratio`, `base_label`, `orbit_label` (the four disagreement tests), `decide_h5` |
| `store.jl` | `read_store`, `find_item`: items found by content |
| `checks.jl` | `check_inputs`, `check_items`, `recheck_inputs`: every refusal, before any output |
| `orbits.jl` | `orbit_record`: all per-orbit quantities |
| `figdata.jl` | `figure_data`: the seven figure tables |
| `pipeline.jl` | `analyze`, the decisions in the order H1 to H5, the control-defect path, the figure hook `FIGURE_DRAWER` |
| `tables.jl` | the tables file and the figure digests |
| `cli.jl` | `cli` |

The plotting is `ext/ChaoslabCairoMakieExt.jl` (CairoMakie as a weak dependency, loaded by `bin/chaoslab.jl` under the verify environment); it sets `FIGURE_DRAWER` and draws `<study>-fig-F1.pdf` to `F7.pdf`.

## What it writes

`<out>/<study>-tables.md`: header (hashes of the four inputs, seed, module version, commit), T1 to T7, the figure data F1 to F7 each with a `digest:` line, the readings the module chose, the deviations it states, and the list of exploratory items (none). With figures, the seven PDF files. After a gated control exceeds its bound: only `<study>-defect.md` and status `:control_defect`.

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| The inputs have the same bytes after the run as before, also after a refusal | hashes taken at the start and checked again before the tables are written |
| Every refusal is an `ArgumentError` naming the cause and happens before anything is written | `checks.jl`; AC1, AC17 |
| The tables file is byte-identical for any thread count and output folder | bootstrap streams keyed by label; sorted iteration; AC16 |
| Decisions come in the order H1 to H5 and use only the numbers of the parameter table (and the constant K = 2 pi, 10 of S) | AC14 |
| A hypothesis whose eligible set has fewer than `min_set` orbits is "not decidable" with the count | `decide_h1` to `decide_h5` |
| Items of the tau sweep with `breakdown_tau` are results, not failures | `check_items`, F7 |

## Second study: p2-indicators (Amendment D89)

`analyze` reads the single header row "Study slug" of the preregistration (`p2_dispatch.jl`) and runs the pipeline of that study: `p1-le-benchmark` is `analyze_p1` (unchanged, `pipeline.jl`), `p2-indicators` is `analyze_p2` (`p2_pipeline.jl`); any other slug is refused. The command gains the options `--p1-store`, `--p1-preregistration`, `--p1-manifest`, `--p1-parameters` (the P1 reference inputs of the labels; a usage error for any other study). Exit code 2 means "figures skipped": the drawers of the six figures are a follow-up to the figures extension (P11). Spec: "Amendment D89" of the analysis spec; parameter table of 58 keys (E22); output `<out>/p2-indicators-tables.md`.

| File | Holds |
|------|-------|
| `p2_rules.jl` | `p2_budget_step`, `p2_nearest_step`, `p2_window_row`, `p2_growth`, `p2_ratio`, `p2_verdict` (R1, R2), `p2_score`, `p2_detection`, `p2_median` |
| `p2_agree.jl` | `p2_contingency`, `p2_unanimity`, `p2_spearman`, `p2_auc`, `p2_regular_range`, `p2_leaving_budget`, `p2_share_below` |
| `p2_stat.jl` | `p2_stat` (H1 to H5), `p2_h6_stats`, `p2_flips`, `p2_decide`, `p2_robust` |
| `p2_boot.jl` | `p2_label`, `p2_resample`, `p2_bootstrap` (joint resampling, streams keyed by label) |
| `p2_sets.jl` | `p2_sets` (E29), `p2_rho`, `p2_grid_label`, `p2_rotation_numbers`, `p2_lift` |
| `p2_store.jl` | `p2_read_store`, `p2_find_item`, `p2_flag`, `p2_grid_starts` |
| `p2_params.jl` | `p2_read_params`, `P2_PARAM_KEYS` (E22) |
| `p2_check.jl` | the input checks of E36 (`p2_check_inputs`, `p2_check_items`) |
| `p2_orbits.jl` | the items each orbit must have, the value records at the recorded steps, the orbit table with the reference class |
| `p2_hyp.jl` | verdict tables, the six decisions, the extension, the R2 version, the gate and the second versions |
| `p2_report.jl`, `p2_items.jl`, `p2_hitems.jl`, `p2_ritems.jl` | AUC, ranges, sensitivities, packages, flags, controls, cost; the items T1 to T17 and the data tables F1 to F6 |
| `p2_render.jl`, `p2_pipeline.jl`, `p2_dispatch.jl` | the tables file, `analyze_p2`, `analyze` and `study_slug` |

## Limits and deferred

| Item | Note |
|------|------|
| K values of S (2 pi and 10) are a constant in `sets.jl` | spec corrections, item 7; a generic form needs a new key and stage 02 changes |
| Chunks that failed twice are refused, not analysed | spec corrections, item 3; the study has none |
| The reader keeps every item in memory and reads files serially | measured on the real store: 54 s, 2.45 GiB |
| Real-store results are in the tables file of study stage 04, not here | `workspaces/study/stages/04-analyze/` |
