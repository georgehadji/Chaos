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

`analyze` reads the single header row "Study slug" of the preregistration (`p2_dispatch.jl`) and runs the pipeline of that study: `p1-le-benchmark` is `analyze_p1` (unchanged, `pipeline.jl`), `p2-indicators` is `analyze_p2` (`p2_pipeline.jl`); any other slug is refused. The command gains the options `--p1-store`, `--p1-preregistration`, `--p1-manifest`, `--p1-parameters` (the P1 reference inputs of the labels; a usage error for any other study). Exit code 2 means "figures skipped" (the plotting package is not loaded) or "figures failed" (a drawer raised; the message names the figure and row). The result carries `status` (`:ok`, `:figures_skipped`, `:figures_failed`), `figures_error` and the `figures` dictionary F1 to F6 handed to the drawers (Amendment D91, `p2_pipeline.jl`, `pipeline.jl`: `_draw_figures`); the drawers are in `ext/p2_figures.jl` and `ext/p2_figures_b.jl`. Spec: "Amendment D89" of the analysis spec; parameter table of 58 keys (E22); output `<out>/p2-indicators-tables.md`.

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
| The reading P11 in the accepted P2 tables file still says that the figure drawers do not exist | owner decision: the file is accepted; changing `P2_READINGS` changes its bytes |

## Third study: p3-abstention (Amendment D93)

The command `calibrate <preregistration> <parameters> [options]` (`calibrate` in Julia; exit code 0 the record was written, also for "no feasible rule" and "sets too small", 1 refusal or usage) computes, from the stored p2 and P1 runs only, the calibration C0 to C10 and the power procedure PP of study `p3-abstention` and writes `<out>/p3-abstention-calibration-record.json` and `<out>/p3-abstention-calibration-tables.md` (table T1 and the data of figure F1 with its digest). `analyze` given the preregistration of this study runs the test side (Amendment D95, below). Spec: "Amendment D93"; parameter table of 71 keys (E47, `p3_params.jl`); nothing under `data/` is written.

Files: `p3_params.jl` (`P3_PARAM_KEYS`, `p3_rat`, `p3_read_params`, E47), `p3_rule.jl` (`p3_stable`, `p3_evidence`, `p3_rule_a`, `p3_rule_ab`, `p3_enumerate`, `p3_space_size`, `p3_rule_cost`, `p3_cost_cap`, `p3_matched`, `p3_d64_costs`, E51 to E53), `p3_sets.jl` (`p3_relabel`, E49 c; `p3_sets`, E50), `p3_select.jl` (`p3_select`, C3 to C8 in exact rationals; `p3_choose`, order of C6), `p3_rng.jl` (the streams of E55: SplitMix64, xoshiro256++, `p3_stream`, `p3_uniform!`, `p3_randint!`, `p3_normal_quantile`), `p3_stat.jl` (exact McNemar tail and critical values, rounding, Hanley and McNeil variance, AUC, binormal shift, votes of HC), `p3_hb.jl`, `p3_pp.jl`, `p3_hc.jl` (PP for HB and HC, E56, E57), `p3_read.jl`, `p3_label.jl`, `p3_verdict.jl`, `p3_c7.jl` (inputs and checks E48, C0 and the change counts E49, verdicts and R_ab, the items C7), `p3_json.jl`, `p3_record.jl`, `p3_tables.jl`, `p3_pipeline.jl` (record, tables file, `calibrate`). Test side (Amendment D95): `count-round <preregistration> <parameters> <test-parameters> [options]` appends one entry to the round log of the stopping rule after verifying the whole log; `analyze <preregistration> <parameters> <test-parameters> [options]` decides HA, HB and HC on the closed log and writes `p3-abstention-test-results.json` and `p3-abstention-tables.md`. Files: `p3t_stats.jl` (exact percentile, ROC, resample, strata, count rule, `p3t_ha`, `p3t_hb_compare`, `p3t_hb`, `p3t_hc`), `p3t_params.jl` (`P3T_PARAM_KEYS`, `p3t_read_params`, E60), `p3t_store.jl` (chunk scan, ref state, counts), `p3t_log.jl` (entry format, verification V1), `p3t_inputs.jl` (checks E61), `p3t_count.jl` (`count_round`), `p3t_data.jl` (items, drops, sets, verdicts, R_ab), `p3t_run.jl` (HB tables, sensitivities, power at achieved sizes), `p3t_figs.jl` (F2, F3), `p3t_json.jl`, `p3t_tables.jl`, `p3t_pipeline.jl` (`analyze_p3t`). Reads datasets of the ref arm only before closure (I39); `count-round` writes only the appended entry (I38), through a temporary file and a rename. Cost on the toy campaign: `count-round` about 18 s, `analyze` about 29 s, at most 1.6 GB. Since Amendment D97 the test side reads the store of the runner (cells `k<i>:<class>:r<round>`, six arms, tg at N_c, failed chunks; verified on a fixture campaign of the runner, O28, O29). Since Amendment D99 `calibrate` without an output folder writes to `stages/04-analyze/output`, the folder `count-round` and `analyze` read the record from; each T5 row of HA item 3 shows the status of its own condition, and the decision line the status of the item. Limits: no real-size test store yet; runner acceptance of the log (O24) on the A2 block is a study gate after P69; the toy does not reach false-regular rules, set-boundary ties or dropped replicates (D95 stage 04 report, V2 to V4).

Invariants: the nine input files and every chunk file listed in the two manifests are hashed before any read and checked again before writing (a change refuses); record and tables file are byte-identical for any thread count and output folder; changes under `src`, `ext`, `bin`, `Project.toml`, `Manifest.toml` (untracked files included) flag the tables file "not a freeze run", no refusal (P59). Cost on the real stores: 105 to 222 s, at most 2.45 GB at 1 to 23 threads, dominated by reading and hashing. Limits: findings V1 to V4 of the stage 04 report (no refusal test of the chunk guard, verdicts not in the record, HC test slightly anti-conservative at AUC 0.5, HB underpowered on the stored data); no byte identity across machines for PP-HC (P58); the analysis folder suite (6 min) is over the D35 ceiling, 24 of 54 refusal cases run only in `verify/analysis-p3.jl` (D38).
