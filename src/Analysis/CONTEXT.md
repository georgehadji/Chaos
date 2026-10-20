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

## Limits and deferred

| Item | Note |
|------|------|
| K values of S (2 pi and 10) are a constant in `sets.jl` | spec corrections, item 7; a generic form needs a new key and stage 02 changes |
| Chunks that failed twice are refused, not analysed | spec corrections, item 3; the study has none |
| The reader keeps every item in memory and reads files serially | measured on the real store: 54 s, 2.45 GiB |
| Real-store results are in the tables file of study stage 04, not here | `workspaces/study/stages/04-analyze/` |
